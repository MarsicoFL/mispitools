#ifndef MISPITOOLS_CORE_LR_DIST_H
#define MISPITOOLS_CORE_LR_DIST_H

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

#include "cpt_engine.h"
#include "result.h"

namespace mispitools {
namespace core {

/// @brief Sparse per-marker LR distribution: (log10 LR, P(.|H1), P(.|H2)).
///
/// Mirrors the output of `R/r_ref_per_marker.R::per_marker_lr_dist_R`.
/// When aggregated, `log10_lr` is sorted ascending with duplicate atoms
/// collapsed (probabilities summed), so two distributions over the same
/// support align by a parallel scan — the layout `lr_dist_compose` (F4.2)
/// folds over.
///
/// Boundary convention (matches R-ref `log10_lr_from_probs` + DESIGN §8.4):
///   * P1 > 0, P2 > 0 → log10_lr = log10(P1) - log10(P2).
///   * P1 > 0, P2 = 0 → log10_lr = +Inf; `has_pos_inf = true`.
///   * P1 = 0, P2 > 0 → log10_lr = -Inf; `has_neg_inf = true`.
///   * P1 = 0, P2 = 0 → skipped (the limit `0 log 0 = 0`; `cpt_marker_joint`
///     already filters these upstream).
/// The (log10_lr, p_h1, p_h2) triple keeps the raw `p_h1` / `p_h2` even on
/// an infinite atom (one of them is exactly 0 there), matching R-ref.
struct LrDist {
    std::vector<double> log10_lr;   ///< sorted ascending when aggregated
    std::vector<double> p_h1;       ///< same length as log10_lr
    std::vector<double> p_h2;       ///< same length as log10_lr
    bool has_pos_inf = false;       ///< a (+Inf) atom: H2 lacks an H1 state
    bool has_neg_inf = false;       ///< a (-Inf) atom: H1 lacks an H2 state

    std::size_t size() const noexcept { return log10_lr.size(); }
    bool empty() const noexcept { return log10_lr.empty(); }
};

/// @brief Relative floor below which two atom keys are the same atom.
///
/// A log10 LR atom is a real number reached by two different arithmetic
/// routes: the C++ engine and the R reference compute the same quantity
/// with different association orders, and a compiler that contracts
/// `a * b + c` into a fused multiply-add (the default on aarch64, and on
/// any target where the FMA instruction is in the baseline ISA) rounds
/// once where another rounds twice. Mathematically equal atoms then differ
/// in the last bits, and grouping by IEEE equality splits them into
/// separate buckets, which makes the size of the support platform
/// dependent.
///
/// The separation is measured, not assumed: over the composed profiles the
/// test suite exercises, consecutive keys either sit within one unit in the
/// last place of each other (the same atom reached twice) or are more than
/// 1e-8 apart in relative terms (genuinely different atoms). The band
/// between those two populations is empty across seven orders of magnitude,
/// and `1e-12` sits in the middle of it: four orders above the rounding
/// noise it has to absorb, four orders below the closest genuine pair it
/// must keep apart.
constexpr double kAtomRelTol = 1e-12;

/// @brief True when key `k` belongs to the group opened by `key`.
///
/// `k` is assumed to come from an ascending sort, so `k >= key`. IEEE
/// equality is tested first, which folds each infinite class into a single
/// bucket (`Inf == Inf`) and leaves `Inf - Inf` (NaN) out of the arithmetic
/// path. Finite keys additionally merge when their gap falls under
/// `extra_tol` (a caller-supplied absolute tolerance, `LrDistComposeOptions
///::merge_tol`) or under the `kAtomRelTol` floor scaled by the magnitude of
/// the key.
inline bool same_atom(double k, double key, double extra_tol = 0.0) {
    if (!(k != key)) return true;
    if (!std::isfinite(k) || !std::isfinite(key)) return false;
    const double tol =
        std::max(extra_tol, kAtomRelTol * std::max(1.0, std::fabs(key)));
    return (k - key) <= tol;
}

/// @brief Per-marker LR distribution from a sparse joint table.
///
/// @param joint joint genotype distribution (output of `cpt_marker_joint`);
///        `p_h1` and `p_h2` are aligned row-by-row.
/// @param aggregate when true (default) the atoms are sorted ascending by
///        `log10_lr` and equal-`log10_lr` atoms are merged (probabilities
///        summed), bit-for-bit with R-ref `aggregate_lr_dist`. When false
///        the atoms keep the joint's row order.
/// @return `LrDist` on the non-null support; errors on length mismatch or
///         a negative probability.
/// @complexity O(n_rows) raw; O(n_rows log n_rows) aggregated (one stable
/// sort + linear collapse).
Result<LrDist> per_marker_lr_dist(const JointTable& joint,
                                  bool aggregate = true);

/// @brief Composition method for `lr_dist_compose` (F4.2).
enum class ComposeMethod : std::uint8_t {
    Exact = 0,   ///< exact sparse convolution (sum-of-keys)
    Grid  = 1     ///< fixed log-LR lattice heuristic (controlled error)
};

/// @brief Options for `lr_dist_compose`.
struct LrDistComposeOptions {
    ComposeMethod method = ComposeMethod::Exact;
    /// Exact mode: collapse finite keys whose ascending gap is `<= merge_tol`
    /// (representative = the group's first key, matching R-ref `aggregate`).
    /// Default 0.0 → grouping by the `kAtomRelTol` floor alone, which merges
    /// rounding-noise duplicates and nothing else, and agrees with the R
    /// reference on every platform. Raise it to bucket genuinely distinct
    /// but near-equal atoms (caps support growth); it is applied as a floor,
    /// so a value under `kAtomRelTol` has no effect.
    double merge_tol = 0.0;
    /// Grid mode: number of lattice points spanning the finite total range
    /// `[Σ min, Σ max]`. Discretisation error is O(span / (grid_points-1)).
    /// Ignored in exact mode. Must be >= 2.
    int grid_points = 512;
};

/// @brief Compose independent per-feature LR distributions.
///
/// The total log10 LR under conditional independence is the sum of the
/// per-feature log10 LRs (ROADMAP §Combinación, mode `independent`), so the
/// composed distribution is the convolution of the per-feature ones.
///
/// `Exact`: sequential fold (`acc = δ₀; acc = conv(acc, dᵢ)`) over a
/// sparse `{log10_lr, p_h1, p_h2}` support — `p_h1` and `p_h2` convolved on
/// the *same* key grid in one pass. Both probabilities are multiplied; a
/// combined atom with both probabilities zero is dropped (the `0 log 0 = 0`
/// limit; this also disposes of the `+Inf` × `-Inf` cross term, whose mass
/// is identically zero). ±Inf atoms (mutation=none) propagate via IEEE
/// arithmetic and aggregate into single buckets, mirroring
/// `per_marker_lr_dist`. Pattern reference: `DNAtools::convolve` (GPL-2+) —
/// re-implemented, not lifted; the integer-offset `+1` shift of that code
/// is structurally inexpressible here (keys are real-valued).
///
/// `Grid`: project every feature onto a common lattice of spacing
/// `delta = span / (grid_points-1)` via mass- and mean-preserving linear
/// splitting, then convolve on the integer lattice (exact index addition,
/// no clamping). The composed mean is preserved exactly; the discretisation
/// error on the shape is O(delta). Requires all-finite supports (±Inf atoms
/// → error; use `Exact`).
///
/// Invariants (SCOUT_DNAtools F4.2): Σp = 1 after each fold;
/// `compose([d]) == d` (idempotence); `compose([d1,d2]) == compose([d2,d1])`
/// (commutativity). An empty input list composes to the identity δ₀
/// (`{0, 1, 1}`).
///
/// @complexity Exact: O(K · |support|²) sparse folds. Grid:
/// O(K · grid_points · |feature|).
Result<LrDist> lr_dist_compose(const std::vector<LrDist>& per_feature,
                               const LrDistComposeOptions& opts = {});

// Placeholder retained for the F0.5 cpp-bootstrap regression test.
int lr_dist_placeholder(int x);

}  // namespace core
}  // namespace mispitools

#endif  // MISPITOOLS_CORE_LR_DIST_H
