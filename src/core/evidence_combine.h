#ifndef MISPITOOLS_CORE_EVIDENCE_COMBINE_H
#define MISPITOOLS_CORE_EVIDENCE_COMBINE_H

#include <cstdint>
#include <vector>

#include "lr_dist.h"
#include "result.h"

namespace mispitools {
namespace core {

/// @brief Evidence-combination mode (DESIGN.md §8.7).
///
/// The numeric values are part of the binding contract (R / JS pass the
/// mode as an int).
enum class CombineMode : std::uint8_t {
    /// Conditional independence given the hypothesis: the total log10 LR is
    /// the sum of the per-feature log10 LRs, so the combined distribution is
    /// the convolution of the per-feature ones. Identical to
    /// `lr_dist_compose` (Exact) — DESIGN.md §8.5 / §8.7.
    Independent = 0,
    /// Egeland-Marsico (2026) chain. The numerator (H1) keeps conditional
    /// independence (Eq. 2 + binary error model); the denominator (H2 /
    /// population) carries a first-order Markov dependence along the feature
    /// ordering, supplied as per-adjacent-pair transition matrices.
    MarkovSE    = 1
};

/// @brief Parameters of the Markov-SE chain (ignored in `Independent`).
///
/// `transition` has exactly `K - 1` entries for `K` features (empty when
/// `K == 1`). `transition[s]` is the row-stochastic transition matrix on
/// the **H2 channel** from feature `s` to feature `s+1`:
///
///   `transition[s][i * n_{s+1} + j] = P(feature s+1 = atom j |
///                                        feature s = atom i, H2)`
///
/// stored row-major, with `n_s = per_feature[s].size()` rows and
/// `n_{s+1} = per_feature[s+1].size()` columns. Each row must be
/// non-negative and sum to 1 (tolerance `1e-9`).
///
/// Setting every row of `transition[s]` equal to the marginal
/// `per_feature[s+1].p_h2` (normalised) makes the H2 joint factorise, so
/// `MarkovSE` then reproduces `Independent` (the independence-equivalence
/// invariant, checked numerically in F6.5 / cross-checked in F6.7 against
/// the Egeland-Marsico reproduction script).
struct MarkovSEParams {
    std::vector<std::vector<double>> transition;  ///< K-1 flattened matrices
};

/// @brief Combine per-feature LR distributions into the case-level one.
///
/// `per_feature` is the vector of per-feature `LrDist` objects (genetic
/// markers and/or non-genetic features — the schema is shared, F6.4). Each
/// atom carries the paired `(p_h1, p_h2)`; `log10_lr` is the per-feature
/// ratio. The combined `LrDist` keeps the genetic-engine boundary
/// convention (`p_h1 > 0, p_h2 = 0` → +Inf; `p_h1 = 0, p_h2 > 0` → -Inf;
/// `0 / 0` dropped) and the same sorted/aggregated layout as
/// `per_marker_lr_dist`.
///
/// * `Independent`: delegates to `lr_dist_compose` (Exact, default
///   options) — bit-for-bit identical to §8.5. `params` is ignored.
/// * `MarkovSE`: enumerates the cartesian product of the per-feature
///   atoms. For a joint outcome `(x_0, …, x_{K-1})`:
///     - H1 mass `= Πₛ p_h1[s][x_s]` (independence under H1);
///     - H2 mass `= p_h2[0][x_0] · Π_{s≥1} transition[s-1][x_{s-1}][x_s]`;
///     - combined `log10_lr = log10(H1) - log10(H2)`.
///   Atoms are then sorted ascending and equal-`log10_lr` atoms collapsed
///   (probabilities summed), exactly as `per_marker_lr_dist`. The
///   per-feature `log10_lr` field is unused in this mode (the combined
///   ratio is recomputed from the joint masses because the H2 channel no
///   longer factorises). Σp_h1 and Σp_h2 are preserved (= 1 when the
///   inputs are normalised).
///
/// Errors (mismatched columns, wrong `transition` shape, non-stochastic
/// rows, an empty `per_feature` in `MarkovSE`, or a path-count blow-up)
/// flow back as `Result::error`; the kernel never throws.
///
/// @complexity Independent: see `lr_dist_compose`. MarkovSE:
/// `O(Πₛ nₛ)` — exponential in the number of features, but the
/// per-feature supports of supplementary evidence are tiny (sex = 2,
/// pigmentation ≤ 3, age bins, …); the fold is capped at
/// `MISPI_MARKOVSE_MAX_PATHS` joint paths.
Result<LrDist> evidence_combine(const std::vector<LrDist>& per_feature,
                                CombineMode mode,
                                const MarkovSEParams& params = {});

// Placeholder retained for the F0.5 cpp-bootstrap regression test.
int evidence_combine_placeholder(int x);

}  // namespace core
}  // namespace mispitools

#endif  // MISPITOOLS_CORE_EVIDENCE_COMBINE_H
