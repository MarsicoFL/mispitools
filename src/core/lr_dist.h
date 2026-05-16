#ifndef MISPITOOLS_CORE_LR_DIST_H
#define MISPITOOLS_CORE_LR_DIST_H

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

// Placeholder retained for the F0.5 cpp-bootstrap regression test.
int lr_dist_placeholder(int x);

}  // namespace core
}  // namespace mispitools

#endif  // MISPITOOLS_CORE_LR_DIST_H
