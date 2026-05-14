#ifndef MISPITOOLS_CORE_KL_ENGINE_H
#define MISPITOOLS_CORE_KL_ENGINE_H

#include <cstdint>

#include "cpt_engine.h"
#include "result.h"

namespace mispitools {
namespace core {

/// @brief Bidirectional KL divergence + expected log10 LR for one marker.
///
/// Mirrors `R/r_ref_per_marker.R::per_marker_kl_R`. Sparse-aware: consumes
/// the joint produced by `cpt_marker_joint` as already filtered (rows with
/// both `p_h1` and `p_h2` zero are dropped upstream).
///
/// Boundary convention (matches R-ref + DESIGN.md §6.2 + SCOUT F2.7):
///   * P1 > 0, P2 > 0 → standard contribution.
///   * P1 > 0, P2 = 0 → log10 LR = +Inf → e_log10_lr_h1 = +Inf,
///                      kl_h1_to_h2 = +Inf, abs_cont_violations_h2 += 1.
///   * P1 = 0, P2 > 0 → log10 LR = -Inf → e_log10_lr_h2 = -Inf,
///                      kl_h2_to_h1 = +Inf, abs_cont_violations_h1 += 1.
///   * P1 = 0, P2 = 0 → skipped (limit `0 log 0 = 0`).
///
/// The `abs_cont_violations_*` and `mass_violations_*` counters expose the
/// KLde-style decomposition (forensIT::KLde) so the binding layer can tell
/// the user *why* a `+Inf` came back. They are diagnostics — finite KL
/// values are unaffected.
struct PerMarkerKL {
    double e_log10_lr_h1 = 0.0;                 ///< E[log10 LR | H1]
    double e_log10_lr_h2 = 0.0;                 ///< E[log10 LR | H2]
    double kl_h1_to_h2 = 0.0;                   ///< KL(P_H1 || P_H2) in nats
    double kl_h2_to_h1 = 0.0;                   ///< KL(P_H2 || P_H1) in nats
    std::int32_t abs_cont_violations_h2 = 0;    ///< rows with P_H1>0, P_H2=0
    std::int32_t abs_cont_violations_h1 = 0;    ///< rows with P_H2>0, P_H1=0
    double mass_violations_h2 = 0.0;            ///< sum P_H1 over abs_cont_h2 rows
    double mass_violations_h1 = 0.0;            ///< sum P_H2 over abs_cont_h1 rows
};

/// @brief Bidirectional KL + expected log10 LR from a sparse joint table.
///
/// @param joint joint genotype distribution (output of `cpt_marker_joint`);
///        `p_h1` and `p_h2` are aligned row-by-row.
/// @return `PerMarkerKL` populated as described in the struct doc; errors
///         are returned when `p_h1` and `p_h2` differ in length or when a
///         negative probability is found.
/// @complexity O(n_rows). One linear pass over the sorted joint.
Result<PerMarkerKL> per_marker_kl(const JointTable& joint);

// Placeholder retained for the F0.5 cpp-bootstrap regression test.
int kl_engine_placeholder(int x);

}  // namespace core
}  // namespace mispitools

#endif  // MISPITOOLS_CORE_KL_ENGINE_H
