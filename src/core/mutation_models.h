#ifndef MISPITOOLS_CORE_MUTATION_MODELS_H
#define MISPITOOLS_CORE_MUTATION_MODELS_H

#include <cstdint>
#include <vector>

#include "marker.h"
#include "result.h"

namespace mispitools {
namespace core {

enum class MutationKind : std::uint8_t {
    None         = 0,
    Equal        = 1,
    Stepwise     = 2,
    Proportional = 3,
    Asymmetric   = 4
};

/// @brief Parameter bundle for a single-marker mutation model.
struct MutationModel {
    MutationKind kind = MutationKind::None;
    double rate  = 0.0;   // R in [0, 1)
    double range = 0.0;   // r for stepwise/asymmetric in (0, 1)
    double rate2 = 0.0;   // out-of-microgroup rate; default 0
    double bias  = 0.5;   // u for asymmetric, in [0, 1]
};

/// @brief Build the K x K identity mutation matrix (None model).
///   M[i,i] = 1, M[i,j] = 0 for j != i.
/// Requires K >= 1.
/// @param K Number of alleles.
Result<std::vector<double>> mutation_matrix_none(AlleleIndex K);

/// @brief Build the K x K mutation matrix for the Equal-rate model.
///   M[i,i] = 1 - R, M[i,j] = R / (K - 1) for j != i.
/// Requires K >= 2 and rate in [0, 1].
/// @param K Number of alleles.
/// @param rate Per-meiosis mutation rate R.
Result<std::vector<double>> mutation_matrix_equal(
    AlleleIndex K, double rate);

/// @brief Build the K x K mutation matrix for the Stepwise model.
/// For each parental allele i:
///   w_j = range ^ |s_j - s_i|  for j != i, with w_i = 0.
///   sw  = sum_j w_j.
///   M[i,j] = (rate / sw) * w_j,  M[i,i] = 1 - rate.
/// Requires K >= 2, rate in [0, 1], `numeric_labels.size() == K`, every
/// label finite, and a strictly positive off-diagonal weight sum on each
/// row (the latter prevents division by zero when range == 0).
/// @param K Number of alleles.
/// @param rate Per-meiosis mutation rate R.
/// @param range Geometric step ratio r in (0, 1].
/// @param numeric_labels Allele labels as doubles (length K).
Result<std::vector<double>> mutation_matrix_stepwise(
    AlleleIndex K, double rate, double range,
    const std::vector<double>& numeric_labels);

/// @brief Dispatch builder. Delegates to the model-specific builder above
/// according to `mut.kind`. Asymmetric is not yet implemented (F5.1).
/// @param mut Parameter bundle.
/// @param n_alleles K.
/// @param numeric_labels Required for Stepwise/Asymmetric (size K, finite
/// entries). Ignored for None / Equal / Proportional.
Result<std::vector<double>> build_mutation_matrix(
    const MutationModel& mut,
    AlleleIndex n_alleles,
    const std::vector<double>& numeric_labels);

// Placeholder retained for the F0.5 cpp-bootstrap regression test.
int mutation_models_placeholder(int x);

}  // namespace core
}  // namespace mispitools

#endif  // MISPITOOLS_CORE_MUTATION_MODELS_H
