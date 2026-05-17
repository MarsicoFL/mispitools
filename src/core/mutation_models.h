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

/// @brief Maximum well-defined `rate` for the Dawid (asymmetric) model
/// given `afreq` and `range`. Beyond this cap the row-stochastic diagonal
/// turns negative and the model is undefined. This is the closed-form
///   1 / max_i sum_{j != i} (1 / afreq[i]) * C * range^|i-j|,
/// with C = (1 - range) / (2 * range * (n - a)), a = (1 - range^n)/(1-range).
/// Mirrors `pedmut::maxRate()` (the `UW` bound). Returns an error if
/// `afreq` is empty, `range` is not in (0, 1), or any frequency is
/// non-finite / non-positive.
/// @param afreq Allele frequencies (length n >= 2, positive, finite).
/// @param range Geometric step ratio in (0, 1).
/// Complexity: O(n^2).
Result<double> dawid_max_rate(
    const std::vector<double>& afreq, double range);

/// @brief Build the K x K mutation matrix for the asymmetric Dawid (2002)
/// model, as implemented by `pedmut::mutationMatrix(model = "dawid")`.
/// With a = (1 - range^K)/(1 - range) and c = rate*(1-range)/(2*range*(K-a)):
///   M[i,j] = c / afreq[i] * range^|i-j|        for j != i,
///   M[i,i] = 1 - sum_{j != i} M[i,j].
/// The step distance uses the allele rank position |i-j| (Dawid is
/// rank-based, not label-based); the matrix is reversible with respect to
/// `afreq`. Requires K >= 2, rate finite in [0, 1], range finite in
/// (0, 1), `afreq.size() == K` with every entry finite and > 0, and a
/// rate not exceeding `dawid_max_rate()` (otherwise the model is
/// undefined and an error carrying the cap is returned).
/// @param K Number of alleles.
/// @param rate Per-meiosis mutation rate R.
/// @param range Geometric step ratio r in (0, 1).
/// @param afreq Allele frequencies (length K, positive, finite).
/// Complexity: O(K^2).
Result<std::vector<double>> mutation_matrix_asymmetric(
    AlleleIndex K, double rate, double range,
    const std::vector<double>& afreq);

/// @brief Dispatch builder. Delegates to the model-specific builder above
/// according to `mut.kind`.
/// @param mut Parameter bundle.
/// @param n_alleles K.
/// @param numeric_labels Required for Stepwise (size K, finite entries).
/// Ignored for None / Equal / Asymmetric.
/// @param afreq Allele frequencies; required for Asymmetric (Dawid),
/// ignored otherwise. Defaults to empty.
Result<std::vector<double>> build_mutation_matrix(
    const MutationModel& mut,
    AlleleIndex n_alleles,
    const std::vector<double>& numeric_labels,
    const std::vector<double>& afreq = {});

// Placeholder retained for the F0.5 cpp-bootstrap regression test.
int mutation_models_placeholder(int x);

}  // namespace core
}  // namespace mispitools

#endif  // MISPITOOLS_CORE_MUTATION_MODELS_H
