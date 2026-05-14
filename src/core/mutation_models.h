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

/// @brief Build the K x K mutation matrix (row-major: M[i,j] = mat[i*K + j])
/// for the given model. Row i is the conditional distribution of the
/// transmitted allele given the parental allele was i. Rows sum to 1.
/// @param mut Parameter bundle.
/// @param n_alleles K.
/// @param numeric_labels Required for Stepwise/Asymmetric (size K, no NaN
/// entries). Ignored for None / Equal / Proportional.
///
/// In F2.2 only MutationKind::None is implemented. Other kinds return
/// Result::error so the binding layer can surface a clear message. The
/// matrices for Equal / Stepwise arrive in F2.3, Asymmetric in F5.1.
Result<std::vector<double>> build_mutation_matrix(
    const MutationModel& mut,
    AlleleIndex n_alleles,
    const std::vector<double>& numeric_labels);

// Placeholder retained for the F0.5 cpp-bootstrap regression test.
int mutation_models_placeholder(int x);

}  // namespace core
}  // namespace mispitools

#endif  // MISPITOOLS_CORE_MUTATION_MODELS_H
