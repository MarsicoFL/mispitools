#include "mutation_models.h"

#include <cstddef>
#include <utility>

namespace mispitools {
namespace core {

Result<std::vector<double>> build_mutation_matrix(
        const MutationModel& mut,
        AlleleIndex K,
        const std::vector<double>& /*numeric_labels*/) {
    if (K <= 0) {
        return err_result<std::vector<double>>(
            "build_mutation_matrix: n_alleles must be positive.");
    }
    std::vector<double> M(static_cast<std::size_t>(K) * K, 0.0);
    if (mut.kind == MutationKind::None) {
        for (AlleleIndex i = 0; i < K; ++i) {
            M[static_cast<std::size_t>(i) * K + i] = 1.0;
        }
        return ok_result(std::move(M));
    }
    return err_result<std::vector<double>>(
        "build_mutation_matrix: only MutationKind::None is implemented "
        "in F2.2. Equal/Stepwise arrive in F2.3, Asymmetric in F5.1.");
}

int mutation_models_placeholder(int x) {
    return x + 1;
}

}  // namespace core
}  // namespace mispitools
