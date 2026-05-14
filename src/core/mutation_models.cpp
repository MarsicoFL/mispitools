#include "mutation_models.h"

#include <cmath>
#include <cstddef>
#include <string>
#include <utility>

namespace mispitools {
namespace core {

namespace {

inline bool rate_in_unit_interval(double r) noexcept {
    return std::isfinite(r) && r >= 0.0 && r <= 1.0;
}

}  // namespace

Result<std::vector<double>> mutation_matrix_none(AlleleIndex K) {
    if (K <= 0) {
        return err_result<std::vector<double>>(
            "mutation_matrix_none: n_alleles must be positive.");
    }
    std::vector<double> M(static_cast<std::size_t>(K) * K, 0.0);
    for (AlleleIndex i = 0; i < K; ++i) {
        M[static_cast<std::size_t>(i) * K + i] = 1.0;
    }
    return ok_result(std::move(M));
}

Result<std::vector<double>> mutation_matrix_equal(AlleleIndex K, double rate) {
    if (K < 2) {
        return err_result<std::vector<double>>(
            "mutation_matrix_equal: equal-rate mutation requires K >= 2.");
    }
    if (!rate_in_unit_interval(rate)) {
        return err_result<std::vector<double>>(
            "mutation_matrix_equal: rate must be finite in [0, 1].");
    }
    const double off = rate / static_cast<double>(K - 1);
    const double on  = 1.0 - rate;
    std::vector<double> M(static_cast<std::size_t>(K) * K, off);
    for (AlleleIndex i = 0; i < K; ++i) {
        M[static_cast<std::size_t>(i) * K + i] = on;
    }
    return ok_result(std::move(M));
}

Result<std::vector<double>> mutation_matrix_stepwise(
        AlleleIndex K, double rate, double range,
        const std::vector<double>& numeric_labels) {
    if (K < 2) {
        return err_result<std::vector<double>>(
            "mutation_matrix_stepwise: stepwise mutation requires K >= 2.");
    }
    if (!rate_in_unit_interval(rate)) {
        return err_result<std::vector<double>>(
            "mutation_matrix_stepwise: rate must be finite in [0, 1].");
    }
    if (!std::isfinite(range) || range < 0.0) {
        return err_result<std::vector<double>>(
            "mutation_matrix_stepwise: range must be finite and >= 0.");
    }
    if (static_cast<AlleleIndex>(numeric_labels.size()) != K) {
        return err_result<std::vector<double>>(
            "mutation_matrix_stepwise: numeric_labels.size() must equal K.");
    }
    for (AlleleIndex i = 0; i < K; ++i) {
        if (!std::isfinite(numeric_labels[static_cast<std::size_t>(i)])) {
            return err_result<std::vector<double>>(
                "mutation_matrix_stepwise: every numeric_label must be "
                "finite (non-numeric allele labels are not supported).");
        }
    }

    std::vector<double> M(static_cast<std::size_t>(K) * K, 0.0);
    std::vector<double> w(static_cast<std::size_t>(K), 0.0);
    for (AlleleIndex i = 0; i < K; ++i) {
        const double s_i = numeric_labels[static_cast<std::size_t>(i)];
        double sw = 0.0;
        for (AlleleIndex j = 0; j < K; ++j) {
            if (j == i) {
                w[static_cast<std::size_t>(j)] = 0.0;
                continue;
            }
            const double s_j = numeric_labels[static_cast<std::size_t>(j)];
            const double steps = std::fabs(s_j - s_i);
            const double wj = std::pow(range, steps);
            w[static_cast<std::size_t>(j)] = wj;
            sw += wj;
        }
        if (!std::isfinite(sw) || sw <= 0.0) {
            return err_result<std::vector<double>>(
                "mutation_matrix_stepwise: row weights sum to zero for "
                "allele index " + std::to_string(static_cast<int>(i)) +
                "; cannot normalize (range == 0 with K >= 2 triggers this).");
        }
        const double scale = rate / sw;
        for (AlleleIndex j = 0; j < K; ++j) {
            M[static_cast<std::size_t>(i) * K + j] =
                scale * w[static_cast<std::size_t>(j)];
        }
        M[static_cast<std::size_t>(i) * K + i] = 1.0 - rate;
    }
    return ok_result(std::move(M));
}

Result<std::vector<double>> build_mutation_matrix(
        const MutationModel& mut,
        AlleleIndex K,
        const std::vector<double>& numeric_labels) {
    switch (mut.kind) {
        case MutationKind::None:
            return mutation_matrix_none(K);
        case MutationKind::Equal:
            return mutation_matrix_equal(K, mut.rate);
        case MutationKind::Stepwise:
            return mutation_matrix_stepwise(K, mut.rate, mut.range,
                                            numeric_labels);
        case MutationKind::Proportional:
            return err_result<std::vector<double>>(
                "build_mutation_matrix: Proportional model is not "
                "implemented (no current hito assigns it).");
        case MutationKind::Asymmetric:
            return err_result<std::vector<double>>(
                "build_mutation_matrix: Asymmetric model arrives in F5.1.");
    }
    return err_result<std::vector<double>>(
        "build_mutation_matrix: unknown MutationKind.");
}

int mutation_models_placeholder(int x) {
    return x + 1;
}

}  // namespace core
}  // namespace mispitools
