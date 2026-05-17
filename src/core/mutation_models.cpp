#include "mutation_models.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
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

namespace {

// Format a double with up to 15 significant digits (no trailing-zero
// padding), matching how R prints the undefined-model cap.
std::string format_g(double x) {
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.15g", x);
    return std::string(buf);
}

// --- Dawid (2002) asymmetric / reversible mutation model -------------
// Adapted from pedmut v0.9.0 (R/mutationMatrix.R, .dawid / maxRate),
// GPL-3, by Magnus D. Vigeland (NMBU). Reference: Dawid A.P. (2002),
// "Properties of diagnostic data distributions", Biometrics 58.
// See package/mispitools/COPYRIGHTS.

// Shared off-diagonal kernel for the Dawid model. Fills `line_sum[i]`
// with sum_{j != i} (1 / afreq[i]) * C * range^|i-j|, i.e. row i of
// (R / rate) with a zero diagonal, where C = (1 - range) /
// (2 * range * (n - a)) and a = (1 - range^n)/(1 - range). The maximum
// of these line sums determines both the undefined-model cap and the
// off-diagonal entries (scaled by `rate`). Callers must have validated
// K >= 2, range in (0, 1), afreq.size() == K and afreq positive/finite.
void dawid_line_sums(AlleleIndex K, double range,
                     const std::vector<double>& afreq,
                     std::vector<double>& line_sum) {
    const int n = static_cast<int>(K);
    const double a = (1.0 - std::pow(range, n)) / (1.0 - range);
    const double c = (1.0 - range) / (2.0 * range * (n - a));
    line_sum.assign(static_cast<std::size_t>(K), 0.0);
    for (int i = 0; i < n; ++i) {
        double s = 0.0;
        for (int j = 0; j < n; ++j) {
            if (i == j) continue;
            s += std::pow(range, std::abs(i - j));
        }
        line_sum[static_cast<std::size_t>(i)] =
            s * c / afreq[static_cast<std::size_t>(i)];
    }
}

// Returns an empty string when inputs are valid, otherwise the reason.
std::string check_dawid_inputs(AlleleIndex K, double range,
                               const std::vector<double>& afreq) {
    if (K < 2) {
        return "Dawid (asymmetric) mutation requires K >= 2.";
    }
    if (!std::isfinite(range) || range <= 0.0 || range >= 1.0) {
        return "Dawid (asymmetric) mutation requires range finite in "
               "(0, 1).";
    }
    if (static_cast<AlleleIndex>(afreq.size()) != K) {
        return "Dawid (asymmetric) mutation requires afreq.size() == K.";
    }
    for (AlleleIndex i = 0; i < K; ++i) {
        const double p = afreq[static_cast<std::size_t>(i)];
        if (!std::isfinite(p) || p <= 0.0) {
            return "Dawid (asymmetric) mutation requires every allele "
                   "frequency to be finite and strictly positive.";
        }
    }
    return std::string{};
}

}  // namespace

Result<double> dawid_max_rate(
        const std::vector<double>& afreq, double range) {
    const AlleleIndex K = static_cast<AlleleIndex>(afreq.size());
    const std::string chk = check_dawid_inputs(K, range, afreq);
    if (!chk.empty()) {
        return err_result<double>("dawid_max_rate: " + chk);
    }
    std::vector<double> line_sum;
    dawid_line_sums(K, range, afreq, line_sum);
    double max_line = 0.0;
    for (double s : line_sum) {
        if (s > max_line) max_line = s;
    }
    if (!(max_line > 0.0)) {
        return err_result<double>(
            "dawid_max_rate: degenerate input (all line sums zero).");
    }
    return ok_result(1.0 / max_line);
}

Result<std::vector<double>> mutation_matrix_asymmetric(
        AlleleIndex K, double rate, double range,
        const std::vector<double>& afreq) {
    const std::string chk = check_dawid_inputs(K, range, afreq);
    if (!chk.empty()) {
        return err_result<std::vector<double>>(
            "mutation_matrix_asymmetric: " + chk);
    }
    if (!rate_in_unit_interval(rate)) {
        return err_result<std::vector<double>>(
            "mutation_matrix_asymmetric: rate must be finite in [0, 1].");
    }

    std::vector<double> line_sum;
    dawid_line_sums(K, range, afreq, line_sum);

    std::vector<double> M(static_cast<std::size_t>(K) * K, 0.0);
    for (AlleleIndex i = 0; i < K; ++i) {
        const double inv_p = 1.0 / afreq[static_cast<std::size_t>(i)];
        const double a = (1.0 - std::pow(range, static_cast<int>(K)))
                         / (1.0 - range);
        const double c = (1.0 - range)
                         / (2.0 * range * (static_cast<int>(K) - a));
        double row_off = 0.0;
        for (AlleleIndex j = 0; j < K; ++j) {
            if (j == i) continue;
            const double m_ij =
                rate * c * inv_p
                * std::pow(range, std::abs(static_cast<int>(i)
                                          - static_cast<int>(j)));
            M[static_cast<std::size_t>(i) * K + j] = m_ij;
            row_off += m_ij;
        }
        const double dg = 1.0 - row_off;
        if (dg < 0.0) {
            const double max_line =
                *std::max_element(line_sum.begin(), line_sum.end());
            return err_result<std::vector<double>>(
                "mutation_matrix_asymmetric: Dawid model undefined; max "
                "rate for the given input is " +
                format_g(1.0 / max_line) + ".");
        }
        M[static_cast<std::size_t>(i) * K + i] = dg;
    }
    return ok_result(std::move(M));
}

Result<std::vector<double>> build_mutation_matrix(
        const MutationModel& mut,
        AlleleIndex K,
        const std::vector<double>& numeric_labels,
        const std::vector<double>& afreq) {
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
            if (afreq.empty()) {
                return err_result<std::vector<double>>(
                    "build_mutation_matrix: Asymmetric (Dawid) model "
                    "requires allele frequencies (afreq).");
            }
            return mutation_matrix_asymmetric(K, mut.rate, mut.range,
                                              afreq);
    }
    return err_result<std::vector<double>>(
        "build_mutation_matrix: unknown MutationKind.");
}

int mutation_models_placeholder(int x) {
    return x + 1;
}

}  // namespace core
}  // namespace mispitools
