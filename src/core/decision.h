#ifndef MISPITOOLS_CORE_DECISION_H
#define MISPITOOLS_CORE_DECISION_H

#include <vector>

#include "lr_dist.h"
#include "result.h"

namespace mispitools {
namespace core {

/// @brief First/second moments of the log10 LR distribution under each
/// hypothesis.
///
/// Computed over the sparse `LrDist` support. `mean_h1` is
/// `E[log10 LR | H1] = Σ p_h1 · log10 LR` (the F1.6/F3.1
/// `e_log10_lr_h1`); `mean_h2` is the H2 analogue. Variances use the
/// same active-weight filter as the mean (`weight > 0`), so a `+Inf`
/// atom (which carries `p_h1 > 0`, `p_h2 = 0`) drives `mean_h1` and
/// `var_h1` to `+Inf` while leaving the H2 moments finite, bit-for-bit
/// with `R/r_ref_per_marker.R::lr_dist_summary_R`. `mass_h1` / `mass_h2`
/// are the total probability (≈ 1 on a well-formed joint); they expose
/// a malformed input rather than silently renormalising.
struct LrDistSummary {
    double mean_h1 = 0.0;   ///< E[log10 LR | H1]
    double mean_h2 = 0.0;   ///< E[log10 LR | H2]
    double var_h1 = 0.0;    ///< Var[log10 LR | H1]
    double var_h2 = 0.0;    ///< Var[log10 LR | H2]
    double sd_h1 = 0.0;     ///< sqrt(var_h1)
    double sd_h2 = 0.0;     ///< sqrt(var_h2)
    double mass_h1 = 0.0;   ///< Σ p_h1
    double mass_h2 = 0.0;   ///< Σ p_h2
    bool has_pos_inf = false;   ///< a (+Inf) atom present
    bool has_neg_inf = false;   ///< a (-Inf) atom present
};

/// @brief Moments of the log10 LR distribution.
/// @param d sparse LR distribution (output of `per_marker_lr_dist` or
///        `lr_dist_compose`); the three columns are row-aligned.
/// @return `LrDistSummary`; errors on a length mismatch or a negative
///         probability.
/// @complexity O(n_atoms), two linear passes (mean then variance).
Result<LrDistSummary> lr_dist_summary(const LrDist& d);

/// @brief Quantiles of the log10 LR distribution under one hypothesis.
///
/// Discrete inverse-CDF (R `quantile` type 1):
/// `Q(p) = inf{ x : F(x) >= p }`, where `F` is the cumulative mass of
/// the atoms with `weight > 0` (`p_h1` if `under_h1`, else `p_h2`),
/// normalised by the total active mass, scanned ascending. `p <= 0`
/// returns the smallest active atom; `p >= 1` the largest. The
/// cumulative comparison carries a `1e-12` slack so a probability that
/// lands exactly on a cumulative breakpoint resolves to that atom,
/// bit-for-bit with `lr_dist_quantile_R`.
///
/// @param d        sparse LR distribution.
/// @param probs    target probabilities (any order; not required sorted).
/// @param under_h1 weight by `p_h1` (true) or `p_h2` (false).
/// @return one quantile per `probs` entry, in input order; errors on a
///         length mismatch, a negative probability, or an empty active
///         support.
/// @complexity O(n_atoms log n_atoms + n_atoms · probs.size()).
Result<std::vector<double>> lr_dist_quantile(const LrDist& d,
                                             const std::vector<double>& probs,
                                             bool under_h1);

/// @brief Error rates of the LR-threshold classifier at one cut point.
///
/// Decision rule: declare a match when `log10 LR > threshold`. Mass at
/// exactly `log10 LR == threshold` is left indeterminate (excluded from
/// both the positive and negative tallies), so `tpr != 1 - fnr` when an
/// atom sits on the threshold — the analytic analogue of
/// `R/threshold_rates.R` (`sum(TPED > t)`, `sum(TPED < t)`, …) over the
/// exact distribution rather than a Monte-Carlo sample.
///   * `fnr = P(log10 LR < threshold | H1)`
///   * `tpr = P(log10 LR > threshold | H1)`
///   * `fpr = P(log10 LR > threshold | H2)`
///   * `tnr = P(log10 LR < threshold | H2)`
///   * `mcc`: Matthews correlation over the (tpr, tnr, fpr, fnr)
///     confusion (each a probability; H1 mass = H2 mass = 1); `0` when
///     the denominator vanishes (mirrors `threshold_rates`).
struct DecisionRates {
    double threshold = 0.0;
    double fpr = 0.0;
    double fnr = 0.0;
    double tpr = 0.0;
    double tnr = 0.0;
    double mcc = 0.0;
};

/// @brief Error rates at a given log10 LR threshold.
/// @param d         sparse LR distribution.
/// @param threshold cut point on `log10 LR`.
/// @return `DecisionRates`; errors on a length mismatch or a negative
///         probability.
/// @complexity O(n_atoms).
Result<DecisionRates> decision_rates(const LrDist& d, double threshold);

/// @brief Analytic ROC curve of the LR-threshold classifier.
///
/// Evaluates `(fpr, tpr)` at every distinct atom value (the natural
/// breakpoints of the step ROC), ordered by ascending threshold. `auc`
/// is the exact rank statistic
/// `P(LR_{H1} > LR_{H2}) + ½ P(LR_{H1} = LR_{H2})` (Mann–Whitney /
/// concordance), which handles ties and ±Inf atoms without quadrature
/// error and equals the trapezoidal area under the returned step curve.
struct RocCurve {
    std::vector<double> threshold;  ///< distinct atom values, ascending
    std::vector<double> fpr;        ///< P(log10 LR > threshold | H2)
    std::vector<double> tpr;        ///< P(log10 LR > threshold | H1)
    std::vector<double> fnr;        ///< 1 - tpr is not assumed; see DecisionRates
    std::vector<double> tnr;        ///< P(log10 LR < threshold | H2)
    double auc = 0.0;               ///< concordance AUC
};

/// @brief ROC curve + AUC of the LR-threshold classifier.
/// @param d sparse LR distribution.
/// @return `RocCurve` with one row per distinct atom value; errors on a
///         length mismatch or a negative probability.
/// @complexity O(n_atoms²) for the concordance AUC (sparse n is small).
Result<RocCurve> roc_curve(const LrDist& d);

/// @brief Weighted-Euclidean optimal threshold.
///
/// Analytic analogue of `R/decision_threshold.R`: over the distinct
/// atom values it minimises `D = sqrt(fnr² + (weight · fpr)²)` and
/// returns the minimising threshold (first one on a tie, scanning
/// ascending) with its rates. `weight > 1` penalises false positives.
struct ThresholdChoice {
    double threshold = 0.0;
    double fpr = 0.0;
    double fnr = 0.0;
    double distance = 0.0;
};

/// @brief Optimal threshold under a weighted FPR/FNR cost.
/// @param d      sparse LR distribution.
/// @param weight relative cost of a false positive (must be > 0).
/// @return `ThresholdChoice`; errors on a length mismatch, a negative
///         probability, or `weight <= 0`.
/// @complexity O(n_atoms²) (rates at each candidate threshold).
Result<ThresholdChoice> choose_threshold_weighted(const LrDist& d,
                                                   double weight);

/// Placeholder retained for the F0.5 cpp-bootstrap regression test.
/// @param x integer input
/// @return x + 1
int decision_placeholder(int x);

}  // namespace core
}  // namespace mispitools

#endif  // MISPITOOLS_CORE_DECISION_H
