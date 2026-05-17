#ifndef MISPITOOLS_CORE_NONGENETIC_LR_H
#define MISPITOOLS_CORE_NONGENETIC_LR_H

#include <vector>

#include "result.h"

namespace mispitools {
namespace core {

/// @brief Internal feature class of a non-genetic feature.
///
/// Mirrors `nongenetic_feature()$feature_class` (R/nongenetic_feature.R):
/// every legacy `lr_*` reduces to one of these three discrete-support
/// shapes. The numeric values are part of the binding contract
/// (R passes the class as an int).
enum class NgFeatureClass {
    Categorical = 0,  ///< sex / region / hair / eyes / pigmentation
    Continuous  = 1,  ///< age
    Date        = 2   ///< birthdate (signed discrepancy bins)
};

/// @brief Plain-old-data view of a non-genetic feature for the kernel.
///
/// This is the non-genetic counterpart of `Marker`: the binding layer
/// flattens an R `nongenetic_feature` object into this POD; the kernel
/// never sees the S3 list. Only the fields relevant to `feature_class`
/// are read (see `nongenetic_cpt()`); the rest are ignored.
struct NongeneticFeature {
    NgFeatureClass feature_class = NgFeatureClass::Categorical;

    // --- categorical ---------------------------------------------------
    int n_categories = 0;                 ///< K (= categories length)
    bool error_is_matrix = false;         ///< confusion matrix vs scalar eps
    std::vector<double> error_matrix;     ///< K*K row-major, positional
    double error_scalar = 0.0;            ///< eps (categorical OR continuous)
    int observed_index = 0;               ///< 0-based TRUE-category index
    bool reference_uniform = false;       ///< model$reference == "uniform"
    std::vector<double> reference_freqs;  ///< categorical p_h2 (len K) OR
                                          ///< closed-date bin freqs (len n_bins)

    // --- continuous ----------------------------------------------------
    double range_lo = 0.0;                ///< model$range[1] (uniform)
    double range_hi = 0.0;                ///< model$range[2] (uniform)
    std::vector<double> sample;           ///< empirical sample (non-uniform)
    double observed_value = 0.0;          ///< recorded numeric observation

    // --- date ----------------------------------------------------------
    int n_bins = 0;                       ///< length(cuts) + 1
    std::vector<double> alpha;            ///< Dirichlet alpha, len n_bins
    bool search_open = true;              ///< model$search == "open"
};

/// @brief (P_H1, P_H2) conditional probability table for a non-genetic
/// feature, plus the numeric support for continuous features.
///
/// `p_h1` / `p_h2` are aligned and share the genetic-engine boundary
/// convention (a state with `p_h1 > 0, p_h2 = 0` yields log10 LR = +Inf
/// downstream; `0 * log(0/x) = 0` in sums). `grid` carries the numeric
/// support of a continuous feature so the binding can label states; it
/// is empty for categorical and date features (their labels are known
/// R-side from `categories` / `cuts`).
struct NongeneticCpt {
    std::vector<double> p_h1;
    std::vector<double> p_h2;
    std::vector<double> grid;
};

/// @brief Build the per-feature CPT under H1 and H2.
///
/// Reproduces R/r_ref_nongenetic.R::ng_cpt_R() bit-for-bit (1e-12 tol in
/// the F6.3 / F6.6 cross-checks):
///
///  - categorical: `p_h1` = row `observed_index` of the K x K confusion
///    matrix (built from `error_scalar` if `error_is_matrix` is false);
///    `p_h2` = uniform(1/K) if `reference_uniform`, else `reference_freqs`
///    (taken as-is, not renormalised, matching the R reference).
///  - continuous: support is `seq(floor(range_lo), ceil(range_hi))` when
///    `reference_uniform`, else `sort(unique(sample))`; `p_h2` is uniform
///    or the empirical histogram; `p_h1` puts mass `1 - eps` on the grid
///    cell nearest `observed_value` (first on ties) and spreads `eps`
///    uniformly over the rest.
///  - date: `p_h1` = `alpha / sum(alpha)` (the deterministic Dirichlet
///    mean); `p_h2` = uniform(1/n_bins) for an open search, else
///    `reference_freqs / sum(reference_freqs)`.
///
/// Errors (invalid K, out-of-range observed index, degenerate closed-date
/// mass) flow back as `Result::error` — never thrown across the boundary.
///
/// @complexity O(card(feature)) plus O(n log n) for the empirical-sample
/// sort in the continuous path.
Result<NongeneticCpt> nongenetic_cpt(const NongeneticFeature& feature);

// Placeholder retained for the F0.5 cpp-bootstrap regression test.
int nongenetic_lr_placeholder(int x);

}  // namespace core
}  // namespace mispitools

#endif  // MISPITOOLS_CORE_NONGENETIC_LR_H
