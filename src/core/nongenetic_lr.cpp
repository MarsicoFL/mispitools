#include "nongenetic_lr.h"

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace mispitools {
namespace core {

namespace {

// Categorical CPT — mirrors ng_cpt_categorical() in r_ref_nongenetic.R.
Result<NongeneticCpt> cpt_categorical(const NongeneticFeature& f) {
    const int K = f.n_categories;
    if (K < 2) {
        return err_result<NongeneticCpt>(
            "a categorical feature needs at least two categories.");
    }
    if (f.observed_index < 0 || f.observed_index >= K) {
        return err_result<NongeneticCpt>(
            "`observed` is not one of the feature categories.");
    }

    NongeneticCpt out;
    out.p_h1.resize(static_cast<std::size_t>(K));

    if (f.error_is_matrix) {
        if (static_cast<int>(f.error_matrix.size()) != K * K) {
            return err_result<NongeneticCpt>(
                "confusion matrix must be K x K (row-major).");
        }
        const std::size_t base =
            static_cast<std::size_t>(f.observed_index) *
            static_cast<std::size_t>(K);
        for (int j = 0; j < K; ++j) {
            out.p_h1[static_cast<std::size_t>(j)] =
                f.error_matrix[base + static_cast<std::size_t>(j)];
        }
    } else {
        const double eps = f.error_scalar;
        const double off = eps / static_cast<double>(K - 1);
        for (int j = 0; j < K; ++j) {
            out.p_h1[static_cast<std::size_t>(j)] =
                (j == f.observed_index) ? (1.0 - eps) : off;
        }
    }

    out.p_h2.resize(static_cast<std::size_t>(K));
    if (f.reference_uniform) {
        const double u = 1.0 / static_cast<double>(K);
        std::fill(out.p_h2.begin(), out.p_h2.end(), u);
    } else {
        if (static_cast<int>(f.reference_freqs.size()) != K) {
            return err_result<NongeneticCpt>(
                "reference frequencies length must equal the number of "
                "categories.");
        }
        out.p_h2 = f.reference_freqs;
    }
    return ok_result(std::move(out));
}

// Continuous CPT — mirrors ng_cpt_continuous() in r_ref_nongenetic.R.
Result<NongeneticCpt> cpt_continuous(const NongeneticFeature& f) {
    NongeneticCpt out;

    if (f.reference_uniform) {
        const double a = std::floor(f.range_lo);
        const double b = std::ceil(f.range_hi);
        // seq(a, b) by 1 (R seq()): from + i, ascending; validation in
        // nongenetic_feature() guarantees range_hi > range_lo so b >= a.
        const long n = static_cast<long>(std::lround(b - a)) + 1L;
        const std::size_t G = (n < 1) ? 0u : static_cast<std::size_t>(n);
        out.grid.resize(G);
        for (std::size_t i = 0; i < G; ++i) {
            out.grid[i] = a + static_cast<double>(i);
        }
        out.p_h2.assign(G, G == 0 ? 0.0 : 1.0 / static_cast<double>(G));
    } else {
        // grid <- sort(unique(sample)); p_h2 = histogram / n.
        std::vector<double> srt = f.sample;
        std::sort(srt.begin(), srt.end());
        std::vector<double> grid;
        grid.reserve(srt.size());
        for (double v : srt) {
            if (grid.empty() || v != grid.back()) grid.push_back(v);
        }
        const std::size_t G = grid.size();
        std::vector<double> counts(G, 0.0);
        for (double v : f.sample) {
            // grid is sorted-unique of the same sample, so an exact
            // match exists (match() semantics).
            const auto it =
                std::lower_bound(grid.begin(), grid.end(), v);
            const std::size_t idx =
                static_cast<std::size_t>(it - grid.begin());
            counts[idx] += 1.0;
        }
        const double total = static_cast<double>(f.sample.size());
        out.grid = std::move(grid);
        out.p_h2.resize(G);
        for (std::size_t i = 0; i < G; ++i) {
            out.p_h2[i] = counts[i] / total;
        }
    }

    const std::size_t G = out.grid.size();
    out.p_h1.resize(G);
    if (G < 2u) {
        if (G == 1u) out.p_h1[0] = 1.0;
        return ok_result(std::move(out));
    }

    // ti <- which.min(abs(grid - observed)): first index of the minimum.
    std::size_t ti = 0;
    double best = std::fabs(out.grid[0] - f.observed_value);
    for (std::size_t i = 1; i < G; ++i) {
        const double d = std::fabs(out.grid[i] - f.observed_value);
        if (d < best) {
            best = d;
            ti = i;
        }
    }
    const double eps = f.error_scalar;
    const double off = eps / static_cast<double>(G - 1u);
    for (std::size_t i = 0; i < G; ++i) out.p_h1[i] = off;
    out.p_h1[ti] = 1.0 - eps;
    return ok_result(std::move(out));
}

// Date CPT — mirrors ng_cpt_date() in r_ref_nongenetic.R.
Result<NongeneticCpt> cpt_date(const NongeneticFeature& f) {
    const int nbins = f.n_bins;
    if (nbins < 1) {
        return err_result<NongeneticCpt>(
            "a date feature needs at least one discrepancy bin.");
    }
    if (static_cast<int>(f.alpha.size()) != nbins) {
        return err_result<NongeneticCpt>(
            "Dirichlet alpha length must equal length(cuts) + 1.");
    }

    NongeneticCpt out;
    out.p_h1.resize(static_cast<std::size_t>(nbins));
    double sum_alpha = 0.0;
    for (double a : f.alpha) sum_alpha += a;
    if (!(sum_alpha > 0.0)) {
        return err_result<NongeneticCpt>(
            "Dirichlet alpha must have positive total mass.");
    }
    for (int i = 0; i < nbins; ++i) {
        out.p_h1[static_cast<std::size_t>(i)] = f.alpha[i] / sum_alpha;
    }

    out.p_h2.resize(static_cast<std::size_t>(nbins));
    if (f.search_open) {
        const double u = 1.0 / static_cast<double>(nbins);
        std::fill(out.p_h2.begin(), out.p_h2.end(), u);
    } else {
        if (static_cast<int>(f.reference_freqs.size()) != nbins) {
            return err_result<NongeneticCpt>(
                "closed-search bin frequencies length must equal "
                "length(cuts) + 1.");
        }
        double s = 0.0;
        for (double v : f.reference_freqs) s += v;
        if (!(s > 0.0)) {
            return err_result<NongeneticCpt>(
                "closed-search bin frequencies must have positive total "
                "mass.");
        }
        for (int i = 0; i < nbins; ++i) {
            out.p_h2[static_cast<std::size_t>(i)] =
                f.reference_freqs[static_cast<std::size_t>(i)] / s;
        }
    }
    return ok_result(std::move(out));
}

// Adapt a one-dimensional non-genetic CPT into the single-member
// JointTable the genetic per-marker kernels consume. One row per
// category / grid cell / discrepancy bin; `states_flat` carries the
// 0-based state index (a single member, so it is the row index). The
// kernels read only `p_h1` / `p_h2`, but the indices keep the layout
// honest and symmetric with the genetic path.
JointTable ng_cpt_to_joint(const NongeneticCpt& c) {
    JointTable jt;
    jt.n_members = 1;
    jt.n_genotypes = static_cast<GenotypeIndex>(c.p_h1.size());
    jt.p_h1 = c.p_h1;
    jt.p_h2 = c.p_h2;
    jt.states_flat.resize(c.p_h1.size());
    for (std::size_t i = 0; i < c.p_h1.size(); ++i) {
        jt.states_flat[i] = static_cast<GenotypeIndex>(i);
    }
    return jt;
}

}  // namespace

Result<NongeneticCpt> nongenetic_cpt(const NongeneticFeature& feature) {
    switch (feature.feature_class) {
        case NgFeatureClass::Categorical:
            return cpt_categorical(feature);
        case NgFeatureClass::Continuous:
            return cpt_continuous(feature);
        case NgFeatureClass::Date:
            return cpt_date(feature);
    }
    return err_result<NongeneticCpt>("unknown non-genetic feature class.");
}

Result<PerMarkerKL> per_feature_kl_nongenetic(
        const NongeneticFeature& feature) {
    auto cpt = nongenetic_cpt(feature);
    if (!cpt.ok()) {
        return err_result<PerMarkerKL>(cpt.error);
    }
    return per_marker_kl(ng_cpt_to_joint(*cpt));
}

Result<LrDist> per_feature_lr_dist_nongenetic(
        const NongeneticFeature& feature, bool aggregate) {
    auto cpt = nongenetic_cpt(feature);
    if (!cpt.ok()) {
        return err_result<LrDist>(cpt.error);
    }
    return per_marker_lr_dist(ng_cpt_to_joint(*cpt), aggregate);
}

int nongenetic_lr_placeholder(int x) {
    return x + 1;
}

}  // namespace core
}  // namespace mispitools
