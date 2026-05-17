// F4.3 — decision-theoretic primitives over a sparse LR distribution.
//
// Every quantity is derived from the `LrDist` columns (the F4.1 /
// F4.2 output) and mirrors a reference in
// `R/r_ref_per_marker.R` (lr_dist_summary_R, lr_dist_quantile_R,
// decision_rates_R, roc_curve_R, choose_threshold_weighted_R)
// bit-for-bit: identical active-weight filtering, identical
// left-to-right (ascending) summation order, identical tie / ±Inf
// conventions. Pure C++17, no Rcpp/R headers; errors flow back through
// `Result<T>` (see DESIGN.md §5 / §8.8).

#include "decision.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace mispitools {
namespace core {

namespace {

// Shared shape / sign validation. The three columns are row-aligned;
// negative probabilities indicate a malformed joint upstream.
std::string validate(const LrDist& d) {
    if (d.log10_lr.size() != d.p_h1.size()
            || d.log10_lr.size() != d.p_h2.size()) {
        return "log10_lr, p_h1 and p_h2 must have the same length.";
    }
    for (std::size_t i = 0; i < d.p_h1.size(); ++i) {
        if (d.p_h1[i] < 0.0 || d.p_h2[i] < 0.0) {
            return "probabilities must be non-negative.";
        }
    }
    return std::string{};
}

// Σ weight·log10_lr over active (weight > 0) atoms, scanned in stored
// (ascending) order — bit-for-bit with weighted_log10_lr_sum() in R.
double weighted_log10_lr_sum(const std::vector<double>& w,
                             const std::vector<double>& lr) {
    bool any_active = false;
    double acc = 0.0;
    for (std::size_t i = 0; i < w.size(); ++i) {
        if (w[i] > 0.0) {
            any_active = true;
            acc += w[i] * lr[i];
        }
    }
    return any_active ? acc : 0.0;
}

// Σ weight·(log10_lr - m)² over active atoms. Matches the R reference:
// no active mass → 0; non-finite mean → +Inf; otherwise the centred
// sum in ascending order (the active log10_lr values are finite here,
// since a ±Inf atom only carries weight under the hypothesis whose
// mean is then already non-finite).
double weighted_centered_sq(const std::vector<double>& w,
                            const std::vector<double>& lr,
                            double m) {
    bool any_active = false;
    for (std::size_t i = 0; i < w.size(); ++i) {
        if (w[i] > 0.0) { any_active = true; break; }
    }
    if (!any_active) return 0.0;
    if (!std::isfinite(m)) return std::numeric_limits<double>::infinity();
    double acc = 0.0;
    for (std::size_t i = 0; i < w.size(); ++i) {
        if (w[i] > 0.0) {
            const double dlt = lr[i] - m;
            acc += w[i] * dlt * dlt;
        }
    }
    return acc;
}

// Ascending-sorted distinct atom values (stable; ±Inf included).
std::vector<double> distinct_sorted(const std::vector<double>& lr) {
    std::vector<double> v(lr);
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
    return v;
}

}  // namespace

Result<LrDistSummary> lr_dist_summary(const LrDist& d) {
    const std::string e = validate(d);
    if (!e.empty()) return err_result<LrDistSummary>(e);

    LrDistSummary s;
    s.mean_h1 = weighted_log10_lr_sum(d.p_h1, d.log10_lr);
    s.mean_h2 = weighted_log10_lr_sum(d.p_h2, d.log10_lr);
    s.var_h1 = weighted_centered_sq(d.p_h1, d.log10_lr, s.mean_h1);
    s.var_h2 = weighted_centered_sq(d.p_h2, d.log10_lr, s.mean_h2);
    s.sd_h1 = std::sqrt(s.var_h1);
    s.sd_h2 = std::sqrt(s.var_h2);

    double m1 = 0.0, m2 = 0.0;
    for (std::size_t i = 0; i < d.log10_lr.size(); ++i) {
        m1 += d.p_h1[i];
        m2 += d.p_h2[i];
        if (std::isinf(d.log10_lr[i])) {
            if (d.log10_lr[i] > 0.0) s.has_pos_inf = true;
            else s.has_neg_inf = true;
        }
    }
    s.mass_h1 = m1;
    s.mass_h2 = m2;
    return ok_result(s);
}

Result<std::vector<double>> lr_dist_quantile(
        const LrDist& d,
        const std::vector<double>& probs,
        bool under_h1) {
    const std::string e = validate(d);
    if (!e.empty()) return err_result<std::vector<double>>(e);

    const std::vector<double>& w = under_h1 ? d.p_h1 : d.p_h2;

    // Active atoms (weight > 0), sorted ascending by log10_lr (stable).
    std::vector<std::size_t> idx;
    idx.reserve(w.size());
    for (std::size_t i = 0; i < w.size(); ++i) {
        if (w[i] > 0.0) idx.push_back(i);
    }
    if (idx.empty()) {
        return err_result<std::vector<double>>(
            "no active support under the requested hypothesis.");
    }
    std::stable_sort(idx.begin(), idx.end(),
                     [&](std::size_t a, std::size_t b) {
                         return d.log10_lr[a] < d.log10_lr[b];
                     });

    double total = 0.0;
    for (std::size_t k : idx) total += w[k];

    std::vector<double> cw(idx.size());
    double run = 0.0;
    for (std::size_t j = 0; j < idx.size(); ++j) {
        run += w[idx[j]];
        cw[j] = run / total;
    }

    const double slack = 1e-12;
    std::vector<double> out;
    out.reserve(probs.size());
    for (double p : probs) {
        if (p <= 0.0) {
            out.push_back(d.log10_lr[idx.front()]);
            continue;
        }
        std::size_t hit = idx.size() - 1;  // p >= 1 → largest atom
        for (std::size_t j = 0; j < cw.size(); ++j) {
            if (cw[j] + slack >= p) { hit = j; break; }
        }
        out.push_back(d.log10_lr[idx[hit]]);
    }
    return ok_result(out);
}

Result<DecisionRates> decision_rates(const LrDist& d, double threshold) {
    const std::string e = validate(d);
    if (!e.empty()) return err_result<DecisionRates>(e);

    double fnr = 0.0, tpr = 0.0, fpr = 0.0, tnr = 0.0;
    for (std::size_t i = 0; i < d.log10_lr.size(); ++i) {
        const double lr = d.log10_lr[i];
        if (lr < threshold) { fnr += d.p_h1[i]; tnr += d.p_h2[i]; }
        else if (lr > threshold) { tpr += d.p_h1[i]; fpr += d.p_h2[i]; }
        // lr == threshold → indeterminate (mirrors threshold_rates).
    }

    DecisionRates r;
    r.threshold = threshold;
    r.fpr = fpr;
    r.fnr = fnr;
    r.tpr = tpr;
    r.tnr = tnr;

    const double TP = tpr, TN = tnr, FP = fpr, FN = fnr;
    const double denom = std::sqrt((TP + FP) * (TP + FN)
                                   * (TN + FP) * (TN + FN));
    r.mcc = (denom == 0.0) ? 0.0
                           : (TP * TN - FP * FN) / denom;
    return ok_result(r);
}

Result<RocCurve> roc_curve(const LrDist& d) {
    const std::string e = validate(d);
    if (!e.empty()) return err_result<RocCurve>(e);

    RocCurve c;
    const std::vector<double> thr = distinct_sorted(d.log10_lr);
    c.threshold = thr;
    c.fpr.resize(thr.size());
    c.tpr.resize(thr.size());
    c.fnr.resize(thr.size());
    c.tnr.resize(thr.size());
    for (std::size_t t = 0; t < thr.size(); ++t) {
        double fpr = 0.0, tpr = 0.0, fnr = 0.0, tnr = 0.0;
        for (std::size_t i = 0; i < d.log10_lr.size(); ++i) {
            const double lr = d.log10_lr[i];
            if (lr < thr[t]) { fnr += d.p_h1[i]; tnr += d.p_h2[i]; }
            else if (lr > thr[t]) { tpr += d.p_h1[i]; fpr += d.p_h2[i]; }
        }
        c.fpr[t] = fpr;
        c.tpr[t] = tpr;
        c.fnr[t] = fnr;
        c.tnr[t] = tnr;
    }

    // Concordance AUC: P(LR_H1 > LR_H2) + ½ P(LR_H1 = LR_H2), the i
    // (outer) / j (inner) nesting matching roc_curve_R.
    double auc = 0.0;
    for (std::size_t i = 0; i < d.log10_lr.size(); ++i) {
        if (d.p_h1[i] <= 0.0) continue;
        for (std::size_t j = 0; j < d.log10_lr.size(); ++j) {
            if (d.p_h2[j] <= 0.0) continue;
            double ind;
            if (d.log10_lr[i] > d.log10_lr[j]) ind = 1.0;
            else if (d.log10_lr[i] == d.log10_lr[j]) ind = 0.5;
            else ind = 0.0;
            auc += d.p_h1[i] * d.p_h2[j] * ind;
        }
    }
    c.auc = auc;
    return ok_result(c);
}

Result<ThresholdChoice> choose_threshold_weighted(const LrDist& d,
                                                  double weight) {
    const std::string e = validate(d);
    if (!e.empty()) return err_result<ThresholdChoice>(e);
    if (!(weight > 0.0)) {
        return err_result<ThresholdChoice>("weight must be positive.");
    }

    const std::vector<double> thr = distinct_sorted(d.log10_lr);
    ThresholdChoice best;
    bool have = false;
    for (double t : thr) {
        double fnr = 0.0, fpr = 0.0;
        for (std::size_t i = 0; i < d.log10_lr.size(); ++i) {
            const double lr = d.log10_lr[i];
            if (lr < t) fnr += d.p_h1[i];
            else if (lr > t) fpr += d.p_h2[i];
        }
        const double wf = weight * fpr;
        const double dist = std::sqrt(fnr * fnr + wf * wf);
        if (!have || dist < best.distance) {
            have = true;
            best.threshold = t;
            best.fpr = fpr;
            best.fnr = fnr;
            best.distance = dist;
        }
    }
    return ok_result(best);
}

int decision_placeholder(int x) {
    return x + 1;
}

}  // namespace core
}  // namespace mispitools
