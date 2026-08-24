#include "lr_dist.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <numeric>
#include <string>
#include <utility>

namespace mispitools {
namespace core {

int lr_dist_placeholder(int x) {
    return x + 1;
}

// Mirrors R/r_ref_per_marker.R::per_marker_lr_dist_R bit-for-bit.
//
// Pattern reference for the F4.2 composition step: DNAtools::convolve
// (GPL-2+) — re-implemented, not lifted. The sparse-sorted (log10_lr,
// p_h1, p_h2) layout below is what `lr_dist_compose` will fold over.
Result<LrDist> per_marker_lr_dist(const JointTable& joint, bool aggregate) {
    if (joint.p_h1.size() != joint.p_h2.size()) {
        return err_result<LrDist>(
            "per_marker_lr_dist: P_H1 and P_H2 vectors have mismatched "
            "lengths.");
    }

    const std::size_t n = joint.p_h1.size();
    const double inf = std::numeric_limits<double>::infinity();

    LrDist out;
    out.log10_lr.reserve(n);
    out.p_h1.reserve(n);
    out.p_h2.reserve(n);

    // Raw atoms in joint row order (R-ref: log10_lr_from_probs then the
    // unaggregated data.frame). Rows with both probabilities zero are
    // skipped (the `0 log 0 = 0` limit; cpt_marker_joint already filters
    // them, so this only guards a hand-built joint).
    for (std::size_t r = 0; r < n; ++r) {
        const double p1 = joint.p_h1[r];
        const double p2 = joint.p_h2[r];
        if (p1 < 0.0 || p2 < 0.0) {
            return err_result<LrDist>(
                "per_marker_lr_dist: negative probability at row "
                + std::to_string(r) + ".");
        }
        const bool pos1 = p1 > 0.0;
        const bool pos2 = p2 > 0.0;

        if (pos1 && pos2) {
            // R-ref: log10(P1) - log10(P2) (not log10(P1/P2)).
            out.log10_lr.push_back(std::log10(p1) - std::log10(p2));
        } else if (pos1) {
            out.log10_lr.push_back(inf);
            out.has_pos_inf = true;
        } else if (pos2) {
            out.log10_lr.push_back(-inf);
            out.has_neg_inf = true;
        } else {
            continue;  // both zero → no atom
        }
        out.p_h1.push_back(p1);
        out.p_h2.push_back(p2);
    }

    if (!aggregate || out.log10_lr.empty()) {
        return ok_result(std::move(out));
    }

    // Aggregate: stable-sort ascending by log10_lr (R `order()` is a
    // radix sort — stable; -Inf first, +Inf last), then collapse equal
    // keys summing the probabilities in the post-sort order (R `tapply`
    // sum). IEEE makes `Inf == Inf` and `-Inf == -Inf`, so the infinite
    // atoms each collapse into a single bucket, matching R-ref.
    const std::size_t m = out.log10_lr.size();
    std::vector<std::size_t> idx(m);
    std::iota(idx.begin(), idx.end(), std::size_t{0});
    std::stable_sort(idx.begin(), idx.end(),
                     [&out](std::size_t a, std::size_t b) {
                         return out.log10_lr[a] < out.log10_lr[b];
                     });

    LrDist agg;
    agg.has_pos_inf = out.has_pos_inf;
    agg.has_neg_inf = out.has_neg_inf;

    std::size_t i = 0;
    while (i < m) {
        const double key = out.log10_lr[idx[i]];
        double s1 = 0.0;
        double s2 = 0.0;
        std::size_t j = i;
        // `!(key != next)` keeps Inf/-Inf grouped via IEEE equality,
        // exactly like R `k[-1] != k[-length(k)]`.
        while (j < m && !(out.log10_lr[idx[j]] != key)) {
            s1 += out.p_h1[idx[j]];
            s2 += out.p_h2[idx[j]];
            ++j;
        }
        agg.log10_lr.push_back(key);  // R-ref uses the group's first key
        agg.p_h1.push_back(s1);
        agg.p_h2.push_back(s2);
        i = j;
    }

    return ok_result(std::move(agg));
}

namespace {

// Sort ascending by key (stable, like R `order()` radix) and collapse
// adjacent keys whose gap is `<= merge_tol` (representative = the group's
// first key, summing the parallel probabilities in post-sort order — the
// same convention as `aggregate_lr_dist`). IEEE makes Inf == Inf, and
// `Inf - Inf` is NaN (any comparison false), so each infinite class
// collapses into its own bucket regardless of `merge_tol`.
LrDist aggregate_sparse(std::vector<double> lr,
                        std::vector<double> p1,
                        std::vector<double> p2,
                        bool has_pos_inf,
                        bool has_neg_inf,
                        double merge_tol) {
    LrDist out;
    out.has_pos_inf = has_pos_inf;
    out.has_neg_inf = has_neg_inf;
    const std::size_t m = lr.size();
    if (m == 0) return out;

    std::vector<std::size_t> idx(m);
    std::iota(idx.begin(), idx.end(), std::size_t{0});
    std::stable_sort(idx.begin(), idx.end(),
                     [&lr](std::size_t a, std::size_t b) {
                         return lr[a] < lr[b];
                     });

    std::size_t i = 0;
    while (i < m) {
        const double key = lr[idx[i]];
        double s1 = 0.0;
        double s2 = 0.0;
        std::size_t j = i;
        while (j < m) {
            const double k = lr[idx[j]];
            const bool same = (k == key) ||
                              (std::isfinite(k) && std::isfinite(key) &&
                               (k - key) <= merge_tol);
            if (!same) break;
            s1 += p1[idx[j]];
            s2 += p2[idx[j]];
            ++j;
        }
        out.log10_lr.push_back(key);
        out.p_h1.push_back(s1);
        out.p_h2.push_back(s2);
        i = j;
    }
    return out;
}

LrDist identity_delta() {
    LrDist d;
    d.log10_lr = {0.0};
    d.p_h1 = {1.0};
    d.p_h2 = {1.0};
    return d;
}

bool dist_has_infinite(const LrDist& d) {
    for (double v : d.log10_lr) {
        if (!std::isfinite(v)) return true;
    }
    return false;
}

// Exact sparse convolution fold. The cartesian product is iterated in the
// same order R-ref `lr_dist_compose_R` uses (feature index outermost,
// accumulator innermost) so the post-sort summation order in
// `aggregate_sparse` matches bit-for-bit.
Result<LrDist> compose_exact(const std::vector<LrDist>& per_feature,
                             double merge_tol) {
    LrDist acc = identity_delta();
    for (const LrDist& d : per_feature) {
        if (d.p_h1.size() != d.log10_lr.size() ||
            d.p_h2.size() != d.log10_lr.size()) {
            return err_result<LrDist>(
                "lr_dist_compose: a per-feature distribution has mismatched "
                "column lengths.");
        }
        std::vector<double> lr;
        std::vector<double> p1;
        std::vector<double> p2;
        const std::size_t cap = acc.size() * d.size();
        lr.reserve(cap);
        p1.reserve(cap);
        p2.reserve(cap);
        bool pos_inf = false;
        bool neg_inf = false;
        for (std::size_t jd = 0; jd < d.size(); ++jd) {
            for (std::size_t ja = 0; ja < acc.size(); ++ja) {
                const double np1 = acc.p_h1[ja] * d.p_h1[jd];
                const double np2 = acc.p_h2[ja] * d.p_h2[jd];
                if (!(np1 > 0.0) && !(np2 > 0.0)) continue;
                const double nlr = acc.log10_lr[ja] + d.log10_lr[jd];
                lr.push_back(nlr);
                p1.push_back(np1);
                p2.push_back(np2);
                if (nlr == std::numeric_limits<double>::infinity()) {
                    pos_inf = true;
                } else if (nlr ==
                           -std::numeric_limits<double>::infinity()) {
                    neg_inf = true;
                }
            }
        }
        acc = aggregate_sparse(std::move(lr), std::move(p1), std::move(p2),
                               pos_inf, neg_inf, merge_tol);
        if (acc.empty()) {
            // Every combined atom carried zero mass — only possible from a
            // degenerate input. Return the empty distribution explicitly.
            return ok_result(std::move(acc));
        }
    }
    return ok_result(std::move(acc));
}

// Fixed-lattice heuristic. Each feature atom is split linearly between the
// two nearest lattice nodes (origin 0, spacing `delta`), which preserves
// total mass and the first moment exactly. Convolution is then exact
// integer-index addition (no offset bug, no clamping). Discretisation
// error on the shape is O(delta).
Result<LrDist> compose_grid(const std::vector<LrDist>& per_feature,
                            int grid_points) {
    if (grid_points < 2) {
        return err_result<LrDist>(
            "lr_dist_compose: grid_points must be >= 2 for the grid method.");
    }
    double g_lo = 0.0;
    double g_hi = 0.0;
    for (const LrDist& d : per_feature) {
        if (d.p_h1.size() != d.log10_lr.size() ||
            d.p_h2.size() != d.log10_lr.size()) {
            return err_result<LrDist>(
                "lr_dist_compose: a per-feature distribution has mismatched "
                "column lengths.");
        }
        if (dist_has_infinite(d)) {
            return err_result<LrDist>(
                "lr_dist_compose: the grid method requires finite supports "
                "(a feature has a +/-Inf atom); use the exact method or set "
                "a positive mutation rate.");
        }
        double lo = std::numeric_limits<double>::infinity();
        double hi = -std::numeric_limits<double>::infinity();
        for (double v : d.log10_lr) {
            lo = std::min(lo, v);
            hi = std::max(hi, v);
        }
        if (!d.empty()) {
            g_lo += lo;
            g_hi += hi;
        }
    }
    const double span = g_hi - g_lo;
    if (!(span > 0.0)) {
        // Degenerate (every feature is a single point): the exact fold is
        // already a single atom; defer to it (merge_tol irrelevant).
        return compose_exact(per_feature, 0.0);
    }
    const double delta = span / static_cast<double>(grid_points - 1);

    // running[i] holds (p_h1, p_h2) mass at lattice value i * delta.
    std::map<long, std::pair<double, double>> running;
    running[0] = {1.0, 1.0};

    for (const LrDist& d : per_feature) {
        std::map<long, std::pair<double, double>> feat;
        for (std::size_t k = 0; k < d.size(); ++k) {
            const double x = d.log10_lr[k] / delta;
            const double fl = std::floor(x);
            const long i0 = static_cast<long>(fl);
            const double fr = x - fl;
            auto& a = feat[i0];
            a.first  += (1.0 - fr) * d.p_h1[k];
            a.second += (1.0 - fr) * d.p_h2[k];
            if (fr > 0.0) {
                auto& b = feat[i0 + 1];
                b.first  += fr * d.p_h1[k];
                b.second += fr * d.p_h2[k];
            }
        }
        std::map<long, std::pair<double, double>> next;
        for (const auto& ra : running) {
            for (const auto& fb : feat) {
                auto& c = next[ra.first + fb.first];
                c.first  += ra.second.first  * fb.second.first;
                c.second += ra.second.second * fb.second.second;
            }
        }
        running.swap(next);
    }

    LrDist out;
    out.log10_lr.reserve(running.size());
    out.p_h1.reserve(running.size());
    out.p_h2.reserve(running.size());
    for (const auto& kv : running) {
        const double p1 = kv.second.first;
        const double p2 = kv.second.second;
        if (!(p1 > 0.0) && !(p2 > 0.0)) continue;
        out.log10_lr.push_back(static_cast<double>(kv.first) * delta);
        out.p_h1.push_back(p1);
        out.p_h2.push_back(p2);
    }
    return ok_result(std::move(out));
}

}  // namespace

Result<LrDist> lr_dist_compose(const std::vector<LrDist>& per_feature,
                               const LrDistComposeOptions& opts) {
    if (opts.merge_tol < 0.0) {
        return err_result<LrDist>(
            "lr_dist_compose: merge_tol must be non-negative.");
    }
    if (per_feature.empty()) {
        return ok_result(identity_delta());
    }
    switch (opts.method) {
        case ComposeMethod::Exact:
            return compose_exact(per_feature, opts.merge_tol);
        case ComposeMethod::Grid:
            return compose_grid(per_feature, opts.grid_points);
    }
    return err_result<LrDist>("lr_dist_compose: unknown method.");
}

}  // namespace core
}  // namespace mispitools
