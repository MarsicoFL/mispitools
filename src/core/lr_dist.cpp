#include "lr_dist.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
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

}  // namespace core
}  // namespace mispitools
