#include "evidence_combine.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <string>
#include <utility>

namespace mispitools {
namespace core {

int evidence_combine_placeholder(int x) {
    return x + 1;
}

namespace {

// Soft cap on the MarkovSE cartesian fold. SE supports are tiny (sex = 2,
// pigmentation ≤ 3, age bins, …) so a realistic case stays far below this;
// the cap only guards a pathological input from running away.
constexpr std::size_t kMarkovSEMaxPaths = 5'000'000;

constexpr double kRowSumTol = 1e-9;  // user-supplied transition rows

// Aggregate sparse atoms with the exact convention of
// `per_marker_lr_dist` (lr_dist.cpp): stable sort ascending by log10_lr,
// collapse the keys that `same_atom` calls one atom (so each ±Inf class
// folds into a single bucket and last-bit duplicates fold with it),
// representative = the group's first key, probabilities summed in post-sort
// order.
LrDist aggregate_atoms(std::vector<double> lr,
                       std::vector<double> p1,
                       std::vector<double> p2,
                       bool has_pos_inf,
                       bool has_neg_inf) {
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
        while (j < m && same_atom(lr[idx[j]], key)) {
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

// One partial path of the MarkovSE fold: cumulative H1 / H2 mass and the
// atom index of the most recent feature (needed to condition the next
// transition row).
struct PartialPath {
    double p1;
    double p2;
    std::size_t last_idx;
};

Result<LrDist> combine_markov_se(const std::vector<LrDist>& per_feature,
                                 const MarkovSEParams& params) {
    const std::size_t K = per_feature.size();
    if (K == 0) {
        return err_result<LrDist>(
            "evidence_combine: MarkovSE requires at least one feature.");
    }
    for (std::size_t s = 0; s < K; ++s) {
        const LrDist& f = per_feature[s];
        if (f.p_h1.size() != f.p_h2.size()) {
            return err_result<LrDist>(
                "evidence_combine: feature " + std::to_string(s)
                + " has mismatched p_h1 / p_h2 lengths.");
        }
        if (f.empty()) {
            return err_result<LrDist>(
                "evidence_combine: feature " + std::to_string(s)
                + " is empty.");
        }
    }
    if (params.transition.size() != K - 1) {
        return err_result<LrDist>(
            "evidence_combine: MarkovSE expects K-1 transition matrices ("
            + std::to_string(K - 1) + "), got "
            + std::to_string(params.transition.size()) + ".");
    }
    for (std::size_t s = 0; s + 1 < K; ++s) {
        const std::size_t ns   = per_feature[s].size();
        const std::size_t nnxt = per_feature[s + 1].size();
        const std::vector<double>& T = params.transition[s];
        if (T.size() != ns * nnxt) {
            return err_result<LrDist>(
                "evidence_combine: transition[" + std::to_string(s)
                + "] must be " + std::to_string(ns) + "x"
                + std::to_string(nnxt) + " (row-major), got "
                + std::to_string(T.size()) + " entries.");
        }
        for (std::size_t i = 0; i < ns; ++i) {
            double rs = 0.0;
            for (std::size_t j = 0; j < nnxt; ++j) {
                const double v = T[i * nnxt + j];
                if (v < 0.0) {
                    return err_result<LrDist>(
                        "evidence_combine: transition[" + std::to_string(s)
                        + "] has a negative entry.");
                }
                rs += v;
            }
            if (std::abs(rs - 1.0) > kRowSumTol) {
                return err_result<LrDist>(
                    "evidence_combine: transition[" + std::to_string(s)
                    + "] row " + std::to_string(i)
                    + " is not stochastic (sums to " + std::to_string(rs)
                    + ").");
            }
        }
    }

    // Path-count blow-up guard (overflow-safe product).
    std::size_t projected = 1;
    for (std::size_t s = 0; s < K; ++s) {
        const std::size_t ns = per_feature[s].size();
        if (ns != 0 && projected > kMarkovSEMaxPaths / ns) {
            return err_result<LrDist>(
                "evidence_combine: MarkovSE joint support exceeds the "
                + std::to_string(kMarkovSEMaxPaths) + "-path cap; reduce "
                "the per-feature support or use Independent mode.");
        }
        projected *= ns;
    }

    // Seed the fold with feature 0 (H2 marginal = its own p_h2).
    std::vector<PartialPath> paths;
    {
        const LrDist& f0 = per_feature[0];
        for (std::size_t j = 0; j < f0.size(); ++j) {
            const double a = f0.p_h1[j];
            const double b = f0.p_h2[j];
            if (a < 0.0 || b < 0.0) {
                return err_result<LrDist>(
                    "evidence_combine: feature 0 has a negative probability.");
            }
            if (!(a > 0.0) && !(b > 0.0)) continue;  // 0/0 → no atom
            paths.push_back(PartialPath{a, b, j});
        }
    }

    // Chain over the remaining features: H1 stays independent, H2 follows
    // the supplied transition row of the previous feature's atom.
    for (std::size_t s = 1; s < K; ++s) {
        const LrDist& fs = per_feature[s];
        const std::size_t ns = fs.size();
        const std::vector<double>& T = params.transition[s - 1];
        std::vector<PartialPath> next;
        next.reserve(paths.size() * ns);
        for (const PartialPath& pp : paths) {
            const std::size_t base = pp.last_idx * ns;
            for (std::size_t j = 0; j < ns; ++j) {
                const double a = fs.p_h1[j];
                if (a < 0.0) {
                    return err_result<LrDist>(
                        "evidence_combine: feature " + std::to_string(s)
                        + " has a negative probability.");
                }
                const double np1 = pp.p1 * a;
                const double np2 = pp.p2 * T[base + j];
                if (!(np1 > 0.0) && !(np2 > 0.0)) continue;
                next.push_back(PartialPath{np1, np2, j});
            }
        }
        paths.swap(next);
    }

    // Materialise the joint atoms (same boundary convention as
    // per_marker_lr_dist: log10(p1) - log10(p2), not log10(p1/p2)).
    const double inf = std::numeric_limits<double>::infinity();
    std::vector<double> lr;
    std::vector<double> p1;
    std::vector<double> p2;
    lr.reserve(paths.size());
    p1.reserve(paths.size());
    p2.reserve(paths.size());
    bool pos_inf = false;
    bool neg_inf = false;
    for (const PartialPath& pp : paths) {
        const bool a = pp.p1 > 0.0;
        const bool b = pp.p2 > 0.0;
        if (a && b) {
            lr.push_back(std::log10(pp.p1) - std::log10(pp.p2));
        } else if (a) {
            lr.push_back(inf);
            pos_inf = true;
        } else if (b) {
            lr.push_back(-inf);
            neg_inf = true;
        } else {
            continue;  // unreachable (filtered above), kept for safety
        }
        p1.push_back(pp.p1);
        p2.push_back(pp.p2);
    }

    return ok_result(aggregate_atoms(std::move(lr), std::move(p1),
                                     std::move(p2), pos_inf, neg_inf));
}

}  // namespace

Result<LrDist> evidence_combine(const std::vector<LrDist>& per_feature,
                                CombineMode mode,
                                const MarkovSEParams& params) {
    switch (mode) {
        case CombineMode::Independent:
            // DESIGN.md §8.7: identical to §8.5 — exact convolution.
            return lr_dist_compose(per_feature, LrDistComposeOptions{});
        case CombineMode::MarkovSE:
            return combine_markov_se(per_feature, params);
    }
    return err_result<LrDist>("evidence_combine: unknown combine mode.");
}

}  // namespace core
}  // namespace mispitools
