#include "kl_engine.h"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <map>
#include <string>
#include <utility>

#include "mutation_models.h"

namespace mispitools {
namespace core {

int kl_engine_placeholder(int x) {
    return x + 1;
}

Result<PerMarkerKL> per_marker_kl(const JointTable& joint) {
    if (joint.p_h1.size() != joint.p_h2.size()) {
        return err_result<PerMarkerKL>(
            "per_marker_kl: P_H1 and P_H2 vectors have mismatched lengths.");
    }

    PerMarkerKL out;
    const std::size_t n = joint.p_h1.size();

    // Match R-ref bit-for-bit: accumulate finite contributions in row order
    // (cpt_marker_joint already lex-sorts), and only collapse to ±Inf at the
    // end if any absolute-continuity violation was seen. Doing it in two
    // stages instead of a single sum-of-Inf-plus-finite avoids relying on
    // IEEE Inf+x semantics (which would also work, but the explicit flag is
    // clearer at the abstraction boundary).
    bool e_h1_inf = false;
    bool e_h2_neg_inf = false;

    for (std::size_t r = 0; r < n; ++r) {
        const double p1 = joint.p_h1[r];
        const double p2 = joint.p_h2[r];
        if (p1 < 0.0 || p2 < 0.0) {
            return err_result<PerMarkerKL>(
                "per_marker_kl: negative probability at row "
                + std::to_string(r) + ".");
        }
        const bool pos1 = p1 > 0.0;
        const bool pos2 = p2 > 0.0;

        if (pos1 && pos2) {
            // R-ref: log10(P1) - log10(P2), then sum P_i * (...) left-to-right.
            const double term = std::log10(p1) - std::log10(p2);
            out.e_log10_lr_h1 += p1 * term;
            out.e_log10_lr_h2 += p2 * term;
        } else if (pos1) {
            // P1 > 0, P2 = 0: H2 lacks the support of H1.
            e_h1_inf = true;
            out.abs_cont_violations_h2 += 1;
            out.mass_violations_h2 += p1;
        } else if (pos2) {
            // P1 = 0, P2 > 0: H1 lacks the support of H2.
            e_h2_neg_inf = true;
            out.abs_cont_violations_h1 += 1;
            out.mass_violations_h1 += p2;
        }
        // (!pos1 && !pos2): cpt_marker_joint already filters these.
    }

    // R uses log(10) (natural log); std::log(10.0) is bit-equivalent.
    const double ln10 = std::log(10.0);
    const double inf  = std::numeric_limits<double>::infinity();

    if (e_h1_inf) {
        out.e_log10_lr_h1 = inf;
        out.kl_h1_to_h2 = inf;
    } else {
        out.kl_h1_to_h2 = out.e_log10_lr_h1 * ln10;
    }

    if (e_h2_neg_inf) {
        out.e_log10_lr_h2 = -inf;
        out.kl_h2_to_h1 = inf;
    } else {
        out.kl_h2_to_h1 = -out.e_log10_lr_h2 * ln10;
    }

    return ok_result(out);
}

Result<std::vector<PerMarkerKL>> per_marker_kl_batch(
        const std::vector<JointTable>& joints) {
    std::vector<PerMarkerKL> out;
    out.reserve(joints.size());

    // One sparse pass per marker, reusing the F3.1 scalar kernel. No
    // pedigree peeling and no mutation-matrix construction happen here:
    // the joints arrive fully built. The first failing marker short-
    // circuits so the caller never sees a partially populated vector.
    for (std::size_t i = 0; i < joints.size(); ++i) {
        auto kl = per_marker_kl(joints[i]);
        if (!kl.ok()) {
            return err_result<std::vector<PerMarkerKL>>(
                "per_marker_kl_batch: marker " + std::to_string(i)
                + ": " + kl.error);
        }
        out.push_back(*kl);
    }

    return ok_result(std::move(out));
}

namespace {

// Serialise the mutation matrix signature into a bit-exact byte string.
// Two markers hit the same cache slot iff the produced K x K matrix is
// identical: (kind, K, rate, range) always participate; numeric_labels
// participate only for Stepwise (where build_mutation_matrix() reads
// them) so an Equal-rate profile with stale label arrays still hits;
// afreq participates only for Asymmetric (Dawid is afreq-dependent).
std::string make_mutation_matrix_key(
        const MutationModel& mut,
        AlleleIndex K,
        const std::vector<double>& labels,
        const std::vector<double>& afreq) {
    std::string s;
    const std::size_t label_bytes =
        (mut.kind == MutationKind::Stepwise)
            ? labels.size() * sizeof(double) : 0;
    const std::size_t afreq_bytes =
        (mut.kind == MutationKind::Asymmetric)
            ? afreq.size() * sizeof(double) : 0;
    s.reserve(1 + sizeof(std::int32_t) + 2 * sizeof(double)
              + label_bytes + afreq_bytes);

    const std::uint8_t kind_byte = static_cast<std::uint8_t>(mut.kind);
    s.append(reinterpret_cast<const char*>(&kind_byte), 1);

    const std::int32_t k_val = K;
    s.append(reinterpret_cast<const char*>(&k_val), sizeof(k_val));

    const double rate_val  = mut.rate;
    const double range_val = mut.range;
    s.append(reinterpret_cast<const char*>(&rate_val),  sizeof(rate_val));
    s.append(reinterpret_cast<const char*>(&range_val), sizeof(range_val));

    if (label_bytes > 0) {
        s.append(reinterpret_cast<const char*>(labels.data()), label_bytes);
    }
    if (afreq_bytes > 0) {
        s.append(reinterpret_cast<const char*>(afreq.data()), afreq_bytes);
    }
    return s;
}

}  // namespace

Result<PerMarkerKLBatch> per_marker_kl_batch(
        const Pedigree& ped,
        const std::vector<Marker>& markers,
        const std::vector<MutationModel>& mutations) {
    if (markers.size() != mutations.size()) {
        return err_result<PerMarkerKLBatch>(
            "per_marker_kl_batch: markers.size() != mutations.size().");
    }

    PerMarkerKLBatch out;

    // std::map keyed on the bit-exact serialisation built above. We use
    // map instead of unordered_map to avoid wiring a custom hash; the
    // expected key count (#distinct mutation signatures across a typical
    // STR panel) is small (<= a few), so the log-factor is negligible.
    std::map<std::string, std::vector<double>> mm_cache;

    // Phase 1: build one sparse joint per marker, reusing the K x K
    // mutation matrix across markers whose (kind, K, rate, range, labels)
    // signature collides. Building from a cached matrix is bit-for-bit
    // identical to rebuilding it: cpt_marker_joint_with_mm() consumes the
    // same row-major matrix regardless of provenance.
    std::vector<JointTable> joints;
    joints.reserve(markers.size());

    for (std::size_t i = 0; i < markers.size(); ++i) {
        const Marker& m = markers[i];
        const MutationModel& mut = mutations[i];

        std::string key = make_mutation_matrix_key(mut, m.n_alleles,
                                                   m.numeric_labels,
                                                   m.freqs);
        auto it = mm_cache.find(key);
        const std::vector<double>* mm_ptr = nullptr;
        if (it != mm_cache.end()) {
            mm_ptr = &it->second;
            out.mutation_matrix_cache_hits += 1;
        } else {
            auto mm = build_mutation_matrix(mut, m.n_alleles,
                                            m.numeric_labels, m.freqs);
            if (!mm.ok()) {
                return err_result<PerMarkerKLBatch>(
                    "per_marker_kl_batch: marker " + std::to_string(i)
                    + ": " + mm.error);
            }
            auto ins = mm_cache.emplace(std::move(key), std::move(*mm));
            mm_ptr = &ins.first->second;
            out.mutation_matrix_cache_misses += 1;
        }

        auto joint = cpt_marker_joint_with_mm(ped, m, *mm_ptr);
        if (!joint.ok()) {
            return err_result<PerMarkerKLBatch>(
                "per_marker_kl_batch: marker " + std::to_string(i)
                + ": " + joint.error);
        }
        joints.push_back(std::move(*joint));
    }

    // Phase 2: delegate the KL pass to the F3.4a pure batch primitive so
    // the cache path and the pre-built-joint path share one kernel. Errors
    // are already marker-indexed by per_marker_kl_batch(joints).
    auto kls = per_marker_kl_batch(joints);
    if (!kls.ok()) {
        return err_result<PerMarkerKLBatch>(kls.error);
    }
    out.entries = std::move(*kls);

    return ok_result(std::move(out));
}

}  // namespace core
}  // namespace mispitools
