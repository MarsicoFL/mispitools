#include "linkage.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <set>
#include <utility>
#include <vector>

namespace mispitools {
namespace core {

namespace {

constexpr double kFreqSumTol = 1e-6;

// Per-pair precomputed tables. Scoped to a single linked_pair_joint call.
struct PairTables {
    AlleleIndex K_A = 0;
    AlleleIndex K_B = 0;
    int HK = 0;                 // K_A * K_B  (haplotypes)
    int D = 0;                  // HK * HK    (ordered diplotypes)
    GenotypeIndex G_A = 0;
    GenotypeIndex G_B = 0;

    std::vector<double> prior_d;   // size D, founder ordered-diplotype prior
    std::vector<double> gamete;    // size D * HK, gamete[d * HK + h]
    std::vector<double> imp;       // size HK, population-imputed gamete

    // d -> (pat_hap, mat_hap); hap -> (aA, aB).
    inline int pat_hap(int d) const noexcept { return d / HK; }
    inline int mat_hap(int d) const noexcept { return d % HK; }
    inline AlleleIndex hap_a(int h) const noexcept {
        return static_cast<AlleleIndex>(h / K_B);
    }
    inline AlleleIndex hap_b(int h) const noexcept {
        return static_cast<AlleleIndex>(h % K_B);
    }
};

Result<PairTables> precompute_pair(const Marker& mA, const Marker& mB,
                                   const std::vector<double>& matA,
                                   const std::vector<double>& matB,
                                   double rho) {
    PairTables t;
    t.K_A = mA.n_alleles;
    t.K_B = mB.n_alleles;
    t.HK = static_cast<int>(t.K_A) * static_cast<int>(t.K_B);
    t.D = t.HK * t.HK;
    t.G_A = n_genotypes(t.K_A);
    t.G_B = n_genotypes(t.K_B);

    if (matA.size() != static_cast<std::size_t>(t.K_A) * t.K_A
            || matB.size() != static_cast<std::size_t>(t.K_B) * t.K_B) {
        return err_result<PairTables>(
            "precompute_pair: mutation matrix size != K * K.");
    }

    // Founder ordered-diplotype prior: HWE per locus, linkage equilibrium
    // between loci. P(d) = fA[hpA]*fB[hpB] * fA[hmA]*fB[hmB].
    t.prior_d.assign(static_cast<std::size_t>(t.D), 0.0);
    for (int d = 0; d < t.D; ++d) {
        const int hp = t.pat_hap(d);
        const int hm = t.mat_hap(d);
        const double v =
            mA.freqs[static_cast<std::size_t>(t.hap_a(hp))]
          * mB.freqs[static_cast<std::size_t>(t.hap_b(hp))]
          * mA.freqs[static_cast<std::size_t>(t.hap_a(hm))]
          * mB.freqs[static_cast<std::size_t>(t.hap_b(hm))];
        t.prior_d[static_cast<std::size_t>(d)] = v;
    }

    // Gamete transmission table: gamete[d][ghap] = P(parent diplotype d
    // transmits haplotype ghap) via the two-locus phased kernel.
    t.gamete.assign(static_cast<std::size_t>(t.D) * t.HK, 0.0);
    for (int d = 0; d < t.D; ++d) {
        const int hp = t.pat_hap(d);
        const int hm = t.mat_hap(d);
        const AlleleIndex p1 = t.hap_a(hp), p2 = t.hap_b(hp);
        const AlleleIndex m1 = t.hap_a(hm), m2 = t.hap_b(hm);
        for (int gh = 0; gh < t.HK; ++gh) {
            const AlleleIndex g1 = t.hap_a(gh), g2 = t.hap_b(gh);
            t.gamete[static_cast<std::size_t>(d) * t.HK + gh] =
                trans_prob_MM(p1, p2, m1, m2, g1, g2, rho,
                              matA, matB, t.K_A, t.K_B);
        }
    }

    // Population-imputed gamete (missing parent ~ founder, integrated
    // against the diplotype prior through the *mutated* transmission;
    // mirrors cpt_engine child_dist_one_missing).
    t.imp.assign(static_cast<std::size_t>(t.HK), 0.0);
    for (int d = 0; d < t.D; ++d) {
        const double pr = t.prior_d[static_cast<std::size_t>(d)];
        if (pr == 0.0) continue;
        for (int gh = 0; gh < t.HK; ++gh) {
            t.imp[static_cast<std::size_t>(gh)] +=
                pr * t.gamete[static_cast<std::size_t>(d) * t.HK + gh];
        }
    }

    return ok_result(std::move(t));
}

// child diplotype dc given parent diplotypes; -1 means parent excluded
// (HWE-imputed). pat-origin haplotype of the child comes from the father,
// mat-origin from the mother.
inline double child_cond(const PairTables& t,
                         int d_fa, int d_mo, int dc) noexcept {
    const int cph = t.pat_hap(dc);
    const int cmh = t.mat_hap(dc);
    const double from_fa = (d_fa >= 0)
        ? t.gamete[static_cast<std::size_t>(d_fa) * t.HK + cph]
        : t.imp[static_cast<std::size_t>(cph)];
    const double from_mo = (d_mo >= 0)
        ? t.gamete[static_cast<std::size_t>(d_mo) * t.HK + cmh]
        : t.imp[static_cast<std::size_t>(cmh)];
    return from_fa * from_mo;
}

// Intermediate diplotype joint. state[row * n_members + i] = diplotype
// index of member i, or -1 if i is unassigned/excluded.
struct PartialJoint {
    std::vector<std::int32_t> states_flat;
    std::vector<double> prob;
    MemberIndex n_members = 0;
    std::size_t n_rows() const noexcept { return prob.size(); }
};

PartialJoint build_joint(const Pedigree& p,
                         const PairTables& t,
                         const std::vector<MemberIndex>& founders_use,
                         const std::vector<MemberIndex>& nf_order,
                         const std::vector<std::uint8_t>& excluded_mask) {
    PartialJoint pj;
    pj.n_members = p.n_members;
    const int D = t.D;
    const std::size_t n_mem = static_cast<std::size_t>(p.n_members);

    if (founders_use.empty()) {
        pj.states_flat.clear();
        pj.prob.clear();
    } else {
        std::size_t total = 1;
        for (std::size_t k = 0; k < founders_use.size(); ++k) {
            total *= static_cast<std::size_t>(D);
        }
        pj.states_flat.assign(total * n_mem, std::int32_t(-1));
        pj.prob.assign(total, 1.0);

        std::size_t period = 1;
        for (std::size_t k = 0; k < founders_use.size(); ++k) {
            const MemberIndex fid = founders_use[k];
            for (std::size_t r = 0; r < total; ++r) {
                const int d = static_cast<int>(
                    (r / period) % static_cast<std::size_t>(D));
                pj.states_flat[r * n_mem + static_cast<std::size_t>(fid)] = d;
                pj.prob[r] *= t.prior_d[static_cast<std::size_t>(d)];
            }
            period *= static_cast<std::size_t>(D);
        }
    }

    for (MemberIndex nf : nf_order) {
        if (pj.n_rows() == 0) break;
        const MemberIndex fa = p.father[static_cast<std::size_t>(nf)];
        const MemberIndex mo = p.mother[static_cast<std::size_t>(nf)];
        const bool fa_known = (fa != kNoParent)
            && (excluded_mask[static_cast<std::size_t>(fa)] == 0);
        const bool mo_known = (mo != kNoParent)
            && (excluded_mask[static_cast<std::size_t>(mo)] == 0);

        const std::size_t cur_rows = pj.n_rows();
        const std::size_t new_rows = cur_rows * static_cast<std::size_t>(D);
        std::vector<std::int32_t> new_states(new_rows * n_mem);
        std::vector<double> new_prob(new_rows);

        for (std::size_t r = 0; r < cur_rows; ++r) {
            const double base_prob = pj.prob[r];
            const int d_fa = fa_known
                ? pj.states_flat[r * n_mem + static_cast<std::size_t>(fa)]
                : -1;
            const int d_mo = mo_known
                ? pj.states_flat[r * n_mem + static_cast<std::size_t>(mo)]
                : -1;
            const bool both_excluded = !fa_known && !mo_known;

            for (int dc = 0; dc < D; ++dc) {
                const std::size_t nr =
                    r * static_cast<std::size_t>(D) + static_cast<std::size_t>(dc);
                std::copy(
                    pj.states_flat.begin() + r * n_mem,
                    pj.states_flat.begin() + (r + 1) * n_mem,
                    new_states.begin() + nr * n_mem);
                new_states[nr * n_mem + static_cast<std::size_t>(nf)] = dc;
                const double cond = both_excluded
                    ? 0.0 : child_cond(t, d_fa, d_mo, dc);
                new_prob[nr] = base_prob * cond;
            }
        }

        std::size_t keep = 0;
        for (std::size_t r = 0; r < new_rows; ++r) {
            if (new_prob[r] > 0.0) {
                if (keep != r) {
                    std::copy(
                        new_states.begin() + r * n_mem,
                        new_states.begin() + (r + 1) * n_mem,
                        new_states.begin() + keep * n_mem);
                    new_prob[keep] = new_prob[r];
                }
                ++keep;
            }
        }
        new_states.resize(keep * n_mem);
        new_prob.resize(keep);
        pj.states_flat = std::move(new_states);
        pj.prob = std::move(new_prob);
    }

    return pj;
}

}  // namespace

Result<LinkedJointTable> linked_pair_joint(
        const Pedigree& p,
        const Marker& m_A,
        const Marker& m_B,
        const MutationModel& mu_A,
        const MutationModel& mu_B,
        double rho) {
    if (!(rho >= 0.0 && rho <= 0.5)) {
        return err_result<LinkedJointTable>(
            "linked_pair_joint: rho must be in [0, 0.5].");
    }
    for (const Marker* m : {&m_A, &m_B}) {
        if (m->n_alleles <= 0
                || m->freqs.size() != static_cast<std::size_t>(m->n_alleles)) {
            return err_result<LinkedJointTable>(
                "linked_pair_joint: marker.freqs size != marker.n_alleles.");
        }
        double s = 0.0;
        for (double f : m->freqs) s += f;
        if (std::fabs(s - 1.0) > kFreqSumTol) {
            return err_result<LinkedJointTable>(
                "linked_pair_joint: marker.freqs do not sum to 1.");
        }
    }
    if (p.n_members <= 0
            || static_cast<MemberIndex>(p.father.size()) != p.n_members
            || static_cast<MemberIndex>(p.mother.size()) != p.n_members) {
        return err_result<LinkedJointTable>(
            "linked_pair_joint: pedigree size invariants violated.");
    }
    if (p.poi != kNoParent && (p.poi < 0 || p.poi >= p.n_members)) {
        return err_result<LinkedJointTable>(
            "linked_pair_joint: poi index out of range.");
    }

    auto mmA = build_mutation_matrix(mu_A, m_A.n_alleles,
                                     m_A.numeric_labels, m_A.freqs);
    if (!mmA.ok()) {
        return err_result<LinkedJointTable>(
            std::string("linked_pair_joint (marker A): ") + mmA.error);
    }
    auto mmB = build_mutation_matrix(mu_B, m_B.n_alleles,
                                     m_B.numeric_labels, m_B.freqs);
    if (!mmB.ok()) {
        return err_result<LinkedJointTable>(
            std::string("linked_pair_joint (marker B): ") + mmB.error);
    }

    auto pt = precompute_pair(m_A, m_B, *mmA, *mmB, rho);
    if (!pt.ok()) {
        return err_result<LinkedJointTable>(
            std::string("linked_pair_joint: ") + pt.error);
    }
    const PairTables& t = *pt;

    const std::vector<MemberIndex> all_founders = founders_of(p);
    const std::vector<MemberIndex> all_nonfounders = nonfounders_of(p);
    const MemberIndex poi = p.poi;
    const bool has_poi = (poi != kNoParent);
    const std::size_t n_mem = static_cast<std::size_t>(p.n_members);

    std::vector<std::uint8_t> empty_mask(n_mem, 0);
    std::vector<MemberIndex> nf_order_h1 =
        ancestral_order(p, all_nonfounders, std::vector<MemberIndex>{});
    PartialJoint h1 = build_joint(p, t, all_founders, nf_order_h1, empty_mask);

    PartialJoint h2;
    if (has_poi) {
        std::vector<MemberIndex> founders_h2 = all_founders;
        founders_h2.erase(
            std::remove(founders_h2.begin(), founders_h2.end(), poi),
            founders_h2.end());
        std::vector<std::uint8_t> h2_mask(n_mem, 0);
        h2_mask[static_cast<std::size_t>(poi)] = 1;
        std::vector<MemberIndex> nf_h2_subset;
        nf_h2_subset.reserve(all_nonfounders.size());
        for (MemberIndex nf : all_nonfounders) {
            if (nf != poi) nf_h2_subset.push_back(nf);
        }
        std::vector<MemberIndex> nf_order_h2 =
            ancestral_order(p, nf_h2_subset, std::vector<MemberIndex>{poi});
        h2 = build_joint(p, t, founders_h2, nf_order_h2, h2_mask);
    }

    // Diplotype-level combined table (states over diplotype indices),
    // built exactly as cpt_marker_joint does over genotype indices.
    std::vector<std::int32_t> dip_states;  // n_rows * n_members
    std::vector<double> dip_p_h1;
    std::vector<double> dip_p_h2;

    if (!has_poi) {
        dip_states = h1.states_flat;
        dip_p_h1 = h1.prob;
        dip_p_h2 = h1.prob;
    } else {
        std::map<std::vector<std::int32_t>, double> h2_lookup;
        for (std::size_t r = 0; r < h2.n_rows(); ++r) {
            std::vector<std::int32_t> k;
            k.reserve(n_mem > 0 ? n_mem - 1 : 0);
            for (MemberIndex i = 0; i < p.n_members; ++i) {
                if (i == poi) continue;
                k.push_back(h2.states_flat[r * n_mem + static_cast<std::size_t>(i)]);
            }
            h2_lookup[std::move(k)] = h2.prob[r];
        }

        const std::size_t n_h1 = h1.n_rows();
        std::vector<double> p_h2_for_h1(n_h1, 0.0);
        for (std::size_t r = 0; r < n_h1; ++r) {
            std::vector<std::int32_t> k;
            k.reserve(n_mem > 0 ? n_mem - 1 : 0);
            for (MemberIndex i = 0; i < p.n_members; ++i) {
                if (i == poi) continue;
                k.push_back(h1.states_flat[r * n_mem + static_cast<std::size_t>(i)]);
            }
            auto it = h2_lookup.find(k);
            const double lookup = (it != h2_lookup.end()) ? it->second : 0.0;
            const int d_poi =
                h1.states_flat[r * n_mem + static_cast<std::size_t>(poi)];
            p_h2_for_h1[r] = t.prior_d[static_cast<std::size_t>(d_poi)] * lookup;
        }

        std::set<std::vector<std::int32_t>> h1_full;
        for (std::size_t r = 0; r < n_h1; ++r) {
            h1_full.insert(std::vector<std::int32_t>(
                h1.states_flat.begin() + r * n_mem,
                h1.states_flat.begin() + (r + 1) * n_mem));
        }

        dip_states = h1.states_flat;
        dip_p_h1 = h1.prob;
        dip_p_h2 = std::move(p_h2_for_h1);

        const int D = t.D;
        for (std::size_t r = 0; r < h2.n_rows(); ++r) {
            for (int d_poi = 0; d_poi < D; ++d_poi) {
                std::vector<std::int32_t> row(n_mem);
                std::copy(
                    h2.states_flat.begin() + r * n_mem,
                    h2.states_flat.begin() + (r + 1) * n_mem,
                    row.begin());
                row[static_cast<std::size_t>(poi)] = d_poi;
                if (h1_full.find(row) != h1_full.end()) continue;
                const double p2 =
                    t.prior_d[static_cast<std::size_t>(d_poi)] * h2.prob[r];
                if (p2 <= 0.0) continue;
                dip_states.insert(dip_states.end(), row.begin(), row.end());
                dip_p_h1.push_back(0.0);
                dip_p_h2.push_back(p2);
            }
        }
    }

    // Collapse latent phase: each member's diplotype -> unordered
    // genotype index at A and at B; aggregate identical observable
    // configurations (key = interleaved A0,B0,A1,B1,...).
    std::map<std::vector<GenotypeIndex>, std::pair<double, double>> agg;
    const std::size_t n_rows = dip_p_h1.size();
    for (std::size_t r = 0; r < n_rows; ++r) {
        std::vector<GenotypeIndex> key(2 * n_mem);
        for (MemberIndex i = 0; i < p.n_members; ++i) {
            const int d = dip_states[r * n_mem + static_cast<std::size_t>(i)];
            const int hp = t.pat_hap(d);
            const int hm = t.mat_hap(d);
            AlleleIndex aA1 = t.hap_a(hp), aA2 = t.hap_a(hm);
            AlleleIndex aB1 = t.hap_b(hp), aB2 = t.hap_b(hm);
            if (aA1 > aA2) std::swap(aA1, aA2);
            if (aB1 > aB2) std::swap(aB1, aB2);
            key[2 * static_cast<std::size_t>(i)] =
                pair_to_genotype_index(aA1, aA2);
            key[2 * static_cast<std::size_t>(i) + 1] =
                pair_to_genotype_index(aB1, aB2);
        }
        auto& cell = agg[key];
        cell.first += dip_p_h1[r];
        cell.second += dip_p_h2[r];
    }

    LinkedJointTable out;
    out.n_members = p.n_members;
    out.n_genotypes_a = t.G_A;
    out.n_genotypes_b = t.G_B;

    // std::map iterates in lex order of the interleaved key, which is
    // exactly the requested (A0,B0,A1,B1,...) row ordering.
    for (const auto& kv : agg) {
        const double ph1 = kv.second.first;
        const double ph2 = kv.second.second;
        if (ph1 <= 0.0 && ph2 <= 0.0) continue;
        const std::vector<GenotypeIndex>& key = kv.first;
        for (MemberIndex i = 0; i < p.n_members; ++i) {
            out.states_a.push_back(key[2 * static_cast<std::size_t>(i)]);
            out.states_b.push_back(key[2 * static_cast<std::size_t>(i) + 1]);
        }
        out.p_h1.push_back(ph1);
        out.p_h2.push_back(ph2);
    }

    return ok_result(std::move(out));
}

int linkage_placeholder(int x) {
    return x + 1;
}

}  // namespace core
}  // namespace mispitools
