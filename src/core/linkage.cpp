#include "linkage.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
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

// ---------------------------------------------------------------------------
// Elston-Stewart online variable elimination (F5.5).
//
// The F5.2 engine materialised the dense phase-conscious joint over *all*
// members; an 8-member first-cousin pedigree blows up to ~(G_A*G_B)^8 rows
// (ESCALATION_20260516_200001.md). Here the joint is instead carried only
// over an active frontier: members are introduced ancestor-first and any
// member that is not in the caller-requested keep-set is summed out the
// moment every one of its informative children has been introduced. This
// is exact variable elimination (correct for any acyclic pedigree; loops
// are broken at the R boundary, SCOUT_pedprobr_linkage.md), so the result
// is bit-for-bit the dense joint when the keep-set is the whole pedigree
// (relevant == {} default), and a proper marginal otherwise.
// ---------------------------------------------------------------------------

// Per-build role of each member. H1: every member plays its pedigree role.
// H2: the POI is detached (an independent HWE founder, prior_d) and is
// excluded as a transmitting parent so its children draw the imputed
// gamete -- exactly the F5.2 dense-merge H2 semantics, expressed locally.
struct PeelSpec {
    std::vector<char> founder_like;  // size n_members
    std::vector<char> excluded;      // size n_members (excluded transmitter)
};

// Sparse joint over a dynamic subset of members. `members[c]` is the
// member id held in column `c`; `state[row * width + c]` is its ordered
// diplotype index. `col_of[m]` is the column of member m, or -1.
struct PeelTable {
    std::vector<MemberIndex> members;
    std::vector<int> col_of;
    std::vector<std::int32_t> state;
    std::vector<double> prob;
    std::size_t width() const noexcept { return members.size(); }
    std::size_t n_rows() const noexcept { return prob.size(); }
};

// Add member `m` as a fresh column, expanding the table by its D ordered
// diplotypes and multiplying in its founder / transmission factor. Rows
// with zero probability are compacted away immediately.
void introduce_member(PeelTable& tb, const PairTables& t, const Pedigree& p,
                      const PeelSpec& sp, MemberIndex m) {
    const int D = t.D;
    const std::size_t mi = static_cast<std::size_t>(m);
    const bool fl = sp.founder_like[mi] != 0;
    const MemberIndex fa = fl ? kNoParent : p.father[mi];
    const MemberIndex mo = fl ? kNoParent : p.mother[mi];

    if (tb.members.empty()) {
        // First member: must be founder-like (recursion introduces
        // ancestors first), so the factor is the diplotype prior.
        tb.members.push_back(m);
        tb.col_of[mi] = 0;
        tb.state.assign(static_cast<std::size_t>(D), 0);
        tb.prob.assign(static_cast<std::size_t>(D), 0.0);
        for (int d = 0; d < D; ++d) {
            tb.state[static_cast<std::size_t>(d)] = d;
            tb.prob[static_cast<std::size_t>(d)] =
                t.prior_d[static_cast<std::size_t>(d)];
        }
        return;
    }

    const std::size_t old_w = tb.width();
    const std::size_t new_w = old_w + 1;
    const int new_col = static_cast<int>(old_w);
    tb.members.push_back(m);
    tb.col_of[mi] = new_col;

    const std::size_t cur_rows = tb.n_rows();
    const std::size_t cap_rows = cur_rows * static_cast<std::size_t>(D);
    std::vector<std::int32_t> ns(cap_rows * new_w);
    std::vector<double> np(cap_rows);

    const bool fa_known = (!fl) && (fa != kNoParent)
        && (sp.excluded[static_cast<std::size_t>(fa)] == 0);
    const bool mo_known = (!fl) && (mo != kNoParent)
        && (sp.excluded[static_cast<std::size_t>(mo)] == 0);
    const int col_fa = fa_known ? tb.col_of[static_cast<std::size_t>(fa)] : -1;
    const int col_mo = mo_known ? tb.col_of[static_cast<std::size_t>(mo)] : -1;
    // A non-founder both of whose parents are excluded contributes no
    // pedigree information and is set to zero, matching the F5.2 dense
    // engine's `both_excluded` branch (never reached with only the POI
    // excluded, retained for parity).
    const bool both_excluded = (!fl) && !fa_known && !mo_known;

    std::size_t keep = 0;
    for (std::size_t r = 0; r < cur_rows; ++r) {
        const double base = tb.prob[r];
        const int d_fa = (col_fa >= 0)
            ? tb.state[r * old_w + static_cast<std::size_t>(col_fa)] : -1;
        const int d_mo = (col_mo >= 0)
            ? tb.state[r * old_w + static_cast<std::size_t>(col_mo)] : -1;
        for (int dc = 0; dc < D; ++dc) {
            double factor;
            if (fl) {
                factor = t.prior_d[static_cast<std::size_t>(dc)];
            } else if (both_excluded) {
                factor = 0.0;
            } else {
                factor = child_cond(t, d_fa, d_mo, dc);
            }
            const double pr = base * factor;
            if (pr <= 0.0) continue;
            std::copy(tb.state.begin() + r * old_w,
                      tb.state.begin() + (r + 1) * old_w,
                      ns.begin() + keep * new_w);
            ns[keep * new_w + old_w] = dc;
            np[keep] = pr;
            ++keep;
        }
    }
    ns.resize(keep * new_w);
    np.resize(keep);
    tb.state = std::move(ns);
    tb.prob = std::move(np);
}

// Sum member `m` out of the table: drop its column and aggregate rows
// that become identical on the remaining members.
void eliminate_member(PeelTable& tb, MemberIndex m) {
    const std::size_t mi = static_cast<std::size_t>(m);
    const int col = tb.col_of[mi];
    if (col < 0) return;
    const std::size_t w = tb.width();
    const std::size_t nw = w - 1;
    const std::size_t rows = tb.n_rows();

    std::map<std::vector<std::int32_t>, double> agg;
    std::vector<std::int32_t> key(nw);
    for (std::size_t r = 0; r < rows; ++r) {
        std::size_t k = 0;
        for (std::size_t c = 0; c < w; ++c) {
            if (c == static_cast<std::size_t>(col)) continue;
            key[k++] = tb.state[r * w + c];
        }
        agg[key] += tb.prob[r];
    }

    tb.members.erase(tb.members.begin() + col);
    tb.col_of[mi] = -1;
    for (MemberIndex mm : tb.members) {
        if (tb.col_of[static_cast<std::size_t>(mm)] > col) {
            --tb.col_of[static_cast<std::size_t>(mm)];
        }
    }

    tb.state.assign(agg.size() * nw, 0);
    tb.prob.assign(agg.size(), 0.0);
    std::size_t r = 0;
    for (const auto& kv : agg) {
        std::copy(kv.first.begin(), kv.first.end(),
                  tb.state.begin() + r * nw);
        tb.prob[r] = kv.second;
        ++r;
    }
}

// Build the joint over `keep` members under `sp`, marginalising every
// other member as soon as it is no longer needed as a parent.
PeelTable build_peeled(const Pedigree& p, const PairTables& t,
                       const PeelSpec& sp,
                       const std::vector<char>& keep) {
    const std::size_t n = static_cast<std::size_t>(p.n_members);

    // children[x] = members whose father or mother is x.
    std::vector<std::vector<MemberIndex>> children(n);
    for (MemberIndex c = 0; c < p.n_members; ++c) {
        const MemberIndex f = p.father[static_cast<std::size_t>(c)];
        const MemberIndex m = p.mother[static_cast<std::size_t>(c)];
        if (f != kNoParent) children[static_cast<std::size_t>(f)].push_back(c);
        if (m != kNoParent) children[static_cast<std::size_t>(m)].push_back(c);
    }

    // added_set = keep U ancestors(keep): exactly the members the
    // ancestor-first recursion will introduce. Children outside it are
    // uninformative subtrees (integrate to 1) and never pin a parent.
    std::vector<char> in_added(n, 0);
    {
        std::vector<MemberIndex> stack;
        for (MemberIndex m = 0; m < p.n_members; ++m) {
            if (keep[static_cast<std::size_t>(m)]) {
                in_added[static_cast<std::size_t>(m)] = 1;
                stack.push_back(m);
            }
        }
        while (!stack.empty()) {
            const MemberIndex x = stack.back();
            stack.pop_back();
            if (sp.founder_like[static_cast<std::size_t>(x)]) continue;
            const MemberIndex par[2] = {
                p.father[static_cast<std::size_t>(x)],
                p.mother[static_cast<std::size_t>(x)]
            };
            for (MemberIndex pp : par) {
                if (pp == kNoParent) continue;
                if (sp.excluded[static_cast<std::size_t>(pp)]) continue;
                if (!in_added[static_cast<std::size_t>(pp)]) {
                    in_added[static_cast<std::size_t>(pp)] = 1;
                    stack.push_back(pp);
                }
            }
        }
    }

    PeelTable tb;
    tb.col_of.assign(n, -1);
    std::vector<char> introduced(n, 0);

    auto all_informative_children_in =
        [&](MemberIndex x) -> bool {
        for (MemberIndex c : children[static_cast<std::size_t>(x)]) {
            if (!in_added[static_cast<std::size_t>(c)]) continue;
            if (!introduced[static_cast<std::size_t>(c)]) return false;
        }
        return true;
    };

    auto eliminate_fixpoint = [&]() {
        bool changed = true;
        while (changed) {
            changed = false;
            for (std::size_t i = 0; i < tb.members.size(); ++i) {
                const MemberIndex x = tb.members[i];
                if (keep[static_cast<std::size_t>(x)]) continue;
                if (!all_informative_children_in(x)) continue;
                eliminate_member(tb, x);
                changed = true;
                break;
            }
        }
    };

    // Iterative ancestor-first introduction (explicit stack: the
    // recursion depth is the pedigree generation count, but a stack
    // keeps the kernel free of deep native recursion).
    std::vector<MemberIndex> work;
    for (MemberIndex m = 0; m < p.n_members; ++m) {
        if (keep[static_cast<std::size_t>(m)]) work.push_back(m);
    }
    std::vector<MemberIndex> stack;
    for (auto it = work.rbegin(); it != work.rend(); ++it) {
        stack.push_back(*it);
    }
    while (!stack.empty()) {
        const MemberIndex m = stack.back();
        const std::size_t mi = static_cast<std::size_t>(m);
        if (introduced[mi]) { stack.pop_back(); continue; }
        const bool fl = sp.founder_like[mi] != 0;
        bool deferred = false;
        if (!fl) {
            const MemberIndex par[2] = {
                p.father[mi], p.mother[mi]
            };
            for (MemberIndex pp : par) {
                if (pp == kNoParent) continue;
                if (sp.excluded[static_cast<std::size_t>(pp)]) continue;
                if (!introduced[static_cast<std::size_t>(pp)]) {
                    stack.push_back(pp);
                    deferred = true;
                }
            }
        }
        if (deferred) continue;  // parents pushed; revisit m later
        stack.pop_back();
        introduced[mi] = 1;
        introduce_member(tb, t, p, sp, m);
        eliminate_fixpoint();
    }
    eliminate_fixpoint();
    return tb;
}

// Permute columns into ascending member-id order so two tables built
// under different specs (H1 vs H2) align by a plain row-key lookup.
void canonicalize(PeelTable& tb) {
    const std::size_t w = tb.width();
    if (w <= 1) return;
    std::vector<std::size_t> perm(w);
    for (std::size_t i = 0; i < w; ++i) perm[i] = i;
    std::sort(perm.begin(), perm.end(),
              [&](std::size_t a, std::size_t b) {
                  return tb.members[a] < tb.members[b];
              });
    bool sorted = true;
    for (std::size_t i = 0; i < w; ++i) {
        if (perm[i] != i) { sorted = false; break; }
    }
    if (sorted) return;
    const std::size_t rows = tb.n_rows();
    std::vector<std::int32_t> ns(rows * w);
    for (std::size_t r = 0; r < rows; ++r) {
        for (std::size_t c = 0; c < w; ++c) {
            ns[r * w + c] = tb.state[r * w + perm[c]];
        }
    }
    std::vector<MemberIndex> nm(w);
    for (std::size_t c = 0; c < w; ++c) nm[c] = tb.members[perm[c]];
    tb.state = std::move(ns);
    tb.members = std::move(nm);
    for (std::size_t c = 0; c < w; ++c) {
        tb.col_of[static_cast<std::size_t>(tb.members[c])] =
            static_cast<int>(c);
    }
}

}  // namespace

Result<LinkedJointTable> linked_pair_joint(
        const Pedigree& p,
        const Marker& m_A,
        const Marker& m_B,
        const MutationModel& mu_A,
        const MutationModel& mu_B,
        double rho,
        const std::vector<MemberIndex>& relevant) {
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
    for (MemberIndex r : relevant) {
        if (r < 0 || r >= p.n_members) {
            return err_result<LinkedJointTable>(
                "linked_pair_joint: relevant member index out of range.");
        }
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

    const MemberIndex poi = p.poi;
    const bool has_poi = (poi != kNoParent);
    const std::size_t n_mem = static_cast<std::size_t>(p.n_members);

    // keep-set: an empty `relevant` means "the whole pedigree" -- this
    // reproduces the F5.2 dense joint exactly (nothing is eliminated).
    // Otherwise only the requested members survive; the POI is forced in
    // so H1 and H2 stay alignable.
    std::vector<char> keep(n_mem, 0);
    if (relevant.empty()) {
        std::fill(keep.begin(), keep.end(), char(1));
    } else {
        for (MemberIndex r : relevant) keep[static_cast<std::size_t>(r)] = 1;
        if (has_poi) keep[static_cast<std::size_t>(poi)] = 1;
    }

    PeelSpec sp_h1;
    sp_h1.founder_like.assign(n_mem, 0);
    sp_h1.excluded.assign(n_mem, 0);
    for (MemberIndex m = 0; m < p.n_members; ++m) {
        sp_h1.founder_like[static_cast<std::size_t>(m)] =
            is_founder(p, m) ? 1 : 0;
    }

    PeelTable h1 = build_peeled(p, t, sp_h1, keep);
    canonicalize(h1);

    PeelTable h2;
    if (has_poi) {
        PeelSpec sp_h2 = sp_h1;
        sp_h2.founder_like[static_cast<std::size_t>(poi)] = 1;
        sp_h2.excluded[static_cast<std::size_t>(poi)] = 1;
        h2 = build_peeled(p, t, sp_h2, keep);
        canonicalize(h2);
    }

    // Outer join of H1 / H2 on the (canonical) keep-member configuration.
    // h1 carries the pedigree joint; h2 carries the same support with the
    // POI replaced by an independent HWE founder -- exactly the F5.2
    // dense-merge semantics, now over the marginalised member set.
    const std::size_t w = h1.width();
    std::vector<std::int32_t> combo_state;  // n_rows * w
    std::vector<double> combo_h1;
    std::vector<double> combo_h2;

    if (!has_poi) {
        combo_state = h1.state;
        combo_h1 = h1.prob;
        combo_h2 = h1.prob;
    } else {
        std::map<std::vector<std::int32_t>, double> h2map;
        for (std::size_t r = 0; r < h2.n_rows(); ++r) {
            h2map[std::vector<std::int32_t>(
                h2.state.begin() + r * w,
                h2.state.begin() + (r + 1) * w)] = h2.prob[r];
        }
        combo_state = h1.state;
        combo_h1 = h1.prob;
        combo_h2.assign(h1.n_rows(), 0.0);
        std::vector<std::int32_t> rkey(w);
        for (std::size_t r = 0; r < h1.n_rows(); ++r) {
            std::copy(h1.state.begin() + r * w,
                      h1.state.begin() + (r + 1) * w, rkey.begin());
            auto it = h2map.find(rkey);
            combo_h2[r] = (it != h2map.end()) ? it->second : 0.0;
        }
        // H2-only configurations (compatible under H2, impossible under
        // the pedigree): p_h1 = 0, kept so KL / the LR distribution see
        // the full support.
        std::map<std::vector<std::int32_t>, std::size_t> h1pos;
        for (std::size_t r = 0; r < h1.n_rows(); ++r) {
            h1pos[std::vector<std::int32_t>(
                h1.state.begin() + r * w,
                h1.state.begin() + (r + 1) * w)] = r;
        }
        for (std::size_t r = 0; r < h2.n_rows(); ++r) {
            std::vector<std::int32_t> k(
                h2.state.begin() + r * w,
                h2.state.begin() + (r + 1) * w);
            if (h1pos.find(k) != h1pos.end()) continue;
            if (h2.prob[r] <= 0.0) continue;
            combo_state.insert(combo_state.end(), k.begin(), k.end());
            combo_h1.push_back(0.0);
            combo_h2.push_back(h2.prob[r]);
        }
    }

    // Collapse the latent phase: each surviving member's ordered
    // diplotype -> unordered genotype index at A and at B. Members that
    // were marginalised out hold a constant placeholder (genotype 0) in
    // both markers -- consumers only ever read kept members. Aggregate
    // identical observable configurations; std::map yields the documented
    // (A0,B0,A1,B1,...) lex row order.
    const std::vector<MemberIndex>& cols = h1.members;
    std::map<std::vector<GenotypeIndex>, std::pair<double, double>> agg;
    const std::size_t n_rows = combo_h1.size();
    for (std::size_t r = 0; r < n_rows; ++r) {
        std::vector<GenotypeIndex> key(2 * n_mem, 0);
        for (std::size_t c = 0; c < w; ++c) {
            const MemberIndex i = cols[c];
            const int d = combo_state[r * w + c];
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
        cell.first += combo_h1[r];
        cell.second += combo_h2[r];
    }

    LinkedJointTable out;
    out.n_members = p.n_members;
    out.n_genotypes_a = t.G_A;
    out.n_genotypes_b = t.G_B;
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
