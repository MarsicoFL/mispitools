#include "cpt_engine.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <cstddef>
#include <map>
#include <set>
#include <utility>
#include <vector>

namespace mispitools {
namespace core {

namespace {

constexpr double kFreqSumTol = 1e-6;

// Precomputed per-marker tables. Lifetimes are scoped to a single
// cpt_marker_joint invocation. All arrays are flat row-major (or flat
// column-major as commented) std::vector<double>.
struct Tables {
    AlleleIndex K = 0;
    GenotypeIndex G = 0;
    std::vector<AlleleIndex> a1;  // size G, allele-1 index per genotype (0-based)
    std::vector<AlleleIndex> a2;  // size G, allele-2 index per genotype
    std::vector<double> hwe;      // size G, HWE genotype prior
    std::vector<double> t_trans;  // size G * K, transmission prob row-major
    std::vector<double> child_dist;            // size G^3, both-parents-known table
    std::vector<double> child_dist_one_missing;// size G^2, one-parent-imputed table

    // Indexing: child_dist is laid out so that (gp, gm, gc) maps to
    //   gp + gm * G + gc * G * G
    // This matches the R reference engine's column-major 3-D array; the
    // C++ side just consumes the same flat indices.
    inline std::size_t child_dist_idx(
            GenotypeIndex gp, GenotypeIndex gm, GenotypeIndex gc) const noexcept {
        return static_cast<std::size_t>(gp)
            + static_cast<std::size_t>(gm) * G
            + static_cast<std::size_t>(gc) * G * G;
    }
    inline std::size_t child_dist_one_missing_idx(
            GenotypeIndex g_known, GenotypeIndex g_child) const noexcept {
        return static_cast<std::size_t>(g_known)
            + static_cast<std::size_t>(g_child) * G;
    }
    inline std::size_t t_trans_idx(
            GenotypeIndex g, AlleleIndex a) const noexcept {
        return static_cast<std::size_t>(g) * K + a;
    }
};

Result<Tables> precompute_tables(const std::vector<double>& freqs,
                                 const std::vector<double>& mut_matrix) {
    Tables t;
    t.K = static_cast<AlleleIndex>(freqs.size());
    t.G = n_genotypes(t.K);

    // Genotype enumeration in column-major-lex order (matches R-ref).
    t.a1.resize(static_cast<std::size_t>(t.G));
    t.a2.resize(static_cast<std::size_t>(t.G));
    for (AlleleIndex j = 0; j < t.K; ++j) {
        for (AlleleIndex i = 0; i <= j; ++i) {
            const GenotypeIndex g = pair_to_genotype_index(i, j);
            t.a1[static_cast<std::size_t>(g)] = i;
            t.a2[static_cast<std::size_t>(g)] = j;
        }
    }

    // HWE genotype priors.
    t.hwe.resize(static_cast<std::size_t>(t.G));
    for (GenotypeIndex g = 0; g < t.G; ++g) {
        const AlleleIndex i = t.a1[static_cast<std::size_t>(g)];
        const AlleleIndex j = t.a2[static_cast<std::size_t>(g)];
        t.hwe[static_cast<std::size_t>(g)] = (i == j)
            ? freqs[static_cast<std::size_t>(i)] * freqs[static_cast<std::size_t>(i)]
            : 2.0 * freqs[static_cast<std::size_t>(i)] * freqs[static_cast<std::size_t>(j)];
    }

    // Meiotic transmission table (pre-mutation). T_meiotic[g, a] = P(child
    // receives allele a | parent has genotype g) without mutation.
    std::vector<double> t_meiotic(static_cast<std::size_t>(t.G) * t.K, 0.0);
    for (GenotypeIndex g = 0; g < t.G; ++g) {
        const AlleleIndex i = t.a1[static_cast<std::size_t>(g)];
        const AlleleIndex j = t.a2[static_cast<std::size_t>(g)];
        if (i == j) {
            t_meiotic[static_cast<std::size_t>(g) * t.K + i] = 1.0;
        } else {
            t_meiotic[static_cast<std::size_t>(g) * t.K + i] = 0.5;
            t_meiotic[static_cast<std::size_t>(g) * t.K + j] = 0.5;
        }
    }

    // Mutation matrix sanity. Empty/wrong-size mut_matrix is an invariant
    // violation: build_mutation_matrix should have returned an error.
    if (mut_matrix.size() != static_cast<std::size_t>(t.K) * t.K) {
        return err_result<Tables>(
            "precompute_tables: mut_matrix size != K * K.");
    }

    // T_trans = T_meiotic %*% M_mut (right-multiplied). Matches the R-ref.
    t.t_trans.assign(static_cast<std::size_t>(t.G) * t.K, 0.0);
    for (GenotypeIndex g = 0; g < t.G; ++g) {
        for (AlleleIndex a = 0; a < t.K; ++a) {
            double s = 0.0;
            for (AlleleIndex b = 0; b < t.K; ++b) {
                s += t_meiotic[static_cast<std::size_t>(g) * t.K + b]
                    * mut_matrix[static_cast<std::size_t>(b) * t.K + a];
            }
            t.t_trans[static_cast<std::size_t>(g) * t.K + a] = s;
        }
    }

    // child_dist[gp, gm, gc] = P(child has genotype gc | parents gp, gm).
    t.child_dist.assign(
        static_cast<std::size_t>(t.G) * t.G * t.G, 0.0);
    for (GenotypeIndex gp = 0; gp < t.G; ++gp) {
        for (GenotypeIndex gm = 0; gm < t.G; ++gm) {
            for (GenotypeIndex gc = 0; gc < t.G; ++gc) {
                const AlleleIndex cA = t.a1[static_cast<std::size_t>(gc)];
                const AlleleIndex cB = t.a2[static_cast<std::size_t>(gc)];
                double v;
                if (cA == cB) {
                    v = t.t_trans[t.t_trans_idx(gp, cA)]
                      * t.t_trans[t.t_trans_idx(gm, cA)];
                } else {
                    v = t.t_trans[t.t_trans_idx(gp, cA)]
                        * t.t_trans[t.t_trans_idx(gm, cB)]
                      + t.t_trans[t.t_trans_idx(gp, cB)]
                        * t.t_trans[t.t_trans_idx(gm, cA)];
                }
                t.child_dist[t.child_dist_idx(gp, gm, gc)] = v;
            }
        }
    }

    // child_dist_one_missing[g_known, g_child]: integrate out the missing
    // parent's genotype against the HWE prior. Used when one parent is
    // excluded under H2.
    t.child_dist_one_missing.assign(
        static_cast<std::size_t>(t.G) * t.G, 0.0);
    for (GenotypeIndex g_K = 0; g_K < t.G; ++g_K) {
        for (GenotypeIndex g_C = 0; g_C < t.G; ++g_C) {
            double s = 0.0;
            for (GenotypeIndex g_other = 0; g_other < t.G; ++g_other) {
                s += t.hwe[static_cast<std::size_t>(g_other)]
                    * t.child_dist[t.child_dist_idx(g_K, g_other, g_C)];
            }
            t.child_dist_one_missing[t.child_dist_one_missing_idx(g_K, g_C)] = s;
        }
    }

    return ok_result(std::move(t));
}

// Intermediate joint produced by the founder-enumeration + peeling phase.
// state[row * n_members + i] is the GenotypeIndex of member i in `row`
// when i has been assigned; otherwise it is the sentinel -1.
struct PartialJoint {
    std::vector<GenotypeIndex> states_flat;
    std::vector<double> prob;
    MemberIndex n_members = 0;

    std::size_t n_rows() const noexcept { return prob.size(); }
};

PartialJoint build_joint(
        const Pedigree& p,
        const Tables& t,
        const std::vector<MemberIndex>& founders_use,
        const std::vector<MemberIndex>& nf_order,
        const std::vector<std::uint8_t>& excluded_mask) {
    PartialJoint pj;
    pj.n_members = p.n_members;
    const GenotypeIndex G = t.G;
    const std::size_t n_mem = static_cast<std::size_t>(p.n_members);

    if (founders_use.empty()) {
        // Mirrors R-ref: 0-row joint propagates as 0 rows everywhere.
        pj.states_flat.clear();
        pj.prob.clear();
    } else {
        // Cartesian product over founders_use, in R `expand.grid` order
        // (first founder varies fastest).
        std::size_t total = 1;
        for (std::size_t k = 0; k < founders_use.size(); ++k) {
            total *= static_cast<std::size_t>(G);
        }
        pj.states_flat.assign(total * n_mem, GenotypeIndex(-1));
        pj.prob.assign(total, 1.0);

        std::size_t period = 1;
        for (std::size_t k = 0; k < founders_use.size(); ++k) {
            const MemberIndex fid = founders_use[k];
            for (std::size_t r = 0; r < total; ++r) {
                const GenotypeIndex g =
                    static_cast<GenotypeIndex>((r / period) % static_cast<std::size_t>(G));
                pj.states_flat[r * n_mem + static_cast<std::size_t>(fid)] = g;
                pj.prob[r] *= t.hwe[static_cast<std::size_t>(g)];
            }
            period *= static_cast<std::size_t>(G);
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
        const std::size_t new_rows = cur_rows * static_cast<std::size_t>(G);
        std::vector<GenotypeIndex> new_states(new_rows * n_mem);
        std::vector<double> new_prob(new_rows);

        for (std::size_t r = 0; r < cur_rows; ++r) {
            const double base_prob = pj.prob[r];
            const GenotypeIndex gp = fa_known
                ? pj.states_flat[r * n_mem + static_cast<std::size_t>(fa)]
                : GenotypeIndex(-1);
            const GenotypeIndex gm = mo_known
                ? pj.states_flat[r * n_mem + static_cast<std::size_t>(mo)]
                : GenotypeIndex(-1);

            for (GenotypeIndex gc = 0; gc < G; ++gc) {
                const std::size_t nr =
                    r * static_cast<std::size_t>(G) + static_cast<std::size_t>(gc);
                std::copy(
                    pj.states_flat.begin() + r * n_mem,
                    pj.states_flat.begin() + (r + 1) * n_mem,
                    new_states.begin() + nr * n_mem);
                new_states[nr * n_mem + static_cast<std::size_t>(nf)] = gc;

                double cond;
                if (fa_known && mo_known) {
                    cond = t.child_dist[t.child_dist_idx(gp, gm, gc)];
                } else if (fa_known) {
                    cond = t.child_dist_one_missing[
                        t.child_dist_one_missing_idx(gp, gc)];
                } else if (mo_known) {
                    cond = t.child_dist_one_missing[
                        t.child_dist_one_missing_idx(gm, gc)];
                } else {
                    // Both parents excluded: structurally impossible if the
                    // caller respected ancestral_order. Treat as zero so the
                    // row drops out silently.
                    cond = 0.0;
                }
                new_prob[nr] = base_prob * cond;
            }
        }

        // Drop zero-probability rows in place to keep the row-set sparse.
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

// ---------------------------------------------------------------------------
// Elston-Stewart online variable elimination (F5.6).
//
// `build_joint` above materialises the dense joint over *every* member;
// a large MP pedigree blows up to ~G^n_members rows (same scalability
// limit the linked engine hit in ESCALATION_20260516_200001.md, base G
// instead of D here). This path instead carries the joint only over an
// active frontier: members are introduced ancestor-first and any member
// not in the caller-requested keep-set is summed out the moment every
// one of its informative children has been introduced. It is exact
// variable elimination on an acyclic pedigree, so it returns a proper
// marginal over the kept members. Mirrors linkage.cpp's F5.5 structure
// (single genotype index here, no latent phase to collapse).
// ---------------------------------------------------------------------------

// Per-build role of each member. H1: every member plays its pedigree
// role. H2: the POI is detached (an independent HWE founder) and is
// excluded as a transmitting parent so its children draw the
// HWE-imputed gamete -- exactly the dense-merge H2 semantics, expressed
// locally.
struct PeelSpec {
    std::vector<char> founder_like;  // size n_members
    std::vector<char> excluded;      // size n_members (excluded transmitter)
};

// Sparse joint over a dynamic subset of members. `members[c]` is the
// member id held in column `c`; `state[row * width + c]` is its
// GenotypeIndex. `col_of[m]` is the column of member m, or -1.
struct PeelTable {
    std::vector<MemberIndex> members;
    std::vector<int> col_of;
    std::vector<GenotypeIndex> state;
    std::vector<double> prob;
    std::size_t width() const noexcept { return members.size(); }
    std::size_t n_rows() const noexcept { return prob.size(); }
};

// child genotype gc given parent genotypes; -1 means parent excluded
// (HWE-imputed). Matches build_joint's cond branch exactly.
inline double child_cond_g(const Tables& t,
                           GenotypeIndex gp, GenotypeIndex gm,
                           GenotypeIndex gc) noexcept {
    if (gp >= 0 && gm >= 0) {
        return t.child_dist[t.child_dist_idx(gp, gm, gc)];
    }
    if (gp >= 0) {
        return t.child_dist_one_missing[
            t.child_dist_one_missing_idx(gp, gc)];
    }
    if (gm >= 0) {
        return t.child_dist_one_missing[
            t.child_dist_one_missing_idx(gm, gc)];
    }
    return 0.0;
}

// Add member `m` as a fresh column, expanding by its G genotypes and
// multiplying in its founder / transmission factor. Zero-probability
// rows are compacted away immediately.

enum class FactorKind : std::uint8_t {
    Owned,
    Hwe,
    OneMissing,
    ChildDist,
    Zero
};

struct Factor {
    FactorKind kind = FactorKind::Owned;
    std::vector<MemberIndex> scope;
    std::vector<double> val;
    std::size_t i_fa = 0;
    std::size_t i_mo = 0;
    std::size_t i_ch = 0;
    std::size_t i_kn = 0;
};

constexpr std::size_t kMaxFactorElems = 50000000;

std::size_t pow_g(GenotypeIndex G, std::size_t k) {
    std::size_t n = 1;
    for (std::size_t i = 0; i < k; ++i) {
        if (n > kMaxFactorElems) return kMaxFactorElems + 1;
        n *= static_cast<std::size_t>(G);
    }
    return n;
}

std::size_t slot_of(const std::vector<MemberIndex>& v, MemberIndex m) {
    return static_cast<std::size_t>(
        std::lower_bound(v.begin(), v.end(), m) - v.begin());
}

std::vector<Factor> build_factors(const Pedigree& p, const Tables& t,
                                  const PeelSpec& sp,
                                  const std::vector<char>& in_added) {
    std::vector<Factor> out;
    for (MemberIndex m = 0; m < p.n_members; ++m) {
        const std::size_t mi = static_cast<std::size_t>(m);
        if (!in_added[mi]) continue;

        Factor f;
        if (sp.founder_like[mi] != 0) {
            f.kind = FactorKind::Hwe;
            f.scope = {m};
            out.push_back(std::move(f));
            continue;
        }

        const MemberIndex fa = p.father[mi];
        const MemberIndex mo = p.mother[mi];
        const bool fa_known = (fa != kNoParent)
            && (sp.excluded[static_cast<std::size_t>(fa)] == 0)
            && (in_added[static_cast<std::size_t>(fa)] != 0);
        const bool mo_known = (mo != kNoParent)
            && (sp.excluded[static_cast<std::size_t>(mo)] == 0)
            && (in_added[static_cast<std::size_t>(mo)] != 0);

        if (fa_known && mo_known) {
            f.kind = FactorKind::ChildDist;
            f.scope = {fa, mo, m};
            std::sort(f.scope.begin(), f.scope.end());
            f.i_fa = slot_of(f.scope, fa);
            f.i_mo = slot_of(f.scope, mo);
            f.i_ch = slot_of(f.scope, m);
        } else if (fa_known || mo_known) {
            const MemberIndex kn = fa_known ? fa : mo;
            f.kind = FactorKind::OneMissing;
            f.scope = {kn, m};
            std::sort(f.scope.begin(), f.scope.end());
            f.i_kn = slot_of(f.scope, kn);
            f.i_ch = slot_of(f.scope, m);
        } else {
            f.kind = FactorKind::Zero;
            f.scope = {m};
        }
        out.push_back(std::move(f));
    }
    return out;
}

Result<Factor> join_factors(const std::vector<Factor>& fs,
                            const std::vector<std::size_t>& which,
                            MemberIndex elim,
                            const Tables& t) {
    const GenotypeIndex G = t.G;

    std::vector<MemberIndex> U;
    for (std::size_t i : which) {
        for (MemberIndex m : fs[i].scope) U.push_back(m);
    }
    std::sort(U.begin(), U.end());
    U.erase(std::unique(U.begin(), U.end()), U.end());
    const std::size_t u = U.size();

    Factor out;
    out.kind = FactorKind::Owned;
    for (MemberIndex m : U) {
        if (m != elim) out.scope.push_back(m);
    }

    const std::size_t out_elems = pow_g(G, out.scope.size());
    if (out_elems > kMaxFactorElems) {
        return err_result<Factor>(
            "cpt_marker_joint: the peeled elimination reached a factor over "
            + std::to_string(out.scope.size()) + " members at "
            + std::to_string(static_cast<long long>(G))
            + " genotypes, past the working-memory cap. Name fewer typed "
              "members, or use a marker with fewer alleles.");
    }
    out.val.assign(out_elems, 0.0);

    const std::size_t nf = which.size();
    std::vector<std::vector<std::size_t>> upos(nf);
    for (std::size_t k = 0; k < nf; ++k) {
        for (MemberIndex m : fs[which[k]].scope) {
            upos[k].push_back(slot_of(U, m));
        }
    }
    std::vector<std::vector<std::size_t>> fstride(
        nf, std::vector<std::size_t>(u, 0));
    for (std::size_t k = 0; k < nf; ++k) {
        if (fs[which[k]].kind != FactorKind::Owned) continue;
        std::size_t st = 1;
        for (std::size_t j = fs[which[k]].scope.size(); j-- > 0; ) {
            fstride[k][upos[k][j]] = st;
            st *= static_cast<std::size_t>(G);
        }
    }
    std::vector<std::size_t> ostride(u, 0);
    {
        std::size_t st = 1;
        for (std::size_t j = out.scope.size(); j-- > 0; ) {
            ostride[slot_of(U, out.scope[j])] = st;
            st *= static_cast<std::size_t>(G);
        }
    }

    std::vector<GenotypeIndex> asg(u, 0);
    std::vector<std::size_t> idxf(nf, 0);
    std::size_t idxo = 0;
    const std::size_t gz = static_cast<std::size_t>(G);

    for (;;) {
        double prod = 1.0;
        for (std::size_t k = 0; k < nf && prod != 0.0; ++k) {
            const Factor& f = fs[which[k]];
            switch (f.kind) {
                case FactorKind::Owned:
                    prod *= f.val[idxf[k]];
                    break;
                case FactorKind::Hwe:
                    prod *= t.hwe[static_cast<std::size_t>(
                        asg[upos[k][0]])];
                    break;
                case FactorKind::OneMissing:
                    prod *= t.child_dist_one_missing[
                        t.child_dist_one_missing_idx(asg[upos[k][f.i_kn]],
                                                     asg[upos[k][f.i_ch]])];
                    break;
                case FactorKind::ChildDist:
                    prod *= t.child_dist[
                        t.child_dist_idx(asg[upos[k][f.i_fa]],
                                         asg[upos[k][f.i_mo]],
                                         asg[upos[k][f.i_ch]])];
                    break;
                case FactorKind::Zero:
                    prod = 0.0;
                    break;
            }
        }
        if (prod != 0.0) out.val[idxo] += prod;

        std::size_t d = u;
        bool carry = true;
        while (carry && d > 0) {
            --d;
            asg[d] = static_cast<GenotypeIndex>(asg[d] + 1);
            for (std::size_t k = 0; k < nf; ++k) idxf[k] += fstride[k][d];
            idxo += ostride[d];
            if (asg[d] < G) {
                carry = false;
            } else {
                asg[d] = 0;
                for (std::size_t k = 0; k < nf; ++k) {
                    idxf[k] -= fstride[k][d] * gz;
                }
                idxo -= ostride[d] * gz;
            }
        }
        if (carry) break;
    }

    return ok_result(std::move(out));
}

Result<PeelTable> build_peeled_factored(const Pedigree& p, const Tables& t,
                                        const PeelSpec& sp,
                                        const std::vector<char>& keep) {
    const std::size_t n = static_cast<std::size_t>(p.n_members);
    const GenotypeIndex G = t.G;

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
            const std::size_t xi = static_cast<std::size_t>(x);
            if (sp.founder_like[xi]) continue;
            const MemberIndex par[2] = { p.father[xi], p.mother[xi] };
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

    std::vector<Factor> fs = build_factors(p, t, sp, in_added);
    std::vector<char> alive(fs.size(), 1);

    for (;;) {
        MemberIndex best = kNoParent;
        std::size_t best_cost = 0;
        std::vector<std::size_t> best_which;
        for (MemberIndex x = 0; x < p.n_members; ++x) {
            const std::size_t xi = static_cast<std::size_t>(x);
            if (!in_added[xi] || keep[xi]) continue;
            std::vector<std::size_t> which;
            std::vector<MemberIndex> U;
            for (std::size_t i = 0; i < fs.size(); ++i) {
                if (!alive[i]) continue;
                if (std::find(fs[i].scope.begin(), fs[i].scope.end(), x)
                        == fs[i].scope.end()) {
                    continue;
                }
                which.push_back(i);
                for (MemberIndex m : fs[i].scope) U.push_back(m);
            }
            if (which.empty()) continue;
            std::sort(U.begin(), U.end());
            U.erase(std::unique(U.begin(), U.end()), U.end());
            const std::size_t cost = pow_g(G, U.size() - 1);
            if (best == kNoParent || cost < best_cost) {
                best = x;
                best_cost = cost;
                best_which = std::move(which);
            }
        }
        if (best == kNoParent) break;

        auto joined = join_factors(fs, best_which, best, t);
        if (!joined.ok()) return err_result<PeelTable>(joined.error);
        for (std::size_t i : best_which) alive[i] = 0;
        fs.push_back(std::move(*joined));
        alive.push_back(1);
    }

    std::vector<std::size_t> rest;
    for (std::size_t i = 0; i < fs.size(); ++i) {
        if (alive[i]) rest.push_back(i);
    }

    PeelTable tb;
    tb.col_of.assign(n, -1);
    if (rest.empty()) return ok_result(std::move(tb));

    auto product = join_factors(fs, rest, kNoParent, t);
    if (!product.ok()) return err_result<PeelTable>(product.error);
    const Factor& f = *product;

    const std::size_t w = f.scope.size();
    tb.members = f.scope;
    for (std::size_t c = 0; c < w; ++c) {
        tb.col_of[static_cast<std::size_t>(f.scope[c])] =
            static_cast<int>(c);
    }
    std::vector<GenotypeIndex> asg(w, 0);
    for (std::size_t flat = 0; flat < f.val.size(); ++flat) {
        if (f.val[flat] > 0.0) {
            tb.state.insert(tb.state.end(), asg.begin(), asg.end());
            tb.prob.push_back(f.val[flat]);
        }
        std::size_t d = w;
        while (d > 0) {
            --d;
            asg[d] = static_cast<GenotypeIndex>(asg[d] + 1);
            if (asg[d] < G) break;
            asg[d] = 0;
        }
    }
    return ok_result(std::move(tb));
}

// F5.6 peeled assembly: marginal joint over the keep-set, H1 / H2 outer
// joined exactly like the dense merge but over the kept members only.
Result<JointTable> cpt_marker_joint_peeled(
        const Pedigree& p, const Tables& t,
        const std::vector<MemberIndex>& relevant) {
    const std::size_t n_mem = static_cast<std::size_t>(p.n_members);
    const MemberIndex poi = p.poi;
    const bool has_poi = (poi != kNoParent);

    std::vector<char> keep(n_mem, 0);
    for (MemberIndex r : relevant) keep[static_cast<std::size_t>(r)] = 1;
    if (has_poi) keep[static_cast<std::size_t>(poi)] = 1;

    PeelSpec sp_h1;
    sp_h1.founder_like.assign(n_mem, 0);
    sp_h1.excluded.assign(n_mem, 0);
    for (MemberIndex m = 0; m < p.n_members; ++m) {
        sp_h1.founder_like[static_cast<std::size_t>(m)] =
            is_founder(p, m) ? 1 : 0;
    }

    auto r_h1 = build_peeled_factored(p, t, sp_h1, keep);
    if (!r_h1.ok()) return err_result<JointTable>(r_h1.error);
    PeelTable h1 = std::move(*r_h1);

    PeelTable h2;
    if (has_poi) {
        PeelSpec sp_h2 = sp_h1;
        sp_h2.founder_like[static_cast<std::size_t>(poi)] = 1;
        sp_h2.excluded[static_cast<std::size_t>(poi)] = 1;
        auto r_h2 = build_peeled_factored(p, t, sp_h2, keep);
        if (!r_h2.ok()) return err_result<JointTable>(r_h2.error);
        h2 = std::move(*r_h2);
    }

    const std::size_t w = h1.width();
    std::vector<GenotypeIndex> combo_state;  // n_rows * w
    std::vector<double> combo_h1;
    std::vector<double> combo_h2;

    if (!has_poi) {
        combo_state = h1.state;
        combo_h1 = h1.prob;
        combo_h2 = h1.prob;
    } else {
        std::map<std::vector<GenotypeIndex>, double> h2map;
        for (std::size_t r = 0; r < h2.n_rows(); ++r) {
            h2map[std::vector<GenotypeIndex>(
                h2.state.begin() + r * w,
                h2.state.begin() + (r + 1) * w)] = h2.prob[r];
        }
        combo_state = h1.state;
        combo_h1 = h1.prob;
        combo_h2.assign(h1.n_rows(), 0.0);
        std::vector<GenotypeIndex> rkey(w);
        for (std::size_t r = 0; r < h1.n_rows(); ++r) {
            std::copy(h1.state.begin() + r * w,
                      h1.state.begin() + (r + 1) * w, rkey.begin());
            auto it = h2map.find(rkey);
            combo_h2[r] = (it != h2map.end()) ? it->second : 0.0;
        }
        // H2-only configurations (compatible under H2, impossible under
        // the pedigree): p_h1 = 0, kept so KL / the LR distribution see
        // the full support. Mirrors the dense engine's phase-2 append.
        std::map<std::vector<GenotypeIndex>, std::size_t> h1pos;
        for (std::size_t r = 0; r < h1.n_rows(); ++r) {
            h1pos[std::vector<GenotypeIndex>(
                h1.state.begin() + r * w,
                h1.state.begin() + (r + 1) * w)] = r;
        }
        for (std::size_t r = 0; r < h2.n_rows(); ++r) {
            std::vector<GenotypeIndex> k(
                h2.state.begin() + r * w,
                h2.state.begin() + (r + 1) * w);
            if (h1pos.find(k) != h1pos.end()) continue;
            if (h2.prob[r] <= 0.0) continue;
            combo_state.insert(combo_state.end(), k.begin(), k.end());
            combo_h1.push_back(0.0);
            combo_h2.push_back(h2.prob[r]);
        }
    }

    // Expand the kept-member configuration to full-pedigree states with
    // a placeholder genotype 0 for marginalised members (consumers only
    // read kept members), then aggregate identical configurations. The
    // std::map yields the documented lex row order over the member-state
    // tuple, matching the dense engine's final sort.
    const std::vector<MemberIndex>& cols = h1.members;
    std::map<std::vector<GenotypeIndex>, std::pair<double, double>> agg;
    const std::size_t n_rows = combo_h1.size();
    for (std::size_t r = 0; r < n_rows; ++r) {
        std::vector<GenotypeIndex> key(n_mem, 0);
        for (std::size_t c = 0; c < w; ++c) {
            key[static_cast<std::size_t>(cols[c])] =
                combo_state[r * w + c];
        }
        auto& cell = agg[key];
        cell.first += combo_h1[r];
        cell.second += combo_h2[r];
    }

    JointTable out;
    out.n_members = p.n_members;
    out.n_genotypes = t.G;
    for (const auto& kv : agg) {
        const double ph1 = kv.second.first;
        const double ph2 = kv.second.second;
        if (ph1 <= 0.0 && ph2 <= 0.0) continue;
        out.states_flat.insert(out.states_flat.end(),
                               kv.first.begin(), kv.first.end());
        out.p_h1.push_back(ph1);
        out.p_h2.push_back(ph2);
    }
    return ok_result(std::move(out));
}

}  // namespace

Result<JointTable> cpt_marker_joint(
        const Pedigree& p,
        const Marker& marker,
        const MutationModel& mut,
        const std::vector<MemberIndex>& relevant) {
    // F3.4: factored out so per_marker_kl_batch() can supply a cached
    // mutation matrix without rebuilding it per marker. This wrapper
    // keeps the historical single-shot entry point: build the K x K
    // matrix from the MutationModel parameters, then delegate.
    auto mm = build_mutation_matrix(mut, marker.n_alleles,
                                    marker.numeric_labels, marker.freqs);
    if (!mm.ok()) {
        return err_result<JointTable>(
            std::string("cpt_marker_joint: ") + mm.error);
    }
    return cpt_marker_joint_with_mm(p, marker, *mm, relevant);
}

Result<JointTable> cpt_marker_joint_with_mm(
        const Pedigree& p,
        const Marker& marker,
        const std::vector<double>& mut_matrix,
        const std::vector<MemberIndex>& relevant) {
    if (marker.n_alleles <= 0
            || marker.freqs.size() != static_cast<std::size_t>(marker.n_alleles)) {
        return err_result<JointTable>(
            "cpt_marker_joint: marker.freqs size != marker.n_alleles.");
    }
    {
        double s = 0.0;
        for (double f : marker.freqs) s += f;
        if (std::fabs(s - 1.0) > kFreqSumTol) {
            return err_result<JointTable>(
                "cpt_marker_joint: marker.freqs do not sum to 1.");
        }
    }
    if (p.n_members <= 0
            || static_cast<MemberIndex>(p.father.size()) != p.n_members
            || static_cast<MemberIndex>(p.mother.size()) != p.n_members) {
        return err_result<JointTable>(
            "cpt_marker_joint: pedigree size invariants violated.");
    }
    if (p.poi != kNoParent && (p.poi < 0 || p.poi >= p.n_members)) {
        return err_result<JointTable>(
            "cpt_marker_joint: poi index out of range.");
    }

    auto pt = precompute_tables(marker.freqs, mut_matrix);
    if (!pt.ok()) {
        return err_result<JointTable>(
            std::string("cpt_marker_joint: ") + pt.error);
    }
    const Tables& t = *pt;

    for (MemberIndex r : relevant) {
        if (r < 0 || r >= p.n_members) {
            return err_result<JointTable>(
                "cpt_marker_joint: relevant member index out of range.");
        }
    }
    // Empty `relevant` => dense joint over every member (the F2 engine,
    // bit-for-bit). A strict subset routes the F5.6 peeled engine, which
    // marginalises every other member by online variable elimination.
    if (!relevant.empty()) {
        return cpt_marker_joint_peeled(p, t, relevant);
    }

    const std::vector<MemberIndex> all_founders = founders_of(p);
    const std::vector<MemberIndex> all_nonfounders = nonfounders_of(p);
    const MemberIndex poi = p.poi;
    const bool has_poi = (poi != kNoParent);

    // H1 build: no member excluded.
    std::vector<std::uint8_t> empty_mask(static_cast<std::size_t>(p.n_members), 0);
    std::vector<MemberIndex> nf_order_h1 =
        ancestral_order(p, all_nonfounders, std::vector<MemberIndex>{});
    PartialJoint h1 = build_joint(p, t, all_founders, nf_order_h1, empty_mask);

    // H2 build: POI excluded (HWE-marginalised).
    PartialJoint h2;
    std::vector<MemberIndex> founders_h2 = all_founders;
    std::vector<std::uint8_t> h2_mask(static_cast<std::size_t>(p.n_members), 0);
    if (has_poi) {
        h2_mask[static_cast<std::size_t>(poi)] = 1;
        founders_h2.erase(
            std::remove(founders_h2.begin(), founders_h2.end(), poi),
            founders_h2.end());
        std::vector<MemberIndex> nf_h2_subset;
        nf_h2_subset.reserve(all_nonfounders.size());
        for (MemberIndex nf : all_nonfounders) {
            if (nf != poi) nf_h2_subset.push_back(nf);
        }
        std::vector<MemberIndex> nf_order_h2 =
            ancestral_order(p, nf_h2_subset, std::vector<MemberIndex>{poi});
        h2 = build_joint(p, t, founders_h2, nf_order_h2, h2_mask);
    }

    const std::size_t n_mem = static_cast<std::size_t>(p.n_members);
    const GenotypeIndex G = t.G;

    JointTable out;
    out.n_members = p.n_members;
    out.n_genotypes = G;

    if (!has_poi) {
        // Degenerate fallback. Not reached when the binding layer always
        // resolves a POI (which it does for marker_model inputs).
        out.states_flat = h1.states_flat;
        out.p_h1 = h1.prob;
        out.p_h2 = h1.prob;
    } else {
        // h2 lookup keyed on the "others" tuple (members minus the POI).
        std::map<std::vector<GenotypeIndex>, double> h2_lookup;
        for (std::size_t r = 0; r < h2.n_rows(); ++r) {
            std::vector<GenotypeIndex> k;
            k.reserve(n_mem > 0 ? n_mem - 1 : 0);
            for (MemberIndex i = 0; i < p.n_members; ++i) {
                if (i == poi) continue;
                k.push_back(h2.states_flat[r * n_mem + static_cast<std::size_t>(i)]);
            }
            h2_lookup[std::move(k)] = h2.prob[r];
        }

        // Phase 1: for each H1 row, P_H2 = hwe[g_poi] * h2_lookup(others).
        const std::size_t n_h1 = h1.n_rows();
        std::vector<double> p_h2_for_h1(n_h1, 0.0);
        for (std::size_t r = 0; r < n_h1; ++r) {
            std::vector<GenotypeIndex> k;
            k.reserve(n_mem > 0 ? n_mem - 1 : 0);
            for (MemberIndex i = 0; i < p.n_members; ++i) {
                if (i == poi) continue;
                k.push_back(h1.states_flat[r * n_mem + static_cast<std::size_t>(i)]);
            }
            auto it = h2_lookup.find(k);
            const double lookup = (it != h2_lookup.end()) ? it->second : 0.0;
            const GenotypeIndex g_poi =
                h1.states_flat[r * n_mem + static_cast<std::size_t>(poi)];
            p_h2_for_h1[r] = t.hwe[static_cast<std::size_t>(g_poi)] * lookup;
        }

        // Phase 2: expand H2 over G choices of POI's genotype; keep rows
        // that are not already in H1's full-tuple support.
        std::set<std::vector<GenotypeIndex>> h1_full;
        for (std::size_t r = 0; r < n_h1; ++r) {
            h1_full.insert(std::vector<GenotypeIndex>(
                h1.states_flat.begin() + r * n_mem,
                h1.states_flat.begin() + (r + 1) * n_mem));
        }

        out.states_flat = h1.states_flat;
        out.p_h1 = h1.prob;
        out.p_h2 = std::move(p_h2_for_h1);

        const std::size_t n_h2 = h2.n_rows();
        for (std::size_t r = 0; r < n_h2; ++r) {
            for (GenotypeIndex g_poi = 0; g_poi < G; ++g_poi) {
                std::vector<GenotypeIndex> row(n_mem);
                std::copy(
                    h2.states_flat.begin() + r * n_mem,
                    h2.states_flat.begin() + (r + 1) * n_mem,
                    row.begin());
                row[static_cast<std::size_t>(poi)] = g_poi;
                if (h1_full.find(row) != h1_full.end()) continue;
                const double p2 =
                    t.hwe[static_cast<std::size_t>(g_poi)] * h2.prob[r];
                out.states_flat.insert(out.states_flat.end(),
                                       row.begin(), row.end());
                out.p_h1.push_back(0.0);
                out.p_h2.push_back(p2);
            }
        }
    }

    // Drop rows with both probabilities zero (mirrors R-ref `keep` filter).
    {
        const std::size_t n = out.p_h1.size();
        std::size_t keep = 0;
        for (std::size_t r = 0; r < n; ++r) {
            if (out.p_h1[r] > 0.0 || out.p_h2[r] > 0.0) {
                if (keep != r) {
                    std::copy(
                        out.states_flat.begin() + r * n_mem,
                        out.states_flat.begin() + (r + 1) * n_mem,
                        out.states_flat.begin() + keep * n_mem);
                    out.p_h1[keep] = out.p_h1[r];
                    out.p_h2[keep] = out.p_h2[r];
                }
                ++keep;
            }
        }
        out.states_flat.resize(keep * n_mem);
        out.p_h1.resize(keep);
        out.p_h2.resize(keep);
    }

    // Lex-sort rows on (state[0], state[1], ..., state[n_members-1]).
    {
        const std::size_t n = out.p_h1.size();
        std::vector<std::size_t> order(n);
        for (std::size_t i = 0; i < n; ++i) order[i] = i;
        std::sort(order.begin(), order.end(),
            [&](std::size_t a, std::size_t b) {
                for (MemberIndex i = 0; i < p.n_members; ++i) {
                    const GenotypeIndex ga = out.states_flat[
                        a * n_mem + static_cast<std::size_t>(i)];
                    const GenotypeIndex gb = out.states_flat[
                        b * n_mem + static_cast<std::size_t>(i)];
                    if (ga != gb) return ga < gb;
                }
                return false;
            });

        std::vector<GenotypeIndex> new_states(out.states_flat.size());
        std::vector<double> new_p_h1(n), new_p_h2(n);
        for (std::size_t i = 0; i < n; ++i) {
            std::copy(
                out.states_flat.begin() + order[i] * n_mem,
                out.states_flat.begin() + (order[i] + 1) * n_mem,
                new_states.begin() + i * n_mem);
            new_p_h1[i] = out.p_h1[order[i]];
            new_p_h2[i] = out.p_h2[order[i]];
        }
        out.states_flat = std::move(new_states);
        out.p_h1 = std::move(new_p_h1);
        out.p_h2 = std::move(new_p_h2);
    }

    return ok_result(std::move(out));
}

int cpt_engine_placeholder(int x) {
    return x + 1;
}

}  // namespace core
}  // namespace mispitools
