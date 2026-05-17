#include "cpt_engine.h"

#include <algorithm>
#include <cmath>
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
void introduce_member_g(PeelTable& tb, const Tables& t, const Pedigree& p,
                        const PeelSpec& sp, MemberIndex m) {
    const GenotypeIndex G = t.G;
    const std::size_t mi = static_cast<std::size_t>(m);
    const bool fl = sp.founder_like[mi] != 0;

    if (tb.members.empty()) {
        // First member: ancestor-first introduction guarantees it is
        // founder-like, so the factor is the HWE genotype prior.
        tb.members.push_back(m);
        tb.col_of[mi] = 0;
        tb.state.assign(static_cast<std::size_t>(G), 0);
        tb.prob.assign(static_cast<std::size_t>(G), 0.0);
        for (GenotypeIndex g = 0; g < G; ++g) {
            tb.state[static_cast<std::size_t>(g)] = g;
            tb.prob[static_cast<std::size_t>(g)] =
                t.hwe[static_cast<std::size_t>(g)];
        }
        return;
    }

    const std::size_t old_w = tb.width();
    const std::size_t new_w = old_w + 1;
    const int new_col = static_cast<int>(old_w);
    tb.members.push_back(m);
    tb.col_of[mi] = new_col;

    const MemberIndex fa = fl ? kNoParent : p.father[mi];
    const MemberIndex mo = fl ? kNoParent : p.mother[mi];
    const bool fa_known = (!fl) && (fa != kNoParent)
        && (sp.excluded[static_cast<std::size_t>(fa)] == 0);
    const bool mo_known = (!fl) && (mo != kNoParent)
        && (sp.excluded[static_cast<std::size_t>(mo)] == 0);
    const int col_fa = fa_known
        ? tb.col_of[static_cast<std::size_t>(fa)] : -1;
    const int col_mo = mo_known
        ? tb.col_of[static_cast<std::size_t>(mo)] : -1;
    // A non-founder with both parents excluded carries no pedigree
    // information; set to zero, matching build_joint's both-missing
    // branch (not reached with only the POI excluded; kept for parity).
    const bool both_excluded = (!fl) && !fa_known && !mo_known;

    const std::size_t cur_rows = tb.n_rows();
    const std::size_t cap_rows = cur_rows * static_cast<std::size_t>(G);
    std::vector<GenotypeIndex> ns(cap_rows * new_w);
    std::vector<double> np(cap_rows);

    std::size_t keep = 0;
    for (std::size_t r = 0; r < cur_rows; ++r) {
        const double base = tb.prob[r];
        const GenotypeIndex gp = (col_fa >= 0)
            ? tb.state[r * old_w + static_cast<std::size_t>(col_fa)]
            : GenotypeIndex(-1);
        const GenotypeIndex gm = (col_mo >= 0)
            ? tb.state[r * old_w + static_cast<std::size_t>(col_mo)]
            : GenotypeIndex(-1);
        for (GenotypeIndex gc = 0; gc < G; ++gc) {
            double factor;
            if (fl) {
                factor = t.hwe[static_cast<std::size_t>(gc)];
            } else if (both_excluded) {
                factor = 0.0;
            } else {
                factor = child_cond_g(t, gp, gm, gc);
            }
            const double pr = base * factor;
            if (pr <= 0.0) continue;
            std::copy(tb.state.begin() + r * old_w,
                      tb.state.begin() + (r + 1) * old_w,
                      ns.begin() + keep * new_w);
            ns[keep * new_w + old_w] = gc;
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
void eliminate_member_g(PeelTable& tb, MemberIndex m) {
    const std::size_t mi = static_cast<std::size_t>(m);
    const int col = tb.col_of[mi];
    if (col < 0) return;
    const std::size_t w = tb.width();
    const std::size_t nw = w - 1;
    const std::size_t rows = tb.n_rows();

    std::map<std::vector<GenotypeIndex>, double> agg;
    std::vector<GenotypeIndex> key(nw);
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
PeelTable build_peeled_g(const Pedigree& p, const Tables& t,
                         const PeelSpec& sp,
                         const std::vector<char>& keep) {
    const std::size_t n = static_cast<std::size_t>(p.n_members);

    std::vector<std::vector<MemberIndex>> children(n);
    for (MemberIndex c = 0; c < p.n_members; ++c) {
        const MemberIndex f = p.father[static_cast<std::size_t>(c)];
        const MemberIndex m = p.mother[static_cast<std::size_t>(c)];
        if (f != kNoParent) children[static_cast<std::size_t>(f)].push_back(c);
        if (m != kNoParent) children[static_cast<std::size_t>(m)].push_back(c);
    }

    // added_set = keep U ancestors(keep): exactly the members the
    // ancestor-first recursion introduces. Children outside it are
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
                eliminate_member_g(tb, x);
                changed = true;
                break;
            }
        }
    };

    // Iterative ancestor-first introduction (explicit stack: the
    // recursion depth is the generation count; a stack keeps the kernel
    // free of deep native recursion).
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
            const MemberIndex par[2] = { p.father[mi], p.mother[mi] };
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
        introduce_member_g(tb, t, p, sp, m);
        eliminate_fixpoint();
    }
    eliminate_fixpoint();
    return tb;
}

// Permute columns into ascending member-id order so the H1 and H2
// tables align by a plain row-key lookup.
void canonicalize_g(PeelTable& tb) {
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
    std::vector<GenotypeIndex> ns(rows * w);
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

// F5.6 peeled assembly: marginal joint over the keep-set, H1 / H2 outer
// joined exactly like the dense merge but over the kept members only.
JointTable cpt_marker_joint_peeled(
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

    PeelTable h1 = build_peeled_g(p, t, sp_h1, keep);
    canonicalize_g(h1);

    PeelTable h2;
    if (has_poi) {
        PeelSpec sp_h2 = sp_h1;
        sp_h2.founder_like[static_cast<std::size_t>(poi)] = 1;
        sp_h2.excluded[static_cast<std::size_t>(poi)] = 1;
        h2 = build_peeled_g(p, t, sp_h2, keep);
        canonicalize_g(h2);
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
    return out;
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
        return ok_result(cpt_marker_joint_peeled(p, t, relevant));
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
