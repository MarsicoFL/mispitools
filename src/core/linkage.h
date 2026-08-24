#ifndef MISPITOOLS_CORE_LINKAGE_H
#define MISPITOOLS_CORE_LINKAGE_H

#include <cmath>
#include <cstddef>
#include <vector>

#include "marker.h"
#include "mutation_models.h"
#include "pedigree.h"
#include "result.h"

namespace mispitools {
namespace core {

/// @brief Sparse joint genotype distribution over a *pair* of linked
/// markers.
///
/// Analogous to `JointTable` (cpt_engine.h) but carries one genotype per
/// marker per member: `states_a[row * n_members + i]` and
/// `states_b[row * n_members + i]` are the 0-based GenotypeIndex of member
/// `i` at marker A and marker B respectively in `row`. `p_h1[row]` /
/// `p_h2[row]` are the matched probabilities under H1 (pedigree topology
/// with recombination `rho`) and H2 (POI HWE-independent at both loci,
/// the two loci themselves at population linkage equilibrium). Rows are
/// the *observable* unordered two-locus genotype configurations: latent
/// phase (which parental haplotype carries which allele) is summed out,
/// then rows are lex-sorted on
/// (A0, B0, A1, B1, ..., A_{n-1}, B_{n-1}) so two tables on the same
/// support align by a parallel scan, exactly like `JointTable`. KL
/// (kl_engine) and the LR distribution (lr_dist) consume `p_h1` / `p_h2`
/// only and apply unchanged (DESIGN.md §14).
struct LinkedJointTable {
    MemberIndex n_members = 0;
    GenotypeIndex n_genotypes_a = 0;
    GenotypeIndex n_genotypes_b = 0;
    std::vector<GenotypeIndex> states_a;  // length n_rows * n_members
    std::vector<GenotypeIndex> states_b;  // length n_rows * n_members
    std::vector<double> p_h1;             // length n_rows
    std::vector<double> p_h2;             // length n_rows

    std::size_t n_rows() const noexcept { return p_h1.size(); }
    bool empty() const noexcept { return p_h1.empty(); }
};

/// @brief Joint genotype distribution for two linked markers under
/// H1 / H2 at recombination fraction `rho`.
///
/// Phase-aware founder enumeration: each member carries an *ordered*
/// two-locus diplotype (a paternal-origin haplotype and a maternal-origin
/// haplotype). Founders are drawn from HWE at each locus and population
/// linkage equilibrium between loci; the ordered enumeration covers both
/// founder phases automatically (including the doubly-heterozygous case,
/// avoiding the classic phase-duplication bug). Non-founders inherit a
/// gamete from each parent through the two-locus transmission kernel
/// `trans_prob_MM` (recombination `rho`, per-locus mutation). Under H2 the
/// POI is removed and HWE-marginalised at both loci (loci independent),
/// mirroring `cpt_marker_joint`'s single-marker H2. The latent phase is
/// summed out and identical observable configurations are aggregated.
///
/// At `rho == 0.5` the loci are unlinked and the joint factorises into
/// the product of the two single-marker joints; the general algorithm is
/// exact in that limit (a fast-path delegation is a deferred F7 perf
/// optimisation, see notes in the .cpp).
///
/// @param p Acyclic pedigree (loops broken at the R boundary). Inbred
///   founders / selfing are not supported for linked markers and must be
///   rejected by the binding layer (SCOUT_pedprobr_linkage.md).
/// @param m_A First marker.
/// @param m_B Second marker.
/// @param mu_A Mutation model for marker A.
/// @param mu_B Mutation model for marker B.
/// @param rho Recombination fraction in [0, 0.5].
/// @param relevant 0-based member indices whose genotypes the caller will
///   read. Empty (default) means "the whole pedigree": every member is
///   kept and the result is bit-for-bit the F5.2 dense joint. When a
///   strict subset is given, every other member is summed out by
///   Elston-Stewart online variable elimination (F5.5) and only the
///   requested members (plus the POI, forced in to keep H1/H2 alignable)
///   carry a meaningful genotype in the returned table; marginalised
///   members hold genotype index 1 (0-based 0) and must not be read.
/// @return The linked-pair joint table, or a Result::error.
/// @complexity Dense path O(D^F * D^N) worst case (D = (K_A*K_B)^2
///   ordered diplotypes). With a strict keep-set the active frontier is
///   the Elston-Stewart treewidth of the relevant sub-pedigree, so
///   MP-typical pedigrees (first-cousin and larger, <=20 founders,
///   <=3 generations) stay tractable (F5.5).
Result<LinkedJointTable> linked_pair_joint(
    const Pedigree& p,
    const Marker& m_A,
    const Marker& m_B,
    const MutationModel& mu_A,
    const MutationModel& mu_B,
    double rho,
    const std::vector<MemberIndex>& relevant = {});

/// @brief Two-locus phased transmission kernel.
///
/// Adapted from pedprobr v1.0.1 (R/peel.R, `.transProbMM`), GPL (>= 2),
/// by M. D. Vigeland (UiO); see COPYRIGHTS. Probability that a parent
/// with paternal-origin haplotype `(p1, p2)` and maternal-origin
/// haplotype `(m1, m2)` transmits the gamete `(g1, g2)` under
/// recombination `rho` and per-locus mutation matrices. `mut_k` is the
/// K_k x K_k row-major matrix with `mut_k[a * K_k + b]` = P(parental
/// allele a -> gamete allele b) at locus k (identity for model none).
/// The four terms are {two non-recombinant phases, two recombinant
/// phases}; their weights sum to 1, and the `0.5` is the Mendelian
/// selection of the grand-parental starting haplotype (already
/// normalised, *not* an extra `rho` factor).
///
/// @param p1,p2 Parent paternal-origin alleles at locus 1, 2.
/// @param m1,m2 Parent maternal-origin alleles at locus 1, 2.
/// @param g1,g2 Gamete alleles at locus 1, 2.
/// @param rho Recombination fraction in [0, 0.5].
/// @param mut1 Locus-1 mutation matrix (K_A x K_A, row-major).
/// @param mut2 Locus-2 mutation matrix (K_B x K_B, row-major).
/// @param K_A Number of alleles at locus 1.
/// @param K_B Number of alleles at locus 2.
inline double trans_prob_MM(
        AlleleIndex p1, AlleleIndex p2,
        AlleleIndex m1, AlleleIndex m2,
        AlleleIndex g1, AlleleIndex g2,
        double rho,
        const std::vector<double>& mut1,
        const std::vector<double>& mut2,
        AlleleIndex K_A, AlleleIndex K_B) noexcept {
    const double l1p = mut1[static_cast<std::size_t>(p1) * K_A + g1];
    const double l1m = mut1[static_cast<std::size_t>(m1) * K_A + g1];
    const double l2p = mut2[static_cast<std::size_t>(p2) * K_B + g2];
    const double l2m = mut2[static_cast<std::size_t>(m2) * K_B + g2];
    return (l1p * l2p * (1.0 - rho) + l1m * l2m * (1.0 - rho)
          + l1p * l2m * rho         + l1m * l2p * rho) * 0.5;
}

/// @brief Map functions adapted from pedprobr v1.0.1 (R/haldane.R),
/// GPL (>= 2), by M. D. Vigeland; see COPYRIGHTS. Convert genetic
/// distance (cM) to / from the recombination fraction `rho`.
inline double haldane_cM_to_rho(double cM) noexcept {
    return 0.5 * (1.0 - std::exp(-cM / 50.0));
}
inline double haldane_rho_to_cM(double rho) noexcept {
    return -50.0 * std::log(1.0 - 2.0 * rho);
}
inline double kosambi_cM_to_rho(double cM) noexcept {
    const double e = std::exp(cM / 25.0);
    return 0.5 * (e - 1.0) / (e + 1.0);
}
inline double kosambi_rho_to_cM(double rho) noexcept {
    return 25.0 * std::log((1.0 + 2.0 * rho) / (1.0 - 2.0 * rho));
}

// Placeholder retained for the F0.5 cpp-bootstrap regression test.
int linkage_placeholder(int x);

}  // namespace core
}  // namespace mispitools

#endif  // MISPITOOLS_CORE_LINKAGE_H
