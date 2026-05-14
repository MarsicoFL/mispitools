#ifndef MISPITOOLS_CORE_CPT_ENGINE_H
#define MISPITOOLS_CORE_CPT_ENGINE_H

#include <cstddef>
#include <vector>

#include "marker.h"
#include "mutation_models.h"
#include "pedigree.h"
#include "result.h"

namespace mispitools {
namespace core {

/// @brief Sparse joint genotype distribution over the pedigree members.
///
/// `states_flat[row * n_members + i]` is the 0-based GenotypeIndex of
/// member `i` in `row`; `p_h1[row]` and `p_h2[row]` are matched
/// probabilities under H1 (pedigree topology) and H2 (POI HWE-independent
/// from the rest). Rows are lex-sorted by the member-state tuple in
/// member-index order, so two `JointTable`s on the same support align
/// by a parallel scan (see DESIGN.md §9).
struct JointTable {
    MemberIndex n_members = 0;
    GenotypeIndex n_genotypes = 0;
    std::vector<GenotypeIndex> states_flat;  // length n_rows * n_members
    std::vector<double> p_h1;                // length n_rows
    std::vector<double> p_h2;                // length n_rows

    std::size_t n_rows() const noexcept { return p_h1.size(); }
    bool empty() const noexcept { return p_h1.empty(); }
};

/// @brief Joint genotype distribution for a single marker under H1 / H2.
///
/// Mirrors R/r_ref_cpt.R::cpt_marker_joint_R() bit-for-bit (1e-12 tol in
/// the cross-check tests). Pedigree-side founder enumeration + Mendelian
/// peeling; under H2 the POI is HWE-marginalised and the missing-parent
/// contribution to its descendants is HWE-imputed.
///
/// Mutation models other than `MutationKind::None` are returned as a
/// `Result::error` in F2.2 and are wired in F2.4.
///
/// @complexity O(G^F * G^N) worst case, where G = K*(K+1)/2, F = number
/// of founders, N = number of non-founders. The row-set shrinks rapidly
/// after each peeling step due to Mendelian sparsity.
Result<JointTable> cpt_marker_joint(
    const Pedigree& ped,
    const Marker& marker,
    const MutationModel& mut);

// Placeholder retained for the F0.5 cpp-bootstrap regression test.
int cpt_engine_placeholder(int x);

}  // namespace core
}  // namespace mispitools

#endif  // MISPITOOLS_CORE_CPT_ENGINE_H
