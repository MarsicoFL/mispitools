#ifndef MISPITOOLS_CORE_PEDIGREE_H
#define MISPITOOLS_CORE_PEDIGREE_H

#include <cstdint>
#include <vector>

namespace mispitools {
namespace core {

using MemberIndex = std::int32_t;
constexpr MemberIndex kNoParent = -1;

/// @brief Plain-old-data pedigree representation.
///
/// All vectors have length n_members. External labels (the R id column)
/// are erased at the binding boundary; the kernel operates on dense
/// 0-based indices only. Members are expected in a topologically valid
/// order: for any non-founder i, father[i] < i and mother[i] < i when
/// they are not kNoParent. The binding layer enforces this; the kernel
/// may rely on it for sorting tie-breaks.
struct Pedigree {
    MemberIndex n_members = 0;
    std::vector<MemberIndex> father;   // size n_members; kNoParent for founders
    std::vector<MemberIndex> mother;   // size n_members; kNoParent for founders
    MemberIndex poi = kNoParent;       // person of interest (excluded under H2)
};

/// @brief True iff `i` has no parents (both father[i] and mother[i] are kNoParent).
bool is_founder(const Pedigree& p, MemberIndex i) noexcept;

/// @brief Founder indices of `p`, in ascending MemberIndex order.
std::vector<MemberIndex> founders_of(const Pedigree& p);

/// @brief Non-founder indices of `p`, in ascending MemberIndex order.
std::vector<MemberIndex> nonfounders_of(const Pedigree& p);

/// @brief Topological order of `nonfounders` treating `excluded` members
/// as already assigned (their genotypes are HWE-imputed for descendants
/// under H2). Returns a permutation of `nonfounders \ excluded` such
/// that each member's parents are scheduled before it (or are excluded).
/// Mirrors R/r_ref_cpt.R::ancestral_order().
std::vector<MemberIndex> ancestral_order(
    const Pedigree& p,
    const std::vector<MemberIndex>& nonfounders,
    const std::vector<MemberIndex>& excluded);

// Placeholder retained for the F0.5 cpp-bootstrap regression test.
int pedigree_placeholder(int x);

}  // namespace core
}  // namespace mispitools

#endif  // MISPITOOLS_CORE_PEDIGREE_H
