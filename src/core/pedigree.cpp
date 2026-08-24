#include "pedigree.h"

#include <algorithm>
#include <cstdint>
#include <utility>

namespace mispitools {
namespace core {

bool is_founder(const Pedigree& p, MemberIndex i) noexcept {
    return p.father[i] == kNoParent && p.mother[i] == kNoParent;
}

std::vector<MemberIndex> founders_of(const Pedigree& p) {
    std::vector<MemberIndex> out;
    out.reserve(static_cast<std::size_t>(p.n_members));
    for (MemberIndex i = 0; i < p.n_members; ++i) {
        if (is_founder(p, i)) out.push_back(i);
    }
    return out;
}

std::vector<MemberIndex> nonfounders_of(const Pedigree& p) {
    std::vector<MemberIndex> out;
    out.reserve(static_cast<std::size_t>(p.n_members));
    for (MemberIndex i = 0; i < p.n_members; ++i) {
        if (!is_founder(p, i)) out.push_back(i);
    }
    return out;
}

std::vector<MemberIndex> ancestral_order(
        const Pedigree& p,
        const std::vector<MemberIndex>& nonfounders,
        const std::vector<MemberIndex>& excluded) {
    std::vector<std::uint8_t> assigned(static_cast<std::size_t>(p.n_members), 0);
    for (MemberIndex i = 0; i < p.n_members; ++i) {
        if (is_founder(p, i)) assigned[static_cast<std::size_t>(i)] = 1;
    }
    for (MemberIndex e : excluded) {
        if (e >= 0 && e < p.n_members) {
            assigned[static_cast<std::size_t>(e)] = 1;
        }
    }

    std::vector<MemberIndex> remaining;
    remaining.reserve(nonfounders.size());
    for (MemberIndex nf : nonfounders) {
        bool is_excluded = false;
        for (MemberIndex e : excluded) {
            if (e == nf) { is_excluded = true; break; }
        }
        if (!is_excluded) remaining.push_back(nf);
    }

    std::vector<MemberIndex> result;
    result.reserve(remaining.size());

    // Topological pass. Mirrors R-ref: scan `remaining` in current order,
    // emit any member whose parents are already assigned. Bounded by
    // remaining.size() outer iterations.
    while (!remaining.empty()) {
        bool progress = false;
        std::vector<MemberIndex> next_remaining;
        next_remaining.reserve(remaining.size());
        for (MemberIndex nf : remaining) {
            const MemberIndex fa = p.father[static_cast<std::size_t>(nf)];
            const MemberIndex mo = p.mother[static_cast<std::size_t>(nf)];
            const bool fa_ok =
                (fa == kNoParent) || (assigned[static_cast<std::size_t>(fa)] != 0);
            const bool mo_ok =
                (mo == kNoParent) || (assigned[static_cast<std::size_t>(mo)] != 0);
            if (fa_ok && mo_ok) {
                result.push_back(nf);
                assigned[static_cast<std::size_t>(nf)] = 1;
                progress = true;
            } else {
                next_remaining.push_back(nf);
            }
        }
        remaining = std::move(next_remaining);
        if (!progress) break;  // unresolved cycle; caller handles
    }
    return result;
}

int pedigree_placeholder(int x) {
    return x + 1;
}

}  // namespace core
}  // namespace mispitools
