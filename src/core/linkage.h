#ifndef MISPITOOLS_CORE_LINKAGE_H
#define MISPITOOLS_CORE_LINKAGE_H

namespace mispitools {
namespace core {

/// Placeholder for the linkage module.
/// @param x integer input
/// @return x + 1
/// Will be replaced in F5.2 with linked_pair_joint() for two markers
/// at recombination fraction theta, used inside cpt_marker_joint when
/// linked markers are declared in the model spec.
int linkage_placeholder(int x);

} // namespace core
} // namespace mispitools

#endif // MISPITOOLS_CORE_LINKAGE_H
