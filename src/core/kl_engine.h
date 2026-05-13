#ifndef MISPITOOLS_CORE_KL_ENGINE_H
#define MISPITOOLS_CORE_KL_ENGINE_H

namespace mispitools {
namespace core {

/// Placeholder for the KL engine module.
/// @param x integer input
/// @return x + 1
/// Will be replaced in F3.1 with per_marker_kl(joint_h1, joint_h2) that
/// returns the bidirectional KL divergence and the H1 / H2 expectations
/// of log10 LR on the shared sparse support.
int kl_engine_placeholder(int x);

} // namespace core
} // namespace mispitools

#endif // MISPITOOLS_CORE_KL_ENGINE_H
