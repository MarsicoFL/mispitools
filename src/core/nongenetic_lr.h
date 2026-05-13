#ifndef MISPITOOLS_CORE_NONGENETIC_LR_H
#define MISPITOOLS_CORE_NONGENETIC_LR_H

namespace mispitools {
namespace core {

/// Placeholder for the non-genetic LR module.
/// @param x integer input
/// @return x + 1
/// Will be replaced in F6.3 / F6.4 with nongenetic_cpt(feature) and the
/// per_feature_kl_nongenetic / per_feature_lr_dist_nongenetic routines
/// for sex / age / region / hair / eyes / pigmentation / birthdate / custom.
int nongenetic_lr_placeholder(int x);

} // namespace core
} // namespace mispitools

#endif // MISPITOOLS_CORE_NONGENETIC_LR_H
