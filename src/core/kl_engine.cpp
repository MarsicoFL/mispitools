#include "kl_engine.h"

#include <cmath>
#include <cstddef>
#include <limits>
#include <string>

namespace mispitools {
namespace core {

int kl_engine_placeholder(int x) {
    return x + 1;
}

Result<PerMarkerKL> per_marker_kl(const JointTable& joint) {
    if (joint.p_h1.size() != joint.p_h2.size()) {
        return err_result<PerMarkerKL>(
            "per_marker_kl: P_H1 and P_H2 vectors have mismatched lengths.");
    }

    PerMarkerKL out;
    const std::size_t n = joint.p_h1.size();

    // Match R-ref bit-for-bit: accumulate finite contributions in row order
    // (cpt_marker_joint already lex-sorts), and only collapse to ±Inf at the
    // end if any absolute-continuity violation was seen. Doing it in two
    // stages instead of a single sum-of-Inf-plus-finite avoids relying on
    // IEEE Inf+x semantics (which would also work, but the explicit flag is
    // clearer at the abstraction boundary).
    bool e_h1_inf = false;
    bool e_h2_neg_inf = false;

    for (std::size_t r = 0; r < n; ++r) {
        const double p1 = joint.p_h1[r];
        const double p2 = joint.p_h2[r];
        if (p1 < 0.0 || p2 < 0.0) {
            return err_result<PerMarkerKL>(
                "per_marker_kl: negative probability at row "
                + std::to_string(r) + ".");
        }
        const bool pos1 = p1 > 0.0;
        const bool pos2 = p2 > 0.0;

        if (pos1 && pos2) {
            // R-ref: log10(P1) - log10(P2), then sum P_i * (...) left-to-right.
            const double term = std::log10(p1) - std::log10(p2);
            out.e_log10_lr_h1 += p1 * term;
            out.e_log10_lr_h2 += p2 * term;
        } else if (pos1) {
            // P1 > 0, P2 = 0: H2 lacks the support of H1.
            e_h1_inf = true;
            out.abs_cont_violations_h2 += 1;
            out.mass_violations_h2 += p1;
        } else if (pos2) {
            // P1 = 0, P2 > 0: H1 lacks the support of H2.
            e_h2_neg_inf = true;
            out.abs_cont_violations_h1 += 1;
            out.mass_violations_h1 += p2;
        }
        // (!pos1 && !pos2): cpt_marker_joint already filters these.
    }

    // R uses log(10) (natural log); std::log(10.0) is bit-equivalent.
    const double ln10 = std::log(10.0);
    const double inf  = std::numeric_limits<double>::infinity();

    if (e_h1_inf) {
        out.e_log10_lr_h1 = inf;
        out.kl_h1_to_h2 = inf;
    } else {
        out.kl_h1_to_h2 = out.e_log10_lr_h1 * ln10;
    }

    if (e_h2_neg_inf) {
        out.e_log10_lr_h2 = -inf;
        out.kl_h2_to_h1 = inf;
    } else {
        out.kl_h2_to_h1 = -out.e_log10_lr_h2 * ln10;
    }

    return ok_result(out);
}

}  // namespace core
}  // namespace mispitools
