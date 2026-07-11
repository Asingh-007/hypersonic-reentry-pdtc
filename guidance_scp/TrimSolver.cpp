#include "TrimSolver.h"
#include <cmath>

namespace guidance_scp {

namespace {

double cmAt(double mach, double alpha_deg, const aero_model::AeroCoefficientTable& aero_table) {
    return aero_table.interpolate(mach, alpha_deg, /*beta_deg=*/0.0,
                                    /*fwd_sym_deg=*/0.0, /*aft_sym_deg=*/0.0, /*aft_diff_deg=*/0.0).Cm;
}

}

double alphaTrimDeg(double mach, const aero_model::AeroCoefficientTable& aero_table,
                     double alpha_min_deg, double alpha_max_deg, double tol_deg, int max_iter) {
    double cm_lo = cmAt(mach, alpha_min_deg, aero_table);
    double cm_hi = cmAt(mach, alpha_max_deg, aero_table);

    if ((cm_lo > 0.0) == (cm_hi > 0.0)) {
        // No sign change -- fall back to the 1-degree-sampled alpha minimizing |Cm|.
        double best_alpha = alpha_min_deg;
        double best_abs_cm = std::abs(cm_lo);
        for (double a = alpha_min_deg; a <= alpha_max_deg; a += 1.0) {
            double abs_cm = std::abs(cmAt(mach, a, aero_table));
            if (abs_cm < best_abs_cm) {
                best_abs_cm = abs_cm;
                best_alpha = a;
            }
        }
        return best_alpha;
    }

    double lo = alpha_min_deg, hi = alpha_max_deg;
    for (int i = 0; i < max_iter && (hi - lo) > tol_deg; ++i) {
        double mid = 0.5 * (lo + hi);
        double cm_mid = cmAt(mach, mid, aero_table);
        if ((cm_mid > 0.0) == (cm_lo > 0.0)) {
            lo = mid;
            cm_lo = cm_mid;
        } else {
            hi = mid;
        }
    }
    return 0.5 * (lo + hi);
}

}  // namespace guidance_scp
