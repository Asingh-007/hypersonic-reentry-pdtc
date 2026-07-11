#ifndef GUIDANCE_SCP_TRIM_SOLVER_H
#define GUIDANCE_SCP_TRIM_SOLVER_H

#include "AeroCoefficientTable.h"


namespace guidance_scp {

// Root-finds alpha_deg such that Cm(mach, alpha_deg, beta_deg=0, 0,0,0) == 0
// (trim: zero pitching moment about the table's moment_ref, symmetric flaps
// neutral) via bisection over [alpha_min_deg, alpha_max_deg].
//
// AeroCoefficientTable::interpolate() clamps queries outside its own grid to
// the nearest boundary rather than extrapolating (see its doc comment), so
// this bisection is always well-defined even if the true trim point lies
// outside the table's fitted envelope, the result is then that boundary's
// trim value.
//
// If Cm doesn't change sign across [alpha_min_deg, alpha_max_deg] (no
// bracketable root possible with a coarse/preliminary aero table), falls
// back to the 1-degree-sampled alpha_deg minimizing |Cm|. This fallback is a
// PLACEHOLDER: a real trimmed vehicle configuration may need per-Mach flap
// bias or a wider search, revisit once a production-quality aero table
// exists.
double alphaTrimDeg(double mach, const aero_model::AeroCoefficientTable& aero_table,
                     double alpha_min_deg = -20.0, double alpha_max_deg = 20.0,
                     double tol_deg = 1e-3, int max_iter = 60);

}  // namespace guidance_scp

#endif  // GUIDANCE_SCP_TRIM_SOLVER_H
