#ifndef GUIDANCE_SCP_FULL_LOOP_NONCONVEX_H
#define GUIDANCE_SCP_FULL_LOOP_NONCONVEX_H

#include "ClarabelSocpSolver.h"
#include <Eigen/Dense>
#include <utility>
#include <vector>

// Reusable, standalone-testable row builders for the full 6DOF SCP loop's
// nonconvex-term handling: lossless-convexified thrust bound, gimbal cone,
// and flap/actuator box constraints.
// 
// Each function returns row content (columns/coefficients/rhs + cone type)
// rather than writing directly into a shared sparse-matrix accumulator
//
// Sign convention (carried forward from Stage A/B, critical -- see the plan's
// context section): a row's coeffs/rhs enforce coeffs.x <= rhs for
// NonnegativeConeT rows (Ax+s=b, s>=0, i.e. s=b-Ax>=0). For SecondOrderConeT
// groups, row 0 is the cone's "scalar" component and the rest are the
// "vector" components, following Stage B's existing SOC row-ordering
namespace guidance_scp {

struct ConeRow {
    std::vector<std::pair<int, double>> coeffs;
    double rhs = 0.0;
};

struct ConeRowGroup {
    std::vector<ConeRow> rows;
    clarabel::SupportedConeT<double> cone;
};

// Lossless-convexified thrust-magnitude bound: ||u|| <= sigma_u (SOC, dim 4),
// reusing Stage B's exact row pattern (reference_stage_b.cpp). The full
// loop's decision variables at these columns are deviations from a
// reference (delta_u, delta_sigma_u), not absolute values -> u_ref/
// sigma_u_ref shift the row rhs so the cone is enforced on the ABSOLUTE
// quantities (u_ref+delta_u, sigma_u_ref+delta_sigma_u), exactly mirroring
// how Stage A's path-constraint rows fold the reference value into rhs
// (e.g. "max_heat_flux - qdot_ref"). Pass a zero reference for a caller
// that has absolute decision variables.
ConeRowGroup buildThrustMagnitudeSocRows(int sigma_u_col, int ux_col, int uy_col, int uz_col,
                                           double sigma_u_ref, const Eigen::Vector3d& u_vec_ref);

// Gimbal cone: cos(deltaE_max)*||u|| <= u . e1_ref (SOC, dim 4). e1_ref is
// evaluated at the reference attitude (fixed within one outer SCP
// iteration, updated between iterations) -- keeps this row linear in u even
// though the cone's axis itself depends on attitude. u_vec_ref shifts the
// row rhs for the same deviation-from-reference reason as above.
ConeRowGroup buildGimbalConeSocRows(int ux_col, int uy_col, int uz_col,
                                      const Eigen::Vector3d& e1_ref, double deltaE_max_rad,
                                      const Eigen::Vector3d& u_vec_ref);

// Symmetric box constraint |var_ref + delta_var| <= bound, as two
// NonnegativeConeT(1) rows (upper and lower bound). Used for flap
// deflection/rate, motor torque, and gimbal-rate box constraints.
std::vector<ConeRowGroup> buildSymmetricBoxRows(int var_col, double var_ref, double bound);

}  // namespace guidance_scp

#endif  // GUIDANCE_SCP_FULL_LOOP_NONCONVEX_H
