#ifndef GUIDANCE_SCP_FULL_LOOP_TRANSITION_H
#define GUIDANCE_SCP_FULL_LOOP_TRANSITION_H

#include "PlanetConfig.h"
#include <Eigen/Dense>

// Phase transition (entry -> flip -> descent), Section 5 of the instructions
// doc: converts Phase 1's terminal translational state (spherical/rotating-
// planet frame) into Phase 2's flat local ENU frame (Stage B's exact
// frame). Attitude/rate/flap state carry over unchanged (same body frame
// convention both phases) -- only the translational states need conversion,
// handled here.
namespace guidance_scp {

// Phase-1-frame coordinates of Phase 2's ENU origin (the intended
// tower-catch point).
struct Phase1ToPhase2Frame {
    double origin_r = 0.0, origin_la = 0.0, origin_lo = 0.0;
};

struct Phase2TranslationalState {
    double rx = 0.0, ry = 0.0, rz = 0.0, vx = 0.0, vy = 0.0, vz = 0.0;
};

// Flat-Earth/local-tangent-plane approximation (PLACEHOLDER, consistent
// with Stage B's own flat-ENU assumption -- valid over the short horizontal
// range covered during flip+descent, no oblateness/curvature correction).
Phase2TranslationalState convertPhase1ToPhase2(double r, double la, double lo,
                                                 double v, double fpa, double v_azi,
                                                 const Phase1ToPhase2Frame& frame,
                                                 const PlanetConfig& planet_config);

// 6x6 Jacobian d(rx,ry,rz,vx,vy,vz)/d(r,la,lo,v,fpa,v_azi) of the map above,
// via the same generalized finite-difference utility used elsewhere in this
// module (FiniteDifferenceJacobian.h), since convertPhase1ToPhase2 is
// nonlinear and must be linearized about the current reference for the
// boundary-matching equality rows in the full subproblem (Milestone 9).
Eigen::MatrixXd computeTransitionJacobian(double r, double la, double lo,
                                            double v, double fpa, double v_azi,
                                            const Phase1ToPhase2Frame& frame,
                                            const PlanetConfig& planet_config);

}  // namespace guidance_scp

#endif  // GUIDANCE_SCP_FULL_LOOP_TRANSITION_H
