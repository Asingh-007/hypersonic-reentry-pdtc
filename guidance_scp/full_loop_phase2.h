#ifndef GUIDANCE_SCP_FULL_LOOP_PHASE2_H
#define GUIDANCE_SCP_FULL_LOOP_PHASE2_H

#include "PlanetConfig.h"
#include "SpacecraftConfig.h"
#include "FlapActuatorConfig.h"
#include "GimbalActuatorConfig.h"
#include <Eigen/Dense>

// Phase 2 (flip + powered descent + tower catch) full 22-state nonlinear
// EOM, flat local ENU frame (Stage B's exact frame/gravity convention, +z
// up, g=[0,0,-g0]). State layout:
//   [m,                              0     mass, kg
//    rx,ry,rz,                       1-3   position, local ENU (+z up)
//    vx,vy,vz,                       4-6   velocity
//    q1,q2,q3,q4,                    7-10  attitude quaternion, scalar-last
//    wx,wy,wz,                       11-13 body angular rate
//    d1,d2,d3,d4,                    14-17 flap deflections, rad
//    ddot1,ddot2,ddot3,ddot4]        18-21 flap deflection RATES, rad/s
// Control: [T, deltaE_dot, phiE_dot, tau_m1..tau_m4] (throttle, gimbal
// rates, flap motor torques).
//
// deltaE/phiE (instantaneous gimbal angle) are not part of the state
// only their rates are controls Instead they are tracked as auxiliary per-node 
// reference bookkeeping outside x2, updated once per outer SCP iteration 
// (exactly how Stage A holds t_scale_s fixed outside its own decision vector)
//
// q here represents body-to-ENU attitude, a different "local" frame than
// Phase 1's body-to-topocentric quaternion (see full_loop_transition.h for
// the boundary-matching implications of this).
namespace guidance_scp {

constexpr int kPhase2StateDim = 22;
constexpr int kPhase2ControlDim = 7;

// Full nonlinear right-hand side. x is 22x1 (mass slot is raw m, NOT
// z=ln(m) -> that convexifying substitution only applies inside the convex
// subproblem, see full_loop_subproblem.h), u is 7x1. Isp_s is passed
// explicitly (a solve-level parameter, like StageBConfig::Isp_s, not a
// vehicle/actuator hardware constant) since no existing config struct holds it.
Eigen::VectorXd phase2Eom(const Eigen::VectorXd& x, const Eigen::VectorXd& u,
                          double deltaE_ref, double phiE_ref, double Isp_s,
                          const PlanetConfig& planet_config,
                          const SpacecraftConfig& spacecraft_config,
                          const FlapActuatorConfig& flap_config,
                          const GimbalActuatorConfig& gimbal_config);

}  // namespace guidance_scp

#endif  // GUIDANCE_SCP_FULL_LOOP_PHASE2_H
