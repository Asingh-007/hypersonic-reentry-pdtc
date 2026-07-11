#ifndef GUIDANCE_SCP_FULL_LOOP_PHASE1_H
#define GUIDANCE_SCP_FULL_LOOP_PHASE1_H

#include "PlanetConfig.h"
#include "SpacecraftConfig.h"
#include "FlapActuatorConfig.h"
#include <Eigen/Dense>

// Phase 1 (entry, engines off) full 21-state nonlinear EOM, per the full
// 6DOF SCP loop plan. State layout:
//   [r,la,lo,V,fpa,v_azi,            0-5   translational (Stage A's exact 6)
//    q1,q2,q3,q4,                    6-9   attitude quaternion, scalar-last
//    wx,wy,wz,                       10-12 body angular rate
//    d1,d2,d3,d4,                    13-16 flap deflections, rad
//    ddot1,ddot2,ddot3,ddot4]        17-20 flap deflection rates, rad/s
// Control: [tau_m1,tau_m2,tau_m3,tau_m4] (flap motor commanded torques).
//
// This is a standalone, purpose-built
// EOM for the full SCP loop's reference/linearization use, following the
// same "reference generators are their own small models" precedent Stage A
// established, reusing DescentDynamics'/AeroAngles'/QuaternionUtils' public
// building blocks rather than DescentDynamics itself (which has no flap
// states and assumes a fixed-mass, flap-neutral vehicle).
namespace guidance_scp {

constexpr int kPhase1StateDim = 21;
constexpr int kPhase1ControlDim = 4;

// Full nonlinear right-hand side. x is 21x1, u is 4x1 (tau_m1..tau_m4),
// returns xdot (21x1). See full_loop_phase1.cpp for the derivation notes
// (bank-angle bridge, hinge-moment-driven flap actuator ODEs).
Eigen::VectorXd phase1Eom(const Eigen::VectorXd& x, const Eigen::VectorXd& u,
                          const PlanetConfig& planet_config,
                          const SpacecraftConfig& spacecraft_config,
                          const FlapActuatorConfig& flap_config);

}  // namespace guidance_scp

#endif  // GUIDANCE_SCP_FULL_LOOP_PHASE1_H
