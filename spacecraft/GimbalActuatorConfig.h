#ifndef SPACECRAFT_GIMBAL_ACTUATOR_CONFIG_H
#define SPACECRAFT_GIMBAL_ACTUATOR_CONFIG_H

#include <Eigen/Dense>

// PLACEHOLDER: no real TVC gimbal hardware spec exists yet. deltaE/phiE (the
// two gimbal-angle axes) are NOT states in this formulation -- only their
// rates are Phase-2 controls (see guidance_scp/full_loop_phase2.h) -- this
// struct only bounds those rates and the total gimbal throw.
// TODO: replace with actual TVC hardware data once available.
struct GimbalActuatorConfig {
    double deltaE_max_rad = 0.13962634015954636;  // +/-8 deg gimbal throw
    double gimbal_rate_max_rad_s = 0.2;           // max gimbal rate, rad/s
    // Aft engine gimbal-point location, body frame, m -- no real vehicle
    // mass-properties/engine-placement data exists yet.
    Eigen::Vector3d engine_gimbal_point_body_m = Eigen::Vector3d(-18.0, 0.0, 0.0);
};

#endif  // SPACECRAFT_GIMBAL_ACTUATOR_CONFIG_H
