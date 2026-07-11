#ifndef AERO_ANGLES_H
#define AERO_ANGLES_H

#include <Eigen/Geometry>

// 3D angle of attack (alpha) and sideslip (beta) are derived from the vehicle's
// translational flight-path state and attitude quaternion.

struct AeroAngles {
    double alpha_rad = 0.0;
    double beta_rad = 0.0;
};

AeroAngles ComputeAeroAngles(const Eigen::Quaterniond& q_body_to_local,
                              double fpa, double v_azi);

// Overload taking the velocity direction (unit vector) directly in
// frame q_body_to_local's rotation is defined against
AeroAngles ComputeAeroAngles(const Eigen::Quaterniond& q_body_to_local,
                              const Eigen::Vector3d& v_hat_local);

#endif // AERO_ANGLES_H
