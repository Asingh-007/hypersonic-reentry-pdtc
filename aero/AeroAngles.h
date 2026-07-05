#ifndef AERO_ANGLES_H
#define AERO_ANGLES_H

#include <Eigen/Geometry>

// Derives 3D angle of attack (alpha) and sideslip (beta) from the vehicle's
// translational flight-path state and attitude quaternion.

struct AeroAngles {
    double alpha_rad = 0.0;
    double beta_rad = 0.0;
};

AeroAngles ComputeAeroAngles(const Eigen::Quaterniond& q_body_to_local,
                              double fpa, double v_azi);

#endif // AERO_ANGLES_H
