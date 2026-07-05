#include "AeroAngles.h"

#include <algorithm>
#include <cmath>

AeroAngles ComputeAeroAngles(const Eigen::Quaterniond& q_body_to_local,
                              double fpa, double v_azi) {
    // Vehicle velocity direction in the local frame (see header comment).
    Eigen::Vector3d v_hat_local(std::sin(fpa),
                                std::cos(fpa) * std::cos(v_azi),
                                std::cos(fpa) * std::sin(v_azi));
    // Wind direction (direction the flow moves past the vehicle) is the
    // negative of the velocity direction.
    Eigen::Vector3d wind_hat_local = -v_hat_local;

    // Rotate into body axes. q is body-to-local, so its conjugate (equal to
    // its inverse for a unit quaternion) rotates local -> body.
    Eigen::Vector3d wind_hat_body = q_body_to_local.conjugate() * wind_hat_local;

    AeroAngles out;
    out.beta_rad = std::asin(std::clamp(wind_hat_body.y(), -1.0, 1.0));
    out.alpha_rad = std::atan2(wind_hat_body.z(), -wind_hat_body.x());
    return out;
}
