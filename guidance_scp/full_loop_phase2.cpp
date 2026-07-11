#include "full_loop_phase2.h"
#include "AeroRegimeDispatch.h"
#include "AeroAngles.h"
#include "DescentDynamics.h"
#include "QuaternionUtils.h"
#include <Eigen/Geometry>
#include <cmath>

namespace guidance_scp {

namespace {

using aero_model::mapGroupDeflectionsToFlapAxes;
using aero_model::FlapAxes;

constexpr double kPi = 3.14159265358979323846;

// Thrust direction in body axes from gimbal angles, mirroring
// freestreamDirectionBody's structure for consistency (nominal thrust along
// -body_x, matching NewtonianAeroModel's nose=+body_x convention).
Eigen::Vector3d thrustDirBody(double deltaE, double phiE) {
    return Eigen::Vector3d(-std::cos(deltaE) * std::cos(phiE),
                            std::sin(phiE),
                            std::sin(deltaE) * std::cos(phiE));
}

}

Eigen::VectorXd phase2Eom(const Eigen::VectorXd& x, const Eigen::VectorXd& u,
                          double deltaE_ref, double phiE_ref, double Isp_s,
                          const PlanetConfig& planet_config,
                          const SpacecraftConfig& spacecraft_config,
                          const FlapActuatorConfig& flap_config,
                          const GimbalActuatorConfig& gimbal_config) {
    const double m = x(0);
    const double rx = x(1), ry = x(2), rz = x(3);
    const double vx = x(4), vy = x(5), vz = x(6);
    const double q1 = x(7), q2 = x(8), q3 = x(9), q4 = x(10);
    const double wx = x(11), wy = x(12), wz = x(13);
    const double d1 = x(14), d2 = x(15), d3 = x(16), d4 = x(17);
    const double ddot1 = x(18), ddot2 = x(19), ddot3 = x(20), ddot4 = x(21);

    const double T = u(0), deltaE_dot = u(1), phiE_dot = u(2);
    const double tau_m[4] = {u(3), u(4), u(5), u(6)};

    Eigen::Quaterniond q(q4, q1, q2, q3);  // body-to-ENU, scalar-last
    Eigen::Vector3d v_enu(vx, vy, vz);

    const double r_eff = planet_config.radius + rz;
    const double rho = DescentDynamics::atmosphereDensity(r_eff, planet_config);
    const double a_sound = DescentDynamics::speedOfSound(r_eff, planet_config);
    const double v_mag = v_enu.norm();
    const double mach = v_mag / a_sound;

    Eigen::Vector3d v_hat_local = (v_mag > 1e-6) ? (v_enu / v_mag) : Eigen::Vector3d(0.0, 0.0, -1.0);
    AeroAngles angles = ComputeAeroAngles(q, v_hat_local);
    const double alpha_deg = angles.alpha_rad * 180.0 / kPi;
    const double beta_deg = angles.beta_rad * 180.0 / kPi;

    FlapAxes axes = mapGroupDeflectionsToFlapAxes(d1, d2, d3, d4);
    const double fwd_sym_deg = axes.fwd_sym_rad * 180.0 / kPi;
    const double aft_sym_deg = axes.aft_sym_rad * 180.0 / kPi;
    const double aft_diff_deg = axes.aft_diff_rad * 180.0 / kPi;

    auto aero = spacecraft_config.aero_table.interpolate(mach, alpha_deg, beta_deg,
                                                           fwd_sym_deg, aft_sym_deg, aft_diff_deg);

    const double qbar = 0.5 * rho * v_mag * v_mag;
    // Body-axis aero force (CA along -body_x, CN along body_z) rotated into ENU via q.
    Eigen::Vector3d F_aero_body = qbar * spacecraft_config.S_ref * Eigen::Vector3d(-aero.CA, 0.0, aero.CN);
    Eigen::Vector3d F_aero_enu = q * F_aero_body;

    // Thrust.
    Eigen::Vector3d thrust_dir_body = thrustDirBody(deltaE_ref, phiE_ref);
    Eigen::Vector3d F_thrust_body = T * thrust_dir_body;
    Eigen::Vector3d F_thrust_enu = q * F_thrust_body;

    const double g_0 = planet_config.g_0;
    Eigen::Vector3d g_enu(0.0, 0.0, -g_0);

    Eigen::Vector3d v_dot_vec = (F_thrust_enu + F_aero_enu) / m + g_enu;
    const double mdot = -T / (Isp_s * g_0);

    // Rotational: TVC + aero torque, both in body axes.
    Eigen::Vector3d tau_TVC = (gimbal_config.engine_gimbal_point_body_m - spacecraft_config.moment_ref)
                                  .cross(F_thrust_body);
    Eigen::Vector3d tau_aero(
        qbar * spacecraft_config.S_ref * spacecraft_config.L_ref * aero.Cl_roll,
        qbar * spacecraft_config.S_ref * spacecraft_config.L_ref * aero.Cm,
        qbar * spacecraft_config.S_ref * spacecraft_config.L_ref * aero.Cn_yaw);
    Eigen::Vector3d tau = tau_TVC + tau_aero;

    Eigen::Quaterniond qdot = QuaternionDerivative(q, wx, wy, wz);
    Eigen::Vector3d w(wx, wy, wz);
    const Eigen::Matrix3d& J = spacecraft_config.inertia;
    Eigen::Vector3d wdot = J.inverse() * (tau - w.cross(J * w));

    // Flap actuator 2nd-order dynamics, identical structure to Phase 1.
    const double ddot_i[4] = {ddot1, ddot2, ddot3, ddot4};
    double ddot_accel[4];
    for (int i = 0; i < 4; ++i) {
        const double H_i = qbar * spacecraft_config.S_ref * spacecraft_config.L_ref * aero.Ch[i];
        ddot_accel[i] = flap_config.Jeff_inv *
            (H_i - flap_config.b_damping_n_m_s * ddot_i[i] + flap_config.N_gear_ratio * tau_m[i]);
    }

    (void)deltaE_dot;  // not used by the EOM itself -- integrated externally into deltaE_ref between iterations
    (void)phiE_dot;

    Eigen::VectorXd xdot(kPhase2StateDim);
    xdot << mdot,
            vx, vy, vz,
            v_dot_vec(0), v_dot_vec(1), v_dot_vec(2),
            qdot.x(), qdot.y(), qdot.z(), qdot.w(),
            wdot(0), wdot(1), wdot(2),
            ddot1, ddot2, ddot3, ddot4,
            ddot_accel[0], ddot_accel[1], ddot_accel[2], ddot_accel[3];
    return xdot;
}

}  // namespace guidance_scp
