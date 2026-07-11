#include "full_loop_phase1.h"
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

// Bridges Phase 1's real 6DOF attitude quaternion back to Stage A's
// bank-angle-parameterized translational EOM (L*cos(sigma)/L*sin(sigma)):
// projects the body-frame "lift plane normal" (0,0,1), rotated into the
// local topocentric frame via q, onto the plane perpendicular to velocity,
// then reads off the angle between that projection and a "vertical" basis
// vector spanning that plane. This only feeds the translational force
// resolution below, rotational dynamics uses the real torques directly,
// so this approximation does not leak into attitude dynamics. Documented,
// deliberate simplification
double extractBankAngleFromAttitude(const Eigen::Quaterniond& q, double fpa, double v_azi) {
    const Eigen::Vector3d n_local = q * Eigen::Vector3d(0.0, 0.0, 1.0);
    const Eigen::Vector3d v_hat_local(std::sin(fpa), std::cos(fpa) * std::cos(v_azi),
                                       std::cos(fpa) * std::sin(v_azi));

    Eigen::Vector3d up_local(1.0, 0.0, 0.0);
    Eigen::Vector3d u_vert = up_local - up_local.dot(v_hat_local) * v_hat_local;
    const double u_vert_norm = u_vert.norm();
    if (u_vert_norm < 1e-9) {
        u_vert = Eigen::Vector3d(0.0, 1.0, 0.0);  // degenerate: velocity nearly vertical
    } else {
        u_vert /= u_vert_norm;
    }
    const Eigen::Vector3d u_lat = v_hat_local.cross(u_vert);

    Eigen::Vector3d n_perp = n_local - n_local.dot(v_hat_local) * v_hat_local;
    const double n_perp_norm = n_perp.norm();
    if (n_perp_norm < 1e-9) return 0.0;
    n_perp /= n_perp_norm;
    return std::atan2(n_perp.dot(u_lat), n_perp.dot(u_vert));
}

}

Eigen::VectorXd phase1Eom(const Eigen::VectorXd& x, const Eigen::VectorXd& u,
                          const PlanetConfig& planet_config,
                          const SpacecraftConfig& spacecraft_config,
                          const FlapActuatorConfig& flap_config) {
    const double r = x(0), la = x(1), v = x(3), fpa = x(4), v_azi = x(5);
    const double q1 = x(6), q2 = x(7), q3 = x(8), q4 = x(9);
    const double wx = x(10), wy = x(11), wz = x(12);
    const double d1 = x(13), d2 = x(14), d3 = x(15), d4 = x(16);
    const double ddot1 = x(17), ddot2 = x(18), ddot3 = x(19), ddot4 = x(20);
    const double tau_m[4] = {u(0), u(1), u(2), u(3)};

    Eigen::Quaterniond q(q4, q1, q2, q3);  // Eigen ctor is (w,x,y,z); q4 is the scalar part

    // Aero closure: density/Mach, real (alpha,beta) from attitude+flight-path
    // (Not a trim assumption -> Stage A's alphaTrimDeg is inapplicable now
    // that attitude/flaps are real states), 4-independent-flap -> 3-axis
    // table query, hinge-moment-extended coefficient lookup.
    const double rho = DescentDynamics::atmosphereDensity(r, planet_config);
    const double a_sound = DescentDynamics::speedOfSound(r, planet_config);
    const double mach = v / a_sound;

    AeroAngles angles = ComputeAeroAngles(q, fpa, v_azi);
    const double alpha_deg = angles.alpha_rad * 180.0 / kPi;
    const double beta_deg = angles.beta_rad * 180.0 / kPi;

    FlapAxes axes = mapGroupDeflectionsToFlapAxes(d1, d2, d3, d4);
    const double fwd_sym_deg = axes.fwd_sym_rad * 180.0 / kPi;
    const double aft_sym_deg = axes.aft_sym_rad * 180.0 / kPi;
    const double aft_diff_deg = axes.aft_diff_rad * 180.0 / kPi;

    auto aero = spacecraft_config.aero_table.interpolate(mach, alpha_deg, beta_deg,
                                                           fwd_sym_deg, aft_sym_deg, aft_diff_deg);

    const double qbar = 0.5 * rho * v * v;
    const double L = qbar * spacecraft_config.S_ref * aero.CL;
    const double D = qbar * spacecraft_config.S_ref * aero.CD;
    const double mass = spacecraft_config.mass;

    // Translational EOM: Stage A's exact structure (reference_stage_a.cpp's
    // stageAEom), with sigma derived from the real attitude above instead of
    // being a free state.
    const double sigma = extractBankAngleFromAttitude(q, fpa, v_azi);

    const double mu = planet_config.mu, j2 = planet_config.j2, r_ref = planet_config.radius;
    const double sin_la = std::sin(la), cos_la = std::cos(la);
    const double r_ratio = r_ref / r;
    const double g_r = (mu / (r * r)) * (1.0 - 3.0 * j2 * (r_ratio * r_ratio) * 0.5 * (3.0 * sin_la * sin_la - 1.0));
    const double g_lambda = (3.0 * mu / (r * r)) * j2 * (r_ratio * r_ratio) * cos_la * sin_la;
    const double omega = planet_config.omega;

    const double cos_fpa = std::cos(fpa), sin_fpa = std::sin(fpa);
    const double cos_vazi = std::cos(v_azi), sin_vazi = std::sin(v_azi);

    const double r_dot = v * sin_fpa;
    const double la_dot = v / r * cos_fpa * cos_vazi;
    const double lo_dot = (v * cos_fpa * sin_vazi) / (r * cos_la);

    const double v_dot = -D / mass - g_r * sin_fpa + g_lambda * cos_fpa * cos_vazi
                        - omega * omega * r * cos_la * (cos_fpa * cos_vazi * cos_la - sin_fpa * sin_la);

    const double fpa_dot_term = mass * v * v / r * cos_fpa + L * std::cos(sigma)
                        - mass * g_r * cos_fpa - mass * g_lambda * sin_fpa * cos_vazi
                        + mass * omega * omega * r * cos_la * (sin_fpa * cos_vazi * sin_la - cos_fpa * cos_la)
                        + 2 * mass * omega * v * sin_vazi * cos_la;
    const double fpa_dot = fpa_dot_term / (mass * v);

    const double v_azi_dot_term = mass * v * v / r * cos_fpa * cos_fpa * sin_vazi * std::tan(la)
                        + L * std::sin(sigma) - mass * g_lambda * sin_vazi
                        + mass * omega * omega * r * sin_vazi * sin_la * cos_la
                        - 2 * mass * omega * v * (sin_fpa * cos_vazi * cos_la - cos_fpa * sin_la);
    const double v_azi_dot = v_azi_dot_term / (mass * v * cos_fpa);

    // Rotational: real hinge-moment/table-derived torque, no gimbal (engines off).
    Eigen::Quaterniond qdot = QuaternionDerivative(q, wx, wy, wz);

    Eigen::Vector3d tau_aero(
        qbar * spacecraft_config.S_ref * spacecraft_config.L_ref * aero.Cl_roll,
        qbar * spacecraft_config.S_ref * spacecraft_config.L_ref * aero.Cm,
        qbar * spacecraft_config.S_ref * spacecraft_config.L_ref * aero.Cn_yaw);
    Eigen::Vector3d w(wx, wy, wz);
    const Eigen::Matrix3d& J = spacecraft_config.inertia;
    Eigen::Vector3d wdot = J.inverse() * (tau_aero - w.cross(J * w));

    // Flap actuator 2nd-order dynamics, per flap: H_i is the hinge moment
    // (physical units) from the table's Ch[i] coefficient
    // converted via the same qbar*S_ref*L_ref normalization as Cm.
    const double ddot_i[4] = {ddot1, ddot2, ddot3, ddot4};
    double ddot_accel[4];
    for (int i = 0; i < 4; ++i) {
        const double H_i = qbar * spacecraft_config.S_ref * spacecraft_config.L_ref * aero.Ch[i];
        ddot_accel[i] = flap_config.Jeff_inv *
            (H_i - flap_config.b_damping_n_m_s * ddot_i[i] + flap_config.N_gear_ratio * tau_m[i]);
    }

    Eigen::VectorXd xdot(kPhase1StateDim);
    xdot << r_dot, la_dot, lo_dot, v_dot, fpa_dot, v_azi_dot,
            qdot.x(), qdot.y(), qdot.z(), qdot.w(),
            wdot(0), wdot(1), wdot(2),
            ddot1, ddot2, ddot3, ddot4,
            ddot_accel[0], ddot_accel[1], ddot_accel[2], ddot_accel[3];
    return xdot;
}

}  // namespace guidance_scp
