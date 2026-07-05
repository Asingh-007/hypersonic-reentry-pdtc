#include "DescentDynamics.h"
#include "SpacecraftConfig.h"
#include "PlanetConfig.h"
#include "ControlInputs.h"
#include "AtmosphereModel.h"
#include <cmath>
#include <algorithm>
#include <iostream>

DescentDynamics::StateVector DescentDynamics::toVector(const DescentState& s) {
    StateVector v;
    v << s.r, s.la, s.lo, s.v, s.fpa, s.v_azi,
         s.q1, s.q2, s.q3, s.q4,
         s.wx, s.wy, s.wz;
    return v;
}

DescentState DescentDynamics::fromVector(const StateVector& v) {
    return DescentState(
        static_cast<float>(v(0)), static_cast<float>(v(1)), static_cast<float>(v(2)),
        static_cast<float>(v(3)), static_cast<float>(v(4)), static_cast<float>(v(5)),
        static_cast<float>(v(6)), static_cast<float>(v(7)), static_cast<float>(v(8)), static_cast<float>(v(9)),
        static_cast<float>(v(10)), static_cast<float>(v(11)), static_cast<float>(v(12)));
}

double DescentDynamics::atmosphereDensity(double r) const {
    double altitude = r - planet_config_.radius;
    if (planet_config_.body == PlanetBody::Mars) {
        return MarsAtmosphereExponential::Compute(altitude).density;
    }
    return EarthAtmosphere1976::Compute(altitude).density;
}

DescentDynamics::StateVector DescentDynamics::derivatives(const StateVector& x,
                                                           const ThrustVectorControlInputs& control_inputs) const {
    double r = x(0);
    double la = x(1);
    double lo = x(2);
    double v = x(3);
    double fpa = x(4);
    double v_azi = x(5);
    double q1 = x(6), q2 = x(7), q3 = x(8), q4 = x(9);
    double wx = x(10), wy = x(11), wz = x(12);

    double thrust = control_inputs.thrust;
    double pitch = control_inputs.pitch;
    double yaw = control_inputs.yaw;

    // J2 (oblateness) gravity perturbation model, superseding the prior
    // spherical/non-oblate assumption (g_lambda is no longer hardcoded to
    // zero). r_ref reuses PlanetConfig::radius as the single reference
    // radius -- this model doesn't distinguish equatorial/mean/polar radius.
    double mu = planet_config_.mu;
    double j2 = planet_config_.j2;
    double r_ref = planet_config_.radius;
    double sin_la = std::sin(la), cos_la = std::cos(la);
    double r_ratio = r_ref / r;
    double g_r = (mu / (r * r)) * (1.0 - 3.0 * j2 * (r_ratio * r_ratio) * 0.5 * (3.0 * sin_la * sin_la - 1.0));
    double g_lambda = (3.0 * mu / (r * r)) * j2 * (r_ratio * r_ratio) * cos_la * sin_la;

    double omega = planet_config_.omega;
    double area = spacecraft_config_.area;
    double c_d = spacecraft_config_.c_d;
    double c_l = spacecraft_config_.c_l;
    double mass = spacecraft_config_.mass;

    // Aerodynamic model: existing constant-Cd/Cl formulation (unchanged),
    // density now from the real atmosphere model instead of a placeholder
    // single exponential.
    // TODO: replace with a panel-based modified-Newtonian model once vehicle
    // geometry/panel data is available.
    double rho = atmosphereDensity(r);
    double lift = 0.5 * rho * v * v * area * c_l;
    double drag = 0.5 * rho * v * v * area * c_d;

    double r_dot = v * std::sin(fpa);
    double la_dot = v / r * std::cos(fpa) * std::cos(v_azi);
    double lo_dot = (v * std::cos(fpa) * std::sin(v_azi)) / (r * std::cos(la));

    double v_dot = (thrust * std::cos(pitch) * std::cos(yaw) - drag) / mass - g_r * std::sin(fpa)
                  + g_lambda * std::cos(fpa) * std::cos(v_azi) - omega * omega * r * std::cos(la)
                  * (std::cos(fpa) * std::cos(v_azi) * std::cos(la) - std::sin(fpa) * std::sin(la));

    double v_azi_dot_term = mass * v * v / r * std::cos(fpa) * std::cos(fpa) * std::sin(v_azi) * std::tan(la)
                        + thrust * std::sin(yaw) - mass * g_lambda * std::sin(v_azi)
                        + mass * omega * omega * r * std::sin(v_azi) * std::sin(la) * std::cos(la)
                        - 2 * mass * omega * v * (std::sin(fpa) * std::cos(v_azi) * std::cos(la) - std::cos(fpa) * std::sin(la));
    double v_azi_dot = v_azi_dot_term / (mass * v * std::cos(fpa));

    double fpa_dot_term = mass * v * v / r * std::cos(fpa) + thrust * std::sin(pitch) * std::cos(yaw)
                        + lift - mass * g_r * std::cos(fpa) - mass * g_lambda * std::sin(fpa) * std::cos(v_azi)
                        + mass * omega * omega * r * std::cos(la) * (std::sin(fpa) * std::cos(v_azi) * std::sin(la) - std::cos(fpa) * std::cos(la))
                        + 2 * mass * omega * v * std::sin(v_azi) * std::cos(la);
    double fpa_dot = fpa_dot_term / (mass * v);

    // Quaternion kinematics.
    Eigen::Quaterniond q(q4, q1, q2, q3); // Eigen ctor is (w,x,y,z); q4 is the scalar part
    Eigen::Quaterniond qdot = QuaternionDerivative(q, wx, wy, wz);

    // Rotational dynamics: torque-free Euler's equation (tau_aero = 0).
    // TODO: tau_aero = 0 (torque-free) -- no panel/geometry data exists yet
    // to compute aerodynamic torque about the center of mass. Once panel
    // geometry and a center-of-pressure model exist, add tau_aero here and
    // change this to J*wdot = tau_aero - w x (J*w).
    Eigen::Vector3d w(wx, wy, wz);
    Eigen::Matrix3d J = Eigen::Vector3d(spacecraft_config_.ixx, spacecraft_config_.iyy, spacecraft_config_.izz).asDiagonal();
    Eigen::Vector3d wdot = J.inverse() * (-w.cross(J * w));

    StateVector xdot;
    xdot << r_dot, la_dot, lo_dot, v_dot, fpa_dot, v_azi_dot,
            qdot.x(), qdot.y(), qdot.z(), qdot.w(),
            wdot(0), wdot(1), wdot(2);
    return xdot;
}

DescentDynamics::TrajectoryHistory DescentDynamics::integrate(
    const DescentState& initial_state,
    const ThrustVectorControlInputs& control_inputs,
    double t0, double tf, double initial_dt, double tol) const
{
    TrajectoryHistory hist;
    StateVector x = toVector(initial_state);
    double t = t0;
    double dt = initial_dt;

    auto record = [&](double t_cur, const StateVector& xv) {
        DescentState s = fromVector(xv);
        hist.t.push_back(t_cur);
        hist.r.push_back(s.r); hist.la.push_back(s.la); hist.lo.push_back(s.lo);
        hist.v.push_back(s.v); hist.fpa.push_back(s.fpa); hist.v_azi.push_back(s.v_azi);
        hist.q1.push_back(s.q1); hist.q2.push_back(s.q2); hist.q3.push_back(s.q3); hist.q4.push_back(s.q4);
        hist.wx.push_back(s.wx); hist.wy.push_back(s.wy); hist.wz.push_back(s.wz);
        double rho = atmosphereDensity(s.r);
        hist.qbar.push_back(0.5 * rho * s.v * s.v);
    };

    record(t, x);

    // Dormand-Prince RK45 (the standard "ode45" method): 5th-order solution
    // advanced, 4th-order embedded solution used only for error estimation.
    // static constexpr double c2 = 1.0 / 5.0, c3 = 3.0 / 10.0, c4 = 4.0 / 5.0, c5 = 8.0 / 9.0;

    static constexpr double a21 = 1.0 / 5.0;
    static constexpr double a31 = 3.0 / 40.0, a32 = 9.0 / 40.0;
    static constexpr double a41 = 44.0 / 45.0, a42 = -56.0 / 15.0, a43 = 32.0 / 9.0;
    static constexpr double a51 = 19372.0 / 6561.0, a52 = -25360.0 / 2187.0, a53 = 64448.0 / 6561.0, a54 = -212.0 / 729.0;
    static constexpr double a61 = 9017.0 / 3168.0, a62 = -355.0 / 33.0, a63 = 46732.0 / 5247.0, a64 = 49.0 / 176.0, a65 = -5103.0 / 18656.0;
    static constexpr double a71 = 35.0 / 384.0, a72 = 0.0, a73 = 500.0 / 1113.0, a74 = 125.0 / 192.0, a75 = -2187.0 / 6784.0, a76 = 11.0 / 84.0;

    // 5th-order weights (propagated solution)
    static constexpr double b1 = 35.0 / 384.0, b2 = 0.0, b3 = 500.0 / 1113.0, b4 = 125.0 / 192.0, b5 = -2187.0 / 6784.0, b6 = 11.0 / 84.0, b7 = 0.0;
    // 4th-order weights (embedded solution, error estimate only)
    static constexpr double bs1 = 5179.0 / 57600.0, bs2 = 0.0, bs3 = 7571.0 / 16695.0, bs4 = 393.0 / 640.0, bs5 = -92097.0 / 339200.0, bs6 = 187.0 / 2100.0, bs7 = 1.0 / 40.0;

    const double dt_min = 1e-3;
    const double dt_max = (tf - t0);
    const double safety = 0.9;
    const double p = 5.0; // order of the propagated solution

    bool warned_floor = false;

    while (t < tf) {
        if (t + dt > tf) dt = tf - t;

        // Floor dt to avoid infinite loops if the dynamics are locally stiff
        bool dt_at_floor = (dt <= dt_min * (1.0 + 1e-9));

        StateVector k1 = derivatives(x, control_inputs);
        StateVector k2 = derivatives(x + dt * (a21 * k1), control_inputs);
        StateVector k3 = derivatives(x + dt * (a31 * k1 + a32 * k2), control_inputs);
        StateVector k4 = derivatives(x + dt * (a41 * k1 + a42 * k2 + a43 * k3), control_inputs);
        StateVector k5 = derivatives(x + dt * (a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4), control_inputs);
        StateVector k6 = derivatives(x + dt * (a61 * k1 + a62 * k2 + a63 * k3 + a64 * k4 + a65 * k5), control_inputs);
        StateVector k7 = derivatives(x + dt * (a71 * k1 + a72 * k2 + a73 * k3 + a74 * k4 + a75 * k5 + a76 * k6), control_inputs);

        StateVector x5 = x + dt * (b1 * k1 + b2 * k2 + b3 * k3 + b4 * k4 + b5 * k5 + b6 * k6 + b7 * k7);
        StateVector x4 = x + dt * (bs1 * k1 + bs2 * k2 + bs3 * k3 + bs4 * k4 + bs5 * k5 + bs6 * k6 + bs7 * k7);

        StateVector err = x5 - x4;

        // RMS-normalized combined absolute/relative error norm (tol doubles
        // as both abs and rel tolerance to match the single-parameter signature).
        double sum_sq = 0.0;
        for (int i = 0; i < 13; ++i) {
            double sc = tol + tol * std::max(std::abs(x(i)), std::abs(x5(i)));
            double e = err(i) / sc;
            sum_sq += e * e;
        }
        double err_norm = std::sqrt(sum_sq / 13.0);

        if (err_norm <= 1.0 || dt_at_floor) {
            if (err_norm > 1.0 && dt_at_floor && !warned_floor) {
                std::cerr << "Warning: RK45 step size hit the floor (" << dt_min
                          << " s) at t=" << t << " without meeting tol; "
                          << "forcing acceptance (dynamics are locally very stiff)."
                          << std::endl;
                warned_floor = true;
            }

            // Accept step.
            t += dt;
            x = x5;

            // Renormalize quaternion sub-block only
            Eigen::Quaterniond q(x(9), x(6), x(7), x(8)); // (w=q4, x=q1, y=q2, z=q3)
            q.normalize();
            x(6) = q.x(); x(7) = q.y(); x(8) = q.z(); x(9) = q.w();

            record(t, x);
        }
        // else: reject step -- do not advance t or x, just shrink dt and retry.

        // Guard against err_norm being NaN/Inf due to a stiff dynamics singularity
        double dt_next;
        if (std::isfinite(err_norm)) {
            double factor = safety * std::pow(err_norm, -1.0 / (p + 1.0));
            factor = std::clamp(factor, 0.2, 5.0);
            dt_next = dt * factor;
        } else {
            dt_next = dt_min;
        }
        dt = std::clamp(dt_next, dt_min, dt_max);
    }

    return hist;
}
