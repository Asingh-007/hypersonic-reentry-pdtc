#ifndef DESCENT_DYNAMICS_H
#define DESCENT_DYNAMICS_H

#include "PlanetConfig.h"
#include "SpacecraftConfig.h"
#include "ControlInputs.h"
#include "QuaternionUtils.h"

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <vector>

struct DescentState{
    // Translational states (spherical coordinates over a rotating planet)
    float r = 0.0f; // radial distance from the center of the Earth
    float la = 0.0f; // latitude
    float lo = 0.0f; // longitude
    float v = 0.0f; // planet-relative speed
    float fpa = 0.0f; // flight path angle
    float v_azi = 0.0f; // velocity azimuth angle

    // Attitude quaternion, scalar-last convention (q1,q2,q3 vector part, q4 scalar part)
    float q1 = 0.0f;
    float q2 = 0.0f;
    float q3 = 0.0f;
    float q4 = 1.0f; // identity quaternion by default

    // Body angular rates in body axes, rad/s
    float wx = 0.0f;
    float wy = 0.0f;
    float wz = 0.0f;

    DescentState(float r, float la, float lo, float v, float fpa, float v_azi,
                 float q1, float q2, float q3, float q4,
                 float wx, float wy, float wz)
        : r(r), la(la), lo(lo), v(v), fpa(fpa), v_azi(v_azi),
          q1(q1), q2(q2), q3(q3), q4(q4), wx(wx), wy(wy), wz(wz) {}
};

class DescentDynamics {
    public:

        // State Vector: 13x1 vector of doubles, representing the full state of the system
        using StateVector = Eigen::Matrix<double, 13, 1>;

        // TrajectoryHistory: structure to hold the history of the trajectory during integration
        struct TrajectoryHistory {
            std::vector<double> t;
            std::vector<double> r, la, lo, v, fpa, v_azi;
            std::vector<double> q1, q2, q3, q4;
            std::vector<double> wx, wy, wz;
            std::vector<double> qbar; // dynamic pressure, derived, for plotting
        };

        DescentDynamics(const PlanetConfig& planet_config, const SpacecraftConfig& spacecraft_config)
            : planet_config_(planet_config), spacecraft_config_(spacecraft_config) {}

        // Full 13-state ODE right-hand side.
        StateVector derivatives(const StateVector& x, const ThrustVectorControlInputs& control_inputs) const;

        // Dormand-Prince RK45 adaptive integrator over [t0, tf]. tol is used
        // as both the absolute and relative error tolerance (RMS-normalized
        // combined error norm), matching MATLAB ode45-style adaptive control.
        TrajectoryHistory integrate(const DescentState& initial_state,
                                     const ThrustVectorControlInputs& control_inputs,
                                     double t0, double tf,
                                     double initial_dt,
                                     double tol) const;

    private:
        PlanetConfig planet_config_;
        SpacecraftConfig spacecraft_config_;

        static StateVector toVector(const DescentState& s);
        static DescentState fromVector(const StateVector& v);

        // Density at radial distance r, delegating to the atmosphere model
        // selected by planet_config_.body (Earth US76 or Mars exponential).
        double atmosphereDensity(double r) const;
};

#endif
