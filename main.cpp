#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <cstdlib>
#include <filesystem>
#include <Eigen/Dense>
#include "DescentDynamics.h"
#include "SpacecraftConfig.h"
#include "PlanetConfig.h"
#include "ControlInputs.h"
#include "HeatingLoadModel.h"
#include "GuidanceConstraints.h"
#include <matplot/matplot.h>

// MSVC only defines M_PI when _USE_MATH_DEFINES is set before the first
// <cmath>/<math.h> include anywhere in the translation unit, which is
// fragile to guarantee across header inclusion order -- use an explicit
// local constant instead.
constexpr double kPi = 3.14159265358979323846;

int main() {
    // Append gnuplot's bin directory to PATH so matplot++ can find it
    {
        const char* existing_path = std::getenv("PATH");
        std::string new_path = "C:\\Program Files\\gnuplot\\bin";
        if (existing_path) {
            new_path += ";";
            new_path += existing_path;
        }
        _putenv_s("PATH", new_path.c_str());
    }

    // Resolved once, up front, since both the aero table path below and the
    // plot output directory later need a location independent of the
    // process' runtime working directory (see the outputs/ note further down).
    std::filesystem::path source_dir = std::filesystem::path(__FILE__).parent_path();

    // Set up initial configurations for the simulation
    PlanetConfig earth_config = PlanetConfig::Earth();
    // Inertia/S_ref/L_ref/moment_ref are PLACEHOLDERS matching
    // aero_model::testutil::makeCylinderNoseFlapBody()'s default geometry --
    // see SpacecraftConfig.h. aero_table.csv is generated offline by
    // aero_table_gen (see aero/GenerateAeroTable.cpp); run that once before
    // running this if the CSV doesn't exist yet.
    Eigen::Matrix3d inertia = Eigen::Vector3d(1000.0, 1000.0, 1500.0).asDiagonal();
    double S_ref = kPi * 4.5 * 4.5;
    double L_ref = 9.0;
    Eigen::Vector3d moment_ref(20.0, 0.0, 0.0);
    double nose_radius_m = 0.85; // PLACEHOLDER: ~0.18x real body_radius=4.633m
                                  // (aero/data/reference_quantities.csv) -- see
                                  // SpacecraftConfig.h nose_radius_m comment.
    std::string aero_table_path = (source_dir / "aero" / "data" / "aero_table.csv").generic_string();
    SpacecraftConfig spacecraft_config(1000.0f, inertia, S_ref, L_ref, moment_ref, nose_radius_m, aero_table_path);
    DescentDynamics descent_dynamics(earth_config, spacecraft_config);
    ThrustVectorControlInputs control_inputs(0.0f, 0.1f, 0.1f);

    // Set initial entry state for the descent simulation
    constexpr float kEntryInterfaceAltitude = 120000.0f; // m
    DescentState initial_state(
        earth_config.radius + kEntryInterfaceAltitude, 0.0f, 0.0f, 7800.0f, -0.1f, 0.0f, // translational
        0.0f, 0.0f, 0.0f, 1.0f,                        // q1,q2,q3,q4 identity
        0.01f, 0.01f, 0.0f                              // wx,wy,wz small placeholder rates
    );

    // Set simulation time parameters
    double t0 = 0.0;
    double tf = 100.0; 
    double initial_dt = 1.0;
    double tol = 1e-6;

    auto history = descent_dynamics.integrate(initial_state, control_inputs, t0, tf, initial_dt, tol);

    std::cout << "Integration complete: " << history.t.size() << " recorded points." << std::endl;
    std::cout << "Final: r = " << history.r.back() << ", v = " << history.v.back() << std::endl;

    GuidanceConstraints constraints;
    ConstraintViolationReport violation_report = checkConstraints(history, constraints);
    if (violation_report.hasViolations()) {
        std::cout << "Constraint check: " << violation_report.violations.size()
                  << " violation(s) found:" << std::endl;
        for (const auto& v : violation_report.violations) {
            std::cout << "  [" << v.constraint_name << "] at t=" << v.time_s
                      << "s (index " << v.index << "): value=" << v.value
                      << ", limit=" << v.limit << ", margin=" << v.margin << std::endl;
        }
    } else {
        std::cout << "Constraint check: no violations (within max_heat_flux="
                  << constraints.max_heat_flux_w_m2 << " W/m^2, max_qbar="
                  << constraints.max_qbar_pa << " Pa, max_load_factor="
                  << constraints.max_load_factor_g << " g)." << std::endl;
    }

    constexpr double kRadToDeg = 180.0 / kPi;
    const size_t n = history.t.size();

    // Set up output directory for plots
    std::filesystem::path output_dir = source_dir / "outputs";
    std::filesystem::create_directories(output_dir);

    try {

    // 1. Altitude vs time
    std::vector<double> altitude(n);
    for (size_t i = 0; i < n; ++i) altitude[i] = history.r[i] - earth_config.radius;
    auto fig1 = matplot::figure(true);
    matplot::plot(history.t, altitude);
    matplot::xlabel("Time (s)");
    matplot::ylabel("Altitude (m)");
    matplot::title("Altitude vs Time");
    matplot::save((output_dir / "altitude_vs_time.png").generic_string());

    // 2. Velocity vs time
    auto fig2 = matplot::figure(true);
    matplot::plot(history.t, history.v);
    matplot::xlabel("Time (s)");
    matplot::ylabel("Velocity (m/s)");
    matplot::title("Velocity vs Time");
    matplot::save((output_dir / "velocity_vs_time.png").generic_string());

    // 3. Flight-path angle (deg) vs time
    std::vector<double> fpa_deg(n);
    for (size_t i = 0; i < n; ++i) fpa_deg[i] = history.fpa[i] * kRadToDeg;
    auto fig3 = matplot::figure(true);
    matplot::plot(history.t, fpa_deg);
    matplot::xlabel("Time (s)");
    matplot::ylabel("Flight Path Angle (deg)");
    matplot::title("Flight Path Angle vs Time");
    matplot::save((output_dir / "fpa_vs_time.png").generic_string());

    // 4. Ground track (longitude vs latitude, deg)
    std::vector<double> la_deg(n), lo_deg(n);
    for (size_t i = 0; i < n; ++i) {
        la_deg[i] = history.la[i] * kRadToDeg;
        lo_deg[i] = history.lo[i] * kRadToDeg;
    }
    auto fig4 = matplot::figure(true);
    matplot::plot(lo_deg, la_deg);
    matplot::xlabel("Longitude (deg)");
    matplot::ylabel("Latitude (deg)");
    matplot::title("Ground Track");
    matplot::save((output_dir / "ground_track.png").generic_string());

    // 5. Dynamic pressure vs time
    auto fig5 = matplot::figure(true);
    matplot::plot(history.t, history.qbar);
    matplot::xlabel("Time (s)");
    matplot::ylabel("Dynamic Pressure (Pa)");
    matplot::title("Dynamic Pressure vs Time");
    matplot::save((output_dir / "dynamic_pressure_vs_time.png").generic_string());

    // 6. Angular velocity components vs time, one figure with legend
    auto fig6 = matplot::figure(true);
    matplot::hold(true);
    matplot::plot(history.t, history.wx)->line_width(2);
    matplot::plot(history.t, history.wy)->line_width(2);
    matplot::plot(history.t, history.wz)->line_width(2);
    matplot::legend(std::vector<std::string>{"wx", "wy", "wz"});
    matplot::xlabel("Time (s)");
    matplot::ylabel("Angular velocity (rad/s)");
    matplot::title("Body Angular Rates vs Time");
    matplot::save((output_dir / "angular_rates_vs_time.png").generic_string());

    // 7. Quaternion norm vs time (numerical-drift sanity check)
    std::vector<double> qnorm(n);
    for (size_t i = 0; i < n; ++i) {
        qnorm[i] = std::sqrt(history.q1[i] * history.q1[i] + history.q2[i] * history.q2[i]
                            + history.q3[i] * history.q3[i] + history.q4[i] * history.q4[i]);
    }
    auto fig7 = matplot::figure(true);
    matplot::plot(history.t, qnorm);
    matplot::xlabel("Time (s)");
    matplot::ylabel("|q|");
    matplot::title("Quaternion Norm vs Time (drift check)");
    matplot::save((output_dir / "quaternion_norm_vs_time.png").generic_string());

    // 8. Heat flux vs time (convective, radiative, total -- one figure with legend)
    auto fig8 = matplot::figure(true);
    matplot::hold(true);
    matplot::plot(history.t, history.heat_flux_conv)->line_width(2);
    matplot::plot(history.t, history.heat_flux_rad)->line_width(2);
    matplot::plot(history.t, history.heat_flux_total)->line_width(2);
    matplot::legend(std::vector<std::string>{"convective (Sutton-Graves)", "radiative (Tauber-Sutton)", "total"});
    matplot::xlabel("Time (s)");
    matplot::ylabel("Stagnation-Point Heat Flux (W/m^2)");
    matplot::title("Heat Flux vs Time");
    matplot::save((output_dir / "heat_flux_vs_time.png").generic_string());

    // 9. Cumulative heat load vs time
    auto fig9 = matplot::figure(true);
    matplot::plot(history.t, history.heat_load);
    matplot::xlabel("Time (s)");
    matplot::ylabel("Cumulative Heat Load (J/m^2)");
    matplot::title("Cumulative Heat Load vs Time");
    matplot::save((output_dir / "heat_load_vs_time.png").generic_string());

    // 10. Load factor vs time
    auto fig10 = matplot::figure(true);
    matplot::plot(history.t, history.load_factor);
    matplot::xlabel("Time (s)");
    matplot::ylabel("Load Factor (g)");
    matplot::title("Load Factor vs Time");
    matplot::save((output_dir / "load_factor_vs_time.png").generic_string());

    } catch (const std::exception& e) {
        // Matplot++ Error Handling: if gnuplot isn't installed or isn't on PATH, matplot++ will throw an exception.
        std::cerr << "Plotting failed: " << e.what() << std::endl;
        std::cerr << "This usually means gnuplot isn't installed or isn't on PATH." << std::endl;
        return 1;
    }

    return 0;
}
