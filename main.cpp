#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <cstdlib>
#include <filesystem>
#include "DescentDynamics.h"
#include "SpacecraftConfig.h"
#include "PlanetConfig.h"
#include "ControlInputs.h"
#include <matplot/matplot.h>

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

    // Set up initial configurations for the simulation
    PlanetConfig earth_config(6371000.0f, 3.986e14f, 1.225f, 8500.0f, 9.80665f, 7.2921159e-5f);
    SpacecraftConfig spacecraft_config(1000.0f, 10.0f, 0.5f, 0.3f, 1000.0f, 1000.0f, 1500.0f);
    DescentDynamics descent_dynamics(earth_config, spacecraft_config);
    ThrustVectorControlInputs control_inputs(0.0f, 0.1f, 0.1f);

    // Set initial entry state for the descent simulation
    constexpr float kEntryInterfaceAltitude = 120000.0f; // m
    DescentState initial_state(
        earth_config.radius + kEntryInterfaceAltitude, 0.0f, 0.0f, 7500.0f, -0.1f, 0.0f, // translational
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

    constexpr double kRadToDeg = 180.0 / 3.14159265358979323846;
    const size_t n = history.t.size();

    // Set up output directory for plots
    std::filesystem::path source_dir = std::filesystem::path(__FILE__).parent_path();
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

    // 8. Heating rate vs time -- TODO: no heating model (e.g. Sutton-Graves)
    // exists in this codebase yet. Add once a stagnation-point heating
    // correlation is implemented.

    } catch (const std::exception& e) {
        // Matplot++ Error Handling: if gnuplot isn't installed or isn't on PATH, matplot++ will throw an exception.
        std::cerr << "Plotting failed: " << e.what() << std::endl;
        std::cerr << "This usually means gnuplot isn't installed or isn't on PATH." << std::endl;
        return 1;
    }

    return 0;
}
