// Standalone Stage A / Stage B demo + plotting executable, analogous to
// aero_table_gen being its own target. Not run automatically, run manually
// to sanity-check Stage A/B's output looks physically reasonable (per the
// SCP guidance planning doc's verification step): altitude/velocity/bank
// angle vs. time with constraint limits overlaid for Stage A; altitude/
// mass/thrust-magnitude/ground-track vs. time for Stage B.
#include "reference_stage_a.h"
#include "reference_stage_b.h"
#include "PlanetConfig.h"
#include "SpacecraftConfig.h"
#include <Eigen/Dense>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <matplot/matplot.h>

constexpr double kPi = 3.14159265358979323846;

int main() {
    {
        const char* existing_path = std::getenv("PATH");
        std::string new_path = "C:\\Program Files\\gnuplot\\bin";
        if (existing_path) { new_path += ";"; new_path += existing_path; }
        _putenv_s("PATH", new_path.c_str());
    }

    std::filesystem::path source_dir = std::filesystem::path(__FILE__).parent_path().parent_path();
    std::filesystem::path output_dir = source_dir / "outputs";
    std::filesystem::create_directories(output_dir);

    PlanetConfig planet_config = PlanetConfig::Earth();
    Eigen::Matrix3d inertia = Eigen::Vector3d(1000.0, 1000.0, 1500.0).asDiagonal();
    double S_ref = kPi * 4.5 * 4.5, L_ref = 9.0;
    Eigen::Vector3d moment_ref(20.0, 0.0, 0.0);
    double nose_radius_m = 0.85;
    std::string aero_table_path = (source_dir / "aero" / "data" / "aero_table.csv").generic_string();
    SpacecraftConfig spacecraft_config(1000.0f, inertia, S_ref, L_ref, moment_ref, nose_radius_m, aero_table_path);

    std::cout << "=== Stage A: entry reference ===" << std::endl;
    // Modest boundary conditions matching the magnitude already validated
    // to converge in GuidanceScpTests.cpp's StageATest as since this demo's
    // purpose is a visual sanity check of a working case, not a convergence-
    // robustness stress test of the simplified single-linearization-per-
    // iteration SCP (a large maneuver like a 30km altitude/3500 m/s change
    // needs either more iterations, a properly-guessed t_scale, or a
    // shrinking trust region to converge
    guidance_scp::StageAState a_initial;
    a_initial.r = planet_config.radius + 60000.0;
    a_initial.v = 4000.0;
    a_initial.fpa = -0.05;

    guidance_scp::StageAState a_terminal;
    a_terminal.r = planet_config.radius + 45000.0;
    a_terminal.la = 0.001; a_terminal.lo = 0.001;
    a_terminal.v = 3000.0;
    a_terminal.fpa = -0.03;
    a_terminal.v_azi = 0.01;

    guidance_scp::StageAConfig a_config;
    a_config.K = 20;
    a_config.t_scale_s = 30.0;
    a_config.max_heat_flux_w_m2 = 5.0e6;
    a_config.max_qbar_pa = 2.0e5;
    a_config.max_load_factor_g = 30.0;

    guidance_scp::StageAResult a_result =
        guidance_scp::solveStageA(a_initial, a_terminal, planet_config, spacecraft_config, a_config);

    std::cout << "converged = " << (a_result.converged ? "true" : "false") << std::endl;
    std::cout << "iterations run = " << a_result.max_nu_per_iter.size() << std::endl;
    for (size_t i = 0; i < a_result.max_nu_per_iter.size(); ++i) {
        std::cout << "  iter " << i << ": max||nu||=" << a_result.max_nu_per_iter[i]
                  << " max(eta)=" << a_result.max_eta_per_iter[i] << std::endl;
    }

    try {
        const auto& h = a_result.history;
        std::vector<double> altitude(h.r.size());
        for (size_t i = 0; i < h.r.size(); ++i) altitude[i] = h.r[i] - planet_config.radius;

        auto figA1 = matplot::figure(true);
        matplot::plot(h.t, altitude);
        matplot::xlabel("Time (s)"); matplot::ylabel("Altitude (m)");
        matplot::title("Stage A: Altitude vs Time");
        matplot::save((output_dir / "stage_a_altitude_vs_time.png").generic_string());

        auto figA2 = matplot::figure(true);
        matplot::plot(h.t, h.v);
        matplot::xlabel("Time (s)"); matplot::ylabel("Velocity (m/s)");
        matplot::title("Stage A: Velocity vs Time");
        matplot::save((output_dir / "stage_a_velocity_vs_time.png").generic_string());

        std::vector<double> sigma_deg(h.sigma.size());
        for (size_t i = 0; i < h.sigma.size(); ++i) sigma_deg[i] = h.sigma[i] * 180.0 / kPi;
        auto figA3 = matplot::figure(true);
        matplot::plot(h.t, sigma_deg);
        matplot::xlabel("Time (s)"); matplot::ylabel("Bank angle (deg)");
        matplot::title("Stage A: Bank Angle vs Time");
        matplot::save((output_dir / "stage_a_bank_angle_vs_time.png").generic_string());

        auto figA4 = matplot::figure(true);
        matplot::hold(true);
        matplot::plot(h.t, h.heat_flux_conv)->line_width(2);
        matplot::plot(h.t, std::vector<double>(h.t.size(), a_config.max_heat_flux_w_m2))->line_style("--");
        matplot::legend(std::vector<std::string>{"heat flux", "limit"});
        matplot::xlabel("Time (s)"); matplot::ylabel("Heat Flux (W/m^2)");
        matplot::title("Stage A: Heat Flux vs Time");
        matplot::save((output_dir / "stage_a_heat_flux_vs_time.png").generic_string());
    } catch (const std::exception& e) {
        std::cerr << "Stage A plotting failed: " << e.what() << std::endl;
    }

    std::cout << "\n=== Stage B: descent reference ===" << std::endl;
    guidance_scp::StageBState b_initial;
    b_initial.m = 8000.0;
    b_initial.rz = 2000.0;
    b_initial.vx = 50.0; b_initial.vz = -80.0;

    guidance_scp::StageBState b_terminal;  // all zero: touchdown at the tower-catch origin

    guidance_scp::StageBConfig b_config;
    b_config.K = 30;
    b_config.Isp_s = 330.0;
    b_config.Tmin_N = 20000.0;
    b_config.Tmax_N = 100000.0;
    b_config.glideslope_deg = 80.0;
    b_config.tf_min_s = 5.0;
    b_config.tf_max_s = 90.0;

    guidance_scp::StageBResult b_result = guidance_scp::solveStageB(b_initial, b_terminal, planet_config, b_config);

    std::cout << "solved = " << (b_result.solved ? "true" : "false") << std::endl;
    if (b_result.solved) {
        std::cout << "tf = " << b_result.tf_s << " s" << std::endl;
        double propellant_frac = (b_initial.m - b_result.history.m.back()) / b_initial.m;
        std::cout << "propellant mass fraction = " << propellant_frac << std::endl;

        try {
            const auto& h = b_result.history;
            auto figB1 = matplot::figure(true);
            matplot::plot(h.t, h.rz);
            matplot::xlabel("Time (s)"); matplot::ylabel("Altitude (m)");
            matplot::title("Stage B: Altitude vs Time");
            matplot::save((output_dir / "stage_b_altitude_vs_time.png").generic_string());

            auto figB2 = matplot::figure(true);
            matplot::plot(h.t, h.m);
            matplot::xlabel("Time (s)"); matplot::ylabel("Mass (kg)");
            matplot::title("Stage B: Mass vs Time");
            matplot::save((output_dir / "stage_b_mass_vs_time.png").generic_string());

            std::vector<double> thrust_mag(h.Tx.size());
            for (size_t i = 0; i < h.Tx.size(); ++i) {
                thrust_mag[i] = std::sqrt(h.Tx[i] * h.Tx[i] + h.Ty[i] * h.Ty[i] + h.Tz[i] * h.Tz[i]);
            }
            auto figB3 = matplot::figure(true);
            matplot::plot(std::vector<double>(h.t.begin(), h.t.end() - 1), thrust_mag);
            matplot::xlabel("Time (s)"); matplot::ylabel("Thrust magnitude (N)");
            matplot::title("Stage B: Thrust Magnitude vs Time");
            matplot::save((output_dir / "stage_b_thrust_vs_time.png").generic_string());

            auto figB4 = matplot::figure(true);
            matplot::plot(h.rx, h.ry);
            matplot::xlabel("East (m)"); matplot::ylabel("North (m)");
            matplot::title("Stage B: Ground Track");
            matplot::save((output_dir / "stage_b_ground_track.png").generic_string());
        } catch (const std::exception& e) {
            std::cerr << "Stage B plotting failed: " << e.what() << std::endl;
        }
    }

    return 0;
}
