// Standalone full 6DOF SCP loop demo + plotting executable, mirroring
// StageDemoMain.cpp's exact structure. Not run automatically -- run
// manually to visually sanity-check the full loop's output looks
// physically reasonable, per the full-loop plan's verification section.
#include "full_scp_loop.h"
#include "full_loop_phase1.h"
#include "full_loop_phase2.h"
#include "PlanetConfig.h"
#include "SpacecraftConfig.h"
#include "FlapActuatorConfig.h"
#include "GimbalActuatorConfig.h"
#include "DescentDynamics.h"
#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cmath>
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
    FlapActuatorConfig flap_config;
    GimbalActuatorConfig gimbal_config;

    std::cout << "=== Stage A: entry reference (prerequisite) ===" << std::endl;
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
    // 30g, not 20g -- see StageDemoMain.cpp's identical comment: 20g makes
    // this exact terminal boundary condition infeasible against the real
    // aero table (natural load factor there is ~23.5g).
    a_config.max_load_factor_g = 30.0;

    guidance_scp::StageAResult a_result =
        guidance_scp::solveStageA(a_initial, a_terminal, planet_config, spacecraft_config, a_config);
    std::cout << "Stage A converged = " << (a_result.converged ? "true" : "false") << std::endl;

    std::cout << "\n=== Stage B: descent reference (prerequisite) ===" << std::endl;
    guidance_scp::StageBState b_initial;
    b_initial.m = 8000.0;
    b_initial.rz = 2000.0;
    b_initial.vx = 50.0; b_initial.vz = -80.0;

    guidance_scp::StageBState b_terminal;  // all zero: touchdown at the tower-catch origin

    guidance_scp::StageBConfig b_config;
    b_config.K = 20;
    b_config.Isp_s = 330.0;
    b_config.Tmin_N = 20000.0;
    b_config.Tmax_N = 100000.0;
    b_config.glideslope_deg = 80.0;
    b_config.tf_min_s = 5.0;
    b_config.tf_max_s = 90.0;

    guidance_scp::StageBResult b_result = guidance_scp::solveStageB(b_initial, b_terminal, planet_config, b_config);
    std::cout << "Stage B solved = " << (b_result.solved ? "true" : "false") << std::endl;

    if (!a_result.converged) {
        std::cerr << "NOTE: Stage A did not converge against the real (richer) aero table within its "
                  << "iteration budget -- a known, already-documented limitation of its simplified "
                  << "SCP (see stage_demo.exe's own prior run). Proceeding anyway with its last "
                  << "attempted reference -- the full loop's own re-linearization/virtual-control "
                  << "machinery is exactly what is being exercised here, on a deliberately imperfect "
                  << "bootstrap." << std::endl;
    }
    if (!b_result.solved) {
        std::cerr << "Stage B prerequisite did not solve -- aborting full-loop demo (Stage B is a "
                  << "single convex solve; without it there is no descent reference to stitch)." << std::endl;
        return 1;
    }

    std::cout << "\n=== Stage C: stitching + full 6DOF SCP loop ===" << std::endl;
    const int K1 = 16, K2 = 16;
    guidance_scp::StageCResult stitched =
        guidance_scp::stitchStageAAndB(a_result, b_result, K1, K2, planet_config, spacecraft_config);

    guidance_scp::Phase1ToPhase2Frame frame;
    // origin_r/la/lo is Phase 2's ENU origin -- the intended TOWER-CATCH
    // point (see full_loop_transition.h's doc comment), i.e. ground level,
    // NOT Stage A's 45000m entry-interface handoff altitude. la/lo are kept
    // at a_terminal's values as a simplifying approximation (tower assumed
    // roughly below the entry corridor's ground track).
    frame.origin_r = planet_config.radius;
    frame.origin_la = a_terminal.la;
    frame.origin_lo = a_terminal.lo;

    Eigen::VectorXd x1_initial = stitched.x1_ref.front();
    Eigen::VectorXd x2_terminal_target = Eigen::VectorXd::Zero(guidance_scp::kPhase2StateDim);
    x2_terminal_target(10) = 1.0;  // identity quaternion at touchdown

    guidance_scp::FullLoopConfig full_loop_config;
    full_loop_config.Isp_s = b_config.Isp_s;
    full_loop_config.Tmin_N = b_config.Tmin_N;
    full_loop_config.Tmax_N = b_config.Tmax_N;
    full_loop_config.glideslope_deg = b_config.glideslope_deg;
    full_loop_config.m_wet_at_handoff_kg = b_initial.m;
    // Stage C's bootstrap Phase-2 attitude and this demo's identity terminal
    // target are both arbitrary defaults, generously separated here for the
    // same reason documented in FullLoopIntegrationTests.cpp.
    full_loop_config.terminal_attitude_error_max_rad = kPi;

    guidance_scp::FullLoopResult result = guidance_scp::solveFullLoop(
        stitched, a_config.t_scale_s, b_result.tf_s, x1_initial, x2_terminal_target, frame,
        planet_config, spacecraft_config, flap_config, gimbal_config, full_loop_config,
        /*max_iters=*/20, /*eps_dyn=*/1.0, /*eps_tr=*/1.0);

    std::cout << "Full loop converged = " << (result.converged ? "true" : "false") << std::endl;
    std::cout << "iterations run = " << result.max_nu_per_iter.size() << std::endl;
    for (size_t i = 0; i < result.max_nu_per_iter.size(); ++i) {
        std::cout << "  iter " << i << ": max||nu||=" << result.max_nu_per_iter[i]
                  << " max(eta)=" << result.max_eta_per_iter[i] << std::endl;
    }

    if (!result.x2_ref.empty()) {
        double propellant_frac = (b_initial.m - result.x2_ref.back()(0)) / b_initial.m;
        std::cout << "propellant mass fraction = " << propellant_frac << std::endl;

        Eigen::Quaterniond q_final(result.x2_ref.back()(10), result.x2_ref.back()(7),
                                     result.x2_ref.back()(8), result.x2_ref.back()(9));
        Eigen::Quaterniond q_target(x2_terminal_target(10), x2_terminal_target(7),
                                     x2_terminal_target(8), x2_terminal_target(9));
        double att_err_deg = 2.0 * std::acos(std::clamp(std::abs(q_final.dot(q_target)), 0.0, 1.0)) * 180.0 / kPi;
        std::cout << "attitude error at tower catch = " << att_err_deg << " deg" << std::endl;
    }

    try {
        auto fig1 = matplot::figure(true);
        matplot::hold(true);
        std::vector<double> iters(result.max_nu_per_iter.size());
        for (size_t i = 0; i < iters.size(); ++i) iters[i] = static_cast<double>(i);
        matplot::semilogy(iters, result.max_nu_per_iter)->line_width(2);
        matplot::semilogy(iters, result.max_eta_per_iter)->line_width(2);
        matplot::legend(std::vector<std::string>{"max||nu||", "max(eta)"});
        matplot::xlabel("Outer iteration"); matplot::ylabel("Convergence diagnostic (log scale)");
        matplot::title("Full Loop: Convergence Diagnostics");
        matplot::save((output_dir / "full_loop_convergence.png").generic_string());

        double dt1 = a_config.t_scale_s / (K1 - 1);
        double dt2 = b_result.tf_s / (K2 - 1);
        std::vector<double> t1(K1), alt1(K1), mach1(K1);
        std::vector<std::array<double, 4>> flap1(K1);
        for (int k = 0; k < K1; ++k) {
            t1[k] = k * dt1;
            alt1[k] = result.x1_ref[k](0) - planet_config.radius;
            mach1[k] = result.x1_ref[k](3) / DescentDynamics::speedOfSound(result.x1_ref[k](0), planet_config);
        }
        std::vector<double> t2(K2), alt2(K2), thrust2(K2 - 1);
        for (int k = 0; k < K2; ++k) {
            t2[k] = a_config.t_scale_s + k * dt2;
            // x2_ref[k](3) is rz in Phase 2's LOCAL ENU frame, whose origin
            // sits at frame.origin_r (the tower-catch point), offset by
            // (frame.origin_r - planet_config.radius) so alt2 shares alt1's
            // "height above the surface" datum in general (this evaluates
            // to 0 now that origin_r is ground level, but stays correct if
            // that ever changes). Without this offset, the two series would
            // disagree by ~frame.origin_r's altitude and produce
            // a spurious cliff at the phase boundary when concatenated below.
            alt2[k] = (frame.origin_r - planet_config.radius) + result.x2_ref[k](3);
        }
        for (int k = 0; k < K2 - 1; ++k) {
            thrust2[k] = result.u2_ref[k](0);
        }

        std::vector<double> t_all = t1;
        t_all.insert(t_all.end(), t2.begin(), t2.end());
        std::vector<double> alt_all = alt1;
        alt_all.insert(alt_all.end(), alt2.begin(), alt2.end());

        auto fig2 = matplot::figure(true);
        matplot::plot(t_all, alt_all)->line_width(2);
        matplot::xlabel("Time (s)"); matplot::ylabel("Altitude (m)");
        matplot::title("Full Loop: Altitude vs Time (Phase 1 + Phase 2)");
        matplot::save((output_dir / "full_loop_altitude_vs_time.png").generic_string());

        auto fig3 = matplot::figure(true);
        matplot::plot(t1, mach1)->line_width(2);
        matplot::xlabel("Time (s)"); matplot::ylabel("Mach");
        matplot::title("Full Loop: Phase 1 Mach vs Time");
        matplot::save((output_dir / "full_loop_mach_vs_time.png").generic_string());

        auto fig4 = matplot::figure(true);
        matplot::hold(true);
        for (int i = 0; i < 4; ++i) {
            std::vector<double> flap_deg(K1);
            for (int k = 0; k < K1; ++k) flap_deg[k] = result.x1_ref[k](13 + i) * 180.0 / kPi;
            matplot::plot(t1, flap_deg)->line_width(1.5);
        }
        matplot::legend(std::vector<std::string>{"d1", "d2", "d3", "d4"});
        matplot::xlabel("Time (s)"); matplot::ylabel("Flap deflection (deg)");
        matplot::title("Full Loop: Phase 1 Flap Deflections vs Time");
        matplot::save((output_dir / "full_loop_flap_deflections_vs_time.png").generic_string());

        auto fig5 = matplot::figure(true);
        matplot::plot(std::vector<double>(t2.begin(), t2.end() - 1), thrust2)->line_width(2);
        matplot::xlabel("Time (s)"); matplot::ylabel("Thrust (N)");
        matplot::title("Full Loop: Phase 2 Thrust Magnitude vs Time");
        matplot::save((output_dir / "full_loop_thrust_vs_time.png").generic_string());
    } catch (const std::exception& e) {
        std::cerr << "Full loop plotting failed: " << e.what() << std::endl;
    }

    return 0;
}
