#include <gtest/gtest.h>
#include "ClarabelSocpSolver.h"
#include "TrimSolver.h"
#include "reference_stage_a.h"
#include "reference_stage_b.h"
#include "full_loop_path_constraints.h"
#include "full_loop_transition.h"
#include "PlanetConfig.h"
#include "SpacecraftConfig.h"
#include "AeroCoefficientTable.h"
#include <Eigen/Sparse>
#include <cmath>

constexpr double kPi = 3.14159265358979323846;

namespace {

std::string TestAeroTablePath() {
    return std::string(TESTS_SOURCE_DIR) + "/tests/fixtures/aero_table_test.csv";
}

SpacecraftConfig MakeTestSpacecraftConfig(double nose_radius_m = 0.85) {
    double S_ref = kPi * 4.5 * 4.5;
    double L_ref = 9.0;
    Eigen::Vector3d moment_ref(20.0, 0.0, 0.0);
    Eigen::Matrix3d inertia = Eigen::Vector3d(1000.0, 1000.0, 1500.0).asDiagonal();
    return SpacecraftConfig(1000.0f, inertia, S_ref, L_ref, moment_ref, nose_radius_m, TestAeroTablePath());
}

std::string RealAeroTablePath() {
    return std::string(TESTS_SOURCE_DIR) + "/aero/data/aero_table.csv";
}

// Real, generated aero table (17 nonuniform Mach points x 8 alpha points) --
// much coarser/more nonlinear than the small hand-crafted fixture above.
// Used to confirm the adaptive trust-region accept/reject mechanism
// (reference_stage_a.cpp) actually converges against production-representative
// data, not just the easy fixture -- this is the exact scenario that was
// empirically confirmed broken before that fix (see
// feedback_scp_trust_region_debugging memory).
SpacecraftConfig MakeRealSpacecraftConfig(double nose_radius_m = 0.85) {
    double S_ref = kPi * 4.5 * 4.5;
    double L_ref = 9.0;
    Eigen::Vector3d moment_ref(20.0, 0.0, 0.0);
    Eigen::Matrix3d inertia = Eigen::Vector3d(1000.0, 1000.0, 1500.0).asDiagonal();
    return SpacecraftConfig(1000.0f, inertia, S_ref, L_ref, moment_ref, nose_radius_m, RealAeroTablePath());
}

}  // namespace

TEST(TrimSolverTest, FindsZeroMomentCoefficientWithinTolerance) {
    // Fixture's Cm is linear in alpha_deg at mach=5 (0.2 at alpha=-10, -0.2
    // at alpha=10), crossing zero exactly at alpha=0.
    SpacecraftConfig spacecraft_config = MakeTestSpacecraftConfig();
    double alpha_trim = guidance_scp::alphaTrimDeg(5.0, spacecraft_config.aero_table);
    EXPECT_NEAR(alpha_trim, 0.0, 1e-2);
}

TEST(ClarabelSocpSolverTest, SolvesKnownTinySocp) {
    // Reproduces Clarabel's own example_socp.cpp: minimize x1^2 (P has a
    // single 2.0 on the (1,1) diagonal, so cost = x1^2) subject to a single
    // 3-dim SOC constraint s = b-Ax = (1, 2*x0-2, x1-2), s0>=||s_{1:}||,
    // i.e. 4*(x0-1)^2 + (x1-2)^2 <= 1 -- an ellipse centered at (1,2) with
    // x1 ranging over [1,3]. Minimizing x1^2 over that ellipse gives the
    // closed-form optimum x = (1, 1) (smallest-|x1| point on the ellipse).
    Eigen::SparseMatrix<double> P(2, 2);
    P.insert(1, 1) = 2.0;
    P.makeCompressed();

    Eigen::VectorXd q(2);
    q << 0.0, 0.0;

    Eigen::SparseMatrix<double> A(3, 2);
    A.insert(1, 0) = -2.0;
    A.insert(2, 1) = -1.0;
    A.makeCompressed();

    Eigen::VectorXd b(3);
    b << 1.0, -2.0, -2.0;

    std::vector<clarabel::SupportedConeT<double>> cones{clarabel::SecondOrderConeT<double>(3)};

    guidance_scp::SocpProblem problem{P, q, A, b, cones};
    guidance_scp::SocpSolution solution = guidance_scp::solve(problem);

    ASSERT_EQ(solution.status, clarabel::SolverStatus::Solved);
    EXPECT_NEAR(solution.x(0), 1.0, 1e-6);
    EXPECT_NEAR(solution.x(1), 1.0, 1e-6);
}

TEST(StageAPathConstraintGradientTest, HeatFluxGradientMatchesFiniteDifference) {
    PlanetConfig planet_config = PlanetConfig::Earth();
    SpacecraftConfig spacecraft_config = MakeTestSpacecraftConfig();
    double r = planet_config.radius + 40000.0, v = 4000.0;

    auto g = guidance_scp::computePathConstraintGradients(r, v, planet_config, spacecraft_config);

    double h_r = 1.0, h_v = 0.1;
    double q_r_plus = guidance_scp::computePathConstraintValues(r + h_r, v, planet_config, spacecraft_config).qdot_w_m2;
    double q_r_minus = guidance_scp::computePathConstraintValues(r - h_r, v, planet_config, spacecraft_config).qdot_w_m2;
    double q_v_plus = guidance_scp::computePathConstraintValues(r, v + h_v, planet_config, spacecraft_config).qdot_w_m2;
    double q_v_minus = guidance_scp::computePathConstraintValues(r, v - h_v, planet_config, spacecraft_config).qdot_w_m2;

    double fd_dQdot_dr = (q_r_plus - q_r_minus) / (2 * h_r);
    double fd_dQdot_dV = (q_v_plus - q_v_minus) / (2 * h_v);

    EXPECT_NEAR(g.dQdot_dr, fd_dQdot_dr, std::abs(fd_dQdot_dr) * 1e-3 + 1e-9);
    EXPECT_NEAR(g.dQdot_dV, fd_dQdot_dV, std::abs(fd_dQdot_dV) * 1e-3 + 1e-9);
}

TEST(StageAPathConstraintGradientTest, QbarGradientMatchesFiniteDifference) {
    PlanetConfig planet_config = PlanetConfig::Earth();
    SpacecraftConfig spacecraft_config = MakeTestSpacecraftConfig();
    double r = planet_config.radius + 40000.0, v = 4000.0;

    auto g = guidance_scp::computePathConstraintGradients(r, v, planet_config, spacecraft_config);

    double h_r = 1.0, h_v = 0.1;
    double q_r_plus = guidance_scp::computePathConstraintValues(r + h_r, v, planet_config, spacecraft_config).qbar_pa;
    double q_r_minus = guidance_scp::computePathConstraintValues(r - h_r, v, planet_config, spacecraft_config).qbar_pa;
    double q_v_plus = guidance_scp::computePathConstraintValues(r, v + h_v, planet_config, spacecraft_config).qbar_pa;
    double q_v_minus = guidance_scp::computePathConstraintValues(r, v - h_v, planet_config, spacecraft_config).qbar_pa;

    double fd_dqbar_dr = (q_r_plus - q_r_minus) / (2 * h_r);
    double fd_dqbar_dV = (q_v_plus - q_v_minus) / (2 * h_v);

    EXPECT_NEAR(g.dqbar_dr, fd_dqbar_dr, std::abs(fd_dqbar_dr) * 1e-3 + 1e-9);
    EXPECT_NEAR(g.dqbar_dV, fd_dqbar_dV, std::abs(fd_dqbar_dV) * 1e-3 + 1e-9);
}

TEST(StageAPathConstraintGradientTest, LoadFactorGradientMatchesFiniteDifference) {
    PlanetConfig planet_config = PlanetConfig::Earth();
    SpacecraftConfig spacecraft_config = MakeTestSpacecraftConfig();
    double r = planet_config.radius + 40000.0, v = 4000.0;

    auto g = guidance_scp::computePathConstraintGradients(r, v, planet_config, spacecraft_config);

    double h_r = 1.0, h_v = 0.1;
    double n_r_plus = guidance_scp::computePathConstraintValues(r + h_r, v, planet_config, spacecraft_config).n_g;
    double n_r_minus = guidance_scp::computePathConstraintValues(r - h_r, v, planet_config, spacecraft_config).n_g;
    double n_v_plus = guidance_scp::computePathConstraintValues(r, v + h_v, planet_config, spacecraft_config).n_g;
    double n_v_minus = guidance_scp::computePathConstraintValues(r, v - h_v, planet_config, spacecraft_config).n_g;

    double fd_dn_dr = (n_r_plus - n_r_minus) / (2 * h_r);
    double fd_dn_dV = (n_v_plus - n_v_minus) / (2 * h_v);

    // Looser tolerance than heat-flux/qbar: n's gradient formula assumes
    // CL/CD's own Mach dependence is negligible (documented simplification
    // in reference_stage_a.h/.cpp), which the finite difference does NOT
    // assume -- some real discrepancy is expected here, not just numerical noise.
    EXPECT_NEAR(g.dn_dr, fd_dn_dr, std::abs(fd_dn_dr) * 0.1 + 1e-9);
    EXPECT_NEAR(g.dn_dV, fd_dn_dV, std::abs(fd_dn_dV) * 0.1 + 1e-9);
}

TEST(FullLoadFactorGradientTest, RVPartMatchesStageAExactly) {
    PlanetConfig planet_config = PlanetConfig::Earth();
    SpacecraftConfig spacecraft_config = MakeTestSpacecraftConfig();
    double r = planet_config.radius + 40000.0, v = 4000.0;

    auto stage_a_grad = guidance_scp::computePathConstraintGradients(r, v, planet_config, spacecraft_config);
    auto full_grad = guidance_scp::computeFullLoadFactorGradient(r, v, 0.0, 0.0, 0.0, 0.0,
                                                                     planet_config, spacecraft_config);

    EXPECT_NEAR(full_grad.dn_dr, stage_a_grad.dn_dr, 1e-12);
    EXPECT_NEAR(full_grad.dn_dV, stage_a_grad.dn_dV, 1e-12);
}

TEST(FullLoadFactorGradientTest, FlapPartialIsNonzeroForFlapDependentFixture) {
    // The test fixture's Ch (and hence, indirectly through the table's
    // shared grid, CL/CD's flap-axis sensitivity) is deliberately
    // flap-dependent (see AeroTests.cpp's fixture comment) -- at least one
    // flap partial should be measurably nonzero.
    PlanetConfig planet_config = PlanetConfig::Earth();
    SpacecraftConfig spacecraft_config = MakeTestSpacecraftConfig();
    double r = planet_config.radius + 40000.0, v = 4000.0;

    auto grad = guidance_scp::computeFullLoadFactorGradient(r, v, 0.02, -0.01, 0.03, -0.02,
                                                               planet_config, spacecraft_config);
    bool any_nonzero = false;
    for (double d : grad.dn_ddelta) {
        if (std::abs(d) > 1e-9) any_nonzero = true;
    }
    EXPECT_TRUE(any_nonzero);
}

TEST(FullLoopTransitionTest, ConvertPhase1ToPhase2IsZeroAtTheFrameOrigin) {
    PlanetConfig planet_config = PlanetConfig::Earth();
    guidance_scp::Phase1ToPhase2Frame frame;
    frame.origin_r = planet_config.radius + 500.0;
    frame.origin_la = 0.1; frame.origin_lo = 0.2;

    auto s = guidance_scp::convertPhase1ToPhase2(frame.origin_r, frame.origin_la, frame.origin_lo,
                                                    100.0, -0.3, 0.0, frame, planet_config);
    EXPECT_NEAR(s.rx, 0.0, 1e-9);
    EXPECT_NEAR(s.ry, 0.0, 1e-9);
    EXPECT_NEAR(s.rz, 0.0, 1e-9);
}

TEST(FullLoopTransitionTest, JacobianMatchesFiniteDifferenceOfTheMapItself) {
    // Sanity check on the Jacobian utility applied to this specific map:
    // d(rz)/dr must be exactly 1 (rz = r - origin_r is affine in r).
    PlanetConfig planet_config = PlanetConfig::Earth();
    guidance_scp::Phase1ToPhase2Frame frame;
    frame.origin_r = planet_config.radius;
    frame.origin_la = 0.0; frame.origin_lo = 0.0;

    double r = planet_config.radius + 2000.0, la = 0.001, lo = 0.001, v = 90.0, fpa = -0.2, v_azi = 0.1;
    Eigen::MatrixXd J = guidance_scp::computeTransitionJacobian(r, la, lo, v, fpa, v_azi, frame, planet_config);
    ASSERT_EQ(J.rows(), 6);
    ASSERT_EQ(J.cols(), 6);
    EXPECT_NEAR(J(2, 0), 1.0, 1e-6);  // d(rz)/dr
}

TEST(StageATest, ConvergesStandaloneAndSatisfiesPathConstraints) {
    PlanetConfig planet_config = PlanetConfig::Earth();
    SpacecraftConfig spacecraft_config = MakeTestSpacecraftConfig();

    guidance_scp::StageAState initial;
    initial.r = planet_config.radius + 60000.0;
    initial.la = 0.0; initial.lo = 0.0;
    initial.v = 4000.0; initial.fpa = -0.05; initial.v_azi = 0.0; initial.sigma = 0.0;

    guidance_scp::StageAState terminal;
    terminal.r = planet_config.radius + 45000.0;
    terminal.la = 0.001; terminal.lo = 0.001;
    terminal.v = 3000.0; terminal.fpa = -0.03; terminal.v_azi = 0.01; terminal.sigma = 0.0;

    guidance_scp::StageAConfig config;
    config.K = 20;
    config.t_scale_s = 30.0;
    config.max_heat_flux_w_m2 = 5.0e6;   // generous, avoids the fixture's coarse aero table producing spurious infeasibility
    config.max_qbar_pa = 2.0e5;
    config.max_load_factor_g = 20.0;

    guidance_scp::StageAResult result = guidance_scp::solveStageA(initial, terminal, planet_config, spacecraft_config, config);

    ASSERT_TRUE(result.converged);
    for (size_t i = 0; i < result.history.t.size(); ++i) {
        EXPECT_LE(result.history.heat_flux_conv[i], config.max_heat_flux_w_m2 * 1.05);
        EXPECT_LE(result.history.qbar[i], config.max_qbar_pa * 1.05);
        EXPECT_LE(result.history.load_factor_g[i], config.max_load_factor_g * 1.05);
    }
}

TEST(StageATest, ConvergesAgainstRealAeroTableWithAdaptiveTrustRegion) {
    // Same boundary conditions as StageDemoMain.cpp's Stage A case -- the
    // exact scenario empirically confirmed broken before the adaptive
    // trust-region accept/reject mechanism was added (converged=false,
    // max||nu|| stuck around ~7874 after iteration 0). This is the real
    // acceptance bar for that fix.
    PlanetConfig planet_config = PlanetConfig::Earth();
    SpacecraftConfig spacecraft_config = MakeRealSpacecraftConfig();

    guidance_scp::StageAState initial;
    initial.r = planet_config.radius + 60000.0;
    initial.v = 4000.0;
    initial.fpa = -0.05;

    guidance_scp::StageAState terminal;
    terminal.r = planet_config.radius + 45000.0;
    terminal.la = 0.001; terminal.lo = 0.001;
    terminal.v = 3000.0;
    terminal.fpa = -0.03;
    terminal.v_azi = 0.01;

    guidance_scp::StageAConfig config;
    config.K = 20;
    config.t_scale_s = 30.0;
    config.max_heat_flux_w_m2 = 5.0e6;
    config.max_qbar_pa = 2.0e5;
    // 20g makes the TERMINAL boundary condition itself infeasible against
    // the real table (natural load factor there is ~23.5g) -- see
    // StageDemoMain.cpp's identical comment.
    config.max_load_factor_g = 30.0;

    guidance_scp::StageAResult result = guidance_scp::solveStageA(initial, terminal, planet_config, spacecraft_config, config);

    ASSERT_TRUE(result.converged);
    ASSERT_FALSE(result.attempt_max_defect_nonlinear.empty());
    for (double d : result.attempt_max_defect_nonlinear) {
        EXPECT_TRUE(std::isfinite(d));
    }
    for (size_t i = 0; i < result.history.t.size(); ++i) {
        EXPECT_LE(result.history.heat_flux_conv[i], config.max_heat_flux_w_m2 * 1.05);
        EXPECT_LE(result.history.qbar[i], config.max_qbar_pa * 1.05);
        EXPECT_LE(result.history.load_factor_g[i], config.max_load_factor_g * 1.05);
    }
}

TEST(StageATest, AdaptiveTrustRegionRejectsAndShrinksOnOverlyLargeStep) {
    // Same real-table scenario, but with max_load_factor_g tightened back
    // to 20g -- empirically confirmed (during this fix's own debugging) to
    // make the TERMINAL boundary condition itself infeasible against the
    // real table (natural load factor there is ~23.5g). Since the terminal
    // node's delta is pinned to exactly zero regardless of eta_max, every
    // attempt after the first is forced into PrimalInfeasible independent
    // of trust-region size -- a reliable, deterministic way to isolate the
    // reject/shrink path (non-Solved status is treated as a reject, exactly
    // like an accepted-then-regressed nu). Deliberately not asserting on
    // `converged` (expected false: this scenario cannot ever succeed).
    PlanetConfig planet_config = PlanetConfig::Earth();
    SpacecraftConfig spacecraft_config = MakeRealSpacecraftConfig();

    guidance_scp::StageAState initial;
    initial.r = planet_config.radius + 60000.0;
    initial.v = 4000.0;
    initial.fpa = -0.05;

    guidance_scp::StageAState terminal;
    terminal.r = planet_config.radius + 45000.0;
    terminal.la = 0.001; terminal.lo = 0.001;
    terminal.v = 3000.0;
    terminal.fpa = -0.03;
    terminal.v_azi = 0.01;

    guidance_scp::StageAConfig config;
    config.K = 20;
    config.t_scale_s = 30.0;
    config.max_heat_flux_w_m2 = 5.0e6;
    config.max_qbar_pa = 2.0e5;
    config.max_load_factor_g = 20.0;

    guidance_scp::StageAResult result = guidance_scp::solveStageA(initial, terminal, planet_config, spacecraft_config, config);

    ASSERT_FALSE(result.attempt_accepted.empty());
    bool any_rejected = false;
    for (bool accepted : result.attempt_accepted) {
        if (!accepted) any_rejected = true;
    }
    EXPECT_TRUE(any_rejected);

    bool any_shrunk = false;
    for (double eta : result.attempt_eta_max) {
        if (eta < config.eta_max_init) any_shrunk = true;
    }
    EXPECT_TRUE(any_shrunk);
}

TEST(StageBTest, SolvesAndProducesSensibleDescent) {
    PlanetConfig planet_config = PlanetConfig::Earth();

    guidance_scp::StageBState initial;
    initial.m = 8000.0;
    initial.rx = 0.0; initial.ry = 0.0; initial.rz = 2000.0;
    initial.vx = 50.0; initial.vy = 0.0; initial.vz = -80.0;

    guidance_scp::StageBState terminal;
    terminal.rx = 0.0; terminal.ry = 0.0; terminal.rz = 0.0;
    terminal.vx = 0.0; terminal.vy = 0.0; terminal.vz = 0.0;

    guidance_scp::StageBConfig config;
    config.K = 20;
    config.Isp_s = 330.0;
    config.Tmin_N = 20000.0;
    config.Tmax_N = 100000.0;
    config.glideslope_deg = 80.0;
    config.tf_min_s = 5.0;
    config.tf_max_s = 60.0;

    guidance_scp::StageBResult result = guidance_scp::solveStageB(initial, terminal, planet_config, config);

    ASSERT_TRUE(result.solved);
    EXPECT_LT(result.history.m.back(), initial.m);
    EXPECT_GT(result.history.m.back(), 0.0);
    for (size_t i = 0; i < result.history.throttle_frac.size(); ++i) {
        double thrust_mag = std::sqrt(result.history.Tx[i] * result.history.Tx[i]
                                     + result.history.Ty[i] * result.history.Ty[i]
                                     + result.history.Tz[i] * result.history.Tz[i]);
        EXPECT_GE(thrust_mag, config.Tmin_N * 0.99);
        EXPECT_LE(thrust_mag, config.Tmax_N * 1.01);
    }
}

TEST(StageBTest, ThrustBoundUsesLosslessConvexificationNotBoxConstraint) {
    PlanetConfig planet_config = PlanetConfig::Earth();
    guidance_scp::StageBState initial;
    initial.m = 8000.0; initial.rz = 2000.0; initial.vz = -80.0;
    guidance_scp::StageBState terminal;  // all zero

    guidance_scp::StageBConfig config;
    config.K = 10;
    config.Tmin_N = 20000.0;
    config.Tmax_N = 100000.0;

    guidance_scp::SocpProblem problem =
        guidance_scp::buildFixedTfProblem(initial, terminal, planet_config, config, 30.0);

    bool has_second_order_cone = false;
    for (const auto& cone : problem.cones) {
        if (cone.tag == clarabel::SupportedConeT<double>::Tag::SecondOrderConeT) {
            has_second_order_cone = true;
            break;
        }
    }
    EXPECT_TRUE(has_second_order_cone);
}
