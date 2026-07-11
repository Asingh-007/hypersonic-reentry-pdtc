#include <gtest/gtest.h>
#include "full_loop_subproblem.h"
#include "full_loop_phase1.h"
#include "full_loop_phase2.h"
#include "PlanetConfig.h"
#include "SpacecraftConfig.h"
#include <set>

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

}  // namespace

TEST(FullLoopVarLayoutTest, OffsetsDoNotOverlapAndTotalMatchesSumOfBlocks) {
    guidance_scp::FullLoopVarLayout layout{5, 5};
    std::set<int> seen;
    auto mark = [&](int col) {
        ASSERT_TRUE(seen.insert(col).second) << "column " << col << " assigned twice";
    };

    for (int k = 0; k < layout.K1; ++k) {
        for (int i = 0; i < guidance_scp::FullLoopVarLayout::kNx1; ++i) mark(layout.dx1_offset(k) + i);
        mark(layout.eta1_offset(k));
    }
    for (int k = 0; k < layout.K1 - 1; ++k) {
        for (int i = 0; i < guidance_scp::FullLoopVarLayout::kNu1; ++i) mark(layout.du1_offset(k) + i);
        for (int i = 0; i < guidance_scp::FullLoopVarLayout::kNx1; ++i) {
            mark(layout.nu1_plus_offset(k) + i);
            mark(layout.nu1_minus_offset(k) + i);
        }
    }
    for (int k = 0; k < layout.K2; ++k) {
        for (int i = 0; i < guidance_scp::FullLoopVarLayout::kNx2; ++i) mark(layout.dx2_offset(k) + i);
        mark(layout.eta2_offset(k));
    }
    for (int k = 0; k < layout.K2 - 1; ++k) {
        for (int i = 0; i < guidance_scp::FullLoopVarLayout::kNu2Vec; ++i) mark(layout.du2_uvec_offset(k) + i);
        mark(layout.du2_sigma_offset(k));
        for (int i = 0; i < guidance_scp::FullLoopVarLayout::kNu2Rest; ++i) mark(layout.du2_rest_offset(k) + i);
        for (int i = 0; i < guidance_scp::FullLoopVarLayout::kNx2; ++i) {
            mark(layout.nu2_plus_offset(k) + i);
            mark(layout.nu2_minus_offset(k) + i);
        }
    }

    EXPECT_EQ(static_cast<int>(seen.size()), layout.total());
    EXPECT_EQ(*seen.begin(), 0);
    EXPECT_EQ(*seen.rbegin(), layout.total() - 1);
}

TEST(FullLoopSubproblemTest, BuildProducesDimensionallyConsistentProblem) {
    // Small, non-propagated (constant) reference -- this test only checks
    // the assembly produces a structurally well-formed SocpProblem (no
    // dimension mismatch/out-of-range column indices), not that it solves
    // to a physically converged answer (that is FullLoopIntegrationTests'
    // job, built on a properly-propagated reference from the outer loop).
    PlanetConfig planet_config = PlanetConfig::Earth();
    SpacecraftConfig spacecraft_config = MakeTestSpacecraftConfig();

    const int K1 = 4, K2 = 4;
    guidance_scp::FullLoopSubproblemInputs in;
    in.planet_config = planet_config;
    in.spacecraft_config = &spacecraft_config;

    Eigen::VectorXd x1_nom(guidance_scp::kPhase1StateDim);
    x1_nom << planet_config.radius + 50000.0, 0.0, 0.0, 4000.0, -0.05, 0.0,
              0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    Eigen::VectorXd u1_nom = Eigen::VectorXd::Zero(guidance_scp::kPhase1ControlDim);
    for (int k = 0; k < K1; ++k) in.x1_ref.push_back(x1_nom);
    for (int k = 0; k < K1 - 1; ++k) in.u1_ref.push_back(u1_nom);
    in.t_scale_1 = 30.0;
    in.x1_initial = x1_nom;

    Eigen::VectorXd x2_nom(guidance_scp::kPhase2StateDim);
    x2_nom << 8000.0, 0.0, 0.0, 2000.0, 20.0, 0.0, -50.0,
              0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    Eigen::VectorXd u2_nom(guidance_scp::kPhase2ControlDim);
    u2_nom << 60000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    for (int k = 0; k < K2; ++k) {
        in.x2_ref.push_back(x2_nom);
        in.deltaE_ref.push_back(0.0);
        in.phiE_ref.push_back(0.0);
    }
    for (int k = 0; k < K2 - 1; ++k) in.u2_ref.push_back(u2_nom);
    in.t_scale_2 = 20.0;

    in.x2_terminal_target = Eigen::VectorXd::Zero(guidance_scp::kPhase2StateDim);
    in.x2_terminal_target(10) = 1.0;  // identity quaternion (scalar part)

    in.transition_frame.origin_r = planet_config.radius;
    in.transition_frame.origin_la = 0.0;
    in.transition_frame.origin_lo = 0.0;

    in.full_loop_config.Isp_s = 330.0;
    in.full_loop_config.Tmin_N = 20000.0;
    in.full_loop_config.Tmax_N = 100000.0;
    in.full_loop_config.m_wet_at_handoff_kg = 8000.0;

    guidance_scp::FullLoopVarLayout layout{K1, K2};
    guidance_scp::SocpProblem problem = guidance_scp::buildFullLoopSubproblem(in);

    EXPECT_EQ(problem.A.cols(), layout.total());
    EXPECT_EQ(problem.A.rows(), static_cast<int>(problem.b.size()));
    EXPECT_EQ(problem.q.size(), layout.total());
    ASSERT_FALSE(problem.cones.empty());

    // Solving is not asserted here (see comment above) -- just confirm the
    // solve CALL itself doesn't throw/crash on this well-formed problem.
    EXPECT_NO_THROW(guidance_scp::solve(problem));
}
