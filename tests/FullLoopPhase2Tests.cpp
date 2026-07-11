#include <gtest/gtest.h>
#include "full_loop_phase2.h"
#include "FiniteDifferenceJacobian.h"
#include "PlanetConfig.h"
#include "SpacecraftConfig.h"
#include "FlapActuatorConfig.h"
#include "GimbalActuatorConfig.h"
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

Eigen::VectorXd MakeNominalPhase2State(double m, double rz, double vx, double vy, double vz) {
    Eigen::VectorXd x(guidance_scp::kPhase2StateDim);
    x << m, 0.0, 0.0, rz,
         vx, vy, vz,
         0.0, 0.0, 0.0, 1.0,  // identity quaternion
         0.0, 0.0, 0.0,       // body rates
         0.0, 0.0, 0.0, 0.0,  // flap deflections
         0.0, 0.0, 0.0, 0.0;  // flap rates
    return x;
}

}  // namespace

TEST(FullLoopPhase2Test, MassDepletionMatchesRocketEquationExactly) {
    PlanetConfig planet_config = PlanetConfig::Earth();
    SpacecraftConfig spacecraft_config = MakeTestSpacecraftConfig();
    FlapActuatorConfig flap_config;
    GimbalActuatorConfig gimbal_config;
    double Isp_s = 330.0;

    Eigen::VectorXd x = MakeNominalPhase2State(8000.0, 2000.0, 50.0, 0.0, -80.0);
    Eigen::VectorXd u(guidance_scp::kPhase2ControlDim);
    double T = 60000.0;
    u << T, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;

    Eigen::VectorXd xdot = guidance_scp::phase2Eom(x, u, /*deltaE_ref=*/0.0, /*phiE_ref=*/0.0, Isp_s,
                                                     planet_config, spacecraft_config, flap_config, gimbal_config);

    double expected_mdot = -T / (Isp_s * planet_config.g_0);
    EXPECT_NEAR(xdot(0), expected_mdot, std::abs(expected_mdot) * 1e-9);
}

TEST(FullLoopPhase2Test, TvcTorqueMatchesHandComputedCrossProduct) {
    // Zero velocity => qbar=0 => tau_aero=F_aero=0 exactly, isolating the
    // TVC-only torque so it can be checked against an independently
    // hand-computed cross product (moment arm x gimbaled thrust force).
    PlanetConfig planet_config = PlanetConfig::Earth();
    SpacecraftConfig spacecraft_config = MakeTestSpacecraftConfig();
    FlapActuatorConfig flap_config;
    GimbalActuatorConfig gimbal_config;  // engine_gimbal_point_body_m = (-18,0,0)
    double Isp_s = 330.0;

    Eigen::VectorXd x = MakeNominalPhase2State(8000.0, 2000.0, 0.0, 0.0, 0.0);
    Eigen::VectorXd u(guidance_scp::kPhase2ControlDim);
    double T = 50000.0;
    u << T, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;

    double deltaE = 0.1, phiE = 0.0;
    Eigen::VectorXd xdot = guidance_scp::phase2Eom(x, u, deltaE, phiE, Isp_s,
                                                     planet_config, spacecraft_config, flap_config, gimbal_config);

    // Independently recompute: thrust_dir_body mirrors freestreamDirectionBody's
    // structure (see full_loop_phase2.cpp), moment arm = gimbal point - moment_ref.
    Eigen::Vector3d thrust_dir_body(-std::cos(deltaE) * std::cos(phiE), std::sin(phiE),
                                     std::sin(deltaE) * std::cos(phiE));
    Eigen::Vector3d F_thrust_body = T * thrust_dir_body;
    Eigen::Vector3d arm = gimbal_config.engine_gimbal_point_body_m - spacecraft_config.moment_ref;
    Eigen::Vector3d expected_tau = arm.cross(F_thrust_body);

    const Eigen::Matrix3d& J = spacecraft_config.inertia;
    Eigen::Vector3d expected_wdot = J.inverse() * expected_tau;  // w=0, so no gyroscopic term

    EXPECT_NEAR(xdot(11), expected_wdot(0), std::abs(expected_wdot(0)) * 1e-6 + 1e-9);
    EXPECT_NEAR(xdot(12), expected_wdot(1), std::abs(expected_wdot(1)) * 1e-6 + 1e-9);
    EXPECT_NEAR(xdot(13), expected_wdot(2), std::abs(expected_wdot(2)) * 1e-6 + 1e-9);
}

TEST(FullLoopPhase2Test, FlapActuatorOdeStructuralCheck) {
    PlanetConfig planet_config = PlanetConfig::Earth();
    SpacecraftConfig spacecraft_config = MakeTestSpacecraftConfig();
    FlapActuatorConfig flap_config;
    GimbalActuatorConfig gimbal_config;

    Eigen::VectorXd x = MakeNominalPhase2State(8000.0, 2000.0, 20.0, 0.0, -50.0);
    x(18) = 0.01; x(19) = -0.02; x(20) = 0.015; x(21) = -0.01;
    Eigen::VectorXd u(guidance_scp::kPhase2ControlDim);
    u << 40000.0, 0.0, 0.0, 1.0, -1.0, 2.0, -2.0;

    Eigen::VectorXd xdot = guidance_scp::phase2Eom(x, u, 0.0, 0.0, 330.0,
                                                     planet_config, spacecraft_config, flap_config, gimbal_config);
    for (int i = 0; i < 4; ++i) {
        EXPECT_NEAR(xdot(14 + i), x(18 + i), 1e-12);  // d(d_i)/dt == ddot_i
    }
}
