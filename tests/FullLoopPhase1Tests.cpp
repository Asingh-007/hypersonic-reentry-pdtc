#include <gtest/gtest.h>
#include "FiniteDifferenceJacobian.h"
#include "full_loop_phase1.h"
#include "PlanetConfig.h"
#include "SpacecraftConfig.h"
#include "FlapActuatorConfig.h"
#include "DescentDynamics.h"
#include "AeroAngles.h"
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

// A nominal, fully-populated 21-state Phase 1 vector (identity attitude,
// zero rates/flaps) -- not aerodynamically trimmed, but well-defined for
// exercising the EOM and its finite-difference Jacobian.
Eigen::VectorXd MakeNominalPhase1State(double r, double v, double fpa) {
    Eigen::VectorXd x(guidance_scp::kPhase1StateDim);
    x << r, 0.0, 0.0, v, fpa, 0.0,
         0.0, 0.0, 0.0, 1.0,   // identity quaternion
         0.0, 0.0, 0.0,        // body rates
         0.0, 0.0, 0.0, 0.0,   // flap deflections
         0.0, 0.0, 0.0, 0.0;   // flap rates
    return x;
}

}  // namespace

TEST(FiniteDifferenceJacobianTest, RecoversKnownLinearToyFunction) {
    // f(x,u) = A*x + B*u for known, arbitrary A (3x3), B (3x2) -- the FD
    // Jacobian of a linear function should recover A, B to near machine
    // precision, confirming the generalized utility works before trusting
    // it on the much larger real EOMs.
    Eigen::Matrix3d A;
    A << 1.0, 2.0, -1.0,
         0.0, -3.0, 4.0,
         5.0, 0.5, -2.0;
    Eigen::Matrix<double, 3, 2> B;
    B << 2.0, -1.0,
         0.0, 3.0,
         -4.0, 1.0;

    guidance_scp::EomFn f = [&](const Eigen::VectorXd& x, const Eigen::VectorXd& u) -> Eigen::VectorXd {
        return A * x + B * u;
    };

    Eigen::VectorXd x_ref(3), u_ref(2);
    x_ref << 1.5, -2.0, 0.5;
    u_ref << 0.3, -0.7;

    guidance_scp::EomJacobian jac = guidance_scp::computeEomJacobianFd(f, x_ref, u_ref);

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            EXPECT_NEAR(jac.Ac(i, j), A(i, j), 1e-6);
        }
        for (int j = 0; j < 2; ++j) {
            EXPECT_NEAR(jac.Bc(i, j), B(i, j), 1e-6);
        }
    }
}

TEST(FullLoopPhase1Test, TranslationalPartialMatchesClosedForm) {
    // r_dot = V*sin(fpa) is elementary and independent of every other
    // state/control -- d(r_dot)/dV = sin(fpa) exactly, regardless of
    // attitude/flaps/aero. A direct transcription-error catcher for the
    // largest, most error-prone new EOM in this plan.
    PlanetConfig planet_config = PlanetConfig::Earth();
    SpacecraftConfig spacecraft_config = MakeTestSpacecraftConfig();
    FlapActuatorConfig flap_config;

    double r = planet_config.radius + 50000.0, v = 4000.0, fpa = -0.07;
    Eigen::VectorXd x_ref = MakeNominalPhase1State(r, v, fpa);
    Eigen::VectorXd u_ref = Eigen::VectorXd::Zero(guidance_scp::kPhase1ControlDim);

    guidance_scp::EomFn f = [&](const Eigen::VectorXd& x, const Eigen::VectorXd& u) {
        return guidance_scp::phase1Eom(x, u, planet_config, spacecraft_config, flap_config);
    };
    guidance_scp::EomJacobian jac = guidance_scp::computeEomJacobianFd(f, x_ref, u_ref);

    EXPECT_NEAR(jac.Ac(0, 3), std::sin(fpa), 1e-6);
}

TEST(FullLoopPhase1Test, FlapActuatorOdeMatchesFormulaDirectly) {
    // Direct transcription check: at an arbitrary state/control, phase1Eom's
    // flap-actuator rows must equal Jeff_inv*(H_i - b*ddot_i + N*tau_m_i)
    // computed independently here from the same aero-table lookup, and
    // d(d_i)/dt must equal the ddot_i state exactly (structural check).
    PlanetConfig planet_config = PlanetConfig::Earth();
    SpacecraftConfig spacecraft_config = MakeTestSpacecraftConfig();
    FlapActuatorConfig flap_config;

    double r = planet_config.radius + 50000.0, v = 4000.0, fpa = -0.05;
    Eigen::VectorXd x = MakeNominalPhase1State(r, v, fpa);
    x(17) = 0.02; x(18) = -0.01; x(19) = 0.03; x(20) = -0.02;  // nonzero flap rates
    Eigen::VectorXd u(guidance_scp::kPhase1ControlDim);
    u << 10.0, -5.0, 8.0, -3.0;

    Eigen::VectorXd xdot = guidance_scp::phase1Eom(x, u, planet_config, spacecraft_config, flap_config);

    for (int i = 0; i < 4; ++i) {
        EXPECT_NEAR(xdot(13 + i), x(17 + i), 1e-12);  // d(d_i)/dt == ddot_i
    }

    // Independently recompute qbar and the table's Ch at this state, using
    // the SAME AeroAngles extraction phase1Eom uses internally (identity
    // attitude does NOT mean alpha=beta=0 -- ComputeAeroAngles' extraction
    // is attitude/flight-path dependent, not just attitude-dependent).
    double a_sound = DescentDynamics::speedOfSound(r, planet_config);
    double mach = v / a_sound;
    double rho = DescentDynamics::atmosphereDensity(r, planet_config);
    double qbar = 0.5 * rho * v * v;
    Eigen::Quaterniond q(1.0, 0.0, 0.0, 0.0);  // identity, matches MakeNominalPhase1State
    AeroAngles angles = ComputeAeroAngles(q, fpa, /*v_azi=*/0.0);
    double alpha_deg = angles.alpha_rad * 180.0 / kPi;
    double beta_deg = angles.beta_rad * 180.0 / kPi;
    auto aero = spacecraft_config.aero_table.interpolate(mach, alpha_deg, beta_deg, 0.0, 0.0, 0.0);
    double ddot_i_state[4] = {0.02, -0.01, 0.03, -0.02};
    double tau_m[4] = {10.0, -5.0, 8.0, -3.0};
    for (int i = 0; i < 4; ++i) {
        double H_i = qbar * spacecraft_config.S_ref * spacecraft_config.L_ref * aero.Ch[i];
        double expected = flap_config.Jeff_inv *
            (H_i - flap_config.b_damping_n_m_s * ddot_i_state[i] + flap_config.N_gear_ratio * tau_m[i]);
        EXPECT_NEAR(xdot(17 + i), expected, std::abs(expected) * 1e-6 + 1e-9);
    }
}

TEST(FullLoopPhase1Test, BFlapMagnitudeScalesWithReferenceQbar) {
    // Section 3c acceptance-checklist requirement: B_flap = d(flap-actuator
    // acceleration)/d(flap deflection), extracted from the FD Jacobian, must
    // be proportional to the reference qbar. Uses two different altitudes
    // (same velocity) so qbar differs substantially; the fixture's Ch is
    // deliberately flap-axis-dependent (see AeroTests.cpp's fixture comment)
    // so this partial is genuinely nonzero.
    PlanetConfig planet_config = PlanetConfig::Earth();
    SpacecraftConfig spacecraft_config = MakeTestSpacecraftConfig();
    FlapActuatorConfig flap_config;
    double v = 4000.0, fpa = -0.05;

    double r_low = planet_config.radius + 20000.0;   // denser air, higher qbar
    double r_high = planet_config.radius + 60000.0;  // thinner air, lower qbar

    auto bflap_and_qbar = [&](double r) {
        Eigen::VectorXd x_ref = MakeNominalPhase1State(r, v, fpa);
        Eigen::VectorXd u_ref = Eigen::VectorXd::Zero(guidance_scp::kPhase1ControlDim);
        guidance_scp::EomFn f = [&](const Eigen::VectorXd& x, const Eigen::VectorXd& u) {
            return guidance_scp::phase1Eom(x, u, planet_config, spacecraft_config, flap_config);
        };
        guidance_scp::EomJacobian jac = guidance_scp::computeEomJacobianFd(f, x_ref, u_ref);
        double b_flap = jac.Ac(17, 13);  // d(ddot1_accel)/d(d1)

        double rho = DescentDynamics::atmosphereDensity(r, planet_config);
        double qbar = 0.5 * rho * v * v;
        return std::make_pair(b_flap, qbar);
    };

    auto [b_flap_low, qbar_low] = bflap_and_qbar(r_low);
    auto [b_flap_high, qbar_high] = bflap_and_qbar(r_high);

    ASSERT_GT(std::abs(b_flap_low), 1e-12);
    ASSERT_GT(std::abs(b_flap_high), 1e-12);

    double ratio_b = b_flap_low / b_flap_high;
    double ratio_qbar = qbar_low / qbar_high;
    EXPECT_NEAR(ratio_b, ratio_qbar, ratio_qbar * 0.05);
}
