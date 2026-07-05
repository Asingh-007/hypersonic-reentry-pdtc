#include <gtest/gtest.h>
#include "DescentDynamics.h"
#include "PlanetConfig.h"
#include "SpacecraftConfig.h"
#include "ControlInputs.h"
#include "QuaternionUtils.h"
#include <cmath>

TEST(QuaternionUtilsTest, KnownValueIdentityQuaternionUnitXRate) {
    Eigen::Quaterniond identity(1.0, 0.0, 0.0, 0.0); // w,x,y,z
    Eigen::Quaterniond qdot = QuaternionDerivative(identity, 1.0, 0.0, 0.0);
    EXPECT_NEAR(qdot.w(), 0.0, 1e-12);
    EXPECT_NEAR(qdot.x(), 0.5, 1e-12);
    EXPECT_NEAR(qdot.y(), 0.0, 1e-12);
    EXPECT_NEAR(qdot.z(), 0.0, 1e-12);
}

namespace {
DescentState MakeNominalState(float q1, float q2, float q3, float q4,
                               float wx, float wy, float wz) {
    return DescentState(
        PlanetConfig::Earth().radius + 100000.0f, 0.0f, 0.0f, 1000.0f, -0.05f, 0.0f,
        q1, q2, q3, q4, wx, wy, wz);
}
}

TEST(DescentDynamicsTest, RoundTripThroughZeroDurationIntegrate) {
    // toVector/fromVector are private, so exercise the round-trip indirectly:
    // a zero-duration integrate() call (t0==tf) records exactly one point --
    // the initial condition, converted to a StateVector and back -- with no
    // integration step taken.
    PlanetConfig planet_config = PlanetConfig::Earth();
    SpacecraftConfig spacecraft_config(1000.0f, 10.0f, 0.5f, 0.3f, 1000.0f, 1000.0f, 1500.0f);
    DescentDynamics dynamics(planet_config, spacecraft_config);
    ThrustVectorControlInputs control_inputs(0.0f, 0.0f, 0.0f);

    DescentState initial = MakeNominalState(0.1f, 0.2f, 0.3f, 0.9274f, 0.01f, 0.02f, 0.03f);

    auto history = dynamics.integrate(initial, control_inputs, 0.0, 0.0, 1.0, 1e-6);

    ASSERT_EQ(history.t.size(), 1u);
    EXPECT_NEAR(history.r[0], initial.r, 1e-2);
    EXPECT_NEAR(history.v[0], initial.v, 1e-4);
    EXPECT_NEAR(history.fpa[0], initial.fpa, 1e-6);
    EXPECT_NEAR(history.wx[0], initial.wx, 1e-6);
    EXPECT_NEAR(history.wy[0], initial.wy, 1e-6);
    EXPECT_NEAR(history.wz[0], initial.wz, 1e-6);
}

TEST(DescentDynamicsTest, TorqueFreeSymmetricBodyWithZeroSpinStaysAtRest) {
    PlanetConfig planet_config = PlanetConfig::Earth();
    SpacecraftConfig spacecraft_config(1000.0f, 10.0f, 0.5f, 0.3f, /*ixx=*/1000.0f, /*iyy=*/1000.0f, /*izz=*/1500.0f);
    DescentDynamics dynamics(planet_config, spacecraft_config);
    ThrustVectorControlInputs control_inputs(0.0f, 0.0f, 0.0f);

    DescentState state = MakeNominalState(0.0f, 0.0f, 0.0f, 1.0f, 0.01f, 0.01f, 0.0f);
    DescentDynamics::StateVector x;
    x << state.r, state.la, state.lo, state.v, state.fpa, state.v_azi,
         state.q1, state.q2, state.q3, state.q4, state.wx, state.wy, state.wz;

    DescentDynamics::StateVector xdot = dynamics.derivatives(x, control_inputs);

    EXPECT_NEAR(xdot(10), 0.0, 1e-9); // wx_dot
    EXPECT_NEAR(xdot(11), 0.0, 1e-9); // wy_dot
    EXPECT_NEAR(xdot(12), 0.0, 1e-9); // wz_dot
}

TEST(DescentDynamicsTest, AngularMomentumMagnitudeConservedUnderTorqueFreeMotion) {
    // Invariant of torque-free rigid-body motion: |J*w| stays
    // constant. Uses asymmetric inertia and nonzero rates on all three axes
    // so the check actually exercises coupled rotational dynamics, not the
    // degenerate zero-motion case above.
    PlanetConfig planet_config = PlanetConfig::Earth();
    SpacecraftConfig spacecraft_config(1000.0f, 10.0f, 0.5f, 0.3f, /*ixx=*/800.0f, /*iyy=*/1000.0f, /*izz=*/1500.0f);
    DescentDynamics dynamics(planet_config, spacecraft_config);
    ThrustVectorControlInputs control_inputs(0.0f, 0.0f, 0.0f);

    DescentState initial = MakeNominalState(0.0f, 0.0f, 0.0f, 1.0f, 0.05f, 0.03f, 0.02f);
    auto history = dynamics.integrate(initial, control_inputs, 0.0, 5.0, 0.1, 1e-8);

    ASSERT_GT(history.t.size(), 1u);
    double ixx = 800.0, iyy = 1000.0, izz = 1500.0;
    double L0 = std::sqrt(std::pow(ixx * history.wx[0], 2) +
                          std::pow(iyy * history.wy[0], 2) +
                          std::pow(izz * history.wz[0], 2));
    for (size_t i = 0; i < history.t.size(); ++i) {
        double L = std::sqrt(std::pow(ixx * history.wx[i], 2) +
                             std::pow(iyy * history.wy[i], 2) +
                             std::pow(izz * history.wz[i], 2));
        EXPECT_NEAR(L, L0, L0 * 1e-3) << "at t=" << history.t[i];
    }
}

TEST(DescentDynamicsTest, QuaternionNormStaysNearOneOverShortRun) {
    PlanetConfig planet_config = PlanetConfig::Earth();
    SpacecraftConfig spacecraft_config(1000.0f, 10.0f, 0.5f, 0.3f, 1000.0f, 1000.0f, 1500.0f);
    DescentDynamics dynamics(planet_config, spacecraft_config);
    ThrustVectorControlInputs control_inputs(0.0f, 0.0f, 0.0f);

    DescentState initial = MakeNominalState(0.05f, 0.05f, 0.05f, std::sqrt(1.0f - 3.0f * 0.05f * 0.05f), 0.02f, 0.01f, 0.03f);
    auto history = dynamics.integrate(initial, control_inputs, 0.0, 10.0, 0.5, 1e-8);

    ASSERT_GT(history.t.size(), 0u);
    for (size_t i = 0; i < history.t.size(); ++i) {
        double norm = std::sqrt(history.q1[i] * history.q1[i] + history.q2[i] * history.q2[i] +
                                history.q3[i] * history.q3[i] + history.q4[i] * history.q4[i]);
        EXPECT_NEAR(norm, 1.0, 1e-6) << "at t=" << history.t[i];
    }
}

TEST(DescentDynamicsTest, NoNaNOverNominalEntryInterfaceRun) {
    PlanetConfig planet_config = PlanetConfig::Earth();
    SpacecraftConfig spacecraft_config(1000.0f, 10.0f, 0.5f, 0.3f, 1000.0f, 1000.0f, 1500.0f);
    DescentDynamics dynamics(planet_config, spacecraft_config);
    ThrustVectorControlInputs control_inputs(0.0f, 0.1f, 0.1f);

    DescentState initial(
        planet_config.radius + 120000.0f, 0.0f, 0.0f, 7500.0f, -0.1f, 0.0f,
        0.0f, 0.0f, 0.0f, 1.0f, 0.01f, 0.01f, 0.0f);
    auto history = dynamics.integrate(initial, control_inputs, 0.0, 20.0, 1.0, 1e-6);

    ASSERT_GT(history.t.size(), 0u);
    for (size_t i = 0; i < history.t.size(); ++i) {
        EXPECT_TRUE(std::isfinite(history.r[i])) << "at t=" << history.t[i];
        EXPECT_TRUE(std::isfinite(history.v[i])) << "at t=" << history.t[i];
        EXPECT_TRUE(std::isfinite(history.fpa[i])) << "at t=" << history.t[i];
        EXPECT_TRUE(std::isfinite(history.wx[i])) << "at t=" << history.t[i];
        EXPECT_TRUE(std::isfinite(history.qbar[i])) << "at t=" << history.t[i];
    }
}
