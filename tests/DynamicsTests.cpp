#include <gtest/gtest.h>
#include "DescentDynamics.h"
#include "PlanetConfig.h"
#include "SpacecraftConfig.h"
#include "ControlInputs.h"
#include "QuaternionUtils.h"
#include <cmath>

// MSVC only defines M_PI when _USE_MATH_DEFINES is set before the first
// <cmath>/<math.h> include anywhere in the translation unit -- use an
// explicit local constant instead (same fix as main.cpp).
constexpr double kPi = 3.14159265358979323846;

TEST(QuaternionUtilsTest, KnownValueIdentityQuaternionUnitXRate) {
    Eigen::Quaterniond identity(1.0, 0.0, 0.0, 0.0); // w,x,y,z
    Eigen::Quaterniond qdot = QuaternionDerivative(identity, 1.0, 0.0, 0.0);
    EXPECT_NEAR(qdot.w(), 0.0, 1e-12);
    EXPECT_NEAR(qdot.x(), 0.5, 1e-12);
    EXPECT_NEAR(qdot.y(), 0.0, 1e-12);
    EXPECT_NEAR(qdot.z(), 0.0, 1e-12);
}

namespace {

std::string TestAeroTablePath() {
    return std::string(TESTS_SOURCE_DIR) + "/tests/fixtures/aero_table_test.csv";
}

SpacecraftConfig MakeTestSpacecraftConfig(const Eigen::Matrix3d& inertia) {
    double S_ref = kPi * 4.5 * 4.5;
    double L_ref = 9.0;
    Eigen::Vector3d moment_ref(20.0, 0.0, 0.0);
    return SpacecraftConfig(1000.0f, inertia, S_ref, L_ref, moment_ref, TestAeroTablePath());
}

// altitude_m defaults to 100 km (used by tests that don't care about the
// magnitude of aerodynamic forces/torques, just general dynamics behavior).
DescentState MakeNominalState(float q1, float q2, float q3, float q4,
                               float wx, float wy, float wz,
                               float altitude_m = 100000.0f) {
    return DescentState(
        PlanetConfig::Earth().radius + altitude_m, 0.0f, 0.0f, 1000.0f, -0.05f, 0.0f,
        q1, q2, q3, q4, wx, wy, wz);
}

}

TEST(DescentDynamicsTest, RoundTripThroughZeroDurationIntegrate) {
    // toVector/fromVector are private, so exercise the round-trip indirectly:
    // a zero-duration integrate() call (t0==tf) records exactly one point --
    // the initial condition, converted to a StateVector and back -- with no
    // integration step taken.
    PlanetConfig planet_config = PlanetConfig::Earth();
    SpacecraftConfig spacecraft_config = MakeTestSpacecraftConfig(Eigen::Vector3d(1000.0, 1000.0, 1500.0).asDiagonal());
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
    // Tested at a very high altitude (above the atmosphere model's upper
    // bound, where density -- and hence qbar/tau_aero -- is exactly zero
    // regardless of the aero table's contents) so this remains a pure test
    // of the torque-free rotational math itself, independent of aero data.
    PlanetConfig planet_config = PlanetConfig::Earth();
    SpacecraftConfig spacecraft_config = MakeTestSpacecraftConfig(Eigen::Vector3d(1000.0, 1000.0, 1500.0).asDiagonal()); // ixx=iyy
    DescentDynamics dynamics(planet_config, spacecraft_config);
    ThrustVectorControlInputs control_inputs(0.0f, 0.0f, 0.0f);

    DescentState state = MakeNominalState(0.0f, 0.0f, 0.0f, 1.0f, 0.01f, 0.01f, 0.0f, /*altitude_m=*/2000000.0f);
    DescentDynamics::StateVector x;
    x << state.r, state.la, state.lo, state.v, state.fpa, state.v_azi,
         state.q1, state.q2, state.q3, state.q4, state.wx, state.wy, state.wz;

    DescentDynamics::StateVector xdot = dynamics.derivatives(x, control_inputs);

    EXPECT_NEAR(xdot(10), 0.0, 1e-9); // wx_dot
    EXPECT_NEAR(xdot(11), 0.0, 1e-9); // wy_dot
    EXPECT_NEAR(xdot(12), 0.0, 1e-9); // wz_dot
}

TEST(DescentDynamicsTest, AngularMomentumMagnitudeConservedUnderTorqueFreeMotion) {
    // Invariant of torque-free rigid-body motion: |J*w| stays constant.
    // Tested at a very high altitude (same reasoning as the test above) so
    // tau_aero is exactly zero and this remains a pure test of the
    // rotational math, using asymmetric inertia and nonzero rates on all
    // three axes so the check actually exercises coupled dynamics.
    PlanetConfig planet_config = PlanetConfig::Earth();
    Eigen::Matrix3d inertia = Eigen::Vector3d(800.0, 1000.0, 1500.0).asDiagonal();
    SpacecraftConfig spacecraft_config = MakeTestSpacecraftConfig(inertia);
    DescentDynamics dynamics(planet_config, spacecraft_config);
    ThrustVectorControlInputs control_inputs(0.0f, 0.0f, 0.0f);

    DescentState initial = MakeNominalState(0.0f, 0.0f, 0.0f, 1.0f, 0.05f, 0.03f, 0.02f, /*altitude_m=*/2000000.0f);
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
    SpacecraftConfig spacecraft_config = MakeTestSpacecraftConfig(Eigen::Vector3d(1000.0, 1000.0, 1500.0).asDiagonal());
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
    SpacecraftConfig spacecraft_config = MakeTestSpacecraftConfig(Eigen::Vector3d(1000.0, 1000.0, 1500.0).asDiagonal());
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

TEST(DescentDynamicsTest, DerivativesProducesNonzeroTauAeroAndSensibleLiftDragAtNonzeroAlpha) {
    // At a low, dense-atmosphere altitude with the vehicle's nose tilted
    // relative to the velocity vector (nonzero alpha), derivatives() should
    // now produce nonzero aero-driven wdot (unlike the old torque-free
    // model) and nonzero, finite lift/drag-influenced translational rates.
    PlanetConfig planet_config = PlanetConfig::Earth();
    SpacecraftConfig spacecraft_config = MakeTestSpacecraftConfig(Eigen::Vector3d(1000.0, 1000.0, 1500.0).asDiagonal());
    DescentDynamics dynamics(planet_config, spacecraft_config);
    ThrustVectorControlInputs control_inputs(0.0f, 0.0f, 0.0f);

    // Identity attitude with fpa=-0.05 rad gives a nonzero angle of attack
    // (body nose not aligned with the velocity direction).
    DescentState state = MakeNominalState(0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, /*altitude_m=*/30000.0f);
    DescentDynamics::StateVector x;
    x << state.r, state.la, state.lo, state.v, state.fpa, state.v_azi,
         state.q1, state.q2, state.q3, state.q4, state.wx, state.wy, state.wz;

    DescentDynamics::StateVector xdot = dynamics.derivatives(x, control_inputs);

    // wdot should no longer be trivially zero now that tau_aero is real.
    bool any_nonzero_wdot = std::abs(xdot(10)) > 1e-12 || std::abs(xdot(11)) > 1e-12 || std::abs(xdot(12)) > 1e-12;
    EXPECT_TRUE(any_nonzero_wdot);
    EXPECT_TRUE(std::isfinite(xdot(3))); // v_dot finite (lift/drag applied)
    EXPECT_TRUE(std::isfinite(xdot(4))); // fpa_dot finite
}
