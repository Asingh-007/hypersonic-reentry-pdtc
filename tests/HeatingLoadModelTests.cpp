#include <gtest/gtest.h>
#include "HeatingLoadModel.h"
#include "GuidanceConstraints.h"
#include "DescentDynamics.h"
#include "PlanetConfig.h"
#include "SpacecraftConfig.h"
#include "ControlInputs.h"
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

DescentState MakeNominalState(float q1, float q2, float q3, float q4,
                               float wx, float wy, float wz,
                               float altitude_m, float speed) {
    return DescentState(
        PlanetConfig::Earth().radius + altitude_m, 0.0f, 0.0f, speed, -0.05f, 0.0f,
        q1, q2, q3, q4, wx, wy, wz);
}

DescentDynamics::StateVector ToStateVector(const DescentState& s) {
    DescentDynamics::StateVector x;
    x << s.r, s.la, s.lo, s.v, s.fpa, s.v_azi,
         s.q1, s.q2, s.q3, s.q4, s.wx, s.wy, s.wz;
    return x;
}

}  // namespace

TEST(HeatingLoadModelTest, SuttonGravesScalesWithVelocityCubed) {
    double rho = 1e-3, nose_radius_m = 1.0;
    double q1 = suttonGravesHeatFlux(rho, 1000.0, nose_radius_m);
    double q2 = suttonGravesHeatFlux(rho, 2000.0, nose_radius_m);
    EXPECT_NEAR(q2 / q1, 8.0, 8.0 * 0.01);
}

TEST(HeatingLoadModelTest, SuttonGravesScalesInverseSqrtNoseRadius) {
    double rho = 1e-3, v = 5000.0;
    double q_small = suttonGravesHeatFlux(rho, v, 1.0);
    double q_large = suttonGravesHeatFlux(rho, v, 4.0);
    // q ~ 1/sqrt(R_n): quadrupling R_n halves heat flux.
    EXPECT_NEAR(q_small / q_large, 2.0, 2.0 * 0.01);
}

TEST(HeatingLoadModelTest, SuttonGravesMatchesKnownMagnitudeAtApolloLikeConditions) {
    // Apollo-like lunar-return entry conditions: rho ~ 4e-4 kg/m^3,
    // V ~ 11000 m/s, nose radius ~ 1m -- expect the well-known ~1-5 MW/m^2
    // peak-heating order of magnitude (loose, factor-of-2-ish tolerance;
    // this is a sanity/scaling cross-check, not an exact-match test).
    double q = suttonGravesHeatFlux(4e-4, 11000.0, 1.0);
    EXPECT_GT(q, 1.0e6);
    EXPECT_LT(q, 1.0e7);
}

// Tauber-Sutton (1991) reference values below were computed independently
// in Python from the exact formula/table transcribed in HeatingLoadModel.h's
// doc comment (Eqs. 1-2, Table 1, air/Earth), not derived from this C++
// implementation -- a real cross-check, mirroring how ShockExpansionAeroTest
// validates against known gas-dynamics-table values elsewhere in this repo.

TEST(HeatingLoadModelTest, TauberSuttonMatchesExactGridPointValue) {
    // V=11000 m/s is an exact Table 1 grid point (f=151), rho=1e-4 kg/m^3,
    // R_n=1m -- independently computed expected value: ~942732.85 W/m^2.
    double q = tauberSuttonRadiativeHeatFlux(1e-4, 11000.0, 1.0);
    EXPECT_NEAR(q, 942732.85, 942732.85 * 1e-3);
}

TEST(HeatingLoadModelTest, TauberSuttonLinearlyInterpolatesBetweenGridPoints) {
    // V=9125 m/s is exactly halfway between Table 1's 9000 (f=1.5) and 9250
    // (f=4.3) grid points -- independently computed expected value using the
    // linearly-interpolated f=2.9: ~18105.47 W/m^2.
    double q = tauberSuttonRadiativeHeatFlux(1e-4, 9125.0, 1.0);
    EXPECT_NEAR(q, 18105.47, 18105.47 * 1e-3);
}

TEST(HeatingLoadModelTest, TauberSuttonIsZeroBelowTableFloor) {
    // Below Table 1's 9000 m/s floor, radiative heating is genuinely
    // negligible for air (the paper's own motivation) -- exact zero, not an
    // extrapolation.
    double q = tauberSuttonRadiativeHeatFlux(1e-4, 8000.0, 1.0);
    EXPECT_NEAR(q, 0.0, 1e-9);
}

TEST(HeatingLoadModelTest, TauberSuttonClampsAboveTableCeiling) {
    // V=17000 m/s exceeds Table 1's 16000 m/s ceiling -- clamped at the last
    // tabulated f=2040 rather than extrapolating further. Independently
    // computed expected value using that clamp: ~12736258.38 W/m^2.
    double q = tauberSuttonRadiativeHeatFlux(1e-4, 17000.0, 1.0);
    EXPECT_NEAR(q, 12736258.38, 12736258.38 * 1e-3);
}

TEST(HeatingLoadModelTest, LoadFactorIsZeroAboveAtmosphereWithNoThrust) {
    PlanetConfig planet_config = PlanetConfig::Earth();
    SpacecraftConfig spacecraft_config = MakeTestSpacecraftConfig();
    ThrustVectorControlInputs control_inputs(0.0f, 0.0f, 0.0f);

    DescentState state = MakeNominalState(0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f,
                                            /*altitude_m=*/2000000.0f, /*speed=*/1000.0f);
    HeatingLoadResult hl = computeHeatingAndLoad(ToStateVector(state), planet_config, spacecraft_config, control_inputs);

    EXPECT_NEAR(hl.load_factor_g, 0.0, 1e-9);
}

TEST(HeatingLoadModelTest, LoadFactorIsNonzeroInDenseAtmosphere) {
    PlanetConfig planet_config = PlanetConfig::Earth();
    SpacecraftConfig spacecraft_config = MakeTestSpacecraftConfig();
    ThrustVectorControlInputs control_inputs(0.0f, 0.0f, 0.0f);

    DescentState state = MakeNominalState(0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f,
                                            /*altitude_m=*/30000.0f, /*speed=*/1000.0f);
    HeatingLoadResult hl = computeHeatingAndLoad(ToStateVector(state), planet_config, spacecraft_config, control_inputs);

    EXPECT_TRUE(std::isfinite(hl.load_factor_g));
    EXPECT_GT(hl.load_factor_g, 1e-6);
}

TEST(HeatingLoadModelTest, ComputeHeatingAndLoadTotalIsConvectivePlusRadiative) {
    PlanetConfig planet_config = PlanetConfig::Earth();
    SpacecraftConfig spacecraft_config = MakeTestSpacecraftConfig();
    ThrustVectorControlInputs control_inputs(0.0f, 0.0f, 0.0f);

    // 30km/1000 m/s is well below Tauber-Sutton's 9000 m/s floor, so
    // heat_flux_rad should be exactly 0 here and total == conv.
    DescentState state = MakeNominalState(0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f,
                                            /*altitude_m=*/30000.0f, /*speed=*/1000.0f);
    HeatingLoadResult hl = computeHeatingAndLoad(ToStateVector(state), planet_config, spacecraft_config, control_inputs);

    EXPECT_NEAR(hl.heat_flux_rad_w_m2, 0.0, 1e-9);
    EXPECT_NEAR(hl.heat_flux_total_w_m2, hl.heat_flux_conv_w_m2 + hl.heat_flux_rad_w_m2, 1e-6);
}

TEST(GuidanceConstraintsTest, CheckConstraintsFlagsClearlyViolatingSyntheticHistory) {
    DescentDynamics::TrajectoryHistory hist;
    hist.t = {0.0, 1.0, 2.0};
    hist.heat_flux_total = {0.0, 2.0e6, 0.0};  // exceeds default max_heat_flux_w_m2 = 1e6
    hist.qbar = {0.0, 60000.0, 0.0};           // exceeds default max_qbar_pa = 50000
    hist.load_factor = {0.0, 8.0, 0.0};        // exceeds default max_load_factor_g = 5.0

    GuidanceConstraints constraints;
    ConstraintViolationReport report = checkConstraints(hist, constraints);

    ASSERT_TRUE(report.hasViolations());
    ASSERT_EQ(report.violations.size(), 3u);
    bool has_heat_flux = false, has_qbar = false, has_load_factor = false;
    for (const auto& v : report.violations) {
        EXPECT_EQ(v.index, 1u);
        EXPECT_NEAR(v.time_s, 1.0, 1e-9);
        if (v.constraint_name == "heat_flux") { has_heat_flux = true; EXPECT_NEAR(v.margin, 1.0e6, 1e-3); }
        if (v.constraint_name == "qbar") { has_qbar = true; EXPECT_NEAR(v.margin, 10000.0, 1e-3); }
        if (v.constraint_name == "load_factor") { has_load_factor = true; EXPECT_NEAR(v.margin, 3.0, 1e-9); }
    }
    EXPECT_TRUE(has_heat_flux);
    EXPECT_TRUE(has_qbar);
    EXPECT_TRUE(has_load_factor);
}

TEST(GuidanceConstraintsTest, CheckConstraintsReportsNoViolationsForCompliantSyntheticHistory) {
    DescentDynamics::TrajectoryHistory hist;
    hist.t = {0.0, 1.0, 2.0};
    hist.heat_flux_total = {0.0, 1.0e5, 0.0};
    hist.qbar = {0.0, 1000.0, 0.0};
    hist.load_factor = {0.0, 1.0, 0.0};

    GuidanceConstraints constraints;
    ConstraintViolationReport report = checkConstraints(hist, constraints);

    EXPECT_FALSE(report.hasViolations());
}

TEST(HeatingLoadModelTest, NoNaNOverNominalEntryInterfaceRunIncludingHeatingFields) {
    PlanetConfig planet_config = PlanetConfig::Earth();
    SpacecraftConfig spacecraft_config = MakeTestSpacecraftConfig();
    DescentDynamics dynamics(planet_config, spacecraft_config);
    ThrustVectorControlInputs control_inputs(0.0f, 0.1f, 0.1f);

    DescentState initial(
        planet_config.radius + 120000.0f, 0.0f, 0.0f, 7500.0f, -0.1f, 0.0f,
        0.0f, 0.0f, 0.0f, 1.0f, 0.01f, 0.01f, 0.0f);
    auto history = dynamics.integrate(initial, control_inputs, 0.0, 20.0, 1.0, 1e-6);

    ASSERT_GT(history.t.size(), 0u);
    double prev_heat_load = -1.0;
    for (size_t i = 0; i < history.t.size(); ++i) {
        EXPECT_TRUE(std::isfinite(history.heat_flux_conv[i])) << "at t=" << history.t[i];
        EXPECT_TRUE(std::isfinite(history.heat_flux_rad[i])) << "at t=" << history.t[i];
        EXPECT_TRUE(std::isfinite(history.heat_flux_total[i])) << "at t=" << history.t[i];
        EXPECT_TRUE(std::isfinite(history.heat_load[i])) << "at t=" << history.t[i];
        EXPECT_TRUE(std::isfinite(history.load_factor[i])) << "at t=" << history.t[i];
        // 7500 m/s entry never reaches Tauber-Sutton's 9000 m/s floor.
        EXPECT_NEAR(history.heat_flux_rad[i], 0.0, 1e-9) << "at t=" << history.t[i];
        EXPECT_GE(history.heat_load[i], prev_heat_load) << "at t=" << history.t[i];
        prev_heat_load = history.heat_load[i];
    }
}
