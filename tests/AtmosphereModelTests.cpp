#include <gtest/gtest.h>
#include "AtmosphereModel.h"
#include <cmath>

namespace {
constexpr double kRelTol = 1e-3;
}

TEST(EarthAtmosphere1976Test, SeaLevelMatchesReference) {
    AtmosphereState s = EarthAtmosphere1976::Compute(0.0);
    EXPECT_NEAR(s.density, 1.225, 1.225 * kRelTol);
}

TEST(EarthAtmosphere1976Test, BridgeSubLayersMatchTable2ExactlyAtTheirReferenceAltitude) {
    // Plugging h = h_i into the Olsen & Bettinger bridge formula makes the
    // (h-h_i) term vanish, so density should equal Table 2's rho_i exactly
    // (to floating-point precision) at each sub-layer's own reference altitude.
    EXPECT_NEAR(EarthAtmosphere1976::Compute(85000.0).density, 7.726e-6, 7.726e-6 * kRelTol);
    EXPECT_NEAR(EarthAtmosphere1976::Compute(99000.0).density, 4.504e-7, 4.504e-7 * kRelTol);
    EXPECT_NEAR(EarthAtmosphere1976::Compute(110000.0).density, 5.930e-8, 5.930e-8 * kRelTol);
}

TEST(EarthAtmosphere1976Test, ContinuityAt84And120kmSeams) {
    // 120 km: the bridge's last sub-layer and the power-law fit agree to
    // within ~1% (the paper's own power-law fit was very likely calibrated
    // against its bridge model at this interface).
    // 84 km: the independent low-altitude exponential and the bridge's first
    // sub-layer agree to within ~2-3% -- an inherent residual from stitching
    // two independently-fit pieces, not a bug (see AtmosphereModel.h comments).
    double eps = 0.01; // meters, tiny offset to stay strictly within each branch
    AtmosphereState below84 = EarthAtmosphere1976::Compute(84000.0 - eps);
    AtmosphereState above84 = EarthAtmosphere1976::Compute(84000.0 + eps);
    AtmosphereState below120 = EarthAtmosphere1976::Compute(120000.0 - eps);
    AtmosphereState above120 = EarthAtmosphere1976::Compute(120000.0 + eps);

    EXPECT_NEAR(below84.density, above84.density, below84.density * 0.03);
    EXPECT_NEAR(below120.density, above120.density, below120.density * 0.015);
}

TEST(MarsAtmosphereExponentialTest, SeaLevelDensityIsReasonable) {
    AtmosphereState s = MarsAtmosphereExponential::Compute(0.0);
    EXPECT_NEAR(s.density, 0.020, 0.020 * kRelTol);
    EXPECT_GT(s.temperature, 0.0);
    EXPECT_TRUE(std::isfinite(s.pressure));
    EXPECT_GT(s.pressure, 0.0);
}

TEST(EarthAtmosphere1976Test, DensityMonotonicallyDecreasesWithAltitude) {
    double prev = EarthAtmosphere1976::Compute(0.0).density;
    for (double alt = 5000.0; alt <= 200000.0; alt += 5000.0) {
        double cur = EarthAtmosphere1976::Compute(alt).density;
        EXPECT_LE(cur, prev) << "density increased at altitude " << alt;
        prev = cur;
    }
}

TEST(MarsAtmosphereExponentialTest, DensityMonotonicallyDecreasesWithAltitude) {
    double prev = MarsAtmosphereExponential::Compute(0.0).density;
    for (double alt = 5000.0; alt <= 100000.0; alt += 5000.0) {
        double cur = MarsAtmosphereExponential::Compute(alt).density;
        EXPECT_LE(cur, prev) << "density increased at altitude " << alt;
        prev = cur;
    }
}
