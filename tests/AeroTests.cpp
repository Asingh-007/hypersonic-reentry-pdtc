#include <gtest/gtest.h>
#include "PanelGeometry.h"
#include "NewtonianAero.h"
#include "TestBodyGenerator.h"
#include "StlMeshLoader.h"
#include "AeroCoefficientTable.h"
#include "AeroAngles.h"
#include "LatinHypercubeSampler.h"
#include "Kriging.h"
#include <cmath>
#include <fstream>
#include <cstdio>

using namespace aero_model;

// MSVC only defines M_PI when _USE_MATH_DEFINES is set before the first
// <cmath>/<math.h> include anywhere in the translation unit -- use an
// explicit local constant instead (same fix as main.cpp).
constexpr double kPi = 3.14159265358979323846;

namespace {
std::string TestAeroTablePath() {
    return std::string(TESTS_SOURCE_DIR) + "/tests/fixtures/aero_table_test.csv";
}
}

TEST(PanelGeometryTest, TriangleAreaNormalCentroidCorrectness) {
    Panel p{Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1, 0, 0), Eigen::Vector3d(0, 1, 0)};
    p.recompute();
    EXPECT_NEAR(p.area, 0.5, 1e-12);
    EXPECT_TRUE(p.centroid.isApprox(Eigen::Vector3d(1.0 / 3.0, 1.0 / 3.0, 0.0), 1e-12));
    EXPECT_TRUE(p.normal.isApprox(Eigen::Vector3d(0, 0, 1), 1e-12));
}

TEST(PanelGeometryTest, RodriguesFlapDeflectionRotatesAboutHingeCorrectly) {
    PanelMesh mesh;
    PanelGroup g;
    g.id = 1;
    g.hinge_point = Eigen::Vector3d(0, 0, 0);
    g.hinge_axis = Eigen::Vector3d::UnitY();
    mesh.addGroup(g);

    Panel p{Eigen::Vector3d(1, 0, 0), Eigen::Vector3d(2, 0, 0), Eigen::Vector3d(1, 1, 0)};
    p.group_id = 1;
    p.recompute();
    mesh.addPanel(p);

    auto deflected = mesh.deflected({{1, kPi / 2.0}}); // 90 degrees about Y axis
    ASSERT_EQ(deflected.size(), 1u);
    // Hand-derived via Rodrigues' formula (also matches the standard
    // rotation-about-Y matrix): (1,0,0)->(0,0,-1), (2,0,0)->(0,0,-2),
    // (1,1,0)->(0,1,-1).
    EXPECT_TRUE(deflected[0].v0.isApprox(Eigen::Vector3d(0, 0, -1), 1e-9));
    EXPECT_TRUE(deflected[0].v1.isApprox(Eigen::Vector3d(0, 0, -2), 1e-9));
    EXPECT_TRUE(deflected[0].v2.isApprox(Eigen::Vector3d(0, 1, -1), 1e-9));
}

TEST(NewtonianAeroModelTest, CpMaxAsymptotesToAndersonLimitAtHighMach) {
    NewtonianAeroModel model(1.4);
    EXPECT_NEAR(model.cpMax(1000.0), 1.839, 1e-2);
}

TEST(NewtonianAeroModelTest, SymmetricBodyAtZeroAlphaZeroBetaZeroFlapGivesNearZeroCL) {
    PanelMesh mesh = testutil::makeCylinderNoseFlapBody();
    NewtonianAeroModel model;
    Eigen::Vector3d moment_ref(20.0, 0.0, 0.0);
    double S_ref = kPi * 4.5 * 4.5, L_ref = 9.0;

    AeroCoefficients c = model.evaluate(mesh, {{1, 0.0}}, 0.0, 0.0, 5.0, moment_ref, S_ref, L_ref);
    // At alpha=beta=0, the flow is tangential to both the axisymmetric
    // body-of-revolution panels and the flat (z=const) flap panels, so
    // every panel's sin(theta) is exactly 0 -- CL should be ~exactly zero,
    // not just approximately.
    EXPECT_NEAR(c.CL, 0.0, 1e-6);
}

TEST(AeroAnglesTest, KnownCaseIdentityQuaternionPureSideslip) {
    // fpa=0, v_azi=0 -> velocity purely along local y-axis; with q=identity
    // (body aligned with local frame), the wind vector is exactly
    // (0,-1,0) in body axes -- pure +/-90deg sideslip. beta is well-defined
    // (-90deg) and checked below. alpha is NOT checked: freestreamDirectionBody
    // collapses to (0, sin(beta), 0) at beta=+/-90deg regardless of alpha (a
    // gimbal-lock-like singularity of this alpha/beta parameterization), so
    // atan2(~0,~0) here is mathematically indeterminate and its sign is at
    // the mercy of floating-point negative-zero propagation through the
    // quaternion rotation -- not a meaningful behavior to pin down.
    Eigen::Quaterniond identity(1.0, 0.0, 0.0, 0.0);
    AeroAngles angles = ComputeAeroAngles(identity, 0.0, 0.0);
    EXPECT_NEAR(angles.beta_rad, -kPi / 2.0, 1e-9);
}

TEST(AeroAnglesTest, ConsistentWithFreestreamDirectionBodyRoundTrip) {
    // ComputeAeroAngles' output, fed back into the same body-wind-direction
    // formula NewtonianAeroModel uses internally, should reproduce the
    // wind vector it was derived from.
    Eigen::Quaterniond q(Eigen::AngleAxisd(0.3, Eigen::Vector3d(0.2, 0.5, -0.1).normalized()));
    double fpa = -0.2, v_azi = 0.7;

    AeroAngles angles = ComputeAeroAngles(q, fpa, v_azi);

    Eigen::Vector3d v_hat_local(std::sin(fpa), std::cos(fpa) * std::cos(v_azi), std::cos(fpa) * std::sin(v_azi));
    Eigen::Vector3d expected_wind_body = q.conjugate() * (-v_hat_local);

    Eigen::Vector3d reconstructed(-std::cos(angles.alpha_rad) * std::cos(angles.beta_rad),
                                  std::sin(angles.beta_rad),
                                  std::sin(angles.alpha_rad) * std::cos(angles.beta_rad));
    EXPECT_TRUE(reconstructed.isApprox(expected_wind_body, 1e-9));
}

namespace {
std::string WriteTempFile(const std::string& contents) {
    std::string path = std::string(TESTS_SOURCE_DIR) + "/tests/fixtures/_tmp_stl_test.stl";
    std::ofstream out(path);
    out << contents;
    out.close();
    return path;
}
}

TEST(StlMeshLoaderTest, RoundTripsATinyHandWrittenAsciiStlFixture) {
    std::string stl_text =
        "solid test\n"
        "facet normal 0 0 1\n"
        "outer loop\n"
        "vertex 0 0 0\n"
        "vertex 1 0 0\n"
        "vertex 0 1 0\n"
        "endloop\n"
        "endfacet\n"
        "endsolid test\n";
    std::string path = WriteTempFile(stl_text);

    PanelMesh mesh = LoadMeshFromStl(path, /*group_id=*/0);
    std::remove(path.c_str());

    ASSERT_EQ(mesh.panels().size(), 1u);
    EXPECT_TRUE(mesh.panels()[0].v0.isApprox(Eigen::Vector3d(0, 0, 0), 1e-12));
    EXPECT_TRUE(mesh.panels()[0].v1.isApprox(Eigen::Vector3d(1, 0, 0), 1e-12));
    EXPECT_TRUE(mesh.panels()[0].v2.isApprox(Eigen::Vector3d(0, 1, 0), 1e-12));
    EXPECT_NEAR(mesh.panels()[0].area, 0.5, 1e-9);
}

TEST(StlMeshLoaderTest, MalformedStlThrows) {
    std::string stl_text =
        "solid test\n"
        "facet normal 0 0 1\n"
        "outer loop\n"
        "vertex 0 0\n"  // missing z coordinate
        "vertex 1 0 0\n"
        "vertex 0 1 0\n"
        "endloop\n"
        "endfacet\n"
        "endsolid test\n";
    std::string path = WriteTempFile(stl_text);
    EXPECT_THROW(LoadMeshFromStl(path, 0), std::runtime_error);
    std::remove(path.c_str());
}

TEST(AeroCoefficientTableTest, ExactMatchAtGridNodes) {
    AeroCoefficientTable table;
    ASSERT_TRUE(table.load(TestAeroTablePath()));
    AeroCoefficients c = table.interpolate(5.0, -10.0, -5.0, -5.0);
    EXPECT_NEAR(c.CL, -0.1, 1e-9);
    EXPECT_NEAR(c.CD, 0.3, 1e-9);
    EXPECT_NEAR(c.Cl_roll, -0.005, 1e-9);
    EXPECT_NEAR(c.Cm, 0.2, 1e-9);
    EXPECT_NEAR(c.Cn_yaw, -0.005, 1e-9);
}

TEST(AeroCoefficientTableTest, SaneInterpolationBetweenNodes) {
    AeroCoefficientTable table;
    ASSERT_TRUE(table.load(TestAeroTablePath()));
    // Midpoint of every axis; fixture is designed so this lands exactly at 0.
    AeroCoefficients c = table.interpolate(7.5, 0.0, 0.0, 0.0);
    EXPECT_NEAR(c.CL, 0.0, 1e-9);
    EXPECT_NEAR(c.Cl_roll, 0.0, 1e-9);
    EXPECT_NEAR(c.Cm, 0.0, 1e-9);
    EXPECT_NEAR(c.Cn_yaw, 0.0, 1e-9);
}

TEST(AeroCoefficientTableTest, ClampsAtGridBoundaryForOutOfRangeQuery) {
    AeroCoefficientTable table;
    ASSERT_TRUE(table.load(TestAeroTablePath()));
    AeroCoefficients below = table.interpolate(1.0, -10.0, -5.0, -5.0); // mach below grid min (5)
    AeroCoefficients at_min = table.interpolate(5.0, -10.0, -5.0, -5.0);
    EXPECT_NEAR(below.CL, at_min.CL, 1e-9);
}

TEST(LatinHypercubeSamplerTest, SampleProducesCorrectShapeAndUnitBounds) {
    LatinHypercubeSampler lhs(42);
    Eigen::MatrixXd design = lhs.sample(10, 3, 500);
    EXPECT_EQ(design.rows(), 10);
    EXPECT_EQ(design.cols(), 3);
    for (int i = 0; i < design.rows(); ++i) {
        for (int j = 0; j < design.cols(); ++j) {
            EXPECT_GE(design(i, j), 0.0);
            EXPECT_LE(design(i, j), 1.0);
        }
    }
}

TEST(LatinHypercubeSamplerTest, ScaleToBoundsMapsCorrectly) {
    Eigen::MatrixXd unit(2, 2);
    unit << 0.0, 1.0,
            1.0, 0.0;
    std::vector<DesignVariable> vars = {{"a", 0.0, 10.0}, {"b", -5.0, 5.0}};
    Eigen::MatrixXd scaled = LatinHypercubeSampler::scaleToBounds(unit, vars);
    EXPECT_NEAR(scaled(0, 0), 0.0, 1e-9);
    EXPECT_NEAR(scaled(0, 1), 5.0, 1e-9);
    EXPECT_NEAR(scaled(1, 0), 10.0, 1e-9);
    EXPECT_NEAR(scaled(1, 1), -5.0, 1e-9);
}

TEST(LatinHypercubeSamplerTest, AugmentWithFixedPointsAppendsCorrectly) {
    Eigen::MatrixXd samples(2, 2);
    samples << 1.0, 2.0, 3.0, 4.0;
    Eigen::MatrixXd fixed(1, 2);
    fixed << 9.0, 9.0;
    Eigen::MatrixXd out = LatinHypercubeSampler::augmentWithFixedPoints(samples, fixed);
    EXPECT_EQ(out.rows(), 3);
    EXPECT_NEAR(out(2, 0), 9.0, 1e-9);
    EXPECT_NEAR(out(2, 1), 9.0, 1e-9);
}

// ---------------------------------------------------------------------
// Kriging
// ---------------------------------------------------------------------

TEST(UniversalKrigingTest, FitsAndPredictsAKnownSyntheticFunctionReasonablyWell) {
    // f(x) = x^2, zero trend function (so the GP has to explain everything).
    UniversalKriging::TrendFn zero_trend = [](const Eigen::VectorXd&) { return 0.0; };
    Eigen::VectorXd length_scales(1);
    length_scales << 2.0;
    UniversalKriging krig(zero_trend, length_scales, /*sigma_f=*/10.0, /*sigma_n=*/1e-3);

    Eigen::MatrixXd X(5, 1);
    X << 0.0, 1.0, 2.0, 3.0, 4.0;
    Eigen::VectorXd y(5);
    y << 0.0, 1.0, 4.0, 9.0, 16.0;
    krig.fit(X, y);

    Eigen::VectorXd query(1);
    query << 2.0;
    auto pred = krig.predict(query);
    EXPECT_NEAR(pred.mean, 4.0, 0.1); // near-exact recovery at a training point

    double rmse = krig.looCV_RMSE();
    EXPECT_TRUE(std::isfinite(rmse));
    EXPECT_GE(rmse, 0.0);
}
