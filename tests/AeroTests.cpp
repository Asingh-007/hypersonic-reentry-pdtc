#include <gtest/gtest.h>
#include "PanelGeometry.h"
#include "NewtonianAero.h"
#include "ShockExpansionAero.h"
#include "AeroRegimeDispatch.h"
#include "TestBodyGenerator.h"
#include "StlMeshLoader.h"
#include "SpacecraftGeometry.h"
#include "FlapHingeData.h"
#include "AeroCoefficientTable.h"
#include "AeroAngles.h"
#include "LatinHypercubeSampler.h"
#include "Kriging.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <cstdio>
#include <string>
#include <utility>

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

    AeroCoefficients c = model.evaluate(mesh, mapFlapAxesToGroupDeflections(0.0, 0.0, 0.0),
                                          0.0, 0.0, 5.0, moment_ref, S_ref, L_ref);
    // At alpha=beta=0 every panel (body + all 4 flat flap patches) is
    // exactly tangential to the flow, so CL should be ~exactly zero.
    EXPECT_NEAR(c.CL, 0.0, 1e-6);
}

// ---------------------------------------------------------------------
// 4-flap-group geometry + flap-axis mapping
// ---------------------------------------------------------------------

TEST(TestBodyGeneratorTest, FourFlapGroupsExistWithExpectedIdsNamesAndHingeAxes) {
    PanelMesh mesh = testutil::makeCylinderNoseFlapBody();
    const auto& groups = mesh.groups();
    ASSERT_EQ(groups.size(), 4u);

    const std::pair<int, std::string> expected[] = {
        {1, "fwd_left"}, {2, "fwd_right"}, {3, "aft_left"}, {4, "aft_right"}};
    for (const auto& [id, name] : expected) {
        auto it = groups.find(id);
        ASSERT_NE(it, groups.end()) << "missing group id " << id;
        EXPECT_EQ(it->second.name, name);
        EXPECT_TRUE(it->second.hinge_axis.isApprox(Eigen::Vector3d::UnitY(), 1e-12));
    }
}

TEST(AeroRegimeDispatchTest, FlapAxisMappingProducesExpectedSignConvention) {
    auto fwd_only = mapFlapAxesToGroupDeflections(0.2, 0.0, 0.0);
    EXPECT_NEAR(fwd_only[1], 0.2, 1e-12);
    EXPECT_NEAR(fwd_only[2], 0.2, 1e-12);
    EXPECT_NEAR(fwd_only[3], 0.0, 1e-12);
    EXPECT_NEAR(fwd_only[4], 0.0, 1e-12);

    auto diff_only = mapFlapAxesToGroupDeflections(0.0, 0.0, 0.1);
    EXPECT_NEAR(diff_only[3] - diff_only[4], 0.1, 1e-12);
    EXPECT_NEAR(diff_only[3] + diff_only[4], 0.0, 1e-12);
}

TEST(AeroRegimeDispatchTest, AftDifferentialFlapProducesNonzeroRollInExpectedDirection) {
    // Rather than hand-derive the expected roll sign (easy to get backwards),
    // this test numerically pins it down; future changes must keep it passing.
    PanelMesh mesh = testutil::makeCylinderNoseFlapBody();
    NewtonianAeroModel model;
    Eigen::Vector3d moment_ref(20.0, 0.0, 0.0);
    double S_ref = kPi * 4.5 * 4.5, L_ref = 9.0;

    auto flap_defl = mapFlapAxesToGroupDeflections(0.0, 0.0, /*aft_diff_rad=*/0.1);
    AeroCoefficients c = model.evaluate(mesh, flap_defl, /*alpha_rad=*/0.1, 0.0, 5.0,
                                          moment_ref, S_ref, L_ref);
    EXPECT_NE(c.Cl_roll, 0.0);

    // Reversing aft_diff_rad must reverse Cl_roll's sign (odd symmetry) --
    // robust regardless of which absolute sign convention was chosen.
    auto flap_defl_neg = mapFlapAxesToGroupDeflections(0.0, 0.0, -0.1);
    AeroCoefficients c_neg = model.evaluate(mesh, flap_defl_neg, 0.1, 0.0, 5.0,
                                              moment_ref, S_ref, L_ref);
    EXPECT_NEAR(c.Cl_roll, -c_neg.Cl_roll, 1e-9);
}

// ---------------------------------------------------------------------
// Real spacecraft geometry (SolidWorks STL import)
// ---------------------------------------------------------------------

TEST(SpacecraftGeometryTest, LoadsRealGeometryWithExpectedGroupsAndSaneReferenceValues) {
    const std::string geometry_dir = std::string(TESTS_SOURCE_DIR) + "/geometry";
    SpacecraftGeometry geo = LoadSpacecraftGeometry(geometry_dir);

    ASSERT_EQ(geo.mesh.groups().size(), 4u);
    for (int id : {1, 2, 3, 4}) {
        EXPECT_NE(geo.mesh.groups().find(id), geo.mesh.groups().end()) << "missing group id " << id;
    }
    EXPECT_GT(geo.mesh.panels().size(), 1000u);  // body (1804 tris) + 4 flaps (794 tris each)

    // Sane magnitude checks (real vehicle is tens of meters long, a few
    // meters in radius) -- not exact-value checks, since these are
    // geometric placeholders derived from the mesh, not fixed constants.
    EXPECT_GT(geo.body_length, 10.0);
    EXPECT_LT(geo.body_length, 200.0);
    EXPECT_GT(geo.body_radius, 0.5);
    EXPECT_LT(geo.body_radius, 50.0);
    EXPECT_GT(geo.S_ref, 0.0);
    EXPECT_NEAR(geo.S_ref, kPi * geo.body_radius * geo.body_radius, 1e-9);
    EXPECT_TRUE(geo.moment_ref.allFinite());
}

// ---------------------------------------------------------------------
// ShockExpansionAeroModel (tangent-wedge + Prandtl-Meyer)
// ---------------------------------------------------------------------

TEST(ShockExpansionAeroTest, PrandtlMeyerNuIsZeroAtMachOne) {
    EXPECT_NEAR(prandtlMeyerNu(1.0, 1.4), 0.0, 1e-9);
}

TEST(ShockExpansionAeroTest, PrandtlMeyerNuMatchesKnownGasDynamicsTableValue) {
    // Classical Prandtl-Meyer table value (e.g. NACA 1135 / Anderson
    // appendix), gamma=1.4: nu(2.0) ~= 26.38 deg.
    double nu_deg = prandtlMeyerNu(2.0, 1.4) * 180.0 / kPi;
    EXPECT_NEAR(nu_deg, 26.38, 0.05);
}

TEST(ShockExpansionAeroTest, InvertPrandtlMeyerRoundTrips) {
    for (double M : {1.2, 1.5, 2.0, 3.0, 5.0, 8.0}) {
        double nu = prandtlMeyerNu(M, 1.4);
        double M_back = invertPrandtlMeyer(nu, 1.4);
        EXPECT_NEAR(M_back, M, 1e-5) << "at M=" << M;
    }
}

TEST(ShockExpansionAeroTest, ObliqueShockBetaMatchesKnownGasDynamicsTableValue) {
    // Classical theta-beta-Mach chart value: M1=2.0, theta=10deg -> weak
    // shock beta ~= 39.3 deg (Anderson / NACA 1135).
    auto beta = solveWeakObliqueShockBeta(2.0, 10.0 * kPi / 180.0, 1.4);
    ASSERT_TRUE(beta.has_value());
    EXPECT_NEAR(*beta * 180.0 / kPi, 39.3, 0.3);
}

TEST(ShockExpansionAeroTest, ObliqueShockDetachesAndReturnsNulloptBeyondMaxDeflection) {
    // M1=1.5's maximum attached-shock deflection is ~12.1deg (theta-beta-M
    // chart) -- 30deg is well beyond detachment.
    auto beta = solveWeakObliqueShockBeta(1.5, 30.0 * kPi / 180.0, 1.4);
    EXPECT_FALSE(beta.has_value());
}

TEST(ShockExpansionAeroTest, SmallAngleObliqueShockCpMatchesLinearizedSupersonicTheory) {
    // For small theta, oblique-shock Cp should roughly agree with linearized
    // (Ackeret) theory: Cp ~= 2*theta/sqrt(M^2-1). 10% tolerance covers
    // linearized theory's own first-order error at theta=2deg.
    const double mach = 3.0;
    const double theta = 2.0 * kPi / 180.0;
    auto beta = solveWeakObliqueShockBeta(mach, theta, 1.4);
    ASSERT_TRUE(beta.has_value());
    const double Mn1 = mach * std::sin(*beta);
    const double Cp_exact = (4.0 / (1.4 + 1.0)) * (Mn1 * Mn1 - 1.0) / (mach * mach);
    const double Cp_linearized = 2.0 * theta / std::sqrt(mach * mach - 1.0);
    EXPECT_NEAR(Cp_exact, Cp_linearized, 0.10 * Cp_linearized);
}

TEST(ShockExpansionAeroTest, SymmetricBodyAtZeroAlphaZeroBetaZeroFlapGivesNearZeroCLAtSupersonicMach) {
    PanelMesh mesh = testutil::makeCylinderNoseFlapBody();
    ShockExpansionAeroModel model;
    Eigen::Vector3d moment_ref(20.0, 0.0, 0.0);
    double S_ref = kPi * 4.5 * 4.5, L_ref = 9.0;

    AeroCoefficients c = model.evaluate(mesh, mapFlapAxesToGroupDeflections(0.0, 0.0, 0.0),
                                          0.0, 0.0, /*mach=*/3.0, moment_ref, S_ref, L_ref);
    EXPECT_NEAR(c.CL, 0.0, 1e-6);
}

// ---------------------------------------------------------------------
// Mach-regime dispatch + transonic/subsonic placeholders
// ---------------------------------------------------------------------

TEST(AeroRegimeDispatchTest, RegimeDispatchPicksExpectedRegimeAtBoundaries) {
    EXPECT_EQ(classifyMachRegime(0.79), MachRegime::kSubsonic);
    EXPECT_EQ(classifyMachRegime(0.8), MachRegime::kTransonic);
    EXPECT_EQ(classifyMachRegime(1.19), MachRegime::kTransonic);
    EXPECT_EQ(classifyMachRegime(1.2), MachRegime::kSupersonic);
    EXPECT_EQ(classifyMachRegime(4.99), MachRegime::kSupersonic);
    EXPECT_EQ(classifyMachRegime(5.0), MachRegime::kHypersonic);
}

TEST(AeroRegimeDispatchTest, SubsonicPlaceholderZeroAlphaGivesHoernerBaselineCD) {
    double S_ref = kPi * 4.5 * 4.5, L_ref = 9.0;
    AeroCoefficients c = subsonicPlaceholderAero(0.0, 0.0, /*body_radius=*/4.5, /*body_length=*/40.0, S_ref, L_ref);
    EXPECT_NEAR(c.CD, 0.8, 1e-9);
    EXPECT_NEAR(c.Cl_roll, 0.0, 1e-12);
    EXPECT_NEAR(c.Cm, 0.0, 1e-12);
    EXPECT_NEAR(c.Cn_yaw, 0.0, 1e-12);
}

TEST(AeroRegimeDispatchTest, SubsonicPlaceholderCrossflowNormalForceGrowsWithAlpha) {
    double S_ref = kPi * 4.5 * 4.5, L_ref = 9.0;
    double prev_CN = -1.0;
    for (double alpha_deg : {10.0, 30.0, 60.0, 89.0}) {
        AeroCoefficients c = subsonicPlaceholderAero(alpha_deg * kPi / 180.0, 0.0, 4.5, 40.0, S_ref, L_ref);
        EXPECT_GT(c.CN, prev_CN);
        prev_CN = c.CN;
    }
}

TEST(AeroRegimeDispatchTest, TransonicPlaceholderMatchesSubsonicAtLowerBoundary) {
    PanelMesh mesh = testutil::makeCylinderNoseFlapBody();
    auto flap_defl = mapFlapAxesToGroupDeflections(0.0, 0.0, 0.0);
    Eigen::Vector3d moment_ref(20.0, 0.0, 0.0);
    double S_ref = kPi * 4.5 * 4.5, L_ref = 9.0;

    AeroCoefficients transonic = transonicPlaceholderAero(mesh, flap_defl, 0.3, 0.0, kMachSubsonicUpper,
                                                            moment_ref, S_ref, L_ref, 4.5, 40.0);
    AeroCoefficients subsonic = subsonicPlaceholderAero(0.3, 0.0, 4.5, 40.0, S_ref, L_ref);
    EXPECT_NEAR(transonic.CL, subsonic.CL, 1e-9);
    EXPECT_NEAR(transonic.CD, subsonic.CD, 1e-9);
}

TEST(AeroRegimeDispatchTest, TransonicPlaceholderMatchesSupersonicAtUpperBoundary) {
    PanelMesh mesh = testutil::makeCylinderNoseFlapBody();
    auto flap_defl = mapFlapAxesToGroupDeflections(0.0, 0.0, 0.0);
    Eigen::Vector3d moment_ref(20.0, 0.0, 0.0);
    double S_ref = kPi * 4.5 * 4.5, L_ref = 9.0;

    AeroCoefficients transonic = transonicPlaceholderAero(mesh, flap_defl, 0.3, 0.0, kMachTransonicUpper,
                                                            moment_ref, S_ref, L_ref, 4.5, 40.0);
    ShockExpansionAeroModel wedge;
    AeroCoefficients supersonic = wedge.evaluate(mesh, flap_defl, 0.3, 0.0, kMachTransonicUpper,
                                                   moment_ref, S_ref, L_ref);
    EXPECT_NEAR(transonic.CL, supersonic.CL, 1e-9);
    EXPECT_NEAR(transonic.CD, supersonic.CD, 1e-9);
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

TEST(StlMeshLoaderTest, WriteMeshToStlRoundTripsThroughLoadMeshFromStl) {
    Panel p{Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(2, 0, 0), Eigen::Vector3d(0, 3, 0)};
    p.recompute();
    std::string path = std::string(TESTS_SOURCE_DIR) + "/tests/fixtures/_tmp_write_test.stl";

    WriteMeshToStl(path, {p}, "round_trip_test");
    PanelMesh reloaded = LoadMeshFromStl(path, /*group_id=*/0);
    std::remove(path.c_str());

    ASSERT_EQ(reloaded.panels().size(), 1u);
    EXPECT_TRUE(reloaded.panels()[0].v0.isApprox(p.v0, 1e-6));
    EXPECT_TRUE(reloaded.panels()[0].v1.isApprox(p.v1, 1e-6));
    EXPECT_TRUE(reloaded.panels()[0].v2.isApprox(p.v2, 1e-6));
    EXPECT_NEAR(reloaded.panels()[0].area, p.area, 1e-6);
}

// ---------------------------------------------------------------------
// FlapHingeData (shared source of truth for real-vehicle flap geometry)
// ---------------------------------------------------------------------

TEST(FlapHingeDataTest, TableHasFourUniqueIdsMatchingFwdAftConvention) {
    const auto& table = FlapHingeTable();
    ASSERT_EQ(table.size(), 4u);
    std::vector<int> ids;
    for (const auto& f : table) ids.push_back(f.group_id);
    std::sort(ids.begin(), ids.end());
    EXPECT_EQ(ids, (std::vector<int>{1, 2, 3, 4}));
}

TEST(FlapHingeDataTest, DeflectingAFlapAboutItsHingeInCadFramePreservesRadiusFromHinge) {
    // Sanity check for GenerateDeflectedGeometry.cpp's use of FlapHingeData
    // directly in the CAD frame (not the model frame SpacecraftGeometry.cpp
    // uses) -- deflecting about (hinge_point_cad_mm, UnitY) should rotate
    // vertices rigidly, preserving their distance from the hinge line.
    const FlapHingeInfo& f = FlapHingeTable()[2];  // aft_top
    const std::string geometry_dir = std::string(TESTS_SOURCE_DIR) + "/geometry";
    PanelMesh mesh = LoadMeshFromStl(geometry_dir + "/" + f.stl_filename, f.group_id);

    PanelGroup g;
    g.id = f.group_id;
    g.hinge_point = f.hinge_point_cad_mm;
    g.hinge_axis = Eigen::Vector3d::UnitY();
    mesh.addGroup(g);

    const double deflection_rad = -9.77187 * kPi / 180.0;
    auto deflected = mesh.deflected({{f.group_id, deflection_rad}});
    ASSERT_EQ(deflected.size(), mesh.panels().size());

    const Eigen::Vector3d& v0_before = mesh.panels()[0].v0;
    const Eigen::Vector3d& v0_after = deflected[0].v0;
    auto radiusFromHinge = [&](const Eigen::Vector3d& v) {
        const Eigen::Vector3d r = v - f.hinge_point_cad_mm;
        return std::sqrt(r.x() * r.x() + r.z() * r.z());  // perpendicular to UnitY hinge axis
    };
    EXPECT_NEAR(radiusFromHinge(v0_before), radiusFromHinge(v0_after), 1e-6);
    EXPECT_NEAR(v0_before.y(), v0_after.y(), 1e-9);  // hinge axis is Y -- unaffected
}

// Fixture (tests/fixtures/aero_table_test.csv): 2^6=64-row complete grid
// with simple linear per-axis formulas (CL=0.01*alpha_deg, CD=0.3,
// Cl_roll=0.001*beta_deg+0.002*aft_diff_deg, Cm=-0.02*alpha_deg,
// Cn_yaw=0.001*beta_deg) chosen so expected values are trivially hand-verifiable.

TEST(AeroCoefficientTableTest, ExactMatchAtGridNodes) {
    AeroCoefficientTable table;
    ASSERT_TRUE(table.load(TestAeroTablePath()));
    AeroCoefficients c = table.interpolate(5.0, -10.0, -5.0, -5.0, -5.0, -5.0);
    EXPECT_NEAR(c.CL, -0.1, 1e-9);
    EXPECT_NEAR(c.CD, 0.3, 1e-9);
    EXPECT_NEAR(c.Cl_roll, -0.015, 1e-9);  // 0.001*(-5) + 0.002*(-5)
    EXPECT_NEAR(c.Cm, 0.2, 1e-9);
    EXPECT_NEAR(c.Cn_yaw, -0.005, 1e-9);
}

TEST(AeroCoefficientTableTest, SaneInterpolationBetweenNodes) {
    AeroCoefficientTable table;
    ASSERT_TRUE(table.load(TestAeroTablePath()));
    // Midpoint of every axis; fixture is designed so this lands exactly at 0.
    AeroCoefficients c = table.interpolate(7.5, 0.0, 0.0, 0.0, 0.0, 0.0);
    EXPECT_NEAR(c.CL, 0.0, 1e-9);
    EXPECT_NEAR(c.Cl_roll, 0.0, 1e-9);
    EXPECT_NEAR(c.Cm, 0.0, 1e-9);
    EXPECT_NEAR(c.Cn_yaw, 0.0, 1e-9);
}

TEST(AeroCoefficientTableTest, ClampsAtGridBoundaryForOutOfRangeQuery) {
    AeroCoefficientTable table;
    ASSERT_TRUE(table.load(TestAeroTablePath()));
    AeroCoefficients below = table.interpolate(1.0, -10.0, -5.0, -5.0, -5.0, -5.0);  // mach below grid min (5)
    AeroCoefficients at_min = table.interpolate(5.0, -10.0, -5.0, -5.0, -5.0, -5.0);
    EXPECT_NEAR(below.CL, at_min.CL, 1e-9);
}

TEST(AeroCoefficientTableTest, InterpolatesAcrossAllSixAxesIndependently) {
    AeroCoefficientTable table;
    ASSERT_TRUE(table.load(TestAeroTablePath()));
    // Varying aft_diff_deg alone (holding beta_deg=0) should change Cl_roll
    // via the 0.002*aft_diff_deg term -- exercises the newly-generalized
    // 6-axis interpolation path (the old 4-axis table had no such axis).
    AeroCoefficients lo = table.interpolate(7.5, 0.0, 0.0, 0.0, 0.0, -5.0);
    AeroCoefficients hi = table.interpolate(7.5, 0.0, 0.0, 0.0, 0.0, 5.0);
    EXPECT_NEAR(hi.Cl_roll - lo.Cl_roll, 0.02, 1e-9);  // 0.002*(5-(-5))
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
