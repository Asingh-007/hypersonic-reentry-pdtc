#include <gtest/gtest.h>
#include "full_loop_nonconvex.h"
#include <cmath>

namespace {
bool HasConeTag(const clarabel::SupportedConeT<double>& cone,
                 clarabel::SupportedConeT<double>::Tag tag) {
    return cone.tag == tag;
}
}  // namespace

TEST(FullLoopNonconvexTest, ThrustMagnitudeRowsAreSecondOrderConeDimFour) {
    auto g = guidance_scp::buildThrustMagnitudeSocRows(0, 1, 2, 3, 0.0, Eigen::Vector3d::Zero());
    EXPECT_TRUE(HasConeTag(g.cone, clarabel::SupportedConeT<double>::Tag::SecondOrderConeT));
    ASSERT_EQ(g.rows.size(), 4u);
    // Row 0 (scalar): -sigma_u <= 0, i.e. coeffs={sigma_u_col:-1}, rhs=0 (zero reference).
    ASSERT_EQ(g.rows[0].coeffs.size(), 1u);
    EXPECT_EQ(g.rows[0].coeffs[0].first, 0);
    EXPECT_NEAR(g.rows[0].coeffs[0].second, -1.0, 1e-12);
    EXPECT_NEAR(g.rows[0].rhs, 0.0, 1e-12);
}

TEST(FullLoopNonconvexTest, ThrustMagnitudeRowsOffsetRhsByReference) {
    // Decision-variable columns are DEVIATIONS from a reference -- a
    // nonzero (sigma_u_ref, u_vec_ref) must shift each row's rhs so the
    // cone is enforced on the ABSOLUTE (ref + delta) quantities.
    Eigen::Vector3d u_vec_ref(1.0, -2.0, 3.0);
    auto g = guidance_scp::buildThrustMagnitudeSocRows(0, 1, 2, 3, 5.0, u_vec_ref);
    EXPECT_NEAR(g.rows[0].rhs, 5.0, 1e-12);
    EXPECT_NEAR(g.rows[1].rhs, 1.0, 1e-12);
    EXPECT_NEAR(g.rows[2].rhs, -2.0, 1e-12);
    EXPECT_NEAR(g.rows[3].rhs, 3.0, 1e-12);
}

TEST(FullLoopNonconvexTest, GimbalConeRowsAreSecondOrderConeDimFour) {
    Eigen::Vector3d e1_ref(-1.0, 0.0, 0.0);
    auto g = guidance_scp::buildGimbalConeSocRows(0, 1, 2, e1_ref, 0.1396, Eigen::Vector3d::Zero());
    EXPECT_TRUE(HasConeTag(g.cone, clarabel::SupportedConeT<double>::Tag::SecondOrderConeT));
    ASSERT_EQ(g.rows.size(), 4u);
}

TEST(FullLoopNonconvexTest, GimbalConeScalarRowUsesCorrectAxisAndScale) {
    // cos(deltaE_max)*||u|| <= u.e1_ref  ==  ||u|| <= (1/cos(deltaE_max)) * u.e1_ref
    // ==  -(1/cos(deltaE_max))*u.e1_ref <= 0 (the SOC's "scalar" row, coeffs*x<=rhs, zero reference).
    Eigen::Vector3d e1_ref(0.0, 1.0, 0.0);
    double deltaE_max = 0.2;
    auto g = guidance_scp::buildGimbalConeSocRows(10, 11, 12, e1_ref, deltaE_max, Eigen::Vector3d::Zero());
    double inv_cos = 1.0 / std::cos(deltaE_max);

    ASSERT_EQ(g.rows[0].coeffs.size(), 3u);
    // e1_ref = (0,1,0) -> only the uy_col (index 11) term should be nonzero.
    for (const auto& [col, val] : g.rows[0].coeffs) {
        if (col == 11) {
            EXPECT_NEAR(val, -inv_cos, 1e-9);
        } else {
            EXPECT_NEAR(val, 0.0, 1e-12);
        }
    }
}

TEST(FullLoopNonconvexTest, SymmetricBoxRowsAreNonnegativeConeWithCorrectSignConvention) {
    // var_ref=2.0, bound=5.0: upper row enforces delta_var <= 3.0,
    // lower row enforces -delta_var <= 7.0 (i.e. delta_var >= -7.0).
    auto rows = guidance_scp::buildSymmetricBoxRows(/*var_col=*/4, /*var_ref=*/2.0, /*bound=*/5.0);
    ASSERT_EQ(rows.size(), 2u);

    for (const auto& g : rows) {
        EXPECT_TRUE(HasConeTag(g.cone, clarabel::SupportedConeT<double>::Tag::NonnegativeConeT));
        ASSERT_EQ(g.rows.size(), 1u);
    }

    const auto& upper = rows[0].rows[0];
    ASSERT_EQ(upper.coeffs.size(), 1u);
    EXPECT_EQ(upper.coeffs[0].first, 4);
    EXPECT_NEAR(upper.coeffs[0].second, 1.0, 1e-12);
    EXPECT_NEAR(upper.rhs, 3.0, 1e-12);

    const auto& lower = rows[1].rows[0];
    ASSERT_EQ(lower.coeffs.size(), 1u);
    EXPECT_EQ(lower.coeffs[0].first, 4);
    EXPECT_NEAR(lower.coeffs[0].second, -1.0, 1e-12);
    EXPECT_NEAR(lower.rhs, 7.0, 1e-12);
}
