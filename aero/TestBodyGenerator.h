#ifndef TEST_BODY_GENERATOR_H
#define TEST_BODY_GENERATOR_H
// TestBodyGenerator.h
//
// Procedurally builds a crude cylinder-plus-nose-cap body with 4 flap
// panel groups, purely so the aero module can be exercised end-to-end
// without a real OpenVSP export. This is the PLACEHOLDER vehicle used by
// aero_table_gen until a real triangulated OML import exists (see StlMeshLoader.h).

#include <cmath>

#include "PanelGeometry.h"

namespace aero_model::testutil {

// MSVC only defines M_PI when _USE_MATH_DEFINES is set before the first
// <cmath>/<math.h> include anywhere in the translation unit, which is
// fragile to guarantee across header inclusion order -- use an explicit
// local constant instead (same fix already applied in main.cpp).
inline constexpr double kPi = 3.14159265358979323846;

inline PanelMesh makeCylinderNoseFlapBody(int n_theta = 24, int n_axial = 10,
                                           double radius = 4.5, double length = 40.0,
                                           double nose_length = 6.0) {
    PanelMesh mesh;
    const double dtheta = 2.0 * kPi / n_theta;

    // Cylindrical body: x in [0, length - nose_length], nose tip at
    // x = length (nose points in +x, matching NewtonianAeroModel's
    // freestream convention).
    const double cyl_length = length - nose_length;
    for (int i = 0; i < n_axial; ++i) {
        const double x0 = cyl_length * static_cast<double>(i) / n_axial;
        const double x1 = cyl_length * static_cast<double>(i + 1) / n_axial;
        for (int j = 0; j < n_theta; ++j) {
            const double th0 = j * dtheta;
            const double th1 = (j + 1) * dtheta;

            Eigen::Vector3d p00(x0, radius * std::cos(th0), radius * std::sin(th0));
            Eigen::Vector3d p01(x0, radius * std::cos(th1), radius * std::sin(th1));
            Eigen::Vector3d p10(x1, radius * std::cos(th0), radius * std::sin(th0));
            Eigen::Vector3d p11(x1, radius * std::cos(th1), radius * std::sin(th1));

            Panel t1{p00, p10, p11};
            t1.group_id = 0;
            t1.recompute();
            Panel t2{p00, p11, p01};
            t2.group_id = 0;
            t2.recompute();
            mesh.addPanel(t1);
            mesh.addPanel(t2);
        }
    }

    // Blunt nose cap: simple hemispherical-ish cap from cyl_length to
    // length via a cosine profile (blunt like Starship's nose, not pointed).
    const int n_nose_rings = 8;
    for (int i = 0; i < n_nose_rings; ++i) {
        const double f0 = static_cast<double>(i) / n_nose_rings;
        const double f1 = static_cast<double>(i + 1) / n_nose_rings;
        const double x0 = cyl_length + nose_length * std::sin(f0 * kPi / 2.0);
        const double x1 = cyl_length + nose_length * std::sin(f1 * kPi / 2.0);
        const double r0 = radius * std::cos(f0 * kPi / 2.0);
        const double r1 = radius * std::cos(f1 * kPi / 2.0);
        for (int j = 0; j < n_theta; ++j) {
            const double th0 = j * dtheta;
            const double th1 = (j + 1) * dtheta;

            Eigen::Vector3d p00(x0, r0 * std::cos(th0), r0 * std::sin(th0));
            Eigen::Vector3d p01(x0, r0 * std::cos(th1), r0 * std::sin(th1));
            Eigen::Vector3d p10(x1, r1 * std::cos(th0), r1 * std::sin(th0));
            Eigen::Vector3d p11(x1, r1 * std::cos(th1), r1 * std::sin(th1));

            Panel t1{p00, p10, p11};
            t1.group_id = 0;
            t1.recompute();
            Panel t2{p00, p11, p01};
            t2.group_id = 0;
            t2.recompute();
            mesh.addPanel(t1);
            mesh.addPanel(t2);
        }
    }

    // Four independent, deflectable flap groups (Starship-like layout): two
    // forward flaps near the nose/cylinder junction on the "top" (+z) side,
    // two aft flaps near the tail on the "bottom" (-z) side, each pair
    // split left(+y)/right(-y). group_id 1=fwd_left, 2=fwd_right,
    // 3=aft_left, 4=aft_right, matching AeroRegimeDispatch.h's fwd_sym/aft_sym/aft_diff naming.
    struct FlapSpec {
        int id;
        const char* name;
        double hinge_x, hinge_z;
        double y_center, y_half_span;
        double chord_sign;  // +1: patch extends toward +x (nose-ward); -1: toward -x (tail-ward)
    };

    const double flap_chord = 6.0;
    const double flap_half_span = 5.0;
    const double y_offset = 6.0;  // left/right patch centers, so they don't overlap at y=0
    const double fwd_hinge_x = cyl_length - 4.0;  // just aft of the nose/cylinder junction
    const double aft_hinge_x = 2.0;               // same position as the original single aft flap

    const FlapSpec flap_specs[4] = {
        {1, "fwd_left",  fwd_hinge_x,  radius, +y_offset, flap_half_span, +1.0},
        {2, "fwd_right", fwd_hinge_x,  radius, -y_offset, flap_half_span, +1.0},
        {3, "aft_left",  aft_hinge_x, -radius, +y_offset, flap_half_span, -1.0},
        {4, "aft_right", aft_hinge_x, -radius, -y_offset, flap_half_span, -1.0},
    };

    for (const FlapSpec& spec : flap_specs) {
        PanelGroup flap_group;
        flap_group.id = spec.id;
        flap_group.name = spec.name;
        flap_group.hinge_point = Eigen::Vector3d(spec.hinge_x, spec.y_center, spec.hinge_z);
        flap_group.hinge_axis = Eigen::Vector3d::UnitY();
        mesh.addGroup(flap_group);

        const int n_flap = 4;
        for (int j = 0; j < n_flap; ++j) {
            const double y0 = spec.y_center - spec.y_half_span + 2.0 * spec.y_half_span * j / n_flap;
            const double y1 = spec.y_center - spec.y_half_span + 2.0 * spec.y_half_span * (j + 1) / n_flap;
            const double x_far = spec.hinge_x + spec.chord_sign * flap_chord;

            Eigen::Vector3d p00(spec.hinge_x, y0, spec.hinge_z);
            Eigen::Vector3d p01(spec.hinge_x, y1, spec.hinge_z);
            Eigen::Vector3d p10(x_far, y0, spec.hinge_z);
            Eigen::Vector3d p11(x_far, y1, spec.hinge_z);

            Panel t1{p00, p10, p11};
            t1.group_id = spec.id;
            t1.recompute();
            Panel t2{p00, p11, p01};
            t2.group_id = spec.id;
            t2.recompute();
            mesh.addPanel(t1);
            mesh.addPanel(t2);
        }
    }

    return mesh;
}

}  // namespace aero_model::testutil

#endif // TEST_BODY_GENERATOR_H
