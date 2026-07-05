#ifndef TEST_BODY_GENERATOR_H
#define TEST_BODY_GENERATOR_H
// TestBodyGenerator.h
//
// Procedurally builds a crude cylinder-plus-nose-cap body with one flap
// panel group, purely so PanelGeometry/NewtonianAero/Kriging/
// LatinHypercubeSampler can be exercised end-to-end without an actual
// OpenVSP export. This is the PLACEHOLDER vehicle used by aero_table_gen
// until a real triangulated OML import exists -- see StlMeshLoader.h for
// that (StlMeshLoader supports OpenVSP/OpenRocket/CAD-exported STL files;
// PanelMesh::addPanel() takes raw triangle vertices either way, so the rest
// of the pipeline is unaffected by the geometry source).

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

    // One flap: a flat rectangular panel patch near the tail, hinged about
    // a y-parallel axis, group_id = 1. Positioned on the -z ("bottom")
    // side of the body, similar to a Starship aft flap.
    PanelGroup flap_group;
    flap_group.id = 1;
    flap_group.name = "aft_flap";
    flap_group.hinge_point = Eigen::Vector3d(2.0, 0.0, -radius);
    flap_group.hinge_axis = Eigen::Vector3d::UnitY();
    mesh.addGroup(flap_group);

    const double flap_chord = 6.0;
    const double flap_span_half = 5.0;
    const int n_flap = 4;
    for (int j = 0; j < n_flap; ++j) {
        const double y0 = -flap_span_half + 2.0 * flap_span_half * j / n_flap;
        const double y1 = -flap_span_half + 2.0 * flap_span_half * (j + 1) / n_flap;

        Eigen::Vector3d p00(2.0, y0, -radius);
        Eigen::Vector3d p01(2.0, y1, -radius);
        Eigen::Vector3d p10(2.0 - flap_chord, y0, -radius);
        Eigen::Vector3d p11(2.0 - flap_chord, y1, -radius);

        Panel t1{p00, p10, p11};
        t1.group_id = 1;
        t1.recompute();
        Panel t2{p00, p11, p01};
        t2.group_id = 1;
        t2.recompute();
        mesh.addPanel(t1);
        mesh.addPanel(t2);
    }

    return mesh;
}

}  // namespace aero_model::testutil

#endif // TEST_BODY_GENERATOR_H
