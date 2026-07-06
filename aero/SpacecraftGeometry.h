#ifndef SPACECRAFT_GEOMETRY_H
#define SPACECRAFT_GEOMETRY_H
// SpacecraftGeometry.h
//
// Loads the real, SolidWorks-exported spacecraft geometry under geometry/
// (5 ASCII STL files: 1 body + 4 independently-hinged flaps) into a single
// PanelMesh with group ids matching AeroRegimeDispatch.h's convention
// (1=fwd_top, 2=fwd_bottom driven by fwd_sym; 3=aft_top, 4=aft_bottom
// driven by aft_sym+/-0.5*aft_diff). This vehicle's flap pairs are
// top/bottom, not left/right like the placeholder TestBodyGenerator body.
//
// The STL export's own axes (X/Z lateral, Y the body's long axis, units
// mm) don't match this codebase's convention (X nose-tail, meters), so
// every vertex is remapped/rescaled on import: model_x=y_cad/1000 (nose
// at +Y_cad -- confirmed by the body tapering to a point over its last
// ~10m), model_y=-x_cad/1000, model_z=z_cad/1000 (chosen to keep the
// frame right-handed).

#include "PanelGeometry.h"
#include <Eigen/Dense>
#include <string>

namespace aero_model {

struct SpacecraftGeometry {
    PanelMesh mesh;
    double S_ref = 0.0;         // frontal area, m^2 -- PLACEHOLDER (pi * body_radius^2)
    double L_ref = 0.0;         // reference length, m -- PLACEHOLDER (2 * body_radius)
    double body_radius = 0.0;   // effective body radius, m -- PLACEHOLDER, for the subsonic Hoerner/Allen-Perkins correlations
    double body_length = 0.0;   // body length, m
    Eigen::Vector3d moment_ref = Eigen::Vector3d::Zero();  // bounding-box-center CG PLACEHOLDER, body frame, m
};

// geometry_dir: path to the folder containing the 5 named STL files
// (unchanged from the SolidWorks export's default names). S_ref/L_ref/
// body_radius/moment_ref are geometric placeholders derived from the body
// mesh's bounding box -- true CG needs mass distribution, not pure geometry.
SpacecraftGeometry LoadSpacecraftGeometry(const std::string& geometry_dir);

}  // namespace aero_model

#endif // SPACECRAFT_GEOMETRY_H
