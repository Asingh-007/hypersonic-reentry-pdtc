#ifndef FLAP_HINGE_DATA_H
#define FLAP_HINGE_DATA_H
// FlapHingeData.h
//
// Shared source of truth for the 4 real spacecraft flaps' file names and
// hinge geometry, in the STL export's own (CAD) frame -- units mm, hinge
// axis along CAD_Y (the body's long axis; see SpacecraftGeometry.h for why).
// Used by SpacecraftGeometry.cpp (which converts to this codebase's model
// frame for aero analysis) and GenerateDeflectedGeometry.cpp (which rotates
// flap STLs directly in the CAD frame for CFD geometry export), so the
// hinge point measurements only live in one place.
//
// group_id matches AeroRegimeDispatch.h's convention: 1/2 = forward pair
// (driven by fwd_sym), 3/4 = aft pair (driven by aft_sym +/- 0.5*aft_diff).
// Hinge points are the measured centroids of each flap's mounting edge (a
// tight, near-planar cluster of vertices found by inspecting the raw STL
// data); the hinge axis is CAD_Y for all 4, since each flap's span runs
// the full length of its mounting edge parallel to the body axis.

#include <array>

#include <Eigen/Dense>

namespace aero_model {

struct FlapHingeInfo {
    int group_id;
    const char* name;
    const char* stl_filename;      // relative to the geometry/ directory
    Eigen::Vector3d hinge_point_cad_mm;
};

inline const std::array<FlapHingeInfo, 4>& FlapHingeTable() {
    static const std::array<FlapHingeInfo, 4> table = {{
        {1, "fwd_top",    "Spacecraft Flap 1.STL", {4549.0, 37600.0, 12230.8}},
        {2, "fwd_bottom", "Spacecraft Flap 2.STL", {4549.0, 37600.0, 3223.6}},
        {3, "aft_top",    "Spacecraft Flap 3.STL", {4546.3, 4000.0, 12230.8}},
        {4, "aft_bottom", "Spacecraft Flap 4.STL", {4551.2, 4000.0, 3223.5}},
    }};
    return table;
}

inline const char* kSpacecraftBodyStlFilename = "Spacecraft Assembly - Spacecraft Body-1.STL";

}  // namespace aero_model

#endif // FLAP_HINGE_DATA_H
