#include "SpacecraftGeometry.h"
#include "StlMeshLoader.h"

#include <filesystem>
#include <limits>

namespace aero_model {

namespace {

// MSVC only defines M_PI when _USE_MATH_DEFINES is set before the first
// <cmath>/<math.h> include anywhere in the translation unit -- use an
// explicit local constant instead (same fix as elsewhere in this repo).
constexpr double kPi = 3.14159265358979323846;

constexpr double kMmToM = 1.0 / 1000.0;  // STL export uses mm; this codebase uses meters throughout.

// model_x=CAD_Y (nose-tail, +model_x=nose), model_y=-CAD_X, model_z=CAD_Z --
// verified right-handed (model_x cross model_y = model_z).
Eigen::Vector3d ToModelFrame(const Eigen::Vector3d& cad) {
    return Eigen::Vector3d(cad.y(), -cad.x(), cad.z()) * kMmToM;
}

void AppendTransformed(const PanelMesh& mesh, PanelMesh& out) {
    for (Panel p : mesh.panels()) {
        p.v0 = ToModelFrame(p.v0);
        p.v1 = ToModelFrame(p.v1);
        p.v2 = ToModelFrame(p.v2);
        p.recompute();
        out.addPanel(p);
    }
}

struct FlapFile {
    int id;
    const char* name;
    const char* filename;
    Eigen::Vector3d hinge_cad_mm;
};

}  // namespace

SpacecraftGeometry LoadSpacecraftGeometry(const std::string& geometry_dir) {
    namespace fs = std::filesystem;
    const fs::path dir(geometry_dir);

    SpacecraftGeometry geo;

    // --- Body (group 0, fixed) ---
    const PanelMesh body_raw = LoadMeshFromStl((dir / "Spacecraft Assembly - Spacecraft Body-1.STL").string(), 0);
    AppendTransformed(body_raw, geo.mesh);

    // --- 4 flaps. Hinge points below are the measured centroids of each
    // flap's mounting edge (a tight, near-planar cluster of vertices found
    // by inspecting the raw STL data); hinge_axis is the body's nose-tail
    // direction for all 4, since each flap's span runs the full length of
    // its mounting edge parallel to the body axis. ---
    const FlapFile flap_files[4] = {
        {1, "fwd_top",    "Spacecraft Flap 1.STL", {4549.0, 37600.0, 12230.8}},
        {2, "fwd_bottom", "Spacecraft Flap 2.STL", {4549.0, 37600.0, 3223.6}},
        {3, "aft_top",    "Spacecraft Flap 3.STL", {4546.3, 4000.0, 12230.8}},
        {4, "aft_bottom", "Spacecraft Flap 4.STL", {4551.2, 4000.0, 3223.5}},
    };
    for (const FlapFile& f : flap_files) {
        const PanelMesh flap_raw = LoadMeshFromStl((dir / f.filename).string(), f.id);
        AppendTransformed(flap_raw, geo.mesh);

        PanelGroup g;
        g.id = f.id;
        g.name = f.name;
        g.hinge_point = ToModelFrame(f.hinge_cad_mm);
        g.hinge_axis = Eigen::Vector3d::UnitX();
        geo.mesh.addGroup(g);
    }

    // --- Reference quantities: PLACEHOLDER geometric proxies (true CG
    // needs mass distribution, not pure geometry), derived from the body
    // mesh's (group 0 only, flaps excluded) bounding box. ---
    Eigen::Vector3d bmin = Eigen::Vector3d::Constant(std::numeric_limits<double>::max());
    Eigen::Vector3d bmax = Eigen::Vector3d::Constant(std::numeric_limits<double>::lowest());
    for (const Panel& p : geo.mesh.panels()) {
        if (p.group_id != 0) continue;
        for (const Eigen::Vector3d& v : {p.v0, p.v1, p.v2}) {
            bmin = bmin.cwiseMin(v);
            bmax = bmax.cwiseMax(v);
        }
    }
    geo.body_length = bmax.x() - bmin.x();
    const double diam_y = bmax.y() - bmin.y();
    const double diam_z = bmax.z() - bmin.z();
    geo.body_radius = 0.25 * (diam_y + diam_z);  // average of the 2 lateral bounding extents, halved
    geo.L_ref = 2.0 * geo.body_radius;
    geo.S_ref = kPi * geo.body_radius * geo.body_radius;
    geo.moment_ref = 0.5 * (bmin + bmax);

    return geo;
}

}  // namespace aero_model
