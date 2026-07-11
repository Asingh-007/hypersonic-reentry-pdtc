#include "SpacecraftGeometry.h"
#include "StlMeshLoader.h"
#include "FlapHingeData.h"

#include <filesystem>
#include <limits>

namespace aero_model {

namespace {

constexpr double kPi = 3.14159265358979323846;

constexpr double kMmToM = 1.0 / 1000.0;  // STL export uses mm; this codebase uses meters throughout.

// model_x=CAD_Y (nose-tail, +model_x=nose), model_y=-CAD_X, model_z=CAD_Z ->
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

}

SpacecraftGeometry LoadSpacecraftGeometry(const std::string& geometry_dir) {
    namespace fs = std::filesystem;
    const fs::path dir(geometry_dir);

    SpacecraftGeometry geo;

    // Body (group 0, fixed)
    const PanelMesh body_raw = LoadMeshFromStl((dir / kSpacecraftBodyStlFilename).string(), 0);
    AppendTransformed(body_raw, geo.mesh);

    // 4 flaps, Hinge axis is the body's nose-tail direction (model +X = CAD_Y) for all 4
    for (const FlapHingeInfo& f : FlapHingeTable()) {
        const PanelMesh flap_raw = LoadMeshFromStl((dir / f.stl_filename).string(), f.group_id);
        AppendTransformed(flap_raw, geo.mesh);

        PanelGroup g;
        g.id = f.group_id;
        g.name = f.name;
        g.hinge_point = ToModelFrame(f.hinge_point_cad_mm);
        g.hinge_axis = Eigen::Vector3d::UnitX();
        geo.mesh.addGroup(g);
    }

    // Reference quantities: PLACEHOLDER geometric proxies (true CG
    // needs mass distribution, not pure geometry), derived from the body
    // mesh's (group 0 only, flaps excluded) bounding box.
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
