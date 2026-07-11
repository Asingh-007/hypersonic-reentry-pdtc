// GenerateDeflectedGeometry.cpp
//
// Offline tool: reads the DOE anchor-point matrix written by aero_table_gen
// (aero/data/doe_points.csv) and, for each point, rotates each flap's
// neutral-position STL about its hinge (reusing PanelMesh::deflected(), the
// same Rodrigues rotation used for the aero-coefficient sweep) by that
// point's flap deflection angle, in the STL export's own CAD frame/units
// (mm, hinge axis = CAD_Y).
//
// Writes both a combined multi-solid assembly.stl (5 "solid...endsolid"
// blocks in one file as kept for reference/debugging) AND 5 separate
// per-part files (spacecraft_body.stl, fwd_top.stl, etc.) as the latter
// are what cfd/union_geometry.py actually consumes, since it needs each
// part as an independent mesh to Boolean-union into one watertight solid
// (the parts overlap at the flap mount points, a normal CAD-assembly
// mating practice that isn't watertight on its own

#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "PanelGeometry.h"
#include "StlMeshLoader.h"
#include "FlapHingeData.h"
#include "AeroRegimeDispatch.h"

using namespace aero_model;

constexpr double kPi = 3.14159265358979323846;

namespace {

struct DoePoint {
    int idx;
    double mach;
    double alpha_deg;
    double aft_sym_deg;
};

std::vector<DoePoint> ReadDoePoints(const std::string& csv_path) {
    std::ifstream file(csv_path);
    if (!file.is_open()) {
        throw std::runtime_error("ReadDoePoints: could not open " + csv_path +
                                  " -- run aero_table_gen first to generate it");
    }
    std::string line;
    std::getline(file, line);  // header, not validated beyond the fixed column order below

    std::vector<DoePoint> points;
    while (std::getline(file, line)) {
        if (line.find_first_not_of(" \t\r\n") == std::string::npos) continue;
        std::istringstream iss(line);
        std::string field;
        std::vector<double> vals;
        while (std::getline(iss, field, ',')) vals.push_back(std::stod(field));
        if (vals.size() != 4) {
            throw std::runtime_error("ReadDoePoints: malformed row in " + csv_path);
        }
        points.push_back({static_cast<int>(vals[0]), vals[1], vals[2], vals[3]});
    }
    return points;
}

}  // namespace

int main() {
    const std::filesystem::path repo_root = std::filesystem::path(__FILE__).parent_path().parent_path();
    const std::filesystem::path geometry_dir = repo_root / "geometry";
    const std::filesystem::path doe_csv_path = repo_root / "aero" / "data" / "doe_points.csv";
    const std::filesystem::path output_root = geometry_dir / "doe_points";

    const std::vector<DoePoint> points = ReadDoePoints(doe_csv_path.string());
    std::cout << "Read " << points.size() << " DOE points from " << doe_csv_path << std::endl;

    for (const DoePoint& pt : points) {
        // Current DOE only samples aft_sym_deg (see GenerateAeroTable.cpp) --
        // fwd_sym/aft_diff are held at 0 until forward-flap/roll CFD
        // correction is in scope.
        const auto flap_defl_rad = mapFlapAxesToGroupDeflections(
            /*fwd_sym_rad=*/0.0, pt.aft_sym_deg * kPi / 180.0, /*aft_diff_rad=*/0.0);

        char point_dir_name[32];
        std::snprintf(point_dir_name, sizeof(point_dir_name), "point_%03d", pt.idx);
        const std::filesystem::path point_dir = output_root / point_dir_name;
        std::filesystem::create_directories(point_dir);

        std::ofstream assembly_out((point_dir / "assembly.stl").string());
        if (!assembly_out.is_open()) {
            throw std::runtime_error("could not open assembly.stl for writing in " + point_dir.string());
        }

        const PanelMesh body_mesh = LoadMeshFromStl((geometry_dir / kSpacecraftBodyStlFilename).string(), 0);
        AppendMeshToStlStream(assembly_out, body_mesh.panels(), "spacecraft_body");
        WriteMeshToStl((point_dir / "spacecraft_body.stl").string(), body_mesh.panels(), "spacecraft_body");

        for (const FlapHingeInfo& f : FlapHingeTable()) {
            PanelMesh flap_mesh = LoadMeshFromStl((geometry_dir / f.stl_filename).string(), f.group_id);
            PanelGroup g;
            g.id = f.group_id;
            g.hinge_point = f.hinge_point_cad_mm;
            g.hinge_axis = Eigen::Vector3d::UnitY();  // CAD frame nose-tail direction
            flap_mesh.addGroup(g);

            const std::vector<Panel> deflected = flap_mesh.deflected(flap_defl_rad);
            AppendMeshToStlStream(assembly_out, deflected, f.name);
            WriteMeshToStl((point_dir / (std::string(f.name) + ".stl")).string(), deflected, f.name);
        }

        // Manifest: the conditions this geometry configuration corresponds
        // to, for the downstream meshing/solve script to read boundary
        // conditions from without re-deriving them.
        std::ofstream manifest((point_dir / "conditions.csv").string());
        manifest << "mach,alpha_deg,aft_sym_deg,fwd_top_deg,fwd_bottom_deg,aft_top_deg,aft_bottom_deg\n";
        manifest << pt.mach << "," << pt.alpha_deg << "," << pt.aft_sym_deg << ","
                  << (flap_defl_rad.at(1) * 180.0 / kPi) << "," << (flap_defl_rad.at(2) * 180.0 / kPi) << ","
                  << (flap_defl_rad.at(3) * 180.0 / kPi) << "," << (flap_defl_rad.at(4) * 180.0 / kPi) << "\n";
    }

    std::cout << "Wrote " << points.size() << " deflected-geometry configurations to " << output_root << std::endl;
    return 0;
}
