#include "StlMeshLoader.h"

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace aero_model {

PanelMesh LoadMeshFromStl(const std::string& path, int group_id) {
    std::ifstream file(path);
    if (!file.is_open()) {
        throw std::runtime_error("LoadMeshFromStl: could not open " + path);
    }

    std::string line;
    int line_no = 0;

    auto nextLine = [&]() -> bool {
        if (!std::getline(file, line)) return false;
        ++line_no;
        return true;
    };

    if (!nextLine()) {
        throw std::runtime_error("LoadMeshFromStl: " + path + " is empty");
    }
    {
        std::istringstream iss(line);
        std::string first_token;
        iss >> first_token;
        if (first_token != "solid") {
            throw std::runtime_error(
                "LoadMeshFromStl: " + path + " does not start with 'solid' -- "
                "binary STL is not supported by this loader");
        }
    }

    PanelMesh mesh;
    int panel_count = 0;

    while (nextLine()) {
        std::istringstream iss(line);
        std::string token;
        iss >> token;

        if (token == "endsolid") {
            break;
        }
        if (token != "facet") {
            // Blank lines / unexpected tokens between facets are tolerated
            // (some exporters add stray whitespace); only "facet" starts a
            // new triangle.
            continue;
        }

        // "facet normal nx ny nz" -- read but not trusted (Panel::recompute()
        // derives the normal from vertex winding, which is authoritative).
        std::string normal_kw;
        double nx, ny, nz;
        iss >> normal_kw >> nx >> ny >> nz;
        if (normal_kw != "normal") {
            throw std::runtime_error("LoadMeshFromStl: expected 'facet normal' at line " +
                                      std::to_string(line_no));
        }

        if (!nextLine()) {
            throw std::runtime_error("LoadMeshFromStl: unexpected end of file after 'facet' at line " +
                                      std::to_string(line_no));
        }
        {
            std::istringstream loop_iss(line);
            std::string loop_kw1, loop_kw2;
            loop_iss >> loop_kw1 >> loop_kw2;
            if (loop_kw1 != "outer" || loop_kw2 != "loop") {
                throw std::runtime_error("LoadMeshFromStl: expected 'outer loop' at line " +
                                          std::to_string(line_no));
            }
        }

        Eigen::Vector3d verts[3];
        for (int i = 0; i < 3; ++i) {
            if (!nextLine()) {
                throw std::runtime_error("LoadMeshFromStl: unexpected end of file inside 'outer loop' at line " +
                                          std::to_string(line_no));
            }
            std::istringstream v_iss(line);
            std::string v_kw;
            double x, y, z;
            v_iss >> v_kw >> x >> y >> z;
            if (v_kw != "vertex" || v_iss.fail()) {
                throw std::runtime_error("LoadMeshFromStl: expected 'vertex x y z' at line " +
                                          std::to_string(line_no));
            }
            verts[i] = Eigen::Vector3d(x, y, z);
        }

        if (!nextLine()) {
            throw std::runtime_error("LoadMeshFromStl: unexpected end of file, expected 'endloop' at line " +
                                      std::to_string(line_no));
        }
        {
            std::istringstream endloop_iss(line);
            std::string endloop_kw;
            endloop_iss >> endloop_kw;
            if (endloop_kw != "endloop") {
                throw std::runtime_error("LoadMeshFromStl: expected 'endloop' at line " +
                                          std::to_string(line_no) + " -- STL triangles must have exactly 3 vertices");
            }
        }

        if (!nextLine()) {
            throw std::runtime_error("LoadMeshFromStl: unexpected end of file, expected 'endfacet' at line " +
                                      std::to_string(line_no));
        }
        {
            std::istringstream endfacet_iss(line);
            std::string endfacet_kw;
            endfacet_iss >> endfacet_kw;
            if (endfacet_kw != "endfacet") {
                throw std::runtime_error("LoadMeshFromStl: expected 'endfacet' at line " +
                                          std::to_string(line_no));
            }
        }

        Panel panel{verts[0], verts[1], verts[2]};
        panel.group_id = group_id;
        panel.recompute();
        mesh.addPanel(panel);
        ++panel_count;
    }

    if (panel_count == 0) {
        throw std::runtime_error("LoadMeshFromStl: no panels parsed from " + path +
                                  " -- check file format");
    }

    return mesh;
}

}  // namespace aero_model
