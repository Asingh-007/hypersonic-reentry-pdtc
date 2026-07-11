#ifndef STL_MESH_LOADER_H
#define STL_MESH_LOADER_H

#include "PanelGeometry.h"
#include <fstream>
#include <string>
#include <vector>

namespace aero_model {

// Loads an ASCII STL file (as commonly exported by OpenVSP, OpenRocket, and
// CAD tools such as SolidWorks) into a PanelMesh. Every triangle in the file
// becomes one Panel with the given group_id (0 = fixed body by default;
// pass a nonzero id and register a matching PanelGroup separately if the
// imported geometry represents a deflectable surface -> STL has no concept
// of hinge groups, so that metadata must be added by hand after loading).
//
// Malformed files raise std::runtime_error with a line-number-anchored
// message rather than silently producing a garbage or empty mesh.
PanelMesh LoadMeshFromStl(const std::string& path, int group_id = 0);

// Writes `panels` out as an ASCII STL file (format symmetric with
// LoadMeshFromStl -- one "facet normal/outer loop/3 vertex/endloop/
// endfacet" block per panel), using each panel's already-computed
// `normal` field (call recompute() first if vertices were modified).
// Throws std::runtime_error if the file can't be opened for writing.
void WriteMeshToStl(const std::string& path, const std::vector<Panel>& panels,
                     const std::string& solid_name);

// Appends one "solid <solid_name> ... endsolid <solid_name>" block to an
// already-open output stream, without opening/closing a file itself --
// call this once per part on the same stream to build a single multi-
// solid ASCII STL file (e.g. a full vehicle assembly as one file, which
// downstream CAD/meshing tools generally import more reliably than a
// list of separate single-part files -- see GenerateDeflectedGeometry.cpp).
void AppendMeshToStlStream(std::ofstream& out, const std::vector<Panel>& panels,
                            const std::string& solid_name);

}  // namespace aero_model

#endif // STL_MESH_LOADER_H
