#include "PanelGeometry.h"
#include <cmath>

namespace aero_model {

namespace {

// Rotates point p about an axis (unit vector) passing through hinge_point,
// by angle_rad radians, using Rodrigues' rotation formula.
Eigen::Vector3d rotateAboutHinge(const Eigen::Vector3d& p,
                                  const Eigen::Vector3d& hinge_point,
                                  const Eigen::Vector3d& hinge_axis,
                                  double angle_rad) {
    const Eigen::Vector3d r = p - hinge_point;
    const double c = std::cos(angle_rad);
    const double s = std::sin(angle_rad);
    const Eigen::Vector3d rotated =
        r * c + hinge_axis.cross(r) * s +
        hinge_axis * (hinge_axis.dot(r)) * (1.0 - c);
    return hinge_point + rotated;
}

}  // namespace

std::vector<Panel> PanelMesh::deflected(
    const std::unordered_map<int, double>& deflections) const {
    std::vector<Panel> out;
    out.reserve(panels_.size());

    for (const Panel& p : panels_) {
        if (p.group_id == 0) {
            out.push_back(p);
            continue;
        }

        auto def_it = deflections.find(p.group_id);
        auto grp_it = groups_.find(p.group_id);
        if (def_it == deflections.end() || grp_it == groups_.end()) {
            // No deflection commanded / group metadata missing: pass through.
            out.push_back(p);
            continue;
        }

        const double angle = def_it->second;
        const PanelGroup& g = grp_it->second;

        Panel rotated = p;
        rotated.v0 = rotateAboutHinge(p.v0, g.hinge_point, g.hinge_axis, angle);
        rotated.v1 = rotateAboutHinge(p.v1, g.hinge_point, g.hinge_axis, angle);
        rotated.v2 = rotateAboutHinge(p.v2, g.hinge_point, g.hinge_axis, angle);
        rotated.recompute();
        out.push_back(rotated);
    }
    return out;
}

}  // namespace aero_model
