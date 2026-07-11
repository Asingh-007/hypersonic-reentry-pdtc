#ifndef NEWTONIAN_AERO_H
#define NEWTONIAN_AERO_H
// NewtonianAero.h
//
// Modified Newtonian panel-method aerodynamics (Lees-modified Newtonian
// impact theory), following the same Cp = Cp_max * sin^2(theta) formulation
// used for the Apollo 10 CM analysis in Olsen & Bettinger (2022). Extended
// here to support flap-deflection sensitivity (Cm vs flap) via PanelMesh's
// deflectable panel groups
//
// This model is intentionally cheap: it is meant to be swept densely over
// (Mach, alpha, beta, flap) OFFLINE by aero/GenerateAeroTable.cpp, and
// corrected against a sparse set of CFD anchor points via UniversalKriging
// once real CFD data exists. DescentDynamics never calls
// this model directly at runtime -- it only queries the resulting
// AeroCoefficientTable.

#include <array>
#include <unordered_map>

#include "PanelGeometry.h"

namespace aero_model {

struct AeroCoefficients {
    double CL = 0.0;       // lift coefficient, wind axes
    double CD = 0.0;       // drag coefficient, wind axes
    double CN = 0.0;       // normal force coefficient, body axes
    double CA = 0.0;       // axial force coefficient, body axes
    double Cl_roll = 0.0;  // rolling moment coefficient about moment_ref, body x-axis
    double Cm = 0.0;       // pitching moment coefficient about moment_ref, body y-axis
    double Cn_yaw = 0.0;   // yawing moment coefficient about moment_ref, body z-axis
    // Per-flap hinge-moment coefficients, Ch[i] for group_id (i+1), normalized
    // by S_ref*L_ref (same convention as Cm) -- moment about that flap's own
    // hinge point/axis, NOT moment_ref. Zero in regimes with no defensible
    // hinge-moment model (subsonic/transonic placeholders, supersonic).
    std::array<double, 4> Ch = {0.0, 0.0, 0.0, 0.0};
};

// Body-axis convention assumed throughout: x points from tail to nose,
// y is the pitch axis, z completes a right-handed triad. Angle of attack
// alpha and sideslip beta are the standard aircraft-body-axes definitions
// (see freestreamDirectionBody). If this does not match vehicle_dynamics'
// body-frame convention, adjust the sign in freestreamDirectionBody()
// accordingly as that is the single place the convention is encoded.
class NewtonianAeroModel {
public:
    explicit NewtonianAeroModel(double gamma = 1.4) : gamma_(gamma) {}

    // Cp behind a normal shock at the stagnation point (Anderson's
    // modified-Newtonian Cp_max), as a function of freestream Mach number.
    double cpMax(double mach_inf) const;

    // Evaluates lift/drag/normal/axial/moment coefficients at a single
    // (alpha, beta, Mach, flap deflection) condition.
    //   moment_ref  : point about which moments are taken, body frame (e.g. CG)
    //   S_ref, L_ref: reference area / length for non-dimensionalization
    //
    AeroCoefficients evaluate(const PanelMesh& mesh,
                               const std::unordered_map<int, double>& flap_deflections_rad,
                               double alpha_rad, double beta_rad, double mach_inf,
                               const Eigen::Vector3d& moment_ref, double S_ref,
                               double L_ref) const;

    // Computes all 4 flaps' hinge-moment coefficients in one pass over the
    // deflected panel set (reuses the same Cp/force-per-panel work evaluate()
    // does, as a separate pass rather than merged into it -- this is a
    // one-time offline sweep, so simplicity beats micro-optimization here).
    // Moment for flap group_id=i+1 is taken about mesh.groups().at(i+1)'s own
    // hinge_point (NOT moment_ref), projected onto that group's hinge_axis,
    // normalized by S_ref*L_ref for dimensional consistency with Cm. Panels
    // with group_id not in [1,4] (the fixed body) do not contribute.
    std::array<double, 4> evaluateHingeMoments(
        const PanelMesh& mesh, const std::unordered_map<int, double>& flap_deflections_rad,
        double alpha_rad, double beta_rad, double mach_inf, double S_ref, double L_ref) const;

private:
    Eigen::Vector3d freestreamDirectionBody(double alpha_rad, double beta_rad) const;

    double gamma_;
};

}  // namespace aero_model

#endif // NEWTONIAN_AERO_H
