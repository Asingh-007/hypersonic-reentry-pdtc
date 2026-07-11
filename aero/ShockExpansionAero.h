#ifndef SHOCK_EXPANSION_AERO_H
#define SHOCK_EXPANSION_AERO_H
// ShockExpansionAero.h
//
// Tangent-Wedge (oblique shock) + Prandtl-Meyer expansion panel method, the
// supersonic (1.2 <= M < ~5) complement to NewtonianAeroModel: each panel's
// local flow-deflection angle is treated as an isolated 2D wedge/expansion
// corner (Anderson's "shock-expansion method", generalized panel-locally
// like NASA's S/HABP code does).
// Offline-only, like NewtonianAeroModel, swept by aero/GenerateAeroTable.cpp
// via AeroRegimeDispatch.h; DescentDynamics never calls this directly.

#include <optional>
#include <unordered_map>

#include "NewtonianAero.h" 
#include "PanelGeometry.h"

namespace aero_model {

// Prandtl-Meyer function nu(M) (radians), valid for M >= 1 (nu(1) = 0).
// Anderson, "Fundamentals of Aerodynamics", eq. 9.42.
double prandtlMeyerNu(double mach, double gamma);

// Maximum Prandtl-Meyer turning angle (M -> infinity, ~130.45deg for
// gamma=1.4); used to clamp the expansion branch defensively against a
// NaN/garbage inverse-PM solve.
double prandtlMeyerNuMax(double gamma);

// Inverts nu(M) = nu_target for M (no closed form) via Newton iteration.
// Returns M >= 1.
double invertPrandtlMeyer(double nu_target, double gamma);

// Solves the theta-beta-Mach relation (Anderson eq. 4.23) for the WEAK
// oblique shock angle beta (radians), seeded near the Mach angle to bias
// toward the physically-realized weak-shock root. Returns std::nullopt if
// the shock is detached as callers should fall back to another model.
std::optional<double> solveWeakObliqueShockBeta(double mach_inf, double theta_rad, double gamma);

// Body-axis convention and evaluate() semantics match NewtonianAeroModel;
// selected instead of it in the supersonic regime by AeroRegimeDispatch.h.
class ShockExpansionAeroModel {
public:
    explicit ShockExpansionAeroModel(double gamma = 1.4) : gamma_(gamma) {}

    AeroCoefficients evaluate(const PanelMesh& mesh,
                               const std::unordered_map<int, double>& flap_deflections_rad,
                               double alpha_rad, double beta_rad, double mach_inf,
                               const Eigen::Vector3d& moment_ref, double S_ref,
                               double L_ref) const;

private:
    Eigen::Vector3d freestreamDirectionBody(double alpha_rad, double beta_rad) const;

    // Cp for a single panel given its signed local flow-deflection angle:
    // theta >= 0 (compressive, flow turns INTO the panel) uses the
    // oblique-shock relation, falling back to local Newtonian
    // cpMax(mach)*sin^2(theta) on shock detachment; theta < 0 (expansive,
    // flow turns away) uses Prandtl-Meyer expansion. Unlike
    // NewtonianAeroModel, every panel contributes (nothing is zeroed out
    // as "leeward") since both branches produce a physically meaningful Cp.
    double panelCp(double mach_inf, double theta_rad) const;

    double gamma_;
};

}  // namespace aero_model

#endif // SHOCK_EXPANSION_AERO_H
