#ifndef AERO_REGIME_DISPATCH_H
#define AERO_REGIME_DISPATCH_H
// AeroRegimeDispatch.h
//
// Dispatches aero-coefficient evaluation to the right method for the
// upstream Mach number
// Regime boundaries:
//   M <  0.8       : subsonic   -- subsonicPlaceholderAero()
//   0.8 <= M < 1.2 : transonic  -- transonicPlaceholderAero()
//   1.2 <= M < 5.0 : supersonic -- ShockExpansionAeroModel (tangent-wedge/Prandtl-Meyer)
//   M >= 5.0       : hypersonic -- NewtonianAeroModel (modified Newtonian)
//

#include <unordered_map>

#include "NewtonianAero.h"
#include "PanelGeometry.h"

namespace aero_model {

inline constexpr double kMachSubsonicUpper   = 0.8;
inline constexpr double kMachTransonicUpper  = 1.2;
inline constexpr double kMachHypersonicLower = 5.0;

enum class MachRegime { kSubsonic, kTransonic, kSupersonic, kHypersonic };

MachRegime classifyMachRegime(double mach_inf);

// Maps the 3 named control axes (radians) to TestBodyGenerator's 4 flap
// groups (1=fwd_left, 2=fwd_right, 3=aft_left, 4=aft_right). Positive
// aft_diff_rad adds to aft_left and subtracts from aft_right; forward-
// differential isn't modeled (both forward flaps share fwd_sym_rad).
std::unordered_map<int, double> mapFlapAxesToGroupDeflections(
    double fwd_sym_rad, double aft_sym_rad, double aft_diff_rad);

// Subsonic (M < 0.8) PLACEHOLDER empirical trend, not CFD-derived (Hoerner
// CD0 + Allen & Perkins crossflow CN) since panel methods are invalid for
// this separated-flow bluff body. Moments are zero-trend -- no defensible
// simple theory for them here.
AeroCoefficients subsonicPlaceholderAero(double alpha_rad, double beta_rad,
                                          double body_radius, double body_length,
                                          double S_ref, double L_ref);

// Transonic (0.8 <= M < 1.2) PLACEHOLDER -- linear fairing in Mach between
// subsonicPlaceholderAero() and ShockExpansionAeroModel at M=1.2, since no
// standalone Cp predictor exists for this regime (Von Karman similarity/
// area rule are scaling guidelines, not predictors). Does not capture real
// transonic effects (shock formation / wave-drag rise) -- real CFD needed.
AeroCoefficients transonicPlaceholderAero(const PanelMesh& mesh,
                                           const std::unordered_map<int, double>& flap_deflections_rad,
                                           double alpha_rad, double beta_rad, double mach_inf,
                                           const Eigen::Vector3d& moment_ref, double S_ref, double L_ref,
                                           double body_radius, double body_length);

// The one function GenerateAeroTable.cpp calls per grid point: classifies
// mach_inf's regime and dispatches to the appropriate model above.
AeroCoefficients evaluateAeroRegime(const PanelMesh& mesh,
                                     const std::unordered_map<int, double>& flap_deflections_rad,
                                     double alpha_rad, double beta_rad, double mach_inf,
                                     const Eigen::Vector3d& moment_ref, double S_ref, double L_ref,
                                     double body_radius, double body_length);

}  // namespace aero_model

#endif // AERO_REGIME_DISPATCH_H
