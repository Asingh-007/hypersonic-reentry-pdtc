#include "AeroRegimeDispatch.h"

#include <algorithm>
#include <cmath>

#include "ShockExpansionAero.h"

namespace aero_model {

namespace {

// PLACEHOLDER: Hoerner CD~0.8 for a
// blunt-nosed cylindrical body at this fineness ratio (40m/9m=4.44)
constexpr double kHoernerCD0Baseline = 0.8;

constexpr double kCdc = 1.2;            // 2D circular-cylinder crossflow Cd, subcritical Re
constexpr double kEtaCrossflow = 0.75;  // Allen & Perkins end-effect factor, L/D~4.4

}

MachRegime classifyMachRegime(double mach_inf) {
    if (mach_inf < kMachSubsonicUpper) return MachRegime::kSubsonic;
    if (mach_inf < kMachTransonicUpper) return MachRegime::kTransonic;
    if (mach_inf < kMachHypersonicLower) return MachRegime::kSupersonic;
    return MachRegime::kHypersonic;
}

std::unordered_map<int, double> mapFlapAxesToGroupDeflections(
    double fwd_sym_rad, double aft_sym_rad, double aft_diff_rad) {
    return {
        {1, fwd_sym_rad},
        {2, fwd_sym_rad},
        {3, aft_sym_rad + 0.5 * aft_diff_rad},
        {4, aft_sym_rad - 0.5 * aft_diff_rad},
    };
}

FlapAxes mapGroupDeflectionsToFlapAxes(double d1, double d2, double d3, double d4) {
    // Exact algebraic inverse of mapFlapAxesToGroupDeflections above: an
    // approximation when d1 != d2 (fwd_sym only has one degree of freedom in
    // the aero table)
    return FlapAxes{(d1 + d2) / 2.0, (d3 + d4) / 2.0, d3 - d4};
}

AeroCoefficients subsonicPlaceholderAero(double alpha_rad, double /*beta_rad*/,
                                          double body_radius, double body_length,
                                          double S_ref, double L_ref) {
    // Allen & Perkins (NACA TR 1048) slender-body crossflow normal force:
    // CN = 2*sin(alpha)*cos(alpha) + eta*Cdc*sin^2(alpha)*(S_planform/S_ref)
    // Beta-dependence and alpha-dependence of the axial force are out of
    // scope for this placeholder
    const double S_planform = 2.0 * body_radius * body_length;
    const double sin_a = std::sin(alpha_rad);
    const double cos_a = std::cos(alpha_rad);
    const double CN = 2.0 * sin_a * cos_a + kEtaCrossflow * kCdc * sin_a * sin_a * (S_planform / S_ref);
    const double CA = kHoernerCD0Baseline;  // roughly Mach-independent at low subsonic (Hoerner)

    AeroCoefficients c;
    c.CN = CN;
    c.CA = CA;
    // Same alpha-only wind-axis rotation used throughout this module.
    c.CL = c.CN * cos_a - c.CA * sin_a;
    c.CD = c.CN * sin_a + c.CA * cos_a;
    c.Cl_roll = 0.0;
    c.Cm = 0.0;
    c.Cn_yaw = 0.0;
    (void)L_ref;  // not needed by this placeholder as moments are zero-trend
    return c;
}

AeroCoefficients transonicPlaceholderAero(const PanelMesh& mesh,
                                           const std::unordered_map<int, double>& flap_deflections_rad,
                                           double alpha_rad, double beta_rad, double mach_inf,
                                           const Eigen::Vector3d& moment_ref, double S_ref, double L_ref,
                                           double body_radius, double body_length) {
    const AeroCoefficients sub =
        subsonicPlaceholderAero(alpha_rad, beta_rad, body_radius, body_length, S_ref, L_ref);
    const ShockExpansionAeroModel wedge;
    const AeroCoefficients sup = wedge.evaluate(mesh, flap_deflections_rad, alpha_rad, beta_rad,
                                                  kMachTransonicUpper, moment_ref, S_ref, L_ref);

    const double t = (mach_inf - kMachSubsonicUpper) / (kMachTransonicUpper - kMachSubsonicUpper);
    const double tc = std::clamp(t, 0.0, 1.0);

    AeroCoefficients c;
    c.CL = sub.CL + tc * (sup.CL - sub.CL);
    c.CD = sub.CD + tc * (sup.CD - sub.CD);
    c.Cl_roll = sub.Cl_roll + tc * (sup.Cl_roll - sub.Cl_roll);
    c.Cm = sub.Cm + tc * (sup.Cm - sub.Cm);
    c.Cn_yaw = sub.Cn_yaw + tc * (sup.Cn_yaw - sub.Cn_yaw);
    return c;
}

AeroCoefficients evaluateAeroRegime(const PanelMesh& mesh,
                                     const std::unordered_map<int, double>& flap_deflections_rad,
                                     double alpha_rad, double beta_rad, double mach_inf,
                                     const Eigen::Vector3d& moment_ref, double S_ref, double L_ref,
                                     double body_radius, double body_length) {
    switch (classifyMachRegime(mach_inf)) {
        case MachRegime::kSubsonic:
            return subsonicPlaceholderAero(alpha_rad, beta_rad, body_radius, body_length, S_ref, L_ref);
        case MachRegime::kTransonic:
            return transonicPlaceholderAero(mesh, flap_deflections_rad, alpha_rad, beta_rad, mach_inf,
                                              moment_ref, S_ref, L_ref, body_radius, body_length);
        case MachRegime::kSupersonic: {
            // TODO: no hinge-moment model for ShockExpansionAeroModel yet
            // Ch stays at its zero-trend default here (acceptable: Phase 1
            // spends the large majority of its trajectory at M>=5).
            const ShockExpansionAeroModel wedge;
            return wedge.evaluate(mesh, flap_deflections_rad, alpha_rad, beta_rad, mach_inf,
                                    moment_ref, S_ref, L_ref);
        }
        case MachRegime::kHypersonic:
        default: {
            const NewtonianAeroModel newton;
            AeroCoefficients c = newton.evaluate(mesh, flap_deflections_rad, alpha_rad, beta_rad, mach_inf,
                                                   moment_ref, S_ref, L_ref);
            // Hinge moments are only modeled in this (Newtonian/hypersonic)
            // regime, where Phase 1 spends the large majority of its
            // trajectory. Subsonic/transonic/supersonic all leave Ch at its zero-trend default
            c.Ch = newton.evaluateHingeMoments(mesh, flap_deflections_rad, alpha_rad, beta_rad,
                                                 mach_inf, S_ref, L_ref);
            return c;
        }
    }
}

}  // namespace aero_model
