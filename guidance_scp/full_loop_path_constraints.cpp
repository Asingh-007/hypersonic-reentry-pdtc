#include "full_loop_path_constraints.h"
#include "reference_stage_a.h"
#include "TrimSolver.h"
#include "AeroRegimeDispatch.h"
#include "DescentDynamics.h"
#include <algorithm>
#include <cmath>

namespace guidance_scp {

namespace {

constexpr double kPi = 3.14159265358979323846;

// Same aero closure convention as Stage A's computeAeroClosure (trim alpha,
// neutral beta) extended with the 4 independent flap states collapsed to
// the table's 3-axis query via mapGroupDeflectionsToFlapAxes.
double nOfFlaps(double r, double v, double d1, double d2, double d3, double d4,
                 const PlanetConfig& planet_config, const SpacecraftConfig& spacecraft_config) {
    double rho = DescentDynamics::atmosphereDensity(r, planet_config);
    double a_sound = DescentDynamics::speedOfSound(r, planet_config);
    double mach = v / a_sound;
    double alpha_deg = alphaTrimDeg(mach, spacecraft_config.aero_table);

    aero_model::FlapAxes axes = aero_model::mapGroupDeflectionsToFlapAxes(d1, d2, d3, d4);
    double fwd_sym_deg = axes.fwd_sym_rad * 180.0 / kPi;
    double aft_sym_deg = axes.aft_sym_rad * 180.0 / kPi;
    double aft_diff_deg = axes.aft_diff_rad * 180.0 / kPi;

    auto aero = spacecraft_config.aero_table.interpolate(mach, alpha_deg, 0.0,
                                                           fwd_sym_deg, aft_sym_deg, aft_diff_deg);
    double qbar = 0.5 * rho * v * v;
    double L = qbar * spacecraft_config.S_ref * aero.CL;
    double D = qbar * spacecraft_config.S_ref * aero.CD;
    double force_mag = std::sqrt(L * L + D * D);
    return force_mag / (spacecraft_config.mass * planet_config.g_0);
}

} 

FullLoadFactorGradient computeFullLoadFactorGradient(double r, double v,
                                                        double d1, double d2, double d3, double d4,
                                                        const PlanetConfig& planet_config,
                                                        const SpacecraftConfig& spacecraft_config) {
    FullLoadFactorGradient g;

    PathConstraintGradients rv_grad = computePathConstraintGradients(r, v, planet_config, spacecraft_config);
    g.dn_dr = rv_grad.dn_dr;
    g.dn_dV = rv_grad.dn_dV;

    double d[4] = {d1, d2, d3, d4};
    for (int i = 0; i < 4; ++i) {
        double h = std::max(std::abs(d[i]), 0.01) * 1e-4;
        double d_plus[4] = {d1, d2, d3, d4};
        double d_minus[4] = {d1, d2, d3, d4};
        d_plus[i] += h;
        d_minus[i] -= h;
        double n_plus = nOfFlaps(r, v, d_plus[0], d_plus[1], d_plus[2], d_plus[3], planet_config, spacecraft_config);
        double n_minus = nOfFlaps(r, v, d_minus[0], d_minus[1], d_minus[2], d_minus[3], planet_config, spacecraft_config);
        g.dn_ddelta[i] = (n_plus - n_minus) / (2.0 * h);
    }

    return g;
}

}  // namespace guidance_scp
