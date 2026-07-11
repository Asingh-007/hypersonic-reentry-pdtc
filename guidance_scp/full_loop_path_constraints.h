#ifndef GUIDANCE_SCP_FULL_LOOP_PATH_CONSTRAINTS_H
#define GUIDANCE_SCP_FULL_LOOP_PATH_CONSTRAINTS_H

#include "PlanetConfig.h"
#include "SpacecraftConfig.h"
#include <array>

// Path constraints for the full 6DOF SCP loop's Phase 1. Reuses Stage A's
// exact closed-form (r,V) gradients for Qdot/qbar/n verbatim (reference_
// stage_a.h's computePathConstraintValues/computePathConstraintGradients,
// unchanged) -> this file only adds the load factor's flap-dependence
// term (n = n(r,V,flap_state) once flaps are real states), 
// via finite difference (per the confirmed decision:
// closed forms are reused where they already exist and are validated;
// only the new partial is finite-differenced).
namespace guidance_scp {

struct FullLoadFactorGradient {
    double dn_dr = 0.0, dn_dV = 0.0;               // Stage A's exact closed form, reused verbatim
    std::array<double, 4> dn_ddelta = {0, 0, 0, 0};  // NEW, finite-differenced
};

// alpha used for the flap-dependence finite difference is the same trim
// assumption Stage A's own (r,V)-only gradients use (computeAeroClosure's
// alphaTrimDeg), consistent aero closure convention throughout, not a new
// attitude-dependent one.
FullLoadFactorGradient computeFullLoadFactorGradient(double r, double v,
                                                        double d1, double d2, double d3, double d4,
                                                        const PlanetConfig& planet_config,
                                                        const SpacecraftConfig& spacecraft_config);

}  // namespace guidance_scp

#endif  // GUIDANCE_SCP_FULL_LOOP_PATH_CONSTRAINTS_H
