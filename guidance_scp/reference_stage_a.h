#ifndef GUIDANCE_SCP_REFERENCE_STAGE_A_H
#define GUIDANCE_SCP_REFERENCE_STAGE_A_H

#include "PlanetConfig.h"
#include "SpacecraftConfig.h"
#include <vector>

// Stage A: standalone entry reference-trajectory generator (3DOF translation
// + bank angle, engines off), per this repo's SCP-guidance planning doc
// Section 8. Deliberately not the full 13-state DescentDynamics -> this is
// its own small, cheap, sequential-SOCP-solved reference generator; the
// full 6DOF coupled SCP loop (flap actuators, mass depletion, gimbal rates)
// is out of scope here
namespace guidance_scp {

// sigma here is ALWAYS bank angle (rad) -> never the classical SCP-
// literature "free time-dilation factor" (that quantity is called t_scale
// throughout this module instead, to avoid the name collision).
struct StageAState {
    double r = 0.0, la = 0.0, lo = 0.0, v = 0.0, fpa = 0.0, v_azi = 0.0, sigma = 0.0;
};

struct StageATrajectoryHistory {
    std::vector<double> t;
    std::vector<double> r, la, lo, v, fpa, v_azi, sigma;
    std::vector<double> sigma_dot;                              // control history (K-1 entries)
    std::vector<double> qbar, heat_flux_conv, load_factor_g;    // for constraint verification/plots
    std::vector<double> alpha_trim_deg, mach;                    // recorded for a future 6DOF lift-up
};

struct StageAConfig {
    int K = 40;                          // discretization node count
    double sigma_max = 1.2;              // rad (~68 deg), PLACEHOLDER bank-angle limit
    double sigma_dot_max = 0.1;          // rad/s, PLACEHOLDER bank-rate limit
    double w_nu = 1.0e4;                 // virtual-control L1 penalty weight, PLACEHOLDER
    double w_trust = 1.0;                // trust-region penalty weight, PLACEHOLDER
    int max_iters = 20;
    double eps_dyn = 1e-3;               // convergence: max virtual-control L1 norm per node
    double eps_tr = 1e-2;                // convergence: max trust-region radius per node
    // t_scale (phase duration) is treated as a FIXED, user-supplied value for
    // this reference generator, not solved as a free dilation variable
    double t_scale_s = 60.0;             // PLACEHOLDER phase-duration guess, s
    double max_heat_flux_w_m2 = 1.0e6;   // PLACEHOLDER, matches GuidanceConstraints' default
    double max_qbar_pa = 50000.0;        // PLACEHOLDER, matches GuidanceConstraints' default
    double max_load_factor_g = 5.0;      // PLACEHOLDER, matches GuidanceConstraints' default

    // --- Adaptive trust-region accept/reject
    double eta_max_init = 50.0;
    double eta_max_floor = 1.0e-3;
    double eta_max_ceiling = 1.0e4;
    double eta_max_shrink = 0.5;
    double eta_max_grow = 2.0;
    // Accept a solved step iff it is the first-ever accepted step (nothing
    // to compare against yet) OR max||nu|| does not regress by more than
    // this factor relative to the best max||nu|| seen so far -> e.g. 1.5
    // allows up to 50% worse than the best-so-far before rejecting.
    double nu_regression_tolerance = 1.5;
    int max_solve_attempts = 60;              // accepted+rejected attempts, bounds total work
    int max_consecutive_floor_rejects = 5;    // abort honestly if stuck at eta_max_floor this many times in a row
};

struct StageAResult {
    StageATrajectoryHistory history;
    bool converged = false;
    std::vector<double> max_nu_per_iter;   // convergence diagnostic, one entry per ACCEPTED iteration
    std::vector<double> max_eta_per_iter;  // convergence diagnostic, one entry per ACCEPTED iteration

    // Adaptive trust-region diagnostics, one entry per solve attempt
    // (accepted or rejected), distinct from the two vectors above, which
    // are accepted-iterations-only (unchanged semantics/length convention).
    std::vector<double> attempt_max_defect_nonlinear;
    std::vector<double> attempt_eta_max;
    std::vector<bool> attempt_accepted;
};

// terminal is treated as a full target state (not a partial/relaxed one)
// every field of `terminal` is matched exactly at the final node. A future
// caller wanting a genuinely partial terminal constraint (e.g. free lo,
// v_azi) would need the boundary-equality block below relaxed per-field.
StageAResult solveStageA(const StageAState& initial, const StageAState& terminal,
                          const PlanetConfig& planet_config, const SpacecraftConfig& spacecraft_config,
                          const StageAConfig& config);

//Exposed for direct unit testing (finite-difference gradient checks)

struct PathConstraintValues {
    double qdot_w_m2 = 0.0;   // Sutton-Graves stagnation-point heat flux (loads/HeatingLoadModel.h, reused exactly)
    double qbar_pa = 0.0;
    double n_g = 0.0;         // load factor, g's (sqrt(L^2+D^2)/(mass*g_0))
};

struct PathConstraintGradients {
    double dQdot_dr = 0.0, dQdot_dV = 0.0;
    double dqbar_dr = 0.0, dqbar_dV = 0.0;
    double dn_dr = 0.0, dn_dV = 0.0;
};

// Aero closure at a given (r, V): density/Mach via DescentDynamics's public
// static atmosphereDensity/speedOfSound, alpha at trim via TrimSolver, lift/
// drag via AeroCoefficientTable::interpolate() at that trim alpha, neutral
// flaps/beta, exactly the closure Stage A's translational EOM and path
// constraints both need.
struct AeroClosure {
    double rho = 0.0, mach = 0.0, alpha_trim_deg = 0.0, qbar_pa = 0.0, lift_n = 0.0, drag_n = 0.0;
};
AeroClosure computeAeroClosure(double r, double v, const PlanetConfig& planet_config,
                                const SpacecraftConfig& spacecraft_config);

PathConstraintValues computePathConstraintValues(double r, double v, const PlanetConfig& planet_config,
                                                    const SpacecraftConfig& spacecraft_config);

// Assumes EarthAtmosphere1976's low-altitude exponential regime (h < 84km,
// see AtmosphereModel.h) -- valid for essentially the entire altitude band
// where these constraints actually bind for a hypersonic entry. See this
// file's .cpp for the closed-form derivation (all three gradients reduce to
// "-(some multiple of kBeta) * the quantity itself").
PathConstraintGradients computePathConstraintGradients(double r, double v, const PlanetConfig& planet_config,
                                                          const SpacecraftConfig& spacecraft_config);

}  // namespace guidance_scp

#endif  // GUIDANCE_SCP_REFERENCE_STAGE_A_H
