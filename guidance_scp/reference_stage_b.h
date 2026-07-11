#ifndef GUIDANCE_SCP_REFERENCE_STAGE_B_H
#define GUIDANCE_SCP_REFERENCE_STAGE_B_H

#include "ClarabelSocpSolver.h"
#include "PlanetConfig.h"
#include <vector>

// Stage B: standalone descent-reference generator (mass-added, lossless-
// convexified minimum-fuel powered descent), per this repo's SCP-guidance
// planning doc Section 8 -> the classical Acikmese & Ploen (2007) convex
// powered-descent-guidance formulation. Uses its own small local flat ENU
// frame anchored at the tower-catch point, independent of the full sim's
// spherical/J2 DescentDynamics frame
namespace guidance_scp {

struct StageBState {
    double m = 0.0;                    // kg
    double rx = 0.0, ry = 0.0, rz = 0.0;  // m, local ENU, +z up
    double vx = 0.0, vy = 0.0, vz = 0.0;  // m/s
};

struct StageBTrajectoryHistory {
    std::vector<double> t;
    std::vector<double> m, rx, ry, rz, vx, vy, vz;
    std::vector<double> Tx, Ty, Tz;       // recovered thrust = m*u, for plotting/verification
    std::vector<double> throttle_frac;    // ||T|| / Tmax_N
};

struct StageBConfig {
    int K = 30;                        // discretization node count
    double Isp_s = 330.0;              // PLACEHOLDER
    double Tmin_N = 0.0, Tmax_N = 0.0;  // PLACEHOLDER: must be set by caller, no sane universal default
    double glideslope_deg = 80.0;      // PLACEHOLDER: half-angle of the approach cone, measured from vertical
    double tf_min_s = 5.0, tf_max_s = 120.0;  // PLACEHOLDER search bracket for the outer free-final-time search
    int tf_search_iters = 15;          // bisection iterations for the outer tf search
};

struct StageBResult {
    StageBTrajectoryHistory history;
    bool solved = false;
    double tf_s = 0.0;   // the t_f the outer search converged to
};

// `initial.m` is the wet mass at the Stage-A/B handoff. `terminal`'s
// position/velocity are matched exactly at the final node (mass is left
// free, maximized, i.e. minimum propellant).

StageBResult solveStageB(const StageBState& initial, const StageBState& terminal,
                          const PlanetConfig& planet_config, const StageBConfig& config);

SocpProblem buildFixedTfProblem(const StageBState& initial, const StageBState& terminal,
                                  const PlanetConfig& planet_config, const StageBConfig& config, double tf);

}  // namespace guidance_scp

#endif  // GUIDANCE_SCP_REFERENCE_STAGE_B_H
