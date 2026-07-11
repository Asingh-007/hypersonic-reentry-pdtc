#ifndef GUIDANCE_SCP_FULL_SCP_LOOP_H
#define GUIDANCE_SCP_FULL_SCP_LOOP_H

#include "reference_stage_a.h"
#include "reference_stage_b.h"
#include "full_loop_subproblem.h"
#include "full_loop_transition.h"
#include <Eigen/Dense>
#include <vector>

// Outer SCP loop: stitches Stage A/B's
// already-converged reduced references into a full 21/22-state initial
// reference (Stage C), then repeatedly linearizes+solves+updates until
// convergence, mirroring Stage A's exact convergence-diagnostic-logging
// pattern (max_nu_per_iter/max_eta_per_iter).
namespace guidance_scp {

struct StageCResult {
    std::vector<Eigen::VectorXd> x1_ref;   // K1 nodes, size kPhase1StateDim
    std::vector<Eigen::VectorXd> u1_ref;   // K1-1 intervals, size kPhase1ControlDim
    std::vector<Eigen::VectorXd> x2_ref;   // K2 nodes, size kPhase2StateDim (raw mass)
    std::vector<Eigen::VectorXd> u2_ref;   // K2-1 intervals, size kPhase2ControlDim (physical)
    std::vector<double> deltaE_ref, phiE_ref;  // K2 nodes
};

// Lifts Stage A's/B's reduced-state histories into the full state space,
// resampled onto K1/K2 nodes (linear interpolation over each stage's own
// time history). New states (attitude, rate, flap position/rate, gimbal
// angle) are held at simple, physically reasonable defaults
StageCResult stitchStageAAndB(const StageAResult& stage_a, const StageBResult& stage_b,
                                int K1, int K2,
                                const PlanetConfig& planet_config,
                                const SpacecraftConfig& spacecraft_config);

struct FullLoopResult {
    bool converged = false;
    std::vector<double> max_nu_per_iter;   // one entry per ACCEPTED iteration
    std::vector<double> max_eta_per_iter;
    std::vector<Eigen::VectorXd> x1_ref, x2_ref;  // final trajectory (converged or last attempted)
    std::vector<Eigen::VectorXd> u1_ref, u2_ref;
    std::vector<double> deltaE_ref, phiE_ref;

    // Adaptive trust-region diagnostics
    std::vector<double> attempt_max_defect_nonlinear;
    std::vector<double> attempt_eta_max;
    std::vector<bool> attempt_accepted;
};

// Runs the inner SCP relinearization loop (Milestone 10.2) starting from
// `initial_reference`, with t_scale_1/t_scale_2 held FIXED throughout (a
// scoped-down simplification vs. the plan's "outer-outer" dilation-factor
// search flagged explicitly, not silently dropped: revisit if empirical
// runs show the fixed dilation factors are inadequate).
FullLoopResult solveFullLoop(const StageCResult& initial_reference,
                              double t_scale_1, double t_scale_2,
                              const Eigen::VectorXd& x1_initial,
                              const Eigen::VectorXd& x2_terminal_target,
                              const Phase1ToPhase2Frame& transition_frame,
                              const PlanetConfig& planet_config,
                              const SpacecraftConfig& spacecraft_config,
                              const FlapActuatorConfig& flap_config,
                              const GimbalActuatorConfig& gimbal_config,
                              const FullLoopConfig& full_loop_config,
                              int max_iters = 20, double eps_dyn = 1.0, double eps_tr = 1.0);

}  // namespace guidance_scp

#endif  // GUIDANCE_SCP_FULL_SCP_LOOP_H
