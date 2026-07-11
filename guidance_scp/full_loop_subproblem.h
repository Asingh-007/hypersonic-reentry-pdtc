#ifndef GUIDANCE_SCP_FULL_LOOP_SUBPROBLEM_H
#define GUIDANCE_SCP_FULL_LOOP_SUBPROBLEM_H

#include "ClarabelSocpSolver.h"
#include "PlanetConfig.h"
#include "SpacecraftConfig.h"
#include "FlapActuatorConfig.h"
#include "GimbalActuatorConfig.h"
#include "full_loop_transition.h"
#include "FiniteDifferenceJacobian.h"
#include <Eigen/Dense>
#include <array>
#include <vector>

// Full 2-phase coupled SCP subproblem assembly
// Combines Phase 1 (guidance_scp/full_loop_phase1.h) and
// Phase 2 (guidance_scp/full_loop_phase2.h) dynamics/constraints into one
// joint convex subproblem per outer SCP iteration, following Stage A/B's
// exact conventions (VarLayout offset-helper pattern, addIneqRow/addEqRow
// coeffs.x<=rhs sign convention, per-component trust-region/virtual-control
// scaling, mandatory solver-status checking).
namespace guidance_scp {

struct FullLoopConfig {
    double Isp_s = 330.0;
    double Tmin_N = 0.0, Tmax_N = 0.0;         // must be set by caller
    double glideslope_deg = 80.0;
    double qbar_relight_max_pa = 5000.0;        // relight condition at the Phase1->2 transition
    double terminal_attitude_error_max_rad = 0.2;  // linearized quaternion-error bound at tower catch
    double w_nu = 1.0e4, w_eta = 1.0;           // SCP penalty weights, PLACEHOLDER
    double m_wet_at_handoff_kg = 0.0;           // must be set by caller (Phase 1 has no mass state)

    // Adaptive trust-region accept/reject, same design as StageAConfig
    double eta_max_init = 50.0;
    double eta_max_floor = 1.0e-3;
    double eta_max_ceiling = 1.0e4;
    double eta_max_shrink = 0.5;
    double eta_max_grow = 2.0;
    // Accept a solved step iff it is the first-ever accepted step OR
    // max||nu|| (combined across both phases) does not regress by more than
    // this factor relative to the best value seen so far.
    double nu_regression_tolerance = 1.5;
    int max_solve_attempts = 80;              // accepted+rejected attempts, bounds total work
    int max_consecutive_floor_rejects = 5;    // abort honestly if stuck at eta_max_floor this many times in a row
};

// All per-outer-iteration inputs: both phases' current reference
// trajectories/controls (deviation-from-reference decision variables are
// built relative to these), the per-phase free dilation factors (held
// fixed within this one convex solve, per the confirmed simplification),
// boundary conditions, and physical/actuator configs.
struct FullLoopSubproblemInputs {
    std::vector<Eigen::VectorXd> x1_ref;   // K1 nodes, each size kPhase1StateDim
    std::vector<Eigen::VectorXd> u1_ref;   // K1-1 intervals, each size kPhase1ControlDim
    double t_scale_1 = 0.0;

    std::vector<Eigen::VectorXd> x2_ref;   // K2 nodes, each size kPhase2StateDim (mass slot RAW m)
    std::vector<Eigen::VectorXd> u2_ref;   // K2-1 intervals, each size kPhase2ControlDim (physical controls)
    std::vector<double> deltaE_ref, phiE_ref;  // K2 nodes, gimbal-angle bookkeeping
    double t_scale_2 = 0.0;

    Eigen::VectorXd x1_initial;            // size 21, fixed Phase-1 boundary
    Eigen::VectorXd x2_terminal_target;    // size 22 (only position/velocity/attitude enforced)

    Phase1ToPhase2Frame transition_frame;

    PlanetConfig planet_config = PlanetConfig::Earth();
    SpacecraftConfig* spacecraft_config = nullptr;  // non-owning, must outlive the call
    FlapActuatorConfig flap_config;
    GimbalActuatorConfig gimbal_config;
    FullLoopConfig full_loop_config;

    double eta_max = 1.0;
};

// Fixed variable-layout offset helper, extending Stage A's VarLayout pattern
// to 2 phases concatenated (phase-major): Phase 1 block
// (delta_x1/delta_u1/nu1_plus/nu1_minus/eta1), then Phase 2 block
// (delta_x2 [mass slot repurposed as delta_z]/delta_u2 [u_vec(3)+sigma_u(1)
// +rest(6): deltaE_dot,phiE_dot,tau_m1..4]/nu2_plus/nu2_minus/eta2).
struct FullLoopVarLayout {
    int K1 = 0, K2 = 0;
    static constexpr int kNx1 = 21, kNu1 = 4;
    static constexpr int kNx2 = 22;
    static constexpr int kNu2Vec = 3, kNu2Sigma = 1, kNu2Rest = 6;
    static constexpr int kNu2 = kNu2Vec + kNu2Sigma + kNu2Rest;  // 10

    int dx1_offset(int k) const { return kNx1 * k; }
    int p1_dx_total() const { return kNx1 * K1; }
    int du1_offset(int k) const { return p1_dx_total() + kNu1 * k; }
    int p1_du_total() const { return kNu1 * (K1 - 1); }
    int nu1_plus_offset(int k) const { return p1_dx_total() + p1_du_total() + kNx1 * k; }
    int nu1_minus_offset(int k) const {
        return p1_dx_total() + p1_du_total() + kNx1 * (K1 - 1) + kNx1 * k;
    }
    int p1_nu_total() const { return 2 * kNx1 * (K1 - 1); }
    int eta1_offset(int k) const { return p1_dx_total() + p1_du_total() + p1_nu_total() + k; }
    int p1_eta_total() const { return K1; }
    int p1_total() const { return p1_dx_total() + p1_du_total() + p1_nu_total() + p1_eta_total(); }

    int dx2_offset(int k) const { return p1_total() + kNx2 * k; }
    int p2_dx_total() const { return kNx2 * K2; }
    int du2_uvec_offset(int k) const { return p1_total() + p2_dx_total() + kNu2 * k; }
    int du2_sigma_offset(int k) const { return du2_uvec_offset(k) + kNu2Vec; }
    int du2_rest_offset(int k) const { return du2_uvec_offset(k) + kNu2Vec + kNu2Sigma; }
    int p2_du_total() const { return kNu2 * (K2 - 1); }
    int nu2_plus_offset(int k) const {
        return p1_total() + p2_dx_total() + p2_du_total() + kNx2 * k;
    }
    int nu2_minus_offset(int k) const {
        return p1_total() + p2_dx_total() + p2_du_total() + kNx2 * (K2 - 1) + kNx2 * k;
    }
    int p2_nu_total() const { return 2 * kNx2 * (K2 - 1); }
    int eta2_offset(int k) const {
        return p1_total() + p2_dx_total() + p2_du_total() + p2_nu_total() + k;
    }
    int p2_eta_total() const { return K2; }

    int total() const {
        return p1_total() + p2_dx_total() + p2_du_total() + p2_nu_total() + p2_eta_total();
    }
};

SocpProblem buildFullLoopSubproblem(const FullLoopSubproblemInputs& in);

// Per-phase characteristic scales for the trust region / virtual-control
// cost (same reasoning as Stage A's kStateScale -> see
// reference_stage_a.cpp). Exposed here (not file-local to
// full_loop_subproblem.cpp) so full_scp_loop.cpp's outer loop can reuse the
// SAME normalization when computing the diagnostic nonlinear defect below --
// a single source of truth, avoiding drift between the two.
inline constexpr std::array<double, 21> kStateScale1 = {
    1000.0, 0.01, 0.01, 100.0, 0.01, 0.01,      // r,la,lo,V,fpa,v_azi (Stage A's exact numbers)
    0.1, 0.1, 0.1, 0.1,                          // quaternion
    0.05, 0.05, 0.05,                            // body rates
    0.05, 0.05, 0.05, 0.05,                      // flap positions
    0.05, 0.05, 0.05, 0.05,                      // flap rates
};
inline constexpr std::array<double, 4> kControlScale1 = {50.0, 50.0, 50.0, 50.0};  // tau_m

inline constexpr std::array<double, 22> kStateScale2 = {
    0.05,                                          // delta_z (ln-mass), NOTE: the diagnostic
                                                     // nonlinear defect below evaluates phase2Eom's
                                                     // RAW mdot (kg/s), not the ln-mass rate, so this
                                                     // component's normalized defect value is only a
                                                     // rough heuristic, not unit-consistent. Harmless:
                                                     // it is diagnostic-only, never gates accept/reject.
    100.0, 100.0, 100.0, 10.0, 10.0, 10.0,         // position/velocity
    0.1, 0.1, 0.1, 0.1,                            // quaternion
    0.05, 0.05, 0.05,                              // body rates
    0.05, 0.05, 0.05, 0.05,                        // flap positions
    0.05, 0.05, 0.05, 0.05,                        // flap rates
};
inline constexpr double kControlScaleVec2 = 5.0;      // u_vec (specific thrust), m/s^2
inline constexpr double kControlScaleSigma2 = 5.0;    // sigma_u, m/s^2
inline constexpr std::array<double, 6> kControlScaleRest2 = {0.01, 0.01, 50.0, 50.0, 50.0, 50.0};

// Diagnostic true-nonlinear one-step defect (Stage A's exact
// computeMaxNonlinearDefectStageA, generalized to arbitrary state dim and an
// arbitrary EomFn). NOT used to gate accept/reject (see
// feedback_scp_trust_region_debugging memory / reference_stage_a.cpp's
// top-of-file comment for why an absolute/eta_max-scaled defect threshold is
// the wrong criterion for this repo's SCP design) -> purely for visibility
// via FullLoopResult::attempt_max_defect_nonlinear. Phase 1 can pass its
// whole trajectory in one call (a single EomFn covers every interval).
// Phase 2 cannot: phase2Eom takes deltaE_ref[k]/phiE_ref[k] as extra
// PER-INTERVAL scalar parameters, so the caller must invoke this once per
// interval with a freshly-rebound EomFn and a 2-node "mini-trajectory"
// {x_ref_candidate[k], x_ref_candidate[k+1]} instead of one whole-trajectory call.
double computeMaxNonlinearDefectGeneric(const EomFn& f, const std::vector<Eigen::VectorXd>& x_ref_candidate,
                                          const std::vector<Eigen::VectorXd>& u_ref_candidate, double dt,
                                          const double* state_scale, int state_dim);

}  // namespace guidance_scp

#endif  // GUIDANCE_SCP_FULL_LOOP_SUBPROBLEM_H
