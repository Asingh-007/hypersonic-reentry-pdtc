#ifndef GUIDANCE_SCP_FINITE_DIFFERENCE_JACOBIAN_H
#define GUIDANCE_SCP_FINITE_DIFFERENCE_JACOBIAN_H

#include <Eigen/Dense>
#include <functional>

// Generalized, dynamic-sized central-finite-difference Jacobian of an
// arbitrary nonlinear vector function f(x,u) -> xdot, generalizing Stage A's
// fixed-7x1 computeEomJacobian() (reference_stage_a.cpp) to arbitrary
// state/control dimension, used by both full-loop phases (up to Phase 2's
// 22-state/7-control) and the phase-transition frame-conversion map
// (Milestone 8, with a zero-length control). Confirmed decision (this
// session): finite differences here, not hand-derived analytical Jacobians
// automatically consistent with whatever the EOM function actually computes.
namespace guidance_scp {

using EomFn = std::function<Eigen::VectorXd(const Eigen::VectorXd& x, const Eigen::VectorXd& u)>;

struct EomJacobian {
    Eigen::MatrixXd Ac;  // n_x x n_x
    Eigen::MatrixXd Bc;  // n_x x n_u
};

// Per-component step h_i = max(|x_i or u_i|, 1.0) * step_frac -> identical
// heuristic to Stage A's computeEomJacobian, already proven numerically
// sound there.
EomJacobian computeEomJacobianFd(const EomFn& f, const Eigen::VectorXd& x_ref,
                                   const Eigen::VectorXd& u_ref, double step_frac = 1e-6);

}  // namespace guidance_scp

#endif  // GUIDANCE_SCP_FINITE_DIFFERENCE_JACOBIAN_H
