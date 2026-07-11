#ifndef CLARABEL_SOCP_SOLVER_H
#define CLARABEL_SOCP_SOLVER_H

#include <clarabel.hpp>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>

// Thin, solver-agnostic wrapper around Clarabel's C++ SOCP interface (already
// vendored/linked in this repo
// Both guidance_scp::solveStageA (per SCP iteration) and
// solveStageB (per free-final-time candidate) go through this.
//
// Problem form (matches Clarabel's own convention exactly): minimize
// 0.5*x'Px + q'x subject to Ax + s = b, s in the direct product of `cones`.
namespace guidance_scp {

struct SocpProblem {
    Eigen::SparseMatrix<double> P; // pass an all-zero (nnz=0) matrix if the cost has no quadratic term
    Eigen::VectorXd q;
    Eigen::SparseMatrix<double> A;
    Eigen::VectorXd b;
    std::vector<clarabel::SupportedConeT<double>> cones;
};

struct SocpSolution {
    Eigen::VectorXd x;
    clarabel::SolverStatus status;
    double obj_val = 0.0;
};

// Builds a clarabel::DefaultSolver from `problem`, solves it, and copies the
// solution out. Clarabel's DefaultSolution::x is an Eigen::Map into the
// solver's internal buffer
SocpSolution solve(const SocpProblem& problem);

}  // namespace guidance_scp

#endif  // CLARABEL_SOCP_SOLVER_H
