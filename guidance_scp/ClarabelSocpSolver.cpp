#include "ClarabelSocpSolver.h"

namespace guidance_scp {

SocpSolution solve(const SocpProblem& problem) {
    clarabel::DefaultSettings<double> settings = clarabel::DefaultSettings<double>::default_settings();

    Eigen::VectorXd q = problem.q;
    Eigen::VectorXd b = problem.b;
    clarabel::DefaultSolver<double> solver(problem.P, q, problem.A, b, problem.cones, settings);
    solver.solve();

    clarabel::DefaultSolution<double> solution = solver.solution();

    SocpSolution result;
    result.x = solution.x;  // copies out of the solver-owned Eigen::Map
    result.status = solution.status;
    result.obj_val = solution.obj_val;
    return result;
}

}  // namespace guidance_scp
