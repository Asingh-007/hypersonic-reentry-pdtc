#include "FiniteDifferenceJacobian.h"

#include <algorithm>
#include <cmath>

namespace guidance_scp {

EomJacobian computeEomJacobianFd(const EomFn& f, const Eigen::VectorXd& x_ref,
                                   const Eigen::VectorXd& u_ref, double step_frac) {
    const int n_x = static_cast<int>(x_ref.size());
    const int n_u = static_cast<int>(u_ref.size());

    EomJacobian jac;
    jac.Ac = Eigen::MatrixXd::Zero(n_x, n_x);
    jac.Bc = Eigen::MatrixXd::Zero(n_x, n_u);

    for (int i = 0; i < n_x; ++i) {
        const double h = std::max(std::abs(x_ref(i)), 1.0) * step_frac;
        Eigen::VectorXd x_plus = x_ref, x_minus = x_ref;
        x_plus(i) += h;
        x_minus(i) -= h;
        Eigen::VectorXd f_plus = f(x_plus, u_ref);
        Eigen::VectorXd f_minus = f(x_minus, u_ref);
        jac.Ac.col(i) = (f_plus - f_minus) / (2.0 * h);
    }

    for (int j = 0; j < n_u; ++j) {
        const double h = std::max(std::abs(u_ref(j)), 1.0) * step_frac;
        Eigen::VectorXd u_plus = u_ref, u_minus = u_ref;
        u_plus(j) += h;
        u_minus(j) -= h;
        Eigen::VectorXd f_plus = f(x_ref, u_plus);
        Eigen::VectorXd f_minus = f(x_ref, u_minus);
        jac.Bc.col(j) = (f_plus - f_minus) / (2.0 * h);
    }

    return jac;
}

}  // namespace guidance_scp
