#include "Kriging.h"

#include <cmath>
#include <stdexcept>

namespace aero_model {

double UniversalKriging::kernel(const Eigen::VectorXd& xi,
                                 const Eigen::VectorXd& xj) const {
    // Anisotropic Matern 5/2:
    //   r^2   = sum_d ((xi_d - xj_d) / theta_d)^2
    //   k(r)  = sigma_f^2 * (1 + sqrt(5) r + 5/3 r^2) * exp(-sqrt(5) r)
    const Eigen::ArrayXd diff = (xi - xj).array() / theta_.array();
    const double r2 = diff.square().sum();
    const double r = std::sqrt(std::max(r2, 0.0));
    const double sqrt5_r = std::sqrt(5.0) * r;
    return sigma_f_ * sigma_f_ * (1.0 + sqrt5_r + (5.0 / 3.0) * r * r) *
           std::exp(-sqrt5_r);
}

void UniversalKriging::fit(const Eigen::MatrixXd& X, const Eigen::VectorXd& y) {
    if (X.rows() != y.size()) {
        throw std::invalid_argument("UniversalKriging::fit: X/y size mismatch");
    }
    if (X.cols() != theta_.size()) {
        throw std::invalid_argument(
            "UniversalKriging::fit: X column count must match length_scales size");
    }

    X_ = X;
    const int n = static_cast<int>(X.rows());

    residual_.resize(n);
    for (int i = 0; i < n; ++i) {
        residual_(i) = y(i) - trend_(X.row(i));
    }

    K_.resize(n, n);
    for (int i = 0; i < n; ++i) {
        for (int j = i; j < n; ++j) {
            const double kij = kernel(X.row(i), X.row(j));
            K_(i, j) = kij;
            K_(j, i) = kij;
        }
        K_(i, i) += sigma_n_ * sigma_n_;  // nugget for noise + conditioning
    }

    K_ldlt_.compute(K_);
    if (K_ldlt_.info() != Eigen::Success) {
        throw std::runtime_error(
            "UniversalKriging::fit: covariance matrix factorization failed "
            "(check length_scales / point spacing / sigma_n)");
    }
    alpha_ = K_ldlt_.solve(residual_);
    K_inv_ = K_ldlt_.solve(Eigen::MatrixXd::Identity(n, n));
}

UniversalKriging::Prediction UniversalKriging::predict(
    const Eigen::VectorXd& x) const {
    const int n = static_cast<int>(X_.rows());
    Eigen::VectorXd k_vec(n);
    for (int i = 0; i < n; ++i) {
        k_vec(i) = kernel(x, X_.row(i));
    }

    const double gp_mean = k_vec.dot(alpha_);
    const double explained_var = k_vec.dot(K_ldlt_.solve(k_vec));
    const double variance =
        std::max(sigma_f_ * sigma_f_ - explained_var, 0.0);  // numerical floor

    Prediction pred;
    pred.mean = trend_(x) + gp_mean;
    pred.variance = variance;
    return pred;
}

double UniversalKriging::looCV_RMSE() const {
    // Rasmussen & Williams (2006) eq. 5.12: the leave-one-out residual at
    // training point i is alpha_i / (K^{-1})_{ii}, computed in closed form
    // from the already-factorized system -- no need to refit N times.
    const int n = static_cast<int>(residual_.size());
    double sse = 0.0;
    for (int i = 0; i < n; ++i) {
        const double loo_error = alpha_(i) / K_inv_(i, i);
        sse += loo_error * loo_error;
    }
    return std::sqrt(sse / static_cast<double>(n));
}

}  // namespace aero_model
