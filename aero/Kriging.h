#ifndef KRIGING_H
#define KRIGING_H
// Kriging.h
//
// Universal Kriging (Gaussian Process regression with a physics-informed
// trend/mean function) used to correct the cheap, dense Newtonian panel
// sweep against a sparse set of expensive CFD anchor points:
//
//     C(x) = trend(x) + Z(x),      Z ~ GP(0, k(x, x'))
//
// trend(x) is evaluated from NewtonianAeroModel (a physics-based prior,
// not just a constant/polynomial), and Z(x) is a zero-mean GP fit to the
// residuals y_i - trend(x_i) at the CFD anchor points x_i. The posterior
// mean reverts to the Newtonian trend far from any CFD anchor and locks
// onto the CFD-corrected value near anchors, with a posterior variance
// that grows with distance from the nearest anchors -- useful for
// adaptive DOE refinement (see LatinHypercubeSampler.h).
//
// Kernel: anisotropic Matern 5/2, one length scale per input dimension
// (e.g. Mach, alpha, flap deflection), which is the standard choice for
// aerodynamic surrogate modeling (smooth but not infinitely differentiable
// like squared-exponential, which tends to over-smooth shock-driven
// nonlinearities).
//
// NOTE: this is an OFFLINE-only utility (see aero/GenerateAeroTable.cpp) --
// DescentDynamics never calls this at runtime, only a precomputed
// AeroCoefficientTable. No real CFD data exists yet, so this is not
// currently wired into the table generator's output; see the "FUTURE WORK"
// comment block in GenerateAeroTable.cpp for the intended hook.

#include <Eigen/Dense>
#include <functional>

namespace aero_model {

class UniversalKriging {
public:
    using TrendFn = std::function<double(const Eigen::VectorXd&)>;

    struct Prediction {
        double mean;
        double variance;
    };

    // trend         : physics-based mean function, e.g. wraps
    //                 NewtonianAeroModel::evaluate(...).CL
    // length_scales : one per input dimension, same units as the inputs
    //                 (e.g. Mach ~1-2, alpha ~10-20 deg, flap ~5-10 deg --
    //                 rescale inputs to comparable ranges before fitting
    //                 if the natural units differ wildly)
    // sigma_f       : GP signal standard deviation (residual scale)
    // sigma_n       : nugget / observation-noise standard deviation,
    //                 also regularizes the linear solve -- keep nonzero
    //                 even for "noise-free" CFD (e.g. 1e-4) for conditioning
    UniversalKriging(TrendFn trend, Eigen::VectorXd length_scales,
                      double sigma_f = 1.0, double sigma_n = 1e-4)
        : trend_(std::move(trend)),
          theta_(std::move(length_scales)),
          sigma_f_(sigma_f),
          sigma_n_(sigma_n) {}

    // X: N x d matrix of CFD anchor input points (rows = points)
    // y: N-vector of CFD-observed coefficient values at those points
    void fit(const Eigen::MatrixXd& X, const Eigen::VectorXd& y);

    // Posterior mean/variance of the blended (trend + GP) prediction at x.
    Prediction predict(const Eigen::VectorXd& x) const;

    // Leave-one-out cross-validation RMSE of the fitted model (Rasmussen &
    // Williams eq. 5.12) -- use this to sanity-check length_scales/sigma_f
    // choices via a small grid search before trusting the final surrogate.
    double looCV_RMSE() const;

private:
    double kernel(const Eigen::VectorXd& xi, const Eigen::VectorXd& xj) const;

    TrendFn trend_;
    Eigen::VectorXd theta_;
    double sigma_f_;
    double sigma_n_;

    Eigen::MatrixXd X_;
    Eigen::VectorXd residual_;
    Eigen::MatrixXd K_;
    Eigen::LDLT<Eigen::MatrixXd> K_ldlt_;
    Eigen::VectorXd alpha_;  // K^{-1} * residual_, cached for predict()
    Eigen::MatrixXd K_inv_;  // cached full inverse, used by looCV_RMSE()
};

}  // namespace aero_model

#endif // KRIGING_H
