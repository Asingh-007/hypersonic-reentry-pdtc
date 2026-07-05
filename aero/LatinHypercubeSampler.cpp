#include "LatinHypercubeSampler.h"

#include <algorithm>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace aero_model {

Eigen::MatrixXd LatinHypercubeSampler::sample(int n_samples, int n_dims,
                                               int n_swap_iters) {
    if (n_samples < 2 || n_dims < 1) {
        throw std::invalid_argument(
            "LatinHypercubeSampler::sample: need n_samples >= 2, n_dims >= 1");
    }

    std::uniform_real_distribution<double> jitter(0.0, 1.0);
    Eigen::MatrixXd design(n_samples, n_dims);

    // Independent stratified permutation per dimension.
    for (int d = 0; d < n_dims; ++d) {
        std::vector<int> bins(n_samples);
        std::iota(bins.begin(), bins.end(), 0);
        std::shuffle(bins.begin(), bins.end(), rng_);
        for (int i = 0; i < n_samples; ++i) {
            const double u = (static_cast<double>(bins[i]) + jitter(rng_)) /
                              static_cast<double>(n_samples);
            design(i, d) = u;
        }
    }

    // Maximin improvement: repeatedly propose swapping two rows' values in
    // a randomly chosen dimension (this preserves the per-dimension
    // stratification -- each column remains a permutation of the same
    // bins) and accept only if it improves the minimum pairwise distance.
    std::uniform_int_distribution<int> dim_pick(0, n_dims - 1);
    std::uniform_int_distribution<int> row_pick(0, n_samples - 1);

    double best_min_dist = minPairwiseDistance(design);
    for (int iter = 0; iter < n_swap_iters; ++iter) {
        const int d = dim_pick(rng_);
        const int i = row_pick(rng_);
        int j = row_pick(rng_);
        while (j == i) {
            j = row_pick(rng_);
        }

        std::swap(design(i, d), design(j, d));
        const double candidate_min_dist = minPairwiseDistance(design);

        if (candidate_min_dist > best_min_dist) {
            best_min_dist = candidate_min_dist;  // keep the swap
        } else {
            std::swap(design(i, d), design(j, d));  // revert
        }
    }

    return design;
}

Eigen::MatrixXd LatinHypercubeSampler::scaleToBounds(
    const Eigen::MatrixXd& unit_samples, const std::vector<DesignVariable>& vars) {
    if (static_cast<int>(vars.size()) != unit_samples.cols()) {
        throw std::invalid_argument(
            "scaleToBounds: vars.size() must match unit_samples.cols()");
    }
    Eigen::MatrixXd out = unit_samples;
    for (int d = 0; d < unit_samples.cols(); ++d) {
        const double lo = vars[d].lower;
        const double hi = vars[d].upper;
        out.col(d) = (lo + unit_samples.col(d).array() * (hi - lo)).matrix();
    }
    return out;
}

Eigen::MatrixXd LatinHypercubeSampler::augmentWithFixedPoints(
    const Eigen::MatrixXd& physical_samples, const Eigen::MatrixXd& fixed_points) {
    if (fixed_points.cols() != physical_samples.cols()) {
        throw std::invalid_argument(
            "augmentWithFixedPoints: column count mismatch");
    }
    Eigen::MatrixXd out(physical_samples.rows() + fixed_points.rows(),
                         physical_samples.cols());
    out.topRows(physical_samples.rows()) = physical_samples;
    out.bottomRows(fixed_points.rows()) = fixed_points;
    return out;
}

double LatinHypercubeSampler::minPairwiseDistance(const Eigen::MatrixXd& pts) const {
    double min_d2 = std::numeric_limits<double>::max();
    const int n = static_cast<int>(pts.rows());
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            const double d2 = (pts.row(i) - pts.row(j)).squaredNorm();
            min_d2 = std::min(min_d2, d2);
        }
    }
    return std::sqrt(min_d2);
}

}  // namespace aero_model
