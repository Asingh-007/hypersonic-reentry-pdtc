#ifndef LATIN_HYPERCUBE_SAMPLER_H
#define LATIN_HYPERCUBE_SAMPLER_H
// LatinHypercubeSampler.h
//
// Space-filling design-of-experiments generator for selecting CFD anchor
// points in the (Mach, alpha, beta, flap-deflection, ...) envelope, used to
// populate the training set for UniversalKriging. A plain grid wastes CFD
// budget (coverage scales as n^d); Latin Hypercube Sampling gives full
// per-dimension stratification for a fixed n, and the maximin swap pass
// below improves the raw LHS draw so points aren't accidentally clustered
// or nearly collinear in the joint space.
//
// NOTE: this is an OFFLINE-only utility (see aero/GenerateAeroTable.cpp) --
// DescentDynamics never calls this at runtime.

#include <Eigen/Dense>
#include <random>
#include <string>
#include <vector>

namespace aero_model {

struct DesignVariable {
    std::string name;
    double lower;
    double upper;
};

class LatinHypercubeSampler {
public:
    explicit LatinHypercubeSampler(unsigned seed = 42) : rng_(seed) {}

    // Draws n_samples in [0,1]^n_dims via LHS (one stratified, jittered
    // sample per bin per dimension, independently permuted across
    // dimensions), then runs n_swap_iters randomized coordinate swaps,
    // each accepted only if it increases the minimum pairwise distance
    // between design points (Morris-Mitchell maximin criterion).
    Eigen::MatrixXd sample(int n_samples, int n_dims, int n_swap_iters = 5000);

    // Maps unit-hypercube samples to physical variable bounds.
    static Eigen::MatrixXd scaleToBounds(const Eigen::MatrixXd& unit_samples,
                                          const std::vector<DesignVariable>& vars);

    // Appends must-include physical-unit points (e.g. known trim AoA at
    // the belly-flop condition, nominal entry-interface Mach) to an LHS
    // design already in physical units. These are not moved by the
    // maximin optimizer since they are added after scaling.
    static Eigen::MatrixXd augmentWithFixedPoints(
        const Eigen::MatrixXd& physical_samples,
        const Eigen::MatrixXd& fixed_points);

private:
    double minPairwiseDistance(const Eigen::MatrixXd& pts) const;
    std::mt19937 rng_;
};

}  // namespace aero_model

#endif // LATIN_HYPERCUBE_SAMPLER_H
