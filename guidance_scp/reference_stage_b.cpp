#include "reference_stage_b.h"
#include "ClarabelSocpSolver.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cmath>

namespace guidance_scp {

namespace {

constexpr double kPi = 3.14159265358979323846;

// Fixed variable layout: [r_0..r_{K-1}, v_0..v_{K-1}, z_0..z_{K-1} (7/node)],
// [u_0..u_{K-2} (3/interval), sigma_u_0..sigma_u_{K-2} (1/interval)].
struct VarLayout {
    int K;
    int r_offset(int k) const { return 7 * k; }
    int v_offset(int k) const { return 7 * k + 3; }
    int z_offset(int k) const { return 7 * k + 6; }
    int u_offset(int k) const { return 7 * K + 4 * k; }
    int sigma_u_offset(int k) const { return 7 * K + 4 * k + 3; }
    int total() const { return 7 * K + 4 * (K - 1); }
};

struct FixedTfResult {
    bool solved = false;
    Eigen::VectorXd x;
};

}

// Exposed (declared in reference_stage_b.h) so GuidanceScpTests.cpp can
// directly inspect the assembled cones, confirming the thrust bound is a
// SecondOrderConeT (lossless convexification), not a plain box constraint.
SocpProblem buildFixedTfProblem(const StageBState& initial, const StageBState& terminal,
                                  const PlanetConfig& planet_config, const StageBConfig& config, double tf) {
    const int K = config.K;
    const double dt = tf / (K - 1);
    const double g0 = planet_config.g_0;
    const double alpha = 1.0 / (config.Isp_s * g0);  // zdot = -alpha*sigma_u
    VarLayout layout{K};
    int n = layout.total();

    std::vector<Eigen::Triplet<double>> A_triplets;
    std::vector<double> b_vec;
    std::vector<clarabel::SupportedConeT<double>> cones;
    int row = 0;
    auto addRow = [&](const std::vector<std::pair<int, double>>& coeffs, double rhs) {
        for (const auto& [col, val] : coeffs) A_triplets.emplace_back(row, col, val);
        b_vec.push_back(rhs);
        ++row;
    };

    // Initial boundary: r_0, v_0, z_0 fixed exactly.
    double z0_initial = std::log(initial.m);
    addRow({{layout.r_offset(0) + 0, 1.0}}, initial.rx); cones.push_back(clarabel::ZeroConeT<double>(1));
    addRow({{layout.r_offset(0) + 1, 1.0}}, initial.ry); cones.push_back(clarabel::ZeroConeT<double>(1));
    addRow({{layout.r_offset(0) + 2, 1.0}}, initial.rz); cones.push_back(clarabel::ZeroConeT<double>(1));
    addRow({{layout.v_offset(0) + 0, 1.0}}, initial.vx); cones.push_back(clarabel::ZeroConeT<double>(1));
    addRow({{layout.v_offset(0) + 1, 1.0}}, initial.vy); cones.push_back(clarabel::ZeroConeT<double>(1));
    addRow({{layout.v_offset(0) + 2, 1.0}}, initial.vz); cones.push_back(clarabel::ZeroConeT<double>(1));
    addRow({{layout.z_offset(0), 1.0}}, z0_initial); cones.push_back(clarabel::ZeroConeT<double>(1));

    // Terminal boundary: r_{K-1}, v_{K-1} matched exactly; z_{K-1} free (maximized).
    int kf = K - 1;
    addRow({{layout.r_offset(kf) + 0, 1.0}}, terminal.rx); cones.push_back(clarabel::ZeroConeT<double>(1));
    addRow({{layout.r_offset(kf) + 1, 1.0}}, terminal.ry); cones.push_back(clarabel::ZeroConeT<double>(1));
    addRow({{layout.r_offset(kf) + 2, 1.0}}, terminal.rz); cones.push_back(clarabel::ZeroConeT<double>(1));
    addRow({{layout.v_offset(kf) + 0, 1.0}}, terminal.vx); cones.push_back(clarabel::ZeroConeT<double>(1));
    addRow({{layout.v_offset(kf) + 1, 1.0}}, terminal.vy); cones.push_back(clarabel::ZeroConeT<double>(1));
    addRow({{layout.v_offset(kf) + 2, 1.0}}, terminal.vz); cones.push_back(clarabel::ZeroConeT<double>(1));

    // Dynamics (exact: rdot=v, vdot=u+g, zdot=-alpha*sigma_u are already
    // linear/affine, so ZOH-per-interval control integrates in closed form
    // with no linearization needed): g = [0,0,-g0].
    double gz = -g0;
    for (int k = 0; k < K - 1; ++k) {
        for (int i = 0; i < 3; ++i) {
            double g_i = (i == 2) ? gz : 0.0;
            addRow({{layout.r_offset(k + 1) + i, 1.0}, {layout.r_offset(k) + i, -1.0},
                    {layout.v_offset(k) + i, -dt}, {layout.u_offset(k) + i, -0.5 * dt * dt}},
                   0.5 * dt * dt * g_i);
            cones.push_back(clarabel::ZeroConeT<double>(1));
        }
        for (int i = 0; i < 3; ++i) {
            double g_i = (i == 2) ? gz : 0.0;
            addRow({{layout.v_offset(k + 1) + i, 1.0}, {layout.v_offset(k) + i, -1.0},
                    {layout.u_offset(k) + i, -dt}},
                   dt * g_i);
            cones.push_back(clarabel::ZeroConeT<double>(1));
        }
        addRow({{layout.z_offset(k + 1), 1.0}, {layout.z_offset(k), -1.0}, {layout.sigma_u_offset(k), dt * alpha}},
               0.0);
        cones.push_back(clarabel::ZeroConeT<double>(1));
    }

    // Lossless-convexified thrust bound: ||u_k|| <= sigma_u_k (SOC, dim 4).
    for (int k = 0; k < K - 1; ++k) {
        addRow({{layout.sigma_u_offset(k), -1.0}}, 0.0);
        addRow({{layout.u_offset(k) + 0, -1.0}}, 0.0);
        addRow({{layout.u_offset(k) + 1, -1.0}}, 0.0);
        addRow({{layout.u_offset(k) + 2, -1.0}}, 0.0);
        cones.push_back(clarabel::SecondOrderConeT<double>(4));
    }

    // Mass-lower-bound-derived linear relaxation on sigma_u_k.
    for (int k = 0; k < K - 1; ++k) {
        double t_k = k * dt;
        double m0_k = initial.m - (config.Tmax_N / (config.Isp_s * g0)) * t_k;
        double z0_k = std::log(m0_k);
        double n0_k = config.Tmin_N / m0_k;
        double n1_k = config.Tmax_N / m0_k;
        // sigma_u_k + n0_k*z_k >= n0_k*(1+z0_k)  =>  -sigma_u_k - n0_k*z_k <= -n0_k*(1+z0_k)
        addRow({{layout.sigma_u_offset(k), -1.0}, {layout.z_offset(k), -n0_k}}, -n0_k * (1.0 + z0_k));
        cones.push_back(clarabel::NonnegativeConeT<double>(1));
        // n1_k*z_k + sigma_u_k <= n1_k*(1+z0_k)
        addRow({{layout.z_offset(k), n1_k}, {layout.sigma_u_offset(k), 1.0}}, n1_k * (1.0 + z0_k));
        cones.push_back(clarabel::NonnegativeConeT<double>(1));
    }

    // Glideslope cone: ||[rx_k-rx_f, ry_k-ry_f]|| <= (rz_k-rz_f)/tan(half-angle from vertical), SOC dim 3.
    double c = 1.0 / std::tan(config.glideslope_deg * kPi / 180.0);
    for (int k = 1; k < K; ++k) {  // skip k=0 (initial point may not itself satisfy the cone)
        addRow({{layout.r_offset(k) + 2, -c}}, -c * terminal.rz);
        addRow({{layout.r_offset(k) + 0, -1.0}}, -terminal.rx);
        addRow({{layout.r_offset(k) + 1, -1.0}}, -terminal.ry);
        cones.push_back(clarabel::SecondOrderConeT<double>(3));
    }

    Eigen::SparseMatrix<double> A(row, n);
    A.setFromTriplets(A_triplets.begin(), A_triplets.end());
    Eigen::VectorXd b = Eigen::Map<Eigen::VectorXd>(b_vec.data(), b_vec.size());
    Eigen::SparseMatrix<double> P(n, n);  // pure linear cost

    Eigen::VectorXd q = Eigen::VectorXd::Zero(n);
    q(layout.z_offset(K - 1)) = -1.0;  // minimize -z_{K-1} = maximize final ln-mass

    return SocpProblem{P, q, A, b, cones};
}

namespace {

FixedTfResult solveFixedTf(const StageBState& initial, const StageBState& terminal,
                             const PlanetConfig& planet_config, const StageBConfig& config, double tf) {
    SocpProblem problem = buildFixedTfProblem(initial, terminal, planet_config, config, tf);
    SocpSolution solution = solve(problem);

    FixedTfResult result;
    result.solved = (solution.status == clarabel::SolverStatus::Solved
                      || solution.status == clarabel::SolverStatus::AlmostSolved);
    result.x = solution.x;
    return result;
}

StageBTrajectoryHistory extractHistory(const FixedTfResult& fixed, const StageBConfig& config, double tf) {
    const int K = config.K;
    const double dt = tf / (K - 1);
    VarLayout layout{K};
    StageBTrajectoryHistory hist;
    for (int k = 0; k < K; ++k) {
        hist.t.push_back(k * dt);
        double z_k = fixed.x(layout.z_offset(k));
        double m_k = std::exp(z_k);
        hist.m.push_back(m_k);
        hist.rx.push_back(fixed.x(layout.r_offset(k) + 0));
        hist.ry.push_back(fixed.x(layout.r_offset(k) + 1));
        hist.rz.push_back(fixed.x(layout.r_offset(k) + 2));
        hist.vx.push_back(fixed.x(layout.v_offset(k) + 0));
        hist.vy.push_back(fixed.x(layout.v_offset(k) + 1));
        hist.vz.push_back(fixed.x(layout.v_offset(k) + 2));
        if (k < K - 1) {
            double ux = fixed.x(layout.u_offset(k) + 0);
            double uy = fixed.x(layout.u_offset(k) + 1);
            double uz = fixed.x(layout.u_offset(k) + 2);
            double Tx = m_k * ux, Ty = m_k * uy, Tz = m_k * uz;
            hist.Tx.push_back(Tx); hist.Ty.push_back(Ty); hist.Tz.push_back(Tz);
            double thrust_mag = std::sqrt(Tx * Tx + Ty * Ty + Tz * Tz);
            hist.throttle_frac.push_back(config.Tmax_N > 0.0 ? thrust_mag / config.Tmax_N : 0.0);
        }
    }
    return hist;
}

}

StageBResult solveStageB(const StageBState& initial, const StageBState& terminal,
                          const PlanetConfig& planet_config, const StageBConfig& config) {
    StageBResult result;

    FixedTfResult hi_attempt = solveFixedTf(initial, terminal, planet_config, config, config.tf_max_s);
    if (!hi_attempt.solved) {
        result.solved = false;
        return result;  // even the most generous t_f in the search bracket is infeasible
    }

    double lo = config.tf_min_s, hi = config.tf_max_s;
    FixedTfResult best = hi_attempt;
    double best_tf = hi;
    for (int i = 0; i < config.tf_search_iters; ++i) {
        double mid = 0.5 * (lo + hi);
        FixedTfResult mid_attempt = solveFixedTf(initial, terminal, planet_config, config, mid);
        if (mid_attempt.solved) {
            hi = mid;
            best = mid_attempt;
            best_tf = mid;
        } else {
            lo = mid;
        }
    }

    result.solved = true;
    result.tf_s = best_tf;
    result.history = extractHistory(best, config, best_tf);
    return result;
}

}  // namespace guidance_scp
