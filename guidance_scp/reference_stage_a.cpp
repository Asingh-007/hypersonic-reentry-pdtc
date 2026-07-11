#include "reference_stage_a.h"
#include "ClarabelSocpSolver.h"
#include "TrimSolver.h"
#include "DescentDynamics.h"
#include "HeatingLoadModel.h"
#include "AtmosphereModel.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>

namespace guidance_scp {

namespace {

using Vec7 = Eigen::Matrix<double, 7, 1>;

Vec7 toVec(const StageAState& s) {
    Vec7 v;
    v << s.r, s.la, s.lo, s.v, s.fpa, s.v_azi, s.sigma;
    return v;
}

StageAState fromVec(const Vec7& v) {
    StageAState s;
    s.r = v(0); s.la = v(1); s.lo = v(2); s.v = v(3); s.fpa = v(4); s.v_azi = v(5); s.sigma = v(6);
    return s;
}

}

AeroClosure computeAeroClosure(double r, double v, const PlanetConfig& planet_config,
                                const SpacecraftConfig& spacecraft_config) {
    AeroClosure c;
    c.rho = DescentDynamics::atmosphereDensity(r, planet_config);
    double a_sound = DescentDynamics::speedOfSound(r, planet_config);
    c.mach = v / a_sound;
    c.alpha_trim_deg = alphaTrimDeg(c.mach, spacecraft_config.aero_table);
    auto aero = spacecraft_config.aero_table.interpolate(c.mach, c.alpha_trim_deg, /*beta_deg=*/0.0,
                                                           0.0, 0.0, 0.0);
    c.qbar_pa = 0.5 * c.rho * v * v;
    c.lift_n = c.qbar_pa * spacecraft_config.S_ref * aero.CL;
    c.drag_n = c.qbar_pa * spacecraft_config.S_ref * aero.CD;
    return c;
}

PathConstraintValues computePathConstraintValues(double r, double v, const PlanetConfig& planet_config,
                                                    const SpacecraftConfig& spacecraft_config) {
    AeroClosure c = computeAeroClosure(r, v, planet_config, spacecraft_config);
    PathConstraintValues pcv;
    pcv.qdot_w_m2 = suttonGravesHeatFlux(c.rho, v, spacecraft_config.nose_radius_m);
    pcv.qbar_pa = c.qbar_pa;
    double force_mag = std::sqrt(c.lift_n * c.lift_n + c.drag_n * c.drag_n);
    pcv.n_g = force_mag / (spacecraft_config.mass * planet_config.g_0);
    return pcv;
}

// Closed-form derivation (all three reduce to "-multiple_of_kBeta * quantity"
// given d(rho)/dr = -kBeta*rho in EarthAtmosphere1976's low-altitude
// exponential regime, h<84km -- see AtmosphereModel.h):
//   Qdot ~ sqrt(rho)*V^3  =>  dQdot/dr = -(kBeta/2)*Qdot,  dQdot/dV = 3*Qdot/V
//   qbar = 0.5*rho*V^2    =>  dqbar/dr = -kBeta*qbar,      dqbar/dV = 2*qbar/V
//   n = sqrt(L^2+D^2)/(m*g0), L,D ~ qbar (CL/CD's own Mach-dependence not
//       separately differentiated -- PLACEHOLDER, see reference_stage_a.h)
//                          =>  dn/dr = -kBeta*n,           dn/dV = 2*n/V
PathConstraintGradients computePathConstraintGradients(double r, double v, const PlanetConfig& planet_config,
                                                          const SpacecraftConfig& spacecraft_config) {
    PathConstraintValues pcv = computePathConstraintValues(r, v, planet_config, spacecraft_config);
    double kBeta = EarthAtmosphere1976::kBeta;

    PathConstraintGradients g;
    g.dQdot_dr = -(kBeta / 2.0) * pcv.qdot_w_m2;
    g.dQdot_dV = 3.0 * pcv.qdot_w_m2 / v;
    g.dqbar_dr = -kBeta * pcv.qbar_pa;
    g.dqbar_dV = 2.0 * pcv.qbar_pa / v;
    g.dn_dr = -kBeta * pcv.n_g;
    g.dn_dV = 2.0 * pcv.n_g / v;
    return g;
}

namespace {

// Stage A's reduced 7-state EOM: DescentDynamics::derivatives()'s existing
// translational terms with thrust=0 (Stage A is engines-off) and the scalar
// `lift` resolved into bank-angle components -- L*cos(sigma) in place of
// `lift` in the vertical-plane (fpa_dot) equation, and a NEW L*sin(sigma)
// lateral term in v_azi_dot (currently zero in DescentDynamics, which has no
// bank angle at all).
Vec7 stageAEom(const Vec7& x, double sigma_dot, const PlanetConfig& planet_config,
               const SpacecraftConfig& spacecraft_config) {
    double r = x(0), la = x(1), v = x(3), fpa = x(4), v_azi = x(5), sigma = x(6);

    AeroClosure aero = computeAeroClosure(r, v, planet_config, spacecraft_config);
    double L = aero.lift_n, D = aero.drag_n;
    double mass = spacecraft_config.mass;

    double mu = planet_config.mu, j2 = planet_config.j2, r_ref = planet_config.radius;
    double sin_la = std::sin(la), cos_la = std::cos(la);
    double r_ratio = r_ref / r;
    double g_r = (mu / (r * r)) * (1.0 - 3.0 * j2 * (r_ratio * r_ratio) * 0.5 * (3.0 * sin_la * sin_la - 1.0));
    double g_lambda = (3.0 * mu / (r * r)) * j2 * (r_ratio * r_ratio) * cos_la * sin_la;
    double omega = planet_config.omega;

    double cos_fpa = std::cos(fpa), sin_fpa = std::sin(fpa);
    double cos_vazi = std::cos(v_azi), sin_vazi = std::sin(v_azi);

    double r_dot = v * sin_fpa;
    double la_dot = v / r * cos_fpa * cos_vazi;
    double lo_dot = (v * cos_fpa * sin_vazi) / (r * cos_la);

    double v_dot = -D / mass - g_r * sin_fpa + g_lambda * cos_fpa * cos_vazi
                  - omega * omega * r * cos_la * (cos_fpa * cos_vazi * cos_la - sin_fpa * sin_la);

    double fpa_dot_term = mass * v * v / r * cos_fpa + L * std::cos(sigma)
                        - mass * g_r * cos_fpa - mass * g_lambda * sin_fpa * cos_vazi
                        + mass * omega * omega * r * cos_la * (sin_fpa * cos_vazi * sin_la - cos_fpa * cos_la)
                        + 2 * mass * omega * v * sin_vazi * cos_la;
    double fpa_dot = fpa_dot_term / (mass * v);

    double v_azi_dot_term = mass * v * v / r * cos_fpa * cos_fpa * sin_vazi * std::tan(la)
                        + L * std::sin(sigma) - mass * g_lambda * sin_vazi
                        + mass * omega * omega * r * sin_vazi * sin_la * cos_la
                        - 2 * mass * omega * v * (sin_fpa * cos_vazi * cos_la - cos_fpa * sin_la);
    double v_azi_dot = v_azi_dot_term / (mass * v * cos_fpa);

    Vec7 xdot;
    xdot << r_dot, la_dot, lo_dot, v_dot, fpa_dot, v_azi_dot, sigma_dot;
    return xdot;
}

// Central-finite-difference Jacobian of stageAEom w.r.t. (x, u) -- see this
// file's top-of-file comment for why finite differences are used here
// rather than an analytically-derived Jacobian.
void computeEomJacobian(const Vec7& x_ref, double u_ref, const PlanetConfig& planet_config,
                         const SpacecraftConfig& spacecraft_config,
                         Eigen::Matrix<double, 7, 7>& Ac, Eigen::Matrix<double, 7, 1>& Bc) {
    constexpr double kStepFrac = 1e-6;
    for (int i = 0; i < 7; ++i) {
        double h = std::max(std::abs(x_ref(i)), 1.0) * kStepFrac;
        Vec7 x_plus = x_ref, x_minus = x_ref;
        x_plus(i) += h;
        x_minus(i) -= h;
        Vec7 f_plus = stageAEom(x_plus, u_ref, planet_config, spacecraft_config);
        Vec7 f_minus = stageAEom(x_minus, u_ref, planet_config, spacecraft_config);
        Ac.col(i) = (f_plus - f_minus) / (2.0 * h);
    }
    double hu = std::max(std::abs(u_ref), 1.0) * kStepFrac;
    Vec7 f_plus = stageAEom(x_ref, u_ref + hu, planet_config, spacecraft_config);
    Vec7 f_minus = stageAEom(x_ref, u_ref - hu, planet_config, spacecraft_config);
    Bc = (f_plus - f_minus) / (2.0 * hu);
}

// Characteristic per-component scales for the trust region and virtual-control cost .
constexpr std::array<double, 7> kStateScale = {
    1000.0,  // r, m
    0.01,    // la, rad
    0.01,    // lo, rad
    100.0,   // v, m/s
    0.01,    // fpa, rad
    0.01,    // v_azi, rad
    0.1,     // sigma, rad
};
constexpr double kControlScale = 0.01;  // sigma_dot, rad/s

// Adaptive-trust-region accept/reject check
double computeMaxNonlinearDefectStageA(const std::vector<Vec7>& x_ref_candidate,
                                        const std::vector<double>& u_ref_candidate, double dt,
                                        const PlanetConfig& planet_config,
                                        const SpacecraftConfig& spacecraft_config) {
    double max_defect = 0.0;
    const int K = static_cast<int>(x_ref_candidate.size());
    for (int k = 0; k + 1 < K; ++k) {
        Vec7 xdot = stageAEom(x_ref_candidate[k], u_ref_candidate[k], planet_config, spacecraft_config);
        Vec7 predicted_next = x_ref_candidate[k] + dt * xdot;
        Vec7 defect = x_ref_candidate[k + 1] - predicted_next;
        for (int i = 0; i < 7; ++i) {
            max_defect = std::max(max_defect, std::abs(defect(i)) / kStateScale[i]);
        }
    }
    return max_defect;
}

// Fixed variable layout within the single stacked decision vector, in the
// order: [delta_x_0..delta_x_{K-1} (7 each)], [delta_u_0..delta_u_{K-2} (1
// each)], [nu_plus_0..nu_plus_{K-2}, nu_minus_0..nu_minus_{K-2} (7 each)],
// [eta_0..eta_{K-1} (1 each)].
struct VarLayout {
    int K;
    int dx_offset(int k) const { return 7 * k; }
    int du_offset(int k) const { return 7 * K + k; }
    int nu_plus_offset(int k) const { return 7 * K + (K - 1) + 7 * k; }
    int nu_minus_offset(int k) const { return 7 * K + (K - 1) + 7 * (K - 1) + 7 * k; }
    int eta_offset(int k) const { return 7 * K + (K - 1) + 14 * (K - 1) + k; }
    int total() const { return 7 * K + (K - 1) + 14 * (K - 1) + K; }
};

}

StageAResult solveStageA(const StageAState& initial, const StageAState& terminal,
                          const PlanetConfig& planet_config, const SpacecraftConfig& spacecraft_config,
                          const StageAConfig& config) {
    const int K = config.K;
    const double dt = config.t_scale_s / (K - 1);
    VarLayout layout{K};

    std::vector<Vec7> x_ref(K);
    std::vector<double> u_ref(K - 1, 0.0);

    // Initial reference guess: propagate the real nonlinear EOM forward from
    // `initial` at zero bank rate (sigma held at initial.sigma throughout).
    x_ref[0] = toVec(initial);
    for (int k = 0; k + 1 < K; ++k) {
        Vec7 xdot = stageAEom(x_ref[k], 0.0, planet_config, spacecraft_config);
        x_ref[k + 1] = x_ref[k] + dt * xdot;
    }

    StageAResult result;

    // Adaptive trust-region state (see StageAConfig::eta_max_* doc comment):
    // eta_max is a HARD cap on every node's eta_k, persisted/updated across
    // solve attempts (both accepted and rejected), distinct from max_iters
    // (accepted iterations only) and total_attempts (accepted+rejected).
    double eta_max = config.eta_max_init;
    int accepted_iters = 0;
    int total_attempts = 0;
    int consecutive_floor_rejects = 0;
    double best_max_nu = std::numeric_limits<double>::infinity();  // best (smallest) max||nu|| across accepted iterations

    while (accepted_iters < config.max_iters && total_attempts < config.max_solve_attempts) {
        ++total_attempts;
        std::vector<Eigen::Matrix<double, 7, 7>> A_k(K - 1);
        std::vector<Vec7> B_k(K - 1);
        for (int k = 0; k < K - 1; ++k) {
            Eigen::Matrix<double, 7, 7> Ac;
            Vec7 Bc;
            computeEomJacobian(x_ref[k], u_ref[k], planet_config, spacecraft_config, Ac, Bc);
            A_k[k] = Eigen::Matrix<double, 7, 7>::Identity() + dt * Ac;
            B_k[k] = dt * Bc;
        }

        std::vector<PathConstraintValues> pcv(K);
        std::vector<PathConstraintGradients> pcg(K);
        for (int k = 0; k < K; ++k) {
            double r_k = x_ref[k](0), v_k = x_ref[k](3);
            pcv[k] = computePathConstraintValues(r_k, v_k, planet_config, spacecraft_config);
            pcg[k] = computePathConstraintGradients(r_k, v_k, planet_config, spacecraft_config);
        }

        int n = layout.total();
        std::vector<Eigen::Triplet<double>> A_triplets;
        std::vector<double> b_vec;
        std::vector<clarabel::SupportedConeT<double>> cones;
        Eigen::VectorXd q = Eigen::VectorXd::Zero(n);

        int row = 0;
        auto addEqRow = [&](const std::vector<std::pair<int, double>>& coeffs, double rhs) {
            for (const auto& [col, val] : coeffs) A_triplets.emplace_back(row, col, val);
            b_vec.push_back(rhs);
            ++row;
        };
        auto addIneqRow = [&](const std::vector<std::pair<int, double>>& coeffs, double rhs) {
            for (const auto& [col, val] : coeffs) A_triplets.emplace_back(row, col, val);
            b_vec.push_back(rhs);
            ++row;
        };

        // Initial boundary: delta_x_0 = 0
        for (int i = 0; i < 7; ++i) {
            addEqRow({{layout.dx_offset(0) + i, 1.0}}, 0.0);
        }
        cones.push_back(clarabel::ZeroConeT<double>(7));

        // Dynamics: delta_x_{k+1} - A_k*delta_x_k - B_k*delta_u_k - nu_plus_k + nu_minus_k = 0
        for (int k = 0; k < K - 1; ++k) {
            for (int i = 0; i < 7; ++i) {
                std::vector<std::pair<int, double>> coeffs;
                coeffs.emplace_back(layout.dx_offset(k + 1) + i, 1.0);
                for (int j = 0; j < 7; ++j) {
                    if (A_k[k](i, j) != 0.0) coeffs.emplace_back(layout.dx_offset(k) + j, -A_k[k](i, j));
                }
                coeffs.emplace_back(layout.du_offset(k), -B_k[k](i));
                coeffs.emplace_back(layout.nu_plus_offset(k) + i, -1.0);
                coeffs.emplace_back(layout.nu_minus_offset(k) + i, 1.0);
                addEqRow(coeffs, 0.0);
            }
            cones.push_back(clarabel::ZeroConeT<double>(7));
        }

        // Terminal boundary: delta_x_{K-1} = terminal - x_ref_{K-1}
        {
            Vec7 target_delta = toVec(terminal) - x_ref[K - 1];
            for (int i = 0; i < 7; ++i) {
                addEqRow({{layout.dx_offset(K - 1) + i, 1.0}}, target_delta(i));
            }
            cones.push_back(clarabel::ZeroConeT<double>(7));
        }

        // Bank-angle box: sigma_max - (sigma_ref_k + delta_sigma_k) >= 0, and the mirrored lower bound.
        for (int k = 0; k < K; ++k) {
            addIneqRow({{layout.dx_offset(k) + 6, 1.0}}, config.sigma_max - x_ref[k](6));
            cones.push_back(clarabel::NonnegativeConeT<double>(1));
            addIneqRow({{layout.dx_offset(k) + 6, -1.0}}, config.sigma_max + x_ref[k](6));
            cones.push_back(clarabel::NonnegativeConeT<double>(1));
        }
        // Bank-rate (control) box.
        for (int k = 0; k < K - 1; ++k) {
            addIneqRow({{layout.du_offset(k), 1.0}}, config.sigma_dot_max - u_ref[k]);
            cones.push_back(clarabel::NonnegativeConeT<double>(1));
            addIneqRow({{layout.du_offset(k), -1.0}}, config.sigma_dot_max + u_ref[k]);
            cones.push_back(clarabel::NonnegativeConeT<double>(1));
        }

        // Path constraints: f_max - f_ref - df/dr*delta_r_k - df/dV*delta_V_k >= 0.
        for (int k = 0; k < K; ++k) {
            addIneqRow({{layout.dx_offset(k) + 0, -pcg[k].dQdot_dr}, {layout.dx_offset(k) + 3, -pcg[k].dQdot_dV}},
                       config.max_heat_flux_w_m2 - pcv[k].qdot_w_m2);
            cones.push_back(clarabel::NonnegativeConeT<double>(1));
            addIneqRow({{layout.dx_offset(k) + 0, -pcg[k].dqbar_dr}, {layout.dx_offset(k) + 3, -pcg[k].dqbar_dV}},
                       config.max_qbar_pa - pcv[k].qbar_pa);
            cones.push_back(clarabel::NonnegativeConeT<double>(1));
            addIneqRow({{layout.dx_offset(k) + 0, -pcg[k].dn_dr}, {layout.dx_offset(k) + 3, -pcg[k].dn_dV}},
                       config.max_load_factor_g - pcv[k].n_g);
            cones.push_back(clarabel::NonnegativeConeT<double>(1));
        }

        // Trust region (L-infinity box, shared per-node DIMENSIONLESS scalar
        // eta_k): |delta_x_k,i| <= kStateScale[i]*eta_k for all i, and
        // |delta_u_k| <= kControlScale*eta_k where a control exists at that
        // node -> per-component scaling is essential here (see kStateScale's
        // doc comment above). addIneqRow enforces (coeffs.x <= rhs), so
        // "delta_x_k,i <= scale*eta_k" is coeffs={dx:1, eta:-scale}, rhs=0,
        // and "-delta_x_k,i <= scale*eta_k" (i.e. delta_x_k,i >= -scale*eta_k)
        // is coeffs={dx:-1, eta:-scale}, rhs=0.
        for (int k = 0; k < K; ++k) {
            for (int i = 0; i < 7; ++i) {
                addIneqRow({{layout.dx_offset(k) + i, 1.0}, {layout.eta_offset(k), -kStateScale[i]}}, 0.0);
                cones.push_back(clarabel::NonnegativeConeT<double>(1));
                addIneqRow({{layout.dx_offset(k) + i, -1.0}, {layout.eta_offset(k), -kStateScale[i]}}, 0.0);
                cones.push_back(clarabel::NonnegativeConeT<double>(1));
            }
            if (k < K - 1) {
                addIneqRow({{layout.du_offset(k), 1.0}, {layout.eta_offset(k), -kControlScale}}, 0.0);
                cones.push_back(clarabel::NonnegativeConeT<double>(1));
                addIneqRow({{layout.du_offset(k), -1.0}, {layout.eta_offset(k), -kControlScale}}, 0.0);
                cones.push_back(clarabel::NonnegativeConeT<double>(1));
            }
            // eta_k >= 0, i.e. -eta_k <= 0.
            addIneqRow({{layout.eta_offset(k), -1.0}}, 0.0);
            cones.push_back(clarabel::NonnegativeConeT<double>(1));
        }

        
        for (int k = 1; k < K - 1; ++k) {
            for (int i = 0; i < 7; ++i) {
                addIneqRow({{layout.dx_offset(k) + i, 1.0}}, kStateScale[i] * eta_max);
                cones.push_back(clarabel::NonnegativeConeT<double>(1));
                addIneqRow({{layout.dx_offset(k) + i, -1.0}}, kStateScale[i] * eta_max);
                cones.push_back(clarabel::NonnegativeConeT<double>(1));
            }
        }
        for (int k = 0; k < K - 1; ++k) {
            addIneqRow({{layout.du_offset(k), 1.0}}, kControlScale * eta_max);
            cones.push_back(clarabel::NonnegativeConeT<double>(1));
            addIneqRow({{layout.du_offset(k), -1.0}}, kControlScale * eta_max);
            cones.push_back(clarabel::NonnegativeConeT<double>(1));
        }
        
        for (int k = 0; k < K - 1; ++k) {
            for (int i = 0; i < 7; ++i) {
                addIneqRow({{layout.nu_plus_offset(k) + i, -1.0}}, 0.0);
                cones.push_back(clarabel::NonnegativeConeT<double>(1));
                addIneqRow({{layout.nu_minus_offset(k) + i, -1.0}}, 0.0);
                cones.push_back(clarabel::NonnegativeConeT<double>(1));
            }
        }

        // Cost: w_nu * sum((nu_plus+nu_minus)/kStateScale[i]) + w_trust * sum(eta_k).
        // Normalizing nu's cost by each component's characteristic scale
        // (same reasoning as the trust region above) so a 1-unit virtual-
        // control defect in r doesn't cost 1e5x more than a 1-unit defect
        // in sigma purely because of unit choice, not physical significance.
        for (int k = 0; k < K - 1; ++k) {
            for (int i = 0; i < 7; ++i) {
                q(layout.nu_plus_offset(k) + i) = config.w_nu / kStateScale[i];
                q(layout.nu_minus_offset(k) + i) = config.w_nu / kStateScale[i];
            }
        }
        for (int k = 0; k < K; ++k) q(layout.eta_offset(k)) = config.w_trust;

        Eigen::SparseMatrix<double> A(row, n);
        A.setFromTriplets(A_triplets.begin(), A_triplets.end());
        Eigen::VectorXd b = Eigen::Map<Eigen::VectorXd>(b_vec.data(), b_vec.size());
        Eigen::SparseMatrix<double> P(n, n);  // all-zero: pure LP cost

        SocpProblem problem{P, q, A, b, cones};
        SocpSolution solution = solve(problem);

        // A non-Solved status (e.g. DualInfeasible/PrimalInfeasible) means
        // solution.x is not a valid minimizer, do NOT apply it as a
        // reference update. Treated exactly like a rejected step: shrink
        // eta_max and retry (classical trust-region philosophy on any
        // failure signal, shrink, don't grow).
        bool non_solved = (solution.status != clarabel::SolverStatus::Solved
                            && solution.status != clarabel::SolverStatus::AlmostSolved);

        double max_nu = 0.0, max_eta = 0.0, max_defect_nonlinear = std::numeric_limits<double>::infinity();
        std::vector<Vec7> x_ref_candidate;
        std::vector<double> u_ref_candidate;
        bool accept = false;

        if (!non_solved) {
            // Build the CANDIDATE updated reference -- do NOT touch
            // x_ref/u_ref yet, since this step may still be rejected below.
            x_ref_candidate = x_ref;
            u_ref_candidate = u_ref;
            for (int k = 0; k < K; ++k) {
                Vec7 delta_x;
                for (int i = 0; i < 7; ++i) delta_x(i) = solution.x(layout.dx_offset(k) + i);
                x_ref_candidate[k] += delta_x;
                max_eta = std::max(max_eta, solution.x(layout.eta_offset(k)));
            }
            for (int k = 0; k < K - 1; ++k) {
                u_ref_candidate[k] += solution.x(layout.du_offset(k));
                double nu_l1 = 0.0;
                for (int i = 0; i < 7; ++i) {
                    nu_l1 += solution.x(layout.nu_plus_offset(k) + i) + solution.x(layout.nu_minus_offset(k) + i);
                }
                max_nu = std::max(max_nu, nu_l1);
            }
            max_defect_nonlinear = computeMaxNonlinearDefectStageA(x_ref_candidate, u_ref_candidate, dt,
                                                                     planet_config, spacecraft_config);

            accept = !std::isfinite(best_max_nu) || (max_nu <= best_max_nu * config.nu_regression_tolerance);
        }

        result.attempt_max_defect_nonlinear.push_back(max_defect_nonlinear);
        result.attempt_eta_max.push_back(eta_max);
        result.attempt_accepted.push_back(accept);

        if (accept) {
            x_ref = x_ref_candidate;
            u_ref = u_ref_candidate;
            ++accepted_iters;
            consecutive_floor_rejects = 0;
            best_max_nu = std::min(best_max_nu, max_nu);

            result.max_nu_per_iter.push_back(max_nu);
            result.max_eta_per_iter.push_back(max_eta);

            if (max_nu < config.eps_dyn && max_eta < config.eps_tr) {
                result.converged = true;
                break;
            }
        } else {
           
            eta_max = std::max(eta_max * config.eta_max_shrink, config.eta_max_floor);

            if (eta_max <= config.eta_max_floor * (1.0 + 1e-9)) {
                ++consecutive_floor_rejects;
            } else {
                consecutive_floor_rejects = 0;
            }

            if (consecutive_floor_rejects >= config.max_consecutive_floor_rejects) {
                std::cerr << "Stage A: eta_max stuck at floor after " << consecutive_floor_rejects
                          << " consecutive rejections -- aborting honestly rather than spinning."
                          << std::endl;
                result.converged = false;
                break;
            }
        }
    }
    if (total_attempts >= config.max_solve_attempts && !result.converged) {
        std::cerr << "Stage A: exhausted max_solve_attempts (" << config.max_solve_attempts
                  << ") without convergence." << std::endl;
    }

    for (int k = 0; k < K; ++k) {
        StageAState s = fromVec(x_ref[k]);
        result.history.t.push_back(k * dt);
        result.history.r.push_back(s.r);
        result.history.la.push_back(s.la);
        result.history.lo.push_back(s.lo);
        result.history.v.push_back(s.v);
        result.history.fpa.push_back(s.fpa);
        result.history.v_azi.push_back(s.v_azi);
        result.history.sigma.push_back(s.sigma);
        if (k < K - 1) result.history.sigma_dot.push_back(u_ref[k]);

        AeroClosure aero = computeAeroClosure(s.r, s.v, planet_config, spacecraft_config);
        PathConstraintValues pcv = computePathConstraintValues(s.r, s.v, planet_config, spacecraft_config);
        result.history.qbar.push_back(pcv.qbar_pa);
        result.history.heat_flux_conv.push_back(pcv.qdot_w_m2);
        result.history.load_factor_g.push_back(pcv.n_g);
        result.history.alpha_trim_deg.push_back(aero.alpha_trim_deg);
        result.history.mach.push_back(aero.mach);
    }

    return result;
}

}  // namespace guidance_scp
