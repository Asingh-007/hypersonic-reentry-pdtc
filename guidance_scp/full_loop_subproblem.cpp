#include "full_loop_subproblem.h"
#include "full_loop_phase1.h"
#include "full_loop_phase2.h"
#include "full_loop_path_constraints.h"
#include "full_loop_nonconvex.h"
#include "reference_stage_a.h"
#include "FiniteDifferenceJacobian.h"
#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <algorithm>
#include <array>
#include <cmath>

namespace guidance_scp {

namespace {

constexpr double kPi = 3.14159265358979323846;

// Same construction as full_loop_phase2.cpp's (private) thrustDirBody --
// duplicated here since that one is file-local; kept in exact sync (both
// mirror freestreamDirectionBody's structure).
Eigen::Vector3d thrustDirBody(double deltaE, double phiE) {
    return Eigen::Vector3d(-std::cos(deltaE) * std::cos(phiE), std::sin(phiE),
                            std::sin(deltaE) * std::cos(phiE));
}

// Angle between two attitude quaternions (unsigned, via |dot product|,
// robust to the q/-q double-cover ambiguity).
double quatErrorAngle(double q1, double q2, double q3, double q4,
                       const Eigen::Quaterniond& q_target) {
    Eigen::Quaterniond q(q4, q1, q2, q3);
    double d = std::clamp(std::abs(q.dot(q_target)), 0.0, 1.0);
    return 2.0 * std::acos(d);
}

// kStateScale1/kControlScale1/kStateScale2/kControlScaleVec2/
// kControlScaleSigma2/kControlScaleRest2 now live in full_loop_subproblem.h
// (exposed, not file-local) so full_scp_loop.cpp's outer loop can reuse the
// exact same normalization for the diagnostic nonlinear defect.

}

double computeMaxNonlinearDefectGeneric(const EomFn& f, const std::vector<Eigen::VectorXd>& x_ref_candidate,
                                          const std::vector<Eigen::VectorXd>& u_ref_candidate, double dt,
                                          const double* state_scale, int state_dim) {
    double max_defect = 0.0;
    const int K = static_cast<int>(x_ref_candidate.size());
    for (int k = 0; k + 1 < K; ++k) {
        Eigen::VectorXd xdot = f(x_ref_candidate[k], u_ref_candidate[k]);
        Eigen::VectorXd predicted_next = x_ref_candidate[k] + dt * xdot;
        Eigen::VectorXd defect = x_ref_candidate[k + 1] - predicted_next;
        for (int i = 0; i < state_dim; ++i) {
            max_defect = std::max(max_defect, std::abs(defect(i)) / state_scale[i]);
        }
    }
    return max_defect;
}

SocpProblem buildFullLoopSubproblem(const FullLoopSubproblemInputs& in) {
    const int K1 = static_cast<int>(in.x1_ref.size());
    const int K2 = static_cast<int>(in.x2_ref.size());
    FullLoopVarLayout layout{K1, K2};
    const int n = layout.total();

    const PlanetConfig& planet_config = in.planet_config;
    const SpacecraftConfig& spacecraft_config = *in.spacecraft_config;
    const FlapActuatorConfig& flap_config = in.flap_config;
    const GimbalActuatorConfig& gimbal_config = in.gimbal_config;
    const FullLoopConfig& cfg = in.full_loop_config;
    const double g0 = planet_config.g_0;
    const double alpha_mdot = 1.0 / (cfg.Isp_s * g0);

    const double dt1 = in.t_scale_1 / (K1 - 1);
    const double dt2 = in.t_scale_2 / (K2 - 1);

    std::vector<Eigen::Triplet<double>> A_triplets;
    std::vector<double> b_vec;
    std::vector<clarabel::SupportedConeT<double>> cones;
    Eigen::VectorXd q = Eigen::VectorXd::Zero(n);
    int row = 0;

    auto addRow = [&](const std::vector<std::pair<int, double>>& coeffs, double rhs) {
        for (const auto& [col, val] : coeffs) A_triplets.emplace_back(row, col, val);
        b_vec.push_back(rhs);
        ++row;
    };
    auto addConeRow = [&](const ConeRow& r) { addRow(r.coeffs, r.rhs); };
    auto addConeRowGroup = [&](const ConeRowGroup& g) {
        for (const auto& r : g.rows) addConeRow(r);
        cones.push_back(g.cone);
    };

    // Phase 1

    // Initial boundary: delta_x1_0 = 0 (x1_ref[0] is maintained == x1_initial
    // by the outer loop, mirroring Stage A's exact convention).
    for (int i = 0; i < 21; ++i) addRow({{layout.dx1_offset(0) + i, 1.0}}, 0.0);
    cones.push_back(clarabel::ZeroConeT<double>(21));

    for (int k = 0; k < K1 - 1; ++k) {
        EomFn f = [&](const Eigen::VectorXd& x, const Eigen::VectorXd& u) {
            return phase1Eom(x, u, planet_config, spacecraft_config, flap_config);
        };
        EomJacobian jac = computeEomJacobianFd(f, in.x1_ref[k], in.u1_ref[k]);
        Eigen::MatrixXd A1 = Eigen::MatrixXd::Identity(21, 21) + dt1 * jac.Ac;
        Eigen::MatrixXd B1 = dt1 * jac.Bc;

        for (int i = 0; i < 21; ++i) {
            std::vector<std::pair<int, double>> coeffs;
            coeffs.emplace_back(layout.dx1_offset(k + 1) + i, 1.0);
            for (int j = 0; j < 21; ++j) {
                if (A1(i, j) != 0.0) coeffs.emplace_back(layout.dx1_offset(k) + j, -A1(i, j));
            }
            for (int c = 0; c < 4; ++c) {
                if (B1(i, c) != 0.0) coeffs.emplace_back(layout.du1_offset(k) + c, -B1(i, c));
            }
            coeffs.emplace_back(layout.nu1_plus_offset(k) + i, -1.0);
            coeffs.emplace_back(layout.nu1_minus_offset(k) + i, 1.0);
            addRow(coeffs, 0.0);
            cones.push_back(clarabel::ZeroConeT<double>(1));
        }
    }

    // Path constraints:  Qdot/qbar (Stage A's (r,V)-only rows,
    // reused verbatim) + load factor (with the new flap-dependence terms).
    for (int k = 0; k < K1; ++k) {
        double r_k = in.x1_ref[k](0), v_k = in.x1_ref[k](3);
        double d1 = in.x1_ref[k](13), d2 = in.x1_ref[k](14), d3 = in.x1_ref[k](15), d4 = in.x1_ref[k](16);

        PathConstraintValues pcv = computePathConstraintValues(r_k, v_k, planet_config, spacecraft_config);
        PathConstraintGradients pcg = computePathConstraintGradients(r_k, v_k, planet_config, spacecraft_config);
        FullLoadFactorGradient ng = computeFullLoadFactorGradient(r_k, v_k, d1, d2, d3, d4,
                                                                     planet_config, spacecraft_config);

        constexpr double kMaxHeatFlux = 5.0e6, kMaxQbar = 2.0e5, kMaxLoadFactor = 30.0;

        addRow({{layout.dx1_offset(k) + 0, -pcg.dQdot_dr}, {layout.dx1_offset(k) + 3, -pcg.dQdot_dV}},
               kMaxHeatFlux - pcv.qdot_w_m2);
        cones.push_back(clarabel::NonnegativeConeT<double>(1));
        addRow({{layout.dx1_offset(k) + 0, -pcg.dqbar_dr}, {layout.dx1_offset(k) + 3, -pcg.dqbar_dV}},
               kMaxQbar - pcv.qbar_pa);
        cones.push_back(clarabel::NonnegativeConeT<double>(1));
        addRow({{layout.dx1_offset(k) + 0, -ng.dn_dr}, {layout.dx1_offset(k) + 3, -ng.dn_dV},
                {layout.dx1_offset(k) + 13, -ng.dn_ddelta[0]}, {layout.dx1_offset(k) + 14, -ng.dn_ddelta[1]},
                {layout.dx1_offset(k) + 15, -ng.dn_ddelta[2]}, {layout.dx1_offset(k) + 16, -ng.dn_ddelta[3]}},
               kMaxLoadFactor - pcv.n_g);
        cones.push_back(clarabel::NonnegativeConeT<double>(1));
    }

    // Box constraints: flap deflection/rate (state, every node), motor
    // torque (control, every interval).
    for (int k = 0; k < K1; ++k) {
        for (int i = 0; i < 4; ++i) {
            for (const auto& g : buildSymmetricBoxRows(layout.dx1_offset(k) + 13 + i, in.x1_ref[k](13 + i),
                                                          flap_config.delta_max_rad)) {
                addConeRowGroup(g);
            }
            for (const auto& g : buildSymmetricBoxRows(layout.dx1_offset(k) + 17 + i, in.x1_ref[k](17 + i),
                                                          flap_config.delta_dot_max_rad_s)) {
                addConeRowGroup(g);
            }
        }
    }
    for (int k = 0; k < K1 - 1; ++k) {
        for (int i = 0; i < 4; ++i) {
            for (const auto& g : buildSymmetricBoxRows(layout.du1_offset(k) + i, in.u1_ref[k](i),
                                                          flap_config.tau_m_max_n_m)) {
                addConeRowGroup(g);
            }
        }
    }

    // Trust region + virtual-control non-negativity + cost (Stage A's exact
    // per-component-scaled pattern, generalized to 21 states/4 controls).
    for (int k = 0; k < K1; ++k) {
        for (int i = 0; i < 21; ++i) {
            addRow({{layout.dx1_offset(k) + i, 1.0}, {layout.eta1_offset(k), -kStateScale1[i]}}, 0.0);
            cones.push_back(clarabel::NonnegativeConeT<double>(1));
            addRow({{layout.dx1_offset(k) + i, -1.0}, {layout.eta1_offset(k), -kStateScale1[i]}}, 0.0);
            cones.push_back(clarabel::NonnegativeConeT<double>(1));
        }
        if (k < K1 - 1) {
            for (int c = 0; c < 4; ++c) {
                addRow({{layout.du1_offset(k) + c, 1.0}, {layout.eta1_offset(k), -kControlScale1[c]}}, 0.0);
                cones.push_back(clarabel::NonnegativeConeT<double>(1));
                addRow({{layout.du1_offset(k) + c, -1.0}, {layout.eta1_offset(k), -kControlScale1[c]}}, 0.0);
                cones.push_back(clarabel::NonnegativeConeT<double>(1));
            }
        }
        addRow({{layout.eta1_offset(k), -1.0}}, 0.0);
        cones.push_back(clarabel::NonnegativeConeT<double>(1));
    }
    for (int k = 0; k < K1 - 1; ++k) {
        for (int i = 0; i < 21; ++i) {
            addRow({{layout.nu1_plus_offset(k) + i, -1.0}}, 0.0);
            cones.push_back(clarabel::NonnegativeConeT<double>(1));
            addRow({{layout.nu1_minus_offset(k) + i, -1.0}}, 0.0);
            cones.push_back(clarabel::NonnegativeConeT<double>(1));
        }
    }
    for (int k = 0; k < K1 - 1; ++k) {
        for (int i = 0; i < 21; ++i) {
            q(layout.nu1_plus_offset(k) + i) = cfg.w_nu / kStateScale1[i];
            q(layout.nu1_minus_offset(k) + i) = cfg.w_nu / kStateScale1[i];
        }
    }
    for (int k = 0; k < K1; ++k) q(layout.eta1_offset(k)) = cfg.w_eta;

    // Hard trust-region cap (adaptive accept/reject)
    for (int k = 1; k < K1; ++k) {
        for (int i = 0; i < 21; ++i) {
            addRow({{layout.dx1_offset(k) + i, 1.0}}, kStateScale1[i] * in.eta_max);
            cones.push_back(clarabel::NonnegativeConeT<double>(1));
            addRow({{layout.dx1_offset(k) + i, -1.0}}, kStateScale1[i] * in.eta_max);
            cones.push_back(clarabel::NonnegativeConeT<double>(1));
        }
    }
    for (int k = 0; k < K1 - 1; ++k) {
        for (int c = 0; c < 4; ++c) {
            addRow({{layout.du1_offset(k) + c, 1.0}}, kControlScale1[c] * in.eta_max);
            cones.push_back(clarabel::NonnegativeConeT<double>(1));
            addRow({{layout.du1_offset(k) + c, -1.0}}, kControlScale1[c] * in.eta_max);
            cones.push_back(clarabel::NonnegativeConeT<double>(1));
        }
    }

    // Phase 2

    // NOTE: Phase 2 has NO independent "delta_x2_0 = 0" boundary block here
    // (unlike Phase 1's fixed x1_initial) as Phase 2's initial state is
    // determined by the Phase 1 -> Phase 2 transition, not an externally
    // fixed constant. The transition equality rows assembled further below
    // (position/velocity via a linearized frame conversion, attitude/rate/
    // flap/mass via direct carry-over) fully constrain all 22 dx2_offset(0)
    // components exactly once

    std::vector<double> sigma_u_ref_per_interval(K2 - 1);
    std::vector<Eigen::Vector3d> u_vec_ref_per_interval(K2 - 1);

    for (int k = 0; k < K2 - 1; ++k) {
        EomFn f = [&](const Eigen::VectorXd& x, const Eigen::VectorXd& u) {
            return phase2Eom(x, u, in.deltaE_ref[k], in.phiE_ref[k], cfg.Isp_s,
                              planet_config, spacecraft_config, flap_config, gimbal_config);
        };
        EomJacobian jac = computeEomJacobianFd(f, in.x2_ref[k], in.u2_ref[k]);  // 22x22, 22x7

        const double m_ref = in.x2_ref[k](0);
        const double T_ref = in.u2_ref[k](0);
        Eigen::Quaterniond q_ref(in.x2_ref[k](10), in.x2_ref[k](7), in.x2_ref[k](8), in.x2_ref[k](9));
        Eigen::Vector3d thrust_dir_ref = thrustDirBody(in.deltaE_ref[k], in.phiE_ref[k]);
        Eigen::Vector3d u_vec_ref = (T_ref / m_ref) * (q_ref * thrust_dir_ref);
        double sigma_u_ref = u_vec_ref.norm();
        if (sigma_u_ref < 1e-9) sigma_u_ref = 1e-9;  // guard: near-zero-thrust reference
        sigma_u_ref_per_interval[k] = sigma_u_ref;
        u_vec_ref_per_interval[k] = u_vec_ref;

        Eigen::MatrixXd Ac2 = Eigen::MatrixXd::Zero(22, 22);
        Eigen::MatrixXd Bc2 = Eigen::MatrixXd::Zero(22, layout.kNu2);  // [uvec(3), sigma_u(1), rest(6)]

        // Row 0 (zdot): Stage B's exact closed form, not FD-derived.
        Bc2(0, 3) = -alpha_mdot;

        // Rows 1..21: chain-rule-reinterpreted FD Jacobian (see top-of-file comment).
        for (int i = 1; i < 22; ++i) {
            for (int j = 1; j < 22; ++j) Ac2(i, j) = jac.Ac(i, j);
            Ac2(i, 0) = jac.Ac(i, 0) * m_ref + jac.Bc(i, 0) * sigma_u_ref * m_ref;
            for (int c = 0; c < 3; ++c) {
                Bc2(i, c) = jac.Bc(i, 0) * m_ref * u_vec_ref(c) / sigma_u_ref;
            }
            Bc2(i, 4) = jac.Bc(i, 1);  // deltaE_dot
            Bc2(i, 5) = jac.Bc(i, 2);  // phiE_dot
            for (int c = 0; c < 4; ++c) Bc2(i, 6 + c) = jac.Bc(i, 3 + c);  // tau_m1..4
        }

        Eigen::MatrixXd A2 = Eigen::MatrixXd::Identity(22, 22) + dt2 * Ac2;
        Eigen::MatrixXd B2 = dt2 * Bc2;

        for (int i = 0; i < 22; ++i) {
            std::vector<std::pair<int, double>> coeffs;
            coeffs.emplace_back(layout.dx2_offset(k + 1) + i, 1.0);
            for (int j = 0; j < 22; ++j) {
                if (A2(i, j) != 0.0) coeffs.emplace_back(layout.dx2_offset(k) + j, -A2(i, j));
            }
            for (int c = 0; c < 3; ++c) {
                if (B2(i, c) != 0.0) coeffs.emplace_back(layout.du2_uvec_offset(k) + c, -B2(i, c));
            }
            if (B2(i, 3) != 0.0) coeffs.emplace_back(layout.du2_sigma_offset(k), -B2(i, 3));
            for (int c = 0; c < 6; ++c) {
                if (B2(i, 4 + c) != 0.0) coeffs.emplace_back(layout.du2_rest_offset(k) + c, -B2(i, 4 + c));
            }
            coeffs.emplace_back(layout.nu2_plus_offset(k) + i, -1.0);
            coeffs.emplace_back(layout.nu2_minus_offset(k) + i, 1.0);
            addRow(coeffs, 0.0);
            cones.push_back(clarabel::ZeroConeT<double>(1));
        }
    }

    // Terminal boundary: position/velocity matched exactly; attitude via a
    // linearized quaternion-error-angle bound; mass free (maximized via cost).
    {
        int kf = K2 - 1;
        for (int i = 1; i <= 6; ++i) {  // rx,ry,rz,vx,vy,vz
            addRow({{layout.dx2_offset(kf) + i, 1.0}}, in.x2_terminal_target(i) - in.x2_ref[kf](i));
            cones.push_back(clarabel::ZeroConeT<double>(1));
        }

        Eigen::Quaterniond q_target(in.x2_terminal_target(10), in.x2_terminal_target(7),
                                     in.x2_terminal_target(8), in.x2_terminal_target(9));
        Eigen::VectorXd x_ref_kf = in.x2_ref[kf];
        double err_ref = quatErrorAngle(x_ref_kf(7), x_ref_kf(8), x_ref_kf(9), x_ref_kf(10), q_target);
        std::array<double, 4> grad{};
        for (int c = 0; c < 4; ++c) {
            double h = 1e-6;
            Eigen::VectorXd xp = x_ref_kf, xm = x_ref_kf;
            xp(7 + c) += h; xm(7 + c) -= h;
            double ep = quatErrorAngle(xp(7), xp(8), xp(9), xp(10), q_target);
            double em = quatErrorAngle(xm(7), xm(8), xm(9), xm(10), q_target);
            grad[c] = (ep - em) / (2.0 * h);
        }
        addRow({{layout.dx2_offset(kf) + 7, -grad[0]}, {layout.dx2_offset(kf) + 8, -grad[1]},
                {layout.dx2_offset(kf) + 9, -grad[2]}, {layout.dx2_offset(kf) + 10, -grad[3]}},
               cfg.terminal_attitude_error_max_rad - err_ref);
        cones.push_back(clarabel::NonnegativeConeT<double>(1));
    }

    // Thrust-magnitude lossless convexification (Milestone 7.1) + mass-
    // lower-bound linear relaxation (Stage B's exact formulas, re-derived
    // here for delta_z/delta_sigma_u decision variables
    for (int k = 0; k < K2 - 1; ++k) {
        addConeRowGroup(buildThrustMagnitudeSocRows(layout.du2_sigma_offset(k), layout.du2_uvec_offset(k) + 0,
                                                       layout.du2_uvec_offset(k) + 1, layout.du2_uvec_offset(k) + 2,
                                                       sigma_u_ref_per_interval[k], u_vec_ref_per_interval[k]));

        double t_k = k * dt2;
        double m0_k = in.x2_ref[0](0) - (cfg.Tmax_N / (cfg.Isp_s * g0)) * t_k;
        double z0_k = std::log(m0_k);
        double n0_k = cfg.Tmin_N / m0_k;
        double n1_k = cfg.Tmax_N / m0_k;
        double z_ref_k = std::log(in.x2_ref[k](0));
        double sigma_u_ref_k = sigma_u_ref_per_interval[k];

        // Absolute constraint: n0_k*(1-(z-z0_k)) <= sigma_u <= n1_k*(1-(z-z0_k)),
        // with z=z_ref_k+delta_z, sigma_u=sigma_u_ref_k+delta_sigma_u.
        // Lower: -delta_sigma_u - n0_k*delta_z <= -n0_k*(1+z0_k) + n0_k*z_ref_k + sigma_u_ref_k
        addRow({{layout.du2_sigma_offset(k), -1.0}, {layout.dx2_offset(k) + 0, -n0_k}},
               -n0_k * (1.0 + z0_k) + n0_k * z_ref_k + sigma_u_ref_k);
        cones.push_back(clarabel::NonnegativeConeT<double>(1));
        // Upper: n1_k*delta_z + delta_sigma_u <= n1_k*(1+z0_k) - n1_k*z_ref_k - sigma_u_ref_k
        addRow({{layout.dx2_offset(k) + 0, n1_k}, {layout.du2_sigma_offset(k), 1.0}},
               n1_k * (1.0 + z0_k) - n1_k * z_ref_k - sigma_u_ref_k);
        cones.push_back(clarabel::NonnegativeConeT<double>(1));

        // Gimbal cone, axis at the reference attitude.
        Eigen::Quaterniond q_ref_k(in.x2_ref[k](10), in.x2_ref[k](7), in.x2_ref[k](8), in.x2_ref[k](9));
        Eigen::Vector3d e1_ref = q_ref_k * Eigen::Vector3d(-1.0, 0.0, 0.0);
        addConeRowGroup(buildGimbalConeSocRows(layout.du2_uvec_offset(k) + 0, layout.du2_uvec_offset(k) + 1,
                                                  layout.du2_uvec_offset(k) + 2, e1_ref,
                                                  gimbal_config.deltaE_max_rad, u_vec_ref_per_interval[k]));
    }

    // Glideslope SOC (Stage B's exact 3-row construction), skipping k=0.
    {
        double c_slope = 1.0 / std::tan(cfg.glideslope_deg * kPi / 180.0);
        double target_rx = in.x2_terminal_target(1), target_ry = in.x2_terminal_target(2),
               target_rz = in.x2_terminal_target(3);
        for (int k = 1; k < K2; ++k) {
            ConeRowGroup g;
            g.rows.push_back({{{layout.dx2_offset(k) + 3, -c_slope}},
                               -c_slope * (target_rz - in.x2_ref[k](3))});
            g.rows.push_back({{{layout.dx2_offset(k) + 1, -1.0}}, -(target_rx - in.x2_ref[k](1))});
            g.rows.push_back({{{layout.dx2_offset(k) + 2, -1.0}}, -(target_ry - in.x2_ref[k](2))});
            g.cone = clarabel::SecondOrderConeT<double>(3);
            addConeRowGroup(g);
        }
    }

    // Box constraints: flap deflection/rate (state), gimbal-rate + motor
    // torque (control).
    for (int k = 0; k < K2; ++k) {
        for (int i = 0; i < 4; ++i) {
            for (const auto& g : buildSymmetricBoxRows(layout.dx2_offset(k) + 14 + i, in.x2_ref[k](14 + i),
                                                          flap_config.delta_max_rad)) {
                addConeRowGroup(g);
            }
            for (const auto& g : buildSymmetricBoxRows(layout.dx2_offset(k) + 18 + i, in.x2_ref[k](18 + i),
                                                          flap_config.delta_dot_max_rad_s)) {
                addConeRowGroup(g);
            }
        }
    }
    for (int k = 0; k < K2 - 1; ++k) {
        for (const auto& g : buildSymmetricBoxRows(layout.du2_rest_offset(k) + 0, in.u2_ref[k](1),
                                                      gimbal_config.gimbal_rate_max_rad_s)) {
            addConeRowGroup(g);
        }
        for (const auto& g : buildSymmetricBoxRows(layout.du2_rest_offset(k) + 1, in.u2_ref[k](2),
                                                      gimbal_config.gimbal_rate_max_rad_s)) {
            addConeRowGroup(g);
        }
        for (int i = 0; i < 4; ++i) {
            for (const auto& g : buildSymmetricBoxRows(layout.du2_rest_offset(k) + 2 + i, in.u2_ref[k](3 + i),
                                                          flap_config.tau_m_max_n_m)) {
                addConeRowGroup(g);
            }
        }
    }

    // Trust region + virtual-control non-negativity + cost.
    for (int k = 0; k < K2; ++k) {
        for (int i = 0; i < 22; ++i) {
            addRow({{layout.dx2_offset(k) + i, 1.0}, {layout.eta2_offset(k), -kStateScale2[i]}}, 0.0);
            cones.push_back(clarabel::NonnegativeConeT<double>(1));
            addRow({{layout.dx2_offset(k) + i, -1.0}, {layout.eta2_offset(k), -kStateScale2[i]}}, 0.0);
            cones.push_back(clarabel::NonnegativeConeT<double>(1));
        }
        if (k < K2 - 1) {
            for (int c = 0; c < 3; ++c) {
                addRow({{layout.du2_uvec_offset(k) + c, 1.0}, {layout.eta2_offset(k), -kControlScaleVec2}}, 0.0);
                cones.push_back(clarabel::NonnegativeConeT<double>(1));
                addRow({{layout.du2_uvec_offset(k) + c, -1.0}, {layout.eta2_offset(k), -kControlScaleVec2}}, 0.0);
                cones.push_back(clarabel::NonnegativeConeT<double>(1));
            }
            addRow({{layout.du2_sigma_offset(k), 1.0}, {layout.eta2_offset(k), -kControlScaleSigma2}}, 0.0);
            cones.push_back(clarabel::NonnegativeConeT<double>(1));
            addRow({{layout.du2_sigma_offset(k), -1.0}, {layout.eta2_offset(k), -kControlScaleSigma2}}, 0.0);
            cones.push_back(clarabel::NonnegativeConeT<double>(1));
            for (int c = 0; c < 6; ++c) {
                addRow({{layout.du2_rest_offset(k) + c, 1.0}, {layout.eta2_offset(k), -kControlScaleRest2[c]}}, 0.0);
                cones.push_back(clarabel::NonnegativeConeT<double>(1));
                addRow({{layout.du2_rest_offset(k) + c, -1.0}, {layout.eta2_offset(k), -kControlScaleRest2[c]}}, 0.0);
                cones.push_back(clarabel::NonnegativeConeT<double>(1));
            }
        }
        addRow({{layout.eta2_offset(k), -1.0}}, 0.0);
        cones.push_back(clarabel::NonnegativeConeT<double>(1));
    }
    for (int k = 0; k < K2 - 1; ++k) {
        for (int i = 0; i < 22; ++i) {
            addRow({{layout.nu2_plus_offset(k) + i, -1.0}}, 0.0);
            cones.push_back(clarabel::NonnegativeConeT<double>(1));
            addRow({{layout.nu2_minus_offset(k) + i, -1.0}}, 0.0);
            cones.push_back(clarabel::NonnegativeConeT<double>(1));
        }
    }
    for (int k = 0; k < K2 - 1; ++k) {
        for (int i = 0; i < 22; ++i) {
            q(layout.nu2_plus_offset(k) + i) = cfg.w_nu / kStateScale2[i];
            q(layout.nu2_minus_offset(k) + i) = cfg.w_nu / kStateScale2[i];
        }
    }
    for (int k = 0; k < K2; ++k) q(layout.eta2_offset(k)) = cfg.w_eta;

    // Hard trust-region cap
    for (int k = 1; k < K2; ++k) {
        bool is_terminal = (k == K2 - 1);
        for (int i = 0; i < 22; ++i) {
            if (is_terminal && i >= 1 && i <= 6) continue;
            addRow({{layout.dx2_offset(k) + i, 1.0}}, kStateScale2[i] * in.eta_max);
            cones.push_back(clarabel::NonnegativeConeT<double>(1));
            addRow({{layout.dx2_offset(k) + i, -1.0}}, kStateScale2[i] * in.eta_max);
            cones.push_back(clarabel::NonnegativeConeT<double>(1));
        }
    }
    for (int k = 0; k < K2 - 1; ++k) {
        for (int c = 0; c < 3; ++c) {
            addRow({{layout.du2_uvec_offset(k) + c, 1.0}}, kControlScaleVec2 * in.eta_max);
            cones.push_back(clarabel::NonnegativeConeT<double>(1));
            addRow({{layout.du2_uvec_offset(k) + c, -1.0}}, kControlScaleVec2 * in.eta_max);
            cones.push_back(clarabel::NonnegativeConeT<double>(1));
        }
        addRow({{layout.du2_sigma_offset(k), 1.0}}, kControlScaleSigma2 * in.eta_max);
        cones.push_back(clarabel::NonnegativeConeT<double>(1));
        addRow({{layout.du2_sigma_offset(k), -1.0}}, kControlScaleSigma2 * in.eta_max);
        cones.push_back(clarabel::NonnegativeConeT<double>(1));
        for (int c = 0; c < 6; ++c) {
            addRow({{layout.du2_rest_offset(k) + c, 1.0}}, kControlScaleRest2[c] * in.eta_max);
            cones.push_back(clarabel::NonnegativeConeT<double>(1));
            addRow({{layout.du2_rest_offset(k) + c, -1.0}}, kControlScaleRest2[c] * in.eta_max);
            cones.push_back(clarabel::NonnegativeConeT<double>(1));
        }
    }

    // Phase transition (Milestone 8) + qbar-relight (Section 5
    {
        int k1f = K1 - 1;
        const Eigen::VectorXd& x1f_ref = in.x1_ref[k1f];
        Phase2TranslationalState conv = convertPhase1ToPhase2(x1f_ref(0), x1f_ref(1), x1f_ref(2), x1f_ref(3),
                                                                  x1f_ref(4), x1f_ref(5), in.transition_frame,
                                                                  planet_config);
        Eigen::MatrixXd J = computeTransitionJacobian(x1f_ref(0), x1f_ref(1), x1f_ref(2), x1f_ref(3),
                                                        x1f_ref(4), x1f_ref(5), in.transition_frame, planet_config);
        std::array<double, 6> conv_vals = {conv.rx, conv.ry, conv.rz, conv.vx, conv.vy, conv.vz};
        std::array<int, 6> x2_idx = {1, 2, 3, 4, 5, 6};  // rx,ry,rz,vx,vy,vz in dx2

        for (int i = 0; i < 6; ++i) {
            std::vector<std::pair<int, double>> coeffs;
            coeffs.emplace_back(layout.dx2_offset(0) + x2_idx[i], 1.0);
            for (int j = 0; j < 6; ++j) {
                if (J(i, j) != 0.0) coeffs.emplace_back(layout.dx1_offset(k1f) + j, -J(i, j));
            }
            addRow(coeffs, conv_vals[i] - in.x2_ref[0](x2_idx[i]));
            cones.push_back(clarabel::ZeroConeT<double>(1));
        }

        // Attitude/rate/flap state carry over unchanged (same body frame convention both phases)
        // x1 indices 6..20 (q1-4,wx-wz,d1-4,ddot1-4) map to x2 indices 7..21.
        for (int i = 0; i < 15; ++i) {
            int x1_idx = 6 + i, x2_idx_i = 7 + i;
            addRow({{layout.dx2_offset(0) + x2_idx_i, 1.0}, {layout.dx1_offset(k1f) + x1_idx, -1.0}},
                   in.x1_ref[k1f](x1_idx) - in.x2_ref[0](x2_idx_i));
            cones.push_back(clarabel::ZeroConeT<double>(1));
        }

        // Mass at handoff: fixed wet-mass config constant (Phase 1 has no
        // mass state), expressed in the z=ln(m) slot.
        double z_target = std::log(cfg.m_wet_at_handoff_kg);
        double z0_ref = std::log(in.x2_ref[0](0));
        addRow({{layout.dx2_offset(0) + 0, 1.0}}, z_target - z0_ref);
        cones.push_back(clarabel::ZeroConeT<double>(1));

        // qbar-relight inequality (Section 5), Stage A's exact closed form
        // reused at the transition node.
        double r_k1f = x1f_ref(0), v_k1f = x1f_ref(3);
        PathConstraintValues pcv = computePathConstraintValues(r_k1f, v_k1f, planet_config, spacecraft_config);
        PathConstraintGradients pcg = computePathConstraintGradients(r_k1f, v_k1f, planet_config, spacecraft_config);
        addRow({{layout.dx1_offset(k1f) + 0, -pcg.dqbar_dr}, {layout.dx1_offset(k1f) + 3, -pcg.dqbar_dV}},
               cfg.qbar_relight_max_pa - pcv.qbar_pa);
        cones.push_back(clarabel::NonnegativeConeT<double>(1));
    }

    Eigen::SparseMatrix<double> A(row, n);
    A.setFromTriplets(A_triplets.begin(), A_triplets.end());
    Eigen::VectorXd b = Eigen::Map<Eigen::VectorXd>(b_vec.data(), b_vec.size());
    Eigen::SparseMatrix<double> P(n, n);  // pure L1/linear cost, no quadratic term

    // Objective: maximize final ln-mass (Stage B's exact pattern) plus the
    // per-phase nu/eta penalty terms already populated into q above.
    q(layout.dx2_offset(K2 - 1) + 0) += -1.0;

    return SocpProblem{P, q, A, b, cones};
}

}  // namespace guidance_scp
