#include "full_scp_loop.h"
#include "full_loop_phase1.h"
#include "full_loop_phase2.h"
#include <Eigen/Geometry>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

namespace guidance_scp {

namespace {

std::vector<double> Resample(const std::vector<double>& src, int new_k) {
    std::vector<double> out(new_k);
    int n = static_cast<int>(src.size());
    for (int k = 0; k < new_k; ++k) {
        double frac = (new_k == 1) ? 0.0 : static_cast<double>(k) / (new_k - 1);
        double pos = frac * (n - 1);
        int lo = static_cast<int>(std::floor(pos));
        int hi = std::min(lo + 1, n - 1);
        double t = pos - lo;
        out[k] = src[lo] + t * (src[hi] - src[lo]);
    }
    return out;
}

Eigen::Vector3d ThrustDirBody(double deltaE, double phiE) {
    return Eigen::Vector3d(-std::cos(deltaE) * std::cos(phiE), std::sin(phiE),
                            std::sin(deltaE) * std::cos(phiE));
}

}

StageCResult stitchStageAAndB(const StageAResult& stage_a, const StageBResult& stage_b,
                                int K1, int K2,
                                const PlanetConfig& planet_config,
                                const SpacecraftConfig& spacecraft_config) {
    (void)spacecraft_config;  // not needed by this simplified bootstrap (see top-of-file comment)
    StageCResult out;
    const auto& ha = stage_a.history;
    const auto& hb = stage_b.history;

    auto r1 = Resample(ha.r, K1), la1 = Resample(ha.la, K1), lo1 = Resample(ha.lo, K1);
    auto v1 = Resample(ha.v, K1), fpa1 = Resample(ha.fpa, K1), vazi1 = Resample(ha.v_azi, K1);

    for (int k = 0; k < K1; ++k) {
        double fpa = fpa1[k], v_azi = vazi1[k];
        Eigen::Vector3d v_hat_local(std::sin(fpa), std::cos(fpa) * std::cos(v_azi), std::cos(fpa) * std::sin(v_azi));
        Eigen::Quaterniond q = Eigen::Quaterniond::FromTwoVectors(Eigen::Vector3d::UnitX(), v_hat_local);
        Eigen::VectorXd x(kPhase1StateDim);
        x << r1[k], la1[k], lo1[k], v1[k], fpa, v_azi,
             q.x(), q.y(), q.z(), q.w(),
             0.0, 0.0, 0.0,
             0.0, 0.0, 0.0, 0.0,
             0.0, 0.0, 0.0, 0.0;
        out.x1_ref.push_back(x);
    }
    for (int k = 0; k < K1 - 1; ++k) out.u1_ref.push_back(Eigen::VectorXd::Zero(kPhase1ControlDim));

    auto m2 = Resample(hb.m, K2), rx2 = Resample(hb.rx, K2), ry2 = Resample(hb.ry, K2), rz2 = Resample(hb.rz, K2);
    auto vx2 = Resample(hb.vx, K2), vy2 = Resample(hb.vy, K2), vz2 = Resample(hb.vz, K2);

    Eigen::Quaterniond q2 = Eigen::Quaterniond::FromTwoVectors(Eigen::Vector3d(-1.0, 0.0, 0.0),
                                                                 Eigen::Vector3d(0.0, 0.0, 1.0));
    for (int k = 0; k < K2; ++k) {
        Eigen::VectorXd x(kPhase2StateDim);
        x << m2[k], rx2[k], ry2[k], rz2[k], vx2[k], vy2[k], vz2[k],
             q2.x(), q2.y(), q2.z(), q2.w(),
             0.0, 0.0, 0.0,
             0.0, 0.0, 0.0, 0.0,
             0.0, 0.0, 0.0, 0.0;
        out.x2_ref.push_back(x);
        out.deltaE_ref.push_back(0.0);
        out.phiE_ref.push_back(0.0);
    }
    for (int k = 0; k < K2 - 1; ++k) {
        Eigen::VectorXd u(kPhase2ControlDim);
        u << m2[k] * planet_config.g_0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
        out.u2_ref.push_back(u);
    }

    return out;
}

FullLoopResult solveFullLoop(const StageCResult& initial_reference,
                              double t_scale_1, double t_scale_2,
                              const Eigen::VectorXd& x1_initial,
                              const Eigen::VectorXd& x2_terminal_target,
                              const Phase1ToPhase2Frame& transition_frame,
                              const PlanetConfig& planet_config,
                              const SpacecraftConfig& spacecraft_config,
                              const FlapActuatorConfig& flap_config,
                              const GimbalActuatorConfig& gimbal_config,
                              const FullLoopConfig& full_loop_config,
                              int max_iters, double eps_dyn, double eps_tr) {
    FullLoopResult result;
    result.x1_ref = initial_reference.x1_ref;
    result.u1_ref = initial_reference.u1_ref;
    result.x2_ref = initial_reference.x2_ref;
    result.u2_ref = initial_reference.u2_ref;
    result.deltaE_ref = initial_reference.deltaE_ref;
    result.phiE_ref = initial_reference.phiE_ref;

    const int K1 = static_cast<int>(result.x1_ref.size());
    const int K2 = static_cast<int>(result.x2_ref.size());
    FullLoopVarLayout layout{K1, K2};
    const double dt1 = t_scale_1 / (K1 - 1);
    const double dt2 = t_scale_2 / (K2 - 1);

    // Adaptive trust-region accept/reject
    double eta_max = full_loop_config.eta_max_init;
    int accepted_iters = 0;
    int total_attempts = 0;
    int consecutive_floor_rejects = 0;
    double best_max_nu = std::numeric_limits<double>::infinity();

    while (accepted_iters < max_iters && total_attempts < full_loop_config.max_solve_attempts) {
        ++total_attempts;

        FullLoopSubproblemInputs in;
        in.x1_ref = result.x1_ref;
        in.u1_ref = result.u1_ref;
        in.t_scale_1 = t_scale_1;
        in.x2_ref = result.x2_ref;
        in.u2_ref = result.u2_ref;
        in.deltaE_ref = result.deltaE_ref;
        in.phiE_ref = result.phiE_ref;
        in.t_scale_2 = t_scale_2;
        in.x1_initial = x1_initial;
        in.x2_terminal_target = x2_terminal_target;
        in.transition_frame = transition_frame;
        in.planet_config = planet_config;
        in.spacecraft_config = const_cast<SpacecraftConfig*>(&spacecraft_config);
        in.flap_config = flap_config;
        in.gimbal_config = gimbal_config;
        in.full_loop_config = full_loop_config;
        in.eta_max = eta_max;

        SocpProblem problem = buildFullLoopSubproblem(in);
        SocpSolution solution = solve(problem);

        bool non_solved = (solution.status != clarabel::SolverStatus::Solved
                            && solution.status != clarabel::SolverStatus::AlmostSolved);

        double max_nu = 0.0, max_eta = 0.0;
        double max_defect_nonlinear = std::numeric_limits<double>::infinity();
        std::vector<Eigen::VectorXd> x1_candidate, u1_candidate, x2_candidate, u2_candidate;
        std::vector<double> deltaE_candidate, phiE_candidate;
        bool accept = false;

        if (!non_solved) {
            // Build candidates from copies, does not mutate result.* until
            // the step is actually accepted below (lesson from Stage A: the
            // previous version mutated result.* unconditionally, which is
            // exactly what an accept/reject mechanism must not do).
            x1_candidate = result.x1_ref;
            u1_candidate = result.u1_ref;
            x2_candidate = result.x2_ref;
            u2_candidate = result.u2_ref;
            deltaE_candidate = result.deltaE_ref;
            phiE_candidate = result.phiE_ref;

            for (int k = 0; k < K1; ++k) {
                Eigen::VectorXd delta(21);
                for (int i = 0; i < 21; ++i) delta(i) = solution.x(layout.dx1_offset(k) + i);
                x1_candidate[k] += delta;
                max_eta = std::max(max_eta, solution.x(layout.eta1_offset(k)));
            }
            // Renormalize Phase 1 quaternions after the additive update.
            for (int k = 0; k < K1; ++k) {
                Eigen::Quaterniond q(x1_candidate[k](9), x1_candidate[k](6), x1_candidate[k](7), x1_candidate[k](8));
                q.normalize();
                x1_candidate[k](6) = q.x(); x1_candidate[k](7) = q.y();
                x1_candidate[k](8) = q.z(); x1_candidate[k](9) = q.w();
            }
            for (int k = 0; k < K1 - 1; ++k) {
                Eigen::VectorXd delta(4);
                for (int i = 0; i < 4; ++i) delta(i) = solution.x(layout.du1_offset(k) + i);
                u1_candidate[k] += delta;
                double nu_l1 = 0.0;
                for (int i = 0; i < 21; ++i) {
                    nu_l1 += solution.x(layout.nu1_plus_offset(k) + i) + solution.x(layout.nu1_minus_offset(k) + i);
                }
                max_nu = std::max(max_nu, nu_l1);
            }

            for (int k = 0; k < K2; ++k) {
                double z_ref_k = std::log(x2_candidate[k](0));
                double delta_z = solution.x(layout.dx2_offset(k) + 0);
                x2_candidate[k](0) = std::exp(z_ref_k + delta_z);
                for (int i = 1; i < 22; ++i) {
                    x2_candidate[k](i) += solution.x(layout.dx2_offset(k) + i);
                }
                max_eta = std::max(max_eta, solution.x(layout.eta2_offset(k)));
            }
            for (int k = 0; k < K2; ++k) {
                Eigen::Quaterniond q(x2_candidate[k](10), x2_candidate[k](7), x2_candidate[k](8), x2_candidate[k](9));
                q.normalize();
                x2_candidate[k](7) = q.x(); x2_candidate[k](8) = q.y();
                x2_candidate[k](9) = q.z(); x2_candidate[k](10) = q.w();
            }

            for (int k = 0; k < K2 - 1; ++k) {
                // Recompute this attempt's (pre-update) reference u_vec/
                // sigma_u exactly as buildFullLoopSubproblem did, using the
                // values that were actually fed into this solve (in.*).
                double m_ref = in.x2_ref[k](0);
                double T_ref = in.u2_ref[k](0);
                Eigen::Quaterniond q_ref(in.x2_ref[k](10), in.x2_ref[k](7), in.x2_ref[k](8), in.x2_ref[k](9));
                Eigen::Vector3d thrust_dir_ref = ThrustDirBody(in.deltaE_ref[k], in.phiE_ref[k]);
                Eigen::Vector3d u_vec_ref = (T_ref / m_ref) * (q_ref * thrust_dir_ref);

                Eigen::Vector3d delta_u_vec(solution.x(layout.du2_uvec_offset(k) + 0),
                                              solution.x(layout.du2_uvec_offset(k) + 1),
                                              solution.x(layout.du2_uvec_offset(k) + 2));
                Eigen::Vector3d u_vec_new = u_vec_ref + delta_u_vec;

                double m_new = x2_candidate[k](0);  // already updated above
                double T_new = m_new * u_vec_new.norm();

                Eigen::Quaterniond q_new(x2_candidate[k](10), x2_candidate[k](7),
                                           x2_candidate[k](8), x2_candidate[k](9));
                Eigen::Vector3d dir_body = q_new.conjugate() * u_vec_new.normalized();
                double phiE_new = std::asin(std::clamp(dir_body.y(), -1.0, 1.0));
                double deltaE_new = std::atan2(dir_body.z(), -dir_body.x());
                deltaE_candidate[k] = deltaE_new;
                phiE_candidate[k] = phiE_new;

                double deltaE_dot_new = in.u2_ref[k](1) + solution.x(layout.du2_rest_offset(k) + 0);
                double phiE_dot_new = in.u2_ref[k](2) + solution.x(layout.du2_rest_offset(k) + 1);
                Eigen::VectorXd u_new(7);
                u_new << T_new, deltaE_dot_new, phiE_dot_new,
                         in.u2_ref[k](3) + solution.x(layout.du2_rest_offset(k) + 2),
                         in.u2_ref[k](4) + solution.x(layout.du2_rest_offset(k) + 3),
                         in.u2_ref[k](5) + solution.x(layout.du2_rest_offset(k) + 4),
                         in.u2_ref[k](6) + solution.x(layout.du2_rest_offset(k) + 5);
                u2_candidate[k] = u_new;

                double nu_l1 = 0.0;
                for (int i = 0; i < 22; ++i) {
                    nu_l1 += solution.x(layout.nu2_plus_offset(k) + i) + solution.x(layout.nu2_minus_offset(k) + i);
                }
                max_nu = std::max(max_nu, nu_l1);
            }

            // Diagnostic true-nonlinear defect
            EomFn f1 = [&](const Eigen::VectorXd& x, const Eigen::VectorXd& u) {
                return phase1Eom(x, u, planet_config, spacecraft_config, flap_config);
            };
            double max_defect_1 = computeMaxNonlinearDefectGeneric(f1, x1_candidate, u1_candidate, dt1,
                                                                      kStateScale1.data(), 21);
            double max_defect_2 = 0.0;
            for (int k = 0; k + 1 < K2; ++k) {
                EomFn f2 = [&](const Eigen::VectorXd& x, const Eigen::VectorXd& u) {
                    return phase2Eom(x, u, deltaE_candidate[k], phiE_candidate[k], full_loop_config.Isp_s,
                                      planet_config, spacecraft_config, flap_config, gimbal_config);
                };
                std::vector<Eigen::VectorXd> pair_x = {x2_candidate[k], x2_candidate[k + 1]};
                std::vector<Eigen::VectorXd> pair_u = {u2_candidate[k]};
                max_defect_2 = std::max(max_defect_2,
                                          computeMaxNonlinearDefectGeneric(f2, pair_x, pair_u, dt2,
                                                                             kStateScale2.data(), 22));
            }
            max_defect_nonlinear = std::max(max_defect_1, max_defect_2);

            // Accept/reject: first-ever step unconditionally, thereafter
            // accept iff max||nu|| (combined across both phases) has not
            // regressed by more than nu_regression_tolerance relative to
            // the best value seen across previously-accepted iterations
            // NOT an absolute/eta_max-scaled defect threshold (Stage A's
            // lesson #4: a large first-correction defect/nu is structural
            // here, not a bad-step signal).
            accept = !std::isfinite(best_max_nu)
                     || (max_nu <= best_max_nu * full_loop_config.nu_regression_tolerance);
        }

        result.attempt_max_defect_nonlinear.push_back(max_defect_nonlinear);
        result.attempt_eta_max.push_back(eta_max);
        result.attempt_accepted.push_back(accept);

        if (accept) {
            result.x1_ref = x1_candidate;
            result.u1_ref = u1_candidate;
            result.x2_ref = x2_candidate;
            result.u2_ref = u2_candidate;
            result.deltaE_ref = deltaE_candidate;
            result.phiE_ref = phiE_candidate;
            ++accepted_iters;
            consecutive_floor_rejects = 0;
            best_max_nu = std::min(best_max_nu, max_nu);

            result.max_nu_per_iter.push_back(max_nu);
            result.max_eta_per_iter.push_back(max_eta);

            if (max_nu < eps_dyn && max_eta < eps_tr) {
                result.converged = true;
                break;
            }
        } else {
            // Reject (whether non-Solved or nu-regressed): shrink the hard
            // cap and retry the SAME outer iteration, never grow on a
            // rejection
            eta_max = std::max(eta_max * full_loop_config.eta_max_shrink, full_loop_config.eta_max_floor);
            if (eta_max <= full_loop_config.eta_max_floor * (1.0 + 1e-9)) {
                ++consecutive_floor_rejects;
            } else {
                consecutive_floor_rejects = 0;
            }
            if (consecutive_floor_rejects >= full_loop_config.max_consecutive_floor_rejects) {
                std::cerr << "Full loop: eta_max stuck at floor after " << consecutive_floor_rejects
                          << " consecutive rejections -- aborting honestly rather than spinning." << std::endl;
                result.converged = false;
                break;
            }
        }
    }

    return result;
}

}  // namespace guidance_scp
