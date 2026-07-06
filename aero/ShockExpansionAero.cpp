#include "ShockExpansionAero.h"

#include <algorithm>
#include <cmath>

namespace aero_model {

namespace {

// MSVC only defines M_PI when _USE_MATH_DEFINES is set before the first
// <cmath>/<math.h> include anywhere in the translation unit -- use an
// explicit local constant instead (same fix as elsewhere in this repo).
constexpr double kPi = 3.14159265358979323846;

constexpr double kNewtonTol = 1e-10;
constexpr int kNewtonMaxIters = 50;

// Residual of the theta-beta-Mach relation (Anderson eq. 4.23), rearranged
// to f(beta) = 0:
//   tan(theta) * [M1^2*(gamma+cos(2*beta)) + 2] = 2*cot(beta)*(M1^2*sin^2(beta) - 1)
double thetaBetaMachResidual(double beta, double mach_inf, double theta_rad, double gamma) {
    const double M2 = mach_inf * mach_inf;
    const double sinb = std::sin(beta);
    const double cosb = std::cos(beta);
    const double cotb = cosb / sinb;
    return std::tan(theta_rad) * (M2 * (gamma + std::cos(2.0 * beta)) + 2.0)
         - 2.0 * cotb * (M2 * sinb * sinb - 1.0);
}

}  // namespace

double prandtlMeyerNu(double mach, double gamma) {
    double M2m1 = mach * mach - 1.0;
    if (M2m1 < 0.0) M2m1 = 0.0;  // guard against mach slightly < 1 from floating point
    const double k = (gamma + 1.0) / (gamma - 1.0);
    return std::sqrt(k) * std::atan(std::sqrt(M2m1 / k)) - std::atan(std::sqrt(M2m1));
}

double prandtlMeyerNuMax(double gamma) {
    const double k = (gamma + 1.0) / (gamma - 1.0);
    return (kPi / 2.0) * (std::sqrt(k) - 1.0);
}

double invertPrandtlMeyer(double nu_target, double gamma) {
    nu_target = std::min(nu_target, prandtlMeyerNuMax(gamma) - 1e-9);
    double M = 1.5;  // generic initial guess; Newton below converges quickly
                      // for the Mach range this table actually uses (~1-15)
    for (int iter = 0; iter < kNewtonMaxIters; ++iter) {
        const double f = prandtlMeyerNu(M, gamma) - nu_target;
        if (std::abs(f) < kNewtonTol) break;
        double M2m1 = M * M - 1.0;
        if (M2m1 < 0.0) M2m1 = 0.0;
        const double dnu_dM = std::sqrt(M2m1) / (M * (1.0 + 0.5 * (gamma - 1.0) * M * M));
        if (dnu_dM < 1e-12) {
            M += 0.05;  // derivative vanishes right at M=1 -- nudge forward
            continue;
        }
        double M_next = M - f / dnu_dM;
        if (M_next < 1.0) M_next = 1.0 + 0.5 * (M - 1.0);  // stay in the valid domain
        M = M_next;
    }
    return M;
}

std::optional<double> solveWeakObliqueShockBeta(double mach_inf, double theta_rad, double gamma) {
    if (mach_inf <= 1.0) return std::nullopt;  // undefined for subsonic upstream flow

    const double mach_angle = std::asin(std::min(1.0, 1.0 / mach_inf));
    const double beta_lo = mach_angle + 1e-6;
    const double beta_hi = kPi / 2.0 - 1e-6;
    if (beta_lo >= beta_hi) return std::nullopt;

    double beta = std::min(beta_hi, mach_angle + 0.05);  // seed just above the Mach angle
    bool converged = false;
    for (int iter = 0; iter < kNewtonMaxIters; ++iter) {
        const double f = thetaBetaMachResidual(beta, mach_inf, theta_rad, gamma);
        if (std::abs(f) < kNewtonTol) {
            converged = true;
            break;
        }
        const double h = 1e-6;
        const double beta_p = std::clamp(beta + h, beta_lo, beta_hi);
        const double beta_m = std::clamp(beta - h, beta_lo, beta_hi);
        const double fp = thetaBetaMachResidual(beta_p, mach_inf, theta_rad, gamma);
        const double fm = thetaBetaMachResidual(beta_m, mach_inf, theta_rad, gamma);
        const double dfdb = (fp - fm) / (beta_p - beta_m);
        if (std::abs(dfdb) < 1e-12) break;
        beta = std::clamp(beta - f / dfdb, beta_lo, beta_hi);
    }
    if (!converged) {
        const double f_final = thetaBetaMachResidual(beta, mach_inf, theta_rad, gamma);
        if (std::abs(f_final) < 1e-6) converged = true;
    }
    if (!converged || beta <= beta_lo || beta >= beta_hi) {
        return std::nullopt;  // detached
    }
    return beta;
}

Eigen::Vector3d ShockExpansionAeroModel::freestreamDirectionBody(double alpha_rad, double beta_rad) const {
    // Identical to NewtonianAeroModel::freestreamDirectionBody -- 3 lines,
    // not worth sharing across a new abstraction.
    return Eigen::Vector3d(-std::cos(alpha_rad) * std::cos(beta_rad),
                            std::sin(beta_rad),
                            std::sin(alpha_rad) * std::cos(beta_rad));
}

double ShockExpansionAeroModel::panelCp(double mach_inf, double theta_rad) const {
    if (theta_rad >= 0.0) {
        const auto beta = solveWeakObliqueShockBeta(mach_inf, theta_rad, gamma_);
        if (!beta) {
            // Detached shock -- fall back to local Newtonian impact theory
            // for this one panel, a bounded degradation that never produces
            // a NaN/undefined Cp.
            const NewtonianAeroModel local_newton(gamma_);
            const double sin_theta = std::sin(theta_rad);
            return local_newton.cpMax(mach_inf) * sin_theta * sin_theta;
        }
        const double Mn1 = mach_inf * std::sin(*beta);
        return (4.0 / (gamma_ + 1.0)) * (Mn1 * Mn1 - 1.0) / (mach_inf * mach_inf);
    } else {
        const double nu1 = prandtlMeyerNu(mach_inf, gamma_);
        const double nu_max = prandtlMeyerNuMax(gamma_);
        const double nu2 = std::min(nu1 + (-theta_rad), nu_max - 1e-9);
        const double M2 = invertPrandtlMeyer(nu2, gamma_);
        const double g = gamma_;
        const double p_ratio = std::pow(
            (1.0 + 0.5 * (g - 1.0) * mach_inf * mach_inf) / (1.0 + 0.5 * (g - 1.0) * M2 * M2),
            g / (g - 1.0));
        return (p_ratio - 1.0) / (0.5 * g * mach_inf * mach_inf);
    }
}

AeroCoefficients ShockExpansionAeroModel::evaluate(
    const PanelMesh& mesh,
    const std::unordered_map<int, double>& flap_deflections_rad, double alpha_rad,
    double beta_rad, double mach_inf, const Eigen::Vector3d& moment_ref, double S_ref,
    double L_ref) const {
    const Eigen::Vector3d V_hat = freestreamDirectionBody(alpha_rad, beta_rad);
    const std::vector<Panel> panels = mesh.deflected(flap_deflections_rad);

    Eigen::Vector3d sum_F_over_q = Eigen::Vector3d::Zero();
    Eigen::Vector3d sum_M_over_q = Eigen::Vector3d::Zero();

    for (const Panel& p : panels) {
        // Signed local flow-deflection angle: theta>0 compressive, theta<0
        // expansive. Unlike NewtonianAeroModel, nothing is skipped here.
        const double sin_theta = std::clamp(-V_hat.dot(p.normal), -1.0, 1.0);
        const double theta = std::asin(sin_theta);
        const double Cp = panelCp(mach_inf, theta);

        const Eigen::Vector3d dF_over_q = -Cp * p.area * p.normal;
        sum_F_over_q += dF_over_q;

        const Eigen::Vector3d r = p.centroid - moment_ref;
        sum_M_over_q += r.cross(dF_over_q);
    }

    AeroCoefficients c;
    const double CX = sum_F_over_q.x() / S_ref;
    const double CZ = sum_F_over_q.z() / S_ref;
    c.CA = -CX;
    c.CN = CZ;

    // Body-to-wind axis transformation. NOTE: alpha-only, an accepted
    // approximation at nonzero beta -- same documented limitation as
    // NewtonianAeroModel::evaluate.
    c.CL = c.CN * std::cos(alpha_rad) - c.CA * std::sin(alpha_rad);
    c.CD = c.CN * std::sin(alpha_rad) + c.CA * std::cos(alpha_rad);

    c.Cl_roll = sum_M_over_q.x() / (S_ref * L_ref);
    c.Cm      = sum_M_over_q.y() / (S_ref * L_ref);
    c.Cn_yaw  = sum_M_over_q.z() / (S_ref * L_ref);

    return c;
}

}  // namespace aero_model
