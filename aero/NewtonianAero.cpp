#include "NewtonianAero.h"

#include <cmath>

namespace aero_model {

double NewtonianAeroModel::cpMax(double mach_inf) const {
    const double g = gamma_;
    const double M2 = mach_inf * mach_inf;

    // Modified-Newtonian Cp_max behind a normal shock (Anderson,
    // "Hypersonic and High-Temperature Gas Dynamics"), as used in
    // Olsen & Bettinger (2022) eq. 23. Valid for M > 1; guard against
    // the denominator vanishing/going negative near M -> 1.
    const double denom = 4.0 * g * M2 - 2.0 * (g - 1.0);
    if (mach_inf < 1.2 || denom <= 1e-9) {
        // Below this, modified-Newtonian theory is not physically meaningful
        // anyway (it is a hypersonic limit theory) -- clamp to the M=1.2
        // value so the sweep stays well-defined for CFD/Kriging blending
        // to correct in the transonic/subsonic regime.
        const double M2c = 1.2 * 1.2;
        const double denomc = 4.0 * g * M2c - 2.0 * (g - 1.0);
        const double term1 = std::pow((g + 1.0) * (g + 1.0) * M2c / denomc,
                                       g / (g - 1.0));
        const double term2 = (1.0 - g + 2.0 * g * M2c) / (g + 1.0);
        return (2.0 / (g * M2c)) * (term1 * term2 - 1.0);
    }

    const double term1 = std::pow((g + 1.0) * (g + 1.0) * M2 / denom, g / (g - 1.0));
    const double term2 = (1.0 - g + 2.0 * g * M2) / (g + 1.0);
    return (2.0 / (g * M2)) * (term1 * term2 - 1.0);
}

Eigen::Vector3d NewtonianAeroModel::freestreamDirectionBody(double alpha_rad, double beta_rad) const {
    // Direction the flow moves, expressed in body axes (x = nose direction,
    // y = pitch axis, z completes right-handed triad). At alpha=beta=0 the
    // flow moves in -x (nose-first flight through still air). Generalizes
    // the original alpha-only formula (-cos(alpha), 0, sin(alpha)) to include
    // sideslip beta; reduces to it exactly at beta=0.
    return Eigen::Vector3d(-std::cos(alpha_rad) * std::cos(beta_rad),
                            std::sin(beta_rad),
                            std::sin(alpha_rad) * std::cos(beta_rad));
}

AeroCoefficients NewtonianAeroModel::evaluate(
    const PanelMesh& mesh,
    const std::unordered_map<int, double>& flap_deflections_rad, double alpha_rad,
    double beta_rad, double mach_inf, const Eigen::Vector3d& moment_ref, double S_ref,
    double L_ref) const {
    const double Cp_max = cpMax(mach_inf);
    const Eigen::Vector3d V_hat = freestreamDirectionBody(alpha_rad, beta_rad);

    const std::vector<Panel> panels = mesh.deflected(flap_deflections_rad);

    Eigen::Vector3d sum_F_over_q = Eigen::Vector3d::Zero();
    Eigen::Vector3d sum_M_over_q = Eigen::Vector3d::Zero();

    for (const Panel& p : panels) {
        // sin(theta): incidence angle between local surface and the flow.
        // Derivation: at a stagnation point the outward normal is
        // antiparallel to the flow direction, giving sin(theta) = 1.
        const double sin_theta = -V_hat.dot(p.normal);
        if (sin_theta <= 0.0) {
            continue;  // leeward panel: Cp = 0 under Newtonian theory
        }
        const double Cp = Cp_max * sin_theta * sin_theta;

        // Pressure force acts into the body, i.e. along -normal.
        const Eigen::Vector3d dF_over_q = -Cp * p.area * p.normal;
        sum_F_over_q += dF_over_q;

        const Eigen::Vector3d r = p.centroid - moment_ref;
        sum_M_over_q += r.cross(dF_over_q);
    }

    AeroCoefficients c;
    const double CX = sum_F_over_q.x() / S_ref;
    const double CZ = sum_F_over_q.z() / S_ref;
    c.CA = -CX;  // positive CA = drag-like, opposing nose direction
    c.CN = CZ;

    // Body-to-wind axis transformation. NOTE: alpha-only, an accepted
    // approximation at nonzero beta -- see the header comment on evaluate().
    c.CL = c.CN * std::cos(alpha_rad) - c.CA * std::sin(alpha_rad);
    c.CD = c.CN * std::sin(alpha_rad) + c.CA * std::cos(alpha_rad);

    // Full 3-axis moment -- sum_M_over_q was already computed above in all
    // 3 components; previously only .y() (Cm) was exposed.
    c.Cl_roll = sum_M_over_q.x() / (S_ref * L_ref);
    c.Cm      = sum_M_over_q.y() / (S_ref * L_ref);
    c.Cn_yaw  = sum_M_over_q.z() / (S_ref * L_ref);

    return c;
}

}  // namespace aero_model
