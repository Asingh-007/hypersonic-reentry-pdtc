#include "full_loop_transition.h"
#include "FiniteDifferenceJacobian.h"
#include <cmath>

namespace guidance_scp {

Phase2TranslationalState convertPhase1ToPhase2(double r, double la, double lo,
                                                 double v, double fpa, double v_azi,
                                                 const Phase1ToPhase2Frame& frame,
                                                 const PlanetConfig& planet_config) {
    (void)planet_config;  // origin_r already carries the planet's reference radius baked in

    Phase2TranslationalState s;
    s.rz = r - frame.origin_r;
    s.rx = (lo - frame.origin_lo) * frame.origin_r * std::cos(frame.origin_la);  // East
    s.ry = (la - frame.origin_la) * frame.origin_r;                              // North

    // Same local-vertical-frame construction AeroAngles.cpp uses internally
    // (up,north,east order there); mapped to ENU (East,North,Up) here.
    s.vx = v * std::cos(fpa) * std::sin(v_azi);
    s.vy = v * std::cos(fpa) * std::cos(v_azi);
    s.vz = v * std::sin(fpa);
    return s;
}

Eigen::MatrixXd computeTransitionJacobian(double r, double la, double lo,
                                            double v, double fpa, double v_azi,
                                            const Phase1ToPhase2Frame& frame,
                                            const PlanetConfig& planet_config) {
    Eigen::VectorXd x_ref(6);
    x_ref << r, la, lo, v, fpa, v_azi;
    Eigen::VectorXd u_ref(0);  // control-free map

    EomFn f = [&](const Eigen::VectorXd& x, const Eigen::VectorXd&) -> Eigen::VectorXd {
        Phase2TranslationalState s = convertPhase1ToPhase2(x(0), x(1), x(2), x(3), x(4), x(5),
                                                              frame, planet_config);
        Eigen::VectorXd out(6);
        out << s.rx, s.ry, s.rz, s.vx, s.vy, s.vz;
        return out;
    };

    EomJacobian jac = computeEomJacobianFd(f, x_ref, u_ref);
    return jac.Ac;
}

}  // namespace guidance_scp
