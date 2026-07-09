#include "HeatingLoadModel.h"
#include "AeroAngles.h"
#include <algorithm>
#include <array>
#include <cmath>

// MSVC only defines M_PI when _USE_MATH_DEFINES is set before the first
// <cmath>/<math.h> include anywhere in the translation unit -- use an
// explicit local constant instead (same fix as main.cpp/DescentDynamics.cpp).
constexpr double kPi = 3.14159265358979323846;

double suttonGravesHeatFlux(double rho, double v, double nose_radius_m) {
    constexpr double k_sg = 1.7415e-4; // kg^0.5 / m, standard Earth/air value
    return k_sg * std::sqrt(rho / nose_radius_m) * v * v * v;
}

namespace {

// Tauber & Sutton (1991) Table 1, air/Earth column (f_E), transcribed
// directly from the paper -- see HeatingLoadModel.h doc comment.
struct TauberSuttonTablePoint { double v_m_s, f; };
constexpr std::array<TauberSuttonTablePoint, 19> kTauberSuttonAirTable = {{
    {9000.0, 1.5},   {9250.0, 4.3},   {9500.0, 9.7},   {9750.0, 19.5},
    {10000.0, 35.0}, {10250.0, 55.0}, {10500.0, 81.0}, {10750.0, 115.0},
    {11000.0, 151.0}, {11500.0, 238.0}, {12000.0, 359.0}, {12500.0, 495.0},
    {13000.0, 660.0}, {13500.0, 850.0}, {14000.0, 1065.0}, {14500.0, 1313.0},
    {15000.0, 1550.0}, {15500.0, 1780.0}, {16000.0, 2040.0},
}};

// Linear interpolation of Table 1, clamped at both ends (below 9000 m/s,
// radiative heating is genuinely negligible for air -- see paper's own
// motivation; above 16000 m/s, holds the last tabulated value as a
// documented mild extrapolation beyond the paper's validated range).
double tauberSuttonAirTableF(double v) {
    if (v <= kTauberSuttonAirTable.front().v_m_s) {
        return (v < kTauberSuttonAirTable.front().v_m_s) ? 0.0 : kTauberSuttonAirTable.front().f;
    }
    if (v >= kTauberSuttonAirTable.back().v_m_s) {
        return kTauberSuttonAirTable.back().f;
    }
    for (size_t i = 0; i + 1 < kTauberSuttonAirTable.size(); ++i) {
        const auto& lo = kTauberSuttonAirTable[i];
        const auto& hi = kTauberSuttonAirTable[i + 1];
        if (v >= lo.v_m_s && v <= hi.v_m_s) {
            double t = (v - lo.v_m_s) / (hi.v_m_s - lo.v_m_s);
            return lo.f + t * (hi.f - lo.f);
        }
    }
    return kTauberSuttonAirTable.back().f; // unreachable given the bounds checks above
}

}  // namespace

double tauberSuttonRadiativeHeatFlux(double rho, double v, double nose_radius_m) {
    if (v < kTauberSuttonAirTable.front().v_m_s) {
        return 0.0; // below the table's floor -- negligible for air (see doc comment)
    }
    constexpr double C = 4.736e4;
    constexpr double b = 1.22;
    double a = std::min(1.072e6 * std::pow(v, -1.88) * std::pow(rho, -0.325), 1.0);
    double f = tauberSuttonAirTableF(v);
    double q_rad_w_cm2 = C * std::pow(nose_radius_m, a) * std::pow(rho, b) * f;
    return q_rad_w_cm2 * 1.0e4; // W/cm^2 -> W/m^2
}

HeatingLoadResult computeHeatingAndLoad(const DescentDynamics::StateVector& x,
                                          const PlanetConfig& planet_config,
                                          const SpacecraftConfig& spacecraft_config,
                                          const ThrustVectorControlInputs& control_inputs) {
    double r = x(0);
    double v = x(3);
    double fpa = x(4);
    double v_azi = x(5);
    double q1 = x(6), q2 = x(7), q3 = x(8), q4 = x(9);

    double rho = DescentDynamics::atmosphereDensity(r, planet_config);
    double a_sound = DescentDynamics::speedOfSound(r, planet_config);
    double mach = v / a_sound;

    Eigen::Quaterniond q(q4, q1, q2, q3);
    AeroAngles angles = ComputeAeroAngles(q, fpa, v_azi);
    double alpha_deg = angles.alpha_rad * 180.0 / kPi;
    double beta_deg = angles.beta_rad * 180.0 / kPi;
    // PLACEHOLDER -- matches DescentDynamics::derivatives()'s hardcoded 0.0
    // flap axes (no flap control state exists yet).
    auto aero = spacecraft_config.aero_table.interpolate(mach, alpha_deg, beta_deg, 0.0, 0.0, 0.0);

    HeatingLoadResult result;
    result.heat_flux_conv_w_m2 = suttonGravesHeatFlux(rho, v, spacecraft_config.nose_radius_m);
    result.heat_flux_rad_w_m2 = tauberSuttonRadiativeHeatFlux(rho, v, spacecraft_config.nose_radius_m);
    result.heat_flux_total_w_m2 = result.heat_flux_conv_w_m2 + result.heat_flux_rad_w_m2;

    double qbar = 0.5 * rho * v * v;
    double lift = qbar * spacecraft_config.S_ref * aero.CL;
    double drag = qbar * spacecraft_config.S_ref * aero.CD;
    double thrust = control_inputs.thrust;
    double pitch = control_inputs.pitch;
    double yaw = control_inputs.yaw;

    // Total non-gravitational specific force (thrust + aero, excluding
    // gravity -- the standard accelerometer-equivalent/sensed-acceleration
    // definition of load factor). NOTE: aero_table only stores CL/CD (no
    // CN/CA columns -- see aero/data/aero_table.csv header), so lift/drag
    // (wind axes) are used instead of body-axis CN/CA. This combines
    // cleanly with thrust's existing pitch/yaw decomposition since both are
    // already in the same velocity-relative frame in
    // DescentDynamics::derivatives(): drag opposes the v_dot-axis thrust
    // component, lift adds to the fpa_dot-axis component (see v_dot_term/
    // fpa_dot_term there) -- no side-force coefficient exists, so side is
    // thrust-only.
    double fx = thrust * std::cos(pitch) * std::cos(yaw) - drag;
    double fy = thrust * std::sin(yaw);
    double fz = thrust * std::sin(pitch) * std::cos(yaw) + lift;
    double force_mag = std::sqrt(fx * fx + fy * fy + fz * fz);

    result.load_factor_g = force_mag / (spacecraft_config.mass * planet_config.g_0);
    return result;
}
