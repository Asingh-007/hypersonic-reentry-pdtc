#ifndef HEATING_LOAD_MODEL_H
#define HEATING_LOAD_MODEL_H

#include "DescentDynamics.h"
#include "PlanetConfig.h"
#include "SpacecraftConfig.h"
#include "ControlInputs.h"

// Sutton-Graves stagnation-point convective heat flux (Sutton & Graves,
// 1971 -- the standard engineering correlation for blunt-body entry heating,
// matching this codebase's existing practice of well-established closed-form
// correlations over full CFD, e.g. Hoerner drag / Allen-Perkins crossflow in
// AeroRegimeDispatch.h). k_sg is the standard Earth/air value; NOT valid for
// Mars/CO2 entry (would need a different constant) -- not handled here since
// main.cpp's active scenario is Earth-only (PlanetConfig::Earth()).
//   q_dot = k_sg * sqrt(rho / R_n) * V^3   [W/m^2]
//   k_sg = 1.7415e-4   [kg^0.5 / m], rho in kg/m^3, V in m/s, R_n in m
double suttonGravesHeatFlux(double rho, double v, double nose_radius_m);

// Tauber-Sutton stagnation-point radiative heat flux for air/Earth entry
// (Tauber, M.E. and Sutton, K., "Stagnation-Point Radiative Heating
// Relations for Earth and Mars Entries," J. Spacecraft and Rockets, Vol.
// 28, No. 1, 1991, pp. 40-42 -- constants/table transcribed directly from
// the original paper, Eqs. 1-2 and Table 1, air/Earth values only, NOT the
// separate Mars/CO2 constants given in the same paper):
//   q_rad [W/cm^2] = C * R_n^a * rho^b * f(V)
//   C = 4.736e4
//   a = min(1.072e6 * V^-1.88 * rho^-0.325, 1.0)   ("a <= 1 must always be
//       met" -- direct quote; note a depends on V and rho, NOT on R_n)
//   b = 1.22
//   f(V): Table 1, linearly interpolated, V in m/s from 9000 to 16000
// Result is converted to W/m^2 (x1e4) to match suttonGravesHeatFlux's units.
// Validity per the paper: V 10-16 km/s for the C/a/b relation (f(V) table
// extends down to 9000 m/s), rho 6.66e-5 to 6.31e-4 kg/m^3, R_n 0.3-3 m.
// Below the table's 9000 m/s floor, radiative heating is genuinely
// negligible for air (the paper's own motivation: it only exceeds
// convective heating at superorbital/lunar-return-class speeds) -- this
// returns 0 there rather than extrapolating. Above 16000 m/s, f(V) is
// clamped at the last tabulated value (mild extrapolation beyond the
// paper's validated range, documented rather than silently assumed exact).
double tauberSuttonRadiativeHeatFlux(double rho, double v, double nose_radius_m);

struct HeatingLoadResult {
    double heat_flux_conv_w_m2 = 0.0;  // Sutton-Graves stagnation-point convective heat flux
    double heat_flux_rad_w_m2 = 0.0;   // Tauber-Sutton stagnation-point radiative heat flux
    double heat_flux_total_w_m2 = 0.0; // conv + rad -- the physically relevant total for TPS/constraints
    double load_factor_g = 0.0;        // total non-gravitational specific force, in g's
};

// Computes stagnation-point heating and load factor at a full 13-state
// StateVector, mirroring (not duplicating -- see DescentDynamics::
// atmosphereDensity/speedOfSound, now public static) the same atmosphere +
// aero-table lookup DescentDynamics::derivatives() does, so it can be called
// independently from integrate()'s record() lambda without changing
// derivatives()'s signature/return type.
HeatingLoadResult computeHeatingAndLoad(const DescentDynamics::StateVector& x,
                                          const PlanetConfig& planet_config,
                                          const SpacecraftConfig& spacecraft_config,
                                          const ThrustVectorControlInputs& control_inputs);

#endif // HEATING_LOAD_MODEL_H
