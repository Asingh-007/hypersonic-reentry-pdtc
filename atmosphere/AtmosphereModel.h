#ifndef ATMOSPHERE_MODEL_H
#define ATMOSPHERE_MODEL_H

#include <array>

// Atmosphere models used by DescentDynamics.
//
//   - Earth: Olsen & Bettinger (2022), "Six Degree-of-Freedom Analysis of the
//            Apollo 10 Atmospheric Reentry," AIAA 2022-2273 -- a 3-piece
//            atmospheric DENSITY model: simple exponential below 84 km, a
//            3-sub-layer polytropic bridge for 84-120 km (Table 2 of that
//            paper), and a power-law curve fit for 120 km-1000 km. This is a
//            density-only model (the paper doesn't give temperature/pressure
//            formulas for it), so AtmosphereState's pressure/temperature
//            fields are left at 0 here -- only .density is meaningful.
//   - Mars:  placeholder exponential model, to be replaced with a Mars-GRAM
//            or MSL/M2020-heritage curve fit once a Mars EDL mission profile
//            is finalized.
//
// All quantities are SI (meters, kelvin, pascals, kg/m^3). Input altitude is
// GEOMETRIC altitude above the reference spherical planet radius, consistent
// with DescentDynamics' rotating spherical-planet EOM.

struct AtmosphereState {
    double density = 0.0;     // kg/m^3
    double pressure = 0.0;    // Pa -> not modeled by EarthAtmosphere1976 (see above); always 0 there
    double temperature = 0.0; // K  -> not modeled by EarthAtmosphere1976 (see above); always 0 there
};

class EarthAtmosphere1976 {
public:
    static AtmosphereState Compute(double geometric_altitude_m);

    static constexpr double kSeaLevelDensity = 1.225; // kg/m^3
    // Exponential decay rate for h < 84 km. Not given directly by Table 2;
    // chosen (matching a prior placeholder value) to be nearly continuous
    // with the bridge's first sub-layer at h=84km (residual discontinuity
    // ~2%, which is an inherent feature of stitching an independent
    // low-altitude fit to Table 2's piecewise-polytropic bridge, not a bug).
    // Public (not private) since guidance_scp's path-constraint gradients
    // need d(rho)/dr = -kBeta*rho directly (see reference_stage_a.cpp).
    static constexpr double kBeta = 0.14e-3;
    static constexpr double kLowAltitudeTopM = 84000.0; // exponential-regime upper bound

private:
    // Reference Earth radius used in the bridge formula's (h-h_i)/r_earth
    // term as the same mean radius used elsewhere in this codebase
    // (PlanetConfig::Earth().radius), not the WGS-84 geopotential constant.
    static constexpr double kEarthRadiusM = 6371000.0;

    static constexpr double kBridgeTopM = 120000.0;
    static constexpr double kCurveFitTopM = 1000000.0;

    // Table 2 (Olsen & Bettinger 2022): one row per bridge sub-layer.
    struct BridgeLayer {
        double h_i_m;   // reference altitude for this sub-layer, m
        double rho_i;   // density at h_i, kg/m^3
        double alpha;
        double delta_h;
    };
    static const std::array<BridgeLayer, 3>& BridgeLayerTable();
    static const BridgeLayer& FindBridgeLayer(double h_m);

    static AtmosphereState LowAltitudeModel(double z_geometric_m);
    static AtmosphereState BridgeModel(double z_geometric_m);
    static AtmosphereState PowerLawModel(double z_geometric_m);
};

class MarsAtmosphereExponential {
public:
    static constexpr double kSeaLevelDensity = 0.020;
    static constexpr double kScaleHeightM = 11100.0;
    static constexpr double kRefTemperatureK = 210.0;
    static constexpr double kRCO2 = 191.84;
    static constexpr double kValidTopM = 40000.0;

    static AtmosphereState Compute(double geometric_altitude_m);
};

#endif // ATMOSPHERE_MODEL_H
