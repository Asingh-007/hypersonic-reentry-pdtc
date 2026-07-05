#include "AtmosphereModel.h"
#include <cmath>

const std::array<EarthAtmosphere1976::BridgeLayer, 3>& EarthAtmosphere1976::BridgeLayerTable() {
    // Table 2, Olsen & Bettinger (2022).
    static const std::array<BridgeLayer, 3> table = {{
        { 85000.0, 7.726e-6, 0.1545455, 197.9740},  // 84-90 km
        { 99000.0, 4.504e-7, 0.1189286, 128.4577},  // 91-106 km
        {110000.0, 5.930e-8, 0.5925240, 432.8484},  // 107-120 km
    }};
    return table;
}

const EarthAtmosphere1976::BridgeLayer& EarthAtmosphere1976::FindBridgeLayer(double h_m) {
    const auto& table = BridgeLayerTable();
    if (h_m <= 90000.0) return table[0];
    if (h_m <= 106000.0) return table[1];
    return table[2];
}

AtmosphereState EarthAtmosphere1976::LowAltitudeModel(double z_geometric_m) {
    AtmosphereState state;
    state.density = kSeaLevelDensity * std::exp(-kBeta * z_geometric_m);
    return state;
}

AtmosphereState EarthAtmosphere1976::BridgeModel(double z_geometric_m) {
    const BridgeLayer& layer = FindBridgeLayer(z_geometric_m);
    const double term = 1.0 + layer.delta_h * (z_geometric_m - layer.h_i_m) / kEarthRadiusM;
    const double exponent = (1.0 + layer.alpha) / layer.alpha;
    const double rho = layer.rho_i * std::pow(1.0 / term, exponent);

    AtmosphereState state;
    state.density = rho;
    return state;
}

AtmosphereState EarthAtmosphere1976::PowerLawModel(double z_geometric_m) {
    const double h_km = z_geometric_m / 1000.0;
    AtmosphereState state;
    state.density = 4.50847623e7 * std::pow(h_km, -7.44605852);
    return state;
}

AtmosphereState EarthAtmosphere1976::Compute(double z_geometric_m) {
    if (z_geometric_m < kLowAltitudeTopM) {
        return LowAltitudeModel(z_geometric_m);
    } else if (z_geometric_m <= kBridgeTopM) {
        return BridgeModel(z_geometric_m);
    } else if (z_geometric_m <= kCurveFitTopM) {
        return PowerLawModel(z_geometric_m);
    }
    return AtmosphereState{};
}

AtmosphereState MarsAtmosphereExponential::Compute(double z_geometric_m) {
    AtmosphereState state;
    state.density = kSeaLevelDensity * std::exp(-z_geometric_m / kScaleHeightM);
    state.temperature = kRefTemperatureK;
    state.pressure = state.density * kRCO2 * state.temperature;
    return state;
}
