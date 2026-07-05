#include "PlanetConfig.h"

// Earth Planet Config
PlanetConfig PlanetConfig::Earth() {
    return PlanetConfig(
        6371000.0f,      // radius, m (mean radius)
        3.986004418e14f, // mu, m^3/s^2 (standard gravitational parameter, WGS-84)
        9.80665f,        // g_0, m/s^2 (standard gravity)
        7.2921159e-5f,   // omega, rad/s (mean sidereal rotation rate)
        1.08263e-3f,     // J2 (WGS-84 / EGM96)
        PlanetBody::Earth
    );
}

// Mars Planet Config
PlanetConfig PlanetConfig::Mars() {
    return PlanetConfig(
        3389500.0f,    // radius, m (mean radius)
        4.282837e13f,  // mu, m^3/s^2
        3.72076f,      // g_0, m/s^2 (surface gravity)
        7.088218e-5f,  // omega, rad/s (rotation rate)
        1.96045e-3f,   // J2
        PlanetBody::Mars
    );
}
