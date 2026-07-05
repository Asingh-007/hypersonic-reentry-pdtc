#ifndef PLANET_CONFIG_H
#define PLANET_CONFIG_H

struct PlanetConfig {
    float radius = 0.0f; // radius of the planet in meters
    float mu = 0.0f; // gravitational parameter of the planet in m^3/s^2
    float rho_0 = 0.0f; // reference atmospheric density at sea level in kg/m^3
    float h_scale = 0.0f; // scale height of the atmosphere in meters
    float g_0 = 9.80665f; // standard gravity at sea level in m/s^2
    float omega = 0.0f; // angular velocity of the planet in rad/s

    PlanetConfig(float radius, float mu, float rho_0, float h_scale, float g_0, float omega)
        : radius(radius), mu(mu), rho_0(rho_0), h_scale(h_scale), g_0(g_0), omega(omega) {}
};


#endif // PLANET_CONFIG_H