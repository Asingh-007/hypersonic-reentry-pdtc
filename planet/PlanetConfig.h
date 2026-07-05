#ifndef PLANET_CONFIG_H
#define PLANET_CONFIG_H

// Selects which atmosphere model DescentDynamics uses for a given PlanetConfig.
enum class PlanetBody { Earth, Mars };

struct PlanetConfig {
    float radius = 0.0f; // reference radius of the planet in meters (also used as r_ref for J2 gravity)
    float mu = 0.0f; // gravitational parameter of the planet in m^3/s^2
    float g_0 = 9.80665f; // standard gravity at sea level in m/s^2
    float omega = 0.0f; // angular velocity of the planet in rad/s
    float j2 = 0.0f; // J2 (oblateness) gravitational perturbation coefficient, dimensionless
    PlanetBody body = PlanetBody::Earth; // selects the atmosphere model DescentDynamics uses

    PlanetConfig(float radius, float mu, float g_0, float omega, float j2, PlanetBody body)
        : radius(radius), mu(mu), g_0(g_0), omega(omega), j2(j2), body(body) {}

    // Earth reference values (WGS-84 / EGM96 standard constants).
    static PlanetConfig Earth();

    // Mars reference values (standard/textbook constants).
    static PlanetConfig Mars();
};

#endif // PLANET_CONFIG_H
