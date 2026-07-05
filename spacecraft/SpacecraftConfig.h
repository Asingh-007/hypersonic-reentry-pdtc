#ifndef SPACECRAFT_CONFIG_H
#define SPACECRAFT_CONFIG_H

struct SpacecraftConfig {
    float mass = 0.0f; // mass of the spacecraft in kg
    float area = 0.0f; // reference area of the spacecraft in m^2
    float c_d = 0.0f; // drag coefficient of the spacecraft
    float c_l = 0.0f; // lift coefficient of the spacecraft

    // PLACEHOLDER: no real vehicle geometry/mass-properties data exists yet.
    // Diagonal (principal-axis) inertia tensor components in kg*m^2; off-diagonal
    // products of inertia are assumed zero (principal-axis alignment).
    // TODO: replace with actual vehicle inertia data once available.
    float ixx = 1000.0f;
    float iyy = 1000.0f;
    float izz = 1500.0f;

    SpacecraftConfig(float mass, float area, float c_d, float c_l,
                      float ixx, float iyy, float izz)
        : mass(mass), area(area), c_d(c_d), c_l(c_l),
          ixx(ixx), iyy(iyy), izz(izz) {}
};

#endif // SPACECRAFT_CONFIG_H