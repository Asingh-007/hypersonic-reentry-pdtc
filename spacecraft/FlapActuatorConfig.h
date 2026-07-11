#ifndef SPACECRAFT_FLAP_ACTUATOR_CONFIG_H
#define SPACECRAFT_FLAP_ACTUATOR_CONFIG_H

// PLACEHOLDER: no real flap-actuator hardware spec exists yet. Values are
// order-of-magnitude plausible for a ~1000 kg entry vehicle's electro-
// mechanical hinge-line actuators (not hydraulic), shared across all 4
// flaps -- no per-flap asymmetry modeled, a reasonable simplification given
// FlapHingeData.h's 4 flaps are geometrically similar pairs.
// TODO: replace with actual actuator hardware data once available.
struct FlapActuatorConfig {
    double Jeff_kg_m2 = 50.0;       // effective actuator+flap inertia about the hinge axis, kg*m^2
    double Jeff_inv = 1.0 / 50.0;   // kept manually in sync with Jeff_kg_m2
    double N_gear_ratio = 100.0;    // motor-to-hinge gear reduction, dimensionless
    double b_damping_n_m_s = 20.0;  // viscous damping coefficient, N*m/(rad/s)
    double tau_m_max_n_m = 500.0;   // max motor torque (pre-gearbox), N*m
    double delta_max_rad = 0.34906585039886590;      // +/-20 deg deflection limit
    double delta_dot_max_rad_s = 0.5;                // max flap deflection RATE (box bound on the ddot_i state)
};

#endif  // SPACECRAFT_FLAP_ACTUATOR_CONFIG_H
