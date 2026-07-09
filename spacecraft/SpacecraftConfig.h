#ifndef SPACECRAFT_CONFIG_H
#define SPACECRAFT_CONFIG_H

#include "AeroCoefficientTable.h"
#include <Eigen/Dense>
#include <stdexcept>
#include <string>

struct SpacecraftConfig {
    float mass = 0.0f; // mass of the spacecraft in kg

    // PLACEHOLDER: no real vehicle mass-properties data exists yet.
    // Full symmetric inertia tensor (kg*m^2), body axes. Off-diagonal
    // products of inertia are supported (no longer assumed zero) -- pass a
    // diagonal matrix if principal-axis alignment is a valid assumption
    // for your vehicle.
    // TODO: replace with actual vehicle inertia data once available.
    Eigen::Matrix3d inertia = Eigen::Vector3d(1000.0, 1000.0, 1500.0).asDiagonal();

    // PLACEHOLDER: reference quantities for aerodynamic non-dimensionalization,
    // matching aero_model::testutil::makeCylinderNoseFlapBody()'s default
    // geometry (radius=4.5 m, length=40 m). Replace once a real OML/CG
    // exists (see aero/StlMeshLoader.h for importing a real OpenVSP/CAD mesh).
    double S_ref = 0.0;                       // reference area, m^2 (pi * r^2 for the placeholder body)
    double L_ref = 0.0;                       // reference length, m
    Eigen::Vector3d moment_ref = Eigen::Vector3d::Zero(); // moment reference point, body frame (approx CG)

    // PLACEHOLDER: nose radius for stagnation-point heating (Sutton-Graves), m.
    // This vehicle is a pointed/rounded nose-cone shape (Starship-like), NOT a
    // blunt Apollo-style capsule, so nose_radius_m is deliberately much smaller
    // than body_radius (aero/data/reference_quantities.csv: body_radius =
    // 4.63296 m) -- picked as ~0.18x body_radius, NOT auto-fit from STL geometry
    // (no nose-tip-curvature extraction exists yet -- see aero/StlMeshLoader.h).
    // TODO: replace with actual vehicle nose-cap radius once available.
    double nose_radius_m = 0.0;

    // Precomputed (mach, alpha_deg, beta_deg, fwd_sym_deg, aft_sym_deg,
    // aft_diff_deg) -> AeroCoefficients lookup table, generated offline by
    // aero/GenerateAeroTable.cpp. DescentDynamics only ever queries this table.
    aero_model::AeroCoefficientTable aero_table;

    // NOTE: unlike every other config struct in this repo, this constructor
    // does real work (loads a file) and can throw -- a deliberate deviation
    // from the usual "plain aggregate, no-op constructor" convention, since
    // an invalid vehicle config (aero table failed to load) shouldn't
    // silently construct into a broken-but-alive state.
    SpacecraftConfig(float mass, const Eigen::Matrix3d& inertia,
                      double S_ref, double L_ref, const Eigen::Vector3d& moment_ref,
                      double nose_radius_m,
                      const std::string& aero_table_csv_path)
        : mass(mass), inertia(inertia), S_ref(S_ref), L_ref(L_ref), moment_ref(moment_ref),
          nose_radius_m(nose_radius_m) {
        if (!aero_table.load(aero_table_csv_path)) {
            throw std::runtime_error("SpacecraftConfig: failed to load aero table from " + aero_table_csv_path);
        }
    }
};

#endif // SPACECRAFT_CONFIG_H
