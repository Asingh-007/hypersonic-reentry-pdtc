// GenerateAeroTable.cpp
//
// Offline tool: sweeps the Newtonian panel model over a (Mach, alpha_deg,
// beta_deg, flap_deg) grid and writes a CSV lookup table that
// AeroCoefficientTable / DescentDynamics consume at runtime. Also prints a
// suggested LHS-DOE CFD anchor-point matrix for future real Fluent runs.
//
// This tool intentionally does NOT fabricate a "Kriging-corrected" table --
// no real CFD data exists yet, so the table written here is the raw
// Newtonian model evaluated on the grid. See the FUTURE WORK comment block
// at the bottom of this file for the intended correction hook once real CFD
// results exist.

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "PanelGeometry.h"
#include "NewtonianAero.h"
#include "TestBodyGenerator.h"
#include "StlMeshLoader.h"
#include "LatinHypercubeSampler.h"

using namespace aero_model;

// MSVC only defines M_PI when _USE_MATH_DEFINES is set before the first
// <cmath>/<math.h> include anywhere in the translation unit -- use an
// explicit local constant instead (same fix as main.cpp).
constexpr double kPi = 3.14159265358979323846;

int main(int argc, char** argv) {
    // --- 1. Build the vehicle mesh ---
    PanelMesh mesh;
    if (argc > 1) {
        std::cout << "Loading mesh from STL: " << argv[1] << std::endl;
        mesh = LoadMeshFromStl(argv[1], /*group_id=*/0);
        std::cerr << "Note: STL import has no concept of flap hinge groups -- "
                     "if this vehicle has deflectable surfaces, call "
                     "mesh.addGroup(...) manually before sweeping flap_deg, "
                     "otherwise flap deflection has no geometric effect and "
                     "the flap_deg column will not vary the result."
                  << std::endl;
    } else {
        std::cout << "Using placeholder TestBodyGenerator body (no STL path given)." << std::endl;
        mesh = testutil::makeCylinderNoseFlapBody();
    }

    NewtonianAeroModel newton;
    const Eigen::Vector3d moment_ref(20.0, 0.0, 0.0);  // approx CG, body frame -- PLACEHOLDER, matches TestBodyGenerator
    const double S_ref = kPi * 4.5 * 4.5;              // reference area, m^2 -- PLACEHOLDER
    const double L_ref = 9.0;                            // reference length, m -- PLACEHOLDER

    // --- 2. Define the grid envelope ---
    // PLACEHOLDER ranges/resolution -- revisit once real mission Mach/alpha/
    // beta/flap envelopes are known. Matches the original demo's envelope
    // for mach/alpha/flap; beta added since DescentDynamics now derives a
    // real 3D angle of attack + sideslip from the attitude quaternion.
    const std::vector<double> mach_grid = {2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
    const std::vector<double> alpha_deg_grid = {0, 10, 20, 30, 40, 50, 60, 70};
    const std::vector<double> beta_deg_grid = {-10, -5, 0, 5, 10};
    const std::vector<double> flap_deg_grid = {-15, -10, -5, 0, 5, 10, 15};

    std::filesystem::path source_dir = std::filesystem::path(__FILE__).parent_path();
    std::filesystem::path data_dir = source_dir / "data";
    std::filesystem::create_directories(data_dir);
    std::filesystem::path csv_path = data_dir / "aero_table.csv";

    std::ofstream out(csv_path);
    if (!out.is_open()) {
        std::cerr << "Failed to open " << csv_path << " for writing." << std::endl;
        return 1;
    }
    out << "mach,alpha_deg,beta_deg,flap_deg,CL,CD,Cl_roll,Cm,Cn_yaw\n";
    out << std::setprecision(17);

    long long row_count = 0;
    for (double mach : mach_grid) {
        for (double alpha_deg : alpha_deg_grid) {
            for (double beta_deg : beta_deg_grid) {
                for (double flap_deg : flap_deg_grid) {
                    double alpha_rad = alpha_deg * kPi / 180.0;
                    double beta_rad = beta_deg * kPi / 180.0;
                    double flap_rad = flap_deg * kPi / 180.0;
                    // Flap group id 1 = TestBodyGenerator's single aft flap
                    // (see the note above for STL-imported meshes).
                    std::unordered_map<int, double> flap_defl{{1, flap_rad}};

                    AeroCoefficients c = newton.evaluate(mesh, flap_defl, alpha_rad, beta_rad,
                                                          mach, moment_ref, S_ref, L_ref);
                    out << mach << "," << alpha_deg << "," << beta_deg << "," << flap_deg << ","
                        << c.CL << "," << c.CD << "," << c.Cl_roll << "," << c.Cm << "," << c.Cn_yaw << "\n";
                    ++row_count;
                }
            }
        }
    }
    out.close();
    std::cout << "Wrote " << row_count << " rows to " << csv_path << std::endl;

    // --- 3. Suggested CFD anchor-point matrix (LHS-DOE), for future real
    // Fluent runs -- NOT used to fabricate any correction here. ---
    std::vector<DesignVariable> vars = {
        {"mach", 2.0, 12.0}, {"alpha_deg", 0.0, 70.0}, {"flap_deg", -15.0, 15.0}};

    LatinHypercubeSampler lhs(/*seed=*/7);
    Eigen::MatrixXd unit_design = lhs.sample(/*n_samples=*/18, /*n_dims=*/3, /*n_swap_iters=*/8000);
    Eigen::MatrixXd doe = LatinHypercubeSampler::scaleToBounds(unit_design, vars);

    // Must-include points: belly-flop trim (M=6, alpha=60deg, flap=0) and
    // entry-interface condition (M=10, alpha=20deg, flap=0).
    Eigen::MatrixXd fixed_points(2, 3);
    fixed_points << 6.0, 60.0, 0.0,
                    10.0, 20.0, 0.0;
    doe = LatinHypercubeSampler::augmentWithFixedPoints(doe, fixed_points);

    std::cout << "\n=== Suggested CFD anchor-point matrix (" << doe.rows() << " points) ===\n";
    std::cout << "  idx      Mach   alpha[deg]   flap[deg]\n";
    std::cout << std::fixed << std::setprecision(4);
    for (int i = 0; i < doe.rows(); ++i) {
        std::cout << "  " << std::setw(3) << i << "   " << std::setw(8) << doe(i, 0)
                   << "   " << std::setw(9) << doe(i, 1) << "   " << std::setw(9) << doe(i, 2) << "\n";
    }

    return 0;
}

// ---------------------------------------------------------------------
// FUTURE WORK (not implemented -- no real CFD data exists yet):
// Once Fluent results exist at the DOE points printed above:
//   1. Load CFD-observed CL/CD/Cl_roll/Cm/Cn_yaw at each DOE point.
//   2. For each coefficient, construct a UniversalKriging instance whose
//      trend function wraps NewtonianAeroModel::evaluate(...).<coef>,
//      fit() on the DOE points + CFD values.
//   3. Re-evaluate over the SAME (mach_grid, alpha_deg_grid, beta_deg_grid,
//      flap_deg_grid) used above, but querying krig.predict(x).mean instead
//      of the raw Newtonian evaluate() call, and overwrite aero_table.csv.
//   4. AeroCoefficientTable and DescentDynamics require ZERO code changes
//      for this -- they only ever see the CSV.
// ---------------------------------------------------------------------
