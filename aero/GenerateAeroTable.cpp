// GenerateAeroTable.cpp
//
// Offline tool: sweeps the Mach-regime-dispatched aero model (see
// AeroRegimeDispatch.h) over a (Mach, alpha_deg, beta_deg, fwd_sym_deg,
// aft_sym_deg, aft_diff_deg) grid, writes a CSV lookup table, and prints a
// suggested LHS-DOE CFD anchor-point matrix for future real Fluent runs.
// No Kriging correction is applied -- no real CFD data exists yet, so this
// is the raw dispatched-model evaluation (see FUTURE WORK at the bottom).

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "PanelGeometry.h"
#include "NewtonianAero.h"
#include "AeroRegimeDispatch.h"
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
                     "call mesh.addGroup(...) manually (ids 1=fwd_left, "
                     "2=fwd_right, 3=aft_left, 4=aft_right) or flap deflection "
                     "will have no geometric effect."
                  << std::endl;
    } else {
        std::cout << "Using placeholder TestBodyGenerator body (no STL path given)." << std::endl;
        mesh = testutil::makeCylinderNoseFlapBody();
    }

    const Eigen::Vector3d moment_ref(20.0, 0.0, 0.0);  // approx CG, body frame -- PLACEHOLDER, matches TestBodyGenerator
    const double S_ref = kPi * 4.5 * 4.5;              // reference area, m^2 -- PLACEHOLDER
    const double L_ref = 9.0;                            // reference length, m -- PLACEHOLDER
    const double body_radius = 4.5;                      // PLACEHOLDER, matches TestBodyGenerator defaults
    const double body_length = 40.0;                     // PLACEHOLDER, matches TestBodyGenerator defaults

    // --- 2. Define the grid envelope ---
    // PLACEHOLDER ranges/resolution -- revisit once real mission envelopes
    // are known. mach_grid spans 0.2-12 (at least one point per regime) so
    // the subsonic/transonic placeholders are actually exercised in the
    // committed CSV, not just dead code -- this vehicle's tower-catch
    // phase is entirely subsonic. Total: 17*8*5*5^3 = 85,000 rows.
    const std::vector<double> mach_grid = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 2, 3, 4,
                                            5, 6, 7, 8, 9, 10, 11, 12};
    const std::vector<double> alpha_deg_grid = {0, 10, 20, 30, 40, 50, 60, 70};
    const std::vector<double> beta_deg_grid = {-10, -5, 0, 5, 10};
    const std::vector<double> fwd_sym_deg_grid = {-15, -7.5, 0, 7.5, 15};
    const std::vector<double> aft_sym_deg_grid = {-15, -7.5, 0, 7.5, 15};
    const std::vector<double> aft_diff_deg_grid = {-15, -7.5, 0, 7.5, 15};

    std::filesystem::path source_dir = std::filesystem::path(__FILE__).parent_path();
    std::filesystem::path data_dir = source_dir / "data";
    std::filesystem::create_directories(data_dir);
    std::filesystem::path csv_path = data_dir / "aero_table.csv";

    std::ofstream out(csv_path);
    if (!out.is_open()) {
        std::cerr << "Failed to open " << csv_path << " for writing." << std::endl;
        return 1;
    }
    out << "mach,alpha_deg,beta_deg,fwd_sym_deg,aft_sym_deg,aft_diff_deg,CL,CD,Cl_roll,Cm,Cn_yaw\n";
    out << std::setprecision(17);

    long long row_count = 0;
    for (double mach : mach_grid) {
        for (double alpha_deg : alpha_deg_grid) {
            for (double beta_deg : beta_deg_grid) {
                for (double fwd_sym_deg : fwd_sym_deg_grid) {
                    for (double aft_sym_deg : aft_sym_deg_grid) {
                        for (double aft_diff_deg : aft_diff_deg_grid) {
                            const double alpha_rad = alpha_deg * kPi / 180.0;
                            const double beta_rad = beta_deg * kPi / 180.0;
                            const double fwd_sym_rad = fwd_sym_deg * kPi / 180.0;
                            const double aft_sym_rad = aft_sym_deg * kPi / 180.0;
                            const double aft_diff_rad = aft_diff_deg * kPi / 180.0;

                            const auto flap_defl = mapFlapAxesToGroupDeflections(
                                fwd_sym_rad, aft_sym_rad, aft_diff_rad);

                            const AeroCoefficients c = evaluateAeroRegime(
                                mesh, flap_defl, alpha_rad, beta_rad, mach, moment_ref, S_ref, L_ref,
                                body_radius, body_length);

                            out << mach << "," << alpha_deg << "," << beta_deg << ","
                                << fwd_sym_deg << "," << aft_sym_deg << "," << aft_diff_deg << ","
                                << c.CL << "," << c.CD << "," << c.Cl_roll << "," << c.Cm << "," << c.Cn_yaw << "\n";
                            ++row_count;
                        }
                    }
                }
            }
        }
    }
    out.close();
    std::cout << "Wrote " << row_count << " rows to " << csv_path << std::endl;

    // --- 3. Suggested CFD anchor-point matrix (LHS-DOE) for future real
    // Fluent runs -- not used to fabricate any correction here. Mach bounds
    // span the full 0.2-12 envelope since CFD is most needed where no
    // theory applies (transonic/subsonic); kept at 3 variables (mach,
    // alpha_deg, aft_sym_deg -- the most operationally significant flap
    // axis) rather than all 6, with fwd_sym_deg/aft_diff_deg left for later. ---
    std::vector<DesignVariable> vars = {
        {"mach", 0.2, 12.0}, {"alpha_deg", 0.0, 70.0}, {"aft_sym_deg", -15.0, 15.0}};

    LatinHypercubeSampler lhs(/*seed=*/7);
    Eigen::MatrixXd unit_design = lhs.sample(/*n_samples=*/18, /*n_dims=*/3, /*n_swap_iters=*/8000);
    Eigen::MatrixXd doe = LatinHypercubeSampler::scaleToBounds(unit_design, vars);

    // Must-include points: belly-flop trim (M=6, alpha=60deg, aft_sym=0)
    // and entry-interface condition (M=10, alpha=20deg, aft_sym=0) --
    // trim-neutral aft_sym assumed for these reference anchors.
    Eigen::MatrixXd fixed_points(2, 3);
    fixed_points << 6.0, 60.0, 0.0,
                    10.0, 20.0, 0.0;
    doe = LatinHypercubeSampler::augmentWithFixedPoints(doe, fixed_points);

    std::cout << "\n=== Suggested CFD anchor-point matrix (" << doe.rows() << " points) ===\n";
    std::cout << "  idx      Mach   alpha[deg]   aft_sym[deg]\n";
    std::cout << std::fixed << std::setprecision(4);
    for (int i = 0; i < doe.rows(); ++i) {
        std::cout << "  " << std::setw(3) << i << "   " << std::setw(8) << doe(i, 0)
                   << "   " << std::setw(9) << doe(i, 1) << "   " << std::setw(12) << doe(i, 2) << "\n";
    }

    return 0;
}

// ---------------------------------------------------------------------
// FUTURE WORK (not implemented -- no real CFD data exists yet): once Fluent
// results exist at the DOE points above, fit a UniversalKriging per
// coefficient (per regime, since the trend function is regime-dependent)
// on the CFD residuals, re-evaluate over the same 6D grid using
// krig.predict() instead of evaluateAeroRegime(), and overwrite this CSV --
// AeroCoefficientTable/DescentDynamics need zero code changes for that.
// ---------------------------------------------------------------------
