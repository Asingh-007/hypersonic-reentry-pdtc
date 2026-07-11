// GenerateAeroTable.cpp
//
// Offline tool: sweeps the Mach-regime-dispatched aero model (see
// AeroRegimeDispatch.h) over a (Mach, alpha_deg, beta_deg, fwd_sym_deg,
// aft_sym_deg, aft_diff_deg) grid, writes a CSV lookup table, and writes/
// prints a suggested LHS-DOE CFD anchor-point matrix for future real
// Fluent runs. No Kriging correction is applied as no real CFD data exists
// yet, so this is the raw dispatched-model evaluation (see FUTURE WORK at
// the bottom).
//
// Pass --doe-only to (re)generate aero/data/doe_points.csv and
// reference_quantities.csv (mesh load only, a few seconds) without
// re-running the full grid sweep (tens of minutes for the real spacecraft
// geometry), useful for downstream tooling like GenerateDeflectedGeometry.cpp
// and cfd/run_doe_point.py that only need those two files.

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
#include "SpacecraftGeometry.h"
#include "LatinHypercubeSampler.h"

using namespace aero_model;

constexpr double kPi = 3.14159265358979323846;

int main(int argc, char** argv) {
    bool doe_only = false;
    std::string stl_override;
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "--doe-only") doe_only = true;
        else stl_override = argv[i];
    }

    const std::filesystem::path source_dir = std::filesystem::path(__FILE__).parent_path();
    const std::filesystem::path data_dir = source_dir / "data";
    std::filesystem::create_directories(data_dir);

    // Suggested CFD anchor-point matrix (LHS-DOE) for future real
    // Fluent runs, not used to fabricate any correction here. Mach bounds
    // span the full 0.2-12 envelope since CFD is most needed where no
    // theory applies (transonic/subsonic); kept at 3 variables (mach,
    // alpha_deg, aft_sym_deg as the most operationally significant flap
    // axis) rather than all 6, with fwd_sym_deg/aft_diff_deg left for later.
    // Runs first (independent of the mesh/grid sweep below) so --doe-only
    // can regenerate it cheaply. ---
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

    const std::filesystem::path doe_csv_path = data_dir / "doe_points.csv";
    std::ofstream doe_out(doe_csv_path);
    if (!doe_out.is_open()) {
        std::cerr << "Failed to open " << doe_csv_path << " for writing." << std::endl;
        return 1;
    }
    doe_out << "idx,mach,alpha_deg,aft_sym_deg\n" << std::setprecision(17);
    for (int i = 0; i < doe.rows(); ++i) {
        doe_out << i << "," << doe(i, 0) << "," << doe(i, 1) << "," << doe(i, 2) << "\n";
    }
    std::cout << "Wrote DOE points to " << doe_csv_path << std::endl;

    // Vehicle Mesh Construction
    // Default: the real SolidWorks-exported geometry under geometry/ (see
    // SpacecraftGeometry.h). A non-flag argv, if given, overrides with a
    // single STL file instead (body-only, no flap groups).
    PanelMesh mesh;
    Eigen::Vector3d moment_ref;
    double S_ref, L_ref, body_radius, body_length;
    if (!stl_override.empty()) {
        std::cout << "Loading mesh from STL: " << stl_override << std::endl;
        mesh = LoadMeshFromStl(stl_override, /*group_id=*/0);
        std::cerr << "Note: STL import has no concept of flap hinge groups -- "
                     "call mesh.addGroup(...) manually (ids 1=fwd_left, "
                     "2=fwd_right, 3=aft_left, 4=aft_right) or flap deflection "
                     "will have no geometric effect."
                  << std::endl;
        moment_ref = Eigen::Vector3d(20.0, 0.0, 0.0);  // PLACEHOLDER, matches TestBodyGenerator
        S_ref = kPi * 4.5 * 4.5;
        L_ref = 9.0;
        body_radius = 4.5;
        body_length = 40.0;
    } else {
        const std::filesystem::path geometry_dir = source_dir.parent_path() / "geometry";
        std::cout << "Loading real spacecraft geometry from " << geometry_dir << std::endl;
        SpacecraftGeometry geo = LoadSpacecraftGeometry(geometry_dir.string());
        mesh = geo.mesh;
        moment_ref = geo.moment_ref;
        S_ref = geo.S_ref;
        L_ref = geo.L_ref;
        body_radius = geo.body_radius;
        body_length = geo.body_length;
        std::cout << "  S_ref=" << S_ref << " L_ref=" << L_ref << " body_radius=" << body_radius
                   << " body_length=" << body_length << " moment_ref=(" << moment_ref.transpose() << ")"
                   << std::endl;
    }

    // Written so downstream tooling (e.g. cfd/run_doe_point.py) reads these
    // once-computed values instead of hardcoding a second copy that could
    // drift out of sync. moment_ref is in this codebase's model frame: W(meters, +X=nose)
    {
        std::ofstream ref_out((data_dir / "reference_quantities.csv").string());
        ref_out << "S_ref,L_ref,body_radius,body_length,moment_ref_x,moment_ref_y,moment_ref_z\n"
                << std::setprecision(17)
                << S_ref << "," << L_ref << "," << body_radius << "," << body_length << ","
                << moment_ref.x() << "," << moment_ref.y() << "," << moment_ref.z() << "\n";
    }

    if (doe_only) {
        return 0;
    }

    // Grid Envelope Definition
    // PLACEHOLDER ranges/resolution, needs to berevisited once real mission envelopes
    // are known. mach_grid spans 0.2-12 (at least one point per regime) so
    // the subsonic/transonic placeholders are actually exercised in the
    // committed CSV, not just dead code. This vehicle's tower-catch
    // phase is entirely subsonic. Total: 17*8*5*5^3 = 85,000 rows.
    const std::vector<double> mach_grid = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 2, 3, 4,
                                            5, 6, 7, 8, 9, 10, 11, 12};
    const std::vector<double> alpha_deg_grid = {0, 10, 20, 30, 40, 50, 60, 70};
    const std::vector<double> beta_deg_grid = {-10, -5, 0, 5, 10};
    const std::vector<double> fwd_sym_deg_grid = {-15, -7.5, 0, 7.5, 15};
    const std::vector<double> aft_sym_deg_grid = {-15, -7.5, 0, 7.5, 15};
    const std::vector<double> aft_diff_deg_grid = {-15, -7.5, 0, 7.5, 15};

    const std::filesystem::path csv_path = data_dir / "aero_table.csv";
    std::ofstream out(csv_path);
    if (!out.is_open()) {
        std::cerr << "Failed to open " << csv_path << " for writing." << std::endl;
        return 1;
    }
    out << "mach,alpha_deg,beta_deg,fwd_sym_deg,aft_sym_deg,aft_diff_deg,CL,CD,Cl_roll,Cm,Cn_yaw,"
           "Ch1,Ch2,Ch3,Ch4\n";
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
                                << c.CL << "," << c.CD << "," << c.Cl_roll << "," << c.Cm << "," << c.Cn_yaw << ","
                                << c.Ch[0] << "," << c.Ch[1] << "," << c.Ch[2] << "," << c.Ch[3] << "\n";
                            ++row_count;
                        }
                    }
                }
            }
        }
    }
    out.close();
    std::cout << "Wrote " << row_count << " rows to " << csv_path << std::endl;

    return 0;
}

// FUTURE WORK (not implemented -- no real CFD data exists yet): once Fluent
// results exist at the DOE points above, fit a UniversalKriging per
// coefficient (per regime, since the trend function is regime-dependent)
// on the CFD residuals, re-evaluate over the same 6D grid using
// krig.predict() instead of evaluateAeroRegime(), and overwrite this CSV --
// AeroCoefficientTable/DescentDynamics need zero code changes for that.

