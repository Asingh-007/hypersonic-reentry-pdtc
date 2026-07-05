#ifndef AERO_COEFFICIENT_TABLE_H
#define AERO_COEFFICIENT_TABLE_H

#include "NewtonianAero.h"
#include <string>
#include <vector>

namespace aero_model {

// Lightweight regular-grid quadrilinear interpolator over a precomputed
// (mach, alpha_deg, beta_deg, flap_deg) -> AeroCoefficients table, generated
// offline by aero/GenerateAeroTable.cpp. This is the ONLY aero query
// DescentDynamics makes at runtime -- it never calls
// NewtonianAeroModel::evaluate() or UniversalKriging::predict() directly
// (RK45 calls derivatives() ~7x per accepted step, and table lookup is O(1)
// regardless of underlying mesh/model complexity; it also means swapping a
// raw-Newtonian table for a future Kriging-corrected one needs zero runtime
// code changes).
class AeroCoefficientTable {
public:
    // Loads a CSV with header row:
    //   mach,alpha_deg,beta_deg,flap_deg,CL,CD,Cl_roll,Cm,Cn_yaw
    // Rows must form a complete regular/rectilinear 4D grid (every
    // combination of the distinct mach/alpha_deg/beta_deg/flap_deg values
    // present must have exactly one row) -- throws std::runtime_error if
    // the grid is incomplete/malformed. Returns false only for I/O failure
    // to open the file.
    bool load(const std::string& csv_path);

    // Quadrilinear interpolation within the grid; query points outside the
    // grid's min/max on any axis are clamped to that axis' nearest boundary
    // rather than extrapolated (extrapolating a Newtonian/CFD surrogate
    // table well outside its fitted envelope is more dangerous than
    // saturating at the edge).
    AeroCoefficients interpolate(double mach, double alpha_deg, double beta_deg, double flap_deg) const;

    bool isLoaded() const { return loaded_; }

private:
    std::vector<double> mach_grid_, alpha_grid_, beta_grid_, flap_grid_;
    // Flattened, row-major over (mach_idx, alpha_idx, beta_idx, flap_idx).
    std::vector<AeroCoefficients> values_;
    bool loaded_ = false;

    // Finds the bracketing grid indices [lo, hi] and interpolation weight
    // t in [0,1] such that value ~= lerp(grid[lo], grid[hi], t), clamping
    // at the grid boundaries.
    static void findBracket(const std::vector<double>& grid, double value, int& lo, int& hi, double& t);

    const AeroCoefficients& at(int im, int ia, int ib, int ifl) const;
};

}  // namespace aero_model

#endif // AERO_COEFFICIENT_TABLE_H
