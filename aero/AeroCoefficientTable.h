#ifndef AERO_COEFFICIENT_TABLE_H
#define AERO_COEFFICIENT_TABLE_H

#include <array>
#include <string>
#include <vector>

#include "NewtonianAero.h"

namespace aero_model {

// Lightweight regular-grid multilinear interpolator over a precomputed
// (mach, alpha_deg, beta_deg, fwd_sym_deg, aft_sym_deg, aft_diff_deg) ->
// AeroCoefficients table generated offline by GenerateAeroTable.cpp. This
// is the ONLY aero query DescentDynamics makes at runtime -- O(1) lookup
// regardless of model complexity, and swapping in a future Kriging-
// corrected table needs zero runtime code changes.
//
// The 6 axes are fixed/named for call-site readability; only the
// interpolation math (a generic 2^6-corner bit-indexed blend) is generic,
// so a future 7th axis only touches kNumAxes/load()/interpolate()'s signature.
class AeroCoefficientTable {
public:
    // Loads a CSV with header row (fixed column order):
    //   mach,alpha_deg,beta_deg,fwd_sym_deg,aft_sym_deg,aft_diff_deg,CL,CD,Cl_roll,Cm,Cn_yaw
    // Rows must form a complete regular/rectilinear 6D grid (every
    // combination of the distinct values present on each axis must have
    // exactly one row) -- throws std::runtime_error if the grid is
    // incomplete/malformed. Returns false only for I/O failure to open
    // the file.
    bool load(const std::string& csv_path);

    // Multilinear interpolation (2^6 = 64-corner blend) within the grid;
    // query points outside the grid's min/max on any axis are clamped to
    // that axis' nearest boundary rather than extrapolated (extrapolating
    // a surrogate table well outside its fitted envelope is more
    // dangerous than saturating at the edge).
    AeroCoefficients interpolate(double mach, double alpha_deg, double beta_deg,
                                  double fwd_sym_deg, double aft_sym_deg, double aft_diff_deg) const;

    bool isLoaded() const { return loaded_; }

private:
    static constexpr int kNumAxes = 6;

    // grids_[axis], axis in the fixed order documented above (0=mach,
    // 1=alpha_deg, 2=beta_deg, 3=fwd_sym_deg, 4=aft_sym_deg, 5=aft_diff_deg).
    std::vector<std::vector<double>> grids_;
    // Flattened, row-major over grids_ in the order above.
    std::vector<AeroCoefficients> values_;
    bool loaded_ = false;

    // Finds the bracketing grid indices [lo, hi] and interpolation weight
    // t in [0,1] such that value ~= lerp(grid[lo], grid[hi], t), clamping
    // at the grid boundaries.
    static void findBracket(const std::vector<double>& grid, double value, int& lo, int& hi, double& t);

    std::size_t flatIndex(const std::array<int, kNumAxes>& idx) const;
};

}  // namespace aero_model

#endif // AERO_COEFFICIENT_TABLE_H
