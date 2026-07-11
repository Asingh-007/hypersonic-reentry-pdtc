#ifndef AERO_COEFFICIENT_TABLE_H
#define AERO_COEFFICIENT_TABLE_H

#include <array>
#include <string>
#include <vector>

#include "NewtonianAero.h"

namespace aero_model {

// Lightweight regular-grid multilinear interpolator over a precomputed
// (mach, alpha_deg, beta_deg, fwd_sym_deg, aft_sym_deg, aft_diff_deg) ->
// AeroCoefficients table generated offline by GenerateAeroTable.cpp. 

class AeroCoefficientTable {
public:
    // Loads a CSV with header row (fixed column order):
    //   mach,alpha_deg,beta_deg,fwd_sym_deg,aft_sym_deg,aft_diff_deg,CL,CD,Cl_roll,Cm,Cn_yaw,Ch1,Ch2,Ch3,Ch4

    bool load(const std::string& csv_path);

    // Multilinear interpolation (2^6 = 64-corner blend) within the grid
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

    // Get flat index given grid coordinates
    std::size_t flatIndex(const std::array<int, kNumAxes>& idx) const;
};

}  // namespace aero_model

#endif // AERO_COEFFICIENT_TABLE_H
