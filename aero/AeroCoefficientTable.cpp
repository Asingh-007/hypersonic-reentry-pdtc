#include "AeroCoefficientTable.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace aero_model {

namespace {

constexpr int kNumCoeffCols = 5;  // CL, CD, Cl_roll, Cm, Cn_yaw

std::vector<double> sortedUnique(std::vector<double> v) {
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end(), [](double a, double b) {
        return std::abs(a - b) < 1e-9;
    }), v.end());
    return v;
}

int indexOf(const std::vector<double>& grid, double value) {
    for (std::size_t i = 0; i < grid.size(); ++i) {
        if (std::abs(grid[i] - value) < 1e-9) return static_cast<int>(i);
    }
    return -1;
}

}  // namespace

bool AeroCoefficientTable::load(const std::string& csv_path) {
    std::ifstream file(csv_path);
    if (!file.is_open()) {
        return false;
    }

    std::string line;
    if (!std::getline(file, line)) {
        throw std::runtime_error("AeroCoefficientTable::load: " + csv_path + " is empty");
    }
    // Header line is read but not otherwise validated beyond existing --
    // column order is fixed by this class' documented format.

    std::vector<std::array<double, kNumAxes>> axis_cols;
    std::vector<AeroCoefficients> coeff_col;

    while (std::getline(file, line)) {
        if (line.find_first_not_of(" \t\r\n") == std::string::npos) continue;  // skip blank lines
        std::istringstream iss(line);
        std::string field;
        std::vector<double> vals;
        while (std::getline(iss, field, ',')) {
            vals.push_back(std::stod(field));
        }
        const std::size_t expected_cols = kNumAxes + kNumCoeffCols;
        if (vals.size() != expected_cols) {
            throw std::runtime_error("AeroCoefficientTable::load: malformed row in " + csv_path +
                                      " (expected " + std::to_string(expected_cols) + " columns, got " +
                                      std::to_string(vals.size()) + ")");
        }
        std::array<double, kNumAxes> axes;
        for (int a = 0; a < kNumAxes; ++a) axes[a] = vals[a];
        axis_cols.push_back(axes);

        AeroCoefficients c;
        c.CL = vals[kNumAxes + 0];
        c.CD = vals[kNumAxes + 1];
        c.Cl_roll = vals[kNumAxes + 2];
        c.Cm = vals[kNumAxes + 3];
        c.Cn_yaw = vals[kNumAxes + 4];
        coeff_col.push_back(c);
    }

    if (axis_cols.empty()) {
        throw std::runtime_error("AeroCoefficientTable::load: no data rows in " + csv_path);
    }

    grids_.assign(kNumAxes, {});
    for (int a = 0; a < kNumAxes; ++a) {
        std::vector<double> col;
        col.reserve(axis_cols.size());
        for (const auto& row : axis_cols) col.push_back(row[a]);
        grids_[a] = sortedUnique(col);
    }

    std::size_t expected = 1;
    for (int a = 0; a < kNumAxes; ++a) expected *= grids_[a].size();
    if (axis_cols.size() != expected) {
        throw std::runtime_error("AeroCoefficientTable::load: " + csv_path +
                                  " does not form a complete rectilinear grid (" +
                                  std::to_string(axis_cols.size()) + " rows, expected " +
                                  std::to_string(expected) + " for a full grid) -- "
                                  "check for missing or duplicate axis-value combinations");
    }

    values_.assign(expected, AeroCoefficients{});
    std::vector<bool> filled(expected, false);
    for (std::size_t row = 0; row < axis_cols.size(); ++row) {
        std::array<int, kNumAxes> idx;
        for (int a = 0; a < kNumAxes; ++a) {
            idx[a] = indexOf(grids_[a], axis_cols[row][a]);
        }
        const std::size_t flat = flatIndex(idx);
        if (filled[flat]) {
            throw std::runtime_error("AeroCoefficientTable::load: " + csv_path +
                                      " has a duplicate axis-value combination row");
        }
        values_[flat] = coeff_col[row];
        filled[flat] = true;
    }

    loaded_ = true;
    return true;
}

void AeroCoefficientTable::findBracket(const std::vector<double>& grid, double value, int& lo, int& hi, double& t) {
    if (grid.size() == 1) {
        lo = hi = 0;
        t = 0.0;
        return;
    }
    if (value <= grid.front()) {
        lo = hi = 0;
        t = 0.0;
        return;
    }
    if (value >= grid.back()) {
        lo = hi = static_cast<int>(grid.size()) - 1;
        t = 0.0;
        return;
    }
    auto it = std::upper_bound(grid.begin(), grid.end(), value);
    hi = static_cast<int>(it - grid.begin());
    lo = hi - 1;
    t = (value - grid[lo]) / (grid[hi] - grid[lo]);
}

std::size_t AeroCoefficientTable::flatIndex(const std::array<int, kNumAxes>& idx) const {
    std::size_t flat = static_cast<std::size_t>(idx[0]);
    for (int a = 1; a < kNumAxes; ++a) {
        flat = flat * grids_[a].size() + static_cast<std::size_t>(idx[a]);
    }
    return flat;
}

AeroCoefficients AeroCoefficientTable::interpolate(double mach, double alpha_deg, double beta_deg,
                                                     double fwd_sym_deg, double aft_sym_deg,
                                                     double aft_diff_deg) const {
    const double query[kNumAxes] = {mach, alpha_deg, beta_deg, fwd_sym_deg, aft_sym_deg, aft_diff_deg};

    std::array<int, kNumAxes> lo{}, hi{};
    std::array<double, kNumAxes> t{};
    for (int a = 0; a < kNumAxes; ++a) {
        findBracket(grids_[a], query[a], lo[a], hi[a], t[a]);
    }

    double CL = 0, CD = 0, Cl_roll = 0, Cm = 0, Cn_yaw = 0;
    const unsigned num_corners = 1u << kNumAxes;
    for (unsigned mask = 0; mask < num_corners; ++mask) {
        double w = 1.0;
        std::array<int, kNumAxes> idx{};
        for (int a = 0; a < kNumAxes; ++a) {
            const bool bit = (mask >> a) & 1u;
            idx[a] = bit ? hi[a] : lo[a];
            w *= bit ? t[a] : (1.0 - t[a]);
        }
        if (w == 0.0) continue;
        const AeroCoefficients& c = values_[flatIndex(idx)];
        CL += w * c.CL;
        CD += w * c.CD;
        Cl_roll += w * c.Cl_roll;
        Cm += w * c.Cm;
        Cn_yaw += w * c.Cn_yaw;
    }

    AeroCoefficients out;
    out.CL = CL;
    out.CD = CD;
    out.Cl_roll = Cl_roll;
    out.Cm = Cm;
    out.Cn_yaw = Cn_yaw;
    return out;
}

}  // namespace aero_model
