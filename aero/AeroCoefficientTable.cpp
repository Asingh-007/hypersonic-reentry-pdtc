#include "AeroCoefficientTable.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace aero_model {

namespace {

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

    std::vector<double> mach_col, alpha_col, beta_col, flap_col;
    std::vector<AeroCoefficients> coeff_col;

    while (std::getline(file, line)) {
        if (line.find_first_not_of(" \t\r\n") == std::string::npos) continue; // skip blank lines
        std::istringstream iss(line);
        std::string field;
        std::vector<double> vals;
        while (std::getline(iss, field, ',')) {
            vals.push_back(std::stod(field));
        }
        if (vals.size() != 9) {
            throw std::runtime_error("AeroCoefficientTable::load: malformed row in " + csv_path +
                                      " (expected 9 columns, got " + std::to_string(vals.size()) + ")");
        }
        mach_col.push_back(vals[0]);
        alpha_col.push_back(vals[1]);
        beta_col.push_back(vals[2]);
        flap_col.push_back(vals[3]);
        AeroCoefficients c;
        c.CL = vals[4];
        c.CD = vals[5];
        c.Cl_roll = vals[6];
        c.Cm = vals[7];
        c.Cn_yaw = vals[8];
        coeff_col.push_back(c);
    }

    if (mach_col.empty()) {
        throw std::runtime_error("AeroCoefficientTable::load: no data rows in " + csv_path);
    }

    mach_grid_ = sortedUnique(mach_col);
    alpha_grid_ = sortedUnique(alpha_col);
    beta_grid_ = sortedUnique(beta_col);
    flap_grid_ = sortedUnique(flap_col);

    const std::size_t expected = mach_grid_.size() * alpha_grid_.size() * beta_grid_.size() * flap_grid_.size();
    if (mach_col.size() != expected) {
        throw std::runtime_error("AeroCoefficientTable::load: " + csv_path +
                                  " does not form a complete rectilinear grid (" +
                                  std::to_string(mach_col.size()) + " rows, expected " +
                                  std::to_string(expected) + " for a full grid) -- "
                                  "check for missing or duplicate (mach,alpha,beta,flap) combinations");
    }

    values_.assign(expected, AeroCoefficients{});
    std::vector<bool> filled(expected, false);
    for (std::size_t row = 0; row < mach_col.size(); ++row) {
        int im = indexOf(mach_grid_, mach_col[row]);
        int ia = indexOf(alpha_grid_, alpha_col[row]);
        int ib = indexOf(beta_grid_, beta_col[row]);
        int ifl = indexOf(flap_grid_, flap_col[row]);
        std::size_t idx = ((static_cast<std::size_t>(im) * alpha_grid_.size() + ia) * beta_grid_.size() + ib) *
                           flap_grid_.size() + ifl;
        if (filled[idx]) {
            throw std::runtime_error("AeroCoefficientTable::load: " + csv_path +
                                      " has a duplicate (mach,alpha,beta,flap) row");
        }
        values_[idx] = coeff_col[row];
        filled[idx] = true;
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

const AeroCoefficients& AeroCoefficientTable::at(int im, int ia, int ib, int ifl) const {
    std::size_t idx = ((static_cast<std::size_t>(im) * alpha_grid_.size() + ia) * beta_grid_.size() + ib) *
                       flap_grid_.size() + ifl;
    return values_[idx];
}

AeroCoefficients AeroCoefficientTable::interpolate(double mach, double alpha_deg, double beta_deg, double flap_deg) const {
    int im0, im1, ia0, ia1, ib0, ib1, ifl0, ifl1;
    double tm, ta, tb, tfl;
    findBracket(mach_grid_, mach, im0, im1, tm);
    findBracket(alpha_grid_, alpha_deg, ia0, ia1, ta);
    findBracket(beta_grid_, beta_deg, ib0, ib1, tb);
    findBracket(flap_grid_, flap_deg, ifl0, ifl1, tfl);

    double CL = 0, CD = 0, Cl_roll = 0, Cm = 0, Cn_yaw = 0;
    const int im_idx[2] = {im0, im1};
    const int ia_idx[2] = {ia0, ia1};
    const int ib_idx[2] = {ib0, ib1};
    const int ifl_idx[2] = {ifl0, ifl1};
    const double m_w[2] = {1.0 - tm, tm};
    const double a_w[2] = {1.0 - ta, ta};
    const double b_w[2] = {1.0 - tb, tb};
    const double f_w[2] = {1.0 - tfl, tfl};

    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            for (int k = 0; k < 2; ++k) {
                for (int l = 0; l < 2; ++l) {
                    double w = m_w[i] * a_w[j] * b_w[k] * f_w[l];
                    if (w == 0.0) continue;
                    const AeroCoefficients& c = at(im_idx[i], ia_idx[j], ib_idx[k], ifl_idx[l]);
                    CL += w * c.CL;
                    CD += w * c.CD;
                    Cl_roll += w * c.Cl_roll;
                    Cm += w * c.Cm;
                    Cn_yaw += w * c.Cn_yaw;
                }
            }
        }
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
