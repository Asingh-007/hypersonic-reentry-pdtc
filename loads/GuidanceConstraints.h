#ifndef GUIDANCE_CONSTRAINTS_H
#define GUIDANCE_CONSTRAINTS_H

#include "DescentDynamics.h"
#include <string>
#include <vector>

// PLACEHOLDER limits for an uncrewed hypersonic entry vehicle -- round,
// physically-plausible numbers pending real structural/TPS margins:
//  - max_heat_flux_w_m2 = 1.0e6 (1 MW/m^2): order-of-magnitude comparable to
//    Apollo-class lunar-return peak stagnation heating (~1-5 MW/m^2);
//    conservative for this vehicle's smaller nose radius (higher heating per
//    Sutton-Graves' 1/sqrt(R_n) scaling vs. a blunt capsule).
//  - max_qbar_pa = 50000.0 (50 kPa): common "max-q" order of magnitude for
//    launch/entry vehicles (cf. Shuttle ascent max-q ~30 kPa, rounded up).
//  - max_load_factor_g = 5.0: uncrewed structural limit, well under typical
//    crewed entry peaks (~4-8g) but high enough not to trivially trip.
// TODO: replace all three with real vehicle structural/TPS margins.
struct GuidanceConstraints {
    double max_heat_flux_w_m2 = 1.0e6;
    double max_qbar_pa = 50000.0;
    double max_load_factor_g = 5.0;
};

// One violation record: which sample index/time exceeded which limit, and by
// how much.
struct ConstraintViolation {
    size_t index = 0;
    double time_s = 0.0;
    std::string constraint_name; // "heat_flux" | "qbar" | "load_factor"
    double value = 0.0;          // actual value at this sample
    double limit = 0.0;          // the limit that was exceeded
    double margin = 0.0;         // value - limit (always > 0 here)
};

struct ConstraintViolationReport {
    std::vector<ConstraintViolation> violations;
    bool hasViolations() const { return !violations.empty(); }
};

// Scans every recorded point in hist against constraints, returning every
// (index, constraint, margin) exceedance -- a trajectory may violate
// multiple constraints at multiple times, all reported (not just the first).
inline ConstraintViolationReport checkConstraints(const DescentDynamics::TrajectoryHistory& hist,
                                                    const GuidanceConstraints& constraints) {
    ConstraintViolationReport report;
    for (size_t i = 0; i < hist.t.size(); ++i) {
        if (hist.heat_flux_total[i] > constraints.max_heat_flux_w_m2) {
            report.violations.push_back({i, hist.t[i], "heat_flux", hist.heat_flux_total[i],
                                          constraints.max_heat_flux_w_m2,
                                          hist.heat_flux_total[i] - constraints.max_heat_flux_w_m2});
        }
        if (hist.qbar[i] > constraints.max_qbar_pa) {
            report.violations.push_back({i, hist.t[i], "qbar", hist.qbar[i],
                                          constraints.max_qbar_pa,
                                          hist.qbar[i] - constraints.max_qbar_pa});
        }
        if (hist.load_factor[i] > constraints.max_load_factor_g) {
            report.violations.push_back({i, hist.t[i], "load_factor", hist.load_factor[i],
                                          constraints.max_load_factor_g,
                                          hist.load_factor[i] - constraints.max_load_factor_g});
        }
    }
    return report;
}

#endif // GUIDANCE_CONSTRAINTS_H
