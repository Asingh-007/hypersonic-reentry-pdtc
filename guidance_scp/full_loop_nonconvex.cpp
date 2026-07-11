#include "full_loop_nonconvex.h"

#include <cmath>

namespace guidance_scp {

ConeRowGroup buildThrustMagnitudeSocRows(int sigma_u_col, int ux_col, int uy_col, int uz_col,
                                           double sigma_u_ref, const Eigen::Vector3d& u_vec_ref) {
    ConeRowGroup g;
    g.rows.push_back({{{sigma_u_col, -1.0}}, sigma_u_ref});
    g.rows.push_back({{{ux_col, -1.0}}, u_vec_ref.x()});
    g.rows.push_back({{{uy_col, -1.0}}, u_vec_ref.y()});
    g.rows.push_back({{{uz_col, -1.0}}, u_vec_ref.z()});
    g.cone = clarabel::SecondOrderConeT<double>(4);
    return g;
}

ConeRowGroup buildGimbalConeSocRows(int ux_col, int uy_col, int uz_col,
                                      const Eigen::Vector3d& e1_ref, double deltaE_max_rad,
                                      const Eigen::Vector3d& u_vec_ref) {
    const double inv_cos = 1.0 / std::cos(deltaE_max_rad);
    ConeRowGroup g;
    g.rows.push_back({{{ux_col, -inv_cos * e1_ref.x()}, {uy_col, -inv_cos * e1_ref.y()},
                        {uz_col, -inv_cos * e1_ref.z()}}, inv_cos * u_vec_ref.dot(e1_ref)});
    g.rows.push_back({{{ux_col, -1.0}}, u_vec_ref.x()});
    g.rows.push_back({{{uy_col, -1.0}}, u_vec_ref.y()});
    g.rows.push_back({{{uz_col, -1.0}}, u_vec_ref.z()});
    g.cone = clarabel::SecondOrderConeT<double>(4);
    return g;
}

std::vector<ConeRowGroup> buildSymmetricBoxRows(int var_col, double var_ref, double bound) {
    // var_ref + delta_var <= bound  =>  delta_var <= bound - var_ref
    ConeRowGroup upper;
    upper.rows.push_back({{{var_col, 1.0}}, bound - var_ref});
    upper.cone = clarabel::NonnegativeConeT<double>(1);

    // var_ref + delta_var >= -bound  =>  -delta_var <= bound + var_ref
    ConeRowGroup lower;
    lower.rows.push_back({{{var_col, -1.0}}, bound + var_ref});
    lower.cone = clarabel::NonnegativeConeT<double>(1);

    return {upper, lower};
}

}  // namespace guidance_scp
