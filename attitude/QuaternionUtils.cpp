#include "QuaternionUtils.h"

Eigen::Quaterniond QuaternionDerivative(const Eigen::Quaterniond& q,
                                        double wx, double wy, double wz) {
    // qdot = 0.5 * q (x) [0, wx, wy, wz]  (Hamilton product form of qdot = 0.5*Omega(w)*q,
    // for a body-to-inertial quaternion with rates expressed in body axes).
    // Order matters here: q * omega, not omega * q.
    Eigen::Quaterniond omega(0.0, wx, wy, wz); // w()=0, x()=wx, y()=wy, z()=wz
    Eigen::Quaterniond qdot = q * omega; // Hamilton product
    qdot.coeffs() *= 0.5;
    return qdot;
}
