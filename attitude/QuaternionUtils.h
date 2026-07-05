#ifndef QUATERNION_UTILS_H
#define QUATERNION_UTILS_H

#include <Eigen/Geometry>

// Attitude quaternion convention used throughout this codebase: scalar-last
// (q1,q2,q3 vector part, q4 scalar part). Eigen::Quaterniond's constructor is
// (w,x,y,z), so callers must construct it as Quaterniond(q4, q1, q2, q3).
Eigen::Quaterniond QuaternionDerivative(const Eigen::Quaterniond& q,
                                        double wx, double wy, double wz);

#endif // QUATERNION_UTILS_H
