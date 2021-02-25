
#ifndef SO3_H_
#define SO3_H_

#include <Eigen/Core>

Eigen::Matrix3d Skew(const Eigen::Vector3d &w);

//Eigen::Matrix3d NormalizeRotation(const Eigen::Matrix3d &R);

Eigen::Matrix3d ExpSO3(const double x, const double y, const double z);
inline Eigen::Matrix3d ExpSO3(const Eigen::Vector3d &w) { return ExpSO3(w[0], w[1], w[2]); }

Eigen::Vector3d LogSO3(const Eigen::Matrix3d &R);

Eigen::Matrix3d RightJacobianSO3(const double x, const double y, const double z);
inline Eigen::Matrix3d RightJacobianSO3(const Eigen::Vector3d &v) { return RightJacobianSO3(v[0], v[1], v[2]); }

Eigen::Matrix3d InverseRightJacobianSO3(const double x, const double y, const double z);
inline Eigen::Matrix3d InverseRightJacobianSO3(const Eigen::Vector3d &v) { return InverseRightJacobianSO3(v[0], v[1], v[2]); }

#endif  // SO3_H_
