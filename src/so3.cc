// SO3 FUNCTIONS

#include "so3.h"

#include <cmath>

#include <Eigen/SVD>

Eigen::Matrix3d Skew(const Eigen::Vector3d &w) {
  Eigen::Matrix3d W;
  W << 0.0, -w.z(), w.y(), w.z(), 0.0, -w.x(), -w.y(),  w.x(), 0.0;
  return W;
}

/*
Eigen::Matrix3d NormalizeRotation(const Eigen::Matrix3d &R) {
  Eigen::JacobiSVD<Eigen::Matrix3d> svd(R, Eigen::ComputeFullU | Eigen::ComputeFullV);
  return svd.matrixU()*svd.matrixV();
}
*/

Eigen::Matrix3d ExpSO3(const double x, const double y, const double z) {
  const double theta2 = x*x+y*y+z*z;
  const double theta  = std::sqrt(theta2);
  Eigen::Matrix3d W;
  W << 0.0, -z, y, z, 0.0, -x, -y,  x, 0.0;
  if (theta < 1e-5)
      return Eigen::Matrix3d::Identity() + W;// + 0.5*W*W;
  else
      return Eigen::Matrix3d::Identity() + W*std::sin(theta)/theta + W*W*(1.0-std::cos(theta))/theta2;
}

Eigen::Vector3d LogSO3(const Eigen::Matrix3d &R) {
  double costheta = 0.5*(R.trace()-1.0);
  if (costheta > +1.0) costheta = +1.0;
  if (costheta < -1.0) costheta = -1.0;
  const double theta = std::acos(costheta);
  const Eigen::Vector3d w(R(2,1)-R(1,2), R(0,2)-R(2,0), R(1,0)-R(0,1));
  if (theta < 1e-5)
      return 0.5*w;
  else
      return 0.5*theta*w/std::sin(theta);
}

Eigen::Matrix3d RightJacobianSO3(const double x, const double y, const double z) {
  const double theta2 = x*x+y*y+z*z;
  const double theta  = std::sqrt(theta2);

  Eigen::Matrix3d W;
  W << 0.0, -z, y, z, 0.0, -x, -y,  x, 0.0;
  if (theta < 1e-5)
      return Eigen::Matrix3d::Identity() - 0.5*W;// + 1.0/6.0*W*W;
  else
      return Eigen::Matrix3d::Identity() - W*(1.0-std::cos(theta))/theta2 + W*W*(theta-std::sin(theta))/(theta2*theta);
}

Eigen::Matrix3d InverseRightJacobianSO3(const double x, const double y, const double z) {
  const double theta2 = x*x+y*y+z*z;
  const double theta  = std::sqrt(theta2);

  Eigen::Matrix3d W;
  W << 0.0, -z, y, z, 0.0, -x, -y,  x, 0.0;
  if (theta < 1e-5)
      return Eigen::Matrix3d::Identity() + 0.5*W;// + 1.0/12.0*W*W;
  else
      return Eigen::Matrix3d::Identity() + 0.5*W + W*W*(1.0/theta2 - (1.0+std::cos(theta))/(2.0*theta*std::sin(theta)));
}
