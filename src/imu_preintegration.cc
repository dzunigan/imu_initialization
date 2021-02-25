/**
* This file is part of ORB-SLAM3
*
* Copyright (C) 2017-2020 Carlos Campos, Richard Elvira, Juan J. Gómez Rodríguez, José M.M. Montiel and Juan D. Tardós, University of Zaragoza.
* Copyright (C) 2014-2016 Raúl Mur-Artal, José M.M. Montiel and Juan D. Tardós, University of Zaragoza.
*
* ORB-SLAM3 is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
* License as published by the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* ORB-SLAM3 is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
* the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License along with ORB-SLAM3.
* If not, see <http://www.gnu.org/licenses/>.
*/

#include "imu_preintegration.h"

#include "so3.h"

#include "util/eigen_defs.h"

namespace IMU {

namespace {
  //Integration of 1 gyro measurement
  class IntegratedRotation {
   public:
    IntegratedRotation() = delete;
    IntegratedRotation(const Eigen::Vector3d& w, const Eigen::Vector3d& bias, const double time)
      : deltaT(time)
    {
      const Eigen::Vector3d dr = time*(w-bias);    
      deltaR = ExpSO3(dr.x(), dr.y(), dr.z());
      rightJ = RightJacobianSO3(dr.x(), dr.y(), dr.z());
    }

  public:
      const double deltaT; //integration time
      Eigen::Matrix3d deltaR; //integrated rotation
      Eigen::Matrix3d rightJ; // right jacobian
  };
}

Preintegrated::Preintegrated(const Eigen::Vector3d& bg, const Eigen::Vector3d& ba) {
  //NgaWalk = Preintegrated::SigmaW;
  Initialize(bg, ba);
}

void Preintegrated::IntegrateNewMeasurement(const Eigen::Vector3d& w,
                                            const Eigen::Vector3d& a,
                                            const double dt) {
  // Position is updated first, as it depends on previously computed velocity and rotation.
  // Velocity is updated secondly, as it depends on previously computed rotation.
  // Rotation is the last to be updated.

  //Matrices to compute covariance
  Eigen::Matrix9d A;
  A.setIdentity();
  
  Eigen::Matrix<double, 9, 6> B;
  B.setZero();

  Eigen::Vector3d acc = a - b.tail<3>();
  //Eigen::Vector3d accW(angVel.x()-b.bwx, angVel.y()-b.bwy, angVel.z()-b.bwz);

  //avgA = (dT*avgA + dR*acc*dt)/(dT+dt);
  //avgW = (dT*avgW + accW*dt)/(dT+dt);

  // Update delta position dP and velocity dV (rely on no-updated delta rotation)
  dP += dV*dt + 0.5*dR*acc*dt*dt;
  dV += dR*acc*dt;

  // Compute velocity and position parts of matrices A and B (rely on non-updated delta rotation)
  Eigen::Matrix3d Wacc = Skew(acc);

  A.block<3, 3>(3, 0) = -dR*dt*Wacc;
  A.block<3, 3>(6, 0) = -0.5*dR*dt*dt*Wacc;
  A.block<3, 3>(6, 3) = Eigen::Matrix3d::Identity()*dt;
  B.block<3, 3>(3, 3) = dR*dt;
  B.block<3, 3>(6, 3) = 0.5*dR*dt*dt;

  // Update position and velocity jacobians wrt bias correction
  JPa = JPa + JVa*dt -0.5*dR*dt*dt;
  JPg = JPg + JVg*dt -0.5*dR*dt*dt*Wacc*JRg;
  JVa = JVa - dR*dt;
  JVg = JVg - dR*dt*Wacc*JRg;

  // Update delta rotation
  IntegratedRotation dRi(w, b.head<3>(), dt);
  dR *= dRi.deltaR;

  // Compute rotation parts of matrices A and B
  A.block<3, 3>(0, 0) = dRi.deltaR.transpose();
  B.block<3, 3>(0, 0) = dRi.rightJ*dt;

  // Update covariance
  C = A*C*A.transpose() + B*Sigma*B.transpose();
  //C.block<6, 6>(9, 9) += NgaWalk;

  // Update rotation jacobian wrt bias correction
  JRg = dRi.deltaR.transpose()*JRg - dRi.rightJ*dt;

  // Total integrated time
  dT += dt;
}

void Preintegrated::SetNewGyroBias(const Eigen::Vector3d& bg) {
  bu.head<3>() = bg;
  db.head<3>() = bg - b.head<3>();
}

void Preintegrated::SetNewAccBias(const Eigen::Vector3d& ba) {
  bu.tail<3>() = ba;
  db.tail<3>() = ba - b.tail<3>();
}

Eigen::Vector3d Preintegrated::GetGyroDeltaBias() const {
  return db.head<3>();
}

Eigen::Vector3d Preintegrated::GetGyroDeltaBias(const Eigen::Vector3d& bg) const {
  return bg-b.head<3>();
}

Eigen::Vector3d Preintegrated::GetGyroOriginalBias() const {
  return b.head<3>();
}

Eigen::Vector3d Preintegrated::GetGyroUpdatedBias() const {
  return bu.head<3>();
}

Eigen::Vector3d Preintegrated::GetAccDeltaBias() const {
  return db.tail<3>();
}

Eigen::Vector3d Preintegrated::GetAccDeltaBias(const Eigen::Vector3d& ba) const {
  return ba-b.tail<3>();
}

Eigen::Vector3d Preintegrated::GetAccOriginalBias() const {
  return b.tail<3>();
}

Eigen::Vector3d Preintegrated::GetAccUpdatedBias() const {
  return bu.tail<3>();
}

Eigen::Matrix3d Preintegrated::GetDeltaRotation(const Eigen::Vector3d& bg) const {
  return dR*ExpSO3(JRg*(bg-b.head<3>()));
}

Eigen::Vector3d Preintegrated::GetDeltaVelocity(const Eigen::Vector3d& bg, const Eigen::Vector3d& ba) const {
    return dV + JVg*(bg-b.head<3>()) + JVa*(ba-b.tail<3>());
}

Eigen::Vector3d Preintegrated::GetDeltaPosition(const Eigen::Vector3d& bg, const Eigen::Vector3d& ba) const {
  return dP + JPg*(bg-b.head<3>()) + JPa*(ba-b.tail<3>());
}

Eigen::Matrix3d Preintegrated::GetUpdatedDeltaRotation() const {
  return dR*ExpSO3(JRg*db.head<3>());
}

Eigen::Vector3d Preintegrated::GetUpdatedDeltaVelocity() const {
  return dV + JVg*db.head<3>() + JVa*db.tail<3>();
}

Eigen::Vector3d Preintegrated::GetUpdatedDeltaPosition() const {
  return dP + JPg*db.head<3>() + JPa*db.tail<3>();
}

Eigen::Matrix3d Preintegrated::GetOriginalDeltaRotation() const {
  return dR;
}

Eigen::Vector3d Preintegrated::GetOriginalDeltaVelocity() const {
  return dV;
}

Eigen::Vector3d Preintegrated::GetOriginalDeltaPosition() const {
  return dP;
}

void Preintegrated::Initialize(const Eigen::Vector3d& bg, const Eigen::Vector3d& ba) {
  dT = 0.0;
  C.setZero();
  
  b.head<3>() = bg;
  b.tail<3>() = ba;
  dR.setIdentity();
  dV.setZero();
  dP.setZero();
  JRg.setZero();
  JVg.setZero();
  JVa.setZero();
  JPg.setZero();
  JPa.setZero();
  //avgA.setZero();
  //avgW.setZero();
  
  bu = b;
  db.setZero();
}

} // namespace IMU
