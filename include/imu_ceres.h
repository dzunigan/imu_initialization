
#ifndef IMU_CERES_H
#define IMU_CERES_H

#include <cmath>
#include <memory>

#include <ceres/local_parameterization.h>
#include <ceres/sized_cost_function.h>

#include "imu_preintegration.h"
#include "so3.h"

#include "util/eigen_defs.h"

class ScaleParameterization : public ceres::LocalParameterization {
 public:
  virtual ~ScaleParameterization() {}
  bool Plus(const double* x,
            const double* delta,
            double* x_plus_delta) const override;
  bool ComputeJacobian(const double* x, double* jacobian) const override;
  int GlobalSize() const override { return 1; }
  int LocalSize() const override { return 1; }
};

// Implements x_plus_delta = x*exp(delta)
bool ScaleParameterization::Plus(const double* x,
                                 const double* delta,
                                 double* x_plus_delta) const {
  x_plus_delta[0] = x[0]*std::exp(delta[0]);
  return true;
}

bool ScaleParameterization::ComputeJacobian(const double* x, double* jacobian) const {
  jacobian[0] = x[0];
  return true;
}

class GravityParameterization : public ceres::LocalParameterization {
 public:
  virtual ~GravityParameterization() {}
  bool Plus(const double* x,
            const double* delta,
            double* x_plus_delta) const override;
  bool ComputeJacobian(const double* x, double* jacobian) const override;
  int GlobalSize() const override { return 9; }
  int LocalSize() const override { return 2; }
};

// Implements SO3 parameterization with delta_z = 0
bool GravityParameterization::Plus(const double* x,
                               const double* delta,
                               double* x_plus_delta) const {
  Eigen::Map<const Eigen::Matrix3d> R(x);
  Eigen::Matrix3d delta_R = ExpSO3(delta[0], delta[1], 0.0);
  Eigen::Map<Eigen::Matrix3d> result(x_plus_delta);
  result = R*delta_R;

  return true;
}

// Implements the derivative of SO3 while keeping the z coordinate
bool GravityParameterization::ComputeJacobian(const double* x, double* jacobian) const {
  Eigen::Map<const Eigen::Matrix3d> R(x);
  Eigen::Map<Eigen::Matrix<double, 9, 2, Eigen::RowMajor>> J(jacobian);
  J.block<3, 1>(0, 0).setZero();
  J.block<3, 1>(3, 0) =  R.block<3, 1>(0, 2);
  J.block<3, 1>(6, 0) = -R.block<3, 1>(0, 1);

  J.block<3, 1>(0, 1) = -R.block<3, 1>(0, 2);
  J.block<3, 1>(3, 1).setZero();
  J.block<3, 1>(6, 1) =  R.block<3, 1>(0, 0);

  return true;
}

// TODO Shouldn't the covariance matrix be (-invJr)*C*(-invJr)T ?
//      where C := pInt->C.block<3, 3>(0, 0)
//      and invJr := InverseRightJacobianSO3(er)
class GyroscopeBiasCostFunction : public ceres::SizedCostFunction<3, 3> {
 public:
  GyroscopeBiasCostFunction(std::shared_ptr<const IMU::Preintegrated> pInt, const Eigen::Matrix3d& Ri, const Eigen::Matrix3d& Rj)
    : pInt(pInt), Ri(Ri), Rj(Rj)
  {
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(pInt->C.block<3, 3>(0, 0));
    SqrtInformation = solver.operatorInverseSqrt();
  }
  virtual ~GyroscopeBiasCostFunction() { }

  bool Evaluate(double const* const* parameters,
                double* residuals,
                double** jacobians) const override {
    Eigen::Map<const Eigen::Vector3d> bg(parameters[0]);

    const Eigen::Matrix3d eR = pInt->GetDeltaRotation(bg).transpose()*Ri.transpose()*Rj;
    const Eigen::Vector3d er = LogSO3(eR);

    Eigen::Map<Eigen::Vector3d> e(residuals);
    e = er;
    e = SqrtInformation*e;

    if (jacobians != nullptr) {
      if (jacobians[0] != nullptr) {
        // wrt gyro bias
        const Eigen::Vector3d dbg = pInt->GetGyroDeltaBias(bg);
        const Eigen::Matrix3d invJr = InverseRightJacobianSO3(er);
        
        Eigen::Map<Eigen::Matrix<double, 3, 3, Eigen::RowMajor>> J(jacobians[0]);
        J = -invJr*eR.transpose()*RightJacobianSO3(pInt->JRg*dbg)*pInt->JRg;
        J = SqrtInformation*J;
      }
    }

    return true;
  }

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

 private:
  std::shared_ptr<const IMU::Preintegrated> pInt;

  const Eigen::Matrix3d Ri, Rj;

  Eigen::Matrix3d SqrtInformation;
};

// velocity1, velocity2, bias_g, bias_a, Rwg, scale
class InertialCostFunction : public ceres::SizedCostFunction<9, 3, 3, 3, 3, 9, 1> {
 public:
  InertialCostFunction(std::shared_ptr<const IMU::Preintegrated> pInt,
                       const Eigen::Matrix3d &R1, const Eigen::Vector3d &p1,
                       const Eigen::Matrix3d &R2, const Eigen::Vector3d &p2,
                       const Eigen::Isometry3d &Tcb = Eigen::Isometry3d::Identity())
    : pInt(pInt), dt(pInt->dT), R1_(R1), R2_(R2), p1_(p1), p2_(p2), Tcb(Tcb)
  {
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix9d> solver(pInt->C);
    SqrtInformation = solver.operatorInverseSqrt();
  }
  virtual ~InertialCostFunction() { }

  bool Evaluate(double const* const* parameters,
                double* residuals,
                double** jacobians) const override {
    Eigen::Map<const Eigen::Vector3d> v1(parameters[0]);
    Eigen::Map<const Eigen::Vector3d> v2(parameters[1]);
    Eigen::Map<const Eigen::Vector3d> bg(parameters[2]);
    Eigen::Map<const Eigen::Vector3d> ba(parameters[3]);
    Eigen::Map<const Eigen::Matrix3d> Rwg(parameters[4]);
    const double s = parameters[5][0];

    const Eigen::Vector3d g = Rwg*IMU::GRAVITY_VECTOR;

    const Eigen::Matrix3d dR = pInt->GetDeltaRotation(bg);
    const Eigen::Vector3d dV = pInt->GetDeltaVelocity(bg, ba);
    const Eigen::Vector3d dP = pInt->GetDeltaPosition(bg, ba);

    const Eigen::Matrix3d R1 = R1_*Tcb.linear();
    const Eigen::Matrix3d R2 = R2_*Tcb.linear();

    const Eigen::Matrix3d eR = dR.transpose()*R1.transpose()*R2;
    const Eigen::Vector3d er = LogSO3(eR);

    Eigen::Map<Eigen::Vector9d> e(residuals);
    e.head<3>()     = er;
    e.segment<3>(3) = R1.transpose()*(s*(v2 - v1) - g*dt) - dV;
    e.tail<3>()     = R1.transpose()*(s*(p2_ - p1_ - v1*dt) + (R2_-R1_)*Tcb.translation() - 0.5*g*dt*dt) - dP;
    
    e = SqrtInformation*e;

    if (jacobians != nullptr) {
      const Eigen::Vector3d dbg = pInt->GetGyroDeltaBias(bg);
      const Eigen::Matrix3d invJr = InverseRightJacobianSO3(er);
      if (jacobians[0] != nullptr) {
        // wrt velocity1
        Eigen::Map<Eigen::Matrix<double, 9, 3, Eigen::RowMajor>> J(jacobians[0]);
        J.block<3, 3>(0, 0).setZero();
        J.block<3, 3>(3, 0) = -s*R1.transpose();
        J.block<3, 3>(6, 0) = -s*dt*R1.transpose();
        J = SqrtInformation*J;
      }
      if (jacobians[1] != nullptr) {
        // wrt velocity2
        Eigen::Map<Eigen::Matrix<double, 9, 3, Eigen::RowMajor>> J(jacobians[1]);
        J.block<3, 3>(0, 0).setZero();
        J.block<3, 3>(3, 0) = s*R1.transpose();
        J.block<3, 3>(6, 0).setZero();
        J = SqrtInformation*J;
      }
      if (jacobians[2] != nullptr) {
        // wrt gyro bias
        Eigen::Map<Eigen::Matrix<double, 9, 3, Eigen::RowMajor>> J(jacobians[2]);
        J.block<3, 3>(0, 0) = -invJr*eR.transpose()*RightJacobianSO3(pInt->JRg*dbg)*pInt->JRg;
        J.block<3, 3>(3, 0) = -pInt->JVg;
        J.block<3, 3>(6, 0) = -pInt->JPg;
        J = SqrtInformation*J;
      }
      if (jacobians[3] != nullptr) {
        // wrt acc bias
        Eigen::Map<Eigen::Matrix<double, 9, 3, Eigen::RowMajor>> J(jacobians[3]);
        J.block<3, 3>(0, 0).setZero();
        J.block<3, 3>(3, 0) = -pInt->JVa;
        J.block<3, 3>(6, 0) = -pInt->JPa;
        J = SqrtInformation*J;
      }
      if (jacobians[4] != nullptr) {
        // wrt Rwg
        Eigen::Map<Eigen::Matrix<double, 9, 9, Eigen::RowMajor>> J(jacobians[4]);
        J.setZero();
        J.block<3, 3>(3, 6) = dt*IMU::GRAVITY_MAGNITUDE*R1.transpose();
        J.block<3, 3>(6, 6) = 0.5*dt*dt*IMU::GRAVITY_MAGNITUDE*R1.transpose();
        J = SqrtInformation*J;
      }
      if (jacobians[5] != nullptr) {
        // wrt scale
        Eigen::Map<Eigen::Vector9d> J(jacobians[5]);
        J.block<3, 1>(0, 0).setZero();
        J.block<3, 1>(3, 0) = R1.transpose()*(v2 - v1);
        J.block<3, 1>(6, 0) = R1.transpose()*(p2_ - p1_ - v1*dt);
        J = SqrtInformation*J;
      }
    }

    return true;
  }

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

 private:
  std::shared_ptr<const IMU::Preintegrated> pInt;
  const double dt;

  const Eigen::Matrix3d R1_, R2_;
  const Eigen::Vector3d p1_, p2_;
  const Eigen::Isometry3d Tcb;
  
  Eigen::Matrix9d SqrtInformation;
};

// bias
class BiasPriorCostFunction : public ceres::SizedCostFunction<3, 3> {
 public:
  BiasPriorCostFunction(const double variance, const Eigen::Vector3d &mean = Eigen::Vector3d::Zero())
    : weight(std::sqrt(variance)), mean(mean) { }
  virtual ~BiasPriorCostFunction() { }

  bool Evaluate(double const* const* parameters,
                double* residuals,
                double** jacobians) const override {
    Eigen::Map<const Eigen::Vector3d> bias(parameters[0]);

    Eigen::Map<Eigen::Vector3d> error(residuals);
    error = weight*(mean - bias);

    if (jacobians != nullptr) {
      if (jacobians[0] != nullptr) {
        Eigen::Map<Eigen::Matrix<double, 3, 3, Eigen::RowMajor>> J(jacobians[0]);
        J = -weight*Eigen::Matrix3d::Identity();
      }
    }

    return true;
  }

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

 private:
  const double weight;
  const Eigen::Vector3d mean;
};

#endif // IMU_CERES_H
