
#ifndef METHODS_H_
#define METHODS_H_

// STL
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <vector>

// Ceres
#include <ceres/ceres.h>

//Eigen
#include <Eigen/Core>
#include <Eigen/Geometry>

// Glog
#include <glog/logging.h>

#include "imu_ceres.h"
#include "imu_preintegration.h"
#include "polynomial.h"
#include "so3.h"

#include "util/svd.h"
#include "util/timer.h"

struct input_t {
  input_t(const Eigen::Isometry3d &T1, const std::uint64_t t1,
          const Eigen::Isometry3d &T2, const std::uint64_t t2,
          std::shared_ptr<IMU::Preintegrated> pInt)
    : t1(t1), t2(t2), T1(T1), T2(T2), pInt(pInt)
  { }

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  const std::uint64_t t1, t2;
  Eigen::Isometry3d T1, T2;
  std::shared_ptr<IMU::Preintegrated> pInt;
  
  std::vector<IMU::Measurement> vMeasurements;
};

struct result_t {

  result_t()
    : success(false)
  { }

  result_t(bool success, std::int64_t solve_ns, double scale,
          const Eigen::Vector3d &bias_g, const Eigen::Vector3d &bias_a,
          const Eigen::Vector3d &gravity)
    : success(success), solve_ns(solve_ns), scale(scale),
      bias_g(bias_g), bias_a(bias_a), gravity(gravity)
  { }

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  bool success;
  std::int64_t solve_ns, velocities_ns;
  double scale;
  Eigen::Vector3d bias_g, bias_a, gravity;
};

using InputType = std::vector<input_t>;
using ResultType = result_t;

void gyroscope_only(const InputType &input,
                    ResultType &result,
                    const Eigen::Matrix3d &Rcb = Eigen::Matrix3d::Identity(), bool use_covarinace = true) {
  double** parameters = new double*[1];
  parameters[0] = new double[3];
  Eigen::Map<Eigen::Vector3d> bias_(parameters[0]);
  bias_.setZero();

  ceres::Problem problem;
  for (unsigned i = 0; i < input.size(); ++i) {
    const Eigen::Isometry3d &T1 = input[i].T1;
    const Eigen::Isometry3d &T2 = input[i].T2;
    const std::shared_ptr<IMU::Preintegrated> pInt = input[i].pInt;

    ceres::CostFunction* cost_function = new GyroscopeBiasCostFunction(pInt, T1.linear()*Rcb, T2.linear()*Rcb, use_covarinace);
    problem.AddResidualBlock(cost_function, nullptr, parameters, 1);
  }

  Timer timer;
  timer.Start();

  ceres::Solver::Options options;
  ceres::Solver::Summary summary;
  ceres::Solve(options, &problem, &summary);

  result.solve_ns = timer.ElapsedNanoSeconds();

  bool converged = (summary.termination_type == ceres::CONVERGENCE);
  if (converged) {
    result.success = true;
    result.bias_g = bias_;
  } else {
    LOG(ERROR) << summary.FullReport();
    result.success = false;
  }

  delete[] parameters[0];
  delete[] parameters;
}

Eigen::VectorXd real_roots(const Eigen::VectorXd &real, const Eigen::VectorXd &imag) {
  CHECK_EQ(real.size(), imag.size());

  Eigen::VectorXd roots(real.size());

	Eigen::VectorXd::Index j = 0;
	for (Eigen::VectorXd::Index i = 0; i < real.size(); ++i) {
	  if (!imag(i)) {
	    roots(j) = real(i);
	    ++j;
	  }
	}

	roots.conservativeResize(j);
	return roots;
}

void analytic_accelerometer(const InputType &input, ResultType &result,
                            const Eigen::Vector3d &bg = Eigen::Vector3d::Zero(),
                            const Eigen::Vector3d &ba = Eigen::Vector3d::Zero(),
                            const Eigen::Isometry3d &Tcb = Eigen::Isometry3d::Identity(),
                            const double prior = 0.0) {
  //LOG(INFO) << "Running proposed method at " << input[0].t1;

  CHECK_GE(prior, 0.0);

  constexpr int n = 7;
  constexpr int q = 4;

  Eigen::MatrixXd M(n, n);
  M.setZero();

  Eigen::VectorXd m(n);
  m.setZero();

  double Q = 0.;
  
  Timer timer;
  timer.Start();

  const Eigen::Vector3d ba_prior_mean = Eigen::Vector3d::Zero();
  
  // accelerometer bias prior
  {    
    Eigen::MatrixXd M_k(3, n);
    M_k.setZero();
    
    M_k.block<3, 3>(0, 1) = Eigen::Matrix3d::Identity();
    
    Eigen::Vector3d pi_k;
    pi_k = ba_prior_mean;
    
    Eigen::Matrix3d Information = prior*Eigen::Matrix3d::Identity();

    M +=  M_k.transpose()*Information*M_k;
    m += -2.*M_k.transpose()*Information*pi_k;
    Q +=  pi_k.transpose()*Information*pi_k;
  }

  for (unsigned i = 1; i < input.size(); ++i) {
    //CHECK_EQ(input[i-1].t2, input[i].t1);
    const Eigen::Isometry3d &T1 = input[i-1].T1;
    const Eigen::Isometry3d &T2 = input[i].T1;
    const Eigen::Isometry3d &T3 = input[i].T2;
    const IMU::Preintegrated &pInt12 = *(input[i-1].pInt);
    const IMU::Preintegrated &pInt23 = *(input[i].pInt);

    Eigen::Matrix3d R1 = T1.linear()*Tcb.linear(); // Rwb
    Eigen::Matrix3d R2 = T2.linear()*Tcb.linear(); // Rwb

    Eigen::Matrix3d A = R1/pInt12.dT;
    Eigen::Matrix3d B = R2/pInt23.dT;

    Eigen::MatrixXd M_k(3, n);
    M_k.setZero();

    M_k.col(0) = (T3.translation() - T2.translation())/pInt23.dT
                  - (T2.translation() - T1.translation())/pInt12.dT;
    M_k.block<3, 3>(0, 1) = A*pInt12.JPa - B*pInt23.JPa - R1*pInt12.JVa;
    M_k.block<3, 3>(0, q) = -0.5*(pInt12.dT + pInt23.dT)*Eigen::Matrix3d::Identity();

    Eigen::Vector3d pi_k;
    pi_k = B*pInt23.GetDeltaPosition(bg, ba) - A*pInt12.GetDeltaPosition(bg, ba) + R1*pInt12.GetDeltaVelocity(bg, ba)
            + (T2.linear() - T1.linear())*Tcb.translation()/pInt12.dT
            - (T3.linear() - T2.linear())*Tcb.translation()/pInt23.dT;

    Eigen::Matrix3d Covariance;
    Covariance  = A*pInt12.C.block<3, 3>(6, 6)*A.transpose();
    Covariance += B*pInt23.C.block<3, 3>(6, 6)*B.transpose();
    Covariance += T1.linear()*pInt12.C.block<3, 3>(3, 3)*T1.linear().transpose();

    Eigen::Matrix3d Information = selfAdjointInverse(Covariance);
    //Eigen::Matrix3d Information = Eigen::Matrix3d::Identity();

    M +=  M_k.transpose()*Information*M_k;
    m += -2.*M_k.transpose()*Information*pi_k;
    Q +=  pi_k.transpose()*Information*pi_k;
  }

  // Solve
  Eigen::Matrix4d A = 2.*M.block<4, 4>(0, 0);
  //LOG(INFO) << StringPrintf("A: %.16f", A);
  
  // TODO Check if A is invertible!!
  //Eigen::Matrix3d A_ = A.block<3, 3>(1, 1);
  //Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> svdA_(A_, Eigen::EigenvaluesOnly);
  //result.svA_ = svdA_.eigenvalues();
  //result.detA_ = A_.determinant();

  Eigen::MatrixXd Bt = 2.*M.block<3, 4>(q, 0);
  Eigen::MatrixXd BtAi = Bt*A.inverse();

  Eigen::Matrix3d D = 2.*M.block<3, 3>(q, q);
  Eigen::Matrix3d S = D - BtAi*Bt.transpose();

  // TODO Check if S is invertible!
  //Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> svdS(S, Eigen::EigenvaluesOnly);
  //result.svS = svdS.eigenvalues();
  //result.detS = S.determinant();
  //LOG(INFO) << StringPrintf("det(S): %.16f", S.determinant());
  //LOG(INFO) << StringPrintf("eigenvalues(S): %.16f %.16f %.16f",
  //                          c[0], svd.eigenvalues()[1], svd.eigenvalues()[2]);

  Eigen::Matrix3d Sa = S.determinant()*S.inverse();
  Eigen::Matrix3d U = S.trace()*Eigen::Matrix3d::Identity() - S;

  Eigen::Vector3d v1 = BtAi*m.head<q>();
  Eigen::Vector3d m2 = m.tail<3>();

  Eigen::Matrix3d X; Eigen::Vector3d Xm2;

  // X = I
  const double c4 = 16.*(v1.dot(  v1) - 2.*v1.dot( m2) + m2.dot( m2));

  X = U; Xm2 = X*m2;
  const double c3 = 16.*(v1.dot(X*v1) - 2.*v1.dot(Xm2) + m2.dot(Xm2));

  X = 2.*Sa + U*U; Xm2 = X*m2;
  const double c2 =  4.*(v1.dot(X*v1) - 2.*v1.dot(Xm2) + m2.dot(Xm2));

  X = Sa*U + U*Sa; Xm2 = X*m2;
  const double c1 =  2.*(v1.dot(X*v1) - 2.*v1.dot(Xm2) + m2.dot(Xm2));

  X = Sa*Sa; Xm2 = X*m2;
  const double c0 =     (v1.dot(X*v1) - 2.*v1.dot(Xm2) + m2.dot(Xm2));

  const double s00 = S(0, 0), s01 = S(0, 1), s02 = S(0, 2);
  const double s11 = S(1, 1), s12 = S(1, 2), s22 = S(2, 2);

  const double t1 = s00 + s11 + s22;
  const double t2 = s00*s11 + s00*s22 + s11*s22
                     - std::pow(s01, 2) - std::pow(s02, 2) - std::pow(s12, 2);
  const double t3 = s00*s11*s22 + 2.*s01*s02*s12
                     - s00*std::pow(s12, 2) - s11*std::pow(s02, 2) - s22*std::pow(s01, 2);

  Eigen::VectorXd coeffs(7);
  coeffs << 64.,
            64.*t1,
            16.*(std::pow(t1, 2) + 2.*t2),
            16.*(t1*t2 + t3),
             4.*(std::pow(t2, 2) + 2.*t1*t3),
             4.*t3*t2,
            std::pow(t3, 2);

  const double G2i = 1. / std::pow(IMU::GRAVITY_MAGNITUDE, 2);

  coeffs(2) -= c4*G2i;
  coeffs(3) -= c3*G2i;
  coeffs(4) -= c2*G2i;
  coeffs(5) -= c1*G2i;
  coeffs(6) -= c0*G2i;

  Eigen::VectorXd real, imag;
  if (!FindPolynomialRootsCompanionMatrix(coeffs, &real, &imag)) {
    LOG(ERROR) << "Failed to find roots\n"
               << StringPrintf("%.16f %.16f %.16f %.16f %.16f %.16f %.16f",
                               coeffs[0], coeffs[1], coeffs[2], coeffs[3],
                               coeffs[4], coeffs[5], coeffs[6]);
    result.success = false;
    result.solve_ns = timer.ElapsedNanoSeconds();
    return;
  }

  Eigen::VectorXd lambdas = real_roots(real, imag);
  if (lambdas.size() == 0) {
    LOG(ERROR) << "No real roots found\n"
               << StringPrintf("%.16f %.16f %.16f %.16f %.16f %.16f %.16f",
                               coeffs[0], coeffs[1], coeffs[2], coeffs[3],
                               coeffs[4], coeffs[5], coeffs[6]);
    result.success = false;
    result.solve_ns = timer.ElapsedNanoSeconds();
    return;
  }

  Eigen::MatrixXd W(n, n);
  W.setZero();
  W.block<3, 3>(q, q) = Eigen::Matrix3d::Identity();

  Eigen::VectorXd solution;
  double min_cost = std::numeric_limits<double>::max();
  for (Eigen::VectorXd::Index i = 0; i < lambdas.size(); ++i) {
    const double lambda = lambdas(i);

    Eigen::FullPivLU<Eigen::MatrixXd> lu(2.*M + 2.*lambda*W);
    Eigen::VectorXd x_ = -lu.inverse()*m;

    double cost = x_.transpose()*M*x_;
    cost += m.transpose()*x_;
    cost += Q;

    if (cost < min_cost) {
      solution = x_;
      min_cost = cost;
    }
  }

  result.solve_ns = timer.ElapsedNanoSeconds();

  const double constraint = solution.transpose()*W*solution;
  if (solution[0] < 1e-3 || constraint < 0.
      || std::abs(std::sqrt(constraint) - IMU::GRAVITY_MAGNITUDE)/IMU::GRAVITY_MAGNITUDE > 1e-3) { // TODO
    LOG(WARNING) << "Discarding bad solution...\n"
                 << StringPrintf("scale: %.16f\n", solution[0])
                 << StringPrintf("constraint: %.16f\n", constraint)
                 << StringPrintf("constraint error: %.2f %", 100.*std::abs(std::sqrt(constraint) - IMU::GRAVITY_MAGNITUDE)/IMU::GRAVITY_MAGNITUDE);
    result.success = false;
    return;
  }
  
  result.success = true;
  result.scale = solution[0];
  result.bias_a = solution.segment<3>(1);
  result.gravity = solution.segment<3>(4);
  
  const int N = input.size() + 1;
  Eigen::VectorXd velocities(3*N);
  
  timer.Start();
  
  // Recover velocities
  for (unsigned int i = 0; i < input.size(); ++i) {
    const Eigen::Isometry3d &T1 = input[i].T1;
    const Eigen::Isometry3d &T2 = input[i].T2;
    const IMU::Preintegrated &pInt12 = *(input[i].pInt);
    
    // -[(g*dt12^2)/2 + R1*dP12 + R1_c*pcb - R2_c*pcb + p1_c*s - p2_c*s]/dt12
    velocities.segment<3>(3*i) = (-0.5*result.gravity*std::pow(pInt12.dT, 2)
                                  - T1.linear()*Tcb.linear()*pInt12.GetDeltaPosition(bg, result.bias_a)
                                  - (T1.linear() - T2.linear())*Tcb.translation()
                                  - (T1.translation() - T2.translation())*result.scale)/pInt12.dT;
  }
  
  // Lask keyframe velocity
  {
    const Eigen::Isometry3d &T1 = input.back().T2;
    const IMU::Preintegrated &pInt12 = *(input.back().pInt);
    
    velocities.tail<3>() = velocities.segment<3>(3*(N-2)) + result.gravity*pInt12.dT
                           + T1.linear()*Tcb.linear()*pInt12.GetDeltaVelocity(bg, result.bias_a);
  }
  
  result.velocities_ns = timer.ElapsedNanoSeconds();
}

void iterative(const InputType &input, ResultType &result, const double initial_scale = 1.,
               const Eigen::Isometry3d &Tcb = Eigen::Isometry3d::Identity(),
               double *cost = nullptr, double prior = 1e5) {
  std::vector<double*> pointers;
  std::vector<double**> pointers2;

  // Global parameters
  double* bg_ptr = new double[3];
  pointers.push_back(bg_ptr);
  Eigen::Map<Eigen::Vector3d> bg(bg_ptr);
  bg.setZero();

  double* ba_ptr = new double[3];
  pointers.push_back(ba_ptr);
  Eigen::Map<Eigen::Vector3d> ba(ba_ptr);
  ba.setZero();

  double* Rwg_ptr = new double[9];
  pointers.push_back(Rwg_ptr);
  Eigen::Map<Eigen::Matrix3d> Rwg(Rwg_ptr);
  Rwg.setIdentity();

  double* s_ptr = new double[1];
  pointers.push_back(s_ptr);
  s_ptr[0] = initial_scale;

  // Local parameters (for each keyframe)
  double* v0_ptr = new double[3];
  pointers.push_back(v0_ptr);
  Eigen::Map<Eigen::Vector3d> v0(v0_ptr);
  v0.setZero();

  double** parameters = new double*[6];
  pointers2.push_back(parameters);
  parameters[0] = v0_ptr;  // v1
  parameters[1] = nullptr; // v2
  parameters[2] = bg_ptr;  // bg
  parameters[3] = ba_ptr;  // ba
  parameters[4] = Rwg_ptr; // Rwg
  parameters[5] = s_ptr;   // scale

  ceres::Problem problem;

  Eigen::Vector3d dirG;
  dirG.setZero();

  for (unsigned i = 0; i < input.size(); ++i) {
    const Eigen::Isometry3d &T1 = input[i].T1;
    const Eigen::Isometry3d &T2 = input[i].T2;
    const std::shared_ptr<IMU::Preintegrated> pInt = input[i].pInt;

    double* v1_ptr = parameters[0];
    Eigen::Map<Eigen::Vector3d> v1(v1_ptr);

    double* v2_ptr = new double[3];
    pointers.push_back(v2_ptr);
    Eigen::Map<Eigen::Vector3d> v2(v2_ptr);

    v2 = (T2.translation() - T1.translation())/pInt->dT;
    v1 = v2;

    parameters[1] = v2_ptr;

    // Rwg initialization
    dirG -= T1.linear()*pInt->GetUpdatedDeltaVelocity();

    ceres::CostFunction* cost_function = new InertialCostFunction(pInt,
                                                                  T1.linear(), T1.translation(),
                                                                  T2.linear(), T2.translation(),
                                                                  Tcb);
    problem.AddResidualBlock(cost_function, nullptr, parameters, 6);

    double** parameters_ = new double*[6];
    pointers2.push_back(parameters_);
    parameters_[0] = parameters[1];
    parameters_[1] = nullptr;
    parameters_[2] = bg_ptr;
    parameters_[3] = ba_ptr;
    parameters_[4] = Rwg_ptr;
    parameters_[5] = s_ptr;

    parameters = parameters_;
  }
  
  ceres::CostFunction* prior_cost_function = new BiasPriorCostFunction(prior);
  problem.AddResidualBlock(prior_cost_function, nullptr, ba_ptr);

  // Initialize Rwg estimate
  dirG = dirG.normalized();
  const Eigen::Vector3d gI = IMU::GRAVITY_VECTOR.normalized();
  const Eigen::Vector3d v = gI.cross(dirG);
  const double cos_theta = gI.dot(dirG);
  const double theta = std::acos(cos_theta);
  Rwg = ExpSO3(v*theta/v.norm());

  // Add local parameterizations
  GravityParameterization* gravity_local_parameterization = new GravityParameterization;
  problem.SetParameterization(Rwg_ptr, gravity_local_parameterization);

  ScaleParameterization* scale_local_parameterization = new ScaleParameterization;
  problem.SetParameterization(s_ptr, scale_local_parameterization);

  ceres::Solver::Options options;
  options.max_num_iterations = 200;
  ceres::Solver::Summary summary;
  
  Timer timer;
  timer.Start();
  
  ceres::Solve(options, &problem, &summary);

  result.solve_ns = timer.ElapsedNanoSeconds();

  bool converged = (summary.termination_type == ceres::CONVERGENCE);
  if (converged) {
    result.success = true;
    result.scale = s_ptr[0];
    result.bias_g = bg;
    result.bias_a = ba;
    result.gravity = Rwg*IMU::GRAVITY_VECTOR;
    
    if (cost) *cost = summary.final_cost;
  } else {
    result.success = false;
  }

  // Free memory
  for (double* ptr : pointers)
    delete[] ptr;
  for (double** ptr : pointers2)
    delete[] ptr;
}

void mqh_accelerometer(const InputType &input, ResultType &result,
                       const Eigen::Vector3d &bg = Eigen::Vector3d::Zero(),
                       const Eigen::Isometry3d &Tcb = Eigen::Isometry3d::Identity()) {

  // L. Huang, S. Pan, S. Wang, P. Zeng and F. Ye, "A fast initialization method of Visual-Inertial Odometry
  //  based on monocular camera," 2018 Ubiquitous Positioning, Indoor Navigation and Location-Based Services
  //  (UPINLBS), 2018, pp. 1-5, doi: 10.1109/UPINLBS.2018.8559929.
  
  double scale;
  
  Eigen::Vector3d ba;
  ba.setZero();
  
  Eigen::Vector3d gravity;
  
  // Gravity Approximation
  // ---------------------
  
  // R. Mur-Artal and J. D. Tardós, "Visual-Inertial Monocular SLAM With Map Reuse,"
  //  in IEEE Robotics and Automation Letters, vol. 2, no. 2, pp. 796-803, April 2017,
  //  doi: 10.1109/LRA.2017.2653359.
  
  const int N = input.size() + 1; // number of keyframes
  
  Timer timer;
  timer.Start();
  
  // Initialize gravity
  {
    const int n = 4;
    Eigen::VectorXd x_(n);

    // Build linear system
    Eigen::MatrixXd A(3*(N-2), n);
    Eigen::VectorXd b(3*(N-2));
    for (unsigned i = 0; i < input.size()-1; ++i) {
      const Eigen::Isometry3d &T1 = input[i].T1;   // camera to world
      const Eigen::Isometry3d &T2 = input[i+1].T1; // camera to world
      const Eigen::Isometry3d &T3 = input[i+1].T2; // camera to world
      const IMU::Preintegrated &pInt12 = *(input[i].pInt);
      const IMU::Preintegrated &pInt23 = *(input[i+1].pInt);

      Eigen::Matrix3d R1 = T1.linear()*Tcb.linear(); // R1wb
      Eigen::Matrix3d R2 = T2.linear()*Tcb.linear(); // R2wb
      
      Eigen::Matrix3d A_ = R1/pInt12.dT;
      Eigen::Matrix3d B_ = R2/pInt23.dT;

      // lambda:
      // (p1_c - p2_c)/dt12 - (p2_c - p3_c)/dt23

      // beta:
      // - dt12/2 - dt23/2

      // gamma:
      // R1*dV12 - (R1*dP12 + R1_c*pcb - R2_c*pcb)/dt12 + (R2*dP23 + R2_c*pcb - R3_c*pcb)/dt23

      // lambda_i
      A.block<3, 1>(3*i, 0) = (T3.translation() - T2.translation())/pInt23.dT
                               - (T2.translation() - T1.translation())/pInt12.dT;

      // beta_i
      A.block<3, 3>(3*i, 1) = -0.5*(pInt12.dT + pInt23.dT)*Eigen::Matrix3d::Identity();

      // gamma_i
      /*
      b.segment<3>(3*i) = R1*pInt12.GetDeltaVelocity(bg, Eigen::Vector3d::Zero())
                          - (R1*pInt12.GetDeltaPosition(bg, Eigen::Vector3d::Zero()) + T1.linear()*Tcb.translation() - T2.linear()*Tcb.translation())/pInt12.dT
                          + (R2*pInt23.GetDeltaPosition(bg, Eigen::Vector3d::Zero()) + T2.linear()*Tcb.translation() - T3.linear()*Tcb.translation())/pInt23.dT;
      */
      b.segment<3>(3*i) = B_*pInt23.GetDeltaPosition(bg, Eigen::Vector3d::Zero()) - A_*pInt12.GetDeltaPosition(bg, Eigen::Vector3d::Zero())
                           + R1*pInt12.GetDeltaVelocity(bg, Eigen::Vector3d::Zero())
                           + (T2.linear() - T1.linear())*Tcb.translation()/pInt12.dT
                           - (T3.linear() - T2.linear())*Tcb.translation()/pInt23.dT;
    }

    // Solution vector
    x_ = pseudoInverse(A, 1e-6)*b;
    
    // Initial gravity solution
    gravity = x_.tail<3>();
  }

  // Refine Gravity and Accelerometer Bias Estimation
  // ------------------------------------------------

  // R. Mur-Artal and J. D. Tardós, "Visual-Inertial Monocular SLAM With Map Reuse,"
  //  in IEEE Robotics and Automation Letters, vol. 2, no. 2, pp. 796-803, April 2017,
  //  doi: 10.1109/LRA.2017.2653359.

  // Initialize Rwg estimate
  Eigen::Vector3d dirG = gravity.tail<3>();
  const Eigen::Vector3d gI = IMU::GRAVITY_VECTOR.normalized();
  const Eigen::Vector3d v = gI.cross(dirG);
  const double theta = std::atan2(v.norm(), gI.dot(dirG));
  Eigen::Matrix3d Rwi = ExpSO3(v.normalized()*theta);

  // Refinement
  {
    const int n = 6;
    Eigen::VectorXd x_(n);

    Eigen::MatrixXd A(3*(N-2), n);
    Eigen::VectorXd b(3*(N-2));
    for (unsigned i = 0; i < input.size()-1; ++i) {
      const Eigen::Isometry3d &T1 = input[i].T1;   // camera to world
      const Eigen::Isometry3d &T2 = input[i+1].T1; // camera to world
      const Eigen::Isometry3d &T3 = input[i+1].T2; // camera to world
      const IMU::Preintegrated &pInt12 = *(input[i].pInt);
      const IMU::Preintegrated &pInt23 = *(input[i+1].pInt);

      Eigen::Matrix3d R1 = T1.linear()*Tcb.linear(); // R1wb
      Eigen::Matrix3d R2 = T2.linear()*Tcb.linear(); // R2wb
      
      Eigen::Matrix3d A_ = R1/pInt12.dT;
      Eigen::Matrix3d B_ = R2/pInt23.dT;

      // lambda:
      // (p1_c - p2_c)/dt12 - (p2_c - p3_c)/dt23

      // phi:
      // G*Rwi*gI_x*(dt12/2 + dt23/2)

      // xi:
      // (R1*JPa12)/dt12 - R1*JVa12 - (R2*JPa23)/dt23

      // psi:
      // R1*dV12 - (R1*dP12 + R1_c*pcb - R2_c*pcb)/dt12 + (R2*dP23 + R2_c*pcb - R3_c*pcb)/dt23 + G*Rwi*gI*(dt12/2 + dt23/2)

      // lambda_i
      A.block<3, 1>(3*i, 0) = (T3.translation() - T2.translation())/pInt23.dT
                               - (T2.translation() - T1.translation())/pInt12.dT;

      // phi_i
      A.block<3, 2>(3*i, 1) = 0.5*IMU::GRAVITY_MAGNITUDE*(pInt12.dT + pInt23.dT)
                               *(Rwi*Skew(gI)).topLeftCorner<3, 2>();

      // xi_i
      A.block<3, 3>(3*i, 3) = A_*pInt12.JPa - B_*pInt23.JPa - R1*pInt12.JVa;

      // psi_i
      /*
      b.segment<3>(3*i) = R1*pInt12.GetDeltaVelocity(bg, Eigen::Vector3d::Zero())
                          - (R1*pInt12.GetDeltaPosition(bg, Eigen::Vector3d::Zero()) + T1.linear()*Tcb.translation() - T2.linear()*Tcb.translation())/pInt12.dT
                          + (R2*pInt23.GetDeltaPosition(bg, Eigen::Vector3d::Zero()) + T2.linear()*Tcb.translation() - T3.linear()*Tcb.translation())/pInt23.dT
                          + 0.5*IMU::GRAVITY_MAGNITUDE*Rwi*gI*(pInt12.dT + pInt23.dT);
      */
      b.segment<3>(3*i) = B_*pInt23.GetDeltaPosition(bg, Eigen::Vector3d::Zero()) - A_*pInt12.GetDeltaPosition(bg, Eigen::Vector3d::Zero())
                           + R1*pInt12.GetDeltaVelocity(bg, Eigen::Vector3d::Zero())
                           + (T2.linear() - T1.linear())*Tcb.translation()/pInt12.dT
                           - (T3.linear() - T2.linear())*Tcb.translation()/pInt23.dT
                           + 0.5*IMU::GRAVITY_MAGNITUDE*(pInt12.dT + pInt23.dT)*Rwi*gI;
    }

    // Solution vector
    x_ = pseudoInverse(A, 1e-6)*b;
    scale = x_(0);
    ba = x_.tail<3>();
    
    // Refine gravity
    Eigen::Vector3d dTheta;
    dTheta.head<2>() = x_.segment<2>(1);
    dTheta(2) = 0.;
    
    gravity = IMU::GRAVITY_MAGNITUDE*Rwi*(gI - Skew(gI)*dTheta);
  }
  
  // The original paper proposes to use the algorithm below to compute scale and velocity.
  // However, according to some benchmarking:
  //  1. It has quadratic complexity
  //  2. The scale factor is less accurate than Mur-Artal's method
  // So we use the original Mur-Artal's formulation for evaluation (the accuracy of the velocities is not tested)
  
  Eigen::VectorXd velocities(3*N);
  
  timer.Start();
  
  // Recover velocities
  for (unsigned int i = 0; i < input.size(); ++i) {
    const Eigen::Isometry3d &T1 = input[i].T1;
    const Eigen::Isometry3d &T2 = input[i].T2;
    const IMU::Preintegrated &pInt12 = *(input[i].pInt);
    
    // -[(g*dt12^2)/2 + R1*dP12 + R1_c*pcb - R2_c*pcb + p1_c*s - p2_c*s]/dt12
    velocities.segment<3>(3*i) = (-0.5*result.gravity*std::pow(pInt12.dT, 2)
                                  - T1.linear()*Tcb.linear()*pInt12.GetDeltaPosition(bg, result.bias_a)
                                  - (T1.linear() - T2.linear())*Tcb.translation()
                                  - (T1.translation() - T2.translation())*result.scale)/pInt12.dT;
  }
  
  // Lask keyframe velocity
  {
    const Eigen::Isometry3d &T1 = input.back().T2;
    const IMU::Preintegrated &pInt12 = *(input.back().pInt);
    
    velocities.tail<3>() = velocities.segment<3>(3*(N-2)) + result.gravity*pInt12.dT
                           + T1.linear()*Tcb.linear()*pInt12.GetDeltaVelocity(bg, result.bias_a);
  }
  
  result.velocities_ns = timer.ElapsedNanoSeconds();

  // Scale and Velocity Estimation
  // -----------------------------

  // T. Qin and S. Shen, "Robust initialization of monocular visual-inertial estimation on aerial robots,"
  //  2017 IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS), 2017, pp. 4225-4232,
  //  doi: 10.1109/IROS.2017.8206284.

/*
  {
    const int n = 3*(N+1) + 1;
    Eigen::VectorXd x_(n);

    // Since Huang et. al. doesn't say how the linear system should be solved, we use Qin and Shen approach,
    // from VINS-Mono. Source code: https://github.com/HKUST-Aerial-Robotics/VINS-Mono
    // File: vins_estimator/src/initial/initial_aligment.cpp
    Eigen::MatrixXd A(n, n);
    A.setZero();

    Eigen::VectorXd b(n);
    b.setZero();

    for (unsigned i = 0; i < input.size(); ++i) {
      const Eigen::Isometry3d &T1 = input[i].T1;   // camera to world
      const Eigen::Isometry3d &T2 = input[i].T2;  // camera to world
      const IMU::Preintegrated &pInt12 = *(input[i].pInt);

      Eigen::Matrix3d R1 = T1.linear()*Tcb.linear(); // R1wb

      // H:
      // [ -dt12, 0, -dt12^2/2, p2_c - p1_c]
      // [    -1, 1,     -dt12,           0]

      // z:
      //  R1*dP12 + R1_c*pcb - R2_c*pcb
      //                        R1*dV12

      Eigen::MatrixXd tmp_A(6, 10);
      tmp_A.setZero();

      tmp_A.block<3, 3>(0, 0) = -pInt12.dT*Eigen::Matrix3d::Identity();
      tmp_A.block<3, 3>(0, 3).setZero();
      tmp_A.block<3, 3>(0, 6) = -0.5*std::pow(pInt12.dT, 2)*Eigen::Matrix3d::Identity();
      tmp_A.block<3, 1>(0, 9) = (T2.translation() - T1.translation()) / 100.0;

      tmp_A.block<3, 3>(3, 0) = -Eigen::Matrix3d::Identity();
      tmp_A.block<3, 3>(3, 3) = Eigen::Matrix3d::Identity();
      tmp_A.block<3, 3>(3, 6) = -pInt12.dT*Eigen::Matrix3d::Identity();
      tmp_A.block<3, 1>(3, 9).setZero();

      Eigen::VectorXd tmp_b(6);
      tmp_b.block<3, 1>(0, 0) = R1*pInt12.GetDeltaPosition(bg, ba) + (T1.linear() - T2.linear())*Tcb.translation();
      tmp_b.block<3, 1>(3, 0) = R1*pInt12.GetDeltaVelocity(bg, ba);

      Eigen::MatrixXd r_A = tmp_A.transpose()*tmp_A;
      Eigen::VectorXd r_b = tmp_A.transpose()*tmp_b;

      A.block<6, 6>(3*i, 3*i) += r_A.topLeftCorner<6, 6>();
      A.bottomRightCorner<4, 4>() += r_A.bottomRightCorner<4, 4>();
      A.block<6, 4>(3*i, n - 4) += r_A.topRightCorner<6, 4>();
      A.block<4, 6>(n - 4, 3*i) += r_A.bottomLeftCorner<4, 6>();

      b.segment<6>(3*i) += r_b.head<6>();
      b.tail<4>() += r_b.tail<4>();
    }

    // Solution vector
    A = A * 1000.0;
    b = b * 1000.0;
    x_ = A.ldlt().solve(b);
    //scale = x_(n - 1) / 100.0;
    //gravity = x_.segment<3>(n - 4);
  }
*/
  
  result.solve_ns = timer.ElapsedNanoSeconds();
  
  if (scale < 0.) {
    result.success = false;
    return;
  }
  
  result.success = true;
  result.scale = scale;
  result.bias_a = ba;
  result.gravity = gravity;
}

#endif // METHODS_H_
