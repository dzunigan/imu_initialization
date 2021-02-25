
#ifndef SIMULATION_H
#define SIMULATION_H

// STL
#include <cmath>

// EIGEN
#include <Eigen/Core>

#include "util/macros.h"
#include "util/random.h" 

// TODO Bias instability (random walk)
class Gyroscope {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  
  Gyroscope(const double sample_rate, const Eigen::Vector3d &bias, const double noise_density)
  {
    RUNTIME_ASSERT(sample_rate > 0.0);
    RUNTIME_ASSERT(noise_density > 0.0);
    
    sample_rate_ = sample_rate;
    sample_interval_ = 1./sample_rate;
    bias_ = bias;
    noise_density2_ = noise_density * noise_density;
    discrete_sqrt_cov_ = noise_density * std::sqrt(sample_rate) * Eigen::Matrix3d::Identity();
  }
  
  Eigen::Vector3d sample(const Eigen::Vector3d &w) const {
    Eigen::Vector3d n = discrete_sqrt_cov_ * Eigen::Vector3d(RandomGaussian(0.0, 1.0),
                                                             RandomGaussian(0.0, 1.0),
                                                             RandomGaussian(0.0, 1.0)
                                             );
    return w + bias_ + n;
  }
  
  double sample_interval() const {
    return sample_interval_;
  }
  
  Eigen::Matrix3d information_matrix() const {
    return 1./(noise_density2_ * sample_rate_) * Eigen::Matrix3d::Identity();
  }
  
  Eigen::Matrix3d covariance_matrix() const {
    return (noise_density2_ * sample_rate_) * Eigen::Matrix3d::Identity();
  }
  
 private:
  double sample_rate_;
  double sample_interval_;
  double noise_density2_;
  Eigen::Vector3d bias_;
  Eigen::Matrix3d discrete_sqrt_cov_;
};

// TODO Bias instability
class Accelerometer {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  
  Accelerometer(const double sample_rate, const Eigen::Vector3d &bias, const Eigen::Vector3d &gravity, const double noise_density)
  {
    RUNTIME_ASSERT(sample_rate > 0.0);
    RUNTIME_ASSERT(noise_density > 0.0);
    
    sample_rate_ = sample_rate;
    sample_interval_ = 1./sample_rate;
    bias_ = bias;
    gravity_ = gravity;
    noise_density2_ = noise_density * noise_density;
    discrete_sqrt_cov_ = noise_density * std::sqrt(sample_rate) * Eigen::Matrix3d::Identity();
  }
  
  Eigen::Vector3d sample(const Eigen::Vector3d &a, const Eigen::Matrix3d &R) const {
    Eigen::Vector3d n = discrete_sqrt_cov_ * Eigen::Vector3d(RandomGaussian(0.0, 1.0),
                                                             RandomGaussian(0.0, 1.0),
                                                             RandomGaussian(0.0, 1.0)
                                             );
    return R.transpose()*(a - gravity_) + bias_ + n;
  }
  
  inline double sample_interval() const {
    return sample_interval_;
  }
  
  inline Eigen::Matrix3d information_matrix() const {
    return 1./(noise_density2_ * sample_rate_) * Eigen::Matrix3d::Identity();
  }
  
  Eigen::Matrix3d covariance_matrix() const {
    return (noise_density2_ * sample_rate_) * Eigen::Matrix3d::Identity();
  }
  
 private:
  double sample_rate_;
  double sample_interval_;
  double noise_density2_;
  Eigen::Vector3d bias_;
  Eigen::Vector3d gravity_;
  Eigen::Matrix3d discrete_sqrt_cov_;
};

#endif // SIMULATION_H
