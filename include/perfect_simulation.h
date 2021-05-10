
#ifndef SIMULATION_H
#define SIMULATION_H

// STL
#include <cmath>

// EIGEN
#include <Eigen/Core>

#include "util/macros.h"
#include "util/random.h" 

class Gyroscope {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  
  Gyroscope(const double sample_rate, const Eigen::Vector3d &bias)
  {
    RUNTIME_ASSERT(sample_rate > 0.0);
    
    sample_rate_ = sample_rate;
    sample_interval_ = 1./sample_rate;
    bias_ = bias;
  }
  
  Eigen::Vector3d sample(const Eigen::Vector3d &w) const {
    return w + bias_;
  }
  
  double sample_interval() const {
    return sample_interval_;
  }
  
 private:
  double sample_rate_;
  double sample_interval_;
  Eigen::Vector3d bias_;
};

class Accelerometer {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  
  Accelerometer(const double sample_rate, const Eigen::Vector3d &bias, const Eigen::Vector3d &gravity)
  {
    RUNTIME_ASSERT(sample_rate > 0.0);
    
    sample_rate_ = sample_rate;
    sample_interval_ = 1./sample_rate;
    bias_ = bias;
    gravity_ = gravity;
  }
  
  Eigen::Vector3d sample(const Eigen::Vector3d &a, const Eigen::Matrix3d &R) const {
    return R.transpose()*(a - gravity_) + bias_;
  }
  
  inline double sample_interval() const {
    return sample_interval_;
  }
  
 private:
  double sample_rate_;
  double sample_interval_;
  Eigen::Vector3d bias_;
  Eigen::Vector3d gravity_;
};

#endif // SIMULATION_H
