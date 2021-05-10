
#define PROGRAM_NAME \
  "experiment02"

#define FLAGS_CASES                                                                                \
  FLAG_CASE(string, logs_dir, "./logs/", "Logs save directory")                                    \
  FLAG_CASE(uint64, nframes, 10, "Number of frames considered for initialization")

#define ARGS_CASES                                                                                 \
  ARG_CASE(dataset_dir)

// STL
#include <algorithm>
#include <cmath>
#include <iterator>
#include <string>
#include <utility>
#include <vector>

// Boost
#include <boost/filesystem.hpp>
#include <boost/range/iterator_range.hpp>

//Eigen
#include <Eigen/Core>
#include <Eigen/Geometry>

// Glog
#include <glog/logging.h>

#include "imu_preintegration.h"
#include "io.h"
#include "methods.h"

#include "util/args.h"
#include "util/csv.h"
#include "util/timer.h"

namespace fs = boost::filesystem;

using Trajectory = std::vector<io::trajectory_t<double>>;
using Groundtruth = std::vector<io::state_t>;
using ImuData = io::ImuData;

// IMU parameters
// EuRoC
const double rate = 200.;
const double dt = 1./rate;
const double ng = 1.7e-4;
const double na = 2.e-3;

Eigen::Isometry3d Tcb;

struct evaluation_t {
  evaluation_t(const std::uint64_t execution_time,
               const double scale_error, const double gyro_bias_error, const double gyro_bias_error2,
               const double acc_bias_error, double acc_bias_error2, const double gravity_error)
    : execution_time(execution_time),
      scale_error(scale_error), gyro_bias_error(gyro_bias_error), gyro_bias_error2(gyro_bias_error2),
      acc_bias_error(acc_bias_error), acc_bias_error2(acc_bias_error2), gravity_error(gravity_error)
  { }

  std::uint64_t execution_time; // nanoseconds
  double scale_error; // percent
  double gyro_bias_error; // percent
  double gyro_bias_error2; // degress
  double acc_bias_error; // percent
  double acc_bias_error2; // degrees
  double gravity_error; // degrees
};

void ValidateArgs() {
  CHECK(fs::is_directory(ARGS_dataset_dir));
}

void ValidateFlags() {
  fs::create_directories(FLAGS_logs_dir);
  CHECK_GT(FLAGS_nframes, 4);
}

Trajectory read_file_TUM(const std::string &path) {
  std::ifstream input(path);

  Trajectory trajectory;
  for (std::string line; std::getline(input, line);) {
    if (line.empty() || line.front() == '#') continue;

    std::istringstream iss(line);
    double timestamp, tx, ty, tz, qw, qx, qy, qz;
    if (iss >> timestamp >> tx >> ty >> tz >> qx >> qy >> qz >> qw) {
      io::pose_t pose(tx, ty, tz, qw, qx, qy, qz);
      trajectory.emplace_back(timestamp, pose);
    }
  }
  return trajectory;
}

Trajectory::const_iterator next(Trajectory::const_iterator i, Trajectory::const_iterator j, const double dt) {
  if (i == j) {
    LOG(WARNING) << "Already at the end...";
    return i;
  }

  const double t = i->timestamp + dt;
  Trajectory::const_iterator it = std::upper_bound(i, j, t,
    [](const double lhs, const io::trajectory_t<double> &rhs) {
      return lhs < rhs.timestamp;
    });
  if (it == i) return i;
  if (it == j) return j;
  Trajectory::const_iterator it_ = std::next(it, -1);
  if ((it->timestamp - t) > (t - it_->timestamp)) return it_;
  else return it;
}

ImuData::const_iterator start_imu(ImuData::const_iterator i, ImuData::const_iterator j, io::timestamp_t t) {
  ImuData::const_iterator it = std::upper_bound(i, j, t,
    [](const io::timestamp_t lhs, const io::imu_data_t &rhs) {
      return lhs < rhs.timestamp;
    });
  if (it == i) return i;
  if (it == j) return j;
  ImuData::const_iterator it_ = std::next(it, -1);
  if ((it->timestamp - t) > (t - it_->timestamp)) return it_;
  else return it;
}

Trajectory::const_iterator start(const Trajectory &trajectory, const io::ImuData &imu_data) {
  Trajectory::const_iterator i = trajectory.cbegin();
  Trajectory::const_iterator i_ = i;
  while (i != trajectory.cend()) {
    Eigen::Vector3d avgA;
    avgA.setZero();

    io::ImuData::const_iterator it = imu_data.cbegin();
    for (unsigned n = 0; n < FLAGS_nframes; ++n) {
      it = start_imu(it, imu_data.cend(), static_cast<io::timestamp_t>(i->timestamp*1e9));
      CHECK(it != imu_data.cend());

      //Trajectory::const_iterator j = next(i, trajectory.cend(), 0.25); // 4 Hz
      Trajectory::const_iterator j = std::next(i, 1);
      CHECK(j != trajectory.cend());

      std::shared_ptr<IMU::Preintegrated> pInt = std::make_shared<IMU::Preintegrated>(Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero());
      while (it != imu_data.cend() && std::abs(it->timestamp*1e-9 - j->timestamp) > 0.0025) {
        const Eigen::Vector3d w(it->w_x, it->w_y, it->w_z);
        const Eigen::Vector3d a(it->a_x, it->a_y, it->a_z);
        pInt->IntegrateNewMeasurement(w, a, dt);
        std::advance(it, 1);
      }
      CHECK(it != imu_data.cend());

      avgA += pInt->dV/pInt->dT;
      i = j;
    }

    avgA /= static_cast<double>(FLAGS_nframes);
    const double avgA_error = std::abs(avgA.norm() - IMU::GRAVITY_MAGNITUDE) / IMU::GRAVITY_MAGNITUDE;
    LOG(INFO) << "Average acceleration: " << 100.*avgA_error;
    if (avgA_error > 5e-3) break;

    i = next(i_, trajectory.cend(), 0.5);
    i_ = i;
  }
  CHECK(i != trajectory.cend());

  return i_;
}

Groundtruth::const_iterator find_closest(Groundtruth::const_iterator i, Groundtruth::const_iterator j, io::timestamp_t t) {
  Groundtruth::const_iterator it = std::upper_bound(i, j, t,
    [](const io::timestamp_t lhs, const io::state_t &rhs) {
      return lhs < rhs.timestamp;
    });
  if (it == i) return i;
  if (it == j) return j;
  Groundtruth::const_iterator it_ = std::next(it, -1);
  if ((it->timestamp - t) > (t - it_->timestamp)) return it_;
  else return it;
}

Eigen::Isometry3d compute_scale(const InputType &input, const Groundtruth &groundtruth, double &scale_factor) {

  io::Trajectory trajectory;
  trajectory.emplace_back(input.front().t1, input.front().T1);
  for (const input_t &entry : input)
    trajectory.emplace_back(entry.t2, entry.T2);

  std::vector<std::pair<unsigned, unsigned>> pairs;
  for (io::Trajectory::const_iterator it = trajectory.cbegin(); it != trajectory.cend(); ++it) {
    Groundtruth::const_iterator jt = std::upper_bound(groundtruth.cbegin(), groundtruth.cend(), it->timestamp,
    [](const io::timestamp_t lhs, const io::state_t &rhs) {
      return lhs < rhs.timestamp;
    });
    if (jt->timestamp - it->timestamp > 2500000) {
      if (jt == groundtruth.cbegin()) continue;
      std::advance(jt, -1);
      if (jt->timestamp - it->timestamp > 2500000) continue;
    }
    pairs.emplace_back(std::distance(trajectory.cbegin(), it), std::distance(groundtruth.cbegin(), jt));
  }

  const int N = pairs.size();
  CHECK_GE(N, 3) << "At least 3 poses are required!";

  Eigen::MatrixXd src(3, N);
  Eigen::MatrixXd dst(3, N);

  int index = 0;
  for (const std::pair<unsigned, unsigned> &match : pairs) {
    const io::pose_t traj = trajectory.at(match.first).pose;
    const io::pose_t ref = groundtruth.at(match.second).pose;

    src.col(index) = Eigen::Vector3d(traj.tx, traj.ty, traj.tz);
    dst.col(index) = Eigen::Vector3d(ref.tx, ref.ty, ref.tz);
    index++;
  }

  Eigen::Matrix4d M = Eigen::umeyama(src, dst, true);

  scale_factor = std::cbrt(M.block<3, 3>(0, 0).determinant());

  Eigen::Isometry3d T;
  T.linear() = M.block<3, 3>(0, 0)/scale_factor;
  T.translation() = M.block<3, 1>(0, 3);
  return T;
}

void save(const std::vector<evaluation_t> &data, const std::string &save_path) {
  const int n = 7;

  Eigen::MatrixXd m(data.size(), n);

  for (unsigned i = 0; i < data.size(); ++i) {
    Eigen::RowVectorXd row(n);
    row << data[i].execution_time,
           data[i].scale_error, data[i].gyro_bias_error, data[i].gyro_bias_error2,
           data[i].acc_bias_error, data[i].acc_bias_error2, data[i].gravity_error;
    m.row(i) = row;
  }

  csv::write(m, save_path);
}

void run(const fs::path &sequence_path) {
  std::string sequence_name = sequence_path.filename().string();
  if (sequence_name == ".")
    sequence_name = sequence_path.parent_path().filename().string();

  LOG(INFO) << "Running experiment: " << sequence_name;
  LOG(INFO) << StringPrintf("With %d keyframes", FLAGS_nframes);

  fs::path trajectory_path = sequence_path / "KeyFrameTrajectory.txt";
  CHECK(fs::is_regular_file(trajectory_path)) << "Path not found: " << trajectory_path.string();

  fs::path groundtruth_path = sequence_path / "state_groundtruth_estimate0" / "data.csv";
  CHECK(fs::is_regular_file(groundtruth_path)) << "Path not found: " << groundtruth_path.string();

  fs::path data_path = sequence_path / "imu0" / "data.csv";
  CHECK(fs::is_regular_file(data_path)) << "Path not found: " << data_path.string();

  Trajectory trajectory_ = read_file_TUM(trajectory_path.string());
  Groundtruth groundtruth = io::read_file<Groundtruth::value_type>(groundtruth_path.string());
  io::ImuData imu_data = io::read_file<io::ImuData::value_type>(data_path.string());

  // Discard first keyframe
  Trajectory trajectory(std::next(trajectory_.cbegin(), 1), trajectory_.cend());

  Trajectory::const_iterator i = start(trajectory, imu_data);
  LOG(INFO) << "Starting at " << static_cast<io::timestamp_t>(i->timestamp*1e9);

  std::vector<evaluation_t> proposed_evaluation;
  std::vector<evaluation_t> proposed_noprior_evaluation;
  std::vector<evaluation_t> iterative_evaluation;
  std::vector<evaluation_t> iterative_noprior_evaluation;
  std::vector<evaluation_t> mqh_evaluation;

  double skipped = 0.;
  Trajectory::const_iterator i_ = i;
  while (i != trajectory.cend()) {
    Groundtruth::const_iterator gt = find_closest(groundtruth.cbegin(), groundtruth.cend(), static_cast<io::timestamp_t>(i->timestamp*1e9));
    CHECK(gt != groundtruth.cend());
    if (gt == groundtruth.cend()) {
      LOG(WARNING) << "Couldn't find groundtruth for " << static_cast<io::timestamp_t>(i->timestamp*1e9);
      break;
    }

    InputType input;

    Eigen::Vector3d avgBg = Eigen::Vector3d(gt->bw_x, gt->bw_y, gt->bw_z);
    Eigen::Vector3d avgBa = Eigen::Vector3d(gt->ba_x, gt->ba_y, gt->ba_z);

    Eigen::Vector3d avgA;
    avgA.setZero();

    std::uint64_t imu_preintegration = 0;

    io::ImuData::const_iterator it = imu_data.cbegin();
    for (unsigned n = 0; n < FLAGS_nframes; ++n) {
      it = start_imu(it, imu_data.cend(), static_cast<io::timestamp_t>(i->timestamp*1e9));
      if (it == imu_data.cend()) {
        LOG(WARNING) << "Couldn't find IMU measurement at " << static_cast<io::timestamp_t>(i->timestamp*1e9);
        break;
      }
      //LOG(INFO) << "Starting IMU at " << it->timestamp;
      //LOG(INFO) << static_cast<io::timestamp_t>(i->timestamp*1e9);

      Trajectory::const_iterator j = std::next(i, 1);
      if (j == trajectory.cend()) {
        LOG(WARNING) << "Couldn't find next frame for " << static_cast<io::timestamp_t>(i->timestamp*1e9);
        break;
      }
      //LOG(INFO) << "Next frame at " << static_cast<io::timestamp_t>(j->timestamp*1e9);

      gt = find_closest(groundtruth.cbegin(), groundtruth.cend(), static_cast<io::timestamp_t>(j->timestamp*1e9));
      if (gt == groundtruth.cend()) {
        LOG(WARNING) << "Couldn't find groundtruth for " << static_cast<io::timestamp_t>(j->timestamp*1e9);
        break;
      }

      avgBg += Eigen::Vector3d(gt->bw_x, gt->bw_y, gt->bw_z);
      avgBa += Eigen::Vector3d(gt->ba_x, gt->ba_y, gt->ba_z);

      Timer timer;
      timer.Start();

      std::shared_ptr<IMU::Preintegrated> pInt = std::make_shared<IMU::Preintegrated>(Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero());
      while (it != imu_data.cend() && std::abs(it->timestamp*1e-9 - j->timestamp) > 0.0025) {
        const Eigen::Vector3d w(it->w_x, it->w_y, it->w_z);
        const Eigen::Vector3d a(it->a_x, it->a_y, it->a_z);
        pInt->IntegrateNewMeasurement(w, a, dt);
        std::advance(it, 1);
      }
      //LOG(INFO) << "IMU stopped at " << it->timestamp;

      imu_preintegration += timer.ElapsedNanoSeconds();

      if (it == imu_data.cend()) {
        LOG(WARNING) << "IMU stream ended!";
        break;
      }

      avgA += pInt->dV/pInt->dT;
      input.emplace_back(i->pose, static_cast<io::timestamp_t>(i->timestamp*1e9),
                         j->pose, static_cast<io::timestamp_t>(j->timestamp*1e9),
                         pInt);
      i = j;
    }

    if (input.size() < FLAGS_nframes) {
      LOG(INFO) << "I don't have " << FLAGS_nframes << " frames. I think dataset ended...";
      break;
    }

    avgBg /= static_cast<double>(FLAGS_nframes+1);
    avgBa /= static_cast<double>(FLAGS_nframes+1);

    avgA /= static_cast<double>(FLAGS_nframes);
    const double avgA_error = std::abs(avgA.norm() - IMU::GRAVITY_MAGNITUDE) / IMU::GRAVITY_MAGNITUDE;
    //LOG(INFO) << "Average acceleration: " << 100.*avgA_error;
    if (avgA_error > 5e-3) {
      //LOG(INFO) << "Average preintegration time: " << imu_integration / count;
      //imu_integration = 0;
      //count = 0;

      std::uint64_t timestamp = input[0].t1;
      //double initialization_time = i->timestamp - i_->timestamp;

      double true_scale = std::numeric_limits<double>::quiet_NaN();
      Eigen::Isometry3d T = compute_scale(input, groundtruth, true_scale);
      //LOG(INFO) << StringPrintf("True scale: %.6f", true_scale);

      {
        ResultType gyroscope_result;
        gyroscope_only(input, gyroscope_result, Tcb.linear());

        ResultType accelerometer_result;
        analytic_accelerometer(input, accelerometer_result, gyroscope_result.bias_g, Eigen::Vector3d::Zero(), Tcb, 1e5);

        ResultType proposed_result;
        proposed_result.success = gyroscope_result.success & accelerometer_result.success;
        proposed_result.solve_ns = gyroscope_result.solve_ns + accelerometer_result.solve_ns + accelerometer_result.velocities_ns;
        proposed_result.scale = accelerometer_result.scale;
        proposed_result.bias_g = gyroscope_result.bias_g;
        proposed_result.bias_a = accelerometer_result.bias_a;
        proposed_result.gravity = accelerometer_result.gravity;

        if (proposed_result.success) {
          const double scale_error = 100.*std::abs(proposed_result.scale - true_scale)/true_scale;
          const double gyro_bias_error = 100.*std::abs(proposed_result.bias_g.norm() - avgBg.norm()) / avgBg.norm();
          const double gyro_bias_error2 = 180.*std::acos(proposed_result.bias_g.normalized().dot(avgBg.normalized()))/EIGEN_PI;
          const double acc_bias_error = 100.*std::abs(proposed_result.bias_a.norm() - avgBa.norm()) / avgBa.norm();
          const double acc_bias_error2 = 180.*std::acos(proposed_result.bias_a.normalized().dot(avgBa.normalized()))/EIGEN_PI;
          const double gravity_error = 180.*std::acos((T.linear()*proposed_result.gravity).normalized().dot(IMU::GRAVITY_VECTOR.normalized()))/EIGEN_PI;
          proposed_evaluation.emplace_back(imu_preintegration + proposed_result.solve_ns,
              scale_error, gyro_bias_error, gyro_bias_error2, acc_bias_error, acc_bias_error2, gravity_error);
        } else
          LOG(ERROR) << "Proposed method failed at " << timestamp;
      }

      {
        ResultType gyroscope_result;
        gyroscope_only(input, gyroscope_result, Tcb.linear());

        ResultType accelerometer_result;
        analytic_accelerometer(input, accelerometer_result, gyroscope_result.bias_g, Eigen::Vector3d::Zero(), Tcb, 0.0);

        ResultType proposed_result;
        proposed_result.success = gyroscope_result.success & accelerometer_result.success;
        proposed_result.solve_ns = gyroscope_result.solve_ns + accelerometer_result.solve_ns + accelerometer_result.velocities_ns;
        proposed_result.scale = accelerometer_result.scale;
        proposed_result.bias_g = gyroscope_result.bias_g;
        proposed_result.bias_a = accelerometer_result.bias_a;
        proposed_result.gravity = accelerometer_result.gravity;

        if (proposed_result.success) {
          const double scale_error = 100.*std::abs(proposed_result.scale - true_scale)/true_scale;
          const double gyro_bias_error = 100.*std::abs(proposed_result.bias_g.norm() - avgBg.norm()) / avgBg.norm();
          const double gyro_bias_error2 = 180.*std::acos(proposed_result.bias_g.normalized().dot(avgBg.normalized()))/EIGEN_PI;
          const double acc_bias_error = 100.*std::abs(proposed_result.bias_a.norm() - avgBa.norm()) / avgBa.norm();
          const double acc_bias_error2 = 180.*std::acos(proposed_result.bias_a.normalized().dot(avgBa.normalized()))/EIGEN_PI;
          const double gravity_error = 180.*std::acos((T.linear()*proposed_result.gravity).normalized().dot(IMU::GRAVITY_VECTOR.normalized()))/EIGEN_PI;
          proposed_noprior_evaluation.emplace_back(imu_preintegration + proposed_result.solve_ns,
              scale_error, gyro_bias_error, gyro_bias_error2, acc_bias_error, acc_bias_error2, gravity_error);
        } else
          LOG(ERROR) << "Proposed w/o prior method failed at " << timestamp;
      }

      {
        ResultType iterative_result;
        double min_cost = std::numeric_limits<double>::max();

        std::int64_t max_solve_time = 0;
        std::vector<double> scale_values = {1., 4., 16.};
        for (const double scale : scale_values) {
          double cost;
          ResultType result;
          //LOG(INFO) << "Initializing with scale " << scale;
          iterative(input, result, scale, Tcb, &cost);
          max_solve_time = std::max(result.solve_ns, max_solve_time);
          if (cost < min_cost) {
            iterative_result = result;
            min_cost = cost;
          }
        }
        iterative_result.solve_ns = max_solve_time;

        if (iterative_result.success) {
          const double scale_error = 100.*std::abs(iterative_result.scale - true_scale)/true_scale;
          const double gyro_bias_error = 100.*std::abs(iterative_result.bias_g.norm() - avgBg.norm()) / avgBg.norm();
          const double gyro_bias_error2 = 180.*std::acos(iterative_result.bias_g.normalized().dot(avgBg.normalized()))/EIGEN_PI;
          const double acc_bias_error = 100.*std::abs(iterative_result.bias_a.norm() - avgBa.norm()) / avgBa.norm();
          const double acc_bias_error2 = 180.*std::acos(iterative_result.bias_a.normalized().dot(avgBa.normalized()))/EIGEN_PI;
          const double gravity_error = 180.*std::acos((T.linear()*iterative_result.gravity).normalized().dot(IMU::GRAVITY_VECTOR.normalized()))/EIGEN_PI;
          iterative_evaluation.emplace_back(imu_preintegration + iterative_result.solve_ns,
              scale_error, gyro_bias_error, gyro_bias_error2, acc_bias_error, acc_bias_error2, gravity_error);
        } else
          LOG(ERROR) << "Iterative method failed at " << timestamp;
      }

      {
        ResultType iterative_result;
        double min_cost = std::numeric_limits<double>::max();

        std::int64_t max_solve_time = 0;
        std::vector<double> scale_values = {1., 4., 16.};
        for (const double scale : scale_values) {
          double cost;
          ResultType result;
          //LOG(INFO) << "Initializing with scale " << scale;
          iterative(input, result, scale, Tcb, &cost, 0.0);
          max_solve_time = std::max(result.solve_ns, max_solve_time);
          if (cost < min_cost) {
            iterative_result = result;
            min_cost = cost;
          }
        }
        iterative_result.solve_ns = max_solve_time;

        if (iterative_result.success) {
          const double scale_error = 100.*std::abs(iterative_result.scale - true_scale)/true_scale;
          const double gyro_bias_error = 100.*std::abs(iterative_result.bias_g.norm() - avgBg.norm()) / avgBg.norm();
          const double gyro_bias_error2 = 180.*std::acos(iterative_result.bias_g.normalized().dot(avgBg.normalized()))/EIGEN_PI;
          const double acc_bias_error = 100.*std::abs(iterative_result.bias_a.norm() - avgBa.norm()) / avgBa.norm();
          const double acc_bias_error2 = 180.*std::acos(iterative_result.bias_a.normalized().dot(avgBa.normalized()))/EIGEN_PI;
          const double gravity_error = 180.*std::acos((T.linear()*iterative_result.gravity).normalized().dot(IMU::GRAVITY_VECTOR.normalized()))/EIGEN_PI;
          iterative_noprior_evaluation.emplace_back(imu_preintegration + iterative_result.solve_ns,
              scale_error, gyro_bias_error, gyro_bias_error2, acc_bias_error, acc_bias_error2, gravity_error);
        } else
          LOG(ERROR) << "Iterative w/o prior method failed at " << timestamp;
      }

      {
        ResultType gyroscope_result;
        gyroscope_only(input, gyroscope_result, Tcb.linear(), false);

        ResultType accelerometer_result;
        mqh_accelerometer(input, accelerometer_result, gyroscope_result.bias_g, Tcb);

        ResultType mqh_result;
        mqh_result.success = gyroscope_result.success & accelerometer_result.success;
        mqh_result.solve_ns = gyroscope_result.solve_ns + accelerometer_result.solve_ns + accelerometer_result.velocities_ns;
        mqh_result.scale = accelerometer_result.scale;
        mqh_result.bias_g = gyroscope_result.bias_g;
        mqh_result.bias_a = accelerometer_result.bias_a;
        mqh_result.gravity = accelerometer_result.gravity;

        if (mqh_result.success) {
          const double scale_error = 100.*std::abs(mqh_result.scale - 1.);
          const double gyro_bias_error = 100.*std::abs(mqh_result.bias_g.norm() - avgBg.norm()) / avgBg.norm();
          const double gyro_bias_error2 = 180.*std::acos(mqh_result.bias_g.normalized().dot(avgBg.normalized()))/EIGEN_PI;
          const double acc_bias_error = 100.*std::abs(mqh_result.bias_a.norm() - avgBa.norm()) / avgBa.norm();
          const double acc_bias_error2 = 180.*std::acos(mqh_result.bias_a.normalized().dot(avgBa.normalized()))/EIGEN_PI;
          const double gravity_error = 180.*std::acos(mqh_result.gravity.normalized().dot(IMU::GRAVITY_VECTOR.normalized()))/EIGEN_PI;
          mqh_evaluation.emplace_back(imu_preintegration + mqh_result.solve_ns,
              scale_error, gyro_bias_error, gyro_bias_error2, acc_bias_error, acc_bias_error2, gravity_error);
        } else
          LOG(ERROR) << "MQH method failed at " << timestamp;
      }

      i = next(i_, trajectory.cend(), 0.5);
      i_ = i;
      skipped = 0.;
    } else {// next attempt
      skipped += 0.5;
      i = next(i_, trajectory.cend(), skipped); // 0.5s
    }
  }

  //LOG(INFO) << StringPrintf("Average preintegration time: %.3f",
  //                          1e-3 * static_cast<double>(imu_integration) / static_cast<double>(count));

  std::string proposed_file = sequence_name + "_ours.csv";
  LOG(INFO) << "Saving evaluation data into " << proposed_file;
  save(proposed_evaluation, proposed_file);

  std::string proposed_noprior_file = sequence_name + "_ours_noprior.csv";
  LOG(INFO) << "Saving evaluation data into " << proposed_noprior_file;
  save(proposed_noprior_evaluation, proposed_noprior_file);

  std::string iterative_file = sequence_name + "_iterative.csv";
  LOG(INFO) << "Saving evaluation data into " << iterative_file;
  save(iterative_evaluation, iterative_file);

  std::string iterative_noprior_file = sequence_name + "_iterative_noprior.csv";
  LOG(INFO) << "Saving evaluation data into " << iterative_noprior_file;
  save(iterative_noprior_evaluation, iterative_noprior_file);

  std::string mqh_file = sequence_name + "_mqh.csv";
  LOG(INFO) << "Saving evaluation data into " << mqh_file;
  save(mqh_evaluation, mqh_file);

  LOG(INFO) << "done." << std::endl;
}

int main(int argc, char* argv[]) {

  Eigen::Matrix3d Rcb;
  Rcb << 0.0148655429818000,	0.999557249008000,	-0.0257744366974000,
         -0.999880929698000,	0.0149672133247000,	0.00375618835797000,
         0.00414029679422000,	0.0257155299480000,	0.999660727178000;
  Eigen::Vector3d tcb;
  tcb << 0.0652229095355085, -0.0207063854927083, -0.00805460246002837;

  Tcb.linear() = Rcb;
  Tcb.translation() = tcb;

  // Handle help flag
  if (args::HelpRequired(argc, argv)) {
    args::ShowHelp();
    return 0;
  }

  // Parse input flags
  args::ParseCommandLineNonHelpFlags(&argc, &argv, true);

  FLAGS_log_dir = FLAGS_logs_dir;
  FLAGS_stderrthreshold = 0;
  google::InitGoogleLogging(argv[0]);

  // Check number of args
  if (argc-1 != args::NumArgs()) {
    args::ShowHelp();
    return -1;
  }

  // Parse input args
  args::ParseCommandLineArgs(argc, argv);

  // Validate input arguments
  ValidateFlags();
  ValidateArgs();

  IMU::Sigma.block<3, 3>(0, 0) = rate*ng*ng * Eigen::Matrix3d::Identity();
  IMU::Sigma.block<3, 3>(3, 3) = rate*na*na * Eigen::Matrix3d::Identity();

  run(ARGS_dataset_dir);

  return 0;
}
