
#define PROGRAM_NAME \
  "experiment01"

#define FLAGS_CASES                                                                                \
  FLAG_CASE(string, logs_dir, "./logs/", "Logs save directory")

#define ARGS_CASES                                                                                 \
  ARG_CASE(dataset_dir)

// STL
#include <algorithm>
#include <cmath>


#include <iterator>
#include <string>
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

namespace fs = boost::filesystem;

using Groundtruth = std::vector<io::state_t>;
using ImuData = io::ImuData;

// IMU parameters
// EuRoC
const double rate = 200.;
const double dt = 1./rate;
const double ng = 1.7e-4;
const double na = 2.e-3;

struct evaluation_t {
  evaluation_t(const std::uint64_t solve_time, const std::uint64_t velocities_time, const std::uint64_t preintegration_time)
    : solve_time(solve_time), velocities_time(velocities_time), preintegration_time(preintegration_time)
  { }

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  std::uint64_t solve_time; // nanoseconds
  std::uint64_t velocities_time; // nanoseconds
  std::uint64_t preintegration_time; // nanoseconds
};

void ValidateArgs() {
  CHECK(fs::is_directory(ARGS_dataset_dir));
}

void ValidateFlags() {
  fs::create_directories(FLAGS_logs_dir);
}

Groundtruth::const_iterator next(Groundtruth::const_iterator i, Groundtruth::const_iterator j, io::timestamp_t dt) {
  if (i == j) {
    LOG(WARNING) << "Already at the end...";
    return i;
  }

  io::timestamp_t t = i->timestamp + dt;
  Groundtruth::const_iterator it = std::lower_bound(i, j, t,
    [](const io::state_t &lhs, const io::timestamp_t rhs) {
      return lhs.timestamp < rhs;
    });
  if (it == i) return i;
  if (it == j) return j;
  Groundtruth::const_iterator it_ = std::next(it, -1);
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

Groundtruth::const_iterator start(const Groundtruth &trajectory, const io::ImuData &imu_data, const unsigned nframes) {
  Groundtruth::const_iterator i = trajectory.cbegin();
  Groundtruth::const_iterator i_ = i;
  while (i != trajectory.cend()) {
    Eigen::Vector3d avgA;
    avgA.setZero();

    io::ImuData::const_iterator it = imu_data.cbegin();
    for (unsigned n = 0; n < nframes; ++n) {
      it = start_imu(it, imu_data.cend(), i->timestamp);
      CHECK(it != imu_data.cend());

      Groundtruth::const_iterator j = next(i, trajectory.cend(), 250000000); // 4 Hz
      CHECK(j != trajectory.cend());

      std::shared_ptr<IMU::Preintegrated> pInt = std::make_shared<IMU::Preintegrated>(Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero());
      while (it != imu_data.cend() && std::llabs(it->timestamp - j->timestamp) > 2500000) {
                                      //it->timestamp < j->timestamp) {
        const Eigen::Vector3d w(it->w_x, it->w_y, it->w_z);
        const Eigen::Vector3d a(it->a_x, it->a_y, it->a_z);
        pInt->IntegrateNewMeasurement(w, a, dt);
        std::advance(it, 1);
      }
      CHECK(it != imu_data.cend());

      avgA += pInt->dV/pInt->dT;
      i = j;
    }

    avgA /= static_cast<double>(nframes);
    const double avgA_error = std::abs(avgA.norm() - IMU::GRAVITY_MAGNITUDE) / IMU::GRAVITY_MAGNITUDE;
    LOG(INFO) << "Average acceleration: " << 100.*avgA_error;
    if (avgA_error > 5e-3) break;

    i = next(i, trajectory.cend(), 500000000); // 0.5s
    i_ = i;
  }
  CHECK(i != trajectory.cend());

  return i_;
}

void save(const std::vector<evaluation_t> &data, const std::string &save_path) {
  const int n = 3;

  Eigen::MatrixXd m(data.size(), n);

  for (unsigned i = 0; i < data.size(); ++i) {
    Eigen::RowVectorXd row(n);
    row << data[i].solve_time, data[i].velocities_time, data[i].preintegration_time;
    m.row(i) = row;
  }

  csv::write(m, save_path);
}

void run(const fs::path &sequence_path) {

  std::string sequence_name = sequence_path.filename().string();
  if (sequence_name == ".")
    sequence_name = sequence_path.parent_path().filename().string();

  LOG(INFO) << "Running experiment: " << sequence_name;

  fs::path trajectory_path = sequence_path / "state_groundtruth_estimate0" / "data.csv";
  CHECK(fs::is_regular_file(trajectory_path)) << "Path not found: " << trajectory_path.string();

  fs::path data_path = sequence_path / "imu0" / "data.csv";
  CHECK(fs::is_regular_file(data_path)) << "Path not found: " << data_path.string();

  Groundtruth trajectory = io::read_file<Groundtruth::value_type>(trajectory_path.string());
  io::ImuData imu_data = io::read_file<io::ImuData::value_type>(data_path.string());

  std::vector<unsigned> possible_nframes = {20, 40, 50, 60, 75, 80, 100, 120, 140, 160, 180, 200, 220, 240};
  for (unsigned nframes : possible_nframes) {
    //Groundtruth::const_iterator i = start(trajectory, imu_data, nframes);
    Groundtruth::const_iterator i = trajectory.cbegin();
    LOG(INFO) << "Starting at " << i->timestamp;
    LOG(INFO) << StringPrintf("With %d frames", nframes);

    std::vector<evaluation_t> proposed_evaluation;
    std::vector<evaluation_t> iterative_evaluation;
    std::vector<evaluation_t> mqh_evaluation;

    //std::uint64_t skipped = 0;
    Groundtruth::const_iterator i_ = i;
    while (i != trajectory.cend()) {
      InputType input;

      Eigen::Vector3d avgBg = Eigen::Vector3d(i->bw_x, i->bw_y, i->bw_z);
      Eigen::Vector3d avgBa = Eigen::Vector3d(i->ba_x, i->ba_y, i->ba_z);

      Eigen::Vector3d avgA;
      avgA.setZero();
      
      std::uint64_t imu_preintegration = 0;

      io::ImuData::const_iterator it = imu_data.cbegin();
      for (unsigned n = 0; n < nframes; ++n) {
        it = start_imu(it, imu_data.cend(), i->timestamp);
        if (it == imu_data.cend()) {
          LOG(WARNING) << "Couldn't find IMU measurement at " << i->timestamp;
          break;
        }

        Groundtruth::const_iterator j = next(i, trajectory.cend(), 250000000); // 4 Hz
        if (j == trajectory.cend()) {
          LOG(WARNING) << "Couldn't find next frame for " << i->timestamp;
          break;
        }

        Timer timer;
        timer.Start();

        std::vector<IMU::Measurement> measurements;
        std::shared_ptr<IMU::Preintegrated> pInt = std::make_shared<IMU::Preintegrated>(Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero());
        while (it != imu_data.cend() && std::llabs(it->timestamp - j->timestamp) > 2500000) {
          const Eigen::Vector3d w(it->w_x, it->w_y, it->w_z);
          const Eigen::Vector3d a(it->a_x, it->a_y, it->a_z);
          pInt->IntegrateNewMeasurement(w, a, dt);
          measurements.emplace_back(w.x(), w.y(), w.z(), a.x(), a.y(), a.z(), dt);
          std::advance(it, 1);
        }
        if (it == imu_data.cend()) {
          LOG(WARNING) << "IMU stream ended!";
          break;
        }

        imu_preintegration += timer.ElapsedNanoSeconds();
        //count++;

        avgBg += Eigen::Vector3d(j->bw_x, j->bw_y, j->bw_z);
        avgBa += Eigen::Vector3d(j->ba_x, j->ba_y, j->ba_z);

        avgA += pInt->dV/pInt->dT;
        input.emplace_back(i->pose, i->timestamp, j->pose, j->timestamp, pInt);
        input.back().vMeasurements = measurements;

        i = j;
      }

      if (input.size() < nframes) {
        LOG(INFO) << StringPrintf("I don't have %d frames. I think dataset ended...", nframes);
        break;
      }

      avgBg /= static_cast<double>(nframes+1);
      avgBa /= static_cast<double>(nframes+1);

      avgA /= static_cast<double>(nframes);
      //const double avgA_error = std::abs(avgA.norm() - IMU::GRAVITY_MAGNITUDE) / IMU::GRAVITY_MAGNITUDE;
      //if (avgA_error > 5e-3) {
        std::uint64_t timestamp = input[0].t1;
        //std::uint64_t initialization_time = i->timestamp - i_->timestamp;

        {
          ResultType gyroscope_result;
          gyroscope_only(input, gyroscope_result);

          ResultType accelerometer_result;
          analytic_accelerometer(input, accelerometer_result, gyroscope_result.bias_g);

          if (gyroscope_result.success & accelerometer_result.success)
            proposed_evaluation.emplace_back(gyroscope_result.solve_ns + accelerometer_result.solve_ns, accelerometer_result.velocities_ns, imu_preintegration);
          else
            LOG(ERROR) << "Proposed method failed at " << timestamp;
        }

        {
          ResultType iterative_result;
          iterative(input, iterative_result, 1., Eigen::Isometry3d::Identity(), nullptr, 0.0);

          if (iterative_result.success)
            iterative_evaluation.emplace_back(iterative_result.solve_ns, std::numeric_limits<double>::quiet_NaN(), imu_preintegration);
          else
            LOG(ERROR) << "Iterative method failed at " << timestamp;
        }

        {
          ResultType gyroscope_result;
          gyroscope_only(input, gyroscope_result, Eigen::Matrix3d::Identity(), false);

          ResultType accelerometer_result;
          mqh_accelerometer(input, accelerometer_result, gyroscope_result.bias_g, Eigen::Isometry3d::Identity());

          if (gyroscope_result.success & accelerometer_result.success)
            mqh_evaluation.emplace_back(gyroscope_result.solve_ns + accelerometer_result.solve_ns, accelerometer_result.velocities_ns, imu_preintegration);
          else
            LOG(ERROR) << "MQH method failed at " << timestamp;
        }

        i = next(i_, trajectory.cend(), 500000000);
        i_ = i;
        //skipped = 0;
      /*
      } else {// next attempt
        skipped += 500000000;
        i = next(i_, trajectory.cend(), skipped); // 0.5s
      }
      */
    }

    std::string proposed_file = StringPrintf("%s_%d_ours.csv", sequence_name.c_str(), nframes);
    LOG(INFO) << "Saving evaluation data into " << proposed_file;
    save(proposed_evaluation, proposed_file);

    std::string iterative_file = StringPrintf("%s_%d_iterative.csv", sequence_name.c_str(), nframes);
    LOG(INFO) << "Saving evaluation data into " << iterative_file;
    save(iterative_evaluation, iterative_file);

    std::string mqh_file = StringPrintf("%s_%d_mqh.csv", sequence_name.c_str(), nframes);
    LOG(INFO) << "Saving evaluation data into " << mqh_file;
    save(mqh_evaluation, mqh_file);
  }

  LOG(INFO) << "done." << std::endl;
}

int main(int argc, char* argv[]) {

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
