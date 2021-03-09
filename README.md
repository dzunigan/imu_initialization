# An Analytical Solution to the IMU Initialization Problem for Visual-Inertial Systems
Implementation of "An Analytical Solution to the IMU Initialization Problem for Visual-Inertial Systems"

**Authors:** [David Zuñiga-Noël](http://mapir.isa.uma.es/dzuniga), [Francisco-Angel Moreno](http://mapir.uma.es/mapirwebsite/index.php/people/199-francisco-moreno-due%C3%B1as) and [Javier Gonzalez-Jimenez](http://mapir.uma.es/mapirwebsite/index.php/people/95-javier-gonzalez-jimenez)

Draft: [arXiv:2103.03389](https://arxiv.org/abs/2103.03389)

## Dependencies (tested)
* CMake (3.10.2):
   ```
   sudo apt install cmake
   ```
* Boost (1.65.1):
   ```
   sudo apt install libboost-all-dev
   ```
* Eigen3 (3.3.4):
   ```
   sudo apt install libeigen3-dev
   ```
* Gflags (2.2.1):
   ```
   sudo apt install libgflags-dev
   ```
* Glog (0.3.5)
  ```
  sudo apt install libgoogle-glog-dev
  ```
* Ceres-solver (2.0.0)

  Install [ceres-solver-2.0.0](http://ceres-solver.org/ceres-solver-2.0.0.tar.gz) following [these](http://ceres-solver.org/installation.html) instructions.
  Requires additionally `libatlas-base-dev` and `libsuitesparse-dev`.

## Build
Make sure all dependencies are correctly installed. To build, just run the provided `build.sh` script:
```
git clone https://github.com/dzunigan/imu_initializaiton
bash build.sh
```
which should build the executables in the `./build` directory.

## Source structure
* The analytical solution is implemented in function `proposed_accelerometer()` in `include/methods.h`

* The non-linear constfunctions for iterative optimzation with ceres can be found in `imu_ceres.h`

* The iterative alternative is implemented in function `iterative()` in `include/methods.h`

* The IMU preintegration code (adapted from [ORB_SLAM3](https://github.com/UZ-SLAMLab/ORB_SLAM3)) is implemented in:
   `include/imu_preintegration.h`
   `src/imu_preintegration.cc`
