cmake_minimum_required(VERSION 3.0)

project(imu-initialization CXX)
set(CMAKE_CXX_STANDARD 11)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake_modules")
# message(STATUS "CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")

## Set CXX optimization flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wpedantic -Wall -Wextra")
# -mtune=native -march=native

## Set build type to Release by default
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
endif(NOT CMAKE_BUILD_TYPE)
message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})

if(CMAKE_BUILD_TYPE STREQUAL "Release")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
else(CMAKE_BUILD_TYPE STREQUAL "Release")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g")
endif(CMAKE_BUILD_TYPE STREQUAL "Release")
message(STATUS "Flags: " ${CMAKE_CXX_FLAGS})

## Dependencies options
set(GFLAGS_PREFER_EXPORTED_GFLAGS_CMAKE_CONFIGURATION OFF CACHE BOOL "Use exported cmake config to find gflags?")

## Dependencies
find_package(Boost REQUIRED COMPONENTS filesystem system)
find_package(Ceres REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Gflags REQUIRED)
find_package(Glog REQUIRED)

include_directories(include ${PROJECT_SOURCE_DIR} ${Boost_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR} ${GFLAGS_INCLUDE_DIR} ${GLOG_INCLUDE_DIRS} ${GTEST_INCLUDE_DIRS})

## Header files
file(GLOB_RECURSE HEADER_FILES include/*.hpp include/*.h)
add_custom_target(header_files SOURCES ${HEADER_FILES})
message(STATUS "Header files: " ${HEADER_FILES})

add_executable(experiment01 src/experiment01.cc src/imu_preintegration.cc src/polynomial.cc src/so3.cc src/string.cc src/timer.cc)
target_link_libraries(experiment01 ${Boost_LIBRARIES} ${CERES_LIBRARIES} ${GFLAGS_LIBRARY} ${GLOG_LIBRARIES})

add_executable(experiment01c src/experiment01c.cc src/imu_preintegration.cc src/polynomial.cc src/so3.cc src/string.cc src/timer.cc)
target_link_libraries(experiment01c ${Boost_LIBRARIES} ${CERES_LIBRARIES} ${GFLAGS_LIBRARY} ${GLOG_LIBRARIES})

add_executable(experiment02 src/experiment02.cc src/imu_preintegration.cc src/polynomial.cc src/so3.cc src/string.cc src/timer.cc)
target_link_libraries(experiment02 ${Boost_LIBRARIES} ${CERES_LIBRARIES} ${GFLAGS_LIBRARY} ${GLOG_LIBRARIES})
