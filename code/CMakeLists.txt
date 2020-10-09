# Source:
# Adapted from TheiaSfM by Kyle Wilson

# ==============================================================================
# Boilerplate
# ==============================================================================

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.3)
IF (COMMAND cmake_policy)
  CMAKE_POLICY(SET CMP0003 NEW)
ENDIF (COMMAND cmake_policy)

set(CMAKE_RELEASE_POSTFIX "")
set(CMAKE_DEBUG_POSTFIX "-debug")

PROJECT(SFMFLEX C CXX)

# ==============================================================================
# Additional cmake find modules and options
# ==============================================================================
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Use a larger inlining threshold for Clang, since it hobbles Eigen,
# resulting in an unreasonably slow version of the blas routines. The
# -Qunused-arguments is needed because CMake passes the inline
# threshold to the linker and clang complains about it and dies.
if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS
    "${CMAKE_CXX_FLAGS} -Qunused-arguments -mllvm -inline-threshold=600")
  # Older versions of Clang (<= 2.9) do not support the 'return-type-c-linkage'
  # option, so check for its presence before adding it to the default flags set.
  include(CheckCXXCompilerFlag)
  check_cxx_compiler_flag("-Wno-return-type-c-linkage"
    HAVE_RETURN_TYPE_C_LINKAGE)
  if (HAVE_RETURN_TYPE_C_LINKAGE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-return-type-c-linkage")
  endif(HAVE_RETURN_TYPE_C_LINKAGE)
endif ()

# ==============================================================================
# External Libraries
# ==============================================================================

# HDF5
SET(MIN_HDF5_VERSION 1.10.5)
FIND_PACKAGE(Hdf5 ${MIN_HDF5_VERSION} REQUIRED)
IF (HDF5_FOUND)
  MESSAGE("-- Found Hdf5 version ${HDF5_VERSION}: ${HDF5_INCLUDE_DIRS}")
ENDIF (HDF5_FOUND)

# Eigen
SET(MIN_EIGEN_VERSION 3.2.2)
FIND_PACKAGE(Eigen ${MIN_EIGEN_VERSION} REQUIRED)
IF (EIGEN_FOUND)
  MESSAGE("-- Found Eigen version ${EIGEN_VERSION}: ${EIGEN_INCLUDE_DIRS}")
ENDIF (EIGEN_FOUND)

# GFlags. The namespace patch is borrow from Ceres Solver (see license in
# FindGflags.cmake)
find_package(Gflags REQUIRED)
if (GFLAGS_FOUND)
  message("-- Found Google Flags: ${GFLAGS_INCLUDE_DIRS} in namespace: ${GFLAGS_NAMESPACE}")
  add_definitions(-DCERES_GFLAGS_NAMESPACE=${GFLAGS_NAMESPACE})
else (GFLAGS_FOUND)
  message(FATAL_ERROR "Can't find Google FLAGS.")
endif (GFLAGS_FOUND)

# Google Logging
MESSAGE("-- Check for Google Log")
FIND_PACKAGE(Glog REQUIRED)
IF (GLOG_FOUND)
  MESSAGE("-- Found Google Log header in: ${GLOG_INCLUDE_DIRS}")
ELSE (GLOG_FOUND)
  MESSAGE(FATAL_ERROR "Can't find Google Log. Please set GLOG_INCLUDE_DIR & "
  "GLOG_LIBRARY")
ENDIF (GLOG_FOUND)

# Ceres
MESSAGE("-- Check for Ceres")
FIND_PACKAGE(Ceres REQUIRED)
IF (CERES_FOUND)
  MESSAGE("-- Found Ceres header in: ${CERES_INCLUDE_DIRS}")
ELSE (CERES_FOUND)
  MESSAGE(FATAL_ERROR "Can't find Ceres. Please set CERES_INCLUDE_DIR & "
  "CERES_LIBRARY")
ENDIF (CERES_FOUND)
MESSAGE("-- Found Ceres library: ${CERES_LIBRARIES}")


INCLUDE_DIRECTORIES(
  src
  "lib/HighFive/include"
  ${CERES_INCLUDE_DIR}
  ${GFLAGS_INCLUDE_DIRS}
  ${GLOG_INCLUDE_DIRS}
  ${EIGEN_INCLUDE_DIRS}
  ${HDF5_INCLUDE_DIRS}
)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
include(OptimizeCompilerFlags)
optimizecompilerflags()

# ==============================================================================
# Targets
# ==============================================================================
add_executable(simple_bundle_adjuster src/simple_bundle_adjuster.cc)
target_link_libraries(simple_bundle_adjuster ceres)

add_executable(compute_jacobian src/compute_jacobian.cc)
target_link_libraries(compute_jacobian ceres ${HDF5_LIBRARIES})

ADD_EXECUTABLE(bundle_adjuster src/bundle_adjuster.cc src/bal_problem.cc)
TARGET_LINK_LIBRARIES(bundle_adjuster ceres ${GFLAGS_LIBRARIES})
