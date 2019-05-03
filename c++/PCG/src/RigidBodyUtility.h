#pragma once

#define EIGEN_USE_MKL
#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <set>

#define RB_M_PI 3.1415926535897932
#define THRESHOLD 1e-10

#define ONLINE_MODE

static Eigen::Vector3d emptyVector3d;

// Globlal Simulation Information
enum simType { PCG, PCG_unopt, Pardiso};
extern std::set<simType> maximalCoordList;
extern std::set<simType> reducedCoordList;
extern std::set<simType> reducedNoMatrixList;

// define the simulation method
extern simType simtype;
