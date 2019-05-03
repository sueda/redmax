#pragma once

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <vector>

// ----------------------------------------------------//
// Stores all global state information needed for      //
// computation of the gradient                         //
// ----------------------------------------------------//
struct StateDeriv
{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	// stores pos & vel state history
	Eigen::VectorXd qBuf;
	// built backwards after simulation completes
	Eigen::VectorXd qhatBuf;

	// track the partials of force and obj_func across time
	std::vector<Eigen::MatrixXd> c_dfdx;
	std::vector<Eigen::MatrixXd> c_dfdv;
	std::vector<Eigen::VectorXd> c_dodx;
	std::vector<Eigen::VectorXd> c_dodv;
};