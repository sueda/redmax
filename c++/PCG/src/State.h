#pragma once

#include "Rigid.h"
#include "Joint.h"

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/PardisoSupport>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/StdVector>

struct Block;

// -----------------------------------------------------//
// Structure of Array style data container for Linkage  //
// System State information								//
// -----------------------------------------------------//
struct State
{
	EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Eigen::Matrix4d)
	EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Matrix6d)

	// time in the simulation
	double &t;
	// Velocity vector - maximal -- phi_ii
	Eigen::VectorXd v;
	// Transformation matrix E -- E_wi
	std::vector<Eigen::Matrix4d> E;
	// Joint positions - reduced
	Eigen::VectorXd q;
	// Velocity vector - reduced
	Eigen::VectorXd qdot;
	// Joint accelerations
	Eigen::VectorXd qddot;

	// Allow for the inertial matricies to have different sizes
	// std::vector<Eigen::Vector3d> size;
	// Mass vector
	std::vector<double> mass;
	// Spatial inertia: M
	std::vector<Matrix6d> M;


	///*********************************************************/
	/////////////////////////////////////////////////////////////
	// 1/29  CHECK ALL ON-THE-FLY ARE NOT SUPERCEDED INCORRECTLY
	/////////////////////////////////////////////////////////////
	///*********************************************************/

	/// Reduced/Maximal PCG SoA
	// SoA constants
	std::vector<int> jindex;
	std::vector<Joint::JType> type;
	std::vector<int> constraintNum;
	std::vector<int> constraint_index;
	std::vector<int> parentIndex;
	std::vector<std::vector<int>> childIndex;
	std::vector<bool> root;
	std::vector<double> damping;
	std::vector<double> d;
	std::vector<double> k;

	// Joints
	std::vector<Vector6d> V;
	std::vector<Eigen::MatrixXd> S;
	std::vector<Eigen::MatrixXd> ST;
	std::vector<Eigen::MatrixXd> Sdot;
	std::vector<Eigen::VectorXd> q0;

	std::vector<Matrix6d> I_j;

	std::vector<Matrix6d> ad_ij;
	std::vector<Matrix6d> ad_ijT;
	std::vector<Matrix6d> ad_jp;
	std::vector<Matrix6d> ad_jpT;

	std::vector<Matrix6d> addot;
	std::vector<Matrix6d> ad_ip;
	std::vector<Matrix6d> ad_ipT;
	std::vector<Matrix6d> add_ip;
	std::vector<Matrix6d> ad_iw;
	std::vector<Matrix6d> ad_wi;

	std::vector<Vector6d> alpha;
	std::vector<Vector6d> beta;

	std::vector<Eigen::MatrixXd> Pr;
	std::vector<Matrix6d> Kmd;
	std::vector<Matrix6d> Dmd;
	std::vector<Matrix6d> D;
	std::vector<Matrix6d> Mhat;
	std::vector<Eigen::MatrixXd> Psi;
	std::vector<Matrix6d> Pi;
	std::vector<Vector6d> Bhat;
	std::vector<Vector6d> Vdot;

	std::vector<Eigen::MatrixXd> adij_S;
	std::vector<Eigen::MatrixXd> ST_adijT;
	std::vector<Eigen::MatrixXd> ST_Mhat;
	std::vector<Eigen::MatrixXd> S_Psi;
	std::vector<Eigen::MatrixXd> Mhat_S_Psi;
	std::vector<Eigen::MatrixXd> ST_Bhat;

	// For non-joint constraints

	// Conjugate Gradient
	struct local_mt
	{
		std::vector<Vector6d> alpha_;
		std::vector<Vector6d> Bhat_;
		std::vector<Vector6d> beta_;
		std::vector<Vector6d> Vdot_;
		std::vector<Eigen::MatrixXd> ST_Bhat_;
	};

	State(double &t) :
		t(t)
	{};
	~State() {};
};