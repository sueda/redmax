#pragma once

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "Rigid.h"

struct StateSolve
{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
		
	enum matrixSolveType { dense, sparse };

	// define the simulation timestep 
	double h;
	// define the damping alpha term
	double alpha;
	// environment gravity
	Eigen::Vector3d grav;
	// baumgarte
	Eigen::Vector3d baumgarte;
	// toggle sparse or dense
	matrixSolveType mst;

	// store the constraint matrix
	std::vector<Eigen::Triplet<double>> LHSlist;
	// stiffness and damping
	std::vector<Eigen::Triplet<double>> Kmlist;	// joint
	std::vector<Eigen::Triplet<double>> Dmlist;	// joint
	std::vector<Eigen::Triplet<double>> bDmlist;	// body
	// fnew vector, holds Mv + hf, RHS of Ax = b for maximal coords
	Eigen::VectorXd fn;
	// only the force acting on each rigid body
	Eigen::VectorXd f;

	// stiffness and damping
	std::vector<Eigen::Triplet<double>> Krlist;
	std::vector<Eigen::Triplet<double>> Drlist;
	// jacobian matrix
	Eigen::MatrixXd J;	// for building 
	std::vector<Eigen::Triplet<double>> Jlist;
	std::vector<Eigen::Triplet<double>> JTlist;
	// and its time derivative
	Eigen::MatrixXd Jd; // for building
	std::vector<Eigen::Triplet<double>> Jdlist;
	// reduced force
	Eigen::VectorXd fr;
	// additional maximal constraints
	std::vector<Eigen::Triplet<double>> Gmlist;
	std::vector<Eigen::Triplet<double>> GmTlist;
	std::vector<Eigen::Triplet<double>> LHSGlist;
	Eigen::VectorXd gm;
	Eigen::MatrixXd Gm; // for CG
	Eigen::MatrixXd GmT; // for CG

	// permanents
	Eigen::SparseMatrix<double> Kr;
	Eigen::SparseMatrix<double> Dr;
	Eigen::SparseMatrix<double> bDm;

	// track the number of constraints
	int row;

	// quick info about the system
	int Mm_dimension;
	int Mr_dimension;
	int Gm_dimension;
	int Gr_dimension;

	// for debugging solves
	//Eigen::MatrixXd LHS;
};