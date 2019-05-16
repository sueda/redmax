#pragma once
#ifndef _Solver_
#define _Solver_

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <fstream>
#include <memory>

#include "State.h"

struct LinkageSystem;
//struct State;
struct StateDeriv;
struct StateSolve;

struct SolverDataTracker
{
	double time_in_solve;
	int num_iterations;
};

class Solver
{
private:
	void loadMSparse(std::unique_ptr<StateSolve> &SS,
		std::unique_ptr<LinkageSystem> &LS,
		std::unique_ptr<State> &S);

	void computeRHS(std::unique_ptr<StateSolve> &SS,
		const std::unique_ptr<LinkageSystem> &LS,
		const std::unique_ptr<State> &S);

	void pcdSaad2003(Eigen::VectorXd &result,
		const Eigen::VectorXd &LHSx0,
		const Eigen::VectorXd &b,
		std::unique_ptr<StateSolve> &SS,
		std::unique_ptr<LinkageSystem> &LS,
		std::unique_ptr<State> &S,
		std::shared_ptr<State::local_mt> lmt = nullptr,
		double tol = 1e-6, int maxit = 1000);

	void pcdSaad2003_unopt(Eigen::VectorXd &result,
		const Eigen::VectorXd &LHSx0,
		const Eigen::VectorXd &b,
		std::unique_ptr<StateSolve> &SS,
		std::unique_ptr<LinkageSystem> &LS,
		std::unique_ptr<State> &S,
		std::shared_ptr<State::local_mt> lmt = nullptr,
		double tol = 1e-6, int maxit = 1000);
    
#ifdef REDMAX_PARDISO
	Eigen::VectorXd solvePardiso(std::unique_ptr<StateSolve> &SS,
                                 std::unique_ptr<LinkageSystem> &LS,
                                 std::unique_ptr<State> &S);
#endif
    
    Eigen::VectorXd solvePCG_unopt(std::unique_ptr<StateSolve> &SS,
                                   std::unique_ptr<LinkageSystem> &LS,
                                   std::unique_ptr<State> &S);
    
	Eigen::VectorXd solvePCG(std::unique_ptr<StateSolve> &SS,
		std::unique_ptr<LinkageSystem> &LS,
		std::unique_ptr<State> &S);

public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	SolverDataTracker trackLastTimestep;

	Solver();

	Eigen::VectorXd solve(std::unique_ptr<StateSolve> &SS,
		std::unique_ptr<LinkageSystem> &LS, 
		std::unique_ptr<State> &S);

	void printDenseToFile(Eigen::MatrixXd & MG,
		Eigen::VectorXd & f,
		Eigen::VectorXd & res,
		std::string id);

	void printSparseToFile(std::vector<Eigen::Triplet<double>>& MG,
		Eigen::VectorXd & f,
		Eigen::VectorXd & res);

	void addSparseToFile(std::ofstream &outfile, std::string name, Eigen::MatrixXd &M);
	void addSparseToFile(std::ofstream &outfile, std::string name, Eigen::SparseMatrix<double> &M);
	void addVectorToFile(std::ofstream &outfile, std::string name, Eigen::VectorXd &M);
};

#endif
