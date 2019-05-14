#pragma once
#ifndef __ConstraintJoint__
#define __ConstraintJoint__

#include "Constraint.h"
#include "RigidBodyUtility.h"
#include "State.h"

#include <string>

struct LinkageSystem;
//struct State;
struct StateDeriv;
struct StateSolve;

// ----------------------------------------------------//
// Rotation axis of a block is defined by the location //
// of the joint connecting it to its parent            //
// ----------------------------------------------------//
class ConstraintJoint : public Constraint
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	static void init(std::unique_ptr<StateSolve> &SS, const std::unique_ptr<LinkageSystem> &LS, std::unique_ptr<State> & S);
	static void loadIntoLinkageSystem(const std::unique_ptr<LinkageSystem> &LS);
	static void loadOptimalJointOrdering(const std::unique_ptr<LinkageSystem> &LS, std::unique_ptr<State> & S);

    static void update(std::unique_ptr<StateSolve> &SS, const std::unique_ptr<LinkageSystem> &LS, const std::unique_ptr<State> &S);
    //TODO:'static' member function overrides a virtual function in a base class
	static void draw(const std::unique_ptr<StateSolve> &SS, const std::unique_ptr<LinkageSystem> & LS, const std::unique_ptr<State> & S);

	static int getConstraintNumMaximal(const std::unique_ptr<LinkageSystem> & LS);
	static int getConstraintNumReduced(const std::unique_ptr<LinkageSystem> & LS);
	static int getConstraintNumAdditnl(const std::unique_ptr<LinkageSystem> & LS);
	static int getConstraintNumAdditnlMaximal(const std::unique_ptr<LinkageSystem> & LS);

	static void computeJ_Jdot_x_unopt(Eigen::VectorXd &J_x, Eigen::VectorXd &Jdot_x, const Eigen::VectorXd &x, const std::unique_ptr<StateSolve> &SS, const std::unique_ptr<LinkageSystem> &LS, const std::unique_ptr<State> &S);
	static void computeJ_x_unopt(Eigen::VectorXd &J_x, const Eigen::VectorXd &x, const std::unique_ptr<StateSolve> &SS, const std::unique_ptr<LinkageSystem> &LS, const std::unique_ptr<State> &S);
	static Eigen::VectorXd computeM_x_unopt(const Eigen::VectorXd &x, const std::unique_ptr<StateSolve> &SS, const std::unique_ptr<LinkageSystem> &LS, const std::unique_ptr<State> &S);
	static Eigen::VectorXd computeLHS_x_unopt(const Eigen::VectorXd & x, const std::unique_ptr<StateSolve>& SS, const std::unique_ptr<LinkageSystem>& LS, const std::unique_ptr<State>& S);
	static void computeStiffnessDampingJoint_unopt(Eigen::VectorXd &LHSx, const Eigen::VectorXd & qdot, const Eigen::VectorXd & q, std::unique_ptr<StateSolve>& SS, const std::unique_ptr<LinkageSystem>& LS, const std::unique_ptr<State>& S);
	static Eigen::VectorXd computeJT_x_unopt(const Eigen::VectorXd &x, const std::unique_ptr<StateSolve> &SS, const std::unique_ptr<LinkageSystem> &LS, const std::unique_ptr<State> &S, std::shared_ptr<State::local_mt> lmt = nullptr);
	static void computeJT_x_parallel_unopt(Eigen::VectorXd &Mrtilde, const Eigen::VectorXd &LHSJqd, Eigen::VectorXd &Mr, const Eigen::VectorXd &MJqd, Eigen::VectorXd &fr, const Eigen::VectorXd &fm_MJdqd, const std::unique_ptr<StateSolve> &SS, const std::unique_ptr<LinkageSystem> &LS, const std::unique_ptr<State> &S);

	static void preprocess_preconditioner_unopt(std::unique_ptr<StateSolve> &SS, std::unique_ptr<LinkageSystem> &LS, const std::unique_ptr<State> &S, bool blkdiag);
	static Eigen::VectorXd computeMinv_x_unopt(const Eigen::VectorXd &x, const std::unique_ptr<StateSolve> &SS, const std::unique_ptr<LinkageSystem> &LS, const std::unique_ptr<State> &S, std::shared_ptr<State::local_mt> lmt = nullptr);

	static void preprocess_PCG(std::unique_ptr<StateSolve> &SS, std::unique_ptr<LinkageSystem> &LS, const std::unique_ptr<State> &S);
	static void computeJ_Jdot_x(Eigen::VectorXd &J_x, Eigen::VectorXd &Jdot_x, const Eigen::VectorXd &x, const std::unique_ptr<StateSolve> &SS, const std::unique_ptr<LinkageSystem> &LS, const std::unique_ptr<State> &S);
	static void computeJ_x(Eigen::VectorXd &J_x, const Eigen::VectorXd &x, const std::unique_ptr<StateSolve> &SS, const std::unique_ptr<LinkageSystem> &LS, const std::unique_ptr<State> &S);
	static Eigen::VectorXd computeM_x(const Eigen::VectorXd &x, const std::unique_ptr<StateSolve> &SS, const std::unique_ptr<LinkageSystem> &LS, const std::unique_ptr<State> &S);
	static Eigen::VectorXd computeLHS_x(const Eigen::VectorXd & x, const std::unique_ptr<StateSolve>& SS, const std::unique_ptr<LinkageSystem>& LS, const std::unique_ptr<State>& S);
	static void computeStiffnessDampingJoint(Eigen::VectorXd &LHSx, const Eigen::VectorXd & qdot, const Eigen::VectorXd & q, std::unique_ptr<StateSolve>& SS, const std::unique_ptr<LinkageSystem>& LS, const std::unique_ptr<State>& S);
	static Eigen::VectorXd computeJT_x(const Eigen::VectorXd &x, const std::unique_ptr<StateSolve> &SS, const std::unique_ptr<LinkageSystem> &LS, const std::unique_ptr<State> &S, std::shared_ptr<State::local_mt> lmt = nullptr);

	static void preprocess_PCG_preconditioner(std::unique_ptr<StateSolve> &SS, std::unique_ptr<LinkageSystem> &LS, const std::unique_ptr<State> &S, bool blkdiag);
	static Eigen::VectorXd computeMinv_x(const Eigen::VectorXd &x, const std::unique_ptr<StateSolve> &SS, const std::unique_ptr<LinkageSystem> &LS, const std::unique_ptr<State> &S, std::shared_ptr<State::local_mt> lmt = nullptr);

private:

	ConstraintJoint() {};

	static void updateJacobian(std::unique_ptr<StateSolve> &SS, const std::unique_ptr<LinkageSystem> &LS, const std::unique_ptr<State> &S);
};

#endif
