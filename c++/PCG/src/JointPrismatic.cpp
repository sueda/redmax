#include "JointHinge.h"

#include "Block.h"
#include "LinkageSystem.h"
#include "State.h"
#include "StateDeriv.h"
#include "StateSolve.h"
#include "JointPrismatic.h"

JointPrismatic::JointPrismatic(std::unique_ptr<State>& S, std::string pn, std::string cn, Eigen::Vector3d pp, Eigen::Vector3d cp, bool r, int i, double angle, Eigen::Vector3d & axis, double k, double d) :
	Joint(S, JType::Prismatic, pn, cn, pp, cp, r, i, 1, k, d)
{
	this->axis = axis;

	this->S = Eigen::MatrixXd::Zero(6, 1);
	this->S.block<3, 1>(3, 0) = axis;
	this->ST = this->S.transpose();
	Sdot = Eigen::MatrixXd::Zero(6, 1);

	q0 = Eigen::VectorXd(constraintNum);
	q0[0] = angle;

	Q = Eigen::Matrix4d::Identity();

	Pr = Eigen::MatrixXd::Zero(1, 1);
	Psi = Eigen::MatrixXd::Zero(1, 1);
	Kmd = Matrix6d::Zero();
	Dmd = Matrix6d::Zero();
}

void JointPrismatic::updateMaxSparse(std::unique_ptr<StateSolve>& SS, const std::unique_ptr<LinkageSystem>& LS, const std::unique_ptr<State>& S, Matrix6d & parentAdj, Matrix6d & childAdj, Eigen::Vector3d & g) const
{
}

void JointPrismatic::mapSelf(int c_i, std::unique_ptr<State>& S)
{
	new (&q) Eigen::Map<Eigen::VectorXd>(&(S->q[c_i]), constraintNum);
	q[0] = q0[0];

	new (&qdot) Eigen::Map<Eigen::VectorXd>(&(S->qdot[c_i]), constraintNum);
	qdot = Eigen::VectorXd::Zero(constraintNum);

	updateSelf();
}

void JointPrismatic::updateSelf()
{
	Eigen::Vector3d t_wj0;
	if (q.size() != 0)
		t_wj0 = axis * q[0];
	else
		t_wj0 = axis * q0[0];
	Q.block<3, 1>(0, 3) = t_wj0;
}

Vector6d JointPrismatic::Sqdot(const std::unique_ptr<State>& S)
{
	return this->S * qdot;
}
