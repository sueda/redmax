#include "JointHinge.h"

#include "Block.h"
#include "LinkageSystem.h"
#include "State.h"
#include "StateDeriv.h"
#include "StateSolve.h"
#include "JointUniversal.h"

JointUniversal::JointUniversal(std::unique_ptr<State>& S, std::string pn, std::string cn, Eigen::Vector3d pp, Eigen::Vector3d cp, bool r, int i, double angleX, double angleY, double k, double d) :
	Joint(S, JType::Universal, pn, cn, pp, cp, r, i, 2, k, d)
{
	this->axis = axis;

	double c1 = std::cos(angleX);
	double c2 = std::cos(angleY);
	double s1 = std::sin(k);
	double s2 = std::sin(d);

	this->S = Eigen::MatrixXd::Zero(6, 2);
	this->S(0,0) = c2;
	this->S(2,0) = s2;
	this->S(1,1) = 1;
	this->ST = this->S.transpose();
	Sdot = Eigen::MatrixXd::Zero(6, 2);
	
	q0[0] = angleX;
	q0[1] = angleY;
	//qdot0[0] = k;
	//qdot0[1] = d;
	//V[0] = 1;
	//V[1] = 1;

	Q = Eigen::Matrix4d::Identity();
	Q(0, 0) = c2;
	Q(0, 1) = s1*c2;
	Q(0, 2) = -c1*s2;
	Q(1, 0) = 0;
	Q(1, 1) = c1;
	Q(1, 2) = s1;
	Q(2, 0) = s2;
	Q(2, 1) = -s1*c2;
	Q(2, 2) = c1*c2;

	Pr = Eigen::MatrixXd::Zero(1, 1);
	Psi = Eigen::MatrixXd::Zero(1, 1);
	Kmd = Matrix6d::Zero();
	Dmd = Matrix6d::Zero();
}

void JointUniversal::updateMaxSparse(std::unique_ptr<StateSolve>& SS, const std::unique_ptr<LinkageSystem>& LS, const std::unique_ptr<State>& S, Matrix6d & parentAdj, Matrix6d & childAdj, Eigen::Vector3d & g) const
{
}

void JointUniversal::mapSelf(int c_i, std::unique_ptr<State>& S)
{
	new (&q) Eigen::Map<Eigen::VectorXd>(&(S->q[c_i]), constraintNum);
	q[0] = q0[0];
	q[1] = q0[1];

	new (&qdot) Eigen::Map<Eigen::VectorXd>(&(S->qdot[c_i]), constraintNum);
	qdot[0] = qdot0[0];
	qdot[1] = qdot0[1];

	updateSelf();
}

void JointUniversal::updateSelf()
{
	double c1, c2, s1, s2;
	if (q.size() != 0)
	{
		c1 = std::cos(q[0]);
		c2 = std::cos(q[1]);
		s1 = std::sin(q[0]);
		s2 = std::sin(q[1]);
	}
	else
	{
		c1 = std::cos(q0[0]);
		c2 = std::cos(q0[1]);
		s1 = std::sin(q0[0]);
		s2 = std::sin(q0[1]);
	}
	double qd1, qd2;
	if (qdot.size() != 0)
	{
		qd1 = qdot[0];
		qd2 = qdot[1];
	}
	else
	{
		qd1 = qdot0[0];
		qd2 = qdot0[1];
	}

	Q(0, 0) = c2;
	Q(1, 0) = s1 * s2;
	Q(2, 0) = -c1 * s2;
	Q(0, 1) = 0;
	Q(1, 1) = c1;
	Q(2, 1) = s1;
	Q(0, 2) = s2;
	Q(1, 2) = -s1 * c2;
	Q(2, 2) = c1 * c2;

	S(0, 0) = c2;
	S(2, 0) = s2;
	Sdot(0, 0) = -s2 * qd1;
	Sdot(2, 0) = c2 * qd2;
	
	ST = S.transpose();
}

Vector6d JointUniversal::Sqdot(const std::unique_ptr<State>& S)
{
	return this->S * qdot;
}
