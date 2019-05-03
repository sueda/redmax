#include "JointHinge.h"

#include "Block.h"
#include "LinkageSystem.h"
#include "State.h"
#include "StateDeriv.h"
#include "StateSolve.h"

JointHinge::JointHinge(std::unique_ptr<State>& S, std::string pn, std::string cn, Eigen::Vector3d pp, Eigen::Vector3d cp, bool r, int i, double angle, Eigen::Vector3d & axis, double k, double d) :
	Joint(S, JType::Hinge, pn, cn, pp, cp, r, i, 1, k, d)
{
	this->axis = axis;

	this->S = Eigen::MatrixXd::Zero(6, 1);
	this->S.block<3, 1>(0, 0) = axis;
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

void JointHinge::updateMaxSparse(std::unique_ptr<StateSolve>& SS, const std::unique_ptr<LinkageSystem>& LS, const std::unique_ptr<State>& S, Matrix6d & parentAdj, Matrix6d & childAdj, Eigen::Vector3d & g) const
{
	int numObj = (int)LS->blocks.size();
	int row = SS->row;

	int parentIndex;
	int childIndex;
	if (!root)
	{
		parentIndex = LS->find(parentName)->joint->jindex;
	}
	childIndex = LS->find(childName)->joint->jindex;

	// Constrain translational movement
	for (int r = 0; r < 3; ++r)
	{
		for (int c = 0; c < 6; ++c)
		{
			if (!root)
			{
				SS->LHSlist.push_back(Eigen::Triplet<double>(numObj * 6 + row + r, parentIndex * 6 + c, parentAdj(r + 3, c)));
				SS->LHSlist.push_back(Eigen::Triplet<double>(parentIndex * 6 + c, numObj * 6 + row + r, parentAdj(r + 3, c)));
			}
			SS->LHSlist.push_back(Eigen::Triplet<double>(numObj * 6 + row + r, childIndex * 6 + c, -1 * childAdj(r + 3, c)));
			SS->LHSlist.push_back(Eigen::Triplet<double>(childIndex * 6 + c, numObj * 6 + row + r, -1 * childAdj(r + 3, c)));
		}
		// Load the corrective factor -1/h*g into f
		SS->fn[numObj * 6 + row + r] = g[r];
	}
	row = row + 3;

	SS->row = row;
}

void JointHinge::mapSelf(int c_i, std::unique_ptr<State>& S)
{
	new (&q) Eigen::Map<Eigen::VectorXd>(&(S->q[c_i]), constraintNum);
	q[0] = q0[0];

	new (&qdot) Eigen::Map<Eigen::VectorXd>(&(S->qdot[c_i]), constraintNum);
	qdot = Eigen::VectorXd::Zero(constraintNum);

	updateSelf();
}

void JointHinge::updateSelf()
{
	Eigen::Matrix3d R_wj0;
	if (q.size() != 0)
		R_wj0 = Eigen::AngleAxisd(q[0], axis);
	else
		R_wj0 = Eigen::AngleAxisd(q0[0], axis);
	Q.block<3, 3>(0, 0) = R_wj0;
}

Vector6d JointHinge::Sqdot(const std::unique_ptr<State>& S)
{
	return this->S * qdot;
}
