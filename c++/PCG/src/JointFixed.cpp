#include "JointFixed.h"

#include "Block.h"
#include "LinkageSystem.h"
#include "State.h"
#include "StateDeriv.h"
#include "StateSolve.h"

JointFixed::JointFixed(std::unique_ptr<State>& S, std::string pn, std::string cn, Eigen::Vector3d pp, Eigen::Vector3d cp, bool r, int i, Eigen::Matrix3d Q0) :
	Joint(S, JType::Fixed, pn, cn, pp, cp, r, i, 0)
{
	this->S = Eigen::MatrixXd::Zero(6, 0);
	this->ST = Eigen::MatrixXd::Zero(0, 6);
	Sdot = Eigen::MatrixXd::Zero(6, 0);

	q0 = Eigen::VectorXd(constraintNum);
	Q = Eigen::Matrix4d::Identity();
	Q.block<3,3>(0,0) = Q0;

	Pr = Eigen::MatrixXd::Zero(0, 0);
	Psi = Eigen::MatrixXd::Zero(0, 0);
	Kmd = Matrix6d::Zero();
	Dmd = Matrix6d::Zero();
	constraint_index = -1;
}

void JointFixed::updateMaxSparse(std::unique_ptr<StateSolve>& SS, const std::unique_ptr<LinkageSystem>& LS, const std::unique_ptr<State>& S, Matrix6d & parentAdj, Matrix6d & childAdj, Eigen::Vector3d & g) const
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
	for (int r = 0; r < 6; ++r)
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
	row = row + 6;

	SS->row = row;
}

void JointFixed::mapSelf(int c_i, std::unique_ptr<State>& S)
{
	//new (&q) Eigen::Map<Eigen::VectorXd>(&(S->q[c_i]), constraintNum);
	q = Eigen::VectorXd::Zero(constraintNum);

	//new (&qdot) Eigen::Map<Eigen::VectorXd>(&(S->qdot[c_i]), constraintNum);
	qdot = Eigen::VectorXd::Zero(constraintNum);

	updateSelf();
}

void JointFixed::updateSelf()
{
	// Q does not change
}

Vector6d JointFixed::Sqdot(const std::unique_ptr<State>& S)
{
	Vector6d temp = Vector6d::Zero();
	return temp;
}
