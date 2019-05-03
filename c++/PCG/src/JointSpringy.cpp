#include "JointSpringy.h"

#include "RigidBodyUtility.h"

#include "Block.h"
#include "LinkageSystem.h"
#include "State.h"
#include "StateDeriv.h"
#include "StateSolve.h"

void JointSpringy::updateMaxSparse(std::unique_ptr<StateSolve>& SS, const std::unique_ptr<LinkageSystem>& LS, const std::unique_ptr<State>& S, Matrix6d & parentAdj, Matrix6d & childAdj, Eigen::Vector3d & g) const
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

	///****************************************************************************************************************************
	// TODO: acute rest angles, negative rest angles, DO NOT ALLOW THE JOINT TO WRAP ITSELF (360 degree rotations)
	// TODO: springy joints IF THEY AFFECT THE ROOT JOINT
	Eigen::Matrix4d bToJ0 = E_jp0 * Rigid::inverse(S->E[parentIndex]) * S->E[childIndex];
	double denom = ((bToJ0.block<3, 3>(0, 0) * childPos).norm() * (E_jp0.block<3, 3>(0, 0) * parentPos).norm());
	Eigen::Vector3d sinthetaN = Rigid::bracket3(bToJ0.block<3, 3>(0, 0) * childPos) * (E_jp0.block<3, 3>(0, 0) * parentPos) /
		((bToJ0.block<3, 3>(0, 0) * childPos).norm() * (E_jp0.block<3, 3>(0, 0) * parentPos).norm());
	double theta = asin(sinthetaN[2]) - (RB_M_PI - angle);
	SS->fn.segment<3>(parentIndex * 6) += E_jp0.inverse().block<3, 3>(0, 0) * Eigen::Vector3d::UnitZ() * k * theta;
	SS->fn.segment<3>(childIndex * 6) += -1 * S->E[childIndex].inverse().block<3, 3>(0, 0) * S->E[parentIndex].block<3, 3>(0, 0) * Eigen::Vector3d::UnitZ() * k *theta;

	SS->row = row;
}

void JointSpringy::mapSelf(int c_i, std::unique_ptr<State>& S)
{
}

void JointSpringy::updateSelf()
{
}
