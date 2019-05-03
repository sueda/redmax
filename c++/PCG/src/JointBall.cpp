#include "JointBall.h"

#include "Block.h"
#include "LinkageSystem.h"
#include "State.h"
#include "StateDeriv.h"
#include "StateSolve.h"

void JointBall::updateMaxSparse(std::unique_ptr<StateSolve>& SS, const std::unique_ptr<LinkageSystem>& LS, const std::unique_ptr<State>& S, Matrix6d & parentAdj, Matrix6d & childAdj, Eigen::Vector3d & g) const
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

void JointBall::mapSelf(int c_i, std::unique_ptr<State>& S)
{
}

void JointBall::updateSelf()
{
}
