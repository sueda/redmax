#include "JointSlider.h"

#include "RigidBodyUtility.h"

#include "Block.h"
#include "LinkageSystem.h"
#include "State.h"
#include "StateDeriv.h"
#include "StateSolve.h"

void JointSlider::updateMaxSparse(std::unique_ptr<StateSolve>& SS, const std::unique_ptr<LinkageSystem>& LS, const std::unique_ptr<State>& S, Matrix6d & parentAdj, Matrix6d & childAdj, Eigen::Vector3d & g) const
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

	// Instead of adding three rows of constraints, we will need to add only two
	for (int r = 0; r < 2; ++r)		// X is the free movement direction
	{
		for (int c = 0; c < 6; ++c)
		{
			if (!root)
			{
				SS->LHSlist.push_back(Eigen::Triplet<double>(numObj * 6 + row + r, parentIndex * 6 + c, parentAdj(r + 4, c)));
				SS->LHSlist.push_back(Eigen::Triplet<double>(parentIndex * 6 + c, numObj * 6 + row + r, parentAdj(r + 4, c)));
			}
			SS->LHSlist.push_back(Eigen::Triplet<double>(numObj * 6 + row + r, childIndex * 6 + c, -1 * childAdj(r + 4, c)));
			SS->LHSlist.push_back(Eigen::Triplet<double>(childIndex * 6 + c, numObj * 6 + row + r, -1 * childAdj(r + 4, c)));
		}
		// Load the corrective factor -1/h*g into f
		SS->fn[numObj * 6 + row + r] = g[r + 1];
	}
	row = row + 2;

	SS->row = row;
}

void JointSlider::mapSelf(int c_i, std::unique_ptr<State>& S)
{
}

void JointSlider::updateSelf()
{
}
