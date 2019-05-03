#define _USE_MATH_DEFINES
#include <ctgmath>
#include <cmath>

#include "ConstraintJoint.h"

#include "RigidBodyMain.h"
#include "ChronoTimer.h"

#include <queue>
#include <memory>

#ifndef EIGEN_NO_STATIC_ASSERT
#define EIGEN_NO_STATIC_ASSERT
#endif
#include <Eigen/SparseCholesky>	
#include <Eigen/OrderingMethods>
#include <Eigen/Core>
#include <Eigen/Geometry>

#ifdef ONLINE_MODE
#include "online/GLSL.h"
#include <glm/gtc/type_ptr.hpp>
#endif

/// Requires OpenMP 4
//#include <omp.h>
//#pragma omp declare reduction (+: Eigen::VectorXd: omp_out=omp_out+omp_in) initializer(omp_priv=Eigen::VectorXd::Zero(omp_orig.size()))

void ConstraintJoint::init(std::unique_ptr<StateSolve> &SS, const std::unique_ptr<LinkageSystem>& LS, std::unique_ptr<State>& S)
{
	SS->Mr_dimension = getConstraintNumReduced(LS);
	SS->Mm_dimension = getConstraintNumMaximal(LS);
	SS->Gr_dimension = getConstraintNumAdditnl(LS);
	SS->Gm_dimension = getConstraintNumAdditnlMaximal(LS);

	int numObj = (int)LS->joints.size();
	if (reducedCoordList.find(simtype) != reducedCoordList.end())
	{
		int constraint_index = SS->Mr_dimension;
		// loop through joints in order
		for (int i = 0; i < numObj; ++i)
		{
			if (simtype == simType::PCG)
				LS->joints[LS->joint_map[i]]->initState_SoA(constraint_index, LS, S);
			else
				LS->joints[LS->joint_map[i]]->initState(constraint_index, LS, S);
			constraint_index -= LS->joints[LS->joint_map[i]]->constraintNum;
		}

		// init extra constraints for joint solver
		for (int i = 0; i < LS->constraints.size(); ++i)
		{
			LS->constraints[i]->initJoint(SS, LS, S);
		}

		// save permanents
		// joint Stiffness and Damping
		int index;
		std::shared_ptr<Joint> j;
		SS->Drlist.clear();
		SS->Krlist.clear();
		for (int i = 0; i < LS->joints.size(); ++i)
		{
			j = LS->joints[i];
			int ci = j->constraint_index;

			for (int cnum = 0; cnum < j->constraintNum; ++cnum)
			{
				if (j->d > THRESHOLD)
					SS->Drlist.push_back(Eigen::Triplet<double>(ci + cnum, ci + cnum, j->d));
				if (j->k > THRESHOLD)
					SS->Krlist.push_back(Eigen::Triplet<double>(ci + cnum, ci + cnum, j->k));
			}
		}
		SS->Kr = Eigen::SparseMatrix<double>(SS->Mr_dimension, SS->Mr_dimension);
		SS->Kr.setFromTriplets(SS->Krlist.begin(), SS->Krlist.end());
		SS->Dr = Eigen::SparseMatrix<double>(SS->Mr_dimension, SS->Mr_dimension);
		SS->Dr.setFromTriplets(SS->Drlist.begin(), SS->Drlist.end());

		// body Damping
		std::shared_ptr<Block> block;
		SS->bDmlist.clear();
		for (int i = 0; i < LS->blocks.size(); ++i)
		{
			block = LS->blocks[i].second;
			index = block->joint->jindex;
			if (block->damping > THRESHOLD)
			{
				for (int z = 0; z < 6; ++z)
				{
					SS->bDmlist.push_back(Eigen::Triplet<double>(index * 6 + z, index * 6 + z, block->damping));
				}
			}
		}
		SS->bDm = Eigen::SparseMatrix<double>(numObj * 6, numObj * 6);
		SS->bDm.setFromTriplets(SS->bDmlist.begin(), SS->bDmlist.end());

		// load Pr
		int num_constraints = 0;
		for (int i = 0; i < numObj; ++i)
		{
			num_constraints = LS->joints[LS->joint_map[i]]->constraintNum;
			LS->joints[LS->joint_map[i]]->Pr = Eigen::MatrixXd::Identity(num_constraints, num_constraints) *
				(SS->h * LS->joints[LS->joint_map[i]]->d + SS->h * SS->h * LS->joints[LS->joint_map[i]]->k);
		}
	}
}

void ConstraintJoint::loadIntoLinkageSystem(const std::unique_ptr<LinkageSystem> & LS)
{
	Eigen::Matrix4d E_wjp;
	Eigen::Matrix4d E_wji;

	std::shared_ptr<Block> parent;
	std::shared_ptr<Block> child;
	std::shared_ptr<Joint> j;
	for (int jindex = 0; jindex < LS->joints.size(); ++jindex)
	{
		j = LS->joints[jindex];

		E_wjp = Eigen::Matrix4d::Zero();
		E_wji = Eigen::Matrix4d::Zero();
		Eigen::Matrix4d temp;
		Eigen::Matrix4d temp_ji;
		if (!j->root)
		{
			parent = LS->find(j->parentName);
			child = LS->find(j->childName);

			// transform from the parent block to the current joint
			E_wjp.block<3, 3>(0, 0) = parent->E_wi0.block<3, 3>(0, 0);
			E_wjp.block<3, 1>(0, 3) = (parent->E_wi0.block<3, 3>(0, 0)*j->parentPos + parent->E_wi0.block<3, 1>(0, 3));
			E_wjp(3, 3) = 1.0;
			//LS->joints[jindex]->E_jp0 = Rigid::inverse(E_wjp)*parent->E_wi0;
			//LS->joints[jindex]->E_pj0 = Rigid::inverse(LS->joints[jindex]->E_jp0);
			//LS->joints[jindex]->ad_jp0 = Rigid::adjoint(LS->joints[jindex]->E_jp0);

			LS->joints[jindex]->E_pj_ij = Eigen::Matrix4d::Identity();
			LS->joints[jindex]->E_pj_ij.block<3, 1>(0, 3) = j->parentPos - parent->joint->childPos;
			LS->joints[jindex]->ad_pj_ij = Rigid::adjoint(LS->joints[jindex]->E_pj_ij);
			LS->joints[jindex]->ad_ij_pj = Rigid::adjoint(Rigid::inverse(LS->joints[jindex]->E_pj_ij));
		}
		else
		{
			// the root's coordinates are defined in world coordinates
			child = LS->find(j->childName);

			////// transform from the current joint to the parent joint
			LS->joints[jindex]->E_pj_ij = Eigen::Matrix4d::Identity();
			LS->joints[jindex]->E_pj_ij.block<3, 1>(0, 3) = j->parentPos;
			LS->joints[jindex]->ad_pj_ij = Rigid::adjoint(LS->joints[jindex]->E_pj_ij);
			LS->joints[jindex]->ad_ij_pj = Rigid::adjoint(Rigid::inverse(LS->joints[jindex]->E_pj_ij));

			E_wjp.block<3, 3>(0, 0) = Eigen::Matrix3d::Identity();
			E_wjp.block<3, 1>(0, 3) = j->parentPos;
			E_wjp(3, 3) = 1.0;
			//LS->joints[jindex]->E_pj0 = E_wjp;
			//LS->joints[jindex]->E_jp0 = Rigid::inverse(E_wjp);
			//LS->joints[jindex]->ad_jp0 = Rigid::adjoint(LS->joints[jindex]->E_jp0);
		}

		if (j->E_ji0_)
		{
			LS->joints[jindex]->E_ji = j->E_ji0;
		}
		else
		{
			E_wji.block<3, 3>(0, 0) = child->E_wi0.block<3, 3>(0, 0);
			E_wji.block<3, 1>(0, 3) = child->E_wi0.block<3, 3>(0, 0)*j->childPos + child->E_wi0.block<3, 1>(0, 3);
			E_wji(3, 3) = 1.0;
			LS->joints[jindex]->E_ji = Rigid::inverse(E_wji)*child->E_wi0;
		}

		LS->joints[jindex]->E_ij = Rigid::inverse(LS->joints[jindex]->E_ji);
		LS->joints[jindex]->ad_ji = Rigid::adjoint(LS->joints[jindex]->E_ji);
		LS->joints[jindex]->ad_ij = Rigid::adjoint(LS->joints[jindex]->E_ij);
		LS->joints[jindex]->ad_ijT = LS->joints[jindex]->ad_ij.transpose();

		if (j->E_pj0_)
		{
			j->E_pj_ij = j->E_pj0;
			j->ad_pj_ij = Rigid::adjoint(j->E_pj_ij);
			j->ad_ij_pj = Rigid::adjoint(Rigid::inverse(j->E_pj_ij));
		}
	}
}

void ConstraintJoint::loadOptimalJointOrdering(const std::unique_ptr<LinkageSystem>& LS, std::unique_ptr<State> & S)
{
	/// perform DFS and reverse
	int root_joint = - 1;
	int numObj = (int)LS->joints.size();
	int numConstraints = 0;
	for (int j = 0; j < numObj; ++j)
	{
		if (LS->joints[j]->root)
		{
			root_joint = j;
		}
		numConstraints += LS->joints[j]->constraintNum;
	}

	// dont handle any free falling stuff for now
	assert(root_joint != -1);
	assert(LS->joints.size() == LS->blocks.size());
	LS->joint_map.clear();

	// assign optimized index values
	int blockID = -1;
	int jointID = 0;

	/// DFS indecies
	LS->joints[root_joint]->DFSIndexSet(blockID, jointID, LS->blocks, LS->joint_map);

	/// BFS indecies
	//std::queue<std::shared_ptr<Joint>> waitlist;
	//std::queue<std::shared_ptr<Joint>> postlist;
	//waitlist.push(LS->joints[root_joint]);
	//blockID = LS->blocks.size();
	//jointID = numConstraints;
	//while (!waitlist.empty())
	//{
	//	std::shared_ptr<Joint> j = waitlist.front();
	//	waitlist.pop();
	//	std::shared_ptr<Block> block = nullptr;
	//	for (int i = 0; i < LS->blocks.size(); ++i)
	//	{
	//		if (LS->blocks[i].first.compare(j->childName) == 0)
	//			block = LS->blocks[i].second;
	//	}
	//	int numChildren = (int)block->c_joints.size();
	//	for (int i = 0; i < numChildren; ++i)
	//	{
	//		waitlist.push(block->c_joints[i]);
	//	}
	//	LS->joint_map.insert(LS->joint_map.end(), j->index0); // looping through joint_map will access all joints such that every parent is in front of it's children
	//	blockID = blockID - 1;
	//	j->jindex = blockID;
	//	jointID = jointID - j->constraintNum;
	//	j->constraint_index = jointID;
	//}

	/// Mirror MATLAB - bridge scene
	//blockID = 0;
	//jointID = 0;
	//std::shared_ptr<Joint> j;
	//for (int i = (int)LS->joints.size() - 1; i >= 0; --i)
	//{
	//	j = LS->joints[i];
	//	LS->joint_map.insert(LS->joint_map.begin(), j->index0); // looping through joint_map will access all joints such that every parent is in front of it's children
	//	j->jindex = blockID;
	//	blockID = blockID + 1;
	//	if (j->constraintNum == 0)
	//		j->constraint_index = jointID;
	//	else
	//		j->constraint_index = jointID;
	//	jointID = jointID + j->constraintNum;
	//}

	// no matrix CG solve
	int jointCI = 0;
	if (reducedNoMatrixList.find(simtype) != reducedNoMatrixList.end())
	{
		jointID = 0;
		//jointCI = 0;
		for (int i = 0; i < LS->constraints.size(); ++i)
		{
			LS->constraints[i]->constraint_index = jointID;
			jointID = jointID + LS->constraints[i]->constraintNum;
			//LS->constraints[i]->constraint_index_m = jointCI;
			//jointCI = jointCI + LS->constraints[i]->constraintNumMax;
		}
	}
	else  // reduced solve with G matrix
	{
		jointID = 0;
		jointCI = 0;
		for (int i = 0; i < LS->constraints.size(); ++i)
		{
			LS->constraints[i]->constraint_index = jointID;
			jointID = jointID + LS->constraints[i]->constraintNum;
			//LS->constraints[i]->constraint_index_m = jointCI;
			//jointCI = jointCI + LS->constraints[i]->constraintNumMax;
		}
	}
}

int ConstraintJoint::getConstraintNumMaximal(const std::unique_ptr<LinkageSystem> & LS)
{
	// Preprocessing - how big does f need to be?
	int constraints = 0;
	for (int i = 0; i < LS->joints.size(); ++i)
	{
		if (LS->joints[i]->type == Joint::JType::Ball || LS->joints[i]->type == Joint::JType::Springy)
			constraints += 3;
		else if (LS->joints[i]->type == Joint::JType::Powered)
			constraints += 4;
		else if (LS->joints[i]->type == Joint::JType::Slider)
			constraints += 2;
		else if (LS->joints[i]->type == Joint::JType::Hinge)
			constraints += 5;
	}
	for (int i = 0; i < LS->constraints.size(); ++i)
	{
		constraints += LS->constraints[i]->constraintNum;
	}
	return constraints;
}

int ConstraintJoint::getConstraintNumReduced(const std::unique_ptr<LinkageSystem>& LS)
{
	// Preprocessing - how big does f need to be?
	//if (ConstraintJoint::red_constraint_num_save)
	//	return ConstraintJoint::red_constraint_num_save;

	int constraints = 0;
	for (int i = 0; i < LS->joints.size(); ++i)
	{
		constraints += LS->joints[i]->constraintNum;
	}
	return constraints;
}

int ConstraintJoint::getConstraintNumAdditnl(const std::unique_ptr<LinkageSystem>& LS)
{
	int constraints = 0;
	for (int i = 0; i < LS->constraints.size(); ++i)
	{
		constraints += LS->constraints[i]->constraintNum;
	}
	return constraints;
}

int ConstraintJoint::getConstraintNumAdditnlMaximal(const std::unique_ptr<LinkageSystem>& LS)
{
	int constraints = 0;
	for (int i = 0; i < LS->constraints.size(); ++i)
	{
		constraints += LS->constraints[i]->constraintNumMax;
	}
	return constraints;
}

void ConstraintJoint::updateJacobian(std::unique_ptr<StateSolve>& SS, const std::unique_ptr<LinkageSystem>& LS, const std::unique_ptr<State>& S)
{
	int row = SS->row;

	std::shared_ptr<Block> parent;
	std::shared_ptr<Block> child;
	int parentIndex;
	int childIndex;
	std::shared_ptr<Joint> j;

	for (int jmapi = 0; jmapi < LS->joints.size(); ++jmapi)
	{
		int list_index = LS->joint_map[jmapi];
		j = LS->joints[list_index];

		child = j->childBlock;
		childIndex = j->jindex;

		//std::cout << j->childName << std::endl;
		//std::cout << j->E_ji << std::endl << std::endl;
		//std::cout << j->E_ij << std::endl << std::endl;
		//std::cout << j->ad_ij << std::endl << std::endl;
		//std::cout << j->E_pj_ij << std::endl << std::endl;
		//std::cout << j->ad_pj_ij << std::endl << std::endl;
		//std::cout << j->ad_ij_pj << std::endl << std::endl;
		//std::cout << "E_wi \n" << S->E[childIndex] << std::endl;

		for (int cnum = 0; cnum < j->constraintNum; ++cnum)
		{
			// load J(i,i) with ifromJiAdS
			if (j->type != Joint::JType::Fixed)
			{
				int rowindex = childIndex * 6;
				int colindex = j->constraint_index + cnum;

				SS->J.block<6, 1>(rowindex, colindex) = j->ad_ij*j->S.block<6, 1>(0, cnum);

				//std::cout << j->ad_ij << std::endl << std::endl;
				//std::cout << j->S << std::endl << std::endl;

				SS->Jd.block<6, 1>(rowindex, colindex) = j->ad_ij*j->Sdot.block<6, 1>(0, cnum);

				// sparse matrix
				for (int z = 0; z < 6; ++z)
				{
					if (std::abs(SS->J(rowindex + z, colindex)) > THRESHOLD)
					{
						SS->Jlist.push_back(Eigen::Triplet<double>(rowindex + z, colindex, SS->J(rowindex + z, colindex)));
						SS->JTlist.push_back(Eigen::Triplet<double>(colindex, rowindex + z, SS->J(rowindex + z, colindex)));
					}
					if (std::abs(SS->Jd(rowindex + z, colindex)) > THRESHOLD)
					{
						SS->Jdlist.push_back(Eigen::Triplet<double>(rowindex + z, colindex, SS->Jd(rowindex + z, colindex)));
					}
				}
			}
		}

		if (!j->root)
		{
			parent = j->parentBlock;
			parentIndex = parent->joint->jindex;

			if (!j->ad_ip_)
			{
				j->ad_ip = Rigid::adjoint(j->E_iw * S->E[parentIndex]);
				j->ad_ipT = j->ad_ip.transpose();
				j->ad_ip_ = true;
			}

			if (!j->ad_iw_)
			{
				j->ad_iw = Rigid::adjoint(j->E_iw);
				j->ad_iw_ = true;
			}

			if (!j->addot_)
			{
				j->addot = Rigid::addot(S->E[childIndex], S->v.segment<6>(childIndex * 6));
				j->addot_ = true;
			}

			if (!parent->joint->ad_wi_)
			{
				parent->joint->ad_wi = Rigid::adjoint(S->E[parentIndex]);
				parent->joint->ad_wi_ = true;
			}

			if (!parent->joint->addot_)
			{
				parent->joint->addot = Rigid::addot(S->E[parentIndex], S->v.segment<6>(parentIndex * 6));
				parent->joint->addot_ = true;
			}

			j->add_ip = -j->ad_iw * (j->addot *	j->ad_iw * parent->joint->ad_wi - parent->joint->addot);

			// add dependencies into the jacobian until the parent block is rooted (the top of the tree)
	 		std::shared_ptr<Block> ancestor = j->parentBlock;
			std::shared_ptr<Joint> p_joint;
			std::shared_ptr<Joint> c_joint;
			int og_child = childIndex;
			int og_parent = parentIndex;
			while (ancestor != nullptr)
			{
				p_joint = ancestor->joint;
				c_joint = ancestor->joint->childBlock->joint;

				for (int cnum = 0; cnum < p_joint->constraintNum; ++cnum)
				{
					int rowindex = og_child * 6;
					int colindex = p_joint->constraint_index + cnum;

					// J(i,a) = pToiAd*J(p,a)
					SS->J.block<6, 1>(rowindex, colindex) =
						j->ad_ip * SS->J.block<6, 1>(og_parent * 6, colindex);

					// Jd(i,a) = pToiAdd*J(p,a) + pToiAd*Jd(p,a)
					SS->Jd.block<6, 1>(rowindex, colindex) =
						j->add_ip * SS->J.block<6, 1>(og_parent * 6, colindex) +
						j->ad_ip * SS->Jd.block<6, 1>(og_parent * 6, colindex);

					// we actually want a sparse matrix
					for (int z = 0; z < 6; ++z)
					{
						if (std::abs(SS->J(rowindex + z, colindex)) > THRESHOLD)
						{
							SS->Jlist.push_back(Eigen::Triplet<double>(rowindex + z, colindex, SS->J(rowindex + z, colindex)));
							SS->JTlist.push_back(Eigen::Triplet<double>(colindex, rowindex + z, SS->J(rowindex + z, colindex)));
						}
						if (std::abs(SS->Jd(rowindex + z, colindex)) > THRESHOLD)
						{
							SS->Jdlist.push_back(Eigen::Triplet<double>(rowindex + z, colindex, SS->Jd(rowindex + z, colindex)));
						}
					}
				}

				// move up in the family tree
				ancestor = ancestor->joint->parentBlock;
			}
		}
		row = row + j->constraintNum;
	}
	SS->row = row;
}

void ConstraintJoint::computeJ_Jdot_x_unopt(Eigen::VectorXd &J_x, Eigen::VectorXd &Jdot_x, const Eigen::VectorXd &x, const std::unique_ptr<StateSolve>& SS, const std::unique_ptr<LinkageSystem>& LS, const std::unique_ptr<State>& S)
{
	J_x = Eigen::VectorXd::Zero(LS->blocks.size() * 6);
	Jdot_x = Eigen::VectorXd::Zero(LS->blocks.size() * 6);

	std::shared_ptr<Block> parent;
	std::shared_ptr<Block> child;
	int parentIndex;
	int childIndex;
	std::shared_ptr<Joint> j;
	for (int jmapi = 0; jmapi < LS->joints.size(); ++jmapi)
	{
		int list_index = LS->joint_map[jmapi];
		j = LS->joints[list_index];

		child = j->childBlock;
		childIndex = j->jindex;

		// y(i) = ad_ij * S * x(j)
		Vector6d Sx = Vector6d::Zero();
		Vector6d Sdotx = Vector6d::Zero();
		for (int row = 0; row < 6; ++row)
		{
			for (int xi = 0; xi < j->constraintNum; ++xi)
			{
				Sx[row] += j->S(row, xi) * x[j->constraint_index + xi];  // j->S*x[childIndex];
				Sdotx[row] += j->Sdot(row, xi) * x[j->constraint_index + xi];
			}
		}
		J_x.segment<6>(childIndex * 6) = j->ad_ij*Sx;
		Jdot_x.segment<6>(childIndex * 6) = j->ad_ij*Sdotx;
		
		if (!j->root)
		{
			parent = j->parentBlock;
			parentIndex = parent->joint->jindex;

			// compute pToi ad
			if (!j->ad_ip_)
			{
				j->ad_ip = Rigid::adjoint(j->E_iw * S->E[parentIndex]);
				j->ad_ipT = j->ad_ip.transpose();
				j->ad_ip_ = true;
			}

			// compute pToi addot
			if (!j->ad_iw_)
			{
				j->ad_iw = Rigid::adjoint(j->E_iw);
				j->ad_iw_ = true;
			}

			if (!j->addot_)
			{
				j->addot = Rigid::addot(S->E[childIndex], S->v.segment<6>(childIndex * 6));
				j->addot_ = true;
			}

			if (!parent->joint->ad_wi_)
			{
				parent->joint->ad_wi = Rigid::adjoint(S->E[parentIndex]);
				parent->joint->ad_wi_ = true;
			}

			if (!parent->joint->addot_)
			{
				parent->joint->addot = Rigid::addot(S->E[parentIndex], S->v.segment<6>(parentIndex * 6));
				parent->joint->addot_ = true;
			}

			if (!j->add_ip_)
			{
				j->add_ip = -j->ad_iw * (j->addot *	j->ad_iw * parent->joint->ad_wi - parent->joint->addot);
				j->add_ip_ = true;
			}

			//add_ip = -ad_iw*(Rigid::addot(S->E[childIndex], S->v.segment<6>(childIndex * 6)) *
			//	ad_iw * Rigid::adjoint(S->E[parentIndex]) - Rigid::addot(S->E[parentIndex], S->v.segment<6>(parentIndex * 6)));

			// y(i) += ad_ip * y(p)
			J_x.segment<6>(childIndex * 6) += j->ad_ip * J_x.segment<6>(parentIndex * 6);

			// z(i) += ad_ip * z(p) + addot_ip * y(p)
			Jdot_x.segment<6>(childIndex * 6) += j->ad_ip * Jdot_x.segment<6>(parentIndex * 6) +
				j->add_ip * J_x.segment<6>(parentIndex * 6);
		}
	}
}

void ConstraintJoint::computeJ_x_unopt(Eigen::VectorXd &J_x, const Eigen::VectorXd &x, const std::unique_ptr<StateSolve>& SS, const std::unique_ptr<LinkageSystem>& LS, const std::unique_ptr<State>& S)
{
	J_x = Eigen::VectorXd::Zero(LS->blocks.size() * 6);

	std::shared_ptr<Block> parent;
	std::shared_ptr<Block> child;
	int parentIndex;
	int childIndex;
	std::shared_ptr<Joint> j;
	for (int jmapi = 0; jmapi < LS->joints.size(); ++jmapi)
	{
		int list_index = LS->joint_map[jmapi];
		j = LS->joints[list_index];

		child = j->childBlock;
		childIndex = j->jindex;

		// y(i) = ad_ij * S * x(j)
		//Vector6d Sx = Vector6d::Zero();
		//for (int row = 0; row < 6; ++row)
		//{
		//	for (int xi = 0; xi < j->constraintNum; ++xi)
		//	{
		//		Sx[row] += j->S(row, xi) * x[j->constraint_index + xi];  // j->S*x[childIndex];
		//	}
		//}
		//J_x.segment<6>(childIndex * 6) = j->ad_ij*Sx;
		for (int row = 0; row < 6; ++row)
		{
			for (int xi = 0; xi < j->constraintNum; ++xi)
			{
				J_x[childIndex * 6 + row] += j->adij_S(row, xi) * x[j->constraint_index + xi];  // j->S*x[childIndex];
			}
		}

		if (!j->root)
		{
			parent = j->parentBlock;
			parentIndex = parent->joint->jindex;

			// compute pToi ad
			if (!j->ad_ip_)
			{
				j->ad_ip = Rigid::adjoint(j->E_iw * S->E[parentIndex]);
				j->ad_ipT = j->ad_ip.transpose();
				j->ad_ip_ = true;
			}

			// y(i) += ad_ip * y(p)
			J_x.segment<6>(childIndex * 6) += j->ad_ip * J_x.segment<6>(parentIndex * 6);
		}
	}
}

Eigen::VectorXd ConstraintJoint::computeM_x_unopt(const Eigen::VectorXd & x, const std::unique_ptr<StateSolve>& SS, const std::unique_ptr<LinkageSystem>& LS, const std::unique_ptr<State>& S)
{
	Eigen::VectorXd M_x = Eigen::VectorXd::Zero(x.size());

	int jindex;
	int childIndex;
	for (int jmapi = 0; jmapi < LS->joints.size(); ++jmapi)
	{
		jindex = LS->joint_map[jmapi];
		childIndex = LS->joints[jindex]->jindex;

		// M * JX
		for (int mi = 0; mi < 6; ++mi)
		{
			M_x[childIndex * 6 + mi] = S->M[childIndex](mi, mi)*x[childIndex * 6 + mi];
		}
	}
	return M_x;
}

Eigen::VectorXd ConstraintJoint::computeLHS_x_unopt(const Eigen::VectorXd & x, const std::unique_ptr<StateSolve>& SS, const std::unique_ptr<LinkageSystem>& LS, const std::unique_ptr<State>& S)
{
	Eigen::VectorXd M_x = Eigen::VectorXd::Zero(x.size());

	int childIndex;
	std::shared_ptr<Block> b;
	for (int bodyi = 0; bodyi < LS->blocks.size(); ++bodyi)
	{
		b = LS->blocks[bodyi].second;
		childIndex = b->joint->jindex;

		// M * JX
		for (int mi = 0; mi < 6; ++mi)
		{
			M_x[childIndex * 6 + mi] = S->M[childIndex](mi, mi)*x[childIndex * 6 + mi];
		}

		// body damping
		if (b->damping > 0)
		{
			// Mr = (Mr + J'*(h * Dm + h * Bd + h * h * Km)*J + (h * Dr - h * h * Kr)) * qdot
			M_x.segment<6>(childIndex * 6) += (SS->h * b->damping * x.segment<6>(childIndex * 6));
			//SS->fr.segment<6>(childIndex * 6) = SS->fr.segment<6>(childIndex * 6) - (b->damping * S->v.segment<6>(childIndex * 6));
		}
	}
//	Eigen::VectorXd temp;
//#pragma omp parallel for private(temp)// reduction(+ : M_x)
//	for (int i = 0; i < LS->constraints.size(); ++i)
//	{
//		temp = Eigen::VectorXd::Zero(M_x.size());
//		LS->constraints[i]->computeCGProd(temp, x, SS, LS, S);  // most expensive in this function
//#pragma omp critical
//		M_x += temp;
//	}
	for (int i = 0; i < LS->constraints.size(); ++i)
	{
		LS->constraints[i]->computeCGProd(M_x, x, SS, LS, S);  // most expensive in this function
	}
	return M_x;
}

void ConstraintJoint::computeStiffnessDampingJoint_unopt(Eigen::VectorXd &LHSx, const Eigen::VectorXd & qdot, const Eigen::VectorXd & q, std::unique_ptr<StateSolve>& SS, const std::unique_ptr<LinkageSystem>& LS, const std::unique_ptr<State>& S)
{
	int childIndex;
	std::shared_ptr<Joint> j;
	for (int jmapi = 0; jmapi < LS->joints.size(); ++jmapi)
	{
		j = LS->joints[jmapi];
		childIndex = j->constraint_index;

		for (int cnum = 0; cnum < j->constraintNum; ++cnum)
		{
			// Mr = Mr + (h * Dr - h * h * Kr) * qdot
			LHSx[childIndex + cnum] += ((SS->h * j->d + SS->h * SS->h * j->k) * qdot[childIndex + cnum]);
		
			// update on SS->fr
			SS->fr[childIndex + cnum] -= j->k * (q[childIndex + cnum] - j->q0[cnum]);
		}
	}
}

Eigen::VectorXd ConstraintJoint::computeJT_x_unopt(const Eigen::VectorXd & x, const std::unique_ptr<StateSolve>& SS, const std::unique_ptr<LinkageSystem>& LS, const std::unique_ptr<State>& S, std::shared_ptr<State::local_mt> lmt)
{
	Eigen::VectorXd JT_x = Eigen::VectorXd::Zero(SS->Mr_dimension);
	Eigen::VectorXd y = Eigen::VectorXd::Zero(x.size());
	Eigen::VectorXd res;

	std::shared_ptr<Block> parent;
	std::shared_ptr<Block> child;
	int parentIndex;
	int childIndex;
	std::shared_ptr<Joint> j;
	for (int jmapi = (int)LS->joints.size() - 1; jmapi >= 0; --jmapi)
	{
		int list_index = (int)LS->joint_map[jmapi];
		j = LS->joints[list_index];

		child = j->childBlock;
		childIndex = j->jindex;

		// initialize y so MJx is not changed out of scope
		y.segment<6>(childIndex * 6) = x.segment<6>(childIndex * 6);

		// yi += child_alpha
		for (int ci = 0; ci < child->c_joints.size(); ++ci)
		{
			if (lmt == nullptr)
				y.segment<6>(childIndex * 6) += child->c_joints[ci]->alpha;
			else
				y.segment<6>(childIndex * 6) += lmt->alpha_[child->c_joints[ci]->jindex];
		}

		if (!j->root)
		{
			parent = j->parentBlock;
			parentIndex = parent->joint->jindex;

			// set alpha for this joint
			if (!j->ad_ip_)
			{
				j->ad_ip = Rigid::adjoint(j->E_iw * S->E[parentIndex]);
				j->ad_ipT = j->ad_ip.transpose();
				j->ad_ip_ = true;
			}
			if(lmt == nullptr)
				j->alpha = j->ad_ipT * y.segment<6>(childIndex * 6);
			else
				lmt->alpha_[childIndex] = j->ad_ipT * y.segment<6>(childIndex * 6);
		}

		res = j->ST_adijT * y.segment<6>(childIndex * 6);  // most expensive in this fucntion

		assert(res.size() == j->constraintNum);
		for (int resi = 0; resi < res.size(); ++resi)
		{
			JT_x[j->constraint_index + resi] = res[resi];
		}
	}
	return JT_x;
}

// ... broken - errors on the order of 1e-1
void ConstraintJoint::computeJT_x_parallel_unopt(Eigen::VectorXd &Mrtilde, const Eigen::VectorXd &LHSJqd, Eigen::VectorXd &Mr, const Eigen::VectorXd &MJqd, Eigen::VectorXd &fr, const Eigen::VectorXd &fm_MJdqd, const std::unique_ptr<StateSolve>& SS, const std::unique_ptr<LinkageSystem>& LS, const std::unique_ptr<State>& S)
{
	Eigen::VectorXd y1 = Eigen::VectorXd::Zero(LHSJqd.size());
	Eigen::VectorXd y2 = Eigen::VectorXd::Zero(MJqd.size());
	Eigen::VectorXd y3 = Eigen::VectorXd::Zero(fm_MJdqd.size());
	Eigen::VectorXd res1;
	Eigen::VectorXd res2;
	Eigen::VectorXd res3;

	Mrtilde = Eigen::VectorXd::Zero(SS->Mr_dimension);
	Mr = Eigen::VectorXd::Zero(SS->Mr_dimension);
	//fr = Eigen::VectorXd::Zero(SS->Mr_dimension);   // already initialized

	std::shared_ptr<Block> parent;
	std::shared_ptr<Block> child;
	int parentIndex;
	int childIndex;
	std::shared_ptr<Joint> j;
	for (int jmapi = (int)LS->joints.size() - 1; jmapi >= 0; --jmapi)
	{
		int list_index = (int)LS->joint_map[jmapi];
		j = LS->joints[list_index];

		child = j->childBlock;
		childIndex = j->jindex;

		// initialize y so MJx is not changed out of scope
		y1.segment<6>(childIndex * 6) = LHSJqd.segment<6>(childIndex * 6);
		y2.segment<6>(childIndex * 6) = MJqd.segment<6>(childIndex * 6);
		y3.segment<6>(childIndex * 6) = fm_MJdqd.segment<6>(childIndex * 6);

		// yi += child_alpha
		for (int ci = 0; ci < child->c_joints.size(); ++ci)
		{
			y1.segment<6>(childIndex * 6) += child->c_joints[ci]->alpha1;
			y2.segment<6>(childIndex * 6) += child->c_joints[ci]->alpha2;
			y3.segment<6>(childIndex * 6) += child->c_joints[ci]->alpha3;
		}

		if (!j->root)
		{
			parent = j->parentBlock;
			parentIndex = parent->joint->jindex;

			// set alpha for this joint
			if (!j->ad_ip_)
			{
				j->ad_ip = Rigid::adjoint(j->E_iw * S->E[parentIndex]);
				j->ad_ipT = j->ad_ip.transpose();
				j->ad_ip_ = true;
			}
			j->alpha1 = j->ad_ipT * y1.segment<6>(childIndex * 6);
			j->alpha2 = j->ad_ipT * y2.segment<6>(childIndex * 6);
			j->alpha3 = j->ad_ipT * y3.segment<6>(childIndex * 6);
		}

		res1 = j->ST * j->ad_ijT * y1.segment<6>(childIndex * 6);
		res2 = j->ST * j->ad_ijT * y2.segment<6>(childIndex * 6);
		res3 = j->ST * j->ad_ijT * y3.segment<6>(childIndex * 6);

		assert(res1.size() == j->constraintNum);
		for (int resi = 0; resi < res1.size(); ++resi)
		{
			Mrtilde[j->constraint_index + resi] = res1[resi];
			Mr[j->constraint_index + resi] = res2[resi];
			fr[j->constraint_index + resi] += res3[resi];
		}
	}
}

void ConstraintJoint::preprocess_preconditioner_unopt(std::unique_ptr<StateSolve>& SS, std::unique_ptr<LinkageSystem>& LS, const std::unique_ptr<State>& S, bool blkdiag)
{
	std::shared_ptr<Block> child;
	int childIndex;
	std::shared_ptr<Joint> j;
	std::shared_ptr<Joint> j_child;

	// backward traversal
	Matrix6d accum_adj;
	for (int jmapi = (int)LS->joints.size() - 1; jmapi >= 0; --jmapi)
	{
		int list_index = (int)LS->joint_map[jmapi];
		j = LS->joints[list_index];

		child = j->childBlock;
		childIndex = j->jindex;

		// Maximal and Reduced damping terms
		// updated outside preconditioner

		accum_adj = Matrix6d::Zero();
		for (int ci = 0; ci < child->c_joints.size(); ++ci)
		{
			j_child = child->c_joints[ci];
			accum_adj = accum_adj + j_child->ad_jpT * j_child->Pi * j_child->ad_jp;
		}

		// toggle use of damping and stiffness blkdiag in Minv 
		if (blkdiag)
		{
			//std::cout << "Dmd:\n" << (-j->Dmd) << std::endl << std::endl;
			//std::cout << "Dd:\n" << (child->D) << std::endl << std::endl;
			//std::cout << "Kmd:\n" << (j->Kmd) << std::endl << std::endl;
			//std::cout << "total:\n" << (-j->Dmd + child->D - SS->h * j->Kmd) << std::endl << std::endl;
			//std::cout << SS->h * j->ad_ijT * (-j->Dmd + child->D - SS->h * j->Kmd) * j->ad_ij << std::endl << std::endl;
			j->Mhat = j->I_j + SS->h * j->ad_ijT * (-j->Dmd + child->D - SS->h * j->Kmd) * j->ad_ij + accum_adj;
		}
		else 
			j->Mhat = j->I_j + accum_adj;

		// fixed joints do not have index in reduced state
		if (j->type != Joint::JType::Fixed)
		{
			j->Psi = (j->ST * j->Mhat * j->S + j->Pr).inverse();
			j->Pi = j->Mhat - j->Mhat * j->S * j->Psi * j->ST * j->Mhat;
		}

		// precomputations
		j->ST_Mhat = j->ST * j->Mhat;
		j->S_Psi = j->S * j->Psi;
		j->Mhat_S_Psi = j->Mhat * j->S_Psi;
	}
}

Eigen::VectorXd ConstraintJoint::computeMinv_x_unopt(const Eigen::VectorXd & x, const std::unique_ptr<StateSolve>& SS, const std::unique_ptr<LinkageSystem>& LS, const std::unique_ptr<State>& S, std::shared_ptr<State::local_mt> lmt)
{
	std::shared_ptr<Block> parent;
	std::shared_ptr<Block> child;
	int childIndex;
	std::shared_ptr<Joint> j;
	std::shared_ptr<Joint> j_child;
	Eigen::VectorXd Minv_x = Eigen::VectorXd::Zero(x.size());

	// backward traversal
	Eigen::VectorXd x_STB;
	for (int jmapi = (int)LS->joints.size() - 1; jmapi >= 0; --jmapi)
	{
		int list_index = (int)LS->joint_map[jmapi];
		j = LS->joints[list_index];

		child = j->childBlock;
		childIndex = j->jindex;

		// Bhat_j
		Vector6d Bhat = Vector6d::Zero();
		for (int ci = 0; ci < child->c_joints.size(); ++ci)
		{
			// beta changes throughout Minv solve
			if(lmt == nullptr)
				Bhat.noalias() += child->c_joints[ci]->ad_jpT * child->c_joints[ci]->beta;
			else
				Bhat.noalias() += child->c_joints[ci]->ad_jpT * lmt->beta_[child->c_joints[ci]->jindex];
		}
		if (lmt == nullptr)
		{
			j->Bhat = Bhat;
			j->ST_Bhat = j->ST * Bhat;
			x_STB = j->ST_Bhat;
		}
		else
		{
			lmt->Bhat_[childIndex] = Bhat;
			lmt->ST_Bhat_[childIndex] = j->ST * Bhat;
			x_STB = lmt->ST_Bhat_[childIndex];
		}

		// beta_j
		for (int i = 0; i < j->constraintNum; ++i)
		{
			x_STB[i] = x[j->constraint_index + i] - x_STB[i];
		}
		Vector6d beta = Bhat;
		if (j->type != Joint::JType::Fixed)
		{
			beta.noalias() += j->Mhat_S_Psi * x_STB;
		}
		if (lmt == nullptr)
			j->beta = beta;
		else
			lmt->beta_[childIndex] = beta;
	}

	// forward traversal
	Eigen::VectorXd y;	// qddot
	for (int jmapi = 0; jmapi < LS->joints.size(); ++jmapi)
	{
		int list_index = LS->joint_map[jmapi];
		j = LS->joints[list_index];
		childIndex = j->jindex;

		if (lmt == nullptr)
			x_STB = j->ST_Bhat;
		else
			x_STB = lmt->ST_Bhat_[childIndex];
		for (int i = 0; i < j->constraintNum; ++i)
		{
			x_STB[i] = x[j->constraint_index + i] - x_STB[i];
		}

		if(lmt == nullptr)
			j->Vdot = Vector6d::Zero();
		else
			lmt->Vdot_[childIndex] = Vector6d::Zero();
		if (!j->root)
		{
			if (lmt == nullptr)
			{
				j->Vdot = j->ad_jp * j->parentBlock->joint->Vdot;
				x_STB -= j->ST_Mhat * j->Vdot;
			}
			else
			{
				lmt->Vdot_[childIndex] = j->ad_jp * lmt->Vdot_[j->parentBlock->joint->jindex];
				x_STB -= j->ST_Mhat * lmt->Vdot_[childIndex];
			}
		}
		
		y = j->Psi * x_STB;
		if (j->type != Joint::JType::Fixed)
		{
			if (lmt == nullptr)
				j->Vdot += j->S*y;
			else
				lmt->Vdot_[childIndex] += j->S*y;
		}
		else
		{
			if (lmt == nullptr)
				j->Vdot = Vector6d::Zero();
			else
				lmt->Vdot_[childIndex] = Vector6d::Zero();
		}

		// load into result vector
		for (int i = 0; i < y.size(); ++i)
		{
			Minv_x[j->constraint_index + i] = y[i];
		}
	}
	return Minv_x;
}

void ConstraintJoint::preprocess_PCG(std::unique_ptr<StateSolve>& SS, std::unique_ptr<LinkageSystem>& LS, const std::unique_ptr<State>& S)
{
	int list_index;
	int jindex;
	std::shared_ptr<Joint> j;
	// Iterate through the joints in order
	for (int jmapi = 0; jmapi < LS->joints.size(); ++jmapi)
	{
		list_index = LS->joint_map[jmapi];
		jindex = S->jindex[list_index];

		j = LS->joints[list_index];

		S->ad_ip[jindex] = Rigid::adjoint(j->E_iw * S->E[S->parentIndex[jindex]]);
		S->ad_ipT[jindex] = S->ad_ip[jindex].transpose();
		S->ad_iw[jindex] = Rigid::adjoint(j->E_iw);
		S->ad_wi[jindex] = Rigid::adjoint(S->E[jindex]);
		S->addot[jindex] = Rigid::addot(S->E[jindex], S->v.segment<6>(jindex * 6));
	}

	// Iterate through the joints in order
	for (int jmapi = 0; jmapi < LS->joints.size(); ++jmapi)
	{
		list_index = LS->joint_map[jmapi];
		jindex = S->jindex[list_index];

		if (!S->root[jindex])
		{
			S->add_ip[jindex] = -S->ad_iw[jindex] * (S->addot[jindex] * S->ad_iw[jindex] * S->ad_wi[S->parentIndex[jindex]] - S->addot[S->parentIndex[jindex]]);
		}
	}
}

void ConstraintJoint::computeJ_Jdot_x(Eigen::VectorXd & J_x, Eigen::VectorXd & Jdot_x, const Eigen::VectorXd & x, const std::unique_ptr<StateSolve>& SS, const std::unique_ptr<LinkageSystem>& LS, const std::unique_ptr<State>& S)
{
	J_x = Eigen::VectorXd::Zero(LS->blocks.size() * 6);
	Jdot_x = Eigen::VectorXd::Zero(LS->blocks.size() * 6);

	int parentIndex;
	int jindex;
	for (int jmapi = 0; jmapi < LS->joints.size(); ++jmapi)
	{
		int list_index = LS->joint_map[jmapi];
		jindex = S->jindex[list_index];

		// y(i) = ad_ij * S * x(j)
		Vector6d Sx = Vector6d::Zero();
		Vector6d Sdotx = Vector6d::Zero();
		for (int row = 0; row < 6; ++row)
		{
			for (int xi = 0; xi < S->constraintNum[jindex]; ++xi)
			{
				Sx[row] += S->S[jindex](row, xi) * x[S->constraint_index[jindex] + xi];
				Sdotx[row] += S->Sdot[jindex](row, xi) * x[S->constraint_index[jindex] + xi];
			}
		}
		J_x.segment<6>(jindex * 6) = S->ad_ij[jindex]*Sx;
		Jdot_x.segment<6>(jindex * 6) = S->ad_ij[jindex]*Sdotx;

		if (!S->root[jindex])
		{
			parentIndex = S->parentIndex[jindex];

			// y(i) += ad_ip * y(p)
			J_x.segment<6>(jindex * 6) += S->ad_ip[jindex] * J_x.segment<6>(parentIndex * 6);

			// z(i) += ad_ip * z(p) + addot_ip * y(p)
			Jdot_x.segment<6>(jindex * 6) += S->ad_ip[jindex] * Jdot_x.segment<6>(parentIndex * 6) +
				S->add_ip[jindex] * J_x.segment<6>(parentIndex * 6);
		}
	}
}

void ConstraintJoint::computeJ_x(Eigen::VectorXd & J_x, const Eigen::VectorXd & x, const std::unique_ptr<StateSolve>& SS, const std::unique_ptr<LinkageSystem>& LS, const std::unique_ptr<State>& S)
{
	J_x = Eigen::VectorXd::Zero(LS->blocks.size() * 6);

	int parentIndex;
	int jindex;
	for (int jmapi = 0; jmapi < LS->joints.size(); ++jmapi)
	{
		int list_index = LS->joint_map[jmapi];
		jindex = S->jindex[list_index];

		for (int row = 0; row < 6; ++row)
		{
			for (int xi = 0; xi < S->constraintNum[jindex]; ++xi)
			{
				J_x[jindex * 6 + row] += S->adij_S[jindex](row, xi) * x[S->constraint_index[jindex] + xi];
			}
		}

		if (!S->root[jindex])
		{
			parentIndex = S->parentIndex[jindex];

			// y(i) += ad_ip * y(p)
			J_x.segment<6>(jindex * 6) += S->ad_ip[jindex] * J_x.segment<6>(parentIndex * 6);
		}
	}
}

Eigen::VectorXd ConstraintJoint::computeM_x(const Eigen::VectorXd & x, const std::unique_ptr<StateSolve>& SS, const std::unique_ptr<LinkageSystem>& LS, const std::unique_ptr<State>& S)
{
	Eigen::VectorXd M_x = Eigen::VectorXd::Zero(x.size());

	int jindex;
	for (int jmapi = 0; jmapi < LS->joints.size(); ++jmapi)
	{
		jindex = S->jindex[LS->joint_map[jmapi]];

		// M * JX
		for (int mi = 0; mi < 6; ++mi)
		{
			M_x[jindex * 6 + mi] = S->M[jindex](mi, mi)*x[jindex * 6 + mi];
		}
	}
	return M_x;
}

Eigen::VectorXd ConstraintJoint::computeLHS_x(const Eigen::VectorXd & x, const std::unique_ptr<StateSolve>& SS, const std::unique_ptr<LinkageSystem>& LS, const std::unique_ptr<State>& S)
{
	Eigen::VectorXd M_x = Eigen::VectorXd::Zero(x.size());

	int jindex;
	for (int jmapi = 0; jmapi < LS->joints.size(); ++jmapi)
	{
		jindex = S->jindex[LS->joint_map[jmapi]];

		// M * JX
		for (int mi = 0; mi < 6; ++mi)
		{
			M_x[jindex * 6 + mi] = S->M[jindex](mi, mi)*x[jindex * 6 + mi];
		}

		// body damping
		if (S->damping[jindex] > 0)
		{
			// Mr = (Mr + J'*(h * Dm + h * Bd + h * h * Km)*J + (h * Dr - h * h * Kr)) * qdot
			M_x.segment<6>(jindex * 6) += (SS->h * S->damping[jindex] * x.segment<6>(jindex * 6));
			//SS->fr.segment<6>(childIndex * 6) = SS->fr.segment<6>(childIndex * 6) - (b->damping * S->v.segment<6>(childIndex * 6));
		}
	}
	for (int i = 0; i < LS->constraints.size(); ++i)
	{
		/// SoA
		LS->constraints[i]->computeCGProd(M_x, x, SS, LS, S);  // most expensive in this function
	}
	return M_x;
}

void ConstraintJoint::computeStiffnessDampingJoint(Eigen::VectorXd & LHSx, const Eigen::VectorXd & qdot, const Eigen::VectorXd & q, std::unique_ptr<StateSolve>& SS, const std::unique_ptr<LinkageSystem>& LS, const std::unique_ptr<State>& S)
{
	int jindex;
	int const_index;
	for (int jmapi = 0; jmapi < LS->joints.size(); ++jmapi)
	{
		jindex = S->jindex[jmapi];
		const_index = S->constraint_index[jindex];

		for (int cnum = 0; cnum < S->constraintNum[jindex]; ++cnum)
		{
			// Mr = Mr + (h * Dr - h * h * Kr) * qdot
			LHSx[const_index + cnum] += ((SS->h * S->d[jindex] + SS->h * SS->h * S->k[jindex]) * qdot[const_index + cnum]);

			// update on SS->fr
			SS->fr[const_index + cnum] -= S->k[jindex] * (q[const_index + cnum] - S->q0[jindex][cnum]);
		}
	}
}

Eigen::VectorXd ConstraintJoint::computeJT_x(const Eigen::VectorXd & x, const std::unique_ptr<StateSolve>& SS, const std::unique_ptr<LinkageSystem>& LS, const std::unique_ptr<State>& S, std::shared_ptr<State::local_mt> lmt)
{
	Eigen::VectorXd JT_x = Eigen::VectorXd::Zero(SS->Mr_dimension);
	Eigen::VectorXd y = Eigen::VectorXd::Zero(x.size());
	Eigen::VectorXd res;

	int parentIndex;
	int jindex;
	for (int jmapi = (int)LS->joints.size() - 1; jmapi >= 0; --jmapi)
	{
		jindex = S->jindex[LS->joint_map[jmapi]];

		// initialize y so MJx is not changed out of scope
		y.segment<6>(jindex * 6) = x.segment<6>(jindex * 6);

		// yi += child_alpha
		int child_jindex;
		for (int ci = 0; ci < S->childIndex[jindex].size(); ++ci)
		{
			child_jindex = S->childIndex[jindex][ci];
			if (lmt == nullptr)
				y.segment<6>(jindex * 6) += S->alpha[child_jindex];
			else
				y.segment<6>(jindex * 6) += lmt->alpha_[child_jindex];
		}

		if (!S->root[jindex])
		{
			parentIndex = S->parentIndex[jindex];

			// set alpha for this joint
			if (lmt == nullptr)
				S->alpha[jindex] = S->ad_ipT[jindex] * y.segment<6>(jindex * 6);
			else
				lmt->alpha_[jindex] = S->ad_ipT[jindex] * y.segment<6>(jindex * 6);
		}

		res = S->ST_adijT[jindex] * y.segment<6>(jindex * 6);  // most expensive in this fucntion

		assert(res.size() == S->constraintNum[jindex]);
		for (int resi = 0; resi < res.size(); ++resi)
		{
			JT_x[S->constraint_index[jindex] + resi] = res[resi];
		}
	}
	return JT_x;
}

void ConstraintJoint::preprocess_PCG_preconditioner(std::unique_ptr<StateSolve>& SS, std::unique_ptr<LinkageSystem>& LS, const std::unique_ptr<State>& S, bool blkdiag)
{
	// backward traversal
	Matrix6d accum_adj;
	int jindex;
	for (int jmapi = (int)LS->joints.size() - 1; jmapi >= 0; --jmapi)
	{
		jindex = S->jindex[LS->joint_map[jmapi]];

		// Maximal and Reduced damping terms
		// updated outside preconditioner

		accum_adj = Matrix6d::Zero();
		int child_jindex;
		for (int ci = 0; ci < S->childIndex[jindex].size(); ++ci)
		{
			child_jindex = S->childIndex[jindex][ci];
			accum_adj = accum_adj + S->ad_jpT[child_jindex] * S->Pi[child_jindex] * S->ad_jp[child_jindex];
		}

		// toggle use of damping and stiffness blkdiag in Minv 
		if (blkdiag)
		{
			S->Mhat[jindex] = S->I_j[jindex] + SS->h * S->ad_ijT[jindex] * (-S->Dmd[jindex] + S->D[jindex] - SS->h * S->Kmd[jindex]) * S->ad_ij[jindex] + accum_adj;
		}
		else
			S->Mhat[jindex] = S->I_j[jindex] + accum_adj;

		// fixed joints do not have index in reduced state
		if (S->type[jindex] != Joint::JType::Fixed)
		{
			S->Psi[jindex] = (S->ST[jindex] * S->Mhat[jindex] * S->S[jindex] + S->Pr[jindex]).inverse();
			S->Pi[jindex] = S->Mhat[jindex] - S->Mhat[jindex] * S->S[jindex] * S->Psi[jindex] * S->ST[jindex] * S->Mhat[jindex];
		}

		// precomputations
		S->ST_Mhat[jindex] = S->ST[jindex] * S->Mhat[jindex];
		S->S_Psi[jindex] = S->S[jindex] * S->Psi[jindex];
		S->Mhat_S_Psi[jindex] = S->Mhat[jindex] * S->S_Psi[jindex];





		//Matrix6d accumsoa = accum_adj;
		//int list_index = (int)LS->joint_map[jmapi];
		//std::shared_ptr<Joint> j = LS->joints[list_index];
		//std::shared_ptr<Block> child = j->childBlock;
		//int childIndex = j->jindex;
		//// Maximal and Reduced damping terms
		//// updated outside preconditioner
		//accum_adj = Matrix6d::Zero();
		//for (int ci = 0; ci < child->c_joints.size(); ++ci)
		//{
		//	std::shared_ptr<Joint> j_child = child->c_joints[ci];
		//	accum_adj = accum_adj + j_child->ad_jpT * j_child->Pi * j_child->ad_jp;
		//}
		//// toggle use of damping and stiffness blkdiag in Minv 
		//if (blkdiag)
		//{
		//	j->Mhat = j->I_j + SS->h * j->ad_ijT * (-j->Dmd + child->D - SS->h * j->Kmd) * j->ad_ij + accum_adj;
		//}
		//else
		//	j->Mhat = j->I_j + accum_adj;
		//// fixed joints do not have index in reduced state
		//if (j->type != Joint::JType::Fixed)
		//{
		//	j->Psi = (j->ST * j->Mhat * j->S + j->Pr).inverse();
		//	j->Pi = j->Mhat - j->Mhat * j->S * j->Psi * j->ST * j->Mhat;
		//}
		//// precomputations
		//j->ST_Mhat = j->ST * j->Mhat;
		//j->S_Psi = j->S * j->Psi;
		//j->Mhat_S_Psi = j->Mhat * j->S_Psi;

		//std::cout << j->childName << "---------------------------------------" << std::endl;
		//std::cout << "Dmd:" << std::endl;
		//for (int i = 0; i < j->Dmd.rows(); ++i)
		//{
		//	for (int k = 0; k < j->Dmd.cols(); ++k)
		//	{
		//		double diff = j->Dmd(i, k) - S->Dmd[jindex](i, k);
		//		if (isnan(diff))
		//			std::cout << "NaN ";
		//		else if (std::abs(diff) > 1e-8)
		//			std::cout << diff << "      ";
		//		//std::cout << Jdot_x[i] << " - " << Jdot_x_SoA[i] << " = " << diff << "      ";
		//	}
		//}
		//std::cout << std::endl;
		//std::cout << "Kmd:" << std::endl;
		//for (int i = 0; i < j->Kmd.rows(); ++i)
		//{
		//	for (int k = 0; k < j->Kmd.cols(); ++k)
		//	{
		//		double diff = j->Kmd(i, k) - S->Kmd[jindex](i, k);
		//		if (isnan(diff))
		//			std::cout << "NaN ";
		//		else if (std::abs(diff) > 1e-8)
		//			std::cout << diff << "      ";
		//		//std::cout << Jdot_x[i] << " - " << Jdot_x_SoA[i] << " = " << diff << "      ";
		//	}
		//}
		//std::cout << std::endl;
		//std::cout << "ad_ij:" << std::endl;
		//for (int i = 0; i < j->ad_ij.rows(); ++i)
		//{
		//	for (int k = 0; k < j->ad_ij.cols(); ++k)
		//	{
		//		double diff = j->ad_ij(i, k) - S->ad_ij[jindex](i, k);
		//		if (isnan(diff))
		//			std::cout << "NaN ";
		//		else if (std::abs(diff) > 1e-8)
		//			std::cout << diff << "      ";
		//		//std::cout << Jdot_x[i] << " - " << Jdot_x_SoA[i] << " = " << diff << "      ";
		//	}
		//}
		//std::cout << std::endl;
		//std::cout << "ad_ijT:" << std::endl;
		//for (int i = 0; i < j->ad_ijT.rows(); ++i)
		//{
		//	for (int k = 0; k < j->ad_ijT.cols(); ++k)
		//	{
		//		double diff = j->ad_ijT(i, k) - S->ad_ijT[jindex](i, k);
		//		if (isnan(diff))
		//			std::cout << "NaN ";
		//		else if (std::abs(diff) > 1e-8)
		//			std::cout << diff << "      ";
		//		//std::cout << Jdot_x[i] << " - " << Jdot_x_SoA[i] << " = " << diff << "      ";
		//	}
		//}
		//std::cout << std::endl;
		//std::cout << "I_j:" << std::endl;
		//for (int i = 0; i < j->I_j.rows(); ++i)
		//{
		//	for (int k = 0; k < j->I_j.cols(); ++k)
		//	{
		//		double diff = j->I_j(i, k) - S->I_j[jindex](i, k);
		//		if (isnan(diff))
		//			std::cout << "NaN ";
		//		else if (std::abs(diff) > 1e-8)
		//			std::cout << diff << "      ";
		//		//std::cout << Jdot_x[i] << " - " << Jdot_x_SoA[i] << " = " << diff << "      ";
		//	}
		//}
		//std::cout << std::endl;
		//std::cout << "accum_adj:" << std::endl;
		//for (int i = 0; i < accum_adj.rows(); ++i)
		//{
		//	for (int k = 0; k < accum_adj.cols(); ++k)
		//	{
		//		double diff = accum_adj(i, k) - accumsoa(i, k);
		//		if (isnan(diff))
		//			std::cout << "NaN ";
		//		else if (std::abs(diff) > 1e-8)
		//			std::cout << diff << "      ";
		//		//std::cout << Jdot_x[i] << " - " << Jdot_x_SoA[i] << " = " << diff << "      ";
		//	}
		//}
		//std::cout << std::endl;
		//std::cout << "Mhat:" << std::endl;
		//for (int i = 0; i < j->Mhat.rows(); ++i)
		//{
		//	for (int k = 0; k < j->Mhat.cols(); ++k)
		//	{
		//		double diff = j->Mhat(i, k) - S->Mhat[jindex](i, k);
		//		if (isnan(diff))
		//			std::cout << "NaN ";
		//		else if (std::abs(diff) > 1e-8)
		//			std::cout << diff << "      ";
		//		//std::cout << Jdot_x[i] << " - " << Jdot_x_SoA[i] << " = " << diff << "      ";
		//	}
		//}
		//std::cout << std::endl;
		//std::cout << "ST_Mhat:" << std::endl;
		//for (int i = 0; i < j->ST_Mhat.rows(); ++i)
		//{
		//	for (int k = 0; k < j->ST_Mhat.cols(); ++k)
		//	{
		//		double diff = j->ST_Mhat(i, k) - S->ST_Mhat[jindex](i, k);
		//		if (isnan(diff))
		//			std::cout << "NaN ";
		//		else if (std::abs(diff) > 1e-8)
		//			std::cout << diff << "      ";
		//		//std::cout << Jdot_x[i] << " - " << Jdot_x_SoA[i] << " = " << diff << "      ";
		//	}
		//}
		//std::cout << std::endl;		
		//std::cout << "S_Psi:" << std::endl;
		//for (int i = 0; i < j->S_Psi.rows(); ++i)
		//{
		//	for (int k = 0; k < j->S_Psi.cols(); ++k)
		//	{
		//		double diff = j->S_Psi(i, k) - S->S_Psi[jindex](i, k);
		//		if (isnan(diff))
		//			std::cout << "NaN ";
		//		else if (std::abs(diff) > 1e-8)
		//			std::cout << diff << "      ";
		//		//std::cout << Jdot_x[i] << " - " << Jdot_x_SoA[i] << " = " << diff << "      ";
		//	}
		//}
		//std::cout << std::endl;
		//std::cout << "Mhat_S_Psi:" << std::endl;
		//for (int i = 0; i < j->Mhat_S_Psi.rows(); ++i)
		//{
		//	for (int k = 0; k < j->Mhat_S_Psi.cols(); ++k)
		//	{
		//		double diff = j->Mhat_S_Psi(i, k) - S->Mhat_S_Psi[jindex](i, k);
		//		if (isnan(diff))
		//			std::cout << "NaN ";
		//		else if (std::abs(diff) > 1e-8)
		//			std::cout << diff << "      ";
		//		//std::cout << Jdot_x[i] << " - " << Jdot_x_SoA[i] << " = " << diff << "      ";
		//	}
		//}
		//std::cout << std::endl;
	}
}

Eigen::VectorXd ConstraintJoint::computeMinv_x(const Eigen::VectorXd & x, const std::unique_ptr<StateSolve>& SS, const std::unique_ptr<LinkageSystem>& LS, const std::unique_ptr<State>& S, std::shared_ptr<State::local_mt> lmt)
{
	int jindex;
	Eigen::VectorXd Minv_x = Eigen::VectorXd::Zero(x.size());

	// backward traversal
	Eigen::VectorXd x_STB;
	for (int jmapi = (int)LS->joints.size() - 1; jmapi >= 0; --jmapi)
	{
		jindex = S->jindex[LS->joint_map[jmapi]];

		// Bhat_j
		Vector6d Bhat = Vector6d::Zero();
		int child_jindex;
		for (int ci = 0; ci < S->childIndex[jindex].size(); ++ci)
		{
			child_jindex = S->childIndex[jindex][ci];
			// beta changes throughout Minv solve
			if (lmt == nullptr)
				Bhat.noalias() += S->ad_jpT[child_jindex] * S->beta[child_jindex];
			else
				Bhat.noalias() += S->ad_jpT[child_jindex] * lmt->beta_[child_jindex];
		}
		if (lmt == nullptr)
		{
			S->Bhat[jindex] = Bhat;
			S->ST_Bhat[jindex] = S->ST[jindex] * Bhat;
			x_STB = S->ST_Bhat[jindex];
		}
		else
		{
			lmt->Bhat_[jindex] = Bhat;
			lmt->ST_Bhat_[jindex] = S->ST[jindex] * Bhat;
			x_STB = lmt->ST_Bhat_[jindex];
		}

		// beta_j
		for (int i = 0; i < S->constraintNum[jindex]; ++i)
		{
			x_STB[i] = x[S->constraint_index[jindex] + i] - x_STB[i];
		}
		Vector6d beta = Bhat;
		if (S->type[jindex] != Joint::JType::Fixed)
		{
			beta.noalias() += S->Mhat_S_Psi[jindex] * x_STB;
		}
		if (lmt == nullptr)
			S->beta[jindex] = beta;
		else
			lmt->beta_[jindex] = beta;




		//lmt = nullptr;

		//int list_index = (int)LS->joint_map[jmapi];
		//std::shared_ptr<Joint> j = LS->joints[list_index];
		//std::shared_ptr<Block> child = j->childBlock;
		//int childIndex = j->jindex;
		//std::cout << j->childName << "---------------------------------------" << std::endl;
		//std::cout << "CHILD CHECK:" << std::endl;
		//// Bhat_j
		//Bhat = Vector6d::Zero();
		//for (int ci = 0; ci < child->c_joints.size(); ++ci)
		//{
		//	// beta changes throughout Minv solve
		//	if (lmt == nullptr)
		//		Bhat.noalias() += child->c_joints[ci]->ad_jpT * child->c_joints[ci]->beta;
		//	else
		//		Bhat.noalias() += child->c_joints[ci]->ad_jpT * lmt->beta_[child->c_joints[ci]->jindex];

		//	int child_jindex = S->childIndex[jindex][ci];
		//	Vector6d tempj = child->c_joints[ci]->ad_jpT * child->c_joints[ci]->beta;
		//	Vector6d tempsoa = S->ad_jpT[child_jindex] * S->beta[child_jindex];
		//	//std::cout << child->c_joints[ci]->ad_jpT << std::endl;
		//	//std::cout << S->ad_jpT[child_jindex] << std::endl;
		//	for (int i = 0; i < j->beta.size(); ++i)
		//	{
		//		double diff = child->c_joints[ci]->beta[i] - S->beta[child_jindex][i];
		//		if (isnan(diff))
		//			std::cout << "NaN ";
		//		else if (std::abs(diff) > 1e-8)
		//			std::cout << ci << ":beta: " << diff << "      ";
		//		//std::cout << Jdot_x[i] << " - " << Jdot_x_SoA[i] << " = " << diff << "      ";
		//	}
		//	std::cout << "   ";
		//	for (int i = 0; i < 6; ++i)
		//	{
		//		for (int k = 0; k < 6; ++k)
		//		{
		//			double diff = child->c_joints[ci]->ad_jpT(i,k) - S->ad_jpT[child_jindex](i,k);
		//			if (isnan(diff))
		//				std::cout << "NaN ";
		//			else if (std::abs(diff) > 1e-8)
		//				std::cout << ci << ":ad_jp: " << diff << "      ";
		//			//std::cout << Jdot_x[i] << " - " << Jdot_x_SoA[i] << " = " << diff << "      ";
		//		}
		//	}
		//	std::cout << "   ";
		//	for (int i = 0; i < j->Bhat.size(); ++i)
		//	{
		//		double diff = tempj[i] - tempsoa[i];
		//		if (isnan(diff))
		//			std::cout << "NaN ";
		//		else if (std::abs(diff) > 1e-8)
		//			std::cout << ci << ":Bhat: " << diff << "      ";
		//		//std::cout << Jdot_x[i] << " - " << Jdot_x_SoA[i] << " = " << diff << "      ";
		//	}
		//	std::cout << "   ";
		//}
		//std::cout << std::endl;
		//if (lmt == nullptr)
		//{
		//	j->Bhat = Bhat;
		//	j->ST_Bhat = j->ST * Bhat;
		//	x_STB = j->ST_Bhat;
		//}
		//else
		//{
		//	lmt->Bhat_[childIndex] = Bhat;
		//	lmt->ST_Bhat_[childIndex] = j->ST * Bhat;
		//	x_STB = lmt->ST_Bhat_[childIndex];
		//}

		//// beta_j
		//for (int i = 0; i < j->constraintNum; ++i)
		//{
		//	x_STB[i] = x[j->constraint_index + i] - x_STB[i];
		//}
		//beta = Bhat;
		//if (j->type != Joint::JType::Fixed)
		//{
		//	beta.noalias() += j->Mhat_S_Psi * x_STB;
		//}
		//if (lmt == nullptr)
		//	j->beta = beta;
		//else
		//	lmt->beta_[childIndex] = beta;

		//std::cout << "Mhat_S_Psi:" << std::endl;
		//for (int i = 0; i < j->Mhat_S_Psi.rows(); ++i)
		//{
		//	for (int k = 0; k < j->Mhat_S_Psi.cols(); ++k)
		//	{
		//		double diff = j->Mhat_S_Psi(i, k) - S->Mhat_S_Psi[jindex](i, k);
		//		if (isnan(diff))
		//			std::cout << "NaN ";
		//		else if (std::abs(diff) > 1e-8)
		//			std::cout << diff << "      ";
		//		//std::cout << Jdot_x[i] << " - " << Jdot_x_SoA[i] << " = " << diff << "      ";
		//	}
		//}
		//std::cout << std::endl;
		//std::cout << "Bhat:" << std::endl;
		//for (int i = 0; i < j->Bhat.size(); ++i)
		//{
		//	double diff = j->Bhat[i] - S->Bhat[jindex][i];
		//	if (isnan(diff))
		//		std::cout << "NaN ";
		//	else if (std::abs(diff) > 1e-8)
		//		std::cout << diff << "      ";
		//	//std::cout << Jdot_x[i] << " - " << Jdot_x_SoA[i] << " = " << diff << "      ";
		//}
		//std::cout << std::endl;

	}

	// forward traversal
	Eigen::VectorXd y;	// qddot
	for (int jmapi = 0; jmapi < LS->joints.size(); ++jmapi)
	{
		jindex = S->jindex[LS->joint_map[jmapi]];

		if (lmt == nullptr)
			x_STB = S->ST_Bhat[jindex];
		else
			x_STB = lmt->ST_Bhat_[jindex];
		for (int i = 0; i < S->constraintNum[jindex]; ++i)
		{
			x_STB[i] = x[S->constraint_index[jindex] + i] - x_STB[i];
		}

		if (lmt == nullptr)
			S->Vdot[jindex] = Vector6d::Zero();
		else
			lmt->Vdot_[jindex] = Vector6d::Zero();
		if (!S->root[jindex])
		{
			if (lmt == nullptr)
			{
				S->Vdot[jindex] = S->ad_jp[jindex] * S->Vdot[S->parentIndex[jindex]];
				x_STB -= S->ST_Mhat[jindex] * S->Vdot[jindex];
			}
			else
			{
				lmt->Vdot_[jindex] = S->ad_jp[jindex] * lmt->Vdot_[S->parentIndex[jindex]];
				x_STB -= S->ST_Mhat[jindex] * lmt->Vdot_[jindex];
			}
		}

		y = S->Psi[jindex] * x_STB;
		if (S->type[jindex] != Joint::JType::Fixed)
		{
			if (lmt == nullptr)
				S->Vdot[jindex] += S->S[jindex]*y;
			else
				lmt->Vdot_[jindex] += S->S[jindex]*y;
		}
		else
		{
			if (lmt == nullptr)
				S->Vdot[jindex] = Vector6d::Zero();
			else
				lmt->Vdot_[jindex] = Vector6d::Zero();
		}

		// load into result vector
		for (int i = 0; i < y.size(); ++i)
		{
			Minv_x[S->constraint_index[jindex] + i] = y[i];
		}
	}
	return Minv_x;
}

void ConstraintJoint::update(std::unique_ptr<StateSolve> &SS, const std::unique_ptr<LinkageSystem> &LS, const std::unique_ptr<State> &S)
{
	switch (simtype)
	{
	case simType::Pardiso:
		updateJacobian(SS, LS, S);
		break;
	default:
		throw std::runtime_error("Solver type has no implementation for joint constraints");
	}
}

void ConstraintJoint::draw(const std::unique_ptr<StateSolve> &SS, const std::unique_ptr<LinkageSystem> & LS, const std::unique_ptr<State> & S)
{
#ifdef ONLINE_MODE
	Eigen::Vector3d ppos, cpos;
	Eigen::Matrix4d E = Eigen::Matrix4d::Identity();
	glm::mat4 temp;
	int parenti, childi;
	std::shared_ptr<Joint> j;

	for (int jindex = 0; jindex < LS->joints.size(); ++jindex)
	{
		j = LS->joints[jindex];

		if (!j->root)
		{
			parenti = LS->find(j->parentName)->joint->jindex;
			E = S->E[parenti];

			ppos = E.block<3, 3>(0, 0)*j->parentPos + E.block<3, 1>(0, 3);
		}
		else
		{
			ppos = j->parentPos;
		}

		glLineWidth(8);
		glBegin(GL_LINES);
		glVertex3f((float)ppos[0], (float)ppos[1], (float)ppos[2]);
		if (j->type == Joint::JType::Hinge)
		{
			ppos = ppos + j->S.block<3, 1>(0, 0) * 0.5;
		}
		else
		{
			ppos = ppos + Eigen::Vector3d(0,0,0.5);
		}
		glVertex3f((float)ppos[0], (float)ppos[1], (float)ppos[2]);
		glEnd();
		glLineWidth(1);

		childi = LS->find(j->childName)->joint->jindex;
		E = S->E[childi];

		cpos = E.block<3, 3>(0, 0)*j->childPos + E.block<3, 1>(0, 3);
		glLineWidth(4);
		glBegin(GL_LINES);
		glVertex3f((float)cpos[0], (float)cpos[1], (float)cpos[2]);
		if (j->type == Joint::JType::Hinge)
		{
			cpos = cpos - j->S.block<3, 1>(0, 0) * 0.5;
		}
		else
		{
			cpos = cpos - Eigen::Vector3d(0, 0, 0.5);
		}
		glVertex3f((float)cpos[0], (float)cpos[1], (float)cpos[2]);
		glEnd();
		glLineWidth(1);
	}
#endif // ONLINE_MODE
}
