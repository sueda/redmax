#include "RigidBodyCreator.h"

#include "Rigid.h"

#include "Block.h"
#include "Constraint.h"
#include "Joint.h"
#include "JointBall.h"
#include "JointFixed.h"
#include "JointHinge.h"
#include "JointPowered.h"
#include "JointPrismatic.h"
#include "JointSlider.h"
#include "JointSpringy.h"
#include "JointUniversal.h"
#include "LinkageSystem.h"
#include "State.h"
#include "StateSolve.h"

#include <iomanip>
#include <fstream>
#include <float.h>


void RigidBodyCreator::loadTree(int n, simType s, std::unique_ptr<LinkageSystem>& LS, std::unique_ptr<State>& S, std::unique_ptr<StateSolve>& SS)
{
	simtype = s;

	double jointStiffness = 1e3;
	double jointDamping = 1e0;
	double springStiffness = 1e5;
	double springDamping = 1e0;
	double bodyDamping = 1e0;

	std::vector<std::shared_ptr<Block>> leaves;

	std::function<void(int level, std::string &parentname, int n, std::unique_ptr<LinkageSystem>& LS, std::unique_ptr<State>& S, std::unique_ptr<StateSolve>& SS)> branch;
	branch = [&branch, &leaves, &jointStiffness, &jointDamping, &bodyDamping, this](int level, std::string &parentname, int n, std::unique_ptr<LinkageSystem>& LS, std::unique_ptr<State>& S, std::unique_ptr<StateSolve>& SS)
	{
		if (level > n)
			return;

		double mass = 1;
		double m_dens = 1;
		double size;

		Eigen::Matrix3d basis;
		std::shared_ptr<Joint> j;
		Joint::JType jointtype = Joint::JType::Hinge;
		Eigen::Vector3d hinge_axis;
		if (level % 2 == 0)
			hinge_axis = Eigen::Vector3d(0, 1, 0);
		else
			hinge_axis = Eigen::Vector3d(0, 0, 1);

		// create the cross block
		double angle;
		std::string crossname = parentname + "x";
		Eigen::Vector3d parent_pos;	
		Eigen::Vector3d pos = Eigen::Vector3d(0, 0, 0);

		// each level should alternate on the x and z axes
		size = 10 * std::pow(0.8, (level - 1));
		//m_dens = mass / (0.5 * 0.5 * size);
		if (level != 1)
		{
			parent_pos = Eigen::Vector3d(1.5 * std::pow(0.95, level - 1), 0, 0);
		}
		else
		{
			parent_pos = Eigen::Vector3d(-1.5 * std::pow(0.95, level - 1), 0, 0);
		}

		angle = 0.5 * RB_M_PI;
		basis = Eigen::AngleAxisd(angle, hinge_axis);
		LS->joints.push_back(std::make_shared<JointHinge>(S, parentname, crossname, parent_pos, pos, false, LS->joints.size(), angle, hinge_axis, jointStiffness, jointDamping));
		addDisplayBlock(LS, S, crossname, m_dens, bodyDamping, size, 0.5, 0.5, parentname, LS->joints[LS->joints.size() - 1], pos, basis, parent_pos);


		// leaves -  a block's name is the binary encoding of it's position
		size = 3 * std::pow(0.95, level);
		std::string parentX = crossname;
		if (level != 1)
		{
			angle = 1.5 * RB_M_PI;
		}
		else
		{
			angle = 0.5 * RB_M_PI;
		}
		pos = Eigen::Vector3d(-1.5 * std::pow(0.95, level), 0, 0);

		// left
		std::string leftname = parentname + "0";
		parent_pos = Eigen::Vector3d(-4.9 * std::pow(0.8, (level - 1)), 0, 0);
		LS->joints.push_back(std::make_shared<JointHinge>(S, parentX, leftname, parent_pos, pos, false, LS->joints.size(), angle, hinge_axis, jointStiffness, jointDamping));
		basis = Eigen::AngleAxisd(angle, hinge_axis);
		addDisplayBlock(LS, S, leftname, m_dens, bodyDamping, size, 0.5, 0.5, parentX, LS->joints[LS->joints.size() - 1], pos, basis, parent_pos);
		if (level == n)
			leaves.push_back(LS->blocks[LS->blocks.size() - 1].second);

		// right
		std::string rightname = parentname + "1";
		parent_pos = Eigen::Vector3d(4.9 * std::pow(0.8, (level - 1)), 0, 0);
		LS->joints.push_back(std::make_shared<JointHinge>(S, parentX, rightname, parent_pos, pos, false, LS->joints.size(), angle, hinge_axis, jointStiffness, jointDamping));
		basis = Eigen::AngleAxisd(angle, hinge_axis);
		addDisplayBlock(LS, S, rightname, m_dens, bodyDamping, size, 0.5, 0.5, parentX, LS->joints[LS->joints.size() - 1], pos, basis, parent_pos);
		if (level == n)
			leaves.push_back(LS->blocks[LS->blocks.size() - 1].second);

		/// removed 1/4/19 - old scene did not exemplify use of preconditioner
		//// elastic "springs"
		//LS->constraints.push_back(std::make_shared<Elastic>(left_block, right_block, 
		//	Eigen::Vector3d(0,0,0), Eigen::Vector3d(0,0,0), 
		//	-1, springStiffness, springDamping));

		//// show the blocks their spring
		//std::shared_ptr<Constraint> constraint = LS->constraints[LS->constraints.size() - 1];
		//left_block->constraints.push_back(constraint);
		//right_block->constraints.push_back(constraint);

		// call recursively
		branch(level + 1, leftname, n, LS, S, SS);

		// call recursively
		branch(level + 1, rightname, n, LS, S, SS);
	};

	Eigen::Matrix3d basis;
	std::string parent_name = "";
	std::string name = "";
	Eigen::Vector3d parent_pos = Eigen::Vector3d::Zero();
	Eigen::Vector3d pos = Eigen::Vector3d(1.5, 0, 0);
	//Eigen::Vector3d pos = Eigen::Vector3d(0, 1.4, 0);
	Eigen::Vector3d hinge_axis = Eigen::Vector3d(0, 0, 1);

	// Make the root joint
	name = "0";
	double size = 3;
	double angle = 0.5*RB_M_PI;
	double m_dens = 1;// / (size * 0.5 * 0.5);
	basis = Eigen::AngleAxisd(angle, hinge_axis);
	LS->joints.push_back(std::make_shared<JointFixed>(S, parent_name, name, parent_pos, pos, true, LS->joints.size(), basis));
	addDisplayBlock(LS, S, name, m_dens, bodyDamping, size, 0.5, 0.5, parent_name, LS->joints[LS->joints.size() - 1], pos, basis, parent_pos);

	// Make the remaining joints
	branch(1, name, n, LS, S, SS);

	// base case 
	if (n == 0)
		leaves.push_back(LS->blocks[LS->blocks.size() - 1].second);

	// add some springs
	std::set<int> idsprings;
	std::shared_ptr<Block> leafi;
	std::shared_ptr<Block> leafj;
	Eigen::Vector3d posi;
	Eigen::Vector3d posj;
	Eigen::Vector3d xl;
	std::shared_ptr<Constraint> constraint;
	for (int i = 0; i < leaves.size(); ++i)
	{
		// for each leaf, add a spring to its closest neighbor on the same X-Z grid line
		// but not to a sibling
		leafi = leaves[i];
		posi = leafi->E_wi0.block<3, 1>(0, 3);
		double minx0 = DBL_MAX;
		double minx0j = 0;
		double minx1 = -DBL_MAX;
		double minx1j = 0;
		double minz0 = DBL_MAX;
		double minz0j = 0;
		double minz1 = -DBL_MAX;
		double minz1j = 0;

		for (int j = 0; j < leaves.size(); ++j)
		{
			if (j == i)
				continue;

			leafj = leaves[j];
			posj = leafj->E_wi0.block<3, 1>(0, 3);
			double dx = posj[0] - posi[0];
			double dz = posj[2] - posi[2];
			if (std::abs(dz) < 1e-2)
			{
				if (dx > 0 && dx < minx0)
				{
					minx0 = dx;
					minx0j = j;
				}
				if (dx < 0 && dx > minx1)
				{
					minx1 = dx;
					minx1j = j;
				}
			}
			if (std::abs(dx) < 1e-2)
			{
				if (dz > 0 && dz < minz0)
				{
					minz0 = dz;
					minz0j = j;
				}
				if (dz < 0 && dz > minz1)
				{
					minz1 = dz;
					minz1j = j;
				}
			}
		}
		xl = Eigen::Vector3d(0.5*leafi->size[0], 0, 0);
		int id;
		if (minx0j != 0)
		{
			leafj = leaves[(int)minx0j];
			if (minx0j < i)
				id = (int)minx0j*1000 + i;
			else
				id = i*1000 + (int)minx0j;

			if (idsprings.find(id) == idsprings.end())
			{
				LS->constraints.push_back(std::make_shared<Elastic>(leafi, leafj,
					xl, xl, -1, springStiffness, springDamping));
				// show the blocks their spring
				constraint = LS->constraints[LS->constraints.size() - 1];
				//leafi->constraints.push_back(constraint);
				//leafj->constraints.push_back(constraint);

				idsprings.insert(idsprings.begin(), id);
			}
		}
		if (minx1j != 0)
		{
			leafj = leaves[(int)minx1j];
			if (minx1j < i)
				id = (int)minx1j * 1000 + i;
			else
				id = i * 1000 + (int)minx1j;

			if (idsprings.find(id) == idsprings.end())
			{
				LS->constraints.push_back(std::make_shared<Elastic>(leafi, leafj,
					xl, xl, -1, springStiffness, springDamping));
				// show the blocks their spring
				constraint = LS->constraints[LS->constraints.size() - 1];
				//leafi->constraints.push_back(constraint);
				//leafj->constraints.push_back(constraint);

				idsprings.insert(idsprings.begin(), id);
			}
		}
		if (minz0j != 0)
		{
			leafj = leaves[(int)minz0j];
			if (minz0j < i)
				id = (int)minz0j * 1000 + i;
			else
				id = i * 1000 + (int)minz0j;

			if (idsprings.find(id) == idsprings.end())
			{
				LS->constraints.push_back(std::make_shared<Elastic>(leafi, leafj,
					xl, xl, -1, springStiffness, springDamping));
				// show the blocks their spring
				constraint = LS->constraints[LS->constraints.size() - 1];
				//leafi->constraints.push_back(constraint);
				//leafj->constraints.push_back(constraint);

				idsprings.insert(idsprings.begin(), id);
			}
		}
		if (minz1j != 0)
		{
			leafj = leaves[(int)minz1j];
			if (minz1j < i)
				id = (int)minz1j * 1000 + i;
			else
				id = i * 1000 + (int)minz1j;

			if (idsprings.find(id) == idsprings.end())
			{
				LS->constraints.push_back(std::make_shared<Elastic>(leafi, leafj,
					xl, xl, -1, springStiffness, springDamping));
				// show the blocks their spring
				constraint = LS->constraints[LS->constraints.size() - 1];
				//leafi->constraints.push_back(constraint);
				//leafj->constraints.push_back(constraint);

				idsprings.insert(idsprings.begin(), id);
			}
		}
	}

	std::cout << LS->constraints.size() << std::endl;

	// add force to corner leaves
	std::pair<double, double> x0z0;
	x0z0 = std::pair<double, double>(DBL_MAX, DBL_MAX);
	std::pair<double, double> x0z1;
	x0z1 = std::pair<double, double>(DBL_MAX, -DBL_MAX);
	std::pair<double, double> x1z0;
	x1z0 = std::pair<double, double>(-DBL_MAX, DBL_MAX);
	std::pair<double, double> x1z1;
	x1z1 = std::pair<double, double>(-DBL_MAX, -DBL_MAX);
	std::shared_ptr<Block> bx0z0 = nullptr;
	std::shared_ptr<Block> bx0z1 = nullptr;
	std::shared_ptr<Block> bx1z0 = nullptr;
	std::shared_ptr<Block> bx1z1 = nullptr;
	for (int i = 0; i < leaves.size(); ++i)
	{
		leafi = leaves[i];
		posi = leafi->E_wi0.block<3, 1>(0, 3);
		double x = posi[0];
		double z = posi[2];
		if (x < x0z0.first || z < x0z0.second)
		{
			x0z0 = std::pair<double, double>(x, z);
			bx0z0 = leafi;
		}
		if (x < x0z1.first || z > x0z1.second)
		{
			x0z1 = std::pair<double, double>(x, z);
			bx0z1 = leafi;
		}
		if (x > x1z0.first || z < x1z0.second)
		{
			x1z0 = std::pair<double, double>(x, z);
			bx1z0 = leafi;
		}
		if (x > x1z1.first || z > x1z1.second)
		{
			x1z1 = std::pair<double, double>(x, z);
			bx1z1 = leafi;
		}
	}
	double mass = 0;
	for (int i = 0; i < LS->blocks.size(); ++i)
	{
		mass = mass + LS->blocks[i].second->mass0;
	}
	// add the point forces
	LS->constraints.push_back(std::make_shared<SpringPoint>(bx0z0, 
		Eigen::Vector3d::Zero(), Eigen::Vector3d(-1, -1, -1), 1e1*mass));
	constraint = LS->constraints[LS->constraints.size() - 1];
		
	LS->constraints.push_back(std::make_shared<SpringPoint>(bx0z1,
		Eigen::Vector3d::Zero(), Eigen::Vector3d(-1, -1, 1), 1e1*mass));
	constraint = LS->constraints[LS->constraints.size() - 1];

	LS->constraints.push_back(std::make_shared<SpringPoint>(bx1z0,
		Eigen::Vector3d::Zero(), Eigen::Vector3d(1, -1, -1), 1e1*mass));
	constraint = LS->constraints[LS->constraints.size() - 1];

	LS->constraints.push_back(std::make_shared<SpringPoint>(bx1z1,
		Eigen::Vector3d::Zero(), Eigen::Vector3d(1, -1, 1),	1e1*mass));
	constraint = LS->constraints[LS->constraints.size() - 1];
}

void RigidBodyCreator::loadBridge(int nbridge, int ntower, simType s, std::unique_ptr<LinkageSystem>& LS, std::unique_ptr<State>& S, std::unique_ptr<StateSolve>& SS)
{
	// change gravity
	SS->grav = Eigen::Vector3d(0, -9.8, 0);

	int numCables = 30;

	// sim params
	//double density = 1e3;
	//double deckStiffness = 1e3;
	//double deckDamping = 1e0;
	//double towerStiffness = 9.2e6;
	//double towerDamping = 1e3;
	//double springStiffness = 1e7;
	//double springDamping = 1e3;
	//double bodyDamping = 1e0;

	// video params
	double density = 1e3;
	double deckStiffness = 1e8;
	double deckDamping = 1e4;
	double towerStiffness = 1e8;
	double towerDamping = 1e3;
	double springStiffness = 2.5e7;
	double springDamping = 1e4;
	double bodyDamping = 1e0;

	//if(ntower > 180)
	//	towerStiffness *= (ntower/180.0);

	double decklength = 24;
	double tower_height = 10;

	// sim dims
	//double ydimd = 0.15;//0.15;
	//double zdimd = 2;// 0.75;

	//double zdimt = 2;// 0.75;
	//double xdimt = 0.5;

	// video dims
	double ydimd = 0.15;//0.15;
	double zdimd = 1.6;// 0.75;

	double zdimt = 1.8;// 0.75;
	double xdimt = 0.5;

	// load the rooted/fixed joints
	Eigen::Matrix3d basis;
	std::shared_ptr<Joint> j;
	Joint::JType jointtype = Joint::JType::Hinge;
	std::string parent_name = "";
	std::string name = "";
	Eigen::Vector3d parent_pos = Eigen::Vector3d(0, 0, 0);
	Eigen::Vector3d pos = Eigen::Vector3d(0, 0, 0);
	Eigen::Vector3d hinge_axis = Eigen::Vector3d(0, 0, 1);

	// Ground
	name = "ground";
	double size = decklength + 1;
	double angle = 0;
	basis = Eigen::AngleAxisd(angle, hinge_axis);
	LS->joints.push_back(std::make_shared<JointFixed>(S, parent_name, name, parent_pos, pos, true, LS->joints.size(), basis));
	addDisplayBlock(LS, S, name, density, bodyDamping, size, 0.5, zdimt, parent_name, LS->joints[0], pos, basis, parent_pos);

	// Anchorage
	angle = 0;
	name = "anchorR";
	parent_name = "ground";
	size = 3;
	pos = Eigen::Vector3d(0, -1.5, 0);
	parent_pos = Eigen::Vector3d(12, 0, 0);
	basis = Eigen::AngleAxisd(angle, hinge_axis);
	LS->joints.push_back(std::make_shared<JointFixed>(S, parent_name, name, parent_pos, pos, false, LS->joints.size(), basis));
	addDisplayBlock(LS, S, name, density, bodyDamping, xdimt, size, zdimt, parent_name, LS->joints[LS->joints.size() - 1], pos, basis, parent_pos);
	std::shared_ptr<Block> anchor_block_right = LS->blocks[LS->blocks.size() - 1].second;

	name = "anchorL";
	parent_name = "ground";
	size = 3;
	pos = Eigen::Vector3d(0, -1.5, 0);
	parent_pos = Eigen::Vector3d(-12, 0, 0);
	basis = Eigen::AngleAxisd(angle, hinge_axis);
	LS->joints.push_back(std::make_shared<JointFixed>(S, parent_name, name, parent_pos, pos, false, LS->joints.size(), basis));
	addDisplayBlock(LS, S, name, density, bodyDamping, xdimt, size, zdimt, parent_name, LS->joints[LS->joints.size() - 1], pos, basis, parent_pos);

	// Piers/towers
	double blocklength = tower_height / ntower;
	Eigen::Vector3d tower_pos = Eigen::Vector3d(0, 0.5*blocklength, 0);
	std::shared_ptr<Block> tower_block_left = LS->blocks[0].second;
	std::shared_ptr<Block> tower_block_right = LS->blocks[0].second;

	angle = 0;
	Eigen::Vector3d right_parent_pos = Eigen::Vector3d(5 * decklength / 24, 0, 0);
	Eigen::Vector3d left_parent_pos = Eigen::Vector3d(-5 * decklength / 24, 0, 0);
	pos = Eigen::Vector3d(0, -0.5*blocklength, 0);
	parent_pos = Eigen::Vector3d(0, 0, 0);
	for (int i = 0; i < ntower; ++i)
	{
		name = "tr" + std::to_string(i);
		LS->joints.push_back(std::make_shared<JointHinge>(S, tower_block_right->joint->childName, name, right_parent_pos, pos, false, LS->joints.size(), angle, hinge_axis, towerStiffness, towerDamping));
		basis = Eigen::AngleAxisd(angle, hinge_axis);
		addDisplayBlock(LS, S, name, density, bodyDamping, xdimt, blocklength, zdimt, tower_block_right->joint->childName, LS->joints[LS->joints.size() - 1], pos, basis, right_parent_pos);
		tower_block_right = LS->blocks[LS->blocks.size() - 1].second;
		right_parent_pos = Eigen::Vector3d(0, 0.5*blocklength, 0);

		name = "tl" + std::to_string(i);
		LS->joints.push_back(std::make_shared<JointHinge>(S, tower_block_left->joint->childName, name, left_parent_pos, pos, false, LS->joints.size(), angle, hinge_axis, towerStiffness, towerDamping));
		basis = Eigen::AngleAxisd(angle, hinge_axis);
		addDisplayBlock(LS, S, name, density, bodyDamping, xdimt, blocklength, zdimt, tower_block_left->joint->childName, LS->joints[LS->joints.size() - 1], pos, basis, left_parent_pos);
		tower_block_left = LS->blocks[LS->blocks.size() - 1].second;
		left_parent_pos = Eigen::Vector3d(0, 0.5*blocklength, 0);
		pos = Eigen::Vector3d(0, -0.5*blocklength, 0);
	}

	// Deck
	double scaleheight = 0.05 * 10 / nbridge;
	blocklength = decklength / (double)nbridge;
	std::vector<double> anglecoeffs(nbridge+1, 0);
	int numangles = (int)std::ceil(nbridge / 2);
	double sumupangles = 0;
	int first_deckpiece = LS->blocks.size();
	for (int i = 0; i < numangles; ++i)
	{
		angle = -std::sinh(scaleheight*((double)(i+1) / (double)(numangles-1)));
		anglecoeffs[numangles - 1 - i] = angle;
		anglecoeffs[nbridge - numangles + 1 + i] = angle;
		sumupangles = sumupangles + angle;
	}
	anglecoeffs[0] = anglecoeffs[0] - sumupangles;
	double arclength = 0;
	angle = 0;
	for (int i = 0; i < nbridge; ++i)
	{
		angle = angle + anglecoeffs[i];
		arclength = arclength + abs(blocklength / cos(angle));
	}
	parent_name = "anchorL";
	blocklength = arclength / nbridge;
	pos = Eigen::Vector3d(-blocklength*0.5, 0, 0);
	parent_pos = Eigen::Vector3d(0, 1.5, 0);
	std::shared_ptr<Block> bridge_block;
	double k = 0;
	double d = 0;
	for (int i = 0; i < nbridge; ++i)
	{
		name = "bp" + std::to_string(i+1);
		angle = anglecoeffs[i];
		if (i != 0)
		{
			k = deckStiffness;
			d = deckDamping;
		}
		LS->joints.push_back(std::make_shared<JointHinge>(S, parent_name, name, parent_pos, pos, false, LS->joints.size(), angle, hinge_axis, k, d));
		basis = Eigen::AngleAxisd(angle, hinge_axis);
		addDisplayBlock(LS, S, name, density, bodyDamping, blocklength, ydimd, zdimd, parent_name, LS->joints[LS->joints.size() - 1], pos, basis, parent_pos);
		bridge_block = LS->blocks[LS->blocks.size() - 1].second;
		parent_name = name;
		pos = Eigen::Vector3d(-blocklength * 0.5, 0, 0);
		parent_pos = Eigen::Vector3d(blocklength*0.5, 0, 0);
		
		//std::shared_ptr<Block> tower_block;
		//if (i < ((int)nbridge / 2))
		//	tower_block = tower_block_left;
		//else
		//	tower_block = tower_block_right;
		//// elastic "springs"-- every block
		//std::shared_ptr<Block> bridge_block = LS->blocks[LS->blocks.size() - 1].second;
		//LS->constraints.push_back(std::make_shared<Elastic>(bridge_block, tower_block,
		//	Eigen::Vector3d(0, 0.5, 0), tower_pos,
		//	-1, springStiffness, springDamping));
	}
	
	// elastic "springs"-- same number of cables
	decklength = blocklength * nbridge;
	double cable_spacing = (decklength)/(numCables+1);
	double next_cable_pos = cable_spacing;
	std::shared_ptr<Block> tower_block;
	for (int i = 0; i < numCables; ++i)
	{
		if (i < ((int)numCables / 2))
			tower_block = tower_block_left;
		else
			tower_block = tower_block_right;

		// find bridge piece at next_cable_pos
		int block_index = floor(next_cable_pos / blocklength);	
		double block_offset = fmod(next_cable_pos, blocklength) - (blocklength / 2.0);
		std::shared_ptr<Block> bridge_block = LS->blocks[first_deckpiece + block_index].second;
		LS->constraints.push_back(std::make_shared<Elastic>(bridge_block, tower_block,
			Eigen::Vector3d(block_offset, 0, 0), tower_pos,
			-1, springStiffness, springDamping));

		next_cable_pos += cable_spacing;
	}

	// Add the loop closing constraint
	LS->constraints.push_back(std::make_shared<CloseHinge>(bridge_block, anchor_block_right,
		Eigen::Vector3d(blocklength * 0.5, 0, 0), Eigen::Vector3d(0, 1.5, 0), Eigen::Vector3d(0, 0, 1), 0, 0, 0));

	//// Add the point force to mimic a car
	//LS->constraints.push_back(std::make_shared<SpringPoint>(LS->blocks[first_deckpiece].second,
	//	Eigen::Vector3d(-blocklength / 2.0, 0, 0), Eigen::Vector3d(0, -1, 0), 4800));

	//// weight of many cars
	//LS->constraints.push_back(std::make_shared<SpringPoint>(LS->blocks[first_deckpiece].second,
	//	Eigen::Vector3d(-blocklength / 2.0, 0, 0), Eigen::Vector3d(0, -1, 0), 4800*500));

	// some weight considerations
	// half-filled steel bridge decking with concrete filling ------  65 lbs per sq.ft.
	// Pont de Normandie bridge (2 cable towers) roughly 7032 ft long and 77 ft wide ~ 90 ratio
	// THIS bridge - 24 / 1.6 ~ 15 ratio -- small bridge
	// call every 1 unit ~ 50 ft
	// THIS bridge - 1200 ft long, 80 ft wide 
	// bridge deck weighs 5200 lbs per foot of length
	// model a larger car weighing 4000 lbs, and 15ft long
	// 15ft of decking ~ 78000lbs, car 4000lbs
	// 15ft of sim'd bridge ~ 0.3units long by 1.6units wide @ 1000 density ~ 480units of weight
	// a car weighing 480 units (mass) would create a force of roughly 4800 units.

	//std::cout << blocklength << std::endl;
}

void RigidBodyCreator::loadSimpleBridge(int nbridge, simType s, std::unique_ptr<LinkageSystem>& LS, std::unique_ptr<State>& S, std::unique_ptr<StateSolve>& SS)
{
	// change gravity
	SS->grav = Eigen::Vector3d(0, -9.8, 0);

	int numCables = 30;

	double density = 1e3;
	double deckStiffness = 1e1;
	double deckDamping = 1e0;
	double springStiffness = 1e7;
	double springDamping = 1e3;
	double bodyDamping = 1e0;

	double decklength = 24;
	double tower_height = 10;

	// load the rooted/fixed joints
	Eigen::Matrix3d basis;
	std::shared_ptr<Joint> j;
	Joint::JType jointtype = Joint::JType::Hinge;
	std::string parent_name = "";
	std::string name = "";
	Eigen::Vector3d parent_pos = Eigen::Vector3d(0, 0, 0);
	Eigen::Vector3d pos = Eigen::Vector3d(0, 0, 0);
	Eigen::Vector3d hinge_axis = Eigen::Vector3d(0, 0, 1);

	// Ground
	name = "ground";
	double size = decklength + 1;
	double angle = 0;
	basis = Eigen::AngleAxisd(angle, hinge_axis);
	LS->joints.push_back(std::make_shared<JointFixed>(S, parent_name, name, parent_pos, pos, true, LS->joints.size(), basis));
	addDisplayBlock(LS, S, name, density, bodyDamping, size, 0.5, 0.5, parent_name, LS->joints[LS->joints.size() - 1], pos, basis, parent_pos);

	// Piers/towers
	double blocklength = tower_height;
	Eigen::Vector3d tower_pos = Eigen::Vector3d(0, 0.5*blocklength, 0);
	std::shared_ptr<Block> tower_block_left = LS->blocks[LS->blocks.size() - 1].second;
	std::shared_ptr<Block> tower_block_right = LS->blocks[LS->blocks.size() - 1].second;

	angle = 0;
	Eigen::Vector3d right_parent_pos = Eigen::Vector3d(5 * decklength / 24, 0, 0);
	Eigen::Vector3d left_parent_pos = Eigen::Vector3d(-5 * decklength / 24, 0, 0);
	pos = Eigen::Vector3d(0, -0.5*blocklength, 0);
	name = "towerright";
	basis = Eigen::AngleAxisd(angle, hinge_axis);
	LS->joints.push_back(std::make_shared<JointFixed>(S, tower_block_right->joint->childName, name, right_parent_pos, pos, false, LS->joints.size(), basis));
	basis = Eigen::AngleAxisd(angle, hinge_axis);
	addDisplayBlock(LS, S, name, density, bodyDamping, 0.5, blocklength, 0.5, tower_block_right->joint->childName, LS->joints[LS->joints.size() - 1], pos, basis, right_parent_pos);
	tower_block_right = LS->blocks[LS->blocks.size() - 1].second;
	right_parent_pos = Eigen::Vector3d(0, 0.5*blocklength, 0);

	name = "towerleft";
	basis = Eigen::AngleAxisd(angle, hinge_axis);
	LS->joints.push_back(std::make_shared<JointFixed>(S, tower_block_left->joint->childName, name, left_parent_pos, pos, false, LS->joints.size(), basis));
	addDisplayBlock(LS, S, name, density, bodyDamping, 0.5, blocklength, 0.5, tower_block_left->joint->childName, LS->joints[LS->joints.size() - 1], pos, basis, left_parent_pos);
	tower_block_left = LS->blocks[LS->blocks.size() - 1].second;
	left_parent_pos = Eigen::Vector3d(0, 0.5*blocklength, 0);
	pos = Eigen::Vector3d(0, -0.5*blocklength, 0);

	// Anchorage
	angle = 0;
	name = "anchorR";
	parent_name = "ground";
	size = 3;
	pos = Eigen::Vector3d(0, -1.5, 0);
	parent_pos = Eigen::Vector3d(12, 0, 0);
	basis = Eigen::AngleAxisd(angle, hinge_axis);
	LS->joints.push_back(std::make_shared<JointFixed>(S, parent_name, name, parent_pos, pos, false, LS->joints.size(), basis));
	addDisplayBlock(LS, S, name, density, bodyDamping, 0.5, size, 0.5, parent_name, LS->joints[LS->joints.size() - 1], pos, basis, parent_pos);
	std::shared_ptr<Block> anchor_block_right = LS->blocks[LS->blocks.size() - 1].second;

	name = "anchorL";
	parent_name = "ground";
	size = 3;
	pos = Eigen::Vector3d(0, -1.5, 0);
	parent_pos = Eigen::Vector3d(-12, 0, 0);
	basis = Eigen::AngleAxisd(angle, hinge_axis);
	LS->joints.push_back(std::make_shared<JointFixed>(S, parent_name, name, parent_pos, pos, false, LS->joints.size(), basis));
	addDisplayBlock(LS, S, name, density, bodyDamping, 0.5, size, 0.5, parent_name, LS->joints[LS->joints.size() - 1], pos, basis, parent_pos);

	// Deck
	double scaleheight = 0.05 * 10 / nbridge;
	blocklength = decklength / (double)nbridge;
	std::vector<double> anglecoeffs(nbridge + 1, 0);
	int numangles = (int)std::ceil(nbridge / 2);
	double sumupangles = 0;
	int first_deckpiece = LS->blocks.size();
	for (int i = 0; i < numangles; ++i)
	{
		angle = -std::sinh(scaleheight*((double)(i + 1) / (double)(numangles - 1)));
		anglecoeffs[numangles - 1 - i] = angle;
		anglecoeffs[nbridge - numangles + 1 + i] = angle;
		sumupangles = sumupangles + angle;
	}
	anglecoeffs[0] = anglecoeffs[0] - sumupangles;
	double arclength = 0;
	angle = 0;
	for (int i = 0; i < nbridge; ++i)
	{
		angle = angle + anglecoeffs[i];
		arclength = arclength + abs(blocklength / cos(angle));
	}
	parent_name = "anchorL";
	blocklength = arclength / nbridge;
	pos = Eigen::Vector3d(-blocklength * 0.5, 0, 0);
	parent_pos = Eigen::Vector3d(0, 1.5, 0);
	std::shared_ptr<Block> bridge_block;
	double k = 0;
	double d = 0;
	for (int i = 0; i < nbridge; ++i)
	{
		name = "bp" + std::to_string(i + 1);
		angle = anglecoeffs[i];
		if (i != 0)
		{
			k = deckStiffness;
			d = deckDamping;
		}
		LS->joints.push_back(std::make_shared<JointHinge>(S, parent_name, name, parent_pos, pos, false, LS->joints.size(), angle, hinge_axis, k, d));
		basis = Eigen::AngleAxisd(angle, hinge_axis);
		addDisplayBlock(LS, S, name, density, bodyDamping, blocklength, 0.5, 0.5, parent_name, LS->joints[LS->joints.size() - 1], pos, basis, parent_pos);
		bridge_block = LS->blocks[LS->blocks.size() - 1].second;
		parent_name = name;
		pos = Eigen::Vector3d(-blocklength * 0.5, 0, 0);
		parent_pos = Eigen::Vector3d(blocklength*0.5, 0, 0);

		//std::shared_ptr<Block> tower_block;
		//if (i < ((int)nbridge / 2))
		//	tower_block = tower_block_left;
		//else
		//	tower_block = tower_block_right;

		//// elastic "springs"
		//std::shared_ptr<Block> bridge_block = LS->blocks[LS->blocks.size() - 1].second;
		//LS->constraints.push_back(std::make_shared<Elastic>(bridge_block, tower_block,
		//	Eigen::Vector3d(0, 0, 0), tower_pos,
		//	-1, springStiffness, springDamping));

		//// show the blocks their spring
		//std::shared_ptr<Constraint> constraint = LS->constraints[LS->constraints.size() - 1];
		//bridge_block->constraints.push_back(constraint);
		//tower_block->constraints.push_back(constraint);
	}

	// elastic "springs"-- same number of cables
	decklength = blocklength * nbridge;
	double cable_spacing = (decklength) / (numCables + 1);
	double next_cable_pos = cable_spacing;
	std::shared_ptr<Block> tower_block;
	for (int i = 0; i < numCables; ++i)
	{
		if (i < ((int)numCables / 2))
			tower_block = tower_block_left;
		else
			tower_block = tower_block_right;

		// find bridge piece at next_cable_pos
		int block_index = floor(next_cable_pos / blocklength);
		double block_offset = fmod(next_cable_pos, blocklength) - (blocklength / 2.0);
		std::shared_ptr<Block> bridge_block = LS->blocks[first_deckpiece + block_index].second;
		LS->constraints.push_back(std::make_shared<Elastic>(bridge_block, tower_block,
			Eigen::Vector3d(block_offset, 0, 0), tower_pos,
			-1, springStiffness, springDamping));

		next_cable_pos += cable_spacing;
	}


	// Add the loop closing constraint
	LS->constraints.push_back(std::make_shared<CloseHinge>(bridge_block, anchor_block_right,
		Eigen::Vector3d(blocklength * 0.5, 0, 0), Eigen::Vector3d(0, 1.5, 0), Eigen::Vector3d(0, 0, 1), 0, 0, 0));
}

void RigidBodyCreator::loadSimpleTree(int n, simType s, std::unique_ptr<LinkageSystem>& LS, std::unique_ptr<State>& S, std::unique_ptr<StateSolve>& SS)
{
	simtype = s;

	double jointStiffness = 1e3;
	double jointDamping = 1e0;
	double springStiffness = 1e5;
	double springDamping = 1e0;
	double bodyDamping = 1e0;

	std::vector<std::shared_ptr<Block>> leaves;

	std::function<void(int level, std::string &parentname, int n, std::unique_ptr<LinkageSystem>& LS, std::unique_ptr<State>& S, std::unique_ptr<StateSolve>& SS)> branch;
	branch = [&branch, &leaves, &jointStiffness, &jointDamping, &bodyDamping, &springStiffness, &springDamping, this](int level, std::string &parentname, int n, std::unique_ptr<LinkageSystem>& LS, std::unique_ptr<State>& S, std::unique_ptr<StateSolve>& SS)
	{
		if (level > n)
			return;

		double mass = 1;
		double m_dens = 1;
		double size;

		double horiz_block_size = 12;
		double horiz_scale = 0.75;

		Eigen::Matrix3d basis;
		std::shared_ptr<Joint> j;
		Joint::JType jointtype = Joint::JType::Hinge;
		Eigen::Vector3d hinge_axis;
		if (level % 2 == 0)
			hinge_axis = Eigen::Vector3d(0, 1, 0);
		else
			hinge_axis = Eigen::Vector3d(0, 0, 1);

		// create the cross block
		double angle;
		std::string crossname = parentname + "x";
		Eigen::Vector3d parent_pos;
		Eigen::Vector3d pos = Eigen::Vector3d(0, 0, 0);

		// each level should alternate on the x and z axes
		size = horiz_block_size * std::pow(horiz_scale, (level - 1));
		//m_dens = mass / (0.5 * 0.5 * size);
		if (level != 1)
		{
			parent_pos = Eigen::Vector3d(1.5 * std::pow(0.95, level - 1), 0, 0);
		}
		else
		{
			parent_pos = Eigen::Vector3d(-1.5 * std::pow(0.95, level - 1), 0, 0);
		}

		angle = 0.5 * RB_M_PI;
		basis = Eigen::AngleAxisd(angle, hinge_axis);
		LS->joints.push_back(std::make_shared<JointHinge>(S, parentname, crossname, parent_pos, pos, false, LS->joints.size(), angle, hinge_axis, jointStiffness, jointDamping));
		addDisplayBlock(LS, S, crossname, m_dens, bodyDamping, size, 0.5, 0.5, parentname, LS->joints[LS->joints.size() - 1], pos, basis, parent_pos);


		// leaves -  a block's name is the binary encoding of it's position
		size = 3 * std::pow(0.95, level);
		std::string parentX = crossname;
		if (level != 1)
		{
			angle = 1.5 * RB_M_PI;
		}
		else
		{
			angle = 0.5 * RB_M_PI;
		}
		pos = Eigen::Vector3d(-1.5 * std::pow(0.95, level), 0, 0);

		// left
		std::string leftname = parentname + "0";
		parent_pos = Eigen::Vector3d(-(horiz_block_size / 2.0 - 0.1) * std::pow(horiz_scale, (level - 1)), 0, 0);
		LS->joints.push_back(std::make_shared<JointHinge>(S, parentX, leftname, parent_pos, pos, false, LS->joints.size(), angle, hinge_axis, jointStiffness, jointDamping));
		basis = Eigen::AngleAxisd(angle, hinge_axis);
		addDisplayBlock(LS, S, leftname, m_dens, bodyDamping, size, 0.5, 0.5, parentX, LS->joints[LS->joints.size() - 1], pos, basis, parent_pos);
		std::shared_ptr<Block> left_block = LS->blocks[LS->blocks.size() - 1].second;
		if (level == n)
		{
			leaves.push_back(left_block);
		}

		// right
		std::string rightname = parentname + "1";
		parent_pos = Eigen::Vector3d((horiz_block_size / 2.0 - 0.1) * std::pow(horiz_scale, (level - 1)), 0, 0);
		LS->joints.push_back(std::make_shared<JointHinge>(S, parentX, rightname, parent_pos, pos, false, LS->joints.size(), angle, hinge_axis, jointStiffness, jointDamping));
		basis = Eigen::AngleAxisd(angle, hinge_axis);
		addDisplayBlock(LS, S, rightname, m_dens, bodyDamping, size, 0.5, 0.5, parentX, LS->joints[LS->joints.size() - 1], pos, basis, parent_pos);
		std::shared_ptr<Block> right_block = LS->blocks[LS->blocks.size() - 1].second;
		if (level == n)
		{
			leaves.push_back(right_block);
		}

		/// removed 1/4/19 - old scene did not exemplify use of preconditioner
		// elastic "springs"
		LS->constraints.push_back(std::make_shared<Elastic>(left_block, right_block, 
			Eigen::Vector3d(0,0,0), Eigen::Vector3d(0,0,0), 
			-1, springStiffness, springDamping));

		// show the blocks their spring
		std::shared_ptr<Constraint> constraint = LS->constraints[LS->constraints.size() - 1];

		// call recursively
		branch(level + 1, leftname, n, LS, S, SS);

		// call recursively
		branch(level + 1, rightname, n, LS, S, SS);
	};

	Eigen::Matrix3d basis;
	std::string parent_name = "";
	std::string name = "";
	Eigen::Vector3d parent_pos = Eigen::Vector3d::Zero();
	Eigen::Vector3d pos = Eigen::Vector3d(1.5, 0, 0);
	//Eigen::Vector3d pos = Eigen::Vector3d(0, 1.4, 0);
	Eigen::Vector3d hinge_axis = Eigen::Vector3d(0, 0, 1);

	// Make the root joint
	name = "0";
	double size = 3;
	double angle = 0.5*RB_M_PI;
	double m_dens = 1;// / (size * 0.5 * 0.5);
	basis = Eigen::AngleAxisd(angle, hinge_axis);
	LS->joints.push_back(std::make_shared<JointFixed>(S, parent_name, name, parent_pos, pos, true, LS->joints.size(), basis));
	addDisplayBlock(LS, S, name, m_dens, bodyDamping, size, 0.5, 0.5, parent_name, LS->joints[LS->joints.size() - 1], pos, basis, parent_pos);

	// Make the remaining joints
	branch(1, name, n, LS, S, SS);

	// base case 
	if (n == 0)
		leaves.push_back(LS->blocks[LS->blocks.size() - 1].second);

	// add some springs
	std::shared_ptr<Block> leafi;
	Eigen::Vector3d posi;
	std::shared_ptr<Constraint> constraint;

	// add force to corner leaves
	std::pair<double, double> x0z0;
	x0z0 = std::pair<double, double>(DBL_MAX, DBL_MAX);
	std::pair<double, double> x0z1;
	x0z1 = std::pair<double, double>(DBL_MAX, -DBL_MAX);
	std::pair<double, double> x1z0;
	x1z0 = std::pair<double, double>(-DBL_MAX, DBL_MAX);
	std::pair<double, double> x1z1;
	x1z1 = std::pair<double, double>(-DBL_MAX, -DBL_MAX);
	std::shared_ptr<Block> bx0z0 = nullptr;
	std::shared_ptr<Block> bx0z1 = nullptr;
	std::shared_ptr<Block> bx1z0 = nullptr;
	std::shared_ptr<Block> bx1z1 = nullptr;
	for (int i = 0; i < leaves.size(); ++i)
	{
		leafi = leaves[i];
		posi = leafi->E_wi0.block<3, 1>(0, 3);
		double x = posi[0];
		double z = posi[2];
		if (x < x0z0.first || z < x0z0.second)
		{
			x0z0 = std::pair<double, double>(x, z);
			bx0z0 = leafi;
		}
		if (x < x0z1.first || z > x0z1.second)
		{
			x0z1 = std::pair<double, double>(x, z);
			bx0z1 = leafi;
		}
		if (x > x1z0.first || z < x1z0.second)
		{
			x1z0 = std::pair<double, double>(x, z);
			bx1z0 = leafi;
		}
		if (x > x1z1.first || z > x1z1.second)
		{
			x1z1 = std::pair<double, double>(x, z);
			bx1z1 = leafi;
		}
	}
	double mass = 0;
	for (int i = 0; i < LS->blocks.size(); ++i)
	{
		mass = mass + LS->blocks[i].second->mass0;
	}
	// add the point forces
	LS->constraints.push_back(std::make_shared<SpringPoint>(bx0z0,
		Eigen::Vector3d::Zero(), Eigen::Vector3d(-1, -1, -1), 1e1*mass));
	constraint = LS->constraints[LS->constraints.size() - 1];

	LS->constraints.push_back(std::make_shared<SpringPoint>(bx0z1,
		Eigen::Vector3d::Zero(), Eigen::Vector3d(-1, -1, 1), 1e1*mass));
	constraint = LS->constraints[LS->constraints.size() - 1];

	LS->constraints.push_back(std::make_shared<SpringPoint>(bx1z0,
		Eigen::Vector3d::Zero(), Eigen::Vector3d(1, -1, -1), 1e1*mass));
	constraint = LS->constraints[LS->constraints.size() - 1];

	LS->constraints.push_back(std::make_shared<SpringPoint>(bx1z1,
		Eigen::Vector3d::Zero(), Eigen::Vector3d(1, -1, 1), 1e1*mass));
	constraint = LS->constraints[LS->constraints.size() - 1];
}

void RigidBodyCreator::loadUmbrella(int n, simType s, std::unique_ptr<LinkageSystem>& LS, std::unique_ptr<State>& S, std::unique_ptr<StateSolve>& SS)
{
	int	nRibBodies = n;
	int nStretcherBodies = (int)std::floor(nRibBodies / 2);
	SS->grav = Eigen::Vector3d(0, 0, 0);
	SS->baumgarte[2] = 0;

	double tubeHeight = 60;
	double tubeRadius = 0.5;
	double runnerHeight = 5;
	double runnerRadius = 1;
	double ribLength = 50;
	double ribRadius = 0.3;
	double stretcherRadius = 0.3;
	double runnerForce = 5e3;
	double ribStiffness = 1e5;
	double ribDamping = 1e0;
	double stretcherStiffness = 1e5;
	double stretcherDamping = 1e0;
	double springStiffness = 1e0;
	double springDamping = 1e0;
	double bodyDamping = 1e0;
	double density = 1e0;

	double windForce = 2;
	//double dropForce = 10;

	// for video
	springStiffness = 6;

	Eigen::Matrix3d basis;
	std::shared_ptr<Joint> j;
	std::string parent_name = "";
	std::string name = "";
	Eigen::Vector3d parent_pos;
	Eigen::Vector3d pos;
	Eigen::Vector3d axis = Eigen::Vector3d(0, 1, 0);

	Eigen::Vector3d body_pos;
	Eigen::Vector3d joint_pos;

	Eigen::Matrix4d ident4 = Eigen::Matrix4d::Identity();
	Eigen::Matrix3d ident3 = Eigen::Matrix3d::Identity();
	Eigen::Matrix4d temp = Eigen::Matrix4d::Identity();

	//std::cout << std::setprecision(16);

	// Tube
	name = "tube";
	double angle = 0.5 * RB_M_PI;
	body_pos = Eigen::Vector3d(0, 0.5*tubeHeight, 0);
	joint_pos = Eigen::Vector3d(0, 0.5*tubeHeight, 0);
	basis = Eigen::AngleAxisd(angle, Eigen::Vector3d::UnitX());
	LS->joints.push_back(std::make_shared<JointFixed>(S, parent_name, name, body_pos, joint_pos, true, LS->joints.size(), ident3));
	//LS->joints.push_back(std::make_shared<JointHinge>(S, parent_name, name, body_pos, joint_pos, true, LS->joints.size(), 0, Eigen::Vector3d(0,0,1), 1e6, 1e2));//ident3));
	temp.block<3, 3>(0, 0) = basis;
	temp.block<3, 1>(0, 3) = body_pos;
	LS->joints[LS->joints.size() - 1]->E_ji0 = temp;
	LS->joints[LS->joints.size() - 1]->E_ji0_ = true;
	LS->joints[LS->joints.size() - 1]->E_pj0 = ident4;
	LS->joints[LS->joints.size() - 1]->E_pj0_ = true;
	addDisplayBlock(LS, S, name, density, bodyDamping, tubeRadius, tubeHeight, 0, parent_name, LS->joints[LS->joints.size() - 1], joint_pos, basis, body_pos, Block::SType::Cylinder);

	// Runner
	name = "runner";
	parent_name = "tube";
	angle = 0.5 * RB_M_PI;
	body_pos = Eigen::Vector3d(0, 0, 0);
	joint_pos = Eigen::Vector3d(0, 10, 0);
	basis = Eigen::AngleAxisd(angle, Eigen::Vector3d::UnitX());
	LS->joints.push_back(std::make_shared<JointPrismatic>(S, parent_name, name, body_pos, joint_pos, false, LS->joints.size(), 0, axis));
	temp.block<3, 3>(0, 0) = basis;
	temp.block<3, 1>(0, 3) = body_pos;
	LS->joints[LS->joints.size() - 1]->E_ji0 = temp;
	LS->joints[LS->joints.size() - 1]->E_ji0_ = true;
	temp.block<3, 3>(0, 0) = ident3;
	temp.block<3, 1>(0, 3) = joint_pos;
	LS->joints[LS->joints.size() - 1]->E_pj0 = temp;
	LS->joints[LS->joints.size() - 1]->E_pj0_ = true;
	addDisplayBlock(LS, S, name, density, bodyDamping, runnerRadius, runnerHeight, 0, parent_name, LS->joints[LS->joints.size() - 1], joint_pos, basis, body_pos, Block::SType::Cylinder);
	std::shared_ptr<Block> runner = LS->blocks[LS->blocks.size() - 1].second;

	// Ribs
	Eigen::Matrix4d R0 = Eigen::Matrix4d::Identity();
	basis = Eigen::AngleAxisd(80 * RB_M_PI / 180, Eigen::Vector3d::UnitX());
	R0.block<3, 3>(0, 0) = basis;
	double ss;
	double theta;
	double angleX = 0;
	double angleY = 0;
	std::vector<int> ribLast;
	ribLast.push_back((int)LS->blocks.size());
	for (int k = 0; k < 8; ++k)
	{
		ss = double(k) / 8.0;
		theta = ss * 2 * RB_M_PI;
		double h = ribLength / nRibBodies;

		name = "rib" + std::to_string(k) + "0";
		parent_name = "tube";
		body_pos = Eigen::Vector3d(0, 0, 0.5*h);
		joint_pos = Eigen::Vector3d(0, tubeHeight, 0);
		basis = Eigen::AngleAxisd(theta, Eigen::Vector3d::UnitY());
		LS->joints.push_back(std::make_shared<JointUniversal>(S, parent_name, name, body_pos, joint_pos, false, LS->joints.size(), 0, 0));
		temp.block<3, 3>(0, 0) = ident3;
		temp.block<3, 1>(0, 3) = body_pos;
		LS->joints[LS->joints.size() - 1]->E_ji0 = temp;
		LS->joints[LS->joints.size() - 1]->E_ji0_ = true;
		temp.block<3, 3>(0, 0) = basis;
		temp.block<3, 1>(0, 3) = joint_pos;
		LS->joints[LS->joints.size() - 1]->E_pj0 = temp * R0;
		LS->joints[LS->joints.size() - 1]->E_pj0_ = true;
		addDisplayBlock(LS, S, name, density, bodyDamping, ribRadius, h, 0, parent_name, LS->joints[LS->joints.size() - 1], joint_pos, basis, body_pos, Block::SType::Cylinder);

		for (int i = 1; i < nRibBodies; ++i)
		{
			name = "rib" + std::to_string(k) + std::to_string(i);
			parent_name = "rib" + std::to_string(k) + std::to_string(i-1);
			body_pos = Eigen::Vector3d(0, 0, 0.5*h);
			joint_pos = Eigen::Vector3d(0, 0, h);
			LS->joints.push_back(std::make_shared<JointUniversal>(S, parent_name, name, body_pos, joint_pos, false, LS->joints.size(), 0, 0, ribStiffness, ribDamping));
			temp.block<3, 3>(0, 0) = ident3;
			temp.block<3, 1>(0, 3) = body_pos;
			LS->joints[LS->joints.size() - 1]->E_ji0 = temp;
			LS->joints[LS->joints.size() - 1]->E_ji0_ = true;
			temp.block<3, 3>(0, 0) = ident3;
			temp.block<3, 1>(0, 3) = joint_pos;
			LS->joints[LS->joints.size() - 1]->E_pj0 = temp;
			LS->joints[LS->joints.size() - 1]->E_pj0_ = true;
			addDisplayBlock(LS, S, name, density, bodyDamping, ribRadius, h, 0, parent_name, LS->joints[LS->joints.size() - 1], joint_pos, basis, body_pos, Block::SType::Cylinder);
		}
		ribLast.push_back((int)LS->blocks.size());
	}
	// Computes body positions
	for (int i = 0; i < LS->joints.size(); ++i)
	{
		LS->joints[i]->update_nostatechange(LS);
		//std::shared_ptr<Block> b = LS->find(LS->joints[i]->childName);
		//std::cout << LS->joints[i]->childName << std::endl;
		//std::cout << b->E_wi0_inc << std::endl;
	}
	// Stretchers
	Eigen::Vector3d x0 = runner->E_wi0_inc.block<3, 1>(0, 3);
	basis = Eigen::AngleAxisd(-RB_M_PI * 0.5, Eigen::Vector3d::UnitX());
	R0.block<3, 3>(0, 0) = basis;

	Eigen::Matrix4d R1 = Eigen::Matrix4d::Identity();
	Eigen::Vector3d x1;
	Eigen::Vector3d dx;
	Eigen::Vector3d ax;
	std::shared_ptr<Block> rib;
	for (int k = 0; k < 8; ++k)
	{
		ss = double(k) / 8.0;
		theta = ss * 2 * RB_M_PI;
		basis = Eigen::AngleAxisd(theta, Eigen::Vector3d::UnitZ());
		R1.block<3,3>(0,0) = basis;
		rib = LS->blocks[(int)(ribLast[k] + std::ceil(nRibBodies/2) - 1)].second;
		//std::cout << rib->joint->childName << std::endl;
		//std::cout << rib->E_wi0_inc << std::endl;
		x1 = rib->E_wi0_inc.block<3, 1>(0, 3);
		dx = x1 - x0;
		double stretcherLength = dx.norm(); 
		dx = dx / stretcherLength;
		ax = dx.cross(Eigen::Vector3d::UnitY());
		ax = ax / ax.norm();
		double ang = -std::acos(dx[1]);
		double h = stretcherLength / nStretcherBodies;

		name = "stretcher" + std::to_string(k) + "c0";
		parent_name = "runner";
		body_pos = Eigen::Vector3d(0, 0, 0.5*h);
		joint_pos = Eigen::Vector3d(0, 0, 0);
		basis = Eigen::AngleAxisd(ang, ax);
		LS->joints.push_back(std::make_shared<JointUniversal>(S, parent_name, name, body_pos, joint_pos, false, LS->joints.size(), 0, 0, 0, 0));
		temp.block<3, 3>(0, 0) = ident3;
		temp.block<3, 1>(0, 3) = body_pos;
		LS->joints[LS->joints.size() - 1]->E_ji0 = temp;
		LS->joints[LS->joints.size() - 1]->E_ji0_ = true;
		temp.block<3, 3>(0, 0) = basis;
		temp.block<3, 1>(0, 3) = joint_pos;
		LS->joints[LS->joints.size() - 1]->E_pj0 = temp * R0 * R1;
		LS->joints[LS->joints.size() - 1]->E_pj0_ = true;
		addDisplayBlock(LS, S, name, density, bodyDamping, stretcherRadius, h, 0, parent_name, LS->joints[LS->joints.size() - 1], joint_pos, basis, body_pos, Block::SType::Cylinder);

		for (int i = 1; i < nStretcherBodies; ++i)
		{
			name = "stretcher" + std::to_string(k) + "c" + std::to_string(i);
			parent_name = "stretcher" + std::to_string(k) + "c" + std::to_string(i - 1);
			body_pos = Eigen::Vector3d(0, 0, 0.5*h);
			joint_pos = Eigen::Vector3d(0, 0, h);
			basis = Eigen::AngleAxisd(theta, Eigen::Vector3d::UnitY());
			LS->joints.push_back(std::make_shared<JointUniversal>(S, parent_name, name, body_pos, joint_pos, false, LS->joints.size(), 0, 0, stretcherStiffness, stretcherDamping));
			temp.block<3, 3>(0, 0) = ident3;
			temp.block<3, 1>(0, 3) = body_pos;
			LS->joints[LS->joints.size() - 1]->E_ji0 = temp;
			LS->joints[LS->joints.size() - 1]->E_ji0_ = true;
			temp.block<3, 3>(0, 0) = ident3;
			temp.block<3, 1>(0, 3) = joint_pos;
			LS->joints[LS->joints.size() - 1]->E_pj0 = temp;
			LS->joints[LS->joints.size() - 1]->E_pj0_ = true;
			addDisplayBlock(LS, S, name, density, bodyDamping, stretcherRadius, h, 0, parent_name, LS->joints[LS->joints.size() - 1], joint_pos, basis, body_pos, Block::SType::Cylinder);
		}

		// Loop closure
		LS->constraints.push_back(std::make_shared<CloseUniversal>(LS->blocks[LS->blocks.size()-1].second, rib,
			Eigen::Vector3d(0, 0, 0.5*h), Eigen::Vector3d(0, 0, 0)));

	}
	// Springs
	int k1;
	for (int k = 0; k < 8; ++k)
	{
		k1 = k + 1;
		if (k1 == 8)
			k1 = 0;

		std::shared_ptr<Block> body0;
		std::shared_ptr<Block> body1;
		for (int i = 0; i < nRibBodies; ++i)
		{
			body0 = LS->blocks[ribLast[k] + i].second;
			body1 = LS->blocks[ribLast[k1] + i].second;
			// elastic "springs"
			std::shared_ptr<Block> bridge_block = LS->blocks[LS->blocks.size() - 1].second;
			LS->constraints.push_back(std::make_shared<Elastic>(body0, body1,
				Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 0, 0),
				-1, springStiffness, springDamping));

		}
	}

	//// OPTIONAL wind force
	//for (int i = 0; i < LS->blocks.size(); ++i)
	//{
	//	LS->constraints.push_back(std::make_shared<SpringPoint>(LS->blocks[i].second, Eigen::Vector3d::Zero(), Eigen::Vector3d(1, 0, 0), windForce));
	//}

	// Interesting force
	LS->constraints.push_back(std::make_shared<SpringPoint>(LS->blocks[ribLast[1]-1].second, Eigen::Vector3d::Zero(), Eigen::Vector3d(0, -1, -1), 0));

	// Opening force
	LS->constraints.push_back(std::make_shared<SpringPoint>(runner, Eigen::Vector3d::Zero(), Eigen::Vector3d(0, 1, 0), runnerForce));
}

void RigidBodyCreator::loadTest(int n, simType s, std::unique_ptr<LinkageSystem>& LS, std::unique_ptr<State>& S, std::unique_ptr<StateSolve>& SS)
{
	if (n == 1)
	{
		Eigen::Matrix3d basis;
		Joint::JType jointtype = Joint::JType::Hinge;
		std::string parent_name = "";
		std::string name = "";
		Eigen::Vector3d parent_pos = Eigen::Vector3d(-3, 0, 0);
		Eigen::Vector3d pos = Eigen::Vector3d(-1.5, 0, 0);
		//Eigen::Vector3d pos = Eigen::Vector3d(0, 1.4, 0);
		Eigen::Vector3d hinge_axis = Eigen::Vector3d(0, 0, 1);
		double angle = 0 * RB_M_PI;

		double bodyDamping = 0;

		// Make the root joint
		name = "root";
		basis = Eigen::AngleAxisd(angle, hinge_axis);
		LS->joints.push_back(std::make_shared<JointFixed>(S, parent_name, name, parent_pos, pos, true, LS->joints.size(), basis));
		addDisplayBlock(LS, S, name, 3, bodyDamping, 3, 0.5, 0.5, parent_name, LS->joints[LS->joints.size() - 1], pos, basis, parent_pos);
		std::shared_ptr<Block> left_block = LS->blocks[LS->blocks.size() - 1].second;

		angle = 1.5 * RB_M_PI;
		parent_pos = Eigen::Vector3d(1.5, 0, 0);
		parent_name = name;
		name = "middle";
		LS->joints.push_back(std::make_shared<JointHinge>(S, parent_name, name, parent_pos, pos, false, LS->joints.size(), angle, hinge_axis));
		basis = Eigen::AngleAxisd(angle, hinge_axis);
		addDisplayBlock(LS, S, name, 3, bodyDamping, 3, 0.5, 0.5, parent_name, LS->joints[LS->joints.size() - 1], pos, basis, parent_pos);
		std::shared_ptr<Block> right_block = LS->blocks[LS->blocks.size() - 1].second;

		angle = 0.5 * RB_M_PI;
		parent_name = name;
		name = "bottom";
		LS->joints.push_back(std::make_shared<JointHinge>(S, parent_name, name, parent_pos, pos, false, LS->joints.size(), angle, hinge_axis, 0, 1e3));
		basis = Eigen::AngleAxisd(angle, hinge_axis);
		addDisplayBlock(LS, S, name, 3, bodyDamping, 3, 0.5, 0.5, parent_name, LS->joints[LS->joints.size() - 1], pos, basis, parent_pos);

		//// elastic "springs"
		std::string namea = "root";
		std::string nameb = "middle";
		Eigen::Vector3d posa = Eigen::Vector3d(-0.6, 0, 0);
		Eigen::Vector3d posb = Eigen::Vector3d(0.6, 0, 0);

		LS->constraints.push_back(std::make_shared<Elastic>(left_block, right_block,
			posa, posb,
			-1, 1e5, 1e3));

		// show the blocks their spring
		std::shared_ptr<Constraint> constraint = LS->constraints[LS->constraints.size() - 1];
		//left_block->constraints.push_back(constraint);
		//right_block->constraints.push_back(constraint);
	}
	else if (n == 2)
	{
		Eigen::Matrix3d basis;
		std::string parent_name = "";
		std::string name = "";
		Eigen::Vector3d parent_pos = Eigen::Vector3d(0, 0, 0);
		Eigen::Vector3d pos = Eigen::Vector3d(0, 5, 0);
		//Eigen::Vector3d pos = Eigen::Vector3d(0, 1.4, 0);
		Eigen::Vector3d hinge_axis = Eigen::Vector3d(0, 0, 1);
		double angle = 0 * RB_M_PI;

		double bodyDamping = 0;
		double density = 1;
		double x_width = 1;
		double y_width = 10;
		double z_width = 1;
		double angleX = 0;
		double angleY = 0;
		bool root = true;

		int n_bodies = 3;
		for (int i = 0; i < 3; ++i)
		{
			if ((i % 2) == 0)
			{
				angleX = RB_M_PI / 8;
				angleY = 0;
			}
			else
			{
				angleX = 0;
				angleY = RB_M_PI / 8;
			}

			name = "body" + std::to_string(i+1);
			double c1 = std::cos(angleX);
			double c2 = std::cos(angleY);
			double s1 = std::sin(angleX);
			double s2 = std::sin(angleY);
			Eigen::Matrix3d Q = Eigen::Matrix3d::Identity();
			Q(0, 0) = c2;
			Q(0, 1) = s1 * c2;
			Q(0, 2) = -c1 * s2;
			Q(1, 0) = 0;
			Q(1, 1) = c1;
			Q(1, 2) = s1;
			Q(2, 0) = s2;
			Q(2, 1) = -s1 * c2;
			Q(2, 2) = c1 * c2;
			basis = Q;

			basis = Eigen::AngleAxisd(angleX, Eigen::Vector3d::UnitX()) * Eigen::AngleAxisd(angleY, Eigen::Vector3d::UnitY());
			std::cout << basis << std::endl << std::endl;
			LS->joints.push_back(std::make_shared<JointUniversal>(S, parent_name, name, parent_pos, pos, root, LS->joints.size(), angleX, angleY));
			addDisplayBlock(LS, S, name, density, bodyDamping, x_width, y_width, z_width, parent_name, LS->joints[LS->joints.size() - 1], pos, basis, parent_pos);
		
			root = false;
			parent_name = name;
			name = "body" + std::to_string(i);
			parent_pos = Eigen::Vector3d(0, -5, 0);

		}
	}
	else if (n == 3)
	{
		Eigen::Matrix3d basis;
		std::string parent_name = "";
		std::string name = "";
		Eigen::Vector3d parent_pos = Eigen::Vector3d(0, 0, 0);
		Eigen::Vector3d pos = Eigen::Vector3d(0, 5, 0);
		//Eigen::Vector3d pos = Eigen::Vector3d(0, 1.4, 0);
		Eigen::Vector3d axis = Eigen::Vector3d(0, 1, 0);
		double angle = 0 * RB_M_PI;

		double bodyDamping = 0;
		double density = 1;
		double x_width = 1;
		double y_width = 10;
		double z_width = 1;
		double angleX = 0;
		double angleY = 0;
		bool root = true;

		int n_bodies = 3;
		for (int i = 0; i < 3; ++i)
		{
			name = "body" + std::to_string(i + 1);
			basis = Eigen::Matrix3d::Identity();
			if (i != 2)
			{
				LS->joints.push_back(std::make_shared<JointUniversal>(S, parent_name, name, parent_pos, pos, root, LS->joints.size(), angleX, angleY));
				addDisplayBlock(LS, S, name, density, bodyDamping, x_width, y_width, z_width, parent_name, LS->joints[LS->joints.size() - 1], pos, basis, parent_pos);
			}
			else
			{
				LS->joints.push_back(std::make_shared<JointPrismatic>(S, parent_name, name, parent_pos, pos, root, LS->joints.size(), angle, axis));
				addDisplayBlock(LS, S, name, density, bodyDamping, x_width, y_width, z_width, parent_name, LS->joints[LS->joints.size() - 1], pos, basis, parent_pos);
			}

			root = false;
			parent_name = name;
			name = "body" + std::to_string(i);
			parent_pos = Eigen::Vector3d(0, -5, 0);
		}
	}
	else if (n == 4)
	{
		Eigen::Matrix3d basis;
		std::string parent_name = "";
		std::string name = "";
		Eigen::Vector3d parent_pos = Eigen::Vector3d(0, 0, 0);
		Eigen::Vector3d pos = Eigen::Vector3d(0, 0, 0);
		//Eigen::Vector3d hinge_axis = Eigen::Vector3d(0, 0, 1);
		double angle = 0 * RB_M_PI;
		double angleX = 0;
		double angleY = 0;
		double bodyDamping = 0;
		double density = 1;

		SS->grav = Eigen::Vector3d(0,0,-980);

		// Make the root joint
		name = "root";
		basis = Eigen::Matrix3d::Identity();
		LS->joints.push_back(std::make_shared<JointFixed>(S, parent_name, name, parent_pos, pos, true, LS->joints.size(), basis));
		addDisplayBlock(LS, S, name, density, bodyDamping, 20, 1, 1, parent_name, LS->joints[LS->joints.size() - 1], pos, basis, parent_pos);

		parent_pos = Eigen::Vector3d(-10, 0, 0);
		pos = Eigen::Vector3d(0, 0, 5);
		parent_name = "root";
		name = "b1";
		LS->joints.push_back(std::make_shared<JointUniversal>(S, parent_name, name, parent_pos, pos, false, LS->joints.size(), angleX, angleY));
		addDisplayBlock(LS, S, name, density, bodyDamping, 1, 1, 10, parent_name, LS->joints[LS->joints.size() - 1], pos, basis, parent_pos);
		
		parent_pos = Eigen::Vector3d(10, 0, 0);
		pos = Eigen::Vector3d(0, 0, 5);
		parent_name = "root";
		name = "b2";
		LS->joints.push_back(std::make_shared<JointUniversal>(S, parent_name, name, parent_pos, pos, false, LS->joints.size(), angleX, angleY));
		addDisplayBlock(LS, S, name, density, bodyDamping, 1, 1, 10, parent_name, LS->joints[LS->joints.size() - 1], pos, basis, parent_pos);
		std::shared_ptr<Block> block3 = LS->blocks[LS->blocks.size() - 1].second;

		parent_pos = Eigen::Vector3d(0, 0, -5);
		pos = Eigen::Vector3d(-10, 0, 0);
		parent_name = "b1";
		name = "b3";
		LS->joints.push_back(std::make_shared<JointUniversal>(S, parent_name, name, parent_pos, pos, false, LS->joints.size(), angleX, angleY));
		addDisplayBlock(LS, S, name, density, bodyDamping, 20, 1, 1, parent_name, LS->joints[LS->joints.size() - 1], pos, basis, parent_pos);
		std::shared_ptr<Block> block4 = LS->blocks[LS->blocks.size() - 1].second;

		parent_pos = Eigen::Vector3d(0, 0, 0);
		pos = Eigen::Vector3d(0, 0, 5);
		parent_name = "b3";
		name = "b4";
		LS->joints.push_back(std::make_shared<JointUniversal>(S, parent_name, name, parent_pos, pos, false, LS->joints.size(), angleX, angleY, 1, 1));
		addDisplayBlock(LS, S, name, density, bodyDamping, 1, 1, 10, parent_name, LS->joints[LS->joints.size() - 1], pos, basis, parent_pos);

		// Add the loop closing constraint
		LS->constraints.push_back(std::make_shared<CloseUniversal>(block3, block4,
			Eigen::Vector3d(0, 0, -5), Eigen::Vector3d(10, 0, 0)));
	}
	else if (n == 5)
	{
		Eigen::Matrix3d basis;
		std::string parent_name = "";
		std::string name = "";
		Eigen::Vector3d parent_pos = Eigen::Vector3d(0, 0, 0);
		Eigen::Vector3d pos = Eigen::Vector3d(-0.5, 0, 1.5);
		Eigen::Vector3d hinge_axis = Eigen::Vector3d(0, 1, 0);
		double angle = 0 * RB_M_PI;
		double bodyDamping = 0;
		double density = 1;

		SS->grav = Eigen::Vector3d(0, 0, -980);

		// Make the root joint
		name = "root";
		basis = Eigen::Matrix3d::Identity();
		LS->joints.push_back(std::make_shared<JointHinge>(S, parent_name, name, parent_pos, pos, true, LS->joints.size(), angle, hinge_axis));
		addDisplayBlock(LS, S, name, density, bodyDamping, 0.5, 3, 1, parent_name, LS->joints[LS->joints.size() - 1], pos, basis, parent_pos, Block::SType::Cylinder);

		parent_pos = Eigen::Vector3d(-0.5, 0, -1.5);
		pos = Eigen::Vector3d(-0.5, 0, 1.5);
		parent_name = "root";
		name = "b1";
		LS->joints.push_back(std::make_shared<JointHinge>(S, parent_name, name, parent_pos, pos, false, LS->joints.size(), angle, hinge_axis));
		addDisplayBlock(LS, S, name, density, bodyDamping, 0.5, 3, 10, parent_name, LS->joints[LS->joints.size() - 1], pos, basis, parent_pos, Block::SType::Cylinder);
	}

}

void RigidBodyCreator::loadChain(int n, simType s, std::unique_ptr<LinkageSystem>& LS, std::unique_ptr<State>& S, std::unique_ptr<StateSolve>& SS)
{
	simtype = s;

	Eigen::Matrix3d basis;
	std::shared_ptr<Joint> j;
	Joint::JType jointtype = Joint::JType::Hinge;
	std::string parent_name = "";
	std::string name = ""; 
	Eigen::Vector3d parent_pos = Eigen::Vector3d::Zero();
	Eigen::Vector3d pos = Eigen::Vector3d(-1.4, 0, 0);
	Eigen::Vector3d hinge_axis = Eigen::Vector3d(0, 0, 1);
	double angle = 1.5 * RB_M_PI;

	double bodyDamping = 0;

	// Make the root joint
	int i = 0;
	name = "block" + std::to_string(i);
	LS->joints.push_back(std::make_shared<JointHinge>(S, parent_name, name, parent_pos, pos, true, LS->joints.size(), angle, hinge_axis));
	//std::cout << name << " " << LS->joints.size() << std::endl;
	basis = Eigen::AngleAxisd(angle, Eigen::Vector3d::UnitZ());
	addDisplayBlock(LS, S, name, 3, bodyDamping, 3, 0.8, 0.5, parent_name, LS->joints[LS->joints.size() - 1], pos, basis, parent_pos);

	// Make the remaining joints
	parent_pos = Eigen::Vector3d(1.4, 0, 0);
	angle = 0;
	for (i = 1; i < n; ++i)
	{
		if (i == n - 1)
			angle = 0.25 * RB_M_PI;

		parent_name = name;
		name = "block" + std::to_string(i);

		LS->joints.push_back(std::make_shared<JointHinge>(S, parent_name, name, parent_pos, pos, false, LS->joints.size(), angle, hinge_axis));
		//std::cout << name << " " << LS->joints.size() << std::endl;
		basis = Eigen::AngleAxisd(angle, Eigen::Vector3d::UnitZ());
		addDisplayBlock(LS, S, name, 3, bodyDamping, 3, 0.8, 0.5, parent_name, LS->joints[LS->joints.size() - 1], pos, basis, parent_pos);
	}
}

void RigidBodyCreator::loadLinkagesfromFile(
	std::unique_ptr<LinkageSystem> & LS,
	std::unique_ptr<State>& S,
	std::unique_ptr<StateSolve>& SS,
	const bool disconnect)
{
	std::ifstream infile;
	std::stringstream buffer;
	infile.open(RESOURCE_DIR + activeLinks);

	// Read the file into the buffer
	buffer << infile.rdbuf();

	// Declarations for link structure
	std::string name = "";
	std::string parent_name = "";
	Eigen::Vector3d pos = Eigen::Vector3d::Zero();
	Eigen::Vector3d parent_pos = Eigen::Vector3d::Zero();
	Eigen::Vector3d physics_pos = Eigen::Vector3d::Zero();
	Eigen::Matrix3d basis;
	double density = 0;
	double xdim = 0;
	double ydim = 0;
	double zdim = 0;
	double angle = 0;
	double physics_density = 0;
	double physics_xdim = 0;
	double physics_ydim = 0;
	double physics_zdim = 0;
	double physics_angle = 0;
	Joint::Joint::JType jointtype = Joint::JType::Ball;
	Eigen::Vector3d hinge_axis = Eigen::Vector3d::Zero();
	double cons1 = 0;
	double cons2 = 0;
	bool invisible = false;
	bool physics = false;
	std::string parent_a_name = "";
	std::string parent_b_name = "";
	Eigen::Vector3d parent_a_pos = Eigen::Vector3d::Zero();
	Eigen::Vector3d parent_b_pos = Eigen::Vector3d::Zero();
	double rest_length = 0;
	double stiffness = 0;
	double damping = 0;
	double body_damping = 0;

	std::string type;
	while (buffer.peek() != EOF)
	{
		if (buffer.peek() == '#')
		{
			buffer.get();
			getline(buffer, type);
			if (type.compare("options") == 0)
			{
				name = "";

				// load the settings
				std::string line;
				while (buffer.peek() != '#' && !buffer.eof())
				{
					getline(buffer, line);

					if (line.substr(0, 4).compare("simt") == 0)
					{
						if (line.length() > 4)
						{
							name = line.substr(5, line.length());

							if (name.compare("pardiso") == 0)
								simtype = simType::Pardiso;
							else if (name.compare("PCG") == 0)
								simtype = simType::PCG;
							else if (name.compare("PCGold") == 0)
								simtype = simType::PCG_unopt;
							/*  ...  */
						}
					}
				}
			}
			else if (type.compare("link") == 0)
			{
				name = "";
				parent_name = "";
				pos = Eigen::Vector3d::Zero();
				parent_pos = Eigen::Vector3d::Zero();
				physics_pos = Eigen::Vector3d::Zero();
				density = 0;
				xdim = 0;
				ydim = 0;
				zdim = 0;
				angle = 0;
				physics_density = 0;
				physics_xdim = 0;
				physics_ydim = 0;
				physics_zdim = 0;
				physics_angle = 0;
				jointtype = Joint::JType::Ball;
				hinge_axis = Eigen::Vector3d::Zero();
				cons1 = 0;
				cons2 = 0;
				invisible = false;
				physics = false;
				stiffness = 0.5;
				damping = 0;
				body_damping = 0;

				// Load a linkage
				std::string line;
				while (buffer.peek() != '#' && !buffer.eof())
				{
					getline(buffer, line);

					if (line.substr(0, 4).compare("name") == 0)
					{
						if (line.length() > 4)
						{
							name = line.substr(5, line.length());
						}
					}
					else if (line.substr(0, 4).compare("prnt") == 0)
					{
						if (line.length() > 4)
						{
							parent_name = line.substr(5, line.length());
						}
					}
					else if (line.substr(0, 4).compare("dens") == 0)
					{
						if (line.length() > 4)
						{
							std::stringstream ss;
							ss.str(line.substr(5, line.length()));
							if (!physics)
							{
								ss >> density;
							}
							else
							{
								ss >> physics_density;
							}
						}
					}
					else if (line.substr(0, 4).compare("size") == 0)
					{
						if (line.length() > 4)
						{
							std::stringstream ss;
							ss.str(line.substr(5, line.length()));
							if (!physics)
							{
								ss >> xdim;
								ss >> ydim;
								ss >> zdim;
							}
							else
							{
								ss >> physics_xdim;
								ss >> physics_ydim;
								ss >> physics_zdim;
							}
						}
					}
					else if (line.substr(0, 4).compare("angl") == 0)
					{
						if (line.length() > 4)
						{
							std::stringstream ss;
							ss.str(line.substr(5, line.length()));
							if (!physics)
							{
								ss >> angle;
							}
							else
							{
								ss >> physics_angle;
							}
						}
					}
					else if (line.substr(0, 4).compare("axis") == 0)
					{
						if (line.length() > 4)
						{
							std::stringstream ss;
							ss.str(line.substr(5, line.length()));
							double x, y, z;
							ss >> x;
							ss >> y;
							ss >> z;
							hinge_axis = Eigen::Vector3d(x, y, z);
						}
					}
					else if (line.substr(0, 4).compare("cpos") == 0)
					{
						if (line.length() > 4)
						{
							std::stringstream ss;
							ss.str(line.substr(5, line.length()));
							double x, y, z;
							ss >> x;
							ss >> y;
							ss >> z;
							if (!physics)
							{
								pos = Eigen::Vector3d(x, y, z);
							}
							else
							{
								physics_pos = Eigen::Vector3d(x, y, z);
							}
						}
					}
					else if (line.substr(0, 4).compare("ppos") == 0)
					{
						if (line.length() > 4)
						{
							std::stringstream ss;
							ss.str(line.substr(5, line.length()));
							double x, y, z;
							ss >> x;
							ss >> y;
							ss >> z;
							parent_pos = Eigen::Vector3d(x, y, z);
						}
					}
					else if (line.substr(0, 4).compare("type") == 0)
					{
						if (line.length() > 4)
						{
							line = line.substr(5, line.length());
							if (line.substr(0, 4).compare("ball") == 0)
							{
								jointtype = Joint::JType::Ball;
							}
							else if (line.substr(0, 5).compare("hinge") == 0)
							{
								jointtype = Joint::JType::Hinge;
								std::stringstream ss;
								ss.str(line.substr(6, line.length()));
								double x, y, z;
								ss >> x;
								ss >> y;
								ss >> z;
								hinge_axis = Eigen::Vector3d(x, y, z);
							}
							else if (line.substr(0, 5).compare("fixed") == 0)
							{
								jointtype = Joint::JType::Fixed;
								std::stringstream ss;
								ss.str(line.substr(6, line.length()));
								double x, y, z;
								ss >> x;
								ss >> y;
								ss >> z;
								hinge_axis = Eigen::Vector3d(x, y, z);
							}
							else if (line.substr(0, 6).compare("slider") == 0)
							{
								jointtype = Joint::JType::Slider;
							}
							else if (line.substr(0, 7).compare("springy") == 0)
							{
								jointtype = Joint::JType::Springy;
								std::stringstream ss;
								ss.str(line.substr(8, line.length()));
								ss >> cons1;
								ss >> cons2;
							}
							else if (line.substr(0, 7).compare("powered") == 0)
							{
								jointtype = Joint::JType::Powered;
								std::stringstream ss;
								ss.str(line.substr(8, line.length()));
								ss >> cons1;
							}
						}
					}
					else if (line.substr(0, 4).compare("stif") == 0)
					{
						if (line.length() > 4)
						{
							std::stringstream ss;
							ss.str(line.substr(5, line.length()));

							ss >> stiffness;
						}
					}
					else if (line.substr(0, 4).compare("damp") == 0)
					{
						if (line.length() > 4)
						{
							std::stringstream ss;
							ss.str(line.substr(5, line.length()));

							ss >> damping;
						}
					}
					else if (line.substr(0, 4).compare("bdmp") == 0)
					{
						if (line.length() > 4)
						{
							std::stringstream ss;
							ss.str(line.substr(5, line.length()));

							ss >> body_damping;
						}
					}
					else if (line.substr(0, 7).compare("physics") == 0)
					{
						physics = true;
					}
					else if (line.substr(0, 9).compare("invisible") == 0)
					{
						invisible = true;
					}
				}

				if (density <= 0 || xdim <= 0 || ydim <= 0 || zdim <= 0)
				{
					std::cout << "Link <" << name << "> is not realistic" << std::endl;
				}
				else
				{
					bool root = false;
					if (parent_name.size() == 0)
						root = true;

					// create and save the joint
					std::shared_ptr<Joint> j;
					basis = Eigen::AngleAxisd(angle, Eigen::Vector3d::UnitZ());
					if (jointtype == Joint::JType::Powered)
					{
                        Eigen::Vector3d cp = Eigen::Vector3d(1.0, 0.0, 0.0);
                        j = std::make_shared<JointPowered>(S, parent_name, name, parent_pos, pos, root, LS->joints.size(), cp); // rotation of joint?
					}
					else if (jointtype == Joint::JType::Springy)
					{
						j = std::make_shared<JointSpringy>(S, parent_name, name, parent_pos, pos, root, LS->joints.size(), cons1, cons2);
					}
					else if (jointtype == Joint::JType::Hinge)
					{
						basis = Eigen::AngleAxisd(angle, hinge_axis);
						j = std::make_shared<JointHinge>(S, parent_name, name, parent_pos, pos, root, LS->joints.size(), angle, hinge_axis, stiffness, damping);
					}
					else if (jointtype == Joint::JType::Fixed)
					{
						basis = Eigen::AngleAxisd(angle, hinge_axis);
						j = std::make_shared<JointFixed>(S, parent_name, name, parent_pos, pos, root, LS->joints.size(), basis);
					}
					else if (jointtype == Joint::JType::Ball)
					{
						j = std::make_shared<JointBall>(S, parent_name, name, parent_pos, pos, root, LS->joints.size());
					}
					else if (jointtype == Joint::JType::Slider)
					{
						j = std::make_shared<JointSlider>(S, parent_name, name, parent_pos, pos, root, LS->joints.size()); // slide axis
					}
					LS->joints.push_back(j);

					// It is good enough for us right now- load the link
					addDisplayBlock(LS, S, name, density, body_damping, xdim, ydim, zdim, parent_name, j, pos, basis, parent_pos);

					if (invisible)
					{
						LS->blocks[LS->blocks.size() - 1].second->invisible = true;
					}

					if (physics && disconnect)
					{
						basis = Eigen::AngleAxisd(physics_angle, Eigen::Vector3d::UnitZ());
						disconnectInertial(LS, S, name, physics_pos, basis, physics_density, physics_xdim, physics_ydim, physics_zdim);
					}
				}
			}
			else if (type.compare("spring") == 0)
			{
				parent_a_name = "";
				parent_b_name = "";
				parent_a_pos = Eigen::Vector3d::Zero();
				parent_b_pos = Eigen::Vector3d::Zero();
				rest_length = 0;
				stiffness = 0;
				damping = 0;

				// Load a spring
				std::string line;
				while (buffer.peek() != '#' && !buffer.eof())
				{
					getline(buffer, line);

					if (line.substr(0, 4).compare("lnka") == 0)
					{
						if (line.length() > 4)
						{
							parent_a_name = line.substr(5, line.length());
						}
					}
					else if (line.substr(0, 4).compare("lnkb") == 0)
					{
						if (line.length() > 4)
						{
							parent_b_name = line.substr(5, line.length());
						}
					}
					else if (line.substr(0, 4).compare("posa") == 0)
					{
						if (line.length() > 4)
						{
							std::stringstream ss;
							ss.str(line.substr(5, line.length()));
							double x, y, z;
							ss >> x;
							ss >> y;
							ss >> z;
							parent_a_pos = Eigen::Vector3d(x, y, z);
						}
					}
					else if (line.substr(0, 4).compare("posb") == 0)
					{
						if (line.length() > 4)
						{
							std::stringstream ss;
							ss.str(line.substr(5, line.length()));
							double x, y, z;
							ss >> x;
							ss >> y;
							ss >> z;
							parent_b_pos = Eigen::Vector3d(x, y, z);
						}
					}
					else if (line.substr(0, 4).compare("rest") == 0)
					{
						if (line.length() > 4)
						{
							std::stringstream ss;
							ss.str(line.substr(5, line.length()));
							ss >> rest_length;
						}
					}
					else if (line.substr(0, 4).compare("stif") == 0)
					{
						if (line.length() > 4)
						{
							std::stringstream ss;
							ss.str(line.substr(5, line.length()));
							ss >> stiffness;
						}
					}
					else if (line.substr(0, 4).compare("damp") == 0)
					{
						if (line.length() > 4)
						{
							std::stringstream ss;
							ss.str(line.substr(5, line.length()));
							ss >> damping;
						}
					}
				}

				std::shared_ptr<Block> pa = LS->find(parent_a_name);
				std::shared_ptr<Block> pb = LS->find(parent_b_name);

				// It is good enough for us to load
				LS->constraints.push_back(std::make_unique<Elastic>(pa, pb, parent_a_pos, parent_b_pos,
					rest_length, stiffness, damping));

				//std::shared_ptr<Constraint> constraint = LS->constraints[LS->constraints.size() - 1];
				//pa->constraints.push_back(constraint);
				//pb->constraints.push_back(constraint);
			}
			else if (type.compare("closejoint") == 0)
			{
				parent_a_name = "";
				parent_b_name = "";
				parent_a_pos = Eigen::Vector3d::Zero();
				parent_b_pos = Eigen::Vector3d::Zero();
				hinge_axis = Eigen::Vector3d::Zero();
				angle = 0;
				stiffness = 0;
				damping = 0;

				// Load a spring
				std::string line;
				while (buffer.peek() != '#' && !buffer.eof())
				{
					getline(buffer, line);

					if (line.substr(0, 4).compare("lnka") == 0)
					{
						if (line.length() > 4)
						{
							parent_a_name = line.substr(5, line.length());
						}
					}
					else if (line.substr(0, 4).compare("lnkb") == 0)
					{
						if (line.length() > 4)
						{
							parent_b_name = line.substr(5, line.length());
						}
					}
					else if (line.substr(0, 4).compare("posa") == 0)
					{
						if (line.length() > 4)
						{
							std::stringstream ss;
							ss.str(line.substr(5, line.length()));
							double x, y, z;
							ss >> x;
							ss >> y;
							ss >> z;
							parent_a_pos = Eigen::Vector3d(x, y, z);
						}
					}
					else if (line.substr(0, 4).compare("posb") == 0)
					{
						if (line.length() > 4)
						{
							std::stringstream ss;
							ss.str(line.substr(5, line.length()));
							double x, y, z;
							ss >> x;
							ss >> y;
							ss >> z;
							parent_b_pos = Eigen::Vector3d(x, y, z);
						}
					}
					else if (line.substr(0, 4).compare("angl") == 0)
					{
						if (line.length() > 4)
						{
							std::stringstream ss;
							ss.str(line.substr(5, line.length()));
							ss >> angle;
						}
					}
					else if (line.substr(0, 4).compare("axis") == 0)
					{
						if (line.length() > 4)
						{
							std::stringstream ss;
							ss.str(line.substr(5, line.length()));
							double x, y, z;
							ss >> x;
							ss >> y;
							ss >> z;
							hinge_axis = Eigen::Vector3d(x, y, z);
						}
					}
					else if (line.substr(0, 4).compare("stif") == 0)
					{
						if (line.length() > 4)
						{
							std::stringstream ss;
							ss.str(line.substr(5, line.length()));
							ss >> stiffness;
						}
					}
					else if (line.substr(0, 4).compare("damp") == 0)
					{
						if (line.length() > 4)
						{
							std::stringstream ss;
							ss.str(line.substr(5, line.length()));
							ss >> damping;
						}
					}
				}

				std::shared_ptr<Block> pa = LS->find(parent_a_name);
				std::shared_ptr<Block> pb = LS->find(parent_b_name);

				// It is good enough for us to load
				LS->constraints.push_back(std::make_unique<CloseHinge>(pa, pb, parent_a_pos, parent_b_pos,
					hinge_axis, angle, stiffness, damping));

				//std::shared_ptr<Constraint> constraint = LS->constraints[LS->constraints.size() - 1];
				//pa->constraints.push_back(constraint);
				//pb->constraints.push_back(constraint);
			}
			else if (type.compare("springpoint") == 0)
			{
				parent_a_name = "";
				parent_a_pos = Eigen::Vector3d::Zero();
				hinge_axis = Eigen::Vector3d::Zero();
				stiffness = 0;

				// Load a spring
				std::string line;
				while (buffer.peek() != '#' && !buffer.eof())
				{
					getline(buffer, line);

					if (line.substr(0, 4).compare("lnka") == 0)
					{
						if (line.length() > 4)
						{
							parent_a_name = line.substr(5, line.length());
						}
					}
					else if (line.substr(0, 4).compare("posa") == 0)
					{
						if (line.length() > 4)
						{
							std::stringstream ss;
							ss.str(line.substr(5, line.length()));
							double x, y, z;
							ss >> x;
							ss >> y;
							ss >> z;
							parent_a_pos = Eigen::Vector3d(x, y, z);
						}
					}
					else if (line.substr(0, 4).compare("dirn") == 0)
					{
						if (line.length() > 4)
						{
							std::stringstream ss;
							ss.str(line.substr(5, line.length()));
							double x, y, z;
							ss >> x;
							ss >> y;
							ss >> z;
							hinge_axis = Eigen::Vector3d(x, y, z);
						}
					}
					else if (line.substr(0, 4).compare("stif") == 0)
					{
						if (line.length() > 4)
						{
							std::stringstream ss;
							ss.str(line.substr(5, line.length()));
							ss >> stiffness;
						}
					}
				}

				std::shared_ptr<Block> pa = LS->find(parent_a_name);
				std::shared_ptr<Block> pb = LS->find(parent_b_name);

				// It is good enough for us to load
				LS->constraints.push_back(std::make_unique<SpringPoint>(pa, parent_a_pos,
					hinge_axis, stiffness));

				//std::shared_ptr<Constraint> constraint = LS->constraints[LS->constraints.size() - 1];
				//pa->constraints.push_back(constraint);
				//pb->constraints.push_back(constraint);
			}
		}
		else
		{
			char s[1024];
			buffer.getline(s, 1024);
		}
	}

	infile.close();
}

void RigidBodyCreator::loadPhysicsfromFile(std::unique_ptr<LinkageSystem> & LS, std::unique_ptr<State>& S)
{
	std::ifstream infile;
	std::stringstream buffer;
	infile.open(RESOURCE_DIR + USER_SOURCE);

	// Read the file into the buffer
	buffer << infile.rdbuf();

	// Declarations for link structure
	std::string name = "";
	std::string parent_name = "";
	Eigen::Vector3d pos = Eigen::Vector3d::Zero();
	Eigen::Vector3d parent_pos = Eigen::Vector3d::Zero();
	Eigen::Vector3d physics_pos = Eigen::Vector3d::Zero();
	Eigen::Matrix3d basis;
	double density = 0;
	double xdim = 0;
	double ydim = 0;
	double zdim = 0;
	double angle = 0;
	double physics_density = 0;
	double physics_xdim = 0;
	double physics_ydim = 0;
	double physics_zdim = 0;
	double physics_angle = 0;
	Joint::JType joint = Joint::JType::Ball;
	double cons1 = 0;
	double cons2 = 0;
	bool invisible = false;
	bool physics = false;

	std::string type;
	while (buffer.peek() != EOF)
	{
		if (buffer.peek() == '#')
		{
			buffer.get();
			getline(buffer, type);
			if (type.compare("link") == 0)
			{
				name = "";
				parent_name = "";
				pos = Eigen::Vector3d::Zero();
				parent_pos = Eigen::Vector3d::Zero();
				physics_pos = Eigen::Vector3d::Zero();
				density = 0;
				xdim = 0;
				ydim = 0;
				zdim = 0;
				angle = 0;
				physics_density = 0;
				physics_xdim = 0;
				physics_ydim = 0;
				physics_zdim = 0;
				physics_angle = 0;
				joint = Joint::JType::Ball;
				cons1 = 0;
				cons2 = 0;
				invisible = false;
				physics = false;

				// Load a linkage
				std::string line;
				while (buffer.peek() != '#' && !buffer.eof())
				{
					getline(buffer, line);

					if (line.substr(0, 4).compare("name") == 0)
					{
						if (line.length() > 4)
						{
							name = line.substr(5, line.length());
						}
					}
					else if (line.substr(0, 4).compare("prnt") == 0)
					{
						if (line.length() > 4)
						{
							parent_name = line.substr(5, line.length());
						}
					}
					else if (line.substr(0, 4).compare("dens") == 0)
					{
						if (line.length() > 4)
						{
							std::stringstream ss;
							ss.str(line.substr(5, line.length()));
							if (!physics)
							{
								ss >> density;
							}
							else
							{
								ss >> physics_density;
							}
						}
					}
					else if (line.substr(0, 4).compare("size") == 0)
					{
						if (line.length() > 4)
						{
							std::stringstream ss;
							ss.str(line.substr(5, line.length()));
							if (!physics)
							{
								ss >> xdim;
								ss >> ydim;
								ss >> zdim;
							}
							else
							{
								ss >> physics_xdim;
								ss >> physics_ydim;
								ss >> physics_zdim;
							}
						}
					}
					else if (line.substr(0, 4).compare("angl") == 0)
					{
						if (line.length() > 4)
						{
							std::stringstream ss;
							ss.str(line.substr(5, line.length()));
							if (!physics)
							{
								ss >> angle;
							}
							else
							{
								ss >> physics_angle;
							}
						}
					}
					else if (line.substr(0, 4).compare("cpos") == 0)
					{
						if (line.length() > 4)
						{
							std::stringstream ss;
							ss.str(line.substr(5, line.length()));
							double x, y, z;
							ss >> x;
							ss >> y;
							ss >> z;
							if (!physics)
							{
								pos = Eigen::Vector3d(x, y, z);
							}
							else
							{
								physics_pos = Eigen::Vector3d(x, y, z);
							}
						}
					}
					else if (line.substr(0, 4).compare("ppos") == 0)
					{
						if (line.length() > 4)
						{
							std::stringstream ss;
							ss.str(line.substr(5, line.length()));
							double x, y, z;
							ss >> x;
							ss >> y;
							ss >> z;
							parent_pos = Eigen::Vector3d(x, y, z);
						}
					}
					else if (line.substr(0, 4).compare("type") == 0)
					{
						if (line.length() > 4)
						{
							line = line.substr(5, line.length());
							if (line.substr(0, 4).compare("ball") == 0)
							{
								joint = Joint::JType::Ball;
							}
							else if (line.substr(0, 6).compare("slider") == 0)
							{
								joint = Joint::JType::Slider;
							}
							else if (line.substr(0, 7).compare("springy") == 0)
							{
								joint = Joint::JType::Springy;
								std::stringstream ss;
								ss.str(line.substr(8, line.length()));
								ss >> cons1;
								ss >> cons2;
							}
							else if (line.substr(0, 7).compare("powered") == 0)
							{
								joint = Joint::JType::Powered;
								std::stringstream ss;
								ss.str(line.substr(8, line.length()));
								ss >> cons1;
							}
						}
					}
					else if (line.substr(0, 7).compare("physics") == 0)
					{
						physics = true;
					}
					else if (line.substr(0, 9).compare("invisible") == 0)
					{
						invisible = true;
					}
				}

				if (density <= 0 || xdim <= 0 || ydim <= 0 || zdim <= 0)
				{
					std::cout << "Link <" << name << "> is not realistic" << std::endl;
				}
				else
				{
					if (physics)
					{
						basis = Eigen::AngleAxisd((physics_angle * RB_M_PI), Eigen::Vector3d::UnitZ());
						disconnectInertial(LS, S, name, physics_pos, basis, physics_density, physics_xdim, physics_ydim, physics_zdim);
					}
				}
			}
		}
		else
		{
			char s[1024];
			buffer.getline(s, 1024);
		}
	}

	infile.close();
}

void RigidBodyCreator::disconnectInertial(
	std::unique_ptr<LinkageSystem>& LS, 
	std::unique_ptr<State>& S, 
	std::string name, 
	Eigen::Vector3d & p, 
	Eigen::Matrix3d & basis, 
	double density, double x, double y, double z)
{
	std::shared_ptr<Block> block = LS->find(name);
	int index = block->joint->jindex;
	Eigen::Matrix4d dToI = Eigen::Matrix4d::Identity();
	dToI.block<3, 3>(0, 0) = basis;
	dToI.block<3, 1>(0, 3) = p;

	if (density != 0)
	{
		block->indensity = density;
	}

	if (x != 0)
	{
		block->insize[0] = x;
	}
	if (y != 0)
	{
		block->insize[1] = y;
	}
	if (z != 0)
	{
		// Do not update based on z - instead change the density
		block->insize[2] = z;
	}

	// Update the inertia to display basis transform matrix
	block->hasInertial = true;
	block->iToD = dToI.inverse();
}

// Eliminate the need to 'resnap' the same setup every pass through the optimizer
void RigidBodyCreator::printLinkages(const std::unique_ptr<LinkageSystem>& LS, const std::unique_ptr<State>& S)
{
	std::ofstream file;
	file.open(RESOURCE_DIR + SNAPPED_LINKAGES);

	std::shared_ptr<Joint> parentJoint;
	std::shared_ptr <Block> block;
	Eigen::Vector3d cpos;
	Eigen::Vector3d ppos;
	Eigen::Vector3d sinthetaN;
	Eigen::AngleAxisd R;
	Eigen::Matrix4d bToJ0;
	double angle = 0;
	double phys_angle = 0;
	std::string type = "ball";
	for (int i = 0; i < LS->blocks.size(); ++i)
	{
		// Currently only works with blocks that have one or no (one root block) parents
		block = LS->blocks[i].second;

		bool hasParent = false;
		int parenti = -1;
		for (int j = 0; j < LS->joints.size(); ++j)
		{
			if (LS->joints[j]->childName.compare(LS->blocks[i].first) == 0)
			{
				hasParent = true;
				parentJoint = LS->joints[j];
				parenti = j;
			}
		}

		int index = block->joint->jindex;

		// move this into the individual joint derived structs
		std::shared_ptr<JointPowered> parentJointPow;
		std::shared_ptr<JointHinge> parentJointHinge;
		switch (parentJoint->type)
		{
		case Joint::JType::Ball:
			type = "ball";
			break;
		case Joint::JType::Springy:
			type = "springy";
			break;
		case Joint::JType::Slider:
			type = "slider";
			break;
		case Joint::JType::Powered:
			parentJointPow = std::dynamic_pointer_cast<JointPowered>(parentJoint);
			type = "powered " + std::to_string(parentJointPow->cpower[0]);
			break;
		case Joint::JType::Hinge:
			parentJointHinge = std::dynamic_pointer_cast<JointHinge>(parentJoint);
			type = "hinge " + std::to_string(parentJointHinge->axis[0]) + " " +
				std::to_string(parentJointHinge->axis[1]) + " " +
				std::to_string(parentJointHinge->axis[2]);
			break;
		default:
			type = "unknown";
			break;
		}

		if (hasParent)
		{
			R = parentJoint->E_ij.block<3, 3>(0, 0);
			double bb = R.angle() / RB_M_PI;

			R = S->E[index].block<3, 3>(0, 0);
			angle = bb * (parentJoint->E_ij(0, 1) < 0 ? -1 : 1);

			R = block->iToD.block<3, 3>(0, 0);
			phys_angle = R.angle() / RB_M_PI;
		}
		else if (!block->hasInertial)
		{
			R = S->E[index].block<3, 3>(0, 0);
			angle = R.angle() / RB_M_PI;
		}
		else
		{
			std::cout << "Flag this for later" << std::endl;
			R = S->E[index].block<3, 3>(0, 0);
			angle = R.angle() / RB_M_PI;

			R = block->iToD.block<3, 3>(0, 0);
			phys_angle = R.angle() / RB_M_PI;
		}

		file << "#link\n";
		file << "name " << LS->blocks[i].first << "\n";
		file << "dens " << block->density << "\n";
		file << "size " << block->size[0] << " " << block->size[1] << " " << block->size[2] << "\n";
		if (hasParent)
		{
			file << "prnt " << parentJoint->parentName << "\n";
			file << "cpos " << parentJoint->childPos[0] << " " << parentJoint->childPos[1] << " " << parentJoint->childPos[2] << "\n";
			file << "ppos " << parentJoint->parentPos[0] << " " << parentJoint->parentPos[1] << " " << parentJoint->parentPos[2] << "\n";
			file << "type " << type << "\n";
		}
		else
		{
			file << "cpos " << S->E[index](0, 3) << " " << S->E[index](1, 3) << " " << S->E[index](2, 3) << "\n";
		}

		if (angle != 0)
			file << "angl " << angle << "\n";

		if (block->invisible)
			file << "invisible\n";

		if (block->hasInertial)
		{
			file << "physics\n";
			file << "angl " << phys_angle << "\n";
			file << "cpos " << block->iToD(0, 3) << " " << block->iToD(1, 3) << " " << block->iToD(2, 3) << "\n";
		}
	}

	file.close();
}

void RigidBodyCreator::addDisplayBlock(
	std::unique_ptr<LinkageSystem>& LS, 
	std::unique_ptr<State>& S,
	std::string name,
	double density, double damping, double xwidth, double ywidth, double zwidth, 
	std::string parent_name, 
	std::shared_ptr<Joint>& j, 
	Eigen::Vector3d & pos, 
	Eigen::Matrix3d & basis, 
	Eigen::Vector3d & parent_jpos,
	Block::SType stype)
{
	bool root = false;
	if (parent_name.size() == 0)
		root = true;
	bool powered = false;
	if (j->type == Joint::JType::Powered)
		powered = true;

	// Load transformation matrix
	Eigen::Matrix4d E0_iToW = Eigen::Matrix4d::Identity();
	Eigen::Matrix4d E0_iTop = Eigen::Matrix4d::Identity();
	std::shared_ptr<Block> parent;
	if (root)
	{
		E0_iToW.block<3, 1>(0, 3) = parent_jpos - basis * pos;
		E0_iToW.block<3, 3>(0, 0) = basis;
		E0_iToW(3, 3) = 1;
		E0_iTop = E0_iToW;
	}
	else
	{
		parent = LS->find(parent_name);

		// Use the global location of the parent to determine the global position of the child
		E0_iToW.block<3, 1>(0, 3) = parent->E_wi0.block<3, 1>(0, 3) + parent->E_wi0.block<3, 3>(0, 0)*parent_jpos - (parent->E_wi0.block<3, 3>(0, 0)*basis) * pos;
		E0_iToW.block<3, 3>(0, 0) = parent->E_wi0.block<3, 3>(0, 0)*basis;
		E0_iToW(3, 3) = 1;
		E0_iTop.block<3, 1>(0, 3) = parent_jpos - basis * pos;
		E0_iTop.block<3, 3>(0, 0) = basis;
		E0_iTop(3, 3) = 1;
	}

	// Load spatial inertia
	Eigen::Vector3d s = Eigen::Vector3d(xwidth, ywidth, zwidth);
	Matrix6d M = Matrix6d::Zero();
	double mass;
	if (stype == Block::SType::Cuboid)
	{
		mass = xwidth*ywidth*zwidth*density;
		M(0, 0) = mass * (ywidth * ywidth + zwidth * zwidth) / 12.0;
		M(1, 1) = mass * (xwidth * xwidth + zwidth * zwidth) / 12.0;
		M(2, 2) = mass * (xwidth * xwidth + ywidth * ywidth) / 12.0;
		M(3, 3) = mass;
		M(4, 4) = mass;
		M(5, 5) = mass;
	}
	else
	{
		// cylinder
		double radius = xwidth;
		double height = ywidth;
		mass = density * RB_M_PI * radius * radius * height;

		M(0, 0) = mass * (3 * radius*radius + height * height) / 12.0;
		M(1, 1) = mass * (3 * radius*radius + height * height) / 12.0;
		M(2, 2) = mass * (radius*radius) * 0.5;
		M(3, 3) = mass;
		M(4, 4) = mass;
		M(5, 5) = mass;
	}


	// create and save the block
	LS->blocks.push_back(std::pair<std::string, std::shared_ptr<Block>>(name,
		std::make_shared < Block > (stype, LS->joints[LS->joints.size() - 1], s, root, powered, density, damping,
			mass, M, E0_iToW, E0_iTop)));
}
