#pragma once

#include "online/Brender/cpp/Brenderable.h"
#include "JSONwrapper.h"
#include "Constraint.h"
#include "Joint.h"
#include "Rigid.h"
#include "State.h"

#include <json/json.h>
#include <iostream>
#include <vector>
#include <memory>

static Eigen::Matrix4d idmat44 = Eigen::Matrix4d::Identity();

struct Block : public Brenderable
{
	enum SType { Cuboid, Cylinder};

	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	// Index of block physical data in State
	//std::string name;
	// shape of block
	SType shape;
	// pointer to its joint
	std::shared_ptr<Joint> joint;
	// a block knows all of its children
	std::vector<std::shared_ptr<Joint>> c_joints;
	// Store the size of the original block
	Eigen::Vector3d size;
	// Store the density of the original block
	double density;
	// Body damping
	double damping;
	Matrix6d D;

	// Store some extra block information (eliminate redundant on-the-go computations)
	Eigen::Matrix4d E_wi0;
	Eigen::Matrix4d E_wi0_inc;
	Eigen::Matrix4d E_pi0;

	// initial values
	Vector6d v0;
	Matrix6d M0;
	double mass0;

	// Toggle immovability of the block
	bool rooted;
	// Track whether or not this block is slave to a powered joint
	bool powered;
	
	// Toggle for disconnected display and inertial blocks
	bool hasInertial;
	// Transformation from the inertial block to the display block
	Eigen::Matrix4d iToD;
	// Toggle the drawing display block
	bool invisible;
	// Store the mass of the inertial block - defaults to the same as display mass
	double indensity;
	// Store the size of the inertial (physics) block - defaults to the same as display size
	Eigen::Vector3d insize;

	// local list of connected constraints
	//std::vector<std::shared_ptr<Constraint>> constraints;

	Block(SType st, std::shared_ptr<Joint> &j, Eigen::Vector3d &s, 
		bool rtd, bool pwd, double d, double da, double mass, 
		Matrix6d &M, Eigen::Matrix4d &E0_iToW_in, 
		Eigen::Matrix4d &E0_iTop_in, bool hi = false, 
		Eigen::Matrix4d &iToD_in = idmat44, bool inv = false) :
		shape(st),
		joint(j),
		size(s),
		density(d),
		damping(da),
		E_wi0(E0_iToW_in),
		E_pi0(E0_iTop_in),
		v0(Vector6d::Zero()),
		M0(M),
		mass0(mass),
		rooted(rtd),
		powered(pwd),
		hasInertial(hi),
		iToD(iToD_in),
		invisible(inv),
		indensity(d),
		insize(s)		
	{
		hasInertial = false;
		invisible = false;
		indensity = 0.0;
		iToD = Eigen::Matrix4d::Identity();
		insize = Eigen::Vector3d::Zero();
		D = damping*Matrix6d::Identity();
	};

	virtual void exportBrender(std::vector<std::shared_ptr<std::ofstream>> outfiles, int frame, double time) const
	{
		Eigen::Quaterniond quat(joint->E_draw.block<3, 3>(0, 0));
		//*(outfiles[0]) << "\"" << joint->childName << "_joint\":\n{\n"
		//	<< "\"scale\": [" << std::to_string(1) << "," << std::to_string(1) << ","
		//	<< std::to_string(1) << "],\n\"location\": [" << std::to_string(joint->E_draw(0, 3)) << ","
		//	<< std::to_string(joint->E_draw(1, 3)) << "," << std::to_string(joint->E_draw(2, 3))
		//	<< "],\n\"quat\": [" << std::to_string(quat.w()) << "," << std::to_string(quat.x()) << ","
		//	<< std::to_string(quat.y()) << "," << std::to_string(quat.z()) << "]\n}";

		Json::Value scale(Json::arrayValue);
		scale[0] = size[0];
		scale[1] = size[1];
		scale[2] = size[2];

		Json::Value location(Json::arrayValue);
		location[0] = joint->E_draw(0, 3);
		location[1] = joint->E_draw(1, 3);
		location[2] = joint->E_draw(2, 3);

		Json::Value q(Json::arrayValue);
		q[0] = quat.w();
		q[1] = quat.x();
		q[2] = quat.y();
		q[3] = quat.z();

		Json::Value v;
		v["scale"] = scale;
		v["location"] = location;
		v["quat"] = q;

		*(outfiles[0]) << v;
	}

	Json::Value exportJson(const std::unique_ptr<State> & S, const std::unique_ptr<JSONwrapper> &J)
	{
		// I used these for th umbrella render
		//Eigen::Matrix4d size4 = Eigen::Matrix4d::Identity();

		//Eigen::Matrix4d E = S->E[joint->jindex] * Eigen::Affine3d(Eigen::AngleAxisd(1.5707963, Eigen::Vector3d(1,0,0))).matrix();

		//// The display and the inertial blocks are the same
		//size4(0, 0) = size[0] / 2;
		//size4(1, 1) = size[0] / 2;
		//size4(2, 2) = 0.05;// size[1] / 2;
		//size4 = E * size4;

		//Eigen::Quaterniond quat(E.block<3, 3>(0, 0));

		//Json::Value scale(Json::arrayValue);
		//scale[0] = size[0]/2;
		//scale[1] = size[1]/2;
		//scale[2] = 0.05;// size[2] / 2;

		//Json::Value location(Json::arrayValue);
		//location[0] = E(0, 3);
		//location[1] = E(1, 3);
		//location[2] = E(2, 3);

		//J->block_scales.push_back(Json::Value(Json::arrayValue));
		//J->block_scales[J->block_scales.size() - 1][0] = size[0];
		//J->block_scales[J->block_scales.size() - 1][1] = size[1];
		//J->block_scales[J->block_scales.size() - 1][2] = size[2];

		//J->block_locations.push_back(Json::Value(Json::arrayValue));
		//J->block_locations[J->block_locations.size() - 1][0] = joint->E_draw(0, 3);
		//J->block_locations[J->block_locations.size() - 1][1] = joint->E_draw(1, 3);
		//J->block_locations[J->block_locations.size() - 1][2] = joint->E_draw(2, 3);

		//J->block_quats.push_back(Json::Value(Json::arrayValue));
		//J->block_quats[J->block_quats.size() - 1][0] = quat.x();
		//J->block_quats[J->block_quats.size() - 1][1] = quat.y();
		//J->block_quats[J->block_quats.size() - 1][2] = quat.z();
		//J->block_quats[J->block_quats.size() - 1][3] = quat.w();

		//J->blocks.push_back(Json::Value());
		//J->blocks[J->blocks.size() - 1]["scale"] = J->block_scales[J->block_scales.size() - 1];
		//J->blocks[J->blocks.size() - 1]["location"] = J->block_locations[J->block_locations.size() - 1];
		//J->blocks[J->blocks.size() - 1]["quat"] = J->block_quats[J->block_quats.size() - 1];

		//std::cout << J->blocks[J->blocks.size() - 1] << std::endl;

		Eigen::Quaterniond quat(joint->E_draw.block<3, 3>(0, 0));

		Json::Value scale(Json::arrayValue);
		scale[0] = size[0]/2.0;
		scale[1] = size[1]/2.0;
		scale[2] = size[2]/2.0;
		//scale[0] = size[0]/2.0;
		//scale[1] = size[0]/2.0;
		//scale[2] = size[1]/2.0;

		Json::Value location(Json::arrayValue);
		location[0] = joint->E_draw(0, 3);
		location[1] = joint->E_draw(1, 3);
		location[2] = joint->E_draw(2, 3);

		Json::Value q(Json::arrayValue);
		q[0] = quat.x();
		q[1] = quat.y();
		q[2] = quat.z();
		q[3] = quat.w();

		Json::Value v;
		v["scale"] = scale;
		v["location"] = location;
		v["quat"] = q;

		return v;
	}

	void updateSpatialIntertia(const std::unique_ptr<State> &S)
	{
		int index = joint->jindex;

		double xdim = 0;
		double ydim = 0;
		double zdim = 0;

		if (shape == SType::Cuboid)
		{
			if (hasInertial)
			{
				S->mass[index] = insize[0] * insize[1] * insize[2] * indensity;
				xdim = insize[0];
				ydim = insize[1];
				zdim = insize[2];
			}
			else
			{
				S->mass[index] = size[0] * size[1] * size[2] * density;
				xdim = size[0];
				ydim = size[1];
				zdim = size[2];
			}

			S->M[index](0, 0) = S->mass[index] * (ydim * ydim + zdim * zdim) / 12.0;
			S->M[index](1, 1) = S->mass[index] * (xdim * xdim + zdim * zdim) / 12.0;
			S->M[index](2, 2) = S->mass[index] * (xdim * xdim + ydim * ydim) / 12.0;
			S->M[index](3, 3) = S->mass[index];
			S->M[index](4, 4) = S->mass[index];
			S->M[index](5, 5) = S->mass[index];
		}
		else
		{
			// cylinder
			double radius = size[0];
			double height = size[1];
			S->mass[index] = radius*radius * height * density * RB_M_PI;
			
			S->M[index](0, 0) = S->mass[index] * (3 * radius*radius + height * height) / 12.0;
			S->M[index](1, 1) = S->mass[index] * (3 * radius*radius + height * height) / 12.0;
			S->M[index](2, 2) = S->mass[index] * (radius*radius) * 0.5;
			S->M[index](3, 3) = S->mass[index];
			S->M[index](4, 4) = S->mass[index];
			S->M[index](5, 5) = S->mass[index];
		}
	}
};