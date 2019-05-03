#pragma once

#include <vector>
#include <map>
#include <memory>
#include <string>

class Constraint;
struct Block;
struct Joint;

// ----------------------------------------------------//
// A collection of blocks(links, rigid bodies) that    //
// are connected by joints to make a linkage system    //
// ----------------------------------------------------//
struct LinkageSystem
{
	// the links, blocks, rigid_bodies
	std::vector< std::pair < std::string, std::shared_ptr<Block> > > blocks;
	// the joints between objects in the linkage system
	std::vector<std::shared_ptr<Joint>> joints;
	// the features that apply constraints to the rigid bodies
	std::vector<std::shared_ptr<Constraint>> constraints;

	// reduced coordinate ordering
	std::vector<int> joint_map;

	std::shared_ptr<Block> find(std::string name)
	{
		for (int i = 0; i < blocks.size(); ++i)
		{
			if (blocks[i].first.compare(name) == 0)
				return blocks[i].second;
		}

		return nullptr;
	}
};