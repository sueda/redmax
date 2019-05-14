#pragma once

#include <memory>

#include "RigidBodyUtility.h"
#include "Block.h"

struct Joint;
struct LinkageSystem;
struct State;
struct StateSolve;

class RigidBodyCreator
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    RigidBodyCreator(double t):
            t(t)
    {};
	~RigidBodyCreator() {};

	enum dynamicType { chain, tree };

	void loadTree(
		int n,
		simType s,
		std::unique_ptr<LinkageSystem>& LS,
		std::unique_ptr<State>& S,
		std::unique_ptr<StateSolve>& SS);
	void loadBridge(
		int nbridge,
		int ntower,
		simType s,
		std::unique_ptr<LinkageSystem>& LS,
		std::unique_ptr<State>& S,
		std::unique_ptr<StateSolve>& SS);
	void loadSimpleBridge(
		int nbridge,
		simType s,
		std::unique_ptr<LinkageSystem>& LS,
		std::unique_ptr<State>& S,
		std::unique_ptr<StateSolve>& SS);
	void loadSimpleTree(
		int n,
		simType s,
		std::unique_ptr<LinkageSystem>& LS,
		std::unique_ptr<State>& S,
		std::unique_ptr<StateSolve>& SS);
	void loadUmbrella(
		int n,
		simType s,
		std::unique_ptr<LinkageSystem>& LS,
		std::unique_ptr<State>& S,
		std::unique_ptr<StateSolve>& SS);
	void loadTest(
		int n,
		simType s,
		std::unique_ptr<LinkageSystem>& LS,
		std::unique_ptr<State>& S,
		std::unique_ptr<StateSolve>& SS);
	void loadChain(
		int n,
		simType s,
		std::unique_ptr<LinkageSystem>& LS,
		std::unique_ptr<State>& S,
		std::unique_ptr<StateSolve>& SS);
	void loadLinkagesfromFile(
		std::unique_ptr<LinkageSystem>& LS,
		std::unique_ptr<State>& S,
		std::unique_ptr<StateSolve>& SS,
		const bool disconnect);
	void loadPhysicsfromFile(
		std::unique_ptr<LinkageSystem>& LS,
		std::unique_ptr<State>& S);

	void disconnectInertial(
		std::unique_ptr<LinkageSystem>& LS,
		std::unique_ptr<State>& S,
		std::string name,
		Eigen::Vector3d &p,
		Eigen::Matrix3d &basis,
		double density = 0.0,
		double x = 0.0,
		double y = 0.0,
		double z = 0.0);

	void setResourceDir(const std::string &RESOURCE_DIR) { this->RESOURCE_DIR = RESOURCE_DIR; };
	void setUserSource(const std::string &FILENAME) { this->USER_SOURCE = FILENAME; };
	void useSourceLinks() { activeLinks = USER_SOURCE; };
	void useSnappedLinks() { activeLinks = SNAPPED_LINKAGES; };
	void useSourceSplines() { activeSplines = USER_SOURCE; };
	void useUserSpline() { activeSplines = USER_SPLINE; };
	void useFitSpline() { activeSplines = FIT_SPLINE; };
	void useResultSpline() { activeSplines = SPLINE_RESULT; };

	std::string getResourceDir() { return RESOURCE_DIR; };

	void printLinkages(
		const std::unique_ptr<LinkageSystem>& LS, 
		const std::unique_ptr<State>& S);
private:
    double &t;		// reference 'global' sim time - from rigid body main

	const std::string SNAPPED_LINKAGES = "linkages_auto_save.txt";
	const std::string USER_SPLINE = "splines_auto_save.txt";
	const std::string FIT_SPLINE = "fit_splines_auto_save.txt";
	const std::string SPLINE_RESULT = "result_spline_auto_save.txt";
	std::string RESOURCE_DIR;
	std::string USER_SOURCE;
	std::string activeLinks;
	std::string activeSplines;

	void addDisplayBlock(
		std::unique_ptr<LinkageSystem> & LS,
		std::unique_ptr<State>& S,
		std::string name,
		double density,
		double damping,
		double xwidth,
		double ywidth,
		double zwidth,
		std::string parent_name,
		std::shared_ptr<Joint> &j,
		Eigen::Vector3d &pos,
		Eigen::Matrix3d &basis,
		Eigen::Vector3d &parent_jpos = emptyVector3d,
		Block::SType stype = Block::SType::Cuboid
	);
};
