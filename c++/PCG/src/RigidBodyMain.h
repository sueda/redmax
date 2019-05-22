#pragma once

#ifndef EIGEN_USE_MKL
#define EIGEN_USE_MKL
#endif
#ifndef EIGEN_DONT_ALIGN_STATICALLY
#define EIGEN_DONT_ALIGN_STATICALLY
#endif
//#ifndef EIGEN_DONT_ALIGN
//#define EIGEN_DONT_ALIGN
//#endif
//#ifndef EIGEN_NO_STATIC_ASSERT
//#define EIGEN_NO_STATIC_ASSERT
//#endif
//#ifndef EIGEN_DISABLE_UNALIGNED_ARRAY_ASSERT
//#define EIGEN_DISABLE_UNALIGNED_ARRAY_ASSERT
//#endif
//#ifndef EIGEN_DONT_VECTORIZE 
//#define EIGEN_DONT_VECTORIZE 
//#endif
#include <Eigen/SparseCholesky>	
#include <Eigen/OrderingMethods>
#include <Eigen/Core>
#include <Eigen/Geometry>

#include "online/Brender/cpp/Brenderable.h"

/** Rigid Body System Files **/
#include "Block.h"
#include "Constraint.h"
#include "ConstraintJoint.h"
#include "Joint.h"
#include "JointBall.h"
#include "JointFixed.h"
#include "JointHinge.h"
#include "JointPowered.h"
#include "JointPrismatic.h"
#include "JointSlider.h"
#include "JointSpringy.h"
#include "JointUniversal.h"
#include "JSONwrapper.h"
#include "LinkageSystem.h"
#include "Rigid.h"
#include "RigidBodyCreator.h"
#include "RigidBodyUtility.h"
#include "Solver.h"
#include "State.h"
#include "StateDeriv.h"
#include "StateSolve.h"
/** Rigid Body System Files **/

#include <vector>
#include <set>
#include <memory>
#include <string>

class Cuboid;
class Program;
class MatrixStack;

//class BrenderWrapper : public Brenderable
//{
//public:
//	BrenderWrapper() :
//		Brenderable()
//	{};
//	~BrenderWrapper() {};
//	
//	virtual void exportBrender(std::vector<std::shared_ptr<std::ofstream>> outfiles, int frame, double time) const;
//
//	void linkRigidBodyMain(std::shared_ptr<RigidBodyMain> &rbm) { this->rbm = rbm; };
//
//	// brender variables
//	std::shared_ptr<RigidBodyMain> rbm;
//	int export_part;
//};
#ifdef REDMAX_JSONCPP
class RigidBodyMain : public Brenderable
#else
class RigidBodyMain
#endif
{
	friend class RigidBodyOptFunc;
//	friend class BrenderWrapper;

public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	enum userAction { selectSimType, editSpline, convergeLS, singleStep, free };
	std::set<userAction> blockingAction { selectSimType, editSpline };
	enum clickAction { press, release };

	RigidBodyMain(const std::string & RESOURCE_DIR, double &t);
	~RigidBodyMain();

	// load from file
	void load(const std::string &FILENAME);
	// create a dynamic linkage system
	void loadTree(int n);
	void loadSimpleTree(int n);
	void loadBridge(int nbridge, int ntower);
	void loadSimpleBridge(int nbridge);
	void loadUmbrella(int n);
	void loadChain(int n);
	void loadTest(int n);
	void reLoad(const bool disconnect = true);
	void updateInitialScene(std::vector<double> &x);

	void init();
	void reset();

	Eigen::VectorXd step();

	void update_user_forces();
	
	void draw(std::unique_ptr<MatrixStack> &MV, const std::unique_ptr<Program> &prog) const;
	void drawJoints(std::unique_ptr<MatrixStack> &MV) const;
	void drawPoints(std::unique_ptr<MatrixStack> &MV, const std::unique_ptr<Program> &prog) const;

	void displayUserActions();
	int handleKeyPress(unsigned int key);
	bool userAskForAction(userAction a);
	bool isUserActionOverride();
	int handleOverrideInput(unsigned int key);

	double getCertNormSum();
	Eigen::Vector3d getLeafCertificate();

	double getTime() const { return t; }
	int getNumLinks() const { return (int)LS->blocks.size(); };
	int getRedmaxVsizes() const { return ConstraintJoint::getConstraintNumReduced(LS); };
	int getRedmaxDOFs() const { return ConstraintJoint::getConstraintNumReduced(LS) + ConstraintJoint::getConstraintNumAdditnl(LS); };
	int getMaximalDOFs() const { return ConstraintJoint::getConstraintNumMaximal(LS) + ConstraintJoint::getConstraintNumAdditnlMaximal(LS); };
	double getSolveTime() const { return solver->trackLastTimestep.time_in_solve; };
	double getCGIterations() const { return solver->trackLastTimestep.num_iterations; };
	void setSolveType(simType s) { simtype = s; };
	void toggleBlockDisplayMode() /*{ showInertial = !showInertial; }*/;
#ifdef REDMAX_JSONCPP
	virtual void exportBrender(std::vector<std::shared_ptr<std::ofstream>> outfiles, int frame, double time) const;
	virtual std::vector<std::string> getBrenderExtensions() const;
	virtual std::vector<std::string> getBrenderNames() const;
	virtual std::vector<int> getBrenderTypes() const;
	virtual int getBrenderCount() const;
#endif
	int export_part;

private:
	bool askForAction(userAction a);
	bool handleSelectSimType(unsigned int key);

	// User interaction
	userAction userActionState;

	// All rigid bodies are drawn with a single shape representation
	std::unique_ptr<Cuboid> cuboid;

	// Drawing and editing variables
	bool verbose = false;
	bool showInertial = false;
	bool fullyConverged = false;
	bool showAllSplines = false;
	int editTrackingParticleID = -1;

	// certificate index
	int leaf_cert_index;

	// Environment variables
	double &t;

	std::unique_ptr<RigidBodyCreator> creator;
	std::unique_ptr<Solver> solver;

	std::unique_ptr<State> S;
	std::unique_ptr<StateDeriv> DS;
	std::unique_ptr<StateSolve> SS;
	std::unique_ptr<LinkageSystem> LS;
#ifdef REDMAX_JSONCPP
	std::unique_ptr<JSONwrapper> J;
#endif
	// used to track sim step in optimization routines
	// used when validating partial derivatives ?? - this may not be nec.
	int ct = 0;
	// used to validate partial derivatives and gradient
	int ctMAX = 11;
	Eigen::VectorXd rand_init_vel;
};
