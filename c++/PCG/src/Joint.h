#pragma once

#include "Rigid.h"
#include "online/Brender/cpp/Brenderable.h"

#include <json/json.h>
#include <memory>
#include <iostream>

struct Block;
struct LinkageSystem;
struct State;
struct StateDeriv;
struct StateSolve;

struct Joint : public Brenderable
{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	enum JType { Hinge, Fixed, Prismatic, Ball, Springy, Slider, Powered, Universal };

	// Joint type
	JType type;
	// Pointer to child block
	std::shared_ptr<Block> parentBlock;
	std::shared_ptr<Block> childBlock;
	// Identity of the two blocks
	std::string parentName;
	std::string childName;
	// Where on the object it connects
	Eigen::Vector3d parentPos;
	Eigen::Vector3d childPos;
	// If it is a grounded (immovable) pin
	bool root;
	// stiffness and damping coefficients
	double k;
	double d;
	// list index
	int index0;
	int jindex;
	int constraint_index;
	int constraintNum;

	// State data
	Eigen::Map<Eigen::VectorXd> qdot;
	Eigen::Map<Eigen::VectorXd> q;

	// maximal joint velocity
	Vector6d V;

	// for reduced coordinates
	Eigen::MatrixXd S;
	Eigen::MatrixXd ST;
	Eigen::MatrixXd Sdot;
	Eigen::VectorXd q0;
	Eigen::VectorXd qdot0;
	Eigen::Matrix4d Q;

	// inertia at joint
	Matrix6d I_j;	// if M_bodycenter is stored in State, this is M_jointcenter

	/// Initial transforms
	// from the parent joint to the current joint - reduced coords
	Eigen::Matrix4d E_pj_ij;
	Matrix6d ad_pj_ij;
	Matrix6d ad_ij_pj;
	// from the parent block to the current joint - maximal
	bool E_pj0_;
	Eigen::Matrix4d E_pj0;
	Eigen::Matrix4d E_jp0;
	Matrix6d ad_jp0;
	// constant joint-body transform
	bool E_ji0_;
	Eigen::Matrix4d E_ji0;
	Eigen::Matrix4d E_ij;
	Eigen::Matrix4d E_ji;
	Matrix6d ad_ij;
	Matrix6d ad_ijT;
	Matrix6d ad_ji;
	/// end - Initial transforms
	// other useful transforms - save from update step
	Eigen::Matrix4d E_pj;
	Matrix6d ad_jp;
	Matrix6d ad_jpT;
	Eigen::Matrix4d E_wj;
	Eigen::Matrix4d E_iw;
	// set to 0 after each step finalization - save computation!
	bool addot_;
	Matrix6d addot;
	bool ad_ip_;
	Matrix6d ad_ip;
	Matrix6d ad_ipT;
	bool add_ip_;
	Matrix6d add_ip;
	bool ad_iw_;
	Matrix6d ad_iw;
	bool ad_wi_;
	Matrix6d ad_wi;
	// variables for no-matrix solve - three JT solve
	Vector6d alpha;/// JT
	Vector6d alpha1;
	Vector6d alpha2;
	Vector6d alpha3;
	// for no-matrix preconditioning
	Eigen::MatrixXd Pr;
	Matrix6d Kmd;
	Matrix6d Dmd;
	Matrix6d Mhat;
	Eigen::MatrixXd Psi;
	Matrix6d Pi;
	Vector6d Bhat;/// minv
	Vector6d beta;/// minv
	Vector6d Vdot;/// minv
	// save computations that are constant within a step
	Eigen::MatrixXd adij_S;
	Eigen::MatrixXd ST_adijT;
	// for minv only - can be handled in preprocess preconditioner
	Eigen::MatrixXd ST_Mhat;
	Eigen::MatrixXd S_Psi;
	Eigen::MatrixXd Mhat_S_Psi;
	Eigen::MatrixXd ST_Bhat;/// minv

	// for draw
	Eigen::Matrix4d E_draw;

	Joint() :
		qdot(NULL, 0),
		q(NULL, 0) 
	{};
	Joint(std::unique_ptr<State> &S,
		JType type,
		std::string parentName,
		std::string childName,
		Eigen::Vector3d parentPos,
		Eigen::Vector3d childPos,
		bool root,
		int index,
		int constraintNum,
		double k = 0,
		double d = 0) :
		type(type),
		parentName(parentName),
		childName(childName),
		parentPos(parentPos),
		childPos(childPos),
		root(root),
		jindex(index),
		index0(index),
		constraint_index(index),
		constraintNum(constraintNum),
		k(k),
		d(d),
		qdot(NULL, 0),
		q(NULL, 0)
	{
		V = Vector6d::Zero();
		q0 = Eigen::VectorXd(constraintNum);
		qdot0 = Eigen::VectorXd::Zero(constraintNum);
		E_pj0_ = false;
		E_ji0_ = false;
		//this->k = 0;
		//this->d = 0;
		//E_pj0 = Eigen::Matrix4d::Identity();
	};

	virtual void updateMaxSparse(std::unique_ptr<StateSolve> &SS, const std::unique_ptr<LinkageSystem> &LS, const std::unique_ptr<State> &S, Matrix6d &adp, Matrix6d &adc, Eigen::Vector3d &g) const = 0;
 
	void DFSIndexSet(int &blockID, int &jointID, 
		std::vector< std::pair < std::string, std::shared_ptr<Block> > > & blocks, 
		std::vector<int>& j_map);
	void initState(int constraint_index, const std::unique_ptr<LinkageSystem>& LS, std::unique_ptr<State> &S);
	void initState_SoA(int constraint_index, const std::unique_ptr<LinkageSystem>& LS, std::unique_ptr<State> &S);
	void update(const std::unique_ptr<LinkageSystem> &LS, std::unique_ptr<State> &S);
	void update_SoA(const std::unique_ptr<LinkageSystem> &LS, std::unique_ptr<State> &S);
	void update_nostatechange(const std::unique_ptr<LinkageSystem> &LS);

	Json::Value exportJson(const std::unique_ptr<LinkageSystem> & LS, const std::unique_ptr<State> & S);
	void exportBrender(std::vector< std::shared_ptr< std::ofstream > > outfiles, int frame, double time) const;

protected:
	virtual void mapSelf(int c_i, std::unique_ptr<State> &S) = 0;
	virtual void updateSelf() = 0;
	virtual Vector6d Sqdot(const std::unique_ptr<State> &S) = 0;
	
	void compute_Ij(std::unique_ptr<State>& S);
};