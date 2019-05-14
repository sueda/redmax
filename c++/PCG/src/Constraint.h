#pragma once

#include "Rigid.h"
#include "RigidBodyUtility.h"
#include "online/Brender/cpp/Brenderable.h"
#include <json/writer.h>
#include <json/value.h>
#include <json/json.h>

#include <memory>

struct Block;
struct Joint;
struct LinkageSystem;
struct State;
struct StateDeriv;
struct StateSolve;

class Constraint
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	enum constraintType { elastic, closehinge, springpoint };

	Constraint();
	Constraint(constraintType type);
	~Constraint() {};

	virtual void initJoint(const std::unique_ptr<StateSolve> &SS, const std::unique_ptr<LinkageSystem> &LS, const std::unique_ptr<State> &S) = 0;

	virtual void computeCGProd(Eigen::VectorXd &M_x, const Eigen::VectorXd &x, const std::unique_ptr<StateSolve>& SS, const std::unique_ptr<LinkageSystem>& LS, const std::unique_ptr<State>& S) = 0;

    // TODO: changed to non virtual 
    void update(std::unique_ptr<StateSolve> &SS, const std::unique_ptr<LinkageSystem> &LS, const std::unique_ptr<State> &S) {}
    
	virtual void update(std::unique_ptr<StateSolve> &SS, const std::unique_ptr<LinkageSystem> &LS, const std::unique_ptr<State> &S, std::unique_ptr<StateDeriv> &DS) = 0;
	virtual void updateJoint(std::unique_ptr<StateSolve> &SS, const std::unique_ptr<LinkageSystem> &LS, const std::unique_ptr<State> &S) = 0;
	virtual void draw(const std::unique_ptr<LinkageSystem> &LS, const std::unique_ptr<State> &S) const = 0;

	virtual Json::Value exportBrender() = 0;
	virtual void exportJson(std::ofstream & outfile, const std::unique_ptr<State> &S) = 0;

	// constraint type
	constraintType type;
	int constraintNum;
	int constraintNumMax;
	// Index for the G matrix
	int constraint_index;
	int constraint_index_m;
};

class Elastic : public Constraint
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	Elastic();
	Elastic(std::shared_ptr<Block> blocka,
		std::shared_ptr<Block> blockb,
		Eigen::Vector3d pa,
		Eigen::Vector3d pb,
		double rest,
		double stiffness,
		double damping);
	~Elastic() {};

	void initJoint(const std::unique_ptr<StateSolve> &SS, const std::unique_ptr<LinkageSystem> &LS, const std::unique_ptr<State> &S);

	void computeCGProd(Eigen::VectorXd &M_x, const Eigen::VectorXd &x, const std::unique_ptr<StateSolve>& SS, const std::unique_ptr<LinkageSystem>& LS, const std::unique_ptr<State>& S);

	void update(std::unique_ptr<StateSolve> &SS, const std::unique_ptr<LinkageSystem> &LS, const std::unique_ptr<State> &S);
	void update(std::unique_ptr<StateSolve> &SS, const std::unique_ptr<LinkageSystem> &LS, const std::unique_ptr<State> &S, std::unique_ptr<StateDeriv> &DS);
	void updateJoint(std::unique_ptr<StateSolve> &SS, const std::unique_ptr<LinkageSystem> &LS, const std::unique_ptr<State> &S);
	void draw(const std::unique_ptr<LinkageSystem> &LS, const std::unique_ptr<State> &S) const;

	Json::Value exportBrender();
	void exportJson(std::ofstream & outfile, const std::unique_ptr<State> &S);

	void initLengthToRest(const std::unique_ptr<LinkageSystem> &LS, const std::unique_ptr<State> &S);

	// parent block a
	std::shared_ptr<Block> blocka;
	// parent block b
	std::shared_ptr<Block> blockb;
	// Position of string on parent block a
	Eigen::Vector3d posa;
	// Position of string on parent block b
	Eigen::Vector3d posb;
	// The rest length of the string
	double L;
	// The stiffness or strength of the string
	double k;
	// Damping factor
	double d;

	Eigen::Matrix<double, 12, 12> K;
	Eigen::Matrix<double, 12, 12> D;
};

class CloseHinge : public Constraint
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	CloseHinge();
	CloseHinge(std::shared_ptr<Block> blocka,
		std::shared_ptr<Block> blockb,
		Eigen::Vector3d pa,
		Eigen::Vector3d pb,
		Eigen::Vector3d axis,
		double rest,
		double stiffness,
		double damping);
	~CloseHinge() {};

	void initJoint(const std::unique_ptr<StateSolve> &SS, const std::unique_ptr<LinkageSystem> &LS, const std::unique_ptr<State> &S);

	void computeCGProd(Eigen::VectorXd &M_x, const Eigen::VectorXd &x, const std::unique_ptr<StateSolve>& SS, const std::unique_ptr<LinkageSystem>& LS, const std::unique_ptr<State>& S);

	void update(std::unique_ptr<StateSolve> &SS, const std::unique_ptr<LinkageSystem> &LS, const std::unique_ptr<State> &S);
	void update(std::unique_ptr<StateSolve> &SS, const std::unique_ptr<LinkageSystem> &LS, const std::unique_ptr<State> &S, std::unique_ptr<StateDeriv> &DS);
	void updateJoint(std::unique_ptr<StateSolve> &SS, const std::unique_ptr<LinkageSystem> &LS, const std::unique_ptr<State> &S);
	void draw(const std::unique_ptr<LinkageSystem> &LS, const std::unique_ptr<State> &S) const;

	Json::Value exportBrender();
	void exportJson(std::ofstream & outfile, const std::unique_ptr<State> &S);

	void initLengthToRest(const std::unique_ptr<LinkageSystem> &LS, const std::unique_ptr<State> &S);

	// parent block a
	std::shared_ptr<Block> blocka;
	// parent block b
	std::shared_ptr<Block> blockb;
	// Position of string on parent block a
	Eigen::Vector3d posa;
	// Position of string on parent block b
	Eigen::Vector3d posb;
	// Axis of hinge
	Eigen::Vector3d axis;
	// Orthonormal vectors
	Eigen::Vector3d v1;
	Eigen::Vector3d v2;
	// The rest angle
	double q0;
	// The stiffness
	double k;
	// Damping factor
	double d;

	Eigen::MatrixXd Gm;

};

class CloseUniversal : public Constraint
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	CloseUniversal();
	CloseUniversal(std::shared_ptr<Block> blocka,
		std::shared_ptr<Block> blockb,
		Eigen::Vector3d pa,
		Eigen::Vector3d pb);
	~CloseUniversal() {};

	void initJoint(const std::unique_ptr<StateSolve> &SS, const std::unique_ptr<LinkageSystem> &LS, const std::unique_ptr<State> &S);

	void computeCGProd(Eigen::VectorXd &M_x, const Eigen::VectorXd &x, const std::unique_ptr<StateSolve>& SS, const std::unique_ptr<LinkageSystem>& LS, const std::unique_ptr<State>& S);

	void update(std::unique_ptr<StateSolve> &SS, const std::unique_ptr<LinkageSystem> &LS, const std::unique_ptr<State> &S);
	void update(std::unique_ptr<StateSolve> &SS, const std::unique_ptr<LinkageSystem> &LS, const std::unique_ptr<State> &S, std::unique_ptr<StateDeriv> &DS);
	void updateJoint(std::unique_ptr<StateSolve> &SS, const std::unique_ptr<LinkageSystem> &LS, const std::unique_ptr<State> &S);
	void draw(const std::unique_ptr<LinkageSystem> &LS, const std::unique_ptr<State> &S) const;

	Json::Value exportBrender();
	void exportJson(std::ofstream & outfile, const std::unique_ptr<State> &S);

	void initLengthToRest(const std::unique_ptr<LinkageSystem> &LS, const std::unique_ptr<State> &S);

	// parent block a
	std::shared_ptr<Block> blocka;
	// parent block b
	std::shared_ptr<Block> blockb;
	// Position of string on parent block a
	Eigen::Vector3d posa;
	// Position of string on parent block b
	Eigen::Vector3d posb;

	// Precomputed
	Matrix3x6d gamma_a;
	Matrix3x6d gamma_b;

	Eigen::MatrixXd Gm;

};

class SpringPoint : public Constraint
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	SpringPoint();
	SpringPoint(std::shared_ptr<Block> blocka,
		Eigen::Vector3d pa,
		Eigen::Vector3d dir,
		double stiffness);
	~SpringPoint() {};

	void initJoint(const std::unique_ptr<StateSolve> &SS, const std::unique_ptr<LinkageSystem> &LS, const std::unique_ptr<State> &S);

	void computeCGProd(Eigen::VectorXd &M_x, const Eigen::VectorXd &x, const std::unique_ptr<StateSolve>& SS, const std::unique_ptr<LinkageSystem>& LS, const std::unique_ptr<State>& S);

	void update(std::unique_ptr<StateSolve> &SS, const std::unique_ptr<LinkageSystem> &LS, const std::unique_ptr<State> &S);
	void update(std::unique_ptr<StateSolve> &SS, const std::unique_ptr<LinkageSystem> &LS, const std::unique_ptr<State> &S, std::unique_ptr<StateDeriv> &DS);
	void updateJoint(std::unique_ptr<StateSolve> &SS, const std::unique_ptr<LinkageSystem> &LS, const std::unique_ptr<State> &S);
	void draw(const std::unique_ptr<LinkageSystem> &LS, const std::unique_ptr<State> &S) const;

	Json::Value exportBrender();
	void exportJson(std::ofstream & outfile, const std::unique_ptr<State> &S);

	// parent block a
	std::shared_ptr<Block> blocka;
	// Position of string on parent block a
	Eigen::Vector3d posa;
	// Direction of the point force
	Eigen::Vector3d dir;
	// The stiffness or strength of the string
	double k;

	Eigen::Matrix<double, 6, 6> K;
};
