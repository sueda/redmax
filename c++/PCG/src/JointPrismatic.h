#pragma once

#include "Joint.h"

struct LinkageSystem;
struct State;
struct StateDeriv;
struct StateSolve;

struct JointPrismatic : public Joint
{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	// axis of rotation
	Eigen::Vector3d axis;

	JointPrismatic(std::unique_ptr<State> &S,
		std::string pn,
		std::string cn,
		Eigen::Vector3d pp,
		Eigen::Vector3d cp,
		bool r,
		int i,
		double angle,
		Eigen::Vector3d &axis,
		double k = 0,
		double d = 0);

	void updateMaxSparse(std::unique_ptr<StateSolve> &SS, const std::unique_ptr<LinkageSystem> &LS, const std::unique_ptr<State> &S, Matrix6d &parentAdj, Matrix6d &childAdj, Eigen::Vector3d &g) const;

private:
	void mapSelf(int c_i, std::unique_ptr<State> &S);
	void updateSelf();
	Vector6d Sqdot(const std::unique_ptr<State> &S);
};