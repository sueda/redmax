#pragma once

#include "Joint.h"

struct LinkageSystem;
struct State;
struct StateDeriv;
struct StateSolve;

struct JointPowered : public Joint
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	// Power constant
	Eigen::Vector3d cpower;
    JointPowered(){}
	JointPowered(std::unique_ptr<State> &S,
		std::string pn,
		std::string cn,
		Eigen::Vector3d pp,
		Eigen::Vector3d cp,
		bool r,
		int i,
		Eigen::Vector3d &cpower) :
	Joint(S, JType::Powered, pn, cn, pp, cp, r, i, 1)
	{
		this->cpower = cpower;
	};
    void updateMaxSparse(std::unique_ptr<StateSolve> &SS, const std::unique_ptr<LinkageSystem> &LS, const std::unique_ptr<State> &S, Matrix6d &parentAdj, Matrix6d &childAdj, Eigen::Vector3d &g) const;
	
private:
	void mapSelf(int c_i, std::unique_ptr<State> &S);
	void updateSelf();
	Vector6d Sqdot(const std::unique_ptr<State> &S)
	{
		return Vector6d::Zero();
	};

};
