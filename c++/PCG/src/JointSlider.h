#pragma once

#include "Joint.h"

struct LinkageSystem;
struct State;
struct StateDeriv;
struct StateSolve;

struct JointSlider : public Joint
{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	JointSlider(std::unique_ptr<State> &S, 
		std::string pn,
		std::string cn,
		Eigen::Vector3d pp,
		Eigen::Vector3d cp,
		bool r,
		int i) :
	Joint(S, JType::Slider, pn, cn, pp, cp, r, i, 2)
	{};

	void updateMaxSparse(std::unique_ptr<StateSolve> &SS, const std::unique_ptr<LinkageSystem> &LS, const std::unique_ptr<State> &S, Matrix6d &parentAdj, Matrix6d &childAdj, Eigen::Vector3d &g) const;
	
private:
	void mapSelf(int c_i, std::unique_ptr<State> &S);
	void updateSelf();
	Vector6d Sqdot(const std::unique_ptr<State> &S)
	{
		return Vector6d::Zero();
	};

};