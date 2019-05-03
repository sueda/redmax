#include "Joint.h"

#include "Block.h"
#include "LinkageSystem.h"
#include "State.h"
#include "StateDeriv.h"
#include "StateSolve.h"

#include "ChronoTimer.h"

void Joint::DFSIndexSet(int &blockID, int &jointID, std::vector< std::pair < std::string, std::shared_ptr<Block> > > & blocks, std::vector<int>& j_map)
{
	// Recursively update the indicies
	std::shared_ptr<Block> block = nullptr;
	for (int i = 0; i < blocks.size(); ++i)
	{
		if (blocks[i].first.compare(childName) == 0)
			block = blocks[i].second;
	}

	int numChildren = (int)block->c_joints.size();
	for (int i = 0; i < numChildren; ++i)
	{
		block->c_joints[i]->DFSIndexSet(blockID, jointID, blocks, j_map);
	}

	j_map.insert(j_map.begin(), index0); // looping through joint_map will access all joints such that every parent is in front of it's children
	blockID = blockID + 1;
	jindex = blockID;
	constraint_index = jointID;
	jointID = jointID + constraintNum;
}

void Joint::initState(int constraint_index, const std::unique_ptr<LinkageSystem>& LS, std::unique_ptr<State>& S)
{
	if (constraintNum == 0)
		this->constraint_index = -1;
	else
		this->constraint_index = constraint_index - constraintNum;
	mapSelf(this->constraint_index, S);
	E_pj = E_pj_ij * Q;
	ad_jp = Rigid::adjoint(Rigid::inverse(E_pj));
	ad_jpT = ad_jp.transpose();

	V = Sqdot(S);
	if (root)
		E_wj = E_pj;
	else
	{
		E_wj = parentBlock->joint->E_wj * E_pj; // E_wj = E_wp
		V = V + ad_jp * (parentBlock->joint->V);
	}
	S->E[jindex] = E_wj * E_ji;
	E_iw = Rigid::inverse(S->E[jindex]);
	S->v.segment<6>(jindex * 6) = ad_ij * V;
	E_draw = S->E[jindex];

	compute_Ij(S);

	Pr = Eigen::MatrixXd::Zero(1, 1);
	Kmd = Matrix6d::Zero();
	Dmd = Matrix6d::Zero();

	// precomputations
	adij_S = ad_ij * this->S;
	ST_adijT = ST * ad_ijT;


	//////////////////////////////////////////////////////////////////////////////////////////////////////
	//S->ad_jp[jindex] = Rigid::adjoint(Rigid::inverse(E_pj)); // for no matrix w/ preconditioner
	//S->ad_jpT[jindex] = S->ad_jp[jindex].transpose();

	//S->V[jindex] = Sqdot(S);
	//if (root)
	//	E_wj = E_pj;
	//else
	//{
	//	E_wj = parentBlock->joint->E_wj * E_pj; // E_wj = E_wp
	//	S->V[jindex] = V + ad_jp * (parentBlock->joint->V);
	//}

	//S->Pr[jindex] = Eigen::MatrixXd::Zero(1, 1);
	//S->Kmd[jindex] = Matrix6d::Zero();
	//S->Dmd[jindex] = Matrix6d::Zero();

	//S->I_j[jindex] = I_j;
	//S->S[jindex] = this->S;
	//S->ST[jindex] = this->S.transpose();
	//S->Sdot[jindex] = this->Sdot;

	//// precomputations
	//S->adij_S[jindex] = ad_ij * this->S;
	//S->ST_adijT[jindex] = S->ST[jindex] * ad_ijT;
	//////////////////////////////////////////////////////////////////////////////////////////////////////


	addot_ = false;
	ad_ip_ = false;
	ad_iw_ = false;
	ad_wi_ = false;
	add_ip_ = false;
}

void Joint::initState_SoA(int constraint_index, const std::unique_ptr<LinkageSystem>& LS, std::unique_ptr<State>& S)
{
	if (constraintNum == 0)
		this->constraint_index = -1;
	else
		this->constraint_index = constraint_index - constraintNum;
	mapSelf(this->constraint_index, S);
	E_pj = E_pj_ij * Q;
	S->ad_jp[jindex] = Rigid::adjoint(Rigid::inverse(E_pj)); // for no matrix w/ preconditioner
	S->ad_jpT[jindex] = S->ad_jp[jindex].transpose();

	S->V[jindex] = Sqdot(S);
	if (root)
		E_wj = E_pj;
	else
	{
		E_wj = parentBlock->joint->E_wj * E_pj; // E_wj = E_wp
		S->V[jindex] = V + S->ad_jp[jindex] * (S->V[parentBlock->joint->jindex]);
	}

	// update maximal state information
	S->E[jindex] = E_wj * E_ji;
	E_iw = Rigid::inverse(S->E[jindex]);
	S->v.segment<6>(jindex * 6) = ad_ij * S->V[jindex];

	E_draw = S->E[jindex];
	compute_Ij(S);

	S->Pr[jindex] = Eigen::MatrixXd::Zero(1, 1);
	S->Kmd[jindex] = Matrix6d::Zero();
	S->Dmd[jindex] = Matrix6d::Zero();

	S->I_j[jindex] = I_j;
	S->S[jindex] = this->S;
	S->ST[jindex] = this->S.transpose();
	S->Sdot[jindex] = this->Sdot;

	// precomputations
	S->adij_S[jindex] = ad_ij * this->S;
	S->ST_adijT[jindex] = S->ST[jindex] * ad_ijT;

	addot_ = false;
	ad_ip_ = false;
	ad_iw_ = false;
	ad_wi_ = false;
	add_ip_ = false;
}

void Joint::update(const std::unique_ptr<LinkageSystem>& LS, std::unique_ptr<State>& S)
{
	updateSelf();
	E_pj = E_pj_ij * Q;
	ad_jp = Rigid::adjoint(Rigid::inverse(E_pj)); // for no matrix w/ preconditioner
	ad_jpT = ad_jp.transpose();

	V = Sqdot(S);
	if (root)
		E_wj = E_pj;
	else
	{
		E_wj = parentBlock->joint->E_wj * E_pj; // E_wj = E_wp
		V = V + ad_jp * (parentBlock->joint->V);
	}
	
	// update maximal state information
	S->E[jindex] = E_wj * E_ji;
	E_iw = Rigid::inverse(S->E[jindex]);
	S->v.segment<6>(jindex * 6) = Rigid::adjoint(Rigid::inverse(E_ji)) * V;

	compute_Ij(S);
	E_draw = S->E[jindex];

	//Pr = Eigen::MatrixXd::Zero(1, 1);  // does not change
	Kmd = Matrix6d::Zero();
	Dmd = Matrix6d::Zero();

	// precomputations
	adij_S = ad_ij * this->S;
	ST_adijT = ST * ad_ijT;


	/////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//S->ad_jp[jindex] = Rigid::adjoint(Rigid::inverse(E_pj)); // for no matrix w/ preconditioner
	//S->ad_jpT[jindex] = S->ad_jp[jindex].transpose();

	//S->V[jindex] = Sqdot(S);
	//if (root)
	//	E_wj = E_pj;
	//else
	//{
	//	E_wj = parentBlock->joint->E_wj * E_pj; // E_wj = E_wp
	//	S->V[jindex] = S->V[jindex] + S->addot[jindex] * (S->V[S->parentIndex[jindex]]);
	//}

	//S->I_j[jindex] = I_j;
	//S->S[jindex] = this->S;
	//S->ST[jindex] = S->S[jindex].transpose();
	//S->Sdot[jindex] = this->Sdot;

	////Pr = Eigen::MatrixXd::Zero(1, 1);  // does not change
	//S->Kmd[jindex] = Matrix6d::Zero();
	//S->Dmd[jindex] = Matrix6d::Zero();

	//// precomputations
	//S->adij_S[jindex] = S->ad_ij[jindex] * this->S;
	//S->ST_adijT[jindex] = S->ST[jindex] * S->ad_ijT[jindex];
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////


	addot_ = false;
	ad_ip_ = false;
	ad_iw_ = false;
	ad_wi_ = false;
	add_ip_ = false;
}

void Joint::update_SoA(const std::unique_ptr<LinkageSystem>& LS, std::unique_ptr<State>& S)
{
	updateSelf();
	E_pj = E_pj_ij * Q;
	S->ad_jp[jindex] = Rigid::adjoint(Rigid::inverse(E_pj)); // for no matrix w/ preconditioner
	S->ad_jpT[jindex] = S->ad_jp[jindex].transpose();

	S->V[jindex] = Sqdot(S);
	if (root)
		E_wj = E_pj;
	else
	{
		E_wj = parentBlock->joint->E_wj * E_pj; // E_wj = E_wp
		S->V[jindex] = S->V[jindex] + S->ad_jp[jindex] * (S->V[S->parentIndex[jindex]]);
	}

	// update maximal state information
	S->E[jindex] = E_wj * E_ji;
	E_iw = Rigid::inverse(S->E[jindex]);
	S->v.segment<6>(jindex * 6) = S->ad_ij[jindex] * S->V[jindex];

	compute_Ij(S);
	E_draw = S->E[jindex];

	S->I_j[jindex] = I_j;
	S->S[jindex] = this->S;
	S->ST[jindex] = S->S[jindex].transpose();
	S->Sdot[jindex] = this->Sdot;

	//Pr = Eigen::MatrixXd::Zero(1, 1);  // does not change
	S->Kmd[jindex] = Matrix6d::Zero();
	S->Dmd[jindex] = Matrix6d::Zero();

	// precomputations
	S->adij_S[jindex] = S->ad_ij[jindex] * this->S;
	S->ST_adijT[jindex] = S->ST[jindex] * S->ad_ijT[jindex];

	addot_ = false;
	ad_ip_ = false;
	ad_iw_ = false;
	ad_wi_ = false;
	add_ip_ = false;
}

void Joint::update_nostatechange(const std::unique_ptr<LinkageSystem>& LS)
{
	updateSelf();
	E_pj = E_pj0 * Q;

	if (root)
		E_wj = E_pj;
	else
	{
		std::shared_ptr<Block> pblock = LS->find(parentName);
		E_wj = pblock->joint->E_wj * E_pj; // E_wj = E_wp
	}
	std::shared_ptr<Block> cblock = LS->find(childName);
	cblock->E_wi0_inc = E_wj * E_ji0;
}

Json::Value Joint::exportJson(const std::unique_ptr<LinkageSystem> & LS, const std::unique_ptr<State> & S)
{
	Eigen::Vector3d pos = E_draw.block<3, 4>(0, 0)*Eigen::Vector4d(childPos[0], childPos[1], childPos[2], 1);
	Eigen::Quaterniond quat(E_draw.block<3, 3>(0, 0));
	//*(outfiles[0]) << "\"" << childName << "_joint\":\n{\n"
	//	<< "\"scale\": [" << std::to_string(1) << "," << std::to_string(1) << ","
	//	<< std::to_string(1) << "],\n\"location\": [" << std::to_string(pos[0]) << ","
	//	<< std::to_string(pos[1]) << "," << std::to_string(pos[2])
	//	<< "],\n\"quat\": [" << std::to_string(quat.w()) << "," << std::to_string(quat.x()) << ","
	//	<< std::to_string(quat.y()) << "," << std::to_string(quat.z()) << "]\n}";

	Eigen::Vector3d ppos, cpos;
	Eigen::Matrix4d E = Eigen::Matrix4d::Identity();
	int parenti, childi;
	if (!root) {
		parenti = LS->find(parentName)->joint->jindex;
		E = S->E[parenti];
		ppos = E.block<3, 3>(0, 0)*parentPos + E.block<3, 1>(0, 3);
	}
	else
	{
		ppos = parentPos;
	}

	Json::Value scale(Json::arrayValue);
	scale[0] = 0.5;
	scale[1] = 0.5;
	scale[2] = 0.5;

	Json::Value location(Json::arrayValue);
	location[0] = ppos[0];
	location[1] = ppos[1];
	location[2] = ppos[2];

	Json::Value q(Json::arrayValue);
	q[0] = quat.w();
	q[1] = quat.x();
	q[2] = quat.y();
	q[3] = quat.z();

	Json::Value v;
	v["scale"] = scale;
	v["location"] = location;
	v["quat"] = q;

	return v;
}

void Joint::exportBrender(std::vector<std::shared_ptr<std::ofstream>> outfiles, int frame, double time) const
{
	
}

void Joint::compute_Ij(std::unique_ptr<State>& S)
{
	//m = this.body.I_i(4);
	//R = this.body.E_ji(1:3, 1 : 3);
	//p = this.body.E_ji(1:3, 4);
	//pBrac = se3.brac(p);
	//Ic = diag(this.body.I_i(1:3));
	//this.I_j = [R*Ic*R'+m*(pBrac'*pBrac), m*pBrac; m*pBrac',m*eye(3)];

	double m = S->M[jindex](4, 4);
	Eigen::Matrix3d R = E_ji.block<3, 3>(0, 0);
	Eigen::Vector3d p = E_ji.block<3, 1>(0, 3);
	Eigen::Matrix3d pBrac = Rigid::bracket3(p);
	Eigen::Matrix3d Ic = Eigen::Matrix3d::Identity();
	Ic(0, 0) = S->M[jindex](0, 0);
	Ic(1, 1) = S->M[jindex](1, 1);
	Ic(2, 2) = S->M[jindex](2, 2);
	
	I_j.block<3, 3>(0,0) = R * Ic * R.transpose() + m * (pBrac.transpose() * pBrac);
	I_j.block<3, 3>(0, 3) = m * pBrac;
	I_j.block<3, 3>(3, 0) = I_j.block<3, 3>(0, 3).transpose();
	I_j.block<3, 3>(3, 3) = Eigen::Matrix3d::Identity() * m;
}
