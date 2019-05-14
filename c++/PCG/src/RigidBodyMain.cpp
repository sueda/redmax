#define _USE_MATH_DEFINES
#include <ctgmath>
#include <cmath>

#include "RigidBodyMain.h"

#include "ChronoTimer.h"
#include "Cuboid.h"

#ifdef ONLINE_MODE
#include "online/MatrixStack.h"
#include "online/Program.h"
#endif

#include <glm/gtc/type_ptr.hpp>
#include <json/writer.h>
#include <json/value.h>
#include <json/json.h>

#include <iostream>
#include <fstream>
#include <memory>

// MAC: reinterpret_cast
template <typename To, typename From>
inline std::shared_ptr<To> reinterpret_pointer_cast(
                                                    std::shared_ptr<From> const & ptr) noexcept
{ return std::shared_ptr<To>(ptr, reinterpret_cast<To *>(ptr.get())); }

static double randDouble(double l, double h)
{
	double r = rand() / (double)RAND_MAX;
	return (1.0 - r) * l + r * h;
}

std::set<simType> reducedCoordList{ PCG, PCG_unopt, Pardiso };
std::set<simType> reducedNoMatrixList{ PCG, PCG_unopt };
simType simtype = PCG;

RigidBodyMain::RigidBodyMain(const std::string & RESOURCE_DIR, double &t) :
	t(t)
{
	userActionState = userAction::free;

	// Create simulation variables
	LS = std::make_unique<LinkageSystem>();
	S = std::make_unique<State>(t);
	DS = std::make_unique<StateDeriv>();
	SS = std::make_unique<StateSolve>();
	J = std::make_unique<JSONwrapper>();

	creator = std::make_unique<RigidBodyCreator>(t);
	solver = std::make_unique<Solver>();

	// save resource location
	creator->setResourceDir(RESOURCE_DIR);

	// All objects have the same base shape
#ifdef ONLINE_MODE
	cuboid = std::make_unique<Cuboid>();
	cuboid->loadMesh(RESOURCE_DIR);
	cuboid->init();
#endif

	/// Default Values
	this->t = 0;

	SS->baumgarte = Eigen::Vector3d(0, 0, 100);

	// Timestep
	SS->h = 1e-2;
	SS->alpha = 0.0;

	// -981 cm/s^2
	SS->grav = Eigen::Vector3d(0.0, -980, 0.0);
	SS->mst = StateSolve::matrixSolveType::sparse;

	DS->c_dfdx.resize(ctMAX);
	DS->c_dfdv.resize(ctMAX);

	int size = 15;
	rand_init_vel.resize(size * 6);
	for (int i = 0; i < size; ++i)
	{
		rand_init_vel.segment<3>(i * 6) = Eigen::Vector3d(randDouble(0.25, 0.95), randDouble(0.25, 0.95), randDouble(0.25, 0.95));
		rand_init_vel.segment<3>(i * 6 + 3) = Eigen::Vector3d(randDouble(0.25, 0.95), randDouble(0.25, 0.95), randDouble(0.25, 0.95));
	}
}

RigidBodyMain::~RigidBodyMain()
{
}

void RigidBodyMain::load(const std::string & FILENAME)
{
	// Save file location
	creator->setUserSource(FILENAME);

	// Load the scenario
	creator->useSourceLinks();
	creator->useSourceSplines();
	creator->loadLinkagesfromFile(LS, S, SS, true);
	init();
}

void RigidBodyMain::loadTree(int n)
{
	creator->loadTree(n, simtype, LS, S, SS);
	init();
}

void RigidBodyMain::loadSimpleTree(int n)
{
	creator->loadSimpleTree(n, simtype, LS, S, SS);
	init();
}

void RigidBodyMain::loadBridge(int nbridge, int ntower)
{
	creator->loadBridge(nbridge, ntower, simtype, LS, S, SS);
	init();
}

void RigidBodyMain::loadSimpleBridge(int nbridge)
{
	creator->loadSimpleBridge(nbridge, simtype, LS, S, SS);
	init();
}

void RigidBodyMain::loadUmbrella(int n)
{
	creator->loadUmbrella(n, simtype, LS, S, SS);
	init();
}

void RigidBodyMain::loadChain(int n)
{
	creator->loadChain(n, simtype, LS, S, SS);
	init();
}

void RigidBodyMain::loadTest(int n)
{
	creator->loadTest(n, simtype, LS, S, SS);
	init();
}

void RigidBodyMain::reLoad(const bool disconnect)
{
	creator->loadLinkagesfromFile(LS, S, SS, disconnect);

	init();
}

void RigidBodyMain::updateInitialScene(std::vector<double>& x)
{
	Eigen::Matrix4d E;
	Eigen::Matrix3d R;

	// Update the rotation and position of each of the blocks
	for (int i = 0; i < LS->blocks.size(); i++)
	{
		E = Eigen::Matrix4d::Identity();
		R = Eigen::AngleAxisd(x[i * 3], Eigen::Vector3d::UnitZ());
		E.block<3, 3>(0, 0) = R;
		E(0, 3) = x[i * 3 + 1];
		E(1, 3) = x[i * 3 + 2];

		// Update the transformation matricies
		S->E[i] = E;
	}

	// Update the joints by reloading them with the new State
	ConstraintJoint::loadIntoLinkageSystem(LS);
}

void RigidBodyMain::init()
{
	int numObj = (int)LS->blocks.size();
	if (simtype == simType::PCG)
	{
		S->jindex.resize(numObj);
		S->type.resize(numObj);
		S->constraintNum.resize(numObj);
		S->constraint_index.resize(numObj);
		S->parentIndex.resize(numObj);
		S->childIndex.resize(numObj);
		S->root.resize(numObj);
		S->damping.resize(numObj);
		S->d.resize(numObj);
		S->k.resize(numObj);

		S->V.resize(numObj);
		S->S.resize(numObj);
		S->ST.resize(numObj);
		S->Sdot.resize(numObj);
		S->q0.resize(numObj);
		S->I_j.resize(numObj);
		S->ad_ij.resize(numObj);
		S->ad_ijT.resize(numObj);
		S->ad_jp.resize(numObj);
		S->ad_jpT.resize(numObj);
		S->addot.resize(numObj);
		S->ad_ip.resize(numObj);
		S->ad_ipT.resize(numObj);
		S->add_ip.resize(numObj);
		S->ad_iw.resize(numObj);
		S->ad_wi.resize(numObj);
		S->alpha.resize(numObj);
		S->beta.resize(numObj);
		S->Pr.resize(numObj);
		S->Kmd.resize(numObj);
		S->Dmd.resize(numObj);
		S->D.resize(numObj);
		S->Mhat.resize(numObj);
		S->Psi.resize(numObj);
		S->Pi.resize(numObj);
		S->Bhat.resize(numObj);
		S->Vdot.resize(numObj);
		S->adij_S.resize(numObj);
		S->ST_adijT.resize(numObj);
		S->ST_Mhat.resize(numObj);
		S->S_Psi.resize(numObj);
		S->Mhat_S_Psi.resize(numObj);
		S->ST_Bhat.resize(numObj);
	}

	// reverse DFS oredering for reduced coords
	int num_constraints;
	std::shared_ptr<Block> c_block;
 	if (reducedCoordList.find(simtype) != reducedCoordList.end())
	{
		// attach backwards joint pointers (child joint/block pairs)
		for (int i = 0; i < LS->joints.size(); ++i)
		{
			if (!LS->joints[i]->root)
			{
				// point to the child joint from parent
				c_block = LS->find(LS->joints[i]->childName);
				LS->find(LS->joints[i]->parentName)->c_joints.push_back(c_block->joint);
			}
			// point every joint to its block
			LS->joints[i]->childBlock = LS->find(LS->joints[i]->childName);
			LS->joints[i]->parentBlock = LS->find(LS->joints[i]->parentName);
		}
		ConstraintJoint::loadOptimalJointOrdering(LS, S);

		//for (int i = 0; i < (int)LS->blocks.size(); ++i)
		//{
		//	// update pos and vel
		//	std::cout << LS->joints[LS->joint_map[i]]->childName << " : " << LS->joints[LS->joint_map[i]]->jindex << " " << LS->joints[LS->joint_map[i]]->constraint_index << std::endl;
		//	//LS->joints[LS->joint_map[i]]->update(LS, S);
		//}
		num_constraints = ConstraintJoint::getConstraintNumReduced(LS);
	}
	else // maximal coordinates dont use joint ordering
	{
		num_constraints = ConstraintJoint::getConstraintNumMaximal(LS);
	}

	ConstraintJoint::loadIntoLinkageSystem(LS);

	// initialize the state containers
	S->qddot = Eigen::VectorXd::Zero(num_constraints);
	S->v = Eigen::VectorXd::Zero(numObj * 6);
	S->qdot = Eigen::VectorXd::Zero(num_constraints);
	S->E.resize(numObj);
	S->q = Eigen::VectorXd::Zero(num_constraints);
	S->mass.resize(numObj);
	S->M.resize(numObj);

	// load the init State
	std::shared_ptr<Block> block;
	for (int i = 0; i < numObj; ++i)
	{
		// link up the block state variables
		block = LS->blocks[i].second;
		int index = block->joint->jindex;

		S->v.segment<6>(index * 6) = block->v0;
		S->E[index] = block->E_wi0;
		S->M[index] = block->M0;
		S->mass[index] = block->mass0;

		// load spatial inertia
		LS->blocks[i].second->updateSpatialIntertia(S);
	}

	// Load joint specific joint variables
	ConstraintJoint::init(SS, LS, S);
	
	for (int i = 0; i < LS->constraints.size(); ++i)
	{
		if (LS->constraints[i]->type == Constraint::constraintType::elastic)
			std::reinterpret_pointer_cast<Elastic>(LS->constraints[i])->initLengthToRest(LS, S);
	}
	//ConstraintJoint::loadIntoState(SS, LS, S);

	// load SoA constants
	if (simtype == simType::PCG)
	{
		// forward traversal
		int jindex;
		std::shared_ptr<Joint> j;
		for (int jmapi = 0; jmapi < LS->joints.size(); ++jmapi) 
		{
			int list_index = LS->joint_map[jmapi];
			j = LS->joints[list_index];
			jindex = j->jindex;

			S->jindex[list_index] = jindex;
			S->type[jindex] = j->type;
			S->constraintNum[jindex] = j->constraintNum;
			S->constraint_index[jindex] = j->constraint_index;
			if (j->parentBlock != nullptr)
				S->parentIndex[jindex] = j->parentBlock->joint->jindex;
			for (int ci = 0; ci < j->childBlock->c_joints.size(); ++ci)
			{
				S->childIndex[jindex].push_back(j->childBlock->c_joints[ci]->jindex);
			}
			S->root[jindex] = j->root;
			S->damping[jindex] = j->childBlock->damping;
			S->d[jindex] = j->d;
			S->k[jindex] = j->k;

			S->ad_ij[jindex] = j->ad_ij;
			S->ad_ijT[jindex] = j->ad_ijT;

			S->Pr[jindex] = j->Pr;
			S->q0[jindex] = j->q0;
			S->D[jindex] = j->childBlock->D;
		}
	}
}

void RigidBodyMain::reset()
{
	t = 0.0;
	LS->blocks.clear();
	LS->joints.clear();
	LS->constraints.clear();
	LS->joint_map.clear();

	S->v = Eigen::VectorXd::Zero(0);
	S->E.clear();
	S->M.clear();
	S->mass.clear();
	S->q = Eigen::VectorXd::Zero(0);
	S->qdot = Eigen::VectorXd::Zero(0);
	S->qddot = Eigen::VectorXd::Zero(0);

	S->jindex.clear();
	S->type.clear();
	S->constraintNum.clear();
	S->constraint_index.clear();
	S->parentIndex.clear();
	S->childIndex.clear();
	S->root.clear();
	S->damping.clear();
	S->d.clear();
	S->k.clear();

	S->V.clear();
	S->S.clear();
	S->ST.clear();
	S->Sdot.clear();
	S->q0.clear();
	S->I_j.clear();
	S->ad_ij.clear();
	S->ad_ijT.clear();
	S->ad_jp.clear();
	S->ad_jpT.clear();
	S->addot.clear();
	S->ad_ip.clear();
	S->ad_ipT.clear();
	S->add_ip.clear();
	S->ad_iw.clear();
	S->ad_wi.clear();
	S->alpha.clear();
	S->beta.clear();
	S->Pr.clear();
	S->Kmd.clear();
	S->Dmd.clear();
	S->D.clear();
	S->Mhat.clear();
	S->Psi.clear();
	S->Pi.clear();
	S->Bhat.clear();
	S->Vdot.clear();
	S->adij_S.clear();
	S->ST_adijT.clear();
	S->ST_Mhat.clear();
	S->S_Psi.clear();
	S->Mhat_S_Psi.clear();
	S->ST_Bhat.clear();

	DS->qBuf = Eigen::VectorXd::Zero(0);

	fullyConverged = false;
	ct = 0;
}

Eigen::VectorXd RigidBodyMain::step()
{
	//for (int i = 0; i < LS->blocks.size(); ++i)
	//{
	//	std::cout << LS->blocks[i].second->joint->jindex << " : " << LS->blocks[i].second->joint->childName << " parent: " << LS->blocks[i].second->joint->parentName << std::endl;
	//}

	// update the time
	t = t + SS->h;
#ifdef ONLINE_MODE
	std::cout << t << std::endl;
#endif

	// solve with timestep SS->h
	Eigen::VectorXd result = solver->solve(SS, LS, S);

	update_user_forces();

	return result;
}

void RigidBodyMain::update_user_forces()
{
	// HARDCODED

	//////BRIDGE
	//if (t < 4.99999)
	//{
	//	//double blocklength = 0.240392;	// for 100

	//	std::shared_ptr<SpringPoint> car = std::reinterpret_pointer_cast<SpringPoint>(LS->constraints[LS->constraints.size() - 1]);
	//	double blocklength = car->blocka->size[0];
	//	double pos_step = 24.0 / (5.0/SS->h); // distance per step
	//	//double pos_on_bridge = (t / 5.0) * 24; //decklength = 24
	//	double pos_new = car->posa[0] + pos_step;

	//	//std::cout << pos_step << " " << pos_new << " ";
	//	if (pos_new > (blocklength / 2.0))	// check that the new pos is still associated with the current block
	//	{
	//		pos_new = fmod(pos_new, (blocklength / 2.0));		
	//		pos_new = -(blocklength / 2.0) + pos_new;			// if not update the pos
	//		int block_index = car->blocka->c_joints[0]->index0; // GET CHILD BLOCK
	//		//std::cout << car->blocka->joint->jindex << " new: " << block_index << " " << LS->blocks[block_index].second->joint->childName << " ";
	//		car->blocka = LS->blocks[block_index].second;		// and update the block associated
	//	}
	//	car->posa = Eigen::Vector3d(pos_new, 0, 0);				// move the car
	//	//std::cout << pos_new << std::endl;
	//}

	// UMBRELLA
	// check if the runner needs to be stopped
	//std::shared_ptr<Block> runner = LS->find("runner");
	//Eigen::Vector3d runner_pos = S->E[runner->joint->jindex].block<3,1>(0,3);
	//// stop runner at 44
	//if (runner_pos[1] >= 36)
	//{
	//	std::reinterpret_pointer_cast<SpringPoint>(LS->constraints[LS->constraints.size() - 1])->k =0;
	//}

	//if (t > 4.4 && t < 5.0)
	//{
	//	//std::reinterpret_pointer_cast<SpringPoint>(LS->constraints[LS->constraints.size() - 2])->k = 200;
	//}
	//else if(t >= 5.0)
	//{
	//	std::reinterpret_pointer_cast<SpringPoint>(LS->constraints[LS->constraints.size() - 2])->k = 0;
	//}

}

void RigidBodyMain::draw(std::unique_ptr<MatrixStack>& MV, const std::unique_ptr<Program>& prog) const
{
#ifdef ONLINE_MODE
	std::shared_ptr<Block> block;
	Eigen::Matrix4d size4;
	Eigen::Matrix4d E;
	glm::mat4 temp;
	for (int i = 0; i < LS->blocks.size(); i++)
	{
		block = LS->blocks[i].second;
		//std::cout << LS->keys[i] << " " << block->joint->jindex << std::endl;
		size4 = Eigen::Matrix4d::Identity();

		E = S->E[block->joint->jindex];
		
		if (!block->invisible)
		{
			glUniform3fv(prog->getUniform("kdFront"), 1, Eigen::Vector3f(1.0, 0.0, 0.0).data());
			glUniform3fv(prog->getUniform("kdBack"), 1, Eigen::Vector3f(1.0, 1.0, 0.0).data());

			if (block->shape == Block::SType::Cuboid)
			{
				// If the inertial block is separate from display
				if (block->hasInertial)
				{
					// Draw display block
					MV->pushMatrix();
					size4(0, 0) = block->size[0] / 2;
					size4(1, 1) = block->size[1] / 2;
					size4(2, 2) = block->size[2] / 2;
					size4 = E * block->iToD * size4;
					for (int r = 0; r < 4; ++r)
					{
						for (int c = 0; c < 4; ++c)
						{
							temp[c][r] = (float)size4(r, c);
						}
					}
					MV->multMatrix(temp);

					glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
					cuboid->draw(prog);

					MV->popMatrix();

					if (showInertial)
					{
						// Draw the inertial block
						glUniform3fv(prog->getUniform("kdFront"), 1, Eigen::Vector3f(0.0, 1.0, 1.0).data());
						glUniform3fv(prog->getUniform("kdBack"), 1, Eigen::Vector3f(0.0, 1.0, 0.0).data());

						MV->pushMatrix();
						size4(0, 0) = block->insize[0] / 2;
						size4(1, 1) = block->insize[1] / 2;
						size4(2, 2) = block->insize[2] / 2;
						size4 = E * size4;
						for (int r = 0; r < 4; ++r)
						{
							for (int c = 0; c < 4; ++c)
							{
								temp[c][r] = (float)size4(r, c);
							}
						}
						MV->multMatrix(temp);

						glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
						cuboid->draw(prog);

						MV->popMatrix();
					}

				}
				else
				{
					MV->pushMatrix();

					glUniform3fv(prog->getUniform("kdFront"), 1, Eigen::Vector3f(1.0, 0.0, 0.0).data());
					glUniform3fv(prog->getUniform("kdBack"), 1, Eigen::Vector3f(1.0, 1.0, 0.0).data());

					// The display and the inertial blocks are the same
					size4(0, 0) = block->size[0] / 2;
					size4(1, 1) = block->size[1] / 2;
					size4(2, 2) = block->size[2] / 2;
					size4 = E * size4;
					for (int r = 0; r < 4; ++r)
					{
						for (int c = 0; c < 4; ++c)
						{
							temp[c][r] = (float)size4(r, c);
						}
					}
					MV->multMatrix(temp);

					glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
					cuboid->draw(prog);

					MV->popMatrix();
				}
			}
			else
			{
				MV->pushMatrix();

				// The display and the inertial blocks are the same
				size4(0, 0) = block->size[0] / 2;
				size4(1, 1) = block->size[0] / 2;
				size4(2, 2) = block->size[1] / 2;
				size4 = E * size4;
				for (int r = 0; r < 4; ++r)
				{
					for (int c = 0; c < 4; ++c)
					{
						temp[c][r] = (float)size4(r, c);
					}
				}
				MV->multMatrix(temp);

				glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
				cuboid->draw(prog);

				MV->popMatrix();
			}
		}
	}
#endif
}

void RigidBodyMain::drawJoints(std::unique_ptr<MatrixStack>& MV) const
{
#ifdef ONLINE_MODE
	MV->pushMatrix();

	// Draw the joint pins for visual help!
	ConstraintJoint::draw(SS, LS, S);

	MV->popMatrix();
#endif
}

void RigidBodyMain::drawPoints(std::unique_ptr<MatrixStack>& MV, const std::unique_ptr<Program>& prog) const
{
#ifdef ONLINE_MODE
	MV->pushMatrix();
	glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));

	// Draw the other constraints!
	for (int i = 0; i < LS->constraints.size(); ++i)
	{
		LS->constraints[i]->draw(LS, S);
	}

	MV->popMatrix();
#endif
}

void RigidBodyMain::displayUserActions()
{
	std::cout << "\nRigidBodyMain User Actions:" << std::endl <<
		"\t : Spacebar will prompt continuous calls to step" << std::endl <<
		"\th: Step the state forward one timestep" << std::endl <<
		"\tj: Select the simulation type, DEFAULT - maximal dense solve" << std::endl <<
		"\tr: Reset the simulation from the saved configuration" << std::endl <<
		std::endl;
}

int RigidBodyMain::handleKeyPress(unsigned int key)
{
	if (isUserActionOverride())
	{
		return handleOverrideInput(key);
	}

	switch (key) {
	case 'h':
		return userAskForAction(RigidBodyMain::userAction::singleStep);
	case 'j':
		return userAskForAction(RigidBodyMain::userAction::selectSimType);
	case 'r':
		reset();
		reLoad(true);
		init();
		break;
	default:
		return false;
	}

	return true;
}

bool RigidBodyMain::askForAction(userAction a)
{
	if (userActionState == userAction::free)
	{
		userActionState = a;
		return true;
	}
	else return false;
}

bool RigidBodyMain::userAskForAction(userAction a)
{
	switch (a) {
	case userAction::selectSimType:
		if (askForAction(a))
		{
			std::cout << "Simulation Type Select:" << std::endl <<
				"\t1: Maximal Coordinates using sparse matricies" << std::endl <<
				"\t2: Maximal Coordinates using dense matricies" << std::endl <<
				"\t3: Sueda's Reduced Coordinates using sparse matricies" << std::endl <<
				"\t4: Sueda's Reduced Coordinates using dense matricies" << std::endl <<
				"\t5: Inverse Kinematics forced solve using sparse matricies" << std::endl << 
				std::endl;
			return true;
		}
		break;
	case userAction::singleStep:
		if (askForAction(a))
		{
			std::cout << t << std::endl;
			step();
			userActionState = userAction::free;
			return true;
		}
		break;
	}

	return false;
}

bool RigidBodyMain::isUserActionOverride()
{
	// if the blocking action list does not contain the action state
	if (blockingAction.find(userActionState) == blockingAction.end())
		return false; // the action is non - blocking
	else return true; // else override user input and send it to the action handler
}

int RigidBodyMain::handleOverrideInput(unsigned int key)
{
	switch (userActionState) {
	case userAction::selectSimType:
		return handleSelectSimType(key);
	default:
		return false;
	}
}

double RigidBodyMain::getCertNormSum()
{
	return S->qdot.norm();
}

Eigen::Vector3d RigidBodyMain::getLeafCertificate()
{
	return S->E[LS->joint_map[LS->blocks.size()-3]].block<3, 1>(0, 3);
}

void RigidBodyMain::toggleBlockDisplayMode()
{
	reset();
	setSolveType(simType::PCG);
	///rigid_body->loadBridge(pow(2,(n+3)), 1);
	loadTree(7);
}

void RigidBodyMain::exportBrender(std::vector<std::shared_ptr<std::ofstream>> outfiles, int frame, double time) const
{
	// print header
	if (export_part == 0)
	{
		Json::Value objs(Json::arrayValue);
		objs[0] = "D:/Git/Collab/RigidBodyMaster/resources/cube.obj";
		objs[1] = "D:/Git/Collab/RigidBodyMaster/resources/sphere2.obj";

		Json::Value states(Json::arrayValue);
		std::vector<Json::Value> v1; //(Json::objectValue);
		for (int i = 0; i < LS->blocks.size(); ++i)
		{
			v1.push_back(Json::Value());
			v1[i]["obj"] = 0;
			v1[i]["name"] = LS->blocks[i].second->joint->childName.c_str();
			v1[i]["group"] = "bridge";
			states.append(v1[i]);
		}
		for (int i = 0; i < LS->joints.size(); ++i)
		{
			v1.push_back(Json::Value());
			v1[v1.size() - 1]["obj"] = 1;
			v1[v1.size() - 1]["name"] = (LS->blocks[i].second->joint->childName + "_joint").c_str();
			v1[v1.size() - 1]["group"] = "joint";
			states.append(v1[v1.size() - 1]);
		}
		for (auto c : LS->constraints) 
		{
			if (c->type == Constraint::constraintType::springpoint)
			{
				v1.push_back(Json::Value());
				v1[v1.size() - 1]["obj"] = 0;
				v1[v1.size() - 1]["name"] = "car";
				v1[v1.size() - 1]["group"] = "car";
				states.append(v1[v1.size() - 1]);
			}
		}

		//for (int i = 0; i < LS->joints.size(); ++i)
		//{
		//	Json::Value v(Json::objectValue);
		//	v["obj"] = 1;
		//	v["name"] = LS->blocks[i].second->joint->childName + "_joint";
		//	states.append(v);
		//}

		//for (int i = 0; i < LS->constraints.size(); ++i)
		//{
		//	// implement
		//}

		//Json::Value jsonscene;
		J->jsonscene["header"]["objs"] = objs;
		J->jsonscene["header"]["states"] = states;

		Json::Value v(Json::objectValue);
		v["frame"] = frame;

		for (int i = 0; i < LS->blocks.size(); ++i)
		{
			v[LS->blocks[i].second->joint->childName.c_str()] = LS->blocks[i].second->exportJson(S, J);
			//LS->blocks[i].second->exportJson(S, J);
			//v[LS->blocks[i].second->joint->childName] = J->blocks[J->blocks.size() - 1];  
		}
		for (int i = 0; i < LS->joints.size(); ++i)
		{
			v[(LS->blocks[i].second->joint->childName + "_joint").c_str()] = LS->joints[i]->exportJson(LS, S);
		}

		//for (int i = 0; i < LS->joints.size(); ++i)
		//{
		//	LS->blocks[i].second->joint->exportJson(v, LS, S);
		//}

		

		for (auto c : LS->constraints) {
			c->exportJson(*outfiles[1], S);

			if (c->type == Constraint::constraintType::springpoint)
				v["car"] = c->exportBrender();
		}

		J->frames.append(v);

		//*(outfiles[0]) << jsonscene;

		//*(outfiles[0]) << "{\n\"header\":\n{\n\"objs\":\n[\n\"D:RigidBodyMaster/resources/cube.obj\",\n"
		//	<< "\"D:RigidBodyMaster/resources/sphere2.obj\"\n],\n\"states\":\n[\n";
		//for (int i = 0; i < LS->blocks.size(); ++i)
		//{
		//	*(outfiles[0]) << "{\n\"obj\": 1,\n\"name\": \"" << LS->blocks[i].second->joint->childName << "\"\n}";
		//	if (!(i == LS->blocks.size() - 1 && LS->joints.size() == 0))
		//		*(outfiles[0]) << ",";
		//	*(outfiles[0]) << "\n";
		//}

		//for (int i = 0; i < LS->joints.size(); ++i)
		//{
		//	*(outfiles[0]) << "{\n\"obj\": 0,\n\"name\": \"" << LS->blocks[i].second->joint->childName << "_joint\"\n}";
		//	//if (!(i == LS->joints.size() - 1 && LS->constraints.size() == 0))
		//	if (i != LS->joints.size() - 1)
		//		*(outfiles[0]) << ",";
		//	*(outfiles[0]) << "\n";
		//}

		//*(outfiles[0]) << "]\n},\n\n\"body\"[\n{\n}\n]";
	}
	// print frame
	else if (export_part == 1)
	{
		Json::Value v(Json::objectValue);
		v["frame"] = frame;

		//Json::Value objs(Json::objectValue);
		//objs["scale"] = "D:RigidBodyMaster/resources/cube.obj";
		//objs["location"] = "D:RigidBodyMaster/resources/cube.obj";
		//v["my_part_name"] = objs;

		for (int i = 0; i < LS->blocks.size(); ++i)
		{
			v[LS->blocks[i].second->joint->childName.c_str()] = LS->blocks[i].second->exportJson(S, J);

		}
		for (int i = 0; i < LS->joints.size(); ++i)
		{
			v[(LS->blocks[i].second->joint->childName + "_joint").c_str()] = LS->joints[i]->exportJson(LS, S);
		}

		//for (int i = 0; i < LS->joints.size(); ++i)
		//{
		//	LS->blocks[i].second->joint->exportJson(v, LS, S);
		//}

		//jsonscene["body"] = frames;

		for (auto c : LS->constraints) {
			c->exportJson(*outfiles[1], S);

			if (c->type == Constraint::constraintType::springpoint)
				v["car"] = c->exportBrender();
		}

		J->frames.append(v);

		//*(outfiles[0]) << ",\n{\n\"frame\": " << std::to_string(frame) << ",\n";
		//for (int i = 0; i < LS->blocks.size(); ++i)
		//{
		//	if (i != 0)
		//		*(outfiles[0]) << ",\n";
		//	*(outfiles[0]) << "\"" << LS->blocks[i].second->joint->childName << "\":\n";
		//	LS->blocks[i].second->exportBrender(outfiles, frame, time);
		//}

		//for (int i = 0; i < LS->joints.size(); ++i)
		//{
		//	*(outfiles[0]) << ",\n\"" << LS->blocks[i].second->joint->childName << "_joint\":\n";
		//	LS->joints[i]->exportBrender(outfiles, frame, time);
		//}

		//*(outfiles[0]) << "\n}\n]\n}";
	}
	else {
		J->jsonscene["body"] = J->frames;
		*(outfiles[0]) << J->jsonscene;
	}
}

std::vector<std::string> RigidBodyMain::getBrenderExtensions() const
{
	//return std::vector<std::string>(1, "json");
	std::vector<std::string> exts;
	exts.push_back("json");
	exts.push_back("strand");
	return exts;
}

std::vector<std::string> RigidBodyMain::getBrenderNames() const
{
	//return std::vector<std::string>(1, "Scene");
	std::vector<std::string> names;
	names.push_back("Scene");
	names.push_back("Elastic");
	return names;
}

std::vector<int> RigidBodyMain::getBrenderTypes() const
{
	//if (export_part == 0) {
	//	return std::vector<int>(1, ResetAppend);
	//}
	//else {
	//	return std::vector<int>(1, Append);
	//}
	//return std::vector<int>(1, ResetAppend);
	std::vector<int> types;
	types.push_back(ResetAppend);
	types.push_back(Truncate);
	return types;
}

int RigidBodyMain::getBrenderCount() const
{
	return 2;
}

bool RigidBodyMain::handleSelectSimType(unsigned int key)
{
	switch (key) {
	case '1':
		simtype = simType::PCG;
		std::cout << "Set to use reduced coord with sparse solve" << std::endl;
		break;
	case '2':
		simtype = simType::Pardiso;
		std::cout << "Set to use reduced coord with dense solve" << std::endl;
		break;
	default:
		return false;
	}

	// selection complete - free user action state
	userActionState = userAction::free;
	return true;
}
