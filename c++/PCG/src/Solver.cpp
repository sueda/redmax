#include "Solver.h"

#include "RigidBodyMain.h"
#include "ChronoTimer.h"

#ifndef EIGEN_NO_STATIC_ASSERT
#define EIGEN_NO_STATIC_ASSERT
#endif

#ifdef REDMAX_PARDISO
#include <Eigen/PardisoSupport>
#endif

#include <Eigen/SparseCholesky>	
#include <Eigen/OrderingMethods>
#include <Eigen/Core>
#include <Eigen/Geometry>

#include <fstream>
#include <iomanip>
#include <chrono>
int count = 0;

void Solver::loadMSparse(std::unique_ptr<StateSolve>& SS, std::unique_ptr<LinkageSystem>& LS, std::unique_ptr<State>& S)
{
	int numObj = (int)LS->blocks.size();

	int index;
	std::shared_ptr<Block> block;
	for (int i = 0; i < LS->blocks.size(); ++i)
	{
		block = LS->blocks[i].second;
		index = block->joint->jindex;

		for (int j = 0; j < 6; ++j)
		{
			for (int k = 0; k < 6; ++k)
			{
				if (S->M[index](j, k) != 0)
				{
					SS->LHSlist.push_back(Eigen::Triplet<double>(index * 6 + j, index * 6 + k, S->M[index](j, k) * (1 + SS->alpha*SS->h)));
				}
			}
		}
	}
}

void Solver::computeRHS(std::unique_ptr<StateSolve>& SS, const std::unique_ptr<LinkageSystem>& LS, const std::unique_ptr<State>& S)
{
	Matrix6d spcp = Matrix6d::Zero();
	Vector6d bodyForces = Vector6d::Zero();
	int i;
	for (int index = 0; index < LS->blocks.size(); ++index)
	{
		if (simtype == simType::PCG)
			i = S->jindex[index];
		else
			i = LS->blocks[index].second->joint->jindex;

		// Define the spatial cross product matrix
		spcp.block<3, 3>(0, 0) << Rigid::bracket3(S->v.segment<3>(i * 6));
		spcp.block<3, 3>(3, 3) << spcp.block<3, 3>(0, 0);

		// Define body forces
		bodyForces.segment<3>(3) << S->E[i].block<3, 3>(0, 0).transpose()*SS->grav*S->mass[i];

		// Compute RHS vector
		if (SS->fn.size() != 0)
		{
			SS->fn.segment<6>(i * 6) += S->M[i] * S->v.segment<6>(i * 6) + SS->h*(spcp.transpose()*S->M[i] * S->v.segment<6>(i * 6) + bodyForces);
		}

		// Compute force vector
		if (SS->f.size() != 0)
		{
			SS->f.segment<6>(i * 6) += spcp.transpose()*S->M[i] * S->v.segment<6>(i * 6) + bodyForces;
		}
	}
}

void Solver::pcdSaad2003(Eigen::VectorXd &result, const Eigen::VectorXd &LHSx0, const Eigen::VectorXd &b, std::unique_ptr<StateSolve> &SS, std::unique_ptr<LinkageSystem> &LS, std::unique_ptr<State> &S, std::shared_ptr<State::local_mt> lmt, double tol, int maxit)
{
	//function[x, iter, xs, rs] = pcgSaad2003(A, b, tol, maxit, M, x0)
	//	% pcgSaad2003 Algorithm 9.1 from[Saad 2003]
	std::vector<double> rs;

	Eigen::VectorXd xj = S->qdot;
	Eigen::VectorXd rj = b - LHSx0;
	Eigen::VectorXd zj = ConstraintJoint::computeMinv_x(rj, SS, LS, S, lmt);

	//if (zj.hasNaN())
	//	std::cout << "minv init" << std::endl;

	Eigen::VectorXd pj = zj;
	double tolsq = tol * tol;
	double res1 = rj.dot(rj);
	double res0 = res1;
	Eigen::VectorXd xs = xj;
	rs.push_back(res1);

	auto start = std::chrono::steady_clock::now();
	Eigen::VectorXd J_p;
	Eigen::VectorXd q;
	Eigen::VectorXd Apj;
	Eigen::VectorXd rj1;
	Eigen::VectorXd zj1;
	int j = 1;
	for (; j < maxit; ++j)
	{
		ConstraintJoint::computeJ_x(J_p, pj, SS, LS, S);
		Apj = ConstraintJoint::computeJT_x(
			ConstraintJoint::computeLHS_x(J_p, SS, LS, S),
			SS, LS, S, lmt);

		//if (J_p.hasNaN())
		//	std::cout << "J * qd" << std::endl;
		//if (Apj.hasNaN())
		//	std::cout << "LHS or JT" << std::endl;

		q = S->q + SS->h * pj;
		// Joint stiffness and damping
		ConstraintJoint::computeStiffnessDampingJoint(Apj, pj, S->q, SS, LS, S);

		//if (Apj.hasNaN())
		//	std::cout << "stiffdampjoint" << std::endl;

		double aj = rj.dot(zj) / (Apj).dot(pj);
		if ((Apj).dot(pj) == 0)
		{
			xj = S->qdot;
			break;
		}
		xj = xj + aj * pj;
		rj1 = rj - aj * Apj;
		res1 = rj1.dot(rj1);
		rs.push_back(res1);
		if (res1 < tolsq*res0)
			break;

		zj1 = ConstraintJoint::computeMinv_x(rj1, SS, LS, S, lmt);

		//if (zj1.hasNaN())
		//	std::cout << "minv" << std::endl;
		//std::cout << j << std::endl;

		double bj = rj1.dot(zj1) / rj.dot(zj);
		pj = zj1 + bj * pj;
		rj = rj1;
		zj = zj1;
	}
	auto end = std::chrono::steady_clock::now();
	auto diff = end - start;
	auto time = std::chrono::duration<double, std::nano>(diff).count();
	trackLastTimestep.time_in_solve = time;
	trackLastTimestep.num_iterations = j;

	if (false)//std::abs(S->t - 1.01)) < 1e-5)
	{
		for (int i = 0; i < rs.size(); ++i)
		{
			std::cout << rs[i] << " ";
		}
		std::cout << std::endl;
		std::cout << "---------" << j << std::endl;
	}
	result = xj;
}

void Solver::pcdSaad2003_unopt(Eigen::VectorXd & result, const Eigen::VectorXd & LHSx0, const Eigen::VectorXd & b, std::unique_ptr<StateSolve>& SS, std::unique_ptr<LinkageSystem>& LS, std::unique_ptr<State>& S, std::shared_ptr<State::local_mt> lmt, double tol, int maxit)
{
	//function[x, iter, xs, rs] = pcgSaad2003(A, b, tol, maxit, M, x0)
	//	% pcgSaad2003 Algorithm 9.1 from[Saad 2003]
	std::vector<double> rs;

	Eigen::VectorXd xj = S->qdot;
	Eigen::VectorXd rj = b - LHSx0;
	Eigen::VectorXd zj = ConstraintJoint::computeMinv_x(rj, SS, LS, S, lmt);
	Eigen::VectorXd pj = zj;
	double tolsq = tol * tol;
	double res1 = rj.dot(rj);
	double res0 = res1;
	Eigen::VectorXd xs = xj;
	rs.push_back(res1);

	auto start = std::chrono::steady_clock::now();
	Eigen::VectorXd J_p;
	Eigen::VectorXd q;
	Eigen::VectorXd Apj;
	Eigen::VectorXd rj1;
	Eigen::VectorXd zj1;
	int j = 1;
	for (; j < maxit; ++j)
	{
		ConstraintJoint::computeJ_x_unopt(J_p, pj, SS, LS, S);
		Apj = ConstraintJoint::computeJT_x_unopt(
			ConstraintJoint::computeLHS_x_unopt(J_p, SS, LS, S),
			SS, LS, S, lmt);

		q = S->q + SS->h * pj;
		// Joint stiffness and damping
		ConstraintJoint::computeStiffnessDampingJoint_unopt(Apj, pj, S->q, SS, LS, S);

		double aj = rj.dot(zj) / (Apj).dot(pj);
		if ((Apj).dot(pj) == 0)
		{
			xj = S->qdot;
			break;
		}
		xj = xj + aj * pj;
		rj1 = rj - aj * Apj;
		res1 = rj1.dot(rj1);
		rs.push_back(res1);
		if (res1 < tolsq*res0)
			break;

		zj1 = ConstraintJoint::computeMinv_x_unopt(rj1, SS, LS, S, lmt);
		double bj = rj1.dot(zj1) / rj.dot(zj);
		pj = zj1 + bj * pj;
		rj = rj1;
		zj = zj1;
	}
	auto end = std::chrono::steady_clock::now();
	auto diff = end - start;
	auto time = std::chrono::duration<double, std::nano>(diff).count();
	trackLastTimestep.time_in_solve = time;
	trackLastTimestep.num_iterations = j;

	if (false)
	{
		for (int i = 0; i < rs.size(); ++i)
		{
			std::cout << rs[i] << " ";
		}
		std::cout << std::endl;
		std::cout << "---------" << j << std::endl;
	}
	result = xj;
}
#ifdef REDMAX_PARDISO
Eigen::VectorXd Solver::solvePardiso(std::unique_ptr<StateSolve>& SS, std::unique_ptr<LinkageSystem>& LS, std::unique_ptr<State>& S)
{
	//ChronoTimer ctime("ldlt", 5);
	//ChronoTimer ltime("load", 7);
	//ctime.tic(0);
	//ctime.tic(1);
	// Preprocessing, eliminate matrix resizing
	//ltime.tic(0);
 	int constraints = ConstraintJoint::getConstraintNumReduced(LS);
	int additional_constraints = ConstraintJoint::getConstraintNumAdditnl(LS);
	int numObj = (int)LS->blocks.size();

	// Initialize containers
	SS->J = Eigen::MatrixXd::Zero(numObj * 6, constraints); // jacobian
	SS->Jd = Eigen::MatrixXd::Zero(numObj * 6, constraints); // jacobian time derivative
	SS->f = Eigen::VectorXd::Zero(numObj * 6); // load with fm
	SS->fr = Eigen::VectorXd::Zero(constraints);
	SS->gm = Eigen::VectorXd::Zero(additional_constraints);
	SS->Dmlist.clear();
	SS->Kmlist.clear();
	SS->LHSlist.clear();
	SS->Jlist.clear();
	SS->Jdlist.clear();
	SS->JTlist.clear();
	SS->LHSGlist.clear();
	SS->Gmlist.clear();
	SS->row = 0;

	// init extra constraints for joint solver
	for (int i = 0; i < LS->constraints.size(); ++i)
	{
		LS->constraints[i]->update(SS, LS, S);
	}
	//if (SS->f.hasNaN())
	//	std::cout << "0 f" << std::endl;

	// Compute fn and M
	SS->alpha = 0;
	int index;
	double value;
	std::shared_ptr<Block> block;
	Matrix6d spcp = Matrix6d::Zero();
	Vector6d bodyForces = Vector6d::Zero();
	for (int i = 0; i < LS->blocks.size(); ++i)
	{
		block = LS->blocks[i].second;
		index = block->joint->jindex;

		for (int r = 0; r < 6; ++r)
		{
			value = S->M[index](r, r) * (1 + SS->alpha*SS->h);
			if (std::abs(value) > THRESHOLD)
				SS->LHSlist.push_back(Eigen::Triplet<double>(index * 6 + r, index * 6 + r, value));
		}

		// Define the spatial cross product matrix
		spcp.block<3, 3>(0, 0) << Rigid::bracket3(S->v.segment<3>(i * 6));
		spcp.block<3, 3>(3, 3) << spcp.block<3, 3>(0, 0);

		// Define body forces
		bodyForces.segment<3>(3) << S->E[i].block<3, 3>(0, 0).transpose()*SS->grav*S->mass[i];

		// Compute force vector				
		SS->f.segment<6>(i * 6) += spcp.transpose()*S->M[i] * S->v.segment<6>(i * 6) + bodyForces;

		// Handle joint stiffness
		// update is in constraintJoint class for CG solver
		int ci = block->joint->constraint_index;
		for (int cnum = 0; cnum < block->joint->constraintNum; ++cnum)
		{
			SS->fr[ci + cnum] -= block->joint->k*(S->q[ci + cnum] - block->joint->q0[cnum]);

			//if(isnan(S->q[ci + cnum]))
			//	std::cout << "q" << std::endl;
		}
	}

	// Load joint constraints
	ConstraintJoint::update(SS, LS, S);

	// solve
	Eigen::SparseMatrix<double> Mr(constraints, constraints);
	Eigen::SparseMatrix<double> Mrtilde(constraints, constraints);
	Eigen::VectorXd frtilde;

	Eigen::SparseMatrix<double> J(numObj * 6, constraints);
	J.setFromTriplets(SS->Jlist.begin(), SS->Jlist.end());
	Eigen::SparseMatrix<double> JT(constraints, numObj * 6);
	JT.setFromTriplets(SS->JTlist.begin(), SS->JTlist.end());
	Eigen::SparseMatrix<double> Jdot(numObj * 6, constraints);
	Jdot.setFromTriplets(SS->Jdlist.begin(), SS->Jdlist.end());

	Eigen::SparseMatrix<double> LHS(numObj * 6, numObj * 6);
	LHS.setFromTriplets(SS->LHSlist.begin(), SS->LHSlist.end());
	Eigen::SparseMatrix<double> Km(numObj * 6, numObj * 6);
	Km.setFromTriplets(SS->Kmlist.begin(), SS->Kmlist.end());
	Eigen::SparseMatrix<double> Dm(numObj * 6, numObj * 6);
	Dm.setFromTriplets(SS->Dmlist.begin(), SS->Dmlist.end());

	Eigen::SparseMatrix<double> KDMm(numObj * 6, numObj * 6);
	KDMm = LHS + SS->h * SS->h * Km + SS->h * Dm + SS->h * SS->bDm;
	Eigen::SparseMatrix<double> KDr(constraints, constraints);
	KDr = SS->h * SS->Dr + SS->h * SS->h * SS->Kr;

	//Mr_sparse = 0.5*(Mr_sparse + Eigen::SparseMatrix<double>(Mr_sparse.transpose()));
	//Mrtilde_sparse = Mr_sparse + KDr_sparse;
	Mr = JT * LHS * J;
	Mrtilde = JT * KDMm * J + KDr;

	SS->fr += JT * (SS->f - LHS * (Jdot * S->qdot));
	frtilde = Mr * S->qdot + SS->h * SS->fr;

	// solve for velocity
	auto start = std::chrono::steady_clock::now();
	if (additional_constraints == 0)
	{
		Eigen::PardisoLDLT< Eigen::SparseMatrix<double>> solver;
		S->qdot = solver.compute(Mrtilde).solve(frtilde);
	}
	else
	{
		Eigen::SparseMatrix<double> Gm(additional_constraints, numObj * 6);
		Gm.setFromTriplets(SS->Gmlist.begin(), SS->Gmlist.end());
		//Eigen::SparseMatrix<double> GmT(numObj * 6, numObj * 6);
		//GmT.setFromTriplets(SS->GmTlist.begin(), SS->GmTlist.end());

		Eigen::SparseMatrix<double> Gr = Gm * J;
		//Eigen::SparseMatrix<double> GrT = JT * GmT;

		std::vector<Eigen::Triplet<double> > tripletList;
		tripletList.reserve(Mrtilde.nonZeros() + Gr.nonZeros() + Gr.nonZeros());
		for (int k = 0; k < Mrtilde.outerSize(); ++k)
		{
			// Mrtilde
			for (Eigen::SparseMatrix<double>::InnerIterator it(Mrtilde, k); it; ++it)
			{
				tripletList.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
			}
		}
		for (int k = 0; k < Gr.outerSize(); ++k)
		{
			for (Eigen::SparseMatrix<double>::InnerIterator it(Gr, k); it; ++it)
			{ 
				// Gr
				tripletList.push_back(Eigen::Triplet<double>(it.row() + Mrtilde.rows(), it.col(), it.value()));
				// GrT
				tripletList.push_back(Eigen::Triplet<double>(it.col(), it.row() + Mrtilde.rows(), it.value()));
			}
		}
		Eigen::SparseMatrix<double> LHSG(constraints + additional_constraints, constraints + additional_constraints);
		LHSG.setFromTriplets(tripletList.begin(), tripletList.end());

		Eigen::VectorXd frG = Eigen::VectorXd::Zero(constraints + additional_constraints);
		frG.segment(0, frtilde.rows()) = frtilde;
		frG.segment(frtilde.rows(), additional_constraints) = -SS->gm * SS->baumgarte[2];

		Eigen::PardisoLU< Eigen::SparseMatrix<double>> solver;
		Eigen::VectorXd result = solver.compute(LHSG).solve(frG);
		S->qdot = result.segment(0,constraints);

		if (false)//std::abs(S->t - 1.01)) < 1e-5)
		{
			// name of file
			std::string decimals = std::to_string((int)(100 * (S->t - (int)std::floor(S->t))));
			if (decimals.size() == 1)
				decimals = "0" + decimals;
			std::string filename = "../matricies/redmax_" + std::to_string(LS->blocks.size()) + "t" + std::to_string((int)std::floor(S->t)) + "_" + decimals + ".m";
			//std::string filename = "../matricies/redmax_" + std::to_string(LS->blocks.size()) + "Preconditioner_noblkdiag_t" + std::to_string((int)std::floor(S->t)) + "_" + decimals + ".m";
			std::ofstream outfile(filename, std::ios::out);
			outfile << std::setprecision(16);

			//ConstraintJoint::preprocess_preconditioner(SS, LS, S, true);
			//Eigen::MatrixXd preconditioner = Eigen::MatrixXd::Zero(constraints, constraints);
			//Eigen::VectorXd temp;
			//Eigen::VectorXd z;
			//for (int k = 0; k < constraints; ++k)
			//{
			//	temp = Eigen::VectorXd::Zero(constraints);
			//	temp[k] = 1.0;
			//	z = ConstraintJoint::computeMinv_x(temp, SS, LS, S);
			//	for (int l = 0; l < constraints; ++l)
			//	{
			//		preconditioner(l, k) = z[l];
			//	}
			//}

			outfile << "h = " << SS->h << ";\n" << std::endl;
			addSparseToFile(outfile, "Mr_c", Mr);
			addSparseToFile(outfile, "Mrtilde_c", Mrtilde);
			addSparseToFile(outfile, "Mm_c", LHS);
			addSparseToFile(outfile, "J_c", J);
			addSparseToFile(outfile, "Jdot_c", Jdot);

			Eigen::SparseMatrix<double> bDm = Dm + SS->bDm;
			addSparseToFile(outfile, "Dm_c", bDm);
			//addSparseToFile(outfile, "MG_c", LHSGt);
			//addSparseToFile(outfile, "Gr_c", Gr);
			addSparseToFile(outfile, "Km_c", Km);
			addSparseToFile(outfile, "LHSG_c", LHSG);
			addSparseToFile(outfile, "Gm_c", Gm);
			addSparseToFile(outfile, "Dr_c", SS->Dr);
			addSparseToFile(outfile, "Kr_c", SS->Kr);
			addSparseToFile(outfile, "Jmr_c", SS->J);
			addSparseToFile(outfile, "Jdotmr_c", SS->Jd);
			addVectorToFile(outfile, "f_c", SS->f);
			addVectorToFile(outfile, "fr_c", SS->fr);
			addVectorToFile(outfile, "frG_c", frG);
			addVectorToFile(outfile, "fm_c", SS->f);
			addVectorToFile(outfile, "gm_c", SS->gm);
			//addVectorToFile(outfile, "frg_c", frGt);
			addVectorToFile(outfile, "frtilde_c", frtilde);
			addVectorToFile(outfile, "rhs_c", frtilde);
			///Eigen::VectorXd qdot0 = S->qdot;
			//addVectorToFile(outfile, "qdot0_c", qdot0);
			addVectorToFile(outfile, "qdot1_c", S->qdot);
			//addSparseToFile(outfile, "P_c", preconditioner);
			
			//SS->fr += JT * (SS->f - LHS * (Jdot * S->qdot));


			outfile.close();
			std::cout << "Printed" << std::endl;
		}

	}
	auto end = std::chrono::steady_clock::now();
	auto diff = end - start;
	auto time = std::chrono::duration<double, std::nano>(diff).count();
	trackLastTimestep.time_in_solve = time;

	// update joint angles
	S->q = S->q + SS->h * S->qdot;

	if (false)//std::abs(S->t - 1.01)) < 1e-5)
	{
		// name of file
		std::string decimals = std::to_string((int)(100 * (S->t - (int)std::floor(S->t))));
		if (decimals.size() == 1)
			decimals = "0" + decimals;
		std::string filename = "../matricies/redmax_" + std::to_string(LS->blocks.size()) + "t" + std::to_string((int)std::floor(S->t)) + "_" + decimals + ".m";
		//std::string filename = "../matricies/redmax_" + std::to_string(LS->blocks.size()) + "Preconditioner_noblkdiag_t" + std::to_string((int)std::floor(S->t)) + "_" + decimals + ".m";
		std::ofstream outfile(filename, std::ios::out);
		outfile << std::setprecision(16);

		//ConstraintJoint::preprocess_preconditioner(SS, LS, S, true);
		//Eigen::MatrixXd preconditioner = Eigen::MatrixXd::Zero(constraints, constraints);
		//Eigen::VectorXd temp;
		//Eigen::VectorXd z;
		//for (int k = 0; k < constraints; ++k)
		//{
		//	temp = Eigen::VectorXd::Zero(constraints);
		//	temp[k] = 1.0;
		//	z = ConstraintJoint::computeMinv_x(temp, SS, LS, S);
		//	for (int l = 0; l < constraints; ++l)
		//	{
		//		preconditioner(l, k) = z[l];
		//	}
		//}

		outfile << "h = " << SS->h << ";\n" << std::endl;
		addSparseToFile(outfile, "Mr_c", Mr);
		addSparseToFile(outfile, "Mrtilde_c", Mrtilde);
		addSparseToFile(outfile, "Mm_c", LHS);

		Eigen::SparseMatrix<double> bDm = Dm + SS->bDm;
		addSparseToFile(outfile, "Dm_c", bDm);
		//addSparseToFile(outfile, "MG_c", LHSGt);
		//addSparseToFile(outfile, "Gr_c", Gr);
		addSparseToFile(outfile, "Km_c", Km);
		addSparseToFile(outfile, "Dr_c", SS->Dr);
		addSparseToFile(outfile, "Kr_c", SS->Kr);
		addSparseToFile(outfile, "Jmr_c", SS->J);
		addSparseToFile(outfile, "Jdotmr_c", SS->Jd);
		addVectorToFile(outfile, "fr_c", SS->fr);
		addVectorToFile(outfile, "fm_c", SS->f);
		//addVectorToFile(outfile, "frg_c", frGt);
		addVectorToFile(outfile, "frtilde_c", frtilde);
		addVectorToFile(outfile, "rhs_c", frtilde);
		///Eigen::VectorXd qdot0 = S->qdot;
		//addVectorToFile(outfile, "qdot0_c", qdot0);
		addVectorToFile(outfile, "qdot1_c", S->qdot);
		//addSparseToFile(outfile, "P_c", preconditioner);

		outfile.close();
		std::cout << "Printed" << std::endl;
	}

	if (false)
	{
		std::string decimals = std::to_string((int)(100 * (S->t - (int)std::floor(S->t))));
		if (decimals.size() == 1)
			decimals = "0" + decimals;
		std::string filename = "../matricies/qdot" + std::to_string((int)std::floor(S->t)) + "_" + decimals + ".m";
		std::ofstream outfile(filename, std::ios::out);
		outfile << std::setprecision(16);
		
		addSparseToFile(outfile, "Mrtilde_c", Mrtilde);
		addSparseToFile(outfile, "Mm_c", LHS);
		addSparseToFile(outfile, "Km_c", Km);
		addSparseToFile(outfile, "Jmr_c", SS->J);
		addVectorToFile(outfile, "qdot1_c", S->qdot);

		outfile.close();

		std::cout << "Printed" << std::endl;
	}

	for (int i = 0; i < (int)LS->blocks.size(); ++i)
	{
		// update pos and vel
		LS->joints[LS->joint_map[i]]->update(LS, S);
	}

	//if (LS->blocks.size() == 603 && count < 3)
	//{
	//	count++;
	//	for (int i = 0; i < S->qdot.size(); ++i)
	//	{
	//		std::cout << S->qdot[i] << " ";
	//	}
	//	std::cout << std::endl;
	//}

	return S->qdot;
}
#endif

Eigen::VectorXd Solver::solvePCG(std::unique_ptr<StateSolve>& SS, std::unique_ptr<LinkageSystem>& LS, std::unique_ptr<State>& S)
{
	// Preprocessing, eliminate matrix resizing
	int constraints = SS->Mr_dimension;
	int additional_constraints = SS->Gr_dimension;
	int numObj = (int)LS->blocks.size();

	// Initialize containers
	SS->f = Eigen::VectorXd::Zero(numObj * 6); // load with fm
	SS->fr = Eigen::VectorXd::Zero(constraints);
	SS->Gm = Eigen::MatrixXd::Zero(additional_constraints, numObj * 6);
	SS->GmT = Eigen::MatrixXd::Zero(numObj * 6, additional_constraints);
	SS->gm = Eigen::VectorXd::Zero(additional_constraints);

	// Precomputations for SoA
	ConstraintJoint::preprocess_PCG(SS, LS, S);

	// load constraint RHS update and CG joint preprocessing
	for (int i = 0; i < LS->constraints.size(); ++i)
	{
		LS->constraints[i]->updateJoint(SS, LS, S);
	}
	// Update maximal body forces (RHS)
	computeRHS(SS, LS, S);

	// preprocessing
	ConstraintJoint::preprocess_PCG_preconditioner(SS, LS, S, true);

	/// preconditioned
	//Eigen::VectorXd z;
	//Eigen::VectorXd p;
	//if (simtype != redCGNoMat_noprec)
	//{
	//	if (simtype == redCGNoMat_noblkdiag)
	//		ConstraintJoint::preprocess_PCG_preconditioner(SS, LS, S, false);
	//	else
	//		ConstraintJoint::preprocess_PCG_preconditioner(SS, LS, S, true);

	//	z = ConstraintJoint::computeMinv_x(r, SS, LS, S);
	//	p = z;
	//}
	//// NOT preconditioned 
	//else
	//{
	//	p = r;
	//}

	// solve
	Eigen::VectorXd qd = S->qdot;		// initial guess is prev. sol.
	Eigen::VectorXd J_x;
	Eigen::VectorXd Jdot_x;

	// J and Jdot in parallel
	ConstraintJoint::computeJ_Jdot_x(J_x, Jdot_x, qd, SS, LS, S);

	Eigen::VectorXd LHSqd = ConstraintJoint::computeJT_x(
		ConstraintJoint::computeLHS_x(J_x, SS, LS, S),
		SS, LS, S);

	// Joint stiffness and damping
	ConstraintJoint::computeStiffnessDampingJoint(LHSqd, qd, S->q, SS, LS, S);

	Eigen::VectorXd Mqd = ConstraintJoint::computeJT_x(
		ConstraintJoint::computeM_x(J_x, SS, LS, S),
		SS, LS, S);
	Eigen::VectorXd fr_save = SS->fr;
	SS->fr += ConstraintJoint::computeJT_x(
		(SS->f - ConstraintJoint::computeM_x(Jdot_x, SS, LS, S)),
		SS, LS, S);
	Eigen::VectorXd frtilde = Mqd + SS->h * SS->fr;

	//////////////////////////////////////////////////////////////////////////
	// Implement Shin's refactoring so JT is only applied once!
	//////////////////////////////////////////////////////////////////////////
	// error accumulation from somewhere
	//// J and Jdot in parallel
	//ConstraintJoint::computeJ_Jdot_x(J_x, Jdot_x, qd, SS, LS, S);
	//Eigen::VectorXd LHSJqd = ConstraintJoint::computeLHS_x(J_x, SS, LS, S);
	//// Joint stiffness and damping
	//ConstraintJoint::computeStiffnessDampingJoint(LHSJqd, qd, S->q, SS, LS, S);
	//Eigen::VectorXd MJqd = ConstraintJoint::computeM_x(J_x, SS, LS, S);
	//Eigen::VectorXd f_MJdqd = SS->f - ConstraintJoint::computeM_x(Jdot_x, SS, LS, S);
	//Eigen::VectorXd Mrtilde_qdot0;
	//Eigen::VectorXd Mr_qdot0;
	//ConstraintJoint::computeJT_x_parallel(Mrtilde_qdot0, LHSJqd, Mr_qdot0, MJqd, SS->fr, f_MJdqd, SS, LS, S);
	//Eigen::VectorXd frtilde = Mr_qdot0 + SS->h * SS->fr;

	if (additional_constraints == 0)
	{
		// Mrtilde\frtilde = qdot
		pcdSaad2003(S->qdot, LHSqd, frtilde, SS, LS, S);
	}
	else if (additional_constraints == 2)
	{
		// Bridge scene has one loop closing constraint
		Eigen::MatrixXd GrT = Eigen::MatrixXd::Zero(constraints, 2);

		Eigen::VectorXd Gmrow1 = SS->Gm.row(0);
		Eigen::VectorXd Gmrow2 = SS->Gm.row(1);
		Eigen::VectorXd GrTcol1 = ConstraintJoint::computeJT_x(Gmrow1, SS, LS, S);
		Eigen::VectorXd GrTcol2 = ConstraintJoint::computeJT_x(Gmrow2, SS, LS, S);
        
        // Run pcg for each row of G
        Eigen::VectorXd MiGt1;
		Eigen::VectorXd MiGt2;
		pcdSaad2003(MiGt1, LHSqd, GrTcol1, SS, LS, S);
		pcdSaad2003(MiGt2, LHSqd, GrTcol2, SS, LS, S);

		// Form LHS and RHS
		Eigen::VectorXd qdot1unc;
		pcdSaad2003(qdot1unc, LHSqd, frtilde, SS, LS, S);
		Eigen::VectorXd RHS = SS->baumgarte[2] * SS->gm;
		Eigen::MatrixXd LHS = Eigen::MatrixXd::Zero(2, 2);
		for (int r = 0; r < constraints; ++r)
		{
			LHS(0, 0) += GrTcol1[r] * MiGt1[r];
			LHS(0, 1) += GrTcol1[r] * MiGt2[r];
			LHS(1, 0) += GrTcol2[r] * MiGt1[r];
			LHS(1, 1) += GrTcol2[r] * MiGt2[r];

			RHS[0] += GrTcol1[r] * qdot1unc[r];
			RHS[1] += GrTcol2[r] * qdot1unc[r];
		}

		// Solve 2-by-2 system for lambda
		Eigen::VectorXd lambda = LHS.ldlt().solve(RHS);

		// Solve!
		Eigen::VectorXd frtilde_GTlambda = frtilde;
		for (int r = 0; r < constraints; ++r)
		{
			frtilde_GTlambda[r] -= GrTcol1[r] * lambda[0] + GrTcol2[r] * lambda[1];
		}
		pcdSaad2003(S->qdot, LHSqd, frtilde_GTlambda, SS, LS, S);
	}
	else
	{
		// Handle larger scenes
		Eigen::SparseMatrix<double> Gr = Eigen::SparseMatrix<double>(additional_constraints, constraints);
		Eigen::SparseMatrix<double> GrT = Eigen::SparseMatrix<double>(constraints, additional_constraints);
		Eigen::SparseMatrix<double> MiGt = Eigen::SparseMatrix<double>(constraints, additional_constraints);
		std::vector< Eigen::Triplet<double> > Grlist;
		std::vector< Eigen::Triplet<double> > GrTlist;
		std::vector< Eigen::Triplet<double> > MiGtlist;

		// Run PCG for each row of G
		int num_joints = (int)LS->joints.size();
#pragma omp parallel for //private(Gmrowx, GrTcolx, MiGtx)
		for (int i = 0; i < additional_constraints; ++i)
		{
			Eigen::VectorXd Gmrowx;
			Eigen::VectorXd GrTcolx;
			Eigen::VectorXd MiGtx;
			std::shared_ptr<State::local_mt> lmt = std::make_shared<State::local_mt>();
			lmt->alpha_.resize(num_joints);
			lmt->Bhat_.resize(num_joints);
			lmt->beta_.resize(num_joints);
			lmt->Vdot_.resize(num_joints);
			lmt->ST_Bhat_.resize(num_joints);

			Gmrowx = SS->Gm.row(i);
			GrTcolx = ConstraintJoint::computeJT_x(Gmrowx, SS, LS, S, lmt);
			pcdSaad2003(MiGtx, LHSqd, GrTcolx, SS, LS, S, lmt);

			for (int j = 0; j < constraints; ++j)
			{
				if (std::abs(GrTcolx[j]) > THRESHOLD)
				{
#pragma omp critical
					Grlist.push_back(Eigen::Triplet<double>(i, j, GrTcolx[j]));
#pragma omp critical
					GrTlist.push_back(Eigen::Triplet<double>(j, i, GrTcolx[j]));
				}
				if (std::abs(MiGtx[j]) > THRESHOLD)
				{
#pragma omp critical
					MiGtlist.push_back(Eigen::Triplet<double>(j, i, MiGtx[j]));
				}
			}
		}

		// Form LHS and RHS
		Eigen::VectorXd qdot1unc;
		pcdSaad2003(qdot1unc, LHSqd, frtilde, SS, LS, S);
		Eigen::VectorXd RHS = SS->baumgarte[2] * SS->gm;
		Eigen::SparseMatrix<double> LHS = Eigen::SparseMatrix<double>(additional_constraints, additional_constraints);

		Gr.setFromTriplets(Grlist.begin(), Grlist.end());
		MiGt.setFromTriplets(MiGtlist.begin(), MiGtlist.end());

		LHS = Gr * MiGt;
		RHS += Gr * qdot1unc;

		// Solve x-by-x system for lambda
#ifdef REDMAX_PARDISO
		Eigen::PardisoLDLT< Eigen::SparseMatrix<double>> solver;
        Eigen::VectorXd lambda = solver.compute(LHS).solve(RHS);
#else
        Eigen::SparseLU< Eigen::SparseMatrix<double>> solver;
        solver.compute(LHS);
        Eigen::VectorXd lambda = solver.solve(RHS);
#endif

		// Solve!
		GrT.setFromTriplets(GrTlist.begin(), GrTlist.end());
		frtilde -= GrT * lambda;
		pcdSaad2003(S->qdot, LHSqd, frtilde, SS, LS, S);


		if (false)//std::abs(S->t - 1.01)) < 1e-5)
		{
			std::string decimals = std::to_string((int)(100 * (S->t - (int)std::floor(S->t))));
			if (decimals.size() == 1)
				decimals = "0" + decimals;
			std::string filename = "../matricies/redcgSoA_" + std::to_string(LS->blocks.size()) + "links_t" + std::to_string((int)std::floor(S->t)) + "_" + decimals + ".m";
			//std::string filename = "../matricies/redmax_" + std::to_string(LS->blocks.size()) + "Preconditioner_noblkdiag_t" + std::to_string((int)std::floor(S->t)) + "_" + decimals + ".m";
			std::ofstream outfile(filename, std::ios::out);
			outfile << std::setprecision(16);

			outfile << "h = " << SS->h << ";\n" << std::endl;
			//addVectorToFile(outfile, "MiGt1_c", MiGt1);
			//addVectorToFile(outfile, "MiGt2_c", MiGt2);
			//addVectorToFile(outfile, "GrTcol1_c", GrTcol1);
			//addVectorToFile(outfile, "GrTcol2_c", GrTcol2);
			addSparseToFile(outfile, "MiGt_p", MiGt);
			addVectorToFile(outfile, "Mrtildeqd_p", LHSqd);
			addSparseToFile(outfile, "LHS_p", LHS);
			addSparseToFile(outfile, "Gr_p", Gr);
			addSparseToFile(outfile, "Gm_p", SS->Gm);
			addSparseToFile(outfile, "GmT_p", SS->GmT);
			addVectorToFile(outfile, "gm_p", SS->gm);
			addVectorToFile(outfile, "qdot1unc_p", qdot1unc);
			addVectorToFile(outfile, "RHS_p", RHS);
			addVectorToFile(outfile, "fr_p", SS->fr);
			addVectorToFile(outfile, "lambda_p", lambda);
			addVectorToFile(outfile, "GTlambda_p", frtilde);
			addVectorToFile(outfile, "fm_p", SS->f);
			addVectorToFile(outfile, "qdot1_p", S->qdot);
			//addVectorToFile(outfile, "oJT_p", onesJT);
			//addVectorToFile(outfile, "oJ_p", onesJ);

			outfile.close();
			std::cout << "Printed" << std::endl;
		}
	}

	// update joint angles
	S->q = S->q + SS->h * S->qdot;

	if (false)
	{
		std::string decimals = std::to_string((int)(100 * (S->t - (int)std::floor(S->t))));
		if (decimals.size() == 1)
			decimals = "0" + decimals;
		std::string filename = "../matricies/pcgSoA_qdot" + std::to_string((int)std::floor(S->t)) + "_" + decimals + ".m";
		std::ofstream outfile(filename, std::ios::out);
		outfile << std::setprecision(16);
		//outfile << "qdot1_p" << " = [\n";
		//for (int i = 0; i < (int)LS->blocks.size(); ++i)
		//{
		//	Eigen::VectorXd tt = LS->joints[LS->joint_map[i]]->qdot;
		//	for (int j = 0; j < tt.size(); j++)
		//	{
		//		outfile << tt[j] << " ";
		//	}
		//}
		//outfile << "]';\n" << std::endl;

		//Eigen::VectorXd rj = frtilde_save - LHSqd;
		//Eigen::VectorXd zj = ConstraintJoint::computeMinv_x(rj, SS, LS, S);
		//Eigen::VectorXd check = ConstraintJoint::computeMinv_x(rj, SS, LS, S);

		Eigen::VectorXd qdc = Eigen::VectorXd::Zero(S->qdot.size());
		for (int i = 0; i < qdc.size(); ++i)
			qdc[i] = 1;// (((double)rand() / (RAND_MAX)) - 0.5) * 4;
		Eigen::VectorXd checkc = ConstraintJoint::computeMinv_x(qdc, SS, LS, S);

		Eigen::VectorXd MJdqd = ConstraintJoint::computeM_x(Jdot_x, SS, LS, S);

		Eigen::VectorXd J_ones;
		Eigen::VectorXd Jdot_ones;
		ConstraintJoint::computeJ_Jdot_x(J_ones, Jdot_ones, qdc, SS, LS, S);

		addVectorToFile(outfile, "Mrtildeqd_l", LHSqd);
		addVectorToFile(outfile, "frtilde_l", frtilde);
		addVectorToFile(outfile, "Mqd_l", Mqd);
		addVectorToFile(outfile, "fr_l", fr_save);
		addVectorToFile(outfile, "fmr_l", SS->fr);
		addVectorToFile(outfile, "fm_l", SS->f);
		//addVectorToFile(outfile, "zj_l", zj);
		addVectorToFile(outfile, "qdot1_l", S->qdot);
		addVectorToFile(outfile, "check_l", checkc);
		addVectorToFile(outfile, "MJdqd_l", MJdqd);
		addVectorToFile(outfile, "Jdqd_l", Jdot_x);
		addVectorToFile(outfile, "Jqd_l", J_x);
		addVectorToFile(outfile, "Jones_l", J_ones);
		addVectorToFile(outfile, "Jdones_l", Jdot_ones);

		outfile.close();
		std::cout << "Printed" << std::endl;
	}

	for (int i = 0; i < (int)LS->blocks.size(); ++i)
	{
		// update pos and vel
		LS->joints[LS->joint_map[i]]->update(LS, S);
	}
	return S->qdot;
}

Eigen::VectorXd Solver::solvePCG_unopt(std::unique_ptr<StateSolve>& SS, std::unique_ptr<LinkageSystem>& LS, std::unique_ptr<State>& S)
{
	// Preprocessing, eliminate matrix resizing
	int constraints = SS->Mr_dimension;
	int additional_constraints = SS->Gr_dimension;
	int numObj = (int)LS->blocks.size();

	// Initialize containers
	SS->f = Eigen::VectorXd::Zero(numObj * 6); // load with fm
	SS->fr = Eigen::VectorXd::Zero(constraints);
	SS->Gm = Eigen::MatrixXd::Zero(additional_constraints, numObj * 6);
	SS->GmT = Eigen::MatrixXd::Zero(numObj * 6, additional_constraints);
	SS->gm = Eigen::VectorXd::Zero(additional_constraints);

	// load constraint RHS update and CG joint preprocessing
	for (int i = 0; i < LS->constraints.size(); ++i)
	{
		LS->constraints[i]->updateJoint(SS, LS, S);
	}
	// Update maximal body forces (RHS)
	computeRHS(SS, LS, S);

	// preprocessing
	ConstraintJoint::preprocess_preconditioner_unopt(SS, LS, S, true);

	//if (simtype == redCGNoMat_noblkdiag)
	//	ConstraintJoint::preprocess_preconditioner_unopt(SS, LS, S, false);
	//else
	//	ConstraintJoint::preprocess_preconditioner_unopt(SS, LS, S, true);


	// solve
	Eigen::VectorXd qd = S->qdot;		// initial guess is prev. sol.
	Eigen::VectorXd J_x;
	Eigen::VectorXd Jdot_x;

	// J and Jdot in parallel
	ConstraintJoint::computeJ_Jdot_x(J_x, Jdot_x, qd, SS, LS, S);

	Eigen::VectorXd LHSqd = ConstraintJoint::computeJT_x(
		ConstraintJoint::computeLHS_x(J_x, SS, LS, S),
		SS, LS, S);
	// Joint stiffness and damping
	ConstraintJoint::computeStiffnessDampingJoint(LHSqd, qd, S->q, SS, LS, S);

	Eigen::VectorXd Mqd = ConstraintJoint::computeJT_x(
		ConstraintJoint::computeM_x(J_x, SS, LS, S),
		SS, LS, S);
	SS->fr += ConstraintJoint::computeJT_x(
		(SS->f - ConstraintJoint::computeM_x(Jdot_x, SS, LS, S)),
		SS, LS, S);
	Eigen::VectorXd frtilde = Mqd + SS->h * SS->fr;
	Eigen::VectorXd frtilde_save = frtilde;

	//////////////////////////////////////////////////////////////////////////
	// Implement Shin's refactoring so JT is only applied once!
	//////////////////////////////////////////////////////////////////////////

	// error accumulation from somewhere
	//// J and Jdot in parallel
	//ConstraintJoint::computeJ_Jdot_x(J_x, Jdot_x, qd, SS, LS, S);
	//Eigen::VectorXd LHSJqd = ConstraintJoint::computeLHS_x(J_x, SS, LS, S);
	//// Joint stiffness and damping
	//ConstraintJoint::computeStiffnessDampingJoint(LHSJqd, qd, S->q, SS, LS, S);
	//Eigen::VectorXd MJqd = ConstraintJoint::computeM_x(J_x, SS, LS, S);
	//Eigen::VectorXd f_MJdqd = SS->f - ConstraintJoint::computeM_x(Jdot_x, SS, LS, S);
	//Eigen::VectorXd Mrtilde_qdot0;
	//Eigen::VectorXd Mr_qdot0;
	//ConstraintJoint::computeJT_x_parallel(Mrtilde_qdot0, LHSJqd, Mr_qdot0, MJqd, SS->fr, f_MJdqd, SS, LS, S);
	//Eigen::VectorXd frtilde = Mr_qdot0 + SS->h * SS->fr;

	if (additional_constraints == 0)
	{
		// Mrtilde\frtilde = qdot
		pcdSaad2003(S->qdot, LHSqd, frtilde, SS, LS, S);
	}
	else if (additional_constraints == 2)
	{
		// Bridge scene has one loop closing constraint
		Eigen::MatrixXd GrT = Eigen::MatrixXd::Zero(constraints, 2);

		Eigen::VectorXd Gmrow1 = SS->Gm.row(0);
		Eigen::VectorXd Gmrow2 = SS->Gm.row(1);
		Eigen::VectorXd GrTcol1 = ConstraintJoint::computeJT_x(Gmrow1, SS, LS, S);
		Eigen::VectorXd GrTcol2 = ConstraintJoint::computeJT_x(Gmrow2, SS, LS, S);

		// Run pcg for each row of G
		Eigen::VectorXd MiGt1;
		Eigen::VectorXd MiGt2;
		pcdSaad2003(MiGt1, LHSqd, GrTcol1, SS, LS, S);
		pcdSaad2003(MiGt2, LHSqd, GrTcol2, SS, LS, S);

		// Form LHS and RHS
		Eigen::VectorXd qdot1unc;
		pcdSaad2003(qdot1unc, LHSqd, frtilde, SS, LS, S);
		Eigen::VectorXd RHS = SS->baumgarte[2] * SS->gm;
		Eigen::MatrixXd LHS = Eigen::MatrixXd::Zero(2, 2);
		for (int r = 0; r < constraints; ++r)
		{
			LHS(0, 0) += GrTcol1[r] * MiGt1[r];
			LHS(0, 1) += GrTcol1[r] * MiGt2[r];
			LHS(1, 0) += GrTcol2[r] * MiGt1[r];
			LHS(1, 1) += GrTcol2[r] * MiGt2[r];

			RHS[0] += GrTcol1[r] * qdot1unc[r];
			RHS[1] += GrTcol2[r] * qdot1unc[r];
		}

		// Solve 2-by-2 system for lambda
		Eigen::VectorXd lambda = LHS.ldlt().solve(RHS);

		// Solve!
		Eigen::VectorXd frtilde_GTlambda = frtilde;
		for (int r = 0; r < constraints; ++r)
		{
			frtilde_GTlambda[r] -= GrTcol1[r] * lambda[0] + GrTcol2[r] * lambda[1];
		}
		pcdSaad2003(S->qdot, LHSqd, frtilde_GTlambda, SS, LS, S);

		if (false)//std::abs(S->t - 1.01)) < 1e-5)
		{
			std::string decimals = std::to_string((int)(100 * (S->t - (int)std::floor(S->t))));
			if (decimals.size() == 1)
				decimals = "0" + decimals;
			std::string filename = "../matricies/redcg_" + std::to_string(LS->blocks.size()) + "links_2case_t" + std::to_string((int)std::floor(S->t)) + "_" + decimals + ".m";
			//std::string filename = "../matricies/redmax_" + std::to_string(LS->blocks.size()) + "Preconditioner_noblkdiag_t" + std::to_string((int)std::floor(S->t)) + "_" + decimals + ".m";
			std::ofstream outfile(filename, std::ios::out);
			outfile << std::setprecision(16);

			outfile << "h = " << SS->h << ";\n" << std::endl;
			addVectorToFile(outfile, "MiGt1_c", MiGt1);
			addVectorToFile(outfile, "MiGt2_c", MiGt2);
			addVectorToFile(outfile, "GrTcol1_c", GrTcol1);
			addVectorToFile(outfile, "GrTcol2_c", GrTcol2);
			addSparseToFile(outfile, "LHS_c", LHS);
			addSparseToFile(outfile, "GrT_c", GrT);
			addSparseToFile(outfile, "Gm_c", SS->Gm);
			addSparseToFile(outfile, "GmT_c", SS->GmT);
			addVectorToFile(outfile, "gm_c", SS->gm);
			addVectorToFile(outfile, "qdot1unc_c", qdot1unc);
			addVectorToFile(outfile, "RHS_c", RHS);
			addVectorToFile(outfile, "fr_c", SS->fr);
			addVectorToFile(outfile, "lambda_c", lambda);
			addVectorToFile(outfile, "GTlambda_c", frtilde_GTlambda);
			addVectorToFile(outfile, "fm_c", SS->f);
			addVectorToFile(outfile, "qdot1_c", S->qdot);

			outfile.close();
			std::cout << "Printed" << std::endl;
		}
	}
	else
	{
		// Handle larger scenes
		Eigen::SparseMatrix<double> Gr = Eigen::SparseMatrix<double>(additional_constraints, constraints);
		Eigen::SparseMatrix<double> GrT = Eigen::SparseMatrix<double>(constraints, additional_constraints);
		Eigen::SparseMatrix<double> MiGt = Eigen::SparseMatrix<double>(constraints, additional_constraints);
		std::vector< Eigen::Triplet<double> > Grlist;
		std::vector< Eigen::Triplet<double> > GrTlist;
		std::vector< Eigen::Triplet<double> > MiGtlist;

		// Run PCG for each row of G
		int num_joints = LS->joints.size();
#pragma omp parallel for //private(Gmrowx, GrTcolx, MiGtx)
		for (int i = 0; i < additional_constraints; ++i)
		{
			Eigen::VectorXd Gmrowx;
			Eigen::VectorXd GrTcolx;
			Eigen::VectorXd MiGtx;
			std::shared_ptr<State::local_mt> lmt = std::make_shared<State::local_mt>();
			lmt->alpha_.resize(num_joints);
			lmt->Bhat_.resize(num_joints);
			lmt->beta_.resize(num_joints);
			lmt->Vdot_.resize(num_joints);
			lmt->ST_Bhat_.resize(num_joints);

			Gmrowx = SS->Gm.row(i);
			GrTcolx = ConstraintJoint::computeJT_x(Gmrowx, SS, LS, S, lmt);
			pcdSaad2003(MiGtx, LHSqd, GrTcolx, SS, LS, S, lmt);

			for (int j = 0; j < constraints; ++j)
			{
				if (std::abs(GrTcolx[j]) > THRESHOLD)
				{
#pragma omp critical
					Grlist.push_back(Eigen::Triplet<double>(i, j, GrTcolx[j]));
#pragma omp critical
					GrTlist.push_back(Eigen::Triplet<double>(j, i, GrTcolx[j]));
				}
				if (std::abs(MiGtx[j]) > THRESHOLD)
				{
#pragma omp critical
					MiGtlist.push_back(Eigen::Triplet<double>(j, i, MiGtx[j]));
				}
			}
		}

		// Form LHS and RHS
		Eigen::VectorXd qdot1unc;
		pcdSaad2003(qdot1unc, LHSqd, frtilde, SS, LS, S);
		Eigen::VectorXd RHS = SS->baumgarte[2] * SS->gm;
		Eigen::SparseMatrix<double> LHS = Eigen::SparseMatrix<double>(additional_constraints, additional_constraints);

		Gr.setFromTriplets(Grlist.begin(), Grlist.end());
		MiGt.setFromTriplets(MiGtlist.begin(), MiGtlist.end());

		LHS = Gr * MiGt;
		RHS += Gr * qdot1unc;

		// Solve x-by-x system for lambda
#ifdef REDMAX_PARDISO
        Eigen::PardisoLDLT< Eigen::SparseMatrix<double>> solver;
        Eigen::VectorXd lambda = solver.compute(LHS).solve(RHS);
#else
        Eigen::SparseLU< Eigen::SparseMatrix<double>> solver;
        solver.compute(LHS);
        Eigen::VectorXd lambda = solver.solve(RHS);
#endif
		// Solve!
		GrT.setFromTriplets(GrTlist.begin(), GrTlist.end());
		frtilde -= GrT * lambda;
		pcdSaad2003(S->qdot, LHSqd, frtilde, SS, LS, S);


		if (false)//std::abs(S->t - 1.01)) < 1e-5)
		{
			std::string decimals = std::to_string((int)(100 * (S->t - (int)std::floor(S->t))));
			if (decimals.size() == 1)
				decimals = "0" + decimals;
			std::string filename = "../matricies/redcg_" + std::to_string(LS->blocks.size()) + "links_t" + std::to_string((int)std::floor(S->t)) + "_" + decimals + ".m";
			//std::string filename = "../matricies/redmax_" + std::to_string(LS->blocks.size()) + "Preconditioner_noblkdiag_t" + std::to_string((int)std::floor(S->t)) + "_" + decimals + ".m";
			std::ofstream outfile(filename, std::ios::out);
			outfile << std::setprecision(16);

			outfile << "h = " << SS->h << ";\n" << std::endl;
			//addVectorToFile(outfile, "MiGt1_c", MiGt1);
			//addVectorToFile(outfile, "MiGt2_c", MiGt2);
			//addVectorToFile(outfile, "GrTcol1_c", GrTcol1);
			//addVectorToFile(outfile, "GrTcol2_c", GrTcol2);
			addSparseToFile(outfile, "MiGt_l", MiGt);
			addVectorToFile(outfile, "Mrtildeqd_l", LHSqd);
			addSparseToFile(outfile, "LHS_l", LHS);
			addSparseToFile(outfile, "Gr_l", Gr);
			addSparseToFile(outfile, "Gm_l", SS->Gm);
			addSparseToFile(outfile, "GmT_l", SS->GmT);
			addVectorToFile(outfile, "gm_l", SS->gm);
			addVectorToFile(outfile, "qdot1unc_l", qdot1unc);
			addVectorToFile(outfile, "RHS_l", RHS);
			addVectorToFile(outfile, "fr_l", SS->fr);
			addVectorToFile(outfile, "lambda_l", lambda);
			addVectorToFile(outfile, "GTlambda_l", frtilde);
			addVectorToFile(outfile, "frtilde_l", frtilde);
			addVectorToFile(outfile, "fm_l", SS->f);
			addVectorToFile(outfile, "qdot1_l", S->qdot);
			//addVectorToFile(outfile, "oJT_p", onesJT);
			//addVectorToFile(outfile, "oJ_p", onesJ);

			outfile.close();
			std::cout << "Printed" << std::endl;
		}

		if (false)//std::abs(S->t - 1.01)) < 1e-5)
		{
			std::string decimals = std::to_string((int)(100 * (S->t - (int)std::floor(S->t))));
			if (decimals.size() == 1)
				decimals = "0" + decimals;
			std::string filename = "../matricies/redcg_" + std::to_string(LS->blocks.size()) + "links_t" + std::to_string((int)std::floor(S->t)) + "_" + decimals + ".m";
			//std::string filename = "../matricies/redmax_" + std::to_string(LS->blocks.size()) + "Preconditioner_noblkdiag_t" + std::to_string((int)std::floor(S->t)) + "_" + decimals + ".m";
			std::ofstream outfile(filename, std::ios::out);
			outfile << std::setprecision(16);

			outfile << "h = " << SS->h << ";\n" << std::endl;
			//addVectorToFile(outfile, "MiGt1_c", MiGt1);
			//addVectorToFile(outfile, "MiGt2_c", MiGt2);
			//addVectorToFile(outfile, "GrTcol1_c", GrTcol1);
			//addVectorToFile(outfile, "GrTcol2_c", GrTcol2);
			addSparseToFile(outfile, "LHS_c", LHS);
			addSparseToFile(outfile, "GrT_c", GrT);
			//addSparseToFile(outfile, "Gm_c", SS->Gm);
			addSparseToFile(outfile, "GmT_c", SS->GmT);
			addVectorToFile(outfile, "gm_c", SS->gm);
			addVectorToFile(outfile, "qdot1unc_c", qdot1unc);
			addVectorToFile(outfile, "RHS_c", RHS);
			addVectorToFile(outfile, "fr_c", SS->fr);
			addVectorToFile(outfile, "lambda_c", lambda);
			addVectorToFile(outfile, "GTlambda_c", frtilde);
			addVectorToFile(outfile, "fm_c", SS->f);
			addVectorToFile(outfile, "qdot1_c", S->qdot);

			outfile.close();
			std::cout << "Printed" << std::endl;
		}
	}

	if (false)
	{
		std::string decimals = std::to_string((int)(100 * (S->t - (int)std::floor(S->t))));
		if (decimals.size() == 1)
			decimals = "0" + decimals;
		std::string filename = "../matricies/pcg_qdot" + std::to_string((int)std::floor(S->t)) + "_" + decimals + ".m";
		std::ofstream outfile(filename, std::ios::out);
		outfile << std::setprecision(16);

		Eigen::VectorXd rj = frtilde_save - LHSqd;
		Eigen::VectorXd zj = ConstraintJoint::computeMinv_x(rj, SS, LS, S);
		Eigen::VectorXd check = ConstraintJoint::computeMinv_x(rj, SS, LS, S);

		Eigen::VectorXd qdc = Eigen::VectorXd::Zero(S->qdot.size());
		for (int i = 0; i < qdc.size(); ++i)
			qdc[i] = 1;// (((double)rand() / (RAND_MAX)) - 0.5) * 4;
		Eigen::VectorXd checkc = ConstraintJoint::computeMinv_x(qdc, SS, LS, S);

		addVectorToFile(outfile, "Mrtildeqd_l", LHSqd);
		addVectorToFile(outfile, "frtilde_l", frtilde_save);
		addVectorToFile(outfile, "Mqd_l", Mqd);
		addVectorToFile(outfile, "fr_l", SS->fr);
		addVectorToFile(outfile, "zj_l", zj);
		addVectorToFile(outfile, "qdot1_l", S->qdot);
		addVectorToFile(outfile, "check_l", checkc);

		outfile.close();
		std::cout << "Printed" << std::endl;
	}

	// update joint angles
	S->q = S->q + SS->h * S->qdot;

	for (int i = 0; i < (int)LS->blocks.size(); ++i)
	{
		// update pos and vel
		LS->joints[LS->joint_map[i]]->update(LS, S);
	}
	return S->qdot;
}

Solver::Solver()
{
	trackLastTimestep = SolverDataTracker();
}

Eigen::VectorXd Solver::solve(std::unique_ptr<StateSolve> &SS, std::unique_ptr<LinkageSystem> &LS, std::unique_ptr<State> &S)
{
	Eigen::VectorXd result;

	trackLastTimestep.num_iterations = 0;
	if (simtype == simType::PCG)
		result = solvePCG(SS, LS, S);
	else if (simtype == simType::PCG_unopt)
		result = solvePCG_unopt(SS, LS, S);
#ifdef REDMAX_PARDISO
	else if (simtype == simType::Pardiso)
		result = solvePardiso(SS, LS, S);
#endif
	else
		std::cerr << "Specified simtype solver has not been set up" << std::endl;

	return result;
}

void Solver::printDenseToFile(Eigen::MatrixXd & MG, Eigen::VectorXd & f, Eigen::VectorXd & res, std::string id)
{
	std::string filename = "../matricies/" + id + ".m";
	std::ofstream outfile(filename, std::ios::out);

	outfile << "MG = [\n";
	for (int r = 0; r < MG.rows(); r++)
	{
		for (int c = 0; c < MG.cols(); c++)
		{
			outfile << MG(r, c) << " ";
		}
		outfile << ";\n";
	}
	outfile << "];\n\nf = [\n";

	for (int i = 0; i < f.size(); i++)
	{
		outfile << f[i] << "; ";
	}
	outfile << "];\n\nres = [";

	for (int i = 0; i < res.size(); i++)
	{
		outfile << res[i] << "; ";
	}
	outfile << "];\n";
}

void Solver::printSparseToFile(std::vector<Eigen::Triplet<double>>& MG, Eigen::VectorXd & f, Eigen::VectorXd & res)
{
	std::string filename = "../matricies/r_sparse.m";
	std::ofstream outfile(filename, std::ios::out);

	outfile << "i = [\n";
	for (int i = 0; i < MG.size(); i++)
	{
		outfile << MG[i].row() + 1 << " ";
	}
	outfile << "]';\n\nj= [\n";
	for (int i = 0; i < MG.size(); i++)
	{
		outfile << MG[i].col() + 1 << " ";
	}
	outfile << "]';\n\nv = [\n";
	for (int i = 0; i < MG.size(); i++)
	{
		outfile << MG[i].value() << " ";
	}
	outfile << "]';\n" << std::endl;

	outfile << "MG = sparse(i,j,v);\n" << std::endl;

	outfile << "f = [\n";
	for (int i = 0; i < f.size(); i++)
	{
		outfile << f[i] << "; ";
	}
	outfile << "];\n\nres = [";
	for (int i = 0; i < res.size(); i++)
	{
		outfile << res[i] << "; ";
	}
	outfile << "];\n";
}

void Solver::addSparseToFile(std::ofstream & outfile, std::string name, Eigen::MatrixXd &M)
{
	std::vector<int> ilist;
	std::vector<int> jlist;
	std::vector<double> vlist;

	for (int r = 0; r < M.rows(); ++r)
	{
		for (int c = 0; c < M.cols(); ++c)
		{
			if (M(r, c) != 0)
			{
				ilist.push_back(r);
				jlist.push_back(c);
				vlist.push_back(M(r,c));
			}
		}
	}

	outfile << "i = [\n";
	for (int i = 0; i < ilist.size(); i++)
	{
		outfile << ilist[i] + 1 << " ";
	}
	outfile << "]';\n\nj= [\n";
	for (int i = 0; i < jlist.size(); i++)
	{
		outfile << jlist[i] + 1 << " ";
	}
	outfile << "]';\n\nv = [\n";
	for (int i = 0; i < vlist.size(); i++)
	{
		outfile << vlist[i] << " ";
	}
	outfile << "]';\n" << std::endl;

	int m = (int)M.rows();
	int n = (int)M.cols();

	outfile << name << " = sparse(i,j,v," + std::to_string(m) + "," + std::to_string(n) + ");\n" << std::endl;
}

void Solver::addSparseToFile(std::ofstream & outfile, std::string name, Eigen::SparseMatrix<double>& M)
{
    // https://stackoverflow.com/questions/28685877/convert-an-eigen-matrix-to-triplet-form-c
    std::string entry = " = sparse(";
    outfile << "ijv = [" << std::endl;
    int e = 1;
    for(int k = 0; k < M.outerSize(); ++k) {
        for(Eigen::SparseMatrix<double>::InnerIterator it(M,k); it; ++it) {
            double v = it.value();
            if(std::abs(v) > 1e-10) {
                outfile << it.row()+1 << " "; // row index
                outfile << it.col()+1 << " "; // col index (here it is equal to k)
                outfile << v << std::endl;
                entry.append("ijv(:," + std::to_string(e) + "),");
                ++e;
            }
        }
    }
    outfile << "];" << std::endl;
    outfile << name <<  entry + std::to_string(M.rows()) + "," + std::to_string(M.cols()) + ");\n" << std::endl;
}

void Solver::addVectorToFile(std::ofstream & outfile, std::string name, Eigen::VectorXd & M)
{
	outfile << name << " = [\n";

	for (int i = 0; i < M.size(); i++)
	{
		outfile << M[i] << " ";
	}
	outfile << "]';\n" << std::endl;
}

