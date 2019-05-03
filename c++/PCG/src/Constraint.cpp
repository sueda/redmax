#include "Constraint.h"

#include "RigidBodyMain.h"
#include "ChronoTimer.h"

#include <iomanip>
#include <fstream>

#ifdef ONLINE_MODE
#include "online/GLSL.h"
#endif

Elastic::Elastic() :
	Constraint(constraintType::elastic)
{
	constraintNum = 0;
	constraintNumMax = 0;
}

Elastic::Elastic(std::shared_ptr<Block> ba, std::shared_ptr<Block> bb,
	Eigen::Vector3d pa, Eigen::Vector3d pb, double rest, double stiffness, double damping) :
	Constraint(constraintType::elastic)
{
	blocka = ba;
	blockb = bb;
	posa = pa;
	posb = pb;
	L = rest;
	k = stiffness;
	d = damping;
	constraintNum = 0;
	constraintNumMax = 0;
}

void Elastic::initJoint(const std::unique_ptr<StateSolve>& SS, const std::unique_ptr<LinkageSystem>& LS, const std::unique_ptr<State>& S)
{

}

void Elastic::computeCGProd(Eigen::VectorXd & M_x, const Eigen::VectorXd & x, const std::unique_ptr<StateSolve>& SS, const std::unique_ptr<LinkageSystem>& LS, const std::unique_ptr<State>& S)
{
	int ia = blocka->joint->jindex;
	int ib = blockb->joint->jindex;
	//if (blocka != nullptr)
	//{
	//	M_x.segment<6>(ia * 6) -= SS->h * SS->h * K.block<6, 6>(0, 0) * x.segment<6>(ia * 6);
	//	M_x.segment<6>(ia * 6) -= SS->h * D.block<6, 6>(0, 0) * x.segment<6>(ia * 6);
	//}
	//if (blockb != nullptr)
	//{
	//	M_x.segment<6>(ib * 6) -= SS->h * SS->h * K.block<6, 6>(6, 6) * x.segment<6>(ib * 6);
	//	M_x.segment<6>(ib * 6) -= SS->h * D.block<6, 6>(6, 6) * x.segment<6>(ib * 6);
	//}
	//if (blocka != nullptr && blockb != nullptr)
	//{
	//	M_x.segment<6>(ia * 6) -= SS->h * SS->h * K.block<6, 6>(0, 6) * x.segment<6>(ib * 6);
	//	M_x.segment<6>(ib * 6) -= SS->h * SS->h * K.block<6, 6>(6, 0) * x.segment<6>(ia * 6);
	//	
	//	M_x.segment<6>(ia * 6) -= SS->h * D.block<6, 6>(0, 6) * x.segment<6>(ib * 6);
	//	M_x.segment<6>(ib * 6) -= SS->h * D.block<6, 6>(6, 0) * x.segment<6>(ia * 6);
	//}
	double h2 = SS->h*SS->h;
	if (blocka != nullptr && blockb != nullptr)
	{
		for (int r = 0; r < 6; ++r)
		{
			for (int c = 0; c < 6; ++c)
			{
				// blocka
				M_x(ia * 6 + r) -= h2 * K(r, c) * x[ia * 6 + c];
				M_x(ia * 6 + r) -= SS->h * D(r, c) * x[ia * 6 + c];
				// blockb
				M_x(ib * 6 + r) -= h2 * K(r + 6, c + 6) * x[ib * 6 + c];
				M_x(ib * 6 + r) -= SS->h * D(r + 6, c + 6) * x[ib * 6 + c];
				// both
				M_x(ia * 6 + r) -= h2 * K(r, c + 6) * x[ib * 6 + c];
				M_x(ib * 6 + r) -= h2 * K(r + 6, c) * x[ia * 6 + c];
				M_x(ia * 6 + r) -= SS->h * D(r, c + 6) * x[ib * 6 + c];
				M_x(ib * 6 + r) -= SS->h * D(r + 6, c) * x[ia * 6 + c];
			}
		}
	}
	else if (blocka != nullptr)
	{
		for (int r = 0; r < 6; ++r)
		{
			for (int c = 0; c < 6; ++c)
			{
				M_x(ia * 6 + r) -= h2 * K(r, c) * x[ia * 6 + c];
				M_x(ia * 6 + r) -= SS->h * D(r, c) * x[ia * 6 + c];
			}
		}
	}
	else if (blockb != nullptr)
	{
		for (int r = 0; r < 6; ++r)
		{
			for (int c = 0; c < 6; ++c)
			{
				M_x(ib * 6 + r) -= h2 * K(r + 6, c + 6) * x[ib * 6 + c];
				M_x(ib * 6 + r) -= SS->h * D(r + 6, c + 6) * x[ib * 6 + c];
			}
		}
	}
}

void Elastic::update(std::unique_ptr<StateSolve> &SS, const std::unique_ptr<LinkageSystem> &LS, const std::unique_ptr<State> &S)
{
	K = Eigen::Matrix<double, 12, 12>::Zero();
	D = Eigen::Matrix<double, 12, 12>::Zero();

	Eigen::Matrix4d Ea = Eigen::Matrix4d::Identity();
	Eigen::Matrix4d Eb = Eigen::Matrix4d::Identity();
	Vector6d phia = Vector6d::Zero();
	Vector6d phib = Vector6d::Zero();

	if (blocka != nullptr)
	{
		Ea = S->E[blocka->joint->jindex];
		phia = S->v.segment<6>(blocka->joint->jindex * 6);
	}
	if (blockb != nullptr)
	{
		Eb = S->E[blockb->joint->jindex];
		phib = S->v.segment<6>(blockb->joint->jindex * 6);
	}
	
	Eigen::Vector4d xa_w = Ea * Eigen::Vector4d(posa[0], posa[1], posa[2], 1);
	Eigen::Vector4d xb_w = Eb * Eigen::Vector4d(posb[0], posb[1], posb[2], 1);
	Eigen::Vector3d dx_w = (xb_w - xa_w).segment<3>(0);
	double l = dx_w.norm();
	double i_l = (1 / l);

	Eigen::Matrix3d Ra = Ea.block<3, 3>(0, 0);
	Eigen::Matrix3d Rb = Eb.block<3, 3>(0, 0);

	Matrix3x6d gamma_a = Rigid::gamma(posa);
	Matrix3x6d gamma_b = Rigid::gamma(posb);

	Eigen::Vector3d va_w = Ra * gamma_a * phia;
	Eigen::Vector3d vb_w = Rb * gamma_b * phib;
	double v = (dx_w.transpose() * i_l) * (vb_w - va_w);

	Vector6d fx_a = -i_l * gamma_a.transpose() * Ra.transpose() * dx_w;
	Vector6d fx_b = i_l * gamma_b.transpose() * Rb.transpose() * dx_w;
	Eigen::Matrix<double, 12, 1> fn;
	fn.block<6,1>(0,0) = fx_a;
	fn.block<6,1>(6,0) = fx_b;

	//std::cout << blocka->joint->childName << " " << blockb->joint->childName << std::endl;
	//std::cout << L << std::endl;

	double fs = k * (l - L) / L - d * v / L;

	Eigen::Matrix<double, 1, 3> dfsdx = -k / L * dx_w.transpose() * i_l;
	Eigen::Matrix<double, 1, 12> dfsdE;
	dfsdE.block<1, 6>(0, 0) = dfsdx * Ra * gamma_a;
	dfsdE.block<1, 6>(0, 6) = -dfsdx * Rb * gamma_b;

	// Stiffness
	Eigen::MatrixXd K_temp = fn * dfsdE;
	K = -0.5 * (K_temp + K_temp.transpose());

	// Damping
	Eigen::Matrix<double, 1, 3> dfmdv = d / L * dx_w.transpose() * i_l;
	Eigen::Matrix<double, 1, 12> dfmdphi;
	dfmdphi.block<1, 6>(0, 0) = dfmdv * Ra * gamma_a;
	dfmdphi.block<1, 6>(0, 6) = -dfmdv * Rb * gamma_b;
	D = fn * dfmdphi;

	/// Simplify stiffness calcs - 1/8/2019
	//Eigen::Matrix3d RaT = Ea.block<3, 3>(0, 0).transpose();
	//Eigen::Matrix3d RbT = Eb.block<3, 3>(0, 0).transpose();
	//Eigen::VectorXd gammaRx = Eigen::VectorXd::Zero(12);
	//gammaRx.segment<6>(0) = gamma_a.transpose()*RaT*dx_w;
	//gammaRx.segment<6>(6) = -gamma_b.transpose()*RbT*dx_w;
	//Eigen::MatrixXd temp = Eigen::MatrixXd::Zero(3, 12);
	//temp.block<3, 6>(0, 0) = Ea.block<3, 3>(0, 0)*gamma_a;
	//temp.block<3, 6>(0, 6) = -Eb.block<3, 3>(0, 0)*gamma_b;
	//Eigen::MatrixXd dfnd = Eigen::MatrixXd::Zero(1, 12);
	//Eigen::MatrixXd d_T = (dx_w / (l*l*l)).transpose();
	//dfnd.block<1, 6>(0, 0) = d_T * temp.block<3, 6>(0, 0);
	//dfnd.block<1, 6>(0, 6) = d_T * temp.block<3, 6>(0, 6);
	//Eigen::MatrixXd Kn1 = -gammaRx * dfnd;
	//Eigen::MatrixXd Kn2 = Eigen::MatrixXd::Zero(12, 12);
	//Kn2.block<3, 3>(0, 3) = Rigid::bracket3(posa);
	//Kn2.block<3, 3>(3, 0) = Rigid::bracket3(RaT * (Ea.block<3,1>(0,3) - xb_w.segment<3>(0)));
	//Kn2.block<3, 3>(0, 0) = Kn2.block<3, 3>(0, 3) * Kn2.block<3, 3>(3, 0);
	//Kn2.block<3, 3>(3, 3) = Eigen::Matrix3d::Identity();
	//Kn2.block<3, 3>(3, 9) = -RaT * Eb.block<3, 3>(0, 0);
	//Kn2.block<3, 3>(0, 9) = Kn2.block<3, 3>(0, 3) * Kn2.block<3, 3>(3, 9);
	//Kn2.block<3, 3>(6, 9) = Rigid::bracket3(posb);
	//Kn2.block<3, 3>(9, 6) = Rigid::bracket3(RbT * (Eb.block<3, 1>(0, 3) - xa_w.segment<3>(0)));
	//Kn2.block<3, 3>(6, 6) = Kn2.block<3, 3>(6, 9) * Kn2.block<3, 3>(9, 6);
	//Kn2.block<3, 3>(9, 9) = Eigen::Matrix3d::Identity();
	//Kn2.block<3, 3>(9, 3) = -RbT * Ea.block<3, 3>(0, 0);
	//Kn2.block<3, 3>(6, 3) = Kn2.block<3, 3>(6, 9) * Kn2.block<3, 3>(9, 3);
	//Kn2.block<3, 3>(3, 6) = -Kn2.block<3, 3>(3, 9) * Kn2.block<3, 3>(6, 9);
	//Kn2.block<3, 3>(0, 6) = Kn2.block<3, 3>(0, 3) * Kn2.block<3, 3>(3, 6);
	//Kn2.block<3, 3>(9, 0) = -Kn2.block<3, 3>(9, 3) * Kn2.block<3, 3>(0, 3);
	//Kn2.block<3, 3>(6, 0) = Kn2.block<3, 3>(6, 9) * Kn2.block<3, 3>(9, 0);
	//Kn2 = i_l * Kn2;
	//d_T = (-k / L / l * dx_w).transpose();
	//Eigen::MatrixXd dfsdE = Eigen::MatrixXd::Zero(1, 12);
	//dfsdE.block<1, 6>(0, 0) = d_T * temp.block<3, 6>(0, 0);
	//dfsdE.block<1, 6>(0, 6) = d_T * temp.block<3, 6>(0, 6);
	//Eigen::Vector3d vadot = RaT.transpose() * gamma_a * phia;
	//Eigen::Vector3d vbdot = RbT.transpose() * gamma_b * phib;
	//double ldot = dx_w.dot(vbdot - vadot) / l;
	//double fs = k * (l - L) / L - d * ldot / L;
	//Eigen::VectorXd fn = -i_l * gammaRx;
	//d_T = (-d / L / l * dx_w).transpose();
	//Eigen::MatrixXd dfsdO = Eigen::MatrixXd::Zero(1, 12);
	//dfsdO.block<1, 6>(0, 0) = d_T * temp.block<3, 6>(0, 0);
	//dfsdO.block<1, 6>(0, 6) = d_T * temp.block<3, 6>(0, 6);
	//D = -fn * dfsdO;
	//Eigen::MatrixXd Ktemp = fn * dfsdE + fs * (Kn1 + Kn2);
	//K = -0.5 * (Ktemp + Ktemp.transpose());

	int ia = blocka->joint->jindex;
	int ib = blockb->joint->jindex;
	//SS->Km.block<6, 6>(ia * 6, ia * 6) = SS->Km.block<6, 6>(ia * 6, ia * 6) - (K.block<6, 6>(0, 0));
	//SS->Km.block<6, 6>(ib * 6, ib * 6) = SS->Km.block<6, 6>(ib * 6, ib * 6) - (K.block<6, 6>(6, 6));
	//SS->Km.block<6, 6>(ia * 6, ib * 6) = SS->Km.block<6, 6>(ia * 6, ib * 6) - (K.block<6, 6>(0, 6));
	//SS->Km.block<6, 6>(ib * 6, ia * 6) = SS->Km.block<6, 6>(ib * 6, ia * 6) - (K.block<6, 6>(6, 0));

	//SS->Dm.block<6, 6>(ia * 6, ia * 6) = SS->Dm.block<6, 6>(ia * 6, ia * 6) - (D.block<6, 6>(0, 0));
	//SS->Dm.block<6, 6>(ib * 6, ib * 6) = SS->Dm.block<6, 6>(ib * 6, ib * 6) - (D.block<6, 6>(6, 6));
	//SS->Dm.block<6, 6>(ia * 6, ib * 6) = SS->Dm.block<6, 6>(ia * 6, ib * 6) - (D.block<6, 6>(0, 6));
	//SS->Dm.block<6, 6>(ib * 6, ia * 6) = SS->Dm.block<6, 6>(ib * 6, ia * 6) - (D.block<6, 6>(6, 0));
	for (int r = 0; r < 6; ++r)
	{
		for (int c = 0; c < 6; ++c)
		{
			if(std::abs(K(r, c)) > THRESHOLD)
				SS->Kmlist.push_back(Eigen::Triplet<double>(ia * 6 + r, ia * 6 + c, -K(r, c)));
			if (std::abs(K(r+6, c+6)) > THRESHOLD)
				SS->Kmlist.push_back(Eigen::Triplet<double>(ib * 6 + r, ib * 6 + c, -K(r+6, c+6)));
			if (std::abs(K(r, c+6)) > THRESHOLD)
				SS->Kmlist.push_back(Eigen::Triplet<double>(ia * 6 + r, ib * 6 + c, -K(r, c+6)));
			if (std::abs(K(r+6, c)) > THRESHOLD)
				SS->Kmlist.push_back(Eigen::Triplet<double>(ib * 6 + r, ia * 6 + c, -K(r+6, c)));

			if (std::abs(D(r, c)) > THRESHOLD)
				SS->Dmlist.push_back(Eigen::Triplet<double>(ia * 6 + r, ia * 6 + c, -D(r, c)));
			if (std::abs(D(r + 6, c + 6)) > THRESHOLD)
				SS->Dmlist.push_back(Eigen::Triplet<double>(ib * 6 + r, ib * 6 + c, -D(r + 6, c + 6)));
			if (std::abs(D(r, c + 6)) > THRESHOLD)
				SS->Dmlist.push_back(Eigen::Triplet<double>(ia * 6 + r, ib * 6 + c, -D(r, c + 6)));
			if (std::abs(D(r + 6, c)) > THRESHOLD)
				SS->Dmlist.push_back(Eigen::Triplet<double>(ib * 6 + r, ia * 6 + c, -D(r + 6, c)));
		}
	}

	//if (isnan(fs))
	//	std::cout << "fs" << std::endl;
	//if(fx_a.hasNaN())
	//	std::cout << "fx_a" << std::endl;
	//if (fx_b.hasNaN())
	//	std::cout << "fx_b" << std::endl;

	SS->f.segment<6>(ia * 6) -= fs * fx_a;
	SS->f.segment<6>(ib * 6) -= fs * fx_b;

	blocka->joint->Kmd = blocka->joint->Kmd + K.block<6, 6>(0, 0);
	blocka->joint->Dmd = blocka->joint->Dmd + D.block<6, 6>(0, 0);
	blockb->joint->Kmd = blockb->joint->Kmd + K.block<6, 6>(6, 6);
	blockb->joint->Dmd = blockb->joint->Dmd + D.block<6, 6>(6, 6);
}

void Elastic::update(std::unique_ptr<StateSolve> &SS, const std::unique_ptr<LinkageSystem> &LS, const std::unique_ptr<State> &S, std::unique_ptr<StateDeriv> &DS)
{
	
}

void Elastic::updateJoint(std::unique_ptr<StateSolve> &SS, const std::unique_ptr<LinkageSystem> &LS, const std::unique_ptr<State>& S)
{
	K = Eigen::Matrix<double, 12, 12>::Zero();
	D = Eigen::Matrix<double, 12, 12>::Zero();

	Eigen::Matrix4d Ea = Eigen::Matrix4d::Identity();
	Eigen::Matrix4d Eb = Eigen::Matrix4d::Identity();
	Vector6d phia = Vector6d::Zero();
	Vector6d phib = Vector6d::Zero();

	if (blocka != nullptr)
	{
		Ea = S->E[blocka->joint->jindex];
		phia = S->v.segment<6>(blocka->joint->jindex * 6);
	}
	if (blockb != nullptr)
	{
		Eb = S->E[blockb->joint->jindex];
		phib = S->v.segment<6>(blockb->joint->jindex * 6);
	}

	Eigen::Vector4d xa_w = Ea * Eigen::Vector4d(posa[0], posa[1], posa[2], 1);
	Eigen::Vector4d xb_w = Eb * Eigen::Vector4d(posb[0], posb[1], posb[2], 1);
	Eigen::Vector3d dx_w = (xb_w - xa_w).segment<3>(0);
	double l = dx_w.norm();
	double i_l = (1 / l);

	Eigen::Matrix3d Ra = Ea.block<3, 3>(0, 0);
	Eigen::Matrix3d Rb = Eb.block<3, 3>(0, 0);

	Matrix3x6d gamma_a = Rigid::gamma(posa);
	Matrix3x6d gamma_b = Rigid::gamma(posb);

	Eigen::Vector3d va_w = Ra * gamma_a * phia;
	Eigen::Vector3d vb_w = Rb * gamma_b * phib;
	double v = (dx_w.transpose() * i_l) * (vb_w - va_w);

	Vector6d fx_a = -i_l * gamma_a.transpose() * Ra.transpose() * dx_w;
	Vector6d fx_b = i_l * gamma_b.transpose() * Rb.transpose() * dx_w;
	Eigen::Matrix<double, 12, 1> fn;
	fn.block<6, 1>(0, 0) = fx_a;
	fn.block<6, 1>(6, 0) = fx_b;

	double fs = k * (l - L) / L - d * v / L;

	Eigen::Matrix<double, 1, 3> dfsdx = -k / L * dx_w.transpose() * i_l;
	Eigen::Matrix<double, 1, 12> dfsdE;
	dfsdE.block<1, 6>(0, 0) = dfsdx * Ra * gamma_a;
	dfsdE.block<1, 6>(0, 6) = -dfsdx * Rb * gamma_b;

	// Stiffness
	Eigen::MatrixXd K_temp = fn * dfsdE;
	K = -0.5 * (K_temp + K_temp.transpose());

	// Damping
	Eigen::Matrix<double, 1, 3> dfmdv = d / L * dx_w.transpose() * i_l;
	Eigen::Matrix<double, 1, 12> dfmdphi;
	dfmdphi.block<1, 6>(0, 0) = dfmdv * Ra * gamma_a;
	dfmdphi.block<1, 6>(0, 6) = -dfmdv * Rb * gamma_b;
	D = fn * dfmdphi;

	int ia = blocka->joint->jindex;
	int ib = blockb->joint->jindex;
	SS->f.segment<6>(ia * 6) -= fs * fx_a;
	SS->f.segment<6>(ib * 6) -= fs * fx_b;

	blocka->joint->Kmd = blocka->joint->Kmd + K.block<6, 6>(0, 0);
	blocka->joint->Dmd = blocka->joint->Dmd + D.block<6, 6>(0, 0);
	blockb->joint->Kmd = blockb->joint->Kmd + K.block<6, 6>(6, 6);
	blockb->joint->Dmd = blockb->joint->Dmd + D.block<6, 6>(6, 6);

	if (simtype == simType::PCG)
	{
		S->Kmd[blocka->joint->jindex] += K.block<6, 6>(0, 0);
		S->Dmd[blocka->joint->jindex] += D.block<6, 6>(0, 0);
		S->Kmd[blockb->joint->jindex] += K.block<6, 6>(6, 6);
		S->Dmd[blockb->joint->jindex] += D.block<6, 6>(6, 6);
	}
}

void Elastic::draw(const std::unique_ptr<LinkageSystem> &LS, const std::unique_ptr<State> &S) const
{
#ifdef ONLINE_MODE
	// Draw the string connection point on block a
	Eigen::Matrix4d displayE = S->E[blocka->joint->jindex] * blocka->iToD;
	Eigen::Vector3d localposa = displayE.block<3, 3>(0, 0)*posa + displayE.block<3, 1>(0, 3);

	glColor3f(0.2f, 0.2f, 0.2f);
	glPointSize(8);
	glBegin(GL_POINTS);
	glVertex3f((float)localposa[0], (float)localposa[1], (float)localposa[2]);
	glEnd();

	// Draw the corresponding point on block b
	displayE = S->E[blockb->joint->jindex] * blockb->iToD;
	Eigen::Vector3d localposb = displayE.block<3, 3>(0, 0)*posb + displayE.block<3, 1>(0, 3);

	glColor3f(0.2f, 0.2f, 0.2f);
	glPointSize(8);
	glBegin(GL_POINTS);
	glVertex3f((float)localposb[0], (float)localposb[1], (float)localposb[2]);
	glEnd();

	// Draw the string itself
	glColor3f(0.4f, 0.4f, 0.4f);
	glLineWidth(4);
	glBegin(GL_LINES);
	glVertex3f((float)localposa[0], (float)localposa[1], (float)localposa[2]);
	glVertex3f((float)localposb[0], (float)localposb[1], (float)localposb[2]);
	glEnd();
#endif
}

Json::Value Elastic::exportBrender()
{
	return Json::Value();
}

void Elastic::exportJson(std::ofstream & outfile, const std::unique_ptr<State> &S)
{
	// Draw the string connection point on block a
	Eigen::Matrix4d displayE = blocka->joint->E_draw;// *blocka->iToD;
	Eigen::Vector3d localposa = displayE.block<3, 3>(0, 0)*posa + displayE.block<3, 1>(0, 3);

	// Draw the corresponding point on block b
	displayE = blockb->joint->E_draw;// *blockb->iToD;
	Eigen::Vector3d localposb = displayE.block<3, 3>(0, 0)*posb + displayE.block<3, 1>(0, 3);

	std::string name = blocka->joint->childName + "_Elastic";

	outfile << "o " << name << std::endl;
	outfile << "g strand" << std::endl;

	outfile << "v " << localposa[0] << " " << localposa[1] << " " << localposa[2] << std::endl;
	outfile << "v " << localposb[0] << " " << localposb[1] << " " << localposb[2] << std::endl;

	outfile << "vt " << 0.0 << std::endl;
	outfile << "vt " << 1.0 << std::endl;
}

void Elastic::initLengthToRest(const std::unique_ptr<LinkageSystem> &LS, const std::unique_ptr<State> &S)
{
	std::shared_ptr<Joint> ja = blocka->joint;
	std::shared_ptr<Joint> jb = blockb->joint;
	Eigen::Matrix4d Ea = S->E[ja->jindex];
	Eigen::Matrix4d Eb = S->E[jb->jindex];

	//Eigen::Vector4d worlda = blocka->E_wi0 * Eigen::Vector4d(posa[0], posa[1], posa[2], 1);
	//Eigen::Vector4d worldb = blockb->E_wi0 * Eigen::Vector4d(posb[0], posb[1], posb[2], 1);
	Eigen::Vector4d worlda = Ea * Eigen::Vector4d(posa[0], posa[1], posa[2], 1);
	Eigen::Vector4d worldb = Eb * Eigen::Vector4d(posb[0], posb[1], posb[2], 1);

	Eigen::Vector3d a = (worldb - worlda).segment<3>(0);
	L = a.norm();
}

Constraint::Constraint()
{
}

Constraint::Constraint(constraintType type) :
	type(type)
{
}

CloseHinge::CloseHinge() :
	Constraint(constraintType::closehinge)
{
	constraintNum = 2;
	constraintNumMax = 12;
}

CloseHinge::CloseHinge(std::shared_ptr<Block> ba, std::shared_ptr<Block> bb, Eigen::Vector3d pa, Eigen::Vector3d pb, Eigen::Vector3d a, 
	double rest, double stiffness, double damping) :
	Constraint(constraintType::closehinge)
{
	blocka = ba;
	blockb = bb;
	posa = pa;
	posb = pb;
	axis = a;
	q0 = rest;
	k = stiffness;
	d = damping;
	constraintNum = 2;
	constraintNumMax = 12;

	double ax = std::abs(axis[0]);
	double ay = std::abs(axis[1]);
	double az = std::abs(axis[2]);
	if (ax <= ay && ax <= az)
		v1 = Eigen::Vector3d(1, 0, 0);
	else if (ay <= az && ay <= ax)
		v1 = Eigen::Vector3d(0, 1, 0);
	else
		v1 = Eigen::Vector3d(0, 0, 1);
	v2 = axis.cross(v1);
	v2 = v2.normalized();
	v1 = v2.cross(axis);
	v1 = v1.normalized();
}

void CloseHinge::initJoint(const std::unique_ptr<StateSolve>& SS, const std::unique_ptr<LinkageSystem>& LS, const std::unique_ptr<State>& S)
{
}

void CloseHinge::computeCGProd(Eigen::VectorXd & M_x, const Eigen::VectorXd & x, const std::unique_ptr<StateSolve>& SS, const std::unique_ptr<LinkageSystem>& LS, const std::unique_ptr<State>& S)
{
	std::shared_ptr<Joint> ja = blocka->joint;
	std::shared_ptr<Joint> jb = blockb->joint;

	int ia = ja->jindex;
	int ib = jb->jindex;
	//unknown
}

void CloseHinge::update(std::unique_ptr<StateSolve>& SS, const std::unique_ptr<LinkageSystem>& LS, const std::unique_ptr<State>& S)
{
	// redmax update
	std::shared_ptr<Joint> ja = blocka->joint;
	std::shared_ptr<Joint> jb = blockb->joint;

	//SS->Gm.block<1, 6>(constraint_index, ja->jindex * 6) = v1.transpose() * S->E[ja->jindex].block<3, 3>(0, 0) * Rigid::gamma(posa);
	//SS->Gm.block<1, 6>(constraint_index + 1, ja->jindex * 6) = v2.transpose() * S->E[ja->jindex].block<3, 3>(0, 0) * Rigid::gamma(posa);

	//SS->Gm.block<1, 6>(constraint_index, jb->jindex * 6) = -v1.transpose() * S->E[jb->jindex].block<3, 3>(0, 0) * Rigid::gamma(posb);
	//SS->Gm.block<1, 6>(constraint_index + 1, jb->jindex * 6) = -v2.transpose() * S->E[jb->jindex].block<3, 3>(0, 0) * Rigid::gamma(posb);

	Vector6d av1 = v1.transpose() * S->E[ja->jindex].block<3, 3>(0, 0) * Rigid::gamma(posa);
	Vector6d av2 = v2.transpose() * S->E[ja->jindex].block<3, 3>(0, 0) * Rigid::gamma(posa);
	Vector6d bv1 = -v1.transpose() * S->E[jb->jindex].block<3, 3>(0, 0) * Rigid::gamma(posb);
	Vector6d bv2 = -v2.transpose() * S->E[jb->jindex].block<3, 3>(0, 0) * Rigid::gamma(posb);

	for (int r = 0; r < 6; ++r)
	{
		if (std::abs(av1[r]) > THRESHOLD)
		{
			SS->Gmlist.push_back(Eigen::Triplet<double>(constraint_index, ja->jindex * 6 + r, av1[r]));
			//SS->GmTlist.push_back(Eigen::Triplet<double>(ja->jindex * 6 + r, constraint_index, av1[r]));
		}
		if (std::abs(av2[r]) > THRESHOLD)
		{
			SS->Gmlist.push_back(Eigen::Triplet<double>(constraint_index + 1, ja->jindex * 6 + r, av2[r]));
			//SS->GmTlist.push_back(Eigen::Triplet<double>(ja->jindex * 6 + r, constraint_index + 1, av2[r]));
		}
		if (std::abs(bv1[r]) > THRESHOLD)
		{
			SS->Gmlist.push_back(Eigen::Triplet<double>(constraint_index, jb->jindex * 6 + r, bv1[r]));
			//SS->GmTlist.push_back(Eigen::Triplet<double>(jb->jindex * 6 + r, constraint_index, bv1[r]));
		}
		if (std::abs(bv2[r]) > THRESHOLD)
		{
			SS->Gmlist.push_back(Eigen::Triplet<double>(constraint_index + 1, jb->jindex * 6 + r, bv2[r]));
			//SS->GmTlist.push_back(Eigen::Triplet<double>(jb->jindex * 6 + r, constraint_index + 1, bv2[r]));
		}
	}


	Eigen::Vector4d dx = S->E[ja->jindex] * Eigen::Vector4d(posa[0], posa[1], posa[2], 1) - S->E[jb->jindex] * Eigen::Vector4d(posb[0], posb[1], posb[2], 1);
	SS->gm[constraint_index] = v1.transpose() * dx.segment<3>(0);
	SS->gm[constraint_index + 1] = v2.transpose() * dx.segment<3>(0);
}

void CloseHinge::update(std::unique_ptr<StateSolve>& SS, const std::unique_ptr<LinkageSystem>& LS, const std::unique_ptr<State>& S, std::unique_ptr<StateDeriv>& DS)
{
}

void CloseHinge::updateJoint(std::unique_ptr<StateSolve>& SS, const std::unique_ptr<LinkageSystem>& LS, const std::unique_ptr<State>& S)
{
	int ia = blocka->joint->jindex;
	int ib = blockb->joint->jindex;

	//Vector6d av1 = v1.transpose() * S->E[ja->jindex].block<3, 3>(0, 0) * Rigid::gamma(posa);
	//Vector6d av2 = v2.transpose() * S->E[ja->jindex].block<3, 3>(0, 0) * Rigid::gamma(posa);
	//Vector6d bv1 = -v1.transpose() * S->E[jb->jindex].block<3, 3>(0, 0) * Rigid::gamma(posb);
	//Vector6d bv2 = -v2.transpose() * S->E[jb->jindex].block<3, 3>(0, 0) * Rigid::gamma(posb);

	//std::cout << av1 << std::endl;
	//std::cout << av2 << std::endl;
	//std::cout << bv1 << std::endl;
	//std::cout << bv2 << std::endl;

	SS->Gm.block<1, 6>(0, ia * 6) = v1.transpose() * S->E[ia].block<3, 3>(0, 0) * Rigid::gamma(posa);
	SS->Gm.block<1, 6>(1, ia * 6) = v2.transpose() * S->E[ia].block<3, 3>(0, 0) * Rigid::gamma(posa);

	SS->Gm.block<1, 6>(0, ib * 6) = -v1.transpose() * S->E[ib].block<3, 3>(0, 0) * Rigid::gamma(posb);
	SS->Gm.block<1, 6>(1, ib * 6) = -v2.transpose() * S->E[ib].block<3, 3>(0, 0) * Rigid::gamma(posb);

	//SS->GmT.block<6, 1>(ia * 6, 0) = SS->Gm.block<1, 6>(0, 0).transpose();
	//SS->GmT.block<6, 1>(ia * 6, 1) = SS->Gm.block<1, 6>(1, 0).transpose();

	//SS->GmT.block<6, 1>(ib * 6, 0) = SS->Gm.block<1, 6>(0, 6).transpose();
	//SS->GmT.block<6, 1>(ib * 6, 1) = SS->Gm.block<1, 6>(1, 6).transpose();

	Eigen::Vector4d temp = S->E[ia] * Eigen::Vector4d(posa[0], posa[1], posa[2], 1) - S->E[ib] * Eigen::Vector4d(posb[0], posb[1], posb[2], 1);
	Eigen::Vector3d dx = Eigen::Vector3d(temp[0], temp[1], temp[2]);
	SS->gm[constraint_index] += v1.dot(dx);
	SS->gm[constraint_index + 1] += v2.dot(dx);
}

void CloseHinge::draw(const std::unique_ptr<LinkageSystem>& LS, const std::unique_ptr<State>& S) const
{
#ifdef ONLINE_MODE
	Eigen::Vector3d apos, bpos;
	Eigen::Matrix4d E = Eigen::Matrix4d::Identity();

	E = S->E[blocka->joint->jindex];
	apos = E.block<3, 3>(0, 0)*posa + E.block<3, 1>(0, 3);

	glLineWidth(8);
	glBegin(GL_LINES);
	glColor3f((float)0.5, (float)1, float(0.5));
	glVertex3f((float)apos[0], (float)apos[1], (float)apos[2]);
	apos = apos + Eigen::Vector3d(0, 0, 0.5);
	glVertex3f((float)apos[0], (float)apos[1], (float)apos[2]);
	glEnd();
	glLineWidth(1);

	E = S->E[blockb->joint->jindex];
	bpos = E.block<3, 3>(0, 0)*posb + E.block<3, 1>(0, 3);

	glLineWidth(4);
	glBegin(GL_LINES);
	glVertex3f((float)bpos[0], (float)bpos[1], (float)bpos[2]);
	bpos = bpos - Eigen::Vector3d(0, 0, 0.5);
	glVertex3f((float)bpos[0], (float)bpos[1], (float)bpos[2]);
	glEnd();
	glLineWidth(1);
#endif
}

Json::Value CloseHinge::exportBrender()
{
	return Json::Value();
}

void CloseHinge::exportJson(std::ofstream & outfile, const std::unique_ptr<State>& S)
{
}

void CloseHinge::initLengthToRest(const std::unique_ptr<LinkageSystem>& LS, const std::unique_ptr<State>& S)
{
}

SpringPoint::SpringPoint() :
	Constraint(constraintType::springpoint)
{
	constraintNum = 0;
	constraintNumMax = 0;
}

SpringPoint::SpringPoint(std::shared_ptr<Block> ba, Eigen::Vector3d pa, Eigen::Vector3d dr, double stiffness) :
	Constraint(constraintType::springpoint)
{
	blocka = ba;
	posa = pa;
	dir = dr;
	k = stiffness;
	constraintNum = 0;
	constraintNumMax = 0;
}

void SpringPoint::initJoint(const std::unique_ptr<StateSolve>& SS, const std::unique_ptr<LinkageSystem>& LS, const std::unique_ptr<State>& S)
{
}

void SpringPoint::computeCGProd(Eigen::VectorXd & M_x, const Eigen::VectorXd & x, const std::unique_ptr<StateSolve>& SS, const std::unique_ptr<LinkageSystem>& LS, const std::unique_ptr<State>& S)
{
	// Stiffness turned off
	//SS->Km.block<6, 6>(blocka->joint->jindex * 6, blocka->joint->jindex * 6) = SS->Km.block<6, 6>(blocka->joint->jindex * 6, blocka->joint->jindex * 6) + K;
}

void SpringPoint::update(std::unique_ptr<StateSolve>& SS, const std::unique_ptr<LinkageSystem>& LS, const std::unique_ptr<State>& S)
{
	K = Eigen::Matrix<double, 6, 6>::Zero();
	Eigen::Vector3d kRTa = k*S->E[blocka->joint->jindex].block<3, 3>(0, 0).transpose()*dir;
	K.block<3, 3>(3, 0) = Rigid::bracket3(kRTa);
	K.block<3, 3>(0, 0) = Rigid::bracket3(posa) * K.block<3, 3>(3, 0);
	K = 0.5 * (K + K.transpose());
	K = 0 * K;

	Matrix3x6d gamma = Rigid::gamma(posa);

	//if (gamma.hasNaN())
	//	std::cout << "gamma" << std::endl;
	//if (kRTa.hasNaN())
	//	std::cout << "kRTa" << std::endl;

	SS->f.segment<6>(blocka->joint->jindex * 6) += gamma.transpose() * kRTa;
}

void SpringPoint::update(std::unique_ptr<StateSolve>& SS, const std::unique_ptr<LinkageSystem>& LS, const std::unique_ptr<State>& S, std::unique_ptr<StateDeriv>& DS)
{
}

void SpringPoint::updateJoint(std::unique_ptr<StateSolve>& SS, const std::unique_ptr<LinkageSystem>& LS, const std::unique_ptr<State>& S)
{
	K = Eigen::Matrix<double, 6, 6>::Zero();
	Eigen::Vector3d kRTa = k * S->E[blocka->joint->jindex].block<3, 3>(0, 0).transpose()*dir;
	K.block<3, 3>(3, 0) = Rigid::bracket3(kRTa);
	K.block<3, 3>(0, 0) = Rigid::bracket3(posa) * K.block<3, 3>(3, 0);
	K = 0.5 * (K + K.transpose());
	K = 0 * K;

	Matrix3x6d gamma = Rigid::gamma(posa);
	SS->f.segment<6>(blocka->joint->jindex * 6) += gamma.transpose() * kRTa; 
}

void SpringPoint::draw(const std::unique_ptr<LinkageSystem>& LS, const std::unique_ptr<State>& S) const
{
#ifdef ONLINE_MODE
	Eigen::Matrix4d E = S->E[blocka->joint->jindex];
	Eigen::Vector3d apos = E.block<3, 3>(0, 0)*posa + E.block<3, 1>(0, 3);

	glLineWidth(8);
	glBegin(GL_LINES);
	glColor3f((float)0.5, (float)1, float(0.5));
	glVertex3f((float)apos[0], (float)apos[1], (float)apos[2]);
	apos = apos + dir * 1;
	glVertex3f((float)apos[0], (float)apos[1], (float)apos[2]);
	glEnd();
	glLineWidth(1);
#endif
}

Json::Value SpringPoint::exportBrender()
{
	Eigen::Quaterniond quat(blocka->joint->E_draw.block<3, 3>(0, 0));

	Json::Value scale(Json::arrayValue);
	scale[0] = 0.3;
	scale[1] = 0.23;
	scale[2] = 0.15;

	Json::Value location(Json::arrayValue);
	location[0] = blocka->joint->E_draw(0, 3);
	location[1] = blocka->joint->E_draw(1, 3);
	location[2] = blocka->joint->E_draw(2, 3);

	Json::Value q(Json::arrayValue);
	q[0] = quat.x();
	q[1] = quat.y();
	q[2] = quat.z();
	q[3] = quat.w();

	Json::Value v;
	v["scale"] = scale;
	v["location"] = location;
	v["quat"] = q;

	return v;
}

void SpringPoint::exportJson(std::ofstream & outfile, const std::unique_ptr<State>& S)
{
}

CloseUniversal::CloseUniversal() :
	Constraint(constraintType::closehinge)
{
	constraintNum = 2;
	constraintNumMax = 12;
}

CloseUniversal::CloseUniversal(std::shared_ptr<Block> ba, std::shared_ptr<Block> bb, Eigen::Vector3d pa, Eigen::Vector3d pb) :
	Constraint(constraintType::closehinge)
{
	blocka = ba;
	blockb = bb;
	posa = pa;
	posb = pb;
	constraintNum = 2;
	constraintNumMax = 12;

	gamma_a = Rigid::gamma(posa);
	gamma_b = Rigid::gamma(posb);
}

void CloseUniversal::initJoint(const std::unique_ptr<StateSolve>& SS, const std::unique_ptr<LinkageSystem>& LS, const std::unique_ptr<State>& S)
{
}

void CloseUniversal::computeCGProd(Eigen::VectorXd & M_x, const Eigen::VectorXd & x, const std::unique_ptr<StateSolve>& SS, const std::unique_ptr<LinkageSystem>& LS, const std::unique_ptr<State>& S)
{
}

void CloseUniversal::update(std::unique_ptr<StateSolve>& SS, const std::unique_ptr<LinkageSystem>& LS, const std::unique_ptr<State>& S)
{
	Eigen::Vector3d v0;
	Eigen::Vector3d v1;
	Eigen::Vector3d v2;

	v0 = blocka->joint->E_wj.block<3,1>(0,0);
	double ax = std::abs(v0[0]);
	double ay = std::abs(v0[1]);
	double az = std::abs(v0[2]);
	if (ax <= ay && ax <= az)
		v1 = Eigen::Vector3d(1, 0, 0);
	else if (ay <= az && ay <= ax)
		v1 = Eigen::Vector3d(0, 1, 0);
	else
		v1 = Eigen::Vector3d(0, 0, 1);
	v2 = v0.cross(v1);
	v2 = v2.normalized();
	v1 = v2.cross(v0);
	v1 = v1.normalized();

	Eigen::Matrix<double, 1, 3> v1T = v1.transpose();
	Eigen::Matrix<double, 1, 3> v2T = v2.transpose();

	// redmax update
	std::shared_ptr<Joint> ja = blocka->joint;
	std::shared_ptr<Joint> jb = blockb->joint;

	Vector6d av1 = v1T * S->E[ja->jindex].block<3, 3>(0, 0) * gamma_a;
	Vector6d av2 = v2T * S->E[ja->jindex].block<3, 3>(0, 0) * gamma_a;
	Vector6d bv1 = -v1T * S->E[jb->jindex].block<3, 3>(0, 0) * gamma_b;
	Vector6d bv2 = -v2T * S->E[jb->jindex].block<3, 3>(0, 0) * gamma_b;

	for (int r = 0; r < 6; ++r)
	{
		if (std::abs(av1[r]) > THRESHOLD)
		{
			SS->Gmlist.push_back(Eigen::Triplet<double>(constraint_index, ja->jindex * 6 + r, av1[r]));
			//SS->GmTlist.push_back(Eigen::Triplet<double>(ja->jindex * 6 + r, constraint_index, av1[r]));
		}
		if (std::abs(av2[r]) > THRESHOLD)
		{
			SS->Gmlist.push_back(Eigen::Triplet<double>(constraint_index + 1, ja->jindex * 6 + r, av2[r]));
			//SS->GmTlist.push_back(Eigen::Triplet<double>(ja->jindex * 6 + r, constraint_index + 1, av2[r]));
		}
		if (std::abs(bv1[r]) > THRESHOLD)
		{
			SS->Gmlist.push_back(Eigen::Triplet<double>(constraint_index, jb->jindex * 6 + r, bv1[r]));
			//SS->GmTlist.push_back(Eigen::Triplet<double>(jb->jindex * 6 + r, constraint_index, bv1[r]));
		}
		if (std::abs(bv2[r]) > THRESHOLD)
		{
			SS->Gmlist.push_back(Eigen::Triplet<double>(constraint_index + 1, jb->jindex * 6 + r, bv2[r]));
			//SS->GmTlist.push_back(Eigen::Triplet<double>(jb->jindex * 6 + r, constraint_index + 1, bv2[r]));
		}
	}

	Eigen::Vector4d dx = S->E[ja->jindex] * Eigen::Vector4d(posa[0], posa[1], posa[2], 1) - S->E[jb->jindex] * Eigen::Vector4d(posb[0], posb[1], posb[2], 1);
	SS->gm[constraint_index] = v1T * dx.segment<3>(0);
	SS->gm[constraint_index + 1] = v2T * dx.segment<3>(0);
}

void CloseUniversal::update(std::unique_ptr<StateSolve>& SS, const std::unique_ptr<LinkageSystem>& LS, const std::unique_ptr<State>& S, std::unique_ptr<StateDeriv>& DS)
{
}

void CloseUniversal::updateJoint(std::unique_ptr<StateSolve>& SS, const std::unique_ptr<LinkageSystem>& LS, const std::unique_ptr<State>& S)
{
	int ia = blocka->joint->jindex;
	int ib = blockb->joint->jindex;

	Eigen::Vector3d v0;
	Eigen::Vector3d v1;
	Eigen::Vector3d v2;

	v0 = blocka->joint->E_wj.block<3, 1>(0, 0);
	double ax = std::abs(v0[0]);
	double ay = std::abs(v0[1]);
	double az = std::abs(v0[2]);
	if (ax <= ay && ax <= az)
		v1 = Eigen::Vector3d(1, 0, 0);
	else if (ay <= az && ay <= ax)
		v1 = Eigen::Vector3d(0, 1, 0);
	else
		v1 = Eigen::Vector3d(0, 0, 1);
	v2 = v0.cross(v1);
	v2 = v2.normalized();
	v1 = v2.cross(v0);
	v1 = v1.normalized();

	Eigen::Matrix<double, 1, 3> v1T = v1.transpose();
	Eigen::Matrix<double, 1, 3> v2T = v2.transpose();

	Vector6d av1 = v1T * S->E[ia].block<3, 3>(0, 0) * gamma_a;
	Vector6d av2 = v2T * S->E[ia].block<3, 3>(0, 0) * gamma_a;
	Vector6d bv1 = -v1T * S->E[ib].block<3, 3>(0, 0) * gamma_b;
	Vector6d bv2 = -v2T * S->E[ib].block<3, 3>(0, 0) * gamma_b;

	SS->Gm.block<1, 6>(constraint_index, ia * 6) = av1;
	SS->Gm.block<1, 6>(constraint_index + 1, ia * 6) = av2;

	SS->Gm.block<1, 6>(constraint_index, ib * 6) = bv1;
	SS->Gm.block<1, 6>(constraint_index + 1, ib * 6) = bv2;

	//SS->GmT.block<6, 1>(ia * 6, constraint_index) = av1.transpose();
	//SS->GmT.block<6, 1>(ia * 6, constraint_index + 1) = av2.transpose();

	//SS->GmT.block<6, 1>(ib * 6, constraint_index) = bv1.transpose();
	//SS->GmT.block<6, 1>(ib * 6, constraint_index + 1) = bv2.transpose();

	Eigen::Vector4d dx = S->E[ia] * Eigen::Vector4d(posa[0], posa[1], posa[2], 1) - S->E[ib] * Eigen::Vector4d(posb[0], posb[1], posb[2], 1);
	SS->gm[constraint_index] = v1T * dx.segment<3>(0);
	SS->gm[constraint_index + 1] = v2T * dx.segment<3>(0);
}

void CloseUniversal::draw(const std::unique_ptr<LinkageSystem>& LS, const std::unique_ptr<State>& S) const
{
#ifdef ONLINE_MODE
	Eigen::Vector3d apos, bpos;
	Eigen::Matrix4d E = Eigen::Matrix4d::Identity();

	E = S->E[blocka->joint->jindex];
	apos = E.block<3, 3>(0, 0)*posa + E.block<3, 1>(0, 3);

	glLineWidth(8);
	glBegin(GL_LINES);
	glColor3f((float)0.5, (float)0.5, float(1.0));
	glVertex3f((float)apos[0], (float)apos[1], (float)apos[2]);
	apos = apos + Eigen::Vector3d(0, 0, 1.5);
	glVertex3f((float)apos[0], (float)apos[1], (float)apos[2]);
	glEnd();
	glLineWidth(1);

	E = S->E[blockb->joint->jindex];
	bpos = E.block<3, 3>(0, 0)*posb + E.block<3, 1>(0, 3);

	glLineWidth(4);
	glBegin(GL_LINES);
	glVertex3f((float)bpos[0], (float)bpos[1], (float)bpos[2]);
	bpos = bpos - Eigen::Vector3d(0, 0, 1.5);
	glVertex3f((float)bpos[0], (float)bpos[1], (float)bpos[2]);
	glEnd();
	glLineWidth(1);
#endif
}

Json::Value CloseUniversal::exportBrender()
{
	return Json::Value();
}

void CloseUniversal::exportJson(std::ofstream & outfile, const std::unique_ptr<State>& S)
{
}

void CloseUniversal::initLengthToRest(const std::unique_ptr<LinkageSystem>& LS, const std::unique_ptr<State>& S)
{
}
