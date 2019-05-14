#define _USE_MATH_DEFINES
#include <ctgmath>
#include <cmath>

#include "Scene.h"

#include "ChronoTimer.h"
#include "Rigid.h"
#include "RigidBodyMain.h"

#include "online/Camera.h"
#include "online/MatrixStack.h"
#include "online/Program.h"

#ifndef _GLIBCXX_USE_NANOSLEEP
#define _GLIBCXX_USE_NANOSLEEP
#endif

#ifndef GLM_FORCE_RADIANS
#define GLM_FORCE_RADIANS
#endif
#include <glm/gtc/type_ptr.hpp>

#ifndef EIGEN_USE_MKL
#define EIGEN_USE_MKL
#endif
#ifndef EIGEN_DONT_ALIGN_STATICALLY
#define EIGEN_DONT_ALIGN_STATICALLY
#endif
#ifndef EIGEN_DONT_ALIGN
#define EIGEN_DONT_ALIGN
#endif
#ifndef EIGEN_NO_STATIC_ASSERT
#define EIGEN_NO_STATIC_ASSERT
#endif
#ifndef EIGEN_DISABLE_UNALIGNED_ARRAY_ASSERT
#define EIGEN_DISABLE_UNALIGNED_ARRAY_ASSERT
#endif
#include <Eigen/SparseCholesky>	
#include <Eigen/OrderingMethods>
#include <Eigen/Core>
#include <Eigen/Geometry>

#include <iostream>
#include <fstream>
#include <string>

static double randDouble(double l, double h)
{
	double r = rand() / (double)RAND_MAX;
	return (1.0 - r) * l + r * h;
}

Scene::Scene(const std::string &RESOURCE_DIR)
{
	t = 0.0;
	rigid_body = std::make_shared<RigidBodyMain>(RESOURCE_DIR, t);
}

Scene::~Scene()
{
}

void Scene::init()
{
	rigid_body->displayUserActions();
}

void Scene::load(const std::string &FILENAME, std::shared_ptr<Camera> &c)
{
	// load the linkage system
	///rigid_body->load(FILENAME);
	rigid_body->setSolveType(simType::PCG);
	//rigid_body->loadTree(5);
	//rigid_body->loadSimpleTree(6);
	//c->setTranslations(0.0f, 14.0f, -34.8f);
   // rigid_body->loadUmbrella(16);
	//rigid_body->loadTest(5);
    rigid_body->loadBridge(20, 20);
	//rigid_body->loadSimpleBridge(30);

	//c->setTranslations(0.0f, -5.0f, -24.8f);
	c->setTranslations(0.0f, -0.0f, -52.8f);
}

int Scene::keyPressed(unsigned int key)
{
	int rbindicator = rigid_body->handleKeyPress(key);

	switch (key) {
	case 'x':
		//rigid_body->load_kf_com();
		std::cout << "keyframing active" << std::endl;
		return false;
	case 'v':
		//rigid_body->validatepartials();
		std::cout << "validating" << std::endl;
		return false;
	}
	
	return rbindicator;
}

int Scene::clickPress(float x, float y, float z)
{
	return 0;
}

int Scene::clickRelease(float x, float y, float z)
{
	return 0;
}

Eigen::VectorXd Scene::step()
{
	Eigen::VectorXd result = rigid_body->step();
	//brender->exportBrender(rigid_body->getTime());
	return result;
}

void Scene::batchTest()
{
	// save to file named...
	std::string filename = "../simdata/pcg_tree_short.m";//pcg_umbrella_optimize.m";

	// choose simulations to run
	std::vector<simType> testSimTypes{
		PCG
	};

	// name every test
	std::vector<std::string> testNames{
		"pcg"
	};

	std::vector<int> testSizes{
		5
		//8, 16, 32, 64, 128
		//10, 20, 50, 100, 200, 500, 1000
		//5, 6, 7, 8, 9, 10//, 11, 12
		//10, 20, 40, 80, 160, 320, 640, 1280
		//6, 12, 24, 48, 96, 192, 384
		//30, 60, 120, 240, 480, 960
	};

	// select test quantities
	bool ground_truth = false;  // prerequisite for using avg step error
	bool avg_step_lstsq = false; // using the first result as ground truth
	bool qdot_cert = false;
	bool leaf_pos_cert = true;
	bool iterations = true;
	// timing enabled

	// rendering?
	bool render = false;
	int torend = 0;
	if (render)
	{
		brender = BrenderManager::getInstance();
		brender->add(rigid_body);
	}

	// run the scene simulation for 10 simulation seconds and collect data
	std::vector<std::vector<Eigen::VectorXd>> raw_result;  // ground truth, EXPENSIVE
	raw_result.resize(testSimTypes.size());

	// save last qdot
	std::vector<std::vector<Eigen::VectorXd>> qdot_certificates;
	qdot_certificates.resize(testSimTypes.size());

	// save location of single leaf at the bottom of the tree
	std::vector<std::vector<Eigen::VectorXd>> leaf_pos_certificates;
	leaf_pos_certificates.resize(testSimTypes.size());

	// average total accum error per step
	std::vector<std::vector<double>> avg_step_error;
	avg_step_error.resize(testSimTypes.size());

	// average num iterations of CG (-1 if N\A)
	std::vector<std::vector<double>> CG_iterations;
	CG_iterations.resize(testSimTypes.size());

	// total time trackers
	std::vector<std::vector<double>> total_time;
	total_time.resize(testSimTypes.size());

	// solver time trackers
	std::vector<std::vector<double>> solve_time;
	solve_time.resize(testSimTypes.size());

	// temporary trackers
	double accum_error;
	double accum_iterations;
	double accum_total_time;
	double accum_solve_time;

	std::vector<int> num_links;
	std::vector<int> r_dofs;
	std::vector<int> r_r;
	std::vector<int> m_dofs;
	for (int n = 0; n < testSizes.size(); ++n)
	{
		int n_save = n;
		n = testSizes[n];

		for (int i = 0; i < testSimTypes.size(); ++i)
		{
			Eigen::VectorXd result;
			int count = 0;
			accum_error = 0;
			accum_iterations = 0;
			accum_total_time = 0;
			accum_solve_time = 0;

			rigid_body->reset();
			rigid_body->setSolveType(testSimTypes[i]);
			//rigid_body->loadUmbrella(n);
			//rigid_body->loadBridge(n, n);
			//rigid_body->loadSimpleBridge(n);
			//rigid_body->loadTree(n);
			rigid_body->loadSimpleTree(n);
			///rigid_body->loadTest(1);
			if (i == 0)
			{
				num_links.push_back(rigid_body->getNumLinks());
				r_dofs.push_back(rigid_body->getRedmaxDOFs());
				m_dofs.push_back(rigid_body->getMaximalDOFs());
				r_r.push_back(rigid_body->getRedmaxVsizes());
			}

			if (render)
			{
				brender->setExportDir("../brenders");
				rigid_body->export_part = 0;
				brender->exportBrender(t);
				rigid_body->export_part = 1;
			}

			if (avg_step_lstsq)
				avg_step_error[i].push_back(0);
			while (t < 0.1)
			{
				if (std::abs(std::fmod(t, 1.0)) < 1e-2)
					std::cout << int(t) << " ";

				auto start = std::chrono::steady_clock::now();

				result = rigid_body->step();

				auto end = std::chrono::steady_clock::now();
				auto diff = end - start;
				auto time = std::chrono::duration<double, std::nano>(diff).count();
				accum_total_time += time;

				if (render && torend % 2 == 0)
				{
					brender->exportBrender(t);
				}
				torend++;

				// save all 'ground truth' results
				if (ground_truth)
					raw_result[i].push_back(rigid_body->step());

				accum_solve_time += rigid_body->getSolveTime();

				if(avg_step_lstsq)
					avg_step_error[i][n] += (raw_result[0][n] - result).norm();

				// number of CG iterations (-1 if no CG performed)
				if (iterations)
					accum_iterations = (accum_iterations*count + rigid_body->getCGIterations()) / (double)(count + 1.0);

				++count;
			}
			if (qdot_cert)
				qdot_certificates[i].push_back(result);
			if (leaf_pos_cert)
				leaf_pos_certificates[i].push_back(rigid_body->getLeafCertificate());
			if (avg_step_lstsq)
				avg_step_error[i][n] = avg_step_error[i][n] / (double)count;
			if (iterations)
				CG_iterations[i].push_back(accum_iterations);
			solve_time[i].push_back(accum_solve_time);
			total_time[i].push_back(accum_total_time);
			std::cout << std::endl << std::endl;
			std::cout << "Simulation  " << std::to_string(n) <<
				" complete for " << testNames[i] << std::endl;
		}

		if (render)
		{
			rigid_body->export_part = 2;
			brender->exportBrender(t);
		}

		std::ofstream outfile(filename, std::ios::out);

		outfile << "n_dofs = [\n";
		for (int r = 0; r < r_dofs.size(); r++)
		{
			outfile << std::to_string(r_dofs[r]) << "; ";
			//outfile << std::to_string((std::pow(2, r - 1) - 2) + (std::pow(2, r) - 1)) << "; ";
		}
		outfile << "];\n";

		outfile << "n_r = [\n";
		for (int r = 0; r < r_r.size(); r++)
		{
			outfile << std::to_string(r_r[r]) << "; ";
			//outfile << std::to_string((std::pow(2, r - 1) - 2) + (std::pow(2, r) - 1)) << "; ";
		}
		outfile << "];\n";

		outfile << "max_dofs = [\n";
		for (int r = 0; r < m_dofs.size(); r++)
		{
			outfile << std::to_string(m_dofs[r]) << "; ";
			//outfile << std::to_string((std::pow(2, r - 1) - 2) + (std::pow(2, r) - 1)) << "; ";
		}
		outfile << "];\n";

		outfile << "n_links = [\n";
		for (int r = 0; r < num_links.size(); r++)
		{
			outfile << std::to_string(num_links[r]) << "; ";
			//outfile << std::to_string((std::pow(2, r - 1) - 2) + (std::pow(2, r) - 1)) << "; ";
		}
		outfile << "];\n";

		// total times
		for (int i = 0; i < testSimTypes.size(); ++i)
		{
			outfile << testNames[i] << "_total_time = [\n";
			for (int r = 0; r < total_time[i].size(); r++)
			{
				outfile << total_time[i][r] << "; ";
			}
			outfile << "];\n";
		}

		// solver times
		for (int i = 0; i < testSimTypes.size(); ++i)
		{
			outfile << testNames[i] << "_solve_time = [\n";
			for (int r = 0; r < solve_time[i].size(); r++)
			{
				outfile << solve_time[i][r] << "; ";
			}
			outfile << "];\n";
		}

		// iterations
		if (iterations)
		{
			for (int i = 0; i < testSimTypes.size(); ++i)
			{
				outfile << testNames[i] << "_iterations = [\n";
				for (int r = 0; r < CG_iterations[i].size(); r++)
				{
					outfile << CG_iterations[i][r] << "; ";
				}
				outfile << "];\n";
			}
		}

		// accumulated squared error
		if (avg_step_lstsq)
		{
			for (int i = 0; i < testSimTypes.size(); ++i)
			{
				outfile << testNames[i] << "_error = [\n";
				for (int r = 0; r < avg_step_error[i].size(); r++)
				{
					outfile << avg_step_error[i][r] << "; ";
				}
				outfile << "];\n";
			}
		}

		// certificates
		if (qdot_cert)
		{
			for (int i = 0; i < testSimTypes.size(); ++i)
			{
				outfile << testNames[i] << "_qdot_certificates = [\n";
				for (int r = 0; r < qdot_certificates[i].size(); r++)
				{
					outfile << qdot_certificates[i][r] << "; ";
				}
				outfile << "];\n";
			}
		}

		if (leaf_pos_cert)
		{
			for (int i = 0; i < testSimTypes.size(); ++i)
			{
				outfile << testNames[i] << "_leaf_pos_certificates = [\n";
				for (int r = 0; r < leaf_pos_certificates[i].size(); r++)
				{
					outfile << leaf_pos_certificates[i][r][0] << " "
						<< leaf_pos_certificates[i][r][1] << " "
						<< leaf_pos_certificates[i][r][2] << "; ";
				}
				outfile << "];\n";
			}
		}
		outfile.close();

		n = n_save;
	}
}

void Scene::displayActions()
{
	rigid_body->displayUserActions();
}

void Scene::draw(std::unique_ptr<MatrixStack> &MV, const std::unique_ptr<Program> &prog) const
{
	rigid_body->draw(MV, prog);
}

void Scene::drawJoints(std::unique_ptr<MatrixStack>& MV) const
{
	rigid_body->drawJoints(MV);
}

void Scene::drawPoints(std::unique_ptr<MatrixStack>& MV, const std::unique_ptr<Program>& prog) const
{
	rigid_body->drawPoints(MV, prog);
}


