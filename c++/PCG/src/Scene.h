#pragma once
#ifndef __Scene__
#define __Scene__

#include <vector>
#include <memory>
#include <string>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <thread>

#include "online/Brender/cpp/BrenderManager.h"

class Camera;
class Program;
class MatrixStack;
class Optimizer;

class RigidBodyMain;

class Scene
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	Scene(const std::string &RESOURCE_DIR);
	~Scene();

	void init();
	void load(const std::string &FILENAME, std::shared_ptr<Camera> &c);

	int keyPressed(unsigned int key);
	int clickPress(float x, float y, float z);
	int clickRelease(float x, float y, float z);

	Eigen::VectorXd step();
	void batchTest();
	void displayActions();

	double getTime() const { return t; };

	// hand-off (main -> rigidbodylinkage) functions
	void draw(std::unique_ptr<MatrixStack> &MV, const std::unique_ptr<Program> &prog) const;
	void drawJoints(std::unique_ptr<MatrixStack> &MV) const;
	void drawPoints(std::unique_ptr<MatrixStack> &MV, const std::unique_ptr<Program> &prog) const;
	
private:
	double t;
	double h;

	std::shared_ptr<RigidBodyMain> rigid_body;
	BrenderManager *brender;
};

#endif
