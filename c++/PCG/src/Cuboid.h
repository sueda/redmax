#pragma once
#ifndef _CUBOID_H_
#define _CUBOID_H_

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <string>
#include <vector>
#include <memory>

class Program;

typedef Eigen::Matrix<double, 6, 6> Matrix6d;

/**
 * A shape defined by a list of triangles
 * - posBuf should be of length 3*ntris
 * - norBuf should be of length 3*ntris (if normals are available)
 * - texBuf should be of length 2*ntris (if texture coords are available)
 * posBufID, norBufID, and texBufID are OpenGL buffer identifiers.
 */
class Cuboid
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	Cuboid();
	virtual ~Cuboid();
	void loadMesh(const std::string &meshName);
	void size(const Eigen::Vector3d &size);
	void init();
	void updatePosBuf(const Eigen::Matrix4d& E);
	void draw(const std::unique_ptr<Program> &prog) const;

//	Matrix6d spatialInertia(double mass);
	std::vector<Eigen::Vector3d> cornerPos(const Eigen::Matrix4d& E);
	
private:
	float Xdim;
	float Ydim;
	float Zdim;

	std::vector<float> posBuf;
	std::vector<float> posBufOrigin;
	std::vector<float> norBuf;
	std::vector<float> texBuf;
	unsigned posBufID;
	unsigned norBufID;
	unsigned texBufID;
};

#endif
