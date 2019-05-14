#include "Cuboid.h"
#include <iostream>

#include "online/GLSL.h"
#include "online/Program.h"

#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"

using namespace std;
using namespace Eigen;

Cuboid::Cuboid() :
	posBufID(0),
	norBufID(0),
	texBufID(0)
{
}

Cuboid::~Cuboid()
{
}

void Cuboid::loadMesh(const string &meshName)
{
	string cubeMesh = meshName + "cube.obj";
	// Load geometry
	tinyobj::attrib_t attrib;
	std::vector<tinyobj::shape_t> shapes;
	std::vector<tinyobj::material_t> materials;
	string errStr;
	bool rc = tinyobj::LoadObj(&attrib, &shapes, &materials, &errStr, cubeMesh.c_str());
	if(!rc) {
		cerr << errStr << endl;
	} else {
		// Some OBJ files have different indices for vertex positions, normals,
		// and texture coordinates. For example, a cube corner vertex may have
		// three different normals. Here, we are going to duplicate all such
		// vertices.
		// Loop over shapes
		for(size_t s = 0; s < shapes.size(); s++) {
			// Loop over faces (polygons)
			size_t index_offset = 0;
			for(size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++) {
				size_t fv = shapes[s].mesh.num_face_vertices[f];
				// Loop over vertices in the face.
				for(size_t v = 0; v < fv; v++) {
					// access to vertex
					tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];
					posBuf.push_back(attrib.vertices[3*idx.vertex_index+0]);
					posBuf.push_back(attrib.vertices[3*idx.vertex_index+1]);
					posBuf.push_back(attrib.vertices[3*idx.vertex_index+2]);
					posBufOrigin.push_back(attrib.vertices[3*idx.vertex_index+0]);
					posBufOrigin.push_back(attrib.vertices[3*idx.vertex_index+1]);
					posBufOrigin.push_back(attrib.vertices[3*idx.vertex_index+2]);
					if(!attrib.normals.empty()) {
						norBuf.push_back(attrib.normals[3*idx.normal_index+0]);
						norBuf.push_back(attrib.normals[3*idx.normal_index+1]);
						norBuf.push_back(attrib.normals[3*idx.normal_index+2]);
					}
					if(!attrib.texcoords.empty()) {
						texBuf.push_back(attrib.texcoords[2*idx.texcoord_index+0]);
						texBuf.push_back(attrib.texcoords[2*idx.texcoord_index+1]);
					}
				}
				index_offset += fv;
				// per-face material (IGNORE)
				shapes[s].mesh.material_ids[f];
			}
		}
	}
	// Cube w, h, d
	Xdim = 2.0;
	Ydim = 2.0;
	Zdim = 2.0;
}

void Cuboid::size(const Eigen::Vector3d &size)
{
	for (unsigned int i = 0; i < posBufOrigin.size() / 3; i++) {
		posBufOrigin[i * 3] = (float)(posBufOrigin[i * 3]* size[0]/Xdim);
		posBufOrigin[i * 3 + 1] = (float)(posBufOrigin[i * 3 + 1]* size[1]/Ydim);
		posBufOrigin[i * 3 + 2] = (float)(posBufOrigin[i * 3 + 2]* size[2]/Zdim);
	}
	Xdim = (float)size[0];
	Ydim = (float)size[1];
	Zdim = (float)size[2];
}

void Cuboid::init()
{
	// Send the position array to the GPU
	if (!posBuf.empty()) {
		glGenBuffers(1, &posBufID);
		glBindBuffer(GL_ARRAY_BUFFER, posBufID);
		glBufferData(GL_ARRAY_BUFFER, posBuf.size()*sizeof(float), &posBuf[0], GL_DYNAMIC_DRAW);
	}
	
	// Send the normal array to the GPU
	if(!norBuf.empty()) {
		glGenBuffers(1, &norBufID);
		glBindBuffer(GL_ARRAY_BUFFER, norBufID);
		glBufferData(GL_ARRAY_BUFFER, norBuf.size()*sizeof(float), &norBuf[0], GL_STATIC_DRAW);
	}
	
	// Send the texture array to the GPU
	if(!texBuf.empty()) {
		glGenBuffers(1, &texBufID);
		glBindBuffer(GL_ARRAY_BUFFER, texBufID);
		glBufferData(GL_ARRAY_BUFFER, texBuf.size()*sizeof(float), &texBuf[0], GL_STATIC_DRAW);
	}
	
	// Unbind the arrays
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	
	GLSL::checkError(GET_FILE_LINE);
}

void Cuboid::updatePosBuf(const Matrix4d & E)
{
	for (unsigned int i = 0; i < posBufOrigin.size()/3; i++) {
		Eigen::Vector4d x = Eigen::Vector4d(posBufOrigin[i*3], posBufOrigin[i*3 + 1], posBufOrigin[i*3 + 2], 1.0);
		Eigen::Vector4d xnew = E*x;
		posBuf[i*3] = (float)xnew[0];
		posBuf[i*3 + 1] = (float)xnew[1];
		posBuf[i*3 + 2] = (float)xnew[2];
	}
 }

void Cuboid::draw(const unique_ptr<Program> &prog) const
{
	GLSL::checkError(GET_FILE_LINE);
	// Bind position buffer
	int h_pos = prog->getAttribute("aPos");
	if (h_pos != -1 && posBufID != 0) {
		glEnableVertexAttribArray(h_pos);
		glBindBuffer(GL_ARRAY_BUFFER, posBufID);
		glBufferData(GL_ARRAY_BUFFER, posBuf.size() * sizeof(float), &posBuf[0], GL_DYNAMIC_DRAW);
		glVertexAttribPointer(h_pos, 3, GL_FLOAT, GL_FALSE, 0, (const void *)0);
	}
	
	// Bind normal buffer
	int h_nor = prog->getAttribute("aNor");
	if(h_nor != -1 && norBufID != 0) {
		glEnableVertexAttribArray(h_nor);
		glBindBuffer(GL_ARRAY_BUFFER, norBufID);
		glVertexAttribPointer(h_nor, 3, GL_FLOAT, GL_FALSE, 0, (const void *)0);
	}
	
	// Bind texcoords buffer
	int h_tex = prog->getAttribute("aTex");
	if(h_tex != -1 && texBufID != 0) {
		glEnableVertexAttribArray(h_tex);
		glBindBuffer(GL_ARRAY_BUFFER, texBufID);
		glVertexAttribPointer(h_tex, 2, GL_FLOAT, GL_FALSE, 0, (const void *)0);
	}
	
	// Draw
	int count = (int)(posBuf.size()/3); // number of indices to be rendered
	glDrawArrays(GL_TRIANGLES, 0, count);
	
	// Disable and unbind
	if(h_tex != -1) {
		glDisableVertexAttribArray(h_tex);
	}
	if(h_nor != -1) {
		glDisableVertexAttribArray(h_nor);
	}
	glDisableVertexAttribArray(h_pos);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	
	GLSL::checkError(GET_FILE_LINE);
}

std::vector<Eigen::Vector3d> Cuboid::cornerPos(const Eigen::Matrix4d& E)
{
	std::vector<Eigen::Vector3d> corners;
	
	for (int x = -1; x < 2; x+=2)
		for (int y = -1; y < 2; y+=2)
			for (int z = -1; z < 2; z+=2) 
			{
				Eigen::Vector4d curCorner = Eigen::Vector4d(x*Xdim/2, y*Ydim/2, z*Zdim/2, 1.0);
				curCorner = E*curCorner;
				corners.push_back(Eigen::Vector3d(curCorner[0], curCorner[1], curCorner[2]));
			}
	return corners;
}
