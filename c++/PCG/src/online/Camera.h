#pragma  once
#ifndef __Camera__
#define __Camera__

#include <memory>

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>

class MatrixStack;

class Camera
{
public:
	enum {
		ROTATE = 0,
		TRANSLATE,
		SCALE
	};
	
	Camera();
	virtual ~Camera();
	void setInitDistance(const float z) { translations.z = -std::abs(z); }
	void setAspect(const float a) { aspect = a; };
	void setFovy(const float f) { fovy = f; };
	void setZnear(const float z) { znear = z; };
	void setZfar(const float z) { zfar = z; };
	void setRotationFactor(const float f) { rfactor = f; };
	void setTranslationFactor(const float f) { tfactor = f; };
	void setScaleFactor(const float f) { sfactor = f; };
	void mouseClicked(const float x, const float y, const bool shift, const bool ctrl, const bool alt);
	void mouseMoved(const float x, const float y);
	void applyProjectionMatrix(const std::unique_ptr<MatrixStack> &P) const;
	void applyViewMatrix(const std::unique_ptr<MatrixStack> &MV) const;

	void setTranslations(const float x, const float y, const float z);

	float getZnear() { return znear; };
	float getZfar() { return zfar; };
	
private:
	float aspect;
	float fovy;
	float znear;
	float zfar;
	glm::vec2 rotations;
	glm::vec3 translations;
	glm::vec2 mousePrev;
	int state;
	float rfactor;
	float tfactor;
	float sfactor;
	float winz;
};

#endif
