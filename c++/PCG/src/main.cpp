#include "online/GLSL.h"
#include "online/Program.h"
#include "online/Camera.h"
#include "online/MatrixStack.h"
#include "Scene.h"
#include "ChronoTimer.h"

#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <fstream>

#ifndef _GLIBCXX_USE_NANOSLEEP
#define _GLIBCXX_USE_NANOSLEEP
#endif
#include <thread>

#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "RigidBodyUtility.h"


using namespace std;
using namespace Eigen;

bool keyToggles[256] = {false}; // only for English keyboards!
bool activateDisplay = true;
bool overrideCamera = false;
bool editingSpline = false;
float savedwinz = 0.0;

GLFWwindow *window; // Main application window
string RESOURCE_DIR = ""; // Where the resources are loaded from

shared_ptr<Camera> camera;
unique_ptr<Program> prog;
unique_ptr<Program> progSimple;
unique_ptr<Scene> scene;

static void error_callback(int error, const char *description)
{
	cerr << description << endl;
}

static void key_callback(GLFWwindow *window, int key, int scancode, int action, int mods)
{
	if(key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) {
		glfwSetWindowShouldClose(window, GL_TRUE);
	}
}

static void char_callback(GLFWwindow *window, unsigned int key)
{
	keyToggles[key] = !keyToggles[key];

	// Only allow the stepper function to call step repeatedly if nothing else is going on
	if (key == (unsigned)' ')
		return;
	else
		keyToggles[(unsigned)' '] = false;

	switch(key) {
	case 'a':
		activateDisplay = !activateDisplay;
		if (activateDisplay)
		{
			std::cout << "Showing display" << std::endl;
		}
		else
		{
			std::cout << "Display turned off" << std::endl;
		}
		break;
	default:
		int effect = scene->keyPressed(key);
		if (effect == 0)
		{
			std::cout << "Key press had no effect" << std::endl;
		}
		else if (effect == -2)
		{
			overrideCamera = false;
		}
	}
}

static void cursor_position_callback(GLFWwindow* window, double xmouse, double ymouse)
{
	int state = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT);
	if(!overrideCamera && state == GLFW_PRESS) {
		camera->mouseMoved((float)xmouse, (float)ymouse);
	}
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
	// Get the current mouse position.
	double xmouse, ymouse;
	glfwGetCursorPos(window, &xmouse, &ymouse);

	// Get current window size.
	int width, height;
	glfwGetWindowSize(window, &width, &height);

	// Handle camera action
	if(!overrideCamera && action == GLFW_PRESS) {
		bool shift = mods & GLFW_MOD_SHIFT;
		bool ctrl  = mods & GLFW_MOD_CONTROL;
		bool alt   = mods & GLFW_MOD_ALT;
		camera->mouseClicked(float(xmouse), float(ymouse), shift, ctrl, alt);
	}

	if (overrideCamera && action == GLFW_PRESS)
	{
		auto P = make_unique<MatrixStack>();
		auto MV = make_unique<MatrixStack>();
		P->pushMatrix();
		camera->applyProjectionMatrix(P);
		MV->pushMatrix();
		camera->applyViewMatrix(MV);
		glm::mat4 p = P->topMatrix();
		glm::mat4 v = MV->topMatrix();
		P->popMatrix();
		MV->popMatrix();

		double xwindow = (xmouse - width / 2.0)/(width/2.0);
		double ywindow = ((height - ymouse) - height / 2.0)/(height/2.0);

		float *data = new float[height*width];
		glReadPixels(0, 0, width, height, GL_DEPTH_COMPONENT, GL_FLOAT, data);
		float winz = data[(int)(xmouse + (height - ymouse)*width)];

		glm::vec4 res;
		glm::mat4 pv = inverse(p*v);

		res[0] = (float)xwindow;
		res[1] = (float)ywindow;
		res[2] = (float)(2.0* winz - 1.0);
		res[3] = 1.0f;
		res = pv*res;

		res[3] = 1.0f / res[3];
		res[0] *= res[3];
		res[1] *= res[3];
		res[2] *= res[3];

		int effect = scene->clickPress(res[0], res[1], res[2]);
		savedwinz = res[2];
	}
	else if (overrideCamera && action == GLFW_RELEASE)
	{
		auto P = make_unique<MatrixStack>();
		auto MV = make_unique<MatrixStack>();
		P->pushMatrix();
		camera->applyProjectionMatrix(P);
		MV->pushMatrix();
		camera->applyViewMatrix(MV);
		glm::mat4 p = P->topMatrix();
		glm::mat4 v = MV->topMatrix();
		P->popMatrix();
		MV->popMatrix();

		double xwindow = (xmouse - width / 2.0) / (width / 2.0);
		double ywindow = ((height - ymouse) - height / 2.0) / (height / 2.0);

		glm::vec4 res;
		glm::mat4 pv = inverse(p*v);

		res[0] = (float)xwindow;
		res[1] = (float)ywindow;
		res[2] = (float)(2.0* savedwinz - 1.0);
		res[3] = 1.0f;
		res = pv*res;

		res[3] = 1.0f / res[3];
		res[0] *= res[3];
		res[1] *= res[3];
		res[2] *= res[3];

		int effect = scene->clickRelease(res[0], res[1], res[2]);
		if (effect == -2)
			overrideCamera = false;
	}
}

static void init()
{
	GLSL::checkVersion();
	
	// Set background color
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	// Enable z-buffer test
	glEnable(GL_DEPTH_TEST);
	// Enable alpha blending
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	
	progSimple = make_unique<Program>();
	progSimple->setShaderNames(RESOURCE_DIR + "simple_vert.glsl", RESOURCE_DIR + "simple_frag.glsl");
	progSimple->setVerbose(false); // Set this to true when debugging.
	progSimple->init();
	progSimple->addUniform("P");
	progSimple->addUniform("MV");

	prog = make_unique<Program>();
	prog->setVerbose(false); // Set this to true when debugging.
	prog->setShaderNames(RESOURCE_DIR + "phong_vert.glsl", RESOURCE_DIR + "phong_frag.glsl");
	prog->init();
	prog->addUniform("P");
	prog->addUniform("MV");
	prog->addUniform("kdFront");
	prog->addUniform("kdBack");
	prog->addAttribute("aPos");
	prog->addAttribute("aNor");

	camera = make_shared<Camera>();
	scene = make_unique<Scene>(RESOURCE_DIR);
	scene->init();
	scene->load("input.txt", camera);

	// If there were any OpenGL errors, this will print something.
	// You can intersperse this line in your code to find the exact location
	// of your OpenGL error.
	GLSL::checkError(GET_FILE_LINE);
}

void render()
{
	// Get current frame buffer size.
	int width, height;
	glfwGetFramebufferSize(window, &width, &height);
	glViewport(0, 0, width, height);
	
	// Use the window size for camera.
	glfwGetWindowSize(window, &width, &height);
	camera->setAspect((float)width/(float)height);
	
	// Clear buffers
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	if(keyToggles[(unsigned)'c']) {
		glEnable(GL_CULL_FACE);
	} else {
		glDisable(GL_CULL_FACE);
	}
	if(keyToggles[(unsigned)'l']) {
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	} else {
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	}
	
	auto P = make_unique<MatrixStack>();
	auto MV = make_unique<MatrixStack>();
	
	// Apply camera transforms
	P->pushMatrix();
	camera->applyProjectionMatrix(P);
	MV->pushMatrix();
	camera->applyViewMatrix(MV);

	// Draw grid
	progSimple->bind();
	glUniformMatrix4fv(progSimple->getUniform("P"), 1, GL_FALSE, glm::value_ptr(P->topMatrix()));
	glUniformMatrix4fv(progSimple->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
	glLineWidth(2.0f);
	float x0 = -10.0f;
	float x1 = 10.0f;
	float z0 = -10.0f;
	float z1 = 10.0f;
	int gridSize = 20;
	glLineWidth(1.0f);
	glBegin(GL_LINES);
	for(int i = 1; i < gridSize; ++i) {
		if(i == gridSize/2) {
			glColor3f(0.1f, 0.1f, 0.1f);
		} else {
			glColor3f(0.8f, 0.8f, 0.8f);
		}
		float x = x0 + i / (float)gridSize * (x1 - x0);
		glVertex3f(x, 0.0f, z0);
		glVertex3f(x, 0.0f, z1);
	}
	for(int i = 1; i < gridSize; ++i) {
		if(i == gridSize/2) {
			glColor3f(0.1f, 0.1f, 0.1f);
		} else {
			glColor3f(0.8f, 0.8f, 0.8f);
		}
		float z = z0 + i / (float)gridSize * (z1 - z0);
		glVertex3f(x0, 0.0f, z);
		glVertex3f(x1, 0.0f, z);
	}
	glEnd();
	glColor3f(0.4f, 0.4f, 0.4f);
	glBegin(GL_LINE_LOOP);
	glVertex3f(x0, 0.0f, z0);
	glVertex3f(x1, 0.0f, z0);
	glVertex3f(x1, 0.0f, z1);
	glVertex3f(x0, 0.0f, z1);
	glEnd();

	scene->drawJoints(MV);
	scene->drawPoints(MV, progSimple);
	progSimple->unbind();

	// Draw scene
	prog->bind();
	glUniformMatrix4fv(prog->getUniform("P"), 1, GL_FALSE, glm::value_ptr(P->topMatrix()));
	MV->pushMatrix();
	scene->draw(MV, prog);
	MV->popMatrix();
	prog->unbind();
	
	// Pop stacks
	MV->popMatrix();
	P->popMatrix();
	
	GLSL::checkError(GET_FILE_LINE);
}

void stepperFunc()
{
	while (true)
	{
		if (keyToggles[(unsigned)' '])
		{
			scene->step();
		}

		this_thread::sleep_for(chrono::microseconds(1));
	}
}

int main(int argc, char **argv)
{
	//if(argc < 2) 
	//{
	//	cout << "Please specify the resource directory." << endl;
	//	return 0;
	//}
	////RESOURCE_DIR = argv[1] + string("/");
	//
	RESOURCE_DIR = "../resources" + string("/");

	// Start timer
	ChronoTimer ctime("Main");

#ifndef ONLINE_MODE
	// Run the selected method on dynamically generated scene with n^2 bodies
	// Initialize scene
	try {
		scene = make_unique<Scene>(RESOURCE_DIR);
		scene->init();
	}
	catch (std::exception &e) {
		cerr << "Exception: " << e.what() << endl;
		return -1;
	}

	scene->batchTest();
#else
	// Set error callback.
	glfwSetErrorCallback(error_callback);
	// Initialize the library.
	if (!glfwInit()) {
		return -1;
	}
	// Create a windowed mode window and its OpenGL context.
	window = glfwCreateWindow(853, 640, "Rigid Bodies", NULL, NULL);
	if (!window) {
		glfwTerminate();
		return -1;
	}
	// Make the window's context current.
	glfwMakeContextCurrent(window);
	// Initialize GLEW.
	glewExperimental = true;
	if (glewInit() != GLEW_OK) {
		cerr << "Failed to initialize GLEW" << endl;
		return -1;
	}
	glGetError(); // A bug in glewInit() causes an error that we can safely ignore.
	cout << "OpenGL version: " << glGetString(GL_VERSION) << endl;
	cout << "GLSL version: " << glGetString(GL_SHADING_LANGUAGE_VERSION) << endl;
	// Initialize scene.
	try {
		init();
	}
	catch (std::exception &e) {
		cerr << "Exception: " << e.what() << endl;
		return -1;
	}
	// Set vsync.
	glfwSwapInterval(1);
	// Set keyboard callback.
	glfwSetKeyCallback(window, key_callback);
	// Set char callback.
	glfwSetCharCallback(window, char_callback);
	// Set cursor position callback.
	glfwSetCursorPosCallback(window, cursor_position_callback);
	// Set mouse button callback.
	glfwSetMouseButtonCallback(window, mouse_button_callback);
	// Start simulation thread.
	thread stepperThread(stepperFunc);
	// Loop until the user closes the window.
	while (!glfwWindowShouldClose(window)) {
		if (activateDisplay)
		{
			// Render scene.
			render();
			// Swap front and back buffers.
			glfwSwapBuffers(window);
		}
		// Poll for and process events.
		glfwPollEvents();
		ctime.toc();
	}
	//// Display timing info
	//ctime.print();
	//scene->printTiming();
	//cout << "Press any key to end..." << endl;
	//char charend;
	//cin >> charend;

	// Quit program.
	stepperThread.detach();
	glfwDestroyWindow(window);
	glfwTerminate();
#endif // ONLINE_MODE

	return 0;
}
