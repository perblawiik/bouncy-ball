#include <glad/glad.h>
#include <glad/glad.c>
#include <GLFW/glfw3.h>

#include "Shader.h"
#include "Matrix.h"
#include "SoftBody.h"
#include "Structs.h"
#include "Mesh.h"

#include <iostream>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <vector>

/** CONSTANTS **/
const unsigned int WINDOW_WIDTH = GLOBAL_CONSTANTS::window::WIDTH;
const unsigned int WINDOW_HEIGHT = GLOBAL_CONSTANTS::window::HEIGHT;
const GLfloat PI = GLOBAL_CONSTANTS::PI;

/** FUNCTION DECLARATIONS **/
// A callback function on the window that gets called each time the window is resized
void framebuffer_size_callback(GLFWwindow* window, int width, int height);

// Checks if a key is currently being pressed
void processInput(GLFWwindow *window);

// Load simulation with specified number of steps and store all particle positions in a single matrix.
void createSimulation(GLFWwindow* window, SoftBody &bouncyBall, Matrix &PARTICLE_POSITION_DATA, Settings &settings);
// Render one cycle of the simulation, if step counter is larger than total steps, the step counter is reset to 1
void renderSimulationStep(int &stepCounter, const int NUM_STEPS, Settings &settings, Matrix &PARTICLE_POSITION_DATA, SoftBody &softBody);
// Save simulation data to a text file
void saveSimulationToTextFile(SoftBody &sb, Matrix &DATA, const Settings &simulationSettings);

// Matrix functions
void mat4Perspective(GLfloat M[], const GLfloat &vertFov, const GLfloat &aspect, const GLfloat &zNear, const GLfloat &zFar);
void mat4Translate(GLfloat M[], const GLfloat &x, const GLfloat &y, const GLfloat &z);
void mat4Identity(GLfloat M[]);

int main() 
{
	// Initialize glfw
	glfwInit();

	// Tell GLFW that we want to use OpenGL version 3.3
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);

	// Tell GLFW we want to explicitly use the core-profile
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

#ifdef __APPLE__
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // Uncomment this statement to fix compilation on OS X
#endif


	// Create window object
	GLFWwindow* window = glfwCreateWindow(WINDOW_WIDTH, WINDOW_HEIGHT, "Bouncy Ball", NULL, NULL);
	if (window == NULL)
	{
		std::cout << "Failed to create GLFW window" << std::endl;
		glfwTerminate();
		return -1;
	}

	glfwMakeContextCurrent(window);
	// Tell GLFW we want to use our resize callback function on every window resize
	glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
	glfwSetWindowAspectRatio(window, 4, 3);

	// Initialize GLAD (load all function pointers for OpenGL)
	if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
	{
		std::cout << "Failed to initialize GLAD" << std::endl;
		return -1;
	}

	// Modelview matrix
	GLfloat MV[16] = {0};
	// Position the object in front of the camera
	mat4Translate(MV, 0.0f, 0.0f, -50.0f);

	// Perspective matrix
	GLfloat P[16] = {0};
	mat4Perspective(P, PI / 3, 1.0f, 0.1f, 1000.0f); // Field of view is set to PI/3 (60 degrees)

	// Shader instance for our bouncy ball
	Shader mainShader("Shaders//vertex.glsl", "Shaders//fragment.glsl");
	// Activate shader
	mainShader.use();
	// Insert Projection Matrix
	mainShader.setFloatMat4("P", P);

	// Settings for the simulation model of the softbody
	Settings settings = {};
	settings.h = 0.01f; // Time step size
	settings.k = 0.1f; // Spring constant
	settings.b = settings.k/50.0f; // Resistance constant
	settings.g = 9.82f; // Gravitation constant
	settings.DIM = 3; // 3-D (x,y,z)
	settings.RADIUS = 10.0f;
	settings.WEIGHT = 1.0f; // Total weight of the system
	settings.TIME_DURATION = 1.0f; // Specifies how long the simulation should be (given in seconds)
	settings.NUM_STEPS = (int)(settings.TIME_DURATION / settings.h); // Specifies how many steps the simulation will be calculated
	std::cout << "Number of simulation steps: " << settings.NUM_STEPS << std::endl;

	// Create a sphere mesh (for the softbody)
	int numHorizontalSegments = 16 ; // Horizontal segments for the sphere (number of vertical segments are always twice the number of horizontal segments)
	Mesh sphere;
	sphere.createSphere(numHorizontalSegments, settings.RADIUS);
	
	// The bouncy ball is a softbody simulation
	SoftBody bouncyBall;
	// Set the mesh of our softbody to the sphere created earlier
	bouncyBall.setMesh(&sphere);
	// Set up all the matrices needed for the simulation
	bouncyBall.setupSimulationModel(settings, settings.RADIUS);

	// Create a XZ-plane as floor
	Mesh floor;
	floor.createPlaneXZ(100.0f, 100.0f);

	// A matrix to store all simulation cycles in
	Matrix PARTICLE_POSITION_DATA(settings.NUM_STEPS, settings.NUM_POINTS*settings.DIM);

	// Compute the entire simulation with specified number of steps
	createSimulation(window, bouncyBall, PARTICLE_POSITION_DATA, settings);

	// Save simulation data to a textfile with given name as parameter
	saveSimulationToTextFile(bouncyBall, PARTICLE_POSITION_DATA, settings);

	// Wireframe mode
	//glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	// Hide the back face of the triangles
	glEnable(GL_CULL_FACE);
	glClearColor(0.25f, 0.25f, 0.25f, 1.0f);

	// Time variables
	GLfloat time = (GLfloat)glfwGetTime();
	GLfloat deltaTime = 0.0f;

	// Counter to keep track of which simulation step to render
	int stepCounter = 1;

	/** RENDER LOOP **/
	while (!glfwWindowShouldClose(window))
	{
		// Read input
		processInput(window);

		// Rendering commands here
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		// Activate shader
		mainShader.use();

		// Place the floor
		mat4Translate(MV, 0.0f, -settings.RADIUS*2.0f, -50.0f);
		// Update with new model view matrix
		mainShader.setFloatMat4("MV", MV);
		// Draw floor
		floor.render();

		// Place sphere infront of the camera
		mat4Translate(MV, 0.0f, 0.0f, -50.0f);
		// Update with new model view matrix
		mainShader.setFloatMat4("MV", MV);	
		// Render (draw) the simulation one step at the time
		renderSimulationStep(
			stepCounter, // Keeps count of which simulation step to draw
			settings.NUM_STEPS, // Number of steps of the entire simulation
			settings, // Struct that contains specifications about the simulation
			PARTICLE_POSITION_DATA, // A matrix that stores all simulation steps
			bouncyBall // SoftBody simulation model (contains all information about the simulation)
		);

		// Update current time
		GLfloat currentTime = (GLfloat)glfwGetTime();
		// Calculate time passed since last render
		deltaTime = currentTime - time;
		// Update time
		time = currentTime;

		// Hold each frame for a specified number of seconds (its based on step size of the simulation right now)
		while ((currentTime - time) < (settings.h - deltaTime)) {
			currentTime = (GLfloat)glfwGetTime();
		}
		
		// Swap buffers and check for keyboard input or mouse movement events
		glfwSwapBuffers(window);
		glfwPollEvents();
	}

	// Clean up all the resources and properly exit the application
	glfwTerminate();

	return 0;
}


void createSimulation(GLFWwindow* window, SoftBody &softBody, Matrix &PARTICLE_POSITION_DATA, Settings &settings)
{
	// Title string to display loading progress
	static char titleString[200];

	// Run the simulation for specified number of steps
	for (int i = 1; i <= settings.NUM_STEPS; ++i) {
		// Clear screen
		glClear(GL_COLOR_BUFFER_BIT);

		// Calculate one cycle of the simulation
		softBody.updateSimulationModel(settings);

		// Store the calculated step as a row in a matrix
		for (int j = 0; j < settings.NUM_POINTS*settings.DIM; ++j) {
			PARTICLE_POSITION_DATA(i, j + 1) = softBody.getParticlePositionArray()[j];
		}

		// Compute the percentage of completion
		float progress = 100.0f * (float)i / (float)settings.NUM_STEPS;
		// Compose a string with loading progress message
		sprintf_s(titleString, "Bouncy Ball (Creating Simulation: %.1f %%)", progress);
		// Display loading progress in the window title
		glfwSetWindowTitle(window, titleString);

		// Swap buffers and check for keyboard input or mouse movement events
		glfwSwapBuffers(window);
		glfwPollEvents();
	}

	// Set title string to loading complete
	sprintf_s(titleString, "Bouncy Ball (Simulation Complete)");
	// Display loading progress in the window title
	glfwSetWindowTitle(window, titleString);
}

void renderSimulationStep(int &stepCounter, const int NUM_STEPS, Settings &settings, Matrix &PARTICLE_POSITION_DATA, SoftBody &softBody )
{
	// Pointer to the vertex array of the mesh
	GLfloat* vertexArray = softBody.getMeshVertexArray();

	// Update the vertex array of the mesh from the simulated position data
	for (int i = 0; i < settings.NUM_POINTS; ++i) {
		vertexArray[i * 8] = PARTICLE_POSITION_DATA(stepCounter, (i * 3) + 1);
		vertexArray[(i * 8) + 1] = PARTICLE_POSITION_DATA(stepCounter, (i * 3) + 2);
		vertexArray[(i * 8) + 2] = PARTICLE_POSITION_DATA(stepCounter, (i * 3) + 3);
	}

	// Update vertex buffer with new vertex array
	softBody.updatePositions();

	// Draw object
	softBody.render();

	// Incremet the counter to the next simulation step, if the counter has reached maximum number of steps, reset counter
	++stepCounter;
	if (stepCounter >= NUM_STEPS) {
		stepCounter = 1;
	}
	vertexArray = nullptr;
}

void saveSimulationToTextFile(SoftBody &sb, Matrix &DATA, const Settings &simulationSettings)
{
	// Out File Stream
	std::ofstream outFile;

	// Vertices
	float* vertices = DATA.getValues();
	std::ostream_iterator<float> float_out_it(outFile, ",");
	outFile.open("SimulationData//vertex_data.txt", std::ofstream::out);
	std::copy(vertices, vertices + DATA.size(), float_out_it);
	/*
	long pos = outFile.tellp();
	outFile.seekp(pos-1);
	outFile.write("\0", 1);
	*/
	outFile.close();

	// Indices
	GLuint* indices = sb.getMeshIndexArray();
	std::ostream_iterator<GLuint> uint_out_it(outFile, ",");
	outFile.open("SimulationData//index_data.txt", std::ofstream::out);
	std::copy(indices, indices + (sb.getNumMeshVertices() * 3), float_out_it);
	outFile.close();

}

void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
	glViewport(0, 0, width, height);
}

void processInput(GLFWwindow* window)
{
	// If user press the escape key, close GLFW
	if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
		glfwSetWindowShouldClose(window, true);
}

// M is the matrix we want to create (an output argument )
// vertFov is the vertical field of view (in the y direction )
// aspect is the aspect ratio of the viewport ( width / height )
// zNear is the distance to the near clip plane ( znear > 0)
// zFar is the distance to the far clip plane ( zfar > znear )
void mat4Perspective(GLfloat M[], const GLfloat &vertFov, const GLfloat &aspect, const GLfloat &zNear, const GLfloat &zFar) {

	GLfloat f = cos(vertFov / 2) / sin(vertFov / 2);

	M[0] = f / aspect; M[4] = 0.0f; M[8] = 0.0f;                              M[12] = 0.0f;
	M[1] = 0.0f;       M[5] = f;    M[9] = 0.0f;                              M[13] = 0.0f;
	M[2] = 0.0f;       M[6] = 0.0f; M[10] = -(zFar + zNear) / (zFar - zNear); M[14] = -(2 * zNear*zFar) / (zFar - zNear);
	M[3] = 0.0f;       M[7] = 0.0f; M[11] = -1.0f;                            M[15] = 0.0f;
}

void mat4Translate(GLfloat M[], const GLfloat &x, const GLfloat &y, const GLfloat &z) {

	M[0] = 1.0f; M[4] = 0.0f; M[8] = 0.0f; M[12] = x;
	M[1] = 0.0f; M[5] = 1.0f; M[9] = 0.0f; M[13] = y;
	M[2] = 0.0f; M[6] = 0.0f; M[10] = 1.0f; M[14] = z;
	M[3] = 0.0f; M[7] = 0.0f; M[11] = 0.0f; M[15] = 1.0f;
}

void mat4Identity(GLfloat M[]) {

	M[0] = 1.0f; M[4] = 0.0f; M[8] = 0.0f;  M[12] = 0.0f;
	M[1] = 0.0f; M[5] = 1.0f; M[9] = 0.0f;  M[13] = 0.0f;
	M[2] = 0.0f; M[6] = 0.0f; M[10] = 1.0f;  M[14] = 0.0f;
	M[3] = 0.0f; M[7] = 0.0f; M[11] = 0.0f;  M[15] = 1.0f;
}