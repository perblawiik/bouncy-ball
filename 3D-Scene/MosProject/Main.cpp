#include <glad/glad.h>
#include <glad/glad.c>
#include <GLFW/glfw3.h>

#include "Shader.h"
#include "Matrix.h"
#include "SoftBody.h"
#include "Structs.h"
#include "Mesh.h"

#include <iostream>

/** CONSTANTS **/
const unsigned int SCR_WIDTH = 800;
const unsigned int SCR_HEIGHT = 600;
const GLfloat PI = 3.14159265359f;

/** FUNCTION DECLARATIONS **/
// A callback function on the window that gets called each time the window is resized
void framebuffer_size_callback(GLFWwindow* window, int width, int height);

// Checks if a key is currently being pressed
void processInput(GLFWwindow *window);

// Perspective matrix
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
	GLFWwindow* window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "Bouncy Ball", NULL, NULL);
	if (window == NULL)
	{
		std::cout << "Failed to create GLFW window" << std::endl;
		glfwTerminate();
		return -1;
	}

	glfwMakeContextCurrent(window);
	// Tell GLFW we want to use our resize callback function on every window resize
	glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

	// Initialize GLAD (load all function pointers for OpenGL)
	if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
	{
		std::cout << "Failed to initialize GLAD" << std::endl;
		return -1;
	}

	// Modelview matrix
	GLfloat MV[16] = {0};
	// Position the object in front of the camera
	mat4Translate(MV, 0.0f, 0.0f, -100.0f);

	// Perspective matrix
	GLfloat P[16] = {0};
	mat4Perspective(P, PI / 3, 1.0f, 0.1f, 1000.0f);

	Shader softBodyShader("Shaders//SoftBody//vertex.glsl", "Shaders//SoftBody//fragment.glsl");

	// Settings for the simulation model of the softbody
	Settings settings = {};
	settings.h = 0.01f; // Step
	settings.k = 0.5f; // Spring constant
	settings.b = 0.005f; // Resistance constant
	settings.g = 9.82f; // Gravitation constant
	settings.DIM = 3; // 2-D
	settings.RADIUS = 20.0f;
	settings.WEIGHT = 5.0f; // Total weight of the system

	// Create a sphere mesh
	int numHorizontalSegments = 8; // Horizontal segments for the sphere (number of vertical segments are always twice the number of horizontal segments)
	Mesh sphere;
	sphere.createSphere(numHorizontalSegments, settings.RADIUS);
	
	// Initiate a softbody instance
	SoftBody bouncyBall;

	// Set the mesh of our softbody to the sphere created earlier
	bouncyBall.setMesh(&sphere);
	
	// Set up all the matrices needed for the simulation
	bouncyBall.setupSimulationModel(settings);

	// Wireframe mode
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	// Hide the back face of the triangles
	//glEnable(GL_CULL_FACE);

	// Time variable
	GLfloat time = (GLfloat)glfwGetTime();

	// Pointer for passing particle position array into the shader
	GLfloat* positions = bouncyBall.getParticlePositionArray();

	/** RENDER LOOP **/
	while (!glfwWindowShouldClose(window))
	{
		// Read input
		processInput(window);

		// Rendering commands here
		glClearColor(0.25f, 0.25f, 0.25f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT);

		// Calculate one cycle of the simulation
		bouncyBall.updateSimulationModel(settings);

		// Activate shader
		softBodyShader.use();

		// Update time and pass in to the shader
		time = (GLfloat)glfwGetTime();
		softBodyShader.setFloat("time", time);

		// Model View Matrix
		softBodyShader.setFloatMat4("MV", MV);
		 // Projection Matrix
		softBodyShader.setFloatMat4("P", P);

		// Insert particle positions in shader
		positions = bouncyBall.getParticlePositionArray();
		softBodyShader.setFloat("positions", positions, settings.DIM * settings.NUM_POINTS);

		// Draw object
		bouncyBall.render();
		
		// Swap buffers and check for keyboard input or mouse movement events
		glfwSwapBuffers(window);
		glfwPollEvents();
	}

	// Clean up all the resources and properly exit the application
	glfwTerminate();

	return 0;
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

	M[0] = f / aspect; M[4] = 0.0f; M[8] = 0.0f;                       M[12] = 0.0f;
	M[1] = 0.0f;     M[5] = f;    M[9] = 0.0f;                       M[13] = 0.0f;
	M[2] = 0.0f;     M[6] = 0.0f; M[10] = -(zFar + zNear) / (zFar - zNear); M[14] = -(2 * zNear*zFar) / (zFar - zNear);
	M[3] = 0.0f;     M[7] = 0.0f; M[11] = -1.0f;                      M[15] = 0.0f;
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