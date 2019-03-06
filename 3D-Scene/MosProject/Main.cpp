#include <glad/glad.h>
#include <glad/glad.c>
#include <GLFW/glfw3.h>

#include "Shader.h"
#include "Matrix.h"
#include "SoftBody.h"
#include "Structs.h"
#include "Mesh.h"
#include "Camera.h"
#include "CameraController.h"
#include "Texture.h"

#include <iostream>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <vector>
#include <string>
#include <sstream>

/** CONSTANTS **/
const GLfloat PI = GLOBAL_CONSTANTS::PI;

/** FUNCTION DECLARATIONS **/
// A callback function on the window that gets called each time the window is resized
void framebufferSizeCallback(GLFWwindow* window, int width, int height);
// Checks if a key is currently being pressed
void processInput(GLFWwindow* window);
// Get IDs for the uniform locations in shader (prints warning message to console if not found)
void getUniformLocationIDs(Shader &shader, GLint &MV, GLint &CV, GLint &P, GLint &color, GLint &cameraPos);

// Load simulation with specified number of steps and store all particle positions in a single matrix.
void createSimulation(GLFWwindow* window, SoftBody &bouncyBall, Matrix &PARTICLE_POSITION_DATA, Settings &settings);
// Render one cycle of the simulation, if step counter is larger than total steps, the step counter is reset to 1
void renderSimulationStep(int &stepCounter, const int NUM_STEPS, Settings &settings, Matrix &PARTICLE_POSITION_DATA, SoftBody &softBody);
// Save simulation data to a text file
void saveSimulationSequence(SoftBody &sb, Matrix &DATA, const Settings &simulationSettings, const std::string &simulationName);


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


	// Activate MSAA (Multisampling Anti Aliasing)
	glfwWindowHint(GLFW_SAMPLES, 4);

	int width = GLOBAL_CONSTANTS::window::WIDTH;
	int height = GLOBAL_CONSTANTS::window::HEIGHT;

	// Create window object
	GLFWwindow* window = glfwCreateWindow(width, height, "Bouncy Ball", NULL, NULL);
	if (window == NULL)
	{
		std::cout << "Failed to create GLFW window" << std::endl;
		glfwTerminate();
		return -1;
	}

	glfwMakeContextCurrent(window);
	// Tell GLFW we want to use our resize callback function on every window resize
	glfwSetFramebufferSizeCallback(window, framebufferSizeCallback);
	glfwSetWindowAspectRatio(window, 16, 9);

	// Initialize GLAD (load all function pointers for OpenGL)
	if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
	{
		std::cout << "Failed to initialize GLAD" << std::endl;
		return -1;
	}

	// Set window to show black background
	glClear(GL_COLOR_BUFFER_BIT);
	glfwSwapBuffers(window);
	
	//glPolygonMode(GL_FRONT_AND_BACK, GL_LINE); // Wireframe mode
	glEnable(GL_CULL_FACE); // Hide the back face of the triangles
	glClearColor(0.25f, 0.25f, 0.25f, 1.0f); // Set the color for every clear call
	glEnable(GL_DEPTH_TEST); // Enable depth so that objects closest to the camera view is visible
	glEnable(GL_MULTISAMPLE);
	
	// Create our main shader (responsible for all the lighting in the scene)
	Shader mainShader("Shaders//vertex.glsl", "Shaders//fragment.glsl");
	// Activate shader
	mainShader.use();

	// Location IDs for shader uniforms
	GLint modelViewMatrixLocationID; // Model View Matrix
	GLint cameraViewMatrixLocationID; // Camera View Matrix
	GLint perspectiveMatrixLocationID; // Perspective Matrix
	GLint colorLocationID; // Surface color for objects
	GLint cameraPositionLocationID; // Camera position
	
	// Get uniform location IDs from the shader (prints warning message to console if not found)
	getUniformLocationIDs(
		mainShader, 
		modelViewMatrixLocationID, 
		cameraViewMatrixLocationID, 
		perspectiveMatrixLocationID, 
		colorLocationID,
		cameraPositionLocationID
	);

	// Modelview matrix
	GLfloat modelViewMatrix[16] = { 0 };
	MATRIX4::identity(modelViewMatrix); // Identity matrix as default
	mainShader.setFloatMat4(modelViewMatrixLocationID, modelViewMatrix); // Insert Model View Matrix

	// Perspective projection matrix
	GLfloat perspectiveMatrix[16] = { 0 };
	MATRIX4::perspective(perspectiveMatrix, PI / 3, (GLfloat)width/(GLfloat)height, 0.1f, 1000.0f); // Field of view is set to PI/3 (60 degrees)
	mainShader.setFloatMat4(perspectiveMatrixLocationID, perspectiveMatrix); // Insert Projection Matrix


	// Lamp shader
	Shader lampShader("Shaders//lampVert.glsl", "Shaders//lampFrag.glsl");
	lampShader.use();
	lampShader.setFloatMat4("perspective", perspectiveMatrix); // Set perspective projection matrix

	// Position of the lamp
	GLfloat lightPosition[] = {
		0.0f, 75.0f, 0.0f
	};

	// Light color of the lamp
	GLfloat lightColor[] = {
		1.0f, 1.0f, 1.0f
	};

	lampShader.setVec3("objectColor", lightColor[0], lightColor[1], lightColor[2]); // Set color

	// Create a sphere mesh for the lamp
	Mesh lamp;
	lamp.createSphere(8, 10.0f);

	// Set light position and color in the main shader
	mainShader.use();
	mainShader.setVec3("lightPosition", lightPosition[0], lightPosition[1], lightPosition[2]);
	mainShader.setVec3("lightColor", lightColor[0], lightColor[1], lightColor[2]);

	// Create a camera for the scene
	Camera camera(&mainShader);
	// Back the camera away from the origin 
	camera.setPosition(0.0f, 0.0f, 100.0f); // (x, y, z)
	// Add lamp shader to the camera view
	camera.addShader(&lampShader);

	// Create a camera controller
	CameraController controller(&camera, window);
	controller.setMouseSensitivity(20.0f);
	
	// Settings for the simulation model of the softbody
	Settings settings = {};
	settings.h = 0.01f; // Time step size
	settings.k = 0.1f; // Spring constant
	settings.b = settings.k/50.0f; // Resistance constant
	settings.g = 9.82f; // Gravitation constant
	settings.DIM = 3; // 3-D (x,y,z)
	settings.RADIUS = 10.0f;
	settings.WEIGHT = 0.75f; // Total weight of the system
	settings.TIME_DURATION = 10.0f; // Specifies how long the simulation should be (given in seconds)
	settings.NUM_STEPS = (int)(settings.TIME_DURATION / settings.h); // Specifies how many steps the simulation will be calculated
	settings.NUM_SEGS = 16; // Horizontal segments for the sphere (number of vertical segments are always twice the number of horizontal segments)
	std::cout << "Number of simulation steps: " << settings.NUM_STEPS << std::endl;

	/*
	// Create a sphere mesh (for the softbody)
	Mesh sphere;
	sphere.createSphere(settings.NUM_SEGS, settings.RADIUS);
	
	// The bouncy ball is a softbody simulation
	SoftBody bouncyBall;
	// Set the mesh of our softbody to the sphere created earlier
	bouncyBall.setMesh(&sphere);
	// Set up all the matrices needed for the simulation
	bouncyBall.setupSimulationModel(settings, settings.RADIUS);

	// A matrix to store all simulation cycles in
	Matrix PARTICLE_POSITION_DATA(settings.NUM_STEPS, settings.NUM_POINTS*6);
	// Compute the entire simulation with specified number of steps
	createSimulation(window, bouncyBall, PARTICLE_POSITION_DATA, settings);
	// Save simulation data to a textfile with given name as parameter
	saveSimulationSequence(bouncyBall, PARTICLE_POSITION_DATA, settings, "BouncyBall_01");
	*/

	// Create a XZ-plane as floor
	Mesh floor;
	floor.createPlaneXZ(600.0f, 600.0f);
	
	Mesh ball;
	ball.loadMeshData("BouncyBall_back_left"); // All animations have the same mesh data
	ball.addAnimation("BouncyBall_back_left"); // ID: 0
	/*
	ball.addAnimation("BouncyBall_back_right"); // ID: 1
	ball.addAnimation("BouncyBall_front_left"); // ID: 2
	ball.addAnimation("BouncyBall_straight_up"); // ID: 3
	ball.addAnimation("BouncyBall_front_right"); // ID: 4
	*/

	// Load textures
	Texture woodenFloorTexture("Files//Textures//wooden_floor.jpg");
	Texture footballTexture("Files//Textures//football.jpg");


	/** SHADOW MAPPING STUFF **/
	/*
	// Framebuffer object for rendering the depth map
	unsigned int depthMapFBO;
	glGenFramebuffers(1, &depthMapFBO);

	const unsigned int SHADOW_WIDTH = 1024, SHADOW_HEIGHT = 1024;

	// Create a 2D texture as depth buffer
	unsigned int depthMap;
	glGenTextures(1, &depthMap);
	glBindTexture(GL_TEXTURE_2D, depthMap);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT, SHADOW_WIDTH, SHADOW_HEIGHT, 0, GL_DEPTH_COMPONENT, GL_FLOAT, NULL);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

	// Attach the depth texture as the framebuffer's depth buffer
	glBindFramebuffer(GL_FRAMEBUFFER, depthMapFBO);
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, depthMap, 0);
	glDrawBuffer(GL_NONE);
	glReadBuffer(GL_NONE);
	glBindFramebuffer(GL_FRAMEBUFFER, 0);

	// Use an orthographic projection matrix for the light source (no perspective needed since the light rays are parallell)
	GLfloat nearPlane = 1.0f, farPlane = 7.5f;
	GLfloat lightProjection[16] = { 0.0f };
	MATRIX4::orthogonal(
		lightProjection, // Matrix
		-10.0f, 10.0f, -10.0f, 10.0f, // Left, right, down, up
		nearPlane, farPlane
	);

	// Create a view matrix to transform each object so they're visible from the light's point of view
	GLfloat lightView[16] = { 0.0f };
	GLfloat up[] = { 0.0f, 1.0f, 0.0f };
	GLfloat center[] = { 0.0f, 0.0f, 0.0f };
	MATRIX4::lookAt(
		lightView,
		lightPosition,
		center,
		up
	);

	// Create a light space matrix that transforms world-space vertices into the space as visible from the light source
	GLfloat lightSpaceMatrix[16] = { 0.0f };
	MATRIX4::multiply(lightProjection, lightView, lightSpaceMatrix);
	
	// Create shader for rendering depth map
	Shader depthMapShader("Shaders//depthMapVert.glsl", "Shaders//depthMapFrag.glsl");
	depthMapShader.use();
	depthMapShader.setFloatMat4("lightSpaceMatrix", lightSpaceMatrix);
	depthMapShader.setFloatMat4("modelView", modelViewMatrix);
	*/

	/*********************/
	
	// Time variables
	GLfloat time = (GLfloat)glfwGetTime();
	GLfloat deltaTime = 0.0f;

	// Counter to keep track of which simulation step to render
	int stepCounter = 1;

	/** RENDER LOOP **/
	while (!glfwWindowShouldClose(window))
	{
		// Update current time
		GLfloat startTime = (GLfloat)glfwGetTime();

		// Calculate time passed since last render
		deltaTime = startTime - time;
		// Update time
		time = startTime;

		// Read inputs
		processInput(window); 
		controller.processInput(deltaTime);


		// Render to depth map
		/*
		depthMapShader.use();
		glViewport(0, 0, SHADOW_WIDTH, SHADOW_HEIGHT);
		glBindFramebuffer(GL_FRAMEBUFFER, depthMapFBO);
		glClear(GL_DEPTH_BUFFER_BIT);

		// Render scene
		// FLOOR
		MATRIX4::translate(modelViewMatrix, 0.0f, -settings.RADIUS*2.0f, 0.0f); // Place the floor
		depthMapShader.setFloatMat4("modelViewMatrix", modelViewMatrix); // Update with new model view matrix
		floor.render(); // Draw

		// ANIMATION 0
		MATRIX4::translate(modelViewMatrix, -settings.RADIUS*2.0f, 0.0f, -(settings.RADIUS*2.0f)); // Position the sphere
		depthMapShader.setFloatMat4("modelViewMatrix", modelViewMatrix); // Update with new model view matrix
		ball.startAnimation(0); // Use animation id 0
		ball.update(); // Update animation
		ball.render(); // Draw
		glBindFramebuffer(GL_FRAMEBUFFER, 0);
		*/

		// Get window size (it will change if the user resizes the window)
		glfwGetWindowSize(window, &width, &height);
		// Set viewport
		glViewport(0, 0, width, height); // The entire window

		/*** Rendering commands here ***/
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		
		// Activate lamp shader
		lampShader.use();
		MATRIX4::translate(modelViewMatrix, lightPosition[0], lightPosition[1], lightPosition[2]); // Place the lamp
		lampShader.setFloatMat4("modelView", modelViewMatrix); // Update with new model view matrix
		lamp.render(); // Draw
		

		// Activate main shader
		mainShader.use();

		// Update camera position in shader
		mainShader.setVec3(
			cameraPositionLocationID, 
			camera.getPosition()[0],
			camera.getPosition()[1],
			camera.getPosition()[2]
		);

		// Use texture
		woodenFloorTexture.use(); 
		// FLOOR
		MATRIX4::translate(modelViewMatrix, 0.0f, -settings.RADIUS*2.0f, 0.0f); // Place the floor
		mainShader.setFloatMat4(modelViewMatrixLocationID, modelViewMatrix); // Update with new model view matrix
		mainShader.setVec3(colorLocationID, 0.75f, 0.75f, 0.75f); // Set color to Light Gray
		floor.render(); // Draw

		
		// Use texture
		footballTexture.use(); 
		// ANIMATION 0
		MATRIX4::translate(modelViewMatrix, -settings.RADIUS*2.0f, 0.0f, -(settings.RADIUS*2.0f)); // Position the sphere
		mainShader.setFloatMat4(modelViewMatrixLocationID, modelViewMatrix); // Update with new model view matrix
		mainShader.setVec3(colorLocationID, 1.0f, 0.25f, 0.25f); // Set color to Red
		ball.startAnimation(0); // Use animation id 0
		ball.update(); // Update animation
		ball.render(); // Draw
		
		
		// ANIMATION 1
		MATRIX4::translate(modelViewMatrix, settings.RADIUS*2.0f, 0.0f, -(settings.RADIUS*2.0f)); // Position the sphere
		mainShader.setFloatMat4(modelViewMatrixLocationID, modelViewMatrix); // Update with new model view matrix
		mainShader.setVec3(colorLocationID, 1.0f, 0.5f, 0.0f); // Set color to Orange
		ball.startAnimation(1); // Use animation id 1
		ball.update(); // Update animation
		ball.render(); // Draw


		// ANIMATION 2
		MATRIX4::translate(modelViewMatrix, -settings.RADIUS*2.0f, 0.0f, 0.0f); // Position the sphere
		mainShader.setFloatMat4(modelViewMatrixLocationID, modelViewMatrix); // Update with new model view matrix
		mainShader.setVec3(colorLocationID, 0.0f, 1.0f, 1.0f); // Set color to Cyan
		ball.startAnimation(2); // Use animation id 2
		ball.update(); // Update animation
		ball.render(); // Draw

		
		// ANIMATION 3
		MATRIX4::translate(modelViewMatrix, 0.0f, 0.0f, 0.0f); // Position the sphere
		mainShader.setFloatMat4(modelViewMatrixLocationID, modelViewMatrix); // Update with new model view matrix
		mainShader.setVec3(colorLocationID, 0.2f, 0.2f, 0.2f); // Set color to Dark Grey
		ball.startAnimation(3); // Use animation id 3
		ball.update(); // Update animation
		ball.render(); // Draw


		// ANIMATION 4
		MATRIX4::translate(modelViewMatrix, settings.RADIUS*2.0f, 0.0f, 0.0f); // Position the sphere
		mainShader.setFloatMat4(modelViewMatrixLocationID, modelViewMatrix); // Update with new model view matrix
		mainShader.setVec3(colorLocationID, 1.0f, 0.0f, 1.0f); // Set color to Magenta
		ball.startAnimation(4); // Use animation id 4
		ball.update(); // Update animation
		ball.render(); // Draw
		
		/*
		footballTexture.use();
		// Place sphere infront of the camera
		MATRIX4::translate(modelViewMatrix, 0.0f, 0.0f, -50.0f);
		// Update with new model view matrix
		mainShader.setFloatMat4(modelViewMatrixLocationID, modelViewMatrix);
		mainShader.setVec3(colorLocationID, 1.0f, 0.5f, 0.0f); // Set color to Orange
		// Render (draw) the simulation one step at the time
		renderSimulationStep(
			stepCounter, // Keeps count of which simulation step to draw
			settings.NUM_STEPS, // Number of steps of the entire simulation
			settings, // Struct that contains specifications about the simulation
			PARTICLE_POSITION_DATA, // A matrix that stores all simulation steps
			bouncyBall // SoftBody simulation model (contains all information about the simulation)
		);
		*/

		// Swap buffers and check for keyboard input or mouse movement events
		glfwSwapBuffers(window);
		glfwPollEvents();

		// Update current time
		GLfloat endTime = (GLfloat)glfwGetTime();

		// Hold each frame for a specified number of seconds (its based on step size of the simulation right now)
		while ((endTime - startTime) < (settings.h)) {

			endTime = (GLfloat)glfwGetTime();
		}
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
		for (int j = 0; j < settings.NUM_POINTS; ++j) {
			// Position Coordinates
			PARTICLE_POSITION_DATA(i, (j * 6) + 1) = softBody.getParticlePositionArray()[j * 3];
			PARTICLE_POSITION_DATA(i, (j * 6) + 2) = softBody.getParticlePositionArray()[j * 3 + 1];
			PARTICLE_POSITION_DATA(i, (j * 6) + 3) = softBody.getParticlePositionArray()[j * 3 + 2];

			// Normals
			PARTICLE_POSITION_DATA(i, (j * 6) + 4) = softBody.getParticleNormalArray()[j * 3];
			PARTICLE_POSITION_DATA(i, (j * 6) + 5) = softBody.getParticleNormalArray()[j * 3 + 1];
			PARTICLE_POSITION_DATA(i, (j * 6) + 6) = softBody.getParticleNormalArray()[j * 3 + 2];
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
		// Position coordinates
		vertexArray[i * 8] = PARTICLE_POSITION_DATA(stepCounter, (i * 6) + 1);
		vertexArray[(i * 8) + 1] = PARTICLE_POSITION_DATA(stepCounter, (i * 6) + 2);
		vertexArray[(i * 8) + 2] = PARTICLE_POSITION_DATA(stepCounter, (i * 6) + 3);
		// Normals
		vertexArray[(i * 8) + 3] = PARTICLE_POSITION_DATA(stepCounter, (i * 6) + 4);
		vertexArray[(i * 8) + 4] = PARTICLE_POSITION_DATA(stepCounter, (i * 6) + 5);
		vertexArray[(i * 8) + 5] = PARTICLE_POSITION_DATA(stepCounter, (i * 6) + 6);
	}

	// Update vertex buffer with new vertex array
	softBody.updatePositions();

	// Draw object
	softBody.render();

	// Incremet the counter to the next simulation step, if the counter has reached maximum number of steps, reset counter
	++stepCounter;
	if (stepCounter > NUM_STEPS) {
		stepCounter = 1;
	}
	vertexArray = nullptr;
}

void saveSimulationSequence(SoftBody &sb, Matrix &DATA, const Settings &simulationSettings, const std::string &simulationName)
{
	// Output File Stream
	std::ofstream outFile;
	if (outFile) {
		// Initiate output stream iterators
		std::ostream_iterator<GLfloat> float_out_it(outFile, "\n"); // Iterator for handling floats
		std::ostream_iterator<GLuint> uint_out_it(outFile, "\n"); // Iterator for handling integers
		std::ostream_iterator<char> out_it_char(outFile, ""); // Iterator for handling characters

		/** Save Simulation Data **/
		std::string filePath = ("Simulations//" + simulationName + ".sim");
		outFile.open(filePath, std::ostream::out);

		// Create a string with time step (h), number of columns and number of rows of the simulation data matrix
		std::string simulationInfo = (
			std::to_string(simulationSettings.h) + "\n" +
			std::to_string(DATA.numRows()) + "\n" +
			std::to_string(DATA.numColumns()) + "\n"
		);
		// Add simulation information at the top of the file
		std::copy(begin(simulationInfo), end(simulationInfo), out_it_char);

		// Create pointer to the simulation data
		float* animationSequence = DATA.getValues();
		// Copy simulation data to the animation sequence file
		std::copy(animationSequence, animationSequence + DATA.size(), float_out_it);
		outFile.close();

		/** Save Mesh Data **/
		filePath = ("Simulations//" + simulationName + ".mesh");
		outFile.open(filePath, std::ofstream::out);

		// Get vertex and triangle info
		int numVertices = sb.getNumMeshVertices();
		int numTriangles = sb.getNumMeshTriangles();

		// Create a string containing vertex and triangle info
		std::string numVerticesAndTriangles = (
			std::to_string(numVertices) + "\n" +
			std::to_string(numTriangles) + "\n"
		);
		// Add vertex and triangle info at the top of the file
		std::copy(begin(numVerticesAndTriangles), end(numVerticesAndTriangles), out_it_char);

		// Create a pointer to the vertex array
		float* vertices = sb.getMeshVertexArray();
		// Copy the vertex array to the mesh data file
		std::copy(vertices, vertices + (numVertices * 8), float_out_it);
		// Create a pointer to the index array
		GLuint* indices = sb.getMeshIndexArray();
		// Copy the index array to the mesh data file
		std::copy(indices, indices + (numTriangles * 3), uint_out_it);
		outFile.close();
	}
}

void framebufferSizeCallback(GLFWwindow* window, int width, int height)
{
	glViewport(0, 0, width, height);
}

void processInput(GLFWwindow* window)
{
	// If user press the escape key, close GLFW
	if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS) {
		glfwSetWindowShouldClose(window, true);
	}
}

void getUniformLocationIDs(Shader &shader, GLint &MV, GLint &CV, GLint &P, GLint &color, GLint &cameraPos)
{
	// Model view matrix
	MV = glGetUniformLocation(shader.ID, "modelView"); 
	if (MV == -1) { // If the variable is not found , -1 is returned
		std::cout << " Unable to locate variable 'modelView' in shader !" << std::endl;
	}

	// Camera view matrix
	CV = glGetUniformLocation(shader.ID, "cameraView");
	if (MV == -1) {
		std::cout << " Unable to locate variable 'cameraView' in shader !" << std::endl;
	}

	// Perspective matrix
	P = glGetUniformLocation(shader.ID, "perspective");
	if (P == -1) { 
		std::cout << " Unable to locate variable 'perspective' in shader !" << std::endl;
	}

	// Object color
	color = glGetUniformLocation(shader.ID, "objectColor");
	if (color == -1) {
		std::cout << " Unable to locate variable 'objectColor' in shader !" << std::endl;
	}

	// Camera position
	cameraPos = glGetUniformLocation(shader.ID, "viewPosition");
	if (cameraPos == -1) {
		std::cout << " Unable to locate variable 'viewPosition' in shader !" << std::endl;
	}
}
