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
#include "Object.h"
#include "Canvas.h"

#include <iostream>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <vector>
#include <string>
#include <sstream>
#include <thread> 

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

// Load animations (used for multithreading)
void loadAnimationsFromFile(Object &bouncyBalls, bool &loadingComplete);


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

	/*************/
	/** SHADERS **/
	/*************/
	// Create our main shader (responsible for all the lighting in the scene)
	Shader mainShader("Shaders//mainShader.vert", "Shaders//mainShader.frag");
	// Activate shader
	mainShader.use();

	// Modelview matrix
	GLfloat modelViewMatrix[16] = { 0 };
	MATRIX4::identity(modelViewMatrix); // Identity matrix as default
	mainShader.setFloatMat4("modelView", modelViewMatrix); // Insert Model View Matrix

	// Perspective projection matrix
	GLfloat perspectiveMatrix[16] = { 0 };
	MATRIX4::perspective(perspectiveMatrix, PI / 3, (GLfloat)width / (GLfloat)height, 0.1f, 1000.0f); // Field of view is set to PI/3 (60 degrees)
	mainShader.setFloatMat4("perspective", perspectiveMatrix); // Insert Projection Matrix

	// Lamp shader
	Shader lampShader("Shaders//lampShader.vert", "Shaders//lampShader.frag");
	lampShader.use();
	lampShader.setFloatMat4("perspective", perspectiveMatrix); // Set perspective projection matrix

	// Position of the lamp
	GLfloat lightPosition[] = {
		0.0f, 80.0f, 0.0f
	};

	// Light color of the lamp
	GLfloat lightColor[] = {
		1.0f, 1.0f, 1.0f
	};

	// Set light position and color in the main shader
	mainShader.use();
	mainShader.setVec3("lightPosition", lightPosition[0], lightPosition[1], lightPosition[2]);
	mainShader.setVec3("lightColor", lightColor[0], lightColor[1], lightColor[2]);


	// Simple loading screen shaders
	Shader loadingShader("Shaders//loadingScreenShader.vert", "Shaders//loadingScreenShader.frag");


	/****************/
	/** SIMULATION **/
	/****************/
	// Settings for the simulation model of the softbody
	Settings settings = {};
	settings.h = 0.01f; // Time step size
	settings.k = 0.1f; // Spring constant
	settings.b = settings.k / 50.0f; // Resistance constant
	settings.g = 9.82f; // Gravitation constant
	settings.DIM = 3; // 3-D (x,y,z)
	settings.RADIUS = 10.0f;
	settings.WEIGHT = 0.75f; // Total weight of the system
	settings.TIME_DURATION = 10.0f; // Specifies how long the simulation should be (given in seconds)
	settings.NUM_STEPS = (int)(settings.TIME_DURATION / settings.h); // Specifies how many steps the simulation will be calculated
	settings.NUM_SEGS = 16; // Horizontal segments for the sphere (number of vertical segments are always twice the number of horizontal segments)
	
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


	/**************/
	/** TEXTURES **/
	/**************/
	Texture woodenFloorTexture("Files//Textures//wooden_floor.jpg");
	Texture stripeTexture("Files//Textures//stripes.png");
	Texture concreteTexture("Files//Textures//concrete_wall.jpg");
	Texture bricksTexture("Files//Textures//bricks.jpg");
	Texture whiteWallTexture("Files//Textures//ceiling.jpg");


	/************/
	/** CAMERA **/
	/************/
	// Create a camera for the scene
	Camera camera(&mainShader);
	// Back the camera away from the origin 
	camera.setPosition(0.0f, 0.0f, 100.0f); // (x, y, z)
	// Add lamp shader to the camera view
	camera.addShader(&lampShader);

	// Create a camera controller
	CameraController cameraController(&camera, window);
	cameraController.setMouseSensitivity(20.0f);


	/*******************/
	/** SCENE OBJECTS **/
	/*******************/
	// Create a sphere mesh for the lamp
	Mesh lampSphere;
	lampSphere.createSphere(8, 10.0f);
	// Create a lamp object
	Object lamp;
	lamp.setMesh(&lampSphere);
	lamp.setColor(lightColor[0], lightColor[1], lightColor[2]);
	lamp.setPosition(lightPosition[0], lightPosition[1], lightPosition[2]);

	// Create a plane mesh for the floor and the ceiling
	GLfloat floorDimX = 500.0f;
	GLfloat floorDimZ = 500.0f;
	GLfloat floorLevel = -settings.RADIUS * 2.0f;
	Mesh plane;
	plane.createPlaneXZ(floorDimX, floorDimZ, 200.0f, floorDimZ);
	// Create floor object
	Object floor;
	floor.setMesh(&plane);
	floor.setTexture(&woodenFloorTexture);
	floor.setPosition(0.0f, floorLevel, 0.0f);

	// Create roof object
	GLfloat ceilingLevel = lightPosition[1] + 10.0f; // Set ceiling level at the lamp position + lamp radius
	Object ceiling;
	ceiling.setMesh(&plane);
	ceiling.setTexture(&whiteWallTexture);
	ceiling.setPosition(0.0f, ceilingLevel, 0.0f);
	ceiling.setRotation(180.0f, 0.0f, 0.0f); // Rotate 180 degrees around x-axis

	// Create cylinder for pillars
	GLfloat pillarHeight = ceilingLevel - floorLevel;
	GLfloat pillarRadius = 10.0f;
	Mesh cylinder;
	cylinder.createCylinder(8, 10.f, pillarHeight); // (Vertical segments, horizontal segments, radius, height)
	// Create cylinder object
	Object pillar;
	pillar.setMesh(&cylinder);
	pillar.setTexture(&concreteTexture);

	// Create and add transformation matrices for different locations to the pillar object
	Transform staticPose;
	staticPose.setPosition(-(floorDimX / 2.0f) + pillarRadius, (pillarHeight / 2.0f) + floorLevel, (floorDimZ / 2.0f) - pillarRadius);
	pillar.addStaticPose(staticPose.matrix4); // Pose 0 (lower left corner)

	staticPose.setPosition((floorDimX / 2.0f) - pillarRadius, (pillarHeight / 2.0f) + floorLevel, (floorDimZ / 2.0f) - pillarRadius);
	pillar.addStaticPose(staticPose.matrix4); // Pose 1 (lower right corner)

	staticPose.setPosition(-(floorDimX / 2.0f) + pillarRadius, (pillarHeight / 2.0f) + floorLevel, -(floorDimZ / 2.0f) + pillarRadius);
	pillar.addStaticPose(staticPose.matrix4); // Pose 2 (upper left corner)

	staticPose.setPosition((floorDimX / 2.0f) - pillarRadius, (pillarHeight / 2.0f) + floorLevel, -(floorDimZ / 2.0f) + pillarRadius);
	pillar.addStaticPose(staticPose.matrix4); // Pose 3 (upper right corner)


	// Create side walls
	Mesh plane2;
	plane2.createPlaneXZ(floorDimZ, pillarHeight, pillarHeight, pillarHeight);
	Object wall;
	wall.setMesh(&plane2);
	wall.setTexture(&bricksTexture);

	// Add several transformation matrices for the same object
	staticPose.setPosition(-(floorDimX / 2.0f), pillarHeight / 2.0f + floorLevel, 0.0f);
	staticPose.setRotation(90.0f, 90.0f, 0.0f); // Rotate 90 degrees around x-axis and 90 degrees around y-axis
	wall.addStaticPose(staticPose.matrix4); // Left side wall

	staticPose.setPosition((floorDimX / 2.0f), pillarHeight / 2.0f + floorLevel, 0.0f);
	staticPose.setRotation(90.0f, -90.0f, 0.0f); // Rotate 90 degrees around x-axis and 90 degrees around y-axis
	wall.addStaticPose(staticPose.matrix4); // Right side wall

	staticPose.setPosition(0.0f, pillarHeight / 2.0f + floorLevel, -(floorDimZ / 2.0f));
	staticPose.setRotation(90.0f, 0.0f, 0.0f); // Rotate 90 degrees around x-axis and 90 degrees around y-axis
	wall.addStaticPose(staticPose.matrix4); // Back wall


	// Create loading screen
	Canvas loadingMessage(&width, &height);
	loadingMessage.createRectangle(500, 40);
	Texture text("Files//Textures//text.png");
	loadingMessage.setTexture(&text);
	loadingMessage.setPosition(0.0f, 100.0f, 0.0f);

	Canvas loadingSymbol(&width, &height);
	loadingSymbol.createRectangle(100, 100);
	Texture gear("Files//Textures//gear.png");
	loadingSymbol.setTexture(&gear);
	loadingSymbol.setPosition(95.0f, -200.0f, 0.0f);
	loadingSymbol.useRotationAnimation(&loadingShader, true);

	Canvas loadingSymbol2(&width, &height);
	loadingSymbol2.createRectangle(80, 80);
	Texture gear2("Files//Textures//gear2.png");
	loadingSymbol2.setTexture(&gear2);
	loadingSymbol2.setPosition(-95.0f, -200.0f, 0.0f);
	loadingSymbol2.useRotationAnimation(&loadingShader, true);

	// Create a mesh and load a prepared mesh for the bouncy ball
	Mesh ball;
	ball.loadMeshData("BouncyBall"); // All animations have the same mesh data
	// Create bouncy ball object
	Object bouncyBalls;
	bouncyBalls.setMesh(&ball);
	bouncyBalls.setTexture(&stripeTexture);
	
	// Load animations from files with a separated thread
	bool loadingComplete = false;
	std::thread loadingThread(loadAnimationsFromFile, std::ref(bouncyBalls), std::ref(loadingComplete));

	GLfloat time = (GLfloat)glfwGetTime();

	// While animations are loading, display a loading screen
	while (!loadingComplete) {

		time = (GLfloat)glfwGetTime();

		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // Clear
		loadingMessage.render(&loadingShader); // Draw

		// Large gear
		loadingShader.setFloat("time", -time);
		loadingSymbol.render(&loadingShader); // Draw

		// Small gear
		loadingShader.setFloat("time", time);
		loadingSymbol2.render(&loadingShader); // Draw

		glfwSwapBuffers(window);// Swap buffers
		glfwPollEvents();
	}

	// Synchronize threads before moving on
	loadingThread.join();
	
	// Time variables
	time = (GLfloat)glfwGetTime();
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
		cameraController.processInput(deltaTime); // Camera navigation

		// Update camera position in shader
		mainShader.use();
		mainShader.setVec3(
			"viewPosition",
			camera.getPosition()[0],
			camera.getPosition()[1],
			camera.getPosition()[2]
		);

		// Get window size (it will change if the user resizes the window)
		glfwGetWindowSize(window, &width, &height);
		// Set viewport
		glViewport(0, 0, width, height); // The entire window

		/*** Rendering commands here ***/
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


		// Draw lamp
		lamp.render(&lampShader);

		// Draw floor
		floor.render(&mainShader);

		// Draw roof
		ceiling.render(&mainShader);

		// Draw pillars
		pillar.renderStaticPoses(&mainShader);

		// Draw walls
		wall.renderStaticPoses(&mainShader);
		
		
		// Draw bouncy balls
		// Animation 0
		bouncyBalls.setColor(1.0f, 0.25f, 0.25f); // Set color to Red
		bouncyBalls.setPosition(-settings.RADIUS*2.0f, 0.0f, -(settings.RADIUS*2.0f));
		bouncyBalls.startAnimation(0);
		bouncyBalls.update(); // Update animation
		bouncyBalls.render(&mainShader); // Draw
		// Animation 1
		bouncyBalls.setColor(1.0f, 0.0f, 1.0f); // Set color to Magenta
		bouncyBalls.setPosition(settings.RADIUS*2.0f, 0.0f, -(settings.RADIUS*2.0f));
		bouncyBalls.startAnimation(1);
		bouncyBalls.update(); // Update animation
		bouncyBalls.render(&mainShader); // Draw
		// Animation 2
		bouncyBalls.setColor(1.0f, 0.5f, 0.0f); // Set color to Orange
		bouncyBalls.setPosition(-settings.RADIUS*2.0f, 0.0f, 0.0f);
		bouncyBalls.startAnimation(2);
		bouncyBalls.update(); // Update animation
		bouncyBalls.render(&mainShader); // Draw
		// Animation 3
		bouncyBalls.setColor(0.0f, 1.0f, 1.0f); // Set color to Cyan
		bouncyBalls.setPosition(0.0f, 0.0f, 0.0f);
		bouncyBalls.startAnimation(3);
		bouncyBalls.update(); // Update animation
		bouncyBalls.render(&mainShader); // Draw
		// Animation 4
		bouncyBalls.setColor(1.0f, 1.0f, 1.0f); // Set color to White
		bouncyBalls.setPosition((floorDimX / 2.0f) - (settings.RADIUS * 4.0f), 0.0f, 0.0f);
		bouncyBalls.startAnimation(4);
		bouncyBalls.update(); // Update animation
		bouncyBalls.render(&mainShader); // Draw
		

		/*
		// DRAW SIMULATION 
		footballTexture.use();
		// Place sphere infront of the camera
		MATRIX4::translate(modelViewMatrix, 0.0f, 0.0f, -50.0f);
		// Update with new model view matrix
		mainShader.setFloatMat4("modelView", modelViewMatrix);
		mainShader.setVec3("objectColor", 1.0f, 0.5f, 0.0f); // Set color to Orange
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

void loadAnimationsFromFile(Object &bouncyBalls, bool &loadingComplete)
{
	bouncyBalls.addAnimation("BouncyBall_back_left"); // ID: 0
	bouncyBalls.addAnimation("BouncyBall_back_right"); // ID: 1
	bouncyBalls.addAnimation("BouncyBall_front_left"); // ID: 2
	bouncyBalls.addAnimation("BouncyBall_straight_up"); // ID: 3
	bouncyBalls.addAnimation("BouncyBall_front_right"); // ID: 4
	loadingComplete = true;
}
