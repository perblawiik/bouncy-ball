#include <glad/glad.h>
#include <glad/glad.c>
#include <GLFW/glfw3.h>

#include "Shader.h"
#include "Matrix.h"

#include <iostream>

/** CONSTANTS **/
const unsigned int SCR_WIDTH = 800;
const unsigned int SCR_HEIGHT = 600;
const GLfloat PI = 3.14159265359f;

/** STRUCTS **/
// Struct with all constant values for the simulation
struct Settings 
{
	GLfloat h; // Step
	GLfloat k; // Spring constant
	GLfloat b; // Resistance constant
	GLfloat g; // Gravitation constant
	GLfloat SPHERE_RADIUS; // Radius of the bouncy ball

	GLint NUM_BONDS;
	GLint NUM_POINTS;
	GLint DIM;
};

/** FUNCTION DECLARATIONS **/
// A callback function on the window that gets called each time the window is resized
void framebuffer_size_callback(GLFWwindow* window, int width, int height);

// Checks if a key is currently being pressed
void processInput(GLFWwindow *window);

// Calculates differential equations of the system
void calculateSimulation(const Settings &s, Matrix &m, Matrix &X, Matrix &I, Matrix &V, Matrix &Vp,
	Matrix &Fk, Matrix &Fkp, Matrix &zeros, Matrix &vec1, Matrix &vec2, Matrix &diff);

// Perspective matrix
void mat4Perspective(GLfloat M[], const GLfloat &vertFov, const GLfloat &aspect, const GLfloat &zNear, const GLfloat &zFar);
void mat4Translate(GLfloat M[], const GLfloat &x, const GLfloat &y, const GLfloat &z);

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
	GLfloat MV[16] = {

		1.0f, 0.0f, 0.0f, 0.0f,
		0.0f, 1.0f, 0.0f, 0.0f,
		0.0f, 0.0f, 1.0f, 0.0f,
		0.0f, 0.0f, 0.0f, 1.0f
	};
	mat4Translate(MV, 0.0f, 0.0f, -30.0f);
	// Perspective matrix
	GLfloat P[16] = {

		1.0f, 0.0f, 0.0f, 0.0f,
		0.0f, 1.0f, 0.0f, 0.0f,
		0.0f, 0.0f, 1.0f, 0.0f,
		0.0f, 0.0f, 0.0f, 1.0f
	};
	mat4Perspective(P, PI / 3, 1.0f, 0.1f, 100.0f);

	Shader myShader("Shaders//vertex.glsl", "Shaders//fragment.glsl");

	// Particle positions per cycle (passed into the shader)
	GLfloat positions[342] = { 0.0f };

	// Columns per row in the vertex array (coordinates + normals)
	int stride = 6;
	
	// Number of horizontal segments of the sphere (2 is minimum)
	int numHorizontalSegments = 8;
	int numVerticalSegments = 2 * numHorizontalSegments;
	int numVertices = 1 + (numHorizontalSegments - 1) * numVerticalSegments + 1; // top + middle + bottom
	std::cout << "numVerts: " << numVertices << std::endl;
	int numTriangles = numVerticalSegments + (numHorizontalSegments - 2) * 4 * numHorizontalSegments + numVerticalSegments; // top + middle + bottom
	std::cout << "numTris: " << numTriangles << std::endl;
	std::cout << "Vertex array size: " << 3 * numVertices << std::endl;
	
	GLfloat *vertices = new GLfloat[numVertices * stride]; // Initialize vertex array
	GLuint *indices = new GLuint[numTriangles * 3]; // Initialize index array
	GLfloat radius = 4.0f;

	/** Generate vertex array **/
	// Bottom vertex
	vertices[0] = 0.0f; vertices[1] = -radius; vertices[2] = 0.0f; // Coordinates
	vertices[3] = 0.0f; vertices[4] = -radius; vertices[5] = 0.0f; // Normal

	GLfloat sampleRate = PI / numHorizontalSegments; // Number of steps 
	GLfloat theta = -PI + sampleRate; // Go from bottom to top (Y € -PI < theta < PI )
	GLfloat phi = 0; // Begin at Z = 0 (Z € 0 < phi < 2PI )

	// Generate middle part vertices with normals
	int index = 5; // Skip first 6 (the bottom vertex with normal already specified)
	for (int i = 0; i < numHorizontalSegments - 1; ++i) {

		float Y = radius*cos(theta); // Y-coordinate
		float R = radius*sin(theta); // radius

		for (int j = 0; j < numVerticalSegments; ++j) {
			// Vertex (x, y, z)
			vertices[++index] = R * sin(phi);
			vertices[++index] = Y;
			vertices[++index] = R * cos(phi);
			// Normal (x, y, z)
			vertices[++index] = R * sin(phi);
			vertices[++index] = Y;
			vertices[++index] = R * cos(phi);

			phi += sampleRate;
		}
		theta += sampleRate;
	}

	// Top vertex
	vertices[++index] = 0.0f; vertices[++index] = radius; vertices[++index] = 0.0f;
	vertices[++index] = 0.0f; vertices[++index] = radius; vertices[++index] = 0.0f;
	
	/** Generate index array */
	// Bottom cap
	index = -1;
	for (int i = 0; i < numVerticalSegments; ++i) {

		indices[++index] = 0;

		if ((i + 2) <= numVerticalSegments) {
			indices[++index] = i + 2;
		}
		else {
			indices[++index] = (i + 2) - numVerticalSegments;
		}
		indices[++index] = i + 1;
	}
	
	
	// Middle part
	int v0 = 1;
	for (int i = 0; i < numHorizontalSegments - 2; i++) {
		for (int j = 0; j < numVerticalSegments-1; ++j) {
			// One rectangle at a time (two triangles)
			indices[++index] = v0;
			indices[++index] = v0 + 1;
			indices[++index] = numVerticalSegments + v0;
			indices[++index] = v0 + 1;
			indices[++index] = numVerticalSegments + v0 + 1;
			indices[++index] = numVerticalSegments + v0;
			++v0;
		}
		indices[++index] = v0;
		indices[++index] = (v0 + 1) - numVerticalSegments;
		indices[++index] = numVerticalSegments + v0;
		indices[++index] = (v0 + 1) - numVerticalSegments;
		indices[++index] = v0 + 1;
		indices[++index] = numVerticalSegments + v0;
		++v0;
	}
	
	
	// Top cap
	int lastVertexIndex = numVertices - 1;
	for (int i = 0; i < numVerticalSegments; ++i) {

		indices[++index] = lastVertexIndex;

		if ((lastVertexIndex - 2 - i) >= lastVertexIndex - numVerticalSegments) {
			indices[++index] = lastVertexIndex - 2 - i;
		}
		else {
			indices[++index] = lastVertexIndex - numVerticalSegments - 1;
		}

		indices[++index] = lastVertexIndex - 1 - i;
	}

	// Vertex Buffer Object, Vertex Array Object, Element Buffer Object
	GLuint VBO, VAO, EBO;
	glGenBuffers(1, &VBO);
	glGenBuffers(1, &EBO);
	glGenVertexArrays(1, &VAO);

	// 1. Bind Vertex Array Object
	glBindVertexArray(VAO);
	// 2. Copy our vertices array in a buffer for OpenGL to use
	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferData(GL_ARRAY_BUFFER, stride * numVertices * sizeof(GLfloat), vertices, GL_STATIC_DRAW);
	// 3. Copy our index array in a element buffer for OpenGL to use
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, 3 * numTriangles * sizeof(GLuint), indices, GL_STATIC_DRAW);

	// Tell OpenGL how it should interpret the vertex data (per vertex attribute)
	// glVertexAttribPointer Parameters :
	// 1. Specifies which vertex attribute we want to configure. 
	//	  This sets the location of the vertex attribute to 0 and since we want to pass data to this vertex attribute, we pass in 0.
	// 2. Specifies the size of the vertex attribute (vec3 is composed of 3 values).
	// 3. Specifies the type of the data (float in this case)
	// 4. Specifies if we want the data to be normalized.
	// 5. Known as "the stride" and tells us the space between consecutive vertex attributes. 
	//    Since the next set of position data is located exactly 3 times the size of a float away we specify that value as the stride.
	// 6. This is the offset of where the position data begins in the buffer. 
	//    Since the position data is at the start of the data array this value is just 0.
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, stride * sizeof(GLfloat), (void*)0); // Vertex coordinates
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, stride * sizeof(GLfloat), (void*)(3 * sizeof(GLfloat))); // normals
	glEnableVertexAttribArray(0); // Vertex coordinates
	glEnableVertexAttribArray(1); // Normals

	// Deactivate (unbind) the VAO and the buffers again.
	// Do NOT unbind the index buffer while the VAO is still bound.
	// The index buffer is an essential part of the VAO state.
	glBindVertexArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

	/**********S I M U L A T I O N**********/

	Settings settings = {};
	settings.h = 0.001f; // Step
	settings.k = 6000.0f; // Spring constant
	settings.b = 10.0f; // Resistance constant
	settings.g = 9.82f; // Gravitation constant
	settings.NUM_BONDS = (2 * numHorizontalSegments * numVerticalSegments) - numVerticalSegments; // Total number of bonds
	settings.NUM_POINTS = numVertices; // Total number of points (particles)
	settings.DIM = 3; // 2-D
	settings.SPHERE_RADIUS = radius;

	// Total weight of the system
	const float WEIGHT = 4.0f;

	//Masses per particle
	Matrix m(settings.NUM_POINTS, 1); // Create Nx1 matrix
	for (int i = 0; i < m.size(); ++i) {
		m[i] = WEIGHT / (float)settings.NUM_POINTS; // All masses divided equally
	}

	//Particle x, y Pos [Xx Xy] / per particle
	Matrix X(settings.NUM_POINTS, settings.DIM);
	for (int i = 0; i < settings.NUM_POINTS; ++i) {
		X[i * settings.DIM] = vertices[i * stride];
		X[i * settings.DIM + 1] = vertices[i * stride + 1];
		X[i * settings.DIM + 2] = vertices[i * stride + 2];
	}

	std::cout << "Number of springs: " << settings.NUM_BONDS << std::endl;
	Matrix I(settings.NUM_BONDS, 2);

	index = 1;
	// Generate horizontal spring bond relations between the particles 
	for (int i = 0; i < numHorizontalSegments - 1; ++i) {

		int n = (i * numVerticalSegments) + 2;
		for (int j = 0; j < numVerticalSegments - 1; ++j) {
			I(index, 1) = (float)(n + j);
			I(index, 2) = (float)(n + j + 1);
			++index;
		}
		I(index, 1) = I((index - 1), 2);
		I(index, 2) = (float)n;
		++index;
	}

	// Generate vertical spring bond relations between the particles 
	for (int i = 1; i <= numVerticalSegments; ++i) {

		int n = 1;
		int m = (n + i);

		// Bottom part
		I(index, 1) = (float)n;
		I(index, 2) = (float)m;
		++index;

		// Middle part
		for (int j = 0; j < numHorizontalSegments -2; ++j) {

			I(index, 1) = (float)m;
			I(index, 2) = (float)m + numVerticalSegments;
			++index;
			m = m + numVerticalSegments;
		}

		// Top part
		I(index, 1) = (float)m;
		I(index, 2) = (float)numVertices;
		++index;
	}

	// Starting velocity [Vx Vy]/ per particle
	Matrix V(settings.NUM_POINTS, settings.DIM); // All set to zero by default

	// Acceleration dV/dt 
	Matrix Vp(settings.NUM_POINTS, settings.DIM); // All set to zero by default

	// Fk spring starting force [F]/ per spring  (applied directionally later, depending on spring orientation) 
	Matrix Fk(settings.NUM_BONDS, settings.DIM); // All set to zero by default

	// Fkp = dFp/dt
	Matrix Fkp(settings.NUM_BONDS, settings.DIM); // All set to zero by default

	// Used to set Vp to zero each cycle
	Matrix zeros(settings.NUM_POINTS, settings.DIM);

	// Dummy vectors (1xDIM matrix)
	Matrix vec1(1, settings.DIM);
	Matrix vec2(1, settings.DIM);

	// Vector from point 1 to point 2;
	Matrix diff(1, settings.DIM);

	/***************************************/

	// Wireframe mode
	//glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	// Hide the back face of the triangles
	glEnable(GL_CULL_FACE);

	/***************************************/

	/** RENDER LOOP **/
	while (!glfwWindowShouldClose(window))
	{
		// Read input
		processInput(window);

		// Rendering commands here
		glClearColor(0.25f, 0.25f, 0.25f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT);

		/**********S I M U L A T I O N**********/
		//-------------------------------------//

		// Calculate one cycle for the system
		calculateSimulation(settings, m, X, I, V, Vp, Fk, Fkp, zeros, vec1, vec2, diff);

		// Update particle positions
		X.copyValues(positions);

		//-------------------------------------//
		/***************************************/

		// Activate shader
		myShader.use();

		// Model View Matrix
		myShader.setFloatMat4("MV", MV);
		 // Projection Matrix
		myShader.setFloatMat4("P", P);
		// Insert particle positions in shader
		myShader.setFloat("positions", positions, numVertices * 3);

		// Draw object
		glBindVertexArray(VAO);
		// (mode, vertex count, type, element array buffer offset)
		glDrawElements(GL_TRIANGLES, numTriangles * 3, GL_UNSIGNED_INT, 0);
		glBindVertexArray(0);
		
		// Swap buffers and check for keyboard input or mouse movement events
		glfwSwapBuffers(window);
		glfwPollEvents();
	}

	// OPTIONAL: de-allocate all resources once they've outlived their purpose:
	glDeleteVertexArrays(1, &VAO);
	glDeleteBuffers(1, &VBO);
	glDeleteBuffers(1, &EBO);

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

void calculateSimulation(const Settings &s, Matrix &m, Matrix &X, Matrix &I, Matrix &V, Matrix &Vp,
	Matrix &Fk, Matrix &Fkp, Matrix &zeros, Matrix &vec1, Matrix &vec2, Matrix &diff)
{
	// Set to zero so the components from each connected spring can be += and added separately
	Vp = zeros; 
	for (int n = 1; n <= s.NUM_BONDS; ++n) // Loop through the springs
	{
		int index1 = (int)I(n, 1); // Index from point 1
		int index2 = (int)I(n, 2); // Index from point 2

		// Get point coordinates based on bond indices I
		X.copyRow(index1, vec1); // Copy position of point 1 to vec1
		X.copyRow(index2, vec2); // Copy position of point 2 to vec2

		// Vector from point 1 to 2
		diff = vec1 - vec2; 

		// Normalise it, used to give the Fk and Fb direction
		diff.normalize();

		// Get deltaV, speed difference between the particles in the spring's direction
		V.copyRow(index1, vec1); // Copy velocity of point 1 to vec1
		V.copyRow(index2, vec2); // Copy velocity of point 2 to vec2
		float dV = diff.dot(vec1 - vec2);

		// Apply spring influence to the connected particles, 
		// first one (Vp1) is added, second is subtracted, in the springs direction since they will be either
		// both pulled towards eachother or drawn away from eachother.
		Vp.copyRow(index1, vec1); // Copy acceleration of point 1 to vec1
		Vp.copyRow(index2, vec2); // Copy acceleration of point 2 to vec2
		float m1 = m(index1, 1); // Mass of point 1
		float m2 = m(index2, 1); // Mass of point 2
		float F = Fk(n, 1); // Spring force
		Vp.replaceRow(index1, vec1 - (diff * (1.0f / m1) * (s.b*dV + F))); 
		Vp.replaceRow(index2, vec2 + (diff * (1.0f / m2) * (s.b*dV + F)));

		// The derivative for Fk
		Fkp(n, 1) = s.k * dV;
	}

	// Add gravity for all points (-g to y coordinate)
	
	for (int r = 1; r <= s.NUM_POINTS; ++r) { 
		Vp(r, 2) = Vp(r, 2) - s.g;
	}
	
	// Approximating the new values using: X_n + 1 = X_n + h * X'_n
	// (they're not supressed for debugging purposes)
	V = V + (Vp * s.h);
	Fk = Fk + (Fkp * s.h);
	X = X + (V * s.h);

	 // Code that flips Y - ward velocity when the particle has Xy < -4.0 
	for (int j = 1; j <= s.NUM_POINTS; ++j) {
		 
		if (X(j, 2) < (-s.SPHERE_RADIUS * 2)) {
			V(j, 2) = 0.0f;
			X(j, 2) = -s.SPHERE_RADIUS * 2;
		}
	}
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