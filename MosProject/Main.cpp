#include <glad/glad.h>
#include <glad/glad.c>
#include <GLFW/glfw3.h>

#include "Shader.h"
#include "Matrix.h"

#include <iostream>

/** CONSTANTS **/
const unsigned int SCR_WIDTH = 800;
const unsigned int SCR_HEIGHT = 600;

/** FUNCTION DECLARATIONS **/
// A callback function on the window that gets called each time the window is resized
void framebuffer_size_callback(GLFWwindow* window, int width, int height);

// Checks if a key is currently being pressed
void processInput(GLFWwindow *window);

void simulationCalculations();

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

	Shader myShader("Shaders//vertex.glsl", "Shaders//fragment.glsl");

	// Rectangle
	float vertices[] = {
		0.5f,  0.5f, 0.0f,  // top right
		0.5f, -0.5f, 0.0f,  // bottom right
		-0.5f, -0.5f, 0.0f,  // bottom left
		-0.5f,  0.5f, 0.0f   // top left 
	};
	unsigned int indices[] = {  
		0, 1, 3,   // first triangle
		1, 2, 3    // second triangle
	};

	// Vertex Buffer Object, Vertex Array Object, Element Buffer Object
	unsigned int VBO, VAO, EBO;
	glGenBuffers(1, &VBO);
	glGenBuffers(1, &EBO);
	glGenVertexArrays(1, &VAO);

	// 1. Bind Vertex Array Object
	glBindVertexArray(VAO);
	// 2. Copy our vertices array in a buffer for OpenGL to use
	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
	// 3. Copy our index array in a element buffer for OpenGL to use
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);

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
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(0);

	// Wireframe mode
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

	// Time variable
	GLfloat time = (GLfloat)glfwGetTime();

	simulationCalculations();

	/** RENDER LOOP **/
	while (!glfwWindowShouldClose(window))
	{
		// Read input
		processInput(window);

		// Rendering commands here
		glClearColor(0.25f, 0.25f, 0.25f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT);

		// Activate shader
		myShader.use();

		// Update time and pass in to the shader
		time = (GLfloat)glfwGetTime();
		myShader.setFloat("time", time);

		// Draw object
		glBindVertexArray(VAO);
		glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
		
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

void simulationCalculations() 
{
	// Step
	const GLfloat h = 0.01f; 
	// Spring constant
	const GLfloat k = 100.0f;
	// Resistance constant
	const GLfloat b = 5.0f; 
	// Gravitation constant
	const GLfloat g = 9.82f;

	const GLint NUM_BONDS = 5;
	const GLint NUM_POINTS = 4;
	const GLint DIM = 2; // 2-D

	//Masses per particle ( m = [1; 1; 1; 1]; )
	Matrix m(NUM_POINTS, 1); // Create 4x1 matrix
	m(1, 1) = 1.0f; m(2, 1) = 1.0f; m(3, 1) = 1.0f; m(4, 1) = 1.0f;
	//Particle x, y Pos [Xx Xy] / per particle ( X = [10 10; 20 10; 15 15; 25 15]; )
	Matrix X(NUM_POINTS, DIM);
	X(1, 1) = 10.0f; X(1, 2) = 10.0f;
	X(2, 1) = 20.0f; X(2, 2) = 10.0f;
	X(3, 1) = 15.0f; X(3, 2) = 15.0f;
	X(4, 1) = 25.0f; X(4, 2) = 15.0f;
	// Particle indices for spring bonds [i1 i2]/ per spring ( I = [1 2; 2 3; 1 3; 3 4; 2 4]; )
	Matrix I(NUM_BONDS, 2);
	I(1, 1) = 1.0f; I(1, 2) = 2.0f;
	I(2, 1) = 2.0f; I(2, 2) = 3.0f;
	I(3, 1) = 1.0f; I(3, 2) = 3.0f;
	I(4, 1) = 3.0f; I(4, 2) = 4.0f;
	I(5, 1) = 2.0f; I(5, 2) = 4.0f;
	// Starting velocity [Vx Vy]/ per particle ( V = [0 0; 0 0; 0 0; 0 0];  )
	Matrix V(NUM_POINTS, DIM);
	// Acceleration dV/dt ( Vp = [0 0; 0 0; 0 0; 0 0]; )
	Matrix Vp(NUM_POINTS, DIM);
	// Fk spring starting force [F]/ per spring  (applied directionally later, depending on spring orientation) ( Fk = [0; 0; 0; 0; 0]; )
	Matrix Fk(NUM_BONDS, 1);
	// Fkp = dFp/dt ( Fkp [0; 0; 0; 0; 0]; )
	Matrix Fkp(NUM_BONDS, 1);

	// Used to set Vp to zero each cycle
	Matrix zeros(NUM_POINTS, DIM);
	// Vector from point 1 to point 2;
	Matrix vec1(1, DIM);
	Matrix vec2(1, DIM);
	Matrix diff(1, DIM);

	int num_steps = 100;
	for (int cycle = 0; cycle < num_steps; ++cycle) {

		// Set to zero so the components from each connected spring can be += and added separately
		Vp = zeros; 
		for (int n = 1; n <= NUM_BONDS; ++n) // Loop through the springs
		{
			// Get point coordinates based on bond indices I
			X.copyRow((int)I(n, 1), vec1);
			X.copyRow((int)I(n, 2), vec2);

			// Get vector from point 1 to 2
			diff = vec1 - vec2; /* MATLAB: dif = X(I(n, 1), :) - X(I(n, 2), :); */

			// Normalise it, used to give the Fk and Fb direction
			diff.normalize(); /* MATLAB: nDif = dif / norm(dif); */

			// Get deltaV, speed difference between the particles in the spring's direction
			V.copyRow((int)I(n, 1), vec1);
			V.copyRow((int)I(n, 2), vec2);
			float dV = diff.dot(vec1 - vec2); /* MATLAB: dV = dot(V(I(n, 1), :) - V(I(n, 2), :), nDif); */

			// Apply spring influence to the connected particles, 
			// first one (Vp1) is added, second is subtracted, in the springs direction since they will be either
			// both pulled towards eachother or drawn away from eachother.
			Vp.copyRow((int)I(n, 1), vec1);
			Vp.copyRow((int)I(n, 2), vec2);
			float m1 = m((int)I(n, 1), 1); // Mass of point 1
			float m2 = m((int)I(n, 2), 1); // Mass of point 2
			float F = Fk(n, 1); // Spring force
			diff = diff * ((1 / m1) * (b*dV + F));
			vec1 = vec1 - diff;
			vec2 = vec2 + diff;
			Vp.replaceRow((int)I(n, 1), vec1); /* MATLAB: Vp(I(n, 1), :) = Vp(I(n, 1), :) - 1 / m(I(n, 1)) * (b*dV + Fk(n))*nDif; */ 
			Vp.replaceRow((int)I(n, 2), vec2); /* MATLAB: Vp(I(n, 2), :) = Vp(I(n, 2), :) + 1 / m(I(n, 2)) * (b*dV + Fk(n))*nDif; */ 

			// The derivative for Fk
			Fkp(n, 1) = k * dV; /* MATLAB: Fkp(n) = k * dV; */
		}

		// Add gravity for all points (-g to y coordinate)
		for (int r = 1; r <= Vp.numRows(); ++r) { /* MATLAB: Vp = Vp - [0 g]; */
			Vp(r, 2) = Vp(r, 2) - g;
		}

		// Approximating the new values using: X_n + 1 = X_n + h * X'_n
		// (they're not supressed for debugging purposes)
		V = V + (Vp * h); /* MATLAB: V = V + h * Vp */
		Fk = Fk + (Fkp * h); /* MATLAB: Fk = Fk + h * Fkp */
		X = X + (V * h); /* MATLAB: X = X + h * V */


		// Print out velocities, spring forces and positions for the particles
		std::cout << "Cycle: " << cycle << std::endl;

		std::cout << "V = " << std::endl;
		for (int i = 1; i <= V.numRows(); ++i) {
			std::cout << V(i, 1) << " " << V(i, 2) << std::endl;
		}
		std::cout << "Fk = " << std::endl;
		for (int i = 1; i <= Fk.numRows(); ++i) {
			std::cout << Fk(i, 1) << std::endl;
		}
		std::cout << "X = " << std::endl;
		for (int i = 1; i <= V.numRows(); ++i) {
			std::cout << X(i, 1) << " " << X(i, 2) << std::endl;
		}

		//% Code that flips Y - ward velocity when the particle has Xy<0
		//V(:, 2) = (X(:, 2)>0).*V(:, 2) - (X(:, 2)<0).*V(:, 2);
	}
}