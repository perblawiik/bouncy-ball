#include <glad/glad.h>
#include <glad/glad.c>
#include <GLFW/glfw3.h>

#include "Shader.h"
#include <iostream>

/** CONSTANTS **/
const unsigned int SCR_WIDTH = 800;
const unsigned int SCR_HEIGHT = 600;

/** FUNCTION DECLARATIONS **/
// A callback function on the window that gets called each time the window is resized
void framebuffer_size_callback(GLFWwindow* window, int width, int height);

// Checks if a key is currently being pressed
void processInput(GLFWwindow *window);

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
	GLfloat time = glfwGetTime();

	/** RENDER LOOP **/
	while (!glfwWindowShouldClose(window))
	{
		// Read input
		processInput(window);

		// Rendering commands here
		glClearColor(0.25f, 0.25f, 0.25f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT);

		// Activate shader
		//glUseProgram(shaderProgram);
		myShader.use();

		// Update time and pass in to the shader
		time = glfwGetTime();
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