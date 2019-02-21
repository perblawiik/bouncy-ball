#pragma once
#ifndef MESH_H
#define MESH_H

#include <glad/glad.h>
#include <math.h>
#include <iostream>

class Mesh 
{
public:

	Mesh()
		: meshIsEmpty(true), VAO(0), VBO(0), EBO(0), vertices(nullptr), indices(nullptr), numVertices(0), numTriangles(0), stride(8)
	{ }

	// Copy constructor
	Mesh(const Mesh &m)
		: meshIsEmpty(m.meshIsEmpty), VAO(m.VAO), VBO(m.VBO), EBO(m.EBO), vertices(new GLfloat[m.numVertices * m.stride]), indices(new GLuint[m.numTriangles * 3]), numVertices(m.numVertices), numTriangles(m.numTriangles), stride(m.stride)
	{
		// Copy vertex array
		for (int i = 0; i < (numVertices * stride); ++i) {
			vertices[i] = m.vertices[i];
		}

		// Copy index array
		for (int i = 0; i < numTriangles * 3; ++i) {
			indices[i] = m.indices[i];
		}
	}

	// Assignment operator
	Mesh &operator=(const Mesh &m)
	{
		// Create a local copy
		Mesh copy(m);

		// Swap member variables
		std::swap(this->meshIsEmpty, copy.meshIsEmpty);
		std::swap(this->VAO, copy.VAO);
		std::swap(this->VBO, copy.VBO);
		std::swap(this->EBO, copy.EBO);
		std::swap(this->vertices, copy.vertices);
		std::swap(this->indices, copy.indices);
		std::swap(this->numVertices, copy.numVertices);
		std::swap(this->numTriangles, copy.numTriangles);
		std::swap(this->stride, copy.stride);

		// Return this instance to allow cascading
		return *this;
	}

	~Mesh()
	{
		this->clean();
	}

	void createPlaneXZ(const GLfloat &WIDTH, const GLfloat &HEIGHT)
	{
		meshIsEmpty = false;

		const GLfloat vertexData[] = {
			// Position Coordinates                    // Normals          // Texture coordinates
			-(WIDTH / 2.0f), 0.0f, -(HEIGHT / 2.0f),   0.0f, 1.0f, 0.0f,   0.0f, 1.0f, // Upper left
			 (WIDTH / 2.0f), 0.0f, -(HEIGHT / 2.0f),   0.0f, 1.0f, 0.0f,   1.0f, 1.0f, // Upper right
			 (WIDTH / 2.0f), 0.0f,  (HEIGHT / 2.0f),   0.0f, 1.0f, 0.0f,   0.0f, 0.0f, // Lower right
			-(WIDTH / 2.0f), 0.0f,  (HEIGHT / 2.0f),   0.0f, 1.0f, 0.0f,   1.0f, 0.0f  // Lower left
		};

		const GLuint indexData[] = {
			0, 2, 1, // First triangle
			0, 3, 2 // Second triangle
		};

		this->numVertices = 4;
		this->numTriangles = 2;
		this->vertices = new GLfloat[numVertices * stride];
		this->indices = new GLuint[numTriangles * 3];

		for (int i = 0; i < numVertices * stride; ++i) {
			vertices[i] = vertexData[i];
		}

		for (int i = 0; i < numTriangles * 3; ++i) {
			indices[i] = indexData[i];
		}

		// Generate ID's for buffers and vertex array to use in the shader
		this->generateBuffersAndVertexArrayObject();

		// Store our vertex array and index array in buffers for OpenGL to use
		this->bindBuffersAndVertexArrayObject();

		// Tell OpenGL how it should interpret the vertex data (per vertex attribute)
		this->configureVertexArrayAttributes();

		// Deactivate (unbind) the VAO and the buffers again.
		this->unbindBuffersAndVertexArrayObject();
	}

	// Generates a sphere mesh with a specified number of horizontal segments and radius as input parameters.
	void createSphere(const int segments, const float &radius)
	{
		meshIsEmpty = false;
		int numHorizontalSegments = segments;

		// Minium amount of horizontal segments is 2
		if (numHorizontalSegments < 2) {
			numHorizontalSegments = 2;
		}

		// Number of vertical segments of the sphere
		int numVerticalSegments = 2 * numHorizontalSegments;
		this->numVertices = 1 + (numHorizontalSegments - 1) * numVerticalSegments + 1; // top + middle + bottom
		this->numTriangles = numVerticalSegments + (numHorizontalSegments - 2) * 4 * numHorizontalSegments + numVerticalSegments; // top + middle + bottom

		// Information about the sphere mesh
		std::cout << "Number of vertices: " << numVertices << std::endl;
		std::cout << "Number of triangles: " << numTriangles << std::endl;
		std::cout << "Sphere radius: " << radius << std::endl;

		this->vertices = new GLfloat[numVertices * stride]; // Initialize vertex array
		this->indices = new GLuint[numTriangles * 3]; // Initialize index array

		/** Generate vertex array **/
		// Bottom vertex
		vertices[0] = 0.0f; vertices[1] = -radius; vertices[2] = 0.0f; // Coordinates
		vertices[3] = 0.0f; vertices[4] = -1.0f; vertices[5] = 0.0f; // Normal
		vertices[6] = 0.5f; vertices[7] = 0.0f;

		const GLfloat PI = GLOBAL_CONSTANTS::PI;
		GLfloat sampleRate = PI / numHorizontalSegments; // Number of steps 
		GLfloat theta = -PI + sampleRate; // Go from bottom to top (Y € -PI < theta < PI )
		GLfloat phi = 0; // Begin at Z = 0 (Z € 0 < phi < 2PI )

		// Generate middle part vertices with normals
		int index = 7; // Skip first 7 (the bottom vertex with normal and texture coordinates already specified)
		for (int i = 0; i < numHorizontalSegments - 1; ++i) {

			float Y = radius * cos(theta); // Y-coordinate
			float R = radius * sin(theta); // radius

			for (int j = 0; j < numVerticalSegments; ++j) {
				// Vertex (x, y, z)
				vertices[++index] = R * sin(phi);
				vertices[++index] = Y;
				vertices[++index] = R * cos(phi);
				// Normal (x, y, z)
				vertices[++index] = R * sin(phi);
				vertices[++index] = Y;
				vertices[++index] = R * cos(phi);
				// Texture Coordinates (s, t)
				vertices[++index] = (float)j/numVerticalSegments;
				vertices[++index] = 1.0f - (float)(i + 1) / numHorizontalSegments;

				phi += sampleRate;
			}
			theta += sampleRate;
		}

		// Top vertex
		vertices[++index] = 0.0f; vertices[++index] = radius; vertices[++index] = 0.0f; // Coordinates
		vertices[++index] = 0.0f; vertices[++index] = radius; vertices[++index] = 0.0f; // Normal
		vertices[++index] = 0.5f; vertices[++index] = 1.0f;

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
			for (int j = 0; j < numVerticalSegments - 1; ++j) {
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

		// Generate ID's for buffers and vertex array to use in the shader
		this->generateBuffersAndVertexArrayObject();

		// Store our vertex array and index array in buffers for OpenGL to use
		this->bindBuffersAndVertexArrayObject();

		// Tell OpenGL how it should interpret the vertex data (per vertex attribute)
		this->configureVertexArrayAttributes();
		
		// Deactivate (unbind) the VAO and the buffers again.
		this->unbindBuffersAndVertexArrayObject();
	}

	void updateVertexBufferData()
	{
		glBindBuffer(GL_ARRAY_BUFFER, VBO);
		glBufferData(GL_ARRAY_BUFFER, stride * numVertices * sizeof(GLfloat), this->vertices, GL_DYNAMIC_DRAW);
		glBindBuffer(GL_ARRAY_BUFFER, 0);
	}

	int getNumVertices()
	{
		return this->numVertices;
	}

	int getStride()
	{
		return this->stride;
	}

	GLfloat* getVertexArray()
	{
		if (meshIsEmpty) {
			std::cout << "WARNING! Mesh::Get Vertex Array::No Vertex Array Exists!" << std::endl;
		}
		return this->vertices;
	}

	GLuint* getIndexArray()
	{
		if (meshIsEmpty) {
			std::cout << "WARNING! Mesh::Get Index Array::No Index Array Exists!" << std::endl;
		}
		return this->indices;
	}

	void render()
	{
		if (!meshIsEmpty) {
			// Bind the VAO
			glBindVertexArray(this->VAO);
			// Draw the mesh (mode, vertex count, type, element array buffer offset)
			glDrawElements(GL_TRIANGLES, this->numTriangles * 3, GL_UNSIGNED_INT, 0);
			// Unbind the VAO
			glBindVertexArray(0);
		}
	}

private:

	bool meshIsEmpty;

	GLuint VAO; // Vertex Array Object ID
	GLuint VBO; // Vertex Buffer Object ID
	GLuint EBO; // Element Buffer Object ID
	GLfloat* vertices; // Vertex array
	GLuint* indices; // Index array

	int stride; // Number of elements per row in the vertex array
	int numTriangles; // Number of triangels of the mesh
	int numVertices; // Number of vertices of the mesh

	void generateBuffersAndVertexArrayObject()
	{
		glGenBuffers(1, &VBO); // Vertex Buffer Object ID
		glGenBuffers(1, &EBO); // Element Buffer Object ID
		glGenVertexArrays(1, &VAO); // Vertex Array Object ID (needed when drawing the object)
	}

	void bindBuffersAndVertexArrayObject()
	{
		// 1. Bind Vertex Array Object
		glBindVertexArray(VAO);
		// 2. Copy our vertex array in a buffer for OpenGL to use
		glBindBuffer(GL_ARRAY_BUFFER, VBO);
		glBufferData(GL_ARRAY_BUFFER, stride * numVertices * sizeof(GLfloat), vertices, GL_STATIC_DRAW);
		// 3. Copy our index array in a element buffer for OpenGL to use
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, 3 * numTriangles * sizeof(GLuint), indices, GL_STATIC_DRAW);
	}

	void unbindBuffersAndVertexArrayObject()
	{
		// Do NOT unbind the index buffer while the VAO is still bound.
		// The index buffer is an essential part of the VAO state.
		glBindVertexArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, 0);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	}

	void configureVertexArrayAttributes()
	{
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
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, stride * sizeof(GLfloat), (void*)(3 * sizeof(GLfloat))); // Normals
		glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, stride * sizeof(GLfloat), (void*)(6 * sizeof(GLfloat))); // Texture coordinates
		glEnableVertexAttribArray(0); // Vertex coordinates
		glEnableVertexAttribArray(1); // Normals
		glEnableVertexAttribArray(2); // Texture coordinates
	}

	void clean()
	{
		// De-allocate vertex array and index array
		if (vertices) {
			delete[] vertices;
			vertices = nullptr;
		}
		if (indices) {
			delete[] indices;
			indices = nullptr;
		}
		
		// OPTIONAL: de-allocate all resources once they've outlived their purpose:
		if (glIsVertexArray(VAO)) {
			glDeleteVertexArrays(1, &VAO);
		}
		this->VAO = 0;

		if (glIsBuffer(VBO)) {
			glDeleteBuffers(1, &VBO);
		}
		this->VBO = 0;

		if (glIsBuffer(EBO)) {
			glDeleteBuffers(1, &EBO);
		}
		this->EBO = 0;

		this->numTriangles = 0;
		this->numVertices = 0;
	}
};

#endif