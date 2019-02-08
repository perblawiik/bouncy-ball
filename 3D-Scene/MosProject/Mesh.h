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
		: VAO(0), vertices(nullptr), indices(nullptr), numVertices(0), numTriangles(0), stride(6)
	{
		meshIsEmpty = true;
	}

	// Copy constructor
	Mesh(const Mesh &m)
		: VAO(m.VAO), vertices(new GLfloat[m.numVertices * m.stride]), indices(new GLuint[m.numTriangles * 3]), numVertices(m.numVertices), numTriangles(m.numTriangles), stride(m.stride)
	{
		meshIsEmpty = false;

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
		std::swap(this->VAO, copy.VAO);
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
		// OPTIONAL: de-allocate all resources once they've outlived their purpose:
		glDeleteVertexArrays(1, &VAO);

		// De-allocate vertex array and index array
		delete[] vertices;
		vertices = nullptr;
		delete[] indices;
		indices = nullptr;
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
		std::cout << "Vertex array size: " << 3 * numVertices << std::endl;
		std::cout << "Radius: " << radius << std::endl;

		// Columns per row in the vertex array (coordinates + normals)
		this->vertices = new GLfloat[numVertices * stride]; // Initialize vertex array
		this->indices = new GLuint[numTriangles * 3]; // Initialize index array

		/** Generate vertex array **/
		// Bottom vertex
		vertices[0] = 0.0f; vertices[1] = -radius; vertices[2] = 0.0f; // Coordinates
		vertices[3] = 0.0f; vertices[4] = -radius; vertices[5] = 0.0f; // Normal

		const GLfloat PI = 3.14159265359f;
		GLfloat sampleRate = PI / numHorizontalSegments; // Number of steps 
		GLfloat theta = -PI + sampleRate; // Go from bottom to top (Y � -PI < theta < PI )
		GLfloat phi = 0; // Begin at Z = 0 (Z � 0 < phi < 2PI )

		// Generate middle part vertices with normals
		int index = 5; // Skip first 6 (the bottom vertex with normal already specified)
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

				phi += sampleRate;
			}
			theta += sampleRate;
		}

		// Top vertex
		vertices[++index] = 0.0f; vertices[++index] = radius; vertices[++index] = 0.0f; // Coordinates
		vertices[++index] = 0.0f; vertices[++index] = radius; vertices[++index] = 0.0f; // Normal

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

		GLuint VBO, EBO;
		glGenBuffers(1, &VBO); // Vertex Buffer Object ID
		glGenBuffers(1, &EBO); // Element Buffer Object ID
		glGenVertexArrays(1, &VAO); // Vertex Array Object ID (needed when drawing the object)

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

	void render()
	{
		if (meshIsEmpty == false) {
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
	GLfloat* vertices; // Vertex array
	GLuint* indices; // Index array

	int stride; // Number of elements per row in the vertex array
	int numTriangles; // Number of triangels of the mesh
	int numVertices; // Number of vertices of the mesh
};

#endif