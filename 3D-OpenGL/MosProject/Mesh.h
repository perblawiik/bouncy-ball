#pragma once
#ifndef MESH_H
#define MESH_H

#include <glad/glad.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

#include "Animation.h"

class Mesh 
{
public:

	Mesh()
		: meshIsEmpty(true), VAO(0), VBO(0), EBO(0), vertices(nullptr), indices(nullptr), numVertices(0), numTriangles(0), stride(8)
	{ }

	~Mesh()
	{
		this->clean();
	}

	void loadMeshData(const std::string &fileName)
	{
		// Clean up any previous mesh data
		this->clean();
		meshIsEmpty = false;

		// Input File Stream
		std::ifstream inFile;
		if (inFile) {
			// Temporary string to store each read line
			std::string line;
			std::string filePath = ("Simulations//Saves//" + fileName + ".mesh");

			// Open file for reading
			inFile.open(filePath, std::ifstream::in);

			// Get number of vertices
			std::getline(inFile, line);
			this->numVertices = std::stoi(line);
			// Get number of triangles
			std::getline(inFile, line);
			this->numTriangles = std::stoi(line);

			// Vertex array
			this->vertices = new GLfloat[this->numVertices * this->stride];
			// Index array
			this->indices = new GLuint[this->numTriangles * 3];

			int indexCounter = 0;
			int vertexCounter = 0;
			while (std::getline(inFile, line))
			{
				if (vertexCounter < (this->numVertices * this->stride)) {
					vertices[vertexCounter] = std::stof(line);
					++vertexCounter;
				}
				else if (indexCounter < (this->numTriangles * 3)) {
					indices[indexCounter] = std::stoi(line);
					++indexCounter;
				}
			}
			inFile.close();

			// Generate ID's for buffers and vertex array to use in the shader
			this->generateBuffersAndVertexArrayObject();

			// Store our vertex array and index array in buffers for OpenGL to use
			this->bindBuffersAndVertexArrayObject();

			// Tell OpenGL how it should interpret the vertex data (per vertex attribute)
			this->configureVertexArrayAttributes();

			// Deactivate (unbind) the VAO and the buffers again.
			this->unbindBuffersAndVertexArrayObject();
		}
	}

	void createPlaneXZ(const GLfloat &WIDTH, const GLfloat &HEIGHT, GLfloat textureWidth = -1.0f, GLfloat textureHeight = -1.0f)
	{
		// Clean up any previous mesh data
		this->clean();
		meshIsEmpty = false;

		GLfloat s, t;
		if (textureWidth > GLOBAL_CONSTANTS::EPSILON) {
			s = WIDTH / textureWidth;
		}
		else {
			s = 1.0f;
		}

		if (textureHeight > GLOBAL_CONSTANTS::EPSILON) {
			t = HEIGHT / textureHeight;
		}
		else {
			t = 1.0f;
		}

		const GLfloat vertexData[] = {
			// Position Coordinates                    // Normals          // Texture coordinates
			-(WIDTH / 2.0f), 0.0f, -(HEIGHT / 2.0f),   0.0f, 1.0f, 0.0f,   0.0f, 0.0f, // Upper left
			 (WIDTH / 2.0f), 0.0f, -(HEIGHT / 2.0f),   0.0f, 1.0f, 0.0f,   s, 0.0f, // Upper right
			 (WIDTH / 2.0f), 0.0f,  (HEIGHT / 2.0f),   0.0f, 1.0f, 0.0f,   s, t, // Lower right
			-(WIDTH / 2.0f), 0.0f,  (HEIGHT / 2.0f),   0.0f, 1.0f, 0.0f,   0.0f, t  // Lower left
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
		// Clean up any previous mesh data
		this->clean();
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

		this->vertices = new GLfloat[numVertices * stride]; // Initialize vertex array
		this->indices = new GLuint[numTriangles * 3]; // Initialize index array

		/** Generate vertex array **/
		// Bottom vertex
		vertices[0] = 0.0f; vertices[1] = -radius; vertices[2] = 0.0f; // Coordinates
		vertices[3] = 0.0f; vertices[4] = -1.0f; vertices[5] = 0.0f; // Normal
		vertices[6] = 0.5f; vertices[7] = 0.0f;

		const GLfloat PI = GLOBAL_CONSTANTS::PI;
		GLfloat sampleRate = PI / numHorizontalSegments; // Number of steps 
		GLfloat theta = -PI + sampleRate; // Go from bottom to top (Y � -PI < theta < PI )
		GLfloat phi = 0.0f; // Begin at Z = 0 (Z � 0 < phi < 2PI )
		
		// Generate middle part vertices with normals
		int index = 7; // Skip first 7 (the bottom vertex with normal and texture coordinates already specified)
		for (int i = 0; i < numHorizontalSegments - 1; ++i) {

			float Y = cos(theta); // Y-coordinate
			float R = sin(theta); // XZ-plane

			phi = 0.0f;
			for (int j = 0; j < numVerticalSegments; ++j) {
				// Vertex (x, y, z)
				vertices[++index] = radius * R * sin(phi);
				vertices[++index] = radius * Y;
				vertices[++index] = radius * R * cos(phi);
				// Normal (x, y, z)
				vertices[++index] = R * sin(phi);
				vertices[++index] = Y;
				vertices[++index] = R * cos(phi);
				// Texture Coordinates (s, t)
				vertices[++index] = phi / (2.0f * PI);
				vertices[++index] = 1.0f + (theta / PI);

				phi += sampleRate;
			}
			theta += sampleRate;
		}

		// Top vertex
		vertices[++index] = 0.0f; vertices[++index] = radius; vertices[++index] = 0.0f; // Coordinates
		vertices[++index] = 0.0f; vertices[++index] = 1.0f; vertices[++index] = 0.0f; // Normal
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
		indices[(numTriangles * 3) - 2] = (lastVertexIndex - 1);

		// Generate ID's for buffers and vertex array to use in the shader
		this->generateBuffersAndVertexArrayObject();

		// Store our vertex array and index array in buffers for OpenGL to use
		this->bindBuffersAndVertexArrayObject();

		// Tell OpenGL how it should interpret the vertex data (per vertex attribute)
		this->configureVertexArrayAttributes();
		
		// Deactivate (unbind) the VAO and the buffers again.
		this->unbindBuffersAndVertexArrayObject();
	}

	void createCylinder2 (int vertSegs, int horizSegs, const GLfloat &radius, const GLfloat &height) {

		this->clean();

		this->meshIsEmpty = false;

		if (horizSegs < 1) {
			horizSegs = 1;
		}
		if (vertSegs < 4) {
			vertSegs = 4;
		}

		this->numVertices = (4 * vertSegs) + 2 + (vertSegs * (horizSegs - 1));
		this->numTriangles = 2 * vertSegs * (1 + horizSegs);

		this->vertices = new GLfloat[numVertices * stride];
		this->indices = new GLuint[3 * numTriangles];

		// Bottom center
		// Vertex coordinates
		vertices[0] = 0.0f; vertices[1] = -(height / 2.0f); vertices[2] = 0.0f;
		// Normal coordinates
		vertices[3] = 0.0f; vertices[4] = -1.0f; vertices[5] = 0.0f;
		// Texture coordinates
		vertices[6] = 0.5f; vertices[7] = 0.5f;


		const GLfloat PI = GLOBAL_CONSTANTS::PI;
		// Go from bottom to top (Y � -PI < theta < PI )
		GLfloat theta = -PI;
		// Begin at Z = 0 (Z � 0 < phi < 2PI )
		GLfloat phi = 0.0f;

		int index = 7;
		// Generate vertices and normals for bottom circle plane (all normals should be (0.0, -1.0, 0.0))
		for (int j = 0; j < vertSegs; ++j) {

			// Vertex (x, y, z)
			vertices[++index] = radius * sin(phi);
			vertices[++index] = -(height / 2.0f); // The bottom circle is on the plane y = -height/2
			vertices[++index] = radius * cos(phi);
			// Normal (x, y, z)
			vertices[++index] = 0.0f;
			vertices[++index] = -1.0f;
			vertices[++index] = 0.0f;
			// Textures (s, t)
			vertices[++index] = cos(phi) * 0.5f + 0.5f; 
			vertices[++index] = sin(phi + PI) * 0.5f + 0.5f;

			phi += (2.0f * PI) / (GLfloat)vertSegs;
		}

		// Begin at Z = 0 (Z � 0 < phi < 2PI )
		phi = 0.0f;
		// Generate middle part vertices with normals (from bottom to top)
		for (int i = 0; i < horizSegs + 1; ++i) {

			GLfloat y = cos(theta);
			for (int j = 0; j < vertSegs; ++j) {

				// Vertex (x, y, z)
				vertices[++index] = radius * sin(phi);
				vertices[++index] = (height / 2.0f) * y;
				vertices[++index] = radius * cos(phi);
				// Normal (x, y, z)
				vertices[++index] = sin(phi);
				vertices[++index] = y;
				vertices[++index] = cos(phi);
				// Textures (s, t)
				vertices[++index] = phi/(2.0f*PI);
				vertices[++index] = abs(y * 0.5f - 0.5f);

				phi += (2.0f * PI) / (GLfloat)vertSegs;
			}
			phi = 0.0f;
			theta += PI / (GLfloat)horizSegs;
		}

		phi = 0.0f;
		// Generate vertices and normals for top circle plane (all normals should be (0.0, 1.0, 0.0))
		for (int j = 0; j < vertSegs; ++j) {

			// Vertex (x, y, z)
			vertices[++index] = radius * sin(phi);
			vertices[++index] = (height / 2.0f);
			vertices[++index] = radius * cos(phi);
			// Normal (x, y, z)
			vertices[++index] = 0.0f;
			vertices[++index] = 1.0f;
			vertices[++index] = 0.0f;
			// Textures (s, t)
			vertices[++index] = cos(phi) * 0.5f + 0.5f;
			vertices[++index] = sin(phi + PI) * 0.5f + 0.5f;

			phi += (2.0f * PI) / (GLfloat)vertSegs;
		}

		// Top center vertex, normal and texture coordinates
		vertices[++index] = 0.0f; vertices[++index] = (height / 2.0f); vertices[++index] = 0.0f;
		vertices[++index] = 0.0f; vertices[++index] = 1.0f; vertices[++index] = 0.0f;
		vertices[++index] = 0.5f; vertices[++index] = 0.5f;

		/* Generate Index Array */
		// Bottom circle plane
		index = -1;
		for (int i = 0; i < vertSegs; ++i) {

			indices[++index] = 0;
			if ((i + 2) <= vertSegs) {
				indices[++index] = i + 2;
			}
			else {
				indices[++index] = (i + 2) - vertSegs;
			}
	
			indices[++index] = i + 1;
		}

		// Middle part
		int v0 = vertSegs + 1;
		for (int i = 0; i < horizSegs; i++) {

			for (int j = 0; j < vertSegs - 1; ++j) {
				// One rectangle at a time (two triangles)
				indices[++index] = v0;
				indices[++index] = v0 + 1;
				indices[++index] = vertSegs + v0;
				indices[++index] = v0 + 1;
				indices[++index] = vertSegs + v0 + 1;
				indices[++index] = vertSegs + v0;
				++v0;
			}
			indices[++index] = v0;
			indices[++index] = (v0 + 1) - vertSegs;
			indices[++index] = vertSegs + v0;
			indices[++index] = (v0 + 1) - vertSegs;
			indices[++index] = v0 + 1;
			indices[++index] = vertSegs + v0;
			++v0;
		}

		// Top circle plane
		int lastVertexIndex = numVertices - 1;
		for (int i = 0; i < vertSegs; ++i) {

			indices[++index] = lastVertexIndex;
			if ((lastVertexIndex - 2 - i) >= lastVertexIndex - vertSegs) {
				indices[++index] = lastVertexIndex - 2 - i;
			}
			else {
				indices[++index] = lastVertexIndex - 1;
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

	void createCylinder(int vertSegs, const GLfloat &radius, const GLfloat &height) {

		this->clean();

		this->meshIsEmpty = false;

		int horizSegs = 1;

		if (vertSegs < 4) {
			vertSegs = 4;
		}

		this->numVertices = 6 * vertSegs + 2;
		this->numTriangles = 2 * vertSegs * (1 + horizSegs);

		this->vertices = new GLfloat[numVertices * stride];
		this->indices = new GLuint[3 * numTriangles];

		// Bottom center
		// Vertex coordinates
		vertices[0] = 0.0f; vertices[1] = -(height / 2.0f); vertices[2] = 0.0f;
		// Normal coordinates
		vertices[3] = 0.0f; vertices[4] = -1.0f; vertices[5] = 0.0f;
		// Texture coordinates
		vertices[6] = 0.5f; vertices[7] = 0.5f;

		const GLfloat PI = GLOBAL_CONSTANTS::PI;
		// Go from bottom to top (Y � -PI < theta < PI )
		GLfloat theta = -PI;
		// Begin at Z = 0 (Z � 0 < phi < 2PI )
		GLfloat phi = 0.0f;

		int index = 7;
		// Generate vertices and normals for bottom circle plane (all normals should be (0.0, -1.0, 0.0))
		for (int j = 0; j < vertSegs; ++j) {

			// Vertex (x, y, z)
			vertices[++index] = radius * sin(phi);
			vertices[++index] = -(height / 2.0f); // The bottom circle is on the plane y = -height/2
			vertices[++index] = radius * cos(phi);
			// Normal (x, y, z)
			vertices[++index] = 0.0f;
			vertices[++index] = -1.0f;
			vertices[++index] = 0.0f;
			// Textures (s, t)
			vertices[++index] = cos(phi) * 0.5f + 0.5f;
			vertices[++index] = sin(phi + PI) * 0.5f + 0.5f;

			phi += (2.0f * PI) / (GLfloat)vertSegs;
		}

		
		// Begin at Z = 0 (Z � 0 < phi < 2PI )
		phi = 0.0f;
		// Generate middle part vertices with normals (from bottom to top)
		for (int i = 0; i < horizSegs + 1; ++i) {

			GLfloat y = cos(theta);
			// Two vertices each iteration which belong to the same face
			for (int j = 0; j < vertSegs; ++j) {

				GLfloat phiNext = phi + (2.0f * PI) / (GLfloat)vertSegs;
				// First Vertex
				// Vertex (x, y, z)
				vertices[++index] = radius * sin(phi);
				vertices[++index] = (height / 2.0f) * y;
				vertices[++index] = radius * cos(phi);
				// Normal (x, y, z)
				vertices[++index] = sin(phi) + sin(phiNext);
				vertices[++index] = y + y;
				vertices[++index] = cos(phi) + cos(phiNext);
				
				// Textures (s, t)
				vertices[++index] = phi / (2.0f*PI);
				vertices[++index] = abs(y * 0.5f - 0.5f);

				// Second Vertex
				// Vertex (x, y, z)
				vertices[++index] = radius * sin(phiNext);
				vertices[++index] = (height / 2.0f) * y;
				vertices[++index] = radius * cos(phiNext);
				// Normal (x, y, z)
				vertices[++index] = sin(phi) + sin(phiNext);;
				vertices[++index] = y + y;
				vertices[++index] = cos(phi) + cos(phiNext);
				
				// Textures (s, t)
				vertices[++index] = phiNext / (2.0f*PI);
				vertices[++index] = abs(y * 0.5f - 0.5f);

				phi = phiNext;
			}
			vertices[index - 1] = 1.0f; // Last vertex s-texture coordinate is always 1

			phi = 0.0f;
			theta += PI / (GLfloat)horizSegs;
		}

		phi = 0.0f;
		// Generate vertices and normals for top circle plane (all normals should be (0.0, 1.0, 0.0))
		for (int j = 0; j < vertSegs; ++j) {

			// Vertex (x, y, z)
			vertices[++index] = radius * sin(phi);
			vertices[++index] = (height / 2.0f);
			vertices[++index] = radius * cos(phi);
			// Normal (x, y, z)
			vertices[++index] = 0.0f;
			vertices[++index] = 1.0f;
			vertices[++index] = 0.0f;
			// Textures (s, t)
			vertices[++index] = cos(phi) * 0.5f + 0.5f;
			vertices[++index] = sin(phi + PI) * 0.5f + 0.5f;

			phi += (2.0f * PI) / (GLfloat)vertSegs;
		}
		
		// Top center vertex, normal and texture coordinates
		vertices[++index] = 0.0f; vertices[++index] = (height / 2.0f); vertices[++index] = 0.0f;
		vertices[++index] = 0.0f; vertices[++index] = 1.0f; vertices[++index] = 0.0f;
		vertices[++index] = 0.5f; vertices[++index] = 0.5f;

		/* Generate Index Array */
		// Bottom circle plane
		index = -1;
		for (int i = 0; i < vertSegs; ++i) {

			indices[++index] = 0;
			if ((i + 2) <= vertSegs) {
				indices[++index] = i + 2;
			}
			else {
				indices[++index] = (i + 2) - vertSegs;
			}

			indices[++index] = i + 1;
		}

		// Middle part
		int v0 = vertSegs + 1;
		for (int i = 0; i < horizSegs; i++) {

			for (int j = 0; j < vertSegs; ++j) {
				// One rectangle at a time (two triangles)
				indices[++index] = v0;
				indices[++index] = v0 + 1;
				indices[++index] = 2 * vertSegs + v0;

				indices[++index] = v0 + 1;
				indices[++index] = 2* vertSegs + v0 + 1;
				indices[++index] = 2 * vertSegs + v0;
				v0 = v0 + 2;
			}
		}

		// Top circle plane
		int lastVertexIndex = numVertices - 1;
		for (int i = 0; i < vertSegs; ++i) {

			indices[++index] = lastVertexIndex;
			if ((lastVertexIndex - 2 - i) >= lastVertexIndex - vertSegs) {
				indices[++index] = lastVertexIndex - 2 - i;
			}
			else {
				indices[++index] = lastVertexIndex - 1;
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

	int getNumTriangles()
	{
		return this->numTriangles;
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