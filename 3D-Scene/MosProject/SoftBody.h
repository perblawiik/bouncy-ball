#pragma once
#ifndef SOFTBODY_H
#define SOFTBODY_H

#include <glad/glad.h>

#include <iostream>
#include <math.h>

#include "Matrix.h"
#include "Structs.h"

class SoftBody 
{
public:

	SoftBody ()
		: VBO(0), VAO(0), EBO(0), vertices(nullptr), indices(nullptr), numVertices(0), numTriangles(0), stride(0)
	{ 
		m = Matrix();
		X = Matrix();
		I = Matrix();
		V = Matrix();
		Vp = Matrix();
		Fk = Matrix();
		Fkp = Matrix();
	}

	~SoftBody()
	{
		delete[] vertices;
		vertices = nullptr;
		delete[] indices;
		indices = nullptr;

		// OPTIONAL: de-allocate all resources once they've outlived their purpose:
		glDeleteVertexArrays(1, &VAO);
		glDeleteBuffers(1, &VBO);
		glDeleteBuffers(1, &EBO);
	}

	// Generates a sphere mesh with a specified number of horizontal segments and radius as input parameters.
	void createSphere(const int segments, const float &radius)
	{
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
		this->stride = 6;
		this->vertices = new GLfloat[numVertices * stride]; // Initialize vertex array
		this->indices = new GLuint[numTriangles * 3]; // Initialize index array
	
		/** Generate vertex array **/
		// Bottom vertex
		vertices[0] = 0.0f; vertices[1] = -radius; vertices[2] = 0.0f; // Coordinates
		vertices[3] = 0.0f; vertices[4] = -radius; vertices[5] = 0.0f; // Normal

		const GLfloat PI = 3.14159265359f;
		GLfloat sampleRate = PI / numHorizontalSegments; // Number of steps 
		GLfloat theta = -PI + sampleRate; // Go from bottom to top (Y € -PI < theta < PI )
		GLfloat phi = 0; // Begin at Z = 0 (Z € 0 < phi < 2PI )

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

		// Vertex Buffer Object, Vertex Array Object, Element Buffer Object
		//GLuint VBO, VAO, EBO;
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
	}

	

	void setupSimulationModel(Settings &s)
	{

		s.NUM_POINTS = numVertices;

		//Masses per particle
		Matrix masses(s.NUM_POINTS, 1); // Create Nx1 matrix
		for (int i = 0; i < masses.size(); ++i) {
			masses[i] = s.WEIGHT / (float)s.NUM_POINTS; // All masses divided equally
		}
		// Store the masses as a member of the class
		this->m = masses;

		//Particle x, y Pos [Xx Xy] / per particle
		Matrix positions(s.NUM_POINTS, s.DIM);
		for (int i = 0; i < s.NUM_POINTS; ++i) {
			positions[i * s.DIM] = vertices[i * stride];
			positions[i * s.DIM + 1] = vertices[i * stride + 1];
			positions[i * s.DIM + 2] = vertices[i * stride + 2];
		}
		// Store the positions as a member of the class
		this->X = positions;

		// Indice table for the spring bonds between particles (Ex. bond between p1 and p2 get connection [1, 2])
		s.NUM_BONDS = 0;
		int NUM_P = 10;
		for (int i = 1; i <= s.NUM_POINTS; ++i) {
			s.NUM_BONDS += s.NUM_POINTS - i;
		}
		std::cout << "Number of springs: " << s.NUM_BONDS << std::endl;

		Matrix bondIndices(s.NUM_BONDS, 2);
		int ROW = 1;
		for (int i = 1; i <= s.NUM_POINTS; ++i) {

			for (int j = i + 1; j <= s.NUM_POINTS; ++j) {

				bondIndices(ROW, 1) = (float)i;
				bondIndices(ROW, 2) = (float)j;
				++ROW;
			}
		}
		// Store the bond indices as a member of the class
		this->I = bondIndices;

		// Starting velocity [Vx Vy]/ per particle
		Matrix velocities(s.NUM_POINTS, s.DIM); // All set to zero by default
		for (int i = 1; i < s.NUM_POINTS; ++i) {
			//velocities(i, 1) = 5.0f;
			//velocities(i, 2) = -50.0f;
		}
		// Store the velocities as a member of the class
		this->V = velocities;

		// Acceleration dV/dt 
		Matrix accelerations(s.NUM_POINTS, s.DIM); // All set to zero by default
		// Store the accelerations as a member of the class
		this->Vp = accelerations;

		// Fk spring starting force [F]/ per spring  (applied directionally later, depending on spring orientation) 
		Matrix springForces(s.NUM_BONDS, s.DIM); // All set to zero by default
		// Store the spring forces as a member of the class
		this->Fk = springForces;

		// Fkp = dFp/dt
		Matrix springForceDerivatives(s.NUM_BONDS, s.DIM); // All set to zero by default
		// Store the derivatives of the spring forces as a member of the class
		this->Fkp = springForceDerivatives;

	}

	void updateSimulationModel(const Settings &s)
	{
		// Used to set Vp to zero each cycle
		Matrix zeros(numVertices, s.DIM);

		// Dummy vectors (1xDIM matrix)
		Matrix vec1(1, s.DIM);
		Matrix vec2(1, s.DIM);

		// Vector from point 1 to point 2;
		Matrix diff(1, s.DIM);

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
		for (int r = 1; r <= Vp.numRows(); ++r) {
			Vp(r, 2) = Vp(r, 2) - s.g;
		}

		// Approximating the new values using: X_n + 1 = X_n + h * X'_n
		// (they're not supressed for debugging purposes)
		V = V + (Vp * s.h);
		Fk = Fk + (Fkp * s.h);
		X = X + (V * s.h);

		// Code that flips Y - ward velocity when the particle has Xy < -4.0 
		for (int j = 1; j <= s.NUM_POINTS; ++j) {

			if (X(j, 2) < (-s.RADIUS * 2)) {
				V(j, 2) = 0.0f;
				X(j, 2) = -s.RADIUS * 2;
			}
		}
	}

	void render()
	{
		// Bind the VAO
		glBindVertexArray(this->VAO);
		// Draw the mesh (mode, vertex count, type, element array buffer offset)
		glDrawElements(GL_TRIANGLES, this->numTriangles * 3, GL_UNSIGNED_INT, 0);
		// Unbind the VAO
		glBindVertexArray(0);
	}

	// Returns an array containing the position coordinates of the particles
	GLfloat* getParticlePositionArray()
	{
		return this->X.getValues();
	}

private:
	GLuint VBO; // Vertex Buffer Object ID
	GLuint VAO; // Vertex Array Object ID
	GLuint EBO; // Element Buffer Object ID

	GLfloat* vertices; // Vertex array
	GLuint* indices; // Index array

	int stride; // Number of elements per row in the vertex array
	int numTriangles; // Number of triangels of the mesh
	int numVertices; // Number of vertices of the mesh

	// Simulation variables
	Matrix m; // The masses for each particle in the simulation
	Matrix X; // The coordinates (x, y, z) for each particle in the simulation
	Matrix I; // The spring bond indices for the particles in the simulation
	Matrix V; // The velocity of each particle in the simulation
	Matrix Vp; // The accelerations of each particle in the simulation
	Matrix Fk; // The spring force for the spring bonds in the simulation
	Matrix Fkp; // The derivative of Fk
};

#endif