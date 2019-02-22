#pragma once
#ifndef SOFTBODY_H
#define SOFTBODY_H

#include <glad/glad.h>

#include <iostream>
#include <math.h>

#include "Matrix.h"
#include "Structs.h"
#include "Mesh.h"

class SoftBody 
{
public:

	SoftBody ()
		: mesh(nullptr)
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
		mesh = nullptr;
	}

	void setMesh(Mesh *newMesh)
	{
		mesh = newMesh;
	}

	void setupSimulationModel(Settings &s, const float &bondRegionRadius)
	{
		// Number of particles
		s.NUM_POINTS = mesh->getNumVertices();

		// Create a matrix containing the position coordinates of the mesh vertices
		this->createParticlePositionMatrix(s);

		// Generate an index table for the spring bonds between particles (Ex. bond between p1 and p2 get connection [1, 2])
		this->createBondIndexMatrix(s, bondRegionRadius);

		//Masses per particle
		Matrix masses(s.NUM_POINTS, 1); // Create Nx1 matrix
		for (int i = 0; i < masses.size(); ++i) {
			masses[i] = s.WEIGHT / (float)s.NUM_POINTS; // The whole mass is distributed equally
		}
		// Store the masses as a member of the class
		this->m = masses;

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
		Matrix zeros(s.NUM_POINTS, s.DIM);

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

			// Normalise it, used to give Fk and Fb direction
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

		// Approximating the new values using: X_n + 1 = X_n + h * X'_n (Eulers)
		// (they're not supressed for debugging purposes)
		V = V + (Vp * s.h);
		Fk = Fk + (Fkp * s.h);
		X = X + (V * s.h);

		// Set velocity to zero if particle hits the ground
		for (int j = 1; j <= s.NUM_POINTS; ++j) {

			// Check if particle position is below the ground
			if (X(j, 2) < (-s.RADIUS * 2.0f)) {
				V(j, 2) = 0.0f; // Set velocity to zero
				X(j, 2) = -s.RADIUS * 2.0f; // Move the position to ground level
			}
		}
	}

	void updatePositions()
	{
		mesh->updateVertexBufferData();
	}

	void render()
	{
		mesh->render();
	}

	// Returns an array containing the position coordinates of the particles
	GLfloat* getParticlePositionArray()
	{
		return this->X.getValues();
	}

	// Return mesh vertex array used with opengl
	GLfloat* getMeshVertexArray()
	{
		return mesh->getVertexArray();
	}

	// Return mesh index array used with opengl
	GLuint* getMeshIndexArray()
	{
		return mesh->getIndexArray();
	}

	int getNumMeshVertices()
	{
		return (mesh->getNumVertices());
	}

	int getNumMeshTriangles()
	{
		return (mesh->getNumTriangles());
	}

private:
	// Mesh for the softbody
	Mesh* mesh;

	// Simulation variables
	Matrix m; // The masses for each particle in the simulation
	Matrix X; // The coordinates (x, y, z) for each particle in the simulation
	Matrix I; // The spring bond indices for the particles in the simulation
	Matrix V; // The velocity of each particle in the simulation
	Matrix Vp; // The accelerations of each particle in the simulation
	Matrix Fk; // The spring force for the spring bonds in the simulation
	Matrix Fkp; // The derivative of Fk

	void createParticlePositionMatrix(Settings &s)
	{
		// Mesh vertex array
		GLfloat* vertices = mesh->getVertexArray();
		int stride = mesh->getStride();

		// Particle position coordinates (x, y, z) per particle
		Matrix positions(s.NUM_POINTS, s.DIM);
		for (int i = 0; i < s.NUM_POINTS; ++i) {
			positions[i * s.DIM] = vertices[i * stride];
			positions[i * s.DIM + 1] = vertices[i *  stride + 1];
			positions[i * s.DIM + 2] = vertices[i *  stride + 2];
		}
		// Store the positions as a member of the class
		this->X = positions;

		// Deactivate pointer
		vertices = nullptr;
	}

	void createBondIndexMatrix(Settings &s, const float &bondRegionRadius)
	{
		// Compute the maximum number of bonds if all particles are attached to each other
		int maxNumBonds = 0;
		for (int i = 1; i <= s.NUM_POINTS; ++i) {
			maxNumBonds += s.NUM_POINTS - i;
		}
		std::cout << "Maximum number of springs: " << maxNumBonds << std::endl;

		Matrix bondIndices(maxNumBonds, 2);
		// A counter for row index
		int row = 1;
		// Generate indice table for the spring bonds between particles (Ex. bond between p1 and p2 get connection [1, 2])
		for (int i = 1; i <= s.NUM_POINTS; ++i) {

			for (int j = i + 1; j <= s.NUM_POINTS; ++j) {

				float dist = 0.0f;
				// Calculate the euclidian distance between the particles
				dist = sqrt(pow(X(j, 1) - X(i, 1), 2) + pow(X(j, 2) - X(i, 2), 2) + pow(X(j, 3) - X(i, 3), 2));
				// Only attach a spring bond if distance is smaller than the given region
				if (dist < bondRegionRadius) {
					bondIndices(row, 1) = (float)i;
					bondIndices(row, 2) = (float)j;
					++row;
				}
			}
		}

		// If number of rows is smaller than the maximum number of spring bonds, create a smaller matrix
		if (row < maxNumBonds) {
			std::cout << "Number of springs used: " << row - 1 << std::endl;

			Matrix bondIndicesFinal(row - 1, 2);
			for (int i = 1; i <= bondIndicesFinal.numRows(); ++i) {
				bondIndicesFinal(i, 1) = bondIndices(i, 1);
				bondIndicesFinal(i, 2) = bondIndices(i, 2);
			}
			// Set number of spring bonds
			s.NUM_BONDS = row - 1;
			// Store the bond indices as a member of the class
			this->I = bondIndicesFinal;
		}
		else {
			// Set number of spring bonds
			s.NUM_BONDS = maxNumBonds;
			// Store the bond indices as a member of the class
			this->I = bondIndices;
		}
	}
};

#endif