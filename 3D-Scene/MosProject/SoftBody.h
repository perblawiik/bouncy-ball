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

	// Default constructor
	SoftBody ()
		: mesh(nullptr), m(0.0f), X(Matrix()), I(Matrix()), V(Matrix()), Vp(Matrix()), Fk(Matrix()), Fkp(Matrix()), normals(Matrix())
	{ }

	~SoftBody()
	{
		mesh = nullptr;
	}

	void setMesh(Mesh *m)
	{
		mesh = m;
	}

	void setupSimulationModel(Settings &s, const float &bondRegionRadius)
	{
		// Number of particles
		s.NUM_POINTS = mesh->getNumVertices();

		// Create a matrix containing the position coordinates of the mesh vertices
		this->createParticlePositionMatrix(s);

		// Create a matrix to store all the normals for each particle
		this->normals = X; // The normals are identical to the initial vertices

		// Generate an index table for the spring bonds between particles (Ex. bond between p1 and p2 get connection [1, 2])
		this->createBondIndexMatrix(s, bondRegionRadius);

		// Set the mass for each particle (total weight divided by number of particles)
		this->m = s.WEIGHT / (float)s.NUM_POINTS;

		// Starting velocity [Vx Vy]/ per particle
		Matrix velocities(s.NUM_POINTS, s.DIM); // All set to zero by default
		for (int i = 1; i < s.NUM_POINTS/2; ++i) {
			velocities(i, 1) = 0.0f; // x
			velocities(i, 2) = 30.0f; // y
			velocities(i, 3) = -30.0f; // z
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
		// Dummy vectors (1xDIM matrix)
		Matrix vec1(1, s.DIM);
		Matrix vec2(1, s.DIM);

		// Vector from point 1 to point 2;
		Matrix diff(1, s.DIM);

		// Set to zero so the components from each connected spring can be += and added separately
		Matrix zeros(s.NUM_POINTS, s.DIM);
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
			float F = Fk(n, 1); // Spring force
			Vp.replaceRow(index1, vec1 - (diff * (1.0f / m) * (s.b*dV + F)));
			Vp.replaceRow(index2, vec2 + (diff * (1.0f / m) * (s.b*dV + F)));

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

		// Update the normals to fit the new particle positions
		this->updateNormals(s);

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

	GLfloat* getParticleNormalArray()
	{
		return this->normals.getValues();
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
	Matrix X; // The coordinates (x, y, z) for each particle
private:
	// Mesh for the softbody
	Mesh* mesh;

	// Simulation variables
	GLfloat m; // The mass for each particle
	
	Matrix I; // The spring bond indices for the particles
	Matrix V; // The velocity of each particle
	Matrix Vp; // The accelerations of each particle
	Matrix Fk; // The spring force for the spring bonds
	Matrix Fkp; // The derivative of Fk
	Matrix normals; // The normals for each particle

	void createParticlePositionMatrix(Settings &s)
	{
		// Mesh vertex array
		GLfloat* vertices = mesh->getVertexArray();
		int stride = mesh->getStride();

		// Particle position coordinates (x, y, z) per particle
		this->X = Matrix(s.NUM_POINTS, s.DIM);
		for (int i = 0; i < s.NUM_POINTS; ++i) {
			X[i * s.DIM] = vertices[i * stride];
			X[i * s.DIM + 1] = vertices[i *  stride + 1];
			X[i * s.DIM + 2] = vertices[i *  stride + 2];
		}

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

	void updateNormals(const Settings &s)
	{
		// Mean position coordinates
		GLfloat Mx = 0.0f, My = 0.0f, Mz = 0.0f;

		// Compute normals for the particles connected to the bottom vertex
		this->computeBottomCapNormals(s.NUM_SEGS, Mx, My, Mz);

		// Compute normals for the middle part
		this->computeMiddleNormals(s.NUM_SEGS, Mx, My, Mz);
		
	}

	void computeBottomCapNormals(const int &H_SEGS, GLfloat &Mx, GLfloat &My, GLfloat &Mz)
	{
		this->computeBottomNormal(H_SEGS, Mx, My, Mz);

		int V_SEGS = H_SEGS * 2;
		GLfloat Fx, Fy, Fz, Lx, Ly, Lz, Bx, By, Bz, Rx, Ry, Rz;
		// 2. The normals for all vertices connected to the bottom vertex is calculated 
		//    by using the mean point from the 4-neighbours.
		int n = 1;
		int i = 2;
		for ( ; i < V_SEGS + 1; ++i) {

			Mx = 0.0f; My = 0.0f; Mz = 0.0f;
			int forward = i + 1;
			int left = i + V_SEGS;
			int right = 1;
			int back = i + V_SEGS - n;
			n = V_SEGS + 1;

			// Forward vector
			Fx = X(forward, 1) - X(i, 1);
			Fy = X(forward, 2) - X(i, 2);
			Fz = X(forward, 3) - X(i, 3);
			// Left vector
			Lx = X(left, 1) - X(i, 1);
			Ly = X(left, 2) - X(i, 2);
			Lz = X(left, 3) - X(i, 3);
			// Backwards vector
			Bx = X(back, 1) - X(i, 1);
			By = X(back, 2) - X(i, 2);
			Bz = X(back, 3) - X(i, 3);
			// Right vector
			Rx = X(right, 1) - X(i, 1);
			Ry = X(right, 2) - X(i, 2);
			Rz = X(right, 3) - X(i, 3);

			// Calculate the cross product (Forward x Left) to get the normal
			Mx += (Fy * Lz) - (Fz * Ly);
			My += -((Fx * Lz) - (Fz * Lx));
			Mz += (Fx * Ly) - (Fy * Lx);
			// Calculate the cross product (Left x Backwards) to get the normal
			Mx += (Ly * Bz) - (Lz * By);
			My += -((Lx * Bz) - (Lz * Bx));
			Mz += (Lx * By) - (Ly * Bx);
			// Calculate the cross product (Backwards x Right) to get the normal
			Mx += (By * Rz) - (Bz * Ry);
			My += -((Bx * Rz) - (Bz * Rx));
			Mz += (Bx * Ry) - (By * Rx);
			// Calculate the cross product (Right x Forward) to get the normal
			Mx += (Ry * Fz) - (Rz * Fy);
			My += -((Rx * Fz) - (Rz * Fx));
			Mz += (Rx * Fy) - (Ry * Fx);

			// The final normal is the sum of all face normals
			normals(i, 1) = Mx;
			normals(i, 2) = My;
			normals(i, 3) = Mz;
		} 

		Mx = 0.0f; My = 0.0f; Mz = 0.0f;
		int forward = i - V_SEGS;
		int left = i + V_SEGS;
		int right = 1;
		int back = i + V_SEGS - n;

		// Forward vector
		Fx = X(forward, 1) - X(i, 1);
		Fy = X(forward, 2) - X(i, 2);
		Fz = X(forward, 3) - X(i, 3);
		// Left vector
		Lx = X(left, 1) - X(i, 1);
		Ly = X(left, 2) - X(i, 2);
		Lz = X(left, 3) - X(i, 3);
		// Backwards vector
		Bx = X(back, 1) - X(i, 1);
		By = X(back, 2) - X(i, 2);
		Bz = X(back, 3) - X(i, 3);
		// Right vector
		Rx = X(right, 1) - X(i, 1);
		Ry = X(right, 2) - X(i, 2);
		Rz = X(right, 3) - X(i, 3);

		// Calculate the cross product (Forward x Left) to get the normal
		Mx += (Fy * Lz) - (Fz * Ly);
		My += -((Fx * Lz) - (Fz * Lx));
		Mz += (Fx * Ly) - (Fy * Lx);
		// Calculate the cross product (Left x Backwards) to get the normal
		Mx += (Ly * Bz) - (Lz * By);
		My += -((Lx * Bz) - (Lz * Bx));
		Mz += (Lx * By) - (Ly * Bx);
		// Calculate the cross product (Backwards x Right) to get the normal
		Mx += (By * Rz) - (Bz * Ry);
		My += -((Bx * Rz) - (Bz * Rx));
		Mz += (Bx * Ry) - (By * Rx);
		// Calculate the cross product (Right x Forward) to get the normal
		Mx += (Ry * Fz) - (Rz * Fy);
		My += -((Rx * Fz) - (Rz * Fx));
		Mz += (Rx * Fy) - (Ry * Fx);

		// The final normal is the sum of all face normals
		normals(i, 1) = Mx;
		normals(i, 2) = My;
		normals(i, 3) = Mz;
	}

	void computeBottomNormal(const int &H_SEGS, GLfloat &Mx, GLfloat &My, GLfloat &Mz)
	{
		int V_SEGS = H_SEGS * 2;
		// 1. The bottom vertex normal is calculated first by using the sum of all face normals
		//    of all connected vertices.
		Mx = My = Mz = 0.0f;
		GLfloat Ax, Ay, Az, Bx, By, Bz;
		int i = 1;
		for (; i < V_SEGS; ++i) {
			// Vector A (vector between the mean point and one of the connected vertices)
			Ax = X(i + 1, 1) - X(1, 1);
			Ay = X(i + 1, 2) - X(1, 2);
			Az = X(i + 1, 3) - X(1, 3);

			// Vector B (vector between the mean point and other one of the connected vertices)
			Bx = X(i + 2, 1) - X(1, 1);
			By = X(i + 2, 2) - X(1, 2);
			Bz = X(i + 2, 3) - X(1, 3);

			// Calculate the cross product (B x A) to get the normal
			Mx += (By * Az) - (Bz * Ay); // x-coordinate
			My += -((Bx * Az) - (Bz * Ax)); // y-coordinate
			Mz += (Bx * Ay) - (By * Ax); // z-coordinate
		}
		// Vector A (vector between the mean point and one of the connected vertices)
		Ax = X(i + 1, 1) - X(1, 1);
		Ay = X(i + 1, 2) - X(1, 2);
		Az = X(i + 1, 3) - X(1, 3);

		// Vector B (vector between the mean point and other one of the connected vertices)
		Bx = X(2, 1) - X(1, 1);
		By = X(2, 2) - X(1, 2);
		Bz = X(2, 3) - X(1, 3);

		// Calculate the cross product (B x A) to get the normal
		Mx += (By * Az) - (Bz * Ay); // x-coordinate
		My += -((Bx * Az) - (Bz * Ax)); // y-coordinate
		Mz += (Bx * Ay) - (By * Ax); // z-coordinate

		// The bottom normal is the sum of all face normals
		normals(1, 1) = Mx; // x-coordinate
		normals(1, 2) = My; // y-coordinate
		normals(1, 3) = Mz; // z-coordinate
	}

	void computeMiddleNormals(const int &H_SEGS, GLfloat &Mx, GLfloat &My, GLfloat &Mz)
	{
		int V_SEGS = H_SEGS * 2;
		GLfloat Fx, Fy, Fz, Lx, Ly, Lz, Bx, By, Bz, Rx, Ry, Rz;

		for (int j = 1; j <= H_SEGS - 3; ++j) {

			int n = 1;
			int i = (V_SEGS+2)*j;
			for (; i < V_SEGS + 1; ++i) {

				Mx = 0.0f; My = 0.0f; Mz = 0.0f;
				int forward = i + 1;
				int left = i + V_SEGS;
				int right = i - V_SEGS;
				int back = i + V_SEGS - n;
				n = V_SEGS + 1;

				// Forward vector
				Fx = X(forward, 1) - X(i, 1);
				Fy = X(forward, 2) - X(i, 2);
				Fz = X(forward, 3) - X(i, 3);
				// Left vector
				Lx = X(left, 1) - X(i, 1);
				Ly = X(left, 2) - X(i, 2);
				Lz = X(left, 3) - X(i, 3);
				// Backwards vector
				Bx = X(back, 1) - X(i, 1);
				By = X(back, 2) - X(i, 2);
				Bz = X(back, 3) - X(i, 3);
				// Right vector
				Rx = X(right, 1) - X(i, 1);
				Ry = X(right, 2) - X(i, 2);
				Rz = X(right, 3) - X(i, 3);

				// Calculate the cross product (Forward x Left) to get the normal
				Mx += (Fy * Lz) - (Fz * Ly);
				My += -((Fx * Lz) - (Fz * Lx));
				Mz += (Fx * Ly) - (Fy * Lx);
				// Calculate the cross product (Left x Backwards) to get the normal
				Mx += (Ly * Bz) - (Lz * By);
				My += -((Lx * Bz) - (Lz * Bx));
				Mz += (Lx * By) - (Ly * Bx);
				// Calculate the cross product (Backwards x Right) to get the normal
				Mx += (By * Rz) - (Bz * Ry);
				My += -((Bx * Rz) - (Bz * Rx));
				Mz += (Bx * Ry) - (By * Rx);
				// Calculate the cross product (Right x Forward) to get the normal
				Mx += (Ry * Fz) - (Rz * Fy);
				My += -((Rx * Fz) - (Rz * Fx));
				Mz += (Rx * Fy) - (Ry * Fx);

				// The final normal is the sum of all face normals
				normals(i, 1) = Mx;
				normals(i, 2) = My;
				normals(i, 3) = Mz;
			}

			Mx = 0.0f; My = 0.0f; Mz = 0.0f;
			int forward = i - V_SEGS;
			int left = i + V_SEGS;
			int right = i - V_SEGS;;
			int back = i + V_SEGS - n;

			// Forward vector
			Fx = X(forward, 1) - X(i, 1);
			Fy = X(forward, 2) - X(i, 2);
			Fz = X(forward, 3) - X(i, 3);
			// Left vector
			Lx = X(left, 1) - X(i, 1);
			Ly = X(left, 2) - X(i, 2);
			Lz = X(left, 3) - X(i, 3);
			// Backwards vector
			Bx = X(back, 1) - X(i, 1);
			By = X(back, 2) - X(i, 2);
			Bz = X(back, 3) - X(i, 3);
			// Right vector
			Rx = X(right, 1) - X(i, 1);
			Ry = X(right, 2) - X(i, 2);
			Rz = X(right, 3) - X(i, 3);

			// Calculate the cross product (Forward x Left) to get the normal
			Mx += (Fy * Lz) - (Fz * Ly);
			My += -((Fx * Lz) - (Fz * Lx));
			Mz += (Fx * Ly) - (Fy * Lx);
			// Calculate the cross product (Left x Backwards) to get the normal
			Mx += (Ly * Bz) - (Lz * By);
			My += -((Lx * Bz) - (Lz * Bx));
			Mz += (Lx * By) - (Ly * Bx);
			// Calculate the cross product (Backwards x Right) to get the normal
			Mx += (By * Rz) - (Bz * Ry);
			My += -((Bx * Rz) - (Bz * Rx));
			Mz += (Bx * Ry) - (By * Rx);
			// Calculate the cross product (Right x Forward) to get the normal
			Mx += (Ry * Fz) - (Rz * Fy);
			My += -((Rx * Fz) - (Rz * Fx));
			Mz += (Rx * Fy) - (Ry * Fx);

			// The final normal is the sum of all face normals
			normals(i, 1) = Mx;
			normals(i, 2) = My;
			normals(i, 3) = Mz;
		}

	}

	void computeTopNormals(const int &segs, GLfloat &Mx, GLfloat &My, GLfloat &Mz)
	{

	}
};

#endif