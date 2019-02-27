#pragma once
#ifndef STRUCTS_H
#define STRUCTS_H

#include <glad/glad.h>


/** STRUCTS **/

// Settings for the softbody simulation
struct Settings
{
	GLfloat h; // Step
	GLfloat k; // Spring constant
	GLfloat b; // Resistance constant
	GLfloat g; // Gravitation constant
	GLfloat RADIUS; // Radius of the sphere
	GLfloat WEIGHT; // Total weight of the sphere
	GLfloat TIME_DURATION; // Determines how long the simulation should be (given in seconds)

	GLint NUM_BONDS;
	GLint NUM_POINTS;
	GLint DIM;
	GLint NUM_STEPS;
};

// Struct with all constant values for the simulation
struct GLOBAL_CONSTANTS 
{
	struct window
	{
		static const int WIDTH;
		static const int HEIGHT;
	};

	static const float PI;
};

const int GLOBAL_CONSTANTS::window::WIDTH = 1280;
const int GLOBAL_CONSTANTS::window::HEIGHT = 720;
const float GLOBAL_CONSTANTS::PI = 3.14159265359f;


// Transformation matrices
struct MATRIX4 
{
	// Creates a translation matrix
	static void translate(GLfloat M[], const GLfloat &x, const GLfloat &y, const GLfloat &z) 
	{
		M[0] = 1.0f; M[4] = 0.0f; M[8] = 0.0f; M[12] = x;
		M[1] = 0.0f; M[5] = 1.0f; M[9] = 0.0f; M[13] = y;
		M[2] = 0.0f; M[6] = 0.0f; M[10] = 1.0f; M[14] = z;
		M[3] = 0.0f; M[7] = 0.0f; M[11] = 0.0f; M[15] = 1.0f;
	}

	// Creates an identity matrix
	static void identity(GLfloat M[]) 
	{
		M[0] = 1.0f; M[4] = 0.0f; M[8] = 0.0f;  M[12] = 0.0f;
		M[1] = 0.0f; M[5] = 1.0f; M[9] = 0.0f;  M[13] = 0.0f;
		M[2] = 0.0f; M[6] = 0.0f; M[10] = 1.0f;  M[14] = 0.0f;
		M[3] = 0.0f; M[7] = 0.0f; M[11] = 0.0f;  M[15] = 1.0f;
	}

	// Creates a scale matrix
	static void scale(GLfloat M[], const GLfloat &scale) 
	{
		M[0] = scale; M[4] = 0.0f;  M[8] = 0.0f;  M[12] = 0.0f;
		M[1] = 0.0f;  M[5] = scale; M[9] = 0.0f;  M[13] = 0.0f;
		M[2] = 0.0f;  M[6] = 0.0f;  M[10] = scale; M[14] = 0.0f;
		M[3] = 0.0f;  M[7] = 0.0f;  M[11] = 0.0f;  M[15] = 1.0f;
	}

	// Creates a rotation matrix in X-direction
	static void rotateX(GLfloat M[], const GLfloat &angle)
	{
		M[0] = 1.0f; M[4] = 0.0f;       M[8] = 0.0f;        M[12] = 0.0f;
		M[1] = 0.0f; M[5] = cos(angle); M[9] = -sin(angle); M[13] = 0.0f;
		M[2] = 0.0f; M[6] = sin(angle); M[10] = cos(angle);  M[14] = 0.0f;
		M[3] = 0.0f; M[7] = 0.0f;       M[11] = 0.0f;        M[15] = 1.0f;
	}

	// Creates a rotation matrix in Y-direction
	static void rotateY(GLfloat M[], const GLfloat &angle)
	{
		M[0] = cos(angle);  M[4] = 0.0f; M[8] = sin(angle); M[12] = 0.0f;
		M[1] = 0.0f;        M[5] = 1.0f; M[9] = 0.0f;       M[13] = 0.0f;
		M[2] = -sin(angle); M[6] = 0.0f; M[10] = cos(angle); M[14] = 0.0f;
		M[3] = 0.0f;        M[7] = 0.0f; M[11] = 0.0f;       M[15] = 1.0f;
	}

	// Creates a rotation matrix in Z-direction
	static void rotateZ(GLfloat M[], const GLfloat &angle)
	{
		M[0] = cos(angle); M[4] = -sin(angle); M[8] = 0.0f; M[12] = 0.0f;
		M[1] = sin(angle); M[5] = cos(angle);  M[9] = 0.0f; M[13] = 0.0f;
		M[2] = 0.0f;       M[6] = 0.0f;        M[10] = 1.0f; M[14] = 0.0f;
		M[3] = 0.0f;       M[7] = 0.0f;        M[11] = 0.0f; M[15] = 1.0f;
	}

	// Inverts a given matrix and returns the matrix combination
	static void invert(GLfloat out[], GLfloat a[]) {
		GLfloat a00 = a[0],
			    a01 = a[1],
			    a02 = a[2],
			    a03 = a[3];

		GLfloat a10 = a[4],
			    a11 = a[5],
			    a12 = a[6],
			    a13 = a[7];

		GLfloat a20 = a[8],
			    a21 = a[9],
			    a22 = a[10],
			    a23 = a[11];

		GLfloat a30 = a[12],
			    a31 = a[13],
			    a32 = a[14],
			    a33 = a[15];

		GLfloat b00 = a00 * a11 - a01 * a10;
		GLfloat b01 = a00 * a12 - a02 * a10;
		GLfloat b02 = a00 * a13 - a03 * a10;
		GLfloat b03 = a01 * a12 - a02 * a11;
		GLfloat b04 = a01 * a13 - a03 * a11;
		GLfloat b05 = a02 * a13 - a03 * a12;
		GLfloat b06 = a20 * a31 - a21 * a30;
		GLfloat b07 = a20 * a32 - a22 * a30;
		GLfloat b08 = a20 * a33 - a23 * a30;
		GLfloat b09 = a21 * a32 - a22 * a31;
		GLfloat b10 = a21 * a33 - a23 * a31;
		GLfloat b11 = a22 * a33 - a23 * a32;

		// Calculate the determinant
		GLfloat det = b00 * b11 - b01 * b10 + b02 * b09 + b03 * b08 - b04 * b07 + b05 * b06;

		if (abs(det) > 0.001f) 
			det = 1.0f / det;

		out[0] = (a11 * b11 - a12 * b10 + a13 * b09) * det;
		out[1] = (a02 * b10 - a01 * b11 - a03 * b09) * det;
		out[2] = (a31 * b05 - a32 * b04 + a33 * b03) * det;
		out[3] = (a22 * b04 - a21 * b05 - a23 * b03) * det;
		out[4] = (a12 * b08 - a10 * b11 - a13 * b07) * det;
		out[5] = (a00 * b11 - a02 * b08 + a03 * b07) * det;
		out[6] = (a32 * b02 - a30 * b05 - a33 * b01) * det;
		out[7] = (a20 * b05 - a22 * b02 + a23 * b01) * det;
		out[8] = (a10 * b10 - a11 * b08 + a13 * b06) * det;
		out[9] = (a01 * b08 - a00 * b10 - a03 * b06) * det;
		out[10] = (a30 * b04 - a31 * b02 + a33 * b00) * det;
		out[11] = (a21 * b02 - a20 * b04 - a23 * b00) * det;
		out[12] = (a11 * b07 - a10 * b09 - a12 * b06) * det;
		out[13] = (a00 * b09 - a01 * b07 + a02 * b06) * det;
		out[14] = (a31 * b01 - a30 * b03 - a32 * b00) * det;
		out[15] = (a20 * b03 - a21 * b01 + a22 * b00) * det;
	}

	// M is the matrix we want to create (an output argument )
	// vertFov is the vertical field of view (in the y direction )
	// aspect is the aspect ratio of the viewport ( width / height )
	// zNear is the distance to the near clip plane ( znear > 0)
	// zFar is the distance to the far clip plane ( zfar > znear )
	static void perspective(GLfloat M[], const GLfloat &vertFov, const GLfloat &aspect, const GLfloat &zNear, const GLfloat &zFar) 
	{
		GLfloat f = cos(vertFov / 2) / sin(vertFov / 2);

		M[0] = f / aspect; M[4] = 0.0f; M[8] = 0.0f;                              M[12] = 0.0f;
		M[1] = 0.0f;       M[5] = f;    M[9] = 0.0f;                              M[13] = 0.0f;
		M[2] = 0.0f;       M[6] = 0.0f; M[10] = -(zFar + zNear) / (zFar - zNear); M[14] = -(2 * zNear*zFar) / (zFar - zNear);
		M[3] = 0.0f;       M[7] = 0.0f; M[11] = -1.0f;                            M[15] = 0.0f;
	}

	// Performs a matrix multiplication
	static void multiply(const GLfloat M1[], const GLfloat M2[], GLfloat result[])
	{
		GLfloat temp[16];
		int index = 0;
		int i = 0;

		while (i < 16) {

			// Each column
			for (int k = 0; k < 4; ++k) {

				temp[index] = (M1[k] * M2[i]) + (M1[k + 4] * M2[i + 1]) + (M1[k + 8] * M2[i + 2]) + (M1[k + 12] * M2[i + 3]);
				++index;
			}
			i = i + 4; // Jump to next row (go pass 4 indexes)
		}

		// Copy the result to one of the matrices
		for (int n = 0; n < 16; ++n) {

			result[n] = temp[n];
		}
	}
};

#endif