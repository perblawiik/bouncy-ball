#pragma once

#ifndef TRANSFORM_H
#define TRANSFORM_H

#include <glad/glad.h>
#include <iostream>

#include "Structs.h"

class Transform
{
public:

	GLfloat* position;
	GLfloat* rotation;
	GLfloat* scale;
	GLfloat* matrix4;

	Transform()
		: position(new GLfloat[3]), rotation(new GLfloat[3]), scale(new GLfloat[3]), matrix4(new GLfloat[16]), parent(nullptr), inverted(false)
	{
		for (int i = 0; i < 3; ++i) {
			position[i] = 0.0f;
			rotation[i] = 0.0f;
			scale[i] = 1.0f;
		}

		MATRIX4::identity(matrix4);
	}

	void setScale(const GLfloat &x, const GLfloat &y, const GLfloat &z)
	{
		scale[0] = x;
		scale[1] = y;
		scale[2] = z;

		this->composeMatrix();
	}

	void setPosition(const GLfloat &x, const GLfloat &y, const GLfloat &z)
	{
		position[0] = x;
		position[1] = y;
		position[2] = z;

		this->composeMatrix();
	}

	void setRotation(const GLfloat &x, const GLfloat &y, const GLfloat &z)
	{
		rotation[0] = x;
		rotation[1] = y;
		rotation[2] = z;

		this->composeMatrix();
	}

	void setParent(Transform* parentTransform)
	{
		this->parent = parentTransform;
		this->composeMatrix();
	}

	void invert()
	{
		if (inverted) {
			inverted = false;
		}
		else {
			inverted = true;
		}
	}

private:

	Transform* parent;

	bool inverted;

	void composeMatrix()
	{
		GLfloat dummy[16];

		/*** First, add scale matrix ***/
		MATRIX4::scale(matrix4, scale[0]);


		/*** Second, add rotation matrix ***/
		// 1. Create the rotation matrix by applying rotation to x, y, z one at the time
		GLfloat rotMat[16];
		MATRIX4::identity(rotMat);
		MATRIX4::rotateZ(dummy, (rotation[2] * GLOBAL_CONSTANTS::PI / 180.0f));
		MATRIX4::multiply(dummy, rotMat, rotMat);
		MATRIX4::rotateX(dummy, (rotation[0] * GLOBAL_CONSTANTS::PI / 180.0f));
		MATRIX4::multiply(dummy, rotMat, rotMat);
		MATRIX4::rotateY(dummy, (rotation[1] * GLOBAL_CONSTANTS::PI / 180.0f));
		// 2. Multiply the final rotation matrix to the scale matrix
		MATRIX4::multiply(dummy, rotMat, rotMat);
		MATRIX4::multiply(rotMat, matrix4, matrix4);

		/*** Third, translation matrix ***/
		MATRIX4::translate(dummy, position[0], position[1], position[2]);
		MATRIX4::multiply(dummy, matrix4, matrix4);

		if (parent) {
			MATRIX4::multiply(parent->matrix4, matrix4, matrix4);
		}

		if (inverted) { // Invert the final transformation matrix
			MATRIX4::invert(matrix4, matrix4);
		}
	}
};


#endif