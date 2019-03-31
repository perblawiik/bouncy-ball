#pragma once

#ifndef CAMERA_H
#define CAMERA_H

#include <glad/glad.h>
#include <vector>

#include "Transform.h"
#include "Shader.h"

class Camera
{
public:
	Camera(Shader* sh)
		: transform(new Transform())
	{ 
		// When we move the camera we actually move the world. 
		// This means that the transform matrix should be inverted.
		shaders.push_back(sh);

		transform->invert();
		this->updateUniformMatrix();
	}

	~Camera()
	{
		delete transform;
		transform = nullptr;

		shaders.clear();
	}

	void addShader(Shader* sh)
	{
		shaders.push_back(sh);
	}

	// Set position vector
	void setPosition(const GLfloat &x, const GLfloat &y, const GLfloat &z)
	{
		this->transform->setPosition(x, y, z);
		this->updateUniformMatrix();
	}

	// Set rotation vector
	void setRotation(const GLfloat &x, const GLfloat &y, const GLfloat &z)
	{
		this->transform->setRotation(x, y, z);
		this->updateUniformMatrix();
	}

	// Rotate camera
	void rotate(const GLfloat &x, const GLfloat &y, const GLfloat &z) 
	{
		this->transform->setRotation(
			this->transform->rotation[0] + x,
			this->transform->rotation[1] + y,
			this->transform->rotation[2] + z
		);
		this->updateUniformMatrix();
	}

	// Move the camera position
	void translate(const GLfloat &x, const GLfloat &y, const GLfloat &z) 
	{
		GLfloat rightVec[3] = {
			this->transform->matrix4[0] * x,
			this->transform->matrix4[4] * x,
			this->transform->matrix4[8] * x
		};

		GLfloat upVec[3] = {
			this->transform->matrix4[1] * y,
			this->transform->matrix4[5] * y,
			this->transform->matrix4[9] * y
		};

		GLfloat forwardVec[3] = {
			this->transform->matrix4[2] * z,
			this->transform->matrix4[6] * z,
			this->transform->matrix4[10] * z
		};

		this->transform->setPosition(
			this->transform->position[0] + rightVec[0] + upVec[0] + forwardVec[0],
			this->transform->position[1] + rightVec[1] + upVec[1] + forwardVec[1],
			this->transform->position[2] + rightVec[2] + upVec[2] + forwardVec[2]
		);

		this->updateUniformMatrix();
	}

	// Update camera view matrix in the shader
	void updateUniformMatrix()
	{
		for (unsigned int i = 0; i < shaders.size(); ++i) {
			shaders[i]->use();
			shaders[i]->setFloatMat4("cameraView", transform->matrix4);
		}
	}

	GLfloat* getPosition()
	{
		return transform->position;
	}

	Transform* getTransform()
	{
		return transform;
	}

private:
	Transform* transform;
	std::vector<Shader*> shaders;
};

#endif