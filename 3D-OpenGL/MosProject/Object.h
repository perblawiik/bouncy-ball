#pragma once

#ifndef OBJECT_H
#define OBJECT_H

#include <glad/glad.h>
#include <vector>

#include "Transform.h"
#include "Mesh.h"
#include "Texture.h"
#include "Shader.h"
#include "Animation.h"

class Object
{
public:

	Object()
		: mesh(nullptr), texture(nullptr), transform(Transform()), color(new GLfloat[3]), animationID(-1), staticPoseID(-1)
	{
		color[0] = 1.0f; color[1] = 1.0f; color[2] = 1.0f; // Set default color to white
	}

	~Object()
	{
		delete[] color;
		color = nullptr;

		if (!animations.empty()) {
			std::for_each(animations.begin(), animations.end(), this->deallocateAnimation);
			animations.clear();
			animationID = -1;
		}

		if (!staticPoses.empty()) {
			std::for_each(staticPoses.begin(), staticPoses.end(), this->deallocatePose);
			staticPoses.clear();
			staticPoseID = -1;
		}
	}

	void setMesh(Mesh* m)
	{
		mesh = m;
	}

	void setTexture(Texture* t)
	{
		texture = t;
	}

	void setColor(const GLfloat &r, const GLfloat &g, const GLfloat &b)
	{
		color[0] = r; color[1] = g; color[2] = b;
	}

	void addAnimation(const std::string &fileName)
	{
		// Input File Stream
		std::ifstream inFile;
		if (inFile) {
			// Temporary string to store each read line
			std::string line;
			std::string filePath = ("Simulations//Saves//" + fileName + ".sim");

			// Open file for reading
			inFile.open(filePath, std::ifstream::in);

			// Get simulation time step
			std::getline(inFile, line);
			float timeStep = std::stof(line);
			// Get number of steps in the simulation
			std::getline(inFile, line);
			int numSteps = std::stoi(line);
			// Get the size each vertex array
			std::getline(inFile, line);
			int numColumns = std::stoi(line);

			// Allocate memory for the data set
			Matrix* simulationData = new Matrix(numSteps, numColumns);

			int counter = 0;
			while (std::getline(inFile, line)) {

				(*simulationData)[counter] = std::stof(line);
				++counter;
			}
			inFile.close();

			this->animations.push_back(new Animation(simulationData, timeStep, numSteps));
			++animationID;
		}
	}

	void startAnimation(const int id)
	{
		if (!animations.empty()) {
			selectAnimation(id);
			animations[id]->startAnimation();
		}
	}

	void pauseAnimation(const int id)
	{
		if (!animations.empty()) {
			selectAnimation(id);
			animations[id]->stopAnimation();
		}
	}

	void setPosition(const GLfloat &x, const GLfloat &y, const GLfloat &z)
	{
		transform.setPosition(x, y, z);
	}

	void setRotation(const GLfloat &x, const GLfloat &y, const GLfloat &z)
	{
		transform.setRotation(x, y, z);
	}

	void setScale(const GLfloat &x, const GLfloat &y, const GLfloat &z)
	{
		transform.setScale(x, y, z);
	}

	int addStaticPose(const GLfloat matrix4[])
	{
		GLfloat *mat = new GLfloat[16];
		for (int i = 0; i < 16; ++i) {
			mat[i] = matrix4[i];
		}
		this->staticPoses.push_back(mat);
		++staticPoseID;
	
		return this->staticPoses.size() - 1;
	}

	void useStaticPose(const int id)
	{
		if (id < (int)staticPoses.size()) {
			this->staticPoseID = id;
		}
		else {
			std::cout << "Warning! Object::Function Use Static Pose::Pose ID doesnt exist!" << std::endl;
		}
		
	}

	void update()
	{
		if (!animations.empty()) {

			animations[animationID]->update();
			// Update the vertex array of the mesh from the animation step data
			for (int i = 0; i < mesh->getNumVertices(); ++i) {
				// Vertex coordinates
				mesh->getVertexArray()[i * 8] = animations[animationID]->getAnimationStepData()[i * 6];
				mesh->getVertexArray()[(i * 8) + 1] = animations[animationID]->getAnimationStepData()[(i * 6) + 1];
				mesh->getVertexArray()[(i * 8) + 2] = animations[animationID]->getAnimationStepData()[(i * 6) + 2];
				// Normal
				mesh->getVertexArray()[(i * 8) + 3] = animations[animationID]->getAnimationStepData()[(i * 6) + 3];
				mesh->getVertexArray()[(i * 8) + 4] = animations[animationID]->getAnimationStepData()[(i * 6) + 4];
				mesh->getVertexArray()[(i * 8) + 5] = animations[animationID]->getAnimationStepData()[(i * 6) + 5];
			}

			mesh->updateVertexBufferData();
		}
	}

	void render(Shader* shader)
	{
		shader->use();

		if (texture)
		{
			glActiveTexture(GL_TEXTURE0);
			texture->use();
		}
		if (mesh)
		{
			shader->setFloat3("objectColor", color);
			if (staticPoseID >= 0 && !staticPoses.empty()) {
				shader->setFloatMat4("modelView", staticPoses[staticPoseID]);
			}
			else {
				shader->setFloatMat4("modelView", transform.matrix4);
			}
			
			mesh->render();
		}
	}

	void renderStaticPoses(Shader* shader)
	{
		if (texture)
		{
			glActiveTexture(GL_TEXTURE0);
			texture->use();
		}
		if (mesh)
		{
			shader->setFloat3("objectColor", color);

			for (unsigned int i = 0; i < staticPoses.size(); ++i) {
				shader->setFloatMat4("modelView", staticPoses[i]);
				mesh->render();
			}
		}
	}

private:

	Mesh* mesh;
	Texture* texture;
	Transform transform;
	GLfloat* color;

	std::vector<Animation*> animations;
	int animationID;

	std::vector<GLfloat*> staticPoses;
	int staticPoseID;

	void selectAnimation(const unsigned int ID)
	{
		// Check if ID exists
		if (ID >= 0 && ID < animations.size()) {
			animationID = ID;
		}
		else {
			std::cout << "Warning! Mesh::Select Animation::Animation ID not found!" << std::endl;
		}
	}

	static void deallocateAnimation(Animation* a)
	{
		delete a;
		a = nullptr;
	}
	static void deallocatePose(GLfloat* array)
	{
		delete[] array;
		array = nullptr;
	}
};

#endif
