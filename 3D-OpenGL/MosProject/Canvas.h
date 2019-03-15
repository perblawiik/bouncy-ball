#pragma once

#ifndef CANVAS_H
#define CANVAS_H

#include <glad/glad.h>
#include "Shader.h"
#include "Texture.h"


class Canvas
{
public:

	Canvas(int* width, int* height)
		: windowWidth(width), windowHeight(height), 
		VAO(0), VBO(0), EBO(0), 
		vertices(nullptr), indices(nullptr), 
		numVertices(0), numTriangles(0), 
		empty(true), rotationIsActive(false),
		texture(nullptr), position(new GLfloat[3])
	{
		position[0] = 0.0f; position[1] = 0.0f; position[2] = 0.0f;
	}

	~Canvas()
	{
		if (vertices) {
			delete[] vertices;
			vertices = nullptr;
		}
		if (indices) {
			delete[] indices;
			indices = nullptr;
		}

		delete[] position;
		position = nullptr;
	}

	void createRectangle(const int xDim, const int yDim)
	{
		empty = false;

		GLfloat rectWidth = (GLfloat)xDim;
		GLfloat rectHeight = (GLfloat)yDim;

		if (rectWidth > (*windowWidth)) {
			rectWidth = 1.0f;
		}
		else {
			rectWidth = rectWidth / (GLfloat)(*windowWidth);
		}
		if (rectHeight > (*windowHeight)) {
			rectHeight = 1.0f;
		}
		else {
			rectHeight = rectHeight / (GLfloat)(*windowHeight);
		}

		// Loading Screen
		GLfloat rectVertices[] = {
			//(x, y, z)                      (s, t)
			-rectWidth,  rectHeight, 0.0f,   0.0f, 0.0f,
			 rectWidth,  rectHeight, 0.0f,	 1.0f, 0.0f,
			 rectWidth, -rectHeight, 0.0f,   1.0f, 1.0f,
			-rectWidth, -rectHeight, 0.0f,   0.0f, 1.0f
		};
		GLuint rectIndices[] = {
			0, 3, 1, // Left triangle
			1, 3, 2 // Right triangle
		};

		this->numVertices = 4;
		this->numTriangles = 2;
		int stride = 5;

		this->vertices = new GLfloat[numVertices * stride];
		this->indices = new GLuint[3 * numTriangles];

		for (int i = 0; i < (numVertices * stride); ++i) {
			vertices[i] = rectVertices[i];
		}
		for (int i = 0; i < (3 * numTriangles); ++i) {
			indices[i] = rectIndices[i];
		}

		// 1. Generate buffers and vertex array object ids
		glGenBuffers(1, &VBO); // Vertex Buffer Object ID
		glGenBuffers(1, &EBO); // Element Buffer Object ID
		glGenVertexArrays(1, &VAO); // Vertex Array Object ID (needed when drawing the object)
		// 2. Bind buffers and vertex array object
		glBindVertexArray(VAO);
		glBindBuffer(GL_ARRAY_BUFFER, VBO);
		glBufferData(GL_ARRAY_BUFFER, stride * numVertices * sizeof(GLfloat), vertices, GL_STATIC_DRAW);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, 3 * numTriangles * sizeof(GLuint), indices, GL_STATIC_DRAW);
		// 3. Configure vertex array attributes
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, stride * sizeof(GLfloat), (void*)0);
		glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, stride * sizeof(GLfloat), (void*)(3 * sizeof(GLfloat))); // Texture coordinates
		glEnableVertexAttribArray(0); // Vertex coordinates
		glEnableVertexAttribArray(1); // Texture coordinates
		// 4. Unbind buffers and vertex array object
		glBindVertexArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, 0);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	}

	void setTexture(Texture* t)
	{
		this->texture = t;
	}

	void setPosition(const GLfloat &x, const GLfloat &y, const GLfloat &z)
	{
		position[0] = x / (GLfloat) (*windowWidth); 
		position[1] = y / (GLfloat) (*windowHeight); 
		position[2] = z;
	}

	void useRotationAnimation(Shader *shader, bool flag)
	{
		rotationIsActive = flag;
	}

	void render(Shader *shader)
	{
		if (!empty) {

			shader->use();

			if (texture) {
				texture->use();
			}

			shader->setBool("quadShouldRotate", rotationIsActive);
			shader->setVec3("translation", position[0], position[1], position[2]);

			glBindVertexArray(VAO); // Bind the VAO
			glDrawElements(GL_TRIANGLES, 3 * numTriangles, GL_UNSIGNED_INT, 0); // Draw
			glBindVertexArray(0); // Unbind the VAO
		}
	}

private:

	int* windowWidth;
	int* windowHeight;

	GLuint VAO;
	GLuint VBO;
	GLuint EBO;

	GLfloat* vertices;
	GLuint* indices;

	int numVertices;
	int numTriangles;

	bool empty;
	bool rotationIsActive;

	Texture* texture;
	GLfloat* position;
};

#endif