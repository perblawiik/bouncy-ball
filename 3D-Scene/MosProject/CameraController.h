#pragma once

#ifndef CAMERACONTROLLER_H
#define CAMERACONTROLLER_H

#include <GLFW/glfw3.h>

#include "Camera.h"

class CameraController
{
public:

	// Constructor
	CameraController(Camera* cam, GLFWwindow* window)
		: camera(cam), window(window), lastX(0.0), lastY(0.0), mouseSensitivity(10.0f), lastLeftClick(GL_FALSE), lastRightClick(GL_FALSE)
	{ 
		glfwGetCursorPos(window, &lastX, &lastY);
	}

	void setMouseSensitivity (const GLfloat &sens)
	{
		mouseSensitivity = sens;
	}

	// Read poll events
	void processInput(const GLfloat &dt)
	{
		if (glfwGetKey(window, GLFW_KEY_W)) {
			camera->translate(0.0f, 0.0f, -dt * 50.0f); // Forward
		}
		if (glfwGetKey(window, GLFW_KEY_S)) {
			camera->translate(0.0f, 0.0f, dt * 50.0f); // Backwards
		}
		if (glfwGetKey(window, GLFW_KEY_A)) {
			camera->translate(-dt * 50.0f, 0.0f, 0.0f); // Left
		}
		if (glfwGetKey(window, GLFW_KEY_D)) {
			camera->translate(dt * 50.0f, 0.0f, 0.0f); // Right
		}

		// Find cursor position
		double currentX, currentY;
		glfwGetCursorPos(window, &currentX, &currentY);
		int currentLeftClick = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT);
		int currentRightClick = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT);

		/* Check if any mouse button is being pressed down */
		if (currentLeftClick && lastLeftClick) {
			// Movement speed is set by deltaPosition
			double moveX = currentX - lastX;
			double moveY = currentY - lastY;

			camera->translate(-mouseSensitivity * dt  * (GLfloat)moveX, 0.0f, 0.0f); // Left/Right
			camera->translate(0.0f, mouseSensitivity * dt * (GLfloat)moveY, 0.0f); // Up/down
		}
		else if (currentRightClick && lastRightClick) {
			// Movement speed is set by deltaPosition
			double moveX = currentX - lastX;
			double moveY = currentY - lastY;

			camera->rotate(-mouseSensitivity * dt  * (GLfloat)moveY, -mouseSensitivity * dt * (GLfloat)moveX, 0.0f);
		}

		lastX = currentX;
		lastY = currentY;
		lastLeftClick = currentLeftClick;
		lastRightClick = currentRightClick;
	}

private:

	Camera* camera;
	GLFWwindow* window;

	double lastX;
	double lastY;
	GLfloat mouseSensitivity;

	int lastLeftClick;
	int lastRightClick;
};


#endif