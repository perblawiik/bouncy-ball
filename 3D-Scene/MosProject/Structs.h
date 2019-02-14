#pragma once
#ifndef STRUCTS_H
#define STRUCTS_H

#include <glad/glad.h>


/** STRUCTS **/
// Struct with all constant values for the simulation

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

struct GLOBAL_CONSTANTS {

	struct window
	{
		static const int WIDTH;
		static const int HEIGHT;
	};

	/*
	enum window
	{
		WIDTH = 800,
		HEIGHT = 600
	};
	*/

	static const float PI;
};

const int GLOBAL_CONSTANTS::window::WIDTH = 800;
const int GLOBAL_CONSTANTS::window::HEIGHT = 600;
const float GLOBAL_CONSTANTS::PI = 3.14159265359f;

#endif