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

	GLint NUM_BONDS;
	GLint NUM_POINTS;
	GLint DIM;
	GLint NUM_STEPS;
};

#endif