#pragma once

#ifndef ANIMATION_H
#define ANIMATION_H

#include "Matrix.h"

class Animation
{
public:
	// Constructor
	Animation(Matrix* dataSet, const float &stepSize, const int &numSteps)
		: stepData(Matrix(1, dataSet->numColumns())), animationDataSet(dataSet), timeStep(stepSize), numSteps(numSteps), stepCounter(1), animationIsActive(false)
	{
		// Set step data to first row in the animation data matrix
		animationDataSet->copyRow(stepCounter, stepData);
	}

	~Animation()
	{
		if (animationDataSet) {
			delete animationDataSet;
			animationDataSet = nullptr;
		}
	}

	void update(const int speed = 1)
	{
		// Check if animation should play
		if (animationIsActive) {
			// Increment step counter to next time step in the simulation
			stepCounter += speed;

			// Check if step counter has reached maximum number of steps
			if (stepCounter > numSteps)
				stepCounter = 1;

			// Update current step data
			animationDataSet->copyRow(stepCounter, stepData);
		}
	}

	void startAnimation()
	{
		animationIsActive = true;
	}

	void stopAnimation()
	{
		animationIsActive = false;
	}

	float *getAnimationStepData()
	{
		return this->stepData.getValues();
	}

private:

	Matrix stepData;
	Matrix* animationDataSet;
	
	float timeStep;

	int numSteps;
	int stepCounter;

	bool animationIsActive;
};

#endif