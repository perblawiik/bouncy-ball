#pragma once
#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>

class Matrix 
{

public:

	// Default constructor
	Matrix() : values(nullptr), rows(0), columns(0), length(0) {}

	// Constructor
	Matrix(const int rows, const int cols)
		: values(new float[rows*cols]), rows(rows), columns(cols), length(rows*cols)
	{
		for (int i = 0; i < length; ++i) {
			values[i] = 0.0f;
		}
	}

	// Copy constructor
	Matrix(const Matrix &m)
		: values(new float[m.length]), rows(m.rows), columns(m.columns), length(m.length)
	{
		for (int i = 0; i < length; ++i) {
			values[i] = m.values[i];
		}
	}

	// Clone function
	Matrix* clone() 
	{
		return new Matrix(*this);
	}

	// Assignment operator
	Matrix &operator=(const Matrix &m)
	{
		// Create a local copy
		Matrix copy(m);

		// Swap member variables
		std::swap(rows, copy.rows);
		std::swap(columns, copy.columns);
		std::swap(length, copy.length);
		std::swap(values, copy.values);

		// Return this instance to allow cascading
		return *this;
	}

	// Function call operator
	float &operator()(const int r, const int c) 
	{
		if (r > rows || rows < 1) {
			std::cout << "WARNING! Matrix row number out of bounds!" << std::endl;
		}
		else if (c > columns || c < 1) {
			std::cout << "WARNING! Matrix column number out of bounds!" << std::endl;
		}

		return values[(r-1)*columns + c-1];
	}

	// Subscript operator
	float &operator[] (const int index) 
	{
		if (index > length || index < 0) {
			std::cout << "WARNING! Matrix index out of bound!" << std::endl;
		}

		return values[index];
	}

	const int size() 
	{
		return length;
	}

	// Destructor
	~Matrix() 
	{
		delete[] values;
		values = nullptr;

		rows = 0;
		columns = 0;
		length = 0;
	}
private:

	float* values;
	int rows, columns, length;
};

#endif