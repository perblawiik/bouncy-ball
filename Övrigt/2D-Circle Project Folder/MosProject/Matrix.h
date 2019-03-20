#pragma once
#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <math.h>

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

	// Assignment operator
	Matrix &operator=(const Matrix &m)
	{
		// Create a local copy
		Matrix copy(m);

		// Swap member variables
		std::swap(this->rows, copy.rows);
		std::swap(this->columns, copy.columns);
		std::swap(this->length, copy.length);
		std::swap(this->values, copy.values);

		// Return this instance to allow cascading
		return *this;
	}

	// Function call operator - Returns the value at position (row, column)
	float &operator()(const int r, const int c) 
	{
		if (r > this->rows || r < 1) {
			std::cout << "WARNING! Matrix::Function Call Operator::Row number out of bounds!" << std::endl;
		}
		else if (c > this->columns || c < 1) {
			std::cout << "WARNING! Matrix::Function Call Operator::Column number out of bounds!" << std::endl;
		}

		return this->values[(r-1)*this->columns + c-1];
	}

	// Subscript operator - Returns the value from given array index 
	float &operator[] (const int index) 
	{
		if (index > this->length || index < 0) {
			std::cout << "WARNING! Matrix::Subscript Operator::Index out of bound!" << std::endl;
		}

		return this->values[index];
	}

	// Subtract operator
	Matrix operator-(const Matrix &rhs)
	{
		if (this->rows != rhs.rows || this->columns != rhs.columns) {
			std::cout << "WARNING! Matrix::Subtraction Operator::Dimensions dont match!" << std::endl;
		}

		Matrix result(this->rows, this->columns);

		for (int i = 0; i < this->length; ++i) {
			result.values[i] = this->values[i] - rhs.values[i];
		}

		return result;
	}

	// Addition operator
	Matrix operator+(const Matrix &rhs)
	{
		if (this->rows != rhs.rows || this->columns != rhs.columns) {
			std::cout << "WARNING! Matrix::Addition Operator::Dimensions dont match!" << std::endl;
		}

		Matrix result(this->rows, this->columns);

		for (int i = 0; i < this->length; ++i) {
			result.values[i] = this->values[i] + rhs.values[i];
		}

		return result;
	}

	// Multiply operator
	Matrix operator*(const float &k)
	{
		Matrix result(this->rows, this->columns);

		for (int i = 0; i < this->length; ++i) {
			result.values[i] = this->values[i] * k;
		}

		return result;
	}

	void copyValues(float out[])
	{
		for (int i = 0; i < this->length; ++i) {
			out[i] = this->values[i];
		} 
	}

	/*  Input parameters
	row_number: Integer that specifies which row to extract.
	out: Matrix (1xN) to copy the specified row into. */
	void copyRow(const int row_number, Matrix &out)
	{
		if (row_number < 1 || row_number > this->rows) {
			std::cout << "WARNING! Matrix::Copy Row Function::Row number out of bounds!" << std::endl;
		}

		for (int i = 1; i <= this->columns; ++i) {
			out.values[i - 1] = this->values[(row_number - 1) * this->columns + (i - 1)];
		}
	}

	void replaceRow(const int row_number, const Matrix &in)
	{
		if (row_number < 1 || row_number > this->rows) {
			std::cout << "WARNING! Matrix::Replace Row Function::Row number out of bounds!" << std::endl;
		}

		for (int i = 1; i <= this->columns; ++i) {
			this->values[(row_number - 1) * this->columns + (i - 1)] = in.values[i - 1];
		}
	}

	// Set magnitude to 1
	void normalize() 
	{
		float sum = 0.0f;
		for (int i = 0; i < this->length; ++i) {
			sum += this->values[i] * this->values[i];
		}
		sum = sqrt(sum);

		// Dont divide by zero
		if (sum > 0.00000001) {
			for (int i = 0; i < this->length; ++i) {
				this->values[i] = this->values[i] / sum;
			}
		}
	}

	// Returns the scalar product of this matrix and input matrix 
	// OBS! Works only for 1xN or Nx1 matrices
	float dot(const Matrix &rhs)
	{
		float dotProduct = 0.0f;

		for (int i = 0; i < this->length; ++i) {
			dotProduct += this->values[i] * rhs.values[i];
		}

		return dotProduct;
	}


	const int size()
	{
		return this->length;
	}

	const int numRows()
	{
		return this->rows;
	}

	const int numColumns()
	{
		return this->columns;
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