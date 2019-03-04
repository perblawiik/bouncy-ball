/***********************************************************************************/
/*** Note that the template of this class is taken from https://learnopengl.com/ ***/
/***********************************************************************************/

#pragma once
#ifndef SHADER_H
#define SHADER_H

#include <glad/glad.h>

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

class Shader
{
public:

	unsigned int ID;

	// Constructor generates the shader on the fly
	Shader(const char* vertexPath, const char* fragmentPath)
	{
		// 1. Retrieve the vertex/fragment source code from filePath
		std::string vertexCode;
		std::string fragmentCode;
		std::ifstream vShaderFile;
		std::ifstream fShaderFile;
		// Ensure ifstream objects can throw exceptions:
		vShaderFile.exceptions(std::ifstream::failbit | std::ifstream::badbit);
		fShaderFile.exceptions(std::ifstream::failbit | std::ifstream::badbit);
		try
		{
			// Open files
			vShaderFile.open(vertexPath);
			fShaderFile.open(fragmentPath);
			std::stringstream vShaderStream, fShaderStream;
			// Read file's buffer contents into streams
			vShaderStream << vShaderFile.rdbuf();
			fShaderStream << fShaderFile.rdbuf();
			// Close file handlers
			vShaderFile.close();
			fShaderFile.close();
			// Convert stream into string
			vertexCode = vShaderStream.str();
			fragmentCode = fShaderStream.str();
		}
		catch (std::ifstream::failure e)
		{
			std::cout << "ERROR::SHADER::FILE_NOT_SUCCESFULLY_READ" << std::endl;
		}
		const char* vShaderCode = vertexCode.c_str();
		const char * fShaderCode = fragmentCode.c_str();

		// 2. Compile shaders
		unsigned int vertex, fragment;
		// Vertex shader
		vertex = glCreateShader(GL_VERTEX_SHADER);
		glShaderSource(vertex, 1, &vShaderCode, NULL);
		glCompileShader(vertex);
		checkCompileErrors(vertex, "VERTEX");
		// Fragment Shader
		fragment = glCreateShader(GL_FRAGMENT_SHADER);
		glShaderSource(fragment, 1, &fShaderCode, NULL);
		glCompileShader(fragment);
		checkCompileErrors(fragment, "FRAGMENT");
		// Shader Program
		ID = glCreateProgram();
		glAttachShader(ID, vertex);
		glAttachShader(ID, fragment);
		glLinkProgram(ID);
		checkCompileErrors(ID, "PROGRAM");
		// Delete the shaders as they're linked into our program now and no longer necessary
		glDeleteShader(vertex);
		glDeleteShader(fragment);
	}

	// Activate the shader
	void use()
	{
		glUseProgram(ID);
	}

	/* Utility uniform functions */
	void setBool(const std::string &name, bool value) const
	{
		glUniform1i(glGetUniformLocation(ID, name.c_str()), (int)value);
	}
	void setBool(const GLint &uniformID, bool value) const
	{
		glUniform1i(uniformID, (int)value);
	}

	void setInt(const std::string &name, int value) const
	{
		glUniform1i(glGetUniformLocation(ID, name.c_str()), value);
	}
	void setInt(const GLint &uniformID, int value) const
	{
		glUniform1i(uniformID, value);
	}

	void setFloat(const std::string &name, float value) const
	{
		glUniform1f(glGetUniformLocation(ID, name.c_str()), value);
	}
	void setFloat(const GLint &uniformID, float value) const
	{
		glUniform1f(uniformID, value);
	}

	// Array
	void setFloat(const std::string &name, float values[], int size) const 
	{
		glUniform1fv(glGetUniformLocation(ID, name.c_str()), size, values);
	}
	void setFloat(const GLint &uniformID, float values[], int size) const
	{
		glUniform1fv(uniformID, size, values);
	}

	// Vec3 float
	void setFloat3(const std::string &name, float values[])
	{
		glUniform3f(glGetUniformLocation(ID, name.c_str()), values[0], values[1], values[2]);
	}
	void setFloat3(const GLint &uniformID, float values[])
	{
		glUniform3f(uniformID, values[0], values[1], values[2]);
	}

	void setVec3(const std::string &name, const float &x, const float &y, const float &z)
	{
		glUniform3f(glGetUniformLocation(ID, name.c_str()), x, y, z);
	}
	void setVec3(const GLint &uniformID, const float &x, const float &y, const float &z)
	{
		glUniform3f(uniformID, x, y, z);
	}

	void setFloatMat4(const std::string &name, float values[]) const
	{
		glUniformMatrix4fv(glGetUniformLocation(ID, name.c_str()), 1, GL_FALSE, values);
	}
	void setFloatMat4(const GLint &uniformID, float values[]) const
	{
		glUniformMatrix4fv(uniformID, 1, GL_FALSE, values);
	}

private:

	// Utility function for checking shader compilation/linking errors.
	void checkCompileErrors(unsigned int shader, std::string type)
	{
		int success;
		char infoLog[1024];
		if (type != "PROGRAM")
		{
			glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
			if (!success)
			{
				glGetShaderInfoLog(shader, 1024, NULL, infoLog);
				std::cout << "ERROR::SHADER_COMPILATION_ERROR of type: " << type << "\n" << infoLog << "\n -- --------------------------------------------------- -- " << std::endl;
			}
		}
		else
		{
			glGetProgramiv(shader, GL_LINK_STATUS, &success);
			if (!success)
			{
				glGetProgramInfoLog(shader, 1024, NULL, infoLog);
				std::cout << "ERROR::PROGRAM_LINKING_ERROR of type: " << type << "\n" << infoLog << "\n -- --------------------------------------------------- -- " << std::endl;
			}
		}
	}
};
#endif