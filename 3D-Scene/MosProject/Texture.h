#pragma once
#ifndef TEXTURE_H
#define TEXTURE_H

#ifndef STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#endif

#include <glad/glad.h>
#include <string>
#include <iostream>


class Texture
{
public:

	GLuint ID;

	Texture(const char* path)
	{
		// Load texture image
		unsigned char *data = stbi_load(path, &width, &height, &numChannels, 0);

		// Create texture ID
		glGenTextures(1, &ID);

		// Bind texture
		glBindTexture(GL_TEXTURE_2D, ID);

		// Set the texture wrapping and filtering options for the currently bound texture
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

		// Generate texture
		// Arg 1. The texture target
		// Arg 2. The mipmap level
		// Arg 3. Kind of color format
		// Arg 4 - 5. Width and height of the texture
		// Arg 6. Some legacy stuff (should always be 0).
		// Arg 7 - 8. Specify the format and datatype of the source image
		// Arg 9. Image data
		if (data) {
			glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, data);
			glGenerateMipmap(GL_TEXTURE_2D);
		}
		else {
			std::cout << "Error when loading texture!" << std::endl;
		}

		// Since we are done using the image data we can deallocate the memory
		stbi_image_free(data);
	}

	void use()
	{
		glBindTexture(GL_TEXTURE_2D, ID);
	}

private:

	int width;
	int height;
	int numChannels;
};


#endif