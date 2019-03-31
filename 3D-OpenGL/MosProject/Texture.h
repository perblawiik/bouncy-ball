#pragma once
#ifndef TEXTURE_H
#define TEXTURE_H

#ifndef STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#endif

#include <glad/glad.h>
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

		// Generate texture
		// Arg 1. The texture target
		// Arg 2. The mipmap level
		// Arg 3. Kind of color format
		// Arg 4 - 5. Width and height of the texture
		// Arg 6. Some legacy stuff (should always be 0).
		// Arg 7 - 8. Specify the format and datatype of the source image
		// Arg 9. Image data
		if (data) {

			GLenum format;
			if (numChannels == 1)
				format = GL_RED;
			else if (numChannels == 3)
				format = GL_RGB;
			else if (numChannels == 4)
				format = GL_RGBA;

			// Bind texture
			glBindTexture(GL_TEXTURE_2D, ID);
			glTexImage2D(GL_TEXTURE_2D, 0, format, width, height, 0, format, GL_UNSIGNED_BYTE, data);
			glGenerateMipmap(GL_TEXTURE_2D);

			// Set the texture wrapping and filtering options for the currently bound texture
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, format == GL_RGBA ? GL_CLAMP_TO_EDGE : GL_REPEAT);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, format == GL_RGBA ? GL_CLAMP_TO_EDGE : GL_REPEAT);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
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