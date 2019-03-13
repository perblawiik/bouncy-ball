#version 330 core
layout (location = 0) in vec3 aPosition;
layout (location = 1) in vec2 aTexCoords;

out vec2 TexCoords;

uniform vec3 translation; // Model View Matrix

void main()
{
	TexCoords = aTexCoords;
	vec3 pos = vec3(aPosition.x + translation.x, aPosition.y + translation.y, aPosition.z + translation.z);
    gl_Position = vec4(pos, 1.0);
} 