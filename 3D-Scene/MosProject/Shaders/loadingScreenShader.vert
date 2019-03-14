#version 330 core
layout (location = 0) in vec3 aPosition;
layout (location = 1) in vec2 aTexCoords;

out vec2 TexCoords;

uniform vec3 translation; 
uniform float time;
uniform bool quadShouldRotate;

void main()
{
	TexCoords = aTexCoords;

	float newX = aPosition.x;
	float newY = aPosition.y;

	if (quadShouldRotate) {
		float aspect = 16.0/9.0;
		newX = (aPosition.x * cos(time) * aspect - aPosition.y * sin(time)) / aspect;
		newY = aPosition.x * sin(time) * aspect + aPosition.y * cos(time);
	}

	vec3 pos = vec3 (
		newX + translation.x,
		newY + translation.y,
		translation.z + aPosition.z
	);

    gl_Position = vec4(pos, 1.0);
} 