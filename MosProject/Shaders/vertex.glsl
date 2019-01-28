#version 330 core
layout (location = 0) in vec3 aPos;

uniform float time;
uniform float positions[34];

void main()
{
	int n = gl_VertexID;

	float x = positions[n* 2];
	float y = positions[n * 2 + 1];
	float z = 0.0f;

	vec3 pos = vec3(x, y, z);

    gl_Position = vec4(0.25*pos, 1.0);
}