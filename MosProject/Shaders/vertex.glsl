#version 330 core
layout (location = 0) in vec3 aPos;

uniform float time;

void main()
{

	float x = aPos.x + 0.2*sin(time + aPos.y*(3.14159265359/4));
	float y = aPos.y + 0.2*sin(time + aPos.x*(3.14159265359/4));
	float z = aPos.z;
    gl_Position = vec4(x, y, z, 1.0);
}