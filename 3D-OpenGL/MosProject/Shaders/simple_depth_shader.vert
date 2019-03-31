#version 330 core
layout (location = 0) in vec3 aPosition;

uniform mat4 lightSpaceMatrix;
uniform mat4 modelView;
uniform vec3 objectColor;

void main()
{
    gl_Position = lightSpaceMatrix * modelView * vec4(aPosition, 1.0);
}  