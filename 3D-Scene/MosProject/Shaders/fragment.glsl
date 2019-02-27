#version 330 core

out vec4 FragColor;
in vec3 shadedColor;

void main()
{
    FragColor = vec4 (shadedColor, 1.0);
} 