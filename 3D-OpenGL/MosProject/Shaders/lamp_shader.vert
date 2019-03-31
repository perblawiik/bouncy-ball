#version 330 core
layout (location = 0) in vec3 aPosition;

// Transformation matrices
uniform mat4 perspective; // Projection Matrix
uniform mat4 modelView; // Model View Matrix
uniform mat4 cameraView; // Camera View Matrix

void main()
{
    gl_Position = perspective * cameraView * modelView * vec4(aPosition, 1.0);
} 