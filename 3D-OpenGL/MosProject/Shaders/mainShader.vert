#version 330 core
layout (location = 0) in vec3 aPosition;
layout (location = 1) in vec3 aNormal;
layout (location = 2) in vec2 aTexCoords;

// Time variable
uniform float time;

// Transformation matrices
uniform mat4 perspective; // Projection Matrix
uniform mat4 modelView; // Model View Matrix
uniform mat4 cameraView; // Camera View Matrix

out vec2 TexCoords;
out vec3 Normal;
out vec3 FragPosition;


void main()
{
	// Final transformation ( Perspective multiplied with the model view )
    mat4 T = perspective * cameraView * modelView;
	// Transform (x,y,z) vertex coordinates with a 4x4 matrix T
    gl_Position = T * vec4(aPosition, 1.0);

	Normal = normalize(mat3(modelView)*aNormal);
	TexCoords = aTexCoords;
	FragPosition = vec3(modelView * vec4(aPosition, 1.0));
}