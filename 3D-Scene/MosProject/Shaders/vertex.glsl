#version 330 core
layout (location = 0) in vec3 Position;
layout (location = 1) in vec3 Normal;

uniform float positions[18];
uniform mat4 P; // Projection Matrix
uniform mat4 MV; // Model View Matrix

out vec3 shadedColor;

vec3 viewDirection = vec3(0.0, 0.0, 1.0); // the view direction - (0 ,0 ,1) in view space
float n = 50.0; // the " shininess " parameter
vec3 ka = vec3(0.0, 0.0, 0.0); // the ambient reflection color
vec3 Ia = vec3(0.6, 0.6, 0.6); // the ambient illumination color
vec3 kd = vec3(1.0, 0.5, 0.0); // the diffuse surface reflection color
vec3 Id = vec3(1.0, 1.0, 1.0); // the diffuse illumination color
vec3 ks = vec3(1.0, 1.0, 1.0); // the specular surface reflection color
vec3 Is = vec3(1.0, 1.0, 1.0); // the specular illumination color
//This assumes that N, L and V are normalized .

void main()
{
	int vertexId = gl_VertexID;
	float x = positions[vertexId * 3];
	float y = positions[vertexId * 3 + 1];
	float z = positions[vertexId * 3 + 2];
	vec3 simulationPos = vec3(x, y, z);

	// Final transformation ( Perspective multiplied with the model view )
    mat4 T = P * MV;

    // To avoid translation of the normals use only 3x3 from our 4x4 transformation matrix
    vec3 interpolatedNormal = normalize(mat3(MV)*Normal);

    vec3 lightDirection = normalize(vec3(1.0, 0.0, 1.0));

    // Reflection direction
    vec3 R = normalize(2.0* dot(interpolatedNormal,lightDirection)*interpolatedNormal - lightDirection); // Could also have used the function reflect ()
    float dotNL = max(dot(interpolatedNormal,lightDirection), 0.0) ; // If negative , set to zero
    float dotRV = max(dot(R,viewDirection), 0.0) ;
    if ( dotNL == 0.0) dotRV = 0.0; // Do not show highlight on the dark side
    shadedColor = Ia*ka + Id*kd* dotNL + Is*ks* pow (dotRV , n);

	// Transform (x,y,z) vertex coordinates with a 4x4 matrix T
    gl_Position = T * vec4(Position, 1.0);
}