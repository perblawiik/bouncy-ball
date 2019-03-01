#version 330 core
layout (location = 0) in vec3 Position;
layout (location = 1) in vec3 Normal;
layout (location = 2) in vec2 TextCoordinates;

uniform float time;
uniform mat4 P; // Projection Matrix
uniform mat4 MV; // Model View Matrix
uniform mat4 CV; // Camera View Matrix
uniform vec3 surfaceColor;
uniform vec3 cameraPosition;

out vec3 shadedColor;
out vec2 texCoords;

out vec3 worldPosition;
out vec3 cameraPos;
out vec3 reflection;
out vec3 ambientDiffuse;
out float isOnDarkSide;

vec3 ka = vec3(0.0, 0.0, 0.0); // the ambient reflection color
vec3 Ia = vec3(1.0, 1.0, 1.0); // the ambient illumination color
vec3 kd = vec3(1.0, 0.5, 0.0); // the diffuse surface reflection color
vec3 Id = vec3(1.0, 1.0, 1.0); // the diffuse illumination color
//This assumes that N, L and V are normalized .

void main()
{
	// Final transformation ( Perspective multiplied with the model view )
    mat4 T = P * CV * MV;

	//vec3 viewDirection = normalize(Position - cameraPosition);
    vec3 normalizedNormal = normalize(mat3(MV)*Normal);
    vec3 lightDirection = normalize(vec3(0.0, 1.0, 1.0));

	kd = surfaceColor;
    // Reflection direction
    reflection = normalize(2.0* dot(normalizedNormal,lightDirection)*normalizedNormal - lightDirection); // Could also have used the function reflect ()
    float dotNL = max(dot(normalizedNormal,lightDirection), 0.0) ; // If negative , set to zero
    //float dotRV = max(dot(R,viewDirection), 0.0) ;
    //if ( dotNL == 0.0) dotRV = 0.0; // Do not show highlight on the dark side


    //shadedColor = Ia*ka + Id*kd* dotNL + Is*ks* pow (dotRV , n);
	ambientDiffuse = Ia*ka + Id*kd* dotNL;

	if (dotNL == 0.0) {
		isOnDarkSide = 0.0;
	}
	else {
		isOnDarkSide = 1.0;
	}
	cameraPos = cameraPosition;

	texCoords = TextCoordinates;

	worldPosition = vec3(MV * vec4(Position, 1.0));

	// Transform (x,y,z) vertex coordinates with a 4x4 matrix T
    gl_Position = T * vec4(Position, 1.0);
}