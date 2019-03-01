#version 330 core

out vec4 FragColor;

layout(pixel_center_integer) in vec4 gl_FragCoord;

//in vec3 shadedColor;
in vec2 texCoords;

// Light stuff
in vec3 cameraPos;
in vec3 reflection;
in vec3 ambientDiffuse;
in float isOnDarkSide;
in vec3 worldPosition;

uniform sampler2D texture;

float n = 15.0; // the " shininess " parameter
vec3 ks = vec3(1.0, 1.0, 1.0); // the specular surface reflection color
vec3 Is = vec3(1.0, 1.0, 1.0); // the specular illumination color

void main()
{

	vec3 viewDirection = normalize(worldPosition - cameraPos);

	float dotRV = max(dot(reflection, viewDirection), 0.0) * isOnDarkSide;
	//if (isOnDarkSide) dotRV = 0.0;

	vec3 shadedColor = ambientDiffuse + Is*ks* pow (dotRV , n);

    FragColor = texture(texture, texCoords) * vec4 (shadedColor, 1.0);
} 