#version 330 core

out vec4 FragColor;

in vec2 TexCoords;
in vec3 Normal;
in vec3 FragPosition;

uniform sampler2D textureImage;
uniform vec3 objectColor;
uniform vec3 viewPosition;
uniform vec3 lightPosition;
uniform vec3 lightColor;

void main()
{
	vec3 viewDirection = normalize(viewPosition - FragPosition);
	vec3 lightDirection = normalize(lightPosition - FragPosition);
	vec3 reflectDirection = reflect(-lightDirection, Normal); 

	// Ambient lighting
	float ambientStrength = 0.1;
    vec3 ambient = ambientStrength * lightColor;

	// Diffuse lighting
	float diff = max(dot(Normal, lightDirection), 0.0);
	vec3 diffuse = (diff * lightColor);

	// Specular lighting
	float specularStrength = 0.5;
	float shininess = 32;
	float spec = pow(max(dot(viewDirection, reflectDirection), 0.0), shininess);
	vec3 specular = specularStrength * spec * lightColor;  

	// Calculate the intensity decrement
	float dist = distance(lightPosition, FragPosition);
	float intensity = min(100.0 / dist, 1.0);

	// Phong lighting model
	vec3 phong = (ambient + diffuse + specular) * objectColor * intensity;

	// Final shaded color (texture * lighting)
    FragColor = texture(textureImage, TexCoords) * vec4 (phong, 1.0);
} 