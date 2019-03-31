#version 330 core

out vec4 FragColor;

in vec2 TexCoords;
in vec3 Normal;
in vec3 FragPosition;
in vec4 FragPosLightSpace;

uniform sampler2D textureImage;
uniform sampler2D shadowMap;

uniform vec3 objectColor;
uniform vec3 viewPosition;
uniform vec3 lightPosition;
uniform vec3 lightColor;


float ShadowCalculation(vec4 fragPosLightSpace, vec3 lightDirection)
{
    // Perspective division
    vec3 projCoords = fragPosLightSpace.xyz / fragPosLightSpace.w;

    // Set to [0,1] range
    projCoords = projCoords * 0.5 + 0.5;

    // Current fragment from the light's perspective
    float currentDepth = projCoords.z;

	// Shadow bias based on the surface angle towards the light (bias removes shadow acne)
    float bias = max(0.0001 * (1.0 - dot(Normal, lightDirection)), 0.00001);

	// Size of a single texel used to offset the texture coordinates
	vec2 texelSize = 1.0 / textureSize(shadowMap, 0);
	float shadow = 0.0;
	// Apply percentage-closer filtering
	for(int x = -1; x <= 1; ++x) {
		for(int y = -1; y <= 1; ++y) {
			float pcfDepth = texture(shadowMap, projCoords.xy + vec2(x, y) * texelSize).r; 
			shadow += currentDepth - bias > pcfDepth ? 1.0 : 0.0;        
		}    
	}
	float sampleSize = 12;
	shadow /= sampleSize;

	// Remove shadows casted outside of the light's frustrum
	if(projCoords.z > 1.0)
        shadow = 0.0;

    return shadow;
} 

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
	float intensity = min(80.0 / dist, 1.0);

	// calculate shadow
    float shadow = ShadowCalculation(FragPosLightSpace, lightDirection);

	// Phong lighting model
	vec3 phong = (ambient + (1.0 - shadow) * (diffuse + specular)) * objectColor * intensity;

	// Final shaded color (texture * lighting)
    FragColor = texture(textureImage, TexCoords) * vec4 (phong, 1.0);
} 