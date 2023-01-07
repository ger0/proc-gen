#version 330
precision mediump float;

vec3     lightPos = vec3(0, 0, -1);
in vec3  iNormal;
in vec3  vertPos;
in vec4  iColor;
in vec3  surfNormal;
in float visibility;

uniform int mode;
uniform vec3 skyColor;

const vec3 lightColor = vec3(0.5, 0.5, 0.4);
const float lightPower = 0.25;
vec3 ambientColor = vec3(0.01, 0.03, 0.08);
const vec3 diffuseColor = vec3(0.33, 0.33, 0.0);
const vec3 specColor = vec3(1.0, 1.0, 1.0);
const float shininess = 2.0;
const float screenGamma = 2.2; // Assume the monitor is calibrated to the sRGB color space

void main() {
	ambientColor = iColor.rgb;
	vec3 normal = normalize(iNormal);

	vec3 up = vec3(0,1,0);
	float classify = dot(surfNormal, up);
	if (classify < 0.2) ambientColor = vec3(0.07,  0.07, 0.07);
	else if (classify < 0.4) ambientColor = vec3(0.07,  0.04, 0.005);

	vec3 lightDir = lightPos - vertPos;
	float distance = length(lightDir);
	distance = distance * distance;
	lightDir = normalize(lightDir);

	float lambertian = max(dot(lightDir, normal), 0.0);
	float specular = 0.0;

	if (lambertian > 0.0) {
		vec3 viewDir = normalize(-vertPos);

		// this is blinn phong
		vec3 halfDir = normalize(lightDir + viewDir);
		float specAngle = max(dot(halfDir, normal), 0.0);
		specular = pow(specAngle, shininess);

		// this is phong (for comparison)
		if (mode == 2) {
    		vec3 reflectDir = reflect(-lightDir, normal);
    		specAngle = max(dot(reflectDir, viewDir), 0.0);
    		// note that the exponent is different here
    		specular = pow(specAngle, shininess/4.0);
		}
	}
	vec3 colorLinear = ambientColor +
		diffuseColor * lambertian * lightColor * lightPower +
		specColor * specular * lightColor * lightPower;
	// apply gamma correction (assume ambientColor, diffuseColor and specColor
	// have been linearized, i.e. have no gamma correction in them)
	vec3 colorGammaCorrected = pow(colorLinear, vec3(1.0 / screenGamma));
	// use the gamma corrected color in the fragment
	// mix with fog
	vec3 outColor = mix(skyColor, colorGammaCorrected, visibility);
	gl_FragColor = vec4(outColor, iColor.a);
}
