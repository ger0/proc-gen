#version 330
layout (location = 0)in vec3 vertex;
layout (location = 1)in vec3 normal;
layout (location = 2)in vec4 color;
//layout (location = 2)in unsigned blockType;

uniform mat4 P, V, M;

out vec3 iNormal;
out vec3 vertPos;
out vec4 iColor;
out float visibility;

const float density  = 0.005;
const float gradient = 7.0;

void main() {
	gl_Position = P * V * M * vec4(vertex, 1.0);
	vec4 vertPos4 = V * M * vec4(vertex, 1.0);
	vertPos = vec3(vertPos4) / vertPos4.w;
	iNormal = vec3(V * M * vec4(normal, 0.0));

	vec4 tempcol = color;

	vec3 surfNormal = vec3(M * vec4(normal, 0.0));
	//lightPos = vec3(0, 0, 0);
	//calculating steepness
	vec3 up = vec3(0.0,1.0,0.0);
	float dong = dot(surfNormal, up);
	if (dong < 0.2) {
		tempcol = vec4(0.07, 0.07, 0.07,  1.0);
	}
	else if (dong < 0.4) {
		tempcol = vec4(0.07, 0.04, 0.005, 1.0);
	}
	iColor = tempcol;

	float distance = length(vertPos4.xyz);
	visibility = exp(-pow((distance * density), gradient));
	visibility = clamp(visibility, 0.0, 1.0);
}
