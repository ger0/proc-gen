#version 330
layout (location = 0)in vec3 vertex;
layout (location = 1)in vec3 normal;
layout (location = 2)in vec4 color;
//layout (location = 2)in unsigned blockType;

uniform mat4 P, V, M;

out vec3 iNormal;
out vec3 surfNormal;
out vec3 vertPos;
out vec4 iColor;
out float visibility;

const float density  = 0.015;
const float gradient = 5.0;

void main() {
	gl_Position = P * V * M * vec4(vertex, 1.0);
	vec4 vertPos4 = V * M * vec4(vertex, 1.0);
	vertPos = vec3(vertPos4) / vertPos4.w;
	iNormal = vec3(V * M * vec4(normal, 0.0));
	surfNormal = vec3(M * vec4(normal, 0.0));
	//lightPos = vec3(0, 0, 0);
	iColor = color;

	float distance = length(vertPos4.xyz);
	visibility = exp(-pow((distance * density), gradient));
	visibility = clamp(visibility, 0.0, 1.0);
}
