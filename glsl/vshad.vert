#version 330
layout (location = 0)in vec3 vertex;
layout (location = 1)in vec3 normal;

uniform mat4 P, V, M;

varying vec3 iNormal;
varying vec3 vertPos;

void main() {
    gl_Position = P * V * vec4(vertex, 1.0);
    vec4 vertPos4 = V * vec4(vertex, 1.0);
    vertPos = vec3(vertPos4) / vertPos4.w;
    iNormal = vec3(M * vec4(normal, 0.0));
}
