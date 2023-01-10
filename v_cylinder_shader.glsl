#version 330
in vec4 vertex;
in vec4 normal;

uniform mat4 P, V, M;

varying vec3 iNormal;
varying vec3 vertPos;

void main() {
    gl_Position = P * V * M * vertex;
    vec4 vertPos4 = V * vertex;
    vertPos = vec3(vertPos4) / vertPos4.w;
    iNormal = vec3(M * normal);
}