#version 330
layout (location = 0)in vec3 vertex;
layout (location = 1)in vec3 normal;

uniform mat4 P, V, M;

out vec3 iNormal;
out vec3 vertPos;
//out vec3 lightPos;

void main() {
    gl_Position = P * V * M * vec4(vertex, 1.0);
    vec4 vertPos4 = V * M * vec4(vertex, 1.0);
    vertPos = vec3(vertPos4) / vertPos4.w;
    iNormal = vec3(V * M * vec4(normal, 0.0));
    //lightPos = vec3(0, 0, 0);
}
