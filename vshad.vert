#version 330
uniform mat4 P;
uniform mat4 V;
uniform mat4 M;

in vec3 vertex;

out   vec3 iColor;

void main() {
   gl_Position = P * V * M * vec4(vertex, 1.f);
   iColor = vec3(1.f, 1.f, 1.f);
}
