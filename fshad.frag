#version 330

in vec3 iColor;

out vec4 pixelColor;

void main(void) {
    pixelColor = vec4(iColor, 1.f);
}
