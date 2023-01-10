// #version 330

// in vec4 vertex; //Vertex coordinates in model space
// in vec4 normal; //in model space

// uniform mat4 M;
// uniform mat4 V;
// uniform mat4 P;

// out vec4 fragNormal;


// void main()
// {

//     fragNormal = normal;
//     gl_Position = P * V * M * vertex;
// }


#version 330
in vec4 vertex;
in vec4 normal;
// in vec3 color;


uniform mat4 P, V, M;
uniform vec3 ambientColor;
uniform vec3 diffuseColor;


varying vec3 iNormal;
varying vec3 iambientColor;
varying vec3 idiffuseColor;
varying vec3 vertPos;


void main() {
    gl_Position = P * V * M * vertex;
    vec4 vertPos4 = V * vertex;
    vertPos = vec3(vertPos4) / vertPos4.w;
    iNormal = vec3(M * normal);
    iambientColor = ambientColor;
    idiffuseColor = diffuseColor;
}