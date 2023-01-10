#version 330

/* vertex shader for objecs with many meshes with different textures – Lambert lighting */

//Uniform variables
uniform mat4 P;
uniform mat4 V;
uniform mat4 M;
uniform vec4 light_position1 = vec4(20,50,-35,0); //in world space
uniform vec4 light_position2 = vec4(-20,50,35,0);


//Attributes
in vec4 vertex; //Vertex coordinates in model space
in vec4 color;
in vec4 normal; //in model space
in vec2 texCoord; //texturing coordinates


//Zmienne interpolowane
out vec4 i_color; // interpolowany kolor
out vec4 n;     // wektor normalny powierzchi w przestrzeni oka
out vec4 l[2];    // znormalizowany wektor do źródła światła w przestrzeni oka
out vec4 v;
out vec2 i_texc;  // współrzędne teksturowania


void main(void) {
    l[0] = normalize(V*light_position1 - V*M*vertex);
    l[1] = normalize(V*light_position2 - V*M*vertex);
    n = normalize(V*M*normal);
    v = normalize(vec4(0,0,0,1) - V*M*vertex); //od powierzchni do obserwatora
    i_color = color;
    i_texc=texCoord;

    gl_Position=P*V*M*vertex;
}
