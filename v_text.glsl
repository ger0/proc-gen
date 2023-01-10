#version 330

/* vertex shader to display the 2D text */

//Attributes
in vec2 vertex; //Vertex coordinates in model space
in vec2 texCoord; //texturing coordinates

//Zmienne interpolowane
out vec2 uv;

void main(){
    vec2 vertexPosition_screenspace = vertex;
    // Output position of the vertex, in clip space
    // map [0..800][0..600] to [-1..1][-1..1]
    vec2 vertexPosition_homoneneousspace = vertexPosition_screenspace - vec2(400,300); // [0..800][0..600] -> [-400..400][-300..300]
    vertexPosition_homoneneousspace /= vec2(400,300);
    gl_Position =  vec4(vertexPosition_homoneneousspace,0,1);

    // UV of the vertex. No special space for this one.
    uv = texCoord;
}
