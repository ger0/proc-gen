#version 330

/* fragment shader for objects with many textures - phong lightening*/

//Uniform variables
uniform mat4 P;
uniform mat4 V;
uniform mat4 M;
uniform vec4 light_position1 = vec4(20,50,-35,0); //in world space
uniform vec4 light_position2 = vec4(-20,50,35,0);


//Attributes
in vec4 vertex; //Vertex coordinates in model space
in vec4 normal; //in model space
in vec2 texCoord; //texturing coordinates
in vec4 tangent; // in model space


//Zmienne interpolowane
out vec4 l[2];    // znormalizowany wektor do źródła światła w przestrzeni oka
out vec4 v;
out vec2 i_texc;  // współrzędne teksturowania


void main(void) {
    /* liczymy kolory*/

    /*lambert*/
    // interpolacted_color = kd*ld*nl
    // kd - kolor powierzchni

    /*phong*/
    //L = ka*la + kd*ld*nl + ks*ls*pow(rv, alfa)
    //ks - kolor materiału dla światła odbitego
    //ld - kolor światła odbitego
    //rv – cos kąta między światłem odbitym a wektorem od powierzchni do obserwatora

    //normal mapping
    // wszystko wyrażamy do przestrzeni styczniej
    //invTBM -> z modelu do stycznej

    vec4 norm_n = normalize(normal);
    vec4 norm_bitangent = normalize(vec4(cross(norm_n.xyz, tangent.xyz),0));
    vec4 norm_tangent = normalize(tangent);
    mat4 TBN = mat4(norm_tangent, norm_bitangent, norm_n, vec4(0,0,0,1));
    mat4 invTBN = inverse(TBN);

    l[0] = normalize(invTBN*inverse(M)*light_position1 - invTBN*vertex);
    l[1] = normalize(invTBN*inverse(M)*light_position2 - invTBN*vertex);
    v = normalize(invTBN*inverse(V*M)*vec4(0,0,0,1) - invTBN*vertex); //od powierzchni do obserwatora

    i_texc=texCoord;

    gl_Position=P*V*M*vertex;
}
