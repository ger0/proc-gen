#include <glm/gtc/type_ptr.hpp>
#include "Sphere.hpp"
#include "shaderprogram.h"
#include <iostream>


Sphere::Sphere(float sectorCount,float stackCount, float radius, float h, glm::vec3 ac, glm::vec3 dc){

  ambientColor = ac;
  diffuseColor = dc;

  generate_indicis_of_sphere(sectorCount, stackCount, radius, h);
}


void Sphere::generate_indicis_of_sphere(float sectorCount,float stackCount, float radius, float h){

// std::vector<float>().swap(verts);
// std::vector<float>().swap(norms);
// std::vector<float>().swap(texCoords);

const float PI = 3.1415926f;

float x, y, z, xy;                              // vertex position
float nx, ny, nz, lengthInv = 1.0f / radius;    // vertex normal
float s, t;                                     // vertex texCoord

float sectorStep = 2 * PI / sectorCount;
float stackStep = PI / stackCount;
float sectorAngle, stackAngle;


for(int i = 0; i <= stackCount; ++i)
{
    stackAngle = PI / 2 - i * stackStep;        // starting from pi/2 to -pi/2
    xy = radius * cosf(stackAngle);             // r * cos(u)
    z = radius * sinf(stackAngle);              // r * sin(u)

    // add (sectorCount+1) vertices per stack
    // the first and last vertices have same position and normal, but different tex coords
    for(int j = 0; j <= sectorCount; ++j)
    {
        sectorAngle = j * sectorStep;           // starting from 0 to 2pi

        // vertex position (x, y, z)
        x = xy * cosf(sectorAngle);             // r * cos(u) * cos(v)
        y = xy * sinf(sectorAngle);             // r * cos(u) * sin(v)
        verts.push_back(glm::vec4(x,y,z,1));
        // verts.push_back(y);
        // verts.push_back(z);

        //normalized vertex normal (nx, ny, nz)
        nx = x * lengthInv;
        ny = y * lengthInv;
        nz = z * lengthInv;
        norms.push_back(glm::vec4(nx,ny,nz,1));


        // // vertex tex coord (s, t) range between [0, 1]
        s = (float)j / sectorCount;
        t = (float)i / stackCount;
        texCoords.push_back(glm::vec2(s, t));
    }


int k1, k2;
for(int i = 0; i < stackCount; ++i)
{
    k1 = i * (sectorCount + 1);     // beginning of current stack
    k2 = k1 + sectorCount + 1;      // beginning of next stack

    for(int j = 0; j < sectorCount; ++j, ++k1, ++k2)
    {
        // 2 triangles per sector excluding first and last stacks
        // k1 => k2 => k1+1
        if(i != 0)
        {
            indices.push_back(k1);
            indices.push_back(k2);
            indices.push_back(k1 + 1);
        }

        // k1+1 => k2 => k2+1
        if(i != (stackCount-1))
        {
            indices.push_back(k1 + 1);
            indices.push_back(k2);
            indices.push_back(k2 + 1);
        }

        // store indices for lines
        // vertical lines for all stacks, k1 => k2
        lineIndices.push_back(k1);
        lineIndices.push_back(k2);
        if(i != 0)  // horizontal lines except 1st stack, k1 => k+1
        {
            lineIndices.push_back(k1);
            lineIndices.push_back(k1 + 1);
        }
    }
}
}
}




void Sphere::draw(GLFWwindow* window, glm::mat4 V, glm::mat4 P, glm::mat4 M){
  sp_sphere->use();  //activate shading program

	//Send parameters to graphics card
  glUniformMatrix4fv(sp_sphere->u("P"),1,false,glm::value_ptr(P));
  glUniformMatrix4fv(sp_sphere->u("V"),1,false,glm::value_ptr(V));
  glUniformMatrix4fv(sp_sphere->u("M"),1,false,glm::value_ptr(M));

  glUniform3fv(sp_sphere->u("ambientColor"),1,glm::value_ptr(ambientColor));
  glUniform3fv(sp_sphere->u("diffuseColor"),1,glm::value_ptr(diffuseColor));



	glEnableVertexAttribArray(sp_sphere->a("vertex")); //Enable sending data to the attribute vertex
  glVertexAttribPointer(sp_sphere->a("vertex"),4,GL_FLOAT,false,0, verts.data()); //Specify source of the data for the attribute vertex

	glEnableVertexAttribArray(sp_sphere->a("normal"));
  glVertexAttribPointer(sp_sphere->a("normal"),4,GL_FLOAT,false,0, norms.data());


  glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, indices.data()); //Draw the object

  glDisableVertexAttribArray(sp_sphere->a("vertex")); //Disable sending data to the attribute vertex
	glDisableVertexAttribArray(sp_sphere->a("normal"));
}
