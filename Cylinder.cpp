#include "Cylinder.hpp"
#include <glm/gtc/type_ptr.hpp>
#include "shaderprogram.h"


Cylinder::Cylinder(float sectorCount, float r, float h){
  generate_indicis_of_cylinder(sectorCount, r, h);
}

// glm::vec4 calculate_normals(glm::vec4 x, glm::vec4 y, glm::vec4 z)
// {
//   glm::vec4 u = y - x;
//   glm::vec4 v = z - x;

//   glm::vec4 normal;

//   normal.x = (u.y * v.z) - (u.z * v.y);
//   normal.y = (u.z * v.x) - (u.x * v.z);
//   normal.z = (u.x * v.y) - (u.y * v.x);
// }

void Cylinder::generate_indicis_of_cylinder(float sectorCount, float r, float h){

  const float PI = 3.1415926f;
  float sectorStep = 2 * PI / sectorCount;
  float sectorAngle;  // radian
  float x,y;

  // top circle
  verts.push_back(glm::vec4(0,0,0,1)); //center
  norms.push_back(glm::vec4(0,0,0,1));

  //norms 
  for(int i = 0; i <= sectorCount; ++i)
  {
      sectorAngle = i * sectorStep;
      x = cos(sectorAngle) * r;
      y = sin(sectorAngle) * r;
      verts.push_back(glm::vec4(x, 0, y, 1));
      norms.push_back(glm::vec4(cos(sectorAngle), x, sin(sectorAngle),1));
  }  
  // bottom circle
  verts.push_back(glm::vec4(0,0,y,1)); //center
  norms.push_back(glm::vec4(0,0,sin(sectorAngle),1));

  for(int i = 0; i <= sectorCount; ++i)
  {
      sectorAngle = i * sectorStep;
      x = cos(sectorAngle) * r;
      y = sin(sectorAngle) * r;
      verts.push_back(glm::vec4(x, h, y, 1));
      norms.push_back(glm::vec4(cos(sectorAngle), h, sin(sectorAngle),1));

  }  



  // top circle
  for (int i = 0; i < sectorCount+1; i++)
  {
    indices.push_back(0);
    indices.push_back(i);
    indices.push_back(i+1);
  }
  
  // walls
  for (int i = 1; i < sectorCount+1; i++)
  {
    indices.push_back(i);
    indices.push_back(i+1);
    indices.push_back(i+sectorCount+2);

    indices.push_back(i);
    indices.push_back(i+sectorCount+1);
    indices.push_back(i+sectorCount+2);
  }

  // bottom circle
  for (int i = sectorCount + 2; i < verts.size()-1; i++)
  {
    indices.push_back(sectorCount);
    indices.push_back(i);
    indices.push_back(i+1);
  }
}




void Cylinder::draw(GLFWwindow* window, glm::mat4 V, glm::mat4 P, glm::mat4 M){
  sp_cylinder->use();  //activate shading program

	//Send parameters to graphics card
  glUniformMatrix4fv(sp_cylinder->u("P"),1,false,glm::value_ptr(P));
  glUniformMatrix4fv(sp_cylinder->u("V"),1,false,glm::value_ptr(V));
  glUniformMatrix4fv(sp_cylinder->u("M"),1,false,glm::value_ptr(M));

	glEnableVertexAttribArray(sp_cylinder->a("vertex")); //Enable sending data to the attribute vertex
  glVertexAttribPointer(sp_cylinder->a("vertex"),4,GL_FLOAT,false,0, verts.data()); //Specify source of the data for the attribute vertex

	glEnableVertexAttribArray(sp_cylinder->a("normal"));
  glVertexAttribPointer(sp_cylinder->a("normal"),4,GL_FLOAT,false,0, norms.data());

  glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, indices.data()); //Draw the object

  glDisableVertexAttribArray(sp_cylinder->a("vertex")); //Disable sending data to the attribute vertex
	glDisableVertexAttribArray(sp_cylinder->a("normal"));
}
