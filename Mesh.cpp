#include "Mesh.hpp"
#include <glm/gtc/type_ptr.hpp>
#include "shaderprogram.h"


Mesh::Mesh(){
  generate_indicis_of_cylinder(20.0f, 0.5f, 3.0f);
}

void Mesh::generate_indicis_of_cylinder(float sectorCount, float r, float h){

  const float PI = 3.1415926f;
  float sectorStep = 2 * PI / sectorCount;
  float sectorAngle;  // radian
  float x,y;

  // top circle
  verts.push_back(glm::vec4(0,0,0,1)); //center
  for(int i = 0; i <= sectorCount; ++i)
  {
      sectorAngle = i * sectorStep;
      x = cos(sectorAngle) * r;
      y = sin(sectorAngle) * r;
      verts.push_back(glm::vec4(x, y, 0, 1));
  }  
  // bottom circle
  verts.push_back(glm::vec4(0,0,h,1)); //center
  for(int i = 0; i <= sectorCount; ++i)
  {
      sectorAngle = i * sectorStep;
      x = cos(sectorAngle) * r;
      y = sin(sectorAngle) * r;
      verts.push_back(glm::vec4(x, y, h, 1));
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




void Mesh::draw(GLFWwindow* window, glm::mat4 V, glm::mat4 P, glm::mat4 M, std::vector<GLuint> textures){
  sp->use();  //activate shading program

	//Send parameters to graphics card
  glUniformMatrix4fv(sp->u("P"),1,false,glm::value_ptr(P));
  glUniformMatrix4fv(sp->u("V"),1,false,glm::value_ptr(V));
  glUniformMatrix4fv(sp->u("M"),1,false,glm::value_ptr(M));

	glEnableVertexAttribArray(sp->a("vertex")); //Enable sending data to the attribute vertex
  glVertexAttribPointer(sp->a("vertex"),4,GL_FLOAT,false,0, verts.data()); //Specify source of the data for the attribute vertex

	// glEnableVertexAttribArray(sp->a("texCoord"));
  // glVertexAttribPointer(sp->a("texCoord"),2,GL_FLOAT,false,0, texCoords.data());

	// glEnableVertexAttribArray(sp->a("normal"));
  // glVertexAttribPointer(sp->a("normal"),4,GL_FLOAT,false,0, norms.data());

  // glEnableVertexAttribArray(sp->a("tangent"));
  // glVertexAttribPointer(sp->a("tangent"),4,GL_FLOAT,false,0, tangents.data());

  // glUniform1i(sp->u("texMapColor"), 0); // powiązanie zmiennej z jednostką teksturującą
	// glActiveTexture(GL_TEXTURE0);
	// glBindTexture(GL_TEXTURE_2D, textures[0]);

  // glUniform1i(sp->u("texMapReflect"), 1);
  // glActiveTexture(GL_TEXTURE1);
	// glBindTexture(GL_TEXTURE_2D, textures[1]);

  // glUniform1i(sp->u("texMapNormal"), 2);
  // glActiveTexture(GL_TEXTURE2);
	// glBindTexture(GL_TEXTURE_2D, textures[2]);

  glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, indices.data()); //Draw the object

  glDisableVertexAttribArray(sp->a("vertex")); //Disable sending data to the attribute vertex
	// glDisableVertexAttribArray(sp->a("texCoord"));
	// glDisableVertexAttribArray(sp->a("normal"));
  // glDisableVertexAttribArray(sp->a("tangent"));
  // glDisableVertexAttribArray(sp->a("texMapColor"));
  // glDisableVertexAttribArray(sp->a("texMapReflect"));
  // glDisableVertexAttribArray(sp->a("texMapNormal"));
}
/////////////////////////////////////////////////////////////////////
