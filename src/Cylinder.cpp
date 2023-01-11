#include "Cylinder.hpp"
#include <glm/gtc/type_ptr.hpp>
#include "main.hpp"

Cylinder::Cylinder(float sectorCount, float R, float r, float h){
  	generate_indicis_of_cylinder(sectorCount, R, r, h);
}

void Cylinder::generate_indicis_of_cylinder(float sectorCount, float R, float r, float h) {
  	const float PI = 3.1415926f;
  	float sectorStep = 2 * PI / sectorCount;
  	float sectorAngle;  // radian
  	float x,y;

  	// top circle
  	verts.push_back(glm::vec3(0,0,0)); //center
  	norms.push_back(glm::vec3(0,0,0));

  	//norms 
  	for(int i = 0; i <= sectorCount; ++i)
  	{
      	sectorAngle = i * sectorStep;
      	x = cos(sectorAngle) * R;
      	y = sin(sectorAngle) * R;
      	verts.push_back(glm::vec3(x, 0, y));
      	norms.push_back(glm::vec3(cos(sectorAngle), x, sin(sectorAngle)));
  	}  
  	// bottom circle
  	verts.push_back(glm::vec3(0,0,y)); //center
  	norms.push_back(glm::vec3(0,0,sin(sectorAngle)));

  	for(int i = 0; i <= sectorCount; ++i)
  	{
      	sectorAngle = i * sectorStep;
      	x = cos(sectorAngle) * r;
      	y = sin(sectorAngle) * r;
      	verts.push_back(glm::vec3(x, h, y));
      	norms.push_back(glm::vec3(cos(sectorAngle), h, sin(sectorAngle)));

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

Mesh Cylinder::genMesh(glm::mat4 M, glm::vec4 color) {
	Mesh mesh;
	for (int i = 0; i < indices.size(); i += 3) {
		for (int j = i + 2; j >= i; j--) {
			auto &idx = indices[j];
			mesh.push_back(Vertex{
					.pos 	= glm::vec3(M * glm::vec4(verts[idx], 1.0)),
					.norm 	= glm::vec3(M * glm::vec4(norms[idx], 0.0)),
					.color 	= color
					});	
		}
	}
	return mesh;
}
