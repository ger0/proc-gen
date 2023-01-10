#include <glm/gtc/type_ptr.hpp>
#include "Sphere.hpp"
#include <iostream>


Sphere::Sphere(float sectorCount,float stackCount, float radius, float h, glm::vec3 ac, glm::vec3 dc){

  	ambientColor = ac;
  	diffuseColor = dc;

  	generate_indicis_of_sphere(sectorCount, stackCount, radius, h);
}

void Sphere::generate_indicis_of_sphere(float sectorCount,float stackCount, float radius, float h){
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

        	//normalized vertex normal (nx, ny, nz)
        	nx = x * lengthInv;
        	ny = y * lengthInv;
        	nz = z * lengthInv;
        	norms.push_back(glm::vec4(nx,ny,nz,1));
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

Mesh Sphere::genMesh(glm::mat4 M) {
	Mesh mesh;
	for (int i = 0; i < indices.size(); i += 3) {
		for (int j = i + 2; j >= i; j--) {
			auto &idx = indices[j];
			mesh.push_back(Vertex{
					.pos 	= glm::vec3(M * glm::vec4(verts[idx], 1.0)),
					.norm 	= glm::vec3(M * glm::vec4(norms[idx], 0.0)),
					.color 	= glm::vec4(ambientColor, 1.0)
					});	
		}
	}
	return mesh;
}
