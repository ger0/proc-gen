#include "model.hpp"

#include <unistd.h>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <glm/gtx/string_cast.hpp>

constexpr float CYLINDER_R = 0.30;
constexpr float CYLINDER_H = 2.5; 

int DEPTH = 3;
std::vector<int> operations;
std::vector<std::vector<float>> sizemod;
std::vector<std::vector<float>> rotations;

Model::Model(glm::vec3 _pos, float _angle_x, float _angle_y){
  	branche = Cylinder(4.0f, CYLINDER_R, CYLINDER_H); //float sectorCount, float r, float h
  	green_sphere = Sphere(5, 5, 4.8, 4.8, glm::vec3(0.0, 0.07, 0.01), glm::vec3(0.5, 0.5, 0.3)); //float sectorCount,float stackCount, float radius, float h
  	brown_sphere = Sphere(5, 5 , 0.4, 0.4, glm::vec3(0.06, 0.02, 0.0), glm::vec3(0.5, 0.5, 0.3)); //float sectorCount,float stackCount, float radius, float h

  	srand (static_cast <unsigned> (time(0)));
  	for (int i=0; i<DEPTH; i++){
    	operations.push_back( rand()%4);
  	}

  	for (int j=0; j<DEPTH + 1; j++){
    	std::vector<float> inner_sizemod;
    	std::vector<float> inner_rotations;
    	for (int i=0; i<10; i++){
      		inner_sizemod.push_back(0.75 + static_cast <float> (rand()) 
      				/( static_cast <float> (RAND_MAX/(1.15-0.75))));
      		inner_rotations.push_back(-20 + static_cast <float> (rand()) 
      				/( static_cast <float> (RAND_MAX/float(20+20))));
    	}
    	sizemod.push_back(inner_sizemod);
    	rotations.push_back(inner_rotations);
  	}

}

void meshPush(Mesh &mesh, Mesh &val) {
    printf("val size: %lu\n", val.size());
    mesh.insert(mesh.end(), val.begin(), val.end());
    printf("size after pushing: %lu\n", mesh.size());
}

Mesh Model::trunk(float size, int depth, glm::vec3 position, glm::mat4 V, glm::mat4 M){
    int operation = operations[depth];
    Mesh mesh;

    if (operation < 1) {
        auto b1 = branch_one(size , depth, position, V, M);          
        meshPush(mesh, b1);
    } 
    if (operation < 2) {
        auto b2 = branch_two(size, depth, position, V, M);          
        meshPush(mesh, b2);
    } else if (operation < 3) {
        auto b3 = branch_three(size, depth, position, V, M);          
        meshPush(mesh, b3);
    } else if (operation < 4) {
        auto b4 = branch_four(size, depth, position, V, M);          
        meshPush(mesh, b4);
    }
    return mesh;
}

Mesh Model::branch_one(float size, int depth, glm::vec3 position, glm::mat4 V, glm::mat4 M){
    Mesh mesh;
  	glm::mat4 P=glm::perspective(glm::radians(50.0f), 1.0f, 1.0f, 120.0f); //compute projection matrix

  	M = glm::translate(M, position);
  	M = glm::scale(M, glm::vec3(size, size, size));
  	auto b0 = branche.genMesh(M, glm::vec4(0.06, 0.02, 0.0, 1.0));
  	meshPush(mesh, b0);

  	M = glm::translate(M, glm::vec3(0, CYLINDER_H, 0));

	if (depth > 0) {   
        //brown_sphere.draw(V, P, M);

        glm::mat4 M_1_rotate;
        M_1_rotate=glm::rotate(M,(float)glm::radians(0.0 + rotations[depth][0]),glm::vec3(1.0f,0.0f,0.0f));
        M_1_rotate=glm::rotate(M_1_rotate,(float)glm::radians(0.0 + rotations[depth][2]),glm::vec3(0.0f,1.0f,0.0f));
        M_1_rotate=glm::rotate(M_1_rotate,(float)glm::radians(0.0 + rotations[depth][1]),glm::vec3(0.0f,0.0f,1.0f));

        auto trnk = trunk(size * 0.9 * sizemod[depth][0], depth-1, position, V, M_1_rotate);
        meshPush(mesh, trnk);
	} else {
      	auto grn_sphre = green_sphere.genMesh(M);
      	meshPush(mesh, grn_sphre);
	}
	return mesh;
}

Mesh Model::branch_two(float size, int depth, glm::vec3 position, glm::mat4 V, glm::mat4 M){
    Mesh mesh;
  	glm::mat4 P=glm::perspective(glm::radians(50.0f), 1.0f, 1.0f, 120.0f); //compute projection matrix

  	M=glm::translate(M, position);
  	M = glm::scale(M, glm::vec3(size, size, size));
  	auto b0 = branche.genMesh(M, glm::vec4(0.06, 0.02, 0.0, 1.0));
  	meshPush(mesh, b0);

  	M = glm::translate(M, glm::vec3(0, CYLINDER_H, 0));

	if (depth > 0) 
	{
        //brown_sphere.draw(V, P, M);

        glm::mat4 M_1_rotate;
        M_1_rotate=glm::rotate(M,(float)glm::radians(0.0 + rotations[depth][0]),glm::vec3(1.0f,0.0f,0.0f));
        M_1_rotate=glm::rotate(M_1_rotate,(float)glm::radians(0.0 + rotations[depth][1]),glm::vec3(0.0f,1.0f,0.0f));
        M_1_rotate=glm::rotate(M_1_rotate,(float)glm::radians(30.0 + rotations[depth][2]),glm::vec3(0.0f,0.0f,1.0f));
        auto t1 = trunk(size*0.9*sizemod[depth][0], depth-1, position, V, M_1_rotate);
        meshPush(mesh, t1);

        glm::mat4 M_2_rotate;
        M_2_rotate=glm::rotate(M,(float)glm::radians(0.0 + rotations[depth][3]),glm::vec3(1.0f,0.0f,0.0f));
        M_2_rotate=glm::rotate(M_2_rotate,(float)glm::radians(180.0 + rotations[depth][4]),glm::vec3(0.0f,1.0f,0.0f));
        M_2_rotate=glm::rotate(M_2_rotate,(float)glm::radians(30.0 + rotations[depth][5]),glm::vec3(0.0f,0.0f,1.0f));
        auto t2 = trunk(size*0.9*sizemod[depth][1], depth-1, position, V, M_2_rotate);
        meshPush(mesh, t2);
	}
	else
	{
      	//green_sphere.draw(V, P, M);
      	auto grn_sphre = green_sphere.genMesh(M);
      	meshPush(mesh, grn_sphre);
	}
	return mesh;
}


Mesh Model::branch_three(float size, int depth, glm::vec3 position, glm::mat4 V, glm::mat4 M){
    Mesh mesh;
  	glm::mat4 P=glm::perspective(glm::radians(50.0f), 1.0f, 1.0f, 120.0f); //compute projection matrix

  	M=glm::translate(M, position);
  	M = glm::scale(M, glm::vec3(size, size, size));
  	//branche.draw(V, P, M);
  	auto b0 = branche.genMesh(M, glm::vec4(0.06, 0.02, 0.0, 1.0));
  	meshPush(mesh, b0);

  	M = glm::translate(M, glm::vec3(0, CYLINDER_H, 0));

	if (depth > 0) 
	{

        //brown_sphere.draw(V, P, M);

        glm::mat4 M_1_rotate;
        M_1_rotate=glm::rotate(M,(float)glm::radians(0.0 + rotations[depth][0]),glm::vec3(1.0f,0.0f,0.0f));
        M_1_rotate=glm::rotate(M_1_rotate,(float)glm::radians(0.0 + rotations[depth][1]),glm::vec3(0.0f,1.0f,0.0f));
        M_1_rotate=glm::rotate(M_1_rotate,(float)glm::radians(30.0 + rotations[depth][2]),glm::vec3(0.0f,0.0f,1.0f));
        auto t1 = trunk(size*0.9*sizemod[depth][0], depth-1, position, V, M_1_rotate);
        meshPush(mesh, t1);

        glm::mat4 M_2_rotate;
        M_2_rotate=glm::rotate(M,(float)glm::radians(0.0f + rotations[depth][3]),glm::vec3(1.0f,0.0f,0.0f));
        M_2_rotate=glm::rotate(M_2_rotate,(float)glm::radians(120.0 + rotations[depth][4]),glm::vec3(0.0f,1.0f,0.0f));
        M_2_rotate=glm::rotate(M_2_rotate,(float)glm::radians(30.0 + rotations[depth][5]),glm::vec3(0.0f,0.0f,1.0f));
        auto t2 = trunk(size*0.9*sizemod[depth][1], depth-1, position, V, M_2_rotate);//sizemod[1]
        meshPush(mesh, t2);

        glm::mat4 M_3_rotate;
        M_3_rotate=glm::rotate(M,(float)glm::radians(0.0f + rotations[depth][6]),glm::vec3(1.0f,0.0f,0.0f));
        M_3_rotate=glm::rotate(M_3_rotate,(float)glm::radians(240.0f + rotations[depth][7]),glm::vec3(0.0f,1.0f,0.0f));
        M_3_rotate=glm::rotate(M_3_rotate,(float)glm::radians(30.0f + rotations[depth][8]),glm::vec3(0.0f,0.0f,1.0f));
        auto t3 = trunk(size*0.9*sizemod[depth][2], depth-1, position, V, M_3_rotate);
        meshPush(mesh, t3);
	}
	else
	{
      	//green_sphere.draw(V, P, M);
      	auto grn_sphre = green_sphere.genMesh(M);
      	meshPush(mesh, grn_sphre);
	}
	return mesh;
}
Mesh Model::branch_four(float size, int depth, glm::vec3 position, glm::mat4 V, glm::mat4 M){
    Mesh mesh;

  	M=glm::translate(M, position);
  	M = glm::scale(M, glm::vec3(size, size, size));
  	auto b0 = branche.genMesh(M, glm::vec4(0.06, 0.02, 0.0, 1.0));
  	meshPush(mesh, b0);

  	M = glm::translate(M, glm::vec3(0, CYLINDER_H, 0));

	if (depth > 0) {
        //brown_sphere.draw(V, P, M);

        glm::mat4 M_1_rotate;
        M_1_rotate=glm::rotate(M,(float)glm::radians(0.0 + rotations[depth][0]),glm::vec3(1.0f,0.0f,0.0f));
        M_1_rotate=glm::rotate(M_1_rotate,(float)glm::radians(0.0 + rotations[depth][1]),glm::vec3(0.0f,1.0f,0.0f));
        M_1_rotate=glm::rotate(M_1_rotate,(float)glm::radians(30.0 + rotations[depth][2]),glm::vec3(0.0f,0.0f,1.0f));
        auto t1 = trunk(size*0.9*sizemod[depth][0], depth-1, position, V, M_1_rotate);
        meshPush(mesh, t1);

        glm::mat4 M_2_rotate;
        M_2_rotate=glm::rotate(M,(float)glm::radians(0.0f + rotations[depth][3]),glm::vec3(1.0f,0.0f,0.0f));
        M_2_rotate=glm::rotate(M_2_rotate,(float)glm::radians(90.0 + rotations[depth][4]),glm::vec3(0.0f,1.0f,0.0f));
        M_2_rotate=glm::rotate(M_2_rotate,(float)glm::radians(30.0 + rotations[depth][5]),glm::vec3(0.0f,0.0f,1.0f));
        auto t2 = trunk(size*0.9*sizemod[depth][1], depth-1, position, V, M_2_rotate);
        meshPush(mesh, t2);

        glm::mat4 M_3_rotate;
        M_3_rotate=glm::rotate(M,(float)glm::radians(0.0f + rotations[depth][6]),glm::vec3(1.0f,0.0f,0.0f));
        M_3_rotate=glm::rotate(M_3_rotate,(float)glm::radians(180.0f + rotations[depth][7]),glm::vec3(0.0f,1.0f,0.0f));
        M_3_rotate=glm::rotate(M_3_rotate,(float)glm::radians(30.0f + rotations[depth][8]),glm::vec3(0.0f,0.0f,1.0f));
        auto t3 = trunk(size*0.9*sizemod[depth][2], depth-1, position, V, M_3_rotate);
        meshPush(mesh, t3);

        glm::mat4 M_4_rotate;
        M_4_rotate=glm::rotate(M,(float)glm::radians(0.0f + rotations[depth][9]),glm::vec3(1.0f,0.0f,0.0f));
        M_4_rotate=glm::rotate(M_4_rotate,(float)glm::radians(270.0f + rotations[depth][10]),glm::vec3(0.0f,1.0f,0.0f));
        M_4_rotate=glm::rotate(M_4_rotate,(float)glm::radians(30.0f  + rotations[depth][11]),glm::vec3(0.0f,0.0f,1.0f));
        auto t4 = trunk(size*0.9*sizemod[depth][3], depth-1, position, V, M_4_rotate);
        meshPush(mesh, t4);
	} else {
      	//green_sphere.draw(V, P, M);
      	auto grn_sphre = green_sphere.genMesh(M);
      	meshPush(mesh, grn_sphre);
	}
	return mesh;
}

Mesh Model::genTree(glm::vec3 position, glm::mat4 V, glm::vec3 scale, glm::mat4 M) {
    M = glm::scale(M, scale);
    return trunk(1.0, DEPTH, position, V, M);   
}
