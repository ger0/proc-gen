#include "model.hpp"

#include <unistd.h>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <glm/gtx/string_cast.hpp>



float CYLINDER_R = 0.30;
float CYLINDER_H = 2.5; 

int DEPTH = 3;
std::vector<int> operations;
std::vector<std::vector<float>> sizemod;
std::vector<std::vector<float>> rotations;


Model::Model(glm::vec3 _pos, float _angle_x, float _angle_y){
  branche = Cylinder(10.0f, CYLINDER_R, CYLINDER_H); //float sectorCount, float r, float h
  green_sphere = Sphere(8, 8, 2.8, 2.8, glm::vec3(0.0, 0.07, 0.01), glm::vec3(0.5, 0.5, 0.3)); //float sectorCount,float stackCount, float radius, float h
  brown_sphere = Sphere(8, 8 , 0.4, 0.4, glm::vec3(0.06, 0.02, 0.0), glm::vec3(0.5, 0.5, 0.3)); //float sectorCount,float stackCount, float radius, float h

  srand (static_cast <unsigned> (time(0)));
  for (int i=0; i<DEPTH; i++){
    operations.push_back( rand()%4);
  }
  
  for (int j=0; j<DEPTH; j++){
    std::vector<float> inner_sizemod;
    std::vector<float> inner_rotations;
    for (int i=0; i<10; i++){
      inner_sizemod.push_back(0.75 + static_cast <float> (rand()) 
      /( static_cast <float> (RAND_MAX/(1.15-0.75))));
      inner_rotations.push_back(-20 + static_cast <float> (rand()) 
      /( static_cast <float> (RAND_MAX/(20+20))));
    }
    sizemod.push_back(inner_sizemod);
    rotations.push_back(inner_rotations);
  }

}

void Model::main_draw(GLFWwindow* window, glm::mat4 V) {
    draw_tree(window, get_position(), V, glm::vec3(0.8, 0.8, 0.8), glm::mat4(1.0f)); // scale = 0.8 
}

void Model::trunk(float size, int depth, GLFWwindow* window, glm::vec3 position, glm::mat4 V, glm::mat4 M){
    
    int operation = operations[depth];

    if (operation < 1)
    {
        branch_one(size , depth, window, position, V, M);          
    }
    if (operation < 2)
    {
        branch_two(size, depth, window, position, V, M);          
    }
    else if (operation < 3)
    {
        branch_three(size, depth, window, position, V, M);          
    }
    else if (operation < 4)
    {
        branch_four(size, depth, window, position, V, M);          
    }
}

void Model::branch_one(float size, int depth, GLFWwindow* window, glm::vec3 position, glm::mat4 V, glm::mat4 M){

  glm::mat4 P=glm::perspective(glm::radians(50.0f), 1.0f, 1.0f, 120.0f); //compute projection matrix
  
  M=glm::translate(M, position);
  M = glm::scale(M, glm::vec3(size, size, size));
  branche.draw(window, V, P, M);
  
  M = glm::translate(M, glm::vec3(0, CYLINDER_H, 0));

		if (depth > 0) 
		{   
        brown_sphere.draw(window, V, P, M);

        glm::mat4 M_1_rotate;
        M_1_rotate=glm::rotate(M,(float)glm::radians(0.0 + rotations[depth][0]),glm::vec3(1.0f,0.0f,0.0f));
        M_1_rotate=glm::rotate(M_1_rotate,(float)glm::radians(0.0 + rotations[depth][2]),glm::vec3(0.0f,1.0f,0.0f));
        M_1_rotate=glm::rotate(M_1_rotate,(float)glm::radians(0.0 + rotations[depth][1]),glm::vec3(0.0f,0.0f,1.0f));

        trunk(size*0.9*sizemod[depth][0], depth-1, window, position, V, M_1_rotate);
		}
		else
		{
      green_sphere.draw(window, V, P, M);
		}
}



void Model::branch_two(float size, int depth, GLFWwindow* window, glm::vec3 position, glm::mat4 V, glm::mat4 M){
  
  glm::mat4 P=glm::perspective(glm::radians(50.0f), 1.0f, 1.0f, 120.0f); //compute projection matrix
  
  M=glm::translate(M, position);
  M = glm::scale(M, glm::vec3(size, size, size));
  branche.draw(window, V, P, M);

  M = glm::translate(M, glm::vec3(0, CYLINDER_H, 0));
		
		if (depth > 0) 
		{
        brown_sphere.draw(window, V, P, M);


        glm::mat4 M_1_rotate;
        M_1_rotate=glm::rotate(M,(float)glm::radians(0.0 + rotations[depth][0]),glm::vec3(1.0f,0.0f,0.0f));
        M_1_rotate=glm::rotate(M_1_rotate,(float)glm::radians(0.0 + rotations[depth][1]),glm::vec3(0.0f,1.0f,0.0f));
        M_1_rotate=glm::rotate(M_1_rotate,(float)glm::radians(30.0 + rotations[depth][2]),glm::vec3(0.0f,0.0f,1.0f));
        trunk(size*0.9*sizemod[depth][0], depth-1, window, position, V, M_1_rotate);
        
        glm::mat4 M_2_rotate;
        M_2_rotate=glm::rotate(M,(float)glm::radians(0.0 + rotations[depth][3]),glm::vec3(1.0f,0.0f,0.0f));
        M_2_rotate=glm::rotate(M_2_rotate,(float)glm::radians(180.0 + rotations[depth][4]),glm::vec3(0.0f,1.0f,0.0f));
        M_2_rotate=glm::rotate(M_2_rotate,(float)glm::radians(30.0 + rotations[depth][5]),glm::vec3(0.0f,0.0f,1.0f));
        trunk(size*0.9*sizemod[depth][1], depth-1, window, position, V, M_2_rotate);
		}
		else
		{
      green_sphere.draw(window, V, P, M);
		}
}


void Model::branch_three(float size, int depth, GLFWwindow* window, glm::vec3 position, glm::mat4 V, glm::mat4 M){
  
  glm::mat4 P=glm::perspective(glm::radians(50.0f), 1.0f, 1.0f, 120.0f); //compute projection matrix
  
  M=glm::translate(M, position);
  M = glm::scale(M, glm::vec3(size, size, size));
  branche.draw(window, V, P, M);

  M = glm::translate(M, glm::vec3(0, CYLINDER_H, 0));

		if (depth > 0) 
		{

        brown_sphere.draw(window, V, P, M);

        glm::mat4 M_1_rotate;
        M_1_rotate=glm::rotate(M,(float)glm::radians(0.0 + rotations[depth][0]),glm::vec3(1.0f,0.0f,0.0f));
        M_1_rotate=glm::rotate(M_1_rotate,(float)glm::radians(0.0 + rotations[depth][1]),glm::vec3(0.0f,1.0f,0.0f));
        M_1_rotate=glm::rotate(M_1_rotate,(float)glm::radians(30.0 + rotations[depth][2]),glm::vec3(0.0f,0.0f,1.0f));
        trunk(size*0.9*sizemod[depth][0], depth-1, window, position, V, M_1_rotate);

        glm::mat4 M_2_rotate;
        M_2_rotate=glm::rotate(M,(float)glm::radians(0.0f + rotations[depth][3]),glm::vec3(1.0f,0.0f,0.0f));
        M_2_rotate=glm::rotate(M_2_rotate,(float)glm::radians(120.0 + rotations[depth][4]),glm::vec3(0.0f,1.0f,0.0f));
        M_2_rotate=glm::rotate(M_2_rotate,(float)glm::radians(30.0 + rotations[depth][5]),glm::vec3(0.0f,0.0f,1.0f));
        trunk(size*0.9*sizemod[depth][1], depth-1, window, position, V, M_2_rotate);//sizemod[1]

        glm::mat4 M_3_rotate;
        M_3_rotate=glm::rotate(M,(float)glm::radians(0.0f + rotations[depth][6]),glm::vec3(1.0f,0.0f,0.0f));
        M_3_rotate=glm::rotate(M_3_rotate,(float)glm::radians(240.0f + rotations[depth][7]),glm::vec3(0.0f,1.0f,0.0f));
        M_3_rotate=glm::rotate(M_3_rotate,(float)glm::radians(30.0f + rotations[depth][8]),glm::vec3(0.0f,0.0f,1.0f));
        trunk(size*0.9*sizemod[depth][2], depth-1, window, position, V, M_3_rotate);
		}
		else
		{
      green_sphere.draw(window, V, P, M);
		}


}
void Model::branch_four(float size, int depth, GLFWwindow* window, glm::vec3 position, glm::mat4 V, glm::mat4 M){

  glm::mat4 P=glm::perspective(glm::radians(50.0f), 1.0f, 1.0f, 120.0f); //compute projection matrix
  
  M=glm::translate(M, position);
  M = glm::scale(M, glm::vec3(size, size, size));
  branche.draw(window, V, P, M);

  M = glm::translate(M, glm::vec3(0, CYLINDER_H, 0));

	if (depth > 0) 
		{
        brown_sphere.draw(window, V, P, M);

        glm::mat4 M_1_rotate;
        M_1_rotate=glm::rotate(M,(float)glm::radians(0.0 + rotations[depth][0]),glm::vec3(1.0f,0.0f,0.0f));
        M_1_rotate=glm::rotate(M_1_rotate,(float)glm::radians(0.0 + rotations[depth][1]),glm::vec3(0.0f,1.0f,0.0f));
        M_1_rotate=glm::rotate(M_1_rotate,(float)glm::radians(30.0 + rotations[depth][2]),glm::vec3(0.0f,0.0f,1.0f));
        trunk(size*0.9*sizemod[depth][0], depth-1, window, position, V, M_1_rotate);

        glm::mat4 M_2_rotate;
        M_2_rotate=glm::rotate(M,(float)glm::radians(0.0f + rotations[depth][3]),glm::vec3(1.0f,0.0f,0.0f));
        M_2_rotate=glm::rotate(M_2_rotate,(float)glm::radians(90.0 + rotations[depth][4]),glm::vec3(0.0f,1.0f,0.0f));
        M_2_rotate=glm::rotate(M_2_rotate,(float)glm::radians(30.0 + rotations[depth][5]),glm::vec3(0.0f,0.0f,1.0f));
        trunk(size*0.9*sizemod[depth][1], depth-1, window, position, V, M_2_rotate);

        glm::mat4 M_3_rotate;
        M_3_rotate=glm::rotate(M,(float)glm::radians(0.0f + rotations[depth][6]),glm::vec3(1.0f,0.0f,0.0f));
        M_3_rotate=glm::rotate(M_3_rotate,(float)glm::radians(180.0f + rotations[depth][7]),glm::vec3(0.0f,1.0f,0.0f));
        M_3_rotate=glm::rotate(M_3_rotate,(float)glm::radians(30.0f + rotations[depth][8]),glm::vec3(0.0f,0.0f,1.0f));
        trunk(size*0.9*sizemod[depth][2], depth-1, window, position, V, M_3_rotate);

        glm::mat4 M_4_rotate;
        M_4_rotate=glm::rotate(M,(float)glm::radians(0.0f + rotations[depth][9]),glm::vec3(1.0f,0.0f,0.0f));
        M_4_rotate=glm::rotate(M_4_rotate,(float)glm::radians(270.0f + rotations[depth][10]),glm::vec3(0.0f,1.0f,0.0f));
        M_4_rotate=glm::rotate(M_4_rotate,(float)glm::radians(30.0f  + rotations[depth][11]),glm::vec3(0.0f,0.0f,1.0f));
        trunk(size*0.9*sizemod[depth][3], depth-1, window, position, V, M_4_rotate);
		}
		else
		{
      green_sphere.draw(window, V, P, M);
		}

}


void Model::draw_tree(GLFWwindow* window, glm::vec3 position, glm::mat4 V, glm::vec3 scale, glm::mat4 M){
    
    M = glm::scale(M, scale);
    trunk(1.0, DEPTH, window, position, V, M);   
    // glm::mat4 P=glm::perspective(glm::radians(50.0f), 1.0f, 1.0f, 120.0f); //compute projection matrix
 
    // std::cout<<glm::to_string(V)<<std::endl;

    // branche = Cylinder(7.0f, CYLINDER_R+ 2, CYLINDER_H+2); //float sectorCount, float r, float h
    // branche.draw(window, V, P, M);
}


//////////////////////////////////////////////////////////////////////


void Model::draw(GLFWwindow* window,float angle_x,float angle_y, glm::vec3 position, glm::mat4 V, glm::vec3 scale){
  
    glm::mat4 P=glm::perspective(glm::radians(50.0f), 1.0f, 1.0f, 120.0f); //compute projection matrix

    glm::mat4 M=glm::mat4(1.0f);
    M=glm::translate(M, position);
    // std::cout << position.x << position.y << position.y << std::endl;
    M=glm::scale(M, scale);
    M=glm::rotate(M,angle_y,glm::vec3(1.0f,0.0f,0.0f)); //Compute model matrix
    M=glm::rotate(M,angle_x,glm::vec3(0.0f,1.0f,0.0f));
    
    M=glm::rotate(M,(float)glm::radians(-90.0f),glm::vec3(1.0f,0.0f,0.0f));

    // draw all meshes
    //for (int i=0; i<meshes.size(); i++){
      branche.draw(window, V, P, M);
    //}
    
}