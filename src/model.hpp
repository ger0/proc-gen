#define GLM_FORCE_RADIANS

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <vector>
#include "Cylinder.hpp"
#include "Sphere.hpp"


constexpr float PI = 3.1415;

struct Model{
  	Cylinder branche;
  	Sphere brown_sphere, green_sphere;

  	Model(glm::vec3 _pos = glm::vec3(0, 0, 0), float _angle_x = 0, float _angle_y = 0);
  	Mesh branch_one(float size, int depth, glm::vec3 position, glm::mat4 V, glm::mat4 M);
  	Mesh branch_two(float size, int depth, glm::vec3 position, glm::mat4 V, glm::mat4 M);
  	Mesh branch_three(float size, int depth, glm::vec3 position, glm::mat4 V, glm::mat4 M);
  	Mesh branch_four(float size, int depth, glm::vec3 position, glm::mat4 V, glm::mat4 M);

  	Mesh trunk(float size, int depth, glm::vec3 position, glm::mat4 V, glm::mat4 M);

  	glm::vec3 pos; // current position

  	Mesh genTree(glm::vec3 position, glm::mat4 V, glm::vec3 scale, glm::mat4 M);

  	float angle_x; //current rotation angle of the object, x axis
  	float angle_y; //current rotation angle of the object, y axis
};

