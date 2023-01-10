#define GLM_FORCE_RADIANS

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include "shaderprogram.h"
#include "lodepng.h"
#include <vector>
#include "Cylinder.hpp"
#include "Sphere.hpp"


const float PI = 3.1415;

class Model{
private:
  
  Cylinder branche;
  Sphere brown_sphere, green_sphere;

  // this atributes will be not needed after marge 
  glm::vec3 pos; // current position
  float angle_x; //current rotation angle of the object, x axis
  float angle_y; //current rotation angle of the object, y axis


public:

  Model(glm::vec3 _pos = glm::vec3(0, 0, 0), float _angle_x = 0, float _angle_y = 0);
  void main_draw(GLFWwindow* window, glm::mat4 V);
  void draw_tree(GLFWwindow* window, glm::vec3 position, glm::mat4 V, glm::vec3 scale, glm::mat4 M);

  void branch_one(float size, int depth, GLFWwindow* window, glm::vec3 position, glm::mat4 V, glm::mat4 M);
  void branch_two(float size, int depth, GLFWwindow* window, glm::vec3 position, glm::mat4 V, glm::mat4 M);
  void branch_three(float size, int depth, GLFWwindow* window, glm::vec3 position, glm::mat4 V, glm::mat4 M);
  void branch_four(float size, int depth, GLFWwindow* window, glm::vec3 position, glm::mat4 V, glm::mat4 M);
  void trunk(float size, int depth, GLFWwindow* window, glm::vec3 position, glm::mat4 V, glm::mat4 M);


  // this function will be not needed after marge 
  void draw(GLFWwindow* window,float angle_x,float angle_y, glm::vec3 position, glm::mat4 V, glm::vec3 scale);
	void set_position(glm::vec3 _pos) { pos = _pos; }
	glm::vec3 get_position() { return pos; }
  float get_angle_x() { return angle_x; }
  float get_angle_y() { return angle_y; }
  void set_angle_x(float _angle_x) { angle_x = _angle_x; }
  void set_angle_y(float _angle_y) { angle_y = _angle_y; }  
};
