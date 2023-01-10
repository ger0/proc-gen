#pragma once
#include <vector>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <iostream>


class Cylinder{
private:
  std::vector<glm::vec4> verts; // vertices
  std::vector<glm::vec4> norms; // normal vectors
  std::vector<unsigned int> indices; // indices to draw correct triangles

public:
  Cylinder(){};
  Cylinder(float sectorCount, float r, float h);
  void generate_indicis_of_cylinder(float sectorCount, float r, float h);
  void draw(GLFWwindow* window, glm::mat4 V, glm::mat4 P, glm::mat4 M);
};