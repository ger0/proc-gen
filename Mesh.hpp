#pragma once
#include <vector>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <iostream>


class Mesh{
private:
  std::vector<glm::vec4> verts; // vertices
  // std::vector<glm::vec4> norms; // normal vectors
  // std::vector<glm::vec2> texCoords; // texturing coordinates
  std::vector<unsigned int> indices; // indices to draw correct triangles
  // std::vector<glm::vec4> tangents; // przestrzeń styczna (for normal mapping)
public:
  Mesh();
  void generate_indicis_of_cylinder(float sectorCount, float r, float h);
  void draw(GLFWwindow* window, glm::mat4 V, glm::mat4 P, glm::mat4 M, std::vector<GLuint> textures);
};