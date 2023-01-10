#pragma once
#include <vector>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <iostream>

#include "main.hpp"


struct Cylinder{
  	std::vector<glm::vec3> verts; // vertices
  	std::vector<glm::vec3> norms; // normal vectors
  	std::vector<glm::vec4> cols;  // color
  	std::vector<unsigned int> indices; // indices to draw correct triangles

  	Cylinder(){};
  	Cylinder(float sectorCount, float r, float h);
  	void generate_indicis_of_cylinder(float sectorCount, float r, float h);
  	void draw(GLFWwindow* window, glm::mat4 V, glm::mat4 P, glm::mat4 M);

	// generates mesh out of the data with predefined colour
	Mesh genMesh(glm::mat4 M, glm::vec4 color);
};
