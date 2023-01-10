#pragma once
#include <vector>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <iostream>

#include "main.hpp"


struct Sphere {
  	std::vector<glm::vec3> verts; // vertices
  	std::vector<glm::vec3> norms; // normal vectors
  	std::vector<glm::vec4> cols;  // colors 
  	//std::vector<glm::vec2> texCoords; // texturing coordinates, not use
  	std::vector<unsigned int> indices; // indices to draw correct triangles

  	glm::vec3 ambientColor, diffuseColor;

  	std::vector<int> lineIndices; // ????????

  	Sphere(){};
  	Sphere(float sectorCount,float stackCount, float radius, float h, glm::vec3 dc, glm::vec3 sc);
  	void generate_indicis_of_sphere(float sectorCount, float stackCount, float r, float h);
  	//void draw(GLFWwindow* window, glm::mat4 V, glm::mat4 P, glm::mat4 M);
	
	// generates mesh out of the data with predefined colour
	Mesh genMesh(glm::mat4 M);
};
