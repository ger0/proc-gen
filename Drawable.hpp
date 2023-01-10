#pragma once
#include "Everything.hpp"
#include "model.hpp"


class Drawable: virtual public Everything {
private:
  Model model;     
  glm::vec3 scale; 
public:
	Drawable(glm::vec3 _scale);
	void draw(GLFWwindow* window, glm::mat4 V);
	glm::mat4 calc_M_matrix();  // model matrix
};