#include "Drawable.hpp"

Drawable::Drawable(glm::vec3 _scale){
    this->scale = _scale;
    //model = Model();
}


void Drawable::draw(GLFWwindow* window, glm::mat4 V) {
    model.draw(window, get_angle_x(), get_angle_y(), get_position(), V, scale);
}

glm::mat4 Drawable::calc_M_matrix(){
  glm::mat4 M=glm::mat4(1.0f);
  M=glm::translate(M, get_position());
  M=glm::scale(M, scale);
  M=glm::rotate(M, get_angle_y(), glm::vec3(1.0f,0.0f,0.0f));
  M=glm::rotate(M, get_angle_x(), glm::vec3(0.0f,1.0f,0.0f));
  return M;
}
