#pragma once
#include <glm/glm.hpp>

// every object in game
class Everything {
private:
	glm::vec3 pos; // current position
  float angle_x; //current rotation angle of the object, x axis
  float angle_y; //current rotation angle of the object, y axis

public:
	Everything(glm::vec3 _pos = glm::vec3(0, 0, 0), float _angle_x = 0, float _angle_y = 0);

	void set_position(glm::vec3 _pos) { pos = _pos; }
	glm::vec3 get_position() { return pos; }
  float get_angle_x() { return angle_x; }
  float get_angle_y() { return angle_y; }
  void set_angle_x(float _angle_x) { angle_x = _angle_x; }
  void set_angle_y(float _angle_y) { angle_y = _angle_y; }
};