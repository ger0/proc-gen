#pragma once
#include "Movable.hpp"
#include "Worm.hpp"
class Worm;
class Camera : public Movable {
private:
	bool walking_mode; // follow the worm or move (in the aiming mode)

public:
	glm::vec3 nose_vector;
	Camera();
	void change_mode(Worm* active_worm);
	void update_pos(glm::vec3 _pos, float delta_angle_x);  // move with worm (walking_mode) or change position after change_mode
	bool get_mode() { return walking_mode; }    //potrzebne żeby kamera poruszała się razem z robaczkiem
	void set_angle_y_restricted(float _angle_y); // block the camera field of view while aiming
};
