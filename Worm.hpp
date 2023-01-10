#pragma once
#include "Drawable.hpp"
#include "Board.hpp"
#include "Movable.hpp"
#include "Camera.hpp"

class Camera;
class Worm : public Movable, public Drawable{
private:
	std::string name;
	int life;
	Board* board;
	Camera* camera;

public:
	Worm(std::string name, Board* board, Camera* camera, const std::string& obj_filename, std::vector<const char*> tex_filenames);
	void update(float speed, float angle_speed, double _time); // update position using speed
	void damage(int how_much); // damage when hit with bullet
	static int count_worms; //number of worms on the board
	int get_life();
};