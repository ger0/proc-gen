#pragma once
#include "Everything.hpp"

class Movable: virtual public Everything {
public:
	void rotate(float angle, float time);
	void turn_right(float angle);
	void move_forward(float amount);
	Movable(){};
};
