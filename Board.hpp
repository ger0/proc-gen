#pragma once
#include "Drawable.hpp"

class Board : public Drawable {
private:
	float x; //size
	float z; //size

public:
	Board();
	float get_x() { return x; }
	float get_z() { return z; }
	float get_height(float x, float z); //heigt in exact point
};
