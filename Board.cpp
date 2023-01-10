#include "Board.hpp"

Board::Board():
  Drawable( glm::vec3(50.0f,50.0f,50.0f)),
  Everything(glm::vec3(0, -26, 0)) {
    x = 29;
    z = 48;
}

float Board::get_height(float x, float z) {
    return 0;
}
