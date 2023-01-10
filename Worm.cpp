#include "Worm.hpp"
int Worm::count_worms = 1;


Worm::Worm(std::string name, Board* board, Camera* camera, const std::string& _model_filename, std::vector<const char*> tex_filenames) :
  Movable(),
  Drawable(glm::vec3(0.5f,0.5f,0.5f)),
  Everything(){
    this->name = name;
    this->life = 100;
    this->board = board;
    this->camera = camera;
    srand(time(NULL));
    int x = 0.0f;
    int z = 0.0f;
    set_position(glm::vec3(x, board->get_height(x, z), z));
    count_worms++;
}


void Worm::update(float speed, float angle_speed, double _time) {
    // movement in world space
    set_angle_x(get_angle_x() + angle_speed * _time);

    float x = get_position()[0] + speed * sin(get_angle_x()) * _time;
    float z = get_position()[2] + speed * cos(get_angle_x()) * _time;
    if (abs(x) < board->get_x() && abs(z) < board->get_z()){
      try {
          float y = board->get_height(x, z);
          set_position(glm::vec3(x, y, z));

          camera->update_pos(get_position(), get_angle_x());
      }
      catch (std::out_of_range) {}
    }
}


void Worm::damage(int how_much) {
    life -= how_much;

    std::cout<< name <<" "<< life<<std::endl;
}

int Worm::get_life() {
    return life;
}
////////////////////////////////////////////////////////////////////

