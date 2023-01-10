#include "Camera.hpp"

Camera::Camera() : Movable(), Everything() {
    walking_mode = true;
    nose_vector = glm::vec3(0.0f, 1.0f, 0.0f);
}

void Camera::change_mode(Worm* active_worm) {
  //switch to different mode
    if (this->walking_mode == true) {
        //opcjonalnie zapisz poprzednie ustawienie kamery wzgl. worma |+deklaracja globalna pos_save/zwracanie
        //pos_save = pos - active_worm->get_position();
        set_position(active_worm->get_position());
    }   // change to aiming
    else {
        // wróć do poprzedniego ustawienia
        // pos = active_worm->get_position() + pos_save;
        set_position(active_worm->get_position() + glm::vec3(2, 10, -15));
    }   // change to walking
    walking_mode = -walking_mode;
}

void Camera::update_pos(glm::vec3 worm_pos, float angle_x) {
    float distance = 10;   //distance between worm and the camera
    set_position(worm_pos + glm::vec3(-sin(angle_x) * distance, 7, -cos(angle_x) * distance));
    set_angle_x(angle_x);
}

void Camera::set_angle_y_restricted(float _angle_y) {
    if (abs(_angle_y) <= 0.6) {
        this->set_angle_y(_angle_y);
    }
}