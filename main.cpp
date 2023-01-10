#include "includes.hpp"

using std::cout;
using std::endl;

float speed = 0;                // current speed
float angle_speed = 0;          // current angle speed
float camera_angle_speed_x = 0; // current camera angle speed
float camera_angle_speed_y = 0;

const float ANGLE_SPEED = PI / 3; // max angle speed of the camera and worm
const float SPEED = 5;            // max speed of the worm

// game window to be closed
class CloseWindowException : public std::exception {
public:
    CloseWindowException() {}
};

// OpenGL error handling
void error_callback(int error, const char* description) {
    fputs(description, stderr);
}

// arrows for moving or aiming, space for shooting or changing mode
void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods) {
    if (action == GLFW_PRESS) {
            if (key == GLFW_KEY_LEFT) { angle_speed = ANGLE_SPEED; }
            if (key == GLFW_KEY_RIGHT) { angle_speed = -ANGLE_SPEED; }
            if (key == GLFW_KEY_UP) { speed = SPEED; }
            if (key == GLFW_KEY_DOWN) { speed = -SPEED; }
    }
    if (action == GLFW_RELEASE) {
            if (key == GLFW_KEY_LEFT) angle_speed = 0;
            if (key == GLFW_KEY_RIGHT) angle_speed = 0;
            if (key == GLFW_KEY_UP) speed = 0;
            if (key == GLFW_KEY_DOWN) speed = 0;
    }
}


void initOpenGLProgram(GLFWwindow* window) {
    glClearColor(0.2, 0.2, 0.9, 1);
    glEnable(GL_DEPTH_TEST);
    glfwSetKeyCallback(window, keyCallback);
    initShaders();     // start shaders for worms, text and other objects
}

void freeOpenGLProgram(GLFWwindow* window) {
    freeShaders();
}



// draw scene of walking worm
void drawSceneWalking(GLFWwindow* window, Camera* camera, std::vector<Model*> objects, Worm* active_worm) {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glm::vec3 observer = camera->get_position();
    glm::vec3 center = active_worm->get_position() + glm::vec3(0, 3, 0); 	// camera is following active worm
    glm::vec3 nose_vector = glm::vec3(0.0f, 1.0f, 0.0f); //(pionowo prostopadły do osi patrzenia)

    // Calculate view matrix
    glm::mat4 V = glm::lookAt(observer, center, nose_vector);

    // draw all objects
    for (int i = 0; i < objects.size(); i++) {
        objects[i]->main_draw(window, V);
    }
    glfwSwapBuffers(window);
}


// init libraries and create window
GLFWwindow* create_window() {
    GLFWwindow* window;
    glfwSetErrorCallback(error_callback);

    if (!glfwInit()) { //Initialize GLFW library
        fprintf(stderr, "Can't initialize GLFW.\n");
        exit(EXIT_FAILURE);
    }

    window = glfwCreateWindow(1000, 1000, "Beautiful tree", NULL, NULL);

    if (!window) //If no window is opened then close the program
    {
        glfwTerminate();
        exit(EXIT_FAILURE);
    }

    glfwMakeContextCurrent(window);
    glfwSwapInterval(1); //During vsync wait for the first refresh

    GLenum err;
    if ((err = glewInit()) != GLEW_OK) { //Initialize GLEW library
        fprintf(stderr, "Can't initialize GLEW: %s\n", glewGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    initOpenGLProgram(window);

    return window;
}



int main(int argc, char** argv)
{
    srand(time(NULL));

    //create objects
    GLFWwindow* window = create_window();

    Board board = Board();
    Camera camera;
    Worm worm1 = Worm("Napoleon", &board, &camera, worm_obj, worm_blue_textures);
    
    
    Model obstacle_orange = Model();

    obstacle_orange.set_position(glm::vec3(0,0,0));
    camera.update_pos(worm1.get_position(), 0);

    // group objects
    std::vector<Model*> objects = {&obstacle_orange};
    std::vector<Worm*> worms = { &worm1 };

    float angle_x, angle_y;

    glfwSetTime(0); //Zero the timer
    //Main loop
    try {
        while (!glfwWindowShouldClose(window))
        {

            worm1.update(speed, angle_speed, glfwGetTime());

            glfwSetTime(0);

            drawSceneWalking(window, &camera, objects, & worm1);
            
            glfwPollEvents();
            if (glfwWindowShouldClose(window)) {
                throw CloseWindowException();
            }            
        }
    }
    catch (CloseWindowException) {
      cout << "Exit\n";
    }
    

    freeOpenGLProgram(window);

    glfwDestroyWindow(window); //Delete OpenGL context and the window.
    glfwTerminate(); //Free GLFW resources
    exit(EXIT_SUCCESS);
}
