#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <cstdio>
#include <cstdlib>
#include <vector>

#include "SimplexNoise/SimplexNoise.h"
#include "shaderprogram.hpp"

constexpr float PI = 3.141592653589793238462f;
constexpr unsigned MAP_W = 128;
constexpr unsigned MAP_H = 128;

ShaderProgram *sp;
GLuint vao;
GLuint vbo;

float deltaTime = 0.f;
float lastFrame = 0.f;

struct Position {
    float x;
    float y;
};

float yaw = -90.f;
float pitch = 0.f;

// mouse position
Position mouseLast;

//constexpr float verts[] = {0.f, 0.f, 0.f, 1.f, 0.f, 0.f, 1.f, 1.f, 0.f};
// 3 coords
std::vector<float> noise((MAP_W + 1) * (MAP_H + 1));
std::vector<glm::vec3> verts(6 * MAP_W * MAP_H);

glm::vec3 cameraPos = glm::vec3(0.f, 0.f, 1.f);
glm::vec3 cameraTarget = glm::vec3(0.f, 0.f, 0.f);
glm::vec3 cameraDir = glm::normalize(cameraPos - cameraTarget);

glm::vec3 upVector = glm::vec3(0.f, 1.f, 0.f);
glm::vec3 cameraRight = glm::normalize(glm::cross(upVector, cameraDir));
glm::vec3 cameraUp = glm::cross(cameraRight, cameraDir);

void drawScene(GLFWwindow *window) {
    // time
    float currentFrame = glfwGetTime();
    deltaTime = currentFrame - lastFrame;
    lastFrame = currentFrame;

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    sp->use();

    glm::mat4 M = glm::mat4(1.f);
    glm::mat4 V = glm::lookAt(cameraPos, cameraPos + cameraTarget, cameraUp);
    glm::mat4 P =
        glm::perspective(90.f * PI / 180.f, 800.f / 600.f, 0.0f, 90.f);
    const float cameraSpeed = 0.05f;

    glUniformMatrix4fv(sp->u("M"), 1, false, glm::value_ptr(M));
    glUniformMatrix4fv(sp->u("V"), 1, false, glm::value_ptr(V));
    glUniformMatrix4fv(sp->u("P"), 1, false, glm::value_ptr(P));

    glBindVertexArray(vao);
    glDrawArrays(GL_TRIANGLES, 0, verts.size() * 3);

    glfwSwapBuffers(window);
}

void errCallback(int error, const char *description) {
    fputs(description, stderr);
}

void mouseCallback(GLFWwindow *window, double xpos, double ypos) {
    float xoffset = xpos - mouseLast.x;
    float yoffset =
        mouseLast.y - ypos; // reversed since y-coordinates range from bottom to top
    mouseLast.x = xpos;
    mouseLast.y = ypos;

    const float sensitivity = 0.1f;
    xoffset *= sensitivity;
    yoffset *= sensitivity;

    yaw += xoffset;
    pitch += yoffset;
    if (pitch > 89.f) {
        pitch = 89.f;
    }
    if (pitch < -89.f) {
        pitch = -89.f;
    }
    glm::vec3 direction;
    direction.x = cos(glm::radians(-yaw)) * cos(glm::radians(-pitch));
    direction.y = sin(glm::radians(-pitch));
    direction.z = sin(glm::radians(-yaw)) * cos(glm::radians(-pitch));
    cameraTarget = glm::normalize(direction);
#ifdef DEBUG
    printf("Yaw: %f, Pitch: %f\n", glm::radians(yaw), glm::radians(pitch));
#endif
}

void keyCallback(GLFWwindow *window, int key, int scancode, int act, int mod) {
    const float cameraSpeed = 100.f * deltaTime;

    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
        cameraPos += cameraSpeed * cameraTarget;
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
        cameraPos -= cameraSpeed * cameraTarget;
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
        cameraPos -=
            glm::normalize(glm::cross(cameraTarget, cameraUp)) * cameraSpeed;
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
        cameraPos +=
            glm::normalize(glm::cross(cameraTarget, cameraUp)) * cameraSpeed;

#ifdef DEBUG
    printf("Pos: x: %f, y: %f, z: %f\n", -cameraPos.x, -cameraPos.y, -cameraPos.z);
#endif
}

void GLAPIENTRY MessageCallback( GLenum source,
                 GLenum type,
                 GLuint id,
                 GLenum severity,
                 GLsizei length,
                 const GLchar* message,
                 const void* userParam ) {
  fprintf( stderr, "GL CALLBACK: %s type = 0x%x, severity = 0x%x, message = %s\n",
           ( type == GL_DEBUG_TYPE_ERROR ? "** GL ERROR **" : "" ),
            type, severity, message );
}

void initProgram(GLFWwindow *window) {
    glEnable(GL_BLEND | GL_DEBUG_OUTPUT);
    glDebugMessageCallback(MessageCallback, 0);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glClearColor(0, 0, 0, 1);

    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
    glfwSetCursorPosCallback(window, mouseCallback);
    glfwSetKeyCallback(window, keyCallback);

    sp = new ShaderProgram("glsl/vshad.vert", "glsl/fshad.frag");

    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

    glGenBuffers(1, &vbo);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, verts.size() * sizeof(float) * 3,
            (void*)verts.data(), GL_STATIC_DRAW);

    glVertexAttribPointer(sp->a("vertex"), 3, GL_FLOAT, GL_FALSE,
                          3 * sizeof(float), 0);
    glEnableVertexAttribArray(sp->a("vertex"));
}

void exitFail(const char *message) {
    fprintf(stderr, message, "\n");
    glfwTerminate();
    exit(EXIT_FAILURE);
}

int main() {
    SimplexNoise noiseGen;
    for (unsigned y = 0; y <= MAP_H; y++) {
        for (unsigned x = 0; x <= MAP_W; x++) {
            noise[x + y * (MAP_W + 1)] = 2.f * noiseGen.fractal(20, x, y);
        }
    }

    constexpr unsigned step = 6;
    for (unsigned i = 0; i < (MAP_H * MAP_W) * step; i += step) {

        unsigned frac = i / step;
        unsigned x = frac % MAP_W;
        unsigned y = frac / MAP_W;

        /*
         * (0)---(1,4)
         *  |   / |
         *  |  /  |
         *  | /   | 
         * (2,3)--(5)
         */

        // first triangle
        for (unsigned off_x: {0, 1}) {
            // (x, y): {(0,0); (1,0); (0,1)} 
            glm::vec3 vert;

            vert.x = float(x + off_x);
            vert.y = noise[(x + off_x) + y * MAP_W];
            vert.z = float(y);

            verts[i + off_x] = vert;

            // y + 1
            if (off_x == 0) {
                auto yy = y + 1.f;
                glm::vec3 vert2;

                vert2.x = float(x);
                vert2.y = noise[(x + off_x) + yy * MAP_W];
                vert2.z = float(yy);

                verts[i + off_x + 2] = vert2;
            }
        }

        // second triangle
        for (unsigned off_x: {0, 1}) {
            // (x, y): {(0,1); (1,0); (1,1)} 
            glm::vec3 vert;

            auto off = 3;
            auto yy = y + 1;

            vert.x = float(x + off_x);
            vert.y = noise[(x + off_x) + yy * MAP_W];
            vert.z = float(yy);

            verts[i + 2 * off_x + off] = vert;

            // y - 1
            if (off_x == 1) {
                glm::vec3 vert2;

                vert2.x = float(x + off_x);
                vert2.y = noise[(x + off_x) + y * MAP_W];
                vert2.z = float(y);

                verts[i + off_x + off] = vert2;
            }
        }
    }

    /*
    for (auto item: verts) {
        printf("x: %f y: %f z: %f\n", item.x, item.y, item.z);
    }
    */

    GLFWwindow *window;
    glfwSetErrorCallback(errCallback);
    if (!glfwInit())
        exitFail("Cannot initialise GLFW.");

    window = glfwCreateWindow(800, 600, "gltry", NULL, NULL);
    if (!window) {
        exitFail("Cannot create window.");
    }

    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);
    glViewport(0, 0, 800, 600);

    if (glewInit() != GLEW_OK)
        exitFail("Cannot initialise GLEW.");
    initProgram(window);

    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();
        drawScene(window);
    }

    // free resources
    delete sp;
    glfwDestroyWindow(window);
    glfwTerminate();

    exit(EXIT_SUCCESS);
}
