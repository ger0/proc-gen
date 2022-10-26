#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp> 
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include <cstdlib>
#include <cstdio>

#include "shaderprogram.hpp"

constexpr float PI = 3.141592653589793238462f;
ShaderProgram *sp;
GLuint vao;
GLuint vbo;

constexpr float verts[] = {
   0.f, 0.f, 0.f,
   1.f, 0.f, 0.f,
   1.f, 1.f, 0.f
};

glm::vec3 cameraPos     = glm::vec3(0.f, 0.f, 3.f);
glm::vec3 cameraTarget  = glm::vec3(0.f, 0.f, 0.f);
glm::vec3 cameraDir     = glm::normalize(cameraPos - cameraTarget);

glm::vec3 upVector      = glm::vec3(0.f, 1.f, 0.f);
glm::vec3 cameraRight   = glm::normalize(glm::cross(upVector, cameraDir));
glm::vec3 cameraUp      = glm::cross(cameraRight, cameraDir);

void drawScene(GLFWwindow *window) {
   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
   sp->use();

   glm::mat4 M = glm::mat4(1.f);
   glm::mat4 V = glm::lookAt(
         glm::vec3(0.f, 0.f, -1.f),
         glm::vec3(0.f, 0.f, 0.f),
         glm::vec3(0.f, 1.f, 0.f)
   );
   glm::mat4 P = glm::perspective(
         90.f * PI / 180.f,
         800.f / 600.f, 0.0f, 90.f
   );
   glUniformMatrix4fv(sp->u("M"), 1, false, glm::value_ptr(M));
   glUniformMatrix4fv(sp->u("V"), 1, false, glm::value_ptr(V));
   glUniformMatrix4fv(sp->u("P"), 1, false, glm::value_ptr(P));

   glBindVertexArray(vao);
   glDrawArrays(GL_TRIANGLES, 0, 3);

   glfwSwapBuffers(window);
}


void errCallback(int error, const char *description) {
   fputs(description, stderr);
}

void mouseCallback(GLFWwindow *window, double xpos, double ypos) {
}

void keyCallback(GLFWwindow *window, int key, int scancode, int act, int mod) {
}

void initProgram(GLFWwindow *window) {
   glEnable(GL_BLEND);
   glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
   glClearColor(0, 0, 0, 1);

   glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
   glfwSetCursorPosCallback(window, mouseCallback);
   glfwSetKeyCallback(window, keyCallback);

   sp = new ShaderProgram("vshad.vert", "fshad.frag");

   glGenVertexArrays(1, &vao);
   glBindVertexArray(vao);

   glGenBuffers(1, &vbo);
   glBindBuffer(GL_ARRAY_BUFFER, vbo);
   glBufferData(GL_ARRAY_BUFFER, 3 * 3 * sizeof(float), verts, GL_STATIC_DRAW);

   glVertexAttribPointer(sp->a("vertex"), 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), 0);
   glEnableVertexAttribArray(sp->a("vertex"));
}

void exitFail(const char *message) {
   fprintf(stderr, message, "\n");
   glfwTerminate();
   exit(EXIT_FAILURE);
}

int main() {
   GLFWwindow *window;
   glfwSetErrorCallback(errCallback);
   if (!glfwInit())  exitFail("Cannot initialise GLFW.");

   window = glfwCreateWindow(800, 600, "gltry", NULL, NULL);
   if (!window) {
      exitFail("Cannot create window.");
   }

   glfwMakeContextCurrent(window);
   glfwSwapInterval(1);
   glViewport(0, 0, 800, 600); 

   if (glewInit() != GLEW_OK) exitFail("Cannot initialise GLEW.");
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
