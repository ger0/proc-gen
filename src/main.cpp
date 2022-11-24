#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <memory>

#include "SimplexNoise/SimplexNoise.h"
#include "shaderprogram.hpp"

constexpr const char* 	WINDOW_TITLE = "GLTRY";
constexpr unsigned 		WINDOW_W = 1360;
constexpr unsigned 		WINDOW_H = 768;

constexpr float PI = 3.141592653589793238462f;

constexpr unsigned MAP_W = 128;
constexpr unsigned MAP_H = 128;

// unique_ptr alias
template<typename... T>
using Uniq_Ptr = std::unique_ptr<T...>;

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

// 3 coords
std::vector<float> noise((MAP_W + 1) * (MAP_H + 1));
std::vector<glm::vec3> verts(6 * MAP_W * MAP_H);

glm::vec3 cameraPos = glm::vec3(0.f, 0.f, 1.f);
glm::vec3 cameraTarget = glm::vec3(0.f, 0.f, 0.f);
glm::vec3 cameraDir = glm::normalize(cameraPos - cameraTarget);

glm::vec3 upVector = glm::vec3(0.f, 1.f, 0.f);
glm::vec3 cameraRight = glm::normalize(glm::cross(upVector, cameraDir));
glm::vec3 cameraUp = glm::cross(cameraRight, cameraDir);

void drawScene(ShaderProgram* sp, GLFWwindow *window) {
	// time
	float currentFrame = glfwGetTime();
	deltaTime = currentFrame - lastFrame;
	lastFrame = currentFrame;

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	sp->use();

	glm::mat4 M = glm::mat4(1.f);
	glm::mat4 V = glm::lookAt(cameraPos, cameraPos + cameraTarget, cameraUp);
	glm::mat4 P = glm::perspective(90.f * PI / 180.f, 800.f / 600.f, 0.0f, 90.f);
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
	float yoffset = mouseLast.y -
    	ypos; // reversed since y-coordinates range from bottom to top
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
	printf("Pos: x: %f, y: %f, z: %f\n", -cameraPos.x, -cameraPos.y,
         	-cameraPos.z);
#endif
}

void GLAPIENTRY MessageCallback(GLenum source, GLenum type, GLuint id,
		GLenum severity, GLsizei length,
		const GLchar *message, const void *userParam) {
  	fprintf(stderr,
          	"GL CALLBACK: %s type = 0x%x, severity = 0x%x, message = %s\n",
          	(type == GL_DEBUG_TYPE_ERROR ? "** GL ERROR **" : ""), type, severity,
          	message);
}

ShaderProgram* initProgram(GLFWwindow *window) {
	glEnable(GL_BLEND | GL_DEBUG_OUTPUT);
	glDebugMessageCallback(MessageCallback, 0);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glClearColor(0, 0, 0, 1);

	glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
	glfwSetCursorPosCallback(window, mouseCallback);
	glfwSetKeyCallback(window, keyCallback);

	ShaderProgram* sp = new ShaderProgram("glsl/vshad.vert", "glsl/fshad.frag");

	glGenVertexArrays(1, &vao);
	glBindVertexArray(vao);

	glGenBuffers(1, &vbo);
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBufferData(GL_ARRAY_BUFFER, verts.size() * sizeof(float) * 3,
        	(void *)verts.data(), GL_STATIC_DRAW);

	glVertexAttribPointer(sp->a("vertex"), 3, GL_FLOAT, GL_FALSE,
        	3 * sizeof(float), 0);
	glEnableVertexAttribArray(sp->a("vertex"));

	return sp;
}

int printFail(const char *message) {
	fprintf(stderr, message, "\n");
	return -1;
}

int main() {
	SimplexNoise noiseGen;
	for (unsigned y = 0; y <= MAP_H; y++) {
    	for (unsigned x = 0; x <= MAP_W; x++) {
      		noise[x + y * (MAP_W + 1)] = 2.f * noiseGen.fractal(20, x, y);
    	}
	}

	constexpr unsigned step = 6;
	for (unsigned i = 0; i < (MAP_H * MAP_W) * step; i+= step) {
		unsigned frac = i / step;
    	unsigned x = frac % MAP_W;
    	unsigned y = frac / MAP_W;
    	/*
     	 * (0,0)---(1,0)
     	 *  |     /  |
     	 *  | 	 /   |
     	 *  |   /    |
     	 * (0,1)---(1,1)
     	 */
     	// vertice offsets to create 2 polygons for one quad
		using t = glm::vec2; 
		auto offsets = {t(0,0), t(0,1), t(1,0), t(1,0), t(1,1), t(0,1)};
		unsigned idx = 0;
		for (auto offset : offsets) {
			glm::vec3 vert;
			vert.x = x + offset.x;
			vert.y = noise[(x + offset.x) + (y + offset.y) * MAP_W];
			vert.z = y + offset.y;
			// here the z coord is mapped to the y coord on the noise vector
			verts[i + idx] = vert;
			idx++;
		}
	}
	/*
   	   for (auto item: verts) {
   	   printf("x: %f y: %f z: %f\n", item.x, item.y, item.z);
   	   }
   	*/

	glfwSetErrorCallback(errCallback);
	if (!glfwInit()) {
		glfwTerminate();
    	return printFail("Cannot initialise GLFW.");
    }

	// RAII, this might be overdone here but it still works 
	auto windowDestroyer = [&](GLFWwindow* window) {
		glfwDestroyWindow(window);
		glfwTerminate();
	};
	Uniq_Ptr<GLFWwindow, decltype(windowDestroyer)> 
		window(
			glfwCreateWindow(WINDOW_W, WINDOW_H, WINDOW_TITLE, NULL, NULL),
			windowDestroyer
			);

	if (!window.get()) return printFail("Cannot create window.");

	glfwMakeContextCurrent(window.get());
	glfwSwapInterval(1);
	glViewport(0, 0, WINDOW_W, WINDOW_H);

	if (glewInit() != GLEW_OK) return printFail("Cannot initialise GLEW.");
    
	Uniq_Ptr<ShaderProgram> sp(initProgram(window.get()));

	while (!glfwWindowShouldClose(window.get())) {
    	glfwPollEvents();
    	drawScene(sp.get(), window.get());
	}
}
