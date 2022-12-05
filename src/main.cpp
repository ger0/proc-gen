#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <array>
#include <memory>

#include "SimplexNoise/SimplexNoise.h"
#include "lodepng/lodepng.h"
#include "shaderprogram.hpp"
#include "marchingcubes.hpp"

using byte = unsigned char;
using uint = unsigned;

#define DEBUG

// unique_ptr alias
template<typename... T>
using Uniq_Ptr = std::unique_ptr<T...>;

using std::array;

constexpr const char* 	WINDOW_TITLE = "GLTRY";
constexpr uint 		WINDOW_W = 800;
constexpr uint 		WINDOW_H = 600;

constexpr uint MAP_W = 256;
constexpr uint MAP_H = 256;

constexpr float ASPECT_RATIO = float(WINDOW_W) / WINDOW_H;
constexpr float FOV = 90.f;

// clipping
constexpr float Z_NEAR = 0.1f;
constexpr float Z_FAR = MAP_W;

constexpr uint NOISE_W = MAP_W + 1;
constexpr uint NOISE_H = MAP_H + 1;
constexpr float noiseScale = 32.f;

template<typename T>
struct Position {
	T x;
	T y;

	Position<T> operator + (glm::vec2 &offset) {
		Position<T> npos;
		npos.x = x + T(offset.x);
		npos.y = y + T(offset.y);
		return npos;
	}

	template<typename Y>
	Position<T> offset(Y off_x, Y off_y) {
		Position<T> npos;
		npos.x = x + T(off_x);
		npos.y = y + T(off_y);
		return npos;
	}
};

GLuint vao;
GLuint vbo;

float deltaTime = 0.f;
float lastFrame = 0.f;

float yaw = -90.f;
float pitch = 0.f;

// mouse position
Position<float> mouseLast;

array<float, NOISE_W * NOISE_H> noise;

array<float, NOISE_W * NOISE_H * 32> noise3D;

// vertex + normal data (6 + 6)
// {pos{x, y, z}, norm{x, y, z}} for every single vertex
array<glm::vec3, (6 + 6) * MAP_W * MAP_H> verts;
array<glm::vec3, (6 + 6) * MAP_W * MAP_W * MAP_H> verts3D;

const float cameraSpeed = 0.05f;

glm::vec3 cameraPos = glm::vec3(0.f, 0.f, 3.f);
glm::vec3 cameraTarget = glm::vec3(0.f, 0.f, 0.f);
glm::vec3 cameraDir = glm::normalize(cameraPos - cameraTarget);

glm::vec3 upVector = glm::vec3(0.f, 1.f, 0.f);
glm::vec3 cameraRight = glm::normalize(glm::cross(upVector, cameraDir));
glm::vec3 cameraUp = glm::cross(cameraDir, cameraRight);

void drawScene(ShaderProgram* sp, GLFWwindow *window) {
	// time
	float currentFrame = glfwGetTime();
	deltaTime = currentFrame - lastFrame;
	lastFrame = currentFrame;

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	sp->use();

	glm::mat4 M = glm::mat4(1.f);
	glm::mat4 V = glm::lookAt(cameraPos, cameraPos + cameraTarget, cameraUp);
	glm::mat4 P = glm::perspective(glm::radians(FOV), ASPECT_RATIO, Z_NEAR, Z_FAR);

	glUniformMatrix4fv(sp->u("M"), 1, false, glm::value_ptr(M));
	glUniformMatrix4fv(sp->u("V"), 1, false, glm::value_ptr(V));
	glUniformMatrix4fv(sp->u("P"), 1, false, glm::value_ptr(P));

	glBindVertexArray(vao);
	glDrawArrays(GL_TRIANGLES, 0, (verts.size() / 2) * 3 );

	glfwSwapBuffers(window);
}

void errCallback(int error, const char *description) {
	fputs(description, stderr);
}

void mouseCallback(GLFWwindow *window, double xpos, double ypos) {
	float xoffset = xpos - mouseLast.x;
	// reversed since y-coordinates range from bottom to top
	float yoffset = mouseLast.y - ypos; 
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
	direction.x = cos(glm::radians(yaw)) * cos(glm::radians(pitch));
	direction.y = sin(glm::radians(pitch));
	direction.z = sin(glm::radians(yaw)) * cos(glm::radians(pitch));
	cameraTarget = glm::normalize(direction);
#ifdef DEBUG
	printf("CameraPointing: %f %f %f\n", 
			cameraTarget.x, cameraTarget.y, cameraTarget.z);
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
	printf("Pos: x: %f, y: %f, z: %f\n", 
			cameraPos.x, cameraPos.y, cameraPos.z);
#endif
}

// DEBUGGING INFO
void GLAPIENTRY MessageCallback(GLenum source, GLenum type, GLuint id,
		GLenum severity, GLsizei length,
		const GLchar *message, const void *userParam) {
  	fprintf(stderr, "GL CALLBACK: %s type = 0x%x, severity = 0x%x, message = %s\n",
			(type == GL_DEBUG_TYPE_ERROR ? "** GL ERROR **" : ""), 
			type, severity, message);
}

void windowResizeCallback(GLFWwindow* window, int width, int height) {
    if (height == 0)
	return;
    glViewport(0, 0, width, height);
}

ShaderProgram* initProgram(GLFWwindow *window) {
#ifdef DEBUG
	glEnable(GL_DEBUG_OUTPUT);
	glDebugMessageCallback(MessageCallback, 0);
#endif

	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);

	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);

	//glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glClearColor(0, 0, 0, 1.f);

	glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
	glfwSetCursorPosCallback(window, mouseCallback);
	glfwSetKeyCallback(window, keyCallback);

    glfwSetWindowSizeCallback(window, windowResizeCallback);

	ShaderProgram* sp = new ShaderProgram("glsl/vshad.vert", "glsl/fshad.frag");

	// generating buffers
	glGenVertexArrays(1, &vao);
	glBindVertexArray(vao);

	glGenBuffers(1, &vbo);
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBufferData(GL_ARRAY_BUFFER, verts.size() * sizeof(float) * 3,
        	(void *)verts.data(), GL_STATIC_DRAW);

	uint count = 3;
	glVertexAttribPointer(sp->a("vertex"), count, GL_FLOAT, GL_FALSE,
        	count * 2 * sizeof(float), 0);
	glEnableVertexAttribArray(sp->a("vertex"));

	glVertexAttribPointer(sp->a("normal"), count, GL_FLOAT, GL_FALSE,
        	count * 2 * sizeof(float), &count);
	glEnableVertexAttribArray(sp->a("normal"));

	return sp;
}

int printFail(const char *message) {
	fprintf(stderr, message, "\n");
	return -1;
}

inline float& noiseAt(Position<uint> pos) {
	// bounds checking
	return noise.at(pos.x + pos.y * NOISE_W);
}

// generates a normal vector for a single vertex by 
// checking noise values around the vertex
void genNormal(uint offset, Position<uint> pos) {
	// ignore vertices on the edge of the map
    if (!(pos.x > 0 && pos.y > 0 && pos.x < (MAP_W) && pos.y < (MAP_H))) {
        verts.at(offset + 1) = glm::vec3(0.f, -1.f, 0.f);
        return;
    }

	// derivatives to figure out normal by using cross product
	glm::vec3 dx(
		2.f, noiseAt(pos.offset(1, 0)) - noiseAt(pos.offset(-1, 0)), 0.f);
	glm::vec3 dz(
		0.f, noiseAt(pos.offset(0, -1)) - noiseAt(pos.offset(0, 1)), 2.f);

    glm::vec3 norm;
	norm = glm::cross(dz, dx);
	norm = glm::normalize(norm);
	verts.at(offset + 1) = norm;
}

// generates mesh with normals out of a noise into verts vector
void genMesh() {
	constexpr uint step = 12; // 6 vertices + 6 normals
	for (uint i = 0; i < (MAP_H * MAP_W - step) * step; i+= step) {
		uint frac = i / step;
		Position<uint> pos = {frac % MAP_W, frac / MAP_W};
		/* vertice offsets to create 2 polygons for one quad
 	 	 * (0,0)---(1,0)
 	 	 *  |     /  |
 	 	 *  |    /   |
 	 	 *  |   /    |
 	 	 * (0,1)---(1,1)
 	 	 */
		using t = glm::vec2; 
		auto offsets = {
			t(0,0), t(0,1), t(1,0), // first triangle
			t(1,0), t(0,1), t(1,1)  // second triangle
		};
		uint idx = 0;
		for (auto offset : offsets) {
			glm::vec3 vert;
			vert.x = pos.x + offset.x;
			vert.y = noiseAt(pos + offset);
			// here the 3D Z coord is mapped to the Y coord on the 2D noise plane
			// the Y coord is the noise value at (X, Z)
			vert.z = pos.y + offset.y;
			verts.at(i + idx) = vert;
			idx += 2;
			genNormal(i + idx, pos);
		}
	}
}

void genMesh3D() {
	constexpr uint step = 12; // 6 vertices + 6 normals

}

void saveNoisePNG(const char* fname, array<float, NOISE_W * NOISE_H>& noise) {
	array<byte, NOISE_W * NOISE_H> picture;
	float max = -99999.f;
	for (uint i = 0; i < noise.size(); i++) {
		auto val = noise.at(i) * 128.f / noiseScale;
		if (val > max) max = val;
		picture.at(i) = noise.at(i) * 128.f / noiseScale;
	}
	printf("MAXXXX: %f\n", max);
    lodepng_encode_file(fname, picture.data(), 
    		NOISE_W, NOISE_H, LCT_GREY, 8);
}

void genWavesMap() {
	auto val = 0.f;
	for (uint y = 0; y <= MAP_H; y++) {
		val += 0.5f;
		if (y % 32 == 0) val = 0.f;
    	for (uint x = 0; x <= MAP_W; x++) {
			noiseAt({x, y}) = val;
		}
	}
#ifdef DEBUG
	saveNoisePNG("waves.png", noise);
#endif
}

void genSimplexMap() {
	SimplexNoise noiseGen(3.f / MAP_W, 1.f, 2.f, 0.5f);
	// noise generation
	for (uint y = 0; y <= MAP_H; y++) {
    	for (uint x = 0; x <= MAP_W; x++) {
    		auto val = noiseGen.fractal(4, x, y) + 1.f;
			noiseAt({x, y}) = noiseScale * val;
		}
	}
#ifdef DEBUG
    saveNoisePNG("noise.png", noise);
#endif
}

void fill3DNoise(array<float, NOISE_W * NOISE_H * 32>& noise) {
	SimplexNoise noiseGen(3.f / MAP_W, 1.f, 2.f, 0.5f);
	// noise generation
	for (uint y = 0; y < 32; y++) {
		array<byte, NOISE_W * NOISE_H> picture;
		for (uint z = 0; z <= MAP_H; z++) {
			for (uint x = 0; x <= MAP_W; x++) {
				auto val = noiseGen.fractal(4, x, y, z);
				noise.at(y * (NOISE_W * NOISE_H) + 
					(x + z * NOISE_W)) = noiseScale * val;
			}
		}
#ifdef DEBUG
		// outputnoise directory must exist
		constexpr uint size = NOISE_W * NOISE_H;
		array<float, size> noiseSlice;
		std::copy(noise.begin() + size * y,
				noise.begin() + size * (y + 1),
				noiseSlice.begin());
		char str[21 + 48 / 8];
		sprintf(str, "outputnoise/noise-%u.png", y);
		saveNoisePNG(str, noiseSlice);
#endif
	}
}

int main() {
	// generate a heightmap and a mesh out of it
	genSimplexMap();
	genMesh();

	fill3DNoise(noise3D);

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
	// pointer to the currently used GLFWwindow
	Uniq_Ptr<GLFWwindow, decltype(windowDestroyer)> window(
				glfwCreateWindow(WINDOW_W, WINDOW_H, WINDOW_TITLE, NULL, NULL),
				windowDestroyer);

	if (!window.get()) return printFail("Cannot create window.");

	glfwMakeContextCurrent(window.get());
	glfwSwapInterval(1);
	glViewport(0, 0, WINDOW_W, WINDOW_H);

	if (glewInit() != GLEW_OK) return printFail("Cannot initialise GLEW.");

	auto shaderDestroyer = [&](ShaderProgram* sp) {
		array<GLuint, 2> buffers({vao, vbo});
		glDeleteBuffers(buffers.size(), buffers.data());
		delete sp;
	};
	// pointer to a ShaderProgram
	Uniq_Ptr<ShaderProgram, decltype(shaderDestroyer)> sp(
			initProgram(window.get()),
			shaderDestroyer);

	while (!glfwWindowShouldClose(window.get())) {
    	glfwPollEvents();
		drawScene(sp.get(), window.get());
	}
}
