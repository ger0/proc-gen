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

// temporary vector to be removed
using std::array, std::vector;

constexpr const char* 	WINDOW_TITLE = "GLTRY";
constexpr uint 		WINDOW_W = 800;
constexpr uint 		WINDOW_H = 600;

constexpr uint MAP_W = 32;
constexpr uint MAP_H = 16;

constexpr float ASPECT_RATIO = float(WINDOW_W) / WINDOW_H;
constexpr float FOV = 90.f;

// clipping
constexpr float Z_NEAR = 0.1f;
constexpr float Z_FAR = 256.f;

constexpr uint NOISE_W = MAP_W + 1;
constexpr uint NOISE_H = MAP_H + 1;
constexpr float noiseScale = 6.f;

template<typename T>
struct Position {
	T x;
	T y;
	T z;

	Position<T> operator + (glm::vec2 &offset) {
		Position<T> npos;
		npos.x = x + T(offset.x);
		npos.y = y + T(offset.y);
		return npos;
	}

	Position<T> operator + (Position<T> &offset) {
		Position<T> npos;
		npos.x = x + offset.x;
		npos.y = y + offset.y;
		npos.z = z + offset.z;
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

struct Vertex {
	glm::vec3 position;
	glm::vec3 normal;
};

GLuint vao;
GLuint vbo;

float deltaTime = 0.f;
float lastFrame = 0.f;

float yaw = -90.f;
float pitch = 0.f;

// mouse position
Position<float> mouseLast;

array<float, NOISE_W * NOISE_W * NOISE_H> noise3D;

// vertex + normal data (6 + 6)
// {pos{x, y, z}, norm{x, y, z}} for every single vertex
vector<glm::vec3> verts3D;

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
	glDrawArrays(GL_TRIANGLES, 0, (verts3D.size() / 2) * 3 );

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

/* #ifdef DEBUG
	printf("CameraPointing: %f %f %f\n", 
			cameraTarget.x, cameraTarget.y, cameraTarget.z);
#endif */
}

void keyCallback(GLFWwindow *window, int key, int scancode, int act, int mod) {
	const float cameraSpeed = 10.f * deltaTime;

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

	// glEnable(GL_CULL_FACE);
	// glCullFace(GL_BACK);

	// glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
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
	// glBufferData(GL_ARRAY_BUFFER, verts.size() * sizeof(float) * 3,
 //        	(void *)verts.data(), GL_STATIC_DRAW);
	glBufferData(GL_ARRAY_BUFFER, verts3D.size() * sizeof(glm::vec3),
         	(void*)verts3D.data(), GL_STATIC_DRAW);

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

// normal generation for 2D heightmap (to be removed)
// checking noise values around the vertex
/* void genNormal(uint offset, Position<uint> pos) {
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
}*/

// converts Position<uint> into a reference to an element in the noise array
float inline& noise3DAt(Position<uint> &pos) {
	return noise3D.at(pos.x + pos.z * NOISE_W + pos.y * (NOISE_W * NOISE_W));
}

void genMesh3D() {
	constexpr uint step = 12; // 6 vertices + 6 normals
	for (uint i = 0; i < (MAP_W * MAP_W * MAP_H) * step; i += step) {
		uint frac = i / step;
		auto xz = frac % (MAP_W * MAP_W); // index on the <x,z> plane
		Position<uint> pos = {
			.x = xz % MAP_W,
			.y = frac / (MAP_W * MAP_W),
			.z = xz / MAP_W
		};
		using P = Position<uint>;

		// cube in binary: (v7, v6, v5, v4, v3, v2, v1, v0)
		byte cube = 0x0;
		// offset in binary: (z, y, x)
		for (byte offset = 0b000; offset <= 0b111; offset++) {
			Position<uint> nPos = {
				.x = pos.x + (offset & 0b001),
				.y = pos.y + ((offset & 0b010) >> 1),
				.z = pos.z + ((offset & 0b100) >> 2),
			};
			float noiseVal = noise3DAt(nPos);
			constexpr float threshold = 0.f;
			if (noiseVal > threshold) {
				cube |= (1 << offset);
			}
		}
		auto mask = march_cubes::EdgeMasks.at(cube);
		// array of triplets of the cube edges, each triplet = 1 triangle
		auto edges = march_cubes::TriangleTable.at(cube);

		// array of vertices for a triangle
		for (uint& edge : edges) {
			// check if theres no more edges left
			if (edge == march_cubes::x) {
				break;
			}

			auto verts = march_cubes::EdgeVertexIndices.at(edge);

			glm::vec3 vertex;

			vertex.x = pos.x + (((verts[0] & 0b001) >> 0) + ((verts[1] & 0b001) >> 0)) / 2.f;
			vertex.y = pos.y + (((verts[0] & 0b010) >> 1) + ((verts[1] & 0b010) >> 1)) / 2.f;
			vertex.z = pos.z + (((verts[0] & 0b100) >> 2) + ((verts[1] & 0b100) >> 2)) / 2.f;

			verts3D.push_back(vertex);
			verts3D.push_back(glm::vec3(0,1,0));
		}
	}
}

void saveNoisePNG(const char* fname, array<float, NOISE_W * NOISE_W>& noise) {
	array<byte, NOISE_W * NOISE_W> picture;
	float max = -99999.f;
	for (uint i = 0; i < noise.size(); i++) {
		auto val = noise.at(i) * 128.f / noiseScale;
		if (val > max) max = val;
		picture.at(i) = noise.at(i) * 128.f / noiseScale;
	}
	lodepng_encode_file(fname, picture.data(), 
			NOISE_W, NOISE_W, LCT_GREY, 8);
}
void fill3DNoise() {
	SimplexNoise noiseGen(0.1f, 1.f, 2.f, 0.5f);
	constexpr uint octaves = 9;
	// noise generation
	for (uint y = 0; y < NOISE_H; y++) {
		array<byte, NOISE_W * NOISE_W> picture;
		for (uint z = 0; z <= MAP_W; z++) {
			for (uint x = 0; x <= MAP_W; x++) {
				auto val = noiseGen.fractal(octaves, x, y, z);
				auto pos = Position<uint>{x, y, z};
				noise3DAt(pos) = val;
			}
		}
#ifdef DEBUG
		// outputnoise directory must exist
		constexpr uint size = NOISE_W * NOISE_W;
		array<float, size> noiseSlice;
		std::copy(noise3D.begin() + size * y,
				noise3D.begin() + size * (y + 1),
				noiseSlice.begin());
		char str[21 + 48 / 8];
		sprintf(str, "outputnoise/noise-%u.png", y);
		saveNoisePNG(str, noiseSlice);
#endif
	}
}

void fillNoiseTest() {
	for (uint y = 0; y < NOISE_H; y++) {
		array<byte, NOISE_W * NOISE_W> picture;
		for (uint z = 0; z <= MAP_W; z++) {
			for (uint x = 0; x <= MAP_W; x++) {
				auto pos = Position<uint>{x, y, z};
				noise3DAt(pos) = -1.f;
			}
		}
		auto pos = Position<uint>{0, 0, 0};
		noise3DAt(pos) = -1.f;
#ifdef DEBUG
		// outputnoise directory must exist
		constexpr uint size = NOISE_W * NOISE_W;
		array<float, size> noiseSlice;
		std::copy(noise3D.begin() + size * y,
				noise3D.begin() + size * (y + 1),
				noiseSlice.begin());
		char str[21 + 48 / 8];
		sprintf(str, "outputnoise/noise-%u.png", y);
		saveNoisePNG(str, noiseSlice);
#endif
	}
}

int main() {
	fill3DNoise();
	//fillNoiseTest();
	genMesh3D();

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
