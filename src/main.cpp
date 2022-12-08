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

constexpr uint NOISE_W = MAP_W + 2;
constexpr uint NOISE_H = MAP_H + 2;

constexpr float noiseScale = 6.f;
constexpr float THRESHOLD = 0.f;

template<typename T>
struct Position {
	T x;
	T y;
	T z;

	Position<T> operator + (Position<T> &offset) {
		Position<T> npos;
		npos.x = x + offset.x;
		npos.y = y + offset.y;
		npos.z = z + offset.z;
		return npos;
	}

	template<typename Y>
	Position<T> offset(Y off_x, Y off_y, Y off_z) {
		Position<T> npos;
		npos.x = x + T(off_x);
		npos.y = y + T(off_y);
		npos.z = z + T(off_z);
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

// mouse position
Position<float> mouseLast;

array<float, NOISE_W * NOISE_W * NOISE_H> noise;

// vertex + normal data (6 + 6)
// {pos{x, y, z}, norm{x, y, z}} for every single vertex
vector<glm::vec3> verts3D;

struct Camera {
	glm::vec3 pos = glm::vec3(0.f, 0.f, 3.f);
	glm::vec3 target = glm::vec3(0.f, 0.f, 0.f);
	glm::vec3 dir = glm::normalize(pos - target);

	const float speed = 0.05f;

	glm::vec3 right = glm::normalize(glm::cross(glm::vec3(0,1,0), dir));
	glm::vec3 up = glm::cross(dir, right);

	float yaw = -90.f;
	float pitch = 0.f;
} camera;


void drawScene(ShaderProgram* sp, GLFWwindow *window) {
	// time
	float currentFrame = glfwGetTime();
	deltaTime = currentFrame - lastFrame;
	lastFrame = currentFrame;

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	sp->use();

	glm::mat4 M = glm::mat4(1.f);
	glm::mat4 V = glm::lookAt(camera.pos, camera.pos + camera.target, camera.up);
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

	camera.yaw += xoffset;
	camera.pitch += yoffset;
	if (camera.pitch > 89.f) {
		camera.pitch = 89.f;
	}
	if (camera.pitch < -89.f) {
		camera.pitch = -89.f;
	}
	glm::vec3 direction;
	direction.x = cos(glm::radians(camera.yaw)) * cos(glm::radians(camera.pitch));
	direction.y = sin(glm::radians(camera.pitch));
	direction.z = sin(glm::radians(camera.yaw)) * cos(glm::radians(camera.pitch));
	camera.target = glm::normalize(direction);
}

void keyCallback(GLFWwindow *window, int key, int scancode, int act, int mod) {
	const float speed = 10.f * deltaTime;

	if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
		camera.pos += speed * camera.target;
	if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
		camera.pos -= speed * camera.target;
	if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
		camera.pos -= 
			glm::normalize(glm::cross(camera.target, camera.up)) * speed;
	if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
    	camera.pos += 
    		glm::normalize(glm::cross(camera.target, camera.up)) * speed;
#ifdef DEBUG
	printf("Pos: x: %f, y: %f, z: %f\n", 
			camera.pos.x, camera.pos.y, camera.pos.z);
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

	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
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

// converts Position<uint> into a reference to an element in the noise array
inline float& noiseAtPos(Position<uint> pos) {
	return noise.at(pos.x + pos.z * NOISE_W + pos.y * (NOISE_W * NOISE_W));
}
// same as above + interpolation for floats
inline float noiseAtPos(glm::vec3 pos) {
	// pivot
	glm::vec3 p0(uint(pos.x), uint(pos.y), uint(pos.z));

	// shifted vectors for sampling
	glm::vec3 px = p0 + glm::vec3(1,0,0);
	glm::vec3 py = p0 + glm::vec3(0,1,0);
	glm::vec3 pz = p0 + glm::vec3(0,0,1);

	// samples
	float sx = noiseAtPos(Position<uint>{uint(px.x), uint(px.y), uint(px.z)});
	float sy = noiseAtPos(Position<uint>{uint(py.x), uint(py.y), uint(py.z)});
	float sz = noiseAtPos(Position<uint>{uint(pz.x), uint(pz.y), uint(pz.z)});

	// distance 
	float dx = glm::distance(px, pos);
	float dy = glm::distance(py, pos);
	float dz = glm::distance(pz, pos);

	// interpolated
	return (dx*sx + dy*sy + dz*sz) / (dx + dy + dz);
}

glm::vec3 interpolateVertex(Position<uint> &v1, float val1, Position<uint> &v2, float val2) {
	// finding a 0 
	glm::vec3 ret;
	// coefficient
	float mu = (THRESHOLD - val1) / (val2 - val1);
	ret.x = mu*((int)v2.x - (int)v1.x) + v1.x;
	ret.y = mu*((int)v2.y - (int)v1.y) + v1.y;
	ret.z = mu*((int)v2.z - (int)v1.z) + v1.z;
	return ret;
}

glm::vec3 calculateNormal(glm::vec3 v1) {
	// normal calculation - 6 samples
	// ignore edges of the map
	if (v1.x <= 0.f || v1.y <= 0.f || v1.z <= 0.f ||
				v1.x >= MAP_W || v1.y >= MAP_H || v1.z >= MAP_W) {
		return glm::vec3(0,0,0);
	}
	// derivatives to figure out normal by using cross product
	// v1 - position of THE VOXEL in world space based on the noise
	float dx = noiseAtPos(v1 + glm::vec3(1, 0, 0)) - noiseAtPos(v1 + glm::vec3(-1, 0, 0));
	float dy = noiseAtPos(v1 + glm::vec3(0,-1, 0)) - noiseAtPos(v1 + glm::vec3( 0, 1, 0));
	float dz = noiseAtPos(v1 + glm::vec3(0, 0, 1)) - noiseAtPos(v1 + glm::vec3( 0, 0,-1));
	return -glm::normalize(glm::vec3(dx, dy, dz));
}

void genMesh() {
	constexpr uint step = 12; // 6 vertices + 6 normals // will be used for an array later on
	for (uint i = 0; i < (MAP_W * MAP_W * MAP_H) * step; i += step) {
		uint frac = i / step;
		auto xz = frac % (MAP_W * MAP_W); // index on the <x,z> plane
		Position<uint> pos = {
			.x = xz % MAP_W,
			.y = frac / (MAP_W * MAP_W),
			.z = xz / MAP_W
		};

		// cube in binary: (v7, v6, v5, v4, v3, v2, v1, v0)
		byte cube = 0x0;
		// offset in binary: (z, y, x)
		for (byte offset = 0b000; offset <= 0b111; offset++) {
			Position<uint> nPos = {
				.x = pos.x + (offset & 0b001),
				.y = pos.y + ((offset & 0b010) >> 1),
				.z = pos.z + ((offset & 0b100) >> 2),
			};
			float noiseVal = noiseAtPos(nPos);
			if (noiseVal > THRESHOLD) {
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

			// position of the first vertex on the edge
			auto v1 = Position<uint>{
				.x = pos.x + ((verts[0] & 0b001) >> 0),
				.y = pos.y + ((verts[0] & 0b010) >> 1),
				.z = pos.z + ((verts[0] & 0b100) >> 2)
			};

			// position of the second vertex on the edge
			auto v2 = Position<uint>{
				.x = pos.x + ((verts[1] & 0b001) >> 0),
				.y = pos.y + ((verts[1] & 0b010) >> 1),
				.z = pos.z + ((verts[1] & 0b100) >> 2)
			};
			// position of the interpolated vertex (based on the noise value at v1 and v2)
			glm::vec3 vertex = interpolateVertex(v1, noiseAtPos(v1), v2, noiseAtPos(v2));
			verts3D.push_back(vertex);
			// interpolating the normal
			verts3D.push_back(calculateNormal(vertex));
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
				noiseAtPos(pos) = val;
			}
		}
#ifdef DEBUG
		// outputnoise directory must exist
		constexpr uint size = NOISE_W * NOISE_W;
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

// debugging map
void fillNoiseTest() {
	for (uint y = 0; y < NOISE_H; y++) {
		array<byte, NOISE_W * NOISE_W> picture;
		for (uint z = 0; z <= MAP_W; z++) {
			for (uint x = 0; x <= MAP_W; x++) {
				auto pos = Position<uint>{x, y, z};
				noiseAtPos(pos) = -1.f;
			}
		}
		auto pos = Position<uint>{0, 0, 0};
		noiseAtPos(pos) = -1.f;
#ifdef DEBUG
		// outputnoise directory must exist
		constexpr uint size = NOISE_W * NOISE_W;
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
	fill3DNoise();
	genMesh();

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
