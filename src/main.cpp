#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <cstdio>
#include <cstdlib>
#include <cstddef>

#include <vector>
#include <array>
#include <memory>
#include <unordered_map>

#include "SimplexNoise/SimplexNoise.h"
#include "lodepng/lodepng.h"
#include "shaderprogram.hpp"
#include "marchingcubes.hpp"
#include "log.hpp"

using byte = unsigned char;
using uint = unsigned;

// unique_ptr alias
template<typename... T>
using Uniq_Ptr = std::unique_ptr<T...>;

// temporary vector to be removed
using std::array;
using std::vector;

constexpr const char* WINDOW_TITLE = "GLTRY";
constexpr uint WINDOW_W = 800;
constexpr uint WINDOW_H = 600;
constexpr float ASPECT_RATIO = float(WINDOW_W) / WINDOW_H;
constexpr float FOV = 90.f;

// clipping distance
constexpr float Z_NEAR = 0.1f;
constexpr float Z_FAR = 256.f;

// size of one chunk
constexpr int CHNK_SIZE = 32;
constexpr int MAP_W = CHNK_SIZE;
constexpr int MAP_H = CHNK_SIZE;


// size of the noise array 
constexpr uint NOISE_W = MAP_W + 2;
constexpr uint NOISE_H = MAP_H + 2;

// iso
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

// world space coordinates
using Pos3i = Position<int>;

struct Vertex {
	glm::vec3 pos;
	glm::vec3 norm;
};

GLuint vao;
//GLuint vbo;

float deltaTime = 0.f;
float lastFrame = 0.f;

// mouse position
Position<float> mouseLast;

using ChnkNoise = array<float, NOISE_W * NOISE_W * NOISE_H>;
//ChnkNoise noise;

// vertex + normal data (6 + 6)
// {pos{x, y, z}, norm{x, y, z}} for every single vertex
//vector<Vertex> mesh;

// ------------------------- chunks ----------------------
struct Chunk {
	// mesh vertex size
	size_t indices = 0;
	// mesh
	GLuint vbo;
	ChnkNoise noise;
};

// used for hashing
using ChnkCoord = uint64_t;
std::unordered_map<ChnkCoord, Chunk> chunkMap;


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

void drawChunk(Chunk &chnk) {
	//LOG_INFO("VBO: %u", chnk.vbo);
	glBindVertexArray(vao);
	//glBindVertexArray(chnk.vbo);
	glBindBuffer(GL_ARRAY_BUFFER, chnk.vbo);
	glDrawArrays(GL_TRIANGLES, 0, chnk.indices);
}

void errCallback(int error, const char *description) {
	LOG_PRINT(description);
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
	LOG_INFO("xyz: %.1f %.1f %.1f", camera.pos.x, camera.pos.y, camera.pos.z);
}

// DEBUGGING INFO
void GLAPIENTRY MessageCallback(GLenum source, GLenum type, GLuint id,
		GLenum severity, GLsizei length,
		const GLchar *message, const void *userParam) {
	LOG_WARN("GL CALLBACK: %s type = 0x%x, severity = 0x%x, message = %s\n",
			(type == GL_DEBUG_TYPE_ERROR ? "** GL ERROR **" : ""), 
			type, severity, message);
}

void windowResizeCallback(GLFWwindow* window, int width, int height) {
    if (height == 0)
	return;
    glViewport(0, 0, width, height);
}

ShaderProgram* initProgram(GLFWwindow *window) {
	LOG_INFO("Initialising openGL...");
#ifdef DEBUG
	//glEnable(GL_DEBUG_OUTPUT);
	//glDebugMessageCallback(MessageCallback, 0);
#endif

	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);

	//glEnable(GL_CULL_FACE);
	//glCullFace(GL_BACK);

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
	return sp;
}

int pFail(const char *message) {
	LOG_PRINT(message);
	return -1;
}

// position must be 32bit intiger otherwise crashiento
// returns 64bit uint containing position of the chunk 
ChnkCoord chunkHasher(Pos3i &pos) {
	constexpr int max_chunks = static_cast<int>(INT32_MAX / CHNK_SIZE);
	static int power = log2(max_chunks);
	auto x = (pos.x + CHNK_SIZE) / CHNK_SIZE;
	auto y = (pos.y + CHNK_SIZE) / CHNK_SIZE;
	auto z = (pos.z + CHNK_SIZE) / CHNK_SIZE;
	return x | z << power | y << power * 2;
};

// converts Position<uint> into a reference to an element in the noise array
inline float& noiseAtPos(Pos3i pos) {
	auto chnkHash = chunkHasher(pos);
	auto &noise = chunkMap[chnkHash].noise;
	auto x = glm::abs(pos.x % CHNK_SIZE);
	auto y = glm::abs(pos.y % CHNK_SIZE);
	auto z = glm::abs(pos.z % CHNK_SIZE);
	return noise.at(x + z * NOISE_W + y * (NOISE_W * NOISE_W));
}

// same as above + interpolation for floats
inline float noiseAtPos(glm::vec3 pos) {
	// pivot
	glm::vec3 p0(int(pos.x), int(pos.y), int(pos.z));

	// shifted vectors for sampling
	glm::vec3 px = p0 + glm::vec3(1,0,0);
	glm::vec3 py = p0 + glm::vec3(0,1,0);
	glm::vec3 pz = p0 + glm::vec3(0,0,1);

	// samples
	float sx = noiseAtPos(Pos3i{int(px.x), int(px.y), int(px.z)});
	float sy = noiseAtPos(Pos3i{int(py.x), int(py.y), int(py.z)});
	float sz = noiseAtPos(Pos3i{int(pz.x), int(pz.y), int(pz.z)});

	// distance 
	float dx = glm::distance(px, pos);
	float dy = glm::distance(py, pos);
	float dz = glm::distance(pz, pos);

	// interpolated
	return (dx*sx + dy*sy + dz*sz) / (dx + dy + dz);
	//return noiseAtPos(Position<uint>{uint(pos.x), uint(pos.y), uint(pos.z)});
}

glm::vec3 interpolateVertex(Pos3i &v1, float val1, Pos3i &v2, float val2) {
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
	return -glm::normalize(glm::vec3(dx, -dy, dz));
	//return glm::vec3(0,0,1);
}

vector<Vertex> genMesh(Pos3i &curr) {
	LOG_INFO("Generating mesh...");
	std::vector<Vertex> mesh;
	constexpr uint step = 12; // 6 vertices + 6 normals // will be used for an array later on
	for (int i = 0; i < (MAP_W * MAP_W * MAP_H) * step; i += step) {
		int frac = i / step;
		auto xz = frac % (MAP_W * MAP_W); // index on the <x,z> plane
		Pos3i pos = {
			.x = curr.x + xz % MAP_W,
			.y = curr.y + frac / (MAP_W * MAP_W),
			.z = curr.z + xz / MAP_W
		};

		// cube in binary: (v7, v6, v5, v4, v3, v2, v1, v0)
		byte cube = 0x0;
		// offset in binary: (z, y, x)
		for (byte offset = 0b000; offset <= 0b111; offset++) {
			Pos3i nPos = {
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
			auto v1 = Pos3i{
				.x = pos.x + int((verts[0] & 0b001) >> 0),
				.y = pos.y + int((verts[0] & 0b010) >> 1),
				.z = pos.z + int((verts[0] & 0b100) >> 2)
			};

			// position of the second vertex on the edge
			auto v2 = Pos3i{
				.x = pos.x + int((verts[1] & 0b001) >> 0),
				.y = pos.y + int((verts[1] & 0b010) >> 1),
				.z = pos.z + int((verts[1] & 0b100) >> 2)
			};

			Vertex vertex;
			// position of the interpolated vertex (based on the noise value at v1 and v2)
			vertex.pos = interpolateVertex(v1, noiseAtPos(v1), v2, noiseAtPos(v2));
			// interpolating the normal
			vertex.norm = calculateNormal(vertex.pos);
			mesh.push_back(vertex);
		}
	}
	return mesh;
}

void genChunk(Pos3i &currPos, ShaderProgram *sp) {
	auto hash = chunkHasher(currPos);
	auto mesh = genMesh(currPos);
	LOG_INFO("Mesh size: %lu", mesh.size());

	GLuint vbo;
	glGenBuffers(1, &vbo);
	LOG_INFO("vbo: %i", vbo);
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBufferData(GL_ARRAY_BUFFER, mesh.size() * sizeof(Vertex),
         	(void*)mesh.data(), GL_STATIC_DRAW);

	uint count = 3;
	glVertexAttribPointer(sp->a("vertex"), count, GL_FLOAT, GL_FALSE,
        	sizeof(Vertex), (GLvoid*)offsetof(Vertex, pos));
	glEnableVertexAttribArray(sp->a("vertex"));

	glVertexAttribPointer(sp->a("normal"), count, GL_FLOAT, GL_FALSE,
        	sizeof(Vertex), (GLvoid*)offsetof(Vertex, norm));
	glEnableVertexAttribArray(sp->a("normal"));

	chunkMap[hash].indices = mesh.size();
	chunkMap[hash].vbo = vbo;
}
// -------------------------------------------------------

void saveNoisePNG(const char* fname, ChnkNoise &noise) {
	LOG_INFO("Saving noise pictures...");
	array<byte, NOISE_W * NOISE_W> picture;
	float max = -99999.f;
	for (uint i = 0; i < noise.size(); i++) {
		picture.at(i) = (1 + noise.at(i)) * 128.f;
	}
	lodepng_encode_file(fname, picture.data(), 
			NOISE_W, NOISE_W, LCT_GREY, 8);
}
ChnkNoise fillNoise(Pos3i &curr) {
	LOG_INFO("Generating noise...");
	ChnkNoise noise;
	SimplexNoise noiseGen(0.1f, 1.f, 2.f, 0.5f);
	constexpr uint octaves = 9;
	// noise generation
	for (int y = 0; y < NOISE_H; y++) {
		//array<byte, NOISE_W * NOISE_W> picture;
		for (int z = 0; z <= MAP_W; z++) {
			for (int x = 0; x <= MAP_W; x++) {
				auto val = noiseGen.fractal(octaves, x, y, z);
				auto pos = Pos3i{x + curr.x, y + curr.y, z + curr.z};
				noiseAtPos(pos) = val;
				noise.at(x + z * NOISE_W + y * (NOISE_W * NOISE_W)) = val;
			}
		}
/* #ifdef DEBUG
		// outputnoise directory must exist
		constexpr uint size = NOISE_W * NOISE_W;
		array<float, size> noiseSlice;
		std::copy(noise.begin() + size * y,
				noise.begin() + size * (y + 1),
				noiseSlice.begin());
		char str[21 + 48 / 8];
		sprintf(str, "outputnoise/noise-%u.png", y);
		saveNoisePNG(str, noiseSlice);
#endif */
	}
	return noise;
}

ChnkNoise fillall(Pos3i &curr) {
	LOG_INFO("Generating noise...");
	ChnkNoise noise;
	// noise generation
	for (int y = 0; y < NOISE_H; y++) {
		//array<byte, NOISE_W * NOISE_W> picture;
		for (int z = 0; z <= MAP_W; z++) {
			for (int x = 0; x <= MAP_W; x++) {
				auto val = 1.f;
				//auto pos = Pos3i{x + curr.x, y + curr.y, z + curr.z};
				//noiseAtPos(pos) = val;
				noise.at(x + z * NOISE_W + y * (NOISE_W * NOISE_W)) = val;
			}
		}
/* #ifdef DEBUG
		// outputnoise directory must exist
		constexpr uint size = NOISE_W * NOISE_W;
		array<float, size> noiseSlice;
		std::copy(noise.begin() + size * y,
				noise.begin() + size * (y + 1),
				noiseSlice.begin());
		char str[21 + 48 / 8];
		sprintf(str, "outputnoise/noise-%u.png", y);
		saveNoisePNG(str, noiseSlice);
#endif */
	}
	return noise;
}

void renderChunks(GLFWwindow *window, ShaderProgram *sp) {
	float currentFrame = glfwGetTime();
	deltaTime = currentFrame - lastFrame;
	lastFrame = currentFrame;

	sp->use();

	glm::mat4 M = glm::mat4(1.f);
	glm::mat4 V = glm::lookAt(camera.pos, camera.pos + camera.target, camera.up);
	glm::mat4 P = glm::perspective(glm::radians(FOV), ASPECT_RATIO, Z_NEAR, Z_FAR);

	glUniformMatrix4fv(sp->u("M"), 1, false, glm::value_ptr(M));
	glUniformMatrix4fv(sp->u("V"), 1, false, glm::value_ptr(V));
	glUniformMatrix4fv(sp->u("P"), 1, false, glm::value_ptr(P));

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	constexpr int CHNK_DIST = 1;

	uint i = 0;
	for (int x = -CHNK_DIST; x <= CHNK_DIST; x++) {
	for (int y = -CHNK_DIST; y <= CHNK_DIST; y++) {
	for (int z = -CHNK_DIST; z <= CHNK_DIST; z++) {
		i++;
		Pos3i currentPos = {
			.x = int(camera.pos.x) + x * int(CHNK_DIST),
			.y = int(camera.pos.y) + y * int(CHNK_DIST),
			.z = int(camera.pos.z) + z * int(CHNK_DIST)
		};
		auto hash = chunkHasher(currentPos);
		Chunk &chnk = chunkMap[hash];
		if (chnk.indices == 0) {
			LOG_INFO("Generating chunk (%u / %.0f)...", i, pow(CHNK_DIST * 2 + 1, 3));
			LOG_INFO("Hash: %lu", hash);
			fillNoise(currentPos);
			genChunk(currentPos, sp);
		}
		drawChunk(chnk);
	}
	}
	}
	glfwSwapBuffers(window);
}

int main() {
	glfwSetErrorCallback(errCallback);
	if (!glfwInit()) {
		glfwTerminate();
		return pFail("Cannot initialise GLFW.");
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

	if (!window.get()) return pFail("Cannot create window.");

	glfwMakeContextCurrent(window.get());
	glfwSwapInterval(1);
	glViewport(0, 0, WINDOW_W, WINDOW_H);

	if (glewInit() != GLEW_OK) return pFail("Cannot initialise GLEW.");

	auto shaderDestroyer = [&](ShaderProgram* sp) {
		array<GLuint, 1> buffers({vao});
		glDeleteBuffers(buffers.size(), buffers.data());
		delete sp;
	};
	// pointer to a ShaderProgram
	Uniq_Ptr<ShaderProgram, decltype(shaderDestroyer)> sp(
			initProgram(window.get()),
			shaderDestroyer);
	
	// main loop
	while (!glfwWindowShouldClose(window.get())) {
    	glfwPollEvents();
    	renderChunks(window.get(), sp.get());
	}
}
