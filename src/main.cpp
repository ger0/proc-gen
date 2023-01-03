#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <cstdio>
#include <cstdlib>
#include <cstddef>

#include <thread>
#include <mutex>

#include <vector>
#include <array>
#include <memory>
#include <unordered_map>

#include <omp.h>

#include "imgui/imgui.h"
#include "imgui/imgui_impl_glfw.h"
#include "imgui/imgui_impl_opengl3.h"

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
constexpr uint NOISE_W = MAP_W + 3;
constexpr uint NOISE_H = MAP_H + 3;

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

// [TODO]: remove later
Pos3i currChnkPos;

GLuint vao;
std::mutex glCriticalSection;

float deltaTime = 0.f;
float lastFrame = 0.f;

// mouse position
Position<float> mouseLast;

struct Vertex {
	glm::vec3 pos;
	glm::vec3 norm;
};

using ChnkNoise = array<float, NOISE_W * NOISE_W * NOISE_H>;

// ------------------------- chunks ----------------------
struct Chunk {
	// mesh vertex size
	size_t indices = 0;
	// when a thread generating a new chunk dies the mesh will be submitted
	// to vbo by the main thread and the vector removed
	GLuint vbo;
	vector<Vertex> mesh;

	ChnkNoise noise;
	std::mutex isUsed;
	bool generated = false;
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

void drawChunk(Chunk &chnk, ShaderProgram* sp) {
	glBindVertexArray(vao);
	glBindBuffer(GL_ARRAY_BUFFER, chnk.vbo);

	uint count = 3;
	// passing position vector to vao
	glVertexAttribPointer(sp->a("vertex"), count, GL_FLOAT, GL_FALSE,
        	sizeof(Vertex), (GLvoid*)offsetof(Vertex, pos));
	glEnableVertexAttribArray(sp->a("vertex"));

	// passing normal vector to vao
	glVertexAttribPointer(sp->a("normal"), count, GL_FLOAT, GL_FALSE,
        	sizeof(Vertex), (GLvoid*)offsetof(Vertex, norm));
	glEnableVertexAttribArray(sp->a("normal"));

	glDrawArrays(GL_TRIANGLES, 0, chnk.indices);

	glBindVertexArray(0);
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
	const float speed = 20.f * deltaTime;
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

	currChnkPos.x = glm::floor(camera.pos.x / CHNK_SIZE);
	currChnkPos.y = glm::floor(camera.pos.y / CHNK_SIZE);
	currChnkPos.z = glm::floor(camera.pos.z / CHNK_SIZE);
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
	glEnable(GL_DEBUG_OUTPUT);
	glDebugMessageCallback(MessageCallback, 0);
#endif

	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);

	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);

	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glClearColor(0.45, 0.716, 0.914, 1.f);

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

// input: worldspace coordinates, output: chunk hash
// position must be 16bit integer otherwise crashiento
// returns 64bit uint containing position of the chunk
ChnkCoord chunkHasher(Pos3i &pos) {
	ChnkCoord x = glm::floor(pos.x / CHNK_SIZE);
	ChnkCoord y = glm::floor(pos.y / CHNK_SIZE);
	ChnkCoord z = glm::floor(pos.z / CHNK_SIZE);
	x = x << 0;
	z = z << 16;
	y = y << 32;
	return x + y + z;
};

// converts Position<uint> into a reference to an element in the noise array
inline float& noiseAtPos(Pos3i pos, ChnkNoise &noise) {
	//auto chnkHash = chunkHasher(pos);
	//auto &chunk = chunkMap[chnkHash];
	auto x = pos.x + 1;
	auto y = pos.y + 1;
	auto z = pos.z + 1;
	return noise.at(x + z * NOISE_W + y * (NOISE_W * NOISE_W));
}

// same as above + interpolation for floats
inline float noiseAtPos(glm::vec3 pos, ChnkNoise &noise) {
	// no interpolation (poormans normal)
	return noiseAtPos(Pos3i{int(pos.x), int(pos.y), int(pos.z)}, noise);
	/*
	// interpolated
	// pivot
	glm::vec3 p0(int(pos.x), int(pos.y), int(pos.z));

	// shifted vectors for sampling
	glm::vec3 px = p0 + glm::vec3(1,0,0);
	glm::vec3 py = p0 + glm::vec3(0,1,0);
	glm::vec3 pz = p0 + glm::vec3(0,0,1);

	// samples
	float sx = noiseAtPos(Pos3i{int(px.x), int(px.y), int(px.z)}, noise );
	float sy = noiseAtPos(Pos3i{int(py.x), int(py.y), int(py.z)}, noise );
	float sz = noiseAtPos(Pos3i{int(pz.x), int(pz.y), int(pz.z)}, noise );

	// distance 
	float dx = glm::distance(px, pos);
	float dy = glm::distance(py, pos);
	float dz = glm::distance(pz, pos);

	return (dx*sx + dy*sy + dz*sz) / (dx + dy + dz); */
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

glm::vec3 calculateNormal(glm::vec3 v1, ChnkNoise &noise) {
	// derivatives to figure out normal by using cross product
	// v1 - position of THE VOXEL in world space based on the noise
	float dx = noiseAtPos(v1 + glm::vec3(1, 0, 0), noise) - noiseAtPos(v1 + glm::vec3(-1, 0, 0), noise);
	float dy = noiseAtPos(v1 + glm::vec3(0,-1, 0), noise) - noiseAtPos(v1 + glm::vec3( 0, 1, 0), noise);
	float dz = noiseAtPos(v1 + glm::vec3(0, 0, 1), noise) - noiseAtPos(v1 + glm::vec3( 0, 0,-1), noise);
	return -glm::normalize(glm::vec3(dx, -dy, dz));
	//return glm::vec3(0,0,1);
}

vector<Vertex> genMesh(Pos3i &curr, ChnkNoise &noise) {
	// cubes containing vertex information for entire chunk
	array<vector<Vertex>, MAP_W * MAP_W * MAP_H> cubes;

	auto timer = glfwGetTime();
#pragma omp parallel for
	for (int iter = 0; iter < (MAP_W * MAP_W * MAP_H); iter++) {
		auto xz = iter % (MAP_W * MAP_W); // index on the <x,z> plane
		Pos3i chunkPos = {
			.x = xz % MAP_W,
			.y = iter / (MAP_W * MAP_W),
			.z = xz / MAP_W
		};
		// cube in binary: (v7, v6, v5, v4, v3, v2, v1, v0)
		byte cube = 0x0;
		// offset in binary: (z, y, x)
		for (byte offset = 0b000; offset <= 0b111; offset++) {
			Pos3i nPos = {
				.x = chunkPos.x + (offset & 0b001),
				.y = chunkPos.y + ((offset & 0b010) >> 1),
				.z = chunkPos.z + ((offset & 0b100) >> 2),
			};
			float noiseVal = noiseAtPos(nPos, noise);
			if (noiseVal > THRESHOLD) {
				cube |= (1 << offset);
			}
		}
		auto mask = march_cubes::EdgeMasks.at(cube);

		// array of triplets of the cube edges, each triplet = 1 triangle
		auto edges = march_cubes::TriangleTable.at(cube);
		
		for (uint& edge : edges) {
			// check if theres no more edges left
			if (edge == march_cubes::x) {
				break;
			}

			auto verts = march_cubes::EdgeVertexIndices.at(edge);

			// position of the first vertex on the edge
			auto v1 = Pos3i{
				.x = chunkPos.x + int((verts[0] & 0b001) >> 0),
				.y = chunkPos.y + int((verts[0] & 0b010) >> 1),
				.z = chunkPos.z + int((verts[0] & 0b100) >> 2)
			};

			// position of the second vertex on the edge
			auto v2 = Pos3i{
				.x = chunkPos.x + int((verts[1] & 0b001) >> 0),
				.y = chunkPos.y + int((verts[1] & 0b010) >> 1),
				.z = chunkPos.z + int((verts[1] & 0b100) >> 2)
			};

			Vertex vertex;
			// position of the interpolated vertex (based on the noise value at v1 and v2)
			vertex.pos = interpolateVertex(v1, noiseAtPos(v1, noise), v2, noiseAtPos(v2, noise));
			// interpolating the normal
			vertex.norm = calculateNormal(vertex.pos, noise);
			// adding worldspace shift
			vertex.pos.x += curr.x;
			vertex.pos.y += curr.y;
			vertex.pos.z += curr.z;

			cubes.at(iter).push_back(vertex);
		}
	}
	vector<Vertex> mesh;
	for (auto &cube : cubes) {
		mesh.insert(mesh.end(), cube.begin(), cube.end());
	}
	// copy elision so it shouldnt return by value
	//LOG_INFO("size: %lu", mesh.size());
	return mesh;
}

void genChunk(Pos3i &currPos, Chunk &chnk) {
	auto mesh = genMesh(currPos, chnk.noise);
	chnk.mesh = mesh;
	chnk.indices = mesh.size();
	chnk.generated = true;
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


ChnkNoise genNoise(Pos3i &curr) {
	ChnkNoise noise;
	static SimplexNoise noiseGen(0.01f, 1.f, 2.f, 0.5f);
	constexpr uint octaves = 8;
	// noise generation
#pragma omp parallel for
	for (int y = 0; y < NOISE_H; y++) {
		for (int z = 0; z < NOISE_W; z++) {
			for (int x = 0; x < NOISE_W; x++) {
				auto pos = Pos3i{x + curr.x, y + curr.y, z + curr.z};
				auto val = noiseGen.fractal(octaves, pos.x -1, pos.y -1, pos.z -1);
				val -= pos.y / 32.f;
				//noiseAtPos(pos) = val;
				noise.at(x + z * NOISE_W + y * (NOISE_W * NOISE_W)) = val;
			}
		}
	}
	return noise;
}

ChnkNoise fillall(Pos3i &curr) {
	LOG_INFO("Generating noise...");
	ChnkNoise noise;
	// noise generation
	for (int y = 0; y < NOISE_H; y++) {
		for (int z = 0; z <= MAP_W; z++) {
			for (int x = 0; x <= MAP_W; x++) {
				auto val = 1.f;
				noise.at(x + z * NOISE_W + y * (NOISE_W * NOISE_W)) = val;
			}
		}
	}
	return noise;
}

void enqChunk(Pos3i cPos, Chunk *chnk) {
	chnk->noise = genNoise(cPos);
	genChunk(cPos, *chnk);
	chnk->isUsed.unlock();
}

void renderChunks(GLFWwindow *window, ShaderProgram *sp) {
	sp->use();

	glm::mat4 M = glm::mat4(1.f);
	glm::mat4 V = glm::lookAt(camera.pos, camera.pos + camera.target, camera.up);
	glm::mat4 P = glm::perspective(glm::radians(FOV), ASPECT_RATIO, Z_NEAR, Z_FAR);

	glUniformMatrix4fv(sp->u("M"), 1, false, glm::value_ptr(M));
	glUniformMatrix4fv(sp->u("V"), 1, false, glm::value_ptr(V));
	glUniformMatrix4fv(sp->u("P"), 1, false, glm::value_ptr(P));


	constexpr int CHNK_DIST = 3;

	ImGui::Begin("DEBUG INFO");
	ImGui::Text("x: %0.1f\n y: %0.1f\n z: %0.1f\n", camera.pos.x, camera.pos.y, camera.pos.z);
	ImGui::Text("[%i %i %i]\n\n", currChnkPos.x, currChnkPos.y, currChnkPos.z);

	for (int x = -CHNK_DIST; x <= CHNK_DIST; x++) {
	for (int y = -CHNK_DIST / 2; y <= CHNK_DIST / 2; y++) {
	for (int z = -CHNK_DIST; z <= CHNK_DIST; z++) {
		Pos3i cPos = {
			.x = (currChnkPos.x + x) * CHNK_SIZE,
			.y = (currChnkPos.y + y) * CHNK_SIZE,
			.z = (currChnkPos.z + z) * CHNK_SIZE
		};
		auto hash = chunkHasher(cPos);

// 		ImGui::Text("Hash[%i %i %i]:\n %lu\n", cPos.x / CHNK_SIZE, 
// 				cPos.y / CHNK_SIZE, cPos.z / CHNK_SIZE, hash); 

		Chunk &chnk = chunkMap[hash];


		if (chnk.isUsed.try_lock()) {
			if (chnk.generated) {
				// submitting the mesh to gpu
				if (chnk.mesh.size() > 0) {
					glBindBuffer(GL_ARRAY_BUFFER, chnk.vbo);
					glBufferData(GL_ARRAY_BUFFER, chnk.indices * sizeof(Vertex),
         					(void*)chnk.mesh.data(), GL_STATIC_DRAW);
         			chnk.mesh.clear();
         			// lag protection
					lastFrame = glfwGetTime();
				}
				drawChunk(chnk, sp);
				chnk.isUsed.unlock();
			} else {
				// vbo for the chunk
				GLuint vbo;
				glGenBuffers(1, &vbo);
				chnk.vbo = vbo;
				LOG_INFO("vbo: %i", vbo);

				LOG_INFO("Generating chunk [%i %i %i]...\nHash: %lu", cPos.x / CHNK_SIZE, 
					cPos.y / CHNK_SIZE, cPos.z / CHNK_SIZE, hash);
				std::thread thread(enqChunk, cPos, &chnk);
				thread.detach();
				//enqChunk(cPos, &chnk);
			}
		}
	}}}
	ImGui::End();
}

int main() {
	glfwSetErrorCallback(errCallback);
	if (!glfwInit()) {
		glfwTerminate();
		return pFail("Cannot initialise GLFW.");
	}

    const char* glsl_version = "#version 130";
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);

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

	// imgui init
	IMGUI_CHECKVERSION();
	ImGui::CreateContext();

	ImGuiIO& io = ImGui::GetIO();

	ImGui::StyleColorsClassic();
	ImGui_ImplGlfw_InitForOpenGL(window.get(), true);

	ImGui_ImplOpenGL3_Init(glsl_version);

	auto shaderDestroyer = [&](ShaderProgram* sp) {
		array<GLuint, 1> buffers({vao});
		glDeleteBuffers(buffers.size(), buffers.data());
		
		// freeing vbos on the gpu
		vector<GLuint> vbos;
		for (auto &chnk : chunkMap) {
			vbos.push_back(chnk.second.vbo);
		}
		glDeleteBuffers(vbos.size(), vbos.data());
		delete sp;
	};
	// pointer to a ShaderProgram
	Uniq_Ptr<ShaderProgram, decltype(shaderDestroyer)> sp(
			initProgram(window.get()),
			shaderDestroyer);

	// main loop
	while (!glfwWindowShouldClose(window.get())) {
		float currentFrame = glfwGetTime();
		deltaTime = currentFrame - lastFrame;
		lastFrame = currentFrame;

    	glfwPollEvents();
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		// imgui
    	ImGui_ImplOpenGL3_NewFrame();
    	ImGui_ImplGlfw_NewFrame();
    	ImGui::NewFrame();

    	renderChunks(window.get(), sp.get());

		ImGui::Render();
		ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

		glfwSwapBuffers(window.get());
	}
}
