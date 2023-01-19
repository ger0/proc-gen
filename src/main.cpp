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

#include <boost/asio/thread_pool.hpp>
#include <boost/asio/post.hpp>

#include "main.hpp"

#include "utils.hpp"
#include "imgui/imgui.h"
#include "imgui/imgui_impl_glfw.h"
#include "imgui/imgui_impl_opengl3.h"

#include "SimplexNoise/SimplexNoise.h"
#include "shaderprogram.hpp"
#include "marchingcubes.hpp"
#include "log.hpp"

#include "model.hpp"
#include "poisson_disk_sampling.h"

uint SEED;
GLuint m_depthTexture;

// shader settings
constexpr auto  SUN_DIRECTION = glm::vec3(0, 1, -0.12);
constexpr float FOG_DENSITY   = 0.0036;
constexpr float FOG_GRADIENT  = 8.0;

// generation parameters
constexpr float LVL_SAND = 2.f;
constexpr float LVL_GRASS= 16.f;
constexpr float LVL_DIRT = 48.f;
constexpr float LVL_ROCK = 96.f;
constexpr float LVL_SNOW = 128.f;

constexpr float MIN_SHARP = 0.35;
constexpr float MAX_SHARP = 0.50;

constexpr float MIN_MOUNT = 0.1;
constexpr float MAX_MOUNT = 2.5f;

constexpr float MIN_HEIGH = 0.1;
constexpr float MAX_HEIGH = 128.f;

constexpr float MIN_HUMID = 0.f;
constexpr float MAX_HUMID = 1.f;

constexpr float MIN_FORES = 0.f;
constexpr float MAX_FORES = 10.f;

// colours
constexpr glm::vec4 waterColor(0.15, 0.216, 0.614, 0.46);
constexpr glm::vec3 skyColor(0.45, 0.716, 0.914);

constexpr const char* WINDOW_TITLE = "Proc-Gen";
constexpr uint WINDOW_W = 1280;
constexpr uint WINDOW_H = 720;
constexpr float ASPECT_RATIO = float(WINDOW_W) / WINDOW_H;
constexpr float FOV = 90.f;

// size of one chunk
constexpr int CHK_SIZE = 32;
constexpr int MAP_W = CHK_SIZE;
constexpr int MAP_H = CHK_SIZE;

// size of the noise array 
constexpr uint NOISE_W = MAP_W + 3;
constexpr uint NOISE_H = MAP_H + 3;

// used for hashing
using ChkHash = uint64_t;

using ChkNoise = array<float, NOISE_W * NOISE_W * NOISE_H>;
// biome map (temporary // more biomes tbd)
using SmoothMap = array<float, NOISE_W * NOISE_W>;

// render distance
constexpr int RENDER_DIST = 8;
constexpr int RENDER_RADIUS = RENDER_DIST * RENDER_DIST;
constexpr float Z_NEAR = 0.1f;
constexpr float Z_FAR = (RENDER_DIST + 3) * CHK_SIZE;

// [TODO]: remove later
Pos3i currChkPos;
boost::asio::thread_pool pool(std::thread::hardware_concurrency());

GLuint vao;
//GLuint fbo;

// used to draw water on y = 0
GLuint waterVBO;

float deltaTime = 0.f;
float lastFrame = 0.f;

// mouse position
glm::vec2 mouseLast;

enum JobState {
	none,
	delegated,
	awaiting,
	generated
};

// critical section - thread synchronization
std::mutex mtx;

// ------------------------- chunks ----------------------
struct TreeSamples {
	using Samples = vector<array<float, 2>>;
	Samples vec;
	std::mutex mtx;
};

struct BiomeChunk {
	std::unordered_map<ChkHash, TreeSamples> treeMap;
	// noises used by 3D noise generator
	JobState state = none;
	std::mutex mtx;
	static constexpr int width = CHK_SIZE + 1;
	using biomeNoise = array<float, width * width>;

	biomeNoise height;
	biomeNoise mountain;
	biomeNoise sharp;
	biomeNoise humidity;
	biomeNoise forest;

	float& at(biomeNoise& noise, int x, int z) {
		return noise.at(x + z * width);
	}

	// noise sampling v[4]:
	// 0 -- 1
	// |    |
	// 2 -- 3
	array<float, 4> sample(biomeNoise &noise, Pos3i &pos) {
		return array<float,4> {
			at(noise, pos.x % (uint)CHK_SIZE + 0, pos.z % (uint)CHK_SIZE + 0),
			at(noise, pos.x % (uint)CHK_SIZE + 1, pos.z % (uint)CHK_SIZE + 0),
			at(noise, pos.x % (uint)CHK_SIZE + 0, pos.z % (uint)CHK_SIZE + 1),
			at(noise, pos.x % (uint)CHK_SIZE + 1, pos.z % (uint)CHK_SIZE + 1),
		};
	}
};

struct TreeModel {
	GLuint vbo;
	uint indices;
};

void meshPush(Mesh &mesh, Mesh &val) {
    mesh.insert(mesh.end(), val.begin(), val.end());
}

struct BiomeSamples {
	array<float, 4> humidity;
	array<float, 4> mountain;
	array<float, 4> sharp;
	array<float, 4> height;
	array<float, 4> forest;
	
	static float avg(array<float, 4> &smp) {
		return (smp[0] + smp[1] + smp[2] + smp[3]) / 4.f;
	}
};

struct Tree {
	glm::vec3 position;
	float rotation;
	float skew;
	TreeModel &model;
};

struct Chunk {
	static constexpr uint MAX_MESH_VERTS = 12 * CHK_SIZE * CHK_SIZE * CHK_SIZE * 3;
	// chunkspace pos
	Pos3i pos;
	BiomeSamples biome;
	TreeSamples  *treeMap;

	// how many faces in the mesh to render
	size_t indices = 0;

	// after a thread generates the mesh, the main thread will submit it to GPU
	GLuint vbo;
	Mesh mesh;
	vector<Tree> trees;

	// 3D noise used for terrain generation
	ChkNoise *noise;

	// status
	std::atomic<JobState> state = ATOMIC_VAR_INIT(JobState::none);
	void clear() {
		mesh.clear();
		mesh.shrink_to_fit();
		delete noise;
	}
	ChkHash getHash() {
		return hasher(pos);
	}
	// input: worldspace coordinates, output: chunk hash
	// position must be 16bit integer otherwise crashiento
	// returns 64bit uint containing position of the chunk
	static ChkHash hasher(Pos3i &pos) {
		ChkHash x = glm::floor(pos.x / CHK_SIZE);
		ChkHash y = glm::floor(pos.y / CHK_SIZE);
		ChkHash z = glm::floor(pos.z / CHK_SIZE);
		x = x << 0;
		z = z << 21;
		y = y << 42;
		return x + y + z;
	};
};

// storage for all chunks
std::unordered_map<ChkHash, Chunk> chunkMap;
// storage for biomes, 1 cell = 1 chunk
std::unordered_map<ChkHash, BiomeChunk> biomeMap;
// pre-generated trees
array<TreeModel, 256> treemodels;

struct Camera {
	glm::vec3 pos = glm::vec3(0.f, 3.f, 3.f);
	glm::vec3 target = glm::vec3(0.f, 3.f, 0.f);
	glm::vec3 dir = glm::normalize(pos - target);

	const float speed = 0.05f;

	bool boost = false;	

	glm::vec3 right = glm::normalize(glm::cross(glm::vec3(0,1,0), dir));
	glm::vec3 up = glm::cross(dir, right);

	float yaw = -90.f;
	float pitch = 0.f;
} camera;

void drawVBO(GLuint vbo, uint indices, ShaderProgram* sp, glm::mat4 M = glm::mat4(1.f)) {
	sp->use();
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, m_depthTexture);

	glm::mat4 V = glm::lookAt(camera.pos, camera.pos + camera.target, camera.up);
	glm::mat4 P = glm::perspective(glm::radians(FOV), ASPECT_RATIO, Z_NEAR, Z_FAR);

	glUniformMatrix4fv(sp->u("M"), 1, false, glm::value_ptr(M));
	glUniform3fvARB(sp->u("cameraPos"), 1, glm::value_ptr(camera.pos));
	glUniformMatrix4fv(sp->u("V"), 1, false, glm::value_ptr(V));
	glUniformMatrix4fv(sp->u("P"), 1, false, glm::value_ptr(P));

	glBindVertexArray(vao);
	glBindBuffer(GL_ARRAY_BUFFER, vbo);

	uint count = 3;
	// passing position vector to vao
	glVertexAttribPointer(sp->a("vertex"), count, GL_FLOAT, GL_FALSE,
        	sizeof(Vertex), (GLvoid*)offsetof(Vertex, pos));
	glEnableVertexAttribArray(sp->a("vertex"));

	// passing normal vector to vao
	glVertexAttribPointer(sp->a("normal"), count, GL_FLOAT, GL_FALSE,
        	sizeof(Vertex), (GLvoid*)offsetof(Vertex, norm));
	glEnableVertexAttribArray(sp->a("normal"));

	glVertexAttribPointer(sp->a("color"), count + 1, GL_FLOAT, GL_FALSE,
        	sizeof(Vertex), (GLvoid*)offsetof(Vertex, color));
	glEnableVertexAttribArray(sp->a("color"));

	glDrawArrays(GL_TRIANGLES, 0, indices);

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
	float speed = 20.f * deltaTime * (camera.boost ? 4.f : 1.f);
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
	if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS) {
		camera.boost = true;
	}
	if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_RELEASE) {
		camera.boost = false;
	}

	currChkPos.x = glm::floor(camera.pos.x / CHK_SIZE);
	currChkPos.y = glm::floor(camera.pos.y / CHK_SIZE);
	currChkPos.z = glm::floor(camera.pos.z / CHK_SIZE);
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

void initProgram(GLFWwindow *window) {
	LOG_INFO("Initialising openGL...");
#ifdef DEBUG
	//glEnable(GL_DEBUG_OUTPUT);
	//glDebugMessageCallback(MessageCallback, 0);
#endif

	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);

	glFrontFace(GL_CW);
	glEnable(GL_BLEND);  

	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glClearColor(skyColor.r, skyColor.g, skyColor.b, 1.f);

	glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
	glfwSetCursorPosCallback(window, mouseCallback);
	glfwSetKeyCallback(window, keyCallback);

	glfwSetWindowSizeCallback(window, windowResizeCallback);


	// depth texture
	glGenTextures(1, &m_depthTexture);

	// generating buffers
	glGenVertexArrays(1, &vao);
	//glBindVertexArray(vao);
	//glGenFramebuffers(1, &fbo);
	//glBindFramebuffer(GL_FRAMEBUFFER, fbo);
}

int pFail(const char *message) {
	LOG_PRINT(message);
	return -1;
}

// converts Position<uint> into a reference to an element in the noise array
inline float noiseAtPos(Pos3i pos, ChkNoise &noise) {
	if (pos.x == INT_MIN) {
		return -1.f;
	}
	auto x = pos.x + 1;
	auto y = pos.y + 1;
	auto z = pos.z + 1;
	return noise.at(x + z * NOISE_W + y * (NOISE_W * NOISE_W));
}

// same as above + interpolation for floats
inline float noiseAtPos(glm::vec3 pos, ChkNoise &noise) {
    Vec3f b = floor(pos);
    Vec3f u = ceil(pos);
    float dx = pos.x - b.x;
    float dy = pos.y - b.y;
    float dz = pos.z - b.z;

    array<float, 8> v;
    v[0] = noiseAtPos(Pos3i{int(b.x), int(b.y), int(b.z)}, noise);
    v[1] = noiseAtPos(Pos3i{int(u.x), int(b.y), int(b.z)}, noise);
    v[2] = noiseAtPos(Pos3i{int(b.x), int(u.y), int(b.z)}, noise);
    v[3] = noiseAtPos(Pos3i{int(u.x), int(u.y), int(b.z)}, noise);
    v[4] = noiseAtPos(Pos3i{int(b.x), int(b.y), int(u.z)}, noise);
    v[5] = noiseAtPos(Pos3i{int(u.x), int(b.y), int(u.z)}, noise);
    v[6] = noiseAtPos(Pos3i{int(b.x), int(u.y), int(u.z)}, noise);
    v[7] = noiseAtPos(Pos3i{int(u.x), int(u.y), int(u.z)}, noise);
    return trilinearInterp(v, dx, dy, dz);
}

// ---------------------------- TREES -------------------------------------
void drawTrees(vector<Tree> trees, ShaderProgram *sp) {
	glDisable(GL_CULL_FACE);
	//glUniform1i(sp->u("isFlat"), GLint(1));
	for (auto &tree : trees) {
		glm::mat4 M = glm::mat4(1.f);
		M = glm::translate(M, tree.position);
		M = glm::rotate(M, tree.rotation, glm::vec3(0, 1, 0));
		M = glm::rotate(M, tree.skew, glm::vec3(1, 0, 0));
		drawVBO(tree.model.vbo, tree.model.indices, sp, M);
	}
}

void assignTreeMap(Pos3i &chkPos, BiomeChunk &bm, Chunk &chk) {
	Pos3i cPos = chkPos * CHK_SIZE;
	chk.treeMap = &bm.treeMap[Chunk::hasher(cPos)];
}

// testing whether there is enough space for a tree
bool testTreePlane(Pos3i &s, ChkNoise &density) {
	for (int x : {-1, 0, 1, 2}) {
	for (int z : {-1, 0, 1, 2}) {
		for (int y : {0, 1}) {
			Pos3i offset = {x, y, z};
			offset = offset + s;
			if (noiseAtPos(offset, density) < THRESHOLD) return false;
		}
		for (int y: {3}) {
			Pos3i offset = {x, y, z};
			offset = offset + s;
			if (noiseAtPos(offset, density) > THRESHOLD) return false;
		}
	}}
	for (int x : {0, 1}) {
	for (int z : {0, 1}) {
		for (int y : {2}) {
			Pos3i offset = {x, y, z};
			offset = offset + s;
			if (noiseAtPos(offset, density) < THRESHOLD) return false;
		}
	}}
	return true;
}

// [TODO:] refactor this shit
void placeTrees(Pos3i &cPos, Chunk *chk) {
	auto &rd = rand;
	auto &map = *chk->treeMap;

	map.mtx.lock();
	auto &vec = map.vec;
	// select trees
	auto iter = vec.begin();
	while(iter != vec.end()) {
		auto &t = *iter;
		bool candidate = false;
		int height = CHK_SIZE - 2;
		Pos3i sample = {
			int(glm::clamp(glm::floor(t[0]), 1.f, float(CHK_SIZE))),
			height, 
			int(glm::clamp(glm::floor(t[1]), 1.f, float(CHK_SIZE))),
		};
		while (height > 0 && candidate == false) {
			height--;
			sample.y = height;
			candidate = testTreePlane(sample, *chk->noise);
		}
		if (candidate) {
			if (height + cPos.y > LVL_SAND) {
				// selecting one of the pregenerated models
				auto &model = treemodels.at(rd() % treemodels.size());
				// y rotation
				float rotation = 2 * glm::pi<float>() * (rd() / float(RAND_MAX));
				// skew : -pi/32 -> +pi/32 angle
				float skew = glm::pi<float>() / 32.f + 
					(glm::pi<float>() / 16.f) * (rd() / float(RAND_MAX));
				auto position  = glm::vec3(
						cPos.x + t[0], 
						height + 1 + cPos.y, 
						cPos.z + t[1]
						);
				Tree tree = {
					.position 	= position,
					.rotation 	= rotation,
					.skew 		= skew,
					.model  	= model
				};
				chk->trees.push_back(tree);
			}
			// erasing the position from the sample vector
			iter = vec.erase(iter);
		}
		else iter++;
	}
	map.mtx.unlock();
}
// ------------------------------------------------------------------------


glm::vec3 calculateNormal(glm::vec3 v1, ChkNoise &noise) {
	// derivatives to figure out normal by using cross product
	// v1 - position of the vertex in the world space based on the noise
	Vec3f norm;
	norm.x = noiseAtPos(v1 + glm::vec3(1, 0, 0), noise) - noiseAtPos(v1 + glm::vec3(-1, 0, 0), noise);
	norm.y = noiseAtPos(v1 + glm::vec3(0, 1, 0), noise) - noiseAtPos(v1 + glm::vec3( 0,-1, 0), noise);
	norm.z = noiseAtPos(v1 + glm::vec3(0, 0, 1), noise) - noiseAtPos(v1 + glm::vec3( 0, 0,-1), noise);
	return -glm::normalize(norm);
}

void genTerrainType(Vertex &vert, BiomeSamples &biome, float x, float z) {
	auto hght = bilinearInterp(biome.height,  x, z);
	auto mntn = bilinearInterp(biome.mountain,x, z);
	auto shrp = bilinearInterp(biome.sharp,   x, z);
	auto humd = bilinearInterp(biome.humidity,x, z);

	//auto hscale = shrp * 2.25f;

	auto mtnH = LVL_ROCK * (0.45 + shrp);
	auto dirtH = LVL_DIRT * (0.45 + mntn);

	glm::vec4 humidScale = glm::vec4(humd, humd, humd, 1.f);

	constexpr glm::vec4 colSandDp  (0.03,  0.03,  0.05, 1.0); // underwater sand
	constexpr glm::vec4 colSand    (0.29,  0.29,  0.08, 1.0); // sand
	constexpr glm::vec4 colGrass   (0.014, 0.084, 0.018, 1.0); // grass
	constexpr glm::vec4 colGrassWet(0.048, 0.16,  0.04,  1.0);// humid grass
	constexpr glm::vec4 colGrassShr(0.066, 0.10,  0.03, 1.0); // sharp grass
	constexpr glm::vec4 colRock    (0.46,  0.46,  0.46, 1.0); // rock
	constexpr glm::vec4 colRock1   (0.32,  0.33,  0.34, 1.0); // dark rock
	constexpr glm::vec4 colSnow    (0.8,   0.8,   0.8,  1.0); // snow
	constexpr glm::vec4 colDirt    (0.21,  0.09,  0.07, 1.0); // dirt
	constexpr glm::vec4 colDirt1   (0.09,  0.05,  0.03, 1.0); // mountain dirt

	// calculating the angle of the surface
	auto dotProd = glm::dot(vert.norm, camera.up);

	// rock - angle
	if (dotProd < 0.2) {
		vert.color = colRock;
	}
	// mountain (rock -> snow)
	else if (vert.pos.y >= mtnH) {
		auto scale = (vert.pos.y - mtnH) / LVL_SNOW;
		vert.color = glm::mix(colRock1, colSnow, scale);
		// grass if its still flat
		if (dotProd > 0.85) {
			vert.color = glm::mix(colGrassShr, vert.color, scale);
		}
	} 
	// dirt mountain (dirt -> rock)
	else if (vert.pos.y >= dirtH) {
		auto scale = (vert.pos.y - dirtH) / (mtnH - dirtH);
		vert.color = glm::mix(colDirt1, colRock1, scale);
		// grass if its still flat
		if (dotProd > 0.8) vert.color = colGrassShr;
	}
	// dirt - angle
	else if (dotProd < 0.4) {
		vert.color = colDirt;
	} 
	// dark sand below water level
	else if (vert.pos.y < -LVL_SAND) { 	
		auto scal = glm::clamp(-(vert.pos.y + LVL_SAND) / 8.f, 0.f, 1.f);
		vert.color = glm::mix(colSand, colSandDp, scal);
	} 
	// sand
	else if (vert.pos.y <= LVL_SAND) {
		vert.color = colSand; 
	} 
	// grass
	else {
		auto sharpScal = (shrp - MIN_SHARP) * (1.f / (MAX_SHARP - MIN_SHARP));
		auto scal = glm::clamp(float(vert.pos.y - LVL_GRASS) / float(dirtH - LVL_GRASS), 0.f, 1.f);
		auto dirtScal = glm::pow(scal, 2);
		vert.color = glm::mix(colGrass, colGrassWet, humidScale);
		vert.color = glm::mix(vert.color, colGrassShr, sharpScal);
		vert.color = glm::mix(vert.color, colDirt1, dirtScal);
	}
}

void genMesh(Pos3i &curr, Chunk* chk) {
	// cubes containing vertex information for entire chunk
	auto &noise = chk->noise;
	auto &mesh  = chk->mesh;
	for (int iter = 0; iter < (MAP_W * MAP_W * MAP_H); iter++) {
		auto xz = iter % (MAP_W * MAP_W); // index on the <x,z> plane
		Pos3i cubePos = {
			.x = xz % MAP_W,
			.y = iter / (MAP_W * MAP_W),
			.z = xz / MAP_W
		};
		array<float, 8> cubevals;
		for (uint i = 0; i < 8; i++) {
			Pos3i nPos = {
				.x = cubePos.x + (int)marching_cubes_offsets[i].x,
				.y = cubePos.y + (int)marching_cubes_offsets[i].y,
				.z = cubePos.z + (int)marching_cubes_offsets[i].z,
			};
			cubevals[i] = noiseAtPos(nPos, *noise);
		}
		Mesh verts;
		polygonise(cubevals, verts);
		uint cntr = 0;
		array<Vertex, 3> face;
		for (auto &vert : verts) {
			// debugging
			if (vert.pos.x < 0){
				vert.pos.x = 0;
				vert.pos.y = 0;
				vert.pos.z = 0;
			}

			// chunkwise shift
			vert.pos.x += cubePos.x;
			vert.pos.y += cubePos.y;
			vert.pos.z += cubePos.z;

			vert.norm = calculateNormal(vert.pos, *noise);

			// adding worldspace shift
			vert.pos.x += curr.x;
			vert.pos.y += curr.y;
			vert.pos.z += curr.z;

			float relx, relz;
			relx = cubePos.x / float(MAP_W);
			relz = cubePos.z / float(MAP_W);

			// set colour
			genTerrainType(vert, chk->biome, relx, relz);

			face[cntr] = vert;
			cntr++;
			if (cntr == 3) {
				mesh.insert(mesh.end(), face.begin(), face.end());
				cntr = 0;
			} 
		}
	}
	chk->indices = mesh.size();
}

// ---------------------------- NOISE GENERATION --------------------------
static SimplexNoise hmapGen;
static SimplexNoise mountGen;
static SimplexNoise sharpGen;
static SimplexNoise humidGen;
static SimplexNoise forestGen;
// checks if a heightmap for the position has been generated, 
// if not then generate it and make other threads trying to 
// access that heightmap wait until its finished
BiomeChunk &setupBiome(Pos3i &curr) {
	// coordinates shifted to chunk origin point
	Pos3i chkCoord;
	chkCoord.x = glm::floor((float)curr.x / CHK_SIZE) * CHK_SIZE;
	chkCoord.z = glm::floor((float)curr.z / CHK_SIZE) * CHK_SIZE;
	chkCoord.y = 0;

	// hashing function alias
	auto key = Chunk::hasher;

	if (biomeMap.count(key(chkCoord)) != 0) {
		return biomeMap[key(chkCoord)];
	}
	// else
	mtx.lock();
	BiomeChunk &bm = biomeMap[key(chkCoord)];
	mtx.unlock();

	bm.mtx.lock();
	if (bm.state == JobState::generated) { 
		bm.mtx.unlock();
		return bm;
	}

	//LOG_INFO("Generating biomemaps...origin: %d %d, hash: %lu", 
			//chkCoord.x, chkCoord.z, key(chkCoord));
	uint octaves = 1;

	// noise generation
	for (int z = 0; z < bm.width; z++) {
	for (int x = 0; x < bm.width; x++) {
		auto offset = Pos3i{chkCoord.x + x, 0, chkCoord.z + z}; 

		// (-1; 2)
		auto mval = mountGen.fractal(octaves, offset.x, offset.z);
		mval = (1.f + mval) / 2;
		bm.at(bm.mountain, x, z) = mval;

		auto sval = sharpGen.fractal(octaves, offset.x, offset.z);
		sval = glm::clamp((sval + 2.f) / 12.f + 0.25f, MIN_SHARP, MAX_SHARP);
		bm.at(bm.sharp, x, z) = sval;

		// (-0.833; 1.833)
		auto hghval = hmapGen.fractal(octaves, offset.x, offset.z);
		hghval = 0.5f + hghval / 1.5f;
		hghval = glm::clamp(float(glm::exp(hghval) - 1), MIN_HEIGH, MAX_HEIGH);
		bm.at(bm.height, x, z) = hghval;

		auto humval = humidGen.fractal(octaves, offset.x, offset.z);
		humval = (humval + 2.f) / 4.f;
		// scaling 
		humval = (humval + MIN_HUMID) * (MAX_HUMID - MIN_HUMID);
		bm.at(bm.humidity, x, z) = humval;

		// forest (0, 10)
		auto forestval = (forestGen.fractal(octaves, offset.x, offset.z) + 2.f) / 4.f;
		forestval = (forestval + MIN_FORES) * (MAX_FORES - MIN_FORES);
		bm.at(bm.forest, x, z) = forestval;
	}}

	// tree placement for each chunk
	for (int z = 0; z < CHK_SIZE; z++) {
	for (int x = 0; x < CHK_SIZE; x++) {
		Pos3i cPos = chkCoord * CHK_SIZE;
		cPos.x += x * CHK_SIZE;
		cPos.z += z * CHK_SIZE;

		array<float, 4> f = bm.sample(bm.forest, cPos);
		auto avg = BiomeSamples::avg(f);
		float r = 32.f / avg;
		auto kxmin = array<float, 2>{0, 0};
		auto kxmax = array<float, 2>{CHK_SIZE - 1, CHK_SIZE - 1};
		auto ran = thinks::poisson_disk_sampling_internal::Hash(rand() % RAND_MAX);
		if (avg < 2.f) avg = 0.f;
		auto map = thinks::PoissonDiskSampling(r, kxmin, kxmax, avg * 2.f, ran);

		auto hash = Chunk::hasher(cPos);
		bm.treeMap[hash].vec = map;
	}
	}
	bm.state = JobState::generated;
	bm.mtx.unlock();
	return bm;
}

// generates 3D noise for use inside the chunk
void genNoise(Pos3i &curr, Chunk *chk) {
	uint octaves = 8;

	ChkNoise &noise = *chk->noise;

	// 2D position of the chunk (chunkspace)
	Pos3i chkPos;
	chkPos.x = glm::floor(float(curr.x) / CHK_SIZE);
	chkPos.z = glm::floor(float(curr.z) / CHK_SIZE);
	chkPos.y = 0;
	
	// will wait until its generated (unless one of the threads fails)
	auto &bm = setupBiome(chkPos);
	assignTreeMap(chkPos, bm, *chk);
	array<float, 4> h = bm.sample(bm.height, 	chkPos);
	array<float, 4> s = bm.sample(bm.sharp, 	chkPos);
	array<float, 4> m = bm.sample(bm.mountain, 	chkPos);
	array<float, 4> w = bm.sample(bm.humidity, 	chkPos);
	array<float, 4> f = bm.sample(bm.forest, 	chkPos);

	chk->biome.height   = h;
	chk->biome.sharp    = s;
	chk->biome.mountain = m;
	chk->biome.humidity = w;
	chk->biome.forest 	= f;

	SimplexNoise noiseGen(1.f, 1.f, 2.f, 0.5f);

	// noise generation loop
	for (int z = 0; z < NOISE_W; z++) {
	for (int x = 0; x < NOISE_W; x++) {
		// x, z as a position relative to the chunk - float [0..1]
		float relx = (x - 1) / float(bm.width);
		float relz = (z - 1) / float(bm.width);
		if (x < 1) relx = 0;
		if (z < 1) relz = 0;
		if (x >= bm.width) relx = 1;
		if (z >= bm.width) relz = 1;

		// roundness
		auto sharp = bilinearInterp(s, relx, relz);
		// height scale 0 - 1
		auto height = bilinearInterp(h, relx, relz);
		// humidity
		auto humid = bilinearInterp(w, relx, relz) / 5;
		// exponential scale of the noise y value, and the additional heightmap
		auto mountn = bilinearInterp(m, relx, relz);
		// turbulence
		auto turbGen = SimplexNoise(0, 0.33f, 1.f, 2.f, 0.5);
		for (int y = 0; y < NOISE_H; y++) {
			//auto yscale = 0.1f + mountn * 1.5f;
			// 3D noise generation
			auto pos = Pos3i{x + curr.x - 1, y + curr.y - 1, z + curr.z - 1};
			float yscale = glm::clamp(float(glm::exp(mountn * 1.15) - 1), MIN_MOUNT, MAX_MOUNT);
			float val;
			val = -humid + height - pos.y / (yscale * CHK_SIZE);
			if (THRESHOLD > 2 + val) {
				val = -1.f;
				noise.at(x + z * NOISE_W + y * (NOISE_W * NOISE_W)) = val;
				continue;
			} else if (THRESHOLD < -2 + val) {
				val = 1.f;
				noise.at(x + z * NOISE_W + y * (NOISE_W * NOISE_W)) = val;
				continue;
			}

			auto turbulence = turbGen.noise(pos.x, pos.z) / 75.f;
			noiseGen = SimplexNoise(SEED, 0.012f, 1.f, 2.f, sharp);
			val += noiseGen.fractal(octaves, pos.x, pos.y / yscale, pos.z) + turbulence;
			noise.at(x + z * NOISE_W + y * (NOISE_W * NOISE_W)) = val;
		}
	}}
}
// ------------------------------------------------------------------------

GLuint bindMesh(Mesh &mesh) {
	GLuint vbo;
	glGenBuffers(1, &vbo);
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBufferData(GL_ARRAY_BUFFER, mesh.size() * sizeof(Vertex),
         	(void*)mesh.data(), GL_STATIC_DRAW);
    return vbo;
}

// generates chunk data for the Chunk passed as an argument
void genChunk(Pos3i cPos, Chunk *chk) {
	chk->mesh.reserve(Chunk::MAX_MESH_VERTS);
	genNoise(cPos, chk);
	genMesh(cPos, chk);
	if (cPos.y >= 0 && cPos.y <= CHK_SIZE) placeTrees(cPos, chk);

	// awaiting retrival
	chk->state = JobState::awaiting;
}

// ---------------------------- PREGENERATION -----------------------------
void pregenChunks(int dist = 8) {
	boost::asio::thread_pool chPool(std::thread::hardware_concurrency());
	for (int x = -dist; x <= dist; x++) {
	for (int y = -1; y <= 3; y++) {
	for (int z = -dist; z <= dist; z++) {
		Pos3i cPos = {
			.x = x * CHK_SIZE,
			.y = y * CHK_SIZE,
			.z = z * CHK_SIZE
		};
		auto hash = Chunk::hasher(cPos);
		Chunk &chk = chunkMap[hash];
		chk.state = JobState::delegated;
		chk.noise = new ChkNoise;
		// delegate chunk generation to a thread
		boost::asio::post(chPool, std::bind(genChunk, cPos, &chk));
	}}}
	chPool.join();
}

void pregenTrees() {
	for (auto &tree : treemodels) {
		Model model = Model();
		auto pos = glm::vec3(0, 0, 0);
		auto M = glm::translate(glm::mat4(1.f), pos);
		auto scale = glm::vec3(1.f, 1.f, 1.f);
		auto mesh = model.genTree(glm::vec3(0,0,0), scale, M);
    	auto vbo = bindMesh(mesh);
    	tree.indices = mesh.size();
    	tree.vbo = vbo;
	}
}
// ------------------------------------------------------------------------

void drawWater(ShaderProgram* sp) {
	glDisable(GL_CULL_FACE);
	glm::mat4 M = glm::mat4(1.f);
	glm::mat4 V = glm::lookAt(camera.pos, camera.pos + camera.target, camera.up);
	glm::mat4 P = glm::perspective(glm::radians(FOV), ASPECT_RATIO, Z_NEAR, Z_FAR);
	M = glm::translate(M, {
			camera.pos.x - (CHK_SIZE * RENDER_DIST), 
			0, 
			camera.pos.z - (CHK_SIZE * RENDER_DIST)
			});

	sp->use();

	glUniform3fvARB(sp->u("skyColor"), 1, glm::value_ptr(skyColor));
	glUniformMatrix4fv(sp->u("V"), 1, false, glm::value_ptr(V));
	glUniformMatrix4fv(sp->u("P"), 1, false, glm::value_ptr(P));

	glUniformMatrix4fv(sp->u("M"), 1, false, glm::value_ptr(M));
	glUniform3fvARB(sp->u("cameraPos"), 1, glm::value_ptr(camera.pos));

	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, m_depthTexture);

	glBindVertexArray(vao);
	glBindBuffer(GL_ARRAY_BUFFER, waterVBO);

	uint count = 3;
	// passing position vector to vao
	glVertexAttribPointer(sp->a("vertex"), count, GL_FLOAT, GL_FALSE,
        	sizeof(Vertex), (GLvoid*)offsetof(Vertex, pos));
	glEnableVertexAttribArray(sp->a("vertex"));

	// passing normal vector to vao
	glVertexAttribPointer(sp->a("normal"), count, GL_FLOAT, GL_FALSE,
        	sizeof(Vertex), (GLvoid*)offsetof(Vertex, norm));
	glEnableVertexAttribArray(sp->a("normal"));

	glVertexAttribPointer(sp->a("color"), count + 1, GL_FLOAT, GL_FALSE,
        	sizeof(Vertex), (GLvoid*)offsetof(Vertex, color));
	glEnableVertexAttribArray(sp->a("color"));

	glDrawArrays(GL_TRIANGLES, 0, 6);

	glBindVertexArray(0);
}

void genWater() {
	glGenBuffers(1, &waterVBO);
	glBindBuffer(GL_ARRAY_BUFFER, waterVBO);
	// visible water surface
	auto X = float(4 * RENDER_DIST * CHK_SIZE);
	const array<Vertex, 3 * 2> verts = {
		Vertex{{-X,0,-X},{0,1,0}, waterColor},
		Vertex{{ X,0, X},{0,1,0}, waterColor},
		Vertex{{-X,0, X},{0,1,0}, waterColor},
		Vertex{{-X,0,-X},{0,1,0}, waterColor},
		Vertex{{ X,0,-X},{0,1,0}, waterColor},
		Vertex{{ X,0, X},{0,1,0}, waterColor}
	};
	glBufferData(GL_ARRAY_BUFFER, verts.size() * sizeof(Vertex),
         	(void*)verts.data(), GL_STATIC_DRAW);
}

// LAZY copies a depth buffer to the texture
void copyDepthBuffer(GLuint texture) {
	glBindTexture(GL_TEXTURE_2D, texture);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT, WINDOW_W, WINDOW_H, 0, GL_DEPTH_COMPONENT, GL_UNSIGNED_BYTE, 0);
	glCopyTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, 0, 0, WINDOW_W, WINDOW_H);
}

void drawChunk(Chunk& chk, ShaderProgram *sp) {
	glEnable(GL_CULL_FACE);
	drawVBO(chk.vbo, chk.indices, sp);
}

void renderChunks(ShaderProgram *sp, ShaderProgram *sptree) {
	sp->use();
	//glBindFramebuffer(GL_FRAMEBUFFER, fbo);

	ImGui::Begin("DEBUG INFO");
	ImGui::Text("x: %0.1f\ny: %0.1f\nz: %0.1f\n", 
			camera.pos.x, camera.pos.y, camera.pos.z);
	ImGui::Text("[%i %i %i]\n\n", currChkPos.x, currChkPos.y, currChkPos.z);

	// rendering or queueing creation of chunks
	for (int x = -RENDER_DIST;   x <= RENDER_DIST;   x++) {
	for (int y = -RENDER_DIST/3; y <= RENDER_DIST/3; y++) {
	for (int z = -RENDER_DIST;   z <= RENDER_DIST;   z++) {
		// skip generation if outside the rendering circle
		if (x*x + y*y + z*z > RENDER_RADIUS) continue;
		Pos3i cPos = {
			.x = (currChkPos.x + x) * CHK_SIZE,
			.y = (currChkPos.y + y) * CHK_SIZE,
			.z = (currChkPos.z + z) * CHK_SIZE
		};
		auto hash = Chunk::hasher(cPos);
		Chunk &chk = chunkMap[hash];
		BiomeChunk bm;

		if (chk.state == JobState::awaiting) {
			// submitting the mesh to gpu
			if (chk.mesh.size() > 0) {
				chk.vbo = bindMesh(chk.mesh);
			}
			chk.clear();
			chk.state = JobState::generated;
		} else if (chk.state == JobState::none){
			chk.state = JobState::delegated;
			// will be removed later
			chk.noise = new ChkNoise;
			// delegate chunk generation to a thread
			boost::asio::post(pool, std::bind(genChunk, cPos, &chk));
			//genChunk(cPos, &chk);
		} 
		if (chk.state == JobState::generated) {
			if (y == 0 && z == 0 && x == 0) {
				ImGui::Text("Biome:\n");
				ImGui::Text("Sharpness: %.2f\n", chk.biome.avg(chk.biome.sharp));
				ImGui::Text("Mountain:  %.2f\n", chk.biome.avg(chk.biome.mountain));
				ImGui::Text("Humidity:  %.2f\n", chk.biome.avg(chk.biome.humidity));
				ImGui::Text("Forest:    %.2f\n", chk.biome.avg(chk.biome.forest));
				ImGui::Text("Height:    %.2f\n", chk.biome.avg(chk.biome.height));
			}
			// drawing a chunk
			drawChunk(chk, sp);
			// drawing the trees
			drawTrees(chk.trees, sptree);
		}
	}}}
	// copy depth to a texture
	copyDepthBuffer(m_depthTexture);
	ImGui::End();
}

int main(int argc, char* argv[]) {
	srand(time(NULL));
	if (argc == 2) {
		SEED = strtoul(argv[1], NULL, 10);
		LOG_INFO("SEED: %u", SEED);
	} else {
		SEED = rand();
	}
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
		// freeing vao
		glDeleteBuffers(1, &vao);
		// freeing vbos on the gpu
		vector<GLuint> vbos;
		for (auto &chk : chunkMap) {
			vbos.push_back(chk.second.vbo);
		}
		for (auto &model : treemodels) {
			vbos.push_back(model.vbo);
		}
		glDeleteBuffers(vbos.size(), vbos.data());
		// freeing fbo
		//glDeleteFramebuffers(1, &fbo);

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
	};

	auto shaderDestoyer = [&](ShaderProgram* sp){delete sp;};
	auto shaderCreate = [](const char* vert, const char* frag) {
		ShaderProgram *sp = new ShaderProgram(vert, frag);
		sp->use();
		glUniform3fvARB(sp->u("skyColor"), 1, glm::value_ptr(skyColor));
		glUniform1f(sp->u("far"), Z_FAR);
		glUniform1f(sp->u("near"), Z_NEAR);
		glUniform3fv(sp->u("sunDir"), 1, glm::value_ptr(
					glm::normalize(SUN_DIRECTION)));
		glUniform1f(sp->u("density"),  FOG_DENSITY);
		glUniform1f(sp->u("gradient"), FOG_GRADIENT);
		return sp;
	};

	// terrain shader
	LOG_INFO("INITIALIZING TERRAIN SHADER");
	Uniq_Ptr<ShaderProgram, decltype(shaderDestroyer)> spterr(
			shaderCreate("glsl/vshad.glsl", "glsl/fshad.glsl"),
			shaderDestroyer);
	// water shader
	LOG_INFO("INITIALIZING WATER SHADER");
	Uniq_Ptr<ShaderProgram, decltype(shaderDestoyer)> spwater(
			shaderCreate("glsl/vwater.glsl", "glsl/fwater.glsl"),
			shaderDestoyer);
	// tree shader
	LOG_INFO("INITIALIZING TREE SHADER");
	Uniq_Ptr<ShaderProgram, decltype(shaderDestoyer)> sptree(
			shaderCreate("glsl/vtree.glsl", "glsl/ftree.glsl"),
			shaderDestoyer);
	initProgram(window.get());

	{// noise generators initialisation
		hmapGen  = SimplexNoise((SEED * 11) % UINT_MAX, 0.02f, 1.0f,2.f, 0.5f);
		mountGen = SimplexNoise((SEED * 17) % UINT_MAX, 0.05f, 1.4f,2.f, 0.5f);
		sharpGen = SimplexNoise((SEED * 23) % UINT_MAX, 0.02f, 1.f, 2.f, 0.5f);
		humidGen = SimplexNoise((SEED * 37) % UINT_MAX, 0.01f, 1.f, 2.f, 0.5f);
		forestGen= SimplexNoise((SEED * 71) % UINT_MAX, 0.1f,  1.f, 2.f, 0.5f);
	}
	{// generating things before rendering
		genWater();
    	pregenTrees();
    	LOG_INFO("Done generating trees!");
    	// waiting for chunk generation to complete
    	LOG_INFO("Map pre-generation...");
    	//pregenChunks(12);
    	//LOG_INFO("Done generating chunks!");
    }
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

    	renderChunks(spterr.get(), sptree.get());
    	drawWater(spwater.get());

		ImGui::Render();
		ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

		glfwSwapBuffers(window.get());
	}
	LOG_INFO("Closing the program...");
	pool.stop();
	pool.join();
}
