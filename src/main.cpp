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
#include <condition_variable>

#include <vector>
#include <array>
#include <memory>
#include <unordered_map>

#include <boost/asio/thread_pool.hpp>
#include <boost/asio/post.hpp>

#include "utils.hpp"
#include "imgui/imgui.h"
#include "imgui/imgui_impl_glfw.h"
#include "imgui/imgui_impl_opengl3.h"

#include "SimplexNoise/SimplexNoise.h"
#include "lodepng/lodepng.h"
#include "shaderprogram.hpp"
#include "marchingcubes.hpp"
#include "log.hpp"

#include "main.hpp"

// colours
constexpr glm::vec4 waterColor(0.25, 0.516, 0.914, 0.46);
constexpr glm::vec3 skyColor(0.45, 0.716, 0.914);

// [TODO]: remove later
Pos3i currChkPos;
boost::asio::thread_pool pool(std::thread::hardware_concurrency() - 1);

GLuint vao;
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
struct BiomeChunk {
	// noises used by 3D noise generator
	JobState state = none;
	std::mutex mtx;
	static constexpr int width = CHK_SIZE + 1;
	using biomeNoise = array<float, width * width>;

	biomeNoise heightmap;
	biomeNoise humidity;
	biomeNoise smoothness;

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

struct Chunk {
	// mesh vertex size
	size_t indices = 0;
	// after a thread generates the mesh, main thread will submit it to GPU
	GLuint vbo;
	vector<Vertex> mesh;

	// 3D noise used for terrain generation
	ChkNoise noise;
	// status
	std::atomic<JobState> state = ATOMIC_VAR_INIT(JobState::none);
};

// storage for all chunks
std::unordered_map<ChkHash, Chunk> chunkMap;
// storage for biomes, 1 cell = 1 chunk
std::unordered_map<ChkHash, BiomeChunk> biomeMap;

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

void drawChunk(Chunk &chk, ShaderProgram* sp) {
	glBindVertexArray(vao);
	glBindBuffer(GL_ARRAY_BUFFER, chk.vbo);

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

	glDrawArrays(GL_TRIANGLES, 0, chk.indices);

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

ShaderProgram* initProgram(GLFWwindow *window) {
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

	ShaderProgram* sp = new ShaderProgram("glsl/vshad.glsl", "glsl/fshad.glsl");
	sp->use();
	glUniform3fvARB(sp->u("skyColor"), 1, glm::value_ptr(skyColor));

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
ChkHash chunkHasher(Pos3i &pos) {
	ChkHash x = glm::floor(pos.x / CHK_SIZE);
	ChkHash y = glm::floor(pos.y / CHK_SIZE);
	ChkHash z = glm::floor(pos.z / CHK_SIZE);
	x = x << 0;
	z = z << 16;
	y = y << 32;
	return x + y + z;
};

// converts Position<uint> into a reference to an element in the noise array
inline float& noiseAtPos(Pos3i pos, ChkNoise &noise) {
	auto x = pos.x + 1;
	auto y = pos.y + 1;
	auto z = pos.z + 1;
	return noise.at(x + z * NOISE_W + y * (NOISE_W * NOISE_W));
}

// same as above + interpolation for floats
inline float noiseAtPos(glm::vec3 pos, ChkNoise &noise) {
	// no interpolation (poormans normal)
	//return noiseAtPos(Pos3i{int(pos.x), int(pos.y), int(pos.z)}, noise);
    Vec3f b = floor(pos);
    Vec3f u = ceil(pos);
    float dx = pos.x - b.x;
    float dy = pos.y - b.y;
    float dz = pos.z - b.z;

    array<float, 8> v;
    v[0] = noiseAtPos(Pos3i{(int)b.x, (int)b.y, (int)b.z}, noise);
    v[1] = noiseAtPos(Pos3i{(int)u.x, (int)b.y, (int)b.z}, noise);
    v[2] = noiseAtPos(Pos3i{(int)b.x, (int)u.y, (int)b.z}, noise);
    v[3] = noiseAtPos(Pos3i{(int)u.x, (int)u.y, (int)b.z}, noise);
    v[4] = noiseAtPos(Pos3i{(int)b.x, (int)b.y, (int)u.z}, noise);
    v[5] = noiseAtPos(Pos3i{(int)u.x, (int)b.y, (int)u.z}, noise);
    v[6] = noiseAtPos(Pos3i{(int)b.x, (int)u.y, (int)u.z}, noise);
    v[7] = noiseAtPos(Pos3i{(int)u.x, (int)u.y, (int)u.z}, noise);
    return trilinearInterp(v, dx, dy, dz);
}

glm::vec3 calculateNormal(glm::vec3 v1, ChkNoise &noise) {
	// derivatives to figure out normal by using cross product
	// v1 - position of the vertex in the world space based on the noise
	Vec3f norm;
	norm.x = noiseAtPos(v1 + glm::vec3(1, 0, 0), noise) - noiseAtPos(v1 + glm::vec3(-1, 0, 0), noise);
	norm.y = noiseAtPos(v1 + glm::vec3(0, 1, 0), noise) - noiseAtPos(v1 + glm::vec3( 0,-1, 0), noise);
	norm.z = noiseAtPos(v1 + glm::vec3(0, 0, 1), noise) - noiseAtPos(v1 + glm::vec3( 0, 0,-1), noise);
	return -glm::normalize(norm);
}

void genMesh(Pos3i &curr, Chunk* chk) {
	// cubes containing vertex information for entire chunk
	auto &noise = chk->noise;
	auto &mesh  = chk->mesh;
	uint constexpr MAX_MESH_VERTS = 12 * CHK_SIZE * CHK_SIZE * CHK_SIZE * 3;
	mesh.reserve(MAX_MESH_VERTS);
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
			cubevals[i] = noiseAtPos(nPos, noise);
		}
		vector<Vertex> verts;
		polygonise(cubevals, verts);

		uint cntr = 0;
		array<Vertex, 3> face;
		for (auto &vert : verts) {
			// chunkwise shift
			vert.pos.x += cubePos.x;
			vert.pos.y += cubePos.y;
			vert.pos.z += cubePos.z;

			vert.norm = calculateNormal(vert.pos, noise);

			// adding worldspace shift
			vert.pos.x += curr.x;
			vert.pos.y += curr.y;
			vert.pos.z += curr.z;

			constexpr glm::vec3 colMask(0.26,  0.26, 0.03);
			constexpr glm::vec4 colSand(0.27,  0.27, 0.07, 1.0);
			constexpr glm::vec4 colGrass(0.01, 0.07, 0.0,  1.0);
			constexpr glm::vec4 colRock(0.07,  0.07, 0.07, 1.0);

			// TODO: remove,  only temporary 
			if (vert.pos.y < -2.f) {
				auto scal = std::min(std::max(-(vert.pos.y + 2.f) / 8.f, 0.f), 1.f);
				vert.color = colSand - glm::vec4(colMask, 0.f) * scal;
			}
			else if (vert.pos.y <= 2.f) vert.color = colSand; 
			//else if (vert.pos.y <= 32.f) vert.color = colGrass;
			else if (vert.pos.y <= 28.f) vert.color = colGrass;
			else vert.color = colRock;

			//vert.color = glm::vec4(0.01, 0.07, 0.0, 1.f);

			//mesh.push_back(vert);
			face[cntr] = vert;
			cntr++;
			if (cntr == 3) {
				mesh.insert(mesh.end(), face.begin(), face.end());
				cntr = 0;
			} 
		}
	}
	mesh.shrink_to_fit();
}
// -------------------------------------------------------

void saveNoisePNG(const char* fname, ChkNoise &noise) {
	LOG_INFO("Saving noise pictures...");
	array<byte, NOISE_W * NOISE_W> picture;
	float max = -99999.f;
	for (uint i = 0; i < noise.size(); i++) {
		picture.at(i) = (1 + noise.at(i)) * 128.f;
	}
	lodepng_encode_file(fname, picture.data(), 
			NOISE_W, NOISE_W, LCT_GREY, 8);
}

// [TODO:] seed consistency concern
SimplexNoise hmapGen(0.03f, 1.f, 2.f, 0.5f);
SimplexNoise humidGen(0.01f, 1.f, 2.f, 0.5f);
SimplexNoise smoothGen(0.02f, 1.f, 2.f, 0.5f);
// checks if a heightmap for the position has been generated, 
// if not then generate it and make other threads trying to 
// access that heightmap wait until its finished
BiomeChunk& updateBiomes(Pos3i &curr) {
	// coordinates shifted to chunk origin point
	Pos3i chkCoord;
	chkCoord.x = glm::floor((float)curr.x / CHK_SIZE) * CHK_SIZE;
	chkCoord.z = glm::floor((float)curr.z / CHK_SIZE) * CHK_SIZE;
	chkCoord.y = 0;

	// hashing function alias
	auto key = chunkHasher;

	if (biomeMap.count(key(chkCoord)) != 0) {
		return biomeMap[key(chkCoord)];
	}
	// else
	LOG_INFO("Entering critical section...");
	mtx.lock();
	BiomeChunk 	&bm = biomeMap[key(chkCoord)];
	LOG_INFO("Leaving critical section...");
	mtx.unlock();

	bm.mtx.lock();
	if (bm.state == JobState::generated) { 
		bm.mtx.unlock();
		return bm;
	}

	LOG_INFO("Generating biomemaps...origin: %d %d, hash: %lu", 
			chkCoord.x, chkCoord.z, key(chkCoord));
	uint octaves = 5;

	// noise generation
	for (int z = 0; z < bm.width; z++) {
	for (int x = 0; x < bm.width; x++) {
		auto offset = Pos3i{chkCoord.x + x, 0, chkCoord.z + z}; 

		// humidity
		auto wval = humidGen.fractal(octaves, offset.x, offset.z);
		wval = (1.f + wval) / 2;
		bm.at(bm.humidity, x, z) = wval;

		// smoothness
		auto sval = smoothGen.fractal(octaves, offset.x, offset.z);
		sval = glm::clamp((sval + 0.75f + wval) / 2.f, 0.35f, 0.5f);
		bm.at(bm.smoothness, x, z) = sval;

		// height
		auto hval = 0.5f + hmapGen.fractal(octaves, offset.x, offset.z) / 3.f;
		hval += -wval;
		bm.at(bm.heightmap, x, z) = hval;
	}}
	bm.state = JobState::generated;
	bm.mtx.unlock();
	return bm;
}

// generates 3D noise for use inside the chunk
void genNoise(Pos3i &curr, Chunk *chk) {
	uint octaves = 8;

	ChkNoise &noise = chk->noise;

	// 2D position of the chunk (chunkspace)
	Pos3i chkPos;
	chkPos.x = glm::floor(float(curr.x) / CHK_SIZE);
	chkPos.z = glm::floor(float(curr.z) / CHK_SIZE);
	chkPos.y = 0;
	
	// will wait until its generated (unless one of the threads fails)
	auto &bm = updateBiomes(chkPos);
	//auto &bm = updateBiomemap(chkCoord);
	array<float, 4> h = bm.sample(bm.heightmap,  chkPos);
	array<float, 4> s = bm.sample(bm.smoothness, chkPos);
	array<float, 4> w = bm.sample(bm.humidity,   chkPos);
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
		auto smooth = bilinearInterp(s, relx, relz);
		// height cutoff
		auto cutoff = bilinearInterp(h, relx, relz);
		// humidity
		auto humid = bilinearInterp(w, relx, relz);
		for (int y = 0; y < NOISE_H; y++) {
			//if (y == 0 && x == 0 && z == 0) LOG_INFO("smooth: %f, octaves: %u", smooth, octaves);
			noiseGen = SimplexNoise(0.008f, 1.f, 2.f, smooth);

			// actual 3D noise generation
			auto pos = Pos3i{x + curr.x - 1, y + curr.y - 1, z + curr.z - 1};
			float val;
			val = noiseGen.fractal(octaves, pos.x, pos.y, pos.z);
			val += cutoff - pos.y / (2.f * CHK_SIZE);
			//val -= (pos.y - 24.f) / 64.f;
			//val = -pos.y + smooth * 100.f;

			noise.at(x + z * NOISE_W + y * (NOISE_W * NOISE_W)) = val;
		}
	}}
}

// generates chunk data for the Chunk passed as an argument
void genChunk(Pos3i cPos, Chunk *chk) {
	genNoise(cPos, chk);
	genMesh(cPos, chk);
	chk->indices   = chk->mesh.size();

	// awaiting retrival
	chk->state = JobState::awaiting;
}

void drawWater(ShaderProgram* sp) {
	glm::mat4 M = glm::mat4(1.f);
	M = glm::translate(M, {camera.pos.x - (CHK_SIZE * RENDER_DIST), 0, camera.pos.z - (CHK_SIZE * RENDER_DIST)});
	glUniformMatrix4fv(sp->u("M"), 1, false, glm::value_ptr(M));

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

	glDisable(GL_CULL_FACE);
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

void renderChunks(GLFWwindow *window, ShaderProgram *sp) {
	sp->use();

	glEnable(GL_CULL_FACE);
	glm::mat4 M = glm::mat4(1.f);
	glm::mat4 V = glm::lookAt(camera.pos, camera.pos + camera.target, camera.up);
	glm::mat4 P = glm::perspective(glm::radians(FOV), ASPECT_RATIO, Z_NEAR, Z_FAR);

	glUniformMatrix4fv(sp->u("M"), 1, false, glm::value_ptr(M));
	glUniformMatrix4fv(sp->u("V"), 1, false, glm::value_ptr(V));
	glUniformMatrix4fv(sp->u("P"), 1, false, glm::value_ptr(P));

	ImGui::Begin("DEBUG INFO");
	ImGui::Text("x: %0.1f\n y: %0.1f\n z: %0.1f\n", camera.pos.x, camera.pos.y, camera.pos.z);
	ImGui::Text("[%i %i %i]\n\n", currChkPos.x, currChkPos.y, currChkPos.z);

	// rendering or queueing creation of chunks
	for (int x = -RENDER_DIST;   x <= RENDER_DIST;   x++) {
	for (int y = -RENDER_DIST/3; y <= RENDER_DIST/3; y++) {
	for (int z = -RENDER_DIST;   z <= RENDER_DIST;   z++) {
		Pos3i cPos = {
			.x = (currChkPos.x + x) * CHK_SIZE,
			.y = (currChkPos.y + y) * CHK_SIZE,
			.z = (currChkPos.z + z) * CHK_SIZE
		};
		auto hash = chunkHasher(cPos);
		Chunk &chk = chunkMap[hash];
		BiomeChunk bm;

		if (chk.state == JobState::awaiting) {
			// submitting the mesh to gpu
			if (chk.mesh.size() > 0) {
				glBindBuffer(GL_ARRAY_BUFFER, chk.vbo);
				glBufferData(GL_ARRAY_BUFFER, chk.indices * sizeof(Vertex),
         				(void*)chk.mesh.data(), GL_STATIC_DRAW);
         		chk.mesh.clear();
			}
			chk.state = JobState::generated;
			drawChunk(chk, sp);
		} else if (chk.state == JobState::none){
			chk.state = JobState::delegated;
			// vbo for the chunk
			GLuint vbo;
			glGenBuffers(1, &vbo);
			chk.vbo = vbo;

			//LOG_INFO("Generating chunk [%i %i %i]...\nHash: %lu", cPos.x / CHK_SIZE, 
				//cPos.y / CHK_SIZE, cPos.z / CHK_SIZE, hash);
			boost::asio::post(pool, std::bind(genChunk, cPos, &chk));
			//genChunk(cPos, &chk);
		} else drawChunk(chk, sp); 
	}}}
    drawWater(sp);
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
		for (auto &chk : chunkMap) {
			vbos.push_back(chk.second.vbo);
		}
		glDeleteBuffers(vbos.size(), vbos.data());
		delete sp;
	};
	// pointer to a ShaderProgram
	Uniq_Ptr<ShaderProgram, decltype(shaderDestroyer)> sp(
			initProgram(window.get()),
			shaderDestroyer);

	genWater();

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
	LOG_INFO("Closing the program...");
	pool.stop();
	pool.join();
}
