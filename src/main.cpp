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

//#include <omp.h>

#include "imgui/imgui.h"
#include "imgui/imgui_impl_glfw.h"
#include "imgui/imgui_impl_opengl3.h"

#include "SimplexNoise/SimplexNoise.h"
#include "lodepng/lodepng.h"
#include "shaderprogram.hpp"
#include "marchingcubes.hpp"
#include "log.hpp"

#include "main.hpp"

// color of the water
glm::vec4 waterColor(0.25, 0.516, 0.914, 0.46);

// unique_ptr alias
template<typename... T>
using Uniq_Ptr = std::unique_ptr<T...>;

// [TODO]: remove later
Pos3i currChkPos;

GLuint vao;

float deltaTime = 0.f;
float lastFrame = 0.f;

// mouse position
Position<float> mouseLast;

using ChkNoise = array<float, NOISE_W * NOISE_W * NOISE_H>;
using ChkHmap = array<float, NOISE_W * NOISE_W>;

// ------------------------- chunks ----------------------
struct Chunk {
	// mesh vertex size
	size_t indices = 0;
	// when a thread generating a new chunk dies the mesh will be submitted
	// to vbo by the main thread and the vector removed
	GLuint vbo;
	vector<Vertex> mesh;

	// 3D noise used for terrain generation
	ChkNoise noise;

	std::mutex isUsed;
	bool generated = false;
};

// used for hashing
using ChkCoord = uint64_t;

std::unordered_map<ChkCoord, Chunk> chunkMap;

struct Camera {
	glm::vec3 pos = glm::vec3(0.f, 0.f, 3.f);
	glm::vec3 target = glm::vec3(0.f, 0.f, 0.f);
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
	glEnable(GL_DEBUG_OUTPUT);
	glDebugMessageCallback(MessageCallback, 0);
#endif

	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);

	glFrontFace(GL_CW);
	glEnable(GL_BLEND);  

	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glClearColor(0.45, 0.716, 0.914, 1.f);

	glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
	glfwSetCursorPosCallback(window, mouseCallback);
	glfwSetKeyCallback(window, keyCallback);

	glfwSetWindowSizeCallback(window, windowResizeCallback);

	ShaderProgram* sp = new ShaderProgram("glsl/vshad.glsl", "glsl/fshad.glsl");

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
ChkCoord chunkHasher(Pos3i &pos) {
	ChkCoord x = glm::floor(pos.x / CHK_SIZE);
	ChkCoord y = glm::floor(pos.y / CHK_SIZE);
	ChkCoord z = glm::floor(pos.z / CHK_SIZE);
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

float trilinearInterp(Vec3f pos, float vals[8]) {
  auto &x = pos.x;
  auto &y = pos.y;
  auto &z = pos.z;

  // Calculate the weights for each of the 8 input values
  float wx = 1 - std::abs(x - std::floor(x));
  float wy = 1 - std::abs(y - std::floor(y));
  float wz = 1 - std::abs(z - std::floor(z));

  // Calculate the weighted average of the 8 input values
  float result = 0;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      for (int k = 0; k < 2; k++) {
        float weight = wx * wy * wz;
        result += vals[4 * i + 2 * j + k] * weight;
        wz = 1 - wz;
      }
      wy = 1 - wy;
    }
    wx = 1 - wx;
  }
  return result;
}


// same as above + interpolation for floats
inline float noiseAtPos(glm::vec3 pos, ChkNoise &noise) {
	// no interpolation (poormans normal)
	//return noiseAtPos(Pos3i{int(pos.x), int(pos.y), int(pos.z)}, noise);
    Vec3f b = floor(pos);
    Vec3f u = ceil(pos);
    float v[8];
    v[0] = noiseAtPos(Pos3i{(int)b.x, (int)b.y, (int)b.z}, noise);
    v[1] = noiseAtPos(Pos3i{(int)u.x, (int)b.y, (int)b.z}, noise);
    v[2] = noiseAtPos(Pos3i{(int)b.x, (int)u.y, (int)b.z}, noise);
    v[3] = noiseAtPos(Pos3i{(int)u.x, (int)u.y, (int)b.z}, noise);
    v[4] = noiseAtPos(Pos3i{(int)b.x, (int)b.y, (int)u.z}, noise);
    v[5] = noiseAtPos(Pos3i{(int)u.x, (int)b.y, (int)u.z}, noise);
    v[6] = noiseAtPos(Pos3i{(int)b.x, (int)u.y, (int)u.z}, noise);
    v[7] = noiseAtPos(Pos3i{(int)u.x, (int)u.y, (int)u.z}, noise);
    return trilinearInterp(pos, v);
}

glm::vec3 calculateNormal(glm::vec3 v1, ChkNoise &noise) {
	// derivatives to figure out normal by using cross product
	// v1 - position of vertex in world space based on the noise
	Vec3f norm;
	norm.x = noiseAtPos(v1 + glm::vec3(1, 0, 0), noise) - noiseAtPos(v1 + glm::vec3(-1, 0, 0), noise);
	norm.y = noiseAtPos(v1 + glm::vec3(0, 1, 0), noise) - noiseAtPos(v1 + glm::vec3( 0,-1, 0), noise);
	norm.z = noiseAtPos(v1 + glm::vec3(0, 0, 1), noise) - noiseAtPos(v1 + glm::vec3( 0, 0,-1), noise);
	return -glm::normalize(norm);
	//return glm::vec3(0,0,1);
}

vector<Vertex> genMesh(Pos3i &curr, ChkNoise &noise) {
	// cubes containing vertex information for entire chunk
	//array<vector<Vertex>, MAP_W * MAP_W * MAP_H> cubes;
	vector<Vertex> mesh;
	uint constexpr MAX_MESH_VERTS = 12 * CHK_SIZE * CHK_SIZE * CHK_SIZE * 3;
	mesh.reserve(MAX_MESH_VERTS);
//#pragma omp parallel for
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

			constexpr glm::vec3 colMask(0.15,  0.15, 0.03);
			constexpr glm::vec4 colSand(0.27,  0.27, 0.07, 1.0);
			constexpr glm::vec4 colGrass(0.01, 0.07, 0.0,  1.0);
			constexpr glm::vec4 colRock(0.07,  0.07, 0.07, 1.0);

			// only temporary 
			if (vert.pos.y < -2.f) vert.color = 
				colSand - glm::vec4(colMask, 0.f) * std::max(glm::abs(vert.pos.y / 4), 1.f);
			else if (vert.pos.y <= 2.f) vert.color = colSand; 
			else if (vert.pos.y <= 32.f) vert.color = colGrass;
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
	return mesh;
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

ChkNoise genNoise(Pos3i &curr) {
	ChkNoise noise;
	//float freq = curr.x == 0 ? 0.1f : 0.1f / curr.x;
	uint octaves = 6;
	SimplexNoise noiseGen(0.004f, 1.f, 2.f, 0.5f);
	// noise generation
//#pragma omp parallel for
	for (int y = 0; y < NOISE_H; y++) {
		for (int z = 0; z < NOISE_W; z++) {
			for (int x = 0; x < NOISE_W; x++) {
				auto pos = Pos3i{x + curr.x, y + curr.y, z + curr.z};
				auto val = noiseGen.fractal(octaves, pos.x -1, pos.y -1, pos.z -1);
				val -= (pos.y - 8.f) / 24.f;
				//noiseAtPos(pos, noise) = val;
				noise.at(x + z * NOISE_W + y * (NOISE_W * NOISE_W)) = val;
			}
		}
	}
	return noise;
}

// for debugging purposes
ChkNoise genNoiseFill(Pos3i &curr) {
	ChkNoise noise;
	// noise generation
//#pragma omp parallel for
	for (int y = 0; y < NOISE_H; y++) {
		for (int z = 0; z < NOISE_W; z++) {
			for (int x = 0; x < NOISE_W; x++) {
				auto val = -1;
				val += (x * x  + y * y + z * z) / float(32 * 32 * 32);
				noise.at(x + z * NOISE_W + y * (NOISE_W * NOISE_W)) = val;
			}
		}
	}
	return noise;
}

// generates chunk data for the Chunk passed as an argument
void genChunk(Pos3i cPos, Chunk *chk) {
	chk->noise     = genNoise(cPos);
	chk->mesh      = genMesh(cPos, chk->noise);
	chk->indices   = chk->mesh.size();
	chk->generated = true;

	// unlock mutex
	chk->isUsed.unlock();
}

GLuint waterVBO;

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
		Vertex{{X,0,X},{0,1,0}, waterColor},
		Vertex{{-X,0,X},{0,1,0}, waterColor},
		Vertex{{-X,0,-X},{0,1,0}, waterColor},
		Vertex{{X,0,-X},{0,1,0}, waterColor},
		Vertex{{X,0,X},{0,1,0}, waterColor}
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

	for (int x = -RENDER_DIST;   x <= RENDER_DIST;   x++) {
	for (int y = -RENDER_DIST/2; y <= RENDER_DIST/2; y++) {
	for (int z = -RENDER_DIST;   z <= RENDER_DIST;   z++) {
		Pos3i cPos = {
			.x = (currChkPos.x + x) * CHK_SIZE,
			.y = (currChkPos.y + y) * CHK_SIZE,
			.z = (currChkPos.z + z) * CHK_SIZE
		};
		auto hash = chunkHasher(cPos);
		Chunk &chk = chunkMap[hash];

		if (chk.isUsed.try_lock()) {
			if (chk.generated) {
				// submitting the mesh to gpu
				if (chk.mesh.size() > 0) {
					glBindBuffer(GL_ARRAY_BUFFER, chk.vbo);
					glBufferData(GL_ARRAY_BUFFER, chk.indices * sizeof(Vertex),
         					(void*)chk.mesh.data(), GL_STATIC_DRAW);
         			chk.mesh.clear();
         			// lag protection
					lastFrame = glfwGetTime();
				}
				drawChunk(chk, sp);
				chk.isUsed.unlock();
			} else {
				// vbo for the chunk
				GLuint vbo;
				glGenBuffers(1, &vbo);
				chk.vbo = vbo;
				lastFrame = glfwGetTime();

				//LOG_INFO("Generating chunk [%i %i %i]...\nHash: %lu", cPos.x / CHK_SIZE, 
					//cPos.y / CHK_SIZE, cPos.z / CHK_SIZE, hash);
				std::thread thread(genChunk, cPos, &chk);
				thread.detach();
				//enqChunk(cPos, &chk);
			}
		}
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
}
