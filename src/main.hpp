#pragma once
#include <array>
#include <cinttypes>
#include <glm/glm.hpp>
#include <memory>
#include <vector>

using std::array;
using std::vector;

using byte = uint8_t;
using uint = unsigned;

// unique_ptr alias
template<typename... T>
using Uniq_Ptr = std::unique_ptr<T...>;

constexpr const char* WINDOW_TITLE = "GLTRY";
constexpr uint WINDOW_W = 800;
constexpr uint WINDOW_H = 600;
constexpr float ASPECT_RATIO = float(WINDOW_W) / WINDOW_H;
constexpr float FOV = 90.f;

constexpr int RENDER_DIST = 5;

// size of one chunk
constexpr int CHK_SIZE = 32;
constexpr int MAP_W = CHK_SIZE;
constexpr int MAP_H = CHK_SIZE;

// render distance
constexpr float Z_NEAR = 0.1f;
constexpr float Z_FAR = (RENDER_DIST + 3) * CHK_SIZE;

// size of the noise array 
constexpr uint NOISE_W = MAP_W + 3;
constexpr uint NOISE_H = MAP_H + 3;

// used for hashing
using ChkHash = uint64_t;

using ChkNoise = array<float, NOISE_W * NOISE_W * NOISE_H>;
// biome map (temporary // more biomes tbd)
using SmoothMap = array<float, NOISE_W * NOISE_W>;

// iso
constexpr float THRESHOLD = 0.f;

// world space coordinates
using Vec3f = glm::vec3;

struct Vertex {
	glm::vec3 pos;
	glm::vec3 norm;
	glm::vec4 color;
};

using Mesh = std::vector<Vertex>;

void meshPush(Mesh &mesh, Mesh &val);

