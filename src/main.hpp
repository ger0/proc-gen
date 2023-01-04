#pragma once
#include <array>
#include <cinttypes>
#include <glm/glm.hpp>
#include <vector>

using std::array;
using byte = uint8_t;


using byte = unsigned char;
using uint = unsigned;

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
constexpr int CHK_SIZE = 32;
constexpr int MAP_W = CHK_SIZE;
constexpr int MAP_H = CHK_SIZE;


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
using Vec3f = glm::vec3;

struct Vertex {
	glm::vec3 pos;
	glm::vec3 norm;
};

