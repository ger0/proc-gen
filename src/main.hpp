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

