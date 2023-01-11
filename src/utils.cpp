#include <cstdio>
#include <array>

#include "log.hpp"
#include "utils.hpp"
#include "main.hpp"

float lerp(float v0, float v1, float t) {
	if (t > 1.f || t < 0.f) {
		LOG_WARN("BILINEAR INTERP PROBLEM! %f", t);
	}
  return (1 - t) * v0 + t * v1;
}

// 4 sampled points, x and y within [0..1]
float bilinearInterp(array<float, 4>& v, float &x, float &y) {
	// sampling 4 edges on the x axis
	auto ex1 = lerp(v[0], v[1], x);
	auto ex2 = lerp(v[2], v[3], x);
	// sampling 2 edges on the y axis
	return lerp(ex1, ex2, y);
}

// x, y, z within [0..1]
float trilinearInterp(array<float, 8> &v, float &x, float &y, float &z) {
  float result = 0;
  auto p1 = array<float, 4>{v[0], v[1], v[2], v[3]};
  auto p2 = array<float, 4>{v[4], v[5], v[6], v[7]};
  float ex1 = bilinearInterp(p1, x, y);
  float ex2 = bilinearInterp(p2, x, y);
  return lerp(ex1, ex2, z);
}
