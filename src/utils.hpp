#pragma once

#include "main.hpp"

template<typename T>
struct Position {
	T x;
	T y;
	T z;

	Position<T> operator + (Position<T> &offset) const;
	Position<T> operator / (int div) const;
	Position<T> operator * (int div) const;
	template<typename Y>
	Position<T> offset(Y off_x, Y off_y, Y off_z) const;
};

using Pos3i = Position<int>;

float trilinearInterp(array<float, 8> &v, float &x, float &y, float &z);
float bilinearInterp(array<float, 4> &v, float &x, float &y);
float lerp(float v0, float v1, float t);
