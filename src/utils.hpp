#pragma once

#include "main.hpp"

template<typename T>
struct Position {
	T x;
	T y;
	T z;

	Position<T> operator + (Position<T> &offset) const {
		Position<T> npos;
		npos.x = x + offset.x;
		npos.y = y + offset.y;
		npos.z = z + offset.z;
		return npos;
	}
	Position<T> operator / (int div) const {
		Position<T> npos;
		npos.x = x / div;
		npos.y = y / div;
		npos.z = z / div;
		return npos;
	}
	Position<T> operator * (int div) const {
		Position<T> npos;
		npos.x = x * div;
		npos.y = y * div;
		npos.z = z * div;
		return npos;
	}
	template<typename Y>
	Position<T> offset(Y off_x, Y off_y, Y off_z) const {
		Position<T> npos;
		npos.x = x + T(off_x);
		npos.y = y + T(off_y);
		npos.z = z + T(off_z);
		return npos;
	}
};

using Pos3i = Position<int>;

float trilinearInterp(array<float, 8> &v, float &x, float &y, float &z);
float bilinearInterp(array<float, 4> &v, float &x, float &y);
float flerp(float v0, float v1, float t);
