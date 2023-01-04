#pragma once
#include "main.hpp"

// Vertex layout:
//
//            7             6
//            +-------------+          
//          / |           / |          
//        /   |         /   |        
//    3 +-----+-------+  2  |          
//      |   4 +-------+-----+ 5        
//      |   /         |   /          
//      | /           | /                  
//    0 +-------------+ 1            
//
constexpr array<Vec3f, 8> marching_cubes_offsets = {
    Vec3f{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
    {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}
};

Vec3f vert_interp(const Vec3f &p1, const Vec3f &p2, float valp1, float valp2);

int polygonise(array<float, 8> &grid, vector<Vertex> &vert_vec);
