#ifndef HELPER_H
#define HELPER_H
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <vector> // for 2D vector
#include <bitset>
#include <cmath>
#include <cstring>
#include <chrono>
#include <fstream>
#include "ParamsInit.h"

using namespace std;

QuatPose getNearestPose(float imageTime, const vector<string>& poses);
void writeToPly(vector<Point> points, const char* fileName);
void writeToPlyNorm(vector<Point> points, vector<Point> normals, const char* fileName);
float interpolation(float x, float y, const float *map, int width);
float interpolation(float x, float y, float z,  const vector<Voxel> &map3d, int depth, int width);

#endif
