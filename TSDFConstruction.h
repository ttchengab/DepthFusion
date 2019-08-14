#ifndef TSDF_H
#define TSDF_H
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
#ifndef EIGEN
#define EIGEN
 // EIGEN
#include "eigenlib/Eigen/Dense"
#include "eigenlib/Eigen/Core"
#include "eigenlib/Eigen/Geometry"
#endif
#include "HelperFunctions.h"
#include "ParamsInit.h"

using namespace Eigen;
using namespace std;

MatrixXf transformation(const Vector3f &translation, float rotx, float roty, float rotz);
MatrixXf qtransformation(const Vector3f &translation, float qx, float qy, float qz, float qw);
void validateTSDF(const vector<Voxel>& voxelsTSDF);
void createTSDF(float* depthMap, MatrixXf tf, vector<Voxel>& voxelsTSDF);

#endif