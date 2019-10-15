#ifndef SURF_H
#define SURF_H

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
#include "HelperFunctions.h"

using namespace std;
void computePoints(float *depthMap, vector<Point>& points);
void computeNormal(const vector<Point>& vertexMap, vector<Point>& normalVectors);
void surfaceMeasurement(float *depthMap, vector<Point>& vertexMap, vector<Point>& normalMap);

#endif
