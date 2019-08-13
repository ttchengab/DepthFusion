#ifndef PREP_H
#define PREP_H
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <vector> // for 2D vector
#include <bitset>
#include <cmath>
#include <cstring>
#include <fstream>
#include "ParamsInit.h"

using namespace std;
vector<string> getGroundTruthPose(const char* fileName);
vector<string> getImageNames(const char* fileName);
float *getDepthMap(const char* fileName);

#endif