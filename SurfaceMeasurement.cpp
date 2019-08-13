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
#include "HelperFunctions.cpp"

using namespace std;

vector<Point> computePoints(float *depthMap){
    vector<Point> points;
    Point currentPt;
    for(int i = 0; i < pixelWidth*pixelHeight; i++){
        float currentDepth = depthMap[i];
        if(currentDepth != 0){
            currentPt.z = currentDepth;
            float x = i%pixelWidth, y = i/pixelWidth;
            currentPt.x = (x-cx)*currentPt.z*fxinv;
            currentPt.y = (y-cy)*currentPt.z*fyinv;
            points.push_back(currentPt);
        }
    }
    return points;
}

vector<Point> computeNormal(vector<Point> vertexMap){
    vector<Point> normalVectors;
    for(int i = 0; i < vertexMap.size(); i++){
        int x = i%pixelWidth, y = i/pixelWidth;
        Point u, v;
        if((x+1)%pixelWidth != 0){
            u.x = vertexMap[i+1].x-vertexMap[i].x;
            u.y = vertexMap[i+1].y-vertexMap[i].y;
            u.z = vertexMap[i+1].z-vertexMap[i].z;
        }
        else{
            u.x = vertexMap[i].x-vertexMap[i-1].x;
            u.y = vertexMap[i].y-vertexMap[i-1].y;
            u.z = vertexMap[i].z-vertexMap[i-1].z;
        }
        if((y+1)<320){
            v.x = vertexMap[i+pixelWidth].x-vertexMap[i].x;
            v.y = vertexMap[i+pixelWidth].y-vertexMap[i].y;
            v.z = vertexMap[i+pixelWidth].z-vertexMap[i].z;    
        }
        else{
            v.x = vertexMap[i].x-vertexMap[i-pixelWidth].x;
            v.y = vertexMap[i].y-vertexMap[i-pixelWidth].y;
            v.z = vertexMap[i].z-vertexMap[i-pixelWidth].z;    
        }
        Point normalVec;
        normalVec.x = u.y*v.z-v.y*u.z;
        normalVec.y = v.x*u.z-u.x*v.z;
        normalVec.z = u.x*v.y-v.x*u.y;
        normalVectors.push_back(normalVec);
    }
    return normalVectors;
}
void surfaceMeasurement(float *depthMap){
    vector<Point> vertexMap = computePoints(depthMap);
    writeToPly(vertexMap, "testMesh.ply");
    vector<Point> normalMap = computeNormal(vertexMap);
    //test normal map
    cout<<"surface measurement finished"<<endl;
}