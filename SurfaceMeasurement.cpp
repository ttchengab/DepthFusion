#include "SurfaceMeasurement.h"

void computePoints(float *depthMap, vector<Point>& points){
    //vector<Point> points;
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
}

void computeNormal(const vector<Point>& vertexMap, vector<Point>& normalVectors){
    //vector<Point> normalVectors;
    int count=0;
    vector<Point> vm, nm;
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
        if((y+1)<480){
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
        //cout<<normalVec.x<<" "<<normalVec.y<<" "<<normalVec.z<<endl;
        if(abs(normalVec.x)>0.1 || abs(normalVec.y)>0.1 || abs(normalVec.z)>0.1){
          count++;
          vm.push_back(vertexMap[i]);
          nm.push_back(normalVec);
        }
    }
    writeToPlyNorm(vm, nm, "justNorm");
    cout<<count<<endl;
}
void surfaceMeasurement(float *depthMap, vector<Point>& vertexMap, vector<Point>& normalMap){
    vertexMap.clear();
    normalMap.clear();
    computePoints(depthMap, vertexMap);
    computeNormal(vertexMap, normalMap);
    writeToPlyNorm(vertexMap, normalMap, "testMeshNorm.ply");
    //computeNormal(vertexMap, );
    //test normal map
    cout<<"surface measurement finished"<<endl;
}
