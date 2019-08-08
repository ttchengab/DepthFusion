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
#include "eigenlib/Eigen/Dense"
#include "eigenlib/Eigen/Core"
#include "eigenlib/Eigen/Geometry"
#include "PreProcessing.cpp"
#include "SurfaceMeasurement.cpp"
#include "TSDFConstruction.cpp"
#include "RayCasting.cpp"

using namespace Eigen;
using namespace std;

vector<float> fuseTSDF(vector<float> voxelsTSDF, vector<float> newTSDF){
    for(int y = 0; y < voxelParams.voxNumy; y++){
        for(int x = 0; x < voxelParams.voxNumx; x++){
            for(int z = 0; z < voxelParams.voxNumz; z++){
                int i = z+x*voxelParams.voxNumz+y*voxelParams.voxNumz*voxelParams.voxNumx;
                voxelsTSDF[i] = (voxelsTSDF[i]*totalWeight+newTSDF[i])/(totalWeight+1);
            }
        }
    }
    totalWeight++;
    cout<<"fusion finished"<<endl;
    return voxelsTSDF;
}

icpPose getICPPose(icpPose prev, float* depthMap){
    //float x = (depthMap-cx)*z00*fxinv;
    MatrixXf Tgkprev(3, 4);
    Tgkprev << 1, prev.alpha, -1*prev.gamma, prev.tx,
               -1*prev.alpha, 1, prev.betta, prev.tx,
               prev.gamma, -1*prev.betta, 1, prev.tz;
    MatrixXf curTransform = Tgkprev.inverse()*Tgkprev;
    vector<Point> points = computePoints(depthMap);
    for(int i = 0; i < points.size(); i++){
        VectorXf originalPt(4);
        originalPt << points[i].x,
                        points[i].y,
                        points[i].z,
                        1;
        VectorXf result = curTransform*originalPt;
        float onScreenx = result(0,0)*fx/result(2,0);
        float onScreeny = result(1,0)*fy/result(2,0);
        if(onScreenx < pixelWidth && onScreeny < pixelHeight && onScreenx >= 0 && onScreeny >= 0){
            float curVoxDepth = interpolation(onScreenx, onScreeny, depthMap, pixelWidth);
            Point uPt;
            uPt.z = curVoxDepth;
            uPt.x = (onScreenx-cx)*uPt.z*fxinv;
            uPt.y = (onScreeny-cy)*uPt.z*fyinv;
            points.push_back(uPt);
        }
    }
    //MatrixXf Tgkprev()
    //worldU = K * Tgkprev * u convert bak to 3d
    //computePoints(worldU), computeNormal(v)
    //prevICPPose
}

int main(){
    vector<string> poses = getGroundTruthPose("groundtruth.txt");
    vector<string> imageNames = getImageNames("depth.txt");
    string firstImg = imageNames[0]+".png.bin";
    float *depthMap;
    bool testmode = 1;
    depthMap = getDepthMap(firstImg.c_str());
    surfaceMeasurement(depthMap);
    QuatPose pose(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    MatrixXf tfo = qtransformation(Vector3f(pose.tx, pose.ty, pose.ty), pose.qx, pose.qy, pose.qz, pose.qw);
    vector<float> voxelsTSDF = createTSDF(depthMap, tfo.inverse());
    totalWeight++;
    if(testmode == 0){
        for(int i = 1; i < 5; i ++){
                cout<<"Image No.";
                cout<<i<<endl;
                string img = imageNames[i]+".png.bin";
                float imgPose = stof(imageNames[i].substr(7,16));
                depthMap = getDepthMap(img.c_str());
                QuatPose newPose = getNearestPose(imgPose, poses);
                MatrixXf calibrationMtx = qtransformation(Vector3f(clibRot.tx, clibRot.ty, clibRot.tz), clibRot.qx, clibRot.qy, clibRot.qz, clibRot.qw);
                MatrixXf ntf = qtransformation(Vector3f(newPose.tx, newPose.ty, newPose.tz), newPose.qx, newPose.qy, newPose.qz, newPose.qw);
                MatrixXf tf = calibrationMtx.inverse()*ntf;
                vector<float> newTSDF = createTSDF(depthMap, tf.inverse());
                voxelsTSDF = fuseTSDF(voxelsTSDF, newTSDF);
    }
    }
    //validateTSDF(voxelsTSDF); //Used to validate TSDF
    MatrixXf viewAngle = transformation(Vector3f(0.0,0.0,0.0), 0.0, 0.0, 0.0);
    rayCasting(voxelsTSDF, viewAngle, 0.0, 0.0, 0.0);
    return 0;
}
