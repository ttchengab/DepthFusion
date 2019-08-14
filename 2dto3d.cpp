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

#include "ParamsInit.h"
#include "PreProcessing.h"
#include "SurfaceMeasurement.h"
#include "TSDFConstruction.h"
//#include "RayCasting.cpp"

using namespace Eigen;
using namespace std;
int totalWeight = 1;
void rayCasting(const vector<Voxel>& voxelsTSDF, MatrixXf rotation, float translationx, float translationy, float translationz){
    float minimumDepth = 0;
    vector<Point>points;
    for(int y = 0; y < pixelHeight; y++){
        for(int x = 0; x < pixelWidth; x++){
            int i = y*pixelWidth+x;
            Point pixPt;
            pixPt.z = 1.0;
            pixPt.x = (x-cx)*fxinv;
            pixPt.y = (y-cy)*fyinv;
            float totalLength = sqrt(pow(pixPt.x, 2)+pow(pixPt.y, 2)+pow(pixPt.z, 2));
            float normx = pixPt.x/totalLength, normy = pixPt.y/totalLength, normz = pixPt.z/totalLength;
            float startx = minimumDepth*normx, starty = minimumDepth*normy, startz = minimumDepth*normz;
            startx = startx - translationx;
            starty = starty - translationy;
            startz = startz - translationz;
            float stepx = normx*voxelParams.voxSize, stepy = normy*voxelParams.voxSize, stepz = normz*voxelParams.voxSize;
            VectorXf poseRotate(4);
            poseRotate << stepx,
                          stepy,
                          stepz,
                              1;
            VectorXf result = rotation*poseRotate;
            stepx = result(0,0);
            stepy = result(1,0);
            stepz = result(2,0);
            float locx = startx, locy = starty, locz = startz, locxv = 0, locyv = 0, loczv = 0;
            locxv = (locx+voxelParams.voxPhysLength/2)/voxelParams.voxSize;
            locyv = (locy+voxelParams.voxPhysWidth/2)/voxelParams.voxSize;
            loczv = (locz)/voxelParams.voxSize;
            
            float curTSDF = interpolation(locxv, locyv, loczv, voxelsTSDF, voxelParams.voxNumz, voxelParams.voxNumx);
            //locxv < voxelParams.voxNumx && locxv >=0 && locyv < voxelParams.voxNumy && locxv >=0 && loczv < voxelParams.voxNumz && locxv >=0
            while(curTSDF > 0 && locxv < voxelParams.voxNumx && locxv >=0 && locyv < voxelParams.voxNumy && locyv >= 0 && loczv < voxelParams.voxNumz && loczv >= 0){
                //cout<<locxv<<" "<<locyv<<" "<<loczv<<endl;
                curTSDF = interpolation(locxv, locyv, loczv, voxelsTSDF, voxelParams.voxNumz, voxelParams.voxNumx);
                locxv = (locx+voxelParams.voxPhysLength/2)/voxelParams.voxSize;
                locyv = (locy+voxelParams.voxPhysWidth/2)/voxelParams.voxSize;
                loczv = (locz)/voxelParams.voxSize;
                locx += stepx;
                locy += stepy;
                locz += stepz;
            }
            if(locxv < voxelParams.voxNumx && locxv >= 0 && locyv < voxelParams.voxNumy && locyv >= 0 && loczv < voxelParams.voxNumz && loczv >= 0){
                Point currentPt;
                //cout<<locx<<" "<<locy<<" "<<locz<<endl;
                currentPt.x = locx;
                currentPt.y = locy;
                currentPt.z = locz;
                points.push_back(currentPt);
            }
        }
    }
    cout<<"Raycast done";
    //writeToPly(points, "rayCastMesh.txt");
    writeToPly(points, "rayCastMesh0813.ply");
}

void fuseTSDF(vector<Voxel>& voxelsTSDF, const vector<Voxel>& newTSDF){
    for(int y = 0; y < voxelParams.voxNumy; y++){
        for(int x = 0; x < voxelParams.voxNumx; x++){
            for(int z = 0; z < voxelParams.voxNumz; z++){
                int i = z+x*voxelParams.voxNumz+y*voxelParams.voxNumz*voxelParams.voxNumx;
                float totalWeight = voxelsTSDF[i].weight+newTSDF[i].weight;
                if(totalWeight == 0) continue;
                voxelsTSDF[i].value = (voxelsTSDF[i].value*voxelsTSDF[i].weight+newTSDF[i].value*newTSDF[i].weight)/totalWeight;
                voxelsTSDF[i].weight = totalWeight;
            }
        }
    }
    cout<<"fusion finished"<<endl;
}

MatrixXf getICPPose(MatrixXf initTransform, float* depthMap, const vector<Voxel> &voxelsTSDF){
    //float x = (depthMap-cx)*z00*fxinv;
    MatrixXf f2fTransform = initTransform.inverse()*initTransform;
    MatrixXf curTransform = initTransform;
    vector<Point> points;
    computePoints(depthMap, points);
    vector<Point> normals;
    computeNormal(points, normals);
    MatrixXf AtA, b;
    MatrixXf Tinc(3, 4);
    for(int z = 1; z < 5; z++){
        vector<Point> prevVPts, globalVkPts, prevNVecs;
        for(int i = 0; i < points.size(); i++){
            VectorXf originalPt(4);
            originalPt << points[i].x,
                        points[i].y,
                        points[i].z,
                        1;
            VectorXf result = f2fTransform*originalPt;
            float ux = result(0,0)*fx/result(2,0);
            float uy = result(1,0)*fy/result(2,0);
            //if(i >= 179817) cout<<i<<endl;
            if(ux < pixelWidth && uy < pixelHeight && ux >= 0 && uy >= 0){
                //cout<<ux<<" "<<uy<<endl;
                Point pixPt;
                pixPt.z = 1.0;
                pixPt.x = (ux-cx)*fxinv;
                pixPt.y = (uy-cy)*fyinv;
                float totalLength = sqrt(pow(pixPt.x, 2)+pow(pixPt.y, 2)+pow(pixPt.z, 2));
                float normx = pixPt.x/totalLength, normy = pixPt.y/totalLength, normz = pixPt.z/totalLength;
                float startx = 0, starty = 0, startz = 0;
                float stepx = normx*voxelParams.voxSize, stepy = normy*voxelParams.voxSize, stepz = normz*voxelParams.voxSize;
                float locx = startx, locy = starty, locz = startz, locxv = 0, locyv = 0, loczv = 0;
                locxv = (locx+voxelParams.voxPhysLength/2)/voxelParams.voxSize;
                locyv = (locy+voxelParams.voxPhysWidth/2)/voxelParams.voxSize;
                loczv = (locz)/voxelParams.voxSize;
                float curTSDF = interpolation(locxv, locyv, loczv, voxelsTSDF, voxelParams.voxNumz, voxelParams.voxNumx);
                while(curTSDF > 0 && locxv < voxelParams.voxNumx && locxv >=0 && locyv < voxelParams.voxNumy && locyv >= 0 && loczv < voxelParams.voxNumz && loczv >= 0){
                    //cout<<locxv<<" "<<locyv<<" "<<loczv<<endl;
                    curTSDF = interpolation(locxv, locyv, loczv, voxelsTSDF, voxelParams.voxNumz, voxelParams.voxNumx);
                    locxv = (locx+voxelParams.voxPhysLength/2)/voxelParams.voxSize;
                    locyv = (locy+voxelParams.voxPhysWidth/2)/voxelParams.voxSize;
                    loczv = (locz)/voxelParams.voxSize;
                    locx += stepx;
                    locy += stepy;
                    locz += stepz;
                }


                if(locxv < voxelParams.voxNumx-1 && locxv > 0 && locyv < voxelParams.voxNumy-1 && locyv > 0 && loczv < voxelParams.voxNumz-1 && loczv > 0){
                    //compute correspondences threshold
                    Point globalPt;
                    globalPt.x = locx;
                    globalPt.y = locy;
                    globalPt.z = locz;
                    //globalPt = Vgk(u^), originalPt = V(u), Tzgk = curTransform
                    VectorXf TgkVku = curTransform*originalPt;
                    float vertexDistance = pow(TgkVku(0,0)-globalPt.x, 2)+pow(TgkVku(1,0)-globalPt.y, 2)+pow(TgkVku(2,0)-globalPt.z, 2);
                    //compute Gradient
                    float fz = voxelsTSDF[int(locyv)*voxelParams.voxNumx*voxelParams.voxNumz+int(locxv)*voxelParams.voxNumz+int(loczv)+1].value;
                    fz-= voxelsTSDF[int(locyv)*voxelParams.voxNumx*voxelParams.voxNumz+int(locxv)*voxelParams.voxNumz+int(loczv)-1].value;
                    float fx = voxelsTSDF[int(locyv)*voxelParams.voxNumx*voxelParams.voxNumz+int(locxv)*(voxelParams.voxNumz+1)+int(loczv)].value;
                    fx -= voxelsTSDF[int(locyv)*voxelParams.voxNumx*voxelParams.voxNumz+int(locxv)*(voxelParams.voxNumz-1)+int(loczv)].value;
                    float fy = voxelsTSDF[int(locyv)*(voxelParams.voxNumx*voxelParams.voxNumz+1)+int(locxv)*voxelParams.voxNumz+int(loczv)].value;
                    fy-= voxelsTSDF[int(locyv)*(voxelParams.voxNumx*voxelParams.voxNumz-1)+int(locxv)*voxelParams.voxNumz+int(loczv)].value;
                    float normLength1 = sqrt(pow(fx, 2)+pow(fy, 2)+pow(fz, 2));
                    Point globalNormal;
                    globalNormal.x = fx/normLength1;
                    globalNormal.y = fy/normLength1;
                    globalNormal.z = fz/normLength1;
                    MatrixXf Rgk = curTransform.block<3,3>(0,0);
                    VectorXf Nku(3);
                    Nku << normals[i].x,
                           normals[i].y,
                           normals[i].z;
                    MatrixXf RgkN = Rgk*Nku;
                    float normLength2 = sqrt(pow(RgkN(0,0), 2)+pow(RgkN(1,0), 2)+pow(RgkN(2,0), 2));
                    float cos = (globalNormal.x*RgkN(0,0)+globalNormal.y*RgkN(1,0)+globalNormal.z*RgkN(2,0))/(normLength1*normLength2);
                    cout<<cos<<" "<<vertexDistance<<endl;
                    if(cos > epislonAng && vertexDistance<epislond){
                        cout<<"hi"<<endl;
                        prevVPts.push_back(globalPt);
                        prevNVecs.push_back(globalNormal);
                        Point globalPt;
                        globalPt.x = result(0,0);
                        globalPt.y = result(1,0);
                        globalPt.z = result(2,0);
                        globalVkPts.push_back(globalPt);
                    }

                }
            }
        }
        cout<<prevVPts.size()<<endl;
        //
        //cout<<"hi"<<endl;
        for(int i = 0; i < globalVkPts.size(); i++){
            MatrixXf G(3,6), Nprev(3,1), Vprev(3,1), Vgk(3,1);
            G << 0, -1*globalVkPts[i].z, globalVkPts[i].y, 1, 0, 0,
                globalVkPts[i].z, 0, -1*globalVkPts[i].x, 0, 1, 0,
                -1*globalVkPts[i].y, globalVkPts[i].x, 0, 0, 0, 1;
            Nprev << prevNVecs[i].x,
                    prevNVecs[i].y,
                    prevNVecs[i].z;
            Vprev << prevVPts[i].x,
                    prevVPts[i].y,
                    prevVPts[i].z;
            Vgk << globalVkPts[i].x,
                globalVkPts[i].y,
                globalVkPts[i].z;
            MatrixXf newAt = G.transpose()*Nprev;
            MatrixXf newAta = newAt*newAt.transpose();
            MatrixXf newB = Nprev.transpose()*(Vprev-Vgk);
            AtA += newAta;
            b += newB;
        }
        VectorXf x = AtA.colPivHouseholderQr().solve(b);

        Tinc << 1, x(2,0), -1*x(1,0), x(3,0),
                -1*x(2,0), 1, x(0,0), x(4,0),
                x(1,0), -1*x(0,0), 1, x(5,0);
        curTransform = Tinc*curTransform.inverse();
        f2fTransform = initTransform*curTransform;
    }
    cout<<"icp finished"<<endl;
    return curTransform;
}

int main(){
    vector<string> poses = getGroundTruthPose("groundtruth.txt");
    vector<string> imageNames = getImageNames("depth.txt");
    string firstImg = imageNames[0]+".png.bin";
    float *depthMap;
    int testmode = 1;
    depthMap = getDepthMap(firstImg.c_str());
    surfaceMeasurement(depthMap);
    QuatPose pose(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    MatrixXf tfo = qtransformation(Vector3f(pose.tx, pose.ty, pose.ty), pose.qx, pose.qy, pose.qz, pose.qw);
    Voxel initVoxel;
    initVoxel.value = 1;
    initVoxel.weight = 0;
    vector<Voxel> voxelsTSDF(voxelParams.voxVolume, initVoxel);
    vector<Voxel> newTSDF(voxelParams.voxVolume, initVoxel);
    createTSDF(depthMap, tfo.inverse(), voxelsTSDF);
    if(testmode == 1){
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
                createTSDF(depthMap, tf.inverse(), newTSDF);
                fuseTSDF(voxelsTSDF, newTSDF);
        }
    }
    else if(testmode == 2){
        MatrixXf tf = tfo;
        for(int i = 1; i < 2; i ++){
                cout<<"Image No.";
                cout<<i<<endl;
                string img = imageNames[i]+".png.bin";
                float imgPose = stof(imageNames[i].substr(7,16));
                depthMap = getDepthMap(img.c_str());
                tf = getICPPose(tf, depthMap, voxelsTSDF);
                createTSDF(depthMap, tf.inverse(), newTSDF);
                fuseTSDF(voxelsTSDF, newTSDF);
        }
    }
    
    validateTSDF(voxelsTSDF); //Used to validate TSDF
    MatrixXf viewAngle = transformation(Vector3f(0.0,0.0,0.0), 0.0, 0.0, 0.0);
    rayCasting(voxelsTSDF, viewAngle, 0.0, 0.0, 0.0);
    return 0;
}
