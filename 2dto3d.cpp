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
using namespace std::chrono;
int totalWeight = 1;

bool pointRayCast(Point & currentPt, const vector<Voxel>& voxelsTSDF, float x, float y, float minimumDepth, MatrixXf rotation, float translationx, float translationy, float translationz){
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
    VectorXf poseRotate(3);
    poseRotate << stepx,
                  stepy,
                  stepz;
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
        if(INTERPOLATE){
            curTSDF = interpolation(locxv, locyv, loczv, voxelsTSDF, voxelParams.voxNumz, voxelParams.voxNumx);
        }
        else{
            curTSDF = voxelsTSDF[int(locyv)*voxelParams.voxNumx*voxelParams.voxNumz+int(locxv)*voxelParams.voxNumz+int(loczv)].value;
        }
        locxv = (locx+voxelParams.voxPhysLength/2)/voxelParams.voxSize;
        locyv = (locy+voxelParams.voxPhysWidth/2)/voxelParams.voxSize;
        loczv = (locz)/voxelParams.voxSize;
        locx += stepx;
        locy += stepy;
        locz += stepz;
    }
    if(curTSDF > -1 && locxv < voxelParams.voxNumx && locxv >= 0 && locyv < voxelParams.voxNumy && locyv >= 0 && loczv < voxelParams.voxNumz && loczv >= 0){
        currentPt.x = locx;
        currentPt.y = locy;
        currentPt.z = locz;
        return true;
        //points.push_back(currentPt);
    }
    return false;
}
void rayCasting(const vector<Voxel>& voxelsTSDF, MatrixXf rotation, float translationx, float translationy, float translationz){
    float minimumDepth = 0;
    vector<Point>points;
    MatrixXf rotationMtx = rotation.block<3,3>(0,0);
    for(int y = 0; y < pixelHeight; y++){
        for(int x = 0; x < pixelWidth; x++){
            int i = y*pixelWidth+x;
            Point currentPt;
            bool pointExist = pointRayCast(currentPt, voxelsTSDF, x, y, 0, rotationMtx, translationx, translationy, translationz);
            if(pointExist) points.push_back(currentPt);
        }
    }
    cout<<"Raycast done";
    //writeToPly(points, "rayCastMesh.txt");
    writeToPly(points, "meshes/FYPRayCastMesh.ply");
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

void computeProjPix(const Point& point, float& ux, float& uy, const MatrixXf& f2fTransform){
  VectorXf originalPt(4);
  originalPt << point.x,
              point.y,
              point.z,
              1;
  VectorXf result = f2fTransform*originalPt;
  ux = result(0,0)*fx/result(2,0)+cx;
  uy = result(1,0)*fy/result(2,0)+cy;
}
void computeNormal(Point& normal, float x, float y, float z, const vector<Voxel> &voxelsTSDF){
    // cout<<"Normals:"<<endl;
    float fx = voxelsTSDF[int(y)*voxelParams.voxNumx*voxelParams.voxNumz+int(x+1)*voxelParams.voxNumz+int(z)].value;
    fx -= voxelsTSDF[floor(y)*voxelParams.voxNumx*voxelParams.voxNumz+int(x)*voxelParams.voxNumz+int(z)].value;
    float fy = voxelsTSDF[int(y+1)*voxelParams.voxNumx*voxelParams.voxNumz+int(x)*voxelParams.voxNumz+int(z)].value;
    fy -= voxelsTSDF[int(y)*voxelParams.voxNumx*voxelParams.voxNumz+int(x)*voxelParams.voxNumz+int(z)].value;
    float fz = voxelsTSDF[int(y)*voxelParams.voxNumx*voxelParams.voxNumz+int(x)*voxelParams.voxNumz+int(z+1)].value;
    fz -= voxelsTSDF[int(y)*voxelParams.voxNumx*voxelParams.voxNumz+int(x)*voxelParams.voxNumz+int(z)].value;
    float normalLength = sqrt(pow(fx, 2)+pow(fy, 2)+pow(fz, 2));
    // if(normalLength == 0){
    //   cout<<"here"<<endl;
    //   cout<<fx<<endl;
    //   cout<<fy<<endl;
    //   cout<<fz<<endl;
    // }
    normal.x = fx/normalLength;
    normal.y = fy/normalLength;
    normal.z = fz/normalLength;
}
bool checkOmega(Point Vcur, Point Vglob, Point Ncur, Point Nglob, MatrixXf& curTransform, float Ed, float Etheta){
    // cout<<"norm global"<<endl;
    // cout<<Nglob.x<<endl;
    // cout<<Nglob.y<<endl;
    // cout<<Nglob.z<<endl;
    // cout<<"norm current"<<endl;
    // cout<<Ncur.x<<endl;
    // cout<<Ncur.y<<endl;
    // cout<<Ncur.z<<endl;
    VectorXf temp(4);
    temp(0,0) = Vcur.x;
    temp(1,0) = Vcur.y;
    temp(2,0) = Vcur.z;
    temp(3,0) = 1;
    VectorXf Vku = curTransform*temp;
    float length = sqrt(pow(Vku(0,0)-Vglob.x, 2)+pow(Vku(1,0)-Vglob.y, 2)+pow(Vku(2,0)-Vglob.z, 2));
    VectorXf Nku(3);
    Nku << Ncur.x,
           Ncur.y,
           Ncur.z;
    MatrixXf Rgk = curTransform.block<3,3>(0,0);
    MatrixXf RgkN = Rgk*Nku;
    float normLength1 = sqrt(pow(RgkN(0,0), 2)+pow(RgkN(1,0), 2)+pow(RgkN(2,0), 2));
    float normLength2 = sqrt(pow(Nglob.x, 2)+pow(Nglob.y, 2)+pow(Nglob.z, 2));
    float angle = acos((Nglob.x*RgkN(0,0)+Nglob.y*RgkN(1,0)+Nglob.z*RgkN(2,0))/(normLength1*normLength2));
    // cout<<"Length: "<<length<<endl;
    // cout<<"Angle: "<<angle<<endl;
    // cout<<"vec1"<<endl;
    // cout<<Vcur.x<<endl;
    // cout<<Vcur.y<<endl;
    // cout<<Vcur.z<<endl;
    //cout<<curTransform;
    return (length <= Ed && angle <= Etheta);

}
void sumAt(const vector<Point>& Vk, const vector<Point>& Nglob, MatrixXf& At){
    for(int i = 0; i < Vk.size(); i++){
        MatrixXf G(3,6), Nprev(3,1);
        G << 0, -1*Vk[i].z, Vk[i].y, 1, 0, 0,
            Vk[i].z, 0, -1*Vk[i].x, 0, 1, 0,
            -1*Vk[i].y, Vk[i].x, 0, 0, 0, 1;
        Nprev << Nglob[i].x,
                Nglob[i].y,
                Nglob[i].z;
        MatrixXf newAt = G.transpose()*Nprev;
        if(i == 0){
          cout<<"G: "<<endl<<G<<endl;
          cout<<"Nprev: "<<endl<<Nprev<<endl;
          cout<<"At: "<<endl<<At<<endl;
          cout<<"newAt: "<<newAt<<endl;
        }
        At = At+newAt;
    }
}
void sumB(const vector<Point>& Vk, const vector<Point>& Nglob, const vector<Point>& Vglob, MatrixXf& B){
    for(int i = 0; i < Vk.size(); i++){
        MatrixXf Nprev(3,1), Vgk(3,1), Vprev(3,1);
        Nprev << Nglob[i].x,
                Nglob[i].y,
                Nglob[i].z;
        Vprev << Vglob[i].x,
                Vglob[i].y,
                Vglob[i].z;
        Vgk << Vk[i].x,
              Vk[i].y,
              Vk[i].z;
        cout<<"Vprev and Vgk"<<endl;
        cout<<Vprev<<endl;
        cout<<Vgk<<endl;
        MatrixXf newB = Nprev.transpose()*(Vprev-Vgk);
        B = B+newB;
    }
}
void updateTinc(VectorXf& x, MatrixXf& curTransform, MatrixXf& f2fTransform, const MatrixXf& initTransform){
    MatrixXf Tinc(3,4);
    Tinc << 1, x(2,0), -1*x(1,0), x(3,0),
            -1*x(2,0), 1, x(0,0), x(4,0),
            x(1,0), -1*x(0,0), 1, x(5,0);
    cout<<"curTransform"<<endl<<curTransform<<endl;
    cout<<"f2fTransform"<<endl<<f2fTransform<<endl;
    cout<<"initTransform"<<endl<<initTransform<<endl;
    curTransform.block<3,4>(0,0) = Tinc*curTransform.inverse();
    curTransform(3,0) = 0;
    curTransform(3,1) = 0;
    curTransform(3,2) = 0;
    curTransform(3,3) = 1;
    f2fTransform = initTransform.inverse()*curTransform;
}
MatrixXf getICPPose(MatrixXf initTransform, const vector<Point> vertexMap, const vector<Point> normalMap, const vector<Voxel> &voxelsTSDF){
    MatrixXf curTransform = initTransform;
    MatrixXf f2fTransform = initTransform.inverse()*curTransform;

    // for(int i = 0; i < normalMap.size(); i++){
    //   cout<<normalMap[i].x<<endl;
    //   cout<<normalMap[i].y<<endl;
    //   cout<<normalMap[i].z<<endl;
    // }
    // cout<<vertexMap.size()<<endl;
    // cout<<normalMap.size()<<endl;
    //return curTransform;
    for(int z = 1; z<3; z++){
      vector<Point> Vk, Nk, Vglob, Nglob;
      int check = 0;
      for(int i = 0; i < vertexMap.size(); i++){
          float projUx = 0;
          float projUy = 0;
          computeProjPix(vertexMap[i], projUx, projUy, f2fTransform);

          if(projUx < pixelWidth-1 && projUy < pixelHeight-1 && projUx >= 0 && projUy >= 0){
            //think about the transform before raycast
            //bool pointExist = pointRayCast()
              MatrixXf rotationMtx = curTransform.block<3,3>(0,0);
              float transX = curTransform(3,0);
              float transY = curTransform(3,1);
              float transZ = curTransform(3,2);
              Point raycastPt, raycastPtNxt;
              bool pointExist = pointRayCast(raycastPt, voxelsTSDF, projUx, projUy, 0, rotationMtx, transX, transY, transZ);
            //  bool nextPtExist = pointRayCast(raycastPt, voxe)
              //cout<<<<endl;

              if(pointExist){
                  // cout<<"x y z"<<endl;
                  // cout<<vertexMap[i].x<<" "<<vertexMap[i].y<<" "<<vertexMap[i].z<<endl;
                  // cout<<"rayx rayy rayz"<<endl;
                  // cout<<raycastPt.x<<" "<<raycastPt.y<<" "<<raycastPt.z<<endl;
                  Point raycastNorm;
                  float rayxv = (raycastPt.x+voxelParams.voxPhysLength/2)/voxelParams.voxSize;
                  float rayyv = (raycastPt.y+voxelParams.voxPhysWidth/2)/voxelParams.voxSize;
                  float rayzv = (raycastPt.z)/voxelParams.voxSize;
                  computeNormal(raycastNorm, rayxv, rayyv, rayzv, voxelsTSDF);
                  // cout<<normalMap[i].x<<" "<<normalMap[i].y<<" "<<normalMap[i].z<<endl;
                  // cout<<raycastNorm.x<<" "<<raycastNorm.y<<" "<<raycastNorm.z<<endl;
                  // cout<<"----------"<<endl;

                  if(checkOmega(vertexMap[i], raycastPt, normalMap[i],raycastNorm, curTransform, 0.1, 3.14/12)){
                      check = 1;
                      // cout<<"threshold passed"<<endl;
                      Vk.push_back(vertexMap[i]);
                      Vglob.push_back(raycastPt);
                      Nk.push_back(normalMap[i]);
                      Nglob.push_back(raycastNorm);
                  }
            }
          }
      }
      if(check == 0){
        continue;
      }
      MatrixXf At(6,1);
      At << 0,0,0,0,0,0;
      MatrixXf B(1,1);
      B << 0;

      sumAt(Vk, Nglob, At);
      sumB(Vk, Nglob, Vglob, B);
      cout<<"At"<<endl;
      cout<<At<<endl;
      MatrixXf AtA = At*At.transpose();
      MatrixXf Atb = At*B;
      VectorXf x = AtA.colPivHouseholderQr().solve(Atb);
      cout<<"B"<<endl<<B<<endl;
      cout<<"AtA"<<endl<<AtA<<endl;
      cout<<"AtB"<<endl<<Atb<<endl;
      cout<<"x"<<endl<<x<<endl;
      updateTinc(x, curTransform, f2fTransform, initTransform);
    }
    cout<<"icp finished"<<endl;
    return curTransform.inverse();
}

int main(){
    vector<string> poses = getGroundTruthPose("groundtruth.txt");
    // rgbd_dataset_freiburg1_teddy/
    // rgbd_dataset_freiburg1_plant/
    // rgbd_dataset_freiburg2_coke
    vector<string> imageNames = getImageNames("depth.txt");
    string root_dir = "binDepths/xyzOriginal/";
    // cout<<imageNames;
    string firstImg = root_dir + imageNames[0] + ".png.bin";
    float *depthMap;
    int testmode = 1;
    cout<<firstImg<<endl;
    depthMap = getDepthMap(firstImg.c_str());
    vector<Point> vertexMap, normalMap;
    QuatPose pose(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

    MatrixXf tfo = qtransformation(Vector3f(pose.tx, pose.ty, pose.ty), pose.qx, pose.qy, pose.qz, pose.qw);
    Voxel initVoxel;
    initVoxel.value = 1;
    initVoxel.weight = 0;
    surfaceMeasurement(depthMap, vertexMap, normalMap);
    vector<Voxel> voxelsTSDF(voxelParams.voxVolume, initVoxel);
    vector<Voxel> newTSDF(voxelParams.voxVolume, initVoxel);
    createTSDF(depthMap, tfo.inverse(), voxelsTSDF);
    if(testmode == 1){
        for(int i = 1; i < 30; i ++){
                cout<<"Image No.";
                cout<<i<<endl;
                string img = root_dir + imageNames[i]+".png.bin";
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
        for(int i = 1; i < 5; i ++){
                cout<<"Image No.";
                cout<<i<<endl;
                string img = root_dir + imageNames[i] + ".png.bin";
                float imgPose = stof(imageNames[i].substr(7,16));
                depthMap = getDepthMap(img.c_str());
                surfaceMeasurement(depthMap, vertexMap, normalMap);
                tf = getICPPose(tf.inverse(), vertexMap, normalMap, voxelsTSDF);
                createTSDF(depthMap, tf.inverse(), newTSDF);
                fuseTSDF(voxelsTSDF, newTSDF);
        }
    }

    validateTSDF(voxelsTSDF); //Used to validate TSDF
    MatrixXf viewAngle = transformation(Vector3f(0.0,0.0,0.0), 0.0, 0.0, 0.0);
    chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
    rayCasting(voxelsTSDF, viewAngle, 0.0, 0.0, 0.0);
    chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double> >(t2 - t1);
    cout << "Raycast took " << time_span.count() << " seconds."<<endl;
    return 0;
}
