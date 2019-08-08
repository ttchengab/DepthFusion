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

using namespace Eigen;
using namespace std;

void rayCasting(const vector<float> &voxelsTSDF, MatrixXf rotation, float translationx, float translationy, float translationz){
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
            if( locxv < voxelParams.voxNumx && locxv >= 0 && locyv < voxelParams.voxNumy && locyv >= 0 && loczv < voxelParams.voxNumz && loczv >= 0){
                Point currentPt;
                currentPt.x = locx;
                currentPt.y = locy;
                currentPt.z = locz;
                points.push_back(currentPt);
            }
        }
    }
    cout<<"Raycast done";
    //writeToPly(points, "rayCastMesh.txt");
    writeToPly(points, "rayCastMesh256.ply");
}