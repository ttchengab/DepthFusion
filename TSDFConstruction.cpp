#include "TSDFConstruction.h"

MatrixXf transformation(const Vector3f &translation, float rotx, float roty, float rotz){
    Transform<float, 3, Affine> t = Transform<float, 3, Affine>::Identity();
    t.translate(translation);
    t.rotate(AngleAxisf(rotx, Vector3f::UnitX()));
    t.rotate(AngleAxisf(roty, Vector3f::UnitY()));
    t.rotate(AngleAxisf(rotz, Vector3f::UnitZ()));
    return t.inverse().matrix();
}

MatrixXf qtransformation(const Vector3f &translation, float qx, float qy, float qz, float qw){
    Transform<float, 3, Affine> t = Transform<float, 3, Affine>::Identity();
    t.translate(translation);
    Quaternionf q;
    q.x() = qx;
    q.y() = qy;
    q.z() = qz;
    q.w() = qw;
    MatrixXf rotMatrix = q.normalized().toRotationMatrix();
    MatrixXf transformMatrix = t.matrix();
    transformMatrix.block<3,3>(0,0) = rotMatrix;
    return transformMatrix;
}

void validateTSDF(const vector<Voxel>& voxelsTSDF){
    vector<Point> points;
    Point currentPt;
    for(int y = 0; y < voxelParams.voxNumy; y++){
        for(int x = 0; x < voxelParams.voxNumx; x++){
            for(int z = 0; z < voxelParams.voxNumz; z++){
                int i = z+x*voxelParams.voxNumz+y*voxelParams.voxNumz*voxelParams.voxNumx;
                if(fabs(voxelsTSDF[i].value)<0.2){
                    currentPt.x = voxelParams.voxSize*(x+1)-voxelParams.voxSize/2;
                    currentPt.y = voxelParams.voxSize*(y+1)-voxelParams.voxSize/2;
                    currentPt.z = voxelParams.voxSize*(z+1)-voxelParams.voxSize/2;
                    points.push_back(currentPt);
                }
            }
        }
    }
    writeToPly(points, "tsdfMeshWTF.ply");
}

void createTSDF(float* depthMap, MatrixXf tf, vector<Voxel>& voxelsTSDF){
    Voxel initVoxel;
    initVoxel.value = 1;
    initVoxel.weight = 0;
    std::fill_n(voxelsTSDF.begin(), voxelsTSDF.size(), initVoxel);

    for(int y = 0; y < voxelParams.voxNumy; y++){
        for(int x = 0; x < voxelParams.voxNumx; x++){
            for(int z = 0; z < voxelParams.voxNumz; z++){
                int i = z+x*voxelParams.voxNumz+y*voxelParams.voxNumz*voxelParams.voxNumx;
                //Computing K[x,y,z] = [x,y]
                float curVoxDepth,curVoxDistance;
                float voxelx = voxelParams.voxSize*x-voxelParams.voxSize/2-voxelParams.voxPhysLength/2;
                float voxely = voxelParams.voxSize*y-voxelParams.voxSize/2-voxelParams.voxPhysLength/2;
                float voxelz = voxelParams.voxSize*z-voxelParams.voxSize/2;
                Vector4f poseTransform(voxelx, voxely, voxelz, 1);
                // VectorXf poseTransform(4);
                // poseTransform << voxelx,
                //                  voxely,
                //                  voxelz,
                //                  1;
                VectorXf result = tf*poseTransform;
                voxelx = result(0,0);
                voxely = result(1,0);
                voxelz = result(2,0);
                float onScreenx = -1;
                float onScreeny = -1;
                if(abs(voxelz-0.0)<0.0001) continue;
                onScreenx = (voxelx*fx)/voxelz+cx-0.5;
                onScreeny = (voxely*fy)/voxelz+cy-0.5;
                // //Computing Voxel Depth
                if(onScreenx < pixelWidth && onScreeny < pixelHeight && onScreenx >= 0 && onScreeny >= 0){
                    curVoxDepth = interpolation(onScreenx, onScreeny, depthMap, pixelWidth);
                    // float z00 = depthMap[int(onScreeny)*pixelWidth+int(onScreenx)];
                    // float x00 = (floor(onScreenx)-cx)*z00*fxinv;
                    // float y00 = (floor(onScreeny)-cy)*z00*fyinv;
                    // curVoxDepth = sqrt(pow(x00, 2)+pow(y00, 2)+pow(z00, 2));
                    if(curVoxDepth < 0) continue; //depth isn't normal
                    curVoxDistance = sqrt(pow(voxelx, 2)+pow(voxely, 2)+pow(voxelz, 2));
                    float SDF = curVoxDepth - curVoxDistance;
                    if(abs(SDF)<voxelParams.truncationThrs){
                        voxelsTSDF[i].value = SDF*voxelParams.truncationThrsInv;
                        voxelsTSDF[i].weight += 1;
                    }
                    else if(SDF>0){
                        voxelsTSDF[i].value = 1;
                    }
                }

            }

        }
    }
    cout<<"TSDF done"<<endl;
}