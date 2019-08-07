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
using namespace Eigen;
using namespace std;
using namespace std::chrono;
const int pixelWidth = 640, pixelHeight = 480;
const float fx = 525.0, fy = 525.0, cx = 319.5, cy = 239.5, depthFactor = 5000.0;
const float fxinv = 1.0/double(fx), fyinv = 1.0/double(fy);
int totalWeight = 0;
struct Point{
    float x, y, z;
};
struct QuatPose{
    float tx, ty, tz;
    float qx, qy, qz, qw;
    QuatPose(float tranx, float trany, float tranz, float qax, float qay, float qaz, float qaw)
        : tx(tranx), ty(trany), tz(tranz), qx(qax), qy(qay), qz(qaz), qw(qaw)
        {}
};
class VoxelParams{
    public:
        float voxPhysLength, voxPhysWidth, voxPhysDepth;
        float voxSize;
        int voxNumx, voxNumy, voxNumz;
        int voxVolume;
        float truncationThrs;
        float truncationThrsInv;
        VoxelParams(float physLength, int numx, int numy, int numz)
            : voxPhysLength(physLength), voxNumx(numx), voxNumy(numy), voxNumz(numz){
            voxPhysWidth = voxPhysLength/voxNumx*voxNumy;
            voxPhysDepth = voxPhysLength/voxNumx*voxNumz;
            voxVolume = voxNumx*voxNumy*voxNumz;
            voxSize = voxPhysLength/voxNumx;
            truncationThrs = voxSize*5;
            truncationThrsInv = 1/truncationThrsInv;
        }
};
//Voxel Parameters Initialization
VoxelParams voxelParams(2.0, 256, 256, 256);

QuatPose clibRot(1.3452, 0.6273, 1.6627, 0.6582, 0.6109, -0.2950, -0.3265);

//miscellaneous
vector<string> getGroundTruthPose(const char* fileName){
    ifstream myReadFile;
    myReadFile.open(fileName);
    char output[100];
    vector<string> poses;
    if (myReadFile.is_open()) {
        int count = 0;
        while (!myReadFile.eof()) {
            myReadFile >> output;
            count++;
            if(count > 16) poses.push_back(output);
        }
    //16
    }
    myReadFile.close();
    return poses;
}
QuatPose getNearestPose(float imageTime, const vector<string>& poses){
    float closestTime = stof(poses[0].substr(7,15));
    float closestDifference = abs(closestTime-imageTime);
    float index = 0;
    for(int i = 0; i+8 < poses.size(); i = i+8){
        float curTime = stof(poses[i].substr(7,15));
        if(abs(curTime-imageTime)<closestDifference){
            closestTime = curTime;
            closestDifference = abs(curTime-imageTime);
            index = i;
        }
    }
    QuatPose nearestPose(stof(poses[index+1]), stof(poses[index+2]), stof(poses[index+3]), stof(poses[index+4]),
                         stof(poses[index+5]), stof(poses[index+6]), stof(poses[index+7]));
    return nearestPose;
}


//Depth Map Retrieval
float *getDepthMap(const char* fileName){
    FILE * pFile;
    long lSize;
    char * buffer;
    size_t result;
    pFile = fopen ( fileName , "rb" );
    errno=0;
    if (pFile==NULL) printf("Error %d \n", errno);
    // obtain file size:
    fseek (pFile , 0 , SEEK_END);
    lSize = ftell (pFile);
    rewind (pFile);
    // allocate memory to contain the whole file:
    buffer = (char*) malloc (sizeof(char)*lSize);
    if (buffer == NULL) fputs ("Memory error",stderr);
    // copy the file into the buffer:
    result = fread (buffer,1,lSize,pFile);
    if (result != lSize) fputs ("Reading error",stderr);
    static float depthMap[pixelWidth*pixelHeight];
    for(int i = 0; i < pixelWidth*pixelHeight; i++){
        if(buffer[2*i]!=0){
            bitset<8> small8(buffer[2*i]);
            bitset<8> big8(buffer[2*i+1]);
            bitset<16> depth(big8.to_string()+small8.to_string());
            depthMap[i]=depth.to_ulong()/depthFactor;
        }
    }
    fclose (pFile);
    free (buffer);
    cout<<"depthMap get"<<endl;
    return depthMap;
}
//Write to Ply file
void writeToPly(vector<Point> points, const char* fileName){
    string s = "ply\nformat ascii 1.0\nelement vertex "+to_string(points.size())+"\nproperty float x\nproperty float y\nproperty float z\nend_header\n";
    FILE * bfile;
    bfile = fopen (fileName, "wb");
    unsigned int N(s.size());
    fwrite(s.c_str(),1, N ,bfile);
    for(int i = 0; i < points.size(); i++){
        string pointS = to_string(points[i].x)+" "+to_string(points[i].y)+" "+to_string(points[i].z);
        if(i != points.size()-1) pointS = pointS+"\n";
        unsigned int Ns(pointS.size());
        fwrite(pointS.c_str(),1,Ns,bfile);
    }
    fflush(bfile);
    fclose (bfile);
}

//Interpolation 2 dimensional and 3 dimensional
float interpolation(float x, float y, const float *map, int width){
    float weightx = x-floor(x);
    float weighty = y-floor(y);
    float z00 = map[int(y)*width+int(x)];
    float z01 = map[(int(y+1))*width+int(x)];
    float z10 = map[int(y)*width+int(x)+1];
    float z11 = map[(int(y+1))*width+int(x)+1];
    float x00 = (floor(x)-cx)*z00*fxinv;
    float x01 = (floor(x)-cx)*z01*fxinv;
    float x10 = (floor(x+1)-cx)*z10*fxinv;
    float x11 = (floor(x+1)-cx)*z11*fxinv;
    float y00 = (floor(y)-cy)*z00*fyinv;
    float y01 = (floor(y)-cy)*z01*fyinv;
    float y10 = (floor(y+1)-cy)*z10*fyinv;
    float y11 = (floor(y+1)-cy)*z11*fyinv;
    float d00 = sqrt(pow(x00, 2)+pow(y00, 2)+pow(z00, 2));
    float d01 = sqrt(pow(x01, 2)+pow(y01, 2)+pow(z01, 2));
    float d10 = sqrt(pow(x10, 2)+pow(y10, 2)+pow(z10, 2));
    float d11 = sqrt(pow(x11, 2)+pow(y11, 2)+pow(z11, 2));
    if(!isnormal(d00) || !isnormal(d01) || !isnormal(d10) || !isnormal(d11)){
        return -1;
    }
    float curVoxDepth = (d00*(1-weightx)+d10*weightx)*(1-weighty)+(d01*(1-weightx)+d11*weightx)*(weighty);
    return curVoxDepth;
}
float interpolation(float x, float y, float z,  const vector<float> &map3d, int depth, int width){
    float interpolatedV;
    //cout<<x<<" "<<y<<" "<<z<<endl;
    if(floor(y+1)<256 && floor(x+1)<256 && floor(z+1)<512){
        float weightx = x-floor(x);
        float weighty = y-floor(y);
        float weightz = z-floor(z);
        float f000 = map3d[int(y)*width*depth+int(x)*depth+int(z)];
        float f001 = map3d[int(y)*width*depth+int(x)*depth+int(z+1)];
        float f010 = map3d[int(y+1)*width*depth+int(x)*depth+int(z)];
        float f011 = map3d[int(y+1)*width*depth+int(x)*depth+int(z+1)];
        float f100 = map3d[int(y)*width*depth+int(x+1)*depth+int(z)];
        float f101 = map3d[int(y)*width*depth+int(x+1)*depth+int(z+1)];
        float f110 = map3d[int(y+1)*width*depth+int(x+1)*depth+int(z)];
        float f111 = map3d[int(y+1)*width*depth+int(x+1)*depth+int(z+1)];
        float xy1 = (f000*(1-weightx)+f100*weightx)*(1-weighty)+(f010*(1-weightx)+f110*weightx)*weighty;
        float xy2 = (f001*(1-weightx)+f101*weightx)*(1-weighty)+(f011*(1-weightx)+f111*weightx)*weighty;
        interpolatedV = xy1*(1-weightz)+xy2*weightz;
    }
    else if(y != 256 && x != 256 && z!=512){
         interpolatedV = map3d[int(y)*width*depth+int(x)*depth+int(z)];
     }
    else{
        interpolatedV = map3d[int(y-1)*width*depth+int(x-1)*depth+int(z-1)];
    }
    return interpolatedV;
}

//Part 1: Surface Measurement Functions
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
}
void surfaceMeasurement(float *depthMap){
    vector<Point> vertexMap = computePoints(depthMap);
    writeToPly(vertexMap, "testMesh.ply");
    vector<Point> normalMap = computeNormal(vertexMap);
    //test normal map
    cout<<"surface measurement finished"<<endl;
}

//Part 2 and 3: Create TSDF and RayCast based on TSDF
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
void validateTSDF(const vector<float> &voxelsTSDF){
    vector<Point> points;
    Point currentPt;
    for(int y = 0; y < voxelParams.voxNumx; y++){
        for(int x = 0; x < voxelParams.voxNumx; x++){
            for(int z = 0; z < voxelParams.voxNumz; z++){
                int i = z+x*voxelParams.voxNumz+y*voxelParams.voxNumz*voxelParams.voxNumx;
                if(fabs(voxelsTSDF[i])<0.2){
                    currentPt.x = voxelParams.voxSize*(x+1)-voxelParams.voxSize/2;
                    currentPt.y = voxelParams.voxSize*(y+1)-voxelParams.voxSize/2;
                    currentPt.z = voxelParams.voxSize*(z+1)-voxelParams.voxSize/2;
                    points.push_back(currentPt);
                }
            }
        }
    }
    writeToPly(points, "tsdfMesh.ply");
}

vector<float> createTSDF(float* depthMap, MatrixXf tf){
    vector<float> voxelsTSDF(voxelParams.voxVolume, 1.0);

    for(int y = 0; y < voxelParams.voxNumy; y++){
        for(int x = 0; x < voxelParams.voxNumx; x++){
            for(int z = 0; z < voxelParams.voxNumz; z++){
                int i = z+x*voxelParams.voxNumz+y*voxelParams.voxNumz*voxelParams.voxNumx;
                //Computing K[x,y,z] = [x,y]
                float curVoxDepth,curVoxDistance;
                float voxelx = voxelParams.voxSize*x-voxelParams.voxSize/2-voxelParams.voxPhysLength/2;
                float voxely = voxelParams.voxSize*y-voxelParams.voxSize/2-voxelParams.voxPhysLength/2;
                float voxelz = voxelParams.voxSize*z-voxelParams.voxSize/2;
                VectorXf poseTransform(4);
                poseTransform << voxelx,
                                 voxely,
                                 voxelz,
                                 1;
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
                        voxelsTSDF[i] = SDF*voxelParams.truncationThrsInv;
                    }
                    else if(SDF>0){
                        voxelsTSDF[i] = 1;
                    }
                }

            }

        }
    }
    cout<<"TSDF done"<<endl;
    return voxelsTSDF;

}
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

vector<string> getImageNames(const char* fileName){
    ifstream myReadFile;
    myReadFile.open(fileName);
    char output[100];
    vector<string> imageNames;
    if (myReadFile.is_open()) {
        int count = 0;
        while (!myReadFile.eof()) {
            myReadFile >> output;
            count++;
            if(count > 9 && count%2 == 0) imageNames.push_back(output);
        }
    }
    myReadFile.close();
    return imageNames;
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
