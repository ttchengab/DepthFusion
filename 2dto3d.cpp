#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <vector> // for 2D vector
#include <bitset>
#include <cmath>
#include "eigenlib/Eigen/Dense"
#include "eigenlib/Eigen/Core"
#include "eigenlib/Eigen/Geometry"
using namespace Eigen;
using namespace std;
const int pixelWidth = 640, pixelHeight = 480;
const float fx = 525.0, fy = 525.0, cx = 319.5, cy = 239.5, depthFactor = 5000.0;
struct Point{
    float x, y, z;
};
class VoxelParams{
    public:
        float voxPhysLength, voxPhysWidth, voxPhysDepth;
        float voxSize;
        int voxNumx, voxNumy, voxNumz;
        int voxVolume;
        float truncationThrs;
        VoxelParams(float physLength, int numx, int numy, int numz)
            : voxPhysLength(physLength), voxNumx(numx), voxNumy(numy), voxNumz(numz){
            voxPhysWidth = voxPhysLength/voxNumx*voxNumy;
            voxPhysDepth = voxPhysLength/voxNumx*voxNumz;
            voxVolume = voxNumx*voxNumy*voxNumz;
            voxSize = voxPhysLength/voxNumx;
            truncationThrs = voxSize*5;
        }
};
//Voxel Parameters Initialization
VoxelParams voxelParams(2.0, 256, 256, 512);

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
    float x00 = (floor(x)-cx)*z00/fx;
    float x01 = (floor(x)-cx)*z01/fx;
    float x10 = (floor(x+1)-cx)*z10/fx;
    float x11 = (floor(x+1)-cx)*z11/fx;
    float y00 = (floor(y)-cy)*z00/fy;
    float y01 = (floor(y)-cy)*z01/fy;
    float y10 = (floor(y+1)-cy)*z10/fy;
    float y11 = (floor(y+1)-cy)*z11/fy;
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
            currentPt.x = (x-cx)*currentPt.z/fx;
            currentPt.y = (y-cy)*currentPt.z/fy;
            points.push_back(currentPt);
        }
    }
    return points;
}
void surfaceMeasurement(float *depthMap){
    vector<Point> points = computePoints(depthMap);
    writeToPly(points, "testMesh.ply");
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
void rayCasting(const vector<float> &voxelsTSDF){
    float minimumDepth = 0;
    vector<Point>points;
    float translationx = 0.0, translationy = 0.0, translationz = 0.0;
    MatrixXf rotation = transformation(Vector3f(0.0,0.0,0.0), M_PI/12.0, 0.0, 0.0);
    for(int y = 0; y < pixelHeight; y++){
        for(int x = 0; x < pixelWidth; x++){
            int i = y*pixelWidth+x;
            Point pixPt;
            pixPt.z = 1.0;
            pixPt.x = (x-cx)/fx;
            pixPt.y = (y-cy)/fy;
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
            while(curTSDF > 0 && locxv < voxelParams.voxNumx && locxv >=0 && locyv < voxelParams.voxNumy && locyv >=0 && loczv < voxelParams.voxNumz && loczv >=0){
                //cout<<locxv<<" "<<locyv<<" "<<loczv<<endl;
                curTSDF = interpolation(locxv, locyv, loczv, voxelsTSDF, voxelParams.voxNumz, voxelParams.voxNumx);
                locxv = (locx+voxelParams.voxPhysLength/2)/voxelParams.voxSize;
                locyv = (locy+voxelParams.voxPhysWidth/2)/voxelParams.voxSize;
                loczv = (locz)/voxelParams.voxSize;
                locx += stepx;
                locy += stepy;
                locz += stepz;
            }
            if( locxv < voxelParams.voxNumx && locxv >=0 && locyv < voxelParams.voxNumy && locyv >=0 && loczv < voxelParams.voxNumz && loczv >=0){
                Point currentPt;
                currentPt.x = locx;
                currentPt.y = locy;
                currentPt.z = locz;
                points.push_back(currentPt);
            }
        }
    }
    cout<<"Raycast done";
    writeToPly(points, "rayCastMesh.txt");
    writeToPly(points, "rayCastMesh.ply");
}
void validateTSDF(const vector<float> &voxelsTSDF){
    vector<Point> points;
    Point currentPt;
    for(int y = 0; y <  voxelParams.voxNumx; y++){
        for(int x = 0; x <  voxelParams.voxNumx; x++){
            for(int z = 0; z < voxelParams.voxNumz; z++){
                int i = z+x*voxelParams.voxNumz+y*voxelParams.voxNumz* voxelParams.voxNumx;
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

void createTSDF(float* depthMap){
    vector<float> voxelsTSDF(voxelParams.voxVolume, 1.0);
    MatrixXf tf = transformation(Vector3f(0.0,0.0,0.0), M_PI/6.0, 0.0, 0.0);
    for(int y = 0; y <  voxelParams.voxNumy; y++){
        for(int x = 0; x <  voxelParams.voxNumx; x++){
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
                float onScreenx = (voxelx*fx)/voxelz+cx-0.5;
                float onScreeny = (voxely*fy)/voxelz+cy-0.5;
                //Computing Voxel Depth
                if(onScreenx < pixelWidth && onScreeny < pixelHeight && onScreenx >= 0 && onScreeny >= 0){
                    curVoxDepth = interpolation(onScreenx, onScreeny, depthMap, pixelWidth);
                    if(curVoxDepth<0) continue; //depth isn't normal
                    curVoxDistance = sqrt(pow(voxelx, 2)+pow(voxely, 2)+pow(voxelz, 2));
                    float SDF = curVoxDepth - curVoxDistance;
                    if(abs(SDF)<voxelParams.truncationThrs){
                        voxelsTSDF[i] = SDF/voxelParams.truncationThrs;
                    }
                    else if(SDF>0){
                        voxelsTSDF[i] = 1;
                    }
                }
            }
        }
    }
    cout<<"TSDF done"<<endl;
    //validateTSDF(voxelsTSDF); //Used to validate TSDF
    rayCasting(voxelsTSDF);
}

int main(){
    float *depthMap;
    depthMap = getDepthMap("1305031102.160407.png.bin");
    cout<<"depthMap get"<<endl;
    surfaceMeasurement(depthMap);
    cout<<"surface measurement done"<<endl;
    createTSDF(depthMap);
    return 0;
}
