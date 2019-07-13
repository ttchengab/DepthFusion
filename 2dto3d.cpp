#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <vector> // for 2D vector
#include <bitset>
#include <cmath>
using namespace std;
const int pixelWidth = 640, pixelHeight = 480;
const float fx = 525.0, fy = 525.0, cx = 319.5, cy = 239.5, depthFactor = 5000.0;
const float voxPhysLength = 2.0;
const int voxelVolume = 256*256*512, voxelNumxy = 256, voxelNumz = 512;
const float voxelSize = 2.0/voxelNumxy;
const float truncationThrs = voxelSize*5;
struct Point{
    float x, y, z;
};
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
    for(int i = 0; i < pixelWidth*pixelHeight; i ++){
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
vector<Point> computePoints(float *depthMap){
    vector<Point> points;
    Point currentPt;
    for(int i = 0; i < pixelWidth*pixelHeight; i ++){
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
void writeToPly(vector<Point> points, const char* fileName){
    string s = "ply\nformat ascii 1.0\nelement vertex "+to_string(points.size())+"\nproperty float x\nproperty float y\nproperty float z\nend_header\n";
    FILE * bfile;
    bfile = fopen (fileName, "wb");
    unsigned int N(s.size());
    fwrite(s.c_str(),1, N ,bfile);
    for(int i = 0; i < points.size(); i ++){
        string pointS = to_string(points[i].x)+" "+to_string(points[i].y)+" "+to_string(points[i].z);
        if(i != points.size()-1) pointS = pointS+"\n";
        unsigned int Ns(pointS.size());
        fwrite(pointS.c_str(),1,Ns,bfile);
    }
    fflush(bfile);
    fclose (bfile);
}
void surfaceMeasurement(float *depthMap){
    //variable initialization
    vector<Point> points = computePoints(depthMap);
    writeToPly(points, "testMesh.ply");
}


void validateTSDF(vector<float> voxelsTSDF, float *depthMap){
    vector<Point> points;
    Point currentPt;
    for(int y = 0; y < voxelNumxy; y++){
        for(int x = 0; x < voxelNumxy; x++){
            for(int z = 0; z < voxelNumz; z++){
                int i = z+x*voxelNumz+y*voxelNumz*voxelNumxy;
                if(fabs(voxelsTSDF[i])<0.2){
                    currentPt.x = voxelSize*(x+1)-voxelSize/2.0;
                    currentPt.y = voxelSize*(y+1)-voxelSize/2.0;
                    currentPt.z = voxelSize*(z+1)-voxelSize/2.0;
                    points.push_back(currentPt);
                }
            }
        }
    }
    /* 
    float minimumDepth = 50.0;
    for(int i = 0; i < pixelWidth*pixelHeight; i ++){
        if(minimumDepth > depthMap[i]) minimumDepth = depthMap[i];
    }
    for(int y = 0; y < pixelHeight; y ++){
        for(int x = 0; x < pixelWidth; x ++){
            int i = y*pixelWidth+x;
        }

    }
    */
    // for(y){
    //     for(x){
    //         float rayx = ;
    //         float rayy = ;
    //         float rayz = ;

    //         float locx = ;
    //         float locy = ;
    //         float locz = ;

    //         float step = voxelSize;
    //         while(tsdf > 0){
    //             locx += step * rayx;

    //             if (distance > ) {
    //                 break;
    //             }
    //         }
    //     }
    // }

    writeToPly(points, "tsdfMesh.ply");
}
void createTSDF(float* depthMap){
    vector<float> voxelsTSDF(voxelVolume, 1.0);
    for(int y = 0; y < voxelNumxy; y++){
        for(int x = 0; x < voxelNumxy; x++){
            for(int z = 0; z < voxelNumz; z++){
                int i = z+x*voxelNumz+y*voxelNumz*voxelNumxy;
                //Computing K[x,y,z] = [x,y]
                float curVoxDepth,curVoxDistance;
                float voxelx = voxelSize*(x+1)-voxelSize/2.0-voxPhysLength/2;
                float voxely = voxelSize*(y+1)-voxelSize/2.0-voxPhysLength/2;
                float voxelz = voxelSize*(z+1)-voxelSize/2.0;
                float onScreenx = (voxelx*fx)/voxelz+cx-0.5;
                float onScreeny = (voxely*fy)/voxelz+cy-0.5;

                //Computing Voxel Depth
                if(onScreenx < 640 && onScreeny < 480 && onScreenx >= 0 && onScreeny >= 0){
                    float weightx = onScreenx - floor(onScreenx);
                    float weighty = onScreeny - floor(onScreeny);
                    float d00 = depthMap[int(onScreeny)*pixelWidth+int(onScreenx)];
                    float d01 = depthMap[int(onScreeny+1)*pixelWidth+int(onScreenx)];
                    float d10 = depthMap[int(onScreeny)*pixelWidth+int(onScreenx)+1];
                    float d11 = depthMap[int(onScreeny+1)*pixelWidth+int(onScreenx)+1];
                    curVoxDepth = (d00*(1-weightx)+d10*weightx)*(1-weighty)+(d01*(1-weightx)+d11*weightx)*(weighty);
                    curVoxDistance = sqrt(pow(voxelx, 2)+pow(voxely, 2)+pow(voxelz, 2));
                    float SDF = curVoxDepth - curVoxDistance;
                    if(abs(SDF)<truncationThrs){
                        voxelsTSDF[i] = SDF/truncationThrs;
                    }
                    else if(SDF>0){
                        voxelsTSDF[i] = 1;
                    }
                }
            }
        }
    }
    validateTSDF(voxelsTSDF, depthMap);
}
int main(){
    float *depthMap;
    depthMap = getDepthMap("1305031102.160407.png.bin");
    surfaceMeasurement(depthMap);
    createTSDF(depthMap);
    /*
    //checking depthMap function
    for(int i = 0; i < pixelWidth*pixelHeight; i ++){
        if(depthMap[i]!=0){
            //cout<<i;
            break;
        }
    }
    */
    return 0;
}
