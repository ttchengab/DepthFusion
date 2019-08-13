#include <iostream>
#include <cmath>

const int pixelWidth = 640, pixelHeight = 480;
const float fx = 525.0, fy = 525.0, cx = 319.5, cy = 239.5, depthFactor = 5000.0;
const float fxinv = 1.0/double(fx), fyinv = 1.0/double(fy);
const float epislond = 0.1, epislonAng = sqrt(3)/2;
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
struct icpPose{
    float betta, gamma,alpha;
    float tx, ty, tz;
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