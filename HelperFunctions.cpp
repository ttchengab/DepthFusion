#include "HelperFunctions.h"


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

void writeToPlyNorm(vector<Point> points, vector<Point> normals, const char* fileName){
    string s = "ply\nformat ascii 1.0\nelement vertex "+to_string(points.size())+"\nproperty float x\nproperty float y\nproperty float z\n";
    s = s+"property uchar red\nproperty uchar green\nproperty uchar blue\nend_header\n";
    FILE * bfile;
    bfile = fopen (fileName, "wb");
    unsigned int N(s.size());
    fwrite(s.c_str(),1, N ,bfile);
    float maxNormx=normals[0].x, maxNormy=normals[0].y, maxNormz=normals[0].z;
    float minNormx=normals[0].x, minNormy=normals[0].y, minNormz=normals[0].z;
    for(int i = 0; i < points.size(); i++){
      if(maxNormx < normals[i].x) maxNormx=normals[i].x;
      if(maxNormy < normals[i].y) maxNormy=normals[i].y;
      if(maxNormz < normals[i].z) maxNormz=normals[i].z;

      if(minNormx > normals[i].x) minNormx=normals[i].x;
      if(minNormy > normals[i].y) minNormy=normals[i].y;
      if(minNormz > normals[i].z) minNormz=normals[i].z;
    }
    // cout<<maxNormx<<" "<<maxNormy<<" "<<maxNormz<<endl;
    // cout<<minNormx<<" "<<minNormy<<" "<<minNormz<<endl;
    for(int i = 0; i < points.size(); i++){
      normals[i].x = (maxNormx-normals[i].x)/(maxNormx-minNormx)*255;
      normals[i].y = (maxNormy-normals[i].y)/(maxNormy-minNormy)*255;
      normals[i].z = (maxNormz-normals[i].z)/(maxNormz-minNormz)*255;
    }
    for(int i = 0; i < points.size(); i++){
        string pointS = to_string(points[i].x)+" "+to_string(points[i].y)+" "+to_string(points[i].z)+" ";
        pointS = pointS+to_string(int(normals[i].x))+" "+to_string(int(normals[i].y))+" "+to_string(int(normals[i].z));
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
float interpolation(float x, float y, float z,  const vector<Voxel> &map3d, int depth, int width){
    float interpolatedV;
    //cout<<x<<" "<<y<<" "<<z<<endl;
    if(floor(y+1)<256 && floor(x+1)<256 && floor(z+1)<512){
        float weightx = x-floor(x);
        float weighty = y-floor(y);
        float weightz = z-floor(z);
        vector<float> f;
        float f000 = map3d[int(y)*width*depth+int(x)*depth+int(z)].value;
        float f001 = map3d[int(y)*width*depth+int(x)*depth+int(z+1)].value;
        float f010 = map3d[int(y+1)*width*depth+int(x)*depth+int(z)].value;
        float f011 = map3d[int(y+1)*width*depth+int(x)*depth+int(z+1)].value;
        float f100 = map3d[int(y)*width*depth+int(x+1)*depth+int(z)].value;
        float f101 = map3d[int(y)*width*depth+int(x+1)*depth+int(z+1)].value;
        float f110 = map3d[int(y+1)*width*depth+int(x+1)*depth+int(z)].value;
        float f111 = map3d[int(y+1)*width*depth+int(x+1)*depth+int(z+1)].value;
        f.push_back(f000);
        f.push_back(f001);
        f.push_back(f010);
        f.push_back(f011);
        f.push_back(f100);
        f.push_back(f101);
        f.push_back(f110);
        f.push_back(f111);
        float maxf = f[0], minf = f[0];
        for(int i = 1; i < f.size(); i++){
            if(f[i]>maxf) maxf = f[i];
            if(f[i]<minf) minf = f[i];
        }
        if(maxf-minf > 1){
            return -2;
        }
        float xy1 = (f000*(1-weightx)+f100*weightx)*(1-weighty)+(f010*(1-weightx)+f110*weightx)*weighty;
        float xy2 = (f001*(1-weightx)+f101*weightx)*(1-weighty)+(f011*(1-weightx)+f111*weightx)*weighty;
        interpolatedV = xy1*(1-weightz)+xy2*weightz;
    }
    else if(y != voxelParams.voxNumy && x != voxelParams.voxNumx && z!= voxelParams.voxNumz){
         interpolatedV = map3d[int(y)*width*depth+int(x)*depth+int(z)].value;
         return -2;
     }
    else{
        interpolatedV = map3d[int(y-1)*width*depth+int(x-1)*depth+int(z-1)].value;
        return -2;
    }
    return interpolatedV;
}
