#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <vector> // for 2D vector
#include <bitset>
#include <cmath>
#include <cstring>
#include <fstream>
#include "ParamsInit.cpp"
using namespace std;

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
