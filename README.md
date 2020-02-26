# KinectFusion
Implementation of KinectFusion to work on TUM datasets.
To run the code you have to download the eigen library and put it in the same directory of the cpp file.

Compile with the command

```g++ 2dto3d.cpp HelperFunctions.cpp PreProcessing.cpp SurfaceMeasurement.cpp TSDFConstruction.cpp -o icptest```

Change testmode between 1 and 2 in 2dto3d.cpp to get the fusion from known data or from icp

