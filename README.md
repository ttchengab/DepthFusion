# KinectFusion
Implementation of KinectFusion to work on TUM datasets.
To run the code you have to download the Eigen library and put it in the same directory with the cpp files.

Compile with the command:

```g++ 2dto3d.cpp HelperFunctions.cpp PreProcessing.cpp SurfaceMeasurement.cpp TSDFConstruction.cpp -o test```

Change testmode between 1 and 2 in 2dto3d.cpp to get the fusion from known data or from icp. A ply file will be generated from running the program, and it can be visualized via tools such as Meshlab.

An explanation of the code can be found here:
https://medium.com/@taying.cheng/understanding-real-time-3d-reconstruction-and-kinectfusion-33d61d1cd402#3fe5-fcb284305ef2

