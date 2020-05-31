# KinectFusion
This is an implementation of a lightweight version of KinectFusion to work on TUM datasets. This program takes in a binary stream instead of the original png depth images. 5 images from the dataset were already given, and if you want to test with more depth maps, please use the ```readImg.ipynb``` converter to convert the files into binary streams first.

To run the code you have to download the Eigen library and put it in the same directory with the cpp files.

Compile with the command:

```g++ 2dto3d.cpp HelperFunctions.cpp PreProcessing.cpp SurfaceMeasurement.cpp TSDFConstruction.cpp -o test```

Change testmode between 1 and 2 in the main function in ```2dto3d.cpp``` to get the fusion from either from groundtruth pose or from icp. A ply file will be generated from running the program, and it can be visualized via tools such as Meshlab.

An explanation of the code can be found here:

https://medium.com/@taying.cheng/understanding-real-time-3d-reconstruction-and-kinectfusion-33d61d1cd402#3fe5-fcb284305ef2

