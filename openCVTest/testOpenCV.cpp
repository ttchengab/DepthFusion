//g++ $(pkg-config --cflags --libs opencv4) -std=c++11  testOpenCV.cpp -o test
#include <iostream>
#include <fstream>
#include <opencv2/imgproc.hpp>
#include <opencv2/calib3d.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/rgbd/kinfu.hpp>

using namespace cv;
using namespace cv::kinfu;
using namespace std;

static vector<string> readDepth(std::string fileList);
static vector<string> readDepth(std::string fileList)
{
    vector<string> v;
    fstream file(fileList);
    if(!file.is_open())
        throw std::runtime_error("Failed to read depth list");
    std::string dir;
    size_t slashIdx = fileList.rfind('/');
    slashIdx = slashIdx != std::string::npos ? slashIdx : fileList.rfind('\\');
    dir = fileList.substr(0, slashIdx);
    while(!file.eof())
    {
        std::string s, imgPath;
        std::getline(file, s);
        if(s.empty() || s[0] == '#') continue;
        std::stringstream ss;
        ss << s;
        double thumb;
        ss >> thumb >> imgPath;
        // cout<<dir<<endl;
        // cout<<imgPath<<endl;
        v.push_back(dir+'/'+imgPath);
        cout<<v[v.size()-1]<<endl;
    }

    return v;
}

struct DepthWriter
{
    DepthWriter(string fileList) :
        file(fileList, ios::out), count(0), dir()
    {
        size_t slashIdx = fileList.rfind('/');
        slashIdx = slashIdx != std::string::npos ? slashIdx : fileList.rfind('\\');
        dir = fileList.substr(0, slashIdx);

        if(!file.is_open())
            throw std::runtime_error("Failed to write depth list");

        file << "# depth maps saved from device" << endl;
        file << "# useless_number filename" << endl;
    }

    void append(InputArray _depth)
    {
        Mat depth = _depth.getMat();
        string depthFname = cv::format("%04d.png", count);
        string fullDepthFname = dir + '/' + depthFname;
        if(!imwrite(fullDepthFname, depth))
            throw std::runtime_error("Failed to write depth to file " + fullDepthFname);
        file << count++ << " " << depthFname << endl;
    }

    fstream file;
    int count;
    string dir;
};

namespace Kinect2Params
{
    static const Size frameSize = Size(512, 424);
    // approximate values, no guarantee to be correct
    static const float focal = 366.1f;
    static const float cx = 258.2f;
    static const float cy = 204.f;
    static const float k1 =  0.12f;
    static const float k2 = -0.34f;
    static const float k3 =  0.12f;
};

struct DepthSource
{
public:
    DepthSource(int cam) :
        DepthSource("/Users/tayingcheng/Desktop/FYP/opencvKinFu/rgbd_dataset_freiburg1_xyz/depth.txt", cam)
    { }

    DepthSource(String fileListName) :
        DepthSource(fileListName, -1)
    { }

    DepthSource(String fileListName, int cam) :
        depthFileList(fileListName.empty() ? vector<string>() : readDepth(fileListName)),
        frameIdx(0),
        vc( cam >= 0 ? VideoCapture(VideoCaptureAPIs::CAP_OPENNI2 + cam) : VideoCapture()),
        undistortMap1(),
        undistortMap2(),
        useKinect2Workarounds(true)
    { }

    UMat getDepth()
    {
        UMat out;
        if (!vc.isOpened())
        {
            if (frameIdx < depthFileList.size())
            {
                Mat f = cv::imread(depthFileList[frameIdx++], IMREAD_ANYDEPTH);
                f.copyTo(out);
            }
            else
            {
                return UMat();
            }
        }
        else
        {
            vc.grab();
            vc.retrieve(out, CAP_OPENNI_DEPTH_MAP);

            // workaround for Kinect 2
            if(useKinect2Workarounds)
            {
                out = out(Rect(Point(), Kinect2Params::frameSize));

                UMat outCopy;
                // linear remap adds gradient between valid and invalid pixels
                // which causes garbage, use nearest instead
                remap(out, outCopy, undistortMap1, undistortMap2, cv::INTER_NEAREST);

                cv::flip(outCopy, out, 1);
            }
        }
        if (out.empty())
            throw std::runtime_error("Matrix is empty");
        return out;
    }

    bool empty()
    {
        return depthFileList.empty() && !(vc.isOpened());
    }

    void updateParams(Params& params)
    {
        if (vc.isOpened())
        {
            // this should be set in according to user's depth sensor
            int w = (int)vc.get(VideoCaptureProperties::CAP_PROP_FRAME_WIDTH);
            int h = (int)vc.get(VideoCaptureProperties::CAP_PROP_FRAME_HEIGHT);

            float focal = (float)vc.get(CAP_OPENNI_DEPTH_GENERATOR | CAP_PROP_OPENNI_FOCAL_LENGTH);

            // it's recommended to calibrate sensor to obtain its intrinsics
            float fx, fy, cx, cy;
            Size frameSize;
            if(useKinect2Workarounds)
            {
                fx = fy = Kinect2Params::focal;
                cx = Kinect2Params::cx;
                cy = Kinect2Params::cy;

                frameSize = Kinect2Params::frameSize;
            }
            else
            {
                fx = fy = focal;
                cx = w/2 - 0.5f;
                cy = h/2 - 0.5f;

                frameSize = Size(w, h);
            }

            Matx33f camMatrix = Matx33f(fx,  0, cx,
                                        0,  fy, cy,
                                        0,   0,  1);

            params.frameSize = frameSize;
            params.intr = camMatrix;
            params.depthFactor = 1000.f;

            Matx<float, 1, 5> distCoeffs;
            distCoeffs(0) = Kinect2Params::k1;
            distCoeffs(1) = Kinect2Params::k2;
            distCoeffs(4) = Kinect2Params::k3;
            if(useKinect2Workarounds)
                initUndistortRectifyMap(camMatrix, distCoeffs, cv::noArray(),
                                        camMatrix, frameSize, CV_16SC2,
                                        undistortMap1, undistortMap2);
        }
    }

    vector<string> depthFileList;
    size_t frameIdx;
    VideoCapture vc;
    UMat undistortMap1, undistortMap2;
    bool useKinect2Workarounds;
};


static const char* keys =
{
    "{help h usage ? | | print this message   }"
    "{depth  | | Path to depth.txt file listing a set of depth images }"
    "{camera |0| Index of depth camera to be used as a depth source }"
    "{coarse | | Run on coarse settings (fast but ugly) or on default (slow but looks better),"
        " in coarse mode points and normals are displayed }"
    "{idle   | | Do not run KinFu, just display depth frames }"
    "{record | | Write depth frames to specified file list"
        " (the same format as for the 'depth' key) }"
};

static const std::string message =
 "\nThis demo uses live depth input or RGB-D dataset taken from"
 "\nhttps://vision.in.tum.de/data/datasets/rgbd-dataset"
 "\nto demonstrate KinectFusion implementation \n";


int main(int argc, char **argv)
{
    //cout<<"hi"<<endl;
    bool coarse = false;
    bool idle = false;
    string recordPath;

    CommandLineParser parser(argc, argv, keys);
    parser.about(message);

    if(!parser.check())
    {
        parser.printMessage();
        parser.printErrors();
        return -1;
    }

    if(parser.has("help"))
    {
        parser.printMessage();
        return 0;
    }
    if(parser.has("coarse"))
    {
        coarse = true;
    }
    // if(parser.has("record"))
    // {
    //     recordPath = parser.get<String>("record");
    // }
    recordPath = "/Users/tayingcheng/Desktop/FYP/opencvKinFu/rgbd_dataset_freiburg1_xyz/depthOut.txt";
    if(parser.has("idle"))
    {
        idle = true;
    }

    Ptr<DepthSource> ds;
    if (parser.has("depth"))
        ds = makePtr<DepthSource>(parser.get<String>("depth"));
    else
        ds = makePtr<DepthSource>(parser.get<int>("camera"));

    if (ds->empty())
    {
        std::cerr << "Failed to open depth source" << std::endl;
        parser.printMessage();
        return -1;
    }

    Ptr<DepthWriter> depthWriter;
    if(!recordPath.empty())
        depthWriter = makePtr<DepthWriter>(recordPath);

    Ptr<Params> params;
    Ptr<KinFu> kf;

    if(coarse)
        params = Params::coarseParams();
    else
        params = Params::defaultParams();
    // These params can be different for each depth sensor
    ds->updateParams(*params);
    // Enables OpenCL explicitly (by default can be switched-off)
    cv::setUseOptimized(true);
    if(!idle)
        kf = KinFu::create(params);
    UMat rendered;
    UMat points;
    UMat normals;
    int64 prevTime = getTickCount();
    for(UMat frame = ds->getDepth(); !frame.empty(); frame = ds->getDepth())
    {
        if(depthWriter)
            depthWriter->append(frame);
        {
            UMat cvt8;
            float depthFactor = params->depthFactor;
            convertScaleAbs(frame, cvt8, 0.25*256. / depthFactor);
            if(!idle)
            {
                imshow("depth", cvt8);

                if(!kf->update(frame))
                {
                    kf->reset();
                    std::cout << "reset" << std::endl;
                }
                kf->render(rendered);
            }
            else
            {
                rendered = cvt8;
            }
        }

        int64 newTime = getTickCount();
        putText(rendered, cv::format("FPS: %2d press R to reset, P to pause, Q to quit",
                                     (int)(getTickFrequency()/(newTime - prevTime))),
                Point(0, rendered.rows-1), FONT_HERSHEY_SIMPLEX, 0.5, Scalar(0, 255, 255));
        prevTime = newTime;

        imshow("render", rendered);

        int c = waitKey(1);
        switch (c)
        {
        case 'r':
            if(!idle)
                kf->reset();
            break;
        case 'q':
            return 0;
        default:
            break;
        }
    }
    return 0;
}
