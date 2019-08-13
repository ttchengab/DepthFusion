#include <iostream>
#include "eigenlib/Eigen/Dense"
#include "eigenlib/Eigen/Core"
#include "eigenlib/Eigen/Geometry"
using namespace Eigen;
using namespace std;


VectorXf transformation(VectorXf cameraPose,Vector3f translation, float rotx, float roty, float rotz){
  Transform<float, 3, Affine> t = Transform<float, 3, Affine>::Identity();
  t.translate(translation);
  t.rotate(AngleAxisf(rotx, Vector3f::UnitX()));
  t.rotate(AngleAxisf(roty, Vector3f::UnitY()));
  t.rotate(AngleAxisf(rotz, Vector3f::UnitZ()));
  VectorXf result;
  result = t.inverse()*cameraPose.homogeneous();
  return result;
}



int main(){
   /*
    VectorXf cameraPose(3);
    cameraPose << 1.0,
                  0.0,
                  0.0;
    VectorXf result = transformation(cameraPose, Vector3f(2,2,2), 0.0, 0.0, 0.0);


    cout<<result;
  */
    float angle = M_PI/4;
    float sinA = std::sin(angle/2);
    float cosA = std::cos(angle/2);

    Quaternionf q;
    q.x() = 0;
    q.y() = 0;
    q.z() = 0;
    q.w() = 0;
    MatrixXf rotMatrix = q.normalized().toRotationMatrix();

    Transform<float, 3, Affine> t = Transform<float, 3, Affine>::Identity();
    t.translate(Vector3f(2,2,2));
    //MatrixXd TransformMatrix = rotMatrix*t.matrix();
    MatrixXf lol = t.matrix();
    lol.block<3,3>(0,0) = rotMatrix;
    cout<<lol<<endl;
}
/*
1305031102.160407 depth/1305031102.160407.png
1305031102.194330 depth/1305031102.194330.png
1305031102.226738 depth/1305031102.226738.png
1305031102.262886 depth/1305031102.262886.png
1305031102.295279 depth/1305031102.295279.png
1305031102.329195 depth/1305031102.329195.png
1305031102.363013 depth/1305031102.363013.png
1305031102.394772 depth/1305031102.394772.png
1305031102.427815 depth/1305031102.427815.png
1305031102.462395 depth/1305031102.462395.png
1305031102.494271 depth/1305031102.494271.png
1305031102.526330 depth/1305031102.526330.png
1305031102.562224 depth/1305031102.562224.png
1305031102.594158 depth/1305031102.594158.png
1305031102.626818 depth/1305031102.626818.png
1305031102.663273 depth/1305031102.663273.png
1305031102.695165 depth/1305031102.695165.png
1305031102.728423 depth/1305031102.728423.png
1305031102.763549 depth/1305031102.763549.png
1305031102.794978 depth/1305031102.794978.png
1305031102.828537 depth/1305031102.828537.png
1305031102.862808 depth/1305031102.862808.png
1305031102.894167 depth/1305031102.894167.png
1305031102.926851 depth/1305031102.926851.png
1305031102.962137 depth/1305031102.962137.png
1305031102.994164 depth/1305031102.994164.png
1305031103.027881 depth/1305031103.027881.png
1305031103.062273 depth/1305031103.062273.png
1305031103.094040 depth/1305031103.094040.png
1305031103.129109 depth/1305031103.129109.png
1305031103.162795 depth/1305031103.162795.png
*/
