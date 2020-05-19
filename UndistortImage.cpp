#include <stdio.h>
#include <opencv2/opencv.hpp>
#include <string>

using std::string;
 
using namespace cv;


/// Function to undistort image
/// Input image: orig
/// Output image as return value
Mat undistort_image( Mat& orig ){
  Mat undistorted = orig.clone();

  //Mat cameraMatrix = Mat::eye( 3, 3, CV_64F ); //K
  //Mat distCoeffs = Mat::zeros( 1, 4, CV_64F ); //D
  //  Here's the constants extracted directly from MATLAB:
  //           FocalLength: [2.7605e+03 2.7670e+03]
  //        PrincipalPoint: [1.9143e+03 1.5964e+03]
  //             ImageSize: [3000 4000]
  //      RadialDistortion: [-0.2398 0.1145]
  //  TangentialDistortion: [0 0]
  //               Skew: 0
  //  
  // The intrinsics fx, fy, cx and cy in your intrinsics matrix are the
  // focal length and principle point above, and are dependant on the
  // image size so best to make sure the images you're using haven't
  // been rescaled.
  // The four values only give the minimum k1, k2, p1, p2 that you need 
  // because calibrating additional parameters weren't enabled by the 
  // students when they did this calibration. Tangential distortion is set to 
  // zero because it was also disabled.  


  Matx33f cameraMatrix( 2.7605e+03, 0,          1.9143e+03, 
                        0,          2.7670e+03, 1.5964e+03, 
                        0,          0,          1 );
  Vec4f distCoeffs( -0.2398, 0.1145, 0, 0 ); 

  undistort(orig, undistorted, cameraMatrix, distCoeffs);

  return undistorted;
}


string get_output_filename(int argc, char**argv ){
    string outputname;
    if ( argc ==2 ) {
      string infilename( argv[1] );
      size_t idx = infilename.find_last_of("/");
      outputname = string("undistorted_") + infilename.substr(idx+1 );
    } else {
      outputname = string( argv[2] );
    }

    return outputname;
}


int main(int argc, char** argv )
{
    if ( argc != 2 && argc != 3 )
    {
        printf("usage: UndistortImage <Input_image_with_path> [Ouptut_image_with_path]\n");
        return -1;
    }

    string outputname = get_output_filename( argc, argv );


    Mat image = imread( argv[1], IMREAD_COLOR );
    if ( !image.data )
    {
        printf("No image data \n");
        return -1;
    }

    // undistort image
    Mat  img_undistorted = undistort_image( image );

    imwrite( outputname, img_undistorted );

    return 0;
}
