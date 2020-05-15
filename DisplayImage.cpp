#include <stdio.h>
#include <opencv2/opencv.hpp>

 
using namespace cv;

Mat& ScanImageAndReduceC(Mat& I, const uchar* const table)
{
  std::cout<<"Image depth is "<<I.depth()<<std::endl;
  int channels = I.channels();
  int nRows = I.rows;
  int nCols = I.cols * channels;

  std::cout<<"Image has "<<channels<<" channels "<<nRows<<" rows x "<<nCols<<" columns"<<std::endl;

  // accept only char type matrices
  //CV_Assert(I.depth() == CV_8U);


  if (I.isContinuous())
    {
      nCols *= nRows;
      nRows = 1;
    }
  int i,j;
  uchar* p;
  for( i = 0; i < nRows; ++i)
    {
      p = I.ptr<uchar>(i);
      for ( j = 0; j < nCols; ++j)
        {
	  p[j] = table[p[j]];
        }
    }
  return I;
}


void MatToBWAverage( Mat& I, Mat*& M ) {
  // new matrix to return
  M = new Mat( I.rows, I.cols / I.channels(), CV_8UC1, Scalar::all(0) );

  std::cout<<"Image depth is "<<I.depth()<<std::endl;
  int channels = I.channels();
  int nRows = I.rows;
  int nCols = I.cols * channels;

  std::cout<<"Image has "<<channels<<" channels "<<nRows<<" rows x "<<nCols<<" columns"<<std::endl;

  //uchar* p;
  uchar * m;
  for( int i = 0; i < nRows; ++i) {
      for ( int  j = 0; j < nCols/channels; ++j) {
	Vec3b rgb = I.at<Vec3b>( i, j*channels );
	M->at<uchar>( i, j ) = (rgb[0]+rgb[1]+rgb[2])/channels; 
      }
  }
}

Mat undistort_image( Mat& orig ){// Mat*& undistorted ) {
  // new matrix to return
  //*undistorted = orig.clone(); 
  //undistorted = new Mat( orig.rows, orig.cols , CV_8UC1, Scalar::all(0) );
  //*/undistorted=orig.clone();
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
  //cameraMatrix.at<double>( 0, 0 ) = 2.7605e+03;
  //cameraMatrix.at<double>( 1, 1 ) = 2.7670e+03;
  //cameraMatrix.at<double>( 0, 2 ) = 1.9143e+03;
  //cameraMatrix.at<double>( 1, 2 ) = 1.5964e+03;

  //distCoeffs.at<double>(0, 0) = -0.2398;
  //distCoeffs.at<double>(0, 1 ) = 0.1145;

  Matx33f cameraMatrix( 2.7605e+03, 0,          1.9143e+03, 
                        0,          2.7670e+03, 1.5964e+03, 
                        0,          0,          1 );
  Vec4f distCoeffs( -0.2398, 0.1145, 0, 0 ); 


  // std::cout<<
 
  //cv::fisheye::undistortImage(orig, undistorted, cameraMatrix, distCoeffs);
  undistort(orig, undistorted, cameraMatrix, distCoeffs);


  return undistorted;
  //
}


int main(int argc, char** argv )
{
    if ( argc != 2 )
    {
        printf("usage: DisplayImage.out <Image_Path>\n");
        return -1;
    }
    Mat image = imread( argv[1], IMREAD_COLOR );
    if ( !image.data )
    {
        printf("No image data \n");
        return -1;
    }
    //**Testing
    String orig_name = (String)argv[1];
    //std::cout << s << std::endl;
    //std::cout << s.String::substr(s.String::rfind('/')+1) << std::endl;
    //**

 // get black and white version
    Mat image_bw = imread( argv[1], IMREAD_GRAYSCALE );


    //uchar table[256];
    //int divideWith=10;
    //for (int i = 0; i < 256; ++i)
    //  table[i] = (uchar)(divideWith * (i/divideWith));

    //Mat img2 = ScanImageAndReduceC(image, table);
    //Mat * img3;
    //MatToBWAverage( image, img3 );

    //namedWindow("Display Image Orig", WINDOW_GUI_EXPANDED );
    //imshow("Display Image Orig", image);

    namedWindow("Display Image", WINDOW_NORMAL );
    imshow("Display Image", image_bw);

    //namedWindow("Display Image My BW", WINDOW_NORMAL );
    //imshow("Display Image My BW", *img3);


    // undistort image
    Mat  img_undistorted = undistort_image( image_bw );

    namedWindow("Display Undistorted Image", WINDOW_GUI_EXPANDED );
    imshow("Display Undistorted Image", img_undistorted );

    imwrite(orig_name.String::substr(orig_name.String::rfind('/')+1) , img_undistorted );

    waitKey(0);
    //delete img3;
    return 0;
}
