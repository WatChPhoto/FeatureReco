#include <stdio.h>
#include <opencv2/opencv.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/imgcodecs.hpp>

#include <string>

using std::string;
 
using namespace cv;


string build_output_filename( const string& in, const string& tag ){
    string outputname;
    size_t idx = in.find_last_of("/");
    outputname = string( tag ) + in.substr(idx+1 );
    return outputname;
}


int main(int argc, char** argv )
{
    if ( argc != 2 && argc != 3 )
    {
        printf("usage: FindBoltLocations <Input_image_with_path> [Ouptut_image_with_path]\n");
        return -1;
    }


    Mat image = imread( argv[1], IMREAD_COLOR );
    if ( !image.data )
    {
        printf("No image data \n");
        return -1;
    }

    string outputname;
    if ( argc == 2 ) outputname = build_output_filename( argv[1], "gausblur" );
    else outputname  = string( argv[2] );


    // Gaussian blur
    Mat img_blur = image.clone();
    int blurpixels = 51;            // size of kernel in pixels (must be odd)
    double blursigma = 10.0;        // sigma of gaussian in pixels
    GaussianBlur( image, img_blur, Size( blurpixels, blurpixels ), blursigma );
    imwrite( outputname, img_blur );
    

    // Bilateral filter
    if ( argc == 2 ) outputname = build_output_filename( argv[1], "bifilter" );
    else outputname  = string( argv[2] );

    Mat img_flt = image.clone();
    int d = 75; // value 5-9 distance around each pixel to filter (must be odd)
    int sigColor = d*2; // range of colours to call the same
    int sigSpace = d/2; // ???
    bilateralFilter ( image, img_flt, d, sigColor, sigSpace );

    imwrite( outputname, img_flt );



    return 0;
}
