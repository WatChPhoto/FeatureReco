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
    if ( argc != 2 )
    {
        printf("usage: FindBoltLocations <Input_image_with_path>\n");
        return -1;
    }


    Mat image = imread( argv[1], IMREAD_COLOR );
    if ( !image.data )
    {
        printf("No image data \n");
        return -1;
    }

    string outputname;
    outputname = build_output_filename( argv[1], "gausblur" );
    
    // Gaussian blur
    Mat img_blur = image.clone();
    int blurpixels = 11;            // size of kernel in pixels (must be odd)
    double blursigma = 10.0;        // sigma of gaussian in pixels
    GaussianBlur( image, img_blur, Size( blurpixels, blurpixels ), blursigma );
    imwrite( outputname, img_blur );
    

    // Bilateral filter
    outputname = build_output_filename( argv[1], "bifilter" );

    Mat img_flt = image.clone();
    int d = 9; // value 5-9 distance around each pixel to filter (must be odd)
    int sigColor = d*2; // range of colours to call the same
    int sigSpace = d/2; // ???
    bilateralFilter ( image, img_flt, d, sigColor, sigSpace );

    imwrite( outputname, img_flt );



    /// Do Sobel edge detection

    int scale = 1;
    int delta = 0;
    int ddepth = CV_16S;

    /// Generate grad_x and grad_y
    Mat grad;
    Mat grad_x, grad_y;
    Mat abs_grad_x, abs_grad_y;

    /// Gradient X
    //Scharr( src_gray, grad_x, ddepth, 1, 0, scale, delta, BORDER_DEFAULT );
    Sobel( img_flt, grad_x, ddepth, 1, 0, 3, scale, delta, BORDER_DEFAULT );
    convertScaleAbs( grad_x, abs_grad_x );

    /// Gradient Y
    //Scharr( src_gray, grad_y, ddepth, 0, 1, scale, delta, BORDER_DEFAULT );
    Sobel( img_flt, grad_y, ddepth, 0, 1, 3, scale, delta, BORDER_DEFAULT );
    convertScaleAbs( grad_y, abs_grad_y );

    /// Total Gradient (approximate)
    addWeighted( abs_grad_x, 0.5, abs_grad_y, 0.5, 0, grad );

    outputname = build_output_filename( argv[1], "sobel" );
    imwrite( outputname, grad );


    return 0;
}
