#include <iostream>
#include <opencv2/imgproc.hpp>
#include <opencv2/imgcodecs.hpp>
#include <TFile.h>
#include "Configuration.hpp"
#include "featureFunctions.hpp"

//Something I want to look into.
//Blob                                                                                                                                      
#include <opencv2/features2d.hpp>
//   

using std::string;
using std::vector;
 
using namespace cv;

int main(int argc, char** argv )
{
    if ( argc != 3 )
    {
        printf("usage: FindBoltLocations <Input_image_with_path> <median-bolt-loc-filename>\n");
        return -1;
    }


    Mat image_color = imread( argv[1],  IMREAD_COLOR ); //IMREAD_GRAYSCALE,
    if ( !image_color.data )
    {
        printf("No image data \n");
        return -1;
    }

    /// build output image
    Mat image;
    cvtColor( image_color, image, COLOR_RGBA2GRAY );

    // Open a root file to put histograms into
    TFile * fout = new TFile("FindBoltLocation.root","RECREATE");


    string outputname;
    
    // Gaussian blur
    Mat img_blur = image.clone();
    try {
    bool do_gaus_blur = (bool)config::Get_int("do_gaus_blur");

    if (do_gaus_blur){
      int blurpixels = config::Get_int("blurpixels");            // size of kernel in pixels (must be odd)
      double blursigma = config::Get_double("blursigma");        // sigma of gaussian in pixels
      GaussianBlur( image, img_blur, Size( blurpixels, blurpixels ), blursigma );
      outputname = build_output_filename( argv[1], "gausblur" );
      imwrite( outputname, img_blur );
    }

    // Bilateral filter
    bool do_bifilter = (bool)config::Get_int("do_bifilter");
    Mat img_flt = img_blur.clone();
    if (do_bifilter) {
      int d = config::Get_int("d"); // value 5-9 distance around each pixel to filter (must be odd)
      int sigColor = config::Get_int("sigColor"); // range of colours to call the same
      int sigSpace = config::Get_int("sigSpace"); // ???
      bilateralFilter ( image, img_flt, d, sigColor, sigSpace );

      outputname = build_output_filename( argv[1], "bifilter" );
      imwrite( outputname, img_flt );
    }


    /// Do Sobel edge detection
    bool do_sobel = (bool)config::Get_int("do_sobel");
    Mat grad = img_flt.clone();
    if ( do_sobel ){
      int scale = config::Get_int("scale");
      int delta = config::Get_int("delta");
      int ddepth = CV_16S;
      
      /// Generate grad_x and grad_y
      
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
    }


    //Canny edge detector.                                                                                                                   
    bool do_canny = (bool)config::Get_int("do_canny");                                              
    int thresh_low = config::Get_int("thresh_low"); //The gradient value below thresh_low will be discarded.                                  
    int thresh_high = config::Get_int("thresh_high"); //The gradient value above thresh_high will be used. The inbetween gradient is kept if \the edge is connected.                                                                                                                        
    Mat img_can = grad.clone();
    if ( do_canny ){
      Canny(grad, img_can, thresh_low, thresh_high );

      outputname = build_output_filename( argv[1], "canny" );

      imwrite(outputname,img_can);

    }

    
    histogram_channel0( img_can, "hbefore_thres_cut" );
    int threshold = config::Get_int("threshold");
    apply_image_threshold( img_can, threshold );
    histogram_channel0( grad, "hafter_thres_cut");
   
    

    //Something I want to look into.
    //Blob detection
    Mat img_blob = img_flt.clone();
      // Setup SimpleBlobDetector parameters.
    SimpleBlobDetector::Params params;

    //detect white
    //params.filterByColor=true;  
    params.blobColor=255;

    // Change thresholds
    params.minThreshold = config::Get_int("blob_minThreshold");
    params.maxThreshold = config::Get_int("blob_maxThreshold");

    // Filter by Area.
     params.filterByArea = config::Get_int("blob_filterByArea");
     params.minArea = config::Get_double("blob_minArea");
     params.maxArea = config::Get_double("blob_maxArea");

    // Filter by Circularity
    params.filterByCircularity = config::Get_int("blob_filterByCircularity");
    params.minCircularity = config::Get_double("blob_minCircularity");

    //Filter by distance
    params.minDistBetweenBlobs = config::Get_double("blob_minDistBetweenBlobs");
    // Filter by Convexity
    params.filterByConvexity = config::Get_int("blob_filterByConvexity");
    params.minConvexity = config::Get_double("blob_minConvexity");

    // Filter by Inertia
    params.filterByInertia = config::Get_int("blob_filterByInertia");
    params.minInertiaRatio = config::Get_double("blob_minInertiaRatio");
 
    // Set up the detector with default parameters.                                                                                         
    Ptr<SimpleBlobDetector> detector = SimpleBlobDetector::create(params);

    // Detect blobs.                                                                                                                        
    std::vector<KeyPoint> keypoints;
    detector->detect( img_blob, keypoints);
    // Draw detected blobs as red circles.                                                                                                  
    // DrawMatchesFlags::DRAW_RICH_KEYPOINTS flag ensures the size of the circle corresponds to the size of blob                            
    //Mat im_with_keypoints;
    //drawKeypoints( img_blob, keypoints, im_with_keypoints, Scalar(0,0,255), DrawMatchesFlags::DRAW_RICH_KEYPOINTS );

    // Show blobs                                                                                                                           
    //imshow("keypoints", im_with_keypoints );                                                                                               
    //imwrite("keypoints.jpg",im_with_keypoints);
    

    /// Make image that just has circle centers from previous Hough Transform on it
    vector<Vec3f> blobs;
    int i=0;
    for(KeyPoint keypoint: keypoints){
      //Point center1 = keypoint.pt;
      int x = keypoint.pt.x;
      int y = keypoint.pt.y;
      float r = ((keypoint.size)+0.0)/2;
      Vec3f temp;
      temp[0]=x;
      temp[1]=y;
      temp[2]=r;
      blobs.push_back(temp);
    
      Point blob_center(x,y);
      circle( img_blob, blob_center, r, Scalar(0,255,0), 1, 8, 0 );
      ++i;
    }
    imwrite("blob.jpg", img_blob);
    
    Mat blob_circles = Mat::zeros( image.size(), image.type() );
    for( size_t i = 0; i < blobs.size(); i++ ) {
      blob_circles.at<uchar>(blobs[i][1], blobs[i][0]) = 255 ;
    }
    imwrite("blobcandidate.jpg", blob_circles);
    // Blobend               
    
    /*Something for later    
    //Contour detection
    Mat img_1 = image.clone();
    Mat dst = Mat::zeros(img_1.rows, img_1.cols, CV_8UC3);
    vector<Vec4i> hierarchy;
    vector<vector<Point> > contours;
    findContours( img_1, contours, hierarchy, RETR_CCOMP, CHAIN_APPROX_NONE );
    int idx = 0;
    // for( ; idx >= 0; idx = hierarchy[idx][0] )
    //{
        Scalar color( 255, 255, 255 );
        drawContours( dst, contours, -2, color );
	
	//}

    //findContours( img_1, contours, RETR_LIST, CHAIN_APPROX_SIMPLE);
    //drawContours(img_1, contours, -1, Scalar(0,255,0));
    imwrite("Contour.jpg", dst);
    //
    */

    /// Hough Transform
    vector<Vec3f> circles;
    int dp = config::Get_int("hough_dp"); // Inverse ratio of the accumulator resolution to the image resolution. For example, if dp=1 , the accumulator has the same resolution as the input image. If dp=2 , the accumulator has half as big width and height. 
    int minDist = config::Get_int("hough_minDist"); // min distance between circles
    int param1 = config::Get_int("hough_param1"); // threshold placed on image
    int param2 = config::Get_int("hough_param2"); // minimum accumulator value to call it a circle
    int minR   = config::Get_int("hough_minR"); //= 3 # minimum radius in pixels
    int maxR   = config::Get_int("hough_maxR"); // = 10 # maximum radius in pixels
    HoughCircles( img_can, circles, HOUGH_GRADIENT, dp, minDist, param1, param2, minR, maxR );
  
    for( size_t i = 0; i < circles.size(); i++ ) {
      Point center(cvRound(circles[i][0]), cvRound(circles[i][1]));
      int radius = cvRound(circles[i][2]);
      // draw the circle center
      //circle( image_color, center, 3, Scalar(0,0,255), 1, 8, 0 );
      // draw the circle outline
      circle( image_color, center, radius, Scalar(0,0,255), 1, 8, 0 );
      std::cout<<"Circle "<<i<<" radius = "<<radius<<" at ( "<<circles[i][0]<<", "<<circles[i][1]<<" )"<<std::endl;
    }
    
    outputname = build_output_filename( argv[1], "circles" );
    imwrite( outputname, image_color );
  
    /// Read in bolt locations
    MedianTextReader *boltreader = MedianTextReader::Get();
    boltreader->set_input_file( string( argv[2] ) );
    const MedianTextData& mtd = boltreader->get_data();
    for ( const MedianTextRecord & rec : mtd ){
      std::cout<< rec;
    }
    
    /// Find matches between circles found and bolts labelled in text file
    vector< IndexMatchDist > bolt_matches = find_closest_matches( circles, mtd );    

    /// Make a histogram of closest distance between one of our circles and one of the text records
    make_bolt_dist_histogram( bolt_matches );
    
    //Make a histogram of closest distance from one of text records to our circles.
    make_bolt_dist_histogram_wrt_txt( circles, mtd, image_color );
    make_bolt_metric_histograms( circles, mtd, bolt_matches, image, image_color );
    
    //Make a histogram for the metric of inbetween points
    histogram_inbetween(circles, mtd, bolt_matches, image);

    //Drawing Michel's point in the picture.
    for ( const MedianTextRecord & rec : mtd ){
      Point center_michel(rec.x(), rec.y());   
      //used radius of 10 pixels and green color.
      circle( image_color, center_michel, 10, Scalar(0,255,0), 1, 8, 0 );
    }
    //end of drawing michel's points.


    /// Make image that just has circle centers from previous Hough Transform on it
    Mat img_circles = Mat::zeros( image.size(), image.type() );
    for( size_t i = 0; i < circles.size(); i++ ) {
      img_circles.at<uchar>(cvRound(circles[i][1]), cvRound(circles[i][0])) = 255 ;
    }

    outputname = build_output_filename( argv[1], "candidatebolts" );
    imwrite( outputname, img_circles );

    
    /// Look for circles of bolts
    vector<Vec3f> circles_of_bolts;
    dp =      1; //if dp=1 , the accum has resolution of input image. If dp=2 , the accumulator has half as big width and height. 
    minDist = 216; //210; // min distance between circles
    param1 = 1; // threshold placed on image
    param2 = 3; //5; // minimum accumulator value to call it a circle
    minR   = 90; //80; //= 3 # minimum radius in pixels
    maxR   = 116; //120; // = 10 # maximum radius in pixels
    if(config::Get_int("do_blob")){
      HoughCircles( blob_circles, circles_of_bolts, HOUGH_GRADIENT, dp, minDist, param1, param2, minR, maxR );
      }
    if(config::Get_int("do_hough")){
	HoughCircles( img_circles, circles_of_bolts, HOUGH_GRADIENT, dp, minDist, param1, param2, minR, maxR );
      }    

    for( size_t i = 0; i < circles_of_bolts.size(); i++ ) {
      Point center(cvRound(circles_of_bolts[i][0]), cvRound(circles_of_bolts[i][1]));
      int radius = cvRound(circles_of_bolts[i][2]);
      // draw the circle center
      //circle( image_color, center, 3, Scalar(0,0,255), 1, 8, 0 );
      // draw the circle outline
      circle( image_color, center, radius, Scalar(255,102,255), 1, 8, 0 );
      std::cout<<"Circle "<<i<<" radius = "<<radius<<" at ( "<<circles[i][0]<<", "<<circles[i][1]<<" )"<<std::endl;
    }
    
    outputname = build_output_filename( argv[1], "michel" );
    imwrite( outputname, image_color );


    } catch ( std::string e ){
      std::cout<<"Error with config file key "<<e<<std::endl;
    }
    
    fout->Write();
    fout->Close();

    return 0;
}
