#include <iostream>
#include <opencv2/imgproc.hpp>
#include <opencv2/imgcodecs.hpp>
#include <TFile.h>
#include "Configuration.hpp"
#include "featureFunctions.hpp"
#include <opencv2/features2d.hpp>     //Blob

using std::string;
using std::vector;
 
using namespace cv;

void draw_text_circles(Mat &img, const MedianTextData& mtd){
for ( const MedianTextRecord & rec : mtd ){
      Point center_michel(rec.x(), rec.y());   
      //used radius of 10 pixels and green color.
      circle( img, center_michel, 10, Scalar(0,255,0), 1, 8, 0 );
    }
}

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
   
    //Blob detection
    Mat img_blob = image.clone();
    Mat img_blob_map = image_color.clone();
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

    /*
    // Draw detected blobs as red circles.                                                                                                  
    // DrawMatchesFlags::DRAW_RICH_KEYPOINTS flag ensures the size of the circle corresponds to the size of blob                            
    Mat im_with_keypoints;
    drawKeypoints( img_blob, keypoints, im_with_keypoints, Scalar(0,0,255), DrawMatchesFlags::DRAW_RICH_KEYPOINTS );

    imwrite("keypoints.jpg",im_with_keypoints);
    */
    
    //blob vector will contain x,y,r
    vector<Vec3f> blobs;
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
      }
      
    //Draws circle from data to the input image
    draw_circle_from_data(blobs, img_blob_map, Scalar(0,0,255));
    outputname = build_output_filename( argv[1], "blob" );
    imwrite( outputname, img_blob_map );        

    // Make image that just has circle centers from blob detection
    Mat blob_circles = Mat::zeros( image.size(), image.type() );
    draw_found_center(blobs, blob_circles);
    outputname = build_output_filename( argv[1], "blobCandidate" );
    imwrite(outputname, blob_circles);
    // Blobend               
    
    /// Hough Transform
    vector<Vec3f> circles;
    int dp = config::Get_int("hough_dp"); // Inverse ratio of the accumulator resolution to the image resolution. For example, if dp=1 , the accumulator has the same resolution as the input image. If dp=2 , the accumulator has half as big width and height. 
    int minDist = config::Get_int("hough_minDist"); // min distance between circles
    int param1 = config::Get_int("hough_param1"); // threshold placed on image
    int param2 = config::Get_int("hough_param2"); // minimum accumulator value to call it a circle
    int minR   = config::Get_int("hough_minR"); //= 3 # minimum radius in pixels
    int maxR   = config::Get_int("hough_maxR"); // = 10 # maximum radius in pixels
    HoughCircles( img_can, circles, HOUGH_GRADIENT, dp, minDist, param1, param2, minR, maxR );
  
    draw_circle_from_data(circles, image_color, Scalar(0,0,255));

    outputname = build_output_filename( argv[1], "circles" );
    imwrite( outputname, image_color );
  
    /// Read in bolt locations
    MedianTextReader *boltreader = MedianTextReader::Get();
    boltreader->set_input_file( string( argv[2] ) );
    const MedianTextData& mtd = boltreader->get_data();
    for ( const MedianTextRecord & rec : mtd ){
      std::cout<< rec;
    }

    //Blob analysis
    TH1D* blob_metric_all  = new TH1D("Blob_Metric_all" , "Bolt metric for all circles ;Inside to outside Intensity ratio; Count",502, -1.5, 255.5);
    TH1D* blob_metric_good = new TH1D("Blob_Metric_good", "Bolt metric for matched circles ;Inside to outside Intensity ratio; Count",502, -1.5, 255.5);
    TH1D* blob_metric_bad  = new TH1D("Blob_Metric_bad" , "Bolt metric for un-matched circles  ;Inside to outside Intensity ratio; Count",502, -1.5, 255.5);

    TH2D* blob_metric_2d   = new TH2D("Blob_Metric_2d", " Bolt metric vs distance to bolt;  Distance to bolt (pixels); Bolt Metric", 100, -0.5, 99.5, 101, -1.5, 20 );
    TH1D *blob_dist = new TH1D("bolt_distance from blob","Distance to closest bolt using blob; distance (pixels); count/bin", 1001, 0.5, 500.5);
    TH1D* text_to_blob_dist = new TH1D("blob_bolt_distance_wrt_text","Distance to closest bolt ; distance (pixels); Count/bin",501, -0.5, 500.5);
    TH1D* blob_metric_inb  = new TH1D("blob_metric_inbetween" , "Bolt metric for inbetween circles ;Inside to outside Intensity ratio; Count",502, -1.5, 255.5); 

    
    //Find matches bewteen blobs found and bolts labelled in the text file
    vector< IndexMatchDist > blob_matches = find_closest_matches( blobs, mtd );    

    //Histogram minimum distance to the bolt from the boob
    make_bolt_dist_histogram( blob_matches, blob_dist );
    
    //Make bolt distance histogram with respect to text
    make_bolt_dist_histogram_wrt_txt( blobs, mtd, text_to_blob_dist );
 
    //Make histograms of all metric, good metric and bad metric. 
    make_bolt_metric_histograms( blobs, mtd, blob_matches, image, img_blob_map, blob_metric_all, blob_metric_good, blob_metric_bad, blob_metric_2d );
    //Make a histogram for the metric of inbetween points
    histogram_inbetween(blobs, mtd, blob_matches, image, blob_metric_inb);
    
    //Hough analysis
    TH1D* hough_metric_all  = new TH1D("Hough_Metric_all" , "Bolt metric for all circles ;Inside to outside Intensity ratio; Count",502, -1.5, 255.5);
    TH1D* hough_metric_good = new TH1D("Hough_Metric_good", "Bolt metric for matched circles ;Inside to outside Intensity ratio; Count",502, -1.5, 255.5);
    TH1D* hough_metric_bad  = new TH1D("Hough_Metric_bad" , "Bolt metric for un-matched circles  ;Inside to outside Intensity ratio; Count",502, -1.5, 255.5);
    TH2D* hough_metric_2d   = new TH2D("Hough_Metric_2d", " Bolt metric vs distance to bolt;  Distance to bolt (pixels); Bolt Metric", 100, -0.5, 99.5, 101, -1.5, 20 );
    TH1D* hough_dist = new TH1D("bolt_distance","Distance to closest bolt Using Hough; distance (pixels); Count/bin",1001, 0.5, 500.5);
    TH1D* text_to_hough_dist = new TH1D("hough_bolt_distance_wrt_text","Distance to closest bolt ; distance (pixels); Count/bin",501, -0.5, 500.5);
    TH1D* hough_metric_inb  = new TH1D("hough_metric_inbetween" , "Bolt metric for inbetween circles ;Inside to outside Intensity ratio; Count",502, -1.5, 255.5); 

    /// Find matches between circles found and bolts labelled in text file
    vector< IndexMatchDist > bolt_matches = find_closest_matches( circles, mtd );    

    /// Make a histogram of closest distance between one of our circles and one of the text records    
    make_bolt_dist_histogram( bolt_matches, hough_dist );

    //Make a histogram of closest distance from one of text records to our circles.
    make_bolt_dist_histogram_wrt_txt( circles, mtd, text_to_hough_dist );

    //Make histograms of all metric, good metric and bad metric. 
    make_bolt_metric_histograms( circles, mtd, bolt_matches, image, image_color, hough_metric_all, hough_metric_good, hough_metric_bad, hough_metric_2d );
    
    //Make a histogram for the metric of inbetween points
    histogram_inbetween(circles, mtd, bolt_matches, image, hough_metric_inb);

   
    /*work on this
    //Make a histogram of closest distance from one of text record to found blob
    make_bolt_dist_histogram_wrt_text(blobs,mdt,img_blob);
    make_bolt_dist_histograms(blobs,mdt,blob_matches,img_blob,coloredimg);
    */

    /// Make image that just has circle centers from previous Hough Transform on it
    Mat img_circles = Mat::zeros( image.size(), image.type() );
    draw_found_center(circles, img_circles);
    
    outputname = build_output_filename( argv[1], "houghCandidate" );
    imwrite( outputname, img_circles );
 
    /// Look for circles of bolts
    vector<Vec3f> circles_of_bolts;
    //Look for circles of bolts from blob
    vector<Vec3f> circles_of_blob;
    dp = config::Get_int("sec_hough_dp"); //if dp=1 , the accum has resolution of input image. If dp=2 , the accumulator has half as big width and height. 
    minDist = config::Get_int("sec_hough_minDist");  // min distance between circles
    param1 = config::Get_int("sec_hough_param1");   // threshold placed on image
    param2 = config::Get_int("sec_hough_param2");  // minimum accumulator value to call it a circle
    minR   = config::Get_int("sec_hough_minR");   //= 3 # minimum radius in pixels
    maxR   = config::Get_int("sec_hough_maxR");  // = 10 # maximum radius in pixels
    
    HoughCircles( blob_circles, circles_of_blob, HOUGH_GRADIENT, dp, minDist, param1, param2, minR, maxR );
    HoughCircles( img_circles, circles_of_bolts, HOUGH_GRADIENT, dp, minDist, param1, param2, minR, maxR );

    //Overlays detected circles from second hough transfrom
    draw_circle_from_data(circles_of_bolts, image_color, Scalar(255,102,255));
    draw_circle_from_data(circles_of_blob, img_blob_map, Scalar(255,102,255));
    
    //Draw text file circles
    draw_text_circles(image_color, mtd);
    draw_text_circles(img_blob_map, mtd);
    
    //save images
    outputname = build_output_filename( argv[1], "hough_text" );
    imwrite( outputname, image_color );
    outputname = build_output_filename( argv[1], "blob_text" );
    imwrite( outputname, img_blob_map );

    } catch ( std::string e ){
      std::cout<<"Error with config file key "<<e<<std::endl;
    }
    
    fout->Write();
    fout->Close();

    return 0;
}
