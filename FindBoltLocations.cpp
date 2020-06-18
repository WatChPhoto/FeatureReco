#include <iostream>
#include <opencv2/imgproc.hpp>
#include <opencv2/imgcodecs.hpp>
#include <TFile.h>
#include "Configuration.hpp"
#include "featureFunctions.hpp"
#include <opencv2/features2d.hpp>     //Blob

#include<cmath>
#define PI 3.141592653589793238
#define RADTODEG(R)((180.0 * R) / PI)

using std::string;
using std::vector;
 
using namespace cv;

//Returns the data from text file
//Returns empty vector if text file is not supplied.
MedianTextData assign(int argc, string argv){
  if(argc==2){MedianTextData a; return a;}
  if(argc==3){
    MedianTextReader *boltreader = MedianTextReader::Get();
    boltreader->set_input_file( string( argv ) );
    return boltreader->get_data();
  }
}

//flags to turn on/off saving images
vector<bool> setup(){
  vector <bool> options;
  try {
    int option = config::Get_int("save_option");
    
    for(int i=0; i<5;i++){
      options.push_back(bool(option%2));
      option /= 10;
    }

  } catch ( std::string e ){
      std::cout<<"Error with config file key "<<e<<std::endl;
    }
  return options;  
}

vector<float> get_angles(vector<Vec3f>final_PMTs, vector<Vec3f>serial_final_bolts){
  vector<float> angle;
  int index=0;
  for( Vec3f pmt: final_PMTs ){
    float a = pmt[0];  //x and y at centre
    float b = pmt[1];
    
    for(index; index<serial_final_bolts.size(); ++index){// Vec3f bolts: serial_final_bolts ){
      float x = serial_final_bolts[index][0];
      float y = serial_final_bolts[index][1];
      //Finds out if this is break point. i.e. bolt belonging to different pmt.
      /*=============testing======================
	if(index>23&&index<48){
	circle(image_final, Point(x,y), 5, Scalar(0,255,255), 3);
	}
	============================================*/
      //Marks the starting index for next pmt as next index.
      if(x==-1&&y==-1){++index; angle.push_back(-5); break;}
      float theta = atan2f((y-b),(x-a));
      
      theta = RADTODEG(theta);
      theta = (theta<0)?(theta+360):theta; //getting angle between 0-360
      //Finding angle wrt to Y-axis ^ (nothing to do with axis direction in OpenCv
      theta = (theta<270)?(theta+90):(theta-270);
      angle.push_back(theta);
      
    }
  }
  
  return angle;
}

int main(int argc, char** argv )
{
    
if ( argc != 2 && argc!=3 )
    {
        printf("usage: FindBoltLocations <Input_image_with_path> [<median-bolt-loc-filename>]\n");
        return -1;
    }
//mode 0 means 2 argument mode and mode 1 means 3 arg mode.
 bool mode = argc-2;

    Mat image_color = imread( argv[1],  IMREAD_COLOR ); //IMREAD_GRAYSCALE,
    if ( !image_color.data )
    {
        printf("No image data \n");
        return -1;
    }
    //option has final, text, candidate, circled, filters
    const vector<bool>& option =  setup();
    
    Mat image_final = image_color.clone();

    /// build output image
    Mat image;
    cvtColor( image_color, image, COLOR_RGBA2GRAY );

    // Open a root file to put histograms into
    TFile * fout = new TFile("FindBoltLocation.root","RECREATE");

    string outputname;
    
    // Gaussian blur
    Mat img_blur = image.clone();
    try {
      //bool debug = config::Get_int("debug");
    bool do_gaus_blur = (bool)config::Get_int("do_gaus_blur");

    if (do_gaus_blur){
      int blurpixels = config::Get_int("blurpixels");            // size of kernel in pixels (must be odd)
      double blursigma = config::Get_double("blursigma");        // sigma of gaussian in pixels
      GaussianBlur( image, img_blur, Size( blurpixels, blurpixels ), blursigma );
      if(option[4]){
      outputname = build_output_filename( argv[1], "gausblur" );
      imwrite( outputname, img_blur );
      }
    }

    // Bilateral filter
    bool do_bifilter = (bool)config::Get_int("do_bifilter");
    Mat img_flt = img_blur.clone();
    if (do_bifilter) {
      int d = config::Get_int("d"); // value 5-9 distance around each pixel to filter (must be odd)
      int sigColor = config::Get_int("sigColor"); // range of colours to call the same
      int sigSpace = config::Get_int("sigSpace"); // ???
      bilateralFilter ( image, img_flt, d, sigColor, sigSpace );
      if(option[4]){
      outputname = build_output_filename( argv[1], "bifilter" );
      imwrite( outputname, img_flt );
      }
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

      if(option[4]){
      outputname = build_output_filename( argv[1], "sobel" );
      imwrite( outputname, grad );
      }
    }

    //Canny edge detector.                                                                                                                   
    bool do_canny = (bool)config::Get_int("do_canny");                                              
    int thresh_low = config::Get_int("thresh_low"); //The gradient value below thresh_low will be discarded.                                  
    int thresh_high = config::Get_int("thresh_high"); //The gradient value above thresh_high will be used. The inbetween gradient is kept if \the edge is connected.                                                                                                                        
    Mat img_can = grad.clone();
    if ( do_canny ){
      Canny(grad, img_can, thresh_low, thresh_high );

      if(option[option.size()-1]){
      outputname = build_output_filename( argv[1], "canny" );
      imwrite(outputname,img_can);
      }
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
    if(option[3]){
    outputname = build_output_filename( argv[1], "blob" );
    imwrite( outputname, img_blob_map );        
    }
    // Make image that just has circle centers from blob detection
    Mat blob_circles = Mat::zeros( image.size(), image.type() );
    draw_found_center(blobs, blob_circles);
    
    if(option[2]){
    outputname = build_output_filename( argv[1], "blobCandidate" );
    imwrite(outputname, blob_circles);
    }
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

    if(option[3]){
    outputname = build_output_filename( argv[1], "hough" );
    imwrite( outputname, image_color );
    }
    /// Read in bolt locations
    //Returns empty vector if text file is not supplied.
    const MedianTextData& mtd = assign(argc, string(argv[argc-1]));
    //Debug information PMt
    //if(debug){
    for ( const MedianTextRecord & rec : mtd ){
      std::cout<< rec;
    }
    //}
    
//Only performs blob and hough analysis if there are three argument.    
if(mode){
    //Blob analysis
    TH1D* blob_metric_all  = new TH1D("Blob_Metric_all" , "Bolt metric for all circles ;Inside to outside Intensity ratio; Count",200, 0.0, 10.0);
    TH1D* blob_metric_good = new TH1D("Blob_Metric_good", "Bolt metric for matched circles ;Inside to outside Intensity ratio; Count",200,  0.0, 10.0);
    TH1D* blob_metric_bad  = new TH1D("Blob_Metric_bad" , "Bolt metric for un-matched circles  ;Inside to outside Intensity ratio; Count",200, 0.0, 10.0);

    TH2D* blob_metric_2d   = new TH2D("Blob_Metric_2d", " Bolt metric vs distance to bolt;  Distance to bolt (pixels); Bolt Metric", 51, -0.5, 49.5, 100, 0., 10 );
    TH1D *blob_dist = new TH1D("bolt_distance from blob","Distance to closest bolt using blob; distance (pixels); count/bin", 51, -0.5, 49.5);
    TH1D* text_to_blob_dist = new TH1D("blob_bolt_distance_wrt_text","Distance to closest bolt ; distance (pixels); Count/bin", 51, -0.5, 49.5);
    TH1D* blob_metric_inb  = new TH1D("blob_metric_inbetween" , "Bolt metric for inbetween circles ;Inside to outside Intensity ratio; Count",200, 0.0, 50.0 ); 
    
    
    //Find matches bewteen blobs found and bolts labelled in the text file
    vector< IndexMatchDist > blob_matches = find_closest_matches( blobs, mtd );    

    //Histogram minimum distance to the bolt from the boob
    make_bolt_dist_histogram( blob_matches, blob_dist );
    
    //Make bolt distance histogram with respect to text
    make_bolt_dist_histogram_wrt_txt( blobs, mtd, text_to_blob_dist );
 
    //Make histograms of all metric, good metric and bad metric. 
    make_bolt_metric_histograms( blobs, blob_matches, image, blob_metric_all, blob_metric_good, blob_metric_bad, blob_metric_2d );
   
    //Draws line between matches in a given image
    draw_line(blobs, blob_matches, mtd, img_blob_map );

    //Make a histogram for the metric of inbetween points
    histogram_inbetween(blobs, mtd, blob_matches, image, blob_metric_inb);
    
    //Hough analysis
    TH1D* hough_metric_all  = new TH1D("Hough_Metric_all" , "Bolt metric for all circles ;Inside to outside Intensity ratio; Count",200, 0., 10.);
    TH1D* hough_metric_good = new TH1D("Hough_Metric_good", "Bolt metric for matched circles ;Inside to outside Intensity ratio; Count",200, 0., 10.);
    TH1D* hough_metric_bad  = new TH1D("Hough_Metric_bad" , "Bolt metric for un-matched circles  ;Inside to outside Intensity ratio; Count",200, 0., 10.);
    TH2D* hough_metric_2d   = new TH2D("Hough_Metric_2d", " Bolt metric vs distance to bolt;  Distance to bolt (pixels); Bolt Metric", 51, -0.5, 49.5, 100, 0, 10 );
    TH1D* hough_dist = new TH1D("hough_bolt_distance","Distance to closest bolt Using Hough; distance (pixels); Count/bin", 51, -0.5, 49.5);
    TH1D* text_to_hough_dist = new TH1D("hough_bolt_distance_wrt_text","Distance to closest bolt ; distance (pixels); Count/bin", 51, -0.5, 49.5);
    TH1D* hough_metric_inb  = new TH1D("hough_metric_inbetween" , "Bolt metric for inbetween circles ;Inside to outside Intensity ratio; Count",200, 0., 10.); 

    /// Find matches between circles found and bolts labelled in text file
    vector< IndexMatchDist > bolt_matches = find_closest_matches( circles, mtd );    

    /// Make a histogram of closest distance between one of our circles and one of the text records    
    make_bolt_dist_histogram( bolt_matches, hough_dist );

    //Make a histogram of closest distance from one of text records to our circles.
    make_bolt_dist_histogram_wrt_txt( circles, mtd, text_to_hough_dist );

    //Make histograms of all metric, good metric and bad metric. 
    make_bolt_metric_histograms( circles, bolt_matches, image, hough_metric_all, hough_metric_good, hough_metric_bad, hough_metric_2d );
    
    //Draws line between matches in a given image
    draw_line(circles, bolt_matches, mtd, image_color );

    //Make a histogram for the metric of inbetween points
    histogram_inbetween(circles, mtd, bolt_matches, image, hough_metric_inb);
    }
   
    /*work on this
    //Make a histogram of closest distance from one of text record to found blob
    make_bolt_dist_histogram_wrt_text(blobs,mdt,img_blob);
    make_bolt_dist_histograms(blobs,mdt,blob_matches,img_blob,coloredimg);
    */

    /// Make image that just has circle centers from previous Hough Transform on it
    Mat img_circles = Mat::zeros( image.size(), image.type() );
    draw_found_center(circles, img_circles);
    
    if(option[2]){
    outputname = build_output_filename( argv[1], "houghCandidate" );
    imwrite( outputname, img_circles );
    }
    /// Look for circles of bolts
    vector<Vec3f> hough_circles_of_bolts;
    //Look for circles of bolts from blob
    vector<Vec3f> blob_circles_of_bolts;
    dp = config::Get_int("sec_hough_dp"); //if dp=1 , the accum has resolution of input image. If dp=2 , the accumulator has half as big width and height. 
    minDist = config::Get_int("sec_hough_minDist");  // min distance between circles
    param1 = config::Get_int("sec_hough_param1");   // threshold placed on image
    param2 = config::Get_int("sec_hough_param2");  // minimum accumulator value to call it a circle
    minR   = config::Get_int("sec_hough_minR");   //= 3 # minimum radius in pixels
    maxR   = config::Get_int("sec_hough_maxR");  // = 10 # maximum radius in pixels
    
    HoughCircles( blob_circles, blob_circles_of_bolts, HOUGH_GRADIENT, dp, minDist, param1, param2, minR, maxR );
    HoughCircles( img_circles, hough_circles_of_bolts, HOUGH_GRADIENT, dp, minDist, param1, param2, minR, maxR );

    //Overlays detected circles from second hough transfrom
    draw_circle_from_data(hough_circles_of_bolts, image_color, Scalar(255,102,255));
    draw_circle_from_data(blob_circles_of_bolts, img_blob_map, Scalar(255,102,255));
    
    //Draw text file circles
    //Only if there are three arguments
    if(mode){
    draw_text_circles(image_color, mtd);
    draw_text_circles(img_blob_map, mtd);
    }
    //save images
    if(option[1]){
    outputname = build_output_filename( argv[1], "hough_text" );
    imwrite( outputname, image_color );
    outputname = build_output_filename( argv[1], "blob_text" );
    imwrite( outputname, img_blob_map );
    }

    /*
      Take circles_of_blob to select which bolts are good blobs
      inputs: 
        vector<Vec3f> blobs;  // (x,y,r)
	vector<Vec3f> circles_of_blob;  // (x,y,r)
	
      output:
        vector<Vec3f> final_bolts; // bolt locations selected
     */
    // only accept PMTs that are more than some number of pixels away from edge of image
    const unsigned trim_pixels = 600;
    int ywidth = image.rows;
    int xwidth = image.cols;
    int xmin = trim_pixels;
    int xmax = xwidth-trim_pixels;
    int ymin = trim_pixels;
    int ymax = ywidth-trim_pixels;
    


    vector<Vec3f> final_bolts; // bolt locations selected
    vector<Vec3f> serial_final_bolts;
    vector<Vec3f> final_PMTs; // PMT circles selected
    vector<float> final_dists;
    TH1D* blob_dist2 = new TH1D("bolt_distance2","Distance to closest bolt using blob + second hough; distance (pixels); Count/bin", 51, -0.5, 49.5);
    
    // loop over circles_of_blob 
    // which is the circles of multiple bolts (multiple blobs making circle around PMT)
    for ( const Vec3f & pmtloc : blob_circles_of_bolts ) {
      // loop over bolts (blobs)0 to see if it is on the pmt circle
      //unsigned nbolts = 0;
      int pmtx = pmtloc[0];
      int pmty = pmtloc[1];
      if ( pmtx < xmin || pmtx > xmax || pmty < ymin || pmty > ymax ) continue;

      vector<Vec3f> bolts_on_this_pmt;
      for ( const Vec3f & boltloc : blobs ) {
	// calculate distance from the PMT circle to the bolt location
	// only add ones with distance less than (2?) pixels to add to bolts_on_this_pmt
	// count/and print them after
       
	float pmtr= pmtloc[2];
        int boltx = boltloc[0];
        int bolty = boltloc[1];
        float dist = std::sqrt(std::pow((pmtx-boltx),2)+std::pow((pmty-bolty),2));
        //if(dist>(pmtr-2) && dist<(pmtr+2)){
	if ( fabs( pmtr - dist ) < 6 ) {
          Vec3f temp;
          temp[0]= boltx;
          temp[1]= bolty;
          temp[2]= boltloc[2];
	  final_dists.push_back( fabs(dist-pmtr) );
          bolts_on_this_pmt.push_back(temp);
	  
	}
      }
      
      // add bolts_on_this_pmt to final_bolts if > some number (5?) of bolts match?
      if( bolts_on_this_pmt.size() > 9 ) {
	final_PMTs.push_back( pmtloc );
	final_bolts.insert(final_bolts.end(), bolts_on_this_pmt.begin(), bolts_on_this_pmt.end());

	//Its a way to know which bolts belong to particular PMT.               
        serial_final_bolts.insert(serial_final_bolts.end(), bolts_on_this_pmt.begin(), bolts_on_this_pmt.end());
	Vec3f temp;
        temp[0]= -1;
        temp[1]= -1;
        temp[2]= -1;
        serial_final_bolts.push_back( temp );
      }
    }
   
    //Fill bolb_dist2 histogram
    for(float final: final_dists){
      blob_dist2->Fill(final);
    }
    
   
    // Make image of final answer
    // add the circles for the PMTs selected
    draw_circle_from_data( final_PMTs, image_final, Scalar(255,102,255), 2);
    draw_circle_from_data( final_bolts, image_final, Scalar( 0, 0, 255), 3);
   
    if(mode){
      draw_text_circles( image_final, mtd );
      //Find bolt match
      vector< IndexMatchDist > final_matches = find_closest_matches( final_bolts, mtd );
      //Draw line
      draw_line(final_bolts, final_matches, mtd, image_final );
      //Distance histogram
      TH1D *final_dist = new TH1D("bolt_distance final","Distance to closest bolt; distance (pixels); count/bin", 51, -0.5, 49.5);
      make_bolt_dist_histogram(final_matches, final_dist);
      /*
    if(option[0]){
      circle(image_final, Point(final_PMTs[0][0],final_PMTs[0][1]), 5, Scalar(0,255,255), 3);
 circle(image_final, Point(serial_final_bolts[0][0],serial_final_bolts[0][1]), 5, Scalar(0,255,255), 3);
    outputname = build_output_filename( argv[1], "final" );
    imwrite( outputname, image_final );
    }
      */

      //template will be {angles of one pmt bolt,-5, angles of next pmt bolts}
      //it's a way to identify which angle belong to which bolt and which pmt.
      //since final_PMTs and seral_final_bolts are in sync this approach work.
      const vector<float>& angles = get_angles( final_PMTs, serial_final_bolts );
    
      for(float ang: angles){
	std::cout<<"angle = "<< ang<<std::endl;
      }

      //testing purpose
    if(option[0]){
      circle(image_final, Point(final_PMTs[1][0],final_PMTs[1][1]), 5, Scalar(0,255,255), 3);
      circle(image_final, Point(serial_final_bolts[36][0], serial_final_bolts[36][1]), 5, Scalar(0,255,255), 3);
    outputname = build_output_filename( argv[1], "final" );
    imwrite( outputname, image_final );
    }     


 std::cout<<"#########################################"<<std::endl;
      std::cout<<"final size "<<final_bolts.size()<<" match size "<<final_matches.size()<<" Mtd size "<<mtd.size()<<std::endl;
      std::cout<<"Pmt0 x "<<final_PMTs[0][0]<<" y "<<final_PMTs[0][1]<<" blob x "<<serial_final_bolts[0][0]<<" y "<<serial_final_bolts[0][1]<<std::endl;
    }
    std::cout<<atan(1)<<atan(-1)<<std::endl;
    } catch ( std::string e ){
      std::cout<<"Error with config file key "<<e<<std::endl;
    }
    
    fout->Write();
    fout->Close();

    return 0;
}
