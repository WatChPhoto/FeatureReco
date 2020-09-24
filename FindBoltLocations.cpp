#include <iostream>
#include <fstream>
#include <opencv2/imgproc.hpp>
#include <opencv2/imgcodecs.hpp>
#include <TFile.h>
#include "Configuration.hpp"
#include "featureFunctions.hpp"
#include <opencv2/features2d.hpp>	//Blob

#include "PMTIdentified.hpp"
#include "bwlabel.hpp"

#include<cmath>

#include "hough_ellipse.hpp"

#include "ellipse_detection_2.1.hpp" //ellipse detection fast


using std::string;
using std::vector;

using namespace cv;

// RGB to CMYK conversion
void rgb2cmyk(const cv::Mat& img, std::vector<cv::Mat>& cmyk ) {
  // Allocate cmyk to store 4f componets
  for (int i = 0; i < 4; i++) {
    cmyk.push_back(cv::Mat(img.size(), CV_8UC1));
  }

  // Get rgb
  std::vector<cv::Mat> rgb;
  cv::split(img, rgb);

  // rgb to cmyk
  for (int i = 0; i < img.rows; i++) {
    for (int j = 0; j < img.cols; j++) {
      float r = (int)rgb[2].at<uchar>(i, j) / 255.;
      float g = (int)rgb[1].at<uchar>(i, j) / 255.;
      float b = (int)rgb[0].at<uchar>(i, j) / 255.;
      float k = std::min(std::min(1- r, 1- g), 1- b);         
 
      cmyk[0].at<uchar>(i, j) = (1 - r - k) / (1 - k) * 255.;
      cmyk[1].at<uchar>(i, j) = (1 - g - k) / (1 - k) * 255.;
      cmyk[2].at<uchar>(i, j) = (1 - b - k) / (1 - k) * 255.;
      cmyk[3].at<uchar>(i, j) = 255 - k * 255.;
    }
  }
}


/// return image of color by index
/// ret_rgb==true, indicies go red, green, blue
/// ret_rgb==false, indices go cyan, magenta, yellow, K
Mat output_image_by_color( const cv::Mat& image, bool ret_rgb, unsigned index, std::string infname, bool write_images=false ){

  std::vector<cv::Mat> cmyk;
  rgb2cmyk( image, cmyk );

  std::vector<cv::Mat> rgb;
  cv::split( image, rgb );

  TH1D* red     = new TH1D("color_red"    ,"Red; intensity; pixel count",256,-0.5,255.5);
  TH1D* green   = new TH1D("color_green"  ,"Green; intensity; pixel count",256,-0.5,255.5);
  TH1D* blue    = new TH1D("color_blue"   ,"Blue; intensity; pixel count",256,-0.5,255.5);
  TH1D* cyan    = new TH1D("color_cyan"   ,"Cyan; intensity; pixel count",256,-0.5,255.5);
  TH1D* magenta = new TH1D("color_magenta","Magenta; intensity; pixel count",256,-0.5,255.5);
  TH1D* yellow  = new TH1D("color_yellow" ,"Yellow; intensity; pixel count",256,-0.5,255.5);
  TH1D* K       = new TH1D("color_K"      ,"K; intensity; pixel count",256,-0.5,255.5);

  std::vector< TH1D* > hrgb = {red, green, blue};
  std::vector< TH1D* > hcmyk = {cyan, magenta, yellow, K};

  for (unsigned icol=0; icol<rgb.size(); ++icol ){
    const Mat& img = rgb[ icol ];
    
    for (int i = 0; i < img.rows; i++) {
      for (int j = 0; j < img.cols; j++) {
	int intensity = img.at<uchar>(i, j);
	hrgb[icol]->Fill( intensity );
      }
    }
  }

  for (unsigned icol=0; icol<cmyk.size(); ++icol ){
    const Mat& img = cmyk[ icol ];
    
    for (int i = 0; i < img.rows; i++) {
      for (int j = 0; j < img.cols; j++) {
	int intensity = img.at<uchar>(i, j);
	hcmyk[icol]->Fill( intensity );
      }
    }
  }
  
  std::vector< string > crgb{"red","green","blue"};
  std::vector< string > ccmyk{"cyan","magenta","yellow","K"};
  
  if (write_images){
    for ( unsigned i=0; i<rgb.size(); ++i ){
      string outputname = build_output_filename (infname, crgb[i]);
      imwrite ( outputname, rgb[i] );
    }
    
    for ( unsigned i=0; i<cmyk.size(); ++i ){
      string outputname = build_output_filename (infname, ccmyk[i]);
      imwrite ( outputname, cmyk[i] );
    }
  }
  if (ret_rgb) return rgb[index];
  return cmyk[index];

}

/// equalize each color and then get K image to return
Mat equalize_by_color( const cv::Mat& image, std::string infname, bool write_images=false ){


  int gridsize = config::Get_int( "clahe_gridsize" );
  int cliplimit = config::Get_int( "clahe_cliplimit" );
  
  Ptr<CLAHE> clahe = createCLAHE();
  clahe->setClipLimit( cliplimit );
  clahe->setTilesGridSize( cv::Size( gridsize, gridsize ) );

  std::vector<cv::Mat> rgb;
  cv::split( image, rgb );

  std::vector<cv::Mat> rgb_eq = rgb;
  for ( unsigned i=0; i<rgb.size(); ++i ){
    clahe->apply( rgb[0], rgb_eq[0] );
  }
  Mat rgb_rejoined;
  cv::merge( &rgb_eq[0], 3, rgb_rejoined );

  std::vector<cv::Mat> cmyk;
  rgb2cmyk( rgb_rejoined, cmyk );


  TH1D* red     = new TH1D("color_redeq"    ,"Red (equalized); intensity; pixel count",256,-0.5,255.5);
  TH1D* green   = new TH1D("color_greeneq"  ,"Green (equalized); intensity; pixel count",256,-0.5,255.5);
  TH1D* blue    = new TH1D("color_blueeq"   ,"Blue (equalized); intensity; pixel count",256,-0.5,255.5);
  TH1D* cyan    = new TH1D("color_cyaneq"   ,"Cyan (equalized); intensity; pixel count",256,-0.5,255.5);
  TH1D* magenta = new TH1D("color_magentaeq","Magenta (equalized); intensity; pixel count",256,-0.5,255.5);
  TH1D* yellow  = new TH1D("color_yelloweq" ,"Yellow (equalized); intensity; pixel count",256,-0.5,255.5);
  TH1D* K       = new TH1D("color_Keq"      ,"K (equalized); intensity; pixel count",256,-0.5,255.5);

  std::vector< TH1D* > hrgb = {red, green, blue};
  std::vector< TH1D* > hcmyk = {cyan, magenta, yellow, K};

  for (unsigned icol=0; icol<rgb.size(); ++icol ){
    const Mat& img = rgb_eq[ icol ];
    
    for (int i = 0; i < img.rows; i++) {
      for (int j = 0; j < img.cols; j++) {
	int intensity = img.at<uchar>(i, j);
	hrgb[icol]->Fill( intensity );
      }
    }
  }

  for (unsigned icol=0; icol<cmyk.size(); ++icol ){
    const Mat& img = cmyk[ icol ];
    
    for (int i = 0; i < img.rows; i++) {
      for (int j = 0; j < img.cols; j++) {
	int intensity = img.at<uchar>(i, j);
	hcmyk[icol]->Fill( intensity );
      }
    }
  }



  std::vector< string > crgb{"redeq","greeneq","blueeq"};
  std::vector< string > ccmyk{"cyaneq","magentaeq","yelloweq","Keq"};
  
  if (write_images){
    {
      string outputname = build_output_filename (infname, "color_equalized");
      imwrite ( outputname, rgb_rejoined );
    }

    for ( unsigned i=0; i<rgb.size(); ++i ){
      string outputname = build_output_filename (infname, crgb[i]);
      imwrite ( outputname, rgb_eq[i] );
    }
    
    for ( unsigned i=0; i<cmyk.size(); ++i ){
      string outputname = build_output_filename (infname, ccmyk[i]);
      imwrite ( outputname, cmyk[i] );
    }
  }

  
  return cmyk[3];
}


/// input is image_clahe (grayscale/bw image)
/// draw blobs on blob_circles (an empty Mat object)
/// writes image if write_images is set to false
vector< Vec3f > blob_detect( const Mat& image_clahe , Mat& blob_circles, bool write_images, const std::string& infname,   
			     std::vector < KeyPoint >& keypoints ){
  //Blob detection
  Mat img_blob = image_clahe.clone ();

  // Setup SimpleBlobDetector parameters.
  SimpleBlobDetector::Params params;
  
  //detect white
  //params.filterByColor=true;  
  params.blobColor = 255;
  
  // Change thresholds
  params.minThreshold = config::Get_int ("blob_minThreshold");
  params.maxThreshold = config::Get_int ("blob_maxThreshold");
  
  // Filter by Area.
  params.filterByArea = config::Get_int ("blob_filterByArea");
  params.minArea = config::Get_double ("blob_minArea");
  params.maxArea = config::Get_double ("blob_maxArea");
  
  // Filter by Circularity
  params.filterByCircularity = config::Get_int ("blob_filterByCircularity");
  params.minCircularity = config::Get_double ("blob_minCircularity");
  
  //Filter by distance
  params.minDistBetweenBlobs = config::Get_double ("blob_minDistBetweenBlobs");
  // Filter by Convexity
  params.filterByConvexity = config::Get_int ("blob_filterByConvexity");
  params.minConvexity = config::Get_double ("blob_minConvexity");
  
  // Filter by Inertia
  params.filterByInertia = config::Get_int ("blob_filterByInertia");
  params.minInertiaRatio = config::Get_double ("blob_minInertiaRatio");
  
  // Set up the detector with set parameters.                                                                                         
  Ptr < SimpleBlobDetector > detector = SimpleBlobDetector::create (params);
  
  // Detect blobs.                                                                                                                        
  detector->detect (img_blob, keypoints);
  
  //blob vector will contain x,y,r
  vector < Vec3f > blobs;
  for (KeyPoint keypoint:keypoints) {
    //Point center1 = keypoint.pt;
    int x = keypoint.pt.x;
    int y = keypoint.pt.y;
    float r = ((keypoint.size) + 0.0) / 2;
    
    Vec3f temp;
    temp[0] = x;
    temp[1] = y;
    temp[2] = r;
    blobs.push_back (temp);
  }
  
  //draw_found_center (blobs, blob_circles);
  blob_circles = Mat::zeros (image_clahe.size (), image_clahe.type ());
  draw_foundblobs( blobs, blob_circles );
 
  if ( write_images ) {
    //Draws circle from data to the input image
    Mat img_blob_map = image_clahe.clone();
    draw_circle_from_data (blobs, img_blob_map, Scalar (0, 0, 255));
    string outputname = build_output_filename (infname, "blob");
    imwrite (outputname, img_blob_map);
  
    // Make image that just has circle centers from blob detection
    outputname = build_output_filename (infname, "blobCandidate");
    imwrite (outputname, blob_circles);
  }

  return blobs;
}


Mat apply_gaussian_blur( const Mat& image, bool write_image, const std::string& infname ){
  Mat img_blur = image.clone ();
  bool verbose = config::Get_int("verbosity");
  bool do_gaus_blur = (bool) config::Get_int ("do_gaus_blur");

  if ( do_gaus_blur ) {
    int blurpixels = config::Get_int ("blurpixels");	// size of kernel in pixels (must be odd)
    double blursigma = config::Get_double ("blursigma");	// sigma of gaussian in pixels
	    
    GaussianBlur( image, img_blur, Size(blurpixels, blurpixels), blursigma);
    if ( write_image ) {
      string outputname = build_output_filename (infname, "gausblur");
      imwrite (outputname, img_blur);
    }
  }
  return img_blur;
}
  

Mat  apply_bilateral_filter( const Mat& img_blur, bool write_image, const std::string& infname ){
  bool do_bifilter = (bool) config::Get_int ("do_bifilter");
  Mat img_flt = img_blur.clone ();
  if ( do_bifilter ) {
    int d = config::Get_int ("d");	// value 5-9 distance around each pixel to filter (must be odd)
    int sigColor = config::Get_int ("sigColor");	// range of colours to call the same
    int sigSpace = config::Get_int ("sigSpace");	// ???
	    
    bilateralFilter (img_blur, img_flt, d, sigColor, sigSpace);
    if (write_image) {
      string outputname = build_output_filename (infname, "bifilter");
      imwrite (outputname, img_flt);
    }
  }
  return img_flt;
}


Mat apply_sobel_edge( const Mat& img_flt, bool write_image, const std::string& infname ){
  bool do_sobel = (bool) config::Get_int ("do_sobel");
  Mat grad = img_flt.clone ();
  if ( do_sobel ) {
    int scale = config::Get_int ("scale");
    int delta = config::Get_int ("delta");
    int ddepth = CV_16S;
    
    /// Generate grad_x and grad_y
    Mat grad_x, grad_y;
    Mat abs_grad_x, abs_grad_y;
    
    /// Gradient X
    //Scharr( src_gray, grad_x, ddepth, 1, 0, scale, delta, BORDER_DEFAULT );
    Sobel (img_flt, grad_x, ddepth, 1, 0, 3, scale, delta, BORDER_DEFAULT);
    convertScaleAbs (grad_x, abs_grad_x);

    /// Gradient Y
    //Scharr( src_gray, grad_y, ddepth, 0, 1, scale, delta, BORDER_DEFAULT );
    Sobel (img_flt, grad_y, ddepth, 0, 1, 3, scale, delta, BORDER_DEFAULT);
    convertScaleAbs (grad_y, abs_grad_y);
    
    /// Total Gradient (approximate)
    addWeighted (abs_grad_x, 0.5, abs_grad_y, 0.5, 0, grad);
    
    if ( write_image ) {
      string outputname = build_output_filename (infname, "sobel");
      imwrite (outputname, grad);
    }
  }
  return grad;
}


Mat apply_canny_edge( const Mat& grad, bool write_image, const std::string& infname ){
  
  bool do_canny = (bool) config::Get_int ("do_canny");
  int thresh_low = config::Get_int ("thresh_low");	//The gradient value below thresh_low will be discarded.                                  
  int thresh_high = config::Get_int ("thresh_high");	//The gradient value above thresh_high will be used. The inbetween gradient is kept if \the edge is connected.                                                                                                                        

  Mat img_can = grad.clone();
  if (do_canny) {
    Canny (grad, img_can, thresh_low, thresh_high);

    if ( write_image ) {
      string outputname = build_output_filename (infname, "canny");
      imwrite (outputname, img_can);
    }
  }
  return img_can;
}



// equalize image using clahe
Mat apply_clahe( const Mat& img_can, bool write_image, const std::string& infname ){
  bool do_clahe = (bool)config::Get_int( "do_clahe" );
  Mat image_clahe;
  if ( do_clahe ){
    std::cout<<"Applying equalization"<<std::endl;
    int gridsize = config::Get_int( "clahe_gridsize" );
    int cliplimit = config::Get_int( "clahe_cliplimit" );
	
    Ptr<CLAHE> clahe = createCLAHE();
    clahe->setClipLimit( cliplimit );
    clahe->setTilesGridSize( cv::Size( gridsize, gridsize ) );
    clahe->apply( img_can, image_clahe );
    std::cout<<"Equalized"<<std::endl;

    if ( write_image ) {
      string outputname = build_output_filename (infname, "clahe");
      imwrite (outputname, image_clahe);
    }
  } else {
    image_clahe = img_can.clone();
  }
  return image_clahe;
}
  

void fast_ellipse_detection( const vector<Vec3f > & blobs, Mat& image_ellipse, bool write_image, const std::string& infname, const MedianTextData & mtd ){
  int de_min_major = config::Get_int( "de_min_major" );//80;
  int de_max_major = config::Get_int( "de_max_major" );//160;
  int de_min_minor = config::Get_int( "de_min_minor" );//80;
  int de_max_minor = config::Get_int( "de_max_minor" );//160;
  int de_threshold = config::Get_int( "de_threshold" );//4;
  std::cout<<"before ellipse"<<std::endl;

  //Fast ellipse detection
  std::vector<ParametricEllipse> f_ellipses= detect_ellipse(blobs, 
							    de_min_major, de_max_major,
							    de_min_minor, de_max_minor, de_threshold );

  //Filling PMTIdentified vector. Data is obtained from Fast ellipse detection.
  std::vector< PMTIdentified > f_ellipse_pmts;
  for( ParametricEllipse ellipses: f_ellipses ) {	    
    // a = b / sqrt( 1-e^2 )
    // b/a = sqrt( 1-e^2 )
    // (b/a)^2 = 1 - e^2
    // e^2 = 1-(b/a)^2
    // e = sqrt( 1- (b/a)^2 )
    float mya = ellipses.a;
    float myb = ellipses.b;
    float myphi = ellipses.alpha;
    if ( myb > mya ){
      float tmp = mya;
      mya = myb;
      myb = tmp;
      myphi += pi/2;
    }
    ellipse_st pmtloc{ myb,
	std::sqrt( 1 - (myb/mya)*(myb/mya) ),
	myphi,
	xypoint( ellipses.centre.x, ellipses.centre.y)
	};
    std::vector< Vec3f > boltlocs = ellipses.bolts;
    std::vector< float > dists = ellipses.dist;
    
    f_ellipse_pmts.push_back( PMTIdentified( pmtloc, boltlocs, dists, 0 ) );
  }
	 

  //Angle made by bolts with vertical(12 O' clock)
  TH1D * fhangboltel = new TH1D ("fhangboltel", "Angles of bolts (fast ellipse); angle (deg)", 360, 0., 360.);
  //Difference in angle from expected angle of the bolt.
  TH1D * fhdangboltel = new TH1D ("fhdangboltel", "Angle of bolt from expected (fast ellipse); #Delta angle (deg)", 60, -15., 15.);
  
  for  (const PMTIdentified & pmt : f_ellipse_pmts) {
    for ( const float &ang : pmt.angles) {
      fhangboltel->Fill (ang);
    }
    for ( const float &dang : pmt.dangs) {
      fhdangboltel->Fill (dang);
    }
  }

  // look for duplicate bolts and keep only best matches
  prune_bolts( f_ellipse_pmts, fhdangboltel->GetMean() );
  // remove pmts below threshold (9 bolts)
  prune_pmts( f_ellipse_pmts, 12, "f_ellipsehough" );
  


  //draw ellipses with line showing shortest distance from the point to the ellipse.
  //draw_ellipses(f_ellipses, image_ellipse );
  for ( const PMTIdentified& her : f_ellipse_pmts ){
    if ( her.bolts.size() < 9 ) continue;
    Size axes(  int(her.circ.get_a()), int(her.circ.get_b()) );
    Point center( int(her.circ.get_xy().x), int(her.circ.get_xy().y) );
    ellipse( image_ellipse, center, axes, RADTODEG( her.circ.get_phi() ), 0., 360,  Scalar (255, 102, 255), 2 );
    
    Scalar my_color( 0, 0, 255 );
    for ( const Vec3f& xyz : her.bolts ){
      circle( image_ellipse, Point( xyz[0], xyz[1] ), 3, my_color, 1, 0 );
      //image_ellipse.at<Scalar>( xy.x, xy.y ) = my_color;
    }
  }  
	  
  // annotate with bolt numbers and angles
  overlay_bolt_angle_boltid( f_ellipse_pmts, image_ellipse );

  if (write_image){
    string outputname = build_output_filename (infname, "fast_ellipses");
    imwrite(outputname, image_ellipse);
  }

  //Fill ellipse_dist histogram
  TH1D * f_ellipse_dist = new TH1D ("f_ellipse_dist",
				    "Distance from bolt to fast PMT ellipse; distance (pixels); Count/bin",
				    51, -0.5, 49.5);
	  
  for (const PMTIdentified & pmt : f_ellipse_pmts) {
    for (const float dist : pmt.dists) {
      f_ellipse_dist->Fill (dist);
    }
  }
       
  //trialend
}


void pmtidentified_histograms( const std::vector< PMTIdentified > & ellipse_pmts, const std::string & tag  ){

  // histogram PMT ellipse parameters
  TH1D * hpmt_locx = new TH1D( (tag+"hpmt_locx").c_str(),"PMT location ; x (pixels); counts/bin",120,0.,4000.);
  TH1D * hpmt_locy = new TH1D( (tag+"hpmt_locy").c_str(),"PMT location ; y (pixels); counts/bin",90,0.,3000.);
  TH1D * hpmt_a    = new TH1D( (tag+"hpmt_a").c_str(),"PMT a ; a (pixels); counts/bin",100,0,300);
  TH1D * hpmt_b    = new TH1D( (tag+"hpmt_b").c_str(),"PMT b ; b (pixels); counts/bin",100,0,300);
  TH1D * hpmt_area = new TH1D( (tag+"hpmt_area").c_str(),"PMT area ; area (px^2); counts/bin",200,0,200000);
  TH1D * hpmt_phi  = new TH1D( (tag+"hpmt_phi").c_str(),"PMT phi ; phi (radians); counts/bin",100, 0, 2*pi );
  TH1D * hpmt_e    = new TH1D( (tag+"hpmt_e").c_str(),"PMT e ; eccentricity; counts/bin",100, 0, 1.0 );
  TH1D * hpmt_pkval= new TH1D( (tag+"hpmt_pkval").c_str(),"PMT peakval ; peakval; counts/bin",256, -0.5, 255.5 );

  // histogram PMT phi as function of locations
  TH2D * hpmt_phixy = new TH2D( (tag+"hpmt_phixy").c_str(),"PMT ellipse angle ; x (pixels); y (pixels)",100,0.,4000.,75,0.,3000.);
  TH2D * hpmt_axy = new TH2D( (tag+"hpmt_axy").c_str(),"PMT ellipse a ; x (pixels); y (pixels)",100,0.,4000.,75,0.,3000.);
  TH2D * hpmt_bxy = new TH2D( (tag+"hpmt_bxy").c_str(),"PMT ellipse b ; x (pixels); y (pixels)",100,0.,4000.,75,0.,3000.);
  TH2D * hpmt_areaxy = new TH2D( (tag+"hpmt_areaxy").c_str(),"PMT ellipse area ; x (px^2); y (pixels)",100,0.,4000.,75,0.,3000.);
  TH2D * hpmt_exy = new TH2D( (tag+"hpmt_exy").c_str(),"PMT ellipse e ; x (pixels); y (pixels)",100,0.,4000.,75,0.,3000.);
  TH2D * hpmt_houghpk = new TH2D( (tag+"hpmt_houghpk").c_str(),"PMT hough peak ; x (pixels); y (pixels)",100,0.,4000.,75,0.,3000.);
  


  // histogram of area as function of radius from center of image
  TH2D * havsr = new TH2D( (tag+"hpmt_avsr").c_str(),"PMT area vs distance from center of image; r (pixels); area (px^2)", 200, 0., 4000., 200, 0., 1000000. );


  


  for ( const PMTIdentified& her : ellipse_pmts ){
    float x = her.circ.get_xy().x ;
    float y = her.circ.get_xy().y ;
    float aa = her.circ.get_a();
    float bb = her.circ.get_b();
    float ee = her.circ.get_e();
    float area = aa*bb*pi;
    float phi = her.circ.get_phi();
    hpmt_pkval->Fill( her.peakval );
    hpmt_houghpk->Fill( x, y, her.peakval );
    hpmt_locx->Fill( x );
    hpmt_locy->Fill( y );
    hpmt_a->Fill( aa );
    hpmt_b->Fill( bb );
    hpmt_area->Fill( area );
    hpmt_phi->Fill( phi );
    hpmt_e->Fill( ee );
    
    hpmt_phixy->Fill( x, y, phi );
    hpmt_axy->Fill( x, y, aa );
    hpmt_bxy->Fill( x, y, bb );
    hpmt_areaxy->Fill( x, y, area );
    hpmt_exy->Fill( x, y, ee );
      
    // distance from center of image?
    float r = std::sqrt(  (x-2000)*(x-2000) + (y-1375)*(y-1375) );
    havsr->Fill( r, area );
  }
}




void slow_ellipse_detection( const std::vector< cv::Vec3f > blobs, Mat& image_houghellipse, 
			     bool write_image, const std::string& infname, const MedianTextData & mtd ){ 

  bool do_ellipse_hough = (bool)config::Get_int( "do_ellipse_hough" );
  if ( do_ellipse_hough ){

    unsigned nbins_bb = (unsigned)config::Get_int("ellipse_hough_nbins_bb");
    unsigned nbins_ee = (unsigned)config::Get_int("ellipse_hough_nbins_ee");
    unsigned nbins_phiphi = (unsigned)config::Get_int("ellipse_hough_nbins_phiphi");
    unsigned nbins_x = (unsigned)config::Get_int("ellipse_hough_nbins_x");
    unsigned nbins_y = (unsigned)config::Get_int("ellipse_hough_nbins_y");
    //# above number of bins multiply as short (2 bytes per bin)
    //# therefore eg. 2 x 40 x 10 x 10 x 2300 x 1300 = 23.8 GB !!!
    float bbmin = (float)config::Get_double("ellipse_hough_bbmin");
    float bbmax = (float)config::Get_double("ellipse_hough_bbmax");
    float eemin = (float)config::Get_double("ellipse_hough_eemin");
    float eemax = (float)config::Get_double("ellipse_hough_eemax");
    float phiphimin = (float)config::Get_double("ellipse_hough_phphimin");
    float phiphimax = (float)config::Get_double("ellipse_hough_phiphimax");
    float xmin = (float)config::Get_double("ellipse_hough_xmin");
    float xmax = (float)config::Get_double("ellipse_hough_xmax");
    float ymin = (float)config::Get_double("ellipse_hough_ymin");
    float ymax = (float)config::Get_double("ellipse_hough_ymax");

    ///===========================================================
    /// Begin ellipse hough transfrom stuff
    EllipseHough h ( nbins_bb, bbmin, bbmax, 
		     nbins_ee, eemin, eemax,
		     nbins_phiphi, phiphimin, phiphimax,
		     nbins_x, xmin, xmax,
		     nbins_y, ymin, ymax );

    h.set_minhits( config::Get_int("ellipse_hough_minhits") );
    h.set_threshold( config::Get_int("ellipse_hough_threshold") );
    h.set_minhits( config::Get_double("ellipse_hough_drscale") );

    std::vector< xypoint > data;
    for ( unsigned i=0 ; i < blobs.size(); ++i ){
      //float radius = blobs[i][2];
      int blobx = blobs[i][0];
      int bloby = blobs[i][1];
      data.push_back( xypoint( blobx , bloby ) );
    }
    HoughEllipseResults hers = h.find_ellipses( data );	  
	  
    /// take hough resutls and fill vector of PMTIdentified info
    std::vector< PMTIdentified > ellipse_pmts;
    for ( const HoughEllipseResult& her : hers ){
      ellipse_st pmtloc{ her.e };
      std::vector< Vec3f > boltlocs;
      std::vector< float > dists;
      for ( const xypoint& xy : her.data ){
	//	if(her.e.dmin( xy )<5){
	boltlocs.push_back( Vec3f( xy.x, xy.y, 3 ) );
	dists.push_back( her.e.dmin( xy ) );
	//	}
      }
      ellipse_pmts.push_back( PMTIdentified( pmtloc, boltlocs, dists, her.peakval ) );
    }

    TH1D * hangboltel = new TH1D ("hangboltel", "Angles of bolts (hough ellipse); angle (deg)", 360, 0., 360.);
    TH1D * hdangboltel = new TH1D ("hdangboltel", "Angle of bolt from expected (hough ellipse); #Delta angle (deg)", 60, -15., 15.);

    for  (const PMTIdentified & pmt : ellipse_pmts) {
      for ( const float &ang : pmt.angles) {
	hangboltel->Fill (ang);
      }
      for ( const float &dang : pmt.dangs) {
	hdangboltel->Fill (dang);
      }
    }

    std::cout<<"========================== Before Pruning PMTS ===================================="<<std::endl;
    for  (const PMTIdentified & pmt : ellipse_pmts) {
      print_pmt_ellipse( std::cout, pmt );
      //std::cout<<pmt<<std::endl;
    }
    if ( write_image ){
      Mat image_before = image_houghellipse.clone();

      for ( const PMTIdentified& her : ellipse_pmts ){
	//if ( her.bolts.size() < 9 ) continue;
	Size axes(  int(her.circ.get_a()), int(her.circ.get_b()) );
	Point center( int(her.circ.get_xy().x), int(her.circ.get_xy().y) );
	ellipse( image_before, center, axes, RADTODEG( her.circ.get_phi() ), 0., 360,  Scalar (255, 102, 255), 2 );
	

	Scalar my_color( 0, 0, 255 );
	for ( const Vec3f& xyz : her.bolts ){
	  circle( image_before, Point( xyz[0], xyz[1] ), 3, my_color, 1, 0 );
	  //image_ellipse.at<Scalar>( xy.x, xy.y ) = my_color;
	}
      }

      // annotate with bolt numbers and angles
      overlay_bolt_angle_boltid( ellipse_pmts, image_before );	  

      string outputname = build_output_filename ( infname , "houghellipse_before");
      std::cout<<"Writing image "<<outputname<<std::endl;
      imwrite (outputname, image_before );
    }



    // histogram PMT locations before pruning	  
    pmtidentified_histograms( ellipse_pmts, "preprune" );


    //   prune_bolts_improved2( ellipse_pmts, hdangboltel->GetMean() );
    prune_bolts_super_improved( ellipse_pmts, hdangboltel->GetMean() );
    //prune_bolts( ellipse_pmts, hdangboltel->GetMean() );
    // look for duplicate bolts and keep only best matches
    //prune_bolts( ellipse_pmts, hdangboltel->GetMean() );
    // remove pmts below threshold (9 bolts)
    //    prune_pmts( ellipse_pmts, 9, "ellipsehough" );
    //prune_pmts( ellipse_pmts, 9, "ellipsehough" );
    prune_pmts_improved( ellipse_pmts, 10, "ellipsehough" );

    std::cout<<"========================== AFTER Pruning PMTS ===================================="<<std::endl;
    for  (const PMTIdentified & pmt : ellipse_pmts) {
      print_pmt_ellipse( std::cout, pmt );
      //std::cout<<pmt<<std::endl;
    }

    //making histogram of angle of ellipse wrt y-axis and size
    //TH2D *ellipse_feature = new TH2D("Ellipse feature", "Ellipse feature; y-value;phi(degree);size; count); 


    //Fill ellipse_dist histogram
    TH1D * ellipse_dist = new TH1D ("ellipse_dist",
				    "Distance from bolt to PMT ellipse; distance (pixels); Count/bin",
				    51, -0.5, 49.5);
	  
    for (const PMTIdentified & pmt : ellipse_pmts) {
      for (const float dist : pmt.dists) {
	if(dist<0){continue;}
	ellipse_dist->Fill (dist);
      }
    }
	

    /// draw all ellipses in her on image_ellipse and write
    ///for ( const HoughEllipseResult& her : hers ){
    for ( const PMTIdentified& her : ellipse_pmts ){
      //if ( her.bolts.size() < 9 ) continue;
      Size axes(  int(her.circ.get_a()), int(her.circ.get_b()) );
      Point center( int(her.circ.get_xy().x), int(her.circ.get_xy().y) );
      ellipse( image_houghellipse, center, axes, RADTODEG( her.circ.get_phi() ), 0., 360,  Scalar (255, 102, 255), 2 );
      
      Scalar my_color( 0, 0, 255 );
      for ( const Vec3f& xyz : her.bolts ){
	circle( image_houghellipse, Point( xyz[0], xyz[1] ), 3, my_color, 2, 0 );
	//image_ellipse.at<Scalar>( xy.x, xy.y ) = my_color;
      }
    }  

    // histogram PMT locations	  
    pmtidentified_histograms( ellipse_pmts, "final" );
    
    // histograms after pruning
    TH1D * hangboltel_cor = new TH1D ("hangboltel_cor", "Angles of bolts (hough ellipse corrected); angle (degrees)", 360, 0., 360.);
    TH1D * hdangboltel_cor = new TH1D ("hdangboltel_cor", "Angle of bolt from expected (hough ellipse corrected); #Delta angle (degrees)", 60, -15., 15.);

    //template will be {angles of one pmt bolt,-5, angles of next pmt bolts}
    //it's a way to identify which angle belong to which bolt and which pmt.
    //since final_PMTs and seral_final_bolts are in sync this approach work.
    //const vector<float>& angles = get_angles( final_PMTs, serial_final_bolts );
    string outfilename = build_output_textfilename( infname, "he_bolts" ); 
    std::ofstream fang_out (outfilename);

    for  (const PMTIdentified & pmt : ellipse_pmts) {
      std::cout << pmt;
      fang_out << pmt;
      for ( const float &ang : pmt.angles) {
	float ang_cor = ang - hdangboltel->GetMean();
	if ( ang_cor < 0 ) ang_cor += 360.0;
	hangboltel_cor->Fill(ang_cor);
      }
      for ( const float &dang : pmt.dangs) {
	hdangboltel_cor->Fill (dang -  hdangboltel->GetMean());
      }
    }

    // add the circles for the bolts after pruning
    for (const PMTIdentified & pmt : ellipse_pmts) {
      draw_circle_from_data (pmt.bolts, image_houghellipse,
			     Scalar (255, 255, 255), 1);
    }


 if ( mtd.size() > 0){ //have truth
      //Find bolt matches between those we found and truth
      find_closest_matches( ellipse_pmts, mtd );
	    
      //Draw line from truth to closest bolt found 
      draw_line (ellipse_pmts, mtd, image_houghellipse);
	    
      //Distance histogram
      TH1D * ellipse_truth_dist = new TH1D ("ellipse_truth_dist",
						  "Distance from true bolt loc to ellipse found bolt; distance (pixels); count/bin",
						  51, -0.5, 49.5);
      for ( const PMTIdentified & pmt : ellipse_pmts ){
	for ( const float dist : pmt.dist_txt ) {
	  if ( dist < bad_dmin ){
	    ellipse_truth_dist->Fill( dist );
	  }
	}
      }
	    
      draw_text_circles (image_houghellipse, mtd);
    }

    // annotate with bolt numbers and angles
    overlay_bolt_angle_boltid( ellipse_pmts, image_houghellipse );	  

    if ( write_image ){
      string outputname = build_output_filename ( infname , "houghellipse");
      imwrite (outputname, image_houghellipse );
    }
  }
}
    

//
// circle_bolt_detection
//
// Inputs:
//   const Mat & img_can -- canvas on which to search for small bolt-like circles
//   bool write_images -- set to true to write images showing found bolts 
//   const string& infname -- input filename used to label output written files
// Output:
//   Mat& image_color -- (non-empty) color image input that circles of found features are drawn onto
//   Mat& img_circles -- (empty) image to draw the found bolts onto 
// Return value:
//   vector< Vec3f > -- vector of found bolt locations
//
vector< Vec3f > circle_bolt_detection( const Mat & img_can, Mat& image_color, Mat& img_circles, bool write_images, const std::string& infname ){

  /// Hough Transform
  vector < Vec3f > circles;
  int dp = config::Get_int ("hough_dp");	// Inverse ratio of the accumulator resolution to the image resolution. For example, if dp=1 , the accumulator has the same resolution as the input image. If dp=2 , the accumulator has half as big width and height. 
  int minDist = config::Get_int ("hough_minDist");	// min distance between circles
  int param1 = config::Get_int ("hough_param1");	// threshold placed on image
  int param2 = config::Get_int ("hough_param2");	// minimum accumulator value to call it a circle
  int minR = config::Get_int ("hough_minR");	//= 3 # minimum radius in pixels
  int maxR = config::Get_int ("hough_maxR");	// = 10 # maximum radius in pixels

  HoughCircles (img_can, circles, HOUGH_GRADIENT, dp, minDist, param1,
		param2, minR, maxR);

  draw_circle_from_data (circles, image_color, Scalar (0, 0, 255));

  if ( write_images ) {
    string outputname = build_output_filename (infname, "hough");
    imwrite (outputname, image_color);
  }


  /// Make image that just has circle centers from previous Hough Transform on it
  img_circles = Mat::zeros (img_can.size (), img_can.type ());
  //draw_found_center (circles, img_circles);
  draw_foundblobs( circles, img_circles );

  if ( write_images ) {
    string outputname = build_output_filename ( infname, "houghCandidate");
    imwrite (outputname, img_circles);
  }
  
  return circles;
}


///
/// pmt_circle_detection
///
/// Look for cirlces of bolts indicating the location of the PMT.  Input is an image with bolt features only as white, and other pixels black.
/// Uses the built in HoughCircle of OpenCV, then apply constraints to reduce the number of falsely found bolts using bolt angle being
/// separted by ~15 degrees, and falesly found PMTs by eliminating overlapping circles.  
///
/// Inputs:
///  const vector< Vec3f >& blobs -- vector of locations of candidate bolts
///  const Mat& image -- image with only the points of interest on it that will be points on the circles found
///  bool write_images -- set to true to save the image produced to a file, otherwise it is not written to disk
///  const string & infname -- input file name used to label the output files (images written to disk, and bolts to text file)
///  const MedianTextData & mtd -- information about true locations of points of interest and which circle they go with (can be empty vector if not available)
///  const string & label -- label to add to histograms produced 
///
/// Output:
///   Mat& image_color -- draw circles found, and color the points used on each circle
///      
void pmt_circle_detection( const std::vector< Vec3f >& blobs, const Mat& image, Mat& image_color, bool write_images, const std::string & infname, const MedianTextData & mtd, const std::string& label  ){

  int dp = config::Get_int ("sec_hough_dp");	//if dp=1 , the accum has resolution of input image. If dp=2 , the accumulator has half as big width and height. 
  int minDist = config::Get_int ("sec_hough_minDist");	// min distance between circles
  int param1 = config::Get_int ("sec_hough_param1");	// threshold placed on image
  int param2 = config::Get_int ("sec_hough_param2");	// minimum accumulator value to call it a circle
  int minR = config::Get_int ("sec_hough_minR");	//= 3 # minimum radius in pixels
  int maxR = config::Get_int ("sec_hough_maxR");	// = 10 # maximum radius in pixels


  /// Look for circles of bolts
  vector < Vec3f > circles_of_bolts;

  HoughCircles (image, circles_of_bolts, HOUGH_GRADIENT, dp,
		minDist, param1, param2, minR, maxR);

  //Overlays detected circles from second hough transfrom
  //draw_circle_from_data ( circles_of_bolts, image_color,
  //			  Scalar (255, 102, 255));

  std::vector< PMTIdentified > final_pmts;
  find_candidate_bolts( blobs, circles_of_bolts, final_pmts, image );

  std::ostringstream osname, ostitle;
  osname << "hangbolt_" << label;
  ostitle << "Angles of bolts " << label << "; Angle clockwise from vertical; counts/bin";
  TH1D * hangbolt = new TH1D (osname.str().c_str(), ostitle.str().c_str(), 360, 0., 360.);

  std::ostringstream osname1, ostitle1;
  osname1 << "hdangbolt_" << label;
  ostitle1 << "Angle of bolt from expected " << label <<" #Delta angle; counts/bin";
  TH1D * hdangbolt = new TH1D (osname1.str().c_str(), ostitle1.str().c_str(), 60, -15., 15.);

  for  (const PMTIdentified & pmt : final_pmts) {
    for ( const float &ang : pmt.angles) {
      hangbolt->Fill (ang);
    }
    for ( const float &dang : pmt.dangs) {
      hdangbolt->Fill (dang);
    }
  }

  // look for duplicate bolts and keep only best matches
  prune_bolts_improved( final_pmts, hdangbolt->GetMean() );
  //prune_bolts( final_pmts, hdangbolt->GetMean() );
  // remove pmts below threshold (12 bolts)
  prune_pmts( final_pmts, 12, label );
	
  std::ostringstream osname2, ostitle2;
  osname2 << "bolt_dist_" << label;
  ostitle2 << "Distance closest bolt to PMT circle " << label <<"; distance (pixels); Counts/bin";

  //Fill blob_dist histogram
  TH1D * blob_dist = new TH1D ( osname2.str().c_str(),
				ostitle2.str().c_str(),
				51, -0.5, 49.5);

  for (const PMTIdentified & pmt : final_pmts) {
    for (const float dist : pmt.dists) {
      blob_dist->Fill (dist);
    }
  }

  // Make image of final answer
  // add the circles for the PMTs selected
  for (const PMTIdentified & pmt : final_pmts) {
    vector < Vec3f > circs;
    circs.push_back( Vec3f( pmt.circ[0], pmt.circ[1], pmt.circ[2] ) );
    draw_circle_from_data (circs, image_color,
			   Scalar (255, 102, 255), 2);
    draw_circle_from_data (pmt.bolts, image_color,
			   Scalar (0, 0, 255), 3);
  }

  if (mtd.size()>0) {
    // draw the true location of the bolts on the final image
    draw_text_circles (image_color, mtd);

    //Find bolt matches between those we found and truth
    find_closest_matches( final_pmts, mtd );

    //Draw line from truth to closest bolt found 
    draw_line (final_pmts, mtd, image_color);

    //Distance to truth histogram
    std::ostringstream osname3, ostitle3;
    osname3 << "bolt_dist_true_" << label;
    ostitle3 << "Distance from true bolt loc to found bolt loc " << label <<"; distance (pixels); Counts/bin";


    TH1D * blob_truth_dist = new TH1D ( osname3.str().c_str(),
					ostitle3.str().c_str(),
					51, -0.5, 49.5);
    for ( const PMTIdentified & pmt : final_pmts ){
      for ( const float dist : pmt.dist_txt ) {
	if ( dist < bad_dmin ){
	  blob_truth_dist->Fill( dist );
	}
      }
    }
  }
	
  // histograms after pruning
  std::ostringstream osname4, ostitle4;
  osname4 << "hangbolt_cor_" << label;
  ostitle4 << "Angles of bolts (corrected) " << label << "; Angle Clockwise from vertical; counts/bin";
  TH1D * hangbolt_cor = new TH1D (osname4.str().c_str(), ostitle4.str().c_str(), 360, 0., 360.);

  std::ostringstream osname5, ostitle5;
  osname5 << "hdangbolt_cor_" << label;
  ostitle5 << "Angle of bolt from expected (corrected) " << label << "; #Delta angle (degrees); counts/bin";
  TH1D * hdangbolt_cor = new TH1D (osname5.str().c_str(), ostitle5.str().c_str(), 60, -15., 15.);


  //template will be {angles of one pmt bolt,-5, angles of next pmt bolts}
  //it's a way to identify which angle belong to which bolt and which pmt.
  //since final_PMTs and seral_final_bolts are in sync this approach work.
  //const vector<float>& angles = get_angles( final_PMTs, serial_final_bolts );
  string outfilename = build_output_textfilename (infname, label);
  std::ofstream fang_out (outfilename);

  for  (const PMTIdentified & pmt : final_pmts) {
    std::cout << pmt;
    fang_out << pmt;
    for ( const float &ang : pmt.angles) {
      float ang_cor = ang - hdangbolt->GetMean();
      if ( ang_cor < 0 ) ang_cor += 360.0;
      hangbolt_cor->Fill(ang_cor);
    }
    for ( const float &dang : pmt.dangs) {
      hdangbolt_cor->Fill (dang -  hdangbolt->GetMean());
    }
  }


  // add the circles for the bolts after pruning
  for (const PMTIdentified & pmt : final_pmts) {
    draw_circle_from_data (pmt.bolts, image_color,
			   Scalar (255, 255, 255), 1);
  }

  if ( write_images ) {
    //overlays the bolt angle and bolt id in input image using final_pmts input
    overlay_bolt_angle_boltid(final_pmts, image_color );
    string outputname = build_output_filename (infname, label );
    imwrite (outputname, image_color);
  }

}


void histogram_blobs_bymtd( const vector< KeyPoint>& keypoints, const MedianTextData & mtd ){

  /* class KeyPoint
    Data structure for salient point detectors.
    Point2f pt
        coordinates of the keypoint
    float size
        diameter of the meaningful keypoint neighborhood
    float angle
        computed orientation of the keypoint
	(-1 if not applicable). Its possible
	values are in a range [0,360)
	degrees. It is measured relative to
	image coordinate system (y-axis is
	directed downward), ie in clockwise.
    float response
        the response by which the most strong keypoints have been selected. Can be used for further sorting or subsampling
    int octave
	octave (pyramid layer) from which the keypoint has been extracted
  */

  TH2D* hblobsizeang = new TH2D("hblobsizeang", "Blob size vs angle; size; angle; counts/bin", 40, -0.5, 39.5, 40, -20.0, 20.0);//0., 360. );
  TH1D* hblobresponse = new TH1D("hblobresponse", "Blob response; response; counts/bin", 200, 0., 400. );
  TH1D* hbloboctave = new TH1D("hbloboctave", "Blob octave; octave; counts/bin", 40, -20., 20. );


  TH2D* hblobsizeangbolt = new TH2D("hblobsizeangbolt", "Blob size vs angle of bolt; size; angle; counts/bin", 40, -0.5, 39.5, 40, -20.,20);//0., 360. );
  TH1D* hblobresponsebolt = new TH1D("hblobresponsebolt", "Blob response of bolt; response; counts/bin", 200, 0., 400. );
  TH1D* hbloboctavebolt = new TH1D("hbloboctavebolt", "Blob octave; octave of bolt; counts/bin", 40, -20., 20. );

 
    
  for( KeyPoint k: keypoints){
    float  mind = 1000000;
    for (const MedianTextRecord & rec : mtd ){
      //for ( const cv::Vec3f & circ : circles ){
      //float dist = std::sqrt( (k.pt.x - rec.x())*(k.pt.x - rec.x()) +
      //                      (k.pt.y - rec.y())*(k.pt.y - rec.y()) );
      float dist = RobustLength( fabs(k.pt.x - rec.x()),fabs(k.pt.y - rec.y()) );
      if ( dist < mind ) {
	bool reverse = true;
	for(unsigned j=0; j<keypoints.size(); ++j){
	  KeyPoint kp = keypoints[j];
	  //float d1 = std::sqrt((kp.pt.x-rec.x())*(kp.pt.x-rec.x())+
	  //                   (kp.pt.y-rec.y())*(kp.pt.y-rec.y()));
          float d1 = RobustLength( fabs(kp.pt.x-rec.x()),fabs(kp.pt.y-rec.y()) );
	  if(d1<dist){ reverse = false; break;}
        }

        if(reverse){ mind = dist; }
	
      }
    }
    /*
  for( KeyPoint k: keypoints){
    float mind = 10000000;
    for(const MedianTextRecord & r: mtd){
      double d = RobustLength(r.x()-k.pt.x, r.y()-k.pt.y);
      if(d<mind){ mind=d; }
    } 
    */
    if(mind<3){
      hblobsizeang->Fill(k.size,k.angle,1.0);
      hblobresponse->Fill(k.response);
      hbloboctave->Fill(k.octave);     
    }
    //    else if(mind==1000000){;}
    else{
      hblobsizeangbolt->Fill(k.size,k.angle,1.0);
      hblobresponsebolt->Fill(k.response);
      hbloboctavebolt->Fill(k.octave);
    }
  }
}


//Removing the noise inside of PMT(Filtering blobs)
void rem_bolts (const vector< Vec3f >& blobs1, vector< Vec3f >& blobs, const cv::Mat img){
  for(int j=0; j<blobs1.size(); ++j){
    Vec3f bl = blobs1[j];
    int x = bl[0];
    int y = bl[1];
    int ra = bl[2];
    
    Vec3f intensity = img.at<Vec3b>(y,x);
    int r = intensity.val[2];
    int g = intensity.val[1];
    int b = intensity.val[0];
    bool near = false;
    int n=0;
    float ang;
    int ind;
    for(int i=0; i<blobs1.size(); ++i){
      if(i==j){continue;}
      float dist = RobustLength(fabs(blobs1[i][0]-x),fabs(blobs1[i][1]-y));
      if (dist<35){//25){//35){
	n++;   //count number of blobs within 35 px
	ind = i;
	ang = atan2f((x-blobs1[i][0]),-(y-blobs1[i][1])); //getting angle with ^ axis wrt image  
	ang = RADTODEG( ang );
	ang = (ang<0)?(ang+360):ang; //getting angle between 0-360
      }
      //if n=1 there is another blob within 300px that forms line with current bolt(bolt in which neighbourhood we are looking at), bolt corresponding to n=1, and current bolt then increase count. 
      if(0 ){//n==1 && dist <300 && i!=ind ){
	float theta = atan2f((x-blobs1[i][0]),-(y-blobs1[i][1]));; //getting angle with ^ axis wrt image  
	theta = RADTODEG(theta);
	theta = (theta<0)?(theta+360):theta; //getting angle between 0-360
	
	float diff = fabs(theta-ang);
	if(fabs(diff-180)<0.5 || diff<0.5 ){n++;
	}
      }
      if(n>=2){
	near = true; break;
      }
    }

    //skip if color in centre is yellowish, or radius is greater than 15 or there are two blobs near current blob or there is line. 
    if((abs(r-g)<70 && g>50 && abs(g-b)>100)||ra>10 || near)continue;
    
    blobs.push_back(bl);
  }

  //now overlaying all bolts and kept bolt in image
  Mat blbs = img.clone();
  draw_foundblobs( blobs1, blbs );
  for ( const Vec3f& xyz : blobs ){
    circle( blbs, Point( xyz[0], xyz[1] ), 3, Scalar(0,0,255), 2, 0 );

  }  

  imwrite("trimmmmm.jpg", blbs);
}



int main (int argc, char **argv) {

      if (argc != 2 && argc != 3) {
	printf("usage: FindBoltLocations <Input_image_with_path> [<median-bolt-loc-filename>]\n");
	return -1;
      }
      //have_truth == 0 means 2 argument mode without true bolt information provided
      //bool have_truth = argc - 2;
      
      Mat image_color = imread (argv[1], IMREAD_COLOR);	//IMREAD_GRAYSCALE,

      Mat image_line = image_color.clone(); //This image will have blobs, mtd and lines
      //trial to find circles of pmt
      /*      Mat tri;
      cvtColor (image_color, tri, COLOR_RGBA2GRAY);
      //     Mat magent =
      cmyk[1].at<uchar>(i, j) = (1 - g - k) / (1 - k) * 255.;
      Mat img_gaus;
      //GaussianBlur( image, img_blur, Size(blurpixels, blurpixels), blursigma)
      GaussianBlur( tri, img_gaus, Size(5, 5), 6);
      
      imwrite("gausblur.jpg",img_gaus);

      Mat image_can = img_gaus.clone();
      Canny (img_gaus, image_can, 200, 230);
      imwrite("canny.jpg",image_can);
      */
      if (!image_color.data) {
	printf ("No image data \n");
	return -1;
      }
      
      //option has final, text, candidate, circled, filters
      const vector < bool > & option = setup_image_saveflags ();
      
      Mat image_final = image_color.clone ();
      Mat image_ellipse = image_color.clone();
      Mat image_houghellipse = image_color.clone();
      Mat img_blob_map = image_color.clone();
      
      /// build output image
      Mat image;
      cvtColor (image_color, image, COLOR_RGBA2GRAY);
      Mat image1 = output_image_by_color( image_color, false, 3, argv[1], true ); // K (intensity)
      Mat image2 = equalize_by_color( image_color, argv[1], true );

      // Open a root file to put histograms into
      TFile * fout = new TFile ("FindBoltLocation.root", "RECREATE");


      try
	{
	  bool verbose = config::Get_int("verbosity");
      
	  // Gaussian blur
	  Mat img_blur = apply_gaussian_blur( image, option[4], argv[1] );

	  // Bilateral filter
	  Mat img_flt = apply_bilateral_filter( img_blur, option[4], argv[1] );

	  // Sobel edge detection
	  Mat grad = apply_sobel_edge( img_flt, option[4], argv[1] );

	  //Canny edge detector.
	  Mat img_can = apply_canny_edge( grad, option[ option.size()-1 ], argv[1] );

	  // equalize image
	  Mat image_clahe = apply_clahe( img_can, (bool)config::Get_int("clahe_write"), argv[1] );
    
	  // threshold cut
	  apply_image_threshold (image_clahe, config::Get_int ("threshold") );

	  // do blob detection
	  Mat blob_circles;
	  std::vector < KeyPoint > keypoints;
	  vector< Vec3f > blobs1 = blob_detect( image_clahe, blob_circles, option[3], argv[1], keypoints );

	  vector < Vec3f > blobs;
	  //Filtering closely packed blobs and yellowish blobs for corner image;
	  bool corner = config::Get_int("corner");
	  if(corner){
	  rem_bolts(blobs1, blobs, image_color);
	  }
	  else{
	    blobs = blobs1;
	  }
	  // Read in truth bolt locations if provided
	  // Returns empty vector if text file is not supplied.
	  const MedianTextData & mtd = assign_data_from_text(argc, string (argv[argc - 1]));
	  //drawing mdt in image
	  draw_text_circles (image_line, mtd);
	  draw_circle_from_data (blobs, image_line,
				 Scalar (0, 0, 255), 1);
	  TH1D *dis = new TH1D("distance of blob from truth","Distance of blob from truth;distance(px);counts/bin",1000,-0.5,999.5);
	  make_bolt_dist_histogram_wrt_txt( blobs, mtd, dis, image_line );
	  string outputname = build_output_filename (argv[1], "mapped");
	  imwrite(outputname, image_line);
	  
	  histogram_blobs_bymtd( keypoints, mtd );

	  //Debug information PMt
	  if(verbose){
	    for (const MedianTextRecord & rec:mtd) {
	      std::cout << rec;
	    }
	  }

	  // fast ellipse detection
	  fast_ellipse_detection( blobs, image_ellipse, true, argv[1], mtd );

	  // slow ellipse detection (hough_ellipse.hpp)
	  slow_ellipse_detection( blobs, image_houghellipse, true, argv[1], mtd ); 
	
	  // circle detection on canny image to find bolts
	  //Mat img_circles; // final image of locatad bolts
	  //vector< Vec3f > circles = circle_bolt_detection( img_can, image_color, img_circles, option[3], argv[1] );

	  // PMT circle detection on bolts found by circle_bolt_detection
	  //pmt_circle_detection( circles, img_circles, img_circles, option[1], argv[1], mtd, "houghbolts" );

	  // PMT circle detection on bolts found by blob detection
	  pmt_circle_detection( blobs, blob_circles, img_blob_map, option[1], argv[1], mtd, "houghblobs" );



	/*
  ///bwlabel trial
	///bwlabel trial
	Mat contour1 = image.clone();
	Mat contour3 = image.clone();
    
	cv::threshold( contour1, contour1, 230, 255,THRESH_BINARY );
	BwLabel b;
	vector<vector<int>> ma = b.find_label(contour1);
	struct dots{
	  std::vector<cv::Point> pts;
	  cv::Point centre;
	  dots(std::vector<cv::Point> pts): pts(pts){}
	};
	std::vector<dots> maa;
	std::vector<int> grps;
	for(unsigned i=0; i<ma.size();i++){
	  //bool brk =false;
	  for(unsigned j=0; j<ma[0].size();j++){
	    if(ma[i][j]!=0){ 
	      bool bel=false;
	      for(int grou : grps){
		if(ma[i][j] == grou){bel=true; break;}
	      }
	      if(!bel){
		
		grps.push_back(ma[i][j]);
	      }
	    }
	      
	  }
	}
	std::cout<<"group size "<<grps.size()<<std::endl;
    
	for(int grn:grps){
	  int ida=-1;
	  int idb=-1;
	  for(int i=0; i<ma.size();i++){
	    bool brk =false;
	    for(int j=0; j<ma[0].size();j++){
	      if(ma[i][j] == grn){ ida = i; idb =j; brk=true; break;} 
	    }
	    if(brk){break;}
	  }
	  
	  int lowy = ((ida-150)>0)?(ida-150):0;
	  int lowx = ((idb-150)>0)?(idb-150):0;
	  int higy = ((ida+150)<ma.size())?(ida+150):ma.size();
	  int higx = ((idb+150)<ma[0].size())?(idb+150):ma[0].size();
	  
	  std::vector<cv::Point> pts;
	  for(int k=lowy; k<higy;k++){
	    for(int l=lowx; l<higx;l++){
	      //for(int k=0; k<ma.size();k++){
	      //for(int l=0; l<ma[k].size();l++){
	      if(ma[k][l]!=0){
		if(ma[k][l]==grn){ pts.push_back(cv::Point(l,k));} 
	      }
	    }
	  }
	  std::cout<<"Pts size "<<pts.size()<<std::endl;
	  if(pts.size()>50 && pts.size()<200){maa.push_back(dots(pts));}
	  // }
	}

	std::vector<cv::Point> ab;
	for(int i=0; i<maa.size();i++){
	  int xmin=maa[i].pts[0].x;
	  int ymin=maa[i].pts[0].y;
	  int xmax=maa[i].pts[0].x;
	  int ymax=maa[i].pts[0].y;
      
	  for(cv::Point xy:maa[i].pts){
	    if(xy.x<xmin){xmin=xy.x;}
	    if(xy.y<ymin){ymin=xy.y;}
	    if(xy.x>xmax){xmax=xy.x;}
	    if(xy.y>ymax){ymax=xy.y;}
	  }
	  maa[i].centre=cv::Point((xmin+xmax)/2.0,(ymin+ymax)/2.0);
	  ab.push_back(cv::Point((xmin+xmax)/2.0,(ymin+ymax)/2.0));
	}
    
	
	for(cv::Point p:ab){
	  int x1 =p.x;
	  int y1 = p.y;
	  for(cv::Point q:ab){
	    int x2 = q.x;
	    int y2 = q.y;
	    
	    float slope = (y2-y1+0.0)/(x2-x1);
	    float dist = sqrt(std::pow(y2-y2,2)+std::pow(x2-x1,2));
	    
	  }
	  
	}
	Mat img3 = Mat::zeros(contour1.rows, contour1.cols, CV_8UC1);

	std::cout<<"size of label "<<maa.size()<<" cols "<<ma[0].size()<<std::endl;
	for(int i=0; i<maa.size(); i++){
	  //for(int j=0; j<maa[i].pts.size(); j++){
	  //if(ma[i][j]>0){
	  //img3.at<uchar>((cv::Point)maa[i].pts[j]) = 255;//30*ma[i][j];
	  img3.at<uchar>(maa[i].centre) = 255;//30*ma[i][j];
	  std::cout<<"x "<<maa[i].centre.x<<" y "<<maa[i].centre.y<<std::endl;
	  //img3.at<Vec3b>((cv::Point)maa[i].pts[j])[1]= 255;//10*ma[i][j];
	  //img3.at<Vec3b>((cv::Point)maa[i].pts[j])[2]= 255;//20*ma[i][j];
	  //}
	  //cout:: 
	}
    
	vector<Vec2f> lines;
	//HoughLines(img3, lines, 1, CV_PI/180, 60, 0, 0,CV_PI/2.2,CV_PI/1.8 );
	HoughLines(img3, lines, 1, CV_PI/180, 30, 0, 0,0,CV_PI );
	//vector<Vec4i> linesP;
	// Draw the lines
	for( size_t i = 0; i < lines.size(); i++ )
	  {
	    float rho = lines[i][0], theta = lines[i][1];
	    Point pt1, pt2;
	    double a = cos(theta), b = sin(theta);
	    double x0 = a*rho, y0 = b*rho;
	    pt1.x = cvRound(x0 + 1000*(-b));
	    pt1.y = cvRound(y0 + 1000*(a));
	    pt2.x = cvRound(x0 - 1000*(-b));
	    pt2.y = cvRound(y0 - 1000*(a));
	    line( img3, pt1, pt2, Scalar(255,255,255), 1, LINE_AA);
	  }
     
	imwrite("threshold.jpg", contour1);
	imwrite("bwlbl.jpg",img3);
	//trialend
	*/
    
    }
    catch (std::string e)
    {
	std::cout << "Error with config file key " << e << std::endl;
    }

    fout->Write ();
    fout->Close ();

    return 0;
}
