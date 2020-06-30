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

#include "ellipse_detection.hpp" //ellipse detection fast
using std::string;
using std::vector;

using namespace cv;

int main (int argc, char **argv) {

  if (argc != 2 && argc != 3) {
    printf("usage: FindBoltLocations <Input_image_with_path> [<median-bolt-loc-filename>]\n");
    return -1;
  }
  //have_truth == 0 means 2 argument mode without true bolt information provided
  bool have_truth = argc - 2;
  
  Mat image_color = imread (argv[1], IMREAD_COLOR);	//IMREAD_GRAYSCALE,
  if (!image_color.data) {
    printf ("No image data \n");
    return -1;
  }

    //option has final, text, candidate, circled, filters
    const vector < bool > & option = setup_image_saveflags ();

    Mat image_final = image_color.clone ();
    Mat image_ellipse = image_color.clone();
    Mat image_houghellipse = image_color.clone();

    /// build output image
    Mat image_orig;
    cvtColor (image_color, image_orig, COLOR_RGBA2GRAY);


    // equalize image
    std::cout<<"Applying equalization"<<std::endl;
    Mat image;
    Ptr<CLAHE> clahe = createCLAHE();
    clahe->setClipLimit( 4.0 );
    //clahe->setTilesGridSize( 16 );
    clahe->apply( image_orig, image );
    std::cout<<"Equalized"<<std::endl;



    /*
    // sharpen image using "unsharp mask" algorithm
    Mat blurred; double sigma = 1, threshold = 5, amount = 1;
    GaussianBlur(image, blurred, Size(), sigma, sigma);
    Mat lowContrastMask = abs(image - blurred) < threshold;
    Mat sharpened = image*(1+amount) + blurred*(-amount);
    image.copyTo(sharpened, lowContrastMask);
    //

    image = sharpened.clone();
    imwrite("sharp.jpg", sharpened);
*/
    // Open a root file to put histograms into
    TFile * fout = new TFile ("FindBoltLocation.root", "RECREATE");

    string outputname;

    // Gaussian blur

    Mat img_blur = image.clone ();
    try
    {
	bool verbose = config::Get_int("verbosity");
	bool do_gaus_blur = (bool) config::Get_int ("do_gaus_blur");

	if ( do_gaus_blur ) {
	      int blurpixels = config::Get_int ("blurpixels");	// size of kernel in pixels (must be odd)
	      double blursigma = config::Get_double ("blursigma");	// sigma of gaussian in pixels

	      GaussianBlur( image, img_blur, Size(blurpixels, blurpixels), blursigma);
	      if ( option[4] ) {
		    outputname = build_output_filename (argv[1], "gausblur");
		    imwrite (outputname, img_blur);
		}
	  }
	// Bilateral filter
	bool do_bifilter = (bool) config::Get_int ("do_bifilter");
	Mat img_flt = img_blur.clone ();
	if ( do_bifilter ) {
	      int d = config::Get_int ("d");	// value 5-9 distance around each pixel to filter (must be odd)
	      int sigColor = config::Get_int ("sigColor");	// range of colours to call the same
	      int sigSpace = config::Get_int ("sigSpace");	// ???

	      bilateralFilter (image, img_flt, d, sigColor, sigSpace);
	      if (option[4]) {
		    outputname = build_output_filename (argv[1], "bifilter");
		    imwrite (outputname, img_flt);
	      }
	}
	/// Do Sobel edge detection
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

	      if (option[4]) {
		    outputname = build_output_filename (argv[1], "sobel");
		    imwrite (outputname, grad);
	      }
	}
	//Canny edge detector.                                                                                                                   
	bool do_canny = (bool) config::Get_int ("do_canny");
	int thresh_low = config::Get_int ("thresh_low");	//The gradient value below thresh_low will be discarded.                                  
	int thresh_high = config::Get_int ("thresh_high");	//The gradient value above thresh_high will be used. The inbetween gradient is kept if \the edge is connected.                                                                                                                        

	Mat img_can = grad.clone();
	if (do_canny) {
	      Canny (grad, img_can, thresh_low, thresh_high);

	      if (option[option.size() - 1]) {
		    outputname = build_output_filename (argv[1], "canny");
		    imwrite (outputname, img_can);
	      }
	}

	histogram_channel0 (img_can, "hbefore_thres_cut");
	int threshold = config::Get_int ("threshold");

	apply_image_threshold (img_can, threshold);
	histogram_channel0 (grad, "hafter_thres_cut");

	//Blob detection
	Mat img_blob = image.clone ();
	Mat img_blob_map = image_color.clone ();
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

	// Set up the detector with default parameters.                                                                                         
	Ptr < SimpleBlobDetector > detector = SimpleBlobDetector::create (params);

	// Detect blobs.                                                                                                                        
	std::vector < KeyPoint > keypoints;
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

	//ellipse trial
	detect_ellipse(blobs, image_ellipse, 90, 120 ,90, 120, 4);
	imwrite("ellipses.jpg",image_ellipse);
	//trialend

	//Draws circle from data to the input image
	draw_circle_from_data (blobs, img_blob_map, Scalar (0, 0, 255));
	if (option[3]) {
	  outputname = build_output_filename (argv[1], "blob");
	  imwrite (outputname, img_blob_map);
	}
	// Make image that just has circle centers from blob detection
	Mat blob_circles = Mat::zeros (image.size (), image.type ());
	draw_found_center (blobs, blob_circles);


	if (option[2]) {
	  outputname = build_output_filename (argv[1], "blobCandidate");
	  imwrite (outputname, blob_circles);
	}
	// Blobend         

	/// Read in truth bolt locations
	//Returns empty vector if text file is not supplied.
	const MedianTextData & mtd = assign_data_from_text(argc, string (argv[argc - 1]));

	//Debug information PMt
	if(verbose){
	  for (const MedianTextRecord & rec:mtd) {
	    std::cout << rec;
	  }
	}

	
	bool do_ellipse_hough = (bool)config::Get_int( "do_ellipse_hough" );
	if ( do_ellipse_hough ){

	  ///===========================================================
	  /// Begin ellipse hough transfrom stuff
	  EllipseHough h;
	  std::vector< xypoint > data;
	  for ( unsigned i=0 ; i < blobs.size(); ++i ){
	    data.push_back( xypoint( blobs[i][0], blobs[i][1] ) );
	  }
	  HoughEllipseResults hers = h.find_ellipses( data );
	  
	  std::cout<< hers <<std::endl;
	  
	  /// draw all ellipses in her on image_ellipse and write
	  for ( const HoughEllipseResult& her : hers ){
	    if ( her.data.size() < 9 ) continue;
	    Size axes(  int(her.e.get_a()), int(her.e.get_b()) );
	    Point center( int(her.e.get_xy().x), int(her.e.get_xy().y) );
	    ellipse( image_houghellipse, center, axes, RADTODEG( her.e.get_phi() ), 0., 360,  Scalar (0, 0, 255) );
	    
	    Scalar my_color( 0, 255, 0 );
	    for ( const xypoint& xy : her.data ){
	      circle( image_houghellipse, Point( xy.x, xy.y ), 3, my_color, 1, 0 );
	      //image_ellipse.at<Scalar>( xy.x, xy.y ) = my_color;
	    }
	  }  
	  
	  /// take hough resutls and fill vector of PMTIdentified info
	  std::vector< PMTIdentified > ellipse_pmts;
	  for ( const HoughEllipseResult& her : hers ){
	    //if ( her.data.size() < 9 ) continue;
	    Vec3f pmtloc{ her.e.get_xy().x, her.e.get_xy().y, her.e.get_b() };
	    std::vector< Vec3f > boltlocs;
	    std::vector< float > dists;
	    for ( const xypoint& xy : her.data ){
	      boltlocs.push_back( Vec3f( xy.x, xy.y, 3 ) );
	      dists.push_back( her.e.dmin( xy ) );
	    }
	    ellipse_pmts.push_back( PMTIdentified( pmtloc, boltlocs, dists ) );
	  }
	  
	  //Fill ellipse_dist histogram
	  TH1D * ellipse_dist = new TH1D ("ellipse_dist",
					  "Distance from bolt to PMT ellipse; distance (pixels); Count/bin",
					  51, -0.5, 49.5);
	  
	  for (const PMTIdentified & pmt : ellipse_pmts) {
	    for (const float dist : pmt.dists) {
	      ellipse_dist->Fill (dist);
	    }
	  }
	
	  if (have_truth){
	    //Find bolt matches between those we found and truth
	    find_closest_matches( ellipse_pmts, mtd );
	    
	    //Draw line from truth to closest bolt found 
	    draw_line (ellipse_pmts, mtd, image_ellipse);
	    
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


	  TH1D * hangboltel = new TH1D ("hangboltel", "Angles of bolts (hough ellipse); angle (deg)", 360, 0., 360.);
	  TH1D * hdangboltel = new TH1D ("hdangboltel", "Angle of bolt from expected (hough ellipse); #Delta angle (deg)", 60, -15., 15.);

	  for  (const PMTIdentified & pmt : ellipse_pmts) {
	    //std::cout << pmt;
	    //fang_out << pmt;
	    for ( const float &ang : pmt.angles) {
	      hangboltel->Fill (ang);
	    }
	    for ( const float &dang : pmt.dangs) {
	      hdangboltel->Fill (dang);
	    }
	  }


	  // look for duplicate bolts and keep only best matches
	  prune_bolts( ellipse_pmts, hdangboltel->GetMean() );
	  // remove pmts below threshold (9 bolts)
	  prune_pmts( ellipse_pmts, 9 );


	  // histogram PMT locations
	  TH1D * hpmt_locx = new TH1D("hpmt_locx","PMT location ; x (pixels); counts/bin",120,0.,4000.);
	  TH1D * hpmt_locy = new TH1D("hpmt_locy","PMT location ; y (pixels); counts/bin",90,0.,3000.);
	  TH1D * hpmt_b    = new TH1D("hpmt_b","PMT b ; b (pixels); counts/bin",80,80,140);
	  // histogram PMT phi as function of locations
	  TH2D * hpmt_phixy = new TH2D("hpmt_phixy","PMT ellipse angle ; x (pixels); y (pixels)",120,0.,4000.,90,0.,3000.);
	  TH2D * hpmt_bxy = new TH2D("hpmt_bxy","PMT ellipse b ; x (pixels); y (pixels)",120,0.,4000.,90,0.,3000.);
	  TH2D * hpmt_exy = new TH2D("hpmt_exy","PMT ellipse e ; x (pixels); y (pixels)",120,0.,4000.,90,0.,3000.);

	  for ( const HoughEllipseResult& her : hers ){
	    float x = her.e.get_xy().x ;
	    float y = her.e.get_xy().y ;
	    hpmt_locx->Fill( x );
	    hpmt_locy->Fill( y );
	    hpmt_b->Fill( her.e.get_b() );
	    hpmt_phixy->Fill( x, y, her.e.get_phi() );
	    hpmt_bxy->Fill( x, y, her.e.get_b() );
	    hpmt_exy->Fill( x, y, her.e.get_e() );
	  }

	  // histograms after pruning
	  TH1D * hangboltel_cor = new TH1D ("hangboltel_cor", "Angles of bolts (hough ellipse corrected); angle (degrees)", 360, 0., 360.);
	  TH1D * hdangboltel_cor = new TH1D ("hdangboltel_cor", "Angle of bolt from expected (hough ellipse corrected); #Delta angle (degrees)", 60, -15., 15.);


	  //template will be {angles of one pmt bolt,-5, angles of next pmt bolts}
	  //it's a way to identify which angle belong to which bolt and which pmt.
	  //since final_PMTs and seral_final_bolts are in sync this approach work.
	  //const vector<float>& angles = get_angles( final_PMTs, serial_final_bolts );
	  string outfilename = build_output_textfilename (argv[1], "he_bolts"); 
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

	  // annotate with bolt numbers and angles
	  overlay_bolt_angle_boltid( ellipse_pmts, image_houghellipse );
	  
	  outputname = build_output_filename (argv[1], "houghellipse");
	  imwrite (outputname, image_houghellipse );
	  


	  /// End ellipse hough transform stuff
	
	}
	

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

	if (option[3]) {
	      outputname = build_output_filename (argv[1], "hough");
	      imwrite (outputname, image_color);
	}


	/// Make image that just has circle centers from previous Hough Transform on it
	Mat img_circles = Mat::zeros (image.size (), image.type ());
	draw_found_center (circles, img_circles);

	if (option[2])
	  {
	      outputname = build_output_filename (argv[1], "houghCandidate");
	      imwrite (outputname, img_circles);
	  }
	/// Look for circles of bolts
	vector < Vec3f > hough_circles_of_bolts;
	//Look for circles of bolts from blob
	vector < Vec3f > blob_circles_of_bolts;
	dp = config::Get_int ("sec_hough_dp");	//if dp=1 , the accum has resolution of input image. If dp=2 , the accumulator has half as big width and height. 
	minDist = config::Get_int ("sec_hough_minDist");	// min distance between circles
	param1 = config::Get_int ("sec_hough_param1");	// threshold placed on image
	param2 = config::Get_int ("sec_hough_param2");	// minimum accumulator value to call it a circle
	minR = config::Get_int ("sec_hough_minR");	//= 3 # minimum radius in pixels
	maxR = config::Get_int ("sec_hough_maxR");	// = 10 # maximum radius in pixels


	HoughCircles (blob_circles, blob_circles_of_bolts, HOUGH_GRADIENT, dp,
		      minDist, param1, param2, minR, maxR);



	HoughCircles (img_circles, hough_circles_of_bolts, HOUGH_GRADIENT, dp,
		      minDist, param1, param2, minR, maxR);

	//Overlays detected circles from second hough transfrom
	draw_circle_from_data (hough_circles_of_bolts, image_color,
			       Scalar (255, 102, 255));
	draw_circle_from_data (blob_circles_of_bolts, img_blob_map,
			       Scalar (255, 102, 255));

	//Draw text file circles
	//Only if there are three arguments
	if (have_truth)
	  {
	      draw_text_circles (image_color, mtd);
	      draw_text_circles (img_blob_map, mtd);
	  }
	//save images if enabled
	if (option[1])
	  {
	      outputname = build_output_filename (argv[1], "hough_text");
	      imwrite (outputname, image_color);
	      outputname = build_output_filename (argv[1], "blob_text");
	      imwrite (outputname, img_blob_map);
	  }


	std::vector< PMTIdentified > final_pmts;
	find_candidate_bolts( blobs, blob_circles_of_bolts, final_pmts, image );

	//Fill bolb_dist2 histogram
	TH1D * blob_dist = new TH1D ("blob_dist",
				     "Distance to closest blob to PMT circle; distance (pixels); Count/bin",
				     51, -0.5, 49.5);

	for (const PMTIdentified & pmt : final_pmts) {
	    for (const float dist : pmt.dists) {
	         blob_dist->Fill (dist);
	    }
	}

	// Make image of final answer
	// add the circles for the PMTs selected
	for (const PMTIdentified & pmt : final_pmts) {
	      vector < Vec3f > circs{ pmt.circ };
	      draw_circle_from_data (circs, image_final,
				     Scalar (255, 102, 255), 2);
	      draw_circle_from_data (pmt.bolts, image_final,
				     Scalar (0, 0, 255), 3);
	}

	if (have_truth) {
	  // draw the true location of the bolts on the final image
	  draw_text_circles (image_final, mtd);

	  //Find bolt matches between those we found and truth
	  find_closest_matches( final_pmts, mtd );

	  //Draw line from truth to closest bolt found 
	  draw_line (final_pmts, mtd, image_final);

	  //Distance histogram
	  TH1D * blob_truth_dist = new TH1D ("blob_truth_dist",
					 "Distance from true bolt loc to blob; distance (pixels); count/bin",
					 51, -0.5, 49.5);
	  for ( const PMTIdentified & pmt : final_pmts ){
	    for ( const float dist : pmt.dist_txt ) {
	      if ( dist < bad_dmin ){
		blob_truth_dist->Fill( dist );
	      }
	    }
	  }
	}
	


	TH1D * hangbolt = new TH1D ("hangbolt", "Angles of bolts", 360, 0., 360.);
	TH1D * hdangbolt = new TH1D ("hdangbolt", "Angle of bolt from expected", 60, -15., 15.);

	for  (const PMTIdentified & pmt : final_pmts) {
	  //std::cout << pmt;
	  //fang_out << pmt;
	  for ( const float &ang : pmt.angles) {
	    hangbolt->Fill (ang);
	  }
	  for ( const float &dang : pmt.dangs) {
	    hdangbolt->Fill (dang);
	  }
	}

	//fang_out.close ();

	// look for duplicate bolts and keep only best matches
	prune_bolts( final_pmts, hdangbolt->GetMean() );


	// histograms after pruning
	TH1D * hangbolt_cor = new TH1D ("hangbolt_cor", "Angles of bolts (corrected); angle (degrees)", 360, 0., 360.);
	TH1D * hdangbolt_cor = new TH1D ("hdangbolt_cor", "Angle of bolt from expected (corrected); #Delta angle (degrees)", 60, -15., 15.);

	//template will be {angles of one pmt bolt,-5, angles of next pmt bolts}
	//it's a way to identify which angle belong to which bolt and which pmt.
	//since final_PMTs and seral_final_bolts are in sync this approach work.
	//const vector<float>& angles = get_angles( final_PMTs, serial_final_bolts );
	string outfilename = build_output_textfilename (argv[1], "bolts");

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
	  draw_circle_from_data (pmt.bolts, image_final,
				 Scalar (255, 255, 255), 1);
	}


	//testing purpose
	if (option[0]) {
	  //overlays the bolt angle and bolt id in input image using final_pmts input
	  overlay_bolt_angle_boltid(final_pmts, image_final);
	  outputname = build_output_filename (argv[1], "final");
	  imwrite (outputname, image_final);
	}


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
