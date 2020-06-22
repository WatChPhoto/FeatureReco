#include <iostream>
#include <fstream>
#include <opencv2/imgproc.hpp>
#include <opencv2/imgcodecs.hpp>
#include <TFile.h>
#include "Configuration.hpp"
#include "featureFunctions.hpp"
#include <opencv2/features2d.hpp>	//Blob

#include "PMTIdentified.hpp"

#include<cmath>

using std::string;
using std::vector;

using namespace cv;

//Returns the data from text file
//Returns empty vector if text file is not supplied.
MedianTextData assign (int argc, string argv) {
    if (argc == 2)
      {
	  MedianTextData a;
	  return a;
      }
    if (argc == 3)
      {
	  MedianTextReader *
	      boltreader = MedianTextReader::Get ();

	  boltreader->set_input_file (string (argv));
	  return boltreader->get_data ();
      }
}

//flags to turn on/off saving images
vector < bool > setup_image_saveflags () {
    vector < bool > options;
    try
    {
	int option = config::Get_int ("save_option");

	for (int i = 0; i < 5; i++)
	  {
	      options.push_back (bool (option % 2));
	      option /= 10;
	  }

    } catch (std::string e)
    {
	std::cout << "Error with config file key " << e << std::endl;
    }
    return options;
}

int
main (int argc, char **argv)
{

    if (argc != 2 && argc != 3)
      {
	  printf
	      ("usage: FindBoltLocations <Input_image_with_path> [<median-bolt-loc-filename>]\n");
	  return -1;
      }
    //have_truth == 0 means 2 argument mode without true bolt information provided
    bool have_truth = argc - 2;

    Mat image_color = imread (argv[1], IMREAD_COLOR);	//IMREAD_GRAYSCALE,
    if (!image_color.data)
      {
	  printf ("No image data \n");
	  return -1;
      }
    //option has final, text, candidate, circled, filters
    const vector < bool > & option = setup_image_saveflags ();

    Mat image_final = image_color.clone ();

    /// build output image
    Mat image;
    cvtColor (image_color, image, COLOR_RGBA2GRAY);

    // Open a root file to put histograms into
    TFile * fout = new TFile ("FindBoltLocation.root", "RECREATE");

    string outputname;

    // Gaussian blur
    Mat img_blur = image.clone ();
    try
    {
	//bool debug = config::Get_int("debug");
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
	/// Read in bolt locations
	//Returns empty vector if text file is not supplied.
	const MedianTextData & mtd = assign (argc, string (argv[argc - 1]));

	//Debug information PMt
	//if(debug){
	for (const MedianTextRecord & rec:mtd) {
	      std::cout << rec;
	}
	//}


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
	
	//template will be {angles of one pmt bolt,-5, angles of next pmt bolts}
	//it's a way to identify which angle belong to which bolt and which pmt.
	//since final_PMTs and seral_final_bolts are in sync this approach work.
	//const vector<float>& angles = get_angles( final_PMTs, serial_final_bolts );
	string outfilename = build_output_textfilename (argv[1], "bolts");

	std::ofstream fang_out (outfilename);

	TH1D * hangbolt = new TH1D ("hangbolt", "Angles of bolts", 360, 0., 360.);
	TH1D * hdangbolt = new TH1D ("hdangbolt", "Angle of bolt from expected", 60, -15., 15.);

	for  (const PMTIdentified & pmt : final_pmts) {
	  std::cout << pmt;
	  fang_out << pmt;
	  for ( const float &ang : pmt.angles) {
	    hangbolt->Fill (ang);
	  }
	  for ( const float &dang : pmt.dangs) {
	    hdangbolt->Fill (dang);
	  }
	}

	fang_out.close ();

	//testing purpose
	if (option[0]) {
	  outputname = build_output_filename (argv[1], "final");
	  imwrite (outputname, image_final);
	}

    
    }
    catch (std::string e)
    {
	std::cout << "Error with config file key " << e << std::endl;
    }

    fout->Write ();
    fout->Close ();

    return 0;
}
