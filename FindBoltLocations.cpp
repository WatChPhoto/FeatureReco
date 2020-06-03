#include <stdio.h>
#include <opencv2/opencv.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/imgcodecs.hpp>
#include <string>
#include <vector>

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>

#include "Configuration.hpp"
#include "MedianTextReader.hpp"

using std::string;
using std::vector;
 
using namespace cv;


string build_output_filename( const string& in, const string& tag ){
    string outputname;
    size_t idx = in.find_last_of("/");
    outputname = string( tag ) + in.substr(idx+1 );
    return outputname;
}



void histogram_channel0( const Mat& img,  std::string hname="hist" ){

  TH1D* hout = new TH1D(hname.c_str(),"Intensity ; Intensity; Count/bin",256, -0.5, 255.5);

  int nRows = img.rows;
  int nCols = img.cols;

  for ( int x = 0; x<nCols; ++x){
    for ( int y = 0; y<nRows; ++y){
      Scalar intensity = img.at<uchar>(y, x);
      hout->Fill( intensity.val[0] );
    }
  }


}

void apply_image_threshold( Mat& img, int threshold=250 ){
  
 
  int nRows = img.rows;
  int nCols = img.cols;
  
  for ( int x = 0; x<nCols; ++x){
    for ( int y = 0; y<nRows; ++y){
      Scalar intensity = img.at<uchar>(y, x);
      if ( intensity.val[0] < threshold ) {
	img.at<uchar>(y, x) = 0;
      }
    }
  }
}



/// Simple structure to hold indices into two vectors of matches with distance between those two indices stored
struct IndexMatchDist {
  unsigned idx_txt;   //index of medianTextReader
  unsigned idx_circ;  //index of circles
  float    dist;
  IndexMatchDist( unsigned itxt=0, unsigned icirc=0, float fdist=-1. ) : idx_txt( itxt ), idx_circ( icirc ), dist( fdist ) { }
};


/// find_closest_matches finds closest match to each entry in mtd to the closest entry in circles
/// without using any circle twice.
vector< IndexMatchDist > find_closest_matches( const vector<Vec3f>& circles, const MedianTextData & mtd ){
  vector< IndexMatchDist > data;
  vector< IndexMatchDist > data121;
  int i=0;

  for (const MedianTextRecord & rec : mtd ){
    float  mindist = 10000;
    int index = 0, j = 0;
    for ( const Vec3f & circ : circles ){
      float dist = std::sqrt( (circ[0] - rec.x())*(circ[0] - rec.x()) +
			      (circ[1] - rec.y())*(circ[1] - rec.y()) );
      if ( dist < mindist ) {mindist = dist; index =j;}  //I am assuming that we will find point closer than d=10000
      ++j;
    }
    IndexMatchDist temp;
    temp.idx_txt= i;   //index of medianTextReader
    temp.idx_circ= index;  //index of circles
    temp.dist = mindist;  //mindist
    data.push_back(temp);
    ++i;
  }

  int k=0;
  for ( IndexMatchDist & rec : data ){
    int mtd_index = rec.idx_txt;
    int index = rec.idx_circ;
    float dist = rec.dist;
    int i=0;
    int new_index=k;
    float min_dist = dist;
    for( IndexMatchDist & ab : data){
      int mtd_index1 = ab.idx_txt;
      int index1 = ab.idx_circ;
      float dist1 = ab.dist;
      
      if( mtd_index!=mtd_index1 && index==index1){   //if point from text file is different but from found circle is same.
	if(dist1<min_dist){
	  min_dist = dist1;
	  new_index = i;
	}
	
      } 
      ++i;
    }

    bool copy = false;
    for( IndexMatchDist& d : data121){
      if( d.idx_txt == data[ new_index ].idx_txt ){
	copy=true;
      }
    }
    if(!copy){
      data121.push_back( data[new_index] );
    }
    ++k;
  }

  return data121;
}




void make_bolt_dist_histogram( const vector< IndexMatchDist > & matches){
  TH1D* hout = new TH1D("bolt_distance","Distance to closest bolt ; distance (pixels); Count/bin",1001, 0.5, 500.5);
  for ( const IndexMatchDist& m : matches ){
    hout->Fill( m.dist );
  }
}


/// Calculate a metric to identify a bolt
/// Ratio of average pixel intensity inside identified bolt radius over
/// pixel intensity outside bolt to 2*radius
/// Returns -1 if no pixels match this criteria.
double calculate_bolt_metric( const Vec3f& circ, const Mat& img ) {

  static int callcount = 0;

  TH1D* hinbolt=nullptr; 
  TH1D* houtbolt=nullptr;
  if ( callcount == 0 ){
    std::cout<<"Make hinbolt and houtbolt histograms"<<std::endl;
    hinbolt = new TH1D("hinbolt","hinbolt",256, -0.5,255.5 ); 
    houtbolt = new TH1D("houtbolt","houtbolt",256, -0.5,255.5 ); 
    ++callcount;
  }
 

  float circ_x = circ[0];
  float circ_y = circ[1];
  float circ_r = circ[2];

  unsigned ninside = 0;
  unsigned noutside =0;
  double avg_inside =0.;
  double avg_outside=0.;

  int nRows = img.rows;
  int nCols = img.cols;
  
  int n = 2; //n factor by which second radius is bigger. 
  int x_low = ((circ_x-n*circ_r-1)>0)?(circ_x-n*circ_r-1):0;
  int x_high = ((circ_x + n*circ_r+1)<nCols)?(circ_x + n*circ_r+1):nCols;
  int y_low = ((circ_y - (n*circ_r)-1)>0)?(circ_y - (n*circ_r)-1):0;
  int y_high = ((circ_y + n*circ_r + 1)<nRows)?(circ_y + n*circ_r + 1):nRows;  


  for ( int x = x_low; x<x_high; ++x){
    for ( int y = y_low; y<y_high; ++y){
      Scalar intensity = img.at<uchar>(y, x);
      double cur_radius = sqrt( (x-circ_x)*(x-circ_x) + (y-circ_y)*(y-circ_y) );
      if ( cur_radius < circ_r ){
	++ninside;
	avg_inside += intensity.val[0]; 
	if (hinbolt) hinbolt->Fill( intensity.val[0] );
      } else if ( cur_radius < n*circ_r) {
	++noutside;
	avg_outside += intensity.val[0];
	if (houtbolt) houtbolt->Fill( intensity.val[0] );

      }
    }
  }

  if ( ninside >0 && noutside >0 ){
    //    return avg_outside/(pow(circ_r+3,2)-pow(circ_r,2));
    //return (avg_inside/ninside)  / (avg_outside/noutside) ;
    return (avg_inside/ninside) / (avg_outside/noutside) ;
  }

    
  return -1.;

}


//edited
void make_bolt_dist_histogram_wrt_txt( const vector<Vec3f>& circles, const MedianTextData& mtd, Mat &img ){
  TH1D* hout1 = new TH1D("bolt_distance_wrt_text","Distance to closest bolt ; distance (pixels); Count/bin",501, -0.5, 500.5);
hout1->SetAxisRange(0,501,"X");

 
  for (const MedianTextRecord & rec : mtd ){
    float  mindist = 10000;
    unsigned x,y;
    for ( const Vec3f & circ : circles ){
      float dist = std::sqrt( (circ[0] - rec.x())*(circ[0] - rec.x()) +
			      (circ[1] - rec.y())*(circ[1] - rec.y()) );
      if ( dist < mindist ) {mindist = dist; x=circ[0]; y = circ[1]; }  //I am assuming that we will find point closer than d=10000
    
    }
    //arrowedLine(Mat& img, Point pt1, Point pt2, const Scalar& color, int thickness=1, int line_type=8, int shift=0, double tipLength=0.1)
    //    if(mindist <100){ arrowedLine(img,Point(rec.x(),rec.y()), Point(x,y),  (0,0,0)); }
    hout1->Fill( mindist );
  }
}




void make_bolt_metric_histograms( const vector<Vec3f>& circles, const MedianTextData& mtd, const vector< IndexMatchDist >& data121,  Mat &imbw,  Mat &imcol){
  
  TH1D* metric_all  = new TH1D("Metric_all" , "Bolt metric for all circles ;Inside to outside Intensity ratio; Count",502, -1.5, 255.5);
  TH1D* metric_good = new TH1D("Metric_good", "Bolt metric for matched circles ;Inside to outside Intensity ratio; Count",502, -1.5, 255.5);
  TH1D* metric_bad  = new TH1D("Metric_bad" , "Bolt metric for un-matched circles  ;Inside to outside Intensity ratio; Count",502, -1.5, 255.5);
 
  TH2D* metric_2d   = new TH2D("Metric_2d", " Bolt metric vs distance to bolt;  Distance to bolt (pixels); Bolt Metric", 100, -0.5, 99.5, 101, -1.5, 255.5 );

  for ( const IndexMatchDist & rec : data121 ){ 
    int index = rec.idx_circ;
    float mindist = rec.dist;
    double metric_val = calculate_bolt_metric( circles[index], imbw ); 
    metric_all->Fill( metric_val );
    if (mindist < 4) { 
      metric_good->Fill( metric_val );
    } else {
      metric_bad->Fill( metric_val );
    }
    metric_2d->Fill( mindist,  metric_val );

    int in = rec.idx_txt;
    int m_x= mtd[in].x();
    int m_y= mtd[in].y();
    int x=circles[index][0];
    int y=circles[index][1];
    // if(mindist <100)
      { arrowedLine(imcol, Point(m_x,m_y), Point(x,y),  (0,0,0)); }
  
  }
  //imwrite("new.jpg",imbw);
  //Edit end



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
    TFile * fout = new TFile("FindBoltLocation.root","recreate");


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

    //Drawing Michel's point in the picture.
    for ( const MedianTextRecord & rec : mtd ){
      Point center_michel(rec.x(), rec.y());   
      //used radius of 10 pixels and green color.
      circle( image_color, center_michel, 10, Scalar(0,255,0), 1, 8, 0 );
    }
    outputname = build_output_filename( argv[1], "michel" );
    imwrite( outputname, image_color );
    //end of drawing michel's points.
    

    } catch ( std::string e ){
      std::cout<<"Error with config file key "<<e<<std::endl;
    }
    
    

    fout->Write();
    fout->Close();

    return 0;
}
