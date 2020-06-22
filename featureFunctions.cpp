#include "featureFunctions.hpp"
#include <iostream>
#include <opencv2/imgproc.hpp>
#include "Configuration.hpp"

std::string build_output_filename( const std::string& in, const std::string& tag ){
  std::string outputname;
  size_t idx = in.find_last_of("/");
  outputname = std::string( tag ) + in.substr(idx+1 );
  return outputname;
}


std::string build_output_textfilename( const std::string& in, const std::string& tag ){
  std::string outputname;
  size_t idx1 = in.find_last_of("/");
  size_t idx2 = in.find_last_of(".");
  outputname = std::string( tag ) + in.substr(idx1+1, idx2-idx1-1 ) + ".txt";
  return outputname;
}

void histogram_channel0( const cv::Mat& img,  std::string hname ){

  TH1D* hout = new TH1D(hname.c_str(),"Intensity ; Intensity; Count/bin",256, -0.5, 255.5);

  int nRows = img.rows;
  int nCols = img.cols;

  for ( int x = 0; x<nCols; ++x){
    for ( int y = 0; y<nRows; ++y){
      cv::Scalar intensity = img.at<uchar>(y, x);
      hout->Fill( intensity.val[0] );
    }
  }
}

void apply_image_threshold( cv::Mat& img, int threshold ){ 
  int nRows = img.rows;
  int nCols = img.cols;
  
  for ( int x = 0; x<nCols; ++x){
    for ( int y = 0; y<nRows; ++y){
      cv::Scalar intensity = img.at<uchar>(y, x);
      if ( intensity.val[0] < threshold ) {
	img.at<uchar>(y, x) = 0;
      }
    }
  }
}


/// Calculate a metric to identify a bolt
/// Ratio of average pixel intensity inside identified bolt radius over
/// pixel intensity outside bolt to 2*radius
/// Returns -1 if no pixels match this criteria.
double calculate_bolt_metric( const cv::Vec3f& circ, const cv::Mat& img ) {
  /*=====================================================================================================
  static int callcount = 0;
  
  TH1D* hinbolt=nullptr; 
  TH1D* houtbolt=nullptr;
  if ( callcount == 0 ){
    std::cout<<"Make hinbolt and houtbolt histograms"<<std::endl;
    hinbolt = new TH1D("hinbolt","hinbolt",256, -0.5,255.5 ); 
    houtbolt = new TH1D("houtbolt","houtbolt",256, -0.5,255.5 ); 
    ++callcount;
  }
  ======================================================================================================*/
  float circ_x = circ[0];
  float circ_y = circ[1];
  float circ_r = circ[2];

  unsigned ninside = 0.;
  unsigned noutside= 0.;
  double avg_inside= 0.;
  double avg_outside=0.;
  double sqr_inside= 0.;
  double sqr_outside=0.;

  int nRows = img.rows;
  int nCols = img.cols;
  
  int n = 2; //n factor by which second radius is bigger. 
  int x_low = ((circ_x-n*circ_r-1)>0)?(circ_x-n*circ_r-1):0;
  int x_high = ((circ_x + n*circ_r+1)<nCols)?(circ_x + n*circ_r+1):nCols;
  int y_low = ((circ_y - (n*circ_r)-1)>0)?(circ_y - (n*circ_r)-1):0;
  int y_high = ((circ_y + n*circ_r + 1)<nRows)?(circ_y + n*circ_r + 1):nRows;  


  for ( int x = x_low; x<x_high; ++x){
    for ( int y = y_low; y<y_high; ++y){
      cv::Scalar intensity = img.at<uchar>(y, x);
      double cur_radius = sqrt( (x-circ_x)*(x-circ_x) + (y-circ_y)*(y-circ_y) );
      if ( cur_radius < circ_r ){
	++ninside;
	avg_inside += intensity.val[0];
	sqr_inside += std::pow(intensity.val[0],2);
	//	if (hinbolt) hinbolt->Fill( intensity.val[0] );
      } else if ( cur_radius < n*circ_r) {
	++noutside;
	avg_outside += intensity.val[0];
	sqr_outside += std::pow(intensity.val[0],2);
	//	if (houtbolt) houtbolt->Fill( intensity.val[0] );

      }
    }
  } 
 
  if ( ninside >0 && noutside >0 ){
    try{
    //    return avg_outside/(pow(circ_r+3,2)-pow(circ_r,2));
    //return (avg_inside/ninside)  / (avg_outside/noutside) ;
    if(config::Get_int("do_avg")){
    return (avg_inside/ninside) / (avg_outside/noutside) ;
    }

    if(config::Get_int("do_rms")){
      return (std::sqrt(sqr_inside)/ninside)/(std::sqrt(sqr_outside)/noutside);
    }  
}
 catch ( std::string e ){
    std::cout<<"Error with config file key "<<e<<std::endl;
  }  
  }
  
  return -1.;
}

void make_bolt_metric_histograms( const std::vector< PMTIdentified > &pmtsfound,  cv::Mat &imbw, TH1D *&metric_all, TH1D *&metric_good, TH1D *&metric_bad, TH2D *&metric_2d){
  /*  
  TH1D* metric_all  = new TH1D("Metric_all" , "Bolt metric for all circles ;Inside to outside Intensity ratio; Count",502, -1.5, 255.5);
  TH1D* metric_good = new TH1D("Metric_good", "Bolt metric for matched circles ;Inside to outside Intensity ratio; Count",502, -1.5, 255.5);
  TH1D* metric_bad  = new TH1D("Metric_bad" , "Bolt metric for un-matched circles  ;Inside to outside Intensity ratio; Count",502, -1.5, 255.5);
 
  TH2D* metric_2d   = new TH2D("Metric_2d", " Bolt metric vs distance to bolt;  Distance to bolt (pixels); Bolt Metric", 100, -0.5, 99.5, 101, -1.5, 255.5 );
  */
  for ( const PMTIdentified & pmtid : pmtsfound ){ 
    for ( unsigned boltidx = 0; boltidx < pmtid.bolts.size(); ++boltidx ){
      float mindist = pmtid.dist_txt[boltidx];      
      if ( mindist >= bad_dmin ) continue;
      double metric_val = calculate_bolt_metric( pmtid.bolts[ boltidx ], imbw ); 
      metric_all->Fill( metric_val );
      if (mindist < 4) { 
	metric_good->Fill( metric_val );
      } else {
	metric_bad->Fill( metric_val );
      }
      metric_2d->Fill( mindist,  metric_val );
    }
  }
}

//Draws line using given data121
//Data 121 maps text data to found circles one-to-one.
void draw_line( const std::vector< PMTIdentified >& pmtsfound, const MedianTextData& mtd, cv::Mat &imcol ){
  for ( const PMTIdentified & pmtid : pmtsfound ){
    for ( unsigned boltidx = 0; boltidx < pmtid.bolts.size(); ++boltidx ){
      if ( pmtid.dist_txt[ boltidx ] >= bad_dmin ) continue;
      int in = pmtid.idx_txt[ boltidx ];
      int m_x= mtd[in].x();
      int m_y= mtd[in].y();

      int x= pmtid.bolts[ boltidx ][0];
      int y= pmtid.bolts[ boltidx ][1];
    //line(Mat& img, Point pt1, Point pt2, const Scalar& color, int thickness=1, int lineType=8, int shift=0)
    line(imcol, cv::Point(m_x,m_y), cv::Point(x,y), cv::Scalar(0,0,0), 2, 8,0);
    }
  }
}

//Make histogram of metric of inbetween points that aren't mapped to the circles from the text file
/*
void histogram_inbetween(const std::vector<cv::Vec3f>& circles, const MedianTextData& mtd, const std::vector< IndexMatchDist >& data121, cv::Mat imbw, TH1D *&metric_inb){
  int x_min = 100000, x_max= -1;
  int y_min = 100000, y_max= -1;
  for(const MedianTextRecord val: mtd){
    int x = val.x();
    int y = val.y();
    if(x<x_min){x_min=x;}
    if(x>x_max){x_max=x;}
    if(y<y_min){y_min=y;}
    if(y>y_max){y_max=y;}
  }
   
  int i = 0;
  for(const cv::Vec3f circ: circles){
    int x = circ[0];
    int y = circ[1];
    if(x>=x_min && x <=x_max && y>=y_min && y<=y_max){
      bool belongs = false;
      for(const IndexMatchDist val: data121){
	int index = val.idx_circ;
	if(index==i){belongs=true; break;}
      }
      if(!belongs){
	double metric_val = calculate_bolt_metric( circ, imbw );
	metric_inb->Fill(metric_val);
      }
    }
    i++;
  }
}
*/

void draw_circle_from_data(const std::vector <cv::Vec3f> data, cv::Mat & image, cv::Scalar color, int line_width  ){
  for( size_t i = 0; i < data.size(); i++ ) {
    cv::Point center(cvRound(data[i][0]), cvRound(data[i][1]));
    int radius = cvRound(data[i][2]);
    // draw the circle center                                                                                                               
    //circle( image_color, center, 3, Scalar(0,0,255), 1, 8, 0 );                                                                           
    // draw the circle outline                                                                                                              
    cv::circle( image, center, radius, color, line_width, 8, 0 );
    std::cout<<"Circle "<<i<<" radius = "<<radius<<" at ( "<<data[i][0]<<", "<<data[i][1]<<" )"<<std::endl;
  }
  std::cout<<"=============================================================================="<<std::endl;
}

void draw_found_center(const std::vector<cv::Vec3f> data, cv::Mat & image){
  for( size_t i = 0; i < data.size(); i++ ) {
    image.at<uchar>(data[i][1], data[i][0]) = 255 ;
  }
}

void draw_text_circles(cv::Mat &img, const MedianTextData& mtd){
  for ( const MedianTextRecord & rec : mtd ){
    cv::Point center_text(rec.x(), rec.y());
    //used radius of 10 pixels and green color.                               
    circle( img, center_text, 10, cv::Scalar(0,255,0), 1, 8, 0 );
  }
}
