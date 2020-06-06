
#include "featureFunctions.hpp"
#include <iostream>                                                                                                                          
#include <TH1D.h>                                                                                                                            
#include <TH1D.h>                                                                                                                            
#include <TH2D.h>                                                                                                                            
#include <opencv2/imgproc.hpp> 

std::string build_output_filename( const std::string& in, const std::string& tag ){
  std::string outputname;
  size_t idx = in.find_last_of("/");
  outputname = std::string( tag ) + in.substr(idx+1 );
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

/// find_closest_matches finds closest match to each entry in mtd to the closest entry in circles
/// without using any circle twice.
std::vector< IndexMatchDist > find_closest_matches( const std::vector<cv::Vec3f>& circles, const MedianTextData & mtd ){
  std::vector< IndexMatchDist > data;
  std::vector< IndexMatchDist > data121;
  int i=0;

  for (const MedianTextRecord & rec : mtd ){
    float  mindist = 10000;
    int index = 0, j = 0;
    for ( const cv::Vec3f & circ : circles ){
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

void make_bolt_dist_histogram( const std::vector< IndexMatchDist > & matches){
  TH1D* hout = new TH1D("bolt_distance","Distance to closest bolt ; distance (pixels); Count/bin",1001, 0.5, 500.5);
  for ( const IndexMatchDist& m : matches ){
    hout->Fill( m.dist );
  }
}


/// Calculate a metric to identify a bolt
/// Ratio of average pixel intensity inside identified bolt radius over
/// pixel intensity outside bolt to 2*radius
/// Returns -1 if no pixels match this criteria.
double calculate_bolt_metric( const cv::Vec3f& circ, const cv::Mat& img ) {

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
      cv::Scalar intensity = img.at<uchar>(y, x);
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
void make_bolt_dist_histogram_wrt_txt( const std::vector<cv::Vec3f>& circles, const MedianTextData& mtd, cv::Mat &img ){
  TH1D* hout1 = new TH1D("bolt_distance_wrt_text","Distance to closest bolt ; distance (pixels); Count/bin",501, -0.5, 500.5);
  hout1->SetAxisRange(0,501,"X");

 
  for (const MedianTextRecord & rec : mtd ){
    float  mindist = 10000;
    unsigned x,y;
    for ( const cv::Vec3f & circ : circles ){
      float dist = std::sqrt( (circ[0] - rec.x())*(circ[0] - rec.x()) +
			      (circ[1] - rec.y())*(circ[1] - rec.y()) );
      if ( dist < mindist ) {mindist = dist; x=circ[0]; y = circ[1]; }  //I am assuming that we will find point closer than d=10000
    
    }
    //arrowedLine(Mat& img, Point pt1, Point pt2, const Scalar& color, int thickness=1, int line_type=8, int shift=0, double tipLength=0.1)
    //    if(mindist <100){ arrowedLine(img,Point(rec.x(),rec.y()), Point(x,y),  (0,0,0)); }
    hout1->Fill( mindist );
  }
}

void make_bolt_metric_histograms( const std::vector<cv::Vec3f>& circles, const MedianTextData& mtd, const std::vector< IndexMatchDist >& data121,  cv::Mat &imbw,  cv::Mat &imcol){
  
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
    { arrowedLine(imcol, cv::Point(m_x,m_y), cv::Point(x,y),  (0,0,0)); }
  
  }
  //imwrite("new.jpg",imbw);
  //Edit end
}

//Make histogram of metric of inbetween points that aren't mapped to the circles from the text file
void histogram_inbetween(const std::vector<cv::Vec3f>& circles, const MedianTextData& mtd, const std::vector< IndexMatchDist >& data121, cv::Mat imbw){
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
  
  TH1D* metric_inb  = new TH1D("Metric_inbetween" , "Bolt metric for inbetween circles ;Inside to outside Intensity ratio; Count",502, -1.5, 255.5); 
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
