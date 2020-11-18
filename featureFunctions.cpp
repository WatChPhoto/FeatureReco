#include "featureFunctions.hpp"
#include <iostream>
#include <opencv2/imgproc.hpp>
#include "Configuration.hpp"


float RobustLength( float v1, float v2){
  double length;
  if(v1>v2){
    length = fabs(v1)*sqrt(1.0+std::pow((v2/v1),2));
  }

  else if(v2>v1){
    length = fabs(v2)*sqrt(1.0+std::pow((v1/v2),2));
  }
  else{
    length = fabs(v2)*sqrt(2.0);
  }

  return length;
}


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

  histogram_channel0 (img, "hbefore_thres_cut");

  for ( int x = 0; x<nCols; ++x){
    for ( int y = 0; y<nRows; ++y){
      cv::Scalar intensity = img.at<uchar>(y, x);
      if ( intensity.val[0] < threshold ) {
	img.at<uchar>(y, x) = 0;
      }
    }
  }

  histogram_channel0 (img, "hafter_thres_cut");

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
      //if we need to hide some bad ones but it doesn't seem like we do.
      //      std::cout<<"dist_txt = "<<pmtid.dist_txt[boltidx]<<std::endl;
      //  if ( pmtid.dist_txt[ boltidx ] >= bad_dmin ) continue;
     
      int in = pmtid.idx_txt[ boltidx ];
      int m_x= mtd[in].x();
      int m_y= mtd[in].y();

      int x= pmtid.bolts[ boltidx ][0];
      int y= pmtid.bolts[ boltidx ][1];
      
      float dist = std::sqrt( (m_x-x)*(m_x-x) + (m_y-y)*(m_y-y) );
      // std::cout<<"dist = "<<dist<<std::endl;
      // if ( dist > 20 ) continue;

      //line(Mat& img, Point pt1, Point pt2, const Scalar& color, int thickness=1, int lineType=8, int shift=0)
      line(imcol, cv::Point(m_x,m_y), cv::Point(x,y), cv::Scalar(0,0,0), 2, 8,0);
    }
  }
}

void make_bolt_dist_histogram_wrt_txt( const std::vector<cv::Vec3f>& circles, const MedianTextData& mtd, TH1D *&hist_dist, cv::Mat& imcol ){
  //TH1D* hout1 = new TH1D("bolt_distance_wrt_text","Distance to closest bolt ; distance (pixels); Count/bin",501, -0.5, 500.5);
  //hout1->SetAxisRange(0,501,"X");
  //goal is if a is closest to b and b is closest to a then they are the map.
  for (const MedianTextRecord & rec : mtd ){
    float  mindist = 1000000;
    int m_x = rec.x();
    int m_y = rec.y();

    int x,y;
    for ( const cv::Vec3f & circ : circles ){
      //float dist = std::sqrt( (circ[0] - rec.x())*(circ[0] - rec.x()) +
      //		      (circ[1] - rec.y())*(circ[1] - rec.y()) );
      float dist = RobustLength( fabs(circ[0] - rec.x()), fabs(circ[1] - rec.y()) );
      if ( dist < mindist ) {
	bool reverse = true;
	for(unsigned j=0; j<mtd.size(); ++j){
	  MedianTextRecord m = mtd[j];
	  // float d1 = std::sqrt((circ[0]-m.x())*(circ[0]-m.x())+
	 //		       (circ[1]-m.y())*(circ[1]-m.y()));
	  float d1 = RobustLength( fabs(circ[0]-m.x()), fabs(circ[1]-m.y()) );
	  if(d1<dist){ reverse = false; break;}
	}
	
	if(reverse){ mindist = dist; x=circ[0]; y = circ[1]; }  
	
      }
    }

    if( mindist!=1000000){
      line(imcol, cv::Point(rec.x(),rec.y()), cv::Point(x,y), cv::Scalar(0,0,0), 2, 8,0);
    
    hist_dist->Fill( mindist );
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

void draw_circle_from_data(const std::vector <cv::Vec3f>& data, cv::Mat & image, cv::Scalar color, int line_width  ){
  for( size_t i = 0; i < data.size(); i++ ) {
    cv::Point center(cvRound(data[i][0]), cvRound(data[i][1]));
    int radius = cvRound(data[i][2]);
    // draw the circle center                                                                                                               
    //circle( image_color, center, 3, Scalar(0,0,255), 1, 8, 0 );                                                                           
    // draw the circle outline                                                                                                              
    cv::circle( image, center, radius, color, line_width, 8, 0 );
    //std::cout<<"Circle "<<i<<" radius = "<<radius<<" at ( "<<data[i][0]<<", "<<data[i][1]<<" )"<<std::endl;
  }
  //  std::cout<<"=============================================================================="<<std::endl;
}

void draw_found_center(const std::vector<cv::Vec3f>& data, cv::Mat & image){
  for( size_t i = 0; i < data.size(); i++ ) {
    image.at<uchar>(data[i][1], data[i][0]) = 255 ;
  }
}


void draw_foundblobs(const std::vector<cv::Vec3f>& data, cv::Mat & image){
  for( size_t i = 0; i < data.size(); i++ ) {
    circle(image, cv::Point( data[i][0], data[i][1] ), data[i][2], cv::Scalar(0,255,255), -1 );
  }
}




void draw_text_circles(cv::Mat &img, const MedianTextData& mtd){
  for ( const MedianTextRecord & rec : mtd ){
    cv::Point center_text(rec.x(), rec.y());
    //used radius of 10 pixels and green color.                               
    circle( img, center_text, 10, cv::Scalar(0,255,0), 1, 8, 0 );
  }
}

//Returns the data from text file
//Returns empty vector if text file is not supplied.
MedianTextData assign_data_from_text(int argc, std::string argv){
  if(argc==2){MedianTextData a; return a;}

  MedianTextReader *boltreader = MedianTextReader::Get();
  boltreader->set_input_file( std::string( argv ) );
  return boltreader->get_data();
}

//Assign data from Michael's /Dan's code
//Returns empty vector if text file is not supplied.
std::vector<cv::Vec3f> fill_bolts_vector(std::string argv){
  //  if(argc==2){MedianTextData a; return a;}
  std::vector <cv::Vec3f> bolts;
  std::string line;
  std::ifstream boltloc(argv);
  if(boltloc.is_open()){
    while(getline(boltloc, line)){
      int indx;
      std::istringstream iss(line);
      cv::Vec3f temp(1);
      iss>>indx>>temp[0]>>temp[1];
      temp[2]=3; //fake radius;
      bolts.push_back(temp);
    }
    boltloc.close();
  }
  else{ std::cout<<"Unable to open "<<argv<<std::endl;}
  return bolts;
}

void write_to_text(std::string argv, const std::vector<cv::Vec3f> &blobs){
  std::ofstream text((argv+".txt").c_str());
  if(text.is_open()){
    for(unsigned i=0; i<blobs.size();i++){
      text<<"-1\t"<<blobs[i][0]<<"\t"<<blobs[i][1]<<std::endl;
    } 
    text.close();
  }
  else{ std::cout<<"Unable to open "<<(argv+".txt").c_str()<<std::endl;}  
}


//flags to turn on/off saving images
std::vector<bool> setup_verbosity(int option){
  std::vector <bool> options;
  for(int i=0; i<5;i++){
    options.push_back(bool(option%2));
    option /= 10;
  }

  return options;
}


//flags to turn on/off saving images                                                                                                          
std::vector < bool > setup_image_saveflags () {
  std::vector < bool > options;
    try
      {
        int option = config::Get_int ("save_option");

	for (unsigned i = 0; i < 5; i++)
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
