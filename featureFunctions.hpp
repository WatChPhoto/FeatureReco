#ifndef FEATUREFUNCTIONS_HPP
#define FEATUREFUNCTIONS_HPP

#include <opencv2/core.hpp>
#include "MedianTextReader.hpp"
#include "PMTIdentified.hpp"
#include <TH1D.h>
#include <TH2D.h>


//Outputs the filename
std::string build_output_filename( const std::string& in, const std::string& tag );
std::string build_output_textfilename( const std::string& in, const std::string& tag );

//Histogram of the intensity.
void histogram_channel0( const cv::Mat& img,  std::string hname="hist" );

//Applies threshold to the image. Threshold is determined using previous histogram.
void apply_image_threshold( cv::Mat& img, int threshold=250 );

//Calculate metric for filtering the real bolt.
double calculate_bolt_metric( const cv::Vec3f& circ, const cv::Mat& img );

//Makes the histogram for bolt metric found vs distance.
void make_bolt_metric_histograms( const std::vector< PMTIdentified >& pmtids,  cv::Mat &imbw, TH1D *&metric_all, TH1D *&metric_good, TH1D *&metric_bad, TH2D *&metric_2d);

void draw_line(const std::vector< PMTIdentified >& pmtids, const MedianTextData& mtd, cv::Mat &imcol );

//Bolt metric for the inbetween points that are not mapped.
//void histogram_inbetween(const std::vector<cv::Vec3f>& circles, const MedianTextData& mtd, const std::vector< IndexMatchDist >& data121, cv::Mat imbw, TH1D *&metric_inb);

//Draws circle from data containing x,y,r info
void draw_circle_from_data(const std::vector <cv::Vec3f> data, cv::Mat & image, cv::Scalar color = (0,255,0), int line_width=1);

//Draws point from data containing x,y,r value at x,y  
void draw_found_center(const std::vector<cv::Vec3f> data, cv::Mat & image);

//Draw circle from text file
void draw_text_circles(cv::Mat &img, const MedianTextData& mtd);

//Returns the data from text file
//Returns empty vector if text file is not supplied.
MedianTextData assign_data_from_text(int argc, std::string argv);

//flags to turn on/off saving images
std::vector<bool> setup_verbosity(int option);
#endif
