#ifndef FEATUREFUNCTIONS_HPP
#define FEATUREFUNCTIONS_HPP

#include <opencv2/core.hpp>
#include "MedianTextReader.hpp"
#include <TH1D.h>
#include <TH2D.h>


//Outputs the filename
std::string build_output_filename( const std::string& in, const std::string& tag );

//Histogram of the intensity.
void histogram_channel0( const cv::Mat& img,  std::string hname="hist" );

//Applies threshold to the image. Threshold is determined using previous histogram.
void apply_image_threshold( cv::Mat& img, int threshold=250 );

//For making one to one match between found circle and text circle.
/// Simple structure to hold indices into two vectors of matches with distance between those two indices stored
struct IndexMatchDist {
  unsigned idx_txt;   //index of medianTextReader
  unsigned idx_circ;  //index of circles
  float    dist;
  IndexMatchDist( unsigned itxt=0, unsigned icirc=0, float fdist=-1. ) : idx_txt( itxt ), idx_circ( icirc ), dist( fdist ) { }
};
std::vector< IndexMatchDist > find_closest_matches( const std::vector<cv::Vec3f>& circles, const MedianTextData & mtd );

//Bolt distance histogram. has one to many mapping.
void make_bolt_dist_histogram( const std::vector< IndexMatchDist > & matches, TH1D *&hout);

//Calculate metric for filtering the real bolt.
double calculate_bolt_metric( const cv::Vec3f& circ, const cv::Mat& img );

//Makes histogram from text file to found circle
void make_bolt_dist_histogram_wrt_txt( const std::vector<cv::Vec3f>& circles, const MedianTextData& mtd, TH1D *&hist_dist );


//Makes the histogram for bolt metric found vs distance.
void make_bolt_metric_histograms( const std::vector<cv::Vec3f>& circles, const MedianTextData& mtd, const std::vector< IndexMatchDist >& data121,  cv::Mat &imbw,  cv::Mat &imcol, TH1D *&metric_all, TH1D *&metric_good, TH1D *&metric_bad, TH2D *&metric_2d);

//Bolt metric for the inbetween points that are not mapped.
void histogram_inbetween(const std::vector<cv::Vec3f>& circles, const MedianTextData& mtd, const std::vector< IndexMatchDist >& data121, cv::Mat imbw, TH1D *&metric_inb);

//Draws circle from data containing x,y,r info
void draw_circle_from_data(const std::vector <cv::Vec3f> data, cv::Mat & image, cv::Scalar color = (0,255,0));

//Draws point from data containing x,y,r value at x,y  
void draw_found_center(const std::vector<cv::Vec3f> data, cv::Mat & image);

//Draw circle from text file
void draw_text_circles(cv::Mat &img, const MedianTextData& mtd);
#endif
