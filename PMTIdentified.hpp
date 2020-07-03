#ifndef _PMTIdentified_hpp_
#define _PMTIdentified_hpp_

#include <iostream>
#include <vector>
#include <opencv2/imgproc.hpp>
#include "MedianTextReader.hpp"
#include <TH1D.h>

const float bad_dmin = 500.; //pixels
const double PI = std::acos(-1);
inline double RADTODEG( double R ){ return (180.0 * R) / PI; }

/// PMTIdentified holds locations and radius of PMTs that were identified
/// For each pmt found it stores a vector of bolts found that are around this PMT
/// Identifies the bolt number by the angle clockwise from 12 o'clock
struct PMTIdentified {
  int                pmtid; // -1 if not identified
  cv::Vec3f          circ; // [0]=xc, [1]=yc, [2]=radius
  std::vector<cv::Vec3f> bolts; // bolts going with this PMT
  std::vector<float> dists; // distance of bolt from PMT circle 
  std::vector<float> angles; // angle of each bolt
  std::vector<float> dangs;  // difference in angle from boltid's angle
  std::vector<int>   boltid; // 1 is at 12 o'clock, 2 ... 24 going around clockwise

  // include comparison to truth if available
  std::vector<int>   idx_txt; // index of medianTextReader
  std::vector<float> dist_txt; // distance to closest matching bolt

  PMTIdentified() : pmtid(-1) { }
  PMTIdentified( const cv::Vec3f& pmtloc, const std::vector< cv::Vec3f >& boltlocs, const std::vector< float > final_dists ) : 
    pmtid(-1), circ( pmtloc ), bolts( boltlocs ), dists( final_dists ) {
    calculate_angles();
    calculate_boltid();
  }
  PMTIdentified( const PMTIdentified& pmt ) :
    pmtid(-1), circ( pmt.circ ), bolts( pmt.bolts ), dists( pmt.dists ),
    angles( pmt.angles ), dangs( pmt.dangs ), boltid( pmt.boltid ),
    idx_txt( pmt.idx_txt ), dist_txt( pmt.dist_txt ) {}

private:

  void calculate_angles();
  void calculate_boltid();
};


	
/// Take circles_of_blob to select which bolts are good blobs
///	   inputs: 
///	   vector<Vec3f> blobs;  // (x,y,r)
///	   vector<Vec3f> circles_of_blob;  // (x,y,r)
///
///	   output:
///	   vector< PMTIdentified > final_bolts; // bolt locations selected
///	 
/// only accept PMTs that are more than some number of pixels away from edge of image
const unsigned trim_pixels = 10;
void find_candidate_bolts( const std::vector< cv::Vec3f >& blobs, 
			   const std::vector< cv::Vec3f >& circles_of_blob,
			   std::vector< PMTIdentified >& pmts_found,
			   const cv::Mat& image );

/// Output operator to output the information about a PMTIdentified
std::ostream& operator<<( std::ostream& os, const PMTIdentified& p );

/// Function to match bolts on PMTIdentfied objects to truth ones in text file
void find_closest_matches( std::vector<PMTIdentified>& final_pmts, const MedianTextData & mtd );

//Bolt distance histogram. has one to many mapping.
void make_bolt_dist_histogram( const std::vector<PMTIdentified > & matches, TH1D *&hout);

//Makes histogram from text file to found circle
//void make_bolt_dist_histogram_wrt_txt( const std::vector<cv::Vec3f>& circles, const MedianTextData& mtd, TH1D *&hist_dist );

//Overlays the boltangle and boltid in the image
void overlay_bolt_angle_boltid(const std::vector< PMTIdentified >& final_pmts, cv::Mat image_final);


/// look for duplicate bolts and keep only best matches
void prune_bolts( std::vector< PMTIdentified >& final_pmts, float ang_offset );
  
/// remove bolts and PMTs from PMTs with fewer than some threshold bolts
void prune_pmts(  std::vector< PMTIdentified >& final_pmts, unsigned numbolts, const std::string& label );


#endif

