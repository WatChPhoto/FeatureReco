#ifndef TO_TEXT_FUNCTIONS_HPP
#define TO_TEXT_FUNCTIONS_HPP

#include<opencv2/core.hpp>
#include<vector>



struct PMTIdentified {
  int           pmtid; // -1 if not identified
  cv::Vec3f         circ; // [0]=xc, [1]=yc, [2]=radius
  std::vector<cv::Vec3f> bolts; // bolts going with this PMT 
  std::vector<float> dists; // distance of bolt from PMT circle
  std::vector<float> angles; // angle of each bolt
  std::vector<float> dangs;  // difference in angle from boltid's angle
  std::vector<int>   boltid; // 1 is at 12 o'clock, 2 ... 24 going around clockwise

  PMTIdentified( const cv::Vec3f& pmtloc, const std::vector< cv::Vec3f >& boltlocs, const std::vector< float > final_dists ) :
    pmtid(-1), circ( pmtloc ), bolts( boltlocs ), dists( final_dists ) {
    calculate_angles();
    calculate_boltid();
  }

private:

  void calculate_angles();
  void calculate_boltid();
};

double RADTODEG( double R ); 

#endif
