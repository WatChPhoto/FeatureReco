#include "to_text_functions.hpp"

double RADTODEG( double R ){const double PI = std::acos(-1); return (180.0 * R) / PI; }


void PMTIdentified::calculate_angles(){
  float a = circ[0];  //x and y at centre
  float b = circ[1];

  for( cv::Vec3f bolt : bolts ) {
    float x = bolt[0];
    float y = bolt[1];
    float theta = atan2f((y-b),(x-a));
    theta = RADTODEG(theta);
    theta = (theta<0)?(theta+360):theta; //getting angle between 0-360
    //Finding angle wrt to Y-axis ^ (nothing to do with axis direction in OpenCv
    theta = (theta<270)?(theta+90):(theta-270);
    angles.push_back( theta );
  }
}

void PMTIdentified::calculate_boltid(){
  float angle_between_bolts = 360.0 / 24; // 24 bolts
  float dang = angle_between_bolts/2;

  for ( float angle : angles ){
    int boltnum = int( (angle+dang) / angle_between_bolts ) + 1;
    if (boltnum==25) boltnum=1;
    boltid.push_back( boltnum );
    // calculate difference in angle from boltid angle
    float boltidang = (boltnum-1) * angle_between_bolts;
    float da = angle - boltidang;
    if ( da > 360.0-dang ) da -= 360.0;
    dangs.push_back( da );
  }
}


