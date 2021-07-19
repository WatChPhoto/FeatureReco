/*#############################################################################
#   Here I changed some functions in "struct ImageDataReader {" to public.    #
#   I commented out private. You can change it back if it is more suitable    #
#   for your application.                                                     #
##############################################################################*/
#ifndef _ImageDataReader_hpp_
#define _ImageDataReader_hpp_

#include <string>
#include <iostream>
#include <vector>
//#include <opencv2/core.hpp>
//#include <opencv2/highgui.hpp>

struct CameraFace {
  float yaw;
  float pitch;
  float roll;
  CameraFace( float y, float p, float r ) : yaw(y), pitch(p), roll(r) { }
  CameraFace() : yaw(0.0), pitch(0.0), roll(0.0){ }
  //cv::Matx33d GetRotationMatrix();
  //cv::Matx13d GetVectorAlongCameraFacing();
  //operator+=
  
};

struct ImageMetaData {
  std::string survey;    // survey name
  unsigned    id;        // survey id
  std::string filename;  // filename
  unsigned    imgnum;    // 3 digit image number (non-unique)--- //not 3 digit in case the name is 045 or 001 ritht? --Tapendra
  unsigned    exposure;
  unsigned    iso;
  float       fstop;
  unsigned    epoch;
  unsigned    year;
  unsigned    month;
  unsigned    day;
  unsigned    hr;
  unsigned    min;
  unsigned    sec;
  unsigned    dt;     // time delta
  float       depth;  // meters?
  float       T;      // temperature
  CameraFace face;
  ImageMetaData() : survey(""), id(0), filename(""), imgnum(0), exposure(0), iso(0),
		    fstop(0.0), epoch(0), year(0), month(0), day(0), hr(0), min(0), sec(0),
		    dt(0), depth(0.0), T(0.0), face(0.0,0.0,0.0) {}
};

std::istream& operator>>( std::istream& istream, ImageMetaData& dt );

// static image metadata holder
struct ImageDataReader {
public:
  std::vector< ImageMetaData > md;
  static ImageDataReader& GetInstance();

  ImageDataReader(ImageDataReader const&) = delete;
  void operator=(ImageDataReader const&)  = delete;
  
  CameraFace GetFace( std::string surveyid, unsigned imgnum ) const;
  ImageMetaData GetMetaData(std::string surveyid, unsigned imgnum) const;
private:// --Just this change.
  ImageDataReader();
  ~ImageDataReader(){ if (instance) delete instance; }
  static ImageDataReader* instance;
};
#endif
