#include "ImageDataReader.hpp"
//#include <opencv2/core.hpp>
//#include <opencv2/highgui.hpp>
#include <fstream>                                                             
#include <cmath>
#include <iostream>
//#include <opencv2/calib3d.hpp>
#include<algorithm>

ImageDataReader * ImageDataReader::instance = nullptr;


long int get_image_num_from_filename( const char* filename ){
  int idx_start;
  int idx_end=-1;
  char * pEnd;
  while ( filename[++idx_end] != '.' ){}
  idx_start=idx_end;
  while ( isdigit( filename[--idx_start] ) ){}
  ++idx_start;
  long int i = strtol( &filename[idx_start], &pEnd, 10 );
  return i;
}

  
std::istream& operator>>( std::istream& is, ImageMetaData& dt ){

  is >> dt.survey;
  is >> dt.id;
  is >> dt.filename;
  long int imgnum = get_image_num_from_filename( dt.filename.c_str() );
  dt.imgnum = unsigned( imgnum/10000 );
  std::string junk;
  is >> junk;
  is >> dt.exposure;
  is >> dt.iso;
  is >> dt.fstop;
  is >> dt.epoch;
  is >> dt.year;
  is >> dt.month;
  is >> dt.day;
  is >> dt.hr;
  is >> dt.min;
  is >> dt.sec;
  is >> dt.dt;
  is >> dt.depth;
  is >> dt.T;
  is >> dt.face.yaw;
  is >> dt.face.pitch;
  is >> dt.face.roll;
  return is;

}




ImageDataReader& ImageDataReader::GetInstance(){
  if ( instance == nullptr ){
    instance = new ImageDataReader();
  }
  return *instance;
}

ImageDataReader::ImageDataReader(){
  std::ifstream in( "ImageMetaData.tsv" );
  ImageMetaData dt;
  while ( in >> dt ){
    md.push_back( dt );
  }
  in.close();
}

CameraFace ImageDataReader::GetFace( std::string surveyid, unsigned imgnum ) const{
  //CameraFace facing;
  unsigned imgcount = 0;
  float yaw=0.;
  float pitch =0.;
  float roll = 0.;
  for ( const ImageMetaData &imd : md ){
    if ( imd.survey == surveyid && imd.imgnum == imgnum ){
      ++imgcount;
      //std::cout<<imd.survey<<", "<<imgnum<<std::endl;
      yaw += imd.face.yaw;
      pitch += imd.face.pitch;
      roll += imd.face.roll;	
    }
  }
  
  if(imgcount!=0){
    yaw /= imgcount;
    pitch /= imgcount;
    roll /= imgcount;
  }

  return CameraFace(yaw,pitch,roll);
}

ImageMetaData ImageDataReader::GetMetaData( std::string surveyid, unsigned imgnum ) const{
  //CameraFace facing;
  ImageMetaData m;
  for ( const ImageMetaData &imd : md ){
    if ( imd.survey == surveyid && imd.imgnum == imgnum ){
      m=imd;
      break;
    }
  } 

  return m;
}


