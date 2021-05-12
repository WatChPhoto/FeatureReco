#ifndef __FEATURERECOTTREE__
#define __FEATURERECOTTREE__

#include "TTree.h"
#include "TFile.h"
#include "TBranch.h"
#include "TObject.h"
#include "TClonesArray.h" 
#include <string>
#include <vector>
#include <iostream>

/// Structures to hold information about one image processed by FeatureReco
/// It is of the form that it can be used to set up a TTree readout


/// BlobData class holds the information about one blob
class BlobData : public TObject {
public:
  BlobData() : 
    entry( 0 ), intensity(0), x(0), y(0), area(0), circularity(0), convexity( 0 ), inertia( 0 ) {}     // default constructor
  BlobData( int aentry, float aintensity, float ax, float ay, float aarea, float acirc, float aconv, float ainert ) :
    entry( aentry), intensity(aintensity), x(ax), y(ay), area(aarea), circularity(acirc), convexity( aconv ), inertia( ainert ) { }
  ~BlobData(){}

  int   entry;         //< entry number      
  float intensity;     //< intensity of blob
  float x;             //< x pixel location
  float y;             //< y pixel location
  float area;          //< area of blob
  float circularity;   //< circularity of blob
  float convexity;     //< convexity of blob
  float inertia;       //< inertia of blob

  ClassDef( BlobData, 1 ); // Data for one blob
};

/// Ellipses are stored in a second TTree
class EllipseData : public TObject {
public:
  EllipseData() : 
    entry(0), xx(0), yy(0), aa(0), bb(0), ee(0), phi(0), peakval(0), chi2(0), ndof(0) {}       // default constructor
  EllipseData(int aentry, float axx, float ayy, float aaa, float abb, float aee, float aphi, float apeakval ) : 
    entry(aentry), xx(axx), yy(ayy), aa(aaa), bb(abb), ee(aee), phi(aphi), peakval(apeakval), chi2(0), ndof(0) {}       // default constructor
  ~EllipseData(){}

  int entry;        //< entry number (pmt number if avail)
  float xx;         //< center of ellipse in pixels in x
  float yy;         //< center of ellipse in pixels in y
  float aa;         //< ellipse minor axis length in pixels
  float bb;         //< ellipse major axis length in pixels
  float ee;         //< ellipse eccentricity (0-1)
  float phi;        //< rotation of ellipse about x axis (radians)
  float peakval;    //< peak in hough space for this ellipse
  float chi2;       //< chi2 with an assumed uncertainty of 1 pixel
  float ndof;       //< number of degrees of freedom
  std::vector< BlobData > blobentry; // copy of bolts that are on ellipse

  int size;
  int SetSize() {
    size = sizeof(EllipseData) + blobentry.size()*sizeof(BlobData);
    return size;
  }
  
  ClassDef( EllipseData, 1 );
};



class ImgProcessSettings : public TObject {
public:
  ImgProcessSettings(){
    imgnum =0;           
    do_clahe =0;         
    clahe_gridsize =0;   
    clahe_cliplimit =0;  
    do_blur =0;          
    blur_pixels =0;      
    blur_sigma =0;       
    do_bilat =0;         
    bilat_sigcolor =0;   
    bilat_sigspace =0;   
    blob_min_thres =0;   
    blob_max_thres =0;   
    blob_filter_area =0; 
    blob_min_area =0;   
    blob_max_area =0;
    blob_filter_circ =0; 
    blob_min_circ =0;
    blob_max_circ =0;
    blob_filter_conv =0 ;
    blob_min_conv =0;
    blob_max_conv =0;
    blob_filter_iner =0; 
    blob_min_iner =0;
    blob_max_iner =0;
    ehough_minhits =0;
    ehough_threshold =0; 
    ehough_drscale =0;   
    ehough_nbins_bb =0;   
    ehough_nbins_ee =0;
    ehough_nbins_phiphi =0; 
    ehough_nbins_x =0;      
    ehough_nbins_y =0;      
    ehough_bbmin =0;        
    ehough_bbmax =0;        
    ehough_eemin =0;        
    ehough_eemax =0;        
    ehough_phimin =0;     
    ehough_phimax =0;    
    ehough_xmin =0;         
    ehough_xmax =0;         
    ehough_ymin =0;         
    ehough_ymax =0;    
  }
  ~ImgProcessSettings(){}
  int   imgnum;           //< image number
  int   do_clahe;         //< true (1) if clahe applied
  int   clahe_gridsize;   //< size of box to equalize in
  int   clahe_cliplimit;  //< clahe limits the amplification by clipping 
  int   do_blur;          //< true (1) if gaussian blur
  int   blur_pixels;      //< kernel size of blur in pixels
  int   blur_sigma;       //< sigma of gaussian blur in pixels
  int   do_bilat;         //< true (1) if bilateral filter applied
  int   bilat_sigcolor;   //< range in color to filter over
  int   bilat_sigspace;   //< range of distance to filter over
  int   blob_min_thres;   //< min threshold for blob
  int   blob_max_thres;   //< max threshold for blob
  int   blob_filter_area; //< true (1) if filter by area
  float blob_min_area;   
  float blob_max_area;
  int   blob_filter_circ; //< true (1) if filter by circularity
  float blob_min_circ;
  float blob_max_circ;
  int   blob_filter_conv; //< true (1) if filter by convexity
  float blob_min_conv;
  float blob_max_conv;
  int   blob_filter_iner; //< true (1) if filter by inertia
  float blob_min_iner;
  float blob_max_iner;
  int   ehough_minhits;
  int   ehough_threshold;  //< hough space count to stop finding ellipses
  float ehough_drscale;    //< multiplier to pixel diagonal length 
  int   ehough_nbins_bb;   
  int   ehough_nbins_ee;
  int   ehough_nbins_phiphi; 
  int   ehough_nbins_x;      
  int   ehough_nbins_y;      
  float ehough_bbmin;        
  float ehough_bbmax;        
  float ehough_eemin;        
  float ehough_eemax;        
  float ehough_phimin;     
  float ehough_phimax;    
  float ehough_xmin;         
  float ehough_xmax;         
  float ehough_ymin;         
  float ehough_ymax;    
  ClassDef(ImgProcessSettings,1);
};


/// Image data TTree branch entry
/// one entry per image and arrays of blobs and ellipses
/// holds meta data about the image, and the processing parameters
class ImageData : public TObject {
public:
  ImageData(){}
  ~ImageData(){}
  ImgProcessSettings         ips;  
  std::vector< BlobData >    fBlobs;
  std::vector< EllipseData > fEllipses; // all ellipses found
  std::vector< EllipseData > fPMTs;     // final ellipses identified as PMTs

  void AddBlob( const BlobData& b ){ fBlobs.push_back( b ); }
  void AddEllipse( const EllipseData& e ){ fEllipses.push_back( e ); }
  void AddPMT( const EllipseData& e ){ fPMTs.push_back( e ); }

  int fEntrySize;
  void SetSize() { 
    fEntrySize = sizeof(ImageData) + 
      fBlobs.size() * sizeof(BlobData);
    for ( EllipseData& ed : fEllipses ){
      fEntrySize += ed.SetSize();
    }
    for ( EllipseData& ed : fPMTs ){
      fEntrySize += ed.SetSize();
    }
    std::cout<<"ImageData::SetSize = "<<fEntrySize<<std::endl;
  }
  void Clear_data(){
    fBlobs.clear();
    fEllipses.clear();
  }
    
  ClassDef( ImageData, 1 ); // Data for one image
};

#endif // __FEATURERECOTTREE__
