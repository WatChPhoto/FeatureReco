/// A simple class to do Hough transform on an ellipse
/// Author:  Blair Jamieson (Jun. 2020)
#ifndef _hough_ellipse_h_
#define _hough_ellipse_h_

#include "xypoint.hpp"
#include "ellipse.hpp"
#include <TH2S.h>
#include <vector>
#include <utility>

//enum   HoughEllipseResultType {HoughEllipse, HoughEllipseUnusedPoints };
const unsigned HoughEllipse = 0;
const unsigned HoughEllipseUnusedPoints = 1;

struct binindices_st {
  unsigned ibb;
  unsigned iee;
  unsigned iphi;
  unsigned ix;
  unsigned iy;
  unsigned pkval;
  binindices_st() : ibb(0), iee(0), iphi(0), ix(0), iy(0), pkval(0) { } 
  binindices_st( unsigned bb, unsigned ee, unsigned phi, unsigned x, unsigned y, unsigned pk) : 
    ibb(bb), iee(ee), iphi(phi), ix(x), iy(y), pkval(pk) { } 
  binindices_st( const binindices_st& bi ) : 
    ibb(bi.ibb), iee(bi.iee), iphi(bi.iphi), ix(bi.ix), iy(bi.iy), pkval(bi.pkval) { } 
  
  bool operator<( const binindices_st& bi ) const{
    return ( pkval < bi.pkval );
  }
  
};



/// Define structure type to store results in
/// A single of Hough transform is a vector of data points
/// and a set of ellipse parameters
struct HoughEllipseResult {
  HoughEllipseResult( const ellipse_st& el, float pkval ) :
    e( el ), peakval(pkval), eltype(HoughEllipseUnusedPoints) { }

  HoughEllipseResult() : e( ellipse_st() ),  peakval(0), eltype(HoughEllipseUnusedPoints) { }
  
  std::vector< xypoint > data;      // the points for this result
  ellipse_st             e;         // ellipse parameters
  float                  peakval;   // peak height in hough space
  unsigned               eltype;      // result type (0 for ellipse, 1 for unused points)
  unsigned               ibb{0};    // b hough array index
  unsigned               iee{0};    // e hough array index
  unsigned               iphi{0};   // phi hough array index
  unsigned               ix{0};     // x hough array index
  unsigned               iy{0};     // y hough array indxe
};

/// Vector of Hough transform results
typedef std::vector< HoughEllipseResult >  HoughEllipseResults ; 

/// number of ellipses found
unsigned num_ellipses( const HoughEllipseResults & hrs );     

/// print results
std::ostream& operator<<( std::ostream& os, const HoughEllipseResults& hrs );
std::ostream& operator<<( std::ostream& os, const HoughEllipseResult& hr );


/// A simple class to do Hough transform on an ellipse
/// Example usage, using default radius and bins range.
/// Data is provided as a std::vector< xypoint > data
///
/// EllipseHough h;
/// HoughEllipseResults res = h.find_ellipses( data );
///
/// cout<<"Number of ellipses found is "<<res.num_ellipses()<<endl;
/// // Loop over ellipses and do a fit to the hits
/// for (int ires=0; ires<res.data.size(); ++ires ){
///     if ( res[ires].eltype == HoughEllipse ){
///         ellipse_fit( res[ires].data );
///
/// Parameters that can be tweaked to tune the transform:
///   minhits:   Minimum number of hits to keep searching for ellipse
///              (default value 5)
///   threshold: Minimum peak height in transform space to call a result
///              (default value 8)
class EllipseHough {
 public:
  /// Constructor sets parameters of the Hough Transform
  /// can set number of bins in:
  ///       bb = shorter axis length
  ///       ee = eccentricity (0=circle, 1=ellipse)
  ///   phiphi = rotation of ellipse about x-axis (radians 0-180 degrees)
  /// center x = center location in x
  /// center y = center location in y
  /// along with range in each of those values
  EllipseHough(
	       //unsigned nbins_bb     = 30,   float bbmin=81.0,   float bbmax=141.0,
	       //unsigned nbins_ee     =  8,   float eemin=0.0,    float eemax=0.32,
	       unsigned nbins_bb     = 40,   float bbmin=101.0,   float bbmax=221.0,
	       unsigned nbins_ee     =  10,   float eemin=0.0,    float eemax=0.6,
		unsigned nbins_phiphi = 10,   float phphimin=0.0, float phiphimax=std::acos(-1),
		unsigned nbins_x      = 2300, float xmin=200,     float xmax=3800,
		unsigned nbins_y      = 1300, float ymin=200,     float ymax=2800);

  ~EllipseHough();

  /// Tunable parameters
  /// Minimum hits to be called a circle
  /// default is 3 hits
  void set_minhits( int n ){ minhits = n; }
  /// Minimum height in hough space to be called a circle
  /// default is 5
  void set_threshold( int n ){ threshold = n; }

  /// Scaling factor for distance from ellipse for hit to be part of ellipse
  /// Value of 2.0 is default, uses sqrt(dx*dx+dy*dy) where dx, dy are
  /// the xc,yc bin sizes.
  void set_distance_factor( float dist_factor ){ drscaling=dist_factor; }
  
  /// Main method to find circles
  const HoughEllipseResults& find_ellipses( const std::vector< xypoint >& data );

  /// Get the histogram of the transformed data
  //std::vector< std::vector< TH2S* > > get_transform() { return fTransformed; }
  //std::vector< std::vector< ellipse_st* > > get_ebins() { return fE; }

  // private:
  // Save binning info
  unsigned fNbb;   float fbbmin;   float fbbmax;
  unsigned fNee;   float feemin;   float feemax;
  unsigned fNphi;  float fphimin;  float fphimax;
  unsigned fNx;    float fxmin;    float fxmax;
  unsigned fNy;    float fymin;    float fymax;

  /// Parameters
  int minhits{9};
  int threshold{27};
  float drscaling{4.0};
  
  /// Store results
  HoughEllipseResults fresults;

  // Histogram used to store the transformed data

  // one XY histogram per bin
  // Store each radius in a separate TDirectory (nbins_bb directories == 80)
  //std::vector< TDirectory* > fDirectories;
  // each TDirectory has nbins_ee * nbins_phiphi histograms == 10*18 = 180
  //std::vector< std::vector< TH2S* > > fTransformed;
  // store ellipse at origin for each xy histogram
  //std::vector< std::vector< ellipse_st* > > fE;


  // Replace histograms with an array
  unsigned short ***** fTransformed; //[fNbb][fNee][fNphi][fNx][fNy]
  ellipse_st *** fE; //[fNbb][fNee][fNphi]


  // top directory to store results
  TDirectory * houghdir;

  // helper methods

  // Fill a HoughResult, pass it a copy of all the hits that
  // and any that don't get used in the HoughResult are returned
  // as unused hits
  std::vector< HoughEllipseResult > find_maximum( std::vector< xypoint >& unused_hits );

  // Reset the hough transformed histograms and fill them
  // with the data passed
  void hough_transform( const std::vector< xypoint >& data );

  // Make a clone of one of hough transformed histograms for
  // a candidate ellipse.  Pass it the candidate number, and
  // histogram pointer
  void save_hough_histo( unsigned num, const HoughEllipseResult & hr );

  // make a histogram of the candidate ellipse points, and draw
  // the hough estimate for ellipse
  void plot_candidate( unsigned num, const HoughEllipseResult & hr );

  // following return -1 if outside bin range
  int get_bbin  ( double bb ) const;
  int get_ebin  ( double ee ) const;
  int get_phibin( double phi ) const;
  int get_xbin  ( double x ) const;
  int get_ybin  ( double y ) const;

  float get_b_frombin( int bin ) const;
  float get_e_frombin( int bin ) const;
  float get_phi_frombin( int bin ) const;
  float get_x_frombin( int bin ) const;
  float get_y_frombin( int bin ) const;


private:
  // helper
  void zero_hough_counts();
  
  //  void hough_transform_phibin( unsigned iphi,  const std::vector< xypoint >& data );
  
  //  void find_maximum_phibin( unsigned iphi, std::vector< binindices_st >& passed_threshold );
    
}; /// end of EllipseHough
  



#endif
