#ifndef _ellipse_h_
#define _ellipse_h_

/// Simple class to hold ellipse parameters and calculated distance from
/// a point to an ellipse.

/// Authors:  Blair Jamieson, Andrew Sikora
/// Dates: Summer 2018

#include <iostream>
#include <cmath>
#include "xypoint.hpp"
#include "opencv2/imgproc.hpp"

using std::ostream;
using std::acos;

const double pi = acos(-1.);


struct ellipse_st {
  // define an ellipse by:
  //   bb == shorter axis length
  //   ee == eccentricity
  //   phiphi == rotation of ellipse about x axis
  //   cc == xypoint of center of ellipse
  ellipse_st( double bb, double ee, double phiphi, xypoint cc );
  ellipse_st() :  b(1.), e(0.), phi(0.), c( xypoint(0.,0.) ), epenalty(0.) { }
  ellipse_st( const ellipse_st& elli ) : 
    b( elli.get_b() ), e( elli.get_e() ), phi( elli.get_phi() ), c( elli.get_xy() ), epenalty(0.)  { }
  ellipse_st( const cv::Vec3f &v ) :
    b( v[2] ), e( 0. ), phi( 0. ), c( xypoint( v[0], v[1] ) ), epenalty(0.)  { }


  void set_bephi( double bb, double ee, double phiphi ){ b=bb; e=ee; phi=phiphi; } 

  // parametric function radius of ellipse 
  // at theta == angle from long axis in radians
  double r_of_theta( double theta ) const;

  // get x, y coordinates on ellipse at
  // theta == angle from long axis in radians
  xypoint xy( double theta ) const;

  // find smallest distance from arbitrary point (px, py) to
  // a point on the ellipse.

  double dmin( xypoint p ) const;

  // faster method of finding smallest distnce-squared from point p 
  // to a point on the ellipse. From:
  // https://www.geometrictools.com/Documentation/DistancePointEllipseEllipsoid.pdf
  double fast_dmin2 (xypoint p) const;

  //function needed for fast_dmin
  double Findtbar (double a, xypoint pprime) const;

  //function needed for fast_dmin
  double F (double a, xypoint pprime, double t0, double t1) const;


  double GetRoot(double r0, double z0, double z1, double g) const;

  double DistancePointEllipse2(xypoint p) const;

  // return penalty term in a chi2 for having a bad eccentricity parameter
  double get_epenalty() const { return epenalty; }

  //return semimajor axis of ellipse
  double get_b() const {return b;}

  //return semiminor axis of ellipse
  double get_a() const {return b/std::sqrt(1-e*e) ;}

  // return angle of ellipse from x-axis
  double get_phi() const {return phi;}

  // return eccentricity of ellipse
  double get_e() const{ return e; }

  // return center of ellipse
  xypoint get_xy() const{ return c; }

  // set xy of center of ellipse
  void set_xy( const xypoint& xy ){ c = xy; }

  // find smallest distance-squared from arbitrary point (px, py) to
  // a point on the ellipse.
  double dmin2( xypoint p ) const;

  // find closest point on ellipse that is smallest distance from 
  // arbitrary point (px, py) to a point on the ellipse.
  xypoint closest_point( xypoint p ) const;


  // treat as Vec3f circle [0]=xc, [1]=yc, [2]=r
  float operator[]( unsigned idx ) const{
    switch (idx){
    case 0: 
      return c.x;
    case 1:
      return c.y;
    case 2:
      return b;
    default:
      return 0.;
    }
  }
      
		 

  friend ostream& operator<<( ostream& os, const ellipse_st &  e );

private:
  double  b; // short axis length
  double  e; // eccentricity 0=circle, 1=parallel lines
  double  phi; // rotation of ellipse from x axis
  xypoint c; // center

  double  epenalty; // penalty term for bad eccentricity
};


std::ostream& operator<<( std::ostream& os, const ellipse_st& e );


#endif
