#ifndef _xypoint_h_
#define _xypoint_h_

/// A minimalist way of storing x,y values
/// that could be added to a std::vector
struct xypoint{
  double x;
  double y;
  xypoint( double xx, double yy ) : x(xx), y(yy) { }; 
  xypoint() :x(0), y(0) { };
};

#endif
