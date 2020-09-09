#ifndef ELLIPSE_INTERSECT
#define ELLIPSE_INTERSECT

#include "ellipse.hpp"
#include "xypoint.hpp"
#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <functional>

enum
  {
    ELLIPSES_SEPARATED,
    ELLIPSES_OVERLAP,
    ELLIPSE0_OUTSIDE_ELLIPSE1_BUT_TANGENT,
    ELLIPSE0_STRICTLY_CONTAINS_ELLIPSE1,
    ELLIPSE0_CONTAINS_ELLIPSE1_BUT_TANGENT,
    ELLIPSE1_STRICTLY_CONTAINS_ELLIPSE0,
    ELLIPSE1_CONTAINS_ELLIPSE0_BUT_TANGENT,
    ELLIPSES_EQUAL
  };

class EllipseIntersect{
public:
  double d0;
  double d1;
  double c0;
  double c1;
   
  EllipseIntersect():d0(0.),d1(0.),c0(0.),c1(0.){}
  int intersect(const ellipse_st e0, const ellipse_st e1);

private:
  double F(double s);
  double DF(double s);
  int Classify(double minSqrDistance, double maxSqrDistance, double z);
  double BisectionRoot(double smin, double smax,double f0, double f1, unsigned int maxIterations, double &s);
  void GetRoots(const double &x1, const double &x2, const double &y1, const double &y2, int &numRoots, double *roots);
  void GetRoots(const double &x1, const double &y1, int &numRoots, double *roots);

};
#endif
