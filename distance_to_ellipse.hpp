#ifndef DISTANCE_TO_ELLIPSE_HPP
#define DISTANCE_TO_ELLIPSE_HPP
#include <opencv2/core.hpp>
//just working with one quadrant but slow method.
double dmin(cv::Point c, double e0, double e1, cv::Point p , cv::Point& q, double phi);

//fast method
//c=centre, e0=a, e1=b, p=point, q=query point/closest point in ellipse, phi=angle made by major axis(a). 
//q is a return value of closest point.
double get_distance(const cv::Point& c, const double& a, const double& b, const cv::Point& p , cv::Point& q, const double& phi);
#endif
