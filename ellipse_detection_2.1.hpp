#ifndef ELLIPSE_DETECTION_HPP
#define ELLIPSE_DETECTION_HPP

#include<opencv2/core.hpp>
#include<opencv2/imgproc.hpp>

double distance(cv::Point p1, cv::Point p2);
class ParametricEllipse{
public:
  std::vector<cv::Vec3f> xypoints;
  cv::Point2i centre;
  int a;
  int b;
  int alpha;
  int freq;
  ParametricEllipse(cv::Point2i centre, int a, int b, int alpha, int freq  );
  void set_xypoints(std::vector<cv::Vec3f> xypoints);
};

bool has_key(const std::vector<cv::Point2i>& bbins, int b, int& index);
void find_max(const std::vector<cv::Point2i>& bbins, int& max_freq, int& max_index);
void detect_ellipse( std::vector<cv::Vec3f> coordinates, cv::Mat& img, const int min_major, const int max_major, const int min_minor, const int max_minor,int min_minor_freq);

#endif
