#ifndef ELLIPSE_DETECTION_HPP
#define ELLIPSE_DETECTION_HPP

#include<opencv2/core.hpp>
#include<opencv2/imgproc.hpp>

double distance(cv::Point p1, cv::Point p2);
class ParametricEllipse{
public:
  cv::Point2i centre;
  double a;
  int b;
  double alpha;
  int freq;
  std::vector<cv::Vec3f> bolts;
  std::vector<cv::Point2i> query;
  std::vector<float> dist;

  ParametricEllipse(cv::Point2i centre, int a, int b, int alpha, int freq  );
};

bool has_key(const std::vector<cv::Point2i>& bbins, int b, int& index);
void find_max(const std::vector<cv::Point2i>& bbins, int& max_freq, int& max_index);
std::vector<ParametricEllipse> detect_ellipse( const std::vector<cv::Vec3f>& input_data, const int min_major, const int max_major, const int min_minor, const int max_minor,int min_minor_freq);
void draw_ellipses(const std::vector<ParametricEllipse>& elData, cv::Mat& img );

#endif
