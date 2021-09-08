#ifndef _camera_hpp_
#define _camera_hpp_

#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <fstream>
#include <cmath>
#include <iostream>
#include <vector>
#include <opencv2/calib3d.hpp>
#include <string>

const double PI=std::acos(-1);

class WorldPoints{
public:
  std::string id;
  double x;
  double y;
  double z;
  
  WorldPoints():id(""),x(0),y(0),z(0){}
  WorldPoints(std::string id, double x, double y, double z):id(id),x(x),y(y),z(z){}
};


//data structure
class Ellipse{
public:
  std::string id; //PMT id
  double x;  //centre of PMt
  double y;
  double b;  //minor axis of ellipse
  double e;  //eccentricity of ellipse
  double phi; //angle with x-axis.

  //TransformationData();
  Ellipse():id(""),x(0),y(0),b(0),e(0),phi(0){};
  Ellipse(std::string id, double x, double y, double b, double e, double phi):id(id),x(x),y(y),b(b),e(e),phi(phi){}
  
  double get_a(){
    return b/std::sqrt(1.-e*e);
  }
};


class cam{
public:
  double r;
  double theta;
  double z;
  double yaw;
  double pitch;
  double roll;
  cv::Matx33d R; //rotation vectro
  cv::Matx31d tv; //translation vector
  cv::Matx31d rv_orig; //original rotation vectro
  cv::Matx31d tv_orig; //original translation vector
  bool button_clk=false;

  cam();
  cam(double r, double theta, double z, double yaw, double pitch, double roll);
  void reset_rvec_tvec();
  void set_camera(double r, double theta, double z, double yaw, double pitch, double roll);
  void increment(std::string var, float step);
  void decrement(std::string var, float step);
  cv::Mat get_scene(const cv::Mat& sc);  
  void print();
  void set_world_points(std::vector<WorldPoints> p);
  bool show_background();
  void flip_background();
private:
  void calculate_R_tv();
  unsigned clip(cv::Point3f n,cv::Matx31d p0,cv::Matx31d p,cv::Matx31d v,cv::Point3f& start,cv::Point3f& end);
  std::vector<WorldPoints> all_pmts;
  bool bkg_opt=false;
  //camera parameters
  const double fx=2.760529621789217e+03;
  const double fy=2.767014510543478e+03;
  const double cx=1.914303537872458e+03;
  const double cy=1.596386868474348e+03;
  
  const cv::Matx33d camera_matrix=cv::Matx33d(fx,0.,cx,
					      0.,fy,cy,
					      0.,0.,1.);
  const cv::Matx14d dist_coeffs=cv::Matx14d(-0.2398, 0.1145, 0, 0);
};

//Icon class
class Icon{
public:
  cv::Point2f pos;
  std::string name;
  cv::Scalar text_color;
  cv::Scalar color;
  cv::Point2f size;
  char sign; //true for + and false for -
  Icon();  
  Icon(std::string name, cv::Point2f pos, cv::Point2f size, cv::Scalar color, cv::Scalar text_color);  
  void show(int x, int y,bool clicked,cv::Mat& img, cam& c1); 
  void set(std::string name, cv::Point2f pos, cv::Point2f size, cv::Scalar color, cv::Scalar text_color);
private:
  bool is_inside(int x, int y);
};

#endif
