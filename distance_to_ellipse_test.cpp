#include<iostream>
#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <time.h>
#include <random>
#include "distance_to_ellipse.hpp"
#include "ellipse_detection_2.1.hpp"

int main(){
  

  const double PI = acos(-1);
  srand (time(NULL));
 
  //a between 50 - 200;
  //centre between 500 to 800 px
  cv::Point c = cv::Point(1000,1000);
  double a = 150;//(rand()%300+50);
  double b =150; //(rand()%300+50);
  double phi = (rand() %361 - 180.0)*(PI+0.0)/180.0;
  std::cout<<"before phi =" <<phi<<std::endl;

  //fill in the data
  ParametricEllipse pa(c, a, b, phi, 1);
  //testEl.push_back(pa);

  //fill in random points to calculate distance from.
  std::vector<cv::Vec3f> points;
  for(int i =0; i<100; i++){
    cv::Vec3f temp {rand()%1500+250, rand()%1500+250, 3};
    points.push_back(temp);
  }
  
  for(int j=0; j<points.size(); ++j){
    cv::Point p = cv::Point(points[j][0], points[j][1]); //bolt's centre.
    cv::Point q;
    double dist =  get_distance(c, a, b, p , q, phi); //get the distance of the bolt from current ellipse
    std::cout<<"distance = "<<dist<<std::endl;
    pa.bolts.push_back(points[j]);
    pa.query.push_back(q);
    pa.dist.push_back(dist);
   }

  cv::Mat test(cv::Size(2000, 2000), CV_64FC3);
  test =  cv::Scalar::all(255) - test;

  
  cv::Size axes( a, b );
  cv::Point center = c;
    
  ellipse( test, pa.centre, axes, 180.0*pa.alpha/PI , 0., 360,  cv::Scalar (0, 0, 0) );
  std::cout<<"after phi =" <<pa.alpha<<std::endl;
    
  //drawing line from bolts to closest point in the ellipse
  for(int i = 0; i<pa.bolts.size(); i++){
    cv::Vec3f bolts = pa.bolts[i];
    //drawing included bolts
    cv::circle( test, cv::Point( bolts[0], bolts[1] ), bolts[2], cv::Scalar(0,0,255), 3, 0 );
    
    //making different colored line for inside and outside
    double value = (std::pow(((bolts[0]-pa.centre.x)*cos(pa.alpha)+(bolts[1]-pa.centre.y)*sin(pa.alpha)),2)/(a*a+0.0))+(std::pow(((bolts[0]-pa.centre.x)*sin(pa.alpha)-(bolts[1]-pa.centre.y)*cos(pa.alpha)),2)/(b*b+0.0))-1.0;
    cv::Scalar col;
    if(value<0){col = cv::Scalar(147,20,255);}
    else{col = cv::Scalar(255,0,0);}

    //drawing line from bolt to closest point in the ellipse
    line(test, cv::Point( cvRound(bolts[0]), cvRound(bolts[1]) ), cv::Point(cvRound(pa.query[i].x),cvRound(pa.query[i].y)), col, 2);

  }

  cv::imwrite("test.jpg", test);
}



