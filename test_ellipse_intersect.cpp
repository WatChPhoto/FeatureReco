#include "ellipse_intersection.hpp"
#include <string>
#include <fstream>
#include <iostream>

const int debug =0;
using namespace cv; 

//Function used to test performance of ellipse_intersection.cpp class
//The program seperates the cases of operlapping ellipse, separated ellipse etc. in different folders [0 1 2 3 4 5 6 7]
//Instructions:
//mkdir runs
//cd run
//mkdir 0 1 2 3 4 5 6 7 txt
//get out of run and make executable.
//make test_ellipse_intersect
//to make file needs ellipse.cpp, ellipse_intersection.cpp as library
//after the executable is made run ./test_ellipse_intersect
//  

int main(){
  /* Testing how eigen() function works in openCv.
  Matx22d M2(2.0, 3.0,
	     4.0, 1.0);
  std::cout<<"sin(30) = "<<sin(30*3.14/180)<<std::endl;
  std::cout<<M2(0,0)<<M2(0,1)<<std::endl;
  std::cout<<M2(1,0)<<M2(1,1)<<std::endl;
  std::vector<double> eigenvalues;
  Matx22d eigenvectors;
  eigenNonSymmetric(M2,eigenvalues, eigenvectors);

  //long double a[2][2]={2,0,0,2};
  //long double b[2][2]={3,0,0,3};
  //long double c[2][2]=a*b;
  
  std::cout<<eigenvalues[0]<<" "<<eigenvalues[1]<<" size "<<eigenvalues.size()<<std::endl;
  std::cout<<"===================="<<std::endl;
  std::cout<<eigenvectors(0,0)<<eigenvectors(0,1)<<std::endl;
  std::cout<<eigenvectors(1,0)<<eigenvectors(1,1)<<std::endl;
 */
  srand (time(NULL));
  
  for(int i=0; i<1; i++){
  Mat m = cv::Mat(cv::Size(500, 500), CV_64FC3, Scalar(0,0,0));
  
  
  double phi0 =(rand() %180);
  double b0=(rand()%60+10);
  double e0=(rand()%10)/10.0;
  double a0=b0/std::sqrt(1-e0*e0);
  
  double x0=(rand()%150+90);
  double y0=(rand()%150+90);
  Size axes(  a0, b0 );
  Point center( x0, y0 );

  ellipse( m, center, axes, phi0, 0., 360,  Scalar (0, 255, 0), 1 );

  double phi1 =(rand() %180);
  double b1=(rand()%60+10);
  double e1=(rand()%10)/10.0;
  double a1=b1/std::sqrt(1-e1*e1);
  double x1=(rand()%150+90);
  double y1=(rand()%150+90);
  
  std::string file_name = "./runs/txt/ellipse_test"+std::to_string(i)+".txt";
  std::ofstream outfile;
  outfile.open(file_name);
  

  if(debug){
    std::cout<<"x0 = "<<x0<<std::endl;
    std::cout<<"y0 = "<<y0<<std::endl;
    std::cout<<"a0 = "<<a0<<std::endl;
    std::cout<<"b0 = "<<b0<<std::endl;
    std::cout<<"e0 = "<<e0<<std::endl;
    std::cout<<"phi0 = "<<phi0<<std::endl;

    std::cout<<"e1 = "<<e1<<std::endl;
    std::cout<<"a1 = "<<a1<<std::endl;
    std::cout<<"b1 = "<<b1<<std::endl;
    std::cout<<"phi1 = "<<phi1<<std::endl;
    std::cout<<"x1 = "<<x1<<std::endl;
    std::cout<<"y1 = "<<y1<<std::endl;
  
    outfile<<"x0 = "<<x0<<std::endl;
    outfile<<"y0 = "<<y0<<std::endl;
    outfile<<"a0 = "<<a0<<std::endl;
    outfile<<"b0 = "<<b0<<std::endl;
    outfile<<"e0 = "<<e0<<std::endl;
    outfile<<"phi0 = "<<phi0<<std::endl;

    outfile<<"e1 = "<<e1<<std::endl;
    outfile<<"a1 = "<<a1<<std::endl;
    outfile<<"b1 = "<<b1<<std::endl;
    outfile<<"phi1 = "<<phi1<<std::endl;
    outfile<<"x1 = "<<x1<<std::endl;
    outfile<<"y1 = "<<y1<<std::endl;
   
  }
  
  Size axes1(  a1, b1 );
  Point center1( x1, y1 );
  
  ellipse( m, center1, axes1, phi1, 0., 360,  Scalar (0, 0, 255), 1 );

  //converting angle to radian
  phi0=phi0*acos(-1)/180.0;
  phi1=phi1*acos(-1)/180.0;
  
  ellipse_st first(b0, e0, phi0, xypoint(x0,y0));
  ellipse_st second(b1, e1, phi1, xypoint(x1,y1));
  
  
  EllipseIntersect E;
  int inters = E.intersect(first, second);
  
  std::string text;
  switch(inters){
  case 0: text = "Ellipse Separated"; break;
  case 1: text="Ellipses Overlap"; break;
  case 2: text="Ellipse0(Green) outside Ellipse1(Red) But Tangent"; break;
  case 3: text="Ellipse0(Green) strictly contains Ellipse1(Red)"; break;
  case 4: text="Ellipse0(Green) strictly contains Ellipse1(Red) but tangent"; break;
  case 5: text="Ellipse1(Red) strictly contains Ellipse0(Green)"; break;
  case 6: text="Ellipse1(Red) contains Ellipse0(Green) but tangent"; break;
  case 7: text="Ellipses Equal"; break;
  }
  

  putText(m, "Green= Ellipse 0; Red = Ellipse 1", cv::Point(2,400),FONT_HERSHEY_SIMPLEX , 0.5, cv::Scalar(0,0,255) );
  
  putText(m, text, cv::Point(2,480),FONT_HERSHEY_SIMPLEX , 0.5, cv::Scalar(0,0,255) );

  putText(m, std::to_string(inters), cv::Point(2,500),FONT_HERSHEY_SIMPLEX , 0.5, cv::Scalar(0,0,255) );

  if(debug){
   std::string name = "./runs/"+std::to_string(inters)+"/ellipse_test"+std::to_string(i)+".jpg";
    imwrite(name, m);
  }
  if(!debug){
  imwrite("ellipse.jpg", m);
  }
  outfile.close();
   }
  return 0;
}
