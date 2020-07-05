#include "distance_to_ellipse.hpp"
#include<iostream>

//first slow method
double dmin(cv::Point c, double e0, double e1, cv::Point p , cv::Point& q, double phi){
  const double PI = std::acos(-1);
  //work in coordinates of ellipse
  //changing coordinate system by rotating angle phi to align with major axis.
  double y_0 = (p.x-c.x)*cos(phi)+(p.y-c.y)*sin(phi); 
  double y_1 = -1.0*(p.x-c.x)*sin(phi)+(p.y-c.y)*cos(phi);
  
  //working in the first quadrant
  double y0 = fabs(y_0);
  double y1 = fabs(y_1);
  //std::cout<<y0<<"   "<<y1<<std::endl; 
  
  //keeping record of the sign
  int f_sign = (y_0<0)?-1:1;
  int s_sign = (y_1<0)?-1:1;
   
  //x0, x1 has closest point and distance is the min distance.
  double x0, x1, min_distance = 100000;
  for(double theta=0; theta<=PI/2; theta += PI/10000){
   
    double x_0 = e0*cos(theta);
    double x_1 = e1*sin(theta);
   
    double distance = sqrt(std::pow((x_1-y1),2)+std::pow((x_0-y0),2));
    
    if(distance < min_distance){
      min_distance = distance;
      x0 = x_0; 
      x1 = x_1;
    }
    
  }
  x0 *= f_sign;
  x1 *= s_sign;
  
  q.x = c.x + x0*cos(-phi) + x1*sin(-phi);
  q.y = c.y - x0*sin(-phi) + x1*cos(-phi);

  return min_distance;
}


//second method.
double RobustLength(const double& v1, const double& v2){
  double length;
  if(v1>v2){
    length = fabs(v1)*sqrt(1.0+std::pow((v2/v1),2));
  }

  if(v2>v1){
    length = fabs(v2)*sqrt(1.0+std::pow((v1/v2),2));
  }
  if(v2==v1){
    length = fabs(v2)*sqrt(2.0);
  }

  return length;
}

//based on https://www.geometrictools.com/Documentation/DistancePointEllipseEllipsoid.pdf
double GetRoot(const double& r0, const double& z0, const double& z1, const double& g){
  double n0 = r0*z0;
  double s0 = z1-1.0; //min limit is z1-1
  double s1 = (g<0)?0.0:(RobustLength(n0,z1)-1.0); //Maxlimit
  double s = 0;
  //double a_g=-1000;
  int maxiterations = 1074; //149 for float and 1074 for double
  for(int i=0; i< maxiterations; ++i){
  //  while(fabs(g)>0.01){  
    s=(s0+s1)/2.0;
    if(s==s0 || s==s1){break;}
    double ratio0 = n0/(s+r0);
    double ratio1 = z1/(s+1);
    double g = std::pow(ratio0,2)+std::pow(ratio1,2)-1.0;
    //    a_g = g;
    if (g>0){s0=s;}
    else if (g<0){s1=s;}
     else {break;}
  }
  //std::cout<<"value of g = "<<a_g<<std::endl;
  return s;

}


//c=centre, e0=major axis, e1=minor axis, p = Given point, q = query point, phi = angle made by major axis with +x-axis.
//phi is 0-180 +phi CCW angle and -phi mean CW angle.
//based on https://www.geometrictools.com/Documentation/DistancePointEllipseEllipsoid.pdf
double get_distance(const cv::Point& c, const double& a, const double& b, const cv::Point& p , cv::Point& q, const double& phi){
  
  if(a<=0 || b<=0){std::cerr<<"axis length less than zero"<<std::endl; exit(1);}
  const double PI = std::acos(-1);
  double e0;
  double e1;
  double rot_angle;
  
  if(a>b){
    //if a is bigger axis alighning the coordinate system of ellipse with longer axis of ellipse
    e0 = a;
    e1 = b;
    rot_angle =phi;
  }
  else if(b>a){
    //if b is a bigger axis then alighning the coordinate system of ellipse with the long axis of ellipse.
    e0 = b; //calling longer axis e0.
    e1 = a;
    rot_angle = PI/2.0+phi; //since a makes phi with x-axis, b should make phi+pi/2 degree angle.
  }
  
  //  std::cout<<cos(PI)<<std::endl;
  double x0, x1, distance;
 
  if(a!=b){
    //work in coordinates of ellipse
    //changing coordinate system by rotating angle rot_angle to align with bigger axis.
    double y_0 = (p.x-c.x)*cos(rot_angle)+(p.y-c.y)*sin(rot_angle); 
    double y_1 = -1.0*(p.x-c.x)*sin(rot_angle)+(p.y-c.y)*cos(rot_angle);
    
    //working in the first quadrant
    double y0 = fabs(y_0);
    double y1 = fabs(y_1);
    //std::cout<<y0<<"   "<<y1<<std::endl; 

    //keeping record of the sign
    int f_sign = (y_0<0)?-1:1;
    int s_sign = (y_1<0)?-1:1;
    //std::cout<<"fsign = "<<f_sign<<" s_sign - "<<s_sign<<std::endl;
    
    if(y1>0){
      if(y0>0){
	double z0 = y0/e0;
	double z1 = y1/e1; 
	double g = std::pow(z0,2)+std::pow(z1,2)-1.0;
	if(g!=0){
	  
	  double r0 = std::pow((e0/e1),2);
	  double sbar = GetRoot(r0, z0, z1, g);

	  x0 = f_sign*r0*y0/(sbar+r0);    
	  x1 = s_sign*y1/(sbar +1.0);

	  //converting back to unrotated axis
	  q.x = c.x + x0*cos(-rot_angle) + x1*sin(-rot_angle);
	  q.y = c.y - x0*sin(-rot_angle) + x1*cos(-rot_angle);

	  distance = sqrt((std::pow((x0-y_0),2)+std::pow((x1-y_1),2)));
	}
	else{
	  x0 = f_sign*y0;
	  x1 = s_sign*y1;
	  
	  //converting back to unrotated axis
	  q.x = c.x + x0*cos(-rot_angle) + x1*sin(-rot_angle);
	  q.y = c.y - x0*sin(-rot_angle) + x1*cos(-rot_angle);
	  
	  distance = 0;
	}
      }
      else{
	x0 = f_sign*0;
	x1 = s_sign*e1;

	//converting back to unrotated axis
	q.x =  c.x + x0*cos(-rot_angle) + x1*sin(-rot_angle);
	q.y = c.y - x0*sin(-rot_angle) + x1*cos(-rot_angle);

	distance = fabs(y1-e1);
      }
    }
    else{
      double numer0 = e0*y0;
      double denom0 = std::pow(e0,2)-std::pow(e1,2);
      if(numer0<denom0){
	double xde0 = numer0/denom0;
	x0 = f_sign*e0*xde0;
	x1 = s_sign*e1*sqrt(1.0-std::pow(xde0,2));
	//converting back to unrotated axis
	q.x = c.x + x0*cos(-rot_angle) + x1*sin(-rot_angle);
	q.y = c.y - x0*sin(-rot_angle) + x1*cos(-rot_angle);

	distance = sqrt(std::pow((x0-y_0),2)+std::pow(x1,2));
      }
      else{
	x0 = f_sign*e0;
	x1 = s_sign*0;
	//converting back to unrotated axis
	q.x = c.x + x0*cos(-rot_angle) + x1*sin(-rot_angle);
	q.y = c.y - x0*sin(-rot_angle) + x1*cos(-rot_angle);

	distance = fabs(y0-e0);
      }
    }
  }
 
  //case of circle.
  else {
    double theta = atan2(p.y-c.y,p.x-c.x);
    
    x0 = c.x + a*cos(theta);
    x1 = c.y + a*sin(theta);
      
    q.x = x0;
    q.y = x1;
	
    distance = sqrt(std::pow((p.x-x0),2)+std::pow((p.y-x1),2));
  }
    
  return distance;
}
