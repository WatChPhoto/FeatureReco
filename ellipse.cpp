
#include "ellipse.hpp"
#include "distance_to_ellipse.hpp"

using std::ostream;

ellipse_st::ellipse_st( double bb, double ee, double phiphi, xypoint cc ) :
  b(bb), e(ee), phi(phiphi), c(cc) { 
  epenalty = 0.0;
  if ( e>1.0 ) { 
    //std::cerr<<"Error eccentricity e="<<e<<" too large, capping at e=1"<<std::endl;
    epenalty = (e-1.0)*(e-1.0)*100;
    e=1.0;
  }
  if ( e<0.0 ) { 
    //std::cerr<<"Error eccentricity e="<<e<<" negative, capping at e=0"<<std::endl;
    epenalty = e*e*100;
    e=0.0;
  }
  if ( bb<0.0 ) bb = fabs( bb );
}

// parametric function radius of ellipse 
// at theta == angle from long axis in radians
double ellipse_st::r_of_theta( double theta ) const {
  double checkfirst = 1 - e*e*cos(theta)*cos(theta);
  if ( checkfirst < 0. ) return 0.;
  return b / sqrt( checkfirst );
}

// get x, y coordinates on ellipse at
// theta == angle from long axis in radians
xypoint ellipse_st::xy( double theta ) const{
  double r = r_of_theta( theta );
  xypoint p = c;
  p.x += r*cos(theta) * cos( phi ) - r*sin(theta) * sin( phi );
  p.y += r*cos(theta) * sin( phi ) + r*sin(theta) * cos( phi );
  return p;
}


// find smallest distance from arbitrary point (px, py) to
// a point on the ellipse.
double ellipse_st::dmin( xypoint p ) const{
  // slow method:  scan points around ellipse


  if (0){ // delete if(0) after testing get_distance
  double curdmin = 999e99;
  for ( double theta=0; theta<2*pi; theta+=2*pi/1000 ){
    xypoint curxy = xy( theta );
    double cur_d = sqrt( (curxy.x-p.x)*(curxy.x-p.x)  +(curxy.y-p.y)*(curxy.y-p.y) );
    if ( cur_d < curdmin ) curdmin = cur_d;
  }
  return curdmin;
  }

  cv::Point q;
  cv::Point center( c.x, c.y );
  cv::Point point( p.x, p.y );
  return get_distance( center, get_a(), b, point, q, phi );
}


// find smallest distance-squared from arbitrary point (px, py) to
// a point on the ellipse.
double ellipse_st::dmin2( xypoint p ) const{
  // slow method:  scan points around ellipse
  if (0){
  double curdmin = 999e99;
  for ( double theta=0; theta<2*pi; theta+=2*pi/1000 ){
    xypoint curxy = xy( theta );
    double cur_d2 = (curxy.x-p.x)*(curxy.x-p.x)  +(curxy.y-p.y)*(curxy.y-p.y) ;
    if ( cur_d2 < curdmin ) curdmin = cur_d2;
  }
  return curdmin;
  }

  double d = dmin( p );
  return d*d;
}


// faster method of finding smallest-squared distnce from point p 
// to a point on the ellipse. From:
// https://www.geometrictools.com/Documentation/DistancePointEllipseEllipsoid.pdf
double ellipse_st::fast_dmin2(xypoint p) const{
  double dmin2;
  bool xflip;
  bool yflip;
  //work in coordinates of ellipse
  xypoint pprime((p.x-c.x)*cos(phi)+(p.y-c.y)*sin(phi), 
		 -1*(p.x-c.x)*sin(phi)+(p.y-c.y)*cos(phi));

  //do calculation as if point is in the first quadrant of the ellipse
  if(pprime.x < 0){
    pprime.x *= -1;
    xflip = true;
  }
  if (pprime.y < 0){
    pprime.y *= -1;
    yflip = true;
  }
  
  double a = b/sqrt(1-e*e);                //long axis of ellipse
  double tbar = Findtbar(a, pprime);       //find roots of F
  double xmin = a*a*pprime.x/(tbar+a*a);   //use tbar to get coordinates
  double ymin = b*b*pprime.y/(tbar+b*b);   //of closest point in ellipse
  //flip back to initial quadrant
  if (xflip) 
    xmin *= -1;
  if (yflip) 
    ymin *= -1;
  
  //calculate minimum distance squared
  dmin2 = (pprime.x-xmin)*(pprime.x-xmin) + (pprime.y-ymin)*(pprime.y-ymin);
  
  return dmin2;

}

// function for which roots are found that give point on ellipse
// with minimum distance to point p
double ellipse_st::Findtbar(double a, xypoint pprime) const{
  
  double t0 = -1*b*b+b*pprime.y; //staring left most point of biseciton
  
  //starting left most point of bisection
  double t1 = -1*b*b+sqrt(a*a*pprime.x*pprime.x + b*b*pprime.y*pprime.y);
  double f = F(a, pprime, t0, t1);
  
  int count =0;

  //do the bisection method to find the roots of F.
  // 0.01 can be reduced for more accurate roots
  while ( sqrt(f*f) > 0.01 ){
    //count++;
    //if(count<20){ 
    //  std::cout<<"t1 = "<<t1<<"\tt0 = "<<t0<<std::endl;
    //  std::cout<<"F = "<<f<<std::endl;
    //}
    
    if(f > 0) 
      t0 =(t1+t0)/2;
    if(f < 0) 
      t1 =(t1+t0)/2;
    f = F(a, pprime, t0, t1);
  }
  double tbar = (t1+t0)/2;
  
  return tbar;
}

double ellipse_st::F(double a, xypoint pprime, double t0, double t1) const{
  double tmp = (t1+t0)/2;

  return (a*pprime.x/(tmp+a*a))*(a*pprime.x/(tmp+a*a)) 
    + (b*pprime.y/(tmp+b*b))*(b*pprime.y/(tmp+b*b)) -1;
}

double ellipse_st::GetRoot(double r0, double z0, double z1, double g) const{
  double n0 = r0*z0;
  double s0 = z1-1;
  double s1 = (g<0 ? 0: sqrt(n0*n0 + z1*z1) -1);
  double s = 0;
  
  while (sqrt(g*g)>0.001){
    s = (s0+s1)/2;
    if(s==s0 || s==s1){break;}
    
    double ratio0 = n0/(s+r0);
    double ratio1 = z1/(s+1);
    g = ratio0*ratio0 + ratio1*ratio1 -1;
    
    if(g>0){
      s0 = s;
    }else if(g<0){
      s1=s;
    }else{
      break;
    }
  }  
  return s;

}

double ellipse_st::DistancePointEllipse2(xypoint p) const{

  bool xflip;
  bool yflip;
  //work in coordinates of ellipse
  xypoint pprime((p.x-c.x)*cos(phi)+(p.y-c.y)*sin(phi), 
		 -1*(p.x-c.x)*sin(phi)+(p.y-c.y)*cos(phi));

  //do calculation as if point is in the first quadrant of the ellipse
  if(pprime.x < 0){
    pprime.x *= -1;
    xflip = true;
  }
  if (pprime.y < 0){
    pprime.y *= -1;
    yflip = true;
  }


  double x0, x1;
  double dmin2;
  double a = b/sqrt(1-e*e);   //long axis of ellipse
  
  if(pprime.y>0){
    if (pprime.x > 0){
      double z0 = pprime.x/a;
      double z1 = pprime.y/b;
      double g = z0*z0 + z1+z1 -1;
      if(g!=0){
	double r0 = (a/b)*(a/b);
	double sbar = GetRoot(r0, z0, z1, g);
	x0 = r0*pprime.x/(sbar+r0);
	x1 = pprime.y/(sbar+1);
      } else{
	x0 = pprime.x;
	x1 = pprime.y;
	dmin2 = 0;
      }
    } else{
      x0 = 0;
      x1 = b;
      dmin2 = (pprime.y-b)*(pprime.y-b);
    }
  }else{
    double number0 = a*pprime.x;
    double denom0 = a*a - b*b;
    if(number0<denom0){
      double xde0 = number0/denom0;
      x0 = a*xde0;
      x1 = b*sqrt(1-xde0*xde0);
      dmin2 = (x0 - pprime.x)*(x0-pprime.x) + (x1-pprime.y)*(x1*pprime.y);
    } else{
      x0 = a;
      x1 = 0;
      dmin2 = (pprime.x-a)*(pprime.x-a);
      
    }
  }
  return dmin2;  
}




// find closest point on ellipse that is smallest distance from 
// arbitrary point (px, py) to a point on the ellipse.
xypoint ellipse_st::closest_point( xypoint p ) const{
  // slow method:  scan points around ellipse
  double curdmin = 999e99;
  xypoint xymin = c;
  for ( double theta=0; theta<2*pi; theta+=2*pi/1000 ){
    xypoint curxy = xy( theta );
    double cur_d = sqrt( (curxy.x-p.x)*(curxy.x-p.x)  +(curxy.y-p.y)*(curxy.y-p.y) );
    if ( cur_d < curdmin ){
      curdmin = cur_d;
      xymin = xy( theta );
    }
  }
  return xymin;
}


ostream& operator<<( ostream& os, const ellipse_st &  e ){

  os << "center (x,y) = (" << e.c.x <<", " << e.c.y << ") "
     << " short radius b = " << e.b
     << " eccentricity = " << e.e
     << " phi = " << e.phi
     <<std::endl;
  return os;
}

