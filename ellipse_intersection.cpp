/*################################################################################
#Code to classify if two ellipse intersect, one is inside other, or separated    #
#Only send angle in radians!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    #
#based on:https://www.geometrictools.com/Documentation/IntersectionOfEllipses.pdf#
#                                                                                #  
##################################################################################*/

#include "ellipse_intersection.hpp"
#include <iostream>
const int debug = 0;

using namespace cv;

//based on article :https://www.geometrictools.com/Documentation/IntersectionOfEllipses.pdf
//Sin(x), cos(x) function only use angle in radian in input.
double degree_to_Rad(double a){
  double PI = acos(-1);
  return (a*PI/180.0);
}

int EllipseIntersect::intersect(const ellipse_st e0, const ellipse_st e1){
  //Getting parameters for first ellipse
  double a0 = e0.get_a();
  double b0 = e0.get_b();
  double ee0 = e0.get_e();
  //double phi0 = (ee0==0.0)?0.0:degree_to_Rad(e0.get_phi()); //when input angle is in degree
  double phi0 = (ee0==0.0)?0.0:e0.get_phi();//when angle is in radian
  xypoint cent0 = e0.get_xy();

  //Getting parameters for Second ellipse
  double a1 = e1.get_a();
  double b1 = e1.get_b();
  double ee1 = e1.get_e();
  //double phi1 = (ee1==0.0)?0.0:degree_to_Rad(e1.get_phi());//when angle is in degree;
  double phi1 = (ee1==0.0)?0.0:e1.get_phi();
  xypoint cent1 = e1.get_xy();
  
  //Ellipse equation is (X-K)^TM(X-K)=1
  //X(vector) is variable; K(vector) = centre of ellipse; 
  //M(matrix) = U0U0/l0^2+U0U1^T/l1^2; where U = unit vector along ellipse axis
  //M is positive definite matrix=> symmetric matrix with all +ve eigenvalue.
  Matx21d K0(cent0.x,
  	     cent0.y);
  
  //M = RDR^T (eigenvalue decomposition; R is orthogonal so R^-1 = R^T)
  //R0=matrix with column vector(unit) pointing in direction of ellipse axis.
  Matx22d R0(cos(phi0), -1.0*sin(phi0),
  	     sin(phi0), cos(phi0));
 
  Matx22d D0(1.0/(a0*a0), 0.,
  	     0., 1.0/(b0*b0));
    
    
  //parameters of ellipse 1
  Matx21d K1(cent1.x,
  	     cent1.y);
  
  Matx22d R1(cos(phi1), -1.0*sin(phi1),
  	     sin(phi1), cos(phi1));
  
  Matx22d D1(1.0/(a1*a1), 0.,
                0., 1.0/(b1*b1));

  //D0^-0.5 
  /*Mat D0NegHalf(2, 2, CV_64FC1);
  m.at<double>(0,0)=a0; m.at<double>(0,1)=0.;
  m.at<double>(1,0)=0.; m.at<double>(1,1)=b0;
  */
  Matx22d D0NegHalf(a0, 0,
  		    0, b0);
    
  //D0Half = D0^0.5
  Matx22d D0Half(1.0/a0, 0,
  		 0, 1.0/b0);
    
  //K= centre of transformed ellipse(2)(variable changed to Y)  
  Matx21d K2 = D0Half*(R0.t())*(K1-K0);
  //Compute M2
  Matx22d D0NegHalfR0TR1 = D0NegHalf*(R0.t())*R1;
  
  //M2 = D^-0.5R0^TR1D1R^TR0D0^-0.5
  Matx22d M2 = D0NegHalfR0TR1*D1*(D0NegHalfR0TR1.t());

  //changing to CV_64FC1 to apply eigen()
  Mat Me(2, 2, CV_64FC1);
  Me.at<double>(0,0)=M2(0,0); Me.at<double>(0,1)=M2(0,1);
  Me.at<double>(1,0)=M2(1,0); Me.at<double>(1,1)=M2(1,1);
  
  //M2 = RDR^T
  //M2 is symmetric matrix(when the shape is ellipse).
  Mat eval(2, 2, CV_64FC1);
  //
  //Matx22d eigenvalues;
  Mat evec(2, 2, CV_64FC1);
  //
  eigen(Me, eval, evec);

  //Filling eigenvalue and eigenvector matrix
  std::vector<double> eigenvalues;
  eigenvalues.push_back(eval.at<double>(0,0));  eigenvalues.push_back(eval.at<double>(1,0));
  Matx22d eigenvector(evec.at<double>(0,0), evec.at<double>(0,1),
		      evec.at<double>(1,0), evec.at<double>(1,1));
  
  Matx22d R = eigenvector.t();
  std::vector<double> D = eigenvalues;
  
  //compute K=R^T*K2
  //K is centre of second ellipse after transforming and aligning axis.
  Matx21d K = R.t()*K2;
  
  if(debug){

    std::cout<<"M200 = "<<M2(0,0)<<std::endl;
    std::cout<<"M201 = "<<M2(0,1)<<std::endl;
    std::cout<<"M210 = "<<M2(1,0)<<std::endl;
    std::cout<<"M211 = "<<M2(1,1)<<std::endl;
    std::cout<<"M2"<<std::endl;
    std::cout<<M2(0,0)<<" "<<M2(0,1)<<std::endl;
    std::cout<<M2(1,0)<<" "<<M2(1,1)<<std::endl;
    std::cout<<"eigen vector"<<std::endl;
    std::cout<<R(0,0)<<" "<<R(0,1)<<std::endl;
    std::cout<<R(1,0)<<" "<<R(1,1)<<std::endl;
    std::cout<<"eigenvalue"<<std::endl;
    std::cout<<"D0 "<<" "<<D[0]<<std::endl;
    std::cout<<"D[1] "<<" "<<D[1]<<std::endl;
    //std::cout<<"eigen values "<<" "<<eigenvalues.at<double>(0,0)<<std::endl;
    //std::cout<<"eigen value "<<" "<<eigenvalues.at<double>(0,1)<<std::endl;
    //std::cout<<"eigen values "<<" "<<eigenvalues.at<double>(1,0)<<std::endl;
    //std::cout<<"eigen value "<<" "<<eigenvalues.at<double>(1,1)<<std::endl;

    //Drawing the transformed and axis aligned version.
    Mat mpa = cv::Mat(cv::Size(500, 500), CV_64FC3, Scalar(0,0,0));
    Point center(  250, 250);
    circle(mpa, center, 30,Scalar (0, 255, 0) , 1);
    Size axes1( 30.0/std::sqrt(D[0]), 30.0/std::sqrt(D[1] ));
    Point center1(  30*K(0,0)+250, 30*K(1,0)+250);
    ellipse( mpa, center1, axes1, 0, 0., 360,  Scalar (0, 0, 255), 1 );
    imwrite("transformedandaxisaligned.jpg", mpa);

    //Drawing the transformed but axis not aligned version.
    Mat mpa2 = cv::Mat(cv::Size(500, 500), CV_64FC3, Scalar(0,0,0));
    Point center0(  250, 250);
    circle(mpa2, center0, 30,Scalar (0, 255, 0) , 1);
    Size axes2( 30.0/std::sqrt(D[0]), 30.0/std::sqrt(D[1]) );
    Point center2(  30*K2(0,0)+250, 30*K2(1,0)+250);
    ellipse( mpa2, center2, axes2, e1.get_phi(), 0., 360,  Scalar (0, 0, 255), 1 );
    imwrite("transformedbutaxisnotaligned.jpg", mpa2);
  }
    
 //Transform ellipse0 is Z^T*Z =1 and transformed ellipse 1 is  (Z-K)^T*D*(Z-K)
    double minSqrDistance = std::numeric_limits<double>::max();
    double maxSqrDistance = 0.;
    
    //The special cacse of common centre
    if(K==Matx21d(0,0)){
      for(int i=0; i<2; ++i){
	double invD = 1.0/D[i];
	if(invD < minSqrDistance){
	  minSqrDistance = invD;
	}
	if(invD>maxSqrDistance){
	  maxSqrDistance = invD;
	}
      }
    
      return Classify(minSqrDistance, maxSqrDistance, 0);
    }
    
    //closest 
    double d0p = D[0], d1p = D[1];
    double c0p = K(0,0)*K(0,0), c1p = K(1,0)*K(1,0);

    if(debug){
    std::cout<<"d0= "<<d0<<std::endl;
    std::cout<<"d1= "<<d1<<std::endl;
    }
    
    //sort to make d0>= d1.which allows to bound root of f(s)= d0*k0^2/(d0*s-1)^2+ d1*k1^2/(d1*s-1)^2-1
    std::vector<cv::Point2d>param(2);
    //if axis(0) bigger than axis 1 for transformed ellipse.
    //making bigger axis as major axis.
    if(d0p>=d1){
      param[0] = cv::Point2d(d0p, c0p); 
      param[1] = cv::Point2d(d1p, c1p);
    }
    else{
      param[0] = cv::Point2d(d1p, c1p);
      param[1] = cv::Point2d(d0p, c0p);
    }

    std::vector<cv::Point2d>valid;
    //if d0>d1
    //if major axis bigger than minor.
    if(param[0].x>param[1].x){
      for(int i=0; i<2; ++i){
	//if the coordinate isn't origin.
	if(param[i].y>0){
	  valid.push_back(param[i]);
	}
      }
        }
    //if d0==d1
    else{
      if(param[0].y+param[1].y>0){ //we want atleast one of them to be >0 else the case will be case of same centre which is handled previously.
	valid.push_back(param[0]);
      }
    }
    
    if(valid.size()==1){	
      d0=valid[0].x;
      c0=valid[0].y;
      d1=param[1].x;
      c1=param[1].y;
    }
    else if(valid.size()==2){
      d0=valid[0].x;
      c0=valid[0].y;
      d1=valid[1].x;
      c1=valid[1].y;
    }

    
    size_t numValid = valid.size();
    int numRoots = 0;
    double roots[4];

    if(numValid==2){
      GetRoots( numRoots, roots);
    }
    else if(numValid ==1){
      GetRootsCircle( numRoots, roots);
    }
    
    //numvalid cannot be zero because we already handled K=0 case.
    for(int i=0; i<numRoots; ++i){
      double s = roots[i];
      double p0 = d0*K(0,0)*s/(d0*s-1.0);
      double p1 = d1*K(1,0)*s/(d1*s-1.0);
      double sqrDistance = p0*p0+p1*p1;
      if(sqrDistance < minSqrDistance){
	minSqrDistance = sqrDistance;
      }
      if(sqrDistance > maxSqrDistance ){
	maxSqrDistance = sqrDistance;
      }
    }

    if(debug){
    std::cout<<"maxSqrdist = "<<maxSqrDistance<<std::endl;
    std::cout<<"minSqrdist = "<<minSqrDistance<<std::endl;
    std::cout<<"d0p*c0p+d1p*c1p = "<<d0p*c0p+d1p*c1p<<std::endl;
    std::cout<<"d0 = "<<d0<<std::endl;
    std::cout<<"d1 = "<<d1<<std::endl;
    std::cout<<"c0 = "<<c0<<std::endl;
    std::cout<<"c1 = "<<c1<<std::endl;
    std::cout<<"d0p = "<<d0p<<std::endl;
    std::cout<<"d1p = "<<d1p<<std::endl;
    std::cout<<"c0p = "<<c0p<<std::endl;
    std::cout<<"c1p = "<<c1p<<std::endl;
    }

    return Classify(minSqrDistance, maxSqrDistance, d0*c0+d1*c1);
}

int EllipseIntersect::Classify(double minSqrDistance, double maxSqrDistance, double z){
  if (maxSqrDistance < 1.0)
    {
      return ELLIPSE0_STRICTLY_CONTAINS_ELLIPSE1;
    }
  else if (maxSqrDistance > 1.0)
    {
      if (minSqrDistance < 1.0)
        {
	  return ELLIPSES_OVERLAP;
        }
      else if (minSqrDistance > 1.0)
        {//circle centre outside of ellipse.
	  //z = (h/a)^2+(k/b)^2 for ellipse 1
	  if (z > 1.0)
            {
	      return ELLIPSES_SEPARATED;
            }
	  else
            {
	      return ELLIPSE1_STRICTLY_CONTAINS_ELLIPSE0;
            }
        }
      else  // minSqrDistance = 1
        {//circle center outside ellipse
	  if (z>1.0)//d0c0pd1c1 > 1.0)
            {
	      return ELLIPSE0_OUTSIDE_ELLIPSE1_BUT_TANGENT;
            }
	  else
            {
	      return ELLIPSE1_CONTAINS_ELLIPSE0_BUT_TANGENT;
            }
        }
    }
  else  // maxSqrDistance = 1
    {
      if (minSqrDistance < 1.0)
        {
	  return ELLIPSE0_CONTAINS_ELLIPSE1_BUT_TANGENT;
        }
      else // minSqrDistance = 1
        {
	  return ELLIPSES_EQUAL;
        }
    }
}

//function f(s)
double EllipseIntersect::F(double s){
  double invN0 =1.0/(d0*s-1.0);
  double invN1 = 1.0/(d1*s-1.0);
  double term0 = d0*c0*invN0*invN0;
  double term1 = d1*c1*invN1*invN1;
  double f = term0+term1-1.0;
  return f;
}

/*
//function f'(s)
double EllipseIntersect::DF(double s){
  double invN0 = 1.0/(d0*s-1.0);
  double invN1 = 1.0/(d1*s-1.0);
  double term0 = d0*d0*c0*invN0*invN0*invN0;
  double term1 = d1*d1*c1*invN1*invN1*invN1;
  double df = -2.0*(term0+term1);
  return df;
}
*/
double EllipseIntersect::BisectionRoot(double smin, double smax,double f0, double f1, unsigned int maxIterations, double &root){
  
  root = smin;
  if(smin<smax){
    
    if(f0==0.){
      root = smin;
      return 1;
    }
   
    if(f1==0.){
      root = smax;
      return 1;
    }
    if(f0*f1>0.){
      return 0;
    }
    unsigned int i;
    for(i=2; i<maxIterations; ++i){
      root = 0.5*(smin+smax);
      if(root == smin|| root == smax){
	break;
      }
      double fm= F(root);
    
      double product = fm*f0;
      if(product<0.){
	smax=root;
	f1 = fm;
      }
      else if (product >0.){
	smin=root;
	f0=fm;
      }
      else{
	break;
      }
    }
    return i;
  }
  else{
    return 0;
  }
}

//d0>d1
void EllipseIntersect::GetRoots( int &numRoots, double *roots){
  
  const unsigned int maxIterations = (unsigned int)(3+std::numeric_limits<double>::digits - std::numeric_limits<double>::min_exponent);
  
  //finding smallest root.
  //unsigned int iterations;
  numRoots = 0;
  double smin=0, smax, s;//, fval;
  
  if(debug){
  std::cout<<"d0 = "<<d0<<std::endl;
  std::cout<<"d1 = "<<d1<<std::endl;
  }
  
  //f(0)
  //if f(0)<=0 then root bounding interval is [0,1/d0]
  //else it is [(1-(1+f(0))^0.5)/d1,0]
  double f0 = d0*c0+d1*c1-1.0;//F(0);
   smax = 1.0/d0;

  if(f0>0.0){
    smin = ((1.0-std::sqrt(d0*c0+d1*c1))/d1);//f(smin) should be <0 because the function is increasing in the interval.
    // smax = 0.0;
    //fval = 
    F(smin);
  }
  else if(f0<0.0){
    smin = 0.0;
   
  }

  if(f0!=0){
    //iterations = 
    BisectionRoot(smin, smax,-1,1, maxIterations, s);
    //fval = 
    F(s);
    roots[numRoots++]=s;
  }
  else{
    roots[numRoots++]=0.0;
  }

  double rho = std::pow((d0*d0*c0/(d1*d1*c1)),(1.0/3.0));
  double smid = (1.0+rho)/(d0+rho*d1);
  double fmid = F(smid);
  if(fmid <0.0){
    smin = 1.0/d0;
    smax = smid;
    //iterations = 
    BisectionRoot(smin, smax,1,-1, maxIterations, s);
    roots[numRoots++] = s;

    smin = smid;
    smax = 1.0/d1;
    //iterations = 
    BisectionRoot(smin, smax,-1,1, maxIterations, s);
    roots[numRoots++] = s;
  }
  else if(fmid==0.0){
    roots[numRoots++]=smid;
  }

  //last root
  smin = 1.0/d1;
  smax = (1.0+std::sqrt(d0*c0+d1*c1))/d1;
  //iterations = 
  BisectionRoot(smin, smax,1,-1, maxIterations, s);
  roots[numRoots++]=s;
}

//case when the second transformed ellipse is also circle.
void EllipseIntersect::GetRootsCircle( int &numRoots, double *roots){
  //here c0 = k0^2; c1 = k1^2. This has been done before sending here so just use c0.
  roots[0] = (1.0-std::sqrt(d0*(c0+c1)))/d0;
  roots[1] = (1.0+std::sqrt(d0*(c0+c1)))/d0;
  numRoots =2; 
}
