#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <fstream>                                                           
#include<TF1.h> 
#include <cmath>
#include <iostream>
#include <vector>
#include <opencv2/calib3d.hpp>
#include <TH2D.h>
#include <TH3D.h>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TFile.h>
#include "ImageDataReader.hpp"
#include "distancelib.hpp"


#include <TMath.h>
#include <TGraph2D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TF2.h>
#include <TH1.h>
#include <Math/Functor.h>
#include <TPolyLine3D.h>
#include <Math/Vector3D.h>
#include <Fit/Fitter.h>

#include <cassert>
#include <TLine.h>

TH1D* h = new TH1D("h","radius Vs error;rad;error",1000,0,0);//ceil(lower_r-10),ceil(upper_r+10));
TH1D* h_min_dist = new TH1D("h_min_dist","Distance to closest;distance;counts/bin",1000,0,0);//ceil(lower_r-10),ceil(upper_r+10));
struct minimize{
  unsigned img_no;
  std::vector<cv::Point2f> ellipses;
  std::vector<Ellipse> ellipses1;
  std::vector<WorldPoints> all_pmts;
  cv::Matx33d camera_matrix;
  cv::Matx14d dist_coeffs; //distortion coefficient
  cv::Matx33d rmat; 

  //for image 239
  double yaw0;
  double pitch0;
  double roll0;

  // drone reading for this image
  double yaw1 ;
  double pitch1 ;
  double roll1 ;  
  
  // drone depth reading for this image
  double z1 ;
  
  cv ::Matx33d rot_AB; // rotation matrix for image 239
  
  cv::Matx33d R;  // rotation matrix for current facing guess (filled by calculate_rotation_matrix)
  cv::Matx31d rvec1; // "
  
  cv::Matx31d dir; //
  double theta; //
  
  std::vector<cv::Point2f>im_points;
  
  minimize( const std::string & filename, const std::string & survey_id ){
    
    std::cout<<"filename = "<<filename<<std::endl;
    float offset = 250; //Offset because images were cropped 250px, but image from projected points is not cropped.
    double PI = std::acos(-1);
    
    //##########################Camera intrinsic parameters.#################################
    double fx=2.760529621789217e+03; //focal length converted to pixel unit in x direction.
    double fy=2.767014510543478e+03; //pixel/cm is different in x and y axis so fx and fy are different.
    double cx=1.914303537872458e+03; //x coordinate of a point in image where light coming from principal axis would appear. 
    double cy=1.596386868474348e+03; //cx~40,000/2 and cy~30,000/2 since image is "40,000 by 30,000" size.

    dist_coeffs = cv::Matx14d(-0.2398, 0.1145, 0, 0); //just four par. radial & tangential.x

    camera_matrix  = cv::Matx33d(fx, 0, cx,
				 0, fy, cy,
				 0,  0,  1);
    //##########################################################################################
    
    //Extinsic parameters for manually labelled image(BarrelSurveyFar, 239)
    cv::Matx31d rvec (1.48220885,-0.85878027,0.75081903);
    cv::Matx31d tvec(-77.19825549,1313.21348846,-794.66238321); 
    cv::Rodrigues(rvec,rmat);

    //Reading image metadata(survey_id, yaw, pitch, roll, temp,...)
    ImageDataReader & idr = ImageDataReader::GetInstance();
    img_no =get_img_num(filename);
    ImageMetaData m = idr.GetMetaData(survey_id,img_no );
    ImageMetaData i239 = idr.GetMetaData("BarrelSurveyFar", 239);

    //inputfile info
    std::cout<<"image "<<filename<<" info"<<std::endl;
    std::cout<<m.face.yaw<<std::endl;
    std::cout<<m.face.pitch<<std::endl;
    std::cout<<m.face.roll<<std::endl;
  
    //manually labelled image info->239BarrelFar was manually labelled to get rvec and tvec
    std::cout<<"image 239 info :"<<std::endl;
    std::cout<<i239.face.yaw<<std::endl;
    std::cout<<i239.face.pitch<<std::endl;
    std::cout<<i239.face.roll<<std::endl;
    std::cout<<"he_pmts"+filename+".txt"<<std::endl;

    ellipses1 = read_ellipses_in_image("he_pmts"+filename+".txt");

    for(int i=0; i<ellipses1.size();i++){
      ellipses.push_back(cv::Point2f(ellipses1[i].x,ellipses1[i].y));
    }


    //for image 239
    yaw0 = (i239.face.yaw)*PI/180.;
    pitch0 = (i239.face.pitch)*PI/180.;
    roll0 = (i239.face.roll)*PI/180.;

    //this is for any image we want to work with.
    //this need to be read for particular image.
    yaw1 = (m.face.yaw)*PI/180.;
    pitch1 = (m.face.pitch)*PI/180.;
    roll1 = m.face.roll*PI/180.;  



    //correcting depth
    double b = 1.02585; //ratio of density of sea water to pure water
    double c = -0.271196;//offset
    z1 = m.depth;
    z1 = c+b*z1;// corrected z-coordinate
    z1 *=100; //converting to cm
    
    std::cout<<"z depth is = "<<z1<<std::endl;


    //rotation about x,y,z axis of camera.
    double thx0=pitch0;
    double thy0=-yaw0;
    double thz0=-roll0;
    
    //given a vector in rotated frame about x-axis, this matrix finds coordinate of vector in unrotated frame.
    cv::Matx33d ABx (1.0,   0.0,      0.0,
		     0.0, cos(thx0), -sin(thx0),
		     0.0, sin(thx0),  cos(thx0));
    cv::Matx33d ABy (cos(thy0), 0.0, sin(thy0),
		     0.0,      1.0,    0.0,
		     -sin(thy0), 0.0, cos(thy0));
    cv::Matx33d ABz (cos(thz0), -sin(thz0), 0.0,
		     sin(thz0),  cos(thz0), 0.0,
		     0.0,        0.0,     1.0 );
    
    rot_AB =  ABz.t()*ABx.t()*ABy.t();    

    calculate_rotation_matrix( yaw1, pitch1, roll1 );
    all_pmts = read_all_world_points();
    
  }

  void calculate_rotation_matrix( double yaw, double pitch, double roll  ){
    double thx1=pitch;
    double thy1=-yaw;
    double thz1=-roll;
    
    cv::Matx33d ACx (1.0,   0.0,      0.0,
		     0.0, cos(thx1), -sin(thx1),
		     0.0, sin(thx1),  cos(thx1));
    cv::Matx33d ACy (cos(thy1), 0.0, sin(thy1),
		     0.0,      1.0,    0.0,
		     -sin(thy1), 0.0, cos(thy1));
    cv::Matx33d ACz (cos(thz1), -sin(thz1), 0.0,
		     sin(thz1),  cos(thz1), 0.0,
		     0.0,        0.0,     1.0 );
    
    cv::Matx33d rot_AC = ACz.t()*ACx.t()*ACy.t();
    //cv::Matx33d rot_AB =  ABy.t()*ABx.t()*ABz.t();
    //cv::Matx33d rot_AC = ACy.t()*ACx.t()*ACz.t();
    
    cv::Matx33d rot_BC = rot_AC*(rot_AB.t()); //given the coordinate in frame B, finds the coordinate in frame C.
    R = rot_BC*rmat;//rmat brings coordinate from world to 045 frame(B frame).x
    // R_o = R;
    cv::Rodrigues(R,rvec1);

    /*
    std::cout<<img_no<<" rotation matrix"<<std::endl;
    std::cout<<R(0,0)<<"\t"<<R(0,1)<<"\t"<<R(0,2)<<std::endl;
    std::cout<<R(1,0)<<"\t"<<R(1,1)<<"\t"<<R(1,2)<<std::endl;
    std::cout<<R(2,0)<<"\t"<<R(2,1)<<"\t"<<R(2,2)<<std::endl;
    
    std::cout<<"Rotation vector"<<std::endl;
    std::cout<<rvec1(0,0)<<"\t"<<rvec1(0,1)<<"\t"<<rvec1(0,2)<<std::endl;
    */
    /*
      for(int i=0; i<v.size(); i++){
      cv::Matx31d camera_coords = v[i].get_camera_coords(rmat, tvec);
      }
    */
    dir =R.t()*cv::Matx31d(0,0,1); //direction of z-axis of camera.
    theta = std::atan2(dir(1,0),dir(0,0))*180./PI; // angle with x-axis.
    theta=(theta>=0)?theta:(360+theta);
    theta=theta*PI/180.;

    //std::cout<<"theta ="<<theta<<std::endl;

  }


  void transform(const double* par){
    //p[] = {x, y, z, yaw, pitch, roll}
    // Clear im_points and refill it for current position and facing
    im_points.clear();
   

    cv::Matx31d tvec1;

    double x = par[0]*cos(par[1]);
    double y = par[0]*sin(par[1]);
    //double x=par[0];
    //double y= par[1];
    double z = par[2];
    double yaw = par[3];
    double pitch = par[4];
    double roll = par[5];

    calculate_rotation_matrix( yaw, pitch, roll );

    cv::Matx31d p(x, y, z);
    cv::Matx31d tv = -R*p;
    std::vector<cv::Point3f> object_in_view;
    std::vector<std::string> id;
    for(int i=0; i<all_pmts.size();i++){
      cv::Matx31d xp = R*cv::Matx31d(all_pmts[i].x,all_pmts[i].y,all_pmts[i].z)+tv;
      
      //Clipping (Removing points behind camera)
      if(xp(2,0)>0){
	object_in_view.push_back(cv::Point3f(all_pmts[i].x,all_pmts[i].y,all_pmts[i].z));
	    
	std::string str = all_pmts[i].id.substr(0,5);
	if(isalpha(str[0])){
	  id.push_back(str);
	}
	else{
	  int num;
	  std::stringstream ss;  
	  ss << str;  
	  ss >> num;
	  id.push_back(std::to_string(num));
	}
      }
	  
    }
    std::vector<cv::Point2f>temp_pt;
    if(object_in_view.size()>0){   
      projectPoints(object_in_view,rvec1,tv,camera_matrix,dist_coeffs,temp_pt);
    }
    
    //clip points if they are outside of image border.
    //fills im_points with clipped points(Points in boundary of image, boundary is second argument).
    clip(temp_pt,cv::Size(40000,30000));
    
  }



  double operator()(const double* par){
    //par[] = {x, y, z, yaw, pitch, roll}
    h_min_dist->Reset();
    transform(par);
    float offset=250; //since images were cropped.
    double sum=0;
    for(int i=0; i<ellipses.size();i++){
      double x0 = ellipses[i].x;
      double y0 = ellipses[i].y;
      double min_dist = 1e+8;
      for(int j=0;j<im_points.size();j++){
	double x1 = im_points[j].x;
	double y1 = im_points[j].y-offset;
	double d = (x1-x0)*(x1-x0)+(y1-y0)*(y1-y0);
	if(d<min_dist){min_dist=d;}
      }
      //      std::cout<<"min_dist = "<<min_dist<<std::endl;
      min_dist/=100;
      h_min_dist->Fill(min_dist);
      sum +=min_dist;
    }

    //double yaw = par[4]; // in degrees???
    // double yaw_err = 5.0; // degrees
    //dsq += ( yaw1 - yaw ) * ( yaw1 - yaw ) / yaw_err / yaw_err; 
    // std::cout<<"Error is = "<<dsq<<std::endl;
    //sum/=35;
    h->Fill(sum);
    return sum;
  }

private:
  void clip(std::vector<cv::Point2f> pts, const cv::Size size){
    for(const cv::Point2f pt:pts){
      if(is_inside(pt,size)){
	im_points.push_back(pt);
      }
    }
  }
  
  bool is_inside(const cv::Point2f & pt,const cv::Size & s){
    // determine if the point is in the boundary of image
    return (pt.x>0 && pt.x<s.width)&&(pt.y>0 && pt.y<s.height);
  }
  
};

void draw_error(std::vector<Ellipse> ellipses, std::vector<cv::Point2f> im_points,float offset, cv::Mat& m){
 
  for(int i=0; i<ellipses.size();i++){
    double x0 = ellipses[i].x;
    double y0 = ellipses[i].y;
    double lmin = 1e+8;
    cv::Point2f p2;
    for(int j=0;j<im_points.size();j++){
      double x1 = im_points[j].x;
      double y1 = im_points[j].y-offset;
      double d = (x1-x0)*(x1-x0)+(y1-y0)*(y1-y0);
      if(d<lmin){lmin=d; p2=cv::Point2f(x1,y1);}
    }
    //draw line from ellipses to points
    cv::line(m,cv::Point2f(x0,y0),p2,cv::Scalar(0,230,0),2);
    cv::circle(m,cv::Point2f(x0,y0),10,cv::Scalar(255, 102, 255),-1);
    cv::Size axes(  int(ellipses[i].get_a()), int(ellipses[i].b) );
    cv::ellipse(m,cv::Point2f(x0,y0),axes,ellipses[i].phi*PI/180.,0,360,cv::Scalar (255, 102, 255),2);
  }
}

void distance_range(double& low, double& high, std::vector<Ellipse> e){
  low =1e+8;
  high=0;
  for(Ellipse ellipses: e){
    double a = ellipses.get_a();
    if(a<low){low=a;}
    if(a>high){high=a;}
  }
  low = 1700.-(86190./low);
  high = 20+1700.-(86190./high);
}

void make_error_histogram(std::vector<cv::Point3f>ev, int z){
  std::string lbl = std::to_string(z);
  TH2D *error = new TH2D(("error for z = "+lbl).c_str(),"Reprojection error;x;y",1500.,0.,1500.,1500.,0.,1500.);
  for(cv::Point3f p:ev){
    error->Fill(p.x,p.y,p.z);
  }
}

void make_2D_histogram(ROOT::Math::Functor& fcn,unsigned idx1, unsigned idx2, const double* pv, float p1_min, float p1_max,float p2_min, float p2_max,const std::string& name){

  double p[6];
  for(int i=0;i<6;i++){
    p[i]=pv[i];
  }
  
  TH2D *er = new TH2D(name.c_str(),"${name};x;y",50,0.,0,50.,0.,0.);
  
  Int_t m = 100;
  Int_t n = 100;
  for (Int_t i=0;i<m;i++) {
    p[idx1]=p1_min+i*((p1_max-p1_min)/m);
    float x = p[idx1];
    for(Int_t j=0;j<n;j++){
      p[idx2]=p2_min+j*((p2_max-p2_min)/n);
      float y = p[idx2];
      double z = fcn(p);
      er->Fill(x,y,z);
    }
  }
}

void make_graph(ROOT::Math::Functor& fcn,unsigned idx, const double* pv, float p_min, float p_max,const std::string& name){
  double p[6];
  for(int i=0;i<6;i++){
    p[i]=pv[i];
  }

  Int_t n = 200;
  Double_t x[n], y[n];
  for (Int_t i=0;i<n;i++) {
    p[idx]=p_min+i*((p_max-p_min)/n);
    x[i] = p[idx];
    y[i] = fcn(p);
    //pv[idx]+=1;
  }
  
  TGraph *gr  = new TGraph(n,x,y);    
  gr->SetName(name.c_str());
  
  TLine *lin = new TLine(pv[idx],gr->GetMinimum(),pv[idx],gr->GetMaximum());
  lin->SetLineWidth(8);
  lin->SetLineColor(2);
  gr->GetListOfFunctions()->Add(lin);
  
  gr->Write();
  //graphends

}



int main(int argc, char **argv){
  
  std::string filename = std::string(argv[argc-1]);
  std::string survey_id = std::string(argv[argc-2]);

  minimize f2min( filename, survey_id );


  ROOT::Fit::Fitter fitter;
  ROOT::Math::Functor fcn( f2min, 6 );
  
  double lower_r, upper_r; 
  distance_range(lower_r,upper_r, f2min.ellipses1 );

  TFile* fout = new TFile("testingdistance.root","recreate");
  
  

  double theta = f2min.theta; // theta in rad or deg????
  double r = (lower_r+upper_r)/2.;
  std::cout<<"lower rad = "<<lower_r<<" higher r = "<<upper_r<<std::endl;
  std::cout<<" radius *** = "<<r<<std::endl;
  std::cout<<"b4 pstart"<<std::endl;
  double pStart[6] = { r,  // make sure this is in radians!!!
		       f2min.theta,
		       f2min.z1,
		       f2min.yaw1,
		       f2min.pitch1,
		       f2min.roll1  };
  


  std::cout<<"after pstart"<<std::endl;
  double pRange[6][2] = {{lower_r-100,upper_r+100},
			 {theta-(15*PI/180.),theta+(15*PI/180.)},
			 {f2min.z1-10,f2min.z1+10},
			 {f2min.yaw1-(15*PI/180.),f2min.yaw1+(15*PI/180.)},
			 {f2min.pitch1-(10*PI/180.),f2min.pitch1+(10*PI/180.)},
			 {f2min.roll1-(3*PI/180.),f2min.roll1+(3*PI/180.)}};
  //just for drawing graph

  std::string pName[6]={"r","theta","z","Yaw","Pitch","Roll"};

  for(unsigned i=0;i<6;i++){
    std::cout<<pName[i]<<" "<<pStart[i]<<" ("<<pRange[i][0]<<" , "<<pRange[i][1]<<")"<<std::endl;
  }
  
  double pv[6];
  for(int i=0;i<6;i++){
    pv[i]=pStart[i];
  }
  for(unsigned i=0;i<6;i++){
    make_graph(fcn,i,pv,pRange[i][0],pRange[i][1],pName[i]);
  }


 // std::cout<<"after prange"<<std::endl;
  //r, theta, z, yaw, pitch, roll
  double pSteps[6];// = {0.1,0.001,0.1,0.001,0.001,0.001};
  for(unsigned i=0;i<6;i++){
    std::cout<<"lower for "<< i<<" is = "<<pRange[i][0]<<std::endl;
    std::cout<<"upper for "<< i<<" is = "<<pRange[i][1]<<std::endl;
    pSteps[i]=0.01*(pRange[i][1]-pRange[i][0]);
    std::cout<<"Psteps is = "<<pSteps[i]<<std::endl;
  }
  
  std::cout<<"affter psteps"<<std::endl;
  //Setting the serarch region for each parameters
  //  
 
  
  
  std::cout<<"after setting range"<<std::endl;
  fitter.SetFCN( fcn, pStart);
  
  
  //for(int i=0; i<6;++i)fitter.Config().ParSettings(i).SetLimits(pRange[i][0],pRange[i][1]);
  for(int i=0; i<6;++i)fitter.Config().ParSettings(i).SetName(pName[i]);
  
  std::cout<<"no of parameters = "<<fitter.Config().NPar()<<std::endl;
  std::cout<<"after setting range"<<std::endl;
  //set minimizer
  std::cout<<"Default max Iterations = "<<fitter.Config().MinimizerOptions().MaxIterations()<<std::endl;
  //  fitter.Config().MinimizerOptions().SetMaxIterations(10000); 
  fitter.Config().MinimizerOptions().SetPrintLevel(10); 
  
  fitter.Config().SetMinimizer("Minuit","Simplex");
  //std::cout<<"minimizer = "<<fitter.GetMinimizer()<<std::endl;//"Minuit");
  std::cout<<"after setting minimizer"<<std::endl;
  // set step sizes different than default ones (0.3 times parameter values)
  for(int i=0; i<6;++i)fitter.Config().ParSettings(i).SetStepSize(0.01);//pSteps[i]);
  fitter.Config().MinimizerOptions().Print(); 
  fitter.Config().MinimizerOptions().SetErrorDef(10); 
  fitter.Config().MinimizerOptions().SetStrategy(2); 
  std::cout<<"before ok"<<std::endl;
  

  bool ok = fitter.FitFCN();
  if (!ok) {
      Error("distancetest","Fit failed");
      return 1;
  }

  std::cout<<"After ok"<<std::endl;

  const ROOT::Fit::FitResult & result = fitter.Result();
  
  
  for(unsigned i=0;i<6;i++){
    make_graph(fcn,i,result.GetParams(),pRange[i][0],pRange[i][1],pName[i]+"final");
    
  }
  
  {
    unsigned i=1;//theta
    unsigned j=3;//yaw
    make_2D_histogram(fcn,i,j, result.GetParams(), pRange[i][0], pRange[i][1],pRange[j][0], pRange[j][1],pName[i]+"vs"+pName[j]);
  }

  f2min(result.GetParams());
  
  
  
  std::cout << "Total final distance square " << result.MinFcnValue() << std::endl;
  result.Print(std::cout);
 

  h_min_dist->Write();
  fout->Write();
  h->Write();
  fout->Close(); 
 


  if(0){
    /*
  double minx=0;
  double miny=0;
  double minz=0;
  std::vector<std::string> final_id;
  std::vector<cv::Point3f> obv;
  std::vector<cv::Point2f> imp;
  
  double range =20;//25;
  double th1 = ((theta-range)>=0)?(theta-range):(360+theta-range);
  double th2 = ((theta+range)>=0)?(theta+range):(360+theta+range);
  th2 = (th2>360)?(th2-360.):th2;
  double high,low;
  //(th2-range>=0)?(th2-range):(360+th2-range)
  if((th2-th1)>0){high=th2; low=th1;}

  else{high=th1;low=th2;}
  std::cout<<"low = "<<low<<std::endl;
  std::cout<<"high = "<<high<<std::endl;
  std::cout<<"arctan of 1,1 is "<<std::atan2(1,1)*180./PI<<std::endl;
  std::cout<<"arctan of -1,1 is"<<std::atan2(1,-1)*180.0/PI<<std::endl;
  std::cout<<"arctan of -1,-1 is"<<std::atan2(-1,-1)*180./PI<<std::endl;
  std::cout<<"arctan of 1,-1 is "<<std::atan2(-1,1)*180./PI<<std::endl;
  
  double lower_r, upper_r;
  distance_range(lower_r,upper_r,ellipses1);
  std::cout<<"lower_r = "<<lower_r<<"\t"<<"upper_r = "<<upper_r<<std::endl;
  double low_r = 1700.-(86190./lower_r);
  double high_r = 20+1700.-(86190./upper_r);
  std::cout<<"low_r = "<<low_r<<"\t"<<"high_r = "<<high_r<<std::endl;

  
  TFile* fout = new TFile("testingdistance.root","recreate");
  //TH3D *error =   new TH3D("reprojection error","Reprojection error;x;y;z",int((h.x-l.x+20)/5),l.x-10,h.x+10,int((h.y-l.y+20)/5),l.y-10,h.y+10,100.,z1-50,z1+50.);

  std::cout<<"z1 = "<<z1<<std::endl;
  
  double ang;
  double rad;
  for(int z=z1-3;z<z1+3; z+=1){	
    //vector to hold value of x,y and the error value.
    std::vector<cv::Point3f>ev; //error value
    //     for(int x=l.x; x<h.x; x+=5){
    //  for(int y=l.y; y<h.y; y+=5){
    for(float r=low_r; r<high_r; r+=1){
      for(float theta=int(low); theta<high; theta+=0.2){
	double x = r*cos(theta*PI/180.);
	double y = r*sin(theta*PI/180.);
	/*
	    std::cout<<"z111 depth is = "<<z<<std::endl;
	  float r=std::sqrt(x*x+y*y);//+z*z);
	double alfa = std::atan2(y,x)*180./PI; //angle x,y make.
	alfa = (alfa>=0)?alfa:(360+alfa);
	if(r<low_r||r>high_r){continue;} //rad of tank=1690cm
	else if((high-low)>(2*range)){ if(alfa>=low && alfa<=high){continue;}}
	else{if(alfa<low||alfa>high){continue;}}
	*/
	
	  
	  //std::cout<<"alfa = "<<alfa<<std::endl;
	  //if(alfa<low ||alfa>high){continue;}
	  //       std::cout<<"passed few the hurdles"<<std::endl;
    /*	cv::Matx31d p(x, y, z);
	cv::Matx31d tv = -R*p;
	std::vector<cv::Point3f> object_in_view;
	std::vector<std::string> id;
	for(int i=0; i<all_pmts.size();i++){
	  cv::Matx31d xp = R*cv::Matx31d(all_pmts[i].x,all_pmts[i].y,all_pmts[i].z)+tv;
	  if(xp(2,0)>0){
	    object_in_view.push_back(cv::Point3f(all_pmts[i].x,all_pmts[i].y,all_pmts[i].z));
	    

	    std::string str = all_pmts[i].id.substr(0,5);
	    if(isalpha(str[0])){
	      id.push_back(str);
	    }
	    else{
	      int num;
	      std::stringstream ss;  
	      ss << str;  
	      ss >> num;
	      id.push_back(std::to_string(num));
	    }
	  }
	  
	}
	if(object_in_view.size()>0){
	  std::vector<cv::Point2f>im_points;
	  projectPoints(object_in_view,rvec1,tv,camera_matrix,dist_coeffs,im_points);
	  
	  double rep_error = reprojection_error(ellipses, im_points,offset);
	  ev.push_back(cv::Point3f(x,y,rep_error)); 
	  //error->Fill(x,y,z,rep_error);
	  //	  std::cout<<"Repo error = "<<rep_error<<std::endl;
	  if(rep_error<min_rep_err){
	    //obv = object_in_view;
	    final_id = id;
	    imp=im_points;
	    min_rep_err = rep_error;
	    minx=x;
	    miny=y;
	    minz=z;
	    tvec1=tv;
	    rad=r;
	    ang=theta;
	  } 
	}
      }
    }
    make_error_histogram(ev,z);
  }
  
  //writing rotation and translation vector to file
  std::ofstream file(filename+"_rot_trans.txt");
  file<<"rot_vec = ["<<rvec1(0,0)<<'\t'<<rvec1(0,1)<<'\t'<<rvec1(0,2)<<"]"<<std::endl;
  file<<"tr_vec = ["<<tvec1(0,0)<<'\t'<<tvec1(0,1)<<'\t'<<tvec1(0,2)<<"]"<<std::endl;
  file.close();
  
  std::cout<<"position is ("<<minx<<" , "<<miny<<" , "<<minz<<" )"<<std::endl;
  //ps=cv::Matx31d(minx,miny,minz);
  
  std::cout<<"rad = "<<rad<<"  ang = "<<ang<<std::endl;
  cv::Matx31d a = R.t()*ACz.t()*ACx.t()*cv::Matx31d(0,0,1);
  std::cout<<"angle of direction is ="<<std::atan2(a(1,0),a(0,0))*180./PI;
  
  /*    //20001: UK B1
    //30001: Korean B1
    std::vector<std::string>label{"20001","UKB2","UKB3","UKB4","UKB5","30001","KB2","KB3","KB4","KB5"};
    std::vector<cv::Point2f>im_points;
    std::vector<cv::Point3f>light_injectors{cv::Point3f(1490.73,768.14,1302.95),cv::Point3f(1490.73,768.14,666.65),cv::Point3f(1490.73,768.14,30.35),cv::Point3f(1490.73,768.14,-676.65),cv::Point3f(1490.73,768.14,-1312.95),cv::Point3f(1490.73,768.14,1232.25),cv::Point3f(1490.73,768.14,595.95),cv::Point3f(1490.73,768.14,-40.35),cv::Point3f(1490.73,768.14,-605.95),cv::Point3f(1490.73,768.14,-1242.25)};
    //li=light_injectors;
    projectPoints(light_injectors,rvec1,tvec1,camera_matrix,dist_coeffs,im_points);
 
    for(int i=0;i<im_points.size();i++){
      int x=im_points[i].x;
      int y=im_points[i].y;
      if(x>0 && x<4000 && y>0 && y-offset<2750){
	cv::circle( img, cv::Point( x, y-offset ), 10, cv::Scalar(0,255,250), -1 );
	cv::putText(img, label[i],cv::Point( x, y-offset ) , cv::FONT_HERSHEY_PLAIN,4, cv::Scalar(0,255,250),3);
      }
    }
   */   
  /*for(int i=0;i<imp.size();i++){
    
    if(imp[i].x>0 && imp[i].x<4000 && imp[i].y>0 && imp[i].y-offset<2700){
      std::string text = final_id[i];
      if(isdigit(text[0])){
	cv::circle( img, cv::Point( imp[i].x, imp[i].y-offset ), 10, cv::Scalar(0,255,250), -1 );
	cv::putText(img, text,cv::Point( imp[i].x, imp[i].y-offset ) , cv::FONT_HERSHEY_PLAIN,4, cv::Scalar(0,255,250),3);
      }
      else{
	cv::circle( img, cv::Point( imp[i].x, imp[i].y-offset ), 10, cv::Scalar(0,255,250), -1 );
	cv::putText(img, text,cv::Point( imp[i].x, imp[i].y-offset ) , cv::FONT_HERSHEY_PLAIN,4, cv::Scalar(0,255,250),3);
      }
    }
  }
  std::string name = std::to_string(img_no)+"best.jpg";
  draw_error(ellipses1,imp,offset,img);
  imwrite(name,img);
 
  //********************************Code done here****************************/


  //relation
  //correcting Z value
  /*Double_t x[] = {15.17,16.06,14.72,13.39,12.29,10.98,9.82,8.32,7.14,5.94,-2.08,2.18,-6.73,-5.29,4.84,2.25,-11.39};
  Double_t y[] = {15.37227,15.89127,14.86923,13.49113,12.01067,10.56916,9.81259,8.305633,6.725209,5.767187,-2.668217,2.070681,-7.142699,-5.689827,4.720788,1.582416,-11.79712};
  TGraph *g = new TGraph((sizeof(x)/sizeof(Double_t)),x,y);
  TF1 *f1 = new TF1("f1","[0]+[1]*x",0,20);
  //TF1 *f = new TF1("f","[0]-[0]*[1]+[1]*x",0,20);
  g->Fit(f1);
  g->Draw("AL");
  g->Write();

  //for minimization output depth
  TF1 *f2 = new TF1("f2","[0]+[1]*x",-20,20);
  Double_t x1[] = {15.17,15.17,16.06,14.72,13.39,12.29,10.98,9.82,8.32,5.94,2.18,4.84,2.25,-2.08,-5.29,-6.73,-11.39};
  Double_t y1[] = {15.26,15.30,16.19,14.81,13.48,12.37,11.00,9.79,8.24,5.81,2.02,4.73,2.03,-2.41,-5.72,-7.17,-11.97};
  TGraph *g1 = new TGraph((sizeof(x1)/sizeof(Double_t)),x1,y1);
  g1->SetName("code output");
  //TF1 *f1 = new TF1("f1","[0]+[1]*x",0,20);
  //TF1 *f = new TF1("f","[0]-[0]*[1]+[1]*x",0,20);
  g1->Fit(f2);
  g1->Write();

  fout->Write();

  //finding the brightest part of the image.
  // img3.rows()
  // img3.cols()
  double r=2*400;
  cv::Point2f maxp;

  double max = -1;
  for(int i=r; i<img3.rows-r; i+=80){
    for(int j=r; j<img3.cols-r; j+=80){ 
      double sum = 0;
      for(int y=i-r;y<i+r;y++){
	for(int x=j-r;x<j+r;x++){
	  if(std::sqrt((x-j)*(x-j)+(y-i)*(y-i))<r){
	    sum+=img3.at<uchar>(y,x);
	  }
	}
      }
      if(sum>max){max=sum; maxp=cv::Point2f(j,i);}
    }
  }
  
  cv::circle( img2, maxp, r, cv::Scalar(0,255,250), 2 );
  cv::line( img2,cv::Point2f(2000,10),cv::Point2f(2000,2700),cv::Scalar( 0, 0, 255 ),2 );
  cv::line( img2,cv::Point2f(maxp.x,10),cv::Point2f(maxp.x,2700),cv::Scalar( 0, 255, 255 ),2 );
  imwrite("maxbright.jpg",img2);

  std::cout<<"Input r theta Z yaw pitch roll :";
  //vector<double>cam2;
  */
}
 
 return 0;
}

