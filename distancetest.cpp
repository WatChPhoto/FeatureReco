#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <fstream>                                                            
#include <cmath>
#include <iostream>
#include <vector>
#include <opencv2/calib3d.hpp>
#include <TH2D.h>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TFile.h>
#include "ImageDataReader.hpp"

//data structure
class TransformationData{
public:
  std::string id; //PMT id
  double el_x;  //centre of PMt
  double el_y;
  double el_b;  //minor axis of ellipse
  double el_e;  //eccentricity of ellipse
  double el_phi; //angle with x-axis.

  //mapped to 3d coordinate.
  double x;
  double y;
  double z;
  
  //TransformationData();
  TransformationData():id(""),el_x(0),el_y(0),el_b(0),el_e(0),el_phi(0),x(0),y(0),z(0){};
  TransformationData(std::string id, double el_x, double el_y, double el_b, double el_e, double el_phi, double x, double y, double z):id(id),el_x(el_x),el_y(el_y),el_b(el_b),el_e(el_e),el_phi(el_phi),x(x),y(y),z(z){}
  
  cv::Matx31d get_camera_coords(cv::Matx33d rmat, cv::Matx31d tvec){
    cv::Matx31d v (x,y,z);
    return rmat*v+tvec;
  }

  double get_a(){
    return el_b/std::sqrt(1.-el_e*el_e);
  }
};

std::istream& operator>>( std::istream& is, TransformationData& t){
  is>>t.id;
  is>>t.el_x;
  is>>t.el_y;
  is>>t.el_b;
  is>>t.el_e;
  is>>t.el_phi;
  is>>t.x;
  is>>t.y;
  is>>t.z;

  return is;
}

class WorldPmt{
public:
  std::string id;
  double x;
  double y;
  double z;
  
  WorldPmt():id(""),x(0),y(0),z(0){}
  WorldPmt(std::string id, double x, double y, double z):id(id),x(x),y(y),z(z){}
};

std::istream& operator>>( std::istream& is, WorldPmt& t){
  is>>t.id;
  is>>t.x;
  is>>t.y;
  is>>t.z;

  return is;
}

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

std::istream& operator>>( std::istream& is, Ellipse& p){
  is>>p.id;
  is>>p.x;
  is>>p.y;
  is>>p.b;
  is>>p.e;
  is>>p.phi;

  return is;
}

//read all pmts in sk detector.
std::vector<WorldPmt> read_all_pmts(std::string filename="SK_all_PMT_locations.txt"){
  std::vector<WorldPmt> all_pmts; 
  WorldPmt t; //temporary container

  std::ifstream myfile(filename);
  while (myfile >> t) {
    all_pmts.push_back(t);
    }
  myfile.close();
  return all_pmts;  
}

//read pmts from an image.
std::vector<Ellipse> read_ellipses_in_image(std::string filename ="045.txt"){
  std::vector<Ellipse>ellipses;
  Ellipse e;
  std::ifstream myfile(filename);
  while(myfile>>e){
    ellipses.push_back(e);
  }
  return ellipses;
}
double reprojection_error(std::vector<cv::Point2f> ellipses, std::vector<cv::Point2f> im_points,float offset){
  double dsq=0;
  for(int i=0; i<ellipses.size();i++){
    double x0 = ellipses[i].x;
    double y0 = ellipses[i].y;
    double lmin = 1000000000000000000000000000;
    for(int j=0;j<im_points.size();j++){
      double x1 = im_points[j].x;
      double y1 = im_points[j].y-offset;
      double d = (x1-x0)*(x1-x0)+(y1-y0)*(y1-y0);
      if(d<lmin){lmin=d;}
    }
    dsq +=lmin;
  }
  return dsq;
}

unsigned get_img_num(std::string filename){
  unsigned img_num;
  std::stringstream ss;  
  ss << filename;  
  ss >> img_num;

  return img_num;
}

class LabelledData{
public:
  std::string imnum;
  std::string id;
  double x; 
  double y;
  std::string by;
  LabelledData():imnum(""),id(""),x(0),y(0),by(""){}
};

std::istream& operator>>( std::istream& is, LabelledData& p){
  is>>p.imnum;
  is>>p.id;
  is>>p.x;
  is>>p.y;
  is>>p.by;
  
  return is;
}
void distance_range(double& low, double& high, std::vector<Ellipse> e){
  low =1000000000000000000000000;
  high=0;
  for(Ellipse ellipses: e){
    double a = ellipses.get_a();
    if(a<low){low=a;}
    if(a>high){high=a;}
  }
}

int main(int argc, char **argv){
  std::string filename = std::string(argv[argc-1]);
  std::string survey_id = std::string(argv[argc-2]);

  std::cout<<"filename = "<<filename<<std::endl;

  float offset = 250;
  float PI = std::acos(-1);
  double fx=2.760529621789217e+03;
  double fy=2.767014510543478e+03;
  double cx=1.914303537872458e+03;
  double cy=1.596386868474348e+03;
  //  cv::Matx31d rvec (1.49286548,-0.67181681,0.61404137);
  //   cv::Matx31d tvec(-161.14340624,1522.59325055,-991.13768429); 
  //  cv::Matx31d rvec (1.4661924,-0.67595185,0.61231617);
  //  cv::Matx31d tvec(-134.48375218,1549.37797742,-1028.98069242); 
  cv::Matx31d rvec (1.48220885,-0.85878027,0.75081903);
  cv::Matx31d tvec(-77.19825549,1313.21348846,-794.66238321); 

  cv::Matx14d dist_coeffs (-0.2398, 0.1145, 0,0); //just four par. radial & tangential.
  cv::Matx33d camera_matrix (fx, 0, cx,
			     0, fy, cy,
			     0,  0,  1);

  cv::Matx33d rmat; 
  cv::Rodrigues(rvec,rmat);


  std::cout<<"here "<<std::endl;
  cv::Mat image = imread ("046.jpg", cv::IMREAD_COLOR);//("045.jpg", cv::IMREAD_COLOR);
  std::cout<<"there "<<std::endl;
  cv::Mat img = image.clone();
  cv::Mat img1 = image.clone();
  cv::Mat img2 = image.clone();
 

  
  if(0){
    std::vector<TransformationData> v;
  std::vector<LabelledData> dat;
    std::ifstream MyReadFile("mapped.txt");
   std::cout<<"dat size = "<<dat.size()<<std::endl;  
  // Use a while loop together with the getline() function to read the file line by line
  
  TransformationData t;
  while (MyReadFile >> t) {
    // Output the text from the file
    //    std::cout << t.id<<std::endl;  
    v.push_back(t);
    //std::cout<<i<<std::endl;
  }
  
  // Close the file
  MyReadFile.close(); 
  
  std::ifstream my1;
  my1.open(filename+"lbl.txt");
  LabelledData d;

  while(my1 >> d){dat.push_back(d); std::cout<<"pushing"<<std::endl;}
  my1.close();

    std::vector<cv::Point3f> object_points;
  
  for(int i=0; i<v.size();i++){
    object_points.push_back( cv::Point3f(v[i].x,v[i].y,v[i].z));
    
    //   std::cout<<v[i].id<<'\t'<<v[i].el_x<<'\t'<<v[i].el_y<<'\t'<<v[i].el_b<<'\t'<<v[i].el_e<<'\t'<<v[i].el_phi<<'\t'<<v[i].x<<'\t'<<v[i].y<<'\t'<<v[i].z<<std::endl;
  }
  std::cout<<"1ere "<<std::endl;    
  std::vector<cv::Point2f>image_points;
  projectPoints(object_points,rvec,tvec,camera_matrix,dist_coeffs,image_points);
    for(int i=0;i<image_points.size();i++){
    cv::circle( image, cv::Point( image_points[i].x, image_points[i].y ), 10, cv::Scalar(0,255,250), -1 );
  }
  imwrite("output.jpg",image);
  
  for(int i=0;i<dat.size();i++){
 
    std::cout<<"(x,y) = "<<"( "<<dat[i].x<<" , "<<dat[i].y<<" )"<<std::endl;
    std::cout<<dat[i].id.substr(6)<<std::endl;
    if(dat[i].id.substr(6)=="00"){
      cv::circle( img2, cv::Point(dat[i].x,dat[i].y), 10, cv::Scalar(0,255,250), -1 );
      std::string text = dat[i].id.substr(0,5);
      text=std::to_string(get_img_num(text));
      cv::putText(img2, text,cv::Point( dat[i].x, dat[i].y ) , cv::FONT_HERSHEY_PLAIN,4, cv::Scalar(0,255,250),3);
    }
  }
  imwrite("manlabelled.jpg",img2);
  //TH1D* zcoord   = new TH1D( "zcoord"  ,"Z coord of features  ; size; distance",256,0.,256.);
  
  //   TCanvas *c1 = new TCanvas("c1","Graph of size vs distance",200,10,500,300);
  
  TFile* fout = new TFile("testingdistance.root","recreate");
  Double_t x[v.size()], y[v.size()];
  Int_t n = v.size();
  for (Int_t i=0;i<n;i++) {
    x[i] = v[i].get_a();
    y[i] = (v[i].get_camera_coords(rmat, tvec))(2,0);
  }
  TGraph* gr = new TGraph(n,x,y);
  gr->SetNameTitle("dist vs a","dist vs a");
  gr->Write();
  
  for (Int_t i=0;i<n;i++) {
    x[i] = v[i].el_b;
    y[i] = (v[i].get_camera_coords(rmat, tvec))(2,0);
  }
  TGraph* gr1 = new TGraph(n,x,y);
  gr1->SetNameTitle("dist vs b","dist vs b");
  gr1->Write();
  
  //camera position is 
  cv::Matx31d p = -rmat.t()*tvec;
  for (Int_t i=0;i<n;i++) {
    x[i] = v[i].el_b;
    y[i] = std::sqrt((p(0,0)-v[i].x)*(p(0,0)-v[i].x)+(p(1,0)-v[i].y)*(p(1,0)-v[i].y)+(p(2,0)-v[i].z)*(p(2,0)-v[i].z));
  }
  TGraph* gr2 = new TGraph(n,x,y);
  gr2->SetNameTitle("b vs dist","b vs dist");
  gr2->Write();
  
  for (Int_t i=0;i<n;i++) {
    x[i] = v[i].get_a();
    y[i] = std::sqrt((p(0,0)-v[i].x)*(p(0,0)-v[i].x)+(p(1,0)-v[i].y)*(p(1,0)-v[i].y)+(p(2,0)-v[i].z)*(p(2,0)-v[i].z));
  }
  TGraph* gr3 = new TGraph(n,x,y);
  gr3->SetNameTitle("a vs dist","a vs dist");
  gr3->Write();
  fout->Write();
 
  }
  
  
  
  /*====================================================================
    this is new part.
    ===================================================================*/
  
  ImageDataReader & idr = ImageDataReader::GetInstance();
  


  ImageMetaData m = idr.GetMetaData(survey_id, get_img_num(filename));
  ImageMetaData i239 = idr.GetMetaData("BarrelSurveyFar", 239);
  std::cout<<m.face.yaw<<std::endl;
  std::cout<<m.face.pitch<<std::endl;
  std::cout<<m.face.roll<<std::endl;

  std::cout<<i239.face.yaw<<std::endl;
  std::cout<<i239.face.pitch<<std::endl;
  std::cout<<i239.face.roll<<std::endl;
  std::cout<<"he_pmts"+filename+".txt"<<std::endl;
  std::vector<Ellipse> ellipses1 = read_ellipses_in_image("he_pmts"+filename+".txt");
  std::vector<cv::Point2f> ellipses;
  for(int i=0; i<ellipses1.size();i++){
    ellipses.push_back(cv::Point2f(ellipses1[i].x,ellipses1[i].y));
  }
  
  for(int i=0;i<ellipses.size();i++){
    cv::circle( img1, cv::Point( ellipses[i].x, ellipses[i].y ), 10, cv::Scalar(0,255,250), -1 );
  }
  imwrite("textreadimg.jpg",img1);

  //0 for 045
  double yaw0 = i239.face.yaw*PI/180.;//64.0*PI/180.;
  double pitch0 = i239.face.pitch*PI/180.;//2.0*PI/180.;
  double roll0 = i239.face.roll*PI/180.;//-1.0*PI/180.;  
  //this is for any image we want to work with.
  //this need to be read for particular image.
  double yaw1 = m.face.yaw*PI/180.;//59.0*PI/180.;
  double pitch1 = m.face.pitch*PI/180.;//2.0*PI/180.;
  double roll1 = m.face.roll*PI/180.;  //-1.0*PI/180.;  
  double z1 = m.depth*100;
  //rotation about x,y,z axis of camera.
  double thx0=-pitch0;
  double thy0=-yaw0;
  double thz0=roll0;
  
  double thx1=-pitch1;
  double thy1=-yaw1;
  double thz1=roll1;
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
  
  cv::Matx33d ACx (1.0,   0.0,      0.0,
		  0.0, cos(thx1), -sin(thx1),
		  0.0, sin(thx1),  cos(thx1));
  cv::Matx33d ACy (cos(thy1), 0.0, sin(thy1),
		  0.0,      1.0,    0.0,
		  -sin(thy1), 0.0, cos(thy1));
  cv::Matx33d ACz (cos(thz1), -sin(thz1), 0.0,
		  sin(thz1),  cos(thz1), 0.0,
		  0.0,        0.0,     1.0 );
  
  //frame A->B->C
  //         ^
  //         |
  //       World
  //given the coordinte in frame A, finds the matrix that will give coordinate in frame B.
  cv::Matx33d rot_AB = ABz.t()*ABx.t()*ABy.t();
  cv::Matx33d rot_AC = ACz.t()*ACx.t()*ACy.t();
  cv::Matx33d rot_BC = rot_AC*(rot_AB.t()); //given the coordinate in frame B, finds the coordinate in frame C.
  cv::Matx33d R = rot_BC*rmat;//rmat brings coordinate from world to 045 frame(B frame).x

  cv::Matx31d rvec1;
  cv::Rodrigues(R,rvec1);
  cv::Matx31d tvec1;

  std::cout<<"046 rotation matrix"<<std::endl;
  std::cout<<R(0,0)<<"\t"<<R(0,1)<<"\t"<<R(0,2)<<std::endl;
  std::cout<<R(1,0)<<"\t"<<R(1,1)<<"\t"<<R(1,2)<<std::endl;
  std::cout<<R(2,0)<<"\t"<<R(2,1)<<"\t"<<R(2,2)<<std::endl;
 
  std::cout<<"Rotation vector"<<std::endl;
  std::cout<<rvec1(0,0)<<"\t"<<rvec1(0,1)<<"\t"<<rvec1(0,2)<<std::endl;

  cv::Matx33d sss;
  cv::Rodrigues(cv::Matx31d(1.48220885,-0.85878027, 0.75081903),sss);
  std::cout<<"239 real"<<std::endl;
  std::cout<<sss(0,0)<<"\t"<<sss(0,1)<<"\t"<<sss(0,2)<<std::endl;
  std::cout<<sss(1,0)<<"\t"<<sss(1,1)<<"\t"<<sss(1,2)<<std::endl;
  std::cout<<sss(2,0)<<"\t"<<sss(2,1)<<"\t"<<sss(2,2)<<std::endl;

  /*
  for(int i=0; i<v.size(); i++){
    cv::Matx31d camera_coords = v[i].get_camera_coords(rmat, tvec);
  }
  */
  
  std::vector<WorldPmt> all_pmts = read_all_pmts();
  
  //std::cout<<"Now printing th read world coordinate of pmts"<<std::endl;
  /*  for(int i=0; i<all_pmts.size();i++){
    // std::cout<<all_pmts[i].id<<"\t"<<all_pmts[i].x<<"\t"<<all_pmts[i].y<<"\t"<<all_pmts[i].z<<std::endl;
  }
  */
  
  //we want to move the camera in x,y plane while keeping the z position same.
  //double z=1485.96;//15.17*100;
  double min_rep_err = 10000000000000000000000000;
  double minx=0;
  double miny=0;
  double minz=0;
  std::vector<int> final_id;
  std::vector<cv::Point3f> obv;
  std::vector<cv::Point2f> imp;
  
  cv::Matx31d dir =R.t()*cv::Matx31d(0,0,1); //direction of z-axis of camera.
  double theta = std::atan2(dir(1,0),dir(0,0))*180./PI; // angle with x-axis.
  theta=(theta>=0)?theta:(360+theta);
  std::cout<<"theta ="<<theta<<std::endl;
  double range = 20;
  double th1 = ((theta-range)>=0)?(theta-range):(360+theta-range);
  double th2 = ((theta+range)>=0)?(theta+range):(360+theta+range);
  th2 = (th2>360)?(th2-360.):th2;
  double high,low;
  
  if((th2-th1)>0){high=th2;low=th1;}

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
  double low_r = 1600.-(86190./lower_r);
  double high_r = 1800.-(86190./upper_r);
  std::cout<<"low_r = "<<low_r<<"\t"<<"high_r = "<<high_r<<std::endl;

  
  for(int x=0; x<1500; x+=10){
    for(int y=0; y<1500; y+=10){
      for(int z=int(z1)-100;z<int(z1)+100; z+=10){
	float r=std::sqrt(x*x+y*y);//+z*z);
	double alfa = std::atan2(y,x)*180./PI; //angle x,y make.
	alfa = (alfa>=0)?alfa:(360+alfa);
	if(r<low_r||r>high_r){continue;} //rad of tank=1690cm
	else if((high-low)>(2*range)){ if(alfa>=low && alfa<=high){continue;}}
	else{if(alfa<low||alfa>high){continue;}}
	
	//std::cout<<"alfa = "<<alfa<<std::endl;
	//if(alfa<low ||alfa>high){continue;}
	//       std::cout<<"passed few the hurdles"<<std::endl;
	cv::Matx31d p(x, y, z);
	cv::Matx31d tv = -R*p;
	std::vector<cv::Point3f> object_in_view;
	std::vector<int> id;
	for(int i=0; i<all_pmts.size();i++){
	  cv::Matx31d xp = R*cv::Matx31d(all_pmts[i].x,all_pmts[i].y,all_pmts[i].z)+tv;
	  if(xp(2,0)>0){
	    object_in_view.push_back(cv::Point3f(all_pmts[i].x,all_pmts[i].y,all_pmts[i].z));
	    int num;
	    std::string str = all_pmts[i].id.substr(0,5);
	    std::stringstream ss;  
	    ss << str;  
	    ss >> num;
	    id.push_back(num);
	    //	   std::cout<<"passed all the hurdles"<<std::endl;
	  }
	  
	}
	if(object_in_view.size()>0){
	  std::vector<cv::Point2f>im_points;
	  projectPoints(object_in_view,rvec1,tv,camera_matrix,dist_coeffs,im_points);
	  
	  double rep_error = reprojection_error(ellipses, im_points,offset);
	  //	  std::cout<<"Repo error = "<<rep_error<<std::endl;
	  if(rep_error<=min_rep_err ){
	    //obv = object_in_view;
	    final_id = id;
	    imp=im_points;
	    min_rep_err = rep_error;
	    minx=x;
	    miny=y;
	    minz=z;
	    tvec1=tv;
	  } 
	}
      }
    }
  }
  
  std::cout<<"position is ("<<minx<<" , "<<miny<<" , "<<minz<<" )"<<std::endl;
  if(1){
    //20001: UK B1
    //30001: Korean B1
    std::vector<std::string>label{"20001","UKB2","UKB3","UKB4","UKB5","30001","KB2","KB3","KB4","KB5"};
    std::vector<cv::Point2f>im_points;
    std::vector<cv::Point3f>object_in_view{cv::Point3f(1490.73,768.14,1302.95),cv::Point3f(1490.73,768.14,666.65),cv::Point3f(1490.73,768.14,30.35),cv::Point3f(1490.73,768.14,-676.65),cv::Point3f(1490.73,768.14,-1312.95),cv::Point3f(1490.73,768.14,1232.25),cv::Point3f(1490.73,768.14,595.95),cv::Point3f(1490.73,768.14,-40.35),cv::Point3f(1490.73,768.14,-605.95),cv::Point3f(1490.73,768.14,-1242.25)};
    projectPoints(object_in_view,rvec1,tvec1,camera_matrix,dist_coeffs,im_points);
 
    for(int i=0;i<im_points.size();i++){
      int x=im_points[i].x;
      int y=im_points[i].y;
      if(x>0 && x<4000 && y>0 && y-offset<2750){
	cv::circle( img, cv::Point( x, y-offset ), 10, cv::Scalar(0,255,250), -1 );
	cv::putText(img, label[i],cv::Point( x, y-offset ) , cv::FONT_HERSHEY_PLAIN,4, cv::Scalar(0,255,250),3);
      }
    }
  }
  /*  std::ofstream my;
    my.open ("answer.txt");
    my << "position is ("<<minx<<" , "<<miny<<" )"<<std::endl;
    my.close();*/
  //std::vector<cv::Point2f> imp;
  //projectPoints(object_in_view,rvec,tv,camera_matrix,dist_coeffs,im_points);
    
    for(int i=0;i<imp.size();i++){
    cv::circle( img, cv::Point( imp[i].x, imp[i].y-offset ), 10, cv::Scalar(0,255,250), -1 );
    if(imp[i].x>0 && imp[i].x<4000 && imp[i].y>0 && imp[i].y-offset<2700){
      std::string text = std::to_string(final_id[i]);
      cv::putText(img, text,cv::Point( imp[i].x, imp[i].y-offset ) , cv::FONT_HERSHEY_PLAIN,4, cv::Scalar(0,255,250),3);
    }
  }
  imwrite("best.jpg",img);
  return 0;
}
