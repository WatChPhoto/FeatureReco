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
 
  /*std::vector<double> d;
    for(int i=0; i<ellipses.size();i++){
    double x0 = ellipses[i].x;
    double y0 = ellipses[i].y;
    double lmin = 1000000000000000000000000000;
    //int indx=-1;
    //double lmin0 = 1000000000000000000000000000;
    for(int j=0;j<im_points.size();j++){
      double x1 = im_points[j].x;
      double y1 = im_points[j].y-offset;
      double d1 = (x1-x0)*(x1-x0)+(y1-y0)*(y1-y0);
      if(d1<lmin){
	bool reverse = true;
	for(int k=0; k<ellipses.size();k++){
	  double x2 = ellipses[k].x;
	  double y2 = ellipses[k].y;
	  double d2 = (x2-x1)*(x2-x1)+(y2-y1)*(y2-y1);
	  if(d2<d1){reverse =false; break;}
	  //if(d2<lmin){indx=k;}
	  
	}
	
	if(reverse){lmin=d1;}
      }
    }
    //dsq +=lmin;
    if(lmin!=1000000000000000000000000000){
      d.push_back(lmin);
    }
  }
  
  //only add the distance that is very small.
  int count=0;
  double sum=0;
  for(int i=0; i<d.size();i++){
    //  std::cout<<"sum entry is  = "<<d[i]<<std::endl;
    if(d[i]<(70*70)){sum += d[i]; count++;}
  }
  
  
  
  if(count>0){
    std::cout<<"sum is = "<<sum<<std::endl;
    std::cout<<"count is = "<<count<<std::endl;
    std::cout<<"error is = "<<sum/count<<std::endl;
    
    return sum; 
    
  }
  else{
    return 10000000000000000000000000000000;
  }
  */
 
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

void make_error_histogram(std::vector<cv::Point3f>ev, int z){
  std::string lbl = std::to_string(z);
  TH2D *error = new TH2D(("error for z = "+lbl).c_str(),"Reprojection error;x;y",1500.,0.,1500.,1500.,0.,1500.);
  for(cv::Point3f p:ev){
    error->Fill(p.x,p.y,p.z);
  }
}

struct cam{
public:
  double r;
  double theta;
  double z;
  double yaw;
  double pitch;
  double roll;
  cam(){
    r=0;
    theta=0;
    z=0;
    yaw=0;
    pitch=0;
    roll=0;
  }
  cam(double r, double theta, double z, double yaw, double pitch, double roll):r(r),theta(theta),z(z),yaw(yaw),pitch(pitch),roll(roll){}
};


bool clk=false; //clicked?
cv::Mat get_scene(std::vector<WorldPmt> all_pmts,std::vector<cv::Point3f>light_injectors,cam cam1,cv::Matx33d R, cv::Matx31d p){
  double PI=std::acos(-1);
  cv::Mat scene(3000, 4000, CV_8UC3, cv::Scalar(255,255,255));
  //camera matrix;
  double fx=2.760529621789217e+03;//1100;//
  double fy=2.767014510543478e+03;//1100;//
  double cx=scene.cols/2.;//1.914303537872458e+03;
  double cy=scene.rows/2.;//1.596386868474348e+03;
  cv::Matx33d camera_matrix(fx,0,cx,
			    0,fy,cy,
			    0,0,1);
  cv::Matx14d dist_coeffs (0, 0, 0, 0); //just four par. radial & tangential.x

  double thx = cam1.pitch*PI/180;
  double thy = cam1.yaw*PI/180;
  double thz = cam1.roll*PI/180;
  
  cv::Matx33d Rx (1.0,   0.0,      0.0,
		   0.0, cos(thx), -sin(thx),
		   0.0, sin(thx),  cos(thx));
  cv::Matx33d Ry (cos(thy), 0.0, sin(thy),
		   0.0,      1.0,    0.0,
		   -sin(thy), 0.0, cos(thy));
  cv::Matx33d Rz (cos(thz), -sin(thz), 0.0,
		   sin(thz),  cos(thz), 0.0,
		   0.0,        0.0,     1.0 );
  //-ve 90 degree rotation in x
  cv::Matx33d Rf0 (1.0,   0.0,      0.0,
  		   0.0, cos(-PI/2.), -sin(-PI/2.),
		   0.0, sin(-PI/2.),  cos(-PI/2.));

  //+ve 90 degree rotation in then y
  cv::Matx33d Rf1 (cos(PI/2.), 0.0, sin(PI/2.),
		   0.0,      1.0,    0.0,
		   -sin(PI/2.), 0.0, cos(PI/2.));

  //given the coordinte in unrotated frame, finds the matrix that will give coordinate in rotated frame.
  //rotation matrix for cam1.
  cv::Matx33d R1 = Rz.t()*Rx.t()*Ry.t()*Rf1.t()*Rf0.t();
  
  //vector along three axes of camera;
  //let's do 30 cm in each dir
  cv::Matx31d x(30,0,0);
  cv::Matx31d y(0,30,0);
  cv::Matx31d z(0,0,30);
  //endpoints of vector
  x = R.t()*x+p;
  y = R.t()*y+p;
  z = R.t()*z+p;
  cv::Matx31d rvec;
  cv::Rodrigues(R1,rvec);
  cv::Matx31d pos(cam1.r*cos(cam1.theta*PI/180),cam1.r*sin(cam1.theta*PI/180.),cam1.z);
  cv::Matx31d tvec = -R1*pos;
  std::vector<cv::Point3f> object_in_view;
  std::vector<cv::Point2f>im_points;
  
  for(int i=0; i<all_pmts.size();i++){
    cv::Matx31d xp = R1*cv::Matx31d(all_pmts[i].x,all_pmts[i].y,all_pmts[i].z)+tvec;
    if(xp(2,0)>0){
      object_in_view.push_back(cv::Point3f(all_pmts[i].x,all_pmts[i].y,all_pmts[i].z));
    } 
  }
  
  projectPoints(object_in_view,rvec,tvec,camera_matrix,dist_coeffs,im_points);  
  
  for(int i=0;i<im_points.size();i++){
    cv::circle( scene, cv::Point( im_points[i].x, im_points[i].y ), 10, cv::Scalar(0,240,240), -1 );
  }
  
  std::vector<cv::Point3f> obj={ cv::Point3f(p(0,0),p(1,0),p(2,0)), cv::Point3f(x(0,0),x(1,0),x(2,0)), cv::Point3f(y(0,0),y(1,0),y(2,0)), cv::Point3f(z(0,0),z(1,0),z(2,0)) };
  std::vector<cv::Point2f>image_points;
  projectPoints(obj,rvec,tvec,camera_matrix,dist_coeffs,image_points);

  
  //Drawing camera direction
  cv::arrowedLine( scene, cv::Point( image_points[0].x, image_points[0].y ),cv::Point( image_points[1].x, image_points[1].y ), cv::Scalar(0,0,250), 2 );
  cv::arrowedLine( scene, cv::Point( image_points[0].x, image_points[0].y ),cv::Point( image_points[2].x, image_points[2].y ), cv::Scalar(0,250,0), 2 );
  cv::arrowedLine( scene, cv::Point( image_points[0].x, image_points[0].y ),cv::Point( image_points[3].x, image_points[3].y ), cv::Scalar(250,0,0), 2 );

  //projecting light injectors
  std::vector<cv::Point2f> injector_points;
  projectPoints(light_injectors,rvec,tvec,camera_matrix,dist_coeffs,injector_points);
  for(int i=0;i<injector_points.size();i++){
      int x=injector_points[i].x;
      int y=injector_points[i].y;
      if(x>0 && x<scene.cols && y>0 && y<scene.rows){
	cv::circle( scene, cv::Point( x, y), 5, cv::Scalar(220,31,237), -1 );
      }
    }

  //drawing sk-coordinate
  std::vector<cv::Point3f> axes;
  std::vector<cv::Point2f>axes_points;
  std::vector<cv::Scalar> color;
  //selecting points that will be in frame.
  double s=-tvec(2,0)/(R1(2,0)*2000);
  
  if(s>=0 && s<1){
    // then add s*<2000,0,0> to the points
    axes.push_back(cv::Point3f(s*2000,0,0));
    axes.push_back(cv::Point3f(2000,0,0));
    color.push_back(cv::Scalar(0,0,255));
  }
  s=-tvec(2,0)/(R1(2,1)*2000);
  if(s>=0 && s<1){
    // then add s*<0,2000,0> to the points
    axes.push_back(cv::Point3f(0,s*2000,0));
    axes.push_back(cv::Point3f(0,2000,0));
    color.push_back(cv::Scalar(0,255,0));
  }
  s=-tvec(2,0)/(R1(2,2)*2000);
  if(s>=0 && s<1){
    // then add s*<0,2000,0> to the points
    axes.push_back(cv::Point3f(0,0,s*2000));
    axes.push_back(cv::Point3f(0,0,2000));
    color.push_back(cv::Scalar(255,0,0));
  }
  
  if(axes.size()>0){
    projectPoints(axes,rvec,tvec,camera_matrix,dist_coeffs,axes_points);
  int j=0;
  for(int i=0; i<axes_points.size(); i+=2){
    cv::arrowedLine( scene, cv::Point( axes_points[i].x, axes_points[i].y ),cv::Point( axes_points[i+1].x, axes_points[i+1].y ), color[j], 2 );
    j++;
  }
}
  
  /*
	cv::Matx31d xp = R1*cv::Matx31d(all_pmts[i].x,all_pmts[i].y,all_pmts[i].z)+tvec;
	if(xp(2,0)>0){
	object_in_view.push_back(cv::Point3f(all_pmts[i].x,all_pmxts[i].y,all_pmts[i].z));
	} 
      */
    //Sk axes direction
  /*
  cv::arrowedLine( scene, cv::Point( axes_points[0].x, axes_points[0].y ),cv::Point( axes_points[1].x, axes_points[1].y ), cv::Scalar(0,0,250), 2 );
  cv::arrowedLine( scene, cv::Point( axes_points[0].x, axes_points[0].y ),cv::Point( axes_points[2].x, axes_points[2].y ), cv::Scalar(0,250,0), 2 );
  cv::arrowedLine( scene, cv::Point( axes_points[0].x, axes_points[0].y ),cv::Point( axes_points[3].x, axes_points[3].y ), cv::Scalar(250,0,0), 2 );
  */
  clk=false;
  return scene;
}




//Icon class
class Icon{
public:
  cv::Mat img;
  cv::Point2f pos;
  std::string text;
  cv::Scalar text_color;
  cv::Scalar color;
  cv::Point2f size;
  int index;
  bool sign; //true for + and false for -
  Icon(){
    pos=cv::Point2f(0,0);
    text="";
    text_color=cv::Scalar(0,0,0);
    color=cv::Scalar(0,0,0);
    size=cv::Point2f(100,50);
    index=-1;
    sign=true;
  }


  

  Icon(cv::Mat &img,cv::Point2f pos, cv::Point2f size, std::string text, cv::Scalar color,int index, bool sign, cv::Scalar text_color=cv::Scalar(255,255,255)):img(img),pos(pos),size(size),text(text),color(color),index(index),sign(sign),text_color(text_color){}
  
  void show(int x, int y,bool clicked,std::vector<float> &c1){
    cv::Scalar fill =color;
    cv::Scalar t_fill=text_color;
    if(is_inside(x,y)){
      fill=cv::Scalar(0,250,0); 
      t_fill=cv::Scalar(0,0,0);
      if(clicked){
	clk=true;
	fill=cv::Scalar(255,255,255); t_fill=cv::Scalar(0,250,0);
	if(index!=-1){
	  if(sign){
	    c1[index]++;
	  }
	  else{
	    c1[index]--;
	  }
	}
	else{
	  // r, theta, z, yaw, pitch, roll.
	  std::cout<<"enter r theta z yaw pitch roll : ";
	  std::cin>>c1[0]>>c1[1]>>c1[2]>>c1[3]>>c1[4]>>c1[4];
	  //getting rid of token in buffer
	  while (std::cin.get() != '\n') 
	    {
	      continue;
	    }

	}
      }
    }
    std::cout<<c1[0]<<" "<<c1[1]<<" "<<c1[2]<<" "<<c1[3]<<" "<<c1[4]<<" "<<c1[5]<<std::endl;
    rectangle(img,pos,pos+size,fill,-1);
    putText(img, text,cv::Point2f(pos.x,pos.y+size.y)+cv::Point2f(10,-10),cv::FONT_HERSHEY_SIMPLEX,2,t_fill,2);
  } 

  bool is_inside(int x, int y){
    if(x>=pos.x && x<=pos.x+size.x && y>=pos.y && y<=pos.y+size.y){
      return true;
    }
    else{
      return false;
    }
  }
};


cv::Mat outs;
cam c1;
cv::Mat image_to_show;
bool left_down=false;
bool left_up=false;
std::vector<float> cam1(6,0);
std::vector<WorldPmt> al_pmt;
std::vector<cv::Point3f> li; 
cv::Matx33d R_o;
cv::Matx31d ps;
//r,theta,z,yaw,pitch,roll;

void onMouse(int event, int x, int y, int flags, void* userdata){
  bool clicked=false;
  if(event==cv::EVENT_LBUTTONDOWN){
    clicked=true;
  }
  image_to_show=outs.clone();
  
  int xmin = image_to_show.cols-900;
  int ymin = image_to_show.rows-800;
  int xmax = image_to_show.cols;
  int ymax = image_to_show.rows;
  if(x>xmin && y>ymin && x<xmax && y<ymax){
    cv::Point2f size = cv::Point2f(200,100);
    int xoffset=100;
    int yoffset=50;
    Icon ri(image_to_show, cv::Point2f(xmin,ymin),size,"Rad+",cv::Scalar(0,0,0),0,true,cv::Scalar(255,255,255));
    Icon rd(image_to_show, cv::Point2f(xmin,ymin+4*(size.y+yoffset)),size,"Rad-",cv::Scalar(0,0,0),0,false,cv::Scalar(255,255,255));
    Icon thi(image_to_show, cv::Point2f(xmin+size.x+xoffset,ymin),size,"th+",cv::Scalar(0,0,0),1,true,cv::Scalar(255,255,255));
    Icon thd(image_to_show, cv::Point2f(xmin+size.x+xoffset,ymin+4*(size.y+yoffset)),size,"th-",cv::Scalar(0,0,0),1,false,cv::Scalar(255,255,255));
    Icon zi(image_to_show, cv::Point2f(xmin+2*(size.x+xoffset),ymin),size,"Z+",cv::Scalar(0,0,0),2,true,cv::Scalar(255,255,255));
    Icon zd(image_to_show, cv::Point2f(xmin+2*(size.x+xoffset),ymin+4*(size.y+yoffset)),size,"Z-",cv::Scalar(0,0,0),2,false,cv::Scalar(255,255,255));

    Icon yawd(image_to_show, cv::Point2f(xmin,ymin+size.y+yoffset),size,"Yaw-",cv::Scalar(0,0,0),3,false,cv::Scalar(255,255,255));
    Icon yawi(image_to_show, cv::Point2f(xmin+2*(size.x+xoffset),ymin+size.y+yoffset),size,"Yaw+",cv::Scalar(0,0,0),3,true,cv::Scalar(255,255,255));
    Icon pitd(image_to_show, cv::Point2f(xmin,ymin+2*(size.y+yoffset)),size,"Pih-",cv::Scalar(0,0,0),4,false,cv::Scalar(255,255,255));
    Icon piti(image_to_show, cv::Point2f(xmin+2*(size.x+xoffset),ymin+2*(size.y+yoffset)),size,"Pih+",cv::Scalar(0,0,0),4,true,cv::Scalar(255,255,255));
    Icon rold(image_to_show, cv::Point2f(xmin,ymin+3*(size.y+yoffset)),size,"Rol-",cv::Scalar(0,0,0),5,false,cv::Scalar(255,255,255));
    Icon roli(image_to_show, cv::Point2f(xmin+2*(size.x+xoffset),ymin+3*(size.y+yoffset)),size,"Rol+",cv::Scalar(0,0,0),5,true,cv::Scalar(255,255,255));

    Icon set_to(image_to_show, cv::Point2f(xmin+size.x+xoffset/2,ymin+size.y+yoffset/2),cv::Point2f(300,450),"<<\nSet",cv::Scalar(0,0,0),-1,true,cv::Scalar(255,255,255));
    
        
    ri.show(x,y,clicked,cam1);
    rd.show(x,y,clicked,cam1);
    thi.show(x,y,clicked,cam1);
    thd.show(x,y,clicked,cam1);
    zi.show(x,y,clicked,cam1);
    zd.show(x,y,clicked,cam1);
    yawi.show(x,y,clicked,cam1);
    yawd.show(x,y,clicked,cam1);
    piti.show(x,y,clicked,cam1);
    pitd.show(x,y,clicked,cam1);
    roli.show(x,y,clicked,cam1);
    rold.show(x,y,clicked,cam1);
    set_to.show(x,y,clicked,cam1);
  }




  if(clk){
    c1=cam(cam1[0],cam1[1],cam1[2],cam1[3],cam1[4],cam1[5]);
    outs = get_scene(al_pmt,li, c1, R_o, ps);
  } 
  
  imshow("Scene",image_to_show);

  //imshow("Scene",outs);
}


int main(int argc, char **argv){
  
  /*
  cv::Matx31d r0 (1.45565547,-0.83269149,0.73368359);
  cv::Matx31d t0(-185.21618339,-211.75779124,-880.66337725); 
  cv::Matx33d R0;
  cv::Rodrigues(r0,R0);
  cv::Matx31d pos = -R0.t()*t0;
  std::cout<<"position is ="<<pos(0,0)<<"\t"<<pos(1,0)<<"\t"<<pos(2,0)<<std::endl;
  */

  std::string filename = std::string(argv[argc-1]);
  std::string survey_id = std::string(argv[argc-2]);

  std::cout<<"filename = "<<filename<<std::endl;

  float offset = 250;
  double PI = std::acos(-1);
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

  cv::Matx14d dist_coeffs (-0.2398, 0.1145, 0, 0); //just four par. radial & tangential.x
  cv::Matx33d camera_matrix (fx, 0, cx,
			     0, fy, cy,
			     0,  0,  1);

  cv::Matx33d rmat; 
  cv::Rodrigues(rvec,rmat);


  std::cout<<"here "<<std::endl;
  cv::Mat image = imread (filename+".jpg", cv::IMREAD_COLOR);//("045.jpg", cv::IMREAD_COLOR);
  
  std::cout<<"there "<<std::endl;
  cv::Mat img = image.clone();
  cv::Mat img1 = image.clone();
  cv::Mat img2 = image.clone();
  
  cv::Mat img3; 
  cvtColor(image, img3, cv::COLOR_BGR2GRAY);

  
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
  double yaw0 = (i239.face.yaw)*PI/180.;//64.0*PI/180.;
  double pitch0 = (i239.face.pitch)*PI/180.;//2.0*PI/180.;
  double roll0 = (i239.face.roll)*PI/180.;//-1.0*PI/180.;  
  //this is for any image we want to work with.
  //this need to be read for particular image.
  double yaw1 = (m.face.yaw)*PI/180.;//59.0*PI/180.;
  double pitch1 = (m.face.pitch)*PI/180.;//2.0*PI/180.;
  double roll1 = m.face.roll*PI/180.;  //-1.0*PI/180.;  
  double b = 1.02585; //ratio of density of sea water to pure water
  double c = -0.271196;//offset
  double z1 = m.depth;
  z1 = c+b*z1;
  z1 *=100; //converting to cm

  std::cout<<"z depth is = "<<z1<<std::endl;
  //rotation about x,y,z axis of camera.
  double thx0=pitch0;
  double thy0=-yaw0;
  double thz0=-roll0;
  
  double thx1=pitch1;
  double thy1=-yaw1;
  double thz1=-roll1;
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
  //cv::Matx33d rot_AB =  ABy.t()*ABx.t()*ABz.t();
  //cv::Matx33d rot_AC = ACy.t()*ACx.t()*ACz.t();

  cv::Matx33d rot_BC = rot_AC*(rot_AB.t()); //given the coordinate in frame B, finds the coordinate in frame C.
  cv::Matx33d R = rot_BC*rmat;//rmat brings coordinate from world to 045 frame(B frame).x
  R_o = R;
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
  al_pmt = all_pmts;
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
  double range =0;//25;
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

  //
  std::vector<cv::Point2f>cord;
  cord.push_back(cv::Point2f(low_r*cos(low),low_r*sin(low)));
  cord.push_back(cv::Point2f( low_r*cos(high),low_r*sin(high)));
  cord.push_back(cv::Point2f(high_r*cos(low),high_r*sin(low)));
  cord.push_back(cv::Point2f(high_r*cos(high),high_r*sin(high)));

  cv::Point2f l=cv::Point2f(2000,2000);
  cv::Point2f h=cv::Point2f(0,0);

  for(int i=0; i<cord.size();i++){
    if(cord[i].x>h.x){
      h.x=cord[i].x;
    }
    if(cord[i].x<l.x){
      l.x=cord[i].x;
    }
    if(cord[i].y>h.y){
      h.y=cord[i].y;
    }
    if(cord[i].y<l.y){
      l.y=cord[i].y;
    }
    std::cout<<"( "<<cord[i].x<<","<<cord[i].y<<" )"<<std::endl;
  }
  std::cout<<"high and low"<<std::endl;
  std::cout<<"( "<<h.x<<","<<h.y<<" )"<<std::endl;
  std::cout<<"( "<<l.x<<","<<l.y<<" )"<<std::endl;
  
  TFile* fout = new TFile("testingdistance.root","recreate");
  //TH3D *error =   new TH3D("reprojection error","Reprojection error;x;y;z",int((h.x-l.x+20)/5),l.x-10,h.x+10,int((h.y-l.y+20)/5),l.y-10,h.y+10,100.,z1-50,z1+50.);

  std::cout<<"z1 = "<<z1<<std::endl;
  
  double ang;
  double rad;
  for(int z=z1-4;z<z1+4; z+=1){	
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
  
  std::cout<<"position is ("<<minx<<" , "<<miny<<" , "<<minz<<" )"<<std::endl;
  ps=cv::Matx31d(minx,miny,minz);
  
  std::cout<<"rad = "<<rad<<"  ang = "<<ang<<std::endl;
  cv::Matx31d a = R.t()*ACz.t()*ACx.t()*cv::Matx31d(0,0,1);
  std::cout<<"angle of direction is ="<<std::atan2(a(1,0),a(0,0))*180./PI;
  
    //20001: UK B1
    //30001: Korean B1
    std::vector<std::string>label{"20001","UKB2","UKB3","UKB4","UKB5","30001","KB2","KB3","KB4","KB5"};
    std::vector<cv::Point2f>im_points;
    std::vector<cv::Point3f>light_injectors{cv::Point3f(1490.73,768.14,1302.95),cv::Point3f(1490.73,768.14,666.65),cv::Point3f(1490.73,768.14,30.35),cv::Point3f(1490.73,768.14,-676.65),cv::Point3f(1490.73,768.14,-1312.95),cv::Point3f(1490.73,768.14,1232.25),cv::Point3f(1490.73,768.14,595.95),cv::Point3f(1490.73,768.14,-40.35),cv::Point3f(1490.73,768.14,-605.95),cv::Point3f(1490.73,768.14,-1242.25)};
    li=light_injectors;
    projectPoints(light_injectors,rvec1,tvec1,camera_matrix,dist_coeffs,im_points);
 
    for(int i=0;i<im_points.size();i++){
      int x=im_points[i].x;
      int y=im_points[i].y;
      if(x>0 && x<4000 && y>0 && y-offset<2750){
	cv::circle( img, cv::Point( x, y-offset ), 10, cv::Scalar(0,255,250), -1 );
	cv::putText(img, label[i],cv::Point( x, y-offset ) , cv::FONT_HERSHEY_PLAIN,4, cv::Scalar(0,255,250),3);
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
 
  //relation
  Double_t x[] = {15.17,16.06,14.72,13.39,12.29,10.98,9.82,8.32,7.14,5.94,-2.08,2.18,-6.73,-5.29,4.84,2.25,-11.39};
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

  outs = get_scene(all_pmts,light_injectors, c1, R, ps);
  //  outs=image_to_show = get_scene(all_pmts,li, c1, R_o, ps);
  while(true){
  cv::namedWindow("Scene", cv::WINDOW_NORMAL);    
  cv::resizeWindow("Scene", 1200, 600);
 
  cv::setMouseCallback("Scene", onMouse);
  
  cv::waitKey(0);
  }
  /*
  //img=
  // while(std::cin>>cam1.r>>cam1.theta>>cam1.z>>cam1.yaw>>cam1.pitch>>cam1.roll){
  while(true){
    cv::Mat scene = get_scene(all_pmts,light_injectors, cam1, R, pos);
    //    outs=scene;
    std::cout<<"r, theta, z, yaw, pitch, roll  :  "<<cam1.r<<" , "<<cam1.theta<<" , "<<cam1.z<<" , "<<cam1.yaw<<" , "<<cam1.pitch<<" , "<<cam1.roll<<std::endl;
    
   //imshow("Scene",scene);
    cv::waitKey(1000);
    std::cout<<"Input r theta Z yaw pitch roll :";
      
    
    //clearing buffer
    while (std::cin.get() != '\n') 
      {
	continue;
      }
    


      }*/
  return 0;
}

