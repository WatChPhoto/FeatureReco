#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <fstream>                
#include <cmath>
#include <iostream>
#include <vector>
#include <opencv2/calib3d.hpp>
#include "distancelib.hpp"
//#include "camera.hpp"


cv::Mat outs;
cv::Mat scene;
cv::Mat bkg_scene;
cam c1;
cv::Mat image_to_show;
int img_no=0;
std::vector<Icon>icons(15);
std::vector<Icon> icons2(2);

void onMouse(int event, int x, int y, int flags, void* userdata){
  bool clicked=false;
  if(event==cv::EVENT_LBUTTONDOWN){
    clicked=true;
  }

  if(c1.button_clk){
    if(!c1.show_background()){
      outs = c1.get_scene(scene);
    }
    else{
      outs = c1.get_scene(bkg_scene);
    }

    std::string nme = std::to_string(img_no)+".jpg";
    std::cout<<"name = "<<nme<<std::endl;
    // imwrite(nme,image_to_show);
    img_no++;
  } 
  
  image_to_show=outs.clone();
  
  int xmin = image_to_show.cols-900;
  int ymin = image_to_show.rows-800;
  int xmax = image_to_show.cols;
  int ymax = image_to_show.rows;
  if(x>xmin && y>ymin && x<xmax && y<ymax){

    for(int i=0;i<icons.size();i++){
      icons[i].show(x,y,clicked,image_to_show,c1);
    }   
  }

 int xmin2 = scene.cols-500;
 int ymin2=0;
 if(x>xmin2 && y>ymin2 && x<xmax && y<100){
   for(int i=0;i<icons2.size();i++){
     icons2[i].show(x,y,clicked,image_to_show,c1);
   }
 }
 
 if(x>0 && y>0 &&x<2500 && y<100){
   std::string text;
   cv::rectangle(image_to_show,cv::Point2f(0,0),cv::Point2f(2500,100),cv::Scalar(255,255,255),-1);
   text ="r="+std::to_string(c1.r).substr(0,5)+" theta="+std::to_string(c1.theta).substr(0,5)+" Z="+std::to_string(c1.z).substr(0,5)+" Yaw="+std::to_string(c1.yaw).substr(0,5)+" Pitch="+std::to_string(c1.pitch).substr(0,5)+" roll="+std::to_string(c1.roll).substr(0,5);
   cv::putText(image_to_show,text,cv::Point(10,90),cv::FONT_HERSHEY_SIMPLEX,2,cv::Scalar(0,0,0),2.5);
 }

  
  imshow("Scene",image_to_show);
}

void read_R_t(std::string filename, cv::Matx31d &r, cv::Matx31d &t){
  std::string line;
  double a,b,c;
  std::ifstream file(filename);
  for(int i=0;i<2;i++){
    getline(file,line);
    unsigned begin=line.find('[');
    unsigned end=line.find(']');
    line=line.substr(begin+1,end-begin);
    std::cout<<"r/t = "<<line<<std::endl;
    std::stringstream ss;
    ss<<line;
    ss>>a>>b>>c;
    if(i==0){r={a,b,c};}
    else{t={a,b,c};}
  }
  file.close();
}

void set_icons(std::vector<Icon>& i,std::vector<Icon>& i2 ){
  int xmin = scene.cols-900;
  int ymin = scene.rows-800;
  int xmax = scene.cols;
  int ymax = scene.rows;
  int xoffset=100;
  int yoffset=50;
  
  cv::Point2f size = cv::Point2f(200,100);
  i[0].set( "Rad+",cv::Point2f(xmin,ymin),size,cv::Scalar(0,0,0),cv::Scalar(255,255,255));
  i[1].set("Rad-",cv::Point2f(xmin,ymin+4*(size.y+yoffset)),size,cv::Scalar(0,0,0),cv::Scalar(255,255,255));
  i[2].set( "th+",cv::Point2f(xmin+size.x+xoffset,ymin),size,cv::Scalar(0,0,0),cv::Scalar(255,255,255));
  i[3].set( "th-",cv::Point2f(xmin+size.x+xoffset,ymin+4*(size.y+yoffset)),size,cv::Scalar(0,0,0),cv::Scalar(255,255,255));
  i[4].set( "Z+",cv::Point2f(xmin+2*(size.x+xoffset),ymin),size,cv::Scalar(0,0,0),cv::Scalar(255,255,255));
  i[5].set( "Z-",cv::Point2f(xmin+2*(size.x+xoffset),ymin+4*(size.y+yoffset)),size,cv::Scalar(0,0,0),cv::Scalar(255,255,255));
  i[6].set( "Yaw-",cv::Point2f(xmin,ymin+size.y+yoffset),size,cv::Scalar(0,0,0),cv::Scalar(255,255,255));
  i[7].set( "Yaw+",cv::Point2f(xmin+2*(size.x+xoffset),ymin+size.y+yoffset),size,cv::Scalar(0,0,0),cv::Scalar(255,255,255));
  i[8].set( "Pih-",cv::Point2f(xmin,ymin+2*(size.y+yoffset)),size,cv::Scalar(0,0,0),cv::Scalar(255,255,255));
  i[9].set( "Pih+",cv::Point2f(xmin+2*(size.x+xoffset),ymin+2*(size.y+yoffset)),size,cv::Scalar(0,0,0),cv::Scalar(255,255,255));
  i[10].set( "Rol-",cv::Point2f(xmin,ymin+3*(size.y+yoffset)),size,cv::Scalar(0,0,0),cv::Scalar(255,255,255));
  i[11].set( "Rol+",cv::Point2f(xmin+2*(size.x+xoffset),ymin+3*(size.y+yoffset)),size,cv::Scalar(0,0,0),cv::Scalar(255,255,255));
  i[12].set( "Set_to:",cv::Point2f(xmin+size.x+xoffset/2,ymin+3*size.y+3*yoffset),cv::Point2f(300,110),cv::Scalar(0,0,0),cv::Scalar(255,255,255));
  i[13].set( "Show_bkg:",cv::Point2f(xmin+size.x+xoffset/2,ymin+size.y+yoffset/2),cv::Point2f(300,110),cv::Scalar(0,0,0),cv::Scalar(255,255,255));
i[14].set( "Org_cam:",cv::Point2f(xmin+size.x+xoffset/2,ymin+2*size.y+2*yoffset),cv::Point2f(300,110),cv::Scalar(0,0,0),cv::Scalar(255,255,255));

  
 int xmin2 = scene.cols-500;
 int ymin2=0;
 i2[0].set( "Save",cv::Point2f(xmin2,ymin2),size,cv::Scalar(0,0,0),cv::Scalar(255,255,255));
 i2[1].set( "Exit",cv::Point2f(xmin2+size.x+xoffset,ymin2),size,cv::Scalar(0,0,0),cv::Scalar(255,255,255));
}

int main(int argc, char **argv){
  //read image.
  std::string filename = std::string(argv[argc-1]);
  //std::string img_num=filename
  cv::Mat image = imread (filename+".jpg", cv::IMREAD_COLOR);//("045.jpg", cv::IMREAD_COLOR);
  
  //read all pmts
  std::vector<WorldPoints> all_pmts = read_all_world_points();

  scene = cv::Mat(3000, 4000, CV_8UC3, cv::Scalar(255,255,255));   
  bkg_scene=cv::Mat(3000, 4000, CV_8UC3, cv::Scalar(255,255,255));
  set_icons(icons,icons2);
  cv::Matx31d r;
  cv::Matx31d t;
  //239_rot_trans.txt 
  std::string text_file =filename+"_rot_trans.txt";
  read_R_t(text_file,r,t);
  //cv::Rodrigues(r,R);
  //p=-R.t()*t;
  c1.rv_orig=r;
  c1.tv_orig=t;
  c1.set_world_points(all_pmts);

  std::cout<<"r_old"<<c1.rv_orig(0,0)<<" , "<<c1.rv_orig(1,0)<<" , "<<c1.rv_orig(2,0)<<std::endl;
  std::cout<<"t_old"<<c1.tv_orig(0,0)<<" , "<<c1.tv_orig(1,0)<<" , "<<c1.tv_orig(2,0)<<std::endl;
  

  
  for(unsigned i=250; i<bkg_scene.rows;i++){
    for(unsigned j=0; j<bkg_scene.cols;j++){
      bkg_scene.at<cv::Vec3b>(i,j)=image.at<cv::Vec3b>(i-250,j);
    }
  }
  
  
  outs = c1.get_scene(scene);
  //  outs=image_to_show = get_scene(all_pmts,li, c1, R_o, ps);
    while(true){
      cv::namedWindow("Scene", cv::WINDOW_NORMAL);    
      cv::resizeWindow("Scene", 1200, 600);
      cv::setMouseCallback("Scene", onMouse);
      cv::waitKey(0);
    }    

    return 0;
}
