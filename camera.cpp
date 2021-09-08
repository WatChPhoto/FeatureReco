#include "camera.hpp"

cam::cam(){
  r=0;
  theta=0;
  z=0;
  yaw=0;
  pitch=0;
  roll=0;
  calculate_R_tv();
}

cam::cam(double r, double theta, double z, double yaw, double pitch, double roll):r(r),theta(theta),z(z),yaw(yaw),pitch(pitch),roll(roll){   
  calculate_R_tv(); 
}

void cam::reset_rvec_tvec(){
  cv::Matx33d R_old;
  cv::Rodrigues(rv_orig,R_old);
  //setting drone orientation 
  yaw=std::atan2(R_old(2,1),R_old(2,0))*180./PI;
  yaw=(yaw>=0)?yaw:360+yaw;
  pitch=std::asin(R_old(2,2))*180./PI;
  roll=std::asin(R_old(0,2)/std::cos(pitch*PI/180.))*180./PI;
  
  //setting drone position
  cv::Matx31d p=-R_old.t()*tv_orig;
  r=std::sqrt(p(0,0)*p(0,0)+p(1,0)*p(1,0));
  theta=std::atan2(p(1,0),p(0,0))*180./PI;
  theta=(theta>=0)?theta:(360+theta);
  z=p(2,0);
  
  calculate_R_tv();
}

void cam::set_camera(double r, double theta, double z, double yaw, double pitch, double roll){
  this->r=r;
  this->theta=theta;
  this->z=z;
  this->yaw=yaw;
  this->pitch=pitch;
  this->roll=roll;
  calculate_R_tv(); 
}
  
void cam::increment(std::string var, float step=1.){
  if(var=="Yaw+"){
    yaw+=step;
  }
  else if(var=="Pih+"){
    pitch+=step;
  }
  else if(var=="Rol+"){
    roll+=step;
  }
  else if(var=="Rad+"){
    r+=step;
  }
  else if(var=="th+"){
    theta+=step;
  }
  else if(var=="Z+"){
      z+=step;
  }
  
  calculate_R_tv();
}

void cam::decrement(std::string var, float step=1.){
  if(var=="Yaw-"){
    yaw-=step;
  }
  else if(var=="Pih-"){
    pitch-=step;
  }
  else if(var=="Rol-"){
    roll-=step;
  }
  else if(var=="Rad-"){
    r-=step;
  }
  else if(var=="th-"){
    theta-=step;
  }
  else if(var=="Z-"){
    z-=step;
  }

  calculate_R_tv();
}

unsigned cam::clip(cv::Point3f n,cv::Matx31d p0,cv::Matx31d p,cv::Matx31d v,cv::Point3f& start,cv::Point3f& end){
  double a = p(0,0)*n.x+p(1,0)*n.y+p(2,0)*n.z;
  double b = (p(0,0)+v(0,0))*n.x+(p(1,0)+v(1,0))*n.y+(p(2,0)+v(2,0))*n.z;
  unsigned count=0;
  if(a>0){count++; start=cv::Point3f(p(0,0),p(1,0),p(2,0));}
  if(b>0){count++; end=cv::Point3f(p(0,0)+v(0,0),p(1,0)+v(1,0),p(2,0)+v(2,0));}
  if(count==1){
    double t=(n.x*(p0(0,0)-p(0,0))+n.y*(p0(1,0)-p(1,0))+n.z*(p0(2,0)-p(2,0)))/(n.x*(v(0,0)-p(0,0))+n.y*(v(1,0)-p(1,0))+n.z*(v(2,0)-p(2,0)));
    cv::Matx31d v1 = (1.0-t)*p+t*v;
    cv::Point3f s = cv::Point3f(v1(0,0),v1(1,0),v1(2,0));
    if(a<=0){end=s;}
    else if(b<=0){start=s;}
  }
  
  return count;
}

cv::Mat cam::get_scene(const cv::Mat& sc){
  cv::Mat scene = sc.clone();
 
  cv::Matx31d rvec;
  cv::Rodrigues(R,rvec);
  cv::Matx31d pos(r*cos(theta*PI/180.),r*sin(theta*PI/180.),z);
  std::vector<cv::Point3f> object_in_view;
  std::vector<cv::Point2f>im_points;
  std::vector<std::string>labels;
  for(int i=0; i<all_pmts.size();i++){
    cv::Matx31d xp = R*cv::Matx31d(all_pmts[i].x,all_pmts[i].y,all_pmts[i].z)+tv;
    if(xp(2,0)>0){
      object_in_view.push_back(cv::Point3f(all_pmts[i].x,all_pmts[i].y,all_pmts[i].z));
      labels.push_back(all_pmts[i].id);
    } 
  }
  
  projectPoints(object_in_view,rvec,tv,camera_matrix,dist_coeffs,im_points);  
  
  for(int i=0;i<im_points.size();i++){
    std::string text=labels[i];
    
    //if the first light injector then use different color
    if(isdigit(text[0])){
      text=text.substr(0,5);
      std::stringstream st(text);
      int x=0;
      st>>x;
      text=std::to_string(x);
      cv::circle( scene, cv::Point( im_points[i].x, im_points[i].y ), 10, cv::Scalar(0,240,240), -1 );
    }
    else{//light injecotr	
      cv::circle( scene, cv::Point( im_points[i].x, im_points[i].y ), 5, cv::Scalar(220,31,237), -1 );
    }
    if(r>800 && r<1750){
      putText(scene, text,cv::Point2f(im_points[i].x,im_points[i].y),cv::FONT_HERSHEY_SIMPLEX,2,cv::Scalar(0,0,240),2);
    }
  }
  //ends here
  
  cv::Matx33d R1;
  cv::Rodrigues(rv_orig,R1);
  cv::Matx31d p = -R1.t()*tv_orig;
  //vector along three axes of camera;
  //let's do 30 cm in each dir
  cv::Matx31d x(30,0,0);
  cv::Matx31d y(0,30,0);
  cv::Matx31d z(0,0,30);
  //endpoints of vector
  x = R1.t()*x+p;
  y = R1.t()*y+p;
  z = R1.t()*z+p;
  
  cv::Point3f start,end;
  std::vector<cv::Scalar> axis_col;
  std::vector<cv::Point3f>obj;
  unsigned counts = clip(cv::Point3f(0,0,1),pos,p,x,start,end);
  if(counts>0){
    obj.push_back(start);
    obj.push_back(end);
    axis_col.push_back(cv::Scalar(0,0,255));
  }

  counts = clip(cv::Point3f(0,0,1),pos,p,y,start,end);
  if(counts>0){
    obj.push_back(start);
    obj.push_back(end);
    axis_col.push_back(cv::Scalar(0,255,0));
  }

  counts = clip(cv::Point3f(0,0,1),pos,p,z,start,end);
  if(counts>0){
    obj.push_back(start);
    obj.push_back(end);
    axis_col.push_back(cv::Scalar(255,0,0));
  }

  std::vector<cv::Point2f>image_points;
  projectPoints(obj,rvec,tv,camera_matrix,dist_coeffs,image_points);    
  
  /* counts = clip(cv::Point3f(1,0,0),cv::Point3f(0,0,0),image_points[0],image_points[1],start,end);
  if(counts>0){
    axis_pts.push_back(start);
    axis_pts.push_back(end);
    final_axis_col.push_back(255,0,0);
  }
  */
  
  //Drawing camera direction
  for(int i=0;i<image_points.size();i+=2){
    cv::line( scene, cv::Point( image_points[i].x, image_points[i].y ),cv::Point( image_points[i+1].x, image_points[i+1].y ), cv::Scalar(0,0,250), 2 );
  }
  /*
  cv::line( scene, cv::Point( image_points[0].x, image_points[0].y ),cv::Point( image_points[1].x, image_points[1].y ), cv::Scalar(0,0,250), 2 );
  cv::line( scene, cv::Point( image_points[0].x, image_points[0].y ),cv::Point( image_points[2].x, image_points[2].y ), cv::Scalar(0,250,0), 2 );
  cv::line( scene, cv::Point( image_points[0].x, image_points[0].y ),cv::Point( image_points[3].x, image_points[3].y ), cv::Scalar(250,0,0), 2 );
  */
  /*
  //drawing sk-coordinate
  std::vector<cv::Point3f> axes={cv::Point3f(0,0,0),cv::Point3f(2000,0,0),cv::Point3f(0,2000,0),cv::Point3f(0,0,2000)};
  
  std::vector<cv::Point2f>axes_points;
  std::vector<cv::Scalar> color;
  
  projectPoints(axes,rvec,tvec,camera_matrix,dist_coeffs,axes_points);
  /*
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
  */
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
  button_clk=false;
  
  return scene;
}
void cam::print(){
    std::cout<<"r = "<<r<<" theta = "<<theta<<" z = "<<z<<" yaw = "<<yaw<<" pitch = "<<pitch<<" roll = "<<roll<<std::endl;
  }

void cam::calculate_R_tv(){
  double thx = pitch*PI/180;
  double thy = -yaw*PI/180;
  double thz = -roll*PI/180;
  
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
  R= Rz.t()*Rx.t()*Ry.t()*Rf1.t()*Rf0.t();
  
  cv::Matx31d pos(r*cos(theta*PI/180.),r*sin(theta*PI/180.),z);
  tv = -R*pos; 
}
 
void cam::set_world_points(std::vector<WorldPoints> p){
  all_pmts=p;
}

bool cam::show_background(){
  return bkg_opt;
}

void cam::flip_background(){
  bkg_opt=!(bkg_opt);
}

//Icon definations 
Icon::Icon(){
  name="";
  pos=cv::Point2f(0,0);
  size=cv::Point2f(100,50);
  text_color=cv::Scalar(0,0,0);
  color=cv::Scalar(0,0,0);
  sign=' ';
}

Icon::Icon(std::string name, cv::Point2f pos, cv::Point2f size, cv::Scalar color, cv::Scalar text_color=cv::Scalar(255,255,255)):name(name),pos(pos),size(size),color(color),text_color(text_color){
  sign=name[name.length()-1];
}

void Icon::show(int x, int y,bool clicked,cv::Mat& img, cam& c1){
  cv::Scalar fill =color;
  cv::Scalar t_fill=text_color;
  if(is_inside(x,y)){
    fill=cv::Scalar(0,250,0); 
    t_fill=cv::Scalar(0,0,0);
    if(clicked){
      c1.button_clk=true;
      fill=cv::Scalar(255,255,255); t_fill=cv::Scalar(0,250,0);
      if(sign=='+'){
	c1.increment(name);
      }
      else if(sign=='-'){
	c1.decrement(name);
      }
      else if(sign==':'){
	if(name=="Set_to:"){
	  // r, theta, z, yaw, pitch, roll.
	  float r, theta, z, yaw, pitch, roll;
	  std::cout<<"enter r theta z yaw pitch roll : ";
	  std::cin>>r>>theta>>z>>yaw>>pitch>>roll;
	  c1.set_camera(r,theta,z,yaw,pitch,roll);
	  //getting rid of token in buffer
	  while (std::cin.get() != '\n') 
	    {
	      continue;
	    }
	}
	if(name=="Show_bkg:"){
	  c1.flip_background();
	}
	if(name=="Org_cam:"){
	  c1.reset_rvec_tvec();
	}
      }
      else{
	if(name=="Save"){
	  std::string filename;
	  std::cout<<"Enter filename to save : ";
	  std::cin>>filename;
	  filename=filename+".jpg";
	  imwrite(filename, img);
	}
	else if(name=="Exit"){
	  std::exit(0);
	  //cv::destroyAllWindows();
	}
	else if(name=="Error"){
	  img;//not complete
	}
      }
    }
  }
  c1.print();
  rectangle(img,pos,pos+size,fill,-1);
  putText(img, name,cv::Point2f(pos.x,pos.y+size.y)+cv::Point2f(10,-10),cv::FONT_HERSHEY_SIMPLEX,2,t_fill,2);
} 

bool Icon::is_inside(int x, int y){
  if(x>=pos.x && x<=pos.x+size.x && y>=pos.y && y<=pos.y+size.y){
    return true;
  }
  else{
    return false;
  }
}

void Icon::set(std::string name, cv::Point2f pos, cv::Point2f size, cv::Scalar color, cv::Scalar text_color=cv::Scalar(255,255,255)){
  this->name=name;
  this->pos=pos;
  this->size=size;
  this->color=color;
  this->text_color=text_color;
  sign=name[name.length()-1];
}
  
