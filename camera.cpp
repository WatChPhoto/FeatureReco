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

void find_intersect(cv::Matx31d n, cv::Matx31d p, cv::Matx31d u, cv::Matx31d v, cv::Matx31d& x){
  double t = n.dot(p-u)/n.dot(v-u);
  x=u+t*(v-u);
}

//If part of object falls outside the image boundary then it givis garbage result. 
//To avoid this clipping is necessary.
unsigned cam::clip(cv::Matx31d n,cv::Matx31d p,cv::Matx31d pt[2], cv::Point3f& start,cv::Point3f& end){ 
  //checking number of points inside the plane and outside the plane.
  //inside if in the direction of normal otherwise outside.
  unsigned in[2]; //point
  int n_in=0;
  for(unsigned i=0;i<2;i++){
    cv::Matx31d u=pt[i]-p;
    float n_dot_u = n(0,0)*u(0,0)+n(1,0)*u(1,0)+n(2,0)*u(2,0);
    if(n_dot_u>=0){in[i]=1; n_in++;}
    else{in[i]=0;}
  }
  
  //If two points are inside no clipping needed. return original points.
  if(n_in==2){
    start=cv::Point3f(pt[0](0,0),pt[0](1,0),pt[0](2,0));
    end=cv::Point3f(pt[1](0,0),pt[1](1,0),pt[1](2,0));
    return n_in;
  }
  //If no points are inside then nothing is returned, since nothing is visible.
  else if(n_in==0){return n_in;}
  //If one point is inside and another is outside then find the intersection with the plane(boundary).
  else if(n_in==1){
    //find the intersection point.
    cv::Matx31d x;
    find_intersect(n,p,pt[0],pt[1],x);
    if(in[0]){
      start=cv::Point3f(pt[0](0,0),pt[0](1,0),pt[0](2,0));
      //set intersection to it to end;
      end=cv::Point3f(x(0,0),x(1,0),x(2,0));
      return n_in;
    }
    else{
      end=cv::Point3f(pt[1](0,0),pt[1](1,0),pt[1](2,0));
      //set intersection to it to start;
      start=cv::Point3f(x(0,0),x(1,0),x(2,0));
      return n_in;
    }
  }

}

//x,y,z axes in world coordinate
//p is origin of axes
//pos is camera position in world coordinate
//axis_col returns color of each axes
//image_points_f is final transformed points
void cam::get_transformed_axes(cv::Matx31d x, cv::Matx31d y, cv::Matx31d z, cv::Matx31d p, cv::Matx31d pos, std::vector<cv::Scalar>& axis_col,  std::vector<cv::Point2f>& image_points_f, cv::Size2i s){
  
  cv::Matx31d rvec;
  cv::Rodrigues(R,rvec);//find rotation vector for our camera.
  cv::Point3f start,end;
  std::vector<cv::Point3f>obj;
  cv::Matx31d lin[2];
  lin[0]=p;
  lin[1]=x;
  cv::Matx31d normal; //z direction of our camera in world coordinate
  normal = R.t()*cv::Matx31d(0,0,1);
  //clipping aginst the z axis of the camera
  unsigned counts = clip(normal,pos+normal,lin,start,end);
  if(counts>0){
    obj.push_back(start);
    obj.push_back(end);
    axis_col.push_back(cv::Scalar(0,0,255));
  }
  
  lin[1]=y;
  counts = clip(normal,pos+normal,lin,start,end);
  if(counts>0){
    obj.push_back(start);
    obj.push_back(end);
    axis_col.push_back(cv::Scalar(0,255,0));
  }

  lin[1]=z;
  counts = clip(normal,pos+normal,lin,start,end);
  if(counts>0){
    obj.push_back(start);
    obj.push_back(end);
    axis_col.push_back(cv::Scalar(255,0,0));
  }

  std::vector<cv::Point2f>image_points;
  if(obj.size()>0){
    projectPoints(obj,rvec,tv,camera_matrix,dist_coeffs,image_points);    
  }

  float w=s.width;
  float h=s.height;
  //clipping in image space aginst left, right, top, and bottom of image.
  for(unsigned i=0;i<image_points.size();i+=2){
    lin[0]=cv::Matx31d(image_points[i].x,image_points[i].y,0);
    lin[1]=cv::Matx31d(image_points[i+1].x,image_points[i+1].y,0);
    counts = clip(cv::Matx31d(1,0,0),cv::Matx31d(0,0,0),lin,start,end); //aginst x=0 plane
    if(counts>0){
      lin[0]=cv::Matx31d(start.x,start.y,0);
      lin[1]=cv::Matx31d(end.x,end.y,0);
      counts = clip(cv::Matx31d(-1,0,0),cv::Matx31d(w,0,0),lin,start,end); //aginst x=image width plane
      if(counts>0){
	lin[0]=cv::Matx31d(start.x,start.y,0);
	lin[1]=cv::Matx31d(end.x,end.y,0);
	counts = clip(cv::Matx31d(0,1,0),cv::Matx31d(0,0,0),lin,start,end); //aginst y=0 plane
	if(counts>0){
	  lin[0]=cv::Matx31d(start.x,start.y,0);
	  lin[1]=cv::Matx31d(end.x,end.y,0);
	  counts = clip(cv::Matx31d(0,-1,0),cv::Matx31d(0,h,0),lin,start,end); //aginst y=image height plane
	  
	  if(counts>0){
	    image_points_f.push_back(cv::Point2f(start.x,start.y));
	    image_points_f.push_back(cv::Point2f(end.x,end.y));
	  }
	}
      }
    }
    
  }
}

void show_error(cv::Mat& m, unsigned err){
  cv::rectangle(m,cv::Point2f(0,120),cv::Point2f(1000,300),cv::Scalar(255,255,255),-1);
  std::string text ="Err = "+ std::to_string(err);
  cv::putText(m,text,cv::Point(10,280),cv::FONT_HERSHEY_SIMPLEX,4,cv::Scalar(0,0,0),10);

}

void cam::draw_rep_error( const std::vector<cv::Point2f>& im_points, cv::Mat& m, float offset=250){
  double sum=0;
  for(int i=0; i<all_ellipses.size();i++){
    double x0 = all_ellipses[i].x;
    double y0 = all_ellipses[i].y+offset;
    double lmin = 1000000000000000000000000000;
    cv::Point2f p2;
    for(int j=0;j<im_points.size();j++){
      double x1 = im_points[j].x;
      double y1 = im_points[j].y;//-offset;
      double d = (x1-x0)*(x1-x0)+(y1-y0)*(y1-y0);
      if(d<lmin){lmin=d; p2=cv::Point2f(x1,y1);}
    }
    sum+=lmin;
    //draw line from ellipses to points 
    cv::line(m,cv::Point2f(x0,y0),p2,cv::Scalar(0,255,0),4);//cv::Scalar(52,128,235),2);
    cv::circle(m,cv::Point2f(x0,y0),10,cv::Scalar(255, 102, 255),-1);
    cv::Size axes(  int(all_ellipses[i].get_a()), int(all_ellipses[i].b) );
    cv::ellipse(m,cv::Point2f(x0,y0),axes,all_ellipses[i].phi*PI/180.,0,360,cv::Scalar (255, 102, 255),1);
  }
  
  show_error(m,round(sum));
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
  //x  draw_error=true;
  //draw error
  if(show_rep_error){
    draw_rep_error( im_points, scene, 250);
  }
  //ends here
  
  //From here on Drawing axes is performed.
  /*###################Drawing Orientation of original drone while the picture was taken############################*/
  cv::Matx33d R1;
  cv::Rodrigues(rv_orig,R1);
  cv::Matx31d p = -R1.t()*tv_orig;
  //vector along three axes of camera;
  //let's do 30 cm in each dir
  cv::Matx31d x(10,0,0);
  cv::Matx31d y(0,10,0);
  cv::Matx31d z(0,0,10);
  //endpoints of vector in world coordinate
  x = R1.t()*x+p;
  y = R1.t()*y+p;
  z = R1.t()*z+p;
  
  cv::Point3f start,end;
  std::vector<cv::Scalar> axis_col;
  std::vector<cv::Point2f>axes_points;
  get_transformed_axes(x, y, z, p, pos, axis_col, axes_points,sc.size());

  //Drawing camera direction
  for(int i=0;i<axes_points.size();i+=2){
    cv::arrowedLine( scene, cv::Point( axes_points[i].x, axes_points[i].y ),cv::Point( axes_points[i+1].x, axes_points[i+1].y ), axis_col[i/2], 2 );
  
  }
  /*$$$$$$$$$$$$$$$$ Done drawing orientation of original drone $$$$$$$$$$$$$$$$*/
 
  
 /*###################Drawing Sk axes############################*/
  x=cv::Matx31d(2000,0,0);
  y=cv::Matx31d(0,2000,0);
  z=cv::Matx31d(0,0,2000);
  p=cv::Matx31d(0,0,0);
  axis_col.clear();
  axes_points.clear();
  get_transformed_axes(x, y, z, p, pos, axis_col, axes_points,sc.size());
  
  for(int i=0;i<axes_points.size();i+=2){
    cv::arrowedLine( scene, cv::Point( axes_points[i].x, axes_points[i].y ),cv::Point( axes_points[i+1].x, axes_points[i+1].y ), axis_col[i/2], 2 );
    
  }
  /*###################### Drawing SK axes ends here ######################################*/
  
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
 
void cam::set_world_points(std::vector<WorldPoints>& p){
  all_pmts=p;
}

void cam::set_ellipses(std::vector<Ellipse>& e){
  all_ellipses=e;
}

bool cam::show_background(){
  return bkg_opt;
}

void cam::flip_background(){
  bkg_opt=!(bkg_opt);
}

void cam::flip_error(){
  show_rep_error=!(show_rep_error);
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
	  c1.flip_error();
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
  
