#include<cmath>
//#include<opencv2/core.hpp>
//#include<opencv2/imgproc.hpp>
#include<iostream>
#include "ellipse_detection.hpp"
double distance(cv::Point p1, cv::Point p2){
  int x1 = p1.x;
  int y1 = p1.y;
  int x2 = p2.x;
  int y2 = p2.y;

  //Pythagorean theorem. a^2+b^2=c^2
  return sqrt(std::pow(((x2-x1)+0.0),2)+std::pow(((y2-y1)+0.0),2));
}

class ParametricEllipse{
public:
  cv::Point2i centre;
  int a;
  int b;
  int alpha;
  
ParametricEllipse(cv::Point2i centre, int a, int b, int alpha  ): centre(centre),a(a),b(b),alpha(alpha){}

};

bool has_key(const std::vector<cv::Point2i>& bbins, int b, int& index){
  for(int i=0; i<bbins.size(); ++i){
    if(bbins[i].x == b){index =i; return true;}
  }
  return false;
}

void find_max(const std::vector<cv::Point2i>& bbins, int& max_freq, int& max_index){ 
  max_freq = -1;
  //max_index = 0;

  for(int i=0; i<bbins.size(); i++){
    if(bbins[i].y>max_freq){max_freq = bbins[i].y; max_index =i;}
  }
}
	
void detect_ellipse(const std::vector<cv::Vec3f>& coordinates, cv::Mat& img, const int min_major, const int max_major, const int min_minor, const int max_minor,int min_minor_freq){
  std::vector<ParametricEllipse> elData;
  for(int i=0;i<coordinates.size(); i++){
    cv::Vec3f ordered_triple = coordinates[i];
    int x1 = ordered_triple[0];
    int y1 = ordered_triple[1];
    //cv::Point centre=cv::Point(-1,-1);
    //double alfa=0;
    //double a=0;

    for(int j=0;j<coordinates.size();j++){
      cv::Vec3f sec_ordered_triple = coordinates[j];
      int x2 = sec_ordered_triple[0];
      int y2 = sec_ordered_triple[1];
      double dist = distance(cv::Point(x1,y1),cv::Point(x2,y2));
      
      if(dist>=2*min_major && dist<=2*max_major){
	std::vector<cv::Point2i> bbins;	
	cv::Point centre = cv::Point(((x1+x2)/2.0),((y1+y2)/2.0));
	double a = (dist+0.0)/2.0; //length of semi major axis
	double alfa = atan2((y2-y1),(x2-x1));

	for(int k=0; k<coordinates.size(); k++){
	  cv::Vec3f third_order_triple = coordinates[k];
	  int rx = third_order_triple[0];
	  int ry = third_order_triple[1];
	  int d = distance(cv::Point(rx, ry),centre);
	  if(d>=min_minor && d<=max_minor){
	    double f = distance(cv::Point(rx,ry), cv::Point(x2,y2));
	    double cost = ((std::pow(a,2)+std::pow(d,2)-std::pow(f,2))+0.0)/(2.0*a*d);//(0.00001+2.0*a*d);
	    //length of minor axis
	    int b = cvRound(sqrt(abs((std::pow((a+0.0),2)*std::pow((d+0.0),2)*(1.0-std::pow(cost,2))/(std::pow((a+0.0),2)-std::pow((d+0.0),2)*std::pow(cost,2))))));//(0.00001+ std::pow((a+0.0),2)-std::pow((d+0.0),2)*cost*2)))));
	    int index;
	    std::cout<<" index "<<index<<std::endl;
	    double essentricity = (a>b)?((b+0.0)/a):((a+0.0)/b);
	    if(has_key(bbins,b,index)){
	      bbins[index].y +=1;
	    }
	
	    else if(essentricity>0.8){
	      bbins.push_back(cv::Point2i(b,1));
	    }
	    
	      int max_index=-1;
	      int max_freq;
	      find_max(bbins, max_freq, max_index);
	      std::cout<<"max_index"<<max_index<<std::endl;	
	      if(max_index>=0){
		int bmax=bbins[max_index].x;
		if(max_freq>min_minor_freq){// && alfa>=0.0 && bmax>=min_minor){
		  elData.push_back(ParametricEllipse(centre, a, bmax, alfa));
		}
		//bbins.clear();
	      }
	  }
	

	} 
      }
  
      
    }
   std::cout<<"im in draw exit"<<std::endl;	
  }
  for(ParametricEllipse ellipses: elData){
    
    cv::RotatedRect rRect = cv::RotatedRect(ellipses.centre, cv::Size2f(ellipses.a,ellipses.b), ellipses.alpha);
      ellipse(img, rRect, cv::Scalar(255,255,255), 3);

  }
}
  



