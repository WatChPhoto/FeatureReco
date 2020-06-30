#include<cmath>
#include<iostream>
#include "ellipse_detection.hpp"
double distance(cv::Point p1, cv::Point p2){
  int x1 = p1.x;
  int y1 = p1.y;
  int x2 = p2.x;
  int y2 = p2.y;

  //Pythagorean theorem. a^2+b^2=c^2
  return sqrt((double)std::pow((x2-x1),2)+std::pow((y2-y1),2));
}


  //constructor 
ParametricEllipse::ParametricEllipse(cv::Point2i centre, int a, int b, int alpha, int freq  ): centre(centre),a(a),b(b),alpha(alpha),freq(freq){}


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
	
void detect_ellipse(std::vector<cv::Vec3f> coordinates, cv::Mat& img, const int min_major, const int max_major, const int min_minor, const int max_minor,int min_minor_freq){
  std::vector<ParametricEllipse> elData;
  
  //Get first point
  for(int i=0;i<coordinates.size(); i++){
    int x1 = coordinates[i][0];//x for first point
    int y1 = coordinates[i][1];//y for first point
    //Get second point
    for(int j=0;j<coordinates.size();j++){
      int x2 = coordinates[j][0]; //x coordinate of second point
      int y2 = coordinates[j][1];  //y coordinate of second point
      //distance between first point and second point.
      double dist = distance(cv::Point(x1,y1),cv::Point(x2,y2));
      
      //if distance in range of length of major axis 
      if(dist>=2*min_major && dist<=2*max_major){
	std::vector<cv::Point2i> bbins;	
	cv::Point centre = cv::Point(((x1+x2)/2.0),((y1+y2)/2.0));//probable centre of ellipse
	double a = dist/2.0; //length of semi major axis
	double alfa = atan2((y2-y1),(x2-x1));
	//third loop to look for point around centre of ellipse where the distance is in range of minor axis.
	std::vector<cv::Vec3f> unused;
	for(int k=0; k<coordinates.size(); k++){
	  int x = coordinates[k][0];
	  int y = coordinates[k][1];
	  int d = distance(cv::Point(x, y),centre);
	  //if the distance from centre to point is in range of semi-minor axis
	  if(d>=min_minor && d<=max_minor){
	    //f=chord length from (x,y) to second point (x2,y2)
	    double f = distance(cv::Point(x,y), cv::Point(x2,y2));
	    //cos(t), t is the angle made by (x,y) with major axis.
	    double cost = (double)((std::pow(a,2)+std::pow(d,2)-std::pow(f,2)))/(2.0*a*d);
	    double sint = sqrt(1.0-std::pow(cost,2));
	    //semi minor axis
	    int b = abs(sint*a*d)/sqrt((std::pow(a,2)-std::pow(d,2)*std::pow(cost,2)));
	    
	    int index;
	    //std::cout<<" index "<<index<<std::endl;
	    double essentricity = (a>b)?((b+0.0)/a):((a+0.0)/b);
	    if(has_key(bbins,b,index)){
	      bbins[index].y +=1;
	    }
	
	    else if(b>=min_minor && b<=max_minor){
	      bbins.push_back(cv::Point2i(b,1));
	    }
	  }
	  else{cv::Vec3f temp;
	    temp[0]=x;
	    temp[1]=y;
	    temp[2]=coordinates[k][2];
	    unused.push_back(temp);
	  }
	}   
	      int max_index=-1;
	      int max_freq;
	      find_max(bbins, max_freq, max_index);
	      //std::cout<<"max_index"<<max_index<<std::endl;	
	      if(max_index>=0){
		int bmax=bbins[max_index].x;
		if(max_freq>min_minor_freq){// && alfa>=0.0 && bmax>=min_minor){
		  int indx=-1;
		  for(int i=0; i<elData.size(); i++){//ParametricEllipse pel: elData){
		    
		    float d_cent = distance(elData[i].centre, centre);
		    //float theta = atan2(centre.y-elData[i].centre.y, centre.x-elData[i].centre.x);
		    //float r_theta = elData[i].alpha-theta;
		    //float disa = sqrt(std::pow(elData[i].b*cos(r_theta)-elData[i].centre.x,2)+std::pow(elData[i].a*sin(r_theta)-elData[i].centre.y,2));
		    float di = (a>elData[i].a)?a:elData[i].a;
		    //float disb =sqrt(std::pow(bmax*cos(r_theta)-centre.x,2) + std::pow(a*sin(r_theta)-centre.y,2) );
		    if(d_cent<2*di && max_freq > elData[i].freq ){indx =i; break;}
		  } 
		  
		  if(indx==-1){
		    elData.push_back(ParametricEllipse(centre, a, bmax, alfa, max_freq ));
		  }
		  else{ elData[indx]= ParametricEllipse(centre, a, bmax, alfa, max_freq );}
		  coordinates = unused; 
		}
		//bbins.clear();
	      }
      }
    }
  
    //std::cout<<"im in draw exit"<<std::endl;	
  }
  for(ParametricEllipse ellipses: elData){
    
    cv::Size axes( ellipses.a, ellipses.b );
    cv::Point center = ellipses.centre;
    const double PI = std::acos(-1);
    ellipse( img, center, axes, 180.0*ellipses.alpha/PI , 0., 360,  cv::Scalar (255, 255, 255) );


  }
}
  



