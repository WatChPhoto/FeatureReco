/*=====================================================================================================================
#                                               Fast Ellipse detection code.                                           #
#    Ellipse detection version 2.1.  Different approach of finalizing the ellipse. Looks at all the ellipse formed     # 
#    that include first point and chooses the one with max peak value(max freq) of b (semi minor axis). If the found   #
#    ellipse having max peak value has min_frequency supplied to consider as ellipse then, the point (x1,y1), (x2,y)   #
#    are changed to (-1000, -1000). This will make sure that those both point will not be consider again because no    #
#    matter which point in the image you will consider the distance to (-500,-500) will be too big to consider for     #
#    semi-major axis or semi-minor axis. Also all other points that are within min_minor to max_minor are removed      #
#    from the data(coordinate). Another reason for changing (x1,y1) and (x2,y2) to be(-1000,-100) instead of removing  #
#    them completely is to avoid skipping data due of removal of data at index.
======================================================================================================================*/

#include<cmath>
#include<iostream>
#include "ellipse_detection_2.1.hpp"
#include "distance_to_ellipse.hpp"

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
    max_index = -1;

  for(int i=0; i<bbins.size(); i++){
    if(bbins[i].y>max_freq){max_freq = bbins[i].y; max_index =i;}
  }
}
	
std::vector<ParametricEllipse> detect_ellipse(const std::vector<cv::Vec3f>& input_data, const int min_major, const int max_major, const int min_minor, const int max_minor,int min_minor_freq){

  std::vector<cv::Vec3f> coordinates = input_data;
  std::vector<ParametricEllipse> elData;

  //Get first point
   for(int i=0;i<coordinates.size(); i++){
     int x1 = coordinates[i][0];//[i][0];//x for first point
     int y1 = coordinates[i][1];//[i][1];//y for first point
     //Get second point
     //lets store all the ellipses related to first point.
     ParametricEllipse probable(cv::Point(-1,-1), -1,-1, -1,-1);
     int max_index;
     int max_freq;
     int bmax;
     int sec_ind;
     std::vector<cv::Vec3f> real_unused;
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
	     //f is chord length from (x,y) to second point (x2,y2)
	     double f = distance(cv::Point(x,y), cv::Point(x2,y2));
	     //cos(t), t is the angle made by (x,y) with major axis.
	     double cost = (double)((std::pow(a,2)+std::pow(d,2)-std::pow(f,2)))/(2.0*a*d);
	     double sint = sqrt(1.0-std::pow(cost,2));
	     //semi minor axis
	     int b = abs(sint*a*d)/sqrt((std::pow(a,2)-std::pow(d,2)*std::pow(cost,2)));
	    
	     int index;
	     if(has_key(bbins,b,index)){
	       bbins[index].y +=1;
	     }
	
	     else if(b>=min_minor && b<=max_minor){
	       bbins.push_back(cv::Point2i(b,1));
	     }
	   }
	 }   
		      
	 find_max(bbins, max_freq, max_index);
	 if(max_index>=0){
	   bmax=bbins[max_index].x;
	 }
	 if( max_freq>probable.freq){
	   probable.centre = centre;
	   probable.a =a;
	   probable.b = bmax;
	   probable.alpha =alfa;
	   probable.freq = max_freq;
	   sec_ind=j;
	 } 
       }
     }

     if(probable.freq > min_minor_freq){
       //Now probable is final ellipse 
       elData.push_back(probable);
       
       //unused will include all the unused points
       std::vector<cv::Vec3f> unused;
       
       //including all the points before index to preserve array structure.
       for(int c=0; c<i; c++){ unused.push_back(coordinates[c]);}
       
       //setting the first point to (-1000,-1000) to preserve array structure.
       coordinates[i][0]=-1000;
       coordinates[i][1]=-1000;
      
       //setting the second point to (-1000,-1000) to preserve array structure.
       coordinates[sec_ind][0]=-1000;
       coordinates[sec_ind][1]=-1000;
             
       //removing all the points that were considered for current ellipse.
       for(int m=i+1; m<coordinates.size();m++){
	 int x = coordinates[m][0];
	 int y = coordinates[m][1];
	 int d = distance(cv::Point(x, y), probable.centre);
	
	 // If the points were considered for current ellipse ignore them
	 if(d>=min_minor && d<=max_minor){ continue;}
      	
	 // else keep those points.
	 cv::Vec3f temp;
	 temp[0]=coordinates[m][0];
	 temp[1]=coordinates[m][1];
	 temp[2]=coordinates[m][2];
	 unused.push_back(temp);
       }  
       
       //Data is modified removing used points for current ellipse.
       coordinates = unused;
     }
   }

   //Filling the closest point, dist  and bolts in the pmt info.
   for(int i=0; i< elData.size(); ++i){

     cv::Point centre =  elData[i].centre;
     double e0 = elData[i].a;
     double e1 = elData[i].b;
     double phi = elData[i].alpha;

     for(int j=0; j<input_data.size(); ++j){
       cv::Point p = cv::Point(input_data[j][0], input_data[j][1]); //bolt's centre.
       double d = distance(centre , p); //distance of current blob from centre.
       
       if(d > 2*(e0+e1)){ continue; } // If the point is more than 2X(major axis+minor axis) away then it is too far.
       cv::Point q;
       double dist =  get_distance(centre, e0, e1, p , q, phi); //get the distance of the bolt from current ellipse
       
       if(dist<6){
	 elData[i].bolts.push_back(input_data[j]);
	 elData[i].query.push_back(q);
	 elData[i].dist.push_back(dist);
       }
     }
   }
 
  return elData;
 }

void draw_ellipses(const std::vector<ParametricEllipse>& elData, cv::Mat& img ){   
  const double PI = std::acos(-1);
  for(const ParametricEllipse ellipses: elData){
    cv::Size axes( ellipses.a, ellipses.b );
    cv::Point center = ellipses.centre;
    
    ellipse( img, center, axes, 180.0*ellipses.alpha/PI , 0., 360,  cv::Scalar (255, 255, 255) );
    
    //drawing line from bolts to closest point in the ellipse
    for(int i = 0; i<ellipses.bolts.size(); i++){
      cv::Vec3f bolts = ellipses.bolts[i];
      //drawing included bolts
      cv::circle( img, cv::Point( bolts[0], bolts[1] ), bolts[2], cv::Scalar(0,0,255), 1, 0 );
    
      //drawing line from bolt to closest point in the ellipse
      line(img, cv::Point( cvRound(bolts[0]), cvRound(bolts[1]) ), cv::Point(cvRound(ellipses.query[i].x),cvRound(ellipses.query[i].y)), cv::Scalar(0,0,0), 1, 8,0);

    }
  }
}
  



