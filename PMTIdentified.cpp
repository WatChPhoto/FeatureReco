#include "PMTIdentified.hpp"
#include <map>
#include <TH1D.h>

using namespace cv;
using std::vector;

/*
  Take circles_of_blob to select which bolts are good blobs
  inputs: 
  vector<Vec3f> blobs;  // (x,y,r)
  vector<Vec3f> circles_of_blob;  // (x,y,r)
  
  output:
  vector<Vec3f> final_bolts; // bolt locations selected
*/
// only accept PMTs that are more than some number of pixels away from edge of image
void find_candidate_bolts( const std::vector< Vec3f >& blobs, 
			   const std::vector< Vec3f >& circles_of_blob,
			   std::vector< PMTIdentified >& pmts_found,
			   const Mat& image ) {

  int ywidth = image.rows;
  int xwidth = image.cols;
  int xmin = trim_pixels;
  int xmax = xwidth - trim_pixels;
  int ymin = trim_pixels;
  int ymax = ywidth - trim_pixels;
  
  vector < float > final_dists;

  
  // loop over circles_of_blob 
  // which is the circles of multiple bolts (multiple blobs making circle around PMT)
  for (const Vec3f & pmtloc : circles_of_blob) {
    // loop over bolts (blobs)0 to see if it is on the pmt circle
    //unsigned nbolts = 0;
    int pmtx = pmtloc[0];
    int pmty = pmtloc[1];

    if (pmtx < xmin || pmtx > xmax || pmty < ymin || pmty > ymax) {
      continue;
    }
    vector < Vec3f > bolts_on_this_pmt;
    for (const Vec3f & boltloc:blobs) {
      // calculate distance from the PMT circle to the bolt location
      // only add ones with distance less than (2?) pixels to add to bolts_on_this_pmt
      // count/and print them after
      
      float pmtr = pmtloc[2];
      int boltx = boltloc[0];
      int bolty = boltloc[1];
      float dist =
	std::sqrt (std::pow ((pmtx - boltx), 2) +
		   std::pow ((pmty - bolty), 2));
      
      if (fabs (pmtr - dist) < 6) {
	Vec3f temp;
	temp[0] = boltx;
	temp[1] = bolty;
	temp[2] = boltloc[2];
	final_dists.push_back (fabs (dist - pmtr));
	bolts_on_this_pmt.push_back (temp);
	
      }
    }
    
    // add bolts_on_this_pmt to final_bolts if > some number (5?) of bolts match?
    if (bolts_on_this_pmt.size () > 9) {
      pmts_found.push_back( PMTIdentified( pmtloc, bolts_on_this_pmt, final_dists ));   
    }
  }
}


std::ostream& operator<<( std::ostream& os, const PMTIdentified& p ){
  for ( int ibolt = 0; ibolt < p.bolts.size(); ++ibolt ){
    os<<p.pmtid<<" ";
    os<<p.circ[0]<<" "<<p.circ[1]<<" "<<p.circ[2]<<" ";
    os<<p.boltid[ibolt]<<" ";
    Vec3f b = p.bolts[ibolt];
    os<<b[0]<<" "<<b[1]<<std::endl;    
  }
  return os;
}


void PMTIdentified::calculate_angles(){
  float a = circ[0];  //x and y at centre
  float b = circ[1];
    
  for( Vec3f bolt : bolts ) {
    float x = bolt[0];
    float y = bolt[1];
    float theta = atan2f((x-a),-(y-b)); //getting angle with ^ axis wrt image  
    theta = RADTODEG(theta);
    theta = (theta<0)?(theta+360):theta; //getting angle between 0-360
    
    angles.push_back( theta );
  }
}

void PMTIdentified::calculate_boltid(){
  float angle_between_bolts = 360.0 / 24; // 24 bolts
  float dang = angle_between_bolts/2;

  for ( float angle : angles ){
    int boltnum = int( (angle+dang) / angle_between_bolts ) + 1;
    if (boltnum==25) boltnum=1;
    boltid.push_back( boltnum );
    // calculate difference in angle from boltid angle
    float boltidang = (boltnum-1) * angle_between_bolts;
    float da = angle - boltidang;
    if ( da > 360.0-dang ) da -= 360.0;
    dangs.push_back( da );
  }
}




/// find_closest_matches finds closest match to each entry in mtd to the closest entry in circles
/// without using any circle twice.  Also maybe reject bolts outside expected angle?
//Returns index of circle, mtd and mindist.
void find_closest_matches( std::vector< PMTIdentified>& final_pmts, const MedianTextData & mtd ){

  // initialize idx_txt and dist_txt in PMTIdentifieds
  for ( PMTIdentified & pmt : final_pmts ){
    pmt.idx_txt.clear();
    pmt.dist_txt.clear();
    for ( unsigned i=0; i<pmt.bolts.size(); ++i ){
      pmt.idx_txt.push_back( 0 );
      pmt.dist_txt.push_back( bad_dmin );
    }
    assert( pmt.idx_txt.size() == pmt.bolts.size() );
    assert( pmt.dist_txt.size() == pmt.bolts.size() );
  } 

  
  for ( unsigned mtr_idx = 0; mtr_idx < mtd.size(); ++ mtr_idx ){
    const MedianTextRecord & rec = mtd[ mtr_idx ];
    unsigned trueboltnum = rec.bolt_num();

    float dmin = bad_dmin;
    size_t dmin_pmtidx = 0;
    size_t dmin_boltidx = 0;
    size_t idx_pmtidx = 0;
    for ( unsigned pmt_idx = 0; pmt_idx < final_pmts.size(); ++pmt_idx ){
      const PMTIdentified& pmtfound = final_pmts[ pmt_idx ];
      for ( unsigned bolt_idx = 0; bolt_idx < pmtfound.bolts.size(); ++bolt_idx ){
	const cv::Vec3f& circ = pmtfound.bolts[ bolt_idx ];
	float dist = std::sqrt( (circ[0] - rec.x())*(circ[0] - rec.x()) +
				(circ[1] - rec.y())*(circ[1] - rec.y()) );
	bool boltanglematch = ( pmtfound.boltid[ bolt_idx ] == trueboltnum );
	if ( dist < dmin && boltanglematch ){
	  dmin = dist;
	  dmin_pmtidx = pmt_idx;
	  dmin_boltidx = bolt_idx;
	}
      }
    }

    // have closest match to MedianTextRecord?
    if ( dmin < 20){//bad_dmin ){
      final_pmts[ dmin_pmtidx ].idx_txt[ dmin_boltidx ] = mtr_idx;
      final_pmts[ dmin_pmtidx ].dist_txt[ dmin_boltidx ] = dmin;
    }
    
  }
  

  // Now get rid of duplicates?
  // for now try with duplicates....
}

//Makes the histogram of minimum distance from the matched bolt from text to found bolt
void make_bolt_dist_histogram( const std::vector< PMTIdentified > & final_pmts, TH1D *&hout){
  for ( const PMTIdentified& pmtmatch : final_pmts ){
    for ( const float& dist : pmtmatch.dists ){
      if ( dist < bad_dmin ){
	hout->Fill( dist );
      }
    }
  }
}



void prune_bolts( std::vector< PMTIdentified >& final_pmts, float ang_offset ){
  float angle_between_bolts = 360.0 / 24; // 24 bolts
  float dang = angle_between_bolts/2;


  for( PMTIdentified& pmt : final_pmts ){
    std::map<  unsigned, std::vector< unsigned > > boltmap; 
    for (unsigned i=0; i<pmt.bolts.size() ; ++i ){
     boltmap[ pmt.boltid[i] ].push_back( i );
    }
    std::vector<unsigned> indices_to_remove;
    for ( const std::pair<const unsigned int, std::vector<unsigned int> >& key : boltmap ){
      unsigned boltnum = key.first;
      if ( key.second.size() > 0 ){
	float mindist = 10000.0;
	unsigned idxkeep = 0;
	float boltidang = (boltnum-1) * angle_between_bolts + ang_offset;
	for ( unsigned idx : key.second ){
	  // compare angle to expected 
	  float da = pmt.angles[ idx ] - boltidang;
	  if ( da > 360.0-dang ) da -= 360.0;
	  if ( da < mindist ){
	    mindist = da;
	    idxkeep = idx;
	  }
	} 
	for ( unsigned idx : key.second ){
	  if ( idx != idxkeep || fabs( mindist )>=4.0 ){ // only keep best if it is within +-4 degrees of expected
	    indices_to_remove.push_back( idx );
	  }
	}

      }	
    }

    // now remove bolts
    std::vector<cv::Vec3f> bolts; // bolts going with this PMT
    std::vector<float> dists; // distance of bolt from PMT circle 
    std::vector<float> angles; // angle of each bolt
    std::vector<float> dangs;  // difference in angle from boltid's angle
    std::vector<int>   boltid; // 1 is at 12 o'clock, 2 ... 24 going around clockwise

    // include comparison to truth if available
    std::vector<int>   idx_txt; // index of medianTextReader
    std::vector<float> dist_txt; // distance to closest matching bolt

    for (unsigned i=0; i<pmt.bolts.size() ; ++i ){
      if ( std::find( indices_to_remove.begin(), 
		      indices_to_remove.end(), i ) == indices_to_remove.end() ){

	bolts.push_back( pmt.bolts[ i ] );
	dists.push_back( pmt.dists[ i ] );
	float cor_ang = pmt.angles[ i ] - ang_offset;
	if ( cor_ang > 360.0 ) cor_ang-=360.0;
	angles.push_back( cor_ang );
	dangs.push_back( pmt.dangs[ i ] - ang_offset );
	boltid.push_back( pmt.boltid[ i ] );
	if ( pmt.idx_txt.size() > 0 ){
	  idx_txt.push_back( pmt.idx_txt[ i ] );
	  dist_txt.push_back( pmt.dist_txt[ i ] );
	}
      }

    }

    pmt.bolts = bolts;
    pmt.dists = dists;
    pmt.angles = angles;
    pmt.boltid = boltid;
    pmt.idx_txt = idx_txt;
    pmt.dist_txt = dist_txt;
  }

  
}


void prune_pmts(  std::vector< PMTIdentified >& final_pmts, unsigned numbolts, const std::string& label ){

  std::vector< PMTIdentified > not_pruned_pmts;


  for ( unsigned i = 0; i < final_pmts.size(); ++i ){ 
    const PMTIdentified& pmt = final_pmts[i];
    bool isinside=false;
    bool hasfewerbolts=false; // only set for intersecting pmts
    for ( unsigned j = i+1; j < final_pmts.size(); ++j ){
      const PMTIdentified& pmtb =  final_pmts[j];
      //if ( i == j ) continue;
      float r1 = pmt.circ[2];
      float r2 = pmtb.circ[2];
      float r1r2 = r1+r2;
      float x1 = pmt.circ[0], y1 = pmt.circ[1];
      float x2 = pmtb.circ[0], y2 = pmtb.circ[1];
      float dist = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) );
      if ( dist < r1 && 
	   r1 < pmtb.circ[2] ) {



	  
	isinside = true;
      } else {
	if ( dist < r1r2 && 
	     pmt.bolts.size() < pmtb.bolts.size() ) {

	  std::cout<<"hasfewerbolts i="<<i<<" j="<<j<<std::endl;
	  std::cout<<"(x1,y1)= ("<<x1<<", "<<y1<< ") r1="<<r1
		   <<" (x2,y2)="<<x2<<", "<<y2<< ") r2="<<r2
		   <<" dist ="<<dist<<" r1r2="<<r1r2<<std::endl;

	  hasfewerbolts = true;
	}
      }
    }
      
    if ( pmt.bolts.size() >= numbolts 
	 //	 || !isinside 
	 	 && !hasfewerbolts 
	 ){
      not_pruned_pmts.push_back( pmt );
    }
  }

  final_pmts = not_pruned_pmts;

  // make some histograms of distances

  std::ostringstream os_name;
  os_name << "prunepmt_closest_" << label;
  TH1D* hdist = new TH1D( os_name.str().c_str(), " ; distance in PMT radii; count/bin", 200, 0., 5. ); 
  
  std::ostringstream os_name2;
  os_name2 << "prunepmt_all_" << label;
  TH1D* hdist_all = new TH1D( os_name2.str().c_str(), " ; distance in PMT radii; count/bin", 200, 0., 5. ); 


  for ( unsigned i = 0; i < final_pmts.size(); ++i ){ 
    const PMTIdentified& pmta = final_pmts[i];
    float closest_dist = 9999999.0;
    for ( unsigned j = i+1; j < final_pmts.size(); ++j ){
      const PMTIdentified& pmtb =  final_pmts[j];
      //if ( i == j ) continue;
      //float r = pmta.circ[2];
      float r1 = pmta.circ[2]; float r2 = pmtb.circ[2];
      float r = (r1+r2)/2;
      float x1 = pmta.circ[0], x2 = pmtb.circ[0];
      float y1 = pmtb.circ[0], y2 = pmtb.circ[1];
      float dist = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) )/r;
      if ( fabs( x1 - x2 ) < r || fabs( y1 - y2 ) < r ) { 
	hdist_all->Fill( dist );
      }
      if ( dist < closest_dist ) closest_dist = dist;
    }
    hdist->Fill( closest_dist );
  }
			 
}




void overlay_bolt_angle_boltid(const std::vector< PMTIdentified >& final_pmts, cv::Mat image_final){
  for( const PMTIdentified& pmt : final_pmts ){
    float a = pmt.circ[0]; //x-coordinate of centre of pmt
    float b = pmt.circ[1];
    for(int i=0; i<pmt.bolts.size();i++){
      std::string txt = std::to_string((int)pmt.angles[i]);
      // A(x1,y1)         P(x,y)            B(x2,y2) 
      //  o-----------------o-----------------o
      //  AP:PB = m:n
      float m = 4;
      float n =1;
      float x = pmt.bolts[i][0];
      float y = pmt.bolts[i][1];
      cv::Point text_at = Point((m*x + n*a)/(m+n),(m*y + n*b)/(m+n));

      //writing bolt angle from ^ in image 
      cv::putText( image_final, txt, text_at, FONT_HERSHEY_DUPLEX, 0.3, cv::Scalar(0,255,0),1);

      txt = std::to_string((int)pmt.boltid[i]);
      m = 3;
      n =2;
      x = pmt.bolts[i][0];
      y = pmt.bolts[i][1];
      text_at = Point((m*x + n*a)/(m+n),(m*y + n*b)/(m+n));

      //writing boltid in image
      cv::putText( image_final, txt, text_at, FONT_HERSHEY_DUPLEX, 0.3, cv::Scalar(0,255,255),1);
    }
  }
}

