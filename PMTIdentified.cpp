#include "PMTIdentified.hpp"

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
    float theta = atan2f((y-b),(x-a));  
    theta = RADTODEG(theta);
    theta = (theta<0)?(theta+360):theta; //getting angle between 0-360
    //Finding angle wrt to Y-axis ^ (nothing to do with axis direction in OpenCv
    theta = (theta<270)?(theta+90):(theta-270);
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
    if ( dmin < bad_dmin ){
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



