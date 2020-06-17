
#include "MedianTextReader.hpp"


MedianTextRecord::MedianTextRecord( unsigned num, unsigned pmt, unsigned bolt, unsigned x, unsigned y, std::string id ) : fNum( num ), fPMT( pmt ), fBolt( bolt ), fX( x) , fY( y ), fID( id ) { }

MedianTextRecord::MedianTextRecord() : fNum( 0 ), fPMT( 0 ), fBolt( 0 ), fX( 0) , fY( 0 ), fID("" ) { }

bool MedianTextRecord::is_bolt(){
  
  /*
    ## Numbering Convention
    
  20001: UK B1
  00: Center of diffuser ball
  
  30001: Korean B1
  00-07: Bolts clockwise from top-left

  <19999: PMTs
   00: Center of dynode (or light reflection for labels PD1/2)
   01-24: Bolts clockwise from top center (+z)
   25: Centroid of light reflection nearest dynode center
   26+: First dark ring around dynode center
  */
  if ( fPMT<19999 && fBolt > 0 && fBolt < 25) { return true; }   //finds if the data is really a bolt.
  return false;
}

std::istream& operator>>( std::istream& in, MedianTextRecord & rec ){
  unsigned num, pmt, bolt, x, y;
  std::string id;
  char dash; 
  in >> num >> pmt >> dash >> bolt >> x >> y >> id;
  rec = MedianTextRecord( num, pmt, bolt, x, y, id );
  return in;
}



MedianTextReader * MedianTextReader::instance = nullptr;

MedianTextReader::MedianTextReader() {}

MedianTextReader* MedianTextReader::Get(){
  if ( instance == nullptr ) instance = new MedianTextReader();
  return instance;
}
  
void MedianTextReader::set_input_file( const std::string& fname ){
  if ( fData.size() > 0 ){
    fData.clear();
  }

  std::ifstream infile( fname );
  MedianTextRecord r;
  while ( infile >> r ){
    if(r.is_bolt()){
      fData.push_back( r );
    }
  }
}
  

std::ostream&  operator<<( std::ostream& out, const MedianTextRecord& r ){

  out << r.photo_num() <<" "
      << r.pmt_num() <<"-"
      << r.bolt_num() <<" "
      << r.x() <<" "
      << r.y() <<" "
      << r.id() <<std::endl;

  return out;
}

