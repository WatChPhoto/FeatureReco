
#ifndef _MedianTextReader_hpp_
#define _MedianTextReader_hpp_


#include <string>
#include <vector>
#include <fstream>

/// Class to store Median Text Record
/// Photo-number  PMT-bolt   x    y  Person
///      045      00555-24 3491 2020 PD3
///
class MedianTextRecord{

public:
  MedianTextRecord( unsigned num, unsigned pmt, unsigned bolt, unsigned x, unsigned y, std::string id);
  MedianTextRecord();

  // Getters
  unsigned photo_num() const { return fNum; } // photo number		   
  unsigned pmt_num() const { return fPMT; }   // PMT number		   
  unsigned bolt_num() const { return fBolt; } // Bolt number		   
  unsigned x() const { return fX; }	      // X pixel location		   
  unsigned y() const { return fY; }	      // Y pixel location		   
  std::string id() const { return fID; }      // ID of person locating feature

  // Setters
  void set_photo_num( unsigned num ) { fNum = num; }  // photo number		   
  void set_pmt_num( unsigned pmt )   { fPMT = pmt; }  // PMT number		   
  void set_bolt_num( unsigned bolt ) { fBolt = bolt; }// Bolt number		   
  void set_x( unsigned x)            { fX = x; }	 // X pixel location		   
  void set_y( unsigned y)            { fY = y; }	 // Y pixel location		   
  void set_id( const std::string& id){ fID = id; }    // ID of person locating feature

private:
  unsigned fNum;   // photo number
  unsigned fPMT;   // PMT number
  unsigned fBolt;  // Bolt number
  unsigned fX;     // X pixel location
  unsigned fY;     // Y pixel location
  std::string fID; // ID of person locating feature

};


typedef std::vector< MedianTextRecord > MedianTextData;

std::istream& operator>>( std::istream& in, MedianTextRecord & rec );

/// Singleton class to read in data and store it
class MedianTextReader {
public:
  static MedianTextReader* Get(); // get instance of clas
  
  void set_input_file( const std::string& fname );
  
  const MedianTextData & get_data() const {return fData; }

private:
  MedianTextReader();
  static MedianTextReader* instance;
  MedianTextData fData;
  

};

std::ostream& operator<<( std::ostream& out, const MedianTextRecord& r );


#endif
