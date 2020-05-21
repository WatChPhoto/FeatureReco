#ifndef _Configuration_hpp_
#define _Configuration_hpp_

#include <string>

namespace config {

  //function to trim the leading and trailing whitespace.
  void trim( std::string &name);
  
  //function to change input string to lowercase.
  void tolower( std::string &name);

  //function to obtain the name of the variable from the line in the data file. 
  std::string find_name(const std::string &line);

  //template function to return template data type.
  //use Get<return_type>("variable_name"); to get the value of variable from Config.txt file
  //template<typename T> T Get(std::string key);

  int   Get_int(std::string);
  double Get_double(std::string);

}



#endif
