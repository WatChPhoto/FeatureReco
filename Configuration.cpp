#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>

#include "Configuration.hpp"

using namespace std;

//function to trim the leading and trailing whitespace.
void config::trim(string &name){
  unsigned low = name.find_first_not_of(' ');
  if(low>0 && low<name.length()){         //aviuds the garbage result.
  name = name.substr(low);
  }
  unsigned high = name.find(' ');
  if(high>0 && high < name.length()){name=name.substr(0,high);}
    
}

//function to change input string to lowercase.
void config::tolower(string &name){
  string temp="";
  for(unsigned i = 0; i< name.length(); i++){
    temp+=::tolower(name[i]) ;
  }
  name = temp;
}

//function to obtain the name of the variable from the line in the data file. 
string config::find_name(const string &line){
  string name;
  unsigned high = line.find_first_of('=');
  name = line.substr(0, high);
  trim(name);

  return name;
}

//template function to return template data type.
//use Get<return_type>("variable_name"); to get the value of variable from Config.txt file
//template<typename T> T config::Get(string key){
double config::Get_double(string key){

  ifstream Load("Config.txt");   //loading the data file
    if ( !Load.is_open() ){
	std::cout<<"Could not open Config.txt"<<std::endl;
	exit(0);
  }
    string line;
    if(key.length()>1){ trim(key);}
    tolower(key);

    while(true){

      if(Load.eof()) {
	cerr <<"Error variable "+key+" not found"<<endl;  
	Load.close(); 
	throw key;
      }
     
      getline(Load,line);            
      
      unsigned low = line.find_first_not_of(' ');
      if(low>0 && low<line.length()){
	line = line.substr(low);              //ignores the leading white space in the line.
      }
      
      if(line.length()==0||line[0]=='#'){    //ignores comment in the Config.txt file. comment start with # but can be changed from here. 
	continue;
      }
      
      string name = find_name(line);
      tolower(name);
      
      if(name==key){  //compares the key to first word of in the line.
	double value;
	int index = line.find('=')+1;  
	line = line.substr(index);
	stringstream ss(line);    //finding the value of the key.
	ss >> value;              //converting the obtained value to T(template) type
	Load.close();            
	return value;
      }
  }
  }





int config::Get_int(string key){

  ifstream Load("Config.txt");   //loading the data file
  if ( !Load.is_open() ){
	std::cout<<"Could not open Config.txt"<<std::endl;
	exit(0);
  }

    string line;
    if(key.length()>1){ trim(key);}
    tolower(key);

    while(true){

      if(Load.eof()) {cerr <<"Error variable "+key+" not found"<<endl;  Load.close(); 
	throw key;
      }
     
      getline(Load,line);            
      
      unsigned low = line.find_first_not_of(' ');
      if(low>0 && low<line.length()){
	line = line.substr(low);              //ignores the leading white space in the line.
      }
      
      if(line.length()==0||line[0]=='#'){    //ignores comment in the Config.txt file. comment start with # but can be changed from here. 
	continue;
      }
      
      string name = find_name(line);
      tolower(name);
      
      if(name==key){  //compares the key to first word of in the line.
	int value;
	int index = line.find('=')+1;  
	line = line.substr(index);
	stringstream ss(line);    //finding the value of the key.
	ss >> value;              //converting the obtained value to T(template) type
	Load.close();            
	return value;
      }
  }
  }


//For testing purpose
/* 
int main(){
  cout << Get<int>("d")<<endl;
  cout <<Get<int>("   sigSPace   ")<<endl;
  cout<< 2.15*Get<double>("blursigma")<<endl;
  //  cout<< Get<int>("dsklf") <<endl;
  cout<<Get<int>(" ")<<endl;

}
*/

