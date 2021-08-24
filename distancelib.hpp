#include"camera.hpp"

std::istream& operator>>( std::istream& is, WorldPoints& t){
  is>>t.id;
  is>>t.x;
  is>>t.y;
  is>>t.z;

  return is;
}

std::istream& operator>>( std::istream& is, Ellipse& p){
  is>>p.id;
  is>>p.x;
  is>>p.y;
  is>>p.b;
  is>>p.e;
  is>>p.phi;

  return is;
}

//read all pmts in sk detector.
std::vector<WorldPoints> read_all_world_points(std::string filename="SK_all_PMT_locations.txt"){
  std::vector<WorldPoints> all_pmts; 
  WorldPoints t; //temporary container

  std::ifstream myfile(filename);
  while (myfile >> t) {
    all_pmts.push_back(t);
  }
  myfile.close();
  return all_pmts;  
}

//read pmts from an image.
std::vector<Ellipse> read_ellipses_in_image(std::string filename ="045.txt"){
  std::vector<Ellipse>ellipses;
  Ellipse e;
  std::ifstream myfile(filename);
  while(myfile>>e){
    ellipses.push_back(e);
  }
  return ellipses;
}

unsigned get_img_num(std::string filename){
  unsigned img_num;
  std::stringstream ss;  
  ss << filename;  
  ss >> img_num;

  return img_num;
}
