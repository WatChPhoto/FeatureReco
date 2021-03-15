#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <numeric>
#include <algorithm>


#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TDirectory.h>
#include <TEllipse.h>
#include <TMarker.h>

struct ellipse{
  int num;
  float x;
  float y;
  float b;
  float e;
  float phi;
  double get_a() const {return b/std::sqrt(1-e*e) ;}
};

std::istream& operator>>( std::istream& is, ellipse& e ){
  is >> e.num;
  is >> e.x;
  is >> e.y;
  is >> e.b;
  is >> e.e;
  is >> e.phi;
  return is;
}


std::ostream& operator<<( std::ostream& os, const ellipse& e ){
  os << "num=" << e.num;
  os << "  x=" << e.x;
  os << "  y=" << e.y;
  os << "  b=" << e.b;
  os << "  e=" << e.e;
  os << " phi="<< e.phi << std::endl;
  return os;
}

struct image_pmtinfo {
  int num;
  std::vector< ellipse > pmts; 
};

const double pi = std::acos(-1);
float RADTODEG( float rad ){ return rad*180/pi; }


TH2S* plot_image_pmtinfo( const image_pmtinfo& info ){
  std::ostringstream os;
  os << "image"<<info.num;
  TH2S * hplot = new TH2S(os.str().c_str(), " ; x pixels; y pixels ", 4000, 0., 4000., 3000, 0., 3000. );

  for ( const ellipse& el : info.pmts ){
    TEllipse *elli = new TEllipse( el.x, el.y, el.b, 
				 el.get_a(), 0., 360., RADTODEG( el.phi ) );
    elli->SetLineWidth(2);
    elli->SetLineColor(kBlue);
    elli->SetFillStyle(0);
    hplot->GetListOfFunctions()->Add( elli );
    
    TMarker *mark = new TMarker( el.x, el.y, 20 );
    mark->SetMarkerColor(kRed);
    hplot->GetListOfFunctions()->Add( mark );
  }
  return hplot;
}

// fill vector of pmts in image_pmtinfo
// have to separately fill image num
std::istream& operator>>( std::istream& is, image_pmtinfo& info ){
  ellipse e;
  while ( is >> e ){
    info.pmts.push_back( e );
  }
  return is;
}


std::ostream& operator<<( std::ostream& os, const image_pmtinfo& info ){
  os << "Image number " << info.num <<std::endl;
  for ( const ellipse & e : info.pmts ){
    os << e;
  }
  return os;
}



float std_average( const std::vector< float > & v ){
  if (v.size()==0) return 0;
  return std::accumulate( v.begin(), v.end(), 0 )/v.size();
}

float std_rms( const std::vector< float > & v ){
  float avg = std_average( v );
  double sq_sum = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
  double stdev = std::sqrt(sq_sum / v.size() - avg * avg);
  return stdev;
}

int main(){

  std::vector< std::string > filenames;
  std::ifstream in( "files_to_summarize.txt" );
  std::string fname;
  while ( in >> fname ){
    filenames.push_back( fname );
  }

  std::vector< image_pmtinfo > images;
  for ( const std::string& fname : filenames ){
    std::cout<<"file = "<< fname<<std::endl;
    std::ifstream infile( fname );
    image_pmtinfo info;
    std::stringstream is( fname.substr(0,3) );
    int fnum;
    is >> fnum;
    info.num = fnum;
    infile >> info;
    images.push_back( info );
  }

  for ( const image_pmtinfo& info : images ){
    std::cout << info;
  }


  

  TFile * fout = new TFile("summarize_pmts_found.root","recreate");

  TDirectory* curdir = gFile->CurrentDirectory();
  TDirectory* td = fout->mkdir( "images" );
  td->cd();
  for ( const image_pmtinfo& info: images ){
    TH2S * hplot = plot_image_pmtinfo( info );
  }
  curdir->cd();
  
  TH1D *  hnpmt = new TH1D("hnpmt"," ; Number of PMTs found; counts/bin", 100, -0.5, 99.5 );
  TH1D *  hb = new TH1D("hb", " ; ellipse b (pixels); counts/bin", 200, 50.,250. );
  TH1D *  he = new TH1D("he", " ; ellipse eccentricity ; counts/bin", 100, 0., 1.);

  
  std::vector< float > imagenums;
  std::vector< float > numpmts;
  std::vector< float > avgb;
  std::vector< float > avge;
  std::vector< float > rmsb;
  std::vector< float > rmse;
  
  for ( const image_pmtinfo& info : images ){
    hnpmt->Fill( info.pmts.size() );
    if ( info.pmts.size() < 25 ) {
      std::cout<<"LESS THAN 25 PMTS for image "<<info.num<<std::endl;
    }
    imagenums.push_back( info.num );
    numpmts.push_back( info.pmts.size() );
    std::vector< float > e_vals;
    std::vector< float > b_vals;
    for ( const ellipse& e : info.pmts ){
      hb->Fill( e.b );
      he->Fill( e.e );

      e_vals.push_back( e.e );
      b_vals.push_back( e.b );
    }
    avge.push_back( std_average( e_vals ) );
    rmse.push_back( std_rms( e_vals ) );
    avgb.push_back( std_average( b_vals ) );
    rmsb.push_back( std_rms( b_vals ) );    
  }

  TGraph * gnpmt = new TGraph( imagenums.size(), &imagenums[0], &numpmts[0] );
  gnpmt->SetName("gnpmt");
  gnpmt->SetTitle(" ; image number; number of PMTs" );
  gnpmt->Write();

  TGraph * gavgb = new TGraph( imagenums.size(), &imagenums[0], &avgb[0] );
  gavgb->SetName("gavgb");
  gavgb->SetTitle(" ; image number; average ellipse b" );
  gavgb->Write();

  TGraph * gavge = new TGraph( imagenums.size(), &imagenums[0], &avge[0] );
  gavge->SetName("gavge");
  gavge->SetTitle(" ; image number; average ellipse e" );
  gavge->Write();

  TGraph * grmsb = new TGraph( imagenums.size(), &imagenums[0], &rmsb[0] );
  grmsb->SetName("grmsb");
  grmsb->SetTitle(" ; image number; RMS ellipse b" );
  grmsb->Write();

  TGraph * grmse = new TGraph( imagenums.size(), &imagenums[0], &rmse[0] );
  grmse->SetName("grmse");
  grmse->SetTitle(" ; image number; RMS ellipse e" );
  grmse->Write();



  fout->Write();
  fout->Close();
	      


  return 0;
}



