#include "FeatureTTree.hpp"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <TH2D.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TMarker.h>
#include <TEllipse.h>
#include <TTree.h>
#include <TSystem.h>
#include <TF1.h>
#include <TGraphErrors.h>

#include "ImageDataReader.hpp"

using std::ostringstream;
using std::string;


const double pi = std::acos(-1);
double RADTODEG( double rad ){ return rad*180/pi; }
double DEGTORAD( double deg ){ return deg*pi/180; }

//class ImageData;
const double Rsk = 1690.0; // radial location of PMTs (cm)
const double th_f = 93.9; // field of view (degrees)

// d_of_theta
// return the in-plane distance to the wall given the camera position (r_K,theta_K) and facing (theta_DC)
// x[0] is angle in image
double d_of_theta( double *x, double *p ){
  // parameters:
  // p[0] = r_K       radial position of camera (cm)
  // p[1] = theta_K   angular position of camera (degrees)
  // p[2] = theta_DC  camera facing direction (degrees)

  double r_K    = p[0]; 
  double th_K   = p[1];
  double th_DC  = p[2];
  
  double thx    = x[0];
  double th_okx = 180.0 - thx + th_K;
  if ( th_okx > 180.0 ) th_okx = 360.0 - th_okx;
  double phx    = RADTODEG( std::asin( r_K / Rsk * std::sin( DEGTORAD( th_okx ) ) ) );
  double th_kox = 180.0 - phx - th_okx;
  double dx     = std::sqrt(  Rsk*Rsk + r_K*r_K - 2*r_K*Rsk*std::cos( DEGTORAD( th_kox ) ) );
  return dx;
}


// d_of_xpixel
// return the in-plane distance to the wall given the camera position (r_K,theta_K) and facing (theta_DC)
// x[0] is xpixel in image
double d_of_xpixel( double *x, double *p ){
  // parameters:
  // p[0] = r_K       radial position of camera (cm)
  // p[1] = theta_K   angular position of camera (degrees)
  // p[2] = theta_DC  camera facing direction (degrees)

  double r_K    = p[0]; 
  double th_K   = p[1];
  double th_DC  = p[2];
  double th_pix0 = th_DC + th_f/2;
  
  double thx    = th_pix0 - x[0] * th_f / 4000.0;
  double th_okx = 180.0 - thx + th_K;
  if ( th_okx > 180.0 ) th_okx = 360.0 - th_okx;
  double phx    = RADTODEG( std::asin( r_K / Rsk * std::sin( DEGTORAD( th_okx ) ) ) );
  double th_kox = 180.0 - phx - th_okx;
  double dx     = std::sqrt(  Rsk*Rsk + r_K*r_K - 2*r_K*Rsk*std::cos( DEGTORAD( th_kox ) ) );
  return dx;
}



void basic_histograms( TFile* fout, const ImageData& idt, const int imgn, bool include2d=false  ){
  
  ostringstream os;
  os<<"image_"<<std::setw(3)<<std::setfill('0')<<imgn;
  TDirectory* td=fout->mkdir( os.str().c_str() );
  td->cd();

  ostringstream os2;
  os2<<std::setw(3)<<std::setfill('0')<<imgn;
  std::string imgnum = os2.str();

  std::vector< TH1* > allhists;

  using std::string;
  TH1D* ha_intensity   = new TH1D( (string("ha_intensity")+imgnum).c_str()  ,"All blob intensity  ; intensity; count/bin",256,0.,256.);
  TH1D* ha_area        = new TH1D( (string("ha_area")+imgnum).c_str()       ,"All blob area       ; area; count/bin",100,0.,100.);
  TH1D* ha_circularity = new TH1D((string("ha_circularity")+imgnum).c_str(),"All blob circularity; circularity; count/bin",100,0.,1.);
  TH1D* ha_convexity   = new TH1D((string("ha_convexity")+imgnum).c_str()  ,"All blob convexity  ; convexity; count/bin",100,0.,1.);
  TH1D* ha_inertia     = new TH1D((string("ha_inertia")+imgnum).c_str()    ,"All blob inertia    ; inertia; count/bin",100,0.,1.);

  allhists.push_back( ha_intensity );
  allhists.push_back( ha_area );
  allhists.push_back( ha_circularity );
  allhists.push_back( ha_convexity );
  allhists.push_back( ha_inertia );
    
  
  TH1D* he_intensity   = new TH1D((string("he_intensity")+imgnum).c_str()  ,"Ellipse blob intensity  ; intensity; count/bin",256,0.,256.);
  TH1D* he_area        = new TH1D((string("he_area"    )+imgnum).c_str()   ,"Ellipse blob area       ; area; count/bin",100,0.,100.);
  TH1D* he_circularity = new TH1D((string("he_circularity")+imgnum).c_str(),"Ellipse blob circularity; circularity; count/bin",100,0.,1.);
  TH1D* he_convexity   = new TH1D((string("he_convexity" )+imgnum).c_str() ,"Ellipse blob convexity  ; convexity; count/bin",100,0.,1.);
  TH1D* he_inertia     = new TH1D((string("he_inertia"   )+imgnum).c_str() ,"Ellipse blob inertia    ; inertia; count/bin",100,0.,1.);

  allhists.push_back( he_intensity );
  allhists.push_back( he_area );
  allhists.push_back( he_circularity );
  allhists.push_back( he_convexity );
  allhists.push_back( he_inertia );
  
  TH2D* hblob_xy=new TH2D((string("hblob_xy")+imgnum).c_str(),"Ellipse bolt locations and found ellipses; x (pixels); y (pixels)",4000,0,4000,3000,0,3000);

  allhists.push_back( hblob_xy );


  TH1D* el_peakval     = new TH1D((string("el_peakval")+imgnum).c_str()    ,"Ellipse hough peak value   ;  peak value; count/bin",100,0.,100.);
  TH1D* el_chi2        = new TH1D((string("el_chi2"   )+imgnum).c_str()    ,"Ellipse #chi^{2}           ;    #chi${2}; count/bin",100,0.,1000.);
  TH1D* el_rchi2       = new TH1D((string("el_rchi2"   )+imgnum).c_str()   ,"Ellipse reduced #chi^{2}   ; reduced #chi${2}; count/bin",100,0.,50.);
  TH1D* el_ndof        = new TH1D((string("el_ndof"   )+imgnum).c_str()    ,"Ellipse number of bolts    ; peak value; count/bin",100,0.,100.);




    
  allhists.push_back( el_peakval );
  allhists.push_back( el_chi2 );
  allhists.push_back( el_rchi2 );

  TH1D* hp_intensity   = new TH1D((string("hp_intensity" )+imgnum).c_str() ,"PMT bolt intensity  ; intensity; count/bin",256,0.,256.);
  TH1D* hp_area        = new TH1D((string("hp_area"     )+imgnum).c_str()  ,"PMT bolt area       ; area; count/bin",100,0.,100.);
  TH1D* hp_circularity = new TH1D((string("hp_circularity")+imgnum).c_str(),"PMT bolt circularity; circularity; count/bin",100,0.,1.);
  TH1D* hp_convexity   = new TH1D((string("hp_convexity")+imgnum).c_str()  ,"PMT bolt convexity  ; convexity; count/bin",100,0.,1.);
  TH1D* hp_inertia     = new TH1D((string("hp_inertia"  )+imgnum).c_str()  ,"PMT bolt inertia    ; inertia; count/bin",100,0.,1.);

  allhists.push_back( hp_intensity );
  allhists.push_back( hp_area );
  allhists.push_back( hp_circularity );
  allhists.push_back( hp_convexity );
  allhists.push_back( hp_inertia );

  
  TH2D* hpblob_xy=new TH2D((string("hpblob_xy")+imgnum).c_str(),"PMT and bolt locations; x (pixels); y (pixels)",4000,0,4000,3000,0,3000);

  allhists.push_back( hpblob_xy );
  
  TH1D* p_peakval     = new TH1D((string("p_peakval" )+imgnum).c_str()   ,"PMT hough peak value   ;       peak value; count/bin",100, 0., 100. );
  TH1D* p_chi2        = new TH1D((string("p_chi2"   )+imgnum).c_str()    ,"PMT #chi^{2}           ;         #chi${2}; count/bin",100, 0., 1000.);
  TH1D* p_rchi2       = new TH1D((string("p_rchi2" )+imgnum).c_str()     ,"PMT reduced #chi^{2}   ; reduced #chi${2}; count/bin",100, 0., 50.  );
  TH1D* p_ndof        = new TH1D((string("p_ndof"   )+imgnum).c_str()    ,"PMT number of bolts    ;       peak value; count/bin",100, 0., 100. );

    
  allhists.push_back( p_peakval );
  allhists.push_back( p_chi2 );
  allhists.push_back( p_rchi2 );

  
  TH2D* hablob_intensity  =nullptr;
  TH2D* hablob_area         =nullptr;
  TH2D* hablob_circularity  =nullptr;
  TH2D* hablob_convexity    =nullptr;
  TH2D* hablob_inertia      =nullptr;

  TH2D* el2_peakval    = nullptr;
  TH2D* el2_rchi2      = nullptr;
  TH2D* el2_ndof       = nullptr;


  TH2D* p2_peakval    = nullptr;
  TH2D* p2_rchi2      = nullptr;
  TH2D* p2_ndof       = nullptr;
  
  if (include2d){
    hablob_intensity  =new TH2D((string("hablob_intensity")+imgnum).c_str()  ,"All blob intensity; x (pixels); y (pixels)",4000,0,4000,3000,0,3000);
    hablob_area       =new TH2D((string("hablob_area" )+imgnum).c_str()      ,"All blob area; x (pixels); y (pixels)",4000,0,4000,3000,0,3000);
    hablob_circularity=new TH2D((string("hablob_circularity")+imgnum).c_str(),"All blob circularity; x (pixels); y (pixels)",4000,0,4000,3000,0,3000);
    hablob_convexity  =new TH2D((string("hablob_convexity"  )+imgnum).c_str(),"All blob convexity; x (pixels); y (pixels)",4000,0,4000,3000,0,3000);
    hablob_inertia    =new TH2D((string("hablob_inertia"  )+imgnum).c_str()  ,"All blob inertia; x (pixels); y (pixels)",4000,0,4000,3000,0,3000);

    allhists.push_back( hablob_intensity );
    allhists.push_back( hablob_area );
    allhists.push_back( hablob_circularity );
    allhists.push_back( hablob_convexity );
    allhists.push_back( hablob_inertia );

    el2_peakval    = new TH2D((string("el2_peakval" )+imgnum).c_str()  ,"Ellipse hough peak value   ; x (pixels); y (pixels)",4000,0,4000,3000,0,3000);
    el2_rchi2      = new TH2D((string("el2_rchi2"  )+imgnum).c_str()   ,"Ellipse reduced #chi^{2}   ; x (pixels); y (pixels)",4000,0,4000,3000,0,3000);
    el2_ndof       = new TH2D((string("el2_ndof"  )+imgnum).c_str()    ,"Ellipse number of bolts    ; x (pixels); y (pixels)",4000,0,4000,3000,0,3000);

    allhists.push_back( el2_peakval );
    allhists.push_back( el2_rchi2 );
    allhists.push_back( el2_ndof );

    p2_peakval    = new TH2D((string("p2_peakval" )+imgnum).c_str()  ,"PMT hough peak value   ; x (pixels); y (pixels)",4000,0,4000,3000,0,3000);
    p2_rchi2      = new TH2D((string("p2_rchi2"  )+imgnum).c_str()   ,"PMT reduced #chi^{2}   ; x (pixels); y (pixels)",4000,0,4000,3000,0,3000);
    p2_ndof       = new TH2D((string("p2_ndof"  )+imgnum).c_str()    ,"PMT number of bolts    ; x (pixels); y (pixels)",4000,0,4000,3000,0,3000);
    allhists.push_back( p2_peakval );
    allhists.push_back( p2_rchi2 );
    allhists.push_back( p2_ndof );
  }
  


  
  for ( const BlobData& b : idt.fBlobs ){
    if (include2d){
      hablob_intensity  ->Fill( b.x, 3000-b.y, b.intensity );
      hablob_area       ->Fill( b.x, 3000-b.y, b.area );
      hablob_circularity->Fill( b.x, 3000-b.y, b.circularity );
      hablob_convexity  ->Fill( b.x, 3000-b.y, b.convexity );
      hablob_inertia    ->Fill( b.x, 3000-b.y, b.inertia );
    }
      
    ha_intensity  ->Fill( b.intensity );
    ha_area       ->Fill( b.area );
    ha_circularity->Fill( b.circularity );
    ha_convexity  ->Fill( b.convexity );
    ha_inertia    ->Fill( b.inertia );
  }

  for ( const EllipseData& el : idt.fEllipses ){
    el_peakval->Fill( el.peakval );
    el_chi2->Fill( el.chi2 );
    el_rchi2->Fill( el.chi2/el.ndof );
    el_ndof->Fill( el.ndof );


    if (include2d){
      el2_peakval->Fill( el.xx, 3000-el.yy, el.peakval );
      el2_rchi2->Fill( el.xx, 3000-el.yy, el.chi2/el.ndof );
      el2_ndof->Fill( el.xx, 3000-el.yy, el.ndof );
    }

    for ( const BlobData& b : el.blobentry ){
      hblob_xy->Fill( b.x, 3000-b.y );
      hblob_xy->SetMarkerStyle(7);
      hblob_xy->SetMarkerColor(kRed);

      he_intensity  ->Fill( b.intensity );
      he_area       ->Fill( b.area );
      he_circularity->Fill( b.circularity );
      he_convexity  ->Fill( b.convexity );
      he_inertia    ->Fill( b.inertia );


    }
    TEllipse* elli = new TEllipse( el.xx, 3000-el.yy, el.bb, el.aa, 0., 360., RADTODEG( el.phi ) );
    elli->SetFillStyle(0);
    elli->SetLineColor(kRed);
    elli->SetLineWidth(2);
    hblob_xy->GetListOfFunctions()->Add( elli );
    TMarker *mark = new TMarker( el.xx, 3000-el.yy, 20 );
    mark->SetMarkerColor(kRed);
    hblob_xy->GetListOfFunctions()->Add( mark );

  }


  for ( const EllipseData& el : idt.fPMTs ){
    p_peakval->Fill( el.peakval );
    p_chi2->Fill( el.chi2 );
    p_rchi2->Fill( el.chi2/el.ndof );
    p_ndof->Fill( el.ndof );

    if (include2d){
      p2_peakval->Fill( el.xx, 3000-el.yy, el.peakval );
      p2_rchi2->Fill( el.xx, 3000-el.yy, el.chi2/el.ndof );
      p2_ndof->Fill( el.xx, 3000-el.yy, el.ndof );
    }

    for ( const BlobData& b : el.blobentry ){
      hpblob_xy->Fill( b.x, 3000-b.y );
      hpblob_xy->SetMarkerStyle(7);
      hpblob_xy->SetMarkerColor(kBlue);

      hp_intensity  ->Fill( b.intensity );
      hp_area       ->Fill( b.area );
      hp_circularity->Fill( b.circularity );
      hp_convexity  ->Fill( b.convexity );
      hp_inertia    ->Fill( b.inertia );


    }
    TEllipse* elli = new TEllipse( el.xx, 3000-el.yy, el.bb, el.aa, 0., 360., RADTODEG( el.phi ) );
    elli->SetFillStyle(0);
    elli->SetLineColor(kBlue);
    elli->SetLineWidth(2);
    hpblob_xy->GetListOfFunctions()->Add( elli );
    TMarker *mark = new TMarker( el.xx, 3000-el.yy, 20 );
    mark->SetMarkerColor(kBlue);
    hpblob_xy->GetListOfFunctions()->Add( mark );

  }

  //for (TH1* hist: allhists){
  //  hist->SetDirectory( td );
  //  hist->Write();
  //}
  //td->Write();
  //td->Close();
}


double piecewiselineh( double * x, double * p ){
  //  y =  ( yc - y0 ) / xc * x + y0   for x < xc
  //       ( y4000 - yc ) /( 4000-xc ) * x + y4000  for x>=xc 
  // parameters:
  // p[0] = xc
  // p[1] = yc
  // p[2] = y0
  // p[2] = y4000

  double xc    = p[0];
  double yc    = p[1];
  double y0    = p[2];
  double y4000 = p[3];

  // force bound on xc:
  if ( xc<=1.0e-6 ) xc=1e-6;
  if ( xc>=3999.999999 ) xc=3999.999999;
  
  double y=0.;
  if ( x[0] < xc ){
    y = (yc-y0)/xc * x[0] + y0;
  } else {
    y = (y4000-yc)/(4000.0-xc) * (x[0]-4000.) + y4000;
  }
  return y;
}



double piecewiselinev( double * x, double * p ){
  //  y =  ( yc - y0 ) / xc * x + y0   for x < xc
  //       ( y3000 - yc ) /( 3000-xc ) * x + y3000  for x>=xc 
  // parameters:
  // p[0] = xc
  // p[1] = yc
  // p[2] = y0
  // p[2] = y3000

  double xc    = p[0];
  double yc    = p[1];
  double y0    = p[2];
  double y3000 = p[3];


  // force bound on xc:
  if ( xc<=1.0e-6 ) xc=1e-6;
  if ( xc>=2999.999999 ) xc=2999.999999;
  
  double y=0.;
  if ( x[0] < xc ){
    y = (yc-y0)/xc * x[0] + y0;
  } else {
    y = (y3000-yc)/(3000.0-xc) * (x[0]-3000.0) + y3000;
  }
  return y;
}

double simplelinev( double * x, double * p ){
  //  y =  ( y3000 - y0 ) / 3000.0 * x + y0   
  //
  // parameters:
  // p[0] = y0
  // p[1] = y3000

  double y0    = p[0];
  double y3000 = p[1];


  return (y3000-y0)/3000.0 * x[0] + y0;
}

double simplelineh( double * x, double * p ){
  //  y =  ( y4000 - y0 ) / 4000.0 * x + y0   
  //
  // parameters:
  // p[0] = y0
  // p[1] = y4000

  double y0    = p[0];
  double y4000 = p[1];


  return (y4000-y0)/4000.0 * x[0] + y0;
}




void ellipse_size_vs_pos( TFile * fout, const ImageData& idt ){

  // get theta_DC for this image
  ImageDataReader & idr = ImageDataReader::GetInstance();
  CameraFace cf = idr.GetFace( "BarrelSurveyFar", idt.ips.imgnum );
  float theta_DC = cf.yaw - 19.0; // 19 degrees offset between yaw and theta_DC
  std::cout<<"Image "<<idt.ips.imgnum<<" yaw is "<<theta_DC<<std::endl;
  
  ostringstream os;
  os<<"esize_"<<std::setw(3)<<std::setfill('0')<<idt.ips.imgnum;
  TDirectory* td=fout->mkdir( os.str().c_str() );
  td->cd();

  ostringstream os2;
  os2<<std::setw(3)<<std::setfill('0')<<idt.ips.imgnum;
  string imgnum = os2.str();


  // make histogram of all ellipse size parameter
  TH1D* hball = new TH1D( (string("hball_")+imgnum).c_str(), " ; ellipse b (pixels); counts/bin", 150, 50.0, 200.0 );
  for ( const EllipseData& el : idt.fPMTs ){
    hball->Fill( el.bb );
  }

  // make histograms of ellipse locations
  float bmean = hball->GetMean();
  float brms  = hball->GetRMS();
  TH1D* hlocx = new TH1D( (string("hlocx_")+imgnum).c_str(), " ; ellipse b (pixels); counts/bin", 200, 0., 4000. );
  TH1D* hlocy = new TH1D( (string("hlocy_")+imgnum).c_str(), " ; ellipse b (pixels); counts/bin", 150, 0., 3000. );
  for ( const EllipseData& el : idt.fPMTs ){
    float x = el.xx;
    float y = 3000-el.yy;
    float b = el.bb;
    // trim any more than 2-sigma from mean
    if ( fabs( b - bmean ) > 2*brms ) continue;
    hlocx->Fill( x );
    hlocy->Fill( y );
  }

  // now make graphs of "distance to wall" versus angle
  float xmin = hlocx->GetMean() - 500.0;
  float xmax = hlocx->GetMean() + 500.0;
  float ymin = hlocy->GetMean() - 500.0;
  float ymax = hlocy->GetMean() + 500.0;

    
  // vectors for a horizontal and vertical band!
  std::vector< float > h_x, hx_e;
  std::vector< float > h_b, hb_e;
  std::vector< float > v_y, vy_e;
  std::vector< float > v_b, vb_e;

  
  // d = -23.9 b + 3700.2
  //(corrected) d = -22.34b + 3513.92
  const double cm_per_pixel = -22.34;//-23.9; // cm of distance per ellipse-b pixel width
  const double d_offset = 3513.92;//3700.2; // cm
  for ( const EllipseData& el : idt.fPMTs ){
    float x = el.xx;
    float y = 3000-el.yy;
    float b = el.bb;

    // trim any more than 2-sigma from mean
    if ( fabs( b - bmean ) > 2*brms ) continue;

    
    // horizontal band (vertical extent from 1200 to 1800?)
    if ( y > ymin && y < ymax ){
      h_x.push_back( x ); hx_e.push_back( 0.0 );
      h_b.push_back( b * cm_per_pixel + d_offset ); hb_e.push_back( -0.5*cm_per_pixel );
    }
    // vertical band (horizontal extent from 1700 to 2300?)
    if ( x > xmin && x < xmax ){
      v_y.push_back( y ); vy_e.push_back( 0.0 );
      v_b.push_back( b * cm_per_pixel + d_offset ); vb_e.push_back( -0.5*cm_per_pixel );
    }
  }

  TGraphErrors * tgh = new TGraphErrors( h_x.size(), &h_x[0], &h_b[0], &hx_e[0], &hb_e[0] );
  tgh->SetName( (string("tgh_") + imgnum).c_str() );
  tgh->SetTitle(" ; x (pixels); distance from ellipse b (cm)");
  tgh->SetMarkerStyle(20);
  tgh->SetMarkerColor(kBlue);
 
  TF1* tfh = new TF1( (string("tfh_") + imgnum).c_str(), d_of_xpixel, 0., 4000., 3 );
  tfh->SetParNames( "r_{K} (cm)       ",
		    "#theta_{K} (deg) ",
		    "#theta_{DC} (deg)" );
  tfh->SetParameter( 0, Rsk - tgh->Eval( 2000.0 ) ); // first guess for r_K
  tfh->SetParameter( 1, theta_DC ); // first guess for theta_K is camera facing
  tfh->SetParameter( 2, theta_DC );
  tfh->FixParameter( 2, theta_DC); // camera facing is fixed
  tgh->Fit( tfh , "Q" );

  if ( tfh->GetParameter(0) < 1000.0 ||
       tfh->GetParameter(0) > 3000.0 ){ 
    // refit with a straight line
    TF1* tfh2 = new TF1( (string("tfh2_") + imgnum).c_str(), simplelineh, 0., 4000., 2 );
    int n = tgh->GetN();
    tfh2->SetLineColor( kGreen );
    tfh2->SetParNames( "y0", "y4000" );
    tfh2->SetParameters( tgh->GetPointY( 0 ), tgh->GetPointY( n-1 ) );
    tgh->Fit( tfh2, "Q" );
  }
  
  tgh->Write();
  
  

  TGraphErrors * tgv = new TGraphErrors( v_y.size(), &v_y[0], &v_b[0], &vy_e[0], &vb_e[0] );
  tgv->SetName( (string("tgv_") + imgnum).c_str() );
  tgv->SetTitle(" ; y (pixels); ellipse b (pixels)");
  tgv->SetMarkerStyle(20);
  tgv->SetMarkerColor(kBlue);

  TF1* tfv = new TF1( (string("tfv_") + imgnum).c_str(), piecewiselinev, 0., 3000., 4 );
  tfv->SetParNames( "xc", "yc", "y0", "y3000" );
  tfv->SetParameters( 1500., tgv->GetMaximum(), tgv->GetMaximum()-10.0, tgv->GetMaximum()-10.0 );
  tgv->Fit( tfv, "Q" );


  if ( tfv->GetParameter(0) < 1000.0 ||
       tfv->GetParameter(0) > 2000.0 ){
    // refit with a straight line
    TF1* tfv2 = new TF1( (string("tfv2_") + imgnum).c_str(), simplelinev, 0., 3000., 2 );
    int n = tgv->GetN();
    tfv2->SetLineColor( kGreen );
    tfv2->SetParNames( "y0", "y3000" );
    tfv2->SetParameters( tgv->GetPointY( 0 ), tgv->GetPointY( n-1 ) );
    tgv->Fit( tfv2, "Q" );
  }


  
  tgv->Write();



}


void analyse_merged_ttree( string filename = "merge_ttrees.root" ){

  gSystem->Load("libFeatureTreeLib.so");

  ImageData * imgdata = new ImageData();
  TFile     * fin     = new TFile( filename.c_str() ,"read");
  TTree     * tree    = (TTree*) fin->Get("ImageTree");
  tree->SetBranchAddress( "ImageData", &imgdata );
  tree->GetEntry(0);

  TFile* fout = new TFile("analyse_ttree.root","recreate");

  for ( unsigned i=0; i<tree->GetEntries(); ++i){
    std::cout<<"Analyse image "<<imgdata->ips.imgnum<<std::endl;
    tree->GetEntry( i );
    basic_histograms( fout, *imgdata, imgdata->ips.imgnum );
    ellipse_size_vs_pos( fout, *imgdata );
  }
  fout->Write();

}

int main(){

  analyse_merged_ttree();
  return 0;
  
}
