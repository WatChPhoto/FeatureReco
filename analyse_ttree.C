//#include "FeatureTTree.hpp"
const double pi = std::acos(-1);
float RADTODEG( float rad ){ return rad*180/pi; }

//class ImageData;

void basic_histograms( const ImageData& idt ){
  TFile* fout = new TFile("analyse_ttree.root","recreate");

  TH1D* ha_intensity   = new TH1D("ha_intensity"  ,"All blob intensity  ; intensity; count/bin",256,0.,256.);
  TH1D* ha_area        = new TH1D("ha_area"       ,"All blob area       ; area; count/bin",100,0.,100.);
  TH1D* ha_circularity = new TH1D("ha_circularity","All blob circularity; circularity; count/bin",100,0.,1.);
  TH1D* ha_convexity   = new TH1D("ha_convexity"  ,"All blob convexity  ; convexity; count/bin",100,0.,1.);
  TH1D* ha_inertia     = new TH1D("ha_inertia"    ,"All blob inertia    ; inertia; count/bin",100,0.,1.);

  TH2D* hablob_intensity  =new TH2D("hablob_intensity"  ,"All blob intensity; x (pixels); y (pixels)",4000,0,4000,3000,0,3000);
  TH2D* hablob_area       =new TH2D("hablob_area"       ,"All blob area; x (pixels); y (pixels)",4000,0,4000,3000,0,3000);
  TH2D* hablob_circularity=new TH2D("hablob_circularity","All blob circularity; x (pixels); y (pixels)",4000,0,4000,3000,0,3000);
  TH2D* hablob_convexity  =new TH2D("hablob_convexity"  ,"All blob convexity; x (pixels); y (pixels)",4000,0,4000,3000,0,3000);
  TH2D* hablob_inertia    =new TH2D("hablob_inertia"    ,"All blob inertia; x (pixels); y (pixels)",4000,0,4000,3000,0,3000);

  TH1D* he_intensity   = new TH1D("he_intensity"  ,"Ellipse blob intensity  ; intensity; count/bin",256,0.,256.);
  TH1D* he_area        = new TH1D("he_area"       ,"Ellipse blob area       ; area; count/bin",100,0.,100.);
  TH1D* he_circularity = new TH1D("he_circularity","Ellipse blob circularity; circularity; count/bin",100,0.,1.);
  TH1D* he_convexity   = new TH1D("he_convexity"  ,"Ellipse blob convexity  ; convexity; count/bin",100,0.,1.);
  TH1D* he_inertia     = new TH1D("he_inertia"    ,"Ellipse blob inertia    ; inertia; count/bin",100,0.,1.);


  TH2D* hblob_xy=new TH2D("hblob_xy","Ellipse bolt locations and found ellipses; x (pixels); y (pixels)",4000,0,4000,3000,0,3000);

  for ( const BlobData& b : idt.fBlobs ){
    hablob_intensity  ->Fill( b.x, 3000-b.y, b.intensity );
    hablob_area       ->Fill( b.x, 3000-b.y, b.area );
    hablob_circularity->Fill( b.x, 3000-b.y, b.circularity );
    hablob_convexity  ->Fill( b.x, 3000-b.y, b.convexity );
    hablob_inertia    ->Fill( b.x, 3000-b.y, b.inertia );

    ha_intensity  ->Fill( b.intensity );
    ha_area       ->Fill( b.area );
    ha_circularity->Fill( b.circularity );
    ha_convexity  ->Fill( b.convexity );
    ha_inertia    ->Fill( b.inertia );
  }

  for ( const EllipseData& el : idt.fEllipses ){
    for ( const BlobData& b : el.blobentry ){
      hblob_xy->Fill( b.x, 3000-b.y );

      he_intensity  ->Fill( b.intensity );
      he_area       ->Fill( b.area );
      he_circularity->Fill( b.circularity );
      he_convexity  ->Fill( b.convexity );
      he_inertia    ->Fill( b.inertia );


    }
    TEllipse* elli = new TEllipse( el.xx, 3000-el.yy, el.bb, el.aa, 0., 360., RADTODEG( el.phi ) );
    elli->SetFillStyle(0);
    hblob_xy->GetListOfFunctions()->Add( elli );
    TMarker *mark = new TMarker( el.xx, 3000-el.yy, 20 );
    hblob_xy->GetListOfFunctions()->Add( mark );

  }

  fout->Write();

}


void analyse_ttree(){

  gSystem->Load("libFeatureTreeLib.so");

  ImageData * imgdata = new ImageData();
  TFile     * fin     = new TFile("FindBoltLocation.root","read");
  TTree     * tree    = (TTree*) fin->Get("ImageTree");
  tree->SetBranchAddress( "ImageData", &imgdata );
  tree->GetEntry(0);

  std::cout<<"Nblobs ="<<imgdata->fBlobs.size()<<std::endl;
  std::cout<<"Nellipses ="<<imgdata->fEllipses.size()<<std::endl;
  for ( const EllipseData& el : imgdata->fEllipses ){
    std::cout<<"Ellipse x,y="<<el.xx<<", "<<el.yy<<std::endl;
    for ( const BlobData& b : el.blobentry ){
      std::cout<<" blob x,y="<<b.x<<", "<<b.y<<std::endl;
    }
  }


  basic_histograms( *imgdata );
}
