#include "hough_ellipse.hpp"
#include <sstream>
#include <fstream>
#include "TDirectory.h"
#include "TEllipse.h"
#include "TMarker.h"


//const double pi = std::acos(-1);
float RADTODEG( float rad ){ return rad*180/pi; }

void process_mem_usage(float& vm_usage, float& resident_set)
{
  vm_usage     = 0.0;
  resident_set = 0.0;

  // the two fields we want
  unsigned long vsize;
  long rss;
  {
    std::string ignore;
    std::ifstream ifs("/proc/self/stat", std::ios_base::in);
    ifs >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore
	>> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore
	>> ignore >> ignore >> vsize >> rss;
  }

  long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
  vm_usage = vsize / 1024.0;
  resident_set = rss * page_size_kb;
}



unsigned num_ellipses( const HoughEllipseResults& hrs) {
  unsigned nc=0;
  for ( const HoughEllipseResult& hr : hrs ){
    if ( hr.type == HoughEllipse ) ++nc;
  }
  return nc;
}

std::ostream& operator<<( std::ostream& os, const HoughEllipseResults& hrs ){
  for ( const HoughEllipseResult& hr : hrs ){
    os << hr;
  }
  return os;
}

std::ostream& operator<<( std::ostream& os, const HoughEllipseResult& hr ){
  if ( hr.type == HoughEllipse ){
    os<<"Hough Ellipse: (xc,yc)= ( "<<hr.e.get_xy().x<<", "<<hr.e.get_xy().y
      <<" )  bb="<<hr.e.get_b()
      <<" eccentricity="<<hr.e.get_e()
      <<" phi="<<hr.e.get_phi()
      <<" hough-peak="<<hr.peakval<<std::endl;
  } else {
    os<<"Hough Unused Hits::"<<std::endl;
  }
  os << "  Nhits="<<hr.data.size()<<std::endl;
  for ( const xypoint& xy : hr.data ){
    os << "    (x,y)= ( "<<xy.x<<", "<<xy.y<<" )"<<std::endl;
  }
  return os;
}


EllipseHough::EllipseHough( unsigned nbins_bb     , float bbmin   , float bbmax     ,
			    unsigned nbins_ee     , float eemin   , float eemax    ,
			    unsigned nbins_phiphi , float phiphimin, float phiphimax,
			    unsigned nbins_x      , float xmin    , float xmax     ,
			    unsigned nbins_y      , float ymin    , float ymax     ) :
  fNbb( nbins_bb )     , fbbmin( bbmin )     , fbbmax( bbmax ),
  fNee( nbins_ee )     , feemin( eemin )     , feemax( eemax ),
  fNphi( nbins_phiphi ), fphimin( phiphimin ), fphimax( phiphimax ), 
  fNx( nbins_x ), fxmin( xmin ), fxmax( xmax ),
  fNy( nbins_y ), fymin( ymin ), fymax( ymax )					   
{
  static unsigned instance_count=0;
  ++instance_count;
  std::ostringstream os;

  TDirectory * curdir = gDirectory;
  houghdir = gDirectory->mkdir( (std::string("ehough_")+std::to_string(instance_count)).c_str() );
  houghdir->cd();

  // loop over bb bins to make directories
  float dbb = (bbmax-bbmin)/nbins_bb;
  
  for ( unsigned ibb=0; ibb<nbins_bb; ++ibb){
    std::cout<<"EllipseHough::EllipseHough building directory "<<ibb<<" / "<<nbins_bb<<std::endl;
    float vmuse =0., memuse=0.;
    process_mem_usage( vmuse, memuse );
    std::cout<<" vmuse = "<<vmuse/1024/1024<<" GB, memuse = "<<memuse/1024/1024<<" GB "<<std::endl;
    float bbm = bbmin + dbb*ibb;
    float bbp = bbmin + dbb*(ibb+1);
    float bb = (bbm+bbp)/2;

    std::string bdirname  = std::string("bb_")+std::to_string(bb);
    std::string bdirtitle = std::string("bb = ")+std::to_string(bb);
    TDirectory* bbdir = houghdir->mkdir( bdirname.c_str(), bdirtitle.c_str() );
    bbdir->cd();
    fDirectories.push_back( bbdir );
    fTransformed.push_back( std::vector< TH2S* >() );
    fE.push_back( std::vector< ellipse_st* >() );

    // loop over ee and phiphi bins to make histograms
    float dee = (eemax-eemin)/nbins_ee;
    for ( unsigned iee=0; iee<nbins_ee; ++iee ){
      float eem = eemin + dee*iee;
      float eep = eemin + dee*(iee+1);
      float ee  = (eem+eep)/2;
      float dphi = (phiphimax-phiphimin)/nbins_phiphi;
      for ( unsigned iphi=0; iphi<nbins_phiphi; ++iphi ){
	float phim = phiphimin + dphi*iphi;
	float phip = phiphimin + dphi*(iphi+1);
	float phi  = (phim+phip)/2;

	fE[ibb].push_back( new ellipse_st( bb, ee, phi, xypoint(0.,0.) ) );
	
	std::string hname = std::string( "ht") + std::to_string(instance_count) 
	  + "_bbin"+std::to_string(ibb) + "_eebin"+std::to_string(iee)
	  + "_phibin"+std::to_string(iphi);
	std::string htitle = 
	  std::string( "bb = ( ") + std::to_string( bbm ) + ", " + std::to_string( bbp ) + ") " +
	  std::string( "e = ( ") + std::to_string( eem ) +  ", " + std::to_string( eep ) + ") " +
	  std::string( "phi = ( ") + std::to_string( phim ) +  ", " + std::to_string( phip ) + ") " + 
	  " ; xc; yc ";

	
	fTransformed[ibb].push_back( new TH2S( hname.c_str(), htitle.c_str(),
					       nbins_x, xmin, xmax,
					       nbins_y, ymin, ymax )  );
    
      }
    }
  }
  curdir->cd();

  float vmuse =0., memuse=0.;
  process_mem_usage( vmuse, memuse );
  std::cout<<"Done. vmuse = "<<vmuse/1024/1024<<" GB, memuse = "<<memuse/1024/1024<<" GB "<<std::endl;
}

EllipseHough::~EllipseHough(){
  for ( std::vector< TH2S* > fT : fTransformed){
    for ( unsigned i=0; i<fT.size(); ++i ){
      //fT[i]->SetDirectory(0);
      //delete fT[i];
    }
    fT.clear();
  }
  fTransformed.clear();
  for ( std::vector< ellipse_st* > est : fE ){
    for ( ellipse_st* e : est ){
      if ( e ) delete e;
    }
    est.clear();
  }
  fE.clear();

}


const HoughEllipseResults& EllipseHough::find_ellipses( const std::vector< xypoint >& data ){
  //float fxbwid = (fxmax-fxmin)/fNx; 
  //float fybwid = (fymax-fymin)/fNy;

  std::cout<<"EllipseHough::find_ellipses on data with "<<data.size()<<" points"<<std::endl;
  fresults.clear();
  std::vector< xypoint > unused_hits = data;
  bool done = false;
  while ( unused_hits.size() > minhits && !done ){
    hough_transform( unused_hits );
    float vmuse =0., memuse=0.;
    process_mem_usage( vmuse, memuse );
    std::cout<<" vmuse = "<<vmuse/1024/1024<<" GB, memuse = "<<memuse/1024/1024<<" GB "<<std::endl;



    HoughEllipseResult hr = find_maximum( unused_hits );
    //HoughEllipseResult hr;
    //do {
    //hr = find_maximum( unused_hits );

    if ( hr.peakval > threshold ) {
      hr.type = HoughEllipse;
    } else {
      done = true;
    }
    
    std::cout<<"Find ellipse "<<fresults.size()+1<<std::endl;
    std::cout<< hr <<std::endl;

    process_mem_usage( vmuse, memuse );
    std::cout<<" vmuse = "<<vmuse/1024/1024<<" GB, memuse = "<<memuse/1024/1024<<" GB "<<std::endl;


    fresults.push_back( hr );
    //if ( done ) break;
    save_hough_histo( fresults.size(), fTransformed[ hr.ibb ][ hr.ihist ] );
    plot_candidate( fresults.size(), hr );
      
      // remove current maximum.
      //for ( unsigned ibb = 0; ibb < fTransformed.size(); ++ibb ) {
      //	for ( unsigned ihist = 0; ihist < fTransformed[ibb].size(); ++ihist ) {
      //	  for ( int dx=-2; dx<3; ++dx ){
      //	    for ( int dy=-2; dy<3; ++dy ){
      //	      int ibin = fTransformed[ ibb ][ ihist ]->FindBin( hr.e.get_xy().x + dx*fxbwid, hr.e.get_xy().y + dy*fybwid ); 
      //	      fTransformed[ ibb ][ ihist ]->SetBinContent( ibin, 0 );
      //	    }
      //	  }
      //	} 
      //}
      
      //} while ( hr.peakval > threshold );
  }

  // cleanup
  //for ( unsigned i=0; i<fTransformed.size(); ++i ){
  //  fTransformed[i]->SetDirectory(0);
    //delete fTransformed[i];
  //}

  return fresults;
}

void EllipseHough::hough_transform( const std::vector< xypoint >& data ){

  std::cout<<"EllipseHough::hough_transform call on data with "<<data.size()<<" points"<<std::endl;

  float bbwid = (fbbmax-fbbmin)/fNbb;
  float fxbwid = (fxmax-fxmin)/fNx; 
  float fybwid = (fymax-fymin)/fNy;
  
  for ( std::vector<TH2S*> hists : fTransformed){
    for ( TH2S* h : hists ){
      h->Reset();
    }
  }
  
  float vmuse =0., memuse=0.;
  process_mem_usage( vmuse, memuse );
  std::cout<<"histos-reset vmuse = "<<vmuse/1024/1024<<" GB, memuse = "<<memuse/1024/1024<<" GB "<<std::endl;
  
  assert( fE.size() == fTransformed.size() );


  for ( const xypoint& xy : data ){
    for ( unsigned ibb=0; ibb < fE.size(); ++ibb ){
      assert( fE[ibb].size() == fTransformed[ibb].size() );
      for ( unsigned ihist=0; ihist < fE[ibb].size(); ++ihist ) {
	ellipse_st * elli = fE[ibb][ihist]; // has set everything but x and y

	elli->set_xy( xy );

	float bb = elli->get_b();
	// pick number of angles based on bb
	unsigned nang = unsigned( 2 * bb / bbwid );
	float   dtheta = 2*pi/nang;
	for ( unsigned itheta = 0; itheta<nang; ++itheta ){
	  xypoint ab = elli->xy( itheta*dtheta );

	  float a = ab.x;
	  float b = ab.y;
	  short weight = 1;
	  //fTransformed[ ibb ][ ihist ]->Fill( a, b);
	  for ( int xx=-1; xx<2; ++xx ){
	    for ( int yy=-1; yy<2; ++yy ){
	      weight = 4-2*abs(xx)-2*abs(yy)+abs(xx*yy);
	      float curx = a+xx*fxbwid;
	      float cury = b+yy*fybwid;
	      if ( curx > fxmin && curx < fxmax &&
		   cury > fymin && cury < fymax ){
		fTransformed[ ibb ][ ihist ] -> Fill( curx, cury, weight );
	      }
	    }
	  }
	}
      }
    }
  }
}


/*
struct binindices_st {
  unsigned ibb;
  unsigned ihist;
  unsigned ix;
  unsigned iy;
  binindices_st() : ibb(0), ihist(0), ix(0), iy(0) { } 
};
*/

//void EllipseHough::find_maximum( std::vector< xypoint >& hits, std::vector< HoughEllipseResult > & result ){
HoughEllipseResult EllipseHough::find_maximum( std::vector< xypoint >& hits ){

  //std::vector< binindices_st > passed_threshold;

  HoughEllipseResult curbest( ellipse_st(), 0 );
  // loop over the hough transform histograms to find the peak
  // and store the "best" circle center and radius
  for ( unsigned ibb=0; ibb < fTransformed.size(); ++ibb ){
    for ( unsigned ihist=0; ihist < fTransformed[ibb].size(); ++ihist ) {
      //float vmuse =0., memuse=0.;
      //process_mem_usage( vmuse, memuse );
      //std::cout<<"find_maximum ibb="<<ibb<<" ihist="<<ihist<<" vmuse = "<<vmuse/1024/1024<<" GB, memuse = "<<memuse/1024/1024<<" GB "<<std::endl;
      TH2S* h = fTransformed[ibb][ihist];

      for ( int ix=1; ix<=h->GetNbinsX(); ++ix ){
	for ( int iy=1; iy<=h->GetNbinsY(); ++iy ){
	  if ( h->GetBinContent( ix, iy ) > curbest.peakval ) {
	    float x = h->GetXaxis()->GetBinCenter( ix );
	    float y = h->GetYaxis()->GetBinCenter( iy );
	    curbest.peakval = h->GetBinContent( ix, iy );
	    curbest.ibb = ibb;
	    curbest.ihist = ihist;
	    ellipse_st * elli = fE[ curbest.ibb ][ curbest.ihist ]; // has set everything but x and y
	    elli->set_xy( xypoint( x, y ) );
	    curbest.e = *elli;
	  }
	}
      }
      //process_mem_usage( vmuse, memuse );
      //std::cout<<"endloop find_maximum ibb="<<ibb<<" ihist="<<ihist<<" vmuse = "<<vmuse/1024/1024<<" GB, memuse = "<<memuse/1024/1024<<" GB "<<std::endl;
    }
  }

  // find the hits that are associated with the ellipse and add them
  // to the result, otherwise add them to list of unused_hits
  // use bin sizes in xc, yc as threshold distance for hit to be from circle
  float bbwid = (fbbmax-fbbmin)/fNbb;
  float fxbwid = (fxmax-fxmin)/fNx; 
  float fybwid = (fymax-fymin)/fNy;
  float rthres = drscaling * std::sqrt( fxbwid*fxbwid + fybwid*fybwid );
  std::cout<<"rthres="<<rthres<<std::endl;
  std::cout<<"curbest.e "<<curbest.e<<std::endl;
  std::vector< xypoint > unused_hits;
  for ( xypoint xy : hits ){
    float dr  = curbest.e.dmin( xy );
    if ( fabs( dr ) < rthres ){
      curbest.data.push_back( xy );
    } else {
      unused_hits.push_back( xy );
    }
  }

  hits = unused_hits;
  return curbest;
}

  
void EllipseHough::save_hough_histo( unsigned num, TH2S* histo ){
  TDirectory* curdir = gDirectory;
  houghdir->cd();

  std::string hname = std::string( histo->GetName() ) + "_cand_" + std::to_string( num );
  TH2S* savehist = (TH2S*)histo->Clone( hname.c_str() );
  savehist->SetName( hname.c_str() );
  savehist->SetDirectory( houghdir );
  curdir->cd();
}

void EllipseHough::plot_candidate( unsigned num, const HoughEllipseResult & hr ){
  TDirectory* curdir = gDirectory;
  houghdir->cd();
  
  std::string hname = std::string("hcircle_")+std::to_string(num);
  TH2S* hplot = (TH2S*)fTransformed[ hr.ibb ][ hr.ihist ]->Clone( hname.c_str() );
  hplot->SetName( hname.c_str() );
  hplot->SetDirectory( houghdir );
  hplot->Reset();

  for ( const xypoint& xy : hr.data ){
    hplot->Fill( xy.x, xy.y );
  }

  TEllipse *el = new TEllipse( hr.e.get_xy().x, hr.e.get_xy().y, hr.e.get_b(), 
			       hr.e.get_a(), 0., 360., RADTODEG( hr.e.get_phi() ) );
  el->SetFillStyle(0);
  hplot->GetListOfFunctions()->Add( el );
  TMarker *mark = new TMarker( hr.e.get_xy().x, hr.e.get_xy().y, 20 );
  hplot->GetListOfFunctions()->Add( mark );

  curdir->cd();
}

