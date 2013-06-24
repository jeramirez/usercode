// -*- C++ -*-
//
// Class:      BTagPAT2dHistos
// 
/**\class BTagPAT2dHistos Analyzers/BTagPAT/src/BTagPAT2dHistos.cc

 Description: <A histogram booker>

 Implementation:
 
 this class books 2d histograms per flavor. 
 default binning nxbin =50, xmin = 0. xmax= 5. 
 default binning nybin =50, ymin = 0. ymax= 5. 
*/
//
// Original Author:  J.E. Ramirez
//
// $Id: BTagPAT2dHistos.cc,v 1.1 2008/07/25 23:25:16 jramirez Exp $

#include "Analyzers/BTagPAT/interface/BTagPAT2dHistos.h"

//
// constructors and destructor
//
BTagPAT2dHistos::BTagPAT2dHistos(const edm::ParameterSet& iConfig):
   nxbin(iConfig.getUntrackedParameter<int>("nxbin",50))
  ,xmin(iConfig.getUntrackedParameter<double>("xmin",0.))
  ,xmax(iConfig.getUntrackedParameter<double>("xmax",5.))
  ,nybin(iConfig.getUntrackedParameter<int>("nybin",50))
  ,ymin(iConfig.getUntrackedParameter<double>("xmin",0.))
  ,ymax(iConfig.getUntrackedParameter<double>("xmax",5.))
  ,dim(iConfig.getUntrackedParameter<int>("xdim",11))
{
   //now do what ever initialization is needed
}

BTagPAT2dHistos::BTagPAT2dHistos(const int pnxbin, const double pxmin, const double pxmax,
                                 const int pnybin, const double pymin, const double pymax
):
   nxbin(pnxbin)
  ,xmin(pxmin)
  ,xmax(pxmax)
  ,nybin(pnybin)
  ,ymin(pymin)
  ,ymax(pymax)
  ,dim(pnxbin)
{
   //now do what ever initialization is needed
}

BTagPAT2dHistos::BTagPAT2dHistos(const int pdim,
                                 const int pnybin, const double pymin, const double pymax
):
   nxbin(0)
  ,xmin(0)
  ,xmax(0)  
  ,nybin(pnybin)
  ,ymin(pymin)
  ,ymax(pymax)
  ,dim(pdim)
{
   //now do what ever initialization is needed
}

BTagPAT2dHistos::~BTagPAT2dHistos()
{ 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//
//---------------------
// Fill histograms
void
BTagPAT2dHistos::Fill( double& x, double& y, const std::string& flavor )
{
   histocontainer_[prefix+flavor]->Fill(x,y);
}
//
// Set histograms
void
BTagPAT2dHistos::Set(const std::string& flavor, TFileDirectory& ptreldir, 
                     const std::string& inprefix, 
                     const std::string& histotitle)
{
  std::string histoid;
//  std::string histotitle;

  prefix  = inprefix;
  histoid = prefix + flavor;        //histotitle = "p_{Trel}"+flavor+" [GeV/c]";
  histocontainer_[histoid]=ptreldir.make<TH2D>(histoid.c_str(),histotitle.c_str(),nxbin,xmin,xmax,nybin,ymin,ymax);
}
//
// Set histograms (input vectors)
void
BTagPAT2dHistos::Set(const std::string& flavor, 
                     TFileDirectory& ptreldir, 
                     double *xbins,
                     const std::string& inprefix, 
                     const std::string& histotitle
                     )
{
  std::string histoid;
//  std::string histotitle;

  prefix  = inprefix;
  histoid = prefix + flavor;        //histotitle = "p_{Trel}"+flavor+" [GeV/c]";
  histocontainer_[histoid]=ptreldir.make<TH2D>(histoid.c_str(),histotitle.c_str(),dim-1,xbins,nybin,ymin,ymax);
}
//
// Save histogram errors
void
BTagPAT2dHistos::Sumw2()
{
  for (std::map<std::string,TH2D*>::const_iterator ih=histocontainer_.begin();
         ih!= histocontainer_.end(); ++ih) {
      TH2D *htemp = ih->second;
      htemp->Sumw2();
  }
}//end method Sumw2

