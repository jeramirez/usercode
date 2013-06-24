// -*- C++ -*-
//
// Class:      BTagPAT1dHistos
// 
/**\class BTagPAT1dHistos Analyzers/BTagPAT/src/BTagPAT1dHistos.cc

 Description: <A histogram booker>

 Implementation:
 
 this class books 1d histograms per flavor. 
 default binning nbin =50, xmin = 0. xmax= 5. 
*/
//
// Original Author:  J.E. Ramirez
//
// $Id: BTagPAT1dHistos.cc,v 1.1 2008/07/25 23:25:24 jramirez Exp $

#include "Analyzers/BTagPAT/interface/BTagPAT1dHistos.h"

//
// constructors and destructor
//
BTagPAT1dHistos::BTagPAT1dHistos(const edm::ParameterSet& iConfig):
   nbins(iConfig.getUntrackedParameter<int>("nbins",50))
  ,xmin(iConfig.getUntrackedParameter<double>("xmin",0.))
  ,xmax(iConfig.getUntrackedParameter<double>("xmax",5.))
{
   //now do what ever initialization is needed
}

BTagPAT1dHistos::BTagPAT1dHistos(const int ptbins, const double ptmin, const double ptmax):
   nbins(ptbins)
  ,xmin(ptmin)
  ,xmax(ptmax)
{
   //now do what ever initialization is needed
}

BTagPAT1dHistos::~BTagPAT1dHistos()
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
BTagPAT1dHistos::Fill( double& ptrel, const std::string& flavor )
{
   histocontainer_[prefix+flavor]->Fill(ptrel);
}
// Fill histograms overload
void
BTagPAT1dHistos::Fill( int& iptrel, const std::string& flavor )
{
   histocontainer_[prefix+flavor]->Fill(iptrel);
}
//
// Set histograms
void
BTagPAT1dHistos::Set(const std::string& flavor, TFileDirectory& ptreldir, 
                        const std::string& inprefix, 
                        const std::string& histotitle)
{
  std::string histoid;
//  std::string histotitle;

  prefix  = inprefix;
  histoid = prefix + flavor;        //histotitle = "p_{Trel}"+flavor+" [GeV/c]";
  histocontainer_[histoid]=ptreldir.make<TH1D>(histoid.c_str(),histotitle.c_str(),nbins,xmin,xmax);
}
//
// Save histogram errors
void
BTagPAT1dHistos::Sumw2()
{
  for (std::map<std::string,TH1D*>::const_iterator ih=histocontainer_.begin();
         ih!= histocontainer_.end(); ++ih) {
      TH1D *htemp = (TH1D*) ih->second;
      htemp->Sumw2();
  }
}//end method Sumw2

