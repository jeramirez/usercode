// -*- C++ -*-
//
// Class:      Histos1d
// 
/**\class Histos1d Analyzers/CascadeProducer/src/Histos1d.cc

 Description: <A histogram booker>

 Implementation:
 
 this class books 1d histograms per flavor. 
 default binning nbin =50, xmin = 0. xmax= 5. 
*/
//
// Original Author:  J.E. Ramirez
//
// $Id: Histos1d.cc,v 1.2 2011/06/15 19:41:00 jramirez Exp $

#include <iostream>
#include "Analyzers/CascadeProducer/interface/Histos1d.h"

//
// constructors and destructor
//
Histos1d::Histos1d(const edm::ParameterSet& iConfig):
   nbins(iConfig.getUntrackedParameter<int>("nbins",50))
  ,xmin(iConfig.getUntrackedParameter<double>("xmin",0.))
  ,xmax(iConfig.getUntrackedParameter<double>("xmax",5.))
{
   //now do what ever initialization is needed
}

Histos1d::Histos1d(int ptbins, double ptmin, double ptmax):
   nbins(ptbins)
  ,xmin(ptmin)
  ,xmax(ptmax)
{
   //now do what ever initialization is needed
}

Histos1d::~Histos1d()
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
Histos1d::Fill( double& dvalue, const std::string& idhisto )
{
   if ( histocontainer_.find(prefix+idhisto) == histocontainer_.end() ){
      std::cout << "histogram="
                << prefix+idhisto
                << " NOT defined"
                << std::endl;
   }
   histocontainer_[prefix+idhisto]->Fill(dvalue);
}
// Fill histograms overload
void
Histos1d::Fill( int& ivalue, const std::string& idhisto )
{
   if ( histocontainer_.find(prefix+idhisto) == histocontainer_.end() ){ 
      std::cout << "histogram=" 
		<< prefix+idhisto
		<< " NOT defined"
		<< std::endl;
   }
   histocontainer_[prefix+idhisto]->Fill(ivalue);
}
//
// Set histograms
void
Histos1d::Set(const std::string& idhisto, TFileDirectory& dir, 
                        const std::string& inprefix, 
                        const std::string& histotitle)
{
  std::string histoid;

  prefix  = inprefix;
  histoid = prefix + idhisto;        //histotitle = "prefix"+idhisto+" [GeV/c]";
  histocontainer_[histoid]=dir.make<TH1D>(histoid.c_str(),histotitle.c_str(),nbins,xmin,xmax);
}
//
// Save histogram errors
void
Histos1d::Sumw2()
{
  for (std::map<std::string,TH1D*>::const_iterator ih=histocontainer_.begin();
         ih!= histocontainer_.end(); ++ih) {
      TH1D *htemp = (TH1D*) ih->second;
      htemp->Sumw2();
  }
}//end method Sumw2

