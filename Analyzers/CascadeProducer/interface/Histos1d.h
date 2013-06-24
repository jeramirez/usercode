#ifndef PatAlgos_Histos1d_H_
#define PatAlgos_Histos1d_H_

// -*- C++ -*-
//
// Package:    CascadeProducer
// Class:      Histos1d
// 
/**\class Histos1d Histos1d.h

 Description: <Define and Fill a set of 1d histograms>

 Implementation:
 
 Create a container of histograms. 
*/
//
// Original Author:  J. E. Ramirez
//
// $Id: Histos1d.h,v 1.4 2011/06/15 19:41:00 jramirez Exp $


// system include files
#include <memory>

// user include files
#include "FWCore/ParameterSet/interface/ParameterSet.h"
//3_3_6 include "PhysicsTools/UtilAlgos/interface/TFileDirectory.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"
#include <string>
#include <map>
#include "TH1D.h"

//
// class declaration
//

class Histos1d {
   public:
      explicit Histos1d(const edm::ParameterSet& iConfig);
      Histos1d(int ptnbins, double ptxmin, double ptxmax);
      ~Histos1d();

      void Set(const std::string& idhisto, TFileDirectory& dir, 
               const std::string& inprefix="prefix_", 
               const std::string& histotitle= "Histo_{title} [GeV/c]");
      void Sumw2();
      void Fill(double&, const std::string& );
      void Fill(int&, const std::string& );
      void SetPrefix(const std::string& inprefix="prefix_"){ prefix= inprefix;};
      bool IsValid(const std::string& idhisto ){
           bool histoisOK=true;
           if ( histocontainer_.find(prefix+idhisto) == histocontainer_.end() ){
              histoisOK=false;
              std::cout << "histogram="
                        << prefix+idhisto
                        << " NOT defined"
                        << std::endl;}
           return histoisOK;};
  private:

      // ----------member data ---------------------------

  const int nbins;
  const double xmin,xmax;
  std::map<std::string,TH1D*> histocontainer_;   // simple map to contain all histograms. Histograms are booked in the beginJob() method
  std::string prefix;

};

#endif
