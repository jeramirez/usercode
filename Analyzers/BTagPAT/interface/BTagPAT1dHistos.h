#ifndef PatAlgos_BTagPAT1dHistos_H_
#define PatAlgos_BTagPAT1dHistos_H_

// -*- C++ -*-
//
// Package:    BTagPAT
// Class:      BTagPAT1dHistos
// 
/**\class BTagPAT1dHistos BTagPAT1dHistos.h

 Description: <Define and Fill a set of 1d histograms depending on flavor>

 Implementation:
 
 Create a container of histograms. 
*/
//
// Original Author:  J. E. Ramirez
//
// $Id: BTagPAT1dHistos.h,v 1.1 2008/07/25 23:23:23 jramirez Exp $


// system include files
#include <memory>

// user include files
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "PhysicsTools/UtilAlgos/interface/TFileDirectory.h"
#include <string>
#include <map>
#include "TH1D.h"

//
// class declaration
//

class BTagPAT1dHistos {
   public:
      explicit BTagPAT1dHistos(const edm::ParameterSet& iConfig);
      BTagPAT1dHistos(const int ptnbins, const double ptxmin, const double ptxmax);
      ~BTagPAT1dHistos();

      void Set(const std::string& flavor, TFileDirectory& ptreldir, 
               const std::string& inprefix="ptrel_", 
               const std::string& histotitle= "p_{Trel} [GeV/c]");
      void Sumw2();
      void Fill(double&, const std::string& );
      void Fill(int&, const std::string& );
      void SetPrefix(const std::string& inprefix="ptrel_"){ prefix= inprefix;};
   private:

      // ----------member data ---------------------------

  const int nbins;
  const double xmin,xmax;
  std::map<std::string,TH1D*> histocontainer_;   // simple map to contain all histograms. Histograms are booked in the beginJob() method
  std::string prefix;

};

#endif
