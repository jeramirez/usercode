#ifndef PatAlgos_BTagPATCommonHistos_H_
#define PatAlgos_BTagPATCommonHistos_H_

// -*- C++ -*-
//
// Package:    BTagPAT
// Class:      BTagPATCommonHistos
// 
/**\class BTagPATCommonHistos BTagPATCommonHistos.h

 Description: <Define and Fill common set of histograms depending on flavor and tagger>

 Implementation:
 
 Create a container of histograms. 
*/
//
// Original Author:  J. E. Ramirez
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Common/interface/View.h"
#include "Analyzers/BTagPAT/interface/PATbJetsSelector.h"

#include "TH1D.h"
#include "TH2D.h"
#include <map>

//
// class declaration
//

class BTagPATCommonHistos {
   public:
      explicit BTagPATCommonHistos(const edm::ParameterSet&);
      ~BTagPATCommonHistos();

      void Set(std::string);
      void Sumw2();
      void Fill(edm::View<pat::Jet>::const_iterator&, std::string);
   private:

      // ----------member data ---------------------------


  std::string  flavor;
  std::map<std::string,TH1D*> histocontainer_;   // simple map to contain all histograms. Histograms are booked in the beginJob() method
  std::map<std::string,TH2D*> h2_; // simple map to contain 2D  histograms. Histograms are booked in the beginJob() method
  std::string  BTagdiscriminator_;
  std::string  BTagpurity_;
  double  BTagdisccut_;
  PATbJetsSelector BTagger;
};

#endif
