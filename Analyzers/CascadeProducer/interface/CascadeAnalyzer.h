// -*- C++ -*-
//
// Package:    CascadeProducer
// Class:      CascadeAnalyzer
// 
/**\class CascadeAnalyzer CascadeAnalyzer.h Analyzers/CascadeProducer/interface/CascadeAnalyzer.h

 Description: analyzer for histograming cascade output 

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Juan Eduardo Ramirez Vargas
//         Created:  Mon Jun 29 17:24:25 CDT 2009
// $Id: CascadeAnalyzer.h,v 1.6 2011/06/23 17:43:13 jramirez Exp $
//
//

#ifndef RECOVERTEX__CASCADE_ANALYZER_H
#define RECOVERTEX__CASCADE_ANALYZER_H
// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "Analyzers/CascadeProducer/interface/Histos1d.h"
//
// class declaration
//

class CascadeAnalyzer : public edm::EDAnalyzer {
   public:
      explicit CascadeAnalyzer(const edm::ParameterSet&);
      ~CascadeAnalyzer();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------
  //Input parameters
  std::string CasAlgo;           //Label of Cascade collection
  std::string VeeAlgo;           //Label of Vee collection
  std::string CasDecayName;      //Name of Cascade Collection (Cascade/Omega)
  edm::InputTag tracksAlgo;      //Label of TrackCollection
  Histos1d  masshisto;           //cascade invariant mass histograms.
  Histos1d  massrejhisto;        //cascade invariant mass histograms rejected because vtx upstream of last track hits 
  Histos1d  veemasshisto;        //Lambda invariant mass histograms.
  Histos1d  chisqhisto;          //chisq histograms.
  Histos1d  clhisto;             //cl histograms.
  const edm::ParameterSet iConfig;
  const double pitrk_pt_cut_;    //pt cut for trk daughter
};

#endif
