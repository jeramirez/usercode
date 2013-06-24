#ifndef PatAlgos_BTagPATClosureAnalyzer_H_
#define PatAlgos_BTagPATClosureAnalyzer_H_

// -*- C++ -*-
//
// Class:      BTagPATClosureAnalyzer
// Original Author:  J.E. Ramirez
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

#include "Analyzers/BTagPAT/interface/BTagPATCommonHistos.h"
#include "Analyzers/BTagPAT/interface/BTagPAT1dHistos.h"
#include "Analyzers/BTagPAT/interface/BTagPAT2dHistos.h"

//
// class declaration
//

class BTagPATClosureAnalyzer : public edm::EDAnalyzer {
   public:
      explicit BTagPATClosureAnalyzer(const edm::ParameterSet&);
      ~BTagPATClosureAnalyzer();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------

  edm::InputTag jetLabel_;              //Input PAT Jet Collection
  edm::InputTag muonLabel_;             //Input PAT Muon Collection
  bool    BTagverbose;                  //Flag for printouts

  std::map <std::string,std::map<std::string,double> > opCut;    //map to associate op and taggers
  std::vector<std::string> taggerslist;                          //vector of taggers names
  BTagPAT1dHistos  PtrelHistos;         //ptrel histograms.
  BTagPAT1dHistos  DeltaRHistos;        //DeltaR histograms.
  BTagPAT1dHistos  Chi2muHistos;        // chi2 distribuiion for muons un jet histograms.
  BTagPAT1dHistos  nmuonhistos;         //number of muons per jet histograms.
  BTagPAT1dHistos  nmuonhits;           //number of hits in muons tracks inside a jet histograms.
  BTagPAT2dHistos  scatterjetpt;        //jet pt  vs ptrel histograms.
  BTagPAT2dHistos  scatterjeteta;       //jet eta vs ptrel histograms.

  std::string BTagtagger_;              //jet away tagger
  std::string BTagpurity_;              //jet away op (Loose,Medium,Tight)
  PATbJetsSelector BTagger;             //helper function (isbtag)

  //jet cut variables
  double MinJetPt_, MaxJetEta_, MaxDeltaR_, MinPtRel_;  
  //muon cut variables
  int MinMuHits_;
  double MinMuPt_,MaxMuEta_;
};

#endif
