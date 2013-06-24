#ifndef HecMesons_h
#define HecMesons_h

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
/*
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
*/
#include "TFile.h"
#include "TList.h"
#include "TH1D.h"
#include "TH1I.h"
#include "TH2D.h"
#include "TTree.h"
#include "TNtuple.h"


// class decleration
//

class HecMesons : public edm::EDAnalyzer {
   public:
      explicit HecMesons(const edm::ParameterSet&);
      ~HecMesons();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------
      edm::InputTag trackTags_;         //used to select what tracks to read from configuration file
      edm::InputTag electronCollection_;
      unsigned int minTracks_;
      bool   check_elC_;
      double Pe_match_;
      double mass_uno_;
      TH1D * histo_trk;
      TH1D * histo_Pti;
      TH1D * histo_ij; 
      TH1D * histo_ijPt; 
      TH1D * histo_ij0;
      TH1D * histo_ij1;
      TH1D * histo_ij2;
      TH1D * histo_elec;
      TH1D * histo_NeTrk;
      TH1D * histo_DNel;
      TH1D * histo_Delp;
      TH1D * histo_DelpC;
      TH1D * histo_Delq;
      TH1D * histo_eiej;
      TH1D * histo_DM;
      TH1D * histo_eEn;
      TH1D * histo_EoP;
      TH1D * histo_Clu;
      TH1D * histo_Eni;
};

//
#endif
