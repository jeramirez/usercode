// -*- C++ -*-
//
// Package:    KinkProducer
// Class:      TestKinkNtuplizer
// 
/**\class TestKinkNtuplizer TestKinkNtuplizer.h Analyzers/KinkProducer/plugins/TestKinkNtplizer.h

 Description: analyzer for ntuple cascade with kinksfinder output 

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Juan Eduardo Ramirez Vargas
//         Created:  Mon Jun 29 17:24:25 CDT 2009
// $Id: TestKinkNtuplizer.h,v 1.1 2012/07/17 21:06:18 jramirez Exp $
//
//

#ifndef KINKTEST___NTUPLE_H
#define KINKTEST___NTUPLE_H
// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TTree.h"
//
// class declaration
//

class TestKinkNtuplizer : public edm::EDAnalyzer {
   public:
      explicit TestKinkNtuplizer(const edm::ParameterSet&);
      ~TestKinkNtuplizer();


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
  const edm::ParameterSet iConfig;
  const double pitrk_pt_cut_;    //pt cut for trk daughter
  const double CasMass_;         //Nominal mass for Xi/Om

  TTree*      tree_;
  double  ximass;
  double  rho_cas;
  double  los;
  double  losnocorr;
  double  probvxi;
  double  probvxi2;
  double  Xiip3d;
  double  bestip3d;
  double  bestip3d2;
  double  bestcl;
  int     hasXiTrk;
  int     hasXiTrk2; //from CascadePixelTrackFinder
//trk
  int trkhits;
  int pixhits;
  int hitsbeforev;
  int hitsafterv;
  int pxlHitsBe;
  int pxlHitsAf;
  double  rho_trk;  
  double  pionip3d;
  double  pionpt;

//vee
  double  lbmass;
  double  ksmass;

  double probvee;
  double pionveeip3d;
  double protonveeip3d; 
//prm  
  double NewPrimaryProb;    //prmCL
  int nprm;                 //number of primaries
};

#endif
