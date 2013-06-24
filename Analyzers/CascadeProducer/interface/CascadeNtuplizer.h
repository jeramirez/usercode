// -*- C++ -*-
//
// Package:    CascadeProducer
// Class:      CascadeNtuplizer
// 
/**\class CascadeNtuplizer CascadeNtuplizer.h Analyzers/CascadeProducer/interface/CascadeNtplizer.h

 Description: analyzer for ntuple cascade output 

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Juan Eduardo Ramirez Vargas
//         Created:  Mon Jun 29 17:24:25 CDT 2009
// $Id: CascadeNtuplizer.h,v 1.5 2011/11/25 03:07:19 jramirez Exp $
//
//

#ifndef RECOVERTEX__CASCADE_NTUPLE_H
#define RECOVERTEX__CASCADE_NTUPLE_H
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

class CascadeNtuplizer : public edm::EDAnalyzer {
   public:
      explicit CascadeNtuplizer(const edm::ParameterSet&);
      ~CascadeNtuplizer();


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

};

#endif
