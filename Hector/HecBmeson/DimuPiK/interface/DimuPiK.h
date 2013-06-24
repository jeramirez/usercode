// -*- C++ -*-
//
// Package:    DimuPiK
// Class:      DimuPiK
// 
/**\class DimuPiK DimuPiK.cc HecMeson/DimuPiK/interface/DimuPiK.h

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Hector Mendez
//         Created:  Tue Oct  2 09:00:11 CDT 2012
// $Id: DimuPiK.h,v 1.1 2013/03/18 19:51:01 mendez Exp $
//
//


#ifndef DimuPiK_h
#define DimuPiK_h

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/PatCandidates/interface/Muon.h"

#include "Analyzers/CascadeProducer/interface/KinematicFitDriver.h"

#include "TFile.h"
#include "TList.h"
#include "TH1D.h"
#include "TH1I.h"
#include "TH2D.h"
#include "TTree.h"
#include "TNtuple.h"


//
// class declaration
//

class DimuPiK : public edm::EDAnalyzer {
   public:
      explicit DimuPiK(const edm::ParameterSet&);
      ~DimuPiK();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   protected:
      void fill_evt(const edm::Event&, int, int, int, int,
                                       int, int, int);
      void fill_IuuC();
      void fill_uuC(double);
      void fill_kPi(int);
      void fill_kKa(int);
      void fill_nPi(int);
      void fill_rPi(int);
      void fill_vtxAll(double, double, double, double, double,
                       double, double, double, double, double,
		       double, double, double, double);
      void fill_k0short(double, double, double, double, double, double, double,
                        double, double, double,
			double, double, double,
                        double, double, double, double, double);  
      void fill_tightMuons();
			   
      void fill_softMuons();
			  
      void fill_tightKMuons();
      
      void fill_softKMuons();
			  
      void init();

      void init(int);
      
      void   HecMuVar(int,std::vector<pat::Muon>::const_iterator, reco::TrackRef);
      bool   HecClosest(int, FreeTrajectoryState, FreeTrajectoryState);
      void   HecHltTrig(const edm::Event&);
      void   HecHltMuTrig(int,const pat::Muon*);
      //bool   HecDimuVtx(double, std::vector<RefCountedKinematicParticle>,reco::Particle::LorentzVector);
      bool   HecDimuVtx(double, KinematicFitDriver,reco::Particle::LorentzVector);
      void   HecTrkVar(int, double, reco::TrackRef);
      double HecRcone(int, double, double);
      
      
      void HecCutMu(int,std::vector<pat::Muon>::const_iterator);
     //Hector
       double myP[11][11]; //--First index: 0=B,     Second index 0=En  
                           //--             1=Psip                1=px
			   //--		    2=iMuon               2=py
			   //--		    3=jMuon               3=pz
			   //--		    4=Pion                4=Inv. Mass
			   //--		    5=Kaon                5=Charge
                           //--             6=K_short             6=pT
			   //--             7=iPion               7=eta
			   //--	            8=jPion.              8=Gen. mass
			   //--                                   9=Id
			   //--                                   10= Mass Diff. 
      void calcMyP(int n, const reco::Candidate *Cand);
   
   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------
      const edm::ParameterSet iConfig;
      edm::InputTag trackTags_; //used to select what tracks to read from configuration file
      edm::InputTag theMuonsLabel_;      //--here I get both global & track Mu [stolen from Heavy Oni2MuMu]
      std::string VeeAlgo_;              //--Label of Lambda & Kshort Vee collection
      bool MyPrint_, ChargedKa_, NeutralKa_, MyPrintMC_, doMC_;
      std::string hlTriggerResults_;
      
      TH1D *histo_trk, *histo_nmu, *histo_muId, *histo_primVtx, *histo_kshort, *histo_lamb;
      TH1D *histo_dimu[3], *histo_dimuPi[6], *histo_dimuPiK[2], *histo_k0s[5];
      TH1D *histo_closest, *histo_MCgen, *histo_NLambda0, *histo_L0daug;
      TH1D *histo_B0Energy, *histo_massB0, *histo_massZ4430, *histo_massKaon, *histo_masspion, *histo_masspsip;
      TH1D *histo_massmuimuj, *histo_massmuijpi, *histo_massBch, *histo_massKsh;
      TH1D *histo_combB,*histo_combpsip,*histo_combpi ,*histo_combKch,*histo_combKsh, *histo_comb3p;
      TH1D *histo_primCL, *histo_primCLHM, *histo_primDz, *histo_primDxy, *histo_primDxyz;
      //--Variables for my functions
      double MuVar[2][13];
      int   nMuVar[2][5];
      double dcaVar[6][4];  //--see function HecClosest() Implementation
      int nTrig[2][2];
      int evtTrig1, evtTrig2, evtTrig3, evtTrig4;
      double uuVtx[27];
      double TrkVar[4][12];
      reco::Particle::LorentzVector uuP4fit; //--(0.0,0.0,0.0,0.0);
      double CutMu[3][7];

//--- Structures for ntupling:
struct Ievt{
  int runNb, eventNb, lumiBlock;   // run number, event number and lumi section
  int allTrk, allMu, allPrim, allK0s;
  int upTrig, umTrig, upLMtrig, umLMtrig;
  int evtTrig1, evtTrig2, evtTrig3, evtTrig4;
  int Comb_uu, Comb_kPi, Comb_kKa;
  void init();
} IevtC_;

struct tightMuC{
  double iNormChi2, iMuMuHits, iMuStations, iMuPxHits, iMuTrkLayer;
  double jNormChi2, jMuMuHits, jMuStations, jMuPxHits, jMuTrkLayer;
  void init();
} tightMuC_;
struct softMuC{
  double iMuInTrkLayer,iMuInNormChi2;
  double jMuInTrkLayer,jMuInNormChi2;
  void init();
} softMuC_;
struct tightKMuC{
  double kNormChi2, kMuMuHits, kMuStations, kMuPxHits, kMuTrkLayer;
  void init();
} tightKMuC_;
struct softKMuC{
  double kMuInTrkLayer,kMuInNormChi2;
  void init();
} softKMuC_;

struct IuuC{
  int iQ,imuId,iNtrkHits;
  int jQ,jmuId,jNtrkHits;
  void init();
} IuuC_;
struct uuC{
  double iEta,iPhi,iPx,iPy,iPz,id0,idz,iCal,iSeg,iIso;
  double jEta,jPhi,jPx,jPy,jPz,jd0,jdz,jCal,jSeg,jIso;
  double Muu, dca, xcpt, ycpt, zcpt, uuDist, uuR;
  double MuuVF, uuVCL, uuVMpx, uuVMpy, uuVMpz;
  double MuuMFit, uuVMCL, Vx, Vy, Vz, eta, phi;
  void init();  
} uuC_;

struct kPi{
  double kPiQ,kPiPx,kPiPy,kPiPz,kPiEn,kPiEta,kPiPhi,kPiNhits,kPid0,kPidz;
  void init();
} kPiC_;
struct nPi{
  double nPiQ,nPiPx,nPiPy,nPiPz,nPiEn,nPiEta,nPiPhi,nPiNhits,nPid0,nPidz;
  void init();
} nPiC_;
struct rPi{
  double rPiQ,rPiPx,rPiPy,rPiPz,rPiEn,rPiEta,rPiPhi,rPiNhits,rPid0,rPidz;
  void init();
} rPiC_;

struct kKa{
  double kKaQ,kKaPx,kKaPy,kKaPz,kKaEn,kKaEta,kKaPhi,kKaNhits,kKad0,kKadz;
  void init();
} kKaC_;

struct vtx{
  double MuukPiFit, MuukPi, MuikPi, MujkPi, Ruu_kPi;
  double Muupk, MkPikKa, Ruu_kKa, MuupkFit, uupkV_CL;
  double pV_CL, uupkLoS, uupkL, cos_alpha;
  void init();
} vtxC_;
struct k0s{
  double Mk0s, Mk0sVFit, k0sV_CL, k0sVM_CL, Mlam, k0sDist, k0sR;
  double cos_uupk0s3D, uuk0sL_L3D, uuk0s_Sigma3D;
  double MuupFit, uupV3_CL, pV_CL;
  double MuupK0s, MkPikK0s, uupL, uupSigma, cos_alpha_uup_prim;
  void init();
} k0sC_;

TTree *hecmu_tree_;
};
#endif
