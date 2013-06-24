// -*- C++ -*-
//
// Package:    HecPsiCascade
// Class:      HecPsiCascade
// 
/**\class HecPsiCascade HecPsiCascade.cc HecBaryons/HecPsiCascade/interface/HecPsiCascade.h

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Hector Mendez
//         Created:  Tue Jun 14 08:50:59 CDT 2011
// $Id: HecPsiCascade.h,v 1.16 2011/10/18 16:26:49 mendez Exp $
//
//
#ifndef HecPsiCascade_h
#define HecPsiCascade_h

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
/*
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
*/
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TFile.h"
#include "TList.h"
#include "TH1D.h"
#include "TH1I.h"
#include "TH2D.h"
#include "TTree.h"
#include "TNtuple.h"

// class declaration  

class HecPsiCascade : public edm::EDAnalyzer {
   public:
      explicit HecPsiCascade(const edm::ParameterSet&);
      ~HecPsiCascade();
      
   protected:
      void fill_evt(const edm::Event&, int, int, int, int, int, int );
      void fill_uuC();
      void fill_uuFit();
      void init();
      
   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------
      const edm::ParameterSet iConfig;
      edm::InputTag trackTags_;          //--used to select what tracks to read from configuration file
      edm::InputTag electronCollection_;
      edm::InputTag theMuonsLabel_;      //--here I get both global & track Mu [stolen from Heavy Oni2MuMu]
      edm::InputTag theCaloMuonsLabel_;
      unsigned int minTracks_;
      bool   check_elC_;
      double Pe_match_;
      std::string mass_constrain_;
      double massC_value;
      float massC_sigma;
      bool onlyDiMu_, MyPrint_;
      std::string CasAlgo_;           //--Label of Cascade collection 
      std::string CasDecayName_;      //--Name of Cascade Collection (Cascade or Omega)
      std::string VeeAlgo_;           //--Label of Lambda collection
    //std::string VeeDecayName_;      //--Name of Cascade Collection (Lambda or Kshort)
                
      TH1D *histo_trk, *histo_Pti, *histo_ij, *histo_ijPt, *histo_nmu, *histo_muId;
      TH1D *histo_dimu, *histo_dimuGG, *histo_dimuGT, *histo_dimuTT, *histo_dimu00, *histo_dimu01;
      TH1D *histo_dimuCan;
      TH1D *histo_LoS;
      TH1D *histo_dimuX0, *histo_dimuX1, *histo_dimuX2, *histo_dimuX3;
      TH1D *histo_cas01, *histo_cas02, *histo_cas03, *histo_cas04;
      TH1D *histo_lam01, *histo_lam02;
      //TH1D * histo_ij0,* histo_ij1,* histo_ij2;
      TH1D * histo_elec;
      //TH1D * histo_NeTrk;
      //TH1D * histo_DNel,* histo_Delp,* histo_DelpC,* histo_Delq;
      TH1D *histo_eiej, *histo_eEn, *histo_EoP;
      //TH1D * histo_Clu,* histo_Eni;
      TH1D *histo_uu_prim_cos;
            
//--- Structures for ntupling:
  struct Ievt{
    int runNb, eventNb, lumiBlock;          // run number, event number and lumi section 
    int Ntk, allMu, allCalMu, allCascade, allLambda,allPrim;  
    void init();
  } Ievt_;
  struct IuuC{
    int iGL,iCharge,imuId,iNtrkHits;
    int jGL,jCharge,jmuId,jNtrkHits;
    void init();
  } IuuC_;
  struct uuC{
    double iEta,iPhi,iPt,iP,id0,idz,iCal,iSeg,iIso;
    double jEta,jPhi,jPt,jP,jd0,jdz,jCal,jSeg,jIso;
    double Muu, dca, xcpt, ycpt, zcpt;
    void init();
  } uuC_;
  struct uuFit{
    double MuuFit, MuuMFit, uuLoS, uuL, uuV_Cl, pV_Cl, uuVM_Cl;
    double uuVx, uuVy, uuVz, xpV, ypV, zpV;
    double uuVM_px, uuVM_py, uuVM_pz;
    double VuuVpCos;
    void init();
  } uuFit_;
  struct ILam{
    int uuPair, xiCan;
    int NPixelXiTracks;
    void init();
  } ILam_;    
  struct Lam{
    double xiM, lamM, uuXiM_0, uuXiM_1, uuXiM_2;
    double MlamVFit, lamV_Cl,  lamVx,    lamVy,    lamVz;
    double MlamMFit, lamVM_Cl, lamVM_px, lamVM_py, lamVM_pz;
    double uuL_L, uuL_LoS;
    void init();
  } Lam_;    
  struct Cas{
    double dca, xcPt, ycPt, zcPt;
    double MxiVFit, xiV_Cl,  xiVx,    xiVy,    xiVz;
    double MxiMFit, xiVM_Cl, xiVM_px, xiVM_py, xiVM_pz;
    double uuXi_L, uuXi_LoS;
    double Xiip3d, VuuVxiCosOp, XiCL;
    double rhoPi, rhoXi, IpPi, IpXiuu;
    void init();
  } Cas_;    
  struct Casb{
    double Mxi_bVFit, xi_bV_Cl, xi_bVx, xi_bVy, xi_bVz;
    double xi_bV_px, xi_bV_py, xi_bV_pz;
    double xi_bL, xi_bLoS, uuxi_bL, uuxi_bLoS;
    double VxibVpCos;
    void init();
  } Casb_;
  TTree *hecmu_tree_;
};
#endif
