// -*- C++ -*-
//
// Package:    Hec2sLambdaPAT
// Class:      Hec2sLambdaPAT
// 
/**\class Hec2sLambdaPAT Hec2sLambdaPAT.cc HecBaryons/Hec2sLambda/src/Hec2sLambdaPAT.h

 Description: [one line class summary]
 Do some analysis to Justify my Salary.
 Reconstruct /\_b using muons and pi+pi- in the final State.

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Hector Mendez
//         Created:  Sun Jan  8 11:56:28 CST 2012
// $Id:
//
//

#ifndef Hec2sLambdaPAT_h
#define Hec2sLambdaPAT_h

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

class Hec2sLambdaPAT : public edm::EDAnalyzer {
   public:
      explicit Hec2sLambdaPAT(const edm::ParameterSet&);
      ~Hec2sLambdaPAT();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      
   protected:
      void fill_evt(const edm::Event&, int, int, int, int, int, int, int,
                    int, int, int, int, int );
      void fill_IuuC(int, int, int, int, int, int );
      void fill_uuC(double,double,double,double,double,double,double,double,double,double,
                    double,double,double,double,double,double,double,double,double,double,
                    double,double,double,double,double,double,double,double,
                    double,double,double,double,double,
                    double,double,double,double,double,double);
      void fill_lambda0(double,double,double,double,double,double,double,double,double,     //--proton
                        double,double,double,double,double,double,double,double,double,     //--pion
                        double,double,double,double,double,double,
                        double,double,double,double,double);    //--Lambda Fit     

      void fill_lambdab(double, double, double, double, double, double, double,
                        double, double, double, double,
                        double, double, double, double, double ); 
      void fill_tightMuons(double, double, double, double, double, double,
                           double, double, double, double, double, double);        
                        
      void fill_pipi(double,double,double,double,double,double,double,double,double,
                     double,double,double,double,double,double,double,double,double,
                     double,double,double,double,
                     double,double,double,
                     double,double);    
      void init();
      
      double SignificanceAbsoluteImpactParameter3D(reco::TransientTrack &, reco::Vertex &,
                           double &, double &,
                           double &, double &, double &)const;

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
      edm::InputTag trackTags_;     //--used to select what tracks to read from configuration file
      edm::InputTag theMuonsLabel_; //--here I get both global & track Mu [stolen from Heavy Oni2MuMu]
      std::string VeeAlgo_;         //--Label of Lambda collection
      bool MyPrint_;
      bool MyPrintMC_;
      bool doMC_;
      std::string hlTriggerResults_;  //--HLT Trigger  (from Keith)

      TH1D *histo_Ntrk, *histo_Nmu, *histo_Nlambda, *histo_Nprim;
      TH1D *histo_dimuM[10], *histo_lam[10], *histo_lamb[10], *histo_diff[2];
      TH1D *histo_mdiff[6], *histo_mpipi[2], *histo_vx[3];
      TH1D *histo_cos3d, *histo_cos2d, *histo_DeltaR, *histo_iPiPhi;
      TH1D *histo_MCgen, *histo_NLambda0, *histo_L0daug;
//      unsigned int hlt_mu3, hlt_mu5, hlt_mu8;
//      unsigned int hlt_2mu0,hlt_2mu3,hlt_2mu3JPsi_v1,hlt_2mu3JPsi_v2,hlt_2mu3_quark_v1,hlt_2mu3_quark_v2; 
//      unsigned int hlt_2mu6p5_dis,hlt_2mu7_dis,hlt_2mu3p5_dis,hlt_2mu4_dis;
//      unsigned int hlt_2mu6p5_pro, hlt_2mu0_pro_v1, hlt_2mu0_pro_v5, hlt_2mu0_pro_v6;
//      unsigned int hlt_2mu0_pro_NoVtx_v2, hlt_2mu0_pro_NoVtx_v3;
//      unsigned int hlt_mu0trk0, hlt_mu3trk0, hlt_mu0trkmu0, hlt_mu3trkmu0, hlt_mu0trkmu0OST, hlt_mu3trkmu0OST;
//      unsigned int hlt_mu0trkmu0OST_tight;
//      unsigned int hlt_L1muOpen, hlt_L12muOpen, hlt_L12muOpenTight, hlt_2mu0L2, hlt_DimuPsi2s, hlt_DimuLowMd;
      
   //--std::vector<int> *mupTrig2mu3;
   //--std::vector<int> *mumTrig2mu3;
//      int upTrig2mu3,upTrig2mu3JPsi,upTrig2mu6p5JPsiDisp,upTrig2mu7JPsiDisp;
//      int umTrig2mu3,umTrig2mu3JPsi,umTrig2mu6p5JPsiDisp,umTrig2mu7JPsiDisp;
      
//      int upTrig2mu4JPsiDisp,upTrig2mu6p5JPsiPro,upTrig2mu0JPsiPro,upTrig2mu0JPsiNoVtxPro;    
//      int umTrig2mu4JPsiDisp,umTrig2mu6p5JPsiPro,umTrig2mu0JPsiPro,umTrig2mu0JPsiNoVtxPro;
//      int evtTrig1, evtTrig2, evtTrig3, evtTrig4;
            
      //--- Structures for ntupling:
  struct Ievt{
    int runNb, eventNb, lumiBlock;          // run number, event number and lumi section 
    int allTrk, allMu, allPrim, allLambda;
    int iComb, upTrig, umTrig, muLMtrig;
    int evtTrig1, evtTrig2, evtTrig3, evtTrig4;
    void init();
  } Ievt_;
  struct IuuC{
    int iQ,imuId,iNtrkHits;
    int jQ,jmuId,jNtrkHits;
    void init();
  } IuuC_;
  struct uuC{
    double iEta,iPhi,iPx,iPy,iPz,id0,idz,iCal,iSeg,iIso;
    double jEta,jPhi,jPx,jPy,jPz,jd0,jdz,jCal,jSeg,jIso;
    double Muu, dca, xcpt, ycpt, zcpt, uuDist, uuR, uuDR;
    double MuuVF, uuVCL, uuVMpx, uuVMpy, uuVMpz;
    double MuuMFit, uuVMCL, cosuul3D, cosuul2D, uuLL3D, uuLSigma3D;
    void init();
  } uuC_;
  struct lambdaC{
    double prQ, prEta, prPhi, prPx, prPy, prPz, prd0, prdz, prNHits;
    double piQ, piEta, piPhi, piPx, piPy, piPz, pid0, pidz, piNHits;
    double Mlambda0;
    double lamCL, Mlambda0VF, lamDist, lamR, Mkshort;
    double MlamMFit, lamVMCL, lamVMpx, lamVMpy, lamVMpz;
    void init();
  } lambdaC_;
  struct lambdabC{
    double lambCL, MlambVF, lambDist, lambR;
    double primCL, primDist, primR, lambdabL3D, lambdabSigma3D, cosalphab3D, cosalphab2D;
    double L_2D, LoS_2D, MuuK0s, primDxy, primDxyz;
    void init();
  } lambdabC_;
  struct tightMuC{
    double iNormChi2, iMuMuHits, iMuStations, iMuPxHits, iMuTrkLayer, iKink;
    double jNormChi2, jMuMuHits, jMuStations, jMuPxHits, jMuTrkLayer, jKink;
    void init();
  } tightMuC_; 
  struct uupipiC{
    double iPiQ,iPiEta,iPiPhi,iPiPx,iPiPy,iPiPz,iPiNtrkHits,iPid0,iPidz;
    double jPiQ,jPiEta,jPiPhi,jPiPx,jPiPy,jPiPz,jPiNtrkHits,jPid0,jPidz;
    double MuupipiVFit,uupipiV_CL,MuupipiMFit,uupipiVM_CL;
    double uupipiVPx,uupipiVPy,uupipiVPz;
    double ipiuuR,jpiuuR;
    void init(); 
  } uupipiC_;  
  TTree *hecmu_tree_;
};
#endif
