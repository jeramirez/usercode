// -*- C++ -*-
//
// Package:    HecPsiLambdaPAT
// Class:      HecPsiLambdaPAT
// 
/**\class HecPsiLambdaPAT HecPsiLambdaPAT.cc HecBaryons/HecPsiLambda/src/HecPsiLambdaPAT.h

 Description: [one line class summary]
 Do some analysis to Justify my Salary.
 Reconstruct /\_b using muons in the final State.

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Hector Mendez
//         Created:  Mon Nov 21 11:34:49 CST 2011
// $Id: HecPsiLambdaPAT.h,v 1.4 2013/01/17 21:25:26 mendez Exp $
//
//

#ifndef HecPsiLambdaPAT_h
#define HecPsiLambdaPAT_h

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

class HecPsiLambdaPAT : public edm::EDAnalyzer {
   public:
      explicit HecPsiLambdaPAT(const edm::ParameterSet&);
      ~HecPsiLambdaPAT();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      
   protected:
      void fill_evt(const edm::Event&, int, int, int, int, int, int, int,
                    int, int, int, int, int );
      void fill_IuuC(int, int, int, int, int, int );
      void fill_uuC(double,double,double,double,double,double,double,double,double,double,
                    double,double,double,double,double,double,double,double,double,double,
                    double,double,double,double,double,double,double,
                    double,double,double,double,double,
                    double,double,double,double,double,double,
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
      edm::InputTag trackTags_;       //--used to select what tracks to read from configuration file
      edm::InputTag theMuonsLabel_;   //--here I get both global & track Mu [stolen from Heavy Oni2MuMu]
      std::string VeeAlgo_;           //--Label of Lambda collection
      bool MyPrint_;
      bool MyPrintMC_;
      bool doMC_;
      std::string hlTriggerResults_;  //--HLT Trigger  (from Keith)
      
      TH1D *histo_Ntrk, *histo_Nmu, *histo_Nlambda, *histo_Nprim;
      TH1D *histo_dimuM[10], *histo_lam[10], *histo_lamb[10], *histo_diff[2];
      TH1D *histo_cos3d, *histo_cos2d;
      TH1D *histo_closest, *histo_sig3d[4], *histo_sig2d[4];
      TH1D *histo_MCgen, *histo_NLambda0, *histo_L0daug;
            
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
    double Muu, dca, xcpt, ycpt, zcpt, uuDist, uuR;
    double MuuVF, uuVCL, uuVMpx, uuVMpy, uuVMpz;
    double MuuMFit, uuVMCL, cosuul3D, cosuul2D, uuLL3D, uuLSigma3D;
    double closestTRK, closestD0, closestSigmaDz, closestDz, closestTRK3d, closestId;
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
  TTree *hecmu_tree_;
};
#endif
