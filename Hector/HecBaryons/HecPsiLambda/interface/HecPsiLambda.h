// -*- C++ -*-
//
// Package:    HecPsiLambda
// Class:      HecPsiLambda
// 
/**\class HecPsiLambda HecPsiLambda.cc HecBaryons/HecPsiLambda/src/HecPsiLambda.h

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Hector Mendez
//         Created:  Mon Nov 21 11:34:49 CST 2011
// $Id: HecPsiLambda.h,v 1.5 2012/06/04 19:44:18 mendez Exp $
//
//

#ifndef HecPsiLambda_h
#define HecPsiLambda_h

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

class HecPsiLambda : public edm::EDAnalyzer {
   public:
      explicit HecPsiLambda(const edm::ParameterSet&);
      ~HecPsiLambda();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      
   protected:
      void fill_evt(const edm::Event&, int, int, int, int, int );
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
                        double, double, double );                    
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
      std::string hlTriggerResults_;  //--HLT Trigger  (from Keith)
      
      TH1D *histo_Ntrk, *histo_Nmu, *histo_Nlambda, *histo_Nprim;
      TH1D *histo_dimuM[10], *histo_lam[10], *histo_lamb[10], *histo_diff[2];
      TH1D *histo_cos3d, *histo_cos2d;
      TH1D *histo_closest, *histo_sig3d[4], *histo_sig2d[4];
       
      //--- Structures for ntupling:
  struct Ievt{
    int runNb, eventNb, lumiBlock;          // run number, event number and lumi section 
    int allTrk, allMu, allPrim, allLambda;
    int iComb; 
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
    double L_2D, LoS_2D, MuuK0s;
    void init();
  } lambdabC_;
  TTree *hecmu_tree_;
};
#endif
