// -*- C++ -*-
//
// Package:    Hec2sLambda
// Class:      Hec2sLambda
// 
/**\class Hec2sLambda Hec2sLambda.cc HecBaryons/Hec2sLambda/src/Hec2sLambda.h

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Hector Mendez
//         Created:  Sun Jan  8 11:56:28 CST 2012
// $Id:
//
//

#ifndef Hec2sLambda_h
#define Hec2sLambda_h

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

class Hec2sLambda : public edm::EDAnalyzer {
   public:
      explicit Hec2sLambda(const edm::ParameterSet&);
      ~Hec2sLambda();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      
   protected:
      void fill_evt(const edm::Event&, int, int, int, int, int );
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
                        double, double, double, double ); 
                        
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
      TH1D *histo_Ntrk, *histo_Nmu, *histo_Nlambda, *histo_Nprim;
      TH1D *histo_dimuM[10], *histo_lam[10], *histo_lamb[10], *histo_diff[2];
      TH1D *histo_mdiff[6], *histo_mpipi[2], *histo_vx[3];
      TH1D *histo_cos3d, *histo_cos2d, *histo_DeltaR, *histo_iPiPhi;
       
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
    void init();
  } lambdabC_;
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
