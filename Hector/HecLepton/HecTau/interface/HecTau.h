//
// class declaration
//

#ifndef HecTau_h
#define HecTau_h

// user include files
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

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
#include "TFile.h"
#include "TList.h"
#include "TH1D.h"
#include "TH1I.h"
#include "TH2D.h"
#include "TTree.h"
#include "TNtuple.h"

class HecTau : public edm::EDAnalyzer {
   public:
      explicit HecTau(const edm::ParameterSet&);
      ~HecTau();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   protected:
      void fill_evt(const edm::Event&, int, int, int, int);

      void fill_IuuC();

      void fill_uuC(double);

      void fill_IuukC(int, int, int, int );

      void fill_uukC(double,double,double,double,double,double,
                     double,double,double,double,double,double,
                     double,double,double,double,double,double,double,double,double,
                     double,double,double,double,double,double,double,double,
                     double,double,double,double,double,double );
		     
      void fill_tightMuons();
			   
      void fill_softMuons();
			  
      void fill_tightKMuons(double,double,double,double,double);
      
      void fill_softKMuons(double,double);
			  
      void init();

      void HecMuVar(int,std::vector<pat::Muon>::const_iterator, reco::TrackRef);
      bool HecClosest(int, FreeTrajectoryState, FreeTrajectoryState);
      void HecHltTrig(const edm::Event&);
      void HecHltMuTrig(int,const pat::Muon*);   
      void HecCutMu(int,std::vector<pat::Muon>::const_iterator);

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
      edm::InputTag trackTags_;          //used to select what tracks to read from configuration file
      edm::InputTag theMuonsLabel_;      //--here I get both global & track Mu [stolen from Heavy Oni2MuMu]
      bool MyPrint_, MyOneKaon_, MyThreeMu_;
      std::string hlTriggerResults_;

      TH1D * histo, *histo_trk, *histo_nmu, *histo_muId;
      TH1D *histo_dimu, *histo_dimuk[5];
      TH1D *histo_DeltaPku, *histo_DeltaQ;
      TH1D *histo_uuu[3];
      TH1D *histo_pV_ClHM;
      int iComb;

      //--> Variables for my functions <--
      double MuVar[2][13];
      int   nMuVar[2][5];
      double dcaVar[6][4];  //--see function HecClosest() Implementation
      int    evtTrig1,evtTrig2,evtTrig3,evtTrig4;
      int    nTrig[3][2];
      double CutMu[2][7];

//--- Structures for ntupling:
  struct Ievt{
    int runNb, eventNb, lumiBlock;// run number, event number and lumi section
    int Ntk, allMu, allPrim,iComb;
    int upTrig,umTrig,umLMtrig,upLMtrig,ukTrig,mukLMtrig;
    int evtTrig1,evtTrig2,evtTrig3,evtTrig4;
    void init();
  } Ievt_;
  struct IuuC{
    int iQ,imuId,iNtrkHits;
    int jQ,jmuId,jNtrkHits;
    void init();
  } IuuC_;
  struct uuC{
    double iEta,iPhi,iPt,iPx,iPy,iPz,id0,idz,iCal,iSeg,iIso,iKink;
    double jEta,jPhi,jPt,jPx,jPy,jPz,jd0,jdz,jCal,jSeg,jIso,jKink;
    double Muu, dca, xcpt, ycpt, zcpt;
    void init();
  } uuC_;
  struct IuukC{
    int kGL,kCharge,kmuId,kNtrkHits;
    void init();
  } IuukC_;
  struct uukC{
    double kEta,kPhi,kPt,kPx,kPy,kPz,kd0,kdz,kCal,kSeg,kIso,kKink;
    double Muuk,Muik,Mujk,MuukFit,uukV_Cl,pV_Cl,primDist,primR,primDz,primDxy,primDxyz;
    double uukL,uukLoS,cos_alpha;
    double iVpt,jVpt,kVpt;
    double uukVx,uukVy,uukVz,PrimVx,PrimVy,PrimVz;
    void init();
  } uukC_;
  struct tightMuC{
    double iNormChi2, iMuMuHits, iMuStations, iMuPxHits, iMuTrkLayer;
    double jNormChi2, jMuMuHits, jMuStations, jMuPxHits, jMuTrkLayer;
    void init();
  } tightMuC_;
  struct softMuC{
    double iMuInNormChi2,iMuInTrkLayer;
    double jMuInNormChi2,jMuInTrkLayer;
    void init();
  } softMuC_;
  struct tightKMuC{
    double kNormChi2, kMuMuHits, kMuStations, kMuPxHits, kMuTrkLayer;
    void init();
  } tightKMuC_;
  struct softKMuC{
    double kMuInNormChi2, kMuInTrkLayer;
    void init();
  } softKMuC_;
  TTree *hecmu_tree_;
            
};
#endif
