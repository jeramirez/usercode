// -*- C++ -*-
//
// Package:    HecPsiLambdaPAT
// Class:      HecPsiLambdaPAT
// 
/**\class HecPsiLambdaPAT HecPsiLambdaPAT.cc HecBaryons/HecPsiLambda/src/HecPsiLambdaPAT.cc

 Description: [one line class summary]
 Do some analysis to Justify my Salary.
 Reconstruct /\_b using muons in the final State.
    
 Implementation:
     [Notes on implementation]
  July 7,2012  Added Maria's MC counting
  Need to be done [Jan 4, 2012]
  - cowboys dimuon
  - L/S in 2D
  - Check primary for pile-up
  - Separate by eta
*/
//
// Original Author:  Hector Mendez
//         Created:  Mon Nov 21 11:34:49 CST 2011
// $Id: HecPsiLambdaPAT.cc,v 1.8 2013/01/18 22:19:07 mendez Exp $
//
//


// system include files
#include <memory>

// user include files
#include "HecBaryons/HecPsiLambda/interface/HecPsiLambdaPAT.h"    //--Hec
int Lb[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};  //--Total number of Lambda_b decays modes
int Fb[] = {0, 0, 0, 0};                    //--Total number of all Lambda_b generated

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonTrackLinks.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/CaloMuon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"       //cascade colletcion
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"       //cascade daughter track

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"          //track   collection
//#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
//#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"  
//#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h" 

#include "Analyzers/CascadeProducer/interface/PrimaryInfo.h"             //--Eduardo Cascade & Primary
#include "Analyzers/CascadeProducer/interface/VertexRefit.h"
#include "Analyzers/CascadeProducer/interface/KinematicFitDriver.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "Analyzers/CascadeProducer/interface/ClosestApproachOnHelixLine.h"

//--Bfield
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"

//--kinemactic fitter
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
//--Trigger (from keith)
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"

#include "Analyzers/CascadeProducer/interface/masses.h"

#include "TH1.h"

#include <vector>
#include <utility>
//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
HecPsiLambdaPAT::HecPsiLambdaPAT(const edm::ParameterSet& iConfig)
:
 trackTags_        (iConfig.getParameter<edm::InputTag>("tracks"    )),
 theMuonsLabel_    (iConfig.getParameter<edm::InputTag>("MuonsLabel")),
 VeeAlgo_          (iConfig.getUntrackedParameter<std::string>("VeeAlgo","generalV0Candidates")),
 MyPrint_          (iConfig.getUntrackedParameter<bool>("MyPrint","false")),
 MyPrintMC_        (iConfig.getUntrackedParameter<bool>("MyPrintMC", "false")),
 doMC_             (iConfig.getUntrackedParameter<bool>("doMC", "false") ),
 hlTriggerResults_ (iConfig.getUntrackedParameter<std::string>("HLTriggerResults",std::string("TriggerResults::HLT")) )
{
  //now do what ever initialization is needed
  edm::Service<TFileService> fs;
  histo_Ntrk     = fs->make<TH1D>("Ntrk"   ,"# Tracks", 200, 0, 200 );
  histo_Nmu      = fs->make<TH1D>("NMuons" ,"# Muons" ,  20, 0,  20 );
  histo_Nlambda  = fs->make<TH1D>("NLambda","# Lambda",  20, 0,  20 );
  histo_Nprim    = fs->make<TH1D>("NPrim"  ,"# Prim"  , 100, 0, 100 );
  
  histo_lam[0]   = fs->make<TH1D>("MassLam0","Mass Lam All" , 200, 1.07, 1.17); //--0.5 MeV/channel
  histo_lam[1]   = fs->make<TH1D>("MassLam1","Mass Vtx"     , 200, 1.07, 1.17);
  histo_lam[2]   = fs->make<TH1D>("MassLam2","Mass VtxMass" , 200, 1.07, 1.17);
  histo_lam[3]   = fs->make<TH1D>("MassLam3","Mass VtxMcons", 200, 1.07, 1.17);
  histo_lam[4]   = fs->make<TH1D>("MassLam4","Mass AllDimu" , 200, 1.07, 1.17);
  histo_lam[5]   = fs->make<TH1D>("MassLam5","repTrk Op"    , 200, 1.07, 1.17);
  histo_lam[6]   = fs->make<TH1D>("MassLam6","repTrk Same"  , 200, 1.07, 1.17);
  histo_lam[7]   = fs->make<TH1D>("MassLam7","Mass uuVtx"   , 200, 1.07, 1.17);
  histo_lam[8]   = fs->make<TH1D>("MassLam8","Mass JLVtx"   , 200, 1.07, 1.17);
  histo_lam[9]   = fs->make<TH1D>("MassLam9","Mass Fin"     , 200, 1.07, 1.17);

//histo_dimuM[0] = fs->make<TH1D>("dimuM0","Muu All"      , 500, 0, 100);   //--200 MeV/channel
//histo_dimuM[1] = fs->make<TH1D>("dimuM1","Muu All"      , 500, 0, 100);
//histo_dimuM[2] = fs->make<TH1D>("dimuM2","Muu All"      , 500, 0, 100);
//histo_dimuM[3] = fs->make<TH1D>("dimuM3","Muu All"      , 500, 0, 100);
  histo_dimuM[4] = fs->make<TH1D>("dimuM4","Muu All"      , 500, 0, 100);
  histo_dimuM[5] = fs->make<TH1D>("dimuM5","repTrk Op"    , 500, 0, 100);
  histo_dimuM[6] = fs->make<TH1D>("dimuM6","repTrk Same"  , 500, 0, 100);
  histo_dimuM[7] = fs->make<TH1D>("dimuM7","Muu uuVtx"    , 500, 0, 100);
  histo_dimuM[8] = fs->make<TH1D>("dimuM8","Muu JLVtx"    , 500, 0, 100);
  histo_dimuM[9] = fs->make<TH1D>("dimuM9","Muu Fin"      , 500, 0, 100);
  
  histo_lamb[8] = fs->make<TH1D>("lamb8","MuuL All" , 200, 0, 20);   //--100 MeV/channel
  histo_lamb[9] = fs->make<TH1D>("lamb9","MuuL Fin" , 200, 0, 20);
  
  float r = 1.0/10000.0;
  histo_diff[0]  = fs->make<TH1D>("Diff0","M recLambda - fitLambda", 100, -r, r);
  histo_diff[1]  = fs->make<TH1D>("Diff1","Muu - MuuFit"           , 100, -r, r);
  
  histo_cos3d = fs->make<TH1D>("cos3d","cos3d uuL", 220,-1.1,1.1);
  histo_cos2d = fs->make<TH1D>("cos2d","cos2d uuL", 220,-1.1,1.1);
  
  histo_sig3d[0] = fs->make<TH1D>("sig3d0","Muon Significance 3d", 100,-5,5);
  histo_sig2d[0] = fs->make<TH1D>("sig2d0","Muon Significance 2d", 100,-5,5);
  histo_sig3d[1] = fs->make<TH1D>("sig3d1","Lam0 Significance 3d", 100,-5,5);
  histo_sig2d[1] = fs->make<TH1D>("sig2d1","Lam0 Significance 2d", 100,-5,5);
  histo_sig3d[2] = fs->make<TH1D>("sig3d2","allT Significance 3d", 100,-5,5);
  histo_sig2d[2] = fs->make<TH1D>("sig2d2","allT Significance 2d", 100,-5,5);
  histo_sig3d[3] = fs->make<TH1D>("sig3d3","Clos Significance 3d", 100,-5,5);
  histo_sig2d[3] = fs->make<TH1D>("sig2d3","Clos Significance 2d", 100,-5,5);
  
  histo_closest  = fs->make<TH1D>("closest","Closest to uu", 300,-3,3);
  histo_MCgen    = fs->make<TH1D>("MCgen","MC generated event", 15,0,15);
  
  histo_NLambda0= fs->make<TH1D>("NLambda0","MC gen Lambda0 daug", 10,0,10);
  histo_L0daug  = fs->make<TH1D>("L0daug","MC gen L0daug id", 4600,-2300,2300);
}
HecPsiLambdaPAT::~HecPsiLambdaPAT()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}
//
// member functions
//

// ------------ method called for each event  ------------
void
HecPsiLambdaPAT::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   //--Tracks Collection
 //edm::Handle<reco::TrackCollection> tracks;
   edm::Handle< std::vector<pat::GenericParticle> >tracks;    //--thePATTrackHandle;
   iEvent.getByLabel(trackTags_,tracks); int allTrk = tracks->size();
   
//--Load Tracks into a vector of Tracks (from Eduardo CascadeFitter.cc)
   //std::vector<reco::Track> theTracks;
   //theTracks.insert( theTracks.end(), 
   //                  tracks->begin(), 
   //                  tracks->end()  );
                        
   //--Muon Collection
 //edm::Handle<reco::MuonCollection> allmuons;
   edm::Handle< std::vector<pat::Muon> >allmuons;             //--thePATMuonHandle;
   iEvent.getByLabel(theMuonsLabel_, allmuons); int allMu = allmuons->size();
   
   //--Handles for Std Vees (Lambdas) and Load them into a vector of composite candidates  
   std::vector<reco::VertexCompositeCandidate> theVees;
   edm::Handle<reco::VertexCompositeCandidateCollection> theVeeHandle;
   iEvent.getByLabel(VeeAlgo_, "Lambda", theVeeHandle);

   theVees.insert( theVees.end(), theVeeHandle->begin(), theVeeHandle->end() );
   int allLambda  = theVees.size();
   
   //--Handles for Primary Vertex  (from Eduardo)
   PrimaryInfo primary(iEvent,iConfig);
  
   //--Handle fron Primary just to get the number of rec Prim vertex (Sep 10, 11)
   edm::Handle<reco::VertexCollection> privtxs;
   //iEvent.getByLabel(thePrimaryVertexLabel, privtxs); int allPrim = privtxs->size();
   //iEvent.getByLabel("PrimaryCollection", privtxs); int allPrim = privtxs->size();
   iEvent.getByLabel("offlinePrimaryVertices", privtxs); int allPrim = privtxs->size();
      
   //--Handles for beamspot used for primary refit
   reco::BeamSpot beamSpot;
   edm::Handle<reco::BeamSpot> beamSpotHandle;
   iEvent.getByLabel("offlineBeamSpot", beamSpotHandle);
   if( beamSpotHandle.isValid() )beamSpot = *beamSpotHandle;
     
   //--Handles for B-field
   edm::ESHandle<MagneticField> bFieldHandle;
   iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);
   //--const MagneticField *BField = bFieldHandle.product();
   
   edm::ESHandle<GlobalTrackingGeometry> globTkGeomHandle;
   iSetup.get<GlobalTrackingGeometryRecord>().get(globTkGeomHandle);
   
   //--Loop over tracks to convert to transient tracks (from Eduardo CascadeFitter.cc)
   std::vector<reco::TrackRef> theTrackRefs;
   std::vector<reco::TransientTrack> theTransTracks;
   //for(unsigned int itrack = 0; itrack < tracks->size(); itrack++) {
         //reco::TrackRef tmpRef( tracks, itrack );
//lets loop through the pat tracks to find the match for this reco::track
       for (  std::vector<pat::GenericParticle>::const_iterator itrack  = tracks->begin();
	                                                        itrack != tracks->end();   ++itrack ) {         
         
         reco::TrackRef tmpRef = itrack->track();
         reco::TransientTrack tmpTk2( *tmpRef, &(*bFieldHandle), globTkGeomHandle );
         theTrackRefs.push_back( tmpRef );
         theTransTracks.push_back( tmpTk2 );
   }   
   
   histo_Ntrk->Fill( allTrk );
   histo_Nmu->Fill( allMu ); 
   histo_Nlambda->Fill( allLambda ); 
   histo_Nprim->Fill( allPrim );
   
   if( MyPrint_ )
       std::cout<<"Hec: Ntrk ="<<allTrk
                  <<" Nmuons ="<<allMu
                  <<" NLambda ="<<allLambda
                  <<" NPrim ="<<allPrim<<std::endl;
   init();     //--clean all structures
   int evtTrig1=0, evtTrig2=0, evtTrig3=0, evtTrig4=0;

//--Do trigger here (code from Keith)   
   unsigned int hlt_mu3=0, hlt_mu5=0, hlt_mu8=0;
   unsigned int hlt_2mu0=0,hlt_2mu3=0,hlt_2mu3JPsi_v1=0,hlt_2mu3JPsi_v2=0,hlt_2mu3_quark_v1=0,hlt_2mu3_quark_v2=0; 
   unsigned int hlt_2mu6p5_dis=0,hlt_2mu7_dis=0,hlt_2mu3p5_dis=0,hlt_2mu4_dis=0;
   unsigned int hlt_2mu6p5_pro=0, hlt_2mu0_pro_v1=0, hlt_2mu0_pro_v5=0, hlt_2mu0_pro_v6=0;
   unsigned int hlt_2mu0_pro_NoVtx_v2=0, hlt_2mu0_pro_NoVtx_v3=0;
   unsigned int hlt_mu0trk0=0, hlt_mu3trk0=0, hlt_mu0trkmu0=0, hlt_mu3trkmu0=0, hlt_mu0trkmu0OST=0, hlt_mu3trkmu0OST=0;
   unsigned int hlt_mu0trkmu0OST_tight=0;
   unsigned int hlt_L1muOpen=0, hlt_L12muOpen=0, hlt_L12muOpenTight=0, hlt_2mu0L2=0, hlt_DimuPsi2s=0, hlt_DimuLowMd=0;
     
   //--first get HLT results
   edm::Handle<edm::TriggerResults> hltresults;
   try {
     std::string const &trig = std::string("TriggerResults::")+hlTriggerResults_;
     iEvent.getByLabel(edm::InputTag(trig),hltresults);
   }
   catch ( ... ) {
     std::cout << "Couldn't get handle on HLT Trigger!" << std::endl;
   }
   if (!hltresults.isValid()) {
     std::cout << "No Trigger Results!" << std::endl;
   }
   else {
     int ntrigs = hltresults->size();
     if (ntrigs == 0){
       std::cout << "No trigger name given in TriggerResults of the input " << std::endl;
   }       
   //--get hold of trigger names - based on TriggerResults object!
   const edm::TriggerNames & triggerNames_ = iEvent.triggerNames(*hltresults);
   
   //HLTConfigProvider hltConfig_;
     
   for(int itrig=0; itrig< ntrigs; itrig++){
      TString trigName = triggerNames_.triggerName(itrig);
      int hltflag = (*hltresults)[itrig].accept();
      
      int trigPrescale = 0;
      //if( hltflag == 1 ){
          //std::string trigName_cp = triggerNames_.triggerName(itrig);
          //trigPrescale = hltConfig_.prescaleValue(itrig, trigName_cp);
      //}
          
      if( MyPrint_ ) 
          std::cout << "Trigger " <<  trigName << " was passed = " <<  hltflag <<" prs_Scale:"<<trigPrescale<<std::endl;
          
      if(trigName=="HLT_Mu3_v3") hlt_mu3 = hltflag;
      if(trigName=="HLT_Mu5_v5") hlt_mu5 = hltflag;
      if(trigName=="HLT_Mu8_v1") hlt_mu8 = hltflag;
      if(trigName=="HLT_Dimuon0_Jpsi_v1"        ) hlt_2mu0 = hltflag;
      if(trigName=="HLT_DoubleMu3_v2"           ) hlt_2mu3 = hltflag;
      if(trigName=="HLT_DoubleMu3_Jpsi_v1"      ) hlt_2mu3JPsi_v1   = hltflag;   // 160329-161176 data
      if(trigName=="HLT_DoubleMu3_Jpsi_v2"      ) hlt_2mu3JPsi_v2   = hltflag;   //  161216-163261 data   also in Brian's MC
      if(trigName=="HLT_DoubleMu3_Quarkonium_v1") hlt_2mu3_quark_v1 = hltflag;   // 160329-161176 data 
      if(trigName=="HLT_DoubleMu3_Quarkonium_v2") hlt_2mu3_quark_v2 = hltflag;   //161216-163261 data    also in Brian's MC
  
      if(trigName=="HLT_Dimuon6p5_Jpsi_Displaced_v1") hlt_2mu6p5_dis = hltflag;  // 163269-163869 data
      if(trigName=="HLT_Dimuon7_Jpsi_Displaced_v1"||
         trigName=="HLT_Dimuon7_Jpsi_Displaced_v3"  ) hlt_2mu7_dis   = hltflag;  // 165088-170248 data
 
      if(trigName=="HLT_DoubleMu3p5_Jpsi_Displaced_v2") hlt_2mu3p5_dis = hltflag; // 170249-173235 data
      if(trigName=="HLT_DoubleMu4_Jpsi_Displaced_v1"  ) hlt_2mu4_dis   = hltflag; // 173236- data
  
      if(trigName=="HLT_Dimuon6p5_Jpsi_v1") hlt_2mu6p5_pro  = hltflag;  // 163269-163869 data  prescaled
      if(trigName=="HLT_Dimuon0_Jpsi_v1"  ) hlt_2mu0_pro_v1 = hltflag;  // 165121-170248 data  prescaled
      if(trigName=="HLT_Dimuon0_Jpsi_v5"  ) hlt_2mu0_pro_v5 = hltflag;  // 170249-173235 data  prescaled
      if(trigName=="HLT_Dimuon0_Jpsi_v6"  ) hlt_2mu0_pro_v6 = hltflag;  // 173263- data  prescaled
 
      if(trigName=="HLT_Dimuon0_Jpsi_NoVertex_v2") hlt_2mu0_pro_NoVtx_v2 = hltflag;  // 170249-173235 data  prescaled
      if(trigName=="HLT_Dimuon0_Jpsi_NoVertex_v3") hlt_2mu0_pro_NoVtx_v3 = hltflag;  // 173263- data  prescaled
      
      if(trigName=="HLT_Mu0_Track0_Jpsi"   ) hlt_mu0trk0      = hltflag;
      if(trigName=="HLT_Mu3_Track0_Jpsi"   ) hlt_mu3trk0      = hltflag;
      if(trigName=="HLT_Mu0_TkMu0_Jpsi"    ) hlt_mu0trkmu0    = hltflag;
      if(trigName=="HLT_Mu3_TkMu0_Jpsi"    ) hlt_mu3trkmu0    = hltflag;
      if(trigName=="HLT_Mu0_TkMu0_OST_Jpsi") hlt_mu0trkmu0OST = hltflag;
      if(trigName=="HLT_Mu3_TkMu0_OST_Jpsi") hlt_mu3trkmu0OST = hltflag;
      
      if(trigName=="HLT_Mu0_TkMu0_OST_Jpsi_Tight_v1") hlt_mu0trkmu0OST_tight = hltflag;
      
      if(trigName=="HLT_L1MuOpen"            ) hlt_L1muOpen       = hltflag;
      if(trigName=="HLT_L1DoubleMuOpen"      ) hlt_L12muOpen      = hltflag;
      if(trigName=="HLT_L1DoubleMuOpen_Tight") hlt_L12muOpenTight = hltflag;
      if(trigName=="HLT_L2DoubleMu0_v4"      ) hlt_2mu0L2         = hltflag;
      if(trigName=="HLT_Dimuon7_PsiPrime_v3"         )hlt_DimuPsi2s  = hltflag;
      if(trigName=="HLT_Dimuon7_LowMass_Displaced_v4")hlt_DimuLowMd  = hltflag;
   } //--for (int itrig=0; itrig< ntrigs; itrig++) {   
   } //--end else for if( !hltresults.isValid() ){
   //--pack event Trigger
   evtTrig1 = hlt_mu3           *1
            + hlt_mu5           *10
            + hlt_mu8           *100
            + hlt_2mu0          *1000
            + hlt_2mu3          *10000
            + hlt_2mu3JPsi_v1   *100000
            + hlt_2mu3JPsi_v2   *1000000
            + hlt_2mu3_quark_v1 *10000000;

   evtTrig2 = hlt_2mu3_quark_v2 *1
            + hlt_2mu6p5_dis    *10
            + hlt_2mu7_dis      *100
            + hlt_2mu3p5_dis    *1000
            + hlt_2mu4_dis      *10000
            + hlt_2mu6p5_pro    *100000
            + hlt_2mu0_pro_v1   *1000000
            + hlt_2mu0_pro_v5   *10000000;

   evtTrig3 = hlt_2mu0_pro_v6       *1
            + hlt_2mu0_pro_NoVtx_v2 *10
            + hlt_2mu0_pro_NoVtx_v3 *100
            + hlt_mu0trk0           *1000
            + hlt_mu3trk0           *10000
            + hlt_mu0trkmu0         *100000
            + hlt_mu3trkmu0         *1000000
            + hlt_mu0trkmu0OST      *10000000;

   evtTrig4 = hlt_mu3trkmu0OST       *1
            + hlt_mu0trkmu0OST_tight *10           
            + hlt_L1muOpen           *100
            + hlt_L12muOpen          *1000
            + hlt_L12muOpenTight     *10000
            + hlt_2mu0L2             *100000
            + hlt_DimuPsi2s          *1000000
            + hlt_DimuLowMd          *10000000;                                     

   if( MyPrint_ ) 
       std::cout << hlt_mu3<<" 1: "<< evtTrig1<<" 2: "<< evtTrig2<< " 3: "<< evtTrig3 <<" 4:"<< evtTrig4<<std::endl;

//------------------------------------MC matching June 27-2012--------------------------------------------
   bool found_JPsi    = false;  //uu
   bool found_Psip1   = false;  //uu
   bool found_Psip2   = false;  //j/psi pi pi
   bool found_LambP_0 = false;  //lambda_0 particle
   bool found_LambA_0 = false;  //lambda_0 antiparticle
   
   if( MyPrintMC_ && doMC_ )
           std::cout << "Calling MC matching " << std::endl;
   if( doMC_ ){
       int Tb[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};    //--All Lambda_b decays mode
       int Lbx[] = {0, 0, 0, 0};                     //--All Lambda_b + X
       std::string genParticles_="genParticlesPlusSim";  // this need edm::HepMCProduct
       //std::string genParticles_="genParticles"; //[Jan 2, 2013]
    
       //--get genParticles  
       edm::Handle<reco::GenParticleCollection> genParticles;
      
       iEvent.getByLabel(genParticles_, genParticles);
       if( MyPrintMC_ )
           std::cout << "Size of genParticle collection is " << genParticles->size() << std::endl;
      
       //--check if any of our signals were generated
       for(size_t k = 0; k < genParticles->size(); k++ ){      
           const reco::Candidate & BCand = (*genParticles)[ k ];            //--Lambda b candidate
                      
           //-------------Lambda b candidate------------------------------------
           if ( BCand.pdgId()== 5122 )Fb[0]++,Lbx[0]++;  //--Count all      /\_b
           if ( BCand.pdgId()==-5122 )Fb[1]++,Lbx[1]++;  //--Count all anti-/\_b 
           if ( abs(BCand.pdgId())==5122 && abs(BCand.daughter(0)->pdgId())!=5122 ){
               if ( BCand.pdgId()== 5122 )Fb[2]++,Lbx[2]++;  //--Count all      /\_b
               if ( BCand.pdgId()==-5122 )Fb[3]++,Lbx[3]++;  //--Count all anti-/\_b  
               if( MyPrintMC_ ){
                   double BCan_En = BCand.energy();
                   double BCan_Px = BCand.px();
                   double BCan_Py = BCand.py();
                   double BCan_Pz = BCand.pz();
                   double BCan_Mass = sqrt( BCan_En*BCan_En - BCan_Px*BCan_Px - BCan_Py*BCan_Py  - BCan_Pz*BCan_Pz);
                   std::cout <<" -->Found lambda_b " << std::endl;
                   std::cout <<" genparticle " << k << " has pdgid = " << BCand.pdgId() 
                             <<" Ndaughter: "<<BCand.numberOfDaughters()
                             <<" B_Energy = "  <<BCan_En
                             <<" B_Mass = "<<BCan_Mass<< std::endl;
                   for(uint i = 0; i < BCand.numberOfDaughters(); i++ )
                       std::cout<<" lamb_b "<<i<<" daughterid: "<<BCand.daughter(i)->pdgId();
                   std::cout <<" "<< std::endl;
               }                               //--end if( MyPrintMC )
              
               if( BCand.numberOfDaughters()==2 ){
                  for( uint i = 0; i < BCand.numberOfDaughters(); i++){
                  
                     //--------JPsi candidate-----------------
                     bool up0=false, un0=false;
                     if( abs(BCand.daughter(i)->pdgId())== 443 ){
                         if( found_JPsi )std::cout <<" double_count: Found J/Psi " << std::endl;
                         if( MyPrintMC_ )
                             std::cout <<" Found J/Psi " << std::endl;
                             
                         const reco::Candidate * genDauJPsi = BCand.daughter(i);  //--JPsi cand
                         
                         if( MyPrintMC_ )
                             std::cout <<" J/Psi " << i << " has pdgid = " << genDauJPsi->pdgId() 
                                       <<" Ndaughter: "<<genDauJPsi->numberOfDaughters()<< std::endl;
                         
                         if( genDauJPsi->numberOfDaughters() == 2 ){
                                          
                            for( uint j = 0; j < genDauJPsi->numberOfDaughters(); j++){
                               if( MyPrintMC_ )
                                   std::cout<<" J/Psi "<<j<<" daughterid: "<< genDauJPsi->daughter(j)->pdgId()<< std::endl;                                                          
                               if( !up0 && genDauJPsi->daughter(j)->pdgId() == 13 )up0 = true;  
                               if( !un0 && genDauJPsi->daughter(j)->pdgId() ==-13 )un0 = true;    
                            }                  //--end for ( uint j = 0; j < genDauJPsi->numberOfDaughters(); j++){
                         }                     //--end if( genDauJPsi->numberOfDaughters() == 2){
                         if( up0 && un0 )found_JPsi = true;                           
                     }                         //--end if( abs(BCand.daughter(i)->pdgId())== 443 ){
                                          
                     //----------Psi prime (uu) candidate---------------
                     up0=false, un0=false;
                     if( abs(BCand.daughter(i)->pdgId())== 100443 ){
                         if( found_Psip1 )std::cout <<" double_count: Found Psi_p1 " << std::endl;
                         if( MyPrintMC_ )
                             std::cout <<" Found Psi_p1 " << std::endl;
                             
                         const reco::Candidate * genDauPsip1 = BCand.daughter(i);
                         
                         if( MyPrintMC_ )
                             std::cout <<" Psi_p1 " << i << " has pdgid = " << genDauPsip1->pdgId() 
                                       <<" Ndaughter: "<<genDauPsip1->numberOfDaughters()<< std::endl;
                                       
                         if( genDauPsip1->numberOfDaughters() == 2 ){
                            for( uint j = 0; j < genDauPsip1->numberOfDaughters(); j++){
                               if( MyPrintMC_ )
                                   std::cout<<" Psi_p1 "<<j<<" daughterid: "<< genDauPsip1->daughter(j)->pdgId()<< std::endl;                              
                               if( !up0 && genDauPsip1->daughter(j)->pdgId() == 13 )up0 = true;   
                               if( !un0 && genDauPsip1->daughter(j)->pdgId() ==-13 )un0 = true;   
                            }                  //--end for ( uint j = 0; j < genDauPsip1->numberOfDaughters(); j++){
                         }                     //--end if( genDauPsip1->numberOfDaughters() == 2){
                         if( up0 && un0 )found_Psip1 = true;
                     }                         //--end if( abs(BCand.daughter(i)->pdgId())== 100  443 ){

                     //-------Psi prime (JPsi pi pi) candidate----------------    
                     up0=false, un0=false;
                     bool pip0=false, pin0=false;
                     if( abs(BCand.daughter(i)->pdgId())== 100443 ){
                         if( found_Psip2 )std::cout <<" double_count: Found Psi_p2 " << std::endl;
                         if( MyPrintMC_ )
                             std::cout <<" Found Psi_p2 " << std::endl;
                             
                         const reco::Candidate * genDauPsip2 = BCand.daughter(i);
                         
                         if( MyPrintMC_ )
                             std::cout <<" Psi_p2 " << i << " has pdgid = " << genDauPsip2->pdgId() 
                                       <<" Ndaughter: "<<genDauPsip2->numberOfDaughters()<< std::endl;
                                       
                         if( genDauPsip2->numberOfDaughters() == 3 ){
                            for( uint j = 0; j < genDauPsip2->numberOfDaughters(); j++){
                               if( MyPrintMC_ )
                                   std::cout<<" Psi_p2 "<<j<<" daughterid: "<< genDauPsip2->daughter(j)->pdgId()<< std::endl;                               
                               if( !pip0 && genDauPsip2->daughter(j)->pdgId() == 211 )pip0 = true;  
                               if( !pin0 && genDauPsip2->daughter(j)->pdgId() ==-211 )pin0 = true;  
                               
                               if( abs(genDauPsip2->daughter(j)->pdgId()) == 443 ){
                                  if( MyPrintMC_ )
                                      std::cout <<" Found J/Psi_p " << std::endl;
                                      
                                  const reco::Candidate * genDauJPsi_p = genDauPsip2->daughter(j);  //--JPsi cand
                                  
                                  if( MyPrintMC_ )
                                      std::cout <<" JPsi_p " << i << " has pdgid = " << genDauJPsi_p->pdgId() 
                                                <<" Ndaughter: "<<genDauJPsi_p->numberOfDaughters()<< std::endl;
                                                
                                   //--------JPsi_p candidate-----------------
                                   if( genDauJPsi_p->numberOfDaughters() == 2){
                                      for( uint n = 0; n < genDauJPsi_p->numberOfDaughters(); n++){
                                         if( MyPrintMC_ )
                                             std::cout<<" JPsi_p "<<n<<" daughterid: "<< genDauJPsi_p->daughter(n)->pdgId()<< std::endl;                                          
                                         if( !up0 && genDauJPsi_p->daughter(n)->pdgId() == 13 )up0 = true;   
                                         if( !un0 && genDauJPsi_p->daughter(n)->pdgId() ==-13 )un0 = true;    
                                      }        //--end for ( uint n = 0; n < genDauJPsi_p->numberOfDaughters(); n++){
                                   }           //--end if( genDauJPsi_p->numberOfDaughters() == 2){          
                               }               //--end if( abs(genDauPsip2->daughter(j)->pdgId())== 443 ){
                            }                  //--end for ( uint j = 0; j < genDauPsip2->numberOfDaughters(); j++){
                         }                     //--end if( genDauPsip2->numberOfDaughters() == 3){
                         if( pip0 && pin0 && up0 && un0 )found_Psip2 = true;
                     }                         //--end if( abs(BCand.daughter(i)->pdgId())== 100  443 ){
                                              
                      //--------Lambda_0 candidate---------------------
                      bool p = false, ap = false;
                      pin0 = false, pip0 = false;
                      if( abs(BCand.daughter(i)->pdgId()) == 3122 ){  
                          if( found_LambP_0 )std::cout<<" double_count: lambdaP_0 "<< std::endl;     
                          if( found_LambA_0 )std::cout<<" double_count: lambdaA_0 "<< std::endl; 
                          if( MyPrintMC_ )
                              std::cout<< " lambdaP_0 "<< std::endl; 
                              
                          const reco::Candidate * genDauLamb_0 = BCand.daughter(i); 
                          
                          if( MyPrintMC_ )
                              std::cout <<" Lamb_0   " << i << " has pdgid = " << genDauLamb_0->pdgId() 
                                        <<" Ndaughter: "<<genDauLamb_0->numberOfDaughters()<<std::endl;
                                                    
                          histo_NLambda0 ->Fill(genDauLamb_0->numberOfDaughters());
                          
                          for( uint j = 0; j < genDauLamb_0->numberOfDaughters(); j++)
                               histo_L0daug   ->Fill(genDauLamb_0->daughter(j)->pdgId());

//------This section below is not working, in the new HLT MC (no g4sim)-------                                           
                           if( genDauLamb_0->numberOfDaughters() == 2){ 
                              //--Lambda 0 particle
                              if( BCand.daughter(i)->pdgId() == 3122 ){    
                                                    
                                  for( uint j = 0; j < genDauLamb_0->numberOfDaughters(); j++){
                                     if( MyPrintMC_ )
                                         std::cout<<" "<<j<<" lambP_0 daughterid: "<< genDauLamb_0->daughter(j)->pdgId()<<std::endl;                                     
                                     if( !p    && genDauLamb_0->daughter(j)->pdgId() == 2212 ) p    = true;   
                                     if( !pin0 && genDauLamb_0->daughter(j)->pdgId() ==-211  ) pin0 = true;
                                  }               //--end for ( uint j = 0; j < genDaulamb_0->numberOfDaughters(); j++){
                              }                  //--end if( BCand.daughter(i)->pdgId() == 3122 ){
                              //--Lambda 0 Antiparticule
                              if( BCand.daughter(i)->pdgId() == -3122 ){
                                 for( uint j = 0; j < genDauLamb_0->numberOfDaughters(); j++){
                                     if( MyPrintMC_ )
                                         std::cout<<" "<<j<<" lambA_0 daughterid: "<< genDauLamb_0->daughter(j)->pdgId()<<std::endl;
                                     if( !ap   && genDauLamb_0->daughter(j)->pdgId() == -2212 ) ap   = true;   
                                     if( !pip0 && genDauLamb_0->daughter(j)->pdgId() ==   211 ) pip0 = true;
                                 }               //--end for ( uint j = 0; j < genDaulamb_0->numberOfDaughters(); j++){
                              }                  //--end if( BCand.daughter(i)->pdgId() == -3122 ){
                           }                     //--end if( genDauLamb_0->numberOfDaughters() == 2){
                           if(  p && pin0 )found_LambP_0 = true;    
                           if( ap && pip0 )found_LambA_0 = true;
//-------Temporaly solution---------------------------
//--                          if( BCand.daughter(i)->pdgId() == 3122 )found_LambP_0 = true;
//--                          if( BCand.daughter(i)->pdgId() ==-3122 )found_LambA_0 = true;
                      }                         //--end if( abs(BCand.daughter(i)->pdgId())== 3122 )
                                                                
                  }                            //--end for( uint i = 0; i < BCand.numberOfDaughters(); i++)
              }                                //--end if(BCand.numberOfDaughters()==2)
           }                                   //--end if ( abs(BCand.pdgId())==5122 )
        }                                      //--end for( size_t k = 0; k < genParticles->size(); k++ )                
        if( MyPrintMC_ ){      
            if( found_JPsi    )std::cout<<found_JPsi  <<" Found J/Psi: "<<std::endl;
            else               std::cout<<found_JPsi  <<" NO Found J/Psi: "<<std::endl;
            
            if( found_Psip1   )std::cout<<found_Psip1 <<" Found Psip1: "<<std::endl;
            else               std::cout<<found_Psip1 <<" NO Found Psip1: "<<std::endl;
            
            if( found_Psip2   )std::cout<<found_Psip2 <<" Found Psip2: "<<std::endl;
            else               std::cout<<found_Psip2 <<" NO Found Psip2: "<<std::endl;
            
            if( found_LambP_0 )std::cout<<found_LambP_0 <<" Found LambP_0: "<<std::endl;
            else               std::cout<<found_LambP_0 <<" NO Found LambP_0: "<<std::endl;
  
            if( found_LambA_0 )std::cout<<found_LambA_0 <<" Found LambA_0: "<<std::endl;
            else               std::cout<<found_LambA_0 <<" NO Found LambA_0: "<<std::endl;
        
        }     //--end if( MyPrintMC )
        
        //-----Counting Generated Lambda_b ------        
        if( found_JPsi  && found_LambP_0 )Lb[1]++,Tb[1]++;
        if( found_JPsi  && found_LambA_0 )Lb[2]++,Tb[2]++;             
        Lb[3] = Lb[1] + Lb[2], Tb[3] = Tb[1] + Tb[2];
        
        if( found_Psip1 && found_LambP_0 )Lb[4]++,Tb[4]++;
        if( found_Psip1 && found_LambA_0 )Lb[5]++,Tb[5]++;
        Lb[6] = Lb[4] + Lb[5], Tb[6] = Tb[4] + Tb[5];
        
        if( found_Psip2 && found_LambP_0 )Lb[7]++,Tb[7]++;
        if( found_Psip2 && found_LambA_0 )Lb[8]++,Tb[8]++;
        Lb[9] = Lb[7] + Lb[8], Tb[9] = Tb[7] + Tb[8];
        
        if( MyPrintMC_ ){
            std::cout<<Tb[1]<<"-J/Psi+L0: "<<Lb[1]<<"  "<<Tb[2]<<"-J/Psi+AL0: "<<Lb[2]<<" Total:-"<<Tb[3]<<" "<<Lb[3]<<std::endl;
            std::cout<<Tb[4]<<"-Psip1+L0: "<<Lb[4]<<"  "<<Tb[5]<<"-Psip1+AL0: "<<Lb[5]<<" Total:-"<<Tb[6]<<" "<<Lb[6]<<std::endl;
            std::cout<<Tb[7]<<"-Psip2+L0: "<<Lb[7]<<"  "<<Tb[8]<<"-Psip2+AL0: "<<Lb[8]<<" Total:-"<<Tb[9]<<" "<<Lb[9]<<std::endl;
        }
        if( (Tb[3]+Tb[6]+Tb[9])>1 ){
            std::cout<<"double_count: "<<Tb[1]<<"-J/Psi+L0: "<<Lb[1]<<"  "<<Tb[2]<<"-J/Psi+AL0: "<<Lb[2]<<" Total:-"<<Tb[3]<<" "<<Lb[3]<<std::endl;
            std::cout<<"double_count: "<<Tb[4]<<"-Psip1+L0: "<<Lb[4]<<"  "<<Tb[5]<<"-Psip1+AL0: "<<Lb[5]<<" Total:-"<<Tb[6]<<" "<<Lb[6]<<std::endl;
            std::cout<<"double_count: "<<Tb[7]<<"-Psip2+L0: "<<Lb[7]<<"  "<<Tb[8]<<"-Psip2+AL0: "<<Lb[8]<<" Total:-"<<Tb[9]<<" "<<Lb[9]<<std::endl;
        }
        for( uint i = 1; i < 10; i++)
             histo_MCgen->Fill( i,Tb[i]);
        int jBin = 10;
        for( uint i = 0; i < 4; i++){
             jBin++;
             histo_MCgen->Fill( jBin,Lbx[i]);
             //std::cout<<" i: "<<i<<",jBin: "<<jBin<<std::endl;
        }
 
       if( MyPrintMC_ )
           std::cout <<"-----------------------"<< std::endl;
   }                                           //--end if( doMC_)
//---------------------------------------------------------------------------------------------

   float muon_sig   = muon_mass_c*1.e-6;
   float pion_sig   = pion_mass_c*1.e-6;
   float proton_sig = proton_mass_c*1.e-6;
   float chi = 0., ndf = 0.;

for( unsigned int veeindex = 0; veeindex < theVees.size(); veeindex++ ){
   double Mlambda0 = theVees[veeindex].mass();
   histo_lam[0]->Fill( Mlambda0 );
       
   //--Tracks from the vee (Lambda_0)
   reco::TrackRef protonveetrk = (dynamic_cast<reco::RecoChargedCandidate *> 
                                 (theVees[veeindex].daughter(0) ) )->track();
   reco::TrackRef pionveetrk  = (dynamic_cast<reco::RecoChargedCandidate *>
                                (theVees[veeindex].daughter(1) ) )->track();
  
   //--Reconstruct Lambda Vertex
   reco::TransientTrack pionTveetrk(     pionveetrk, &(*bFieldHandle) );
   reco::TransientTrack protonTveetrk( protonveetrk, &(*bFieldHandle) );
      
   KinematicParticleFactoryFromTransientTrack pFactory;
 
   std::vector<RefCountedKinematicParticle> LamParticle;
   LamParticle.push_back(pFactory.particle( pionTveetrk,   pion_mass_c,   chi, ndf, pion_sig ));
   LamParticle.push_back(pFactory.particle( protonTveetrk, proton_mass_c, chi, ndf, proton_sig ));               
   KinematicFitDriver LamRec( LamParticle, "Lambda" );
   if( LamRec.isValid() ){

      //--Output lambda Vtx fit variables
      double MlamVFit   = LamRec.mass();  //--valid vertex Lambda Mass
      double lamV_CL    = ChiSquaredProbability( LamRec.chi2(), LamRec.ndof() );
     
      histo_lam[1] ->Fill( Mlambda0 );
      histo_diff[0]->Fill( Mlambda0 - MlamVFit );       //--fitted Lambda mass
  
      //--Here do Lambda0 Mass constrains Fit. (proton + pion) to Lambda
      bool lamMcut = std::abs( MlamVFit - lambda_mass_c ) < 0.01;      //--10 MeV  around Nominal Lambda
      LamRec.AddMassConstraint( lambda_mass_c, lambda_mass_c*1.e-6 );
      if( LamRec.isValid() ){

      //--Output lambda mass constrain fit variables
      double MlamMFit = LamRec.mass();                                         //--Valid vertex Lam Mass 
      double lamVM_CL = ChiSquaredProbability( LamRec.chi2(), LamRec.ndof() ); //--Confidence Level Lam Vtx 
      reco::Particle::LorentzVector lamVM_P4 = LamRec.P4();            //--lambda0 4-Momentum with Mass constrains
      reco::Particle::LorentzVector piVM_P4  = LamRec.P4FirstChild();  //--pion Momentum with Vtx constrains
      reco::Particle::LorentzVector prVM_P4  = LamRec.P4NextChild();   //--proton      
      reco::Particle::Point lamVM_vx = LamRec.VertexDecay();           //--Vx,Vy,Vz
      reco::Vertex::CovarianceMatrix lamVM_Cov = LamRec.VertexCov();
      double lamVMx = lamVM_vx.x();
      double lamVMy = lamVM_vx.y();
      double lamVMz = lamVM_vx.z();
      double lamDist= sqrt( lamVMx*lamVMx + lamVMy*lamVMy + lamVMz*lamVMz );
      double lamR   = sqrt( lamVMx*lamVMx + lamVMy*lamVMy );

      histo_lam[2]->Fill( MlamMFit );
      histo_lam[3]->Fill( Mlambda0 );
      
      //--Calculate Kaon_short invariant mass
    //double lamE  = prVM_P4.E()  + piVM_P4.E();
    //double lamEpr = sqrt( proton_mass_c*proton_mass_c + prVM_P4.P()*prVM_P4.P() );
      double lamEpr = sqrt( pion_mass_c*pion_mass_c + prVM_P4.P()*prVM_P4.P() );  //--Assuming PION for this Proton
      double lamEpi = sqrt( pion_mass_c*pion_mass_c + piVM_P4.P()*piVM_P4.P() );
      double lamE   = lamEpr + lamEpi;      
      double lamPx  = prVM_P4.Px() + piVM_P4.Px();    //--Momentum at Vtx
      double lamPy  = prVM_P4.Py() + piVM_P4.Py();
      double lamPz  = prVM_P4.Pz() + piVM_P4.Pz();
      
      reco::Particle::LorentzVector K0s(0.0,0.0,0.0,0.0);  //--Fill Lorentz K0short P4                  
      K0s.SetPxPyPzE( lamPx, lamPy, lamPz, lamE );
      //double MK0s = K0s.M();                    //--same as myLamM
      //double myLamM  = sqrt( lamE*lamE - lamPx*lamPx - lamPy*lamPy - lamPz*lamPz );    //--Kshort when Proton is a Pion
      double myLamM  = K0s.M();
      double myfLamM = sqrt( lamVM_P4.E()*lamVM_P4.E() - lamVM_P4.Px()*lamVM_P4.Px()   //--Mass constrains Momentum
                                                       - lamVM_P4.Py()*lamVM_P4.Py()
                                                       - lamVM_P4.Pz()*lamVM_P4.Pz() );                    
      if( MyPrint_ )                       
          std::cout<<MlamMFit<<" M:"<<myfLamM<<" Class:"<<Mlambda0-myLamM
          <<" Fit:"<<MlamVFit-myLamM<<" <--"<<myLamM<<std::endl;
      
      //--Save proton & pion variables      
      double prQ        = protonveetrk->charge();
      double prEta      = protonveetrk->eta();
      double prPhi      = protonveetrk->phi();
      double prPx       = prVM_P4.Px();  //protonveetrk->p4().Pt();  Transverse
      double prPy       = prVM_P4.Py();  //protonveetrk->p4().P();
      double prPz       = prVM_P4.Pz();  //protonveetrk->p4().Pz();  Longitudinal 
      double prd0       = protonveetrk->d0();
      double prdz       = protonveetrk->dz();   
      double prNtrkHits = protonveetrk->numberOfValidHits();
      
      double piQ        = pionveetrk->charge();
      double piEta      = pionveetrk->eta();
      double piPhi      = pionveetrk->phi();
      double piPx       = piVM_P4.Px();    //pionveetrk->p4().Pt();   Transverse
      double piPy       = piVM_P4.Py();    //pionveetrk->p4().P();
      double piPz       = piVM_P4.Pz();    //pionveetrk->p4().Pz();  Longitudinal
      double pid0       = pionveetrk->d0();
      double pidz       = pionveetrk->dz();   
      double piNtrkHits = pionveetrk->numberOfValidHits();
      
      //--std::cout<<prPt<<" prPt:"<<sqrt(prVM_P4.Px()*prVM_P4.Px()+prVM_P4.Py()*prVM_P4.Py())
      //--         <<piPt<<" piPt:"<<sqrt(piVM_P4.Px()*piVM_P4.Px()+piVM_P4.Py()*piVM_P4.Py())<<std::endl;
      
      fill_lambda0(prQ, prEta, prPhi, prPx, prPy, prPz, prd0, prdz, prNtrkHits,
                   piQ, piEta, piPhi, piPx, piPy, piPz, pid0, pidz, piNtrkHits,
                   Mlambda0, lamV_CL, MlamVFit, lamDist, lamR, myLamM,
                   MlamMFit, lamVM_CL, lamVM_P4.px(), lamVM_P4.py(), lamVM_P4.pz() );
      
      //--Here Dimuon Loop
      int iComb = 0;
    //for(reco::MuonCollection::const_iterator iMuon=allmuons->begin();iMuon!= allmuons->end();++iMuon){ 
      for(std::vector<pat::Muon>::const_iterator iMuon=allmuons->begin();iMuon!= allmuons->end();++iMuon){   //--iMuonP
         int iQ      = iMuon->charge();
         double iEta = iMuon->eta();
         double iPhi = iMuon->phi();
         //double iPx  = iMuon->p4().Px();
         //double iPy  = iMuon->p4().Py();
         //double iPz  = iMuon->p4().Pz();
         double iCal = iMuon->caloCompatibility();
         double iSeg = muon::segmentCompatibility(*iMuon);
         double iIso = iMuon->isolationR03().sumPt;
         reco::TrackRef iMuTrackRef = iMuon->innerTrack(); //--.index();
        
         int muSA=0, muGL=0, muTK=0;          //--Muon Id
         if( iMuon->isStandAloneMuon() )muSA = 1;   //--2^0
         if( iMuon->isGlobalMuon()     )muGL = 10;  //--2^1
         if( iMuon->isTrackerMuon()    )muTK = 100; //--2^2  
         int imuId = muTK + muGL + muSA;
        
       //for(reco::MuonCollection::const_iterator jMuon=iMuon+1;jMuon!=allmuons->end();++jMuon){ 
         for(std::vector<pat::Muon>::const_iterator jMuon=iMuon+1;jMuon!=allmuons->end();++jMuon){  //--iMuonN
            int jQ      = jMuon->charge();
            double jEta = jMuon->eta();
            double jPhi = jMuon->phi();
            //double jPx  = jMuon->p4().Px();
            //double jPy  = jMuon->p4().Py();
            //double jPz  = jMuon->p4().Pz();
            double jCal = jMuon->caloCompatibility();
            double jSeg = muon::segmentCompatibility(*jMuon);
            double jIso = jMuon->isolationR03().sumPt;
            reco::TrackRef jMuTrackRef = jMuon->innerTrack();   //--.index();
            
            muSA=0, muGL=0, muTK=0;            //--Muon Id
            if( jMuon->isStandAloneMuon() )muSA = 1;   //--2^0
            if( jMuon->isGlobalMuon()     )muGL = 10;  //--2^1
            if( jMuon->isTrackerMuon()    )muTK = 100; //--2^2
            int jmuId = muTK + muGL + muSA;
            
            bool goodDimuon = (imuId==111 && jmuId>=11)    //--1 G and 1 (SG or T or ST or G)
                            ||(jmuId==111 && imuId>=11);
            if( ( iQ*jQ < 0 )&& goodDimuon ){  //--Opp. charge & muon Id & unique Trk
                                               //--1  : S
                                               //--10 : G
                                               //--11 : G-S
                                               //--100: T
                                               //--101: T-S
                                               //--110: T-G  doesn't exist
                                               //--111: T-G-S
              double Muu = (iMuon->p4() + jMuon->p4() ).M();
              histo_dimuM[4]->Fill( Muu );
              histo_lam[4]  ->Fill( MlamVFit );
              
              bool repTrk = false;     //--if true I reject Same tracks
              if( protonveetrk == iMuTrackRef || protonveetrk == jMuTrackRef ||
                  pionveetrk   == iMuTrackRef || pionveetrk   == jMuTrackRef )repTrk = true;
              if( repTrk ){
                  //std::cout<<"Hec: repeated Trk. "<<std::endl;
                  if( iQ*jQ < 0 ){
                      histo_dimuM[5]->Fill( Muu );
                      histo_lam[5]  ->Fill( MlamVFit );
                  } else {
                      histo_dimuM[6]->Fill( Muu );
                      histo_lam[6]  ->Fill( MlamVFit );
                  }
              }
              if( !repTrk ){
              
              int iNtrkHits = iMuTrackRef->numberOfValidHits();
              double id0    = iMuTrackRef->d0();
              double idz    = iMuTrackRef->dz();
              double iKink  = iMuon->combinedQuality().trkKink;//--Chi squared
          
              int jNtrkHits = jMuTrackRef->numberOfValidHits();
              double jd0    = jMuTrackRef->d0();
              double jdz    = jMuTrackRef->dz();
              double jKink  = jMuon->combinedQuality().trkKink;//--Chi squared
              
              //--Do here the Tight Muon Variables  Jan 2013
              //--These 2 are not working in pat::
              //double iMuPFmu     = iMuon->isPFMuon();                 
              //double iMuBestTrk  = iMuon->muonBestTrack()->dxy(primary.BestVertex()->position());
              double iNormChi2   = 99;
              double iMuMuHits   = -1;
              if( iMuon->isGlobalMuon() ){
                  iNormChi2 = iMuon->globalTrack()->normalizedChi2();
                  iMuMuHits = iMuon->globalTrack()->hitPattern().numberOfValidMuonHits();
              }
              double iMuStations = iMuon->numberOfMatchedStations();
              double iMuPxHits   = iMuon->innerTrack()->hitPattern().numberOfValidPixelHits();
              double iMuTrkLayer = iMuon->track()->hitPattern().trackerLayersWithMeasurement();
              
              double jNormChi2   = 99;
              double jMuMuHits   = -1;
              if( jMuon->isGlobalMuon() ){
                  jNormChi2 = jMuon->globalTrack()->normalizedChi2();
                  jMuMuHits = jMuon->globalTrack()->hitPattern().numberOfValidMuonHits();
              }
              double jMuStations = jMuon->numberOfMatchedStations();
              double jMuPxHits   = jMuon->innerTrack()->hitPattern().numberOfValidPixelHits();
              double jMuTrkLayer = jMuon->track()->hitPattern().trackerLayersWithMeasurement();
              
              fill_tightMuons(iNormChi2,iMuMuHits,iMuStations,iMuPxHits,iMuTrkLayer,iKink
                             ,jNormChi2,jMuMuHits,jMuStations,jMuPxHits,jMuTrkLayer,jKink);
                              
              if( MyPrint_ ){
              std::cout<<"Tight Muon Variables: "<<iNormChi2    <<"  "<<jNormChi2  
                       <<"  "                    <<iMuMuHits    <<"  "<<jMuMuHits  
                       <<"  "                    <<iMuStations  <<"  "<<jMuStations
                       <<"  "                    <<iMuPxHits    <<"  "<<jMuPxHits  
                       <<"  "                    <<iMuTrkLayer  <<"  "<<jMuTrkLayer
                       <<std::endl;
               }      
              //--Re-fit dimuons to a common Vertex
              reco::TransientTrack iMuTT( iMuTrackRef, &(*bFieldHandle) );
              reco::TransientTrack jMuTT( jMuTrackRef, &(*bFieldHandle) );
           
              KinematicParticleFactoryFromTransientTrack pFactory;
           
              std::vector<RefCountedKinematicParticle> uuParticle;
              uuParticle.push_back(pFactory.particle( iMuTT, muon_mass_c, chi, ndf, muon_sig ));
              uuParticle.push_back(pFactory.particle( jMuTT, muon_mass_c, chi, ndf, muon_sig ));
          
              // Trajectory states to calculate DCA for the 2 tracks
              FreeTrajectoryState iState = pFactory.particle( iMuTT, muon_mass_c, chi, ndf, muon_sig )->currentState().freeTrajectoryState();
              FreeTrajectoryState jState = pFactory.particle( jMuTT, muon_mass_c, chi, ndf, muon_sig )->currentState().freeTrajectoryState();

              //--Measure distance between tracks at their closest approach
              ClosestApproachInRPhi cApp;           
              double dca_uu = -9999;              //--Flag  Closest Approach distance between muons
              GlobalPoint cxPt;                   //--cxPt.x(): Crosing point in X
              cApp.calculate(iState, jState);     //--cxPt.y(): Crosing point in Y 
              if( !cApp.status() ) {              //--cxPt.z(): Crosing point in Z
                std::cout<<"Hec: uu bad status capp"<<std::endl;
                //--continue;    Need to check where this continue will go ??
              } else {
                dca_uu = std::abs( cApp.distance() );
                cxPt = cApp.crossingPoint();
              }
              if( MyPrint_ ){
                  std::cout<<"uu dca = "<<dca_uu<<" CroossPoint ="<<cxPt<<std::endl;
                  std::cout<<" mui= "<<imuId<<" muj= "<<jmuId<<" mass ="<<Muu<<std::endl;
              } 
              //--r = sqrt( x^2 + y^2 )
              //--Cut on fiducial tracking volume [-120 cm < r < 120 cm  and -260 cm < Z < 260 cm]
              //--Pixel [-20 cm < r < 20 cm  and -60 cm < Z < 60 cm] to run vertex fit
              if( ( dca_uu >= 0. && dca_uu <= 5 ) && 
                  ( sqrt( cxPt.x()*cxPt.x() + cxPt.y()*cxPt.y() ) < 120.) &&
                    std::abs(cxPt.z()) < 300. ){
	 
              //--Here Dimuon Vertex
              KinematicFitDriver uuRec( uuParticle, "MuMu" );
              if( !uuRec.isValid() ){
                  if( MyPrint_ )
                      std::cout<<"Hec: Mass_uu: "<<Muu<<" DiMuV BadFit:"<<std::endl;
                  //--continue;   //--I don't know where this goes
              } else {              
              //--Output Vertex fit variables 
              double MuuVFit = uuRec.mass();
              double uuV_CL  = ChiSquaredProbability( uuRec.chi2(), uuRec.ndof() );
              
              //--Here Dimu Mass constraint to massC_value (= jpsi)
              double massC_value = 0;
              bool uu1sMcut = std::abs( MuuVFit - jpsi_mass_c )<0.2;    //--200 MeV
              bool uu2sMcut = std::abs( MuuVFit - psi2s_mass_c)<0.2;    //--200 MeV
              if( uu1sMcut )massC_value = jpsi_mass_c;
              if( uu2sMcut )massC_value = psi2s_mass_c;
              
              double MuuMFit = -1;   //--Mass Fitted Dimuon mass 
              double uuVM_CL = -1;   //--Confidence Level uu Vtx
              if( uu1sMcut || uu2sMcut ){       
                uuRec.AddMassConstraint( massC_value, massC_value*1.e-6 );          
                if( uuRec.isValid() ){ 
                  MuuMFit = uuRec.mass();
                  uuVM_CL = ChiSquaredProbability( uuRec.chi2(), uuRec.ndof() );
                }  //--End good J/Psi or Psi(2s) mass constrain            
              }    //--End uu1sMcut || uu2sMcut
              
              //--Dimuon fitted variables
              reco::Particle::LorentzVector uuV_P4  = uuRec.P4();            //--DiMu 4-Momentum
              reco::Particle::LorentzVector iMuV_P4 = uuRec.P4FirstChild();  //--Mui
              reco::Particle::LorentzVector jMuV_P4 = uuRec.P4NextChild();   //--Muj
              reco::Particle::Point uuV_vx = uuRec.VertexDecay();            //--Vx,Vy,Vz 
              reco::Vertex::CovarianceMatrix uuV_Cov = uuRec.VertexCov();
              double uuVx     = uuV_vx.x();
              double uuVy     = uuV_vx.y();
              double uuVz     = uuV_vx.z();
              double uuDist   = sqrt( uuVx*uuVx + uuVy*uuVy + uuVz*uuVz );
              double uuR      = sqrt( uuVx*uuVx + uuVy*uuVy );

              if( MyPrint_ ){
                  std::cout<<"Hec: Collection uu M= "<<Muu<<" Vtx ReFit: "<<MuuVFit<<std::endl;
                  std::cout<<"Hec: CL = "<<uuV_CL
                           <<" Vtx (x,y,z): ("<<uuVx<<","<<uuVy<<","<<uuVz<<")"<<std::endl;
              }     
              //--Check dimuon invariant mass
              //double uuE  = iMuV_P4.E()  + jMuV_P4.E();
              //double iuE = sqrt( muon_mass_c*muon_mass_c + iMuV_P4.P()*iMuV_P4.P() );
              //double juE = sqrt( muon_mass_c*muon_mass_c + jMuV_P4.P()*jMuV_P4.P() );
              //double uuE = iuE + juE;      
              //double uuPx  = iMuV_P4.Px() + jMuV_P4.Px();    //--Momentum at Vtx
              //double uuPy  = iMuV_P4.Py() + jMuV_P4.Py();
              //double uuPz  = iMuV_P4.Pz() + jMuV_P4.Pz();
      
              double iPx  = iMuV_P4.Px();
              double iPy  = iMuV_P4.Py();
              double iPz  = iMuV_P4.Pz();
              
              double jPx  = jMuV_P4.Px();
              double jPy  = jMuV_P4.Py();
              double jPz  = jMuV_P4.Pz();
         
              //double myuuM  = sqrt( uuE*uuE - uuPx*uuPx - uuPy*uuPy - uuPz*uuPz );    //--dimuon at Vertex
              //double myfuuM = sqrt( uuV_P4.E()*uuV_P4.E() - uuV_P4.Px()*uuV_P4.Px()   //--Mass constrains Momentum
              //                                            - uuV_P4.Py()*uuV_P4.Py()
              //                                            - uuV_P4.Pz()*uuV_P4.Pz() );
              //if( MyPrint_ )
              //    std::cout<<MuuMFit<<" M:"<<myfuuM<<" Class:"<<Muu-myuuM<<" Fit:"<<MuuVFit-myuuM<<std::endl;
                  
              histo_dimuM[7]->Fill( Muu );
              histo_lam[7]  ->Fill( MlamVFit );
              histo_diff[1] ->Fill( Muu - MuuVFit );
                      
              //--Dimuon Vertex Isolation [Feb 10, 2012]
              reco::Vertex uuVTX = *uuRec.RKVertex();
              if( MyPrint_ )
                  std::cout<<"-----Vertex Isolation---"<<uuV_vx.x()
                                                  <<" "<<uuV_vx.y()
                                                  <<" "<<uuV_vx.z()<<std::endl;
              double closestTRK   =  100, closestSigmaD0 = -100, closestD0 = -100,
                     closestTRK3d =  100, closestSigmaDz = -100, closestDz = -100,
                     closestId = 0;
              for(unsigned int trkindex = 0; trkindex < theTransTracks.size(); trkindex++) {
                 double SigmaDz = 0, SigmaD0 = 0, TrakDz, TrakD0, trkip2d;
                 double trkip3d = SignificanceAbsoluteImpactParameter3D( theTransTracks[trkindex], uuVTX,
                                  SigmaDz, TrakDz, SigmaD0, TrakD0,trkip2d );

                 if( iMuTrackRef == theTrackRefs[trkindex] ||
                     jMuTrackRef == theTrackRefs[trkindex] ){
                     if( MyPrint_ )
                         std::cout<<trkindex<<" muons trkip3d ="<<trkip3d
                                            <<" SigmaDz ="      <<SigmaDz
                                            <<" dz ="           <<TrakDz
                                            <<" SigmaD0 ="      <<SigmaD0
                                            <<" d0 ="           <<TrakD0
                                            <<" Significanc2d ="<<trkip2d
                                            <<" uuVtx CL ="     <<uuV_CL<<std::endl;                                        
                     histo_sig3d[0]->Fill( trkip3d );
                     histo_sig2d[0]->Fill( trkip2d );
                 } else if( pionveetrk   == theTrackRefs[trkindex]
                         || protonveetrk == theTrackRefs[trkindex] ){                                        
                     histo_sig3d[1]->Fill( trkip3d );
                     histo_sig2d[1]->Fill( trkip2d );
                 } else {                                       
                     histo_sig3d[2]->Fill( trkip3d );
                     histo_sig2d[2]->Fill( trkip2d );
                 }
                 //--I don't remove the lambda_0 products
                 double lamId_pi = 0, lamId_pr = 0;
                 if( pionveetrk == theTrackRefs[trkindex] )
                     lamId_pi = 1;
                 if( protonveetrk == theTrackRefs[trkindex] )
                     lamId_pr = 10;
                 if( iMuTrackRef== theTrackRefs[trkindex] || jMuTrackRef  == theTrackRefs[trkindex] )continue;   //--remove muons

                 if( std::abs(trkip2d) < std::abs(closestTRK) ){    //--cut in 2d significance
                     closestTRK     = trkip2d;
                     closestSigmaD0 = SigmaD0;
                     closestD0      = TrakD0;
                     closestSigmaDz = SigmaDz;
                     closestDz      = TrakDz;
                     closestTRK3d   = trkip3d;
                     closestId      = lamId_pr + lamId_pi;
                 }
          
              }     //--end  for(unsigned int trkindex                                      
              histo_sig3d[3]->Fill( closestTRK3d );
              histo_sig2d[3]->Fill( closestTRK );
              
              histo_closest->Fill( closestTRK );
              if( MyPrint_ ){
                  std::cout<<"-closestTRk ="  <<closestTRK
                           <<" closestSigmaD0 ="<<closestSigmaD0
                           <<" closestD0 ="   <<closestD0<<std::endl;
                  std::cout<<"---------------------------------------------------------"<<std::endl;
              }            
              //--Dimuon and Lambda-0 Vertex [ Vertexing 2 neutral is bad ? ]
              //-->std::vector<RefCountedKinematicParticle> Lambda_bParticle;
              //Lambda_bParticle.push_back( uuRec.RKParent() );      
              //-->Lambda_bParticle.push_back(pFactory.particle( iMuTT, muon_mass_c, chi, ndf, muon_sig ));
              //-->Lambda_bParticle.push_back(pFactory.particle( jMuTT, muon_mass_c, chi, ndf, muon_sig ));
              //-->Lambda_bParticle.push_back( LamRec.RKParent() );

              //--Check if Lambda0 comes from Dimuon Vertex (calculate cosine in 3D and 2D)
              GlobalVector Vector_uuV_to_lamVM( lamVM_vx.x()- uuV_vx.x(),
                                                lamVM_vx.y()- uuV_vx.y(),
                                                lamVM_vx.z()- uuV_vx.z() );                           
              GlobalVector Lambda_P3( lamVM_P4.px(), lamVM_P4.py(), lamVM_P4.pz() ); 
              double cos_uul3D = ( Lambda_P3.dot( Vector_uuV_to_lamVM ) )/
                                 ( Vector_uuV_to_lamVM.mag() * Lambda_P3.mag() );
              
              double cos_uul2D = ( (lamVM_vx.x()- uuV_vx.x())*lamVM_P4.px() + 
                                   (lamVM_vx.y()- uuV_vx.y())*lamVM_P4.py() )
                                 /sqrt( (lamVM_vx.x()- uuV_vx.x())*(lamVM_vx.x()- uuV_vx.x()) +
                                        (lamVM_vx.y()- uuV_vx.y())*(lamVM_vx.y()- uuV_vx.y()) )
                                 /sqrt( lamVM_P4.px()*lamVM_P4.px() +
                                        lamVM_P4.py()*lamVM_P4.py() );
                                        
              histo_cos3d->Fill( cos_uul3D );
              histo_cos2d->Fill( cos_uul2D );
              
              //--this cos_new3d and 2d doesn't work here because this function use a vertex as Primary
              //--which is ok but also need that my vertex has MotherDecayVertex
              //reco::Vertex uuVTX = *uuRec.RKVertex();
              //double cos_new3d = uuRec.CosAlpha(  uuVTX );
              //double cos_new2d = uuRec.CosAlpha2D(uuVTX );
              //std::cout<<"Hec: Cosine= 3D: "<<cos_uul3D<<" Ed: "<<cos_new3d<<std::endl;
              //std::cout<<"Hec: Cosine= 2D: "<<cos_uul2D<<" Ed: "<<cos_new2d<<std::endl;
              
              //double cos_new3d = uuRec.CosAlpha(  refitVertexPrim );  
              //double cos_new2d = uuRec.CosAlpha2d(refitVertexPrim );
              
              if( MyPrint_ ){
                  std::cout<<"Hec: Pointing uu_lambda0 cos(uul)= 3D: "<<cos_uul3D<<" 2D: "<<cos_uul2D<<std::endl;
                  std::cout<<"lambda: "<<LamRec.chi2()<<" Ndof: "<<LamRec.ndof()<<" CL: "<<lamVM_CL<<std::endl; 
                  std::cout<<"dimuon: chi**2: "<<uuRec.chi2()<<" Ndof: "<<uuRec.ndof()<<" CL: "<<uuV_CL<<std::endl;
                  std::cout<<" mui= "<<imuId<<" muj= "<<jmuId<<std::endl;
                  std::cout<<"Dimuon Vertex: Vx= "<<uuVx <<" Vy= "<<uuVy <<" Vz= "<<uuVz <<" Dist="<<uuDist <<" R="<<uuR<<std::endl;
                  std::cout<<"Lambda Vertex: Vx= "<<lamVMx<<" Vy= "<<lamVMy<<" Vz= "<<lamVMz<<" Dist="<<lamDist<<" R="<<lamR<<std::endl;
                  //std::cout<<"Puu dot Plam: "<<uu_lam<<" "<<cos_uu_lam<<std::endl;
              }

              //--Skip Events here:   ---***-----***-----***-----***-----***---
              //--if( uuV_CL>0.1 && lamV_CL>0.1 && 
              //--   MuuFit>2.8 && MuuFit<3.3 && 
              //--  MlamVFit>1.11 && MlamVFit<1.125 ){   //--dimuon CL + other cuts to avoid vertex crash
              if( cos_uul3D > 0.75 ){
              
              //--Check separation Lambda from Dimuon (L over Sigma) in 3D
              typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D_uuL;
              typedef ROOT::Math::SVector<double, 3> SVector3_uuL;
 
              SMatrixSym3D_uuL totalCov_uuL = lamVM_Cov + uuV_Cov;   //--3D Separation Dimuon - Lambda
              SVector3_uuL uuL_PrimarySep3D( uuV_vx.x() - lamVM_vx.x(),
                                             uuV_vx.y() - lamVM_vx.y(),
                                             uuV_vx.z() - lamVM_vx.z() ); 
              double uuL_L3D     = ROOT::Math::Mag( uuL_PrimarySep3D );
              double uuL_Sigma3D = sqrt(ROOT::Math::Similarity( totalCov_uuL, uuL_PrimarySep3D )) / uuL_L3D;
              
              //--Check separation Lambda from Dimuon (L over Sigma) in 2D  [need to check] ? ? ?
              //typedef ROOT::Math::SMatrix<double, 2, 2, ROOT::Math::MatRepSym<double, 2> > SMatrixSym2D_uuL;
              //typedef ROOT::Math::SVector<double, 2> SVector2_uuL;
              
              //SMatrixSym2D_uuL totalCov2D_uuL = lamVM_Cov[1][1] + uuV_Cov[1][1];   //--2D Separation Dimuon - Lambda
              //SVector2_uuL uuL_PrimarySep2D( uuV_vx.x() - lamVM_vx.x(), uuV_vx.y() - lamVM_vx.y() );
              //double uuL_L2D     = ROOT::Math::Mag( uuL_PrimarySep2D );
              
              //--Try here new Eduardo's Functions 
              //reco::Vertex uuVTX = *uuRec.RKVertex();  defined above     
              //double uuLoS_new3d = uuRec.SignificanceSeparation( uuVTX );
              //double L_new3d     = uuRec.Separation( refitVertexPrim );
              //double LoS_new2d   = uuRec.SignificanceSeparation2D( refitVertexPrim );
              //double L_new2d     = uuRec.Separation2D( refitVertexPrim );
              //std::cout<<"uuLoS ="<<uuL_L3D/uuL_Sigma3D<<" Ed:"<<uuLoS_new3d<<std::endl;

              if( MyPrint_ ){ 
                std::cout<<"Dimuon Vertex: Vx= "<<uuVx <<" Vy= "<<uuVy <<" Vz= "<<uuVz 
                         <<" Dist="<<uuDist <<" R="<<uuR<<std::endl;
                std::cout<<"Lambda Vertex: Vx= "<<lamVMx<<" Vy= "<<lamVMy<<" Vz= "<<lamVMz
                         <<" Dist="<<lamDist<<" R="<<lamR<<std::endl;
                std::cout<<"------------------"<<std::endl;
              } 

              //-->KinematicFitDriver Lambda_bRec( Lambda_bParticle, "Lambda_b" ); //--Vertex Lambda_b (Lambda + dimuon)
              //-->if( Lambda_bRec.isValid() ){

              //--Output Lambda_b Vtx fit variables
              //-->//double Mlam_bVFit = Lambda_bRec.mass(); 
              //-->//double lam_bV_CL  = ChiSquaredProbability( Lambda_bRec.chi2(), Lambda_bRec.ndof() );
              //-->//reco::Particle::LorentzVector lam_bV_P4 = Lambda_bRec.P4();  //--Lambda_b 4-momentum
              //-->//reco::Particle::Point lam_bV_vx = Lambda_bRec.VertexDecay();   //--Vx,Vy,Vz
              //-->//reco::Vertex::CovarianceMatrix lam_bV_Cov = Lambda_bRec.VertexCov();
              //-->//double lam_bVx  = lam_bV_vx.x();
              //-->//double lam_bVy  = lam_bV_vx.y();
              //-->//double lam_bVz  = lam_bV_vx.z();
              //-->//double lambDist = sqrt( lam_bVx*lam_bVx + lam_bVy*lam_bVy +lam_bVz*lam_bVz );
              //-->//double lambR    = sqrt( lam_bVx*lam_bVx + lam_bVy*lam_bVy );

              double lam_bV_CL =  1;
              double lambDist  = -1;
              double lambR     = -1;
              double Mlam_bVFit = ( uuRec.P4() + LamRec.P4() ).M();                   //--Mass Constrains Fitted Momentum for Both (uu & Lambda0)
              //double Mlam_buFit = ( iMuV_P4->P4() + jMuV_P4->P4() + LamRec.P4() ).M();  //--unFitted Mass for uu(unfitted) + Lambda0(fitted)
              
              double MuuK0s = ( uuRec.P4() + K0s ).M(); 
              
              histo_dimuM[8]->Fill( Muu );
              histo_lam[8]  ->Fill( MlamVFit );
              histo_lamb[8] ->Fill( Mlam_bVFit );

              //--Primary Vertex (excluding dimuon and Lambda tracks)
              std::vector<reco::TrackRef> ExclusionList;      
              ExclusionList.push_back( pionveetrk );       
              ExclusionList.push_back( protonveetrk );
              ExclusionList.push_back( iMuTrackRef );
              ExclusionList.push_back( jMuTrackRef ); 
              VertexRefit RefitPrimaryHM(primary.BestVertex() ,ExclusionList ,bFieldHandle ,beamSpot );    //--Highest Multiplicity Primary 
              reco::Vertex refitVertexPrimHM = RefitPrimaryHM.Refitted(); 
              double primCLHM   = ChiSquaredProbability( refitVertexPrimHM.chi2(), refitVertexPrimHM.ndof());       

              GlobalPoint  secver = uuRec.RKVertex()->position();
            //-GlobalVector secmom = uuRec.RKParent()->currentState().globalMomentum();       //--only dimuon       
              GlobalVector secmomTot = GlobalVector(uuRec.P4().Px()+LamRec.P4().Px()          //--this include dimuon and lambda_0
                                                   ,uuRec.P4().Py()+LamRec.P4().Py()
                                                   ,uuRec.P4().Pz()+LamRec.P4().Pz());
              if( MyPrint_ ){
                //std::cout<<"Lambda_b: vector P="<<secmom<<std::endl;
                  std::cout<<"        : uuRec   ="<<uuRec.P4().Px()<<std::endl;
                  std::cout<<"        : ---->   ="<<secmomTot<<std::endl;
              }
            //-VertexRefit RefitPrimary(primary.HigherCosAlpha(secver,secmom) ,ExclusionList ,bFieldHandle ,beamSpot ); //--pointing
              VertexRefit RefitPrimary(primary.HigherCosAlpha(secver,secmomTot) ,ExclusionList ,bFieldHandle ,beamSpot ); //--pointing
              reco::Vertex refitVertexPrim = RefitPrimary.Refitted();
              double primCL   = ChiSquaredProbability( refitVertexPrim.chi2(), refitVertexPrim.ndof());
                 
              //--std::cout<<" CL: High Multiplicity ="<<primCLHM<<" Pointing Angle ="<<primCL<<" ->"<<primCLHM-primCL<<std::endl;
              
              //--GlobalPoint PosVertex = GlobalPoint(refitVertexPrim.x(), refitVertexPrim.y(), refitVertexPrim.z() );
              double primDist = sqrt( refitVertexPrim.x()*refitVertexPrim.x()
                                    + refitVertexPrim.y()*refitVertexPrim.y()
                                    + refitVertexPrim.z()*refitVertexPrim.z() );
              double primR    = sqrt( refitVertexPrim.x()*refitVertexPrim.x()
                                    + refitVertexPrim.y()*refitVertexPrim.y() );
 
              double primDz  = refitVertexPrim.z() - refitVertexPrimHM.z();
              double primDxy = sqrt( (refitVertexPrim.x() - refitVertexPrimHM.x())*(refitVertexPrim.x() - refitVertexPrimHM.x())
                                   + (refitVertexPrim.y() - refitVertexPrimHM.y())*(refitVertexPrim.y() - refitVertexPrimHM.y()) );
              double primDxyz= sqrt( (refitVertexPrim.x() - refitVertexPrimHM.x())*(refitVertexPrim.x() - refitVertexPrimHM.x())
                                   + (refitVertexPrim.y() - refitVertexPrimHM.y())*(refitVertexPrim.y() - refitVertexPrimHM.y())
                                   + (refitVertexPrim.z() - refitVertexPrimHM.z())*(refitVertexPrim.z() - refitVertexPrimHM.z()) );
              
              if( MyPrint_ ){
                  //std::cout<<"Lambda_b: CL="<<lam_bV_CL<<std::endl; 
                  //std::cout<<"Lambda_b: Vtx = ("<<lam_bVx<<","<<lam_bVy<<","<<lam_bVz<<")"<<std::endl;
                  //std::cout<<"Lambda_b: P = ("<<lam_bV_px<<","<<lam_bV_py<<","<<lam_bV_pz<<")"<<std::endl;
                  std::cout<<"Lambda_b: primCL ="<<primCL<<" High Mult: "<<primCLHM<<std::endl;
                  std::cout<<"Lambda_b: primDz ="<<primDz<<" primDxy ="<<primDxy<<" primDxyz ="<<primDxyz<<std::endl;
              }

              //--Calculate Lambda_b Mass, L/Sigma for Primary, L/Sigma Dimuon-Lambda

              //--Check separation Lambda_b from primary (L over Sigma)
              typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D_Lambda_b;
              typedef ROOT::Math::SVector<double, 3> SVector3_Lambda_b;

              //-->SMatrixSym3D_Lambda_b totalCov_Lambda_b = refitVertexPrim.covariance() + lam_bV_Cov;
              //-->SVector3_Lambda_b Lambda_bPrimarySep3D( lam_bV_vx.x() - refitVertexPrim.x(),
              //-->                                        lam_bV_vx.y() - refitVertexPrim.y(),
              //-->                                        lam_bV_vx.z() - refitVertexPrim.z() );
              
              //--Take the dimuon vertex as the Secondary and measure L/Sigma from there (Dec 27, 11)
              SMatrixSym3D_Lambda_b totalCov_Lambda_b = refitVertexPrim.covariance() + uuV_Cov;
              SVector3_Lambda_b Lambda_bPrimarySep3D( uuV_vx.x() - refitVertexPrim.x(),
                                                      uuV_vx.y() - refitVertexPrim.y(),
                                                      uuV_vx.z() - refitVertexPrim.z() );

              double Lambda_bL3D     = ROOT::Math::Mag( Lambda_bPrimarySep3D );
              double Lambda_bSigma3D = sqrt(ROOT::Math::Similarity( totalCov_Lambda_b, Lambda_bPrimarySep3D )) / Lambda_bL3D;
              double Lambda_bLoS3D   = Lambda_bL3D / Lambda_bSigma3D;

              //--Here Pointing Secondary Vertex  (dimuon + lambda0) to Primary  [Sep 9, 2011]

              //-->GlobalVector Lambda_b_PrimSec( lam_bV_vx.x() - refitVertexPrim.x(),
              //-->                               lam_bV_vx.y() - refitVertexPrim.y(),
              //-->                               lam_bV_vx.z() - refitVertexPrim.z() );                            
              //-->GlobalVector Lambda_b_P3( lam_bV_P4.px(), lam_bV_P4.py(), lam_bV_P4.pz() );

              GlobalVector Lambda_b_PrimSec( uuV_vx.x() - refitVertexPrim.x(),
                                             uuV_vx.y() - refitVertexPrim.y(),
                                             uuV_vx.z() - refitVertexPrim.z() );                            
              GlobalVector Lambda_b_P3( (uuV_P4+lamVM_P4).px(),
                                        (uuV_P4+lamVM_P4).py(),
                                        (uuV_P4+lamVM_P4).pz() );
              double V_dot_Plamb = Lambda_b_P3.dot( Lambda_b_PrimSec );

              double cos_alphab3D = V_dot_Plamb /( Lambda_b_PrimSec.mag() * Lambda_b_P3.mag() );
              
              //--Cosine in 2D
              double cos_alphab2D = ( (uuV_vx.x() - refitVertexPrim.x())*(uuV_P4+lamVM_P4).px() + 
                                      (uuV_vx.y() - refitVertexPrim.y())*(uuV_P4+lamVM_P4).py() )
                               /sqrt( (uuV_vx.x() - refitVertexPrim.x())*(uuV_vx.x() - refitVertexPrim.x()) +
                                      (uuV_vx.y() - refitVertexPrim.y())*(uuV_vx.y() - refitVertexPrim.y()) )
                               /sqrt( (uuV_P4+lamVM_P4).px()*(uuV_P4+lamVM_P4).px() +
                                      (uuV_P4+lamVM_P4).py()*(uuV_P4+lamVM_P4).py() );
                                      
              //--NEW [FEB 6, 2012]
              double LoS_new3d   = uuRec.SignificanceSeparation( refitVertexPrim );
              double L_new3d     = uuRec.Separation( refitVertexPrim );
              double LoS_new2d   = uuRec.SignificanceSeparation2D( refitVertexPrim );
              double L_new2d     = uuRec.Separation2D( refitVertexPrim );           
                                      
              if( MyPrint_ ){
                  std::cout<<"Hec: Pointing Xi_b cos alpha= "<<cos_alphab3D<<" 2D:"<<cos_alphab2D<<std::endl;
                  std::cout<<"Hec: L/S ="<<Lambda_bLoS3D<<std::endl;  
                  std::cout<<"Hec: L/S= "<<Lambda_bLoS3D<<" Ed:"<<LoS_new3d<<std::endl;     
                  std::cout<<"Hec: L  = "<<Lambda_bL3D  <<" Ed:"<<L_new3d  <<std::endl;    
                  std::cout<<"2D: L/S= "<<LoS_new2d<<" L = "<<L_new2d  <<std::endl;
              }              
              if( Lambda_bLoS3D>3 && cos_alphab3D>0.95 && uu1sMcut && lamMcut ){
                histo_dimuM[9]->Fill( Muu );
                histo_lam[9]  ->Fill( MlamVFit );
                histo_lamb[9] ->Fill( Mlam_bVFit );
              }
              
              fill_lambdab(Mlam_bVFit, lam_bV_CL, lambDist, lambR, primCL, primDist, primR,
                           Lambda_bL3D, Lambda_bSigma3D, cos_alphab3D, cos_alphab2D,
                           L_new2d, LoS_new2d, MuuK0s, primDxy, primDxyz );
                            
              //--Fill Ntuple and re-initialize all variables
              iComb++;
              fill_IuuC(iQ, imuId, iNtrkHits, jQ, jmuId, jNtrkHits );
              fill_uuC(iEta,iPhi,iPx,iPy,iPz,id0,idz,iCal,iSeg,iIso,
                       jEta,jPhi,jPx,jPy,jPz,jd0,jdz,jCal,jSeg,jIso,
                       Muu, dca_uu, cxPt.x(), cxPt.y(), cxPt.z(), uuDist, uuR,
                       MuuVFit, uuV_CL,uuV_P4.Px(),uuV_P4.Py(),uuV_P4.Pz(),
                       MuuMFit, uuVM_CL, cos_uul3D, cos_uul2D, uuL_L3D, uuL_Sigma3D,
                       closestTRK, closestD0, closestSigmaDz, closestDz, closestTRK3d,closestId);

              //--Trigger from Keith [Apr 2, 12]
              const pat::Muon *patMuonP = &(*iMuon);
              const pat::Muon *patMuonM = &(*jMuon);
              if( iQ<0 ){
                  patMuonM = &(*iMuon);
                  patMuonP = &(*jMuon);
              }                            
              //--unpacks J/Psi trigger bit
              int upTrig2mu3=!patMuonP->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered3").empty();
              int umTrig2mu3=!patMuonM->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered3").empty();
              
              int upTrig2mu3JPsi=!patMuonP->triggerObjectMatchesByFilter("hltDoubleMu3JpsiL3Filtered").empty();
              int umTrig2mu3JPsi=!patMuonM->triggerObjectMatchesByFilter("hltDoubleMu3JpsiL3Filtered").empty();
              
              int upTrig2mu6p5JPsiDisp=!patMuonP->triggerObjectMatchesByFilter("hltDimuon6p5JpsiDisplacedL3Filtered").empty();
              int umTrig2mu6p5JPsiDisp=!patMuonM->triggerObjectMatchesByFilter("hltDimuon6p5JpsiDisplacedL3Filtered").empty();
                            
              int upTrig2mu7JPsiDisp=!patMuonP->triggerObjectMatchesByFilter("hltDisplacedmumuFilterJpsi").empty(); // also the same for HLT_DoubleMu3p5_Jpsi_Displaced_v2
              int umTrig2mu7JPsiDisp=!patMuonM->triggerObjectMatchesByFilter("hltDisplacedmumuFilterJpsi").empty(); // also the same for HLT_DoubleMu3p5_Jpsi_Displaced_v2
              
              int upTrig2mu4JPsiDisp=!patMuonP->triggerObjectMatchesByFilter("hltDisplacedmumuFilterDoubleMu4Jpsi").empty();
              int umTrig2mu4JPsiDisp=!patMuonM->triggerObjectMatchesByFilter("hltDisplacedmumuFilterDoubleMu4Jpsi").empty();             

              int upTrig2mu6p5JPsiPro=!patMuonP->triggerObjectMatchesByFilter("hltDimuon6p5JpsiL3Filtered").empty();
              int umTrig2mu6p5JPsiPro=!patMuonM->triggerObjectMatchesByFilter("hltDimuon6p5JpsiL3Filtered").empty();

              int upTrig2mu0JPsiPro=!patMuonP->triggerObjectMatchesByFilter("hltVertexmumuFilterJpsi").empty();
              int umTrig2mu0JPsiPro=!patMuonM->triggerObjectMatchesByFilter("hltVertexmumuFilterJpsi").empty();
 
              int upTrig2mu0JPsiNoVtxPro=!patMuonP->triggerObjectMatchesByFilter("hltJpsiNoVertexingL3Filtered").empty();
              int umTrig2mu0JPsiNoVtxPro=!patMuonM->triggerObjectMatchesByFilter("hltJpsiNoVertexingL3Filtered").empty();

              if( MyPrint_ ){
                  std::cout << "Mu Trigger mu3 "         <<upTrig2mu3            << "  " <<umTrig2mu3 <<std::endl;
                  std::cout << "Mu TriggerJPsi "         <<upTrig2mu3JPsi        << "  " <<umTrig2mu3JPsi <<std::endl;
                  std::cout << "Mu TriggerJPsiDisp "     <<upTrig2mu6p5JPsiDisp  << "  " <<umTrig2mu6p5JPsiDisp <<std::endl;
                  std::cout << "Mu Trigger7JPsiDisp "    <<upTrig2mu7JPsiDisp    << "  " <<umTrig2mu7JPsiDisp <<std::endl;
                  std::cout << "Mu Trigger4JPsiDisp "    <<upTrig2mu4JPsiDisp    << "  " <<upTrig2mu4JPsiDisp <<std::endl;
                  std::cout << "Mu Trigger5JPsiDisp "    <<upTrig2mu6p5JPsiPro   << "  " <<umTrig2mu6p5JPsiPro <<std::endl;
                  std::cout << "Mu Trigger0JPsiDisp "    <<upTrig2mu0JPsiPro     << "  " <<umTrig2mu0JPsiPro <<std::endl;
                  std::cout << "Mu TriggerNoVtxJPsiDisp "<<upTrig2mu0JPsiNoVtxPro<< "  " <<umTrig2mu0JPsiNoVtxPro <<std::endl;
              }
              
              int upTrig  = upTrig2mu3             *1
                          + upTrig2mu3JPsi         *10
                          + upTrig2mu6p5JPsiDisp   *100
                          + upTrig2mu7JPsiDisp     *1000
                          + upTrig2mu4JPsiDisp     *10000
                          + upTrig2mu6p5JPsiPro    *100000
                          + upTrig2mu0JPsiPro      *1000000
                          + upTrig2mu0JPsiNoVtxPro *10000000;
              int umTrig  = umTrig2mu3             *1
                          + umTrig2mu3JPsi         *10
                          + umTrig2mu6p5JPsiDisp   *100
                          + umTrig2mu7JPsiDisp     *1000
                          + umTrig2mu4JPsiDisp     *10000
                          + umTrig2mu6p5JPsiPro    *100000
                          + umTrig2mu0JPsiPro      *1000000
                          + umTrig2mu0JPsiNoVtxPro *10000000;
                          
              //--Unpacks Displaced LowMass_Dimuon trigger bit
              int umLMDimuon6p5=!patMuonM->triggerObjectMatchesByPath("HLT_Dimuon6p5_LowMass_Displaced_v*").empty();
              int upLMDimuon6p5=!patMuonP->triggerObjectMatchesByPath("HLT_Dimuon6p5_LowMass_Displaced_v*").empty();
              
              int umLMDimuon7  =!patMuonM->triggerObjectMatchesByPath("HLT_Dimuon7_LowMass_Displaced_v*").empty();
              int upLMDimuon7  =!patMuonP->triggerObjectMatchesByPath("HLT_Dimuon7_LowMass_Displaced_v*").empty();
              
              int umLMDoubleMu4=!patMuonM->triggerObjectMatchesByPath("HLT_DoubleMu4_LowMass_Displaced_v*").empty();
              int upLMDoubleMu4=!patMuonP->triggerObjectMatchesByPath("HLT_DoubleMu4_LowMass_Displaced_v*").empty();
              
              int umLMDoubleMu4p5=!patMuonM->triggerObjectMatchesByPath("HLT_DoubleMu4p5_LowMass_Displaced_v*").empty();
              int upLMDoubleMu4p5=!patMuonP->triggerObjectMatchesByPath("HLT_DoubleMu4p5_LowMass_Displaced_v*").empty();
              
              int umLMDoubleMu5=!patMuonM->triggerObjectMatchesByPath("HLT_DoubleMu5_LowMass_Displaced_v*").empty();             
              int upLMDoubleMu5=!patMuonP->triggerObjectMatchesByPath("HLT_DoubleMu5_LowMass_Displaced_v*").empty();
              
              if( MyPrint_ ){
                  std::cout << "Dimuon6p5 " <<  umLMDimuon6p5 << "  " <<  upLMDimuon6p5 <<std::endl;
                  std::cout << "Dimuon7 " <<  umLMDimuon7 << "  " <<  upLMDimuon7 <<std::endl;
                  std::cout << "DoubleMu4 " <<  umLMDoubleMu4 << "  " <<  upLMDoubleMu4 <<std::endl;
                  std::cout << "DoubleMu4p5 " <<  umLMDoubleMu4p5 << "  " <<  upLMDoubleMu4p5 <<std::endl;
                  std::cout << "DoubleMu5 " <<  umLMDoubleMu5 << "  " <<  upLMDoubleMu5 <<std::endl;         
                  std::cout << "-----------------------------"<<std::endl;
              }
              int muLMtrig = umLMDimuon6p5   *1
                           + upLMDimuon6p5   *10
                           + umLMDimuon7     *100
                           + upLMDimuon7     *1000
                           + umLMDoubleMu4   *10000
                           + upLMDoubleMu4   *100000
                           + umLMDoubleMu4p5 *1000000
                           + upLMDoubleMu4p5 *10000000     
                           + umLMDoubleMu5   *100000000            
                           + upLMDoubleMu5   *1000000000;
 
              
              
              fill_evt( iEvent, allTrk, allMu, allPrim, allLambda, iComb,
                        upTrig, umTrig, muLMtrig, evtTrig1, evtTrig2, evtTrig3, evtTrig4 );
              
              hecmu_tree_->Fill();

              //mupTrig2mu3->clear();
              //}   //--if( Lambda_bRec.isValid() ){ //--Dimuon+Lambda Vertex
              //}   //--if( uuV_CL > 0.01 ){   //--dimuon CL
              }   //--if( cos_uul3D > 0.75 ){
              }   //--else{  valid Dimuon Vertex
              }   //--Fiducial cut for dimuon
              }   //--no repeated tracks between muons and Lambda
            }     //--Opp charge & Muon Id       
         }        //--end muonCollection for j              
      }           //--end muonCollection for i
   }              //--end if(  LamRec.isValid()   Good Lambda Mass Constrain
   }              //--end if(  LamRec.isValid()   Good Lambda Vertex
}                 //--end for( unsigned int veeindex )    

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
HecPsiLambdaPAT::beginJob()
{  
  std::cout << " HecPsiLambdaPAT::beginJob" << std::endl;
  hecmu_tree_ = new TTree("Lambda","Dettach Ntuple");
  int bufsize = 64000;

  hecmu_tree_->Branch("Dimuon",  &uuC_,    "iEta/D:iPhi:iPx:iPy:iPz:id0:idz:iCal:iSeg:iIso:jEta:jPhi:jPx:jPy:jPz:jd0:jdz:jCal:jSeg:jIso:Muu:dca:xcpt:ycpt:zcpt:uuDist:uuR:MuuVF:uuVCL:uuVMpx:uuVMpy:uuVMpz:MuuMFit:uuVMCL:cosuul3D:cosuul2D:uuLL3D:uuLSigma3D:closestTRK:closestD0:closestSigmaDz:closestDz:closestTRK3d:closestId", bufsize);
  hecmu_tree_->Branch("Lambda0", &lambdaC_,"prQ/D:prEta:prPhi:prPx:prPy:prPz:prd0:prdz:prNHits:piQ:piEta:piPhi:piPx:piPy:piPz:pid0:pidz:piNHits:Mlambda0:lamCL:Mlambda0VF:lamDist:lamR:Mkshort:MlamMFit:lamVMCL:lamVMpx:lamVMpy:lamVMpz", bufsize);  
  hecmu_tree_->Branch("Lambdab", &lambdabC_,"lambCL/D:MlambVF:lambDist:lambR:primCL:primDist:primR:lambdabL3D:lambdabSigma3D:cosAlphab3D:cosAlphab2D:L_2D:LoS_2D:MuuK0s:primDxy:primDxyz", bufsize);  
  hecmu_tree_->Branch("TightMu", &tightMuC_,"iNormChi2/D:iMuMuHits:iMuStations:iMuPxHits:iMuTrkLayer:iKink:jNormChi2:jMuMuHits:jMuStations:jMuPxHits:jMuTrkLayer:jKink", bufsize);  
 
  hecmu_tree_->Branch("iHeader", &Ievt_, "run/I:evtnum:lumib:allTrk:allMu:allPrim:allLambda:iComb:upTrig:umTrig:muLMtrig:evtTrig1:evtTrig2:evtTrig3:evtTrig4", bufsize);
  hecmu_tree_->Branch("iDimuon", &IuuC_, "iQ/I:imuId:iNtrkHits:jQ:jmuId:jNtrkHits", bufsize);

}

// ------------ method called once each job just after ending the event loop  ------------
void 
HecPsiLambdaPAT::endJob() 
{
  std::cout << " HecPsiLambdaPAT::endJob" << std::endl; 
  //-------------------Maria--------------------
  std::cout << " --My Final MC number-- " << std::endl;
  std::cout << " --JPsi_Lamb0P---  "<< Lb[1] 
            << " --JPsi_Lamb0A---  "<< Lb[2] 
            << " --JPsi_Lamb0PA--- "<< Lb[3] <<std::endl;

  std::cout << " --Psi1_Lamb0P---  "<< Lb[4]
            << " --Psi1_Lamb0A---  "<< Lb[5] 
            << " --Psi1_Lamb0PA--- "<< Lb[6] <<std::endl;
            
  std::cout << " --Psi2_Lamb0P---  "<< Lb[7] 
            << " --Psi2_Lamb0A---  "<< Lb[8] 
            << " --Psi2_Lamb0PA--- "<< Lb[9] <<std::endl;
  std::cout  <<std::endl;
  double All_Lambda_b = Fb[2] + Fb[3];
  std::cout << " --All_Lamb0P---  "<< Fb[0] 
            << " --All_Lamb0A---  "<< Fb[1] <<std::endl;
  std::cout << " --Can_Lamb0P---  "<< Fb[2] 
            << " --Can_Lamb0A---  "<< Fb[3] <<" All Lambda_b = "<<All_Lambda_b<<std::endl; 
  std::cout  <<std::endl;
            
  std::cout << " --MC Rates numbers-- " <<std::endl;
  double JpsiLambdabRatio   = 0;    //--J/Psi/Lambda_b
  double Psi1LambdabRatio   = 0;    //--Psi(2s)/Lambda_b for Psi(2s)-->u+u-
  double Psi2LambdabRatio   = 0;    //--Psi(2s)/Lambda_b for Psi(2s)-->u+u-Pi+Pi-
  double Psi1JpsiRatio = 0;         //--Psi(2s)-->u+u-
  double Psi2JpsiRatio = 0;         //--Psi(2s)-->Psi(1s)Pi+Pi-
  double Psi1Psi2Ratio = 0;         //--Psi(2s)-->u+u- / Psi(2s)-->Psi(1s)Pi+Pi-
  //--error of f=a/b is sigma_f = f sqrt[ (f+1)/a ]  (if a->0 => sigma_f ->inf, if a->inf => sigma_f -> 0)
  double sJpsiLambdabRatio   = 0;    //--J/Psi/Lambda_b
  double sPsi1LambdabRatio   = 0;    //--Psi(2s)/Lambda_b for Psi(2s)-->u+u-
  double sPsi2LambdabRatio   = 0;    //--Psi(2s)/Lambda_b for Psi(2s)-->u+u-Pi+Pi-
  double sPsi1JpsiRatio = 0;         //--Psi(2s)-->u+u-
  double sPsi2JpsiRatio = 0;         //--Psi(2s)-->Psi(1s)Pi+Pi-
  double sPsi1Psi2Ratio = 0;         //--Psi(2s)-->u+u- / Psi(2s)-->Psi(1s)Pi+Pi-
  
  if( All_Lambda_b!= 0 ){
      JpsiLambdabRatio = double( Lb[3] )/ All_Lambda_b;
      Psi1LambdabRatio = double( Lb[6] )/ All_Lambda_b;
      Psi2LambdabRatio = double( Lb[9] )/ All_Lambda_b;
      
      sJpsiLambdabRatio = JpsiLambdabRatio*sqrt( (JpsiLambdabRatio+1)/ All_Lambda_b);
      sPsi1LambdabRatio = Psi1LambdabRatio*sqrt( (Psi1LambdabRatio+1)/ All_Lambda_b);
      sPsi2LambdabRatio = Psi2LambdabRatio*sqrt( (Psi2LambdabRatio+1)/ All_Lambda_b);
  }
  if( Lb[3] != 0 ){
      Psi1JpsiRatio = double( Lb[6] )/ double( Lb[3] );
      Psi2JpsiRatio = double( Lb[9] )/ double( Lb[3] );
      
      sPsi1JpsiRatio = Psi1JpsiRatio*sqrt( (Psi1JpsiRatio+1)/ double( Lb[3] ));
      sPsi2JpsiRatio = Psi2JpsiRatio*sqrt( (Psi2JpsiRatio+1)/ double( Lb[3] ));
  }
  if( Lb[9] != 0 ){
      Psi1Psi2Ratio  = double( Lb[6] )/ double( Lb[9] );
      sPsi1Psi2Ratio = Psi1Psi2Ratio*sqrt(( Psi1Psi2Ratio+1)/ double( Lb[9] ));
  }
  std::cout << "Lb->L0Jpsi[p+pi-u+u-]/Lb:                  JpsiLambdabRatio = "<< JpsiLambdabRatio <<" +/- "<<sJpsiLambdabRatio<<std::endl;
  std::cout << "Lb->L0Psi(2s)[p+pi-u+u-]/Lambda_b:         JpsiLambdabRatio = "<< Psi1LambdabRatio <<" +/- "<<sPsi1LambdabRatio<<std::endl;
  std::cout << "Lb->L0Psi(2s)[p+pi-Jpsi pi+pi-]/Lambda_b:  JpsiLambdabRatio = "<< Psi2LambdabRatio <<" +/- "<<sPsi2LambdabRatio<<std::endl;
  std::cout <<std::endl;
  std::cout << "Lb->L0Psi(2s)[u+u-]/Lb->L0Jpsi:                Psi1JpsiRatio = "<< Psi1JpsiRatio<<" +/- "<<sPsi1JpsiRatio <<"  Hec: Expected: 0.055"<<std::endl;
  std::cout << "Lb->L0Psi(2s)[Psi(1s)Pi+Pi-]/Lb->L0Jpsi:       Psi2JpsiRatio = "<< Psi2JpsiRatio<<" +/- "<<sPsi2JpsiRatio <<"  Hec: Expected: 0.150"<<std::endl;
  std::cout << "Lb->L0Psi(2s)[u+u-]/Lb->L0Psi(2s)[JpsiPi+Pi-]: Psi1psi2Ratio = "<< Psi1Psi2Ratio<<" +/- "<<sPsi1Psi2Ratio <<"  Hec: Expected: 0.366"<<std::endl;
  std::cout <<std::endl;
  std::cout <<"----MC Br Table----"<<std::endl;
  std::cout <<"Lambda_b->Lambda J/psi->(p pi) (mu mu):       0.168*1     *0.639 = 10.74%"<<std::endl;
  std::cout <<"Lambda_b->Lambda psi2s->(p pi) (mu mu):       0.053*0.1741*0.639 =  0.59%"<<std::endl;
  std::cout <<"Lambda_b->Lambda psi2s->(p pi) (mu mu pi pi): 0.053*0.4762*0.639 =  1.62%"<<std::endl;      
}

// ------------ method called when starting to processes a run  ------------
void 
HecPsiLambdaPAT::beginRun(edm::Run const&, edm::EventSetup const&)
{
  std::cout << " HecPsiLambdaPAT::beginRun" << std::endl;
}

// ------------ method called when ending the processing of a run  ------------
void 
HecPsiLambdaPAT::endRun(edm::Run const&, edm::EventSetup const&)
{
  std::cout << " HecPsiLambdaPAT::endRun" << std::endl;
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
HecPsiLambdaPAT::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
  std::cout << " HecPsiLambdaPAT::beginLuminosityBlock" << std::endl;
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
HecPsiLambdaPAT::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
  std::cout << " HecPsiLambdaPAT::endLuminosityBlock" << std::endl;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HecPsiLambdaPAT::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

 //Specify that only 'tracks' is allowed
 //To use, remove the default given above and uncomment below
 //ParameterSetDescription desc;
 //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
 //descriptions.addDefault(desc);
}
//--from here below are my functions
void HecPsiLambdaPAT::init()
{
  Ievt_.init();
  IuuC_.init();
  uuC_.init();
  lambdaC_.init();
  lambdabC_.init();
  tightMuC_.init();
}
void HecPsiLambdaPAT::Ievt::init()
{
  runNb      = -1; eventNb    = -1; lumiBlock  = -1;
  allTrk     = -1; allMu      = -1; allPrim    = -1; allLambda  = -1;
  iComb      =  0; upTrig     = -1; umTrig     = -1; muLMtrig   = -1;
  evtTrig1   =  0; evtTrig2   =  0; evtTrig3   =  0; evtTrig4   =  0;
}
void HecPsiLambdaPAT::fill_evt(const edm::Event& iEvent, int allTrk, int allMu, int allPrim, int allLambda,
                              int iComb, int upTrig, int umTrig, int muLMtrig,
                              int evtTrig1, int evtTrig2, int evtTrig3, int evtTrig4)
{
     Ievt_.runNb      = iEvent.id().run(); 
     Ievt_.eventNb    = iEvent.id().event();
     Ievt_.lumiBlock  = iEvent.luminosityBlock();
     Ievt_.allTrk     = allTrk;         //--tracks->size();
     Ievt_.allMu      = allMu;          //--->allmuons->size();
     Ievt_.allPrim    = allPrim;
     Ievt_.allLambda  = allLambda;
     Ievt_.iComb      = iComb;
     Ievt_.upTrig     = upTrig;
     Ievt_.umTrig     = umTrig;
     Ievt_.muLMtrig   = muLMtrig;
     Ievt_.evtTrig1   = evtTrig1;
     Ievt_.evtTrig2   = evtTrig2;
     Ievt_.evtTrig3   = evtTrig3;
     Ievt_.evtTrig4   = evtTrig4;
     if( MyPrint_ )
       std::cout<<"Hec: Run: "<<Ievt_.runNb<<" Evt: "<<Ievt_.eventNb<<" Lumi:"<<Ievt_.lumiBlock
                <<" Ntk: "<<Ievt_.allTrk<<" Mu:"<<Ievt_.allMu
                <<" NPrimVtx:"<<Ievt_.allPrim<<" #Lambda:"<<allLambda
                <<std::endl;
}
void HecPsiLambdaPAT::IuuC::init()
{
  iQ = 0; imuId = -1; iNtrkHits = -1;
  jQ = 0; jmuId = -1; jNtrkHits = -1;
}
void HecPsiLambdaPAT::fill_IuuC(int iQ, int imuId, int iNtrkHits,
                             int jQ, int jmuId, int jNtrkHits)
{
  IuuC_.iQ = iQ; IuuC_.imuId = imuId; IuuC_.iNtrkHits = iNtrkHits;
  IuuC_.jQ = jQ; IuuC_.jmuId = jmuId; IuuC_.jNtrkHits = jNtrkHits;
}
void HecPsiLambdaPAT::uuC::init()
{
  iEta =9.9; iPhi =-1.0; iPx =-999; iPy =-999; iPz=-999; id0 =-1.0; idz =-1.0; iCal =-1.0; iSeg =-1.0; iIso =-1.0;
  jEta =9.9; jPhi =-1.0; jPx =-999; jPy =-999; jPz=-999; jd0 =-1.0; jdz =-1.0; jCal =-1.0; jSeg =-1.0; jIso =-1.0;
  Muu = -1.0; dca = -1.0; xcpt = -1.0; ycpt = -1.0; zcpt = -1.0; uuDist = -1.0; uuR = -1.0; 
  MuuVF=-1.0; uuVCL=-1.0; uuVMpx=-1;uuVMpy=-1;uuVMpz=-1;
  MuuMFit=-1; uuVMCL=-1; cosuul3D=-999; cosuul2D=-999; uuLL3D=-1; uuLSigma3D=-1;
  closestTRK=99; closestD0=99; closestSigmaDz=99; closestDz=99; closestTRK3d=99; closestId=0;
}
void HecPsiLambdaPAT::fill_uuC(double iEta,double iPhi,double iPx,double iPy,double iPz,double id0,double idz,double iCal,double iSeg,double iIso,
                            double jEta,double jPhi,double jPx,double jPy,double jPz,double jd0,double jdz,double jCal,double jSeg,double jIso,
                            double Muu,double dca,double  xcpt,double ycpt,double  zcpt, double dist, double R2dim,
                            double MuuFit, double uuV_CL, double uuVMpx, double uuVMpy, double uuVMpz,
                            double MuuMFit, double uuVM_CL, double cos_uul3D, double cos_uul2D, double uuL_L3D, double uuL_Sigma3D,
                            double closestTRK, double closestD0, double closestSigmaDz, double closestDz, double closestTRK3d, double closestId)
{
  uuC_.iEta    = iEta;
  uuC_.iPhi    = iPhi;
  uuC_.iPx     = iPx ;
  uuC_.iPy     = iPy ;
  uuC_.iPz     = iPz ;
  uuC_.id0     = id0 ;
  uuC_.idz     = idz ;
  uuC_.iCal    = iCal;
  uuC_.iSeg    = iSeg;
  uuC_.iIso    = iIso;
  uuC_.jEta    = jEta;
  uuC_.jPhi    = jPhi;
  uuC_.jPx     = jPx ;
  uuC_.jPy     = jPy ;
  uuC_.jPz     = jPz ;
  uuC_.jd0     = jd0 ;
  uuC_.jdz     = jdz ;
  uuC_.jCal    = jCal;
  uuC_.jSeg    = jSeg;
  uuC_.jIso    = jIso;
  uuC_.Muu     = Muu ;
  uuC_.dca     = dca ;
  uuC_.xcpt    = xcpt;
  uuC_.ycpt    = ycpt;
  uuC_.zcpt    = zcpt;
  uuC_.uuDist  = dist;
  uuC_.uuR     = R2dim;
  uuC_.MuuVF   = MuuFit;
  uuC_.uuVCL   = uuV_CL;
  uuC_.uuVMpx  = uuVMpx;
  uuC_.uuVMpy  = uuVMpy;
  uuC_.uuVMpz  = uuVMpz;
  uuC_.MuuMFit    = MuuMFit;
  uuC_.uuVMCL     = uuVM_CL;
  uuC_.cosuul3D   = cos_uul3D;
  uuC_.cosuul2D   = cos_uul2D;
  uuC_.uuLL3D     = uuL_L3D;
  uuC_.uuLSigma3D = uuL_Sigma3D;
  uuC_.closestTRK     = closestTRK;
  uuC_.closestD0      = closestD0;
  uuC_.closestSigmaDz = closestSigmaDz;
  uuC_.closestDz      = closestDz;
  uuC_.closestTRK3d   = closestTRK3d;
  uuC_.closestId      = closestId;
}
void HecPsiLambdaPAT::lambdaC::init()
{
  prQ = 0; prEta = -10; prPhi = -10; prPx = -9999; prPy = -9999; prPz = -9999; prd0 = -1; prdz = -1; prNHits = -1;
  piQ = 0; piEta = -10; piPhi = -10; piPx = -9999; piPy = -9999; piPz = -9999; pid0 = -1; pidz = -1; piNHits = -1;
  Mlambda0 = -1; lamCL      = -1; Mlambda0VF = -1; lamDist    = -1; lamR       = -1; Mkshort = -1;
  MlamMFit = -1; lamVMCL   = -1; lamVMpx = -999; lamVMpy = -999; lamVMpz = -999;
}
void HecPsiLambdaPAT::fill_lambda0(double prQ, double prEta, double prPhi, double prPx, double prPy, double prPz,
                                double prd0, double prdz, double prNtrkHits,
                                double piQ, double piEta, double piPhi, double piPx, double piPy, double piPz,
                                double pid0, double pidz, double piNtrkHits,
                                double Mlambda0, double lamVCL, double Mlambda0Vfit, double dist, double R2dim, double myKshort,
                                double MlamMFit, double lamVM_CL, double lamVMpx, double lamVMpy, double lamVMpz)
{ 
   lambdaC_.prQ      = prQ;
   lambdaC_.prEta    = prEta;
   lambdaC_.prPhi    = prPhi;
   lambdaC_.prPx     = prPx; 
   lambdaC_.prPy     = prPy; 
   lambdaC_.prPz     = prPz; 
   lambdaC_.prd0     = prd0; 
   lambdaC_.prdz     = prdz;
   lambdaC_.prNHits  = prNtrkHits;
   lambdaC_.piQ      = piQ;
   lambdaC_.piEta    = piEta;
   lambdaC_.piPhi    = piPhi;
   lambdaC_.piPx     = piPx; 
   lambdaC_.piPy     = piPy; 
   lambdaC_.piPz     = piPz; 
   lambdaC_.pid0     = pid0; 
   lambdaC_.pidz     = pidz;
   lambdaC_.piNHits  = piNtrkHits;
   lambdaC_.Mlambda0   = Mlambda0;
   lambdaC_.lamCL      = lamVCL;
   lambdaC_.Mlambda0VF = Mlambda0Vfit;
   lambdaC_.lamDist    = dist;
   lambdaC_.lamR       = R2dim;
   lambdaC_.Mkshort    = myKshort;
   lambdaC_.MlamMFit   = MlamMFit;
   lambdaC_.lamVMCL    = lamVM_CL;
   lambdaC_.lamVMpx    = lamVMpx;
   lambdaC_.lamVMpy    = lamVMpy;
   lambdaC_.lamVMpz    = lamVMpz;
}
void HecPsiLambdaPAT::lambdabC::init()
{
  lambCL = -1; MlambVF = -1; lambDist = -1; lambR = -1;
  primCL = -1; primDist = -1; primR  = -1;
  lambdabL3D = -1; lambdabSigma3D = -1; cosalphab3D = -10; cosalphab2D = -10;
  L_2D = -1; LoS_2D = -1; MuuK0s = -1; primDxy=99; primDxyz=99;
}
void HecPsiLambdaPAT::fill_lambdab(double Mlam_bVFit, double lam_bV_CL, double lambDist, double lambR, 
                                double primCL, double primDist, double primR,
                                double Lambda_bL3D, double Lambda_bSigma3D, double cos_alphab3D, double cos_alphab2D,
                                double L_new2d, double LoS_new2d, double MuuK0s, double primDxy, double primDxyz)
{
  lambdabC_.lambCL      = lam_bV_CL;
  lambdabC_.MlambVF     = Mlam_bVFit;
  lambdabC_.lambDist    = lambDist;
  lambdabC_.lambR       = lambR;
  lambdabC_.primCL      = primCL;
  lambdabC_.primDist    = primDist;
  lambdabC_.primR       = primR;
  lambdabC_.lambdabL3D     = Lambda_bL3D;
  lambdabC_.lambdabSigma3D = Lambda_bSigma3D;
  lambdabC_.cosalphab3D    = cos_alphab3D;
  lambdabC_.cosalphab2D    = cos_alphab2D;
  lambdabC_.L_2D           = L_new2d;
  lambdabC_.LoS_2D         = LoS_new2d;
  lambdabC_.MuuK0s         = MuuK0s;
  lambdabC_.primDxy        = primDxy;
  lambdabC_.primDxyz       = primDxyz;
}
void HecPsiLambdaPAT::tightMuC::init()
{
  iNormChi2 = 0, iMuMuHits = -1, iMuStations = -1, iMuPxHits = -1, iMuTrkLayer = -1, iKink = 99;
  jNormChi2 = 0, jMuMuHits = -1, jMuStations = -1, jMuPxHits = -1, jMuTrkLayer = -1, jKink = 99;
}            
void HecPsiLambdaPAT::fill_tightMuons(double iNormChi2,double iMuMuHits,double iMuStations,double iMuPxHits,double iMuTrkLayer,double iKink
                                     ,double jNormChi2,double jMuMuHits,double jMuStations,double jMuPxHits,double jMuTrkLayer,double jKink)
{
 tightMuC_.iNormChi2     = iNormChi2;
 tightMuC_.iMuMuHits     = iMuMuHits;
 tightMuC_.iMuStations   = iMuStations;
 tightMuC_.iMuPxHits     = iMuPxHits;
 tightMuC_.iMuTrkLayer   = iMuTrkLayer;
 tightMuC_.iKink         = iKink;
 tightMuC_.jNormChi2     = jNormChi2;
 tightMuC_.jMuMuHits     = jMuMuHits;
 tightMuC_.jMuStations   = jMuStations;
 tightMuC_.jMuPxHits     = jMuPxHits;
 tightMuC_.jMuTrkLayer   = jMuTrkLayer;
 tightMuC_.jKink         = jKink;
}
//Method to compute 3D absolute impact parameter
//based on information of Transversal Impact Parameter (2D) D0 and
//Longitudinal Impact Parameter (along Z)
// IP3d = sqrt(D0^2 + Dz^2)
//The error is propagated adding in quadrature:
// IP3de^2 = (d(IP3d)/d(D0))^2 * D0e^2 + (d(IP3d)/d(Dz))^2 * Dze^2
//simplified to:
// IP3de = sqrt(D0e^2*D0^2 + Dze^2*Dz^2)/IP3d

double HecPsiLambdaPAT::SignificanceAbsoluteImpactParameter3D(
                               reco::TransientTrack &TTrack, reco::Vertex &Vertex,
                               double &TrakDzE, double &TrakDz,  
                               double &TrakD0E, double &TrakD0, double &trkip2d) const{
                               
   GlobalPoint PosVertex = GlobalPoint(Vertex.x(),Vertex.y(),Vertex.z());
   TrajectoryStateClosestToPoint TTraj = TTrack.trajectoryStateClosestToPoint(PosVertex);
          TrakD0  = TTraj.isValid()?TTraj.perigeeParameters().transverseImpactParameter()
                                   :-1000.;
          TrakD0E = TTraj.isValid()?TTraj.perigeeError().transverseImpactParameterError()
                                   :-1000.;
          TrakDz  = TTraj.isValid()?TTraj.perigeeParameters().longitudinalImpactParameter()
                                   :-1000.;
          TrakDzE = TTraj.isValid()?TTraj.perigeeError().longitudinalImpactParameterError()
                                   :-1000.;
   double TrakD3  = TTraj.isValid()?sqrt(TrakD0*TrakD0 + TrakDz*TrakDz):-10000.;
   double TrakD3E = TTraj.isValid()?sqrt(TrakD0E*TrakD0E*TrakD0*TrakD0
                                       + TrakDzE*TrakDzE*TrakDz*TrakDz)/TrakD3:-10000.;
          trkip2d = TrakD0E>0 ? TrakD0/TrakD0E :9999;    //--significance in 2d L/Sigma
        //trkip3d = TrakD3E>0 ? TrakD3/TrakD3E :9999;    //--significance in 3d L/Sigma

   return ( TrakD3E>0 ? TrakD3/TrakD3E :9999 );   
}
//define this as a plug-in
DEFINE_FWK_MODULE(HecPsiLambdaPAT);
