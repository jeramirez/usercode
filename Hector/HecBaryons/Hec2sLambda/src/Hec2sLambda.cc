// -*- C++ -*-
//
// Package:    Hec2sLambda
// Class:      Hec2sLambda
// 
/**\class Hec2sLambda Hec2sLambda.cc HecBaryons/Hec2sLambda/src/Hec2sLambda.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Hector Mendez
//         Created:  Sun Jan  8 11:56:28 CST 2012
// $Id: Hec2sLambda.cc,v 1.2 2012/06/04 22:11:56 mendez Exp $
//
//

// system include files
#include <memory>

// user include files
#include "HecBaryons/Hec2sLambda/interface/Hec2sLambda.h"    //--Hec

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

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

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"       //track   collection
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"//cascade colletcion
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"//cascade daughter track

#include "Analyzers/CascadeProducer/interface/PrimaryInfo.h"             //--Eduardo Cascade & Primary
#include "Analyzers/CascadeProducer/interface/VertexRefit.h"
#include "Analyzers/CascadeProducer/interface/KinematicFitDriver.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "Analyzers/CascadeProducer/interface/ClosestApproachOnHelixLine.h"

//--kinemactic fitter
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

//--Bfield
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "Analyzers/CascadeProducer/interface/masses.h"

#include "TH1.h"

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
Hec2sLambda::Hec2sLambda(const edm::ParameterSet& iConfig)
:
 trackTags_     (iConfig.getParameter<edm::InputTag>("tracks"    )),
 theMuonsLabel_ (iConfig.getParameter<edm::InputTag>("MuonsLabel")),
 VeeAlgo_       (iConfig.getUntrackedParameter<std::string>("VeeAlgo","generalV0Candidates")),
 MyPrint_       (iConfig.getUntrackedParameter<bool>("MyPrint","False"))
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
//histo_lam[5]   = fs->make<TH1D>("MassLam5","repTrk Op"    , 200, 1.07, 1.17);
//histo_lam[6]   = fs->make<TH1D>("MassLam6","repTrk Same"  , 200, 1.07, 1.17);
  histo_lam[7]   = fs->make<TH1D>("MassLam7","Mass uuVtx"   , 200, 1.07, 1.17);
  histo_lam[8]   = fs->make<TH1D>("MassLam8","Mass JLVtx"   , 200, 1.07, 1.17);
  histo_lam[9]   = fs->make<TH1D>("MassLam9","Mass Fin"     , 200, 1.07, 1.17);

//histo_dimuM[0] = fs->make<TH1D>("dimuM0","Muu All"      , 500, 0, 100);   //--200 MeV/channel
//histo_dimuM[1] = fs->make<TH1D>("dimuM1","Muu All"      , 500, 0, 100);
  histo_dimuM[2] = fs->make<TH1D>("dimuM2","Muu All Os"   , 500, 0, 100);
  histo_dimuM[3] = fs->make<TH1D>("dimuM3","Muu All Ss"   , 500, 0, 100);
  histo_dimuM[4] = fs->make<TH1D>("dimuM4","Muu All"      , 500, 0, 100);
//histo_dimuM[5] = fs->make<TH1D>("dimuM5","repTrk Op"    , 500, 0, 100);
//histo_dimuM[6] = fs->make<TH1D>("dimuM6","repTrk Same"  , 500, 0, 100);
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
  
  histo_mdiff[0] = fs->make<TH1D>("mdiff0","Muupipi - Muu 0", 1000, 0, 10);
  histo_mdiff[1] = fs->make<TH1D>("mdiff1","Muupipi - Muu 1", 1000, 0, 10);
  histo_mdiff[2] = fs->make<TH1D>("mdiff2","Muupipi - Muu 2", 1000, 0, 10);
  histo_mdiff[3] = fs->make<TH1D>("mdiff3","Muupipi - Muu 3", 1000, 0, 10);
  histo_mdiff[4] = fs->make<TH1D>("mdiff4","Muupipi - Muu 4", 1000, 0, 10);
  histo_mdiff[5] = fs->make<TH1D>("mdiff5","Muupipi - Muu 5", 1000, 0, 10);
  
  histo_mpipi[0] = fs->make<TH1D>("mpipi0","Mpipi 0", 200, 0, 5);
  histo_mpipi[1] = fs->make<TH1D>("mpipi1","Mpipi 1", 200, 0, 5);
  
  histo_vx[0] = fs->make<TH1D>("vx0","Delta Vx", 100,-5, 5);
  histo_vx[1] = fs->make<TH1D>("vx1","Delta Vy", 100,-5, 5);
  histo_vx[2] = fs->make<TH1D>("vx2","Delta Vz", 100,-5, 5);
  
  histo_iPiPhi = fs->make<TH1D>("iPiPhi", "Delta Phi",  66, -3.3,3.3);
  histo_DeltaR = fs->make<TH1D>("DeltaR", "Delta R (uu-pi)" , 100,0,10);
}
Hec2sLambda::~Hec2sLambda()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}
//
// member functions
//

// ------------ method called for each event  ------------
void
Hec2sLambda::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   //--Tracks Collection
   edm::Handle<reco::TrackCollection> tracks;
   iEvent.getByLabel(trackTags_,tracks); int allTrk = tracks->size();
   
   //--Muon Collection
   edm::Handle<reco::MuonCollection> allmuons;
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
   edm::ESHandle<GlobalTrackingGeometry> globTkGeomHandle;
   iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);
   //const MagneticField *BField = bFieldHandle.product();
   iSetup.get<GlobalTrackingGeometryRecord>().get(globTkGeomHandle);
   
   //--Load Tracks into a vector of Tracks (from Eduardo CascadeFitter.cc)
   std::vector<reco::Track> theTracks;
   theTracks.insert( theTracks.end(), 
                     tracks->begin(), 
		     tracks->end()  );

   //--Loop over tracks to convert to transient tracks (from Eduardo CascadeFitter.cc)
   std::vector<reco::TrackRef> theTrackRefs;
   std::vector<reco::TransientTrack> theTransTracks;
   for(unsigned int itrack = 0; itrack < tracks->size(); itrack++) {
         reco::TrackRef tmpRef( tracks, itrack );
         reco::TransientTrack tmpTk2( *tmpRef, &(*bFieldHandle), globTkGeomHandle );
         theTrackRefs.push_back( tmpRef );
	 theTransTracks.push_back( tmpTk2 );
   }   
   histo_Ntrk   ->Fill( allTrk );
   histo_Nmu    ->Fill( allMu ); 
   histo_Nlambda->Fill( allLambda ); 
   histo_Nprim  ->Fill( allPrim );
   
   if( MyPrint_ )std::cout<<"Hec: Ntrk ="<<allTrk
                          <<" Nmuons ="<<allMu
                          <<" NLambda ="<<allLambda
                          <<" NPrim ="<<allPrim<<std::endl;
   init();     //--clean all structures
   
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
      reco::Particle::LorentzVector lamVM_P4 = LamRec.P4();            //--lambda0 P4 with Mass constrains
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
      double lamEpr = sqrt( pion_mass_c*pion_mass_c + prVM_P4.P()*prVM_P4.P() ); //--Assuming Pi for the Pr
      double lamEpi = sqrt( pion_mass_c*pion_mass_c + piVM_P4.P()*piVM_P4.P() );
      double lamE   = lamEpr + lamEpi;      
      double lamPx  = prVM_P4.Px() + piVM_P4.Px();    //--Momentum at Vtx
      double lamPy  = prVM_P4.Py() + piVM_P4.Py();
      double lamPz  = prVM_P4.Pz() + piVM_P4.Pz();
      
      reco::Particle::LorentzVector K0s(0.0,0.0,0.0,0.0);  //--Fill Lorentz K0short P4                  
      K0s.SetPxPyPzE( lamPx, lamPy, lamPz, lamE );
      //double MK0s = K0s.M();                    //--same as myLamM
      //double myLamM  = sqrt( lamE*lamE - lamPx*lamPx - lamPy*lamPy - lamPz*lamPz );
      double myLamM  = K0s.M();  //--Kshort when Pr is a Pi
      double myfLamM = sqrt( lamVM_P4.E()*lamVM_P4.E() - lamVM_P4.Px()*lamVM_P4.Px()   //--Mass constrains P4
                                                       - lamVM_P4.Py()*lamVM_P4.Py()
                                                       - lamVM_P4.Pz()*lamVM_P4.Pz() );
      if( MyPrint_ )
          std::cout<<MlamMFit<<" M:"<<myfLamM<<" Class:"<<Mlambda0-myLamM<<" Fit:"<<MlamVFit-myLamM<<" <--"<<myLamM
          <<std::endl;
      
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
      for(reco::MuonCollection::const_iterator iMuon=allmuons->begin();iMuon!= allmuons->end();++iMuon){ 
         int iQ      = iMuon->charge();
         double iEta = iMuon->eta();
         double iPhi = iMuon->phi();           
         double iuE = sqrt( muon_mass_c*muon_mass_c + iMuon->p4().P()*iMuon->p4().P() );
         reco::Particle::LorentzVector iMup4(0.0,0.0,0.0,0.0);                      //--Fill Mu P4 Vector                 
         iMup4.SetPxPyPzE(iMuon->p4().Px(),iMuon->p4().Py(),iMuon->p4().Pz(),iuE);  //--TLorentzVector a(0,0,0,0)
         
         double iCal = iMuon->caloCompatibility();
         double iSeg = muon::segmentCompatibility(*iMuon);
         double iIso = iMuon->isolationR03().sumPt;
         reco::TrackRef iMuTrackRef = iMuon->innerTrack(); //--.index();
        
         int muSA=0, muGL=0, muTK=0;          //--Muon Id
         if( iMuon->isStandAloneMuon() )muSA = 1;   //--2^0
         if( iMuon->isGlobalMuon()     )muGL = 10;  //--2^1
         if( iMuon->isTrackerMuon()    )muTK = 100; //--2^2  
         int imuId = muTK + muGL + muSA;
        
         for(reco::MuonCollection::const_iterator jMuon=iMuon+1;jMuon!=allmuons->end();++jMuon){ 
            int jQ      = jMuon->charge();
            double jEta = jMuon->eta();
            double jPhi = jMuon->phi();         
            double juE = sqrt( muon_mass_c*muon_mass_c + jMuon->p4().P()*jMuon->p4().P() );
            reco::Particle::LorentzVector jMup4(0.0,0.0,0.0,0.0);                          
            jMup4.SetPxPyPzE(jMuon->p4().Px(),jMuon->p4().Py(),jMuon->p4().Pz(),juE);

            double jCal = jMuon->caloCompatibility();
            double jSeg = muon::segmentCompatibility(*jMuon);
            double jIso = jMuon->isolationR03().sumPt;
            reco::TrackRef jMuTrackRef = jMuon->innerTrack();   //--.index();
            
            muSA=0, muGL=0, muTK=0;            //--Muon Id
            if( jMuon->isStandAloneMuon() )muSA = 1;   //--2^0
            if( jMuon->isGlobalMuon()     )muGL = 10;  //--2^1
            if( jMuon->isTrackerMuon()    )muTK = 100; //--2^2
            int jmuId = muTK + muGL + muSA;
            
            double Muu = (iMuon->p4() + jMuon->p4() ).M();
            
            bool goodDimuon = (imuId==111 && jmuId>=11)    //--1 G and 1 (SG or T or ST or G)
                            ||(jmuId==111 && imuId>=11);

          //bool mySelDimu = std::abs( Muu - jpsi_mass_c )<0.2;    //--200 MeV
            bool mySelDimu = Muu > 2.8 && Muu < 4.0;
               //mySelDimu = true;
            if( iQ*jQ < 0 )
              histo_dimuM[2]->Fill( Muu );
            else
              histo_dimuM[3]->Fill( Muu );
            if( ( iQ*jQ < 0 )&& goodDimuon && mySelDimu ){ //--Opp. charge & muon Id & unique Trk

              histo_dimuM[4]->Fill( Muu );
              histo_lam[4]  ->Fill( MlamVFit );

              bool repTrk = false;     //--if true I reject Same tracks
              if( protonveetrk == iMuTrackRef || protonveetrk == jMuTrackRef ||
                  pionveetrk   == iMuTrackRef || pionveetrk   == jMuTrackRef )repTrk = true;
              if( !repTrk ){
              
              int iNtrkHits = iMuTrackRef->numberOfValidHits();
              double id0    = iMuTrackRef->d0();
              double idz    = iMuTrackRef->dz();
           
              int jNtrkHits = jMuTrackRef->numberOfValidHits();
              double jd0    = jMuTrackRef->d0();
              double jdz    = jMuTrackRef->dz();
              
              //--Save Dimuon P4 in a LorentzVector
              double uuPx = iMuon->p4().Px() + jMuon->p4().Px();
              double uuPy = iMuon->p4().Py() + jMuon->p4().Py();
              double uuPz = iMuon->p4().Pz() + jMuon->p4().Pz();
              double uuEn = iuE + juE;
              reco::Particle::LorentzVector uuP4(0.0,0.0,0.0,0.0);                          
              uuP4.SetPxPyPzE(uuPx,uuPy,uuPz,uuEn);
              double uuEta = uuP4.eta();
              double uuPhi = uuP4.phi();
              double uuDR  = sqrt( uuEta*uuEta + uuPhi*uuPhi );

              //double uuPhi2 = atan2(uuPy,uuPx);  //--This is the same as uuPhi
              //if( uuPhi2<0 )uuPhi2 = uuPhi2 + 6.2832;          
              //std::cout<<"uuPhi="<<uuPhi<<" uuPhi2="<<uuPhi2<<std::endl;
                             
              //--Re-fit dimuons to a common Vertex
              reco::TransientTrack iMuTT( iMuTrackRef, &(*bFieldHandle) );
              reco::TransientTrack jMuTT( jMuTrackRef, &(*bFieldHandle) );

              //--Trajectory states to calculate DCA for the 2 tracks
              FreeTrajectoryState iState = pFactory.particle( iMuTT, muon_mass_c, chi, ndf, muon_sig )
                                         ->currentState().freeTrajectoryState();
              FreeTrajectoryState jState = pFactory.particle( jMuTT, muon_mass_c, chi, ndf, muon_sig )
                                         ->currentState().freeTrajectoryState();

              //--Measure distance between tracks at their closest approach
              ClosestApproachInRPhi cApp;           
              double dca_uu = 9999;               //--Flag  Closest Approach distance between muons
              GlobalPoint cxPt(999,999,999);      //--cxPt.x(): Crosing point in X
              cApp.calculate(iState, jState);     //--cxPt.y(): Crosing point in Y 
              if( cApp.status() ) {               //--cxPt.z(): Crosing point in Z
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
              if( dca_uu<=1 && 
                  sqrt( cxPt.x()*cxPt.x() + cxPt.y()*cxPt.y() )< 120  &&
                  std::abs( cxPt.z() )<300 ){
   
              histo_dimuM[7]->Fill( Muu );
              histo_lam[7]  ->Fill( MlamVFit );

              //--Here I do Di-pion
              for(unsigned int iPi=0; iPi<theTrackRefs.size(); iPi++){          //--pion i
                 reco::TrackRef  iPiTrack = theTrackRefs[iPi];
                 bool iPiTrkOk = iPiTrack!=protonveetrk && iPiTrack!=pionveetrk &&
                                 iPiTrack!=iMuTrackRef  && iPiTrack!=jMuTrackRef;
                 double iPiQ   = iPiTrack->charge();
                 double iPiEta = iPiTrack->eta();
                 double iPiPhi = iPiTrack->phi();
                 double iPiEn  = sqrt( pion_mass_c*pion_mass_c + iPiTrack->p()*iPiTrack->p() );
                 
                 //--Calculate delta phi between dimuon and the pion
                 double phi1 = uuPhi, phi2 = iPiPhi;
                 if( phi1<0 )phi1 = phi1 + 6.2832;
                 if( phi2<0 )phi2 = phi2 + 6.2832;
                 double iD_Phi = phi1 - phi2;
                 if( iD_Phi>3.1416 )
                     iD_Phi = iD_Phi - 6.2832;
                 else if( iD_Phi<-3.1416 )
                          iD_Phi = iD_Phi + 6.2832;                 
                 //double iD_Phi  = uuPhi-iPiPhi;
                 //if( iD_Phi>3.1416 ) iD_Phi = 6.2832 - iD_Phi;
                 double ipiuuR = sqrt( (uuEta-iPiEta)*(uuEta-iPiEta) + iD_Phi*iD_Phi );
                 histo_iPiPhi->Fill(iD_Phi);
                 histo_DeltaR->Fill( ipiuuR );
                 
                 double iPiNtrkHits = iPiTrack->numberOfValidHits();
                 double iPid0       = iPiTrack->d0();
                 double iPidz       = iPiTrack->dz();
           
                 reco::Particle::LorentzVector iPip4(0.0,0.0,0.0,0.0);  //--Fill Pion P4                  
                 iPip4.SetPxPyPzE( iPiTrack->px(), iPiTrack->py(), iPiTrack->pz(), iPiEn );
                      
                 for(unsigned int jPi=iPi+1; jPi<theTrackRefs.size(); jPi++){   //--pion j
                    reco::TrackRef  jPiTrack = theTrackRefs[jPi];
                    bool jPiTrkOk = jPiTrack!=protonveetrk && jPiTrack!=pionveetrk &&
                                    jPiTrack!=iMuTrackRef  && jPiTrack!=jMuTrackRef;
                    double jPiQ   = jPiTrack->charge();
                    double jPiEta = jPiTrack->eta();
                    double jPiPhi = jPiTrack->phi();
                    double jPiEn  = sqrt( pion_mass_c*pion_mass_c + jPiTrack->p()*jPiTrack->p() );
                    
                    //--Calculate delta phi between dimuon and the pion
                    phi1 = uuPhi, phi2 = jPiPhi;
                    if( phi1<0 )phi1 = phi1 + 6.2832;
                    if( phi2<0 )phi2 = phi2 + 6.2832;
                    double jD_Phi = phi1 - phi2;
                    if( jD_Phi>3.1416 )
                        jD_Phi = jD_Phi - 6.2832;
                    else if( jD_Phi<-3.1416 )
                             jD_Phi = jD_Phi + 6.2832;
                    double jpiuuR = sqrt( (uuEta-jPiEta)*(uuEta-jPiEta) + jD_Phi*jD_Phi );
                                      
                    double jPiNtrkHits = jPiTrack->numberOfValidHits();
                    double jPid0       = jPiTrack->d0();
                    double jPidz       = jPiTrack->dz();
           
                    reco::Particle::LorentzVector jPip4(0.0,0.0,0.0,0.0);  //--Fill Pion P4                 
                    jPip4.SetPxPyPzE( jPiTrack->px(), jPiTrack->py(), jPiTrack->pz(), jPiEn );

                    double Muupipi = (iMup4 + jMup4 + iPip4 + jPip4).M();  //--4 body invariant mass
                    histo_mdiff[0]->Fill( Muupipi - Muu );

                    //--Check Invariant Mass Calculation for Mpipi
                    double Enpipi = iPiEn + jPiEn;
                    double Pxpipi = iPiTrack->px() + jPiTrack->px();
                    double Pypipi = iPiTrack->py() + jPiTrack->py();
                    double Pzpipi = iPiTrack->pz() + jPiTrack->pz();
                    double Ptpipi = sqrt( Pxpipi*Pxpipi + Pypipi*Pypipi );
                    double Mpipi  = sqrt( Enpipi*Enpipi - Pxpipi*Pxpipi
                                                        - Pypipi*Pypipi
                                                        - Pzpipi*Pzpipi );
                    histo_mpipi[0]->Fill( Mpipi );

                    //--Select good diPions to avoid so many combinations                  
                    bool goodPions = iPiQ*jPiQ<0 && iPiTrkOk && jPiTrkOk
                                  && iPiTrack->pt()>0.25 && jPiTrack->pt()>0.25 && Ptpipi>0.75
                                  && ipiuuR<1 && jpiuuR<1;
                    if( goodPions ){
                    histo_mdiff[1]->Fill( Muupipi - Muu );

                    //--dimuon + Pion Vertex
                    reco::TransientTrack iPiTT( *iPiTrack    , &(*bFieldHandle) );
                    reco::TransientTrack jPiTT( *jPiTrack    , &(*bFieldHandle) );
           
                    std::vector<RefCountedKinematicParticle> uupipiParticle;
                    uupipiParticle.push_back(pFactory.particle( iMuTT, muon_mass_c, chi, ndf, muon_sig ));
                    uupipiParticle.push_back(pFactory.particle( jMuTT, muon_mass_c, chi, ndf, muon_sig ));
                    uupipiParticle.push_back(pFactory.particle( iPiTT, pion_mass_c, chi, ndf, pion_sig ));
                    uupipiParticle.push_back(pFactory.particle( jPiTT, pion_mass_c, chi, ndf, pion_sig ));
       
                    //--Make sure that all tracks are close together so we can try to find a common vertex
                    FreeTrajectoryState iPiState = pFactory.particle(iPiTT,pion_mass_c,chi,ndf,pion_sig)
                                                  ->currentState().freeTrajectoryState();
                    FreeTrajectoryState jPiState = pFactory.particle(jPiTT,pion_mass_c,chi,ndf,pion_sig)
                                                  ->currentState().freeTrajectoryState();
                                                  
                    //--Measure distance between tracks at their closest approach
                    //--muon-i with pion-i and with pion-j
                    //--if muon-i is close enough to pion-i, then
                    //--pion-i is close also to muon-j because muon-i & muon-j are close together (see dca_uu)
                    //--but I'm doing all to make sure they are all close together
                    double dcaii= 9999, dcaij= 9999,
                           dcaji= 9999, dcajj= 9999;
                    GlobalPoint cxPtii(999,999,999),cxPtij(99,99,99),
                                cxPtji(999,999,999),cxPtjj(99,99,99);
                    cApp.calculate(iState, iPiState);
                    if( cApp.status() ) { 
                      dcaii  = fabs( cApp.distance() );
                      cxPtii = cApp.crossingPoint();
                    }
                    cApp.calculate(iState, jPiState);
                    if( cApp.status() ) {
                      dcaij  = fabs( cApp.distance() );
                      cxPtij = cApp.crossingPoint();
                    }
                    cApp.calculate(jState, iPiState);
                    if( cApp.status() ) { 
                      dcaji  = fabs( cApp.distance() );
                      cxPtji = cApp.crossingPoint();
                    }
                    cApp.calculate(jState, jPiState);
                    if( cApp.status() ) {
                      dcajj  = fabs( cApp.distance() );
                      cxPtjj = cApp.crossingPoint();
                    }
                    //--r = sqrt( x^2 + y^2 )
                    //--Cut on fiducial tracking volume [-120 cm < r < 120 cm  and -260 cm < Z < 260 cm]
                    //--Pixel [-20 cm < r < 20 cm  and -60 cm < Z < 60 cm] to run vertex fit
                    bool allClose = dcaii<1 && dcaij<1 && dcaji<1 && dcajj<1;
                    bool allPoint = sqrt( cxPtii.x()*cxPtii.x() + 
                                          cxPtii.y()*cxPtii.y() )<120 && std::abs( cxPtii.z() )<300
                                 && sqrt( cxPtij.x()*cxPtij.x() + 
                                          cxPtij.y()*cxPtij.y() )<120 && std::abs( cxPtij.z() )<300
                                 && sqrt( cxPtji.x()*cxPtji.x() + 
                                          cxPtji.y()*cxPtji.y() )<120 && std::abs( cxPtji.z() )<300
                                 && sqrt( cxPtjj.x()*cxPtjj.x() + 
                                          cxPtjj.y()*cxPtjj.y() )<120 && std::abs( cxPtjj.z() )<300;
                    if( allClose && allPoint ){

                    histo_mdiff[2]->Fill( Muupipi - Muu );
                    
                    KinematicFitDriver uupipiRec( uupipiParticle, "MuuMpipi" );  //--Dimuon-PiPi fit
                    if( uupipiRec.isValid() ){
                    
                    //--Dimuon-PiPi fitted variables
                    double MuupipiVFit = uupipiRec.mass();
                    double uupipiV_CL  = ChiSquaredProbability( uupipiRec.chi2(), uupipiRec.ndof() );
                    
                    //--Here Dimu + Dipi Mass constrains to massC_value psi(2s)
                    double massC_value = 0;
                    bool uupipi2sMcut = std::abs( MuupipiVFit - psi2s_mass_c)<0.2; //--200 MeV
                    if(  uupipi2sMcut )massC_value = psi2s_mass_c;
              
                    double MuupipiMFit = -1;   //--Mass Fitted Dimuon mass 
                    double uupipiVM_CL = -1;   //--Confidence Level uu Vtx
                    if( uupipi2sMcut && mySelDimu ){       
                      uupipiRec.AddMassConstraint( massC_value, massC_value*1.e-6 );          
                      if( uupipiRec.isValid() ){ 
                        MuupipiMFit = uupipiRec.mass();
                        uupipiVM_CL = ChiSquaredProbability( uupipiRec.chi2(), uupipiRec.ndof() );
                      }  //--End good Psi(2s) mass constrain            
                    }    //--End  uupipi2sMcut              
                    //--Dimuon fitted variables
                    reco::Particle::LorentzVector uupipiV_P4  = uupipiRec.P4();         //--DiMu-PiPi P4
                    reco::Particle::LorentzVector iMu2V_P4 = uupipiRec.P4FirstChild();  //--Mu-i
                    reco::Particle::LorentzVector jMu2V_P4 = uupipiRec.P4NextChild();   //--Mu-j
                    reco::Particle::LorentzVector iPi2V_P4 = uupipiRec.P4NextChild();   //--Pi-i
                    reco::Particle::LorentzVector jPi2V_P4 = uupipiRec.P4NextChild();   //--Pi-i
                    reco::Particle::Point uupipiV_vx = uupipiRec.VertexDecay();         //--Vx,Vy,Vz 
                    reco::Vertex::CovarianceMatrix uupipiV_Cov = uupipiRec.VertexCov();

                    //--dimuon parameters
                    double iPx  = iMu2V_P4.Px();
                    double iPy  = iMu2V_P4.Py();
                    double iPz  = iMu2V_P4.Pz();
                    double iEn  = sqrt( muon_mass_c*muon_mass_c + iMu2V_P4.P()*iMu2V_P4.P() );
              
                    double jPx  = jMu2V_P4.Px();
                    double jPy  = jMu2V_P4.Py();
                    double jPz  = jMu2V_P4.Pz();
                    double jEn  = sqrt( muon_mass_c*muon_mass_c + jMu2V_P4.P()*jMu2V_P4.P() );
                    
                    reco::Particle::LorentzVector uuV_P4;                
                    uuV_P4.SetPxPyPzE( iPx+jPx, iPy+jPy, iPz+jPz, iEn+jEn );
                    
                    double MuuVFit = (iMu2V_P4 + jMu2V_P4).M();
                    double uuV_CL  = uupipiV_CL;
                    double uuDist  = -1;
                    double uuR     = -1;
                    double MuuMFit = MuupipiMFit;  //--these 2 tmp var are for uupipi instead than uu
                    double uuVM_CL = uupipiVM_CL;
                     
                    bool mySel2Dimu = std::abs( MuuVFit - jpsi_mass_c )<0.2;    //--200 MeV
                            
                    //if( uupipiV_CL>0.01 ){
                    histo_diff[1] ->Fill( Muu - MuuVFit );  
                    histo_mdiff[3]->Fill( MuupipiVFit - MuuVFit );
                    
                    //--Dimuon and Lambda-0 Vertex [ Vertexing 2 neutral is bad ? ]                   
                    //--Alternative Vertexing for a Neutral
                    //--Check if Lambda0 comes from DimuonDiPion Vertex and calculate cosine in 3D and 2D
                    GlobalVector Vector_uupipiV_to_lamVM( lamVM_vx.x()- uupipiV_vx.x(),
                                                          lamVM_vx.y()- uupipiV_vx.y(),
                                                          lamVM_vx.z()- uupipiV_vx.z() );                           
                    GlobalVector Lambda_P3( lamVM_P4.px(), lamVM_P4.py(), lamVM_P4.pz() ); 
                    double cos_uupipil3D = ( Lambda_P3.dot( Vector_uupipiV_to_lamVM ) )/
                                       ( Vector_uupipiV_to_lamVM.mag() * Lambda_P3.mag() );

                    double cos_uupipil2D = ( (lamVM_vx.x()- uupipiV_vx.x())*lamVM_P4.px() + 
                                           (lamVM_vx.y()- uupipiV_vx.y())*lamVM_P4.py() )
                                    /sqrt( (lamVM_vx.x()- uupipiV_vx.x())*(lamVM_vx.x()- uupipiV_vx.x()) +
                                           (lamVM_vx.y()- uupipiV_vx.y())*(lamVM_vx.y()- uupipiV_vx.y()) )
                                    /sqrt(  lamVM_P4.px()*lamVM_P4.px() +
                                            lamVM_P4.py()*lamVM_P4.py() );

                    histo_cos3d->Fill( cos_uupipil3D );
                    histo_cos2d->Fill( cos_uupipil2D );
                    //--Skip Events here:
                    if( cos_uupipil3D > 0.9 ){
                                               
                    histo_mdiff[4]->Fill( MuupipiVFit - MuuVFit );

                    //--Check separation Lambda from Dimuon (L over Sigma) in 3D
                    typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > 
                                       SMatrixSym3D_uupipiL;
                    typedef ROOT::Math::SVector<double, 3> 
                                       SVector3_uupipiL;

                    SMatrixSym3D_uupipiL totalCov_uupipiL = lamVM_Cov + uupipiV_Cov; //--3D Separation 2u2Pi-Lambda
                    SVector3_uupipiL uupipiL_PrimarySep3D( uupipiV_vx.x() - lamVM_vx.x(),
                                                           uupipiV_vx.y() - lamVM_vx.y(),
                                                           uupipiV_vx.z() - lamVM_vx.z() ); 
                    double uupipiL_L3D     = ROOT::Math::Mag( uupipiL_PrimarySep3D );
                    double uupipiL_Sigma3D = sqrt(ROOT::Math::Similarity( totalCov_uupipiL, uupipiL_PrimarySep3D ))
                                           / uupipiL_L3D;

                    //--Check separation Lambda from Dimuon (L over Sigma) in 2D  [need to check] ? ? ?
                    double lam_bV_CL =  1;
                    double lambDist  = -1;
                    double lambR     = -1;
                    double Mlam_bVFit = ( uupipiRec.P4() + LamRec.P4() ).M();
                    //double Mlam_buFit = ( iMuV_P4->P4() + jMuV_P4->P4() + LamRec.P4() ).M();  //--unFitted Mass for uu(unfitted) + Lambda0(fitted)

                    double MuuK0s = ( uupipiRec.P4() + K0s ).M(); 

                    histo_dimuM[8]->Fill( MuuVFit );
                    histo_lam[8]  ->Fill( MlamVFit );
                    histo_lamb[8] ->Fill( Mlam_bVFit );

                    //--Primary Vertex (excluding dimuon and Lambda tracks)
                    std::vector<reco::TrackRef> ExclusionList;      
                    ExclusionList.push_back( pionveetrk );       
                    ExclusionList.push_back( protonveetrk );
                    ExclusionList.push_back( iMuTrackRef );
                    ExclusionList.push_back( jMuTrackRef ); 
                    ExclusionList.push_back( iPiTrack );  
                    ExclusionList.push_back( jPiTrack ); 
                    VertexRefit RefitPrimary(primary.BestVertex(),ExclusionList,bFieldHandle,beamSpot);
               
                    reco::Vertex refitVertexPrim = RefitPrimary.Refitted();
                    GlobalPoint PosVertex = GlobalPoint(refitVertexPrim.x(),
                                                        refitVertexPrim.y(),
                                                        refitVertexPrim.z() );
 
                    double primCL   = ChiSquaredProbability( refitVertexPrim.chi2(), refitVertexPrim.ndof());
                    double primDist = sqrt( refitVertexPrim.x()*refitVertexPrim.x()
                                          + refitVertexPrim.y()*refitVertexPrim.y()
                                          + refitVertexPrim.z()*refitVertexPrim.z() );
                    double primR    = sqrt( refitVertexPrim.x()*refitVertexPrim.x()
                                          + refitVertexPrim.y()*refitVertexPrim.y() );
                                          
                    //--Calculate Lambda_b Mass, L/Sigma for Primary, L/Sigma Dimuon-Lambda

                    //--Check separation Lambda_b from primary (L over Sigma)
                    typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> >
                                       SMatrixSym3D_Lambda_b;
                    typedef ROOT::Math::SVector<double, 3> SVector3_Lambda_b;

                    //-->SMatrixSym3D_Lambda_b totalCov_Lambda_b = refitVertexPrim.covariance() + lam_bV_Cov;
                    //-->SVector3_Lambda_b Lambda_bPrimarySep3D( lam_bV_vx.x() - refitVertexPrim.x(),
                    //-->                                        lam_bV_vx.y() - refitVertexPrim.y(),
                    //-->                                        lam_bV_vx.z() - refitVertexPrim.z() );
               
                    //--Take the dimuon vertex as the Secondary and measure L/Sigma from there (Dec 27, 11)
                    SMatrixSym3D_Lambda_b totalCov_Lambda_b = refitVertexPrim.covariance() + uupipiV_Cov;
                    SVector3_Lambda_b Lambda_bPrimarySep3D( uupipiV_vx.x() - refitVertexPrim.x(),
                                                            uupipiV_vx.y() - refitVertexPrim.y(),
                                                            uupipiV_vx.z() - refitVertexPrim.z() );
 
                    double Lambda_bL3D     = ROOT::Math::Mag( Lambda_bPrimarySep3D );
                    double Lambda_bSigma3D = sqrt(ROOT::Math::Similarity( totalCov_Lambda_b, Lambda_bPrimarySep3D ))
                                           / Lambda_bL3D;
                    double Lambda_bLoS3D   = Lambda_bL3D / Lambda_bSigma3D;
 
                    //--Here Pointing Secondary Vertex  (dimuon + lambda0) to Primary  [Sep 9, 2011]
 
                    //-->GlobalVector Lambda_b_PrimSec( lam_bV_vx.x() - refitVertexPrim.x(),
                    //-->                               lam_bV_vx.y() - refitVertexPrim.y(),
                    //-->                               lam_bV_vx.z() - refitVertexPrim.z() );                            
                    //-->GlobalVector Lambda_b_P3( lam_bV_P4.px(), lam_bV_P4.py(), lam_bV_P4.pz() );
 
                    GlobalVector Lambda_b_PrimSec( uupipiV_vx.x() - refitVertexPrim.x(),
                                                   uupipiV_vx.y() - refitVertexPrim.y(),
                                                   uupipiV_vx.z() - refitVertexPrim.z() );                            
                    GlobalVector Lambda_b_P3( (uupipiV_P4+lamVM_P4).px(),
                                              (uupipiV_P4+lamVM_P4).py(),
                                              (uupipiV_P4+lamVM_P4).pz() );
                    double V_dot_Plamb = Lambda_b_P3.dot( Lambda_b_PrimSec );
 
                    double cos_alphab3D = V_dot_Plamb /( Lambda_b_PrimSec.mag() * Lambda_b_P3.mag() );
               
                    //--Cosine in 2D
                    double cos_alphab2D = ( (uupipiV_vx.x() - refitVertexPrim.x())*(uupipiV_P4+lamVM_P4).px() + 
                                            (uupipiV_vx.y() - refitVertexPrim.y())*(uupipiV_P4+lamVM_P4).py() )
                                     /sqrt( (uupipiV_vx.x() - refitVertexPrim.x())*(uupipiV_vx.x() - refitVertexPrim.x()) +
                                            (uupipiV_vx.y() - refitVertexPrim.y())*(uupipiV_vx.y() - refitVertexPrim.y()) )
                                     /sqrt( (uupipiV_P4+lamVM_P4).px()*(uupipiV_P4+lamVM_P4).px() +
                                            (uupipiV_P4+lamVM_P4).py()*(uupipiV_P4+lamVM_P4).py() );
                                            
                    histo_vx[0]->Fill( refitVertexPrim.x() - uupipiV_vx.x() );
                    histo_vx[1]->Fill( refitVertexPrim.y() - uupipiV_vx.y() );
                    histo_vx[2]->Fill( refitVertexPrim.z() - uupipiV_vx.z() );

                    if( Lambda_bLoS3D>3 && primCL>0.01 && uupipi2sMcut && uupipiVM_CL>0.01
                     && lamMcut && mySel2Dimu){
                        histo_dimuM[9]->Fill( MuuVFit );
                        histo_lam[9]  ->Fill( MlamVFit );
                        histo_lamb[9] ->Fill( Mlam_bVFit );
                        histo_mpipi[1]->Fill( Mpipi );
                        histo_mdiff[5]->Fill( MuupipiVFit - MuuVFit );

                    }
                    fill_pipi(iPiQ,iPiEta,iPiPhi,iPiTrack->px(),iPiTrack->py(),iPiTrack->pz(),iPiNtrkHits,iPid0,iPidz,
                              jPiQ,jPiEta,jPiPhi,jPiTrack->px(),jPiTrack->py(),jPiTrack->pz(),jPiNtrkHits,jPid0,jPidz,
                              MuupipiVFit,uupipiV_CL,MuupipiMFit,uupipiVM_CL,
                              uupipiV_P4.Px(),uupipiV_P4.Py(),uupipiV_P4.Pz(),
                              ipiuuR,jpiuuR);                              
                    fill_lambdab(Mlam_bVFit, lam_bV_CL, lambDist, lambR, primCL, primDist, primR,
                                 Lambda_bL3D, Lambda_bSigma3D, cos_alphab3D, cos_alphab2D);
                             
                    //--Fill Ntuple and re-initialize all variables
                    iComb++;
                    fill_evt( iEvent, allTrk, allMu, allPrim, allLambda, iComb );
                    fill_IuuC(iQ, imuId, iNtrkHits, jQ, jmuId, jNtrkHits );
                    fill_uuC(iEta,iPhi,iPx,iPy,iPz,id0,idz,iCal,iSeg,iIso,
                             jEta,jPhi,jPx,jPy,jPz,jd0,jdz,jCal,jSeg,jIso,
                             Muu, dca_uu, cxPt.x(), cxPt.y(), cxPt.z(), uuDist, uuR,uuDR,
                             MuuVFit, uuV_CL,uuV_P4.Px(),uuV_P4.Py(),uuV_P4.Pz(),
                             MuuMFit, uuVM_CL, cos_uupipil3D, cos_uupipil2D, uupipiL_L3D, uupipiL_Sigma3D);
                         
                    hecmu_tree_->Fill();
                    //}   //--if( Lambda_bRec.isValid() ){ //--Dimuon+Lambda Vertex
                    //}   //--if( uuV_CL > 0.01 ){   //--dimuon CL                                           
                    } //--if( cos_uul2D > 0.9 ){
//                     } //--End if( uupipiV_CL>0.1 )
                    } //--End if( uupipiRec.isValid() ){
                    } //--End if( allClose && allPoint){
                    } //--End goodPions
                 }    //--End for(unsigned int jPi=0                 
              }       //--End for(unsigned int iPi=0
//               }       //--else{  valid Dimuon Vertex
              }       //--Fiducial cut for dimuon
              }       //--no repeated tracks between muons and Lambda
            }         //--Opp charge & Muon Id       
         }            //--end muonCollection for j              
      }               //--end muonCollection for i
   }                  //--end if(  LamRec.isValid()   Good Lambda Mass Constrain
   }                  //--end if(  LamRec.isValid()   Good Lambda Vertex
}                     //--end for( unsigned int veeindex )    

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
Hec2sLambda::beginJob()
{  
  std::cout << " Hec2sLambda::beginJob" << std::endl;
  hecmu_tree_ = new TTree("Lambda","Dettach Ntuple");
  int bufsize = 64000;

  hecmu_tree_->Branch("Dimuon",  &uuC_,    "iEta/D:iPhi:iPx:iPy:iPz:id0:idz:iCal:iSeg:iIso:jEta:jPhi:jPx:jPy:jPz:jd0:jdz:jCal:jSeg:jIso:Muu:dca:xcpt:ycpt:zcpt:uuDist:uuR:uuDR:MuuVF:uuVCL:uuVMpx:uuVMpy:uuVMpz:MuuMFit:uuVMCL:cosuul3D:cosuul2D:uuLL3D:uuLSigma3D", bufsize);
  hecmu_tree_->Branch("Lambda0", &lambdaC_,"prQ/D:prEta:prPhi:prPx:prPy:prPz:prd0:prdz:prNHits:piQ:piEta:piPhi:piPx:piPy:piPz:pid0:pidz:piNHits:Mlambda0:lamCL:Mlambda0VF:lamDist:lamR:Mkshort:MlamMFit:lamVMCL:lamVMpx:lamVMpy:lamVMpz", bufsize);  
  hecmu_tree_->Branch("Lambdab", &lambdabC_,"lambCL/D:MlambVF:lambDist:lambR:primCL:primDist:primR:lambdabL3D:lambdabSigma3D:cosAlphab3D:cosAlphab2D", bufsize);  
  hecmu_tree_->Branch("uupipi",  &uupipiC_, "iPiQ/D:iPiEta:iPiPhi:iPiPx:iPiPy:iPiPz:iPiNtrkHits:iPid0:iPidz:jPiQ:jPiEta:jPiPhi:jPiPx:jPiPy:jPiPz:jPiNtrkHits:jPid0:jPidz:MuupipiVFit:uupipiV_CL:MuupipiMFit:uupipiVM_CL:uupipiVPx:uupipiVPy:uupipiVPz:ipiuuR:jpiuuR", bufsize); 

  hecmu_tree_->Branch("iHeader", &Ievt_, "run/I:evtnum:lumib:allTrk:allMu:allPrim:allLambda:iComb", bufsize);
  hecmu_tree_->Branch("iDimuon", &IuuC_, "iQ/I:imuId:iNtrkHits:jQ:jmuId:jNtrkHits", bufsize);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
Hec2sLambda::endJob() 
{
  std::cout << " Hec2sLambda::endJob" << std::endl;
}

// ------------ method called when starting to processes a run  ------------
void 
Hec2sLambda::beginRun(edm::Run const&, edm::EventSetup const&)
{
  std::cout << " Hec2sLambda::beginRun" << std::endl;
}

// ------------ method called when ending the processing of a run  ------------
void 
Hec2sLambda::endRun(edm::Run const&, edm::EventSetup const&)
{
  std::cout << " Hec2sLambda::endRun" << std::endl;
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
Hec2sLambda::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
  std::cout << " Hec2sLambda::beginLuminosityBlock" << std::endl;
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
Hec2sLambda::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
  std::cout << " Hec2sLambda::endLuminosityBlock" << std::endl;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Hec2sLambda::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
void Hec2sLambda::init()
{
  Ievt_.init();
  IuuC_.init();
  uuC_.init();
  lambdaC_.init();
  lambdabC_.init();
  uupipiC_.init();
}
void Hec2sLambda::Ievt::init()
{
  runNb      = -1; eventNb    = -1; lumiBlock  = -1;
  allTrk     = -1; allMu      = -1; allPrim    = -1; allLambda  = -1;
  iComb      =  0;
}
void Hec2sLambda::fill_evt(const edm::Event& iEvent, int allTrk, int allMu, int allPrim, int allLambda, int iComb )
{
     Ievt_.runNb      = iEvent.id().run(); 
     Ievt_.eventNb    = iEvent.id().event();
     Ievt_.lumiBlock  = iEvent.luminosityBlock();
     Ievt_.allTrk     = allTrk;         //--tracks->size();
     Ievt_.allMu      = allMu;          //--->allmuons->size();
     Ievt_.allPrim    = allPrim;
     Ievt_.allLambda  = allLambda;
     Ievt_.iComb      = iComb;
     if( MyPrint_ )
       std::cout<<"Hec: Run: "<<Ievt_.runNb<<" Evt: "<<Ievt_.eventNb<<" Lumi:"<<Ievt_.lumiBlock
                <<" Ntk: "<<Ievt_.allTrk<<" Mu:"<<Ievt_.allMu
                <<" NPrimVtx:"<<Ievt_.allPrim<<" #Lambda:"<<allLambda
                <<std::endl;
}
void Hec2sLambda::IuuC::init()
{
  iQ = 0; imuId = -1; iNtrkHits = -1;
  jQ = 0; jmuId = -1; jNtrkHits = -1;
}
void Hec2sLambda::fill_IuuC(int iQ, int imuId, int iNtrkHits,
                             int jQ, int jmuId, int jNtrkHits)
{
  IuuC_.iQ = iQ; IuuC_.imuId = imuId; IuuC_.iNtrkHits = iNtrkHits;
  IuuC_.jQ = jQ; IuuC_.jmuId = jmuId; IuuC_.jNtrkHits = jNtrkHits;
}
void Hec2sLambda::uuC::init()
{
  iEta =9.9; iPhi =-1.0; iPx =-999; iPy =-999; iPz=-999; id0 =-1.0; idz =-1.0; iCal =-1.0; iSeg =-1.0; iIso =-1.0;
  jEta =9.9; jPhi =-1.0; jPx =-999; jPy =-999; jPz=-999; jd0 =-1.0; jdz =-1.0; jCal =-1.0; jSeg =-1.0; jIso =-1.0;
  Muu = -1.0; dca = -1.0; xcpt = -1.0; ycpt = -1.0; zcpt = -1.0; uuDist = -1.0; uuR = -1.0; uuDR = 999; 
  MuuVF=-1.0; uuVCL=-1.0; uuVMpx=-1;uuVMpy=-1;uuVMpz=-1;
  MuuMFit=-1; uuVMCL=-1; cosuul3D=-999; cosuul2D=-999; uuLL3D=-1; uuLSigma3D=-1;
}
void Hec2sLambda::fill_uuC(double iEta,double iPhi,double iPx,double iPy,double iPz,double id0,double idz,double iCal,double iSeg,double iIso,
                            double jEta,double jPhi,double jPx,double jPy,double jPz,double jd0,double jdz,double jCal,double jSeg,double jIso,
                            double Muu,double dca,double  xcpt,double ycpt,double  zcpt, double dist, double R2dim, double uuDR,
                            double MuuFit, double uuV_CL, double uuVMpx, double uuVMpy, double uuVMpz,
                            double MuuMFit, double uuVM_CL, double cos_uul3D, double cos_uul2D, double uuL_L3D, double uuL_Sigma3D)
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
  uuC_.uuDR    = uuDR;
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
}
void Hec2sLambda::lambdaC::init()
{
  prQ = 0; prEta = -10; prPhi = -10; prPx = -9999; prPy = -9999; prPz = -9999; prd0 = -1; prdz = -1; prNHits = -1;
  piQ = 0; piEta = -10; piPhi = -10; piPx = -9999; piPy = -9999; piPz = -9999; pid0 = -1; pidz = -1; piNHits = -1;
  Mlambda0 = -1; lamCL      = -1; Mlambda0VF = -1; lamDist    = -1; lamR       = -1; Mkshort = -1;
  MlamMFit = -1; lamVMCL   = -1; lamVMpx = -999; lamVMpy = -999; lamVMpz = -999;
}
void Hec2sLambda::fill_lambda0(double prQ, double prEta, double prPhi, double prPx, double prPy, double prPz,
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
void Hec2sLambda::lambdabC::init()
{
  lambCL = -1; MlambVF = -1; lambDist = -1; lambR = -1;
  primCL = -1; primDist = -1; primR  = -1;
  lambdabL3D = -1; lambdabSigma3D = -1; cosalphab3D = -10; cosalphab2D = -10;
}
void Hec2sLambda::fill_lambdab(double Mlam_bVFit, double lam_bV_CL, double lambDist, double lambR, 
                                double primCL, double primDist, double primR,
                                double Lambda_bL3D, double Lambda_bSigma3D, double cos_alphab3D, double cos_alphab2D)
{
  lambdabC_.lambCL      = lam_bV_CL;
  lambdabC_.MlambVF     = Mlam_bVFit;
  lambdabC_.lambDist    = lambDist;
  lambdabC_.lambR       = lambR;
  lambdabC_.primCL      = primCL;
  lambdabC_.primDist    = primDist;
  lambdabC_.primR       = primR;
  lambdabC_.lambdabL3D      = Lambda_bL3D;
  lambdabC_.lambdabSigma3D  = Lambda_bSigma3D;
  lambdabC_.cosalphab3D     = cos_alphab3D;
  lambdabC_.cosalphab2D     = cos_alphab2D;
}
void Hec2sLambda::uupipiC::init()
{
 iPiQ=0;iPiEta=-9;iPiPhi=-999;iPiPx=-999;iPiPy=-999;iPiPz=-999;iPiNtrkHits=-999;iPid0=-999;iPidz=-999;
 jPiQ=0;jPiEta=-9;jPiPhi=-999;jPiPx=-999;jPiPy=-999;jPiPz=-999;jPiNtrkHits=-999;jPid0=-999;jPidz=-999;
 MuupipiVFit=-1;uupipiV_CL=-1;MuupipiMFit=-1;uupipiVM_CL =-1;
 uupipiVPx=-1;uupipiVPy=-1;uupipiVPz=-1;
 ipiuuR=999; jpiuuR=999; 
}
void Hec2sLambda::fill_pipi(double iQ,double iEta,double iPhi,double ipx,double ipy,double ipz,double iNtrkHits,double id0,double idz,
                            double jQ,double jEta,double jPhi,double jpx,double jpy,double jpz,double jNtrkHits,double jd0,double jdz,
                            double MuupipiVFit,double uupipiV_CL,double MuupipiMFit,double uupipiVM_CL,
                            double uupipiVPx,double uupipiVPy,double uupipiVPz,
                            double ipiuuR,double jpiuuR)
{
  uupipiC_.iPiQ        = iQ;
  uupipiC_.iPiEta      = iEta;
  uupipiC_.iPiPhi      = iPhi;
  uupipiC_.iPiPx       = ipx;
  uupipiC_.iPiPy       = ipy;
  uupipiC_.iPiPz       = ipz;
  uupipiC_.iPiNtrkHits = iNtrkHits;
  uupipiC_.iPid0       = id0;
  uupipiC_.iPidz       = idz;
  uupipiC_.jPiQ        = jQ;
  uupipiC_.jPiEta      = jEta;
  uupipiC_.jPiPhi      = jPhi;
  uupipiC_.jPiPx       = jpx;
  uupipiC_.jPiPy       = jpy;
  uupipiC_.jPiPz       = jpz;
  uupipiC_.jPiNtrkHits = jNtrkHits;
  uupipiC_.jPid0       = jd0;
  uupipiC_.jPidz       = jdz;
  uupipiC_.MuupipiVFit = MuupipiVFit;
  uupipiC_.uupipiV_CL  = uupipiV_CL;
  uupipiC_.MuupipiMFit = MuupipiMFit;
  uupipiC_.uupipiVM_CL = uupipiVM_CL;
  uupipiC_.uupipiVPx   = uupipiVPx;
  uupipiC_.uupipiVPy   = uupipiVPy;
  uupipiC_.uupipiVPz   = uupipiVPz; 
  uupipiC_.ipiuuR      = ipiuuR;
  uupipiC_.jpiuuR      = jpiuuR;
}
//Method to compute 3D absolute impact parameter
//based on information of Transversal Impact Parameter (2D) D0 and
//Longitudinal Impact Parameter (along Z)
// IP3d = sqrt(D0^2 + Dz^2)
//The error is propagated adding in quadrature:
// IP3de^2 = (d(IP3d)/d(D0))^2 * D0e^2 + (d(IP3d)/d(Dz))^2 * Dze^2
//simplified to:
// IP3de = sqrt(D0e^2*D0^2 + Dze^2*Dz^2)/IP3d

double Hec2sLambda::SignificanceAbsoluteImpactParameter3D(
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
DEFINE_FWK_MODULE(Hec2sLambda);
