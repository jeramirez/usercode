#ifndef DIMUONTOOLS_DIMUCANDIDATES_H
#define DIMUONTOOLS_DIMUCANDIDATES_H
// -*- C++ -*-
//
// Package:    HecBaryons
// Class:      DiMuCandidates
//
/**\class DiMuCandidates DiMuCandidates.h Analyzers/DiMuonTools/interface/DiMuCandidates.h

 Description: Select DiMu candidates based on configurable cuts

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Juan Eduardo Ramirez Vargas
//         Created:  Mon Jun 29 17:24:25 CDT 2009
// $Id: DiMuCandidates.h,v 1.3 2013/05/17 17:23:24 jramirez Exp $
//
//

// system include files
#include <memory>
//references
#include "FWCore/Framework/interface/ESHandle.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
// user include files
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicVertex.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
//tools
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

//bfield
#include "MagneticField/Engine/interface/MagneticField.h"

//primary
#include "Analyzers/CascadeProducer/interface/PrimaryInfo.h"

//
// class declaration
//

class DiMuCandidates {
   public:
      explicit DiMuCandidates(const edm::Event& iEvent);
      explicit DiMuCandidates(const edm::ParameterSet& iConfig
                                ,const edm::Event& iEvent
		                ,const edm::EventSetup& iSetup);
      ~DiMuCandidates();      
      bool IsValid(){return isValid;};
      reco::TrackRef Muontrack(unsigned int DiMuon, int leg){
        return leg==1?Muon1[DiMuon].innerTrack()
                    :leg==2?Muon2[DiMuon].innerTrack()
                           : reco::TrackRef();
      }
      const GlobalPoint DiMuVertexDecay(unsigned int DiMuon){
        return (MuuVertex.find(DiMuon)!=MuuVertex.end()?
                GlobalPoint(MuuVertex[DiMuon].x(),
                            MuuVertex[DiMuon].y(),
                            MuuVertex[DiMuon].z()
                           ):
                GlobalPoint(0,0,0)
               ); 
      }
      GlobalVector MomentumAtDecay(unsigned int DiMuon){
        return (MuuVertexP4.find(DiMuon)!=MuuVertexP4.end()?
                GlobalVector(MuuVertexP4[DiMuon].x(),
                             MuuVertexP4[DiMuon].y(),
                             MuuVertexP4[DiMuon].z()
                            ):
                GlobalVector(0,0,0) 
               );
      }
      double mass(unsigned int DiMuon){
        return MuuMassFit.find(DiMuon)!=MuuMassFit.end()?
                                      MuuMassFit[DiMuon]:0;
      };
      bool hasKinematicFit(unsigned int DiMuon){
        return DiMuon<IsKinematicFitReady.size()?
                     IsKinematicFitReady[DiMuon]:false;
      };
      std::vector<RefCountedKinematicParticle> Candidate(unsigned int DiMuon){
        return (hasKinematicFit(DiMuon)&&uuParticle.find(DiMuon)!=uuParticle.end()?
                     uuParticle[DiMuon]: std::vector<RefCountedKinematicParticle>() ) ;
      };
      reco::Vertex Vertex(unsigned int DiMuon){
        return (hasKinematicFit(DiMuon)&&MuuRecoVertex.find(DiMuon)!=MuuRecoVertex.end()?
                     MuuRecoVertex[DiMuon] : reco::Vertex() );
      };
      //degrees of Freedom of Particle Decay Vertex
      double ndof(unsigned int DiMuon){
        return (hasKinematicFit(DiMuon)?MuuVertexInfo[DiMuon]->degreesOfFreedom():-1.); 
      };
      //Chi Squared of Particle Decay Vertex
      double chi2(unsigned int DiMuon){
        return (hasKinematicFit(DiMuon)?MuuVertexInfo[DiMuon]->chiSquared():-1.); 
      };
      //Chi2 Probability
      double CL(unsigned int DiMuon){
        double  idof  =  ndof(DiMuon);
        return (idof>0?ChiSquaredProbability(chi2(DiMuon),idof):-1.);
      };
      unsigned int size(){
         return Muon1.size()==Muon2.size()?Muon1.size():0;
      };
      unsigned int totalmuons(){
         return theMuons.size();
      };
      reco::Vertex PrimaryRefitted(unsigned int DiMuon);
      double delsig(unsigned int DiMuon);
//      double SignificanceImpactParameter3DAtPrimary(unsigned int casindex);
//      double SignificanceImpactParameter3DAtVertex(unsigned int casindex,reco::Vertex &PrmVtx);
      bool dcaisvalid(FreeTrajectoryState &iState, FreeTrajectoryState &jState, double &dca_uu, GlobalPoint &cxPt);
      int muGL(unsigned int DiMuon, int leg);
      void init(const edm::Event& iEvent,const edm::EventSetup& iSetup);
      void init(edm::InputTag
               ,const edm::Event& iEvent
	       ,const edm::EventSetup& iSetup);
      void SelectDiMu(const edm::Event& iEvent, const edm::EventSetup& iSetup);
      void analyze(const edm::Event& iEvent,const edm::EventSetup& iSetup);
      void LoadMuons(const edm::Event&,edm::InputTag);

   private:
      // ----------member data ---------------------------
      bool isValid;                                            //flag of validity
      bool MuonsAreLoaded;                                     //flag of tracks
      PrimaryInfo primary;                                     //Primary Vertex
      reco::BeamSpot beamSpot;                                 //BeamSpot
      edm::ESHandle<MagneticField> bFieldHandle;               //BField 
      edm::ESHandle<GlobalTrackingGeometry> globTkGeomHandle;  //Global tracker geometry
      edm::ESHandle<TransientTrackBuilder> theTTBHandle;       //Transient Track Builder
      std::vector<reco::Muon> theMuons;            //Input Muon Collection
      std::vector<reco::Muon> Muon1;               //First Leg Dimuon Candidate
      std::vector<reco::Muon> Muon2;               //Second Leg Dimuon Candidate
      std::vector<bool> IsKinematicFitReady;       //Flag for Dimuon Candidate
      std::map< unsigned int,  std::vector<RefCountedKinematicParticle> > uuParticle;
      std::map< unsigned int,  double > MuuMassFit;
      std::map< unsigned int,  reco::Particle::Point > MuuVertex;
      std::map< unsigned int,  reco::Vertex::CovarianceMatrix > MuuVertex_Cov;
      std::map< unsigned int,  reco::Particle::LorentzVector > MuuVertexP4;
      std::map< unsigned int,  reco::Vertex > MuuRecoVertex;
      std::map< unsigned int,  RefCountedKinematicVertex > MuuVertexInfo;
      std::map< unsigned int,  RefCountedKinematicParticle > MuuCandidate;
     
};
#endif
