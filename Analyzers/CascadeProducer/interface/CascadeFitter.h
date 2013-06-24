// -*- C++ -*-
//
// Package:    CascadeProducer
// Class:      CascadeFitter
//
/**\class CascadeFitter CascadeFitter.h Analyzers/CascadeProducer/interface/CascadeFitter.h

 Description: Algorithm to find Cascades/Omegas

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Juan Eduardo Ramirez Vargas
//         Created:  Tue Jun 23 17:52:47 CDT 2009
// $Id: CascadeFitter.h,v 1.8 2011/06/16 19:40:04 jramirez Exp $
//
//

#ifndef RECOVERTEX__CASCADE_FITTER_H
#define RECOVERTEX__CASCADE_FITTER_H

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "RecoVertex/KinematicFitPrimitives/interface/KinematicRefittedTrackState.h"
#include "RecoVertex/VertexPrimitives/interface/CachingVertex.h"
#include "RecoVertex/KinematicFitPrimitives/interface/Matrices.h"


#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

//
// class declaration
//

class CascadeFitter {
  public:
     CascadeFitter(const edm::ParameterSet& iConfig,
	           const edm::Event& iEvent, 
		   const edm::EventSetup& iSetup);
    ~CascadeFitter();

  const reco::VertexCompositeCandidateCollection& getCascades() const;
  const reco::VertexCompositeCandidateCollection& getOmegas() const;

 private:
  //Input parameters
  std::string veesAlgo;         //Label of V0 collection
  std::string V0DecayName;      //Name of VO Collection Decay (Lambda/Kshort)
  edm::InputTag tracksAlgo;     //Label of TrackCollection
  bool fVerbose;                //verbose switch
  bool vtxFromfit;              //Use cas vertex from fit otherwise crosspoint
  bool massFromfit;             //Use cas inv mass from fit otherwise from geometry reconstruction
  double XiDCACut;              //Cut DCA tk and vee (Xi decay vertex)
  double LpiCut;                //Cut invariant mass Lambda-pi max value
  double LkaCut;                //Cut invariant mass Lambda-k  max value
  double LMassWidthCut;         //Cut Lambda mass window size around nominal value
  bool refitPrimary;            //flag to refit primary after excluding tracks
  double vtxSigCut2D;           //Cut on transversal separation of Xi vertex and primary
  double vtxSigCut3D;           //Cut on 3D separation of Xi vertex and primary
  double IP3DCut;               //Cut on 3D IP separation of Xi track daugthers and primary
	  
  // STL vector of VertexCompositeCandidate that will be filled with VertexCompositeCandidates by DoTheFit()
  reco::VertexCompositeCandidateCollection theXis;
  reco::VertexCompositeCandidateCollection theOmegas;

  // Method containing the fitter algorithm
  void DoTheFit(const edm::Event& iEvent, const edm::EventSetup& iSetup);

  // Method to select proton/pion traks from Lambda sorted first proton, second pion
  const std::pair<reco::TrackRef,reco::TrackRef> GetLambdaDaughters(
                          reco::VertexCompositeCandidate& TheLambda) const;
  const double SignificanceAbsoluteImpactParameter3D(
		                 reco::TransientTrack &TTrack,
		                 reco::Vertex &Vertex) const;

};
 
#endif
