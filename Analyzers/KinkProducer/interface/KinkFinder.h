#ifndef KINKPRODUCER_FINDER_H
#define KINKPRODUCER_FINDER_H
// -*- C++ -*-
//
// Package:    KinkProducer
// Class:      KinkFinder
//
/**\class KinkFinder KinkFinder.h Analyzers/KinkProducer/interface/KinkFinder.h

 Description: Gives a list of pair tracks (pixel,pixeless) matching with good vertex

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Juan Eduardo Ramirez Vargas
//         Created:  Mon Jun 29 17:24:25 CDT 2009
// $Id: KinkFinder.h,v 1.1 2012/07/17 21:06:17 jramirez Exp $
//
//

// system include files
#include <memory>
//references
#include "FWCore/Framework/interface/ESHandle.h"
// user include files
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"       //track   collection
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"    //cascade collection


//
// class declaration
//

class KinkFinder {
   public:
      explicit KinkFinder();
      explicit KinkFinder(const edm::ParameterSet& iConfig
		         ,const edm::EventSetup& iSetup
		         ,std::vector<reco::TrackRef> TrkRefs
		         ,std::vector<reco::TrackRef> VetoTrks
			 );
      explicit KinkFinder(double myvtxcut
		         ,const edm::EventSetup& iSetup
		         ,std::vector<reco::TrackRef> TrkRefs
		         ,std::vector<reco::TrackRef> VetoTrks
			 );
      ~KinkFinder();
      bool IsValid(){return !NumberOfParents.empty();};
      
      bool VertexIsValid(reco::TransientTrack PxTrack,reco::TransientTrack PxlessTrack, TransientVertex& theRecoVertex);
      bool VertexIsValid(std::vector<reco::TransientTrack> ListTracks, TransientVertex& theRecoVertex);
      bool VertexIsGood(TransientVertex& theRecoVertex);
      bool VertexIsDownStreamTrk(reco::TrackRef InTrack,TransientVertex& theRecoVertex);
      bool VertexIsUpStreamTrk(reco::TrackRef InTrack,TransientVertex& theRecoVertex);
      double RhoOuterPosition(reco::TrackRef InTrack);
      double RhoInnerPosition(reco::TrackRef InTrack);
      std::vector<reco::TrackRef> ListOfPixelTracks(std::vector<reco::TrackRef> InputTracks);
      std::vector<reco::TrackRef> ListOfPixellessTracks(std::vector<reco::TrackRef> InputTracks);
      std::vector<reco::TrackRef> ListOfTracks(std::vector<reco::TrackRef> InputTracks,
                                               std::vector<reco::TrackRef> VetoTracks  );
      std::vector<reco::TrackRef> ListOfParents(reco::TrackRef  TrkRefs);
      void MakeListOfTTracks(const edm::EventSetup& iSetup, std::vector<reco::TrackRef> TrkRefs );				       
      void analyze(const edm::EventSetup& iSetup, std::vector<reco::TrackRef> TrkRefs);

   private:
      // ----------member data ---------------------------
      bool isValid;                                            //flag of validity
      double vtxCUT;                                 //CL cut default is 1%
      std::map<reco::TrackRef,int > NumberOfParents; //Number of parent pix tracks linked to pixless track
      std::map<reco::TrackRef,std::vector<reco::TrackRef> > KinkCand; //list of parent pix tracks linked to given pixless track
      std::map<reco::TrackRef,int > NumberOfParentsElegibles;//Number of parent pix tracks linked to pixless track no overlap
      std::map<reco::TrackRef,std::vector<reco::TrackRef> > KinkCandElegibles; //list of parent pix tracks linked to given pixless track no overlap
      std::map<reco::TrackRef,reco::TransientTrack > ListOfTTracks; //Transient Tracks  corresponding to Track Ref
};
#endif
