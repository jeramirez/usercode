// -*- C++ -*-
//
// Package:    KinkProducer
// Class:      KinkFinder
// 
/**\class KinkFinder KinkFinder.cc Analyzers/KinkProducer/src/KinkFinder.cc

 Description: find pair track with pixel and track pixeless consistent with a common vertex  

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Juan Eduardo Ramirez Vargas
//         Created:  Mon Jun 29 17:24:25 CDT 2009
// $Id: KinkFinder.cc,v 1.1 2012/07/17 21:06:17 jramirez Exp $
//
//


// system include files
#include <memory>

// user include files
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
//references
#include "FWCore/Framework/interface/ESHandle.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
//tools
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateClosestToPoint.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
//bfield
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "Analyzers/KinkProducer/interface/KinkFinder.h"
//kinemactic fitter
//#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
//#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
//#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
//#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
//dummy constructor
KinkFinder::KinkFinder():
isValid(false)
,vtxCUT(0.01)
{
}
//input in configurable parameters
KinkFinder::KinkFinder(const edm::ParameterSet& iConfig
                      ,const edm::EventSetup& iSetup
		      ,std::vector<reco::TrackRef> TrkRefs
		      ,std::vector<reco::TrackRef> VetoTrks
		      ):
isValid(false)
,vtxCUT(iConfig.getParameter<double>("kinkvtxcut"))
{
   //now do what ever initialization is needed

   analyze(iSetup,ListOfTracks(TrkRefs,VetoTrks));

}//end constructor
//input directly as variables
KinkFinder::KinkFinder(double myvtxcut
                      ,const edm::EventSetup& iSetup
		      ,std::vector<reco::TrackRef> TrkRefs
		      ,std::vector<reco::TrackRef> VetoTrks
		      ):
isValid(false)
,vtxCUT(myvtxcut)
{
   analyze(iSetup,ListOfTracks(TrkRefs,VetoTrks));

}//end constructor2


KinkFinder::~KinkFinder()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
KinkFinder::analyze(const edm::EventSetup& iSetup, std::vector<reco::TrackRef> TrkRefs)
{
   using namespace edm;


//KinkFinder start
   std::vector<reco::TrackRef> PxTracks    = ListOfPixelTracks(TrkRefs);
   std::vector<reco::TrackRef> PxlessTracks= ListOfPixellessTracks(TrkRefs);
   MakeListOfTTracks(iSetup,TrkRefs);
   //Loop over pixeless tracks (espected to be less)
   for( std::vector<reco::TrackRef>::iterator iter=PxlessTracks.begin();
     iter != PxlessTracks.end();
     iter++       
   ){
      reco::TrackRef PxlessTrack      = *iter;
      reco::TransientTrack pixlesstt  = ListOfTTracks[PxlessTrack];
      std::vector<reco::TrackRef> KinkParentCand;
      std::vector<reco::TrackRef> KinkParentCandElegible;
      for( std::vector<reco::TrackRef>::iterator itpx=PxTracks.begin();
        itpx != PxTracks.end();
        itpx++       
      ){
	 reco::TrackRef PxTrack      = *itpx;
	 //reject if overlaps outerpos with inner poss
	 bool isElegible=RhoInnerPosition(PxlessTrack) > RhoOuterPosition(PxTrack);
	 reco::TransientTrack pixtt = ListOfTTracks[PxTrack];
	 TransientVertex theRecoVertex;
	 if ( VertexIsValid(pixtt,pixlesstt,theRecoVertex)
	   && VertexIsGood(theRecoVertex)
	   && VertexIsDownStreamTrk(PxTrack,theRecoVertex)
	   && VertexIsUpStreamTrk(PxlessTrack,theRecoVertex)
	  ){
	     KinkParentCand.push_back(PxTrack);
	     if (isElegible)
	        KinkParentCandElegible.push_back(PxTrack);
	  }//if pass cond
	  
	  
      }//px trk loop
      NumberOfParents[PxlessTrack]=KinkParentCand.size();
      KinkCand[PxlessTrack]=KinkParentCand;
      NumberOfParentsElegibles[PxlessTrack]=KinkParentCandElegible.size();
      KinkCandElegibles[PxlessTrack]=KinkParentCandElegible;

   }//pixless trk loop
    
//KinkFinder end



}//end analyzer

std::vector<reco::TrackRef>
KinkFinder::ListOfPixelTracks(std::vector<reco::TrackRef> InputTracks)
{
  std::vector<reco::TrackRef> tmplist;
  for(unsigned int trkindex = 0;
       trkindex < InputTracks.size();
       trkindex++) {
       if (InputTracks[trkindex]->hitPattern().numberOfValidPixelHits()>0) 
          tmplist.push_back(InputTracks[trkindex]);
  }
  return tmplist;
}
std::vector<reco::TrackRef>
KinkFinder::ListOfPixellessTracks(std::vector<reco::TrackRef> InputTracks)
{
  std::vector<reco::TrackRef> tmplist;
  for(unsigned int trkindex = 0;
       trkindex < InputTracks.size();
       trkindex++) {
       if ( InputTracks[trkindex]->hitPattern().numberOfValidPixelHits()==0 ) 
          tmplist.push_back(InputTracks[trkindex]);
  }
  return tmplist;
}
std::vector<reco::TrackRef>
KinkFinder::ListOfTracks(std::vector<reco::TrackRef> InputTracks,
                               std::vector<reco::TrackRef> VetoTracks  )
{
  std::vector<reco::TrackRef> tmplist;
  for( std::vector<reco::TrackRef>::iterator iter=InputTracks.begin();
       iter != InputTracks.end();
       iter++       
     ){
     bool veto=false;
     
     for (std::vector<reco::TrackRef>::iterator vetoer=VetoTracks.begin();
          vetoer != VetoTracks.end();
	  vetoer++
         ){
	 if (*vetoer == *iter ) {
	   veto =true;
	   break;
	 } 
     }//loop veto track list

     tmplist.push_back(*iter);

  }//loop input track list
  return tmplist;
}
double
KinkFinder::RhoInnerPosition(reco::TrackRef InTrack)
{
 return InTrack->innerPosition().rho();
}
double
KinkFinder::RhoOuterPosition(reco::TrackRef InTrack)
{
 return InTrack->outerPosition().rho();
}
bool
KinkFinder::VertexIsValid(reco::TransientTrack PxTrack,reco::TransientTrack PxlessTrack,TransientVertex& theRecoVertex)
{
   std::vector<reco::TransientTrack> ListTracks;  //two tracks,  pix track and pixless trak as daughter
   ListTracks.push_back(PxTrack);
   ListTracks.push_back(PxlessTrack);
   return VertexIsValid(ListTracks,theRecoVertex);   
}
bool
KinkFinder::VertexIsValid(std::vector<reco::TransientTrack> ListTracks,TransientVertex& theRecoVertex)
{
   bool valid=false;
   // Check if pixtrk and pixlesstrk make a vertex
   
   KalmanVertexFitter theKalmanFitter(true);
   theRecoVertex = theKalmanFitter.vertex(ListTracks);
   
   //Valid vertex
   if( !theRecoVertex.isValid()|| theRecoVertex.totalChiSquared() < 0. ) {
      valid = false;
   }else{
      valid =true;
   }
   
   
   return valid;
}
bool
KinkFinder::VertexIsGood(TransientVertex& theRecoVertex)
{
  bool IsGood=false;
  double vtxchi2 = theRecoVertex.totalChiSquared();
  double vtxndof = theRecoVertex.degreesOfFreedom();
  double vtxcl   = ChiSquaredProbability(vtxchi2,vtxndof);
  //at least 1 % CL
  if (vtxcl > vtxCUT){
      IsGood=true;
  }
  return IsGood;
   
}
bool
KinkFinder::VertexIsDownStreamTrk(reco::TrackRef InTrack,TransientVertex& theRecoVertex)
{
  //vertex should be after last hit of pixel trak
  return theRecoVertex.position().perp() > RhoOuterPosition(InTrack);
}
bool
KinkFinder::VertexIsUpStreamTrk(reco::TrackRef InTrack,TransientVertex& theRecoVertex)
{
  //vertex should be before first hit of pixless trak
  return theRecoVertex.position().perp() < RhoInnerPosition(InTrack);
}
void
KinkFinder::MakeListOfTTracks(const edm::EventSetup& iSetup, std::vector<reco::TrackRef> TrkRefs )
{
   // Handles for B-field
   edm::ESHandle<MagneticField> bFieldHandle;
   iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);
   
   // Handles for Tracker Geometry
   edm::ESHandle<GlobalTrackingGeometry> globTkGeomHandle;
   iSetup.get<GlobalTrackingGeometryRecord>().get(globTkGeomHandle);

   ListOfTTracks.clear(); //Clean container 

   for( std::vector<reco::TrackRef>::iterator iter=TrkRefs.begin();
        iter != TrkRefs.end();
        iter++       
      ){
        reco::TrackRef Track      = *iter;
	reco::TransientTrack tt(Track, &(*bFieldHandle), globTkGeomHandle );
	ListOfTTracks[Track] = tt;
   }//loop trk ref
   
   return;

}
std::vector<reco::TrackRef> KinkFinder::ListOfParents(reco::TrackRef  TrkRefs){
   return KinkCand.find(TrkRefs)!=KinkCand.end()?
                               KinkCand[TrkRefs]:
			       std::vector<reco::TrackRef> ();
}
