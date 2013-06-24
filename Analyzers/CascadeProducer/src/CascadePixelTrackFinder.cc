// -*- C++ -*-
//
// Package:    CascadeProducer
// Class:      CascadePixelTrackFinder
// 
/**\class CascadePixelTrackFinder CascadePixelTrackFinder.cc Analyzers/CascadeProducer/src/CascadePixelTrackFinder.cc

 Description: find track with pixel information which is identified as  Xi/Omega particle 

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Juan Eduardo Ramirez Vargas
//         Created:  Mon Jun 29 17:24:25 CDT 2009
// $Id: CascadePixelTrackFinder.cc,v 1.3 2011/10/25 18:38:34 jramirez Exp $
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
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"       //track   collection
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"    //cascade collection
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"    //cascade daughter track
#include "DataFormats/Common/interface/TriggerResults.h"                 //triggers
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerRecord.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
//#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
//tools
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateClosestToPoint.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
//bfield
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "Analyzers/CascadeProducer/interface/CascadePixelTrackFinder.h"
#include "Analyzers/CascadeProducer/interface/TrackerInfo.h"
#include "Analyzers/CascadeProducer/interface/VertexRefit.h"
//MC
#include "Analyzers/CascadeProducer/interface/CascadeSimAnalyzer.h"
//kinemactic fitter
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
//
// constants, enums and typedefs
//
#include "Analyzers/CascadeProducer/interface/masses.h"

//
// static data member definitions
//

//
// constructors and destructor
//
//dummy constructor
CascadePixelTrackFinder::CascadePixelTrackFinder():
isValid(false),TracksAreLoaded(false)
,CascadesAreLoaded(false),hasPixelTrack(0)
{
}
//input in configurable parameters
CascadePixelTrackFinder::CascadePixelTrackFinder(const edm::ParameterSet& iConfig
                                                ,const edm::Event& iEvent
                                                ,const edm::EventSetup& iSetup):
isValid(false),TracksAreLoaded(false)
,CascadesAreLoaded(false),hasPixelTrack(0)
{
   //now do what ever initialization is needed

   // Get the track reco algorithm from the ParameterSet
   edm::InputTag tracksAlgo = iConfig.getParameter<edm::InputTag>("trackingAlgo");

   // Get the cas reco algorithm from the ParameterSet
   std::string CasAlgo = iConfig.getParameter<std::string>("CasAlgo");

   // Cas Collection Label
   std::string CasDecayName = iConfig.getUntrackedParameter<std::string>("CasDecayName","Cascade");
   
   LoadTracks(iEvent,tracksAlgo);
   if( !isValid ) {
    return;
   }

   LoadCascades(iEvent,CasAlgo,CasDecayName);
   if( !isValid ) {
    return;
   }

   isValid = true;
   analyze(iEvent,iSetup);

}//end constructor
//input directly as variables
CascadePixelTrackFinder::CascadePixelTrackFinder(edm::InputTag tracksAlgo
                                                ,std::string CasAlgo
                                                ,std::string CasDecayName
                                                ,const edm::Event& iEvent
                                                ,const edm::EventSetup& iSetup):
isValid(false),TracksAreLoaded(false)
,CascadesAreLoaded(false),hasPixelTrack(0)
{
   LoadTracks(iEvent,tracksAlgo);
   if( !isValid ) {
    return;
   }

   LoadCascades(iEvent,CasAlgo,CasDecayName);
   if( !isValid ) {
    return;
   }

   analyze(iEvent,iSetup);

}//end constructor2


CascadePixelTrackFinder::~CascadePixelTrackFinder()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
CascadePixelTrackFinder::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;



  //Load Tracker Geometry
   TrackerInfo tracker(iEvent,iSetup); 

   // Handles for B-field
   edm::ESHandle<MagneticField> bFieldHandle;
   iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);
//   const MagneticField* BField = bFieldHandle.product();
   
   // Handles for Tracker Geometry
   edm::ESHandle<GlobalTrackingGeometry> globTkGeomHandle;
   iSetup.get<GlobalTrackingGeometryRecord>().get(globTkGeomHandle);

   //Handles for Transient Track Builder 
   edm::ESHandle<TransientTrackBuilder> theTTBHandle;
   iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBHandle);
   
   
   //We have tracks and Cascades and are valids 
   if ( TracksAreLoaded 
     && CascadesAreLoaded
     && isValid){
   
   //Loop over cascades
   for(unsigned int casindex = 0; 
       casindex < theCascades.size(); 
       casindex++) {
     
     //get Xi vertex
     const GlobalPoint XiVtxPos(XiVertexDecay(casindex));

//pi trk
     reco::TrackRef piontrk      = Piontrack(casindex);

     reco::TransientTrack pionTtrk(piontrk, &(*bFieldHandle), globTkGeomHandle );
     math::XYZPoint innerhitPos= piontrk->innerPosition();
     double rho_trk = piontrk->innerPosition().rho();
     int pixhits    = piontrk->hitPattern().numberOfValidPixelHits();

     //Number of pixel hits on XI/Omega daughter track
     pionpixhits[casindex]=pixhits;

//vee
     reco::TrackRef protonveetrk =  ProtonVeetrack(casindex);
     reco::TrackRef pionveetrk   =  PionVeetrack(casindex);


// loop over track (pixel tracks)
// to find real Xi track
     std::vector<reco::TrackRef> XiTrkCand;
     std::vector<reco::TrackRef> XiTrkReject; //incompatible with vertex found
     std::vector<reco::TrackRef> XiTrkFail;   //pitrk and xitrk don't make a vertex
     std::vector<reco::TrackRef> XiTrkBad;    //pitrk and xitrk bad vertex
     std::vector<reco::TrackRef> XiTrkOverlapMother;      //overlap parent track
     std::vector<reco::TrackRef> XiTrkOverlapDaughter;    //overlap daughter track
     reco::TrackRef BestXiTrk;
     double bestcl=0;                                //for sorting tracks
     double bestip3d=-999;                           //for sorting tracks -
     std::vector<double> XiTrkip3d;                  //list of ip3d
     std::vector<double> XiTrkcl;                    //list of cl
     reco::TransientTrack pixtrk;

     //init vars
     XiTracks[casindex]=XiTrkCand;
     NumberXiTracks[casindex]=0;

     if (pixhits == 0) {
       bestip3d=999; // pixelless +
       for(unsigned int trkindex = 0; 
       trkindex < theTracks.size(); 
       trkindex++) {
	  if (piontrk      == TrkRefs[trkindex] ||
	      protonveetrk == TrkRefs[trkindex] ||
	      pionveetrk   == TrkRefs[trkindex]  ) continue;
	  if (TrkRefs[trkindex]->hitPattern().numberOfValidPixelHits()==0) continue;
          //get most outer hit position of pixel track
          const GlobalPoint trkPosOuter( TrkRefs[trkindex]->outerPosition().x(), 
                                         TrkRefs[trkindex]->outerPosition().y(), 
				         TrkRefs[trkindex]->outerPosition().z()
				       );	      
	  double rho_trkouter = trkPosOuter.perp();

          //Reject if overlaps with most inner pixeless track
	  if (rho_trkouter > rho_trk) continue;
          //
          //make transient trak (xi track candidate)
          reco::TransientTrack tmptt(TrkRefs[trkindex], &(*bFieldHandle), globTkGeomHandle);

          // Check if xitrk and piontrk make a vertex
          TransientVertex theXiRecoVertex;
          KalmanVertexFitter theKalmanFitter(true);
          std::vector<reco::TransientTrack> ListTracks;  //two tracks,  Xi track and pion trak Xi daughter
          ListTracks.push_back(tmptt);
          ListTracks.push_back(pionTtrk);
          theXiRecoVertex = theKalmanFitter.vertex(ListTracks);
          //Valid vertex
          if( !theXiRecoVertex.isValid()|| theXiRecoVertex.totalChiSquared() < 0. ) {
            XiTrkFail.push_back(TrkRefs[trkindex]);
            continue;
          }
          double vtxchi2 = theXiRecoVertex.totalChiSquared();
          double vtxndof = theXiRecoVertex.degreesOfFreedom();
          double vtxcl  = ChiSquaredProbability(vtxchi2,vtxndof);
          //at least 1 % CL
          if (vtxcl < 0.01){
            XiTrkBad.push_back(TrkRefs[trkindex]);
            continue;
          }
          //vertex should be after last hit of pixel trak
          //and before first hit of pion xi daughter trak
	  if (theXiRecoVertex.position().perp() < rho_trkouter){
	    XiTrkOverlapMother.push_back(TrkRefs[trkindex]);
	    continue;
	  }
	  if (theXiRecoVertex.position().perp() > rho_trk){
	    XiTrkOverlapDaughter.push_back(TrkRefs[trkindex]);
	    continue;
	  }

          //Impact parameter near Xi vtx decay
          TrajectoryStateClosestToPoint XiTraj =
                     tmptt.trajectoryStateClosestToPoint(XiVtxPos);
          if (!XiTraj.isValid()){
            XiTrkReject.push_back(TrkRefs[trkindex]);
            continue;
          }
	  //significance of dca for vertex
          double xitrkD0  = XiTraj.perigeeParameters().transverseImpactParameter();
          double xitrkD0E = XiTraj.perigeeError().transverseImpactParameterError();
          double xitrkDz  = XiTraj.perigeeParameters().longitudinalImpactParameter();
          double xitrkDzE = XiTraj.perigeeError().longitudinalImpactParameterError();
          double xitrkD3  = sqrt(xitrkD0*xitrkD0 + xitrkDz*xitrkDz);
          double xitrkD3E = sqrt(xitrkD0E*xitrkD0E*xitrkD0*xitrkD0 
	                  +      xitrkDzE*xitrkDzE*xitrkDz*xitrkDz)/xitrkD3;
          double xitrkip3d = xitrkD3/xitrkD3E;

          //best cl at begining
          if (vtxcl > bestcl){
             BestXiTrk = TrkRefs[trkindex];
             pixtrk   = tmptt;
             bestip3d = xitrkip3d;
	     bestcl   = vtxcl;
          }
          //list of impact parameters and the associated track
 	  XiTrkCand.push_back(TrkRefs[trkindex]);
 	  XiTrkip3d.push_back(xitrkip3d);
 	  XiTrkcl.push_back(vtxcl);
	  

       }//end traks loop
     }//endif pixhits == 0
     NumberXiTracks[casindex]=XiTrkCand.size();
     XiTracks[casindex]=XiTrkCand;
     XiImpactParameterList[casindex]=XiTrkip3d;
     XiCLList[casindex]=XiTrkcl;
     XiBestImpactParameter[casindex]=bestip3d;
     XiBestCL[casindex]=bestcl;
     BestXiTracks[casindex]=BestXiTrk;
     if (XiTrkCand.size()>0){
      std::cout << "CascadePixelTrackFinder Found " 
                << XiTrkCand.size()
                << " Xi trk cand, rejected  "
                <<  XiTrkReject.size()
                << " ,failed fit "
                <<  XiTrkFail.size()
                << " ,bad vertex (cl<1%) "
                <<  XiTrkBad.size()
                << " ,overlap parent "
                <<  XiTrkOverlapMother.size()
                << " ,overlap daughter "
                <<  XiTrkOverlapDaughter.size()
                << " ip3d = "
                <<  bestip3d
                << " cl = "
                <<  bestcl
		<< std::endl;
     }

 
   }//loop cascades
   }//We have tracks and Cascades valids 

}//end analyzer
void
CascadePixelTrackFinder::LoadTracks(const edm::Event& iEvent,
                                    edm::InputTag tracksLabel)
{
  TracksAreLoaded=true; //flag of tracks means that LoadTracks  runned
  // Handles for tracks
  edm::Handle<reco::TrackCollection> theTrackHandle;
  // Get the tracks from the event
  iEvent.getByLabel(tracksLabel, theTrackHandle);
  if( !theTrackHandle->size() ) {
    isValid = false;
    hasPixelTrack = 0;    
    return;
  }
  //Load Tracks into a vector of Tracks
  theTracks.insert( theTracks.end(), 
                    theTrackHandle->begin(), 
		    theTrackHandle->end() 
		    );
  //Load Refs		    
  for(unsigned int trkindex = 0;
    trkindex < theTracks.size();
    trkindex++) {
          reco::TrackRef tmptrk( theTrackHandle, trkindex );
          TrkRefs[trkindex]=tmptrk;
  }
  isValid=true;   

}
void
CascadePixelTrackFinder::LoadCascades(const edm::Event& iEvent, 
                                            std::string CasAlgo,
					    std::string CasDecayName)
{
  CascadesAreLoaded=true;
  //Handles for Cascades
  edm::Handle<reco::VertexCompositeCandidateCollection> theCasHandle;
  iEvent.getByLabel(CasAlgo, CasDecayName, theCasHandle);
  if( !theCasHandle->size() ) {
    isValid = false;
    hasPixelTrack = 0;
    return;
  }
  //Load  Cascades into a vector of composite candidates
  theCascades.insert( theCascades.end(), 
                     theCasHandle->begin(),
		     theCasHandle->end() 
		     );
  isValid=true;  
}
