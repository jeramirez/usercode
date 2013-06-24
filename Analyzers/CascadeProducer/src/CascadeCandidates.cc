// -*- C++ -*-
//
// Package:    CascadeProducer
// Class:      CascadeCandidates
// 
/**\class CascadeCleaner CascadeCandidates.cc Analyzers/CascadeProducer/src/CascadeCandidates.cc

 Description: Skim Xi/Omega after fitter 

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Juan Eduardo Ramirez Vargas
//         Created:  Mon Jun 29 17:24:25 CDT 2009
// $Id: CascadeCandidates.cc,v 1.7 2012/07/18 18:46:09 jramirez Exp $
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
//#include "FWCore/Framework/interface/ESHandle.h"
//#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
//#include "TrackingTools/TransientTrack/interface/TransientTrack.h"       //track   collection
//#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"    //cascade collection
//#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"    //cascade daughter track
#include "DataFormats/Common/interface/TriggerResults.h"                 //triggers
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerRecord.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
//#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
//tools
//#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateClosestToPoint.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
//bfield
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
//#include "DataFormats/VertexReco/interface/VertexFwd.h"
//#include "DataFormats/VertexReco/interface/Vertex.h"

#include "Analyzers/CascadeProducer/interface/CascadeCandidates.h"
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
//constructor based on config file
CascadeCandidates::CascadeCandidates(const edm::ParameterSet& iConfig
                                                ,const edm::Event& iEvent
                                                ,const edm::EventSetup& iSetup):
isValid(false),TracksAreLoaded(false)
,CascadesAreLoaded(false),hasPixelTrack(0)
,primary(iEvent,iConfig)
{
   //now do what ever initialization is needed

   // Get the track reco algorithm from the ParameterSet
   edm::InputTag tracksAlgo = iConfig.getParameter<edm::InputTag>("trackingAlgo");

   // Get the cas reco algorithm from the ParameterSet
   std::string CasAlgo = iConfig.getParameter<std::string>("CasAlgo");

   // Cas Collection Label
   std::string CasDecayName = iConfig.getUntrackedParameter<std::string>("CasDecayName","Cascade");
   

   //get primary
   //primary(iEvent,iConfig);

   //basic init
   init(iEvent,iSetup);

   
   //get tracks
   LoadTracks(iEvent,tracksAlgo);
   if( !isValid ) {
    return;
   }

   //get cascades
   LoadCascades(iEvent,CasAlgo,CasDecayName);
   if( !isValid ) {
    return;
   }

   isValid = true;
   analyze(iEvent,iSetup);

}//end constructor
//constructor based on input tracks and cascades
CascadeCandidates::CascadeCandidates(edm::InputTag tracksAlgo
                                                ,std::string CasAlgo
                                                ,std::string CasDecayName
                                                ,const edm::Event& iEvent
                                                ,const edm::EventSetup& iSetup):
isValid(false),TracksAreLoaded(false)
,CascadesAreLoaded(false),hasPixelTrack(0)
,primary(iEvent,"offlinePrimaryVertices")
{

   //primary(iEvent,iConfig);  //get primary
   //basic init
   init(iEvent,iSetup);
   //get tracks
   LoadTracks(iEvent,tracksAlgo);
   if( !isValid ) {
    return;
   }

   //get cascades
   LoadCascades(iEvent,CasAlgo,CasDecayName);
   if( !isValid ) {
    return;
   }

   analyze(iEvent,iSetup);

}//end constructor2
//constructor dummy based on iConfig and iEvent
CascadeCandidates::CascadeCandidates(const edm::ParameterSet& iConfig
                                    ,const edm::Event& iEvent):
isValid(false),TracksAreLoaded(false)
,CascadesAreLoaded(false),hasPixelTrack(0)
,primary(iEvent,iConfig)
{
}//end constructor 3 dummy 1
//constructor dummy based iEvent needed for primary 
CascadeCandidates::CascadeCandidates(const edm::Event& iEvent):
isValid(false),TracksAreLoaded(false)
,CascadesAreLoaded(false),hasPixelTrack(0)
,primary(iEvent,"offlinePrimaryVertices")
{
}//end constructor 4 dummy 2
CascadeCandidates::~CascadeCandidates()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
//Basic initialization
void
CascadeCandidates::init(const edm::Event& iEvent
                       ,const edm::EventSetup& iSetup)
{
  //get Beamspot
    edm::Handle<reco::BeamSpot> beamSpotHandle;
   iEvent.getByLabel("offlineBeamSpot", beamSpotHandle);
   if ( beamSpotHandle.isValid() ) beamSpot = *beamSpotHandle;

   //get BField
   iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);
//   const MagneticField* BField = bFieldHandle.product();
   
   //get global geometry
   iSetup.get<GlobalTrackingGeometryRecord>().get(globTkGeomHandle);
   
   //get Transient Track Builder
   iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBHandle);
}
//Wrapper initialization to decouple from unpacking tracks and Cascades
void
CascadeCandidates::init(edm::InputTag tracksAlgo
                       ,std::string CasAlgo
                       ,std::string CasDecayName
                       ,const edm::Event& iEvent
                       ,const edm::EventSetup& iSetup)
{
  //Basic Initalization
  init(iEvent,iSetup);
  //get tracks
  if (!TracksAreLoaded) {
     LoadTracks(iEvent,tracksAlgo);
     if( !isValid ) {
        return;
     }
  }
  //get cascades
  if (!CascadesAreLoaded) {
     LoadCascades(iEvent,CasAlgo,CasDecayName);
     if( !isValid ) {
        return;
     }
  }
  //Analyze Xi tracks cases
  if (TracksAreLoaded&&CascadesAreLoaded)
     analyze(iEvent,iSetup);

  
} 
void
CascadeCandidates::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;



  //Load Tracker Geometry
   TrackerInfo tracker(iEvent,iSetup); 

   //We have tracks and Cascades and are valids 
   if ( TracksAreLoaded 
     && CascadesAreLoaded
     && isValid){
   
   //Loop over cascades
   for(unsigned int casindex = 0; 
       casindex < theCascades.size(); 
       casindex++) {
     
     //get Xi vertex
     const GlobalPoint XiVtxPos(theCascades[casindex].vx(), 
                                theCascades[casindex].vy(), 
				theCascades[casindex].vz()
				);	      
//pi trk
    //get Xi track daughter
    reco::TrackRef piontrk      = Piontrack(casindex);

    reco::TransientTrack pionTtrk(piontrk, &(*bFieldHandle), globTkGeomHandle );
    double rho_trk = RhoTrkDaughter(casindex);
    int pixhits    = piontrk->hitPattern().numberOfValidPixelHits();

     //Number of pixel hits on XI/Omega daughter track
     pionpixhits[casindex]      = pixhits;

//vee
    //get Xi grandchildren tracks of Lambda daughter
    reco::TrackRef pionveetrk   = PionVeetrack(casindex);
    reco::TrackRef protonveetrk = ProtonVeetrack(casindex);


// loop over track (pixel tracks)
// to find real Xi track
     std::vector<reco::TrackRef> XiTrkCand;
     std::vector<reco::TrackRef> XiTrkReject; //incompatible with vertex found
     std::vector<reco::TrackRef> XiTrkFail;   //pitrk and xitrk don't make a vertex
     std::vector<reco::TrackRef> XiTrkBad;    //pitrk and xitrk bad vertex
     std::vector<reco::TrackRef> XiTrkVtxs;   // /\pi vtx incompatible with Xi pi(daughter) vtx (separated)
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
          //make transient trak
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
	  //Vertices should be consistent to each other within 3sigma
	  reco::Vertex NewXiVertex = theXiRecoVertex;
          double XiVtxSigSep = delsig(casindex,NewXiVertex);
	  if (XiVtxSigSep > 3){
	    XiTrkVtxs.push_back(TrkRefs[trkindex]);
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
      LogInfo("CascadeCandidates") << "CascadeCandidates Found " 
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
                << " ,vtx separated "
                <<  XiTrkVtxs.size()
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
CascadeCandidates::LoadTracks(const edm::Event& iEvent,
                                    edm::InputTag tracksLabel)
{
  // Handles for tracks
  edm::Handle<reco::TrackCollection> theTrackHandle;
  // Get the tracks from the event
  TracksAreLoaded = iEvent.getByLabel(tracksLabel, theTrackHandle);
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
CascadeCandidates::LoadCascades(const edm::Event& iEvent, 
                                            std::string CasAlgo,
					    std::string CasDecayName)
{
  //Handles for Cascades
  edm::Handle<reco::VertexCompositeCandidateCollection> theCasHandle;
  CascadesAreLoaded=iEvent.getByLabel(CasAlgo, CasDecayName, theCasHandle);
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
reco::TransientTrack 
CascadeCandidates::GetXiTransientTrack(unsigned int casindex){
//get Xi vertex
    const GlobalPoint XiVtxPos(XiVertexDecay(casindex)
                              );

//get Xi momentum
    GlobalVector XiVectorP(MomentumAtDecay(casindex)
                           );
//get Xi track daughter
    reco::TrackRef piontrk      = Piontrack(casindex);

//get Bfield
    const MagneticField* BField = bFieldHandle.product();

//get pion transient track
    reco::TransientTrack pionTtrk(piontrk, &(*bFieldHandle), globTkGeomHandle );

//Build Xi free trajectory state
     FreeTrajectoryState XiFTS
        =FreeTrajectoryState(XiVtxPos,
	                     XiVectorP,
			     piontrk->charge(),
			     BField);
//use pion track at Xi decay to estimate errors
     //to be included in  Xi Transient Track
     TrajectoryStateClosestToPoint PionTrajAtXiDecay =
                     pionTtrk.trajectoryStateClosestToPoint(XiVtxPos);
     FreeTrajectoryState piFTS
         = PionTrajAtXiDecay.theState();
     XiFTS.setCartesianError(piFTS.cartesianError());
     
     return (*theTTBHandle).build(XiFTS);

}
double 
CascadeCandidates::SignificanceImpactParameter3DAtVertex(unsigned int casindex,reco::Vertex &PrmVtx){
     reco::TransientTrack XiTTrack = GetXiTransientTrack(casindex);
     return SignificanceAbsoluteImpactParameter3D(XiTTrack,PrmVtx);
}
double 
CascadeCandidates::SignificanceImpactParameter3DAtPrimary(unsigned int casindex){

     reco::TransientTrack XiTTrack = GetXiTransientTrack(casindex);
     reco::Vertex prmrefitted      = PrimaryRefitted(casindex);

     return SignificanceAbsoluteImpactParameter3D(XiTTrack,prmrefitted);
}
//Method to compute 3D absolute impact parameter
//based on information of Transversal Impact Parameter (2D) D0 and
//Longitudinal Impact Parameter (along Z)
// IP3d = sqrt(D0^2 + Dz^2)
//The error is propagated adding in quadrature:
// IP3de^2 = (d(IP3d)/d(D0))^2 * D0e^2 + (d(IP3d)/d(Dz))^2 * Dze^2
//simplified to:
// IP3de = sqrt(D0e^2*D0^2 + Dze^2*Dz^2)/IP3d
double 
CascadeCandidates::SignificanceAbsoluteImpactParameter3D(
                               reco::TransientTrack &TTrack,
                               reco::Vertex &Vertex) const{
   GlobalPoint PosVertex = GlobalPoint(Vertex.x(),Vertex.y(),Vertex.z());
   TrajectoryStateClosestToPoint TTraj =
                    TTrack.trajectoryStateClosestToPoint(PosVertex);
   double TrakD0  = TTraj.isValid()?TTraj.perigeeParameters().transverseImpactParameter()
                                   :-1000.;
   double TrakD0E = TTraj.isValid()?TTraj.perigeeError().transverseImpactParameterError()
                                   :-1000.;
   double TrakDz  = TTraj.isValid()?TTraj.perigeeParameters().longitudinalImpactParameter()
                                   :-1000.;
   double TrakDzE = TTraj.isValid()?TTraj.perigeeError().longitudinalImpactParameterError()
                                   :-1000.;
   double TrakD3  = sqrt(TrakD0*TrakD0 + TrakDz*TrakDz);
   double TrakD3E = sqrt(TrakD0E*TrakD0E*TrakD0*TrakD0
                  +      TrakDzE*TrakDzE*TrakDz*TrakDz)/TrakD3;

   return (TrakD3E>0?TrakD3/TrakD3E:0);
}
//Significance separation
//Xi respect its primary
double 
CascadeCandidates::delsig(unsigned int casindex){
//get refitted primary
     reco::Vertex refitVertexPrim = PrimaryRefitted(casindex);
     return delsig(casindex,refitVertexPrim);
}
//Significance separation
//Xi respect any given vertex
double 
CascadeCandidates::delsig(unsigned int casindex, reco::Vertex& refitVertexPrim){
//L/s
     typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
     typedef ROOT::Math::SVector<double, 3> SVector3;

//get Xi vertex
    const GlobalPoint XiVtxPos(XiVertexDecay(casindex)
                              );


     SMatrixSym3D totalCov = refitVertexPrim.covariance()
                           + theCascades[casindex].vertexCovariance();

//     GlobalVector Xiprmsep = XiVtxPos - PosVertex;

     SVector3 Xi_PrimarySep3D( XiVtxPos.x() - refitVertexPrim.x(),
                               XiVtxPos.y() - refitVertexPrim.y(),
                               XiVtxPos.z() - refitVertexPrim.z()
                             );

//     return ROOT::Math::Similarity(totalCov, Xi_PrimarySep3D)>0?Xiprmsep.mag2()/sqrt(ROOT::Math::Similarity(totalCov, Xi_PrimarySep3D))
     return ROOT::Math::Similarity(totalCov, Xi_PrimarySep3D)>0?
            ROOT::Math::Dot(Xi_PrimarySep3D,Xi_PrimarySep3D)/sqrt(ROOT::Math::Similarity(totalCov, Xi_PrimarySep3D))
                                                               :-999;


}
//get Primary Refitted
reco::Vertex 
CascadeCandidates::PrimaryRefitted(unsigned int casindex){
//exclude Xi pion and /\0 tracks
     std::vector<reco::TrackRef> ExclusionList;
     ExclusionList.push_back(Piontrack(casindex));
     ExclusionList.push_back(PionVeetrack(casindex));
     ExclusionList.push_back(ProtonVeetrack(casindex));
     VertexRefit RefitPrimary(primary.BestVertex()
                                      ,ExclusionList
                                      ,bFieldHandle
                                      ,beamSpot
                                      );
     return RefitPrimary.Refitted();

}
reco::VertexCollection::const_iterator
CascadeCandidates::ChoosePrimaryByCosAlpha(unsigned int casindex){
//get Xi vertex
    GlobalPoint XiVtxPos(XiVertexDecay(casindex)
                        );
//get Xi momentum
    GlobalVector XiVectorP(MomentumAtDecay(casindex)
                          );
//get primary
    return primary.HigherCosAlpha(XiVtxPos,XiVectorP);
}
