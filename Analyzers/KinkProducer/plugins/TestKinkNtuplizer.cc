// -*- C++ -*-
//
// Package:    KinkProducer
// Class:      TestKinkNtuplizer
// 
/**\class TestKinkNtuplizer TestKinkNutplizer.cc Analyzers/KinkProducer/plugins/TestKinkNtplizer.cc

 Description: analyzer for ntuple cascade compared with kink finder output 

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Juan Eduardo Ramirez Vargas
//         Created:  Mon Jun 29 17:24:25 CDT 2009
// $Id: TestKinkNtuplizer.cc,v 1.1 2012/07/17 21:06:18 jramirez Exp $
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
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"    //cascade colletcion
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

#include "Analyzers/KinkProducer/plugins/TestKinkNtuplizer.h"
#include "Analyzers/KinkProducer/interface/KinkFinder.h"
//#include "Analyzers/CascadeProducer/interface/CascadePixelTrackFinder.h"
#include "Analyzers/CascadeProducer/interface/CascadeCandidates.h"
#include "Analyzers/CascadeProducer/interface/KinematicFitDriver.h"
#include "Analyzers/CascadeProducer/interface/PrimaryInfo.h"
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
TestKinkNtuplizer::TestKinkNtuplizer(const edm::ParameterSet& iConfig):
pitrk_pt_cut_(iConfig.getUntrackedParameter<double>("pitrk_pt_cut",0))
,CasMass_(iConfig.getUntrackedParameter<double>("CasMass",1.32171))
,ximass(0),rho_cas(-1),los(-1),losnocorr(-1),probvxi(-1)
,rho_trk(-1),pionip3d(-1)
,probvee(-1),pionveeip3d(-1),protonveeip3d(-1)
,NewPrimaryProb(-1)
{
   //now do what ever initialization is needed

   // Get the track reco algorithm from the ParameterSet
   tracksAlgo = iConfig.getParameter<edm::InputTag>("trackingAlgo");

   // Get the vee reco algorithm from the ParameterSet
   VeeAlgo = iConfig.getParameter<std::string>("VeeAlgo");

   // Get the cas reco algorithm from the ParameterSet
   CasAlgo = iConfig.getParameter<std::string>("CasAlgo");

   // Cas Collection Label
   CasDecayName = iConfig.getUntrackedParameter<std::string>("CasDecayName","Cascade");
   

}


TestKinkNtuplizer::~TestKinkNtuplizer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
TestKinkNtuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;



  // Handles for tracks
  edm::Handle<reco::TrackCollection> theTrackHandle;
  // Get the tracks from the event
  iEvent.getByLabel(tracksAlgo, theTrackHandle);
  if( !theTrackHandle->size() ) {
    std::cout << "Warning<--EmptyTrackingCollection-->" 
                              << "No tracks in event!"
			      << std::endl;
    return;
  }
  //Load Tracks into a vector of Tracks
  std::vector<reco::Track> theTracks;
  theTracks.insert( theTracks.end(), 
                    theTrackHandle->begin(), 
		    theTrackHandle->end() 
		    );
   
 //Handles for Primary Vertex
  PrimaryInfo primary(iEvent,iConfig);
  PrimaryInfo primarybs(iEvent,"offlinePrimaryVerticesWithBS");
  nprm = primary.size();
  
 //Load Tracker Geometry
  TrackerInfo tracker(iEvent,iSetup); 

  //Handles for Cascades
  edm::Handle<reco::VertexCompositeCandidateCollection> theCasHandle;
  iEvent.getByLabel(CasAlgo, CasDecayName, theCasHandle);
  if( !theCasHandle->size() ) {
    return;
  }
  //Load  Cascades into a vector of composite candidates
  std::vector<reco::VertexCompositeCandidate> theCascades;
  theCascades.insert( theCascades.end(), 
                     theCasHandle->begin(),
		     theCasHandle->end() 
		     );
  
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
   
   //Handles for beamspot
   reco::BeamSpot beamSpot;
   edm::Handle<reco::BeamSpot> beamSpotHandle;
   iEvent.getByLabel("offlineBeamSpot", beamSpotHandle);
   if ( beamSpotHandle.isValid() ) beamSpot = *beamSpotHandle; 


   
    std::cout << "NXi="<<theCascades.size()
              <<" Ntrak="<< theTracks.size() 
              << " ";
//ed++ MC MC MC  MC MC MC
  CascadeSimAnalyzer SimData(iEvent,iSetup);
  if (SimData.SimTrackerIsValid())
    SimData.print(); 
//ed-- MC MC MC  MC MC MC


   CascadeCandidates MyXiTrack(tracksAlgo,CasAlgo,CasDecayName,iEvent,iSetup);
//kink test start
   std::vector<reco::TrackRef> TrkRefs;
   std::vector<reco::TrackRef> TrkVeto;
   for(unsigned int trkindex = 0; 
       trkindex < theTracks.size(); 
       trkindex++) {
          reco::TrackRef tmptrk( theTrackHandle, trkindex );
          TrkRefs.push_back(tmptrk);
   }
   KinkFinder KinkList(0.01,iSetup,TrkRefs,TrkVeto);
//kink test end

   //Loop over cascades
   for(unsigned int casindex = 0; 
       casindex < theCascades.size(); 
       casindex++) {
     ximass=theCascades[casindex].mass();
     std::cout << "#Ximass=" <<  ximass << std::endl;
     std::cout << "Number of Pixel Tracks =" << MyXiTrack.NumPixelXiTracks(casindex); //<< std::endl;     

     rho_cas = MyXiTrack.RhoVertexDecay(casindex);

     double  chi2 =  theCascades[casindex].vertexChi2();
     double  ndof =  theCascades[casindex].vertexNdof();
     probvxi  = ChiSquaredProbability(chi2,ndof);
     probvxi2 = ndof>1?ChiSquaredProbability(chi2,ndof-1):probvxi;
     //get Xi vertex
     const GlobalPoint XiVtxPos(theCascades[casindex].vx(), 
                                theCascades[casindex].vy(), 
				theCascades[casindex].vz()
				);	      
     //get Xi momentum
     GlobalVector XiVectorP(theCascades[casindex].momentum().x(),
                            theCascades[casindex].momentum().y(),
			    theCascades[casindex].momentum().z()
                           );
     //Xi vertex error
     //cov matrix to be used as Xi vtx error
     GlobalError XiCovError(
        theCascades[casindex].vertexCovariance(0,0),
        theCascades[casindex].vertexCovariance(1,0),theCascades[casindex].vertexCovariance(1,1),
        theCascades[casindex].vertexCovariance(2,0),theCascades[casindex].vertexCovariance(2,1),theCascades[casindex].vertexCovariance(2,2)
                           );
     
//pi trk
     reco::TrackRef piontrk =  (dynamic_cast<reco::RecoChargedCandidate *>
                                ( theCascades[casindex].daughter(0) )
				)->track();

     // cut on histogram  pion trk pt cut from conf file (def off=0)
     double piontrk_pt = piontrk->pt();
// pi trk pt_filter
     pionpt = piontrk_pt;
     if (piontrk_pt < pitrk_pt_cut_ ) {
       std::cout << " filtered" << std::endl;
       continue;
     }else{
       std::cout << " " << std::endl;
     }
     rho_trk = MyXiTrack.RhoTrkDaughter(casindex);
     reco::TransientTrack pionTtrk(piontrk, &(*bFieldHandle), globTkGeomHandle );

     math::XYZPoint innerhitPos= piontrk->innerPosition();
     trkhits = piontrk->numberOfValidHits();
     pixhits = piontrk->hitPattern().numberOfValidPixelHits();
     //check for hits before/after Xi Vertex
     const reco::HitPattern& hp = piontrk->hitPattern();
     hitsbeforev=0;
     hitsafterv=0;
     pxlHitsBe=0;
     pxlHitsAf=0;
     for (int i = 0; i < hp.numberOfHits(); i++) {
         uint32_t hit = hp.getHitPattern(i);
	 if (hp.trackerHitFilter(hit) && hp.validHitFilter(hit)) {
	    uint32_t subDet = hp.getSubStructure(hit);
            uint32_t layer = hp.getLayer(hit);
            double maxRZ = tracker.rangeRZmax(subDet,layer);
	    double rhoZ  = tracker.IsBarrel(subDet)?
	                                    rho_cas:
 		        std::abs(theCascades[casindex].vz());
	    if (rhoZ > maxRZ){
	      hitsbeforev++; 
	      if(tracker.IsPixel(subDet))pxlHitsBe++;
	    }else{
	      hitsafterv++;
	      if(tracker.IsPixel(subDet))pxlHitsAf++;
	    }


	 
         }//end if hit filter and valid

     }//for 
     std::cout << "Number of  Hits=" << piontrk->hitPattern().numberOfHits()
               << "\n Number of Valid Hits="<< trkhits
               << "\n Number of Valid Pixel Hits="<< piontrk->hitPattern().numberOfValidPixelHits()
               << "\n Number of Pixel Hits Before="<< pxlHitsBe
               << "\n Number of Pixel Hits After="<< pxlHitsAf
               << std::endl ;
 
//vee
     double chi2vee = theCascades[casindex].daughter(1)->vertexChi2();
     double ndofvee = theCascades[casindex].daughter(1)->vertexNdof();
     probvee = ChiSquaredProbability(chi2vee,ndofvee);
     
     reco::TrackRef protonveetrk =  MyXiTrack.ProtonVeetrack(casindex); 
     reco::TrackRef pionveetrk   =  MyXiTrack.PionVeetrack(casindex);
     reco::TransientTrack pionTveetrk(pionveetrk, &(*bFieldHandle) );
     reco::TransientTrack protonTveetrk(protonveetrk, &(*bFieldHandle) );
     
// loop over tracks (pixel tracks)
//  to find real Xi track
     std::list<reco::TrackRef> XiTrkCand;
     std::vector<reco::TrackRef> XiTrkReject; //incompatible with vertex found
     std::vector<reco::TrackRef> XiTrkRejTraj; //incompatible with vertex found
     std::vector<reco::TrackRef> XiTrkFail;   //pitrk and xitrk don't make a vertex
     std::vector<reco::TrackRef> XiTrkBad;    //pitrk and xitrk bad vertex
     std::vector<reco::TrackRef> XiTrkOverlap;//pitrk and xitrk overlap vertex track hit
     bestcl=0;                                //for sorting tracks
     bestip3d=-999;                           //for sorting tracks -
     reco::TransientTrack pixtrk;
     if (pixhits == 0) {
       bestip3d=999; // pixelless +
       for(unsigned int trkindex = 0; 
       trkindex < theTracks.size(); 
       trkindex++) {
          reco::TrackRef tmptrk( theTrackHandle, trkindex );
	  if (piontrk      == tmptrk ||
	      protonveetrk == tmptrk ||
	      pionveetrk   == tmptrk  ) continue;
	  if (tmptrk->hitPattern().numberOfValidPixelHits()==0) continue;
          const GlobalPoint trkPosOuter( tmptrk->outerPosition().x(), 
                                         tmptrk->outerPosition().y(), 
				         tmptrk->outerPosition().z()
				       );	      
	  double rho_trkouter = trkPosOuter.perp();
	  if (rho_trkouter > rho_trk) continue;
          //
          //test if pixel track is compatible with vertex already found
          reco::TransientTrack tmptt(tmptrk, &(*bFieldHandle), globTkGeomHandle);
	  TransientVertex theXiRecoVertex;
	  std::vector<reco::TransientTrack> FakeListTracks;  //is just 1 Xi track
	  FakeListTracks.push_back(tmptt);
	  KalmanVertexFitter theKalmanFitter(true);
	  theXiRecoVertex = theKalmanFitter.vertex(FakeListTracks,XiVtxPos,XiCovError);
          // Check Xi vertex validity 
          // and compatibility with pixel trk
          if( !theXiRecoVertex.isValid()|| theXiRecoVertex.totalChiSquared() < 0. ) {
            XiTrkReject.push_back(tmptrk);
            continue;
          }
          // Reject also if trk if far from current vertex
//          double ttchi2 = theXiRecoVertex.totalChiSquared();
//          double ttndof = theXiRecoVertex.degreesOfFreedom();
//          double ttcl  = ChiSquaredProbability(ttchi2,ttndof);
//          if (ttcl > 0.01){
//            XiTrkReject.push_back(tmptrk);
//            continue;
//          }
          //Impact parameter near Xi vtx decay
          TrajectoryStateClosestToPoint XiTraj =
		     tmptt.trajectoryStateClosestToPoint(XiVtxPos);
          if (!XiTraj.isValid()){
            XiTrkRejTraj.push_back(tmptrk);
            continue;
          }
          // Check if xitrk and piontrk make a vertex
          FakeListTracks.push_back(pionTtrk);
          theXiRecoVertex = theKalmanFitter.vertex(FakeListTracks);
          if( !theXiRecoVertex.isValid()|| theXiRecoVertex.totalChiSquared() < 0. ) {
            XiTrkFail.push_back(tmptrk);
            continue;
          }
          double vtxchi2 = theXiRecoVertex.totalChiSquared();
          double vtxndof = theXiRecoVertex.degreesOfFreedom();
          double vtxcl  = ChiSquaredProbability(vtxchi2,vtxndof);
          if (vtxcl < 0.01){
            XiTrkBad.push_back(tmptrk);
            continue;
          }
          //GlobalPoint theXiRecoVertex.position()
	  //vertex should be  between pixtrk end and pion trk beginning
	  if (theXiRecoVertex.position().perp() < rho_trkouter){
	    XiTrkOverlap.push_back(tmptrk);
	    continue;
	  }
	  if (theXiRecoVertex.position().perp() > rho_trk){
	    XiTrkOverlap.push_back(tmptrk);
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
//           //track within 30 sigma of Xi vertex decay
//           if (xitrkip3d > 30){
//             XiTrkReject.push_back(tmptrk);
//             continue;
//           }

          //best cl at begining
          if (vtxcl > bestcl){
             XiTrkCand.push_front(tmptrk);
             pixtrk   = tmptt;
             bestip3d = xitrkip3d;
	     bestcl   = vtxcl;
          }else{
	     XiTrkCand.push_back(tmptrk);
          }
//	  if (tmptrk->numberOfValidHits() == 
//	      tmptrk->hitPattern().numberOfValidPixelHits()){
//	      
//	  }//only pixel trks
	  

       }//end traks loop
     }//endif pixhits == 0
     hasXiTrk=XiTrkCand.size();
     hasXiTrk2=MyXiTrack.NumPixelXiTracks(casindex);
     bestip3d2=MyXiTrack.BestSignificanceImpactParameterAtPixelTrack(casindex);
     if (XiTrkCand.size()>0){
      std::vector<reco::TrackRef> KinkParents=KinkList.ListOfParents(piontrk);
      std::cout << "Kink Parents =" << KinkParents.size()  << std::endl;
      std::cout << "Found " 
                << XiTrkCand.size()
                << " Xi trk candidates, rejected incompatible with vertex "
                <<  XiTrkReject.size()
                << " ,rejected  because trajectory "
                <<  XiTrkRejTraj.size()
                << " ,fail vertex "
                <<  XiTrkFail.size()
                << " ,bad vertex "
                <<  XiTrkBad.size()
                << " ,overlap vertex "
                <<  XiTrkOverlap.size()
                << " ip3d = "
                <<  bestip3d
                << " cl = "
                <<  bestcl
		<< std::endl;
      if (bestip3d2 != bestip3d ){
          for (unsigned int iXi=0;iXi<MyXiTrack.size_i(casindex);iXi++){
             std::cout << "bestip3d[" << iXi <<"] again =" 
//	               << MyXiTrack.ImpactParameter_i(casindex,iXi) 
	               << MyXiTrack.ImpactParameterPixelTrack_i(casindex,iXi) 
		       << " cl[" << iXi << "]="
//		       << MyXiTrack.CL_i(casindex,iXi) 
		       << MyXiTrack.CLpixelTrack_i(casindex,iXi) 
		       << std::endl;
	  }
	  hasXiTrk=-XiTrkCand.size(); //tag as negative when mismatch
      }
     }

//lambda mass
     ParticleMass pion_mass     = pion_mass_c;
     ParticleMass proton_mass   = proton_mass_c;
     float pion_sigma   = pion_mass*1.e-6;       //to avoid singularities in the covariance matrix.
     float proton_sigma   = proton_mass*1.e-6;   //to avoid singularities in the covariance matrix.
     float chi = 0.;
     float ndf = 0.;
     KinematicParticleFactoryFromTransientTrack pFactory;     //Creating a KinematicParticleFactory
     std::vector<RefCountedKinematicParticle> LambdaParticle;
     LambdaParticle.push_back(pFactory.particle(pionTveetrk,pion_mass,chi,ndf,pion_sigma));
     LambdaParticle.push_back(pFactory.particle(protonTveetrk,proton_mass,chi,ndf,proton_sigma));
     KinematicFitDriver LambaRec(LambdaParticle,"Lambda");
     lbmass=LambaRec.mass();
     
//veto Ks
     //Check if Ks is reconstructed and veto it
     bool IsKs =false;
     ksmass=0.;
     std::vector<RefCountedKinematicParticle> KsParticle;
     KsParticle.push_back(pFactory.particle(pionTveetrk,pion_mass,chi,ndf,pion_sigma));
     KsParticle.push_back(pFactory.particle(protonTveetrk,pion_mass,chi,ndf,pion_sigma));
     KinematicParticleVertexFitter Ksfitter;
     RefCountedKinematicTree KsVertexFitTree;
     KsVertexFitTree = Ksfitter.fit(KsParticle);
     if (KsVertexFitTree->isValid()) {
        KsVertexFitTree->movePointerToTheTop();
        RefCountedKinematicParticle Ks_Particle = KsVertexFitTree->currentParticle();
        ksmass = Ks_Particle->currentState().mass();
        IsKs = std::abs(ksmass  - 0.497648) < 0.020; //Declare Ks a window with 20MeV of PDG
     }

//primary vertex
     //refit primary vertex (exclude Xi and /\0 tracks)	      
     reco::Vertex refitVertexPrim = MyXiTrack.PrimaryRefitted(casindex);
     GlobalPoint PosVertex = GlobalPoint(refitVertexPrim.x(),refitVertexPrim.y(),refitVertexPrim.z());

     
     NewPrimaryProb = ChiSquaredProbability(refitVertexPrim.chi2(),
		     refitVertexPrim.ndof());
		     
		     
		     
//Impact parameter 		     
     TrajectoryStateClosestToPoint PionTraj =
		     pionTtrk.trajectoryStateClosestToPoint(PosVertex);
     double PionD0  = PionTraj.perigeeParameters().transverseImpactParameter();
     double PionD0E = PionTraj.perigeeError().transverseImpactParameterError();
     double PionDz  = PionTraj.perigeeParameters().longitudinalImpactParameter();
     double PionDzE = PionTraj.perigeeError().longitudinalImpactParameterError();
     double PionD3  = sqrt(PionD0*PionD0 + PionDz*PionDz);
     double PionD3E = sqrt(PionD0E*PionD0E*PionD0*PionD0 
	            +      PionDzE*PionDzE*PionDz*PionDz)/PionD3;
     pionip3d = PionD3/PionD3E;
     std::pair<bool,Measurement1D> Pion3DIpPair = IPTools::absoluteImpactParameter3D(pionTtrk, refitVertexPrim);
     double PionD3Tool = -1000;
     if(Pion3DIpPair.first){
	     PionD3Tool = Pion3DIpPair.second.significance();
     }

 
     TrajectoryStateClosestToPoint PionVeeTraj =
		     pionTveetrk.trajectoryStateClosestToPoint(PosVertex);
     double PionVeeD0  = PionVeeTraj.perigeeParameters().transverseImpactParameter();
     double PionVeeD0E = PionVeeTraj.perigeeError().transverseImpactParameterError();
     double PionVeeDz  = PionVeeTraj.perigeeParameters().longitudinalImpactParameter();
     double PionVeeDzE = PionVeeTraj.perigeeError().longitudinalImpactParameterError();
     double PionVeeD3  = sqrt(PionVeeD0*PionVeeD0 + PionVeeDz*PionVeeDz);
     double PionVeeD3E = sqrt(PionVeeD0E*PionVeeD0E*PionVeeD0*PionVeeD0 
	               + PionVeeDzE*PionVeeDzE*PionVeeDz*PionVeeDz)/PionVeeD3;
     pionveeip3d = PionVeeD3/PionVeeD3E;		    
     std::pair<bool,Measurement1D> PionVee3DIpPair =
		     IPTools::absoluteImpactParameter3D(pionTveetrk, refitVertexPrim);
     double PionVeeD3Tool = -1000;
     if(PionVee3DIpPair.first){
	     PionVeeD3Tool = PionVee3DIpPair.second.significance();
     }

     TrajectoryStateClosestToPoint ProtonVeeTraj =
		     protonTveetrk.trajectoryStateClosestToPoint(PosVertex);
     double ProtonVeeD0  = ProtonVeeTraj.perigeeParameters().transverseImpactParameter();
     double ProtonVeeD0E = ProtonVeeTraj.perigeeError().transverseImpactParameterError();
     double ProtonVeeDz  = ProtonVeeTraj.perigeeParameters().longitudinalImpactParameter();
     double ProtonVeeDzE = ProtonVeeTraj.perigeeError().longitudinalImpactParameterError();
     double ProtonVeeD3  = sqrt(ProtonVeeD0*ProtonVeeD0 + ProtonVeeDz*ProtonVeeDz);
     double ProtonVeeD3E = sqrt(ProtonVeeD0E*ProtonVeeD0E*ProtonVeeD0*ProtonVeeD0 
	                 + ProtonVeeDzE*ProtonVeeDzE*ProtonVeeDz*ProtonVeeDz)/ProtonVeeD3;
     protonveeip3d = ProtonVeeD3/ProtonVeeD3E;		    
     std::pair<bool,Measurement1D> ProtonVee3DIpPair =
		     IPTools::absoluteImpactParameter3D(protonTveetrk, refitVertexPrim);
     double ProtonVeeD3Tool = -1000;          
     if(ProtonVee3DIpPair.first){
	     ProtonVeeD3Tool = ProtonVee3DIpPair.second.significance();
     }

//Xi transient track object (not Xitrack)
     reco::TransientTrack XiTTrack = MyXiTrack.GetXiTransientTrack(casindex);
//Significance Impact Parameter
     Xiip3d = MyXiTrack.SignificanceImpactParameter3DAtPrimary(casindex);		    

     

//L/s no corr                     
     typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
     typedef ROOT::Math::SVector<double, 3> SVector3;
     
     SMatrixSym3D totalCov = refitVertexPrim.covariance()
		           + theCascades[casindex].vertexCovariance();
     GlobalPoint XiVertexDecay  = GlobalPoint(theCascades[casindex].vx(),
		                              theCascades[casindex].vy(),
		                              theCascades[casindex].vz()
		                 );
     GlobalVector Xiprmsep = XiVertexDecay - PosVertex;
     double XiprmsepE = sqrt( totalCov(0,0) + totalCov(1,1) + totalCov(2,2) );

     losnocorr=Xiprmsep.mag()/XiprmsepE;  //L/sigma no correlation
     
//L/s
     los= MyXiTrack.delsig(casindex);

       tree_->Fill();
	      

   }//loop cascades

}//end analyzer


// ------------ method called once each job just before starting event loop  ------------
void 
TestKinkNtuplizer::beginJob()
{
  edm::Service<TFileService> fs;

  tree_ = fs->make<TTree>("CascadeTree","Xi/Omega ntuple");
  tree_->Branch("ximass",&ximass,"ximass/d");
  tree_->Branch("rho_cas",&rho_cas);
  tree_->Branch("rho_trk",&rho_trk);
  tree_->Branch("LoSnocorr",&losnocorr);
  tree_->Branch("LoS",&los);
  tree_->Branch("xiCL",&probvxi);
  tree_->Branch("xiCL2",&probvxi2);
  tree_->Branch("hasXiTrk",&hasXiTrk);
  tree_->Branch("hasXiTrk2",&hasXiTrk2);
  tree_->Branch("veemass",&lbmass);
  tree_->Branch("ksmass",&ksmass);
  tree_->Branch("veeCL",&probvee);
  tree_->Branch("prmCL",&NewPrimaryProb);
  tree_->Branch("pionip3d",&pionip3d);
  tree_->Branch("pionpt",&pionpt);
  tree_->Branch("piveeip3d",&pionveeip3d);
  tree_->Branch("prveeip3d",&protonveeip3d);
  tree_->Branch("xiip3d",&Xiip3d);
  tree_->Branch("bestip3d",&bestip3d);   //ip3d at decay of XiTrk
  tree_->Branch("bestip3d2",&bestip3d2); //ip3d at decay of XiTrk PixelTrack
  tree_->Branch("bestcl",&bestcl);       //vtx cl at decay of XiTrk
  tree_->Branch("trkhits",&trkhits);
  tree_->Branch("pixhits",&pixhits);
  tree_->Branch("hitsbef",&hitsbeforev);
  tree_->Branch("hitsaft",&hitsafterv);
  tree_->Branch("pixbefor",&pxlHitsBe);
  tree_->Branch("pixafter",&pxlHitsAf);
  tree_->Branch("nprm",&nprm);

}//begin job

// ------------ method called once each job just after ending the event loop  ------------
void 
TestKinkNtuplizer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(TestKinkNtuplizer);
