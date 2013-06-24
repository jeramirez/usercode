// -*- C++ -*-
//
// Package:    CascadeProducer
// Class:      CascadeFitter
// 
/**\class CascadeFitter CascadeFitter.cc Analyzers/CascadeProducer/src/CascadeFitter.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Juan Eduardo Ramirez Vargas
//         Created:  Tue Jun 23 17:52:47 CDT 2009
// $Id: CascadeFitter.cc,v 1.14 2011/06/24 22:02:49 jramirez Exp $
//
//
//
//references
#include "FWCore/Framework/interface/ESHandle.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
//bfield and geometry
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
//MC
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
//kinemactic fitter
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

#include "Analyzers/CascadeProducer/interface/CascadeFitter.h"
#include "Analyzers/CascadeProducer/interface/ClosestApproachOnHelixLine.h"
#include "Analyzers/CascadeProducer/interface/PrimaryInfo.h"
#include "Analyzers/CascadeProducer/interface/VertexRefit.h"
#include "Analyzers/CascadeProducer/interface/KinematicFitDriver.h"
#include "Analyzers/CascadeProducer/interface/masses.h"

// Constructor
CascadeFitter::CascadeFitter(const edm::ParameterSet& iConfig,
		             const edm::Event& iEvent, 
			     const edm::EventSetup& iSetup) {

  // verbose swich
  fVerbose      = iConfig.getUntrackedParameter<bool>("verbose", false);
  // Get the vee reco algorithm from the ParameterSet
  veesAlgo      = iConfig.getParameter<std::string>("v0Algo");

  // Vee Collection Label
  V0DecayName   = iConfig.getUntrackedParameter<std::string>("v0DecayName","Lambda");

  // Get the track reco algorithm from the ParameterSet
  tracksAlgo    = iConfig.getParameter<edm::InputTag>("trackingAlgo");

  //bools for cuts on Kinematic fit
  vtxFromfit    = iConfig.getParameter<bool>(std::string("VtxFromFit"));
  massFromfit   = iConfig.getParameter<bool>(std::string("MassFromFit"));
  //cuts for producer
  XiDCACut      = iConfig.getParameter<double>(std::string("XiDCACut"));
  LpiCut        = iConfig.getParameter<double>(std::string("mLambdaPiCut"));
  LkaCut        = iConfig.getParameter<double>(std::string("mLambdaKCut"));
  LMassWidthCut = iConfig.getParameter<double>(std::string("LambdaMassWidthCut"));
  //Xi decay vertex away from primary
  vtxSigCut2D   = iConfig.getParameter<double>(std::string("vtxSignificance2DCut"));
  vtxSigCut3D   = iConfig.getParameter<double>(std::string("vtxSignificance3DCut"));
  IP3DCut       = iConfig.getParameter<double>(std::string("IP3DCut"));
  refitPrimary  = iConfig.getParameter<bool>(std::string("refitPrimary"));


  DoTheFit(iEvent, iSetup);


}//end constructor

//destructor (empty)
CascadeFitter::~CascadeFitter() {
}

// Method containing the fitter algorithm
void CascadeFitter::DoTheFit(const edm::Event& iEvent, 
                             const edm::EventSetup& iSetup) {

  // Handles for tracks
  edm::Handle<reco::TrackCollection> theTrackHandle;
  // Get the tracks from the event
  iEvent.getByLabel(tracksAlgo, theTrackHandle);
  if( !theTrackHandle->size() ) {
    if (fVerbose)
       std::cout <<"Warning-->EmptyTrackingCollection-->" 
                 << "No tracks in event!"
	         << std::endl;
    return;//nothing to do
  }
  //Load Tracks into a vector of Tracks
  std::vector<reco::Track> theTracks;
  theTracks.insert( theTracks.end(), 
                    theTrackHandle->begin(), 
		    theTrackHandle->end() 
		    );

  // Handles for B-field, and tracker geometry
  edm::ESHandle<MagneticField> bFieldHandle;
  edm::ESHandle<GlobalTrackingGeometry> globTkGeomHandle;
  //Get B-field and tracker geometry from event setup
  iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);
  iSetup.get<GlobalTrackingGeometryRecord>().get(globTkGeomHandle);

  const MagneticField* BField = bFieldHandle.product();

  //Handles for beamspot used for primary refit
  reco::BeamSpot beamSpot;
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByLabel("offlineBeamSpot", beamSpotHandle);
  if ( beamSpotHandle.isValid() ) beamSpot = *beamSpotHandle; 
   
  //Handles for Lambdas
  edm::Handle<reco::VertexCompositeCandidateCollection> theVeeHandle;
  iEvent.getByLabel(veesAlgo, V0DecayName, theVeeHandle);
  if( !theVeeHandle->size() ) {
//    throw cms::Exception("EmptyVeesCollection") 
    if (fVerbose)
       std::cout << "Warning-->EmptyVeesCollection-->"
                               << "No vees in event or not reconstructed but needed for Cascade Reconstruction"
  			       << std::endl;
    return;
  }

  //Load  Lambdas into a vector of composite candidates
  std::vector<reco::VertexCompositeCandidate> theLambdas;
  theLambdas.insert( theLambdas.end(), 
                     theVeeHandle->begin(),
		     theVeeHandle->end() 
		     );


   //Loop over tracks to convert to transient tracks
   std::vector<reco::TrackRef> theTrackRefs;
   std::vector<reco::TransientTrack> theTransTracks;
   for(unsigned int itrack = 0; 
       itrack < theTrackHandle->size(); 
       itrack++) {
         reco::TrackRef tmpRef( theTrackHandle, itrack );
         reco::TransientTrack tmpTk2( *tmpRef, &(*bFieldHandle), globTkGeomHandle );
         theTrackRefs.push_back( tmpRef );
	 theTransTracks.push_back( tmpTk2 );
   }
   
   //Get primary
   PrimaryInfo primary(iEvent,"offlinePrimaryVertices");
   
   //constants used for kinematic fit
   ParticleMass pion_mass   = pion_mass_c;
   ParticleMass kaon_mass   = kaon_mass_c;
   ParticleMass proton_mass = proton_mass_c;
   ParticleMass lambda_mass = lambda_mass_c;
   //to avoid singularities in the covariance matrix.
   float pion_sigma   = pion_mass*1.e-6;
   float kaon_sigma   = kaon_mass*1.e-6;
   float proton_sigma = proton_mass*1.e-6;
   float lambda_sigma = 0.000006;
   float chi = 0.;
   float ndf = 0.;

   //Loop over vees
   for(unsigned int veeindex = 0; 
       veeindex < theLambdas.size(); 
       veeindex++) {

      //basic rejection
      if ( theLambdas[veeindex].mass() < lambda_mass_c - LMassWidthCut 
        || theLambdas[veeindex].mass() > lambda_mass_c + LMassWidthCut
	 ) continue;                    //mass window cut;
      if ( theLambdas[veeindex].pdgId() !=  3122 
        && theLambdas[veeindex].pdgId() != -3122) continue; //make sure are /\0

      const std::pair<reco::TrackRef,reco::TrackRef> 
               theDaughterTracks      = GetLambdaDaughters(theLambdas[veeindex]);
      reco::TrackRef ProtonTrackRef   = theDaughterTracks.first;
      reco::TrackRef PionTrackRef     = theDaughterTracks.second;


      //get vee vertex
      const GlobalPoint VeeVtxPos(theLambdas[veeindex].vx(), 
                                  theLambdas[veeindex].vy(), 
				  theLambdas[veeindex].vz()
				 );	      
      //get vee momentum
      GlobalVector veeVectorP(theLambdas[veeindex].momentum().x(),
                              theLambdas[veeindex].momentum().y(),
			      theLambdas[veeindex].momentum().z()
                             );

      //get vertex error
      //Vee vertex cov matrix to be used as Xi vtx error
      GlobalError veecoverror(
        theLambdas[veeindex].vertexCovariance(0,0),
        theLambdas[veeindex].vertexCovariance(1,0),theLambdas[veeindex].vertexCovariance(1,1),
        theLambdas[veeindex].vertexCovariance(2,0),theLambdas[veeindex].vertexCovariance(2,1),theLambdas[veeindex].vertexCovariance(2,2)
                             );
      
      reco::TransientTrack protonTT(ProtonTrackRef, &(*bFieldHandle),globTkGeomHandle );
      reco::TransientTrack pionTT(PionTrackRef, &(*bFieldHandle),globTkGeomHandle );

      //Build vee trajetory state
      FreeTrajectoryState veeFTStrack =FreeTrajectoryState(VeeVtxPos,veeVectorP,0,BField);
	      
      //Creating a KinematicParticleFactory
      KinematicParticleFactoryFromTransientTrack pFactory;
      std::vector<RefCountedKinematicParticle> LambdaParticles;
      LambdaParticles.push_back(pFactory.particle(protonTT,proton_mass,chi,ndf,proton_sigma));
      LambdaParticles.push_back(pFactory.particle(pionTT,pion_mass,chi,ndf,pion_sigma));

      KinematicFitDriver LambdaRec(LambdaParticles,"vee");
	    if (vtxFromfit){
         // add mass constraint to the vee fit
         // otherwise cascade kinematic will fit fail in a buggy way
         LambdaRec.AddMassConstraint(lambda_mass,lambda_sigma);
         if (!LambdaRec.isValid()) continue;
         veecoverror = LambdaRec.error();
      }
      //Vee vertex cov matrix to be used as Xi vtx error
      AlgebraicSymMatrix33 veecov = veecoverror.matrix_new();

	
      //Now Loop over single tracks
      for(unsigned int itrack = 0; 
	  itrack < theTrackRefs.size();
	  itrack++) {
         reco::TrackRef tmpRef = theTrackRefs[itrack];
	 reco::TransientTrack pionTT2(tmpRef, &(*bFieldHandle),globTkGeomHandle );

         //trajectory state closest to point
         FreeTrajectoryState PionFTStrack = pionTT2.impactPointTSCP().theState();

	 if(!pionTT2.impactPointTSCP().isValid() ) continue;

	 GlobalVector trackVectorP = PionFTStrack.momentum();
	 GlobalPoint FTSposition   = PionFTStrack.position();
	 
	 ClosestApproachOnHelixLine cApp;
         GlobalPoint crosspoint;
	 std::pair<GlobalPoint, GlobalPoint> points;
	 double dca=-1;
	 double dl=-999;
         if (cApp.calculate(PionFTStrack.charge(),
	                    PionFTStrack.momentum(),
			    PionFTStrack.position(),
			    veeFTStrack.momentum(),
			    veeFTStrack.position(),
	                    PionFTStrack.parameters().magneticField()
	                   )
	     ){
	    crosspoint=cApp.crossingPoint();
	    points=cApp.points();
	    dca = cApp.distance();
	    dl = cApp.d_l();
	  }//cApp
          //Basic Selection of Tracks
	  if (
	      !(tmpRef == ProtonTrackRef)                      //track diff from vee daughter
	   && !(tmpRef == PionTrackRef)
	   &&  (tmpRef->charge()*ProtonTrackRef->charge() < 0) // /\0 pi- or anti/\0pi+
//This cut loose 1/3 of signal ?
//           &&  (dl > 0)                                        // /\ and pi point same hemisphere
           &&  (dca < XiDCACut)                                // dca < cut
           &&  (sqrt( crosspoint.x()*crosspoint.x() 
		    + crosspoint.y()*crosspoint.y() ) < 120.)   // fidutial volume
           &&  (abs(crosspoint.z()) < 300.)                     // fidutial volume
			  
  	     ){
	     
	      TrajectoryStateClosestToPoint trackTSCP = 
	          pionTT2.trajectoryStateClosestToPoint( crosspoint );
	      if( !trackTSCP.isValid()  ){
	         std::cout << "couln't get a pion state closest to crosspoint" << std::endl; 
	         continue;
	      }

	      //Preliminary Rejection way off candidates
	      GlobalVector ptrack = trackTSCP.momentum();
	      double pionESq  = ptrack.mag2() + piMassSquared;
	      double kaonESq  = ptrack.mag2() + kMassSquared;
	      double veeESq   = theLambdas[veeindex].momentum().mag2() + lambdaMassSquared;
	      double pionE    = sqrt(pionESq);
	      double kaonE    = sqrt(kaonESq);
	      double veeE     = sqrt(veeESq);
	      double totalEXi = pionE + veeE;
	      double totalEOm = kaonE + veeE;
	      double totalEXiSq = totalEXi*totalEXi;
	      double totalEOmSq = totalEOm*totalEOm;
	      double totalPSq   = ( ptrack + veeVectorP ).mag2();
	      double Xi_mass    = sqrt( totalEXiSq - totalPSq);
	      double Om_mass    = sqrt( totalEOmSq - totalPSq);

              		       
              if (Xi_mass > LpiCut && Om_mass > LkaCut)    //skip if mass to  big
	          continue;

	      //fit track  
	      //fit Xi pion daughter to vertex at crosspoint
	      //using /\0 vertex error 
	      TransientVertex theXiRecoVertex;
	      std::vector<reco::TransientTrack> FakeListTracks;  //is just 1 Xi track
	      FakeListTracks.push_back(pionTT2);
	      KalmanVertexFitter theKalmanFitter(true);
	      GlobalError VertexError(veecov);              //use Vee vertex error 
	      theXiRecoVertex = theKalmanFitter.vertex(FakeListTracks,crosspoint,VertexError);

              // Check Xi vertex validity
              if( !theXiRecoVertex.isValid()|| theXiRecoVertex.totalChiSquared() < 0. ) {
	        continue;
              }
	      
              // Create reco::Vertex object for use in creating the Candidate
              reco::Vertex theXiVertex = theXiRecoVertex;
	      

	      
	      reco::Vertex refitPrimVertex = *primary.BestVertex();	      
	      //refit primary vertex (exclude Xi and /\0 tracks)
	      if (refitPrimary){
	        std::vector<reco::TrackRef> ExclusionList;
	        ExclusionList.push_back(PionTrackRef);
	        ExclusionList.push_back(ProtonTrackRef);
	        ExclusionList.push_back(tmpRef);
	        VertexRefit RefitPrimary(primary.BestVertex()
			                ,ExclusionList
			                ,bFieldHandle
			                ,beamSpot
			                );
	        refitPrimVertex = RefitPrimary.Refitted();
              }


	      //Check separation from primary
	      typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
	      typedef ROOT::Math::SVector<double, 3> SVector3;
	      
	      SMatrixSym3D totalCov = refitPrimVertex.covariance() 
			            + theXiVertex.covariance();
	      SVector3 Xi_PrimarySep2D( theXiVertex.x() - refitPrimVertex.x(),
			                theXiVertex.y() - refitPrimVertex.y(),
			            0 );
	      SVector3 Xi_PrimarySep3D( theXiVertex.x() - refitPrimVertex.x(),
			                theXiVertex.y() - refitPrimVertex.y(),
			                theXiVertex.z() - refitPrimVertex.z()
			              );


              double rVtxMag      = ROOT::Math::Mag(Xi_PrimarySep2D);
              double sigmaRvtxMag = sqrt(ROOT::Math::Similarity(totalCov, Xi_PrimarySep2D)) / rVtxMag;
              double rVtxMag3D    = ROOT::Math::Mag(Xi_PrimarySep3D);
              double sigma3DvtxMag= sqrt(ROOT::Math::Similarity(totalCov,Xi_PrimarySep3D)) / rVtxMag3D;
	      double sigIP3dPion  = SignificanceAbsoluteImpactParameter3D(pionTT2,refitPrimVertex);
	      double sigIP3dPionVee  = SignificanceAbsoluteImpactParameter3D(pionTT,refitPrimVertex);
	      double sigIP3dProtonVee= SignificanceAbsoluteImpactParameter3D(protonTT,refitPrimVertex);
			     

	      //reject if significance2D < vtxSigCut2D
	      //reject if significance3D < vtxSigCut3D
	      //reject if IP3Dsignificance < IP3DCut
	      if ( rVtxMag / sigmaRvtxMag < vtxSigCut2D ) continue;
	      if ( rVtxMag3D / sigma3DvtxMag < vtxSigCut3D ) continue;
	      if ( sigIP3dPion < IP3DCut ) continue;
	      if ( sigIP3dPionVee < IP3DCut ) continue;
	      if ( sigIP3dProtonVee < IP3DCut ) continue;
	      
	      GlobalPoint XiVertexDecay(theXiVertex.x(), theXiVertex.y(), theXiVertex.z());
	      if( theXiRecoVertex.hasRefittedTracks() ) {
	         trackTSCP = 
	         theXiRecoVertex.refittedTracks().begin()->
			      trajectoryStateClosestToPoint(
			            XiVertexDecay );
		 
	      }else{
	         trackTSCP = pionTT2.trajectoryStateClosestToPoint( 
			            XiVertexDecay );
		      
	      }
	      if( !trackTSCP.isValid()  ){
	         std::cout << "couln't get a pion state closest to XiVertex fit" << std::endl; 
	         continue;
	      }
	      
	      //Start saving Xi info
	      ptrack = trackTSCP.momentum();
	      pionESq  = ptrack.mag2() + piMassSquared;
	      kaonESq  = ptrack.mag2() + kMassSquared;
	      pionE    = sqrt(pionESq);
	      kaonE    = sqrt(kaonESq);
	      totalEXi = pionE + veeE;
	      totalEOm = kaonE + veeE;
	      totalEXiSq = totalEXi*totalEXi;
	      totalEOmSq = totalEOm*totalEOm;
	      totalPSq   = ( ptrack + veeVectorP ).mag2();
	      Xi_mass    = sqrt( totalEXiSq - totalPSq);
	      Om_mass    = sqrt( totalEOmSq - totalPSq);
	      
              if (Xi_mass > LpiCut && Om_mass > LkaCut)    //again skip if mass to  big
	          continue;
              // Create momentum 4-vectors for the 2 candidate types
	      reco::Particle::LorentzVector pionP4(ptrack.x(), 
					           ptrack.y(), 
						   ptrack.z(), 
					           pionE);
	      reco::Particle::LorentzVector kaonP4(ptrack.x(), 
					           ptrack.y(), 
						   ptrack.z(), 
					           kaonE);
	      reco::Particle::LorentzVector xiP4(ptrack.x()+veeVectorP.x(), 
					           ptrack.y()+veeVectorP.y(), 
						   ptrack.z()+veeVectorP.z(), 
					           totalEXi);

	      reco::Particle::LorentzVector omP4(ptrack.x()+veeVectorP.x(), 
					         ptrack.y()+veeVectorP.y(), 
						 ptrack.z()+veeVectorP.z(), 
					         totalEOm);
	      
	      reco::Particle::Point Xivtx(XiVertexDecay.x(), XiVertexDecay.y(), XiVertexDecay.z());
	      reco::Particle::Point Omvtx(XiVertexDecay.x(), XiVertexDecay.y(), XiVertexDecay.z());

	      double XivtxChi2(theXiVertex.chi2());
	      double OmvtxChi2(theXiVertex.chi2());
	      
	      double XivtxNdof(theXiVertex.ndof());
	      double OmvtxNdof(theXiVertex.ndof());
	      
	      reco::Vertex::CovarianceMatrix XivtxCov(theXiVertex.covariance());	      
	      reco::Vertex::CovarianceMatrix OmvtxCov(theXiVertex.covariance());
	      
	      if (massFromfit || vtxFromfit){
	      //----Now Kinematic Fit Process----- 
	      //Do Cascade vertex fit (uses Lambda mass constraint)
              std::vector<RefCountedKinematicParticle> CascadeParticle;
	      CascadeParticle.push_back(pFactory.particle(pionTT2,pion_mass,chi,ndf,pion_sigma));
	      CascadeParticle.push_back(LambdaRec.RKParent());
	      KinematicFitDriver XiRec(CascadeParticle,"Xi");

	      //Do Omega vertex fit (uses Lambda mass constraint)
              std::vector<RefCountedKinematicParticle> OmegaParticle;
	      OmegaParticle.push_back(pFactory.particle(pionTT2,kaon_mass,chi,ndf,kaon_sigma));
	      OmegaParticle.push_back(LambdaRec.RKParent());
	      KinematicFitDriver OmegaRec(OmegaParticle,"Omega");
	      
	      //Reject if kinematic fit fail
	      if (!XiRec.isValid() && !OmegaRec.isValid()){ 
      		 continue;
	      }
	      

 	      if (massFromfit){ //replace P4 parameters from kinematic fit
					Xi_mass = XiRec.mass();
					xiP4    = XiRec.P4();
	        pionP4  = XiRec.P4FirstChild();
	        Om_mass = OmegaRec.mass();
					omP4    = OmegaRec.P4();
	        kaonP4  = OmegaRec.P4FirstChild();
	      }
	      if (vtxFromfit){//replace vertex parameters from kinematic fit
				  Xivtx     = XiRec.VertexDecay(); //from kinematic fit					    
	        XivtxChi2 = XiRec.chi2();        //from kinematic fit
	        XivtxNdof = XiRec.ndof();        //from kinematic fit
	        XivtxCov  = XiRec.VertexCov();   //from kinematic fit
					Omvtx     = OmegaRec.VertexDecay(); //from kinematic fit					    
	        OmvtxChi2 = OmegaRec.chi2();        //from kinematic fit
	        OmvtxNdof = OmegaRec.ndof();        //from kinematic fit
				  OmvtxCov  = OmegaRec.VertexCov();   //from kinematic fit
        }

              }//endif massFromfit or vtxFromfit
	      
	      // Create the VertexCompositeCandidate object that will be stored in the Event
	      reco::VertexCompositeCandidate* theXi    = 0;
	      reco::VertexCompositeCandidate* theOmega = 0;
	      theXi    = new reco::VertexCompositeCandidate(tmpRef->charge(), xiP4, 
			      Xivtx, XivtxCov, XivtxChi2, XivtxNdof);
	      theOmega = new reco::VertexCompositeCandidate(tmpRef->charge(), omP4,
			      Omvtx, OmvtxCov, OmvtxChi2, OmvtxNdof);
	      
	      // Create daughter candidates for the VertexCompositeCandidates
	      reco::RecoChargedCandidate thePionCand(tmpRef->charge(), pionP4, Xivtx);
	      thePionCand.setTrack(tmpRef);
	      reco::RecoChargedCandidate theKaonCand(tmpRef->charge(), kaonP4, Omvtx);
	      theKaonCand.setTrack(tmpRef);
	      
	      theXi->addDaughter(thePionCand);
	      theXi->addDaughter(theLambdas[veeindex]);
	      theXi->setPdgId(-(tmpRef->charge())*3312);
	      theOmega->addDaughter(theKaonCand);
	      theOmega->addDaughter(theLambdas[veeindex]);
	      theOmega->setPdgId(-(tmpRef->charge())*3334);
	      
	      if ((Xi_mass > 0) && (Xi_mass < LpiCut ) ) 
	        theXis.push_back( *theXi );
	      if ( (Om_mass > 0) && (Om_mass <LkaCut ) )
	        theOmegas.push_back( *theOmega );
	      
	      delete theXi;
	      delete theOmega;
	      theXi=theOmega=0;
          

	  }//Basic Track Selection
	}//loop over tracks  
   }//loop over vees


}//Do the fit

// Get methods
const reco::VertexCompositeCandidateCollection& CascadeFitter::getCascades() const {
  return theXis;
}

const reco::VertexCompositeCandidateCollection& CascadeFitter::getOmegas() const {
  return theOmegas;
}

const std::pair<reco::TrackRef,reco::TrackRef> CascadeFitter::GetLambdaDaughters(
                              reco::VertexCompositeCandidate& TheLambda) const {
   //unpack daughter track info
   std::vector<reco::RecoChargedCandidate> v0daughters;
   std::vector<reco::TrackRef> theDaughterTracks;
   for (unsigned int i = 0; 
        i < TheLambda.numberOfDaughters(); 
	i++) {
         v0daughters.push_back( *(dynamic_cast<reco::RecoChargedCandidate *> 
			       (TheLambda.daughter(i))) );
   }
   for(unsigned int j = 0; 
       j < v0daughters.size(); 
       j++) {
         theDaughterTracks.push_back(v0daughters[j].track());
   }

   bool isparticle(TheLambda.pdgId()>0);
   reco::TrackRef ProtonTrackRef;
   reco::TrackRef PionTrackRef;
   if ( (theDaughterTracks[0]->charge() > 0 && isparticle)     //is particle     and trk 0 + (proton +)
     || (theDaughterTracks[0]->charge() < 0 && !isparticle)    //is antiparticle and trk 0 - (proton -) 
         ){
         ProtonTrackRef = theDaughterTracks[0];
         PionTrackRef   = theDaughterTracks[1];              
   }else{
         ProtonTrackRef = theDaughterTracks[1];
         PionTrackRef   = theDaughterTracks[0];              
   }
   
   return std::pair<reco::TrackRef,reco::TrackRef>(ProtonTrackRef,PionTrackRef);
}
//Method to compute 3D absolute impact parameter 
//based on information of Transversal Impact Parameter (2D) D0 and
//Longitudinal Impact Parameter (along Z)
// IP3d = sqrt(D0^2 + Dz^2)
//The error is propagated adding in quadrature:
// IP3de^2 = (d(IP3d)/d(D0))^2 * D0e^2 + (d(IP3d)/d(Dz))^2 * Dze^2 
//simplified to:
// IP3de = sqrt(D0e^2*D0^2 + Dze^2*Dz^2)/IP3d
const double CascadeFitter::SignificanceAbsoluteImpactParameter3D(
		                 reco::TransientTrack &TTrack,
		                 reco::Vertex &Vertex) const{
   GlobalPoint PosVertex = GlobalPoint(Vertex.x(),Vertex.y(),Vertex.z());
   TrajectoryStateClosestToPoint TTraj =
		    TTrack.trajectoryStateClosestToPoint(PosVertex);
   double TrakD0  = TTraj.perigeeParameters().transverseImpactParameter();
   double TrakD0E = TTraj.perigeeError().transverseImpactParameterError();
   double TrakDz  = TTraj.perigeeParameters().longitudinalImpactParameter();
   double TrakDzE = TTraj.perigeeError().longitudinalImpactParameterError();
   double TrakD3  = sqrt(TrakD0*TrakD0 + TrakDz*TrakDz);
   double TrakD3E = sqrt(TrakD0E*TrakD0E*TrakD0*TrakD0 
	          +      TrakDzE*TrakDzE*TrakDz*TrakDz)/TrakD3;

   return (TrakD3E>0?TrakD3/TrakD3E:0);
}
