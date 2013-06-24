// -*- C++ -*-
//
// Package:    DiMuonTools
// Class:      DiMuCandidates
//
/**\class DiMuCandidates DiMuCandidates.cc Analyzers/DiMuonTools/src/DiMuCandidates.cc

 Description: J/Psi Candidates selector

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Juan Eduardo Ramirez Vargas
//         Created:  Mon Jun 29 17:24:25 CDT 2009
// $Id: DiMuCandidates.cc,v 1.3 2013/05/17 17:23:23 jramirez Exp $
//
//

//user code 
#include "Analyzers/DiMuonTools/interface/DiMuCandidates.h"
#include "Analyzers/CascadeProducer/interface/KinematicFitDriver.h"
#include "Analyzers/CascadeProducer/interface/VertexRefit.h"
#include "Analyzers/CascadeProducer/interface/masses.h"

//bfield
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
//kinemactic fitter
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"


//
// constructors and destructor
//
//constructor dummy
DiMuCandidates::DiMuCandidates(const edm::Event& iEvent):
isValid(false)
,MuonsAreLoaded(false)
,primary(iEvent,"offlinePrimaryVertices")
{
}//end constructor dummy
//constructor dummy based on config file
DiMuCandidates::DiMuCandidates(const edm::ParameterSet& iConfig
                              ,const edm::Event& iEvent
                              ,const edm::EventSetup& iSetup):
isValid(false)
,MuonsAreLoaded(false)
,primary(iEvent,iConfig)
//,theMuonsLabel_(iConfig.getParameter<edm::InputTag>("MuonsLabel"    ))
{
}//end constructor dummy based on config file
//destructor
DiMuCandidates::~DiMuCandidates()
{
}//end destructor


//
// member functions
//

// ------------ method called to for each event  ------------
//Basic initialization
void
DiMuCandidates::init(const edm::Event& iEvent
                    ,const edm::EventSetup& iSetup)
{
   //get Beamspot
   edm::Handle<reco::BeamSpot> beamSpotHandle;
   iEvent.getByLabel("offlineBeamSpot", beamSpotHandle);
   if ( beamSpotHandle.isValid() ) beamSpot = *beamSpotHandle;

   //get BField
   iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);

}//end init
//Wrapper initialization to decouple from unpacking Bfield and Muons
void
DiMuCandidates::init(edm::InputTag theMuonsLabel_
                    ,const edm::Event& iEvent
                    ,const edm::EventSetup& iSetup)
{
  //Basic Initalization
  init(iEvent,iSetup);
  //get muons
  if (!MuonsAreLoaded) {
     LoadMuons(iEvent,theMuonsLabel_);
     if( !isValid ) {
        return;
     }
  }//endif Muons Are Loaded

  //Analyze Di Muon cases
  if (MuonsAreLoaded&&isValid){
     SelectDiMu(iEvent,iSetup);
     analyze(iEvent,iSetup);
  }
}//end  init2 the wrapper
void
DiMuCandidates::LoadMuons(const edm::Event& iEvent,
                                    edm::InputTag theMuonsLabel_)
{
 //--Get Muon Collection
  edm::Handle<reco::MuonCollection> allmuons;
  MuonsAreLoaded = iEvent.getByLabel(theMuonsLabel_,allmuons);
  if( !MuonsAreLoaded ) {
    isValid = false;
    return;
  }//end if


  if ( allmuons->size() > 0 ) {
    //Load Muons into a vector of Muons
    theMuons.insert( theMuons.end(),
                     allmuons->begin(),
                     allmuons->end()
                   );
    isValid = true;
  }else{
    isValid = false;//empty muon collection
  }
  
}//end LoadMuons
void
DiMuCandidates::SelectDiMu(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   //Loop over Muon1
   for(unsigned int iMuon = 0;
       iMuon < theMuons.size()-1;
       iMuon++) {
      //Loop over Muon2
      for(unsigned int jMuon = iMuon+1;
        jMuon < theMuons.size();
        jMuon++) {
        //require oposite muons
        if( theMuons[iMuon].charge()*theMuons[jMuon].charge() < 0 ){
          //save selection
          Muon1.push_back(theMuons[iMuon]);
          Muon2.push_back(theMuons[jMuon]);
          IsKinematicFitReady.push_back(false);
        }
      }//end for Muon2==jMuon
  }//end for Muon1==iMuon
}//end SelectDiMu
void
DiMuCandidates::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
//make sure Muon1 and Muon2 are same size
if ( Muon1.size() != Muon2.size() ){
  std::cout << "DiMuCandidates went wrong Muon1 and Muon2 diferent sizes";
  return;
}
for(unsigned int DiMuon = 0;
       DiMuon < Muon1.size();
       DiMuon++) {
  //require at least T-S or T-S-G
  if ( muGL(DiMuon,1)>=101 
    && muGL(DiMuon,2)>=101 ){
    
    //get ref tracks
    reco::TrackRef    iMuTrackRef = Muon1[DiMuon].innerTrack();
    reco::TrackRef    jMuTrackRef = Muon2[DiMuon].innerTrack();

    //--Re-fit dimuons to Vertex
    reco::TransientTrack iMuTT( iMuTrackRef, &(*bFieldHandle) );
    reco::TransientTrack jMuTT( jMuTrackRef, &(*bFieldHandle) );

    // Trajectory states to calculate DCA for the 2 tracks
    KinematicParticleFactoryFromTransientTrack pFactory;
    //FreeTrajectoryState iState = pFactory.particle( iMuTT, muon_mass_c, chi, ndf, muon_sig )->currentState().freeTrajectoryState();
    //FreeTrajectoryState jState = pFactory.particle( jMuTT, muon_mass_c, chi, ndf, muon_sig )->currentState().freeTrajectoryState();
    FreeTrajectoryState iState = iMuTT.impactPointTSCP().theState();
    FreeTrajectoryState jState = jMuTT.impactPointTSCP().theState();

    //--Measure distance between tracks at their closest approach
    double dca_uu;
    GlobalPoint cxPt;
    if (dcaisvalid(iState,jState,dca_uu,cxPt) ){   
       //--Cut on fiducial tracking volume [-120 cm < r < 120 cm  and -260 cm < Z < 260 cm]   r = sqrt( x^2 + y^2 )
       //--Pixel [-20 cm < r < 20 cm  and -60 cm < Z < 60 cm] to run vertex fit
       if( ( dca_uu >= 0. && dca_uu <= 5 ) && 
           ( sqrt( cxPt.x()*cxPt.x() + cxPt.y()*cxPt.y() ) < 120.) &&
                std::abs(cxPt.z()) < 300. ){

          float muon_sig   = muon_mass_c*1.e-6;
          float chi = 0., ndf = 0.;

          uuParticle[DiMuon].push_back(pFactory.particle( iMuTT, muon_mass_c, chi, ndf, muon_sig ));
          uuParticle[DiMuon].push_back(pFactory.particle( jMuTT, muon_mass_c, chi, ndf, muon_sig ));
          KinematicFitDriver uuRec( uuParticle[DiMuon], "MuMu" );
          IsKinematicFitReady[DiMuon]= uuRec.isValid();
          MuuMassFit[DiMuon]         = uuRec.mass();
          MuuVertex_Cov[DiMuon]      = uuRec.VertexCov();
          MuuVertex[DiMuon]          = uuRec.VertexDecay();
          MuuVertexP4[DiMuon]        = uuRec.P4();
          MuuRecoVertex[DiMuon]      = uuRec.Vertex();
          MuuVertexInfo[DiMuon]      = uuRec.RKVertex();
          MuuCandidate[DiMuon]       = uuRec.RKParent();
       }//if fidutial cut
    }//if dca is valid
  }//if at least global muons
}//end for DiMu
}//analyze
//Compute DistanceClosestApproach
bool
DiMuCandidates::dcaisvalid(FreeTrajectoryState &iState, FreeTrajectoryState &jState, double &dca_uu, GlobalPoint &cxPt)
{

//--Measure distance between tracks at their closest approach
ClosestApproachInRPhi cApp;           
dca_uu = -9999;
cApp.calculate(iState, jState);

if( !cApp.status() ) {
     std::cout << "uu bad status capp" << std::endl;
} else {
     dca_uu = std::abs( cApp.distance() );
     cxPt   = cApp.crossingPoint();
}
return cApp.status();

}//end of dcaisvalid
//muon id flag for Hector Control Histograms
int
DiMuCandidates::muGL(unsigned int DiMuon, int leg){

if (leg <1 || leg >2) return 0;
if (Muon1.size()-DiMuon <= 0) return 0;
reco::Muon iMuon;
if (leg == 1){
  iMuon = Muon1[DiMuon];
}else if (leg == 2){
  iMuon = Muon2[DiMuon];
}
int muSA=0, muGL=0, muTK=0;   
if( iMuon.isStandAloneMuon() )muSA = 1;     //--2^0
if( iMuon.isGlobalMuon()     )muGL = 10;    //--2^1
if( iMuon.isTrackerMuon()    )muTK = 100;   //--2^2
return  (muTK + muGL + muSA);

}//end muGL
//get Primary Refitted
reco::Vertex
DiMuCandidates::PrimaryRefitted(unsigned int DiMuon){
//exclude Muon1 and Muon2
     std::vector<reco::TrackRef> ExclusionList;
     ExclusionList.push_back(Muontrack(DiMuon,1) );
     ExclusionList.push_back(Muontrack(DiMuon,2) );
     VertexRefit RefitPrimary(primary.BestVertex()
                                      ,ExclusionList
                                      ,bFieldHandle
                                      ,beamSpot
                                      );
     return RefitPrimary.Refitted();

}
double
DiMuCandidates::delsig(unsigned int DiMuon){
//L/s
     typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
     typedef ROOT::Math::SVector<double, 3> SVector3;

//get refitted primary
     reco::Vertex refitVertexPrim = PrimaryRefitted(DiMuon);
//get DiMu vertex
    const GlobalPoint DiMuVtxPos(DiMuVertexDecay(DiMuon)
                              );


     SMatrixSym3D totalCov = refitVertexPrim.covariance()
                           + MuuVertex_Cov[DiMuon];

     SVector3 Xi_PrimarySep3D( DiMuVtxPos.x() - refitVertexPrim.x(),
                               DiMuVtxPos.y() - refitVertexPrim.y(),
                               DiMuVtxPos.z() - refitVertexPrim.z()
                             );

     return ROOT::Math::Similarity(totalCov, Xi_PrimarySep3D)>0?
            ROOT::Math::Dot(Xi_PrimarySep3D,Xi_PrimarySep3D)/sqrt(ROOT::Math::Similarity(totalCov, Xi_PrimarySep3D))
                                                               :-999;


}

