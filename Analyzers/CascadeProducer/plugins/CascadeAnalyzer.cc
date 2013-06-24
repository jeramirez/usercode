// -*- C++ -*-
//
// Package:    CascadeProducer
// Class:      CascadeAnalyzer
// 
/**\class CascadeAnalyzer CascadeAnalyzer.cc Analyzers/CascadeProducer/src/CascadeAnalyzer.cc

 Description: analyzer for histograming cascade output 

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Juan Eduardo Ramirez Vargas
//         Created:  Mon Jun 29 17:24:25 CDT 2009
// $Id: CascadeAnalyzer.cc,v 1.1 2011/06/29 22:03:56 jramirez Exp $
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
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"       //track   collection
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"    //cascade colletcion
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"    //cascade daughter track
#include "DataFormats/Common/interface/TriggerResults.h"                 //triggers
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerRecord.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
//tools
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateClosestToPoint.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
//bfield
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "Analyzers/CascadeProducer/interface/CascadeAnalyzer.h"
#include "Analyzers/CascadeProducer/interface/PrimaryInfo.h"
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
CascadeAnalyzer::CascadeAnalyzer(const edm::ParameterSet& iConfig):
masshisto(iConfig.getUntrackedParameter<int>("nbins",100),
          iConfig.getUntrackedParameter<double>("xmin",1.27),
	  iConfig.getUntrackedParameter<double>("xmax",1.37))
,massrejhisto(iConfig.getUntrackedParameter<int>("nbins",100),
          iConfig.getUntrackedParameter<double>("xmin",1.27),
	  iConfig.getUntrackedParameter<double>("xmax",1.37))
,veemasshisto(100,1.06,1.16)
,chisqhisto(101,0,100)
,clhisto(100,0.0,1.0)
,pitrk_pt_cut_(iConfig.getUntrackedParameter<double>("pitrk_pt_cut",0))
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


CascadeAnalyzer::~CascadeAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
CascadeAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;



  // Handles for tracks
  edm::Handle<reco::TrackCollection> theTrackHandle;
  // Get the tracks from the event
  iEvent.getByLabel(tracksAlgo, theTrackHandle);
  if( !theTrackHandle->size() ) {
//    throw cms::Exception("EmptyTrackingCollection") 
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
   
  //Handles for Std Vees (Lambdas)
  bool ThereAreVees =false; 
  edm::Handle<reco::VertexCompositeCandidateCollection> theVeeHandle;
  iEvent.getByLabel(VeeAlgo, "Lambda", theVeeHandle);
  if( theVeeHandle->size() ) {
    ThereAreVees =true;
  }
  std::vector<reco::VertexCompositeCandidate> theVees;
  if (ThereAreVees){
    //Load  Vees into a vector of composite candidates
  theVees.insert( theVees.end(), 
                  theVeeHandle->begin(),
		  theVeeHandle->end() 
		  );
  }
 //Handles for Primary Vertex
  PrimaryInfo primary(iEvent,iConfig);
  PrimaryInfo primarybs(iEvent,"offlinePrimaryVerticesWithBS");

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
   
   //Handles for beamspot
   reco::BeamSpot beamSpot;
   edm::Handle<reco::BeamSpot> beamSpotHandle;
   iEvent.getByLabel("offlineBeamSpot", beamSpotHandle);
   if ( beamSpotHandle.isValid() ) beamSpot = *beamSpotHandle; 


   
     std::cout << theCascades.size()<< " "
            << theTracks.size() << " ";
//ed++ MC MC MC  MC MC MC
  CascadeSimAnalyzer SimData(iEvent,iSetup);
  if (SimData.SimTrackerIsValid())
    SimData.print(); 
//ed-- MC MC MC  MC MC MC


   if (ThereAreVees){
     //Loop over std. vees
     for(unsigned int veeindex = 0; 
         veeindex < theVees.size(); 
         veeindex++) {
        double veemass=theVees[veeindex].mass();
        veemasshisto.Fill(veemass,"#Lambda inv. mass");
       //refit mass using Kinematic Particle Factory
       reco::TrackRef trkvee0 = (dynamic_cast<reco::RecoChargedCandidate *>
		     (theVees[veeindex].daughter(0)) )->track();
       reco::TrackRef trkvee1 = (dynamic_cast<reco::RecoChargedCandidate *>
		     (theVees[veeindex].daughter(1)) )->track();
     
       bool isparticle(theVees[veeindex].pdgId()>0);
       reco::TrackRef ProtonTrackRef;
       reco::TrackRef PionTrackRef;
       if ( (trkvee0->charge() > 0 && isparticle)     //is particle     and trk 0 + (proton +)
       ||   (trkvee0->charge() < 0 && !isparticle)    //is antiparticle and trk 0 - (proton -)
          ){
         ProtonTrackRef = trkvee0;
         PionTrackRef   = trkvee1;
       }else{
         ProtonTrackRef = trkvee1;
         PionTrackRef   = trkvee0;
       }
       reco::TransientTrack protonTT(ProtonTrackRef, &(*bFieldHandle) );
       reco::TransientTrack pionTT(PionTrackRef, &(*bFieldHandle) );
	 
       //Creating a KinematicParticleFactory
       KinematicParticleFactoryFromTransientTrack pFactory;
       ParticleMass pion_mass   = pion_mass_c;
       ParticleMass proton_mass = proton_mass_c;
       //to avoid singularities in the covariance matrix.
       float pion_sigma   = pion_mass*1.e-6;
       float proton_sigma = proton_mass*1.e-6;
       float chi = 0.;
       float ndf = 0.;
       std::vector<RefCountedKinematicParticle> LambdaParticle;
       LambdaParticle.push_back(pFactory.particle(protonTT,proton_mass,chi,ndf,proton_sigma));
       LambdaParticle.push_back(pFactory.particle(pionTT,pion_mass,chi,ndf,pion_sigma));
       KinematicParticleVertexFitter veefitter;
       RefCountedKinematicTree VeeVertexFitTree;
       VeeVertexFitTree = veefitter.fit(LambdaParticle);
       if (VeeVertexFitTree->isValid()) {
	  VeeVertexFitTree->movePointerToTheTop();
	  RefCountedKinematicParticle Vee_Particle = VeeVertexFitTree->currentParticle();
	  double refitmass = Vee_Particle->currentState().mass();
	  double deltamass = refitmass  - veemass + 1.1;
          veemasshisto.Fill(refitmass,"#Lambda mass refit");
          veemasshisto.Fill(deltamass,"#Delta#Lambda mass");
       }
       //Check if Ks is reconstructed
       std::vector<RefCountedKinematicParticle> KsParticle;
       KsParticle.push_back(pFactory.particle(protonTT,pion_mass,chi,ndf,pion_sigma));
       KsParticle.push_back(pFactory.particle(pionTT,pion_mass,chi,ndf,pion_sigma));
       KinematicParticleVertexFitter Ksfitter;
       RefCountedKinematicTree KsVertexFitTree;
       KsVertexFitTree = Ksfitter.fit(KsParticle);
       if (KsVertexFitTree->isValid()) {
          KsVertexFitTree->movePointerToTheTop();
          RefCountedKinematicParticle Ks_Particle = KsVertexFitTree->currentParticle();
          double refitmass = Ks_Particle->currentState().mass();
          double deltamass = refitmass  - 0.497648 + 1.1;
          bool IsKs= std::abs(refitmass  - 0.497648) < 0.020;
          if (IsKs){veemasshisto.Fill(deltamass,"Ks mass");
	  }else{veemasshisto.Fill(deltamass,"Ks mass veto");
          }
          veemasshisto.Fill(deltamass,"#Delta Ks mass");
       }

     }//for over vees
   }//if There are vees
   //Loop over cascades
   for(unsigned int casindex = 0; 
       casindex < theCascades.size(); 
       casindex++) {
     double ximass=theCascades[casindex].mass();
     std::cout << "#Ximass=" <<  ximass << std::endl;
     
     double r_cas = sqrt(theCascades[casindex].vx()*theCascades[casindex].vx() +
                       theCascades[casindex].vy()*theCascades[casindex].vy() );     
     double chisq = theCascades[casindex].vertexNormalizedChi2();    
     double  chi2 =  theCascades[casindex].vertexChi2();
     double  ndof =  theCascades[casindex].vertexNdof();
     double probvxi= ChiSquaredProbability(chi2,ndof);
//pi trk
     reco::TrackRef piontrk =  (dynamic_cast<reco::RecoChargedCandidate *>
                                ( theCascades[casindex].daughter(0) )
				)->track();

     // cut on histogram  pion trk pt cut from conf file (def off=0)
     double piontrk_pt = piontrk->pt();
     if (piontrk_pt < pitrk_pt_cut_ ) continue;
     reco::TransientTrack pionTtrk(piontrk, &(*bFieldHandle) );
     math::XYZPoint innerhitPos= piontrk->innerPosition();
     double r_trk = sqrt(innerhitPos.x()*innerhitPos.x() +
                         innerhitPos.y()*innerhitPos.y() );
     double trkchi2 = piontrk->chi2();
     double trkndof = piontrk->ndof();
     int trkhits = piontrk->numberOfValidHits();
 
//vee
     double chi2vee = theCascades[casindex].daughter(1)->vertexChi2();
     double ndofvee = theCascades[casindex].daughter(1)->vertexNdof();
     double probvee = ChiSquaredProbability(chi2vee,ndofvee);
     reco::VertexCompositeCandidate* CasVee =
		     dynamic_cast<reco::VertexCompositeCandidate*>
		     (theCascades[casindex].daughter(1)
		     );
     
     reco::TrackRef protonveetrk =  (dynamic_cast<reco::RecoChargedCandidate *> 
		                     (CasVee->daughter(0) ) 
		                    )->track();
     reco::TrackRef pionveetrk =  (dynamic_cast<reco::RecoChargedCandidate *> 
		                     (CasVee->daughter(1) ) 
		                    )->track();
     reco::TransientTrack pionTveetrk(pionveetrk, &(*bFieldHandle) );
     reco::TransientTrack protonTveetrk(protonveetrk, &(*bFieldHandle) );
     

//veto Ks
     //Check if Ks is reconstructed and veto it
     bool IsKs =false;
     KinematicParticleFactoryFromTransientTrack pFactory;     //Creating a KinematicParticleFactory
     ParticleMass pion_mass   = pion_mass_c;
     float pion_sigma   = pion_mass*1.e-6;       //to avoid singularities in the covariance matrix.
     float chi = 0.;
     float ndf = 0.;
     std::vector<RefCountedKinematicParticle> KsParticle;
     KsParticle.push_back(pFactory.particle(pionTveetrk,pion_mass,chi,ndf,pion_sigma));
     KsParticle.push_back(pFactory.particle(protonTveetrk,pion_mass,chi,ndf,pion_sigma));
     KinematicParticleVertexFitter Ksfitter;
     RefCountedKinematicTree KsVertexFitTree;
     KsVertexFitTree = Ksfitter.fit(KsParticle);
     if (KsVertexFitTree->isValid()) {
        KsVertexFitTree->movePointerToTheTop();
        RefCountedKinematicParticle Ks_Particle = KsVertexFitTree->currentParticle();
        double ksmass = Ks_Particle->currentState().mass();
        IsKs = std::abs(ksmass  - 0.497648) < 0.020; //Declare Ks a window with 20MeV of PDG
     }

//primary vertex
     //refit primary vertex (exclude Xi and /\0 tracks)	      
     std::vector<reco::TrackRef> ExclusionList;	      
     ExclusionList.push_back(piontrk);	      
     ExclusionList.push_back(pionveetrk);	      
     ExclusionList.push_back(protonveetrk);	      
     VertexRefit RefitPrimary(primary.BestVertex()
			              ,ExclusionList
			              ,bFieldHandle
			              ,beamSpot
			              );
     reco::Vertex refitVertexPrim = RefitPrimary.Refitted();
     GlobalPoint PosVertex = GlobalPoint(refitVertexPrim.x(),refitVertexPrim.y(),refitVertexPrim.z());

     double NewPrimaryProb=-1;
     NewPrimaryProb = ChiSquaredProbability(refitVertexPrim.chi2(),
		     refitVertexPrim.ndof());
     TrajectoryStateClosestToPoint PionTraj =
		     pionTtrk.trajectoryStateClosestToPoint(PosVertex);
     double PionD0  = PionTraj.perigeeParameters().transverseImpactParameter();
     double PionD0E = PionTraj.perigeeError().transverseImpactParameterError();
     double PionDz  = PionTraj.perigeeParameters().longitudinalImpactParameter();
     double PionDzE = PionTraj.perigeeError().longitudinalImpactParameterError();
     double PionD3  = sqrt(PionD0*PionD0 + PionDz*PionDz);
     double PionD3E = sqrt(PionD0E*PionD0E*PionD0*PionD0 
	            +      PionDzE*PionDzE*PionDz*PionDz)/PionD3;
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
     std::pair<bool,Measurement1D> ProtonVee3DIpPair =
		     IPTools::absoluteImpactParameter3D(protonTveetrk, refitVertexPrim);
     double ProtonVeeD3Tool = -1000;          
     if(ProtonVee3DIpPair.first){
	     ProtonVeeD3Tool = ProtonVee3DIpPair.second.significance();
     }
                    
     typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
     SMatrixSym3D totalCov = refitVertexPrim.covariance()
		           + theCascades[casindex].vertexCovariance();
     GlobalPoint XiVertexDecay  = GlobalPoint(theCascades[casindex].vx(),
		                              theCascades[casindex].vy(),
		                              theCascades[casindex].vz()
		                 );
     GlobalVector Xiprmsep = XiVertexDecay - PosVertex;
     double XiprmsepE = sqrt( totalCov(0,0) + totalCov(1,1) + totalCov(2,2) );

          std::cout << "r_trk=" <<  r_trk
               << " r_cas="<< r_cas
               << " chi2="<< chi2
               << " ndof="<< ndof
               << " prob="<< probvxi
               << " trk_chisq="<< trkchi2
               << " trk_ndof="<< trkndof
               << " trk_nhits="<< trkhits
               << " prob_vee="<< probvee
	       << '\n'
	       << " prminfo chi2="<<primary.chi2()
	       << " prminfo ndof="<< primary.ndof()
	       << " prminfo prb="<< primary.cl()
	       << '\n'
	       << " Xi vertex cxx="<< theCascades[casindex].vertexCovariance(0,0)
	       << " Xi vertex cyy="<< theCascades[casindex].vertexCovariance(1,1)
	       << " Xi vertex czz="<< theCascades[casindex].vertexCovariance(2,2)
	       << '\n'
	       << " Prm vertex cxx="<< refitVertexPrim.covariance(0,0)
	       << " Prm vertex cyy="<< refitVertexPrim.covariance(1,1)
	       << " Prm vertex czz="<< refitVertexPrim.covariance(2,2)
	       << '\n'
	       << " prmbs chi2="<<primarybs.chi2()
	       << " prmbs ndof="<< primarybs.ndof()
	       << " prmbs prb="<< primarybs.cl()
	       << '\n'
	       << " prmvtx refit chi2="<< refitVertexPrim.chi2()
	       << " prmvtx refit ndof="<< refitVertexPrim.ndof()
	       << " prmvtx refit prb="<< NewPrimaryProb
              << std::endl;
	      
// loop over transversal R (i->R)
     for (int i=0;i<16;i++){
       std::stringstream ss;
       ss << i;
       std::string Xi_massID("#Xi_" + ss.str());              //all           r_cas > i
       std::string Xi_massID1("#Xi1_" + ss.str());            //all           r_cas > i
       std::string Xi_massIDpass0("pass0_#Xi_" + ss.str());   //all           r_cas > i
       std::string Xi_massID1cas("all1cas_#Xi1_" + ss.str());   //all         r_cas > i
       std::string Xi_massIDtype1("type1_#Xi_" + ss.str());   //              r_cas <  r_trk 
				// (transversal Xi vtx < pi trk most inner hits trans distance)
       std::string Xi_massIDtype2("type2_#Xi_" + ss.str());   //              r_cas >  r_trk 
       std::string Xi_massIDprob0("prob0_#Xi_" + ss.str());   //   XiCL>0.01 
       std::string Xi_massIDprob1("prob1_#Xi_" + ss.str());   //   XiCL>0.01 && r_cas < r_trk (type 1)
       std::string Xi_massIDprob2("prob2_#Xi_" + ss.str());   //   XiCL>0.01 && r_cas > r_trk (type 2)
       std::string Xi_massIDprob3("probv_#Xi_" + ss.str());   //   /\0CL>0.01 
       std::string Xi_massIDhits0("hits6_#Xi_" + ss.str());   //   trkhits > 5
       std::string Xi_massIDhits1("hits5_#Xi_" + ss.str());   //   trkhits > 4
       std::string Xi_massIDhits2("hits4_#Xi_" + ss.str());   //   trkhits > 3
       std::string Xi_massIDIP3D1("IP3D(#pi)_#Xi_" + ss.str()); //   IP3D_pi > 3sigma
       std::string Xi_massIDIP2D1("IP2D(#pi)_#Xi_" + ss.str()); //   IP2D_pi > 3sigma
       std::string Xi_massIDIP3D2("IP3D(#Lambda)_#Xi_" + ss.str());//IP3D_/\0 > 3sigma
       std::string Xi_massIDIP2D2("IP2D(#Lambda)_#Xi_" + ss.str());//IP2D_/\0 > 3sigma
       std::string Xi_massIDIP3D("IP3D_#Xi_" + ss.str());     //IP3D_Xi trks > 3sigma
       std::string Xi_massIDprim0("prm0_#Xi_" + ss.str());    //   prmCL>0.01
       std::string Xi_massIDprim1("prm1_#Xi_" + ss.str());    //   prmCL < prmrefitCL && prmrefitCL>0.01
       std::string Xi_massIDprim2("prmre_#Xi_" + ss.str());   //   prmrefitCL>0.01
       std::string Xi_massIDprim3("prmreIP3D(#pi)_#Xi_" + ss.str());//prmrefitCL>0.01 && IP3D_pi > 3sigma
       std::string Xi_massIDprim4("prmreIP2D(#pi)_#Xi_" + ss.str());//prmrefitCL>0.01 && IP2D_pi > 3sigma
       std::string Xi_massIDprim5("prmreIP3D(#Lambda)_#Xi_" + ss.str());//prmrefitCL>0.01 && IP3D_/\0 > 3sigma
       std::string Xi_massIDprim6("prmreIP2D(#Lambda)_#Xi_" + ss.str());//prmrefitCL>0.01 && IP3D_/\0 > 3sigma
       std::string Xi_massIDprim7("prmreIP3D_#Xi_" + ss.str());         //prmrefitCL>0.01 && IP3D_Xi trks > 3sigma
       std::string Xi_massIDprim8("prmrehits4_#Xi_" + ss.str());        //prmrefitCL>0.01 && trkhits > 3
       std::string Xi_massIDprim9("prmrehits5_#Xi_" + ss.str());        //prmrefitCL>0.01 && trkhits > 4
       std::string Xi_massIDprim10("prmrehits6_#Xi_" + ss.str());       //prmrefitCL>0.01 && trkhits > 5
       std::string Xi_massIDprim11("prmrehits4IP3D_#Xi_" + ss.str());   //prmrefitCL>0.01 && trkhits > 3 && IP3D_Xi trks > 3sigma
       std::string Xi_massIDprim12("prmrehits5IP3D_#Xi_" + ss.str());   //prmrefitCL>0.01 && trkhits > 4 && IP3D_Xi trks > 3sigma
       std::string Xi_massIDprim13("prmrehits6IP3D_#Xi_" + ss.str());   //prmrefitCL>0.01 && trkhits > 5 && IP3D_Xi trks > 3sigma
       std::string Xi_massIDprim14("prmrevxi_#Xi_" + ss.str());         //prmrefitCL>0.01 && XiCL>0.01
       std::string Xi_massIDprim15("prmrevxivee_#Xi_" + ss.str());      //prmrefitCL>0.01 && XiCL>0.01 && /\0CL>0.01
       std::string Xi_massIDprim16("prmrevxiveeIP3D_#Xi_" + ss.str());  //prmrefitCL>0.01 && XiCL>0.01 && /\0CL>0.01&& IP3D_Xi trks > 3sigma
       std::string Xi_massIDprim17("prmrevxiveels3_#Xi_" + ss.str());  //prmrefitCL>0.01 && XiCL>0.01 && /\0CL>0.01&& l > 3sigma
       std::string Xi_massIDprim18("prmrevxiveels4_#Xi_" + ss.str());  //prmrefitCL>0.01 && XiCL>0.01 && /\0CL>0.01&& l > 4sigma
       std::string Xi_massIDprim19("prmrevxiveels3hit3_#Xi_" + ss.str());  //prmrefitCL>0.01 && XiCL>0.01 && /\0CL>0.01&& l > 3sigma trkhits > 3
       std::string Xi_massIDprim20("prmrevxiveels3hit4_#Xi_" + ss.str());  //prmrefitCL>0.01 && XiCL>0.01 && /\0CL>0.01&& l > 3sigma trkhits > 4
       std::string Xi_massIDprim21("prmrevxiveels3hit5_#Xi_" + ss.str());  //prmrefitCL>0.01 && XiCL>0.01 && /\0CL>0.01&& l > 3sigma trkhits > 5
       std::string Xi_massIDprim22("prmrevxiveels3Rvtx_#Xi_" + ss.str());  //prmrefitCL>0.01 && XiCL>0.01 && /\0CL>0.01&& l > 3sigma r_cas < r_trk
       std::string Xi_massIDprim23("prmrevxiveeIPTool3_#Xi_" + ss.str());  //prmrefitCL>0.01 && XiCL>0.01 && /\0CL>0.01&& IP3DTools_Xi trks > 3sigma
       std::string Xi_massIDprim24("prmrevxiveels3IP3D_#Xi_" + ss.str());  //prmrefitCL>0.01 && XiCL>0.01 && /\0CL>0.01&& IP3D_Xi trks > 3sigma&& l > 3sigma
       std::string Xi_massIDprim25("prmrevxiveels3IPTool3_#Xi_" + ss.str());  //prmrefitCL>0.01 && XiCL>0.01 && /\0CL>0.01&& IP3DTools_Xi trks > 3sigma&& l > 3sigma
       std::string Xi_massIDprim26("prmrevxiveeIsKsls3_#Xi_" + ss.str());  //prmrefitCL>0.01 && XiCL>0.01 && /\0CL>0.01&& l > 3sigma
       std::string Xi_massIDprim27("prmrevxiveeIsKsls4_#Xi_" + ss.str());  //prmrefitCL>0.01 && XiCL>0.01 && /\0CL>0.01&& l > 4sigma
       std::string Xi_massIDprim28("prmrevxiveeIsKsls3hit3_#Xi_" + ss.str());  //prmrefitCL>0.01 && XiCL>0.01 && /\0CL>0.01&& l > 3sigma trkhits > 3
       std::string Xi_massIDprim29("prmrevxiveeIsKsls3hit4_#Xi_" + ss.str());  //prmrefitCL>0.01 && XiCL>0.01 && /\0CL>0.01&& l > 3sigma trkhits > 4
       std::string Xi_massIDprim30("prmrevxiveeIsKsls3hit5_#Xi_" + ss.str());  //prmrefitCL>0.01 && XiCL>0.01 && /\0CL>0.01&& l > 3sigma trkhits > 5
       std::string Xi_massIDprim31("prmrevxiveeIsKsls3Rvtx_#Xi_" + ss.str());  //prmrefitCL>0.01 && XiCL>0.01 && /\0CL>0.01&& l > 3sigma r_cas < r_trk
       std::string Xi_massIDprim32("prmrevxiveeIsKsIPTool3_#Xi_" + ss.str());  //prmrefitCL>0.01 && XiCL>0.01 && /\0CL>0.01&& IP3DTools_Xi trks > 3sigma
//NoKs
       std::string Xi_massIDprim33("prmrevxiveeNoKsls3_#Xi_" + ss.str());  //prmrefitCL>0.01 && XiCL>0.01 && /\0CL>0.01&& l > 3sigma
       std::string Xi_massIDprim34("prmrevxiveeNoKsls4_#Xi_" + ss.str());  //prmrefitCL>0.01 && XiCL>0.01 && /\0CL>0.01&& l > 4sigma
       std::string Xi_massIDprim35("prmrevxiveeNoKsls3hit3_#Xi_" + ss.str());  //prmrefitCL>0.01 && XiCL>0.01 && /\0CL>0.01&& l > 3sigma trkhits > 3
       std::string Xi_massIDprim36("prmrevxiveeNoKsls3hit4_#Xi_" + ss.str());  //prmrefitCL>0.01 && XiCL>0.01 && /\0CL>0.01&& l > 3sigma trkhits > 4
       std::string Xi_massIDprim37("prmrevxiveeNoKsls3hit5_#Xi_" + ss.str());  //prmrefitCL>0.01 && XiCL>0.01 && /\0CL>0.01&& l > 3sigma trkhits > 5
       std::string Xi_massIDprim38("prmrevxiveeNoKsls3Rvtx_#Xi_" + ss.str());  //prmrefitCL>0.01 && XiCL>0.01 && /\0CL>0.01&& l > 3sigma r_cas < r_trk
       std::string Xi_massIDprim39("prmrevxiveeNoKsIPTool3_#Xi_" + ss.str());  //prmrefitCL>0.01 && XiCL>0.01 && /\0CL>0.01&& IP3DTools_Xi trks > 3sigma
//NoksEnd
//CL xvi > 5%
       std::string Xi_massIDprim40("prmrevxi05veels3_#Xi_" + ss.str());  //prmrefitCL>0.01 && XiCL>0.01 && /\0CL>0.01&& l > 3sigma
       std::string Xi_massIDprim41("prmrevxi05veels4_#Xi_" + ss.str());  //prmrefitCL>0.01 && XiCL>0.01 && /\0CL>0.01&& l > 4sigma
       std::string Xi_massIDprim42("prmrevxi05veels3hit3_#Xi_" + ss.str());  //prmrefitCL>0.01 && XiCL>0.01 && /\0CL>0.01&& l > 3sigma trkhits > 3
       std::string Xi_massIDprim43("prmrevxi05veels3hit4_#Xi_" + ss.str());  //prmrefitCL>0.01 && XiCL>0.01 && /\0CL>0.01&& l > 3sigma trkhits > 4
       std::string Xi_massIDprim44("prmrevxi05veels3hit5_#Xi_" + ss.str());  //prmrefitCL>0.01 && XiCL>0.01 && /\0CL>0.01&& l > 3sigma trkhits > 5
       std::string Xi_massIDprim45("prmrevxi05veels3Rvtx_#Xi_" + ss.str());  //prmrefitCL>0.01 && XiCL>0.01 && /\0CL>0.01&& l > 3sigma r_cas < r_trk
       std::string Xi_massIDprim46("prmrevxi05veels3IPTool3_#Xi_" + ss.str());  //prmrefitCL>0.01 && XiCL>0.01 && /\0CL>0.01&& IP3DTools_Xi trks > 3sigma
//end CL xvi > 5%
//Pixelless Xi daughter trk
       std::string Xi_massIDprim47("prmrevxiveeNoPixells3_#Xi_" + ss.str());  //prmrefitCL>0.01 && XiCL>0.01 && /\0CL>0.01&& l > 3sigma
       std::string Xi_massIDprim48("prmrevxiveeNoPixells4_#Xi_" + ss.str());  //prmrefitCL>0.01 && XiCL>0.01 && /\0CL>0.01&& l > 4sigma
       std::string Xi_massIDprim49("prmrevxiveeNoPixells3hit3_#Xi_" + ss.str());  //prmrefitCL>0.01 && XiCL>0.01 && /\0CL>0.01&& l > 3sigma trkhits > 3
       std::string Xi_massIDprim50("prmrevxiveeNoPixells3hit4_#Xi_" + ss.str());  //prmrefitCL>0.01 && XiCL>0.01 && /\0CL>0.01&& l > 3sigma trkhits > 4
       std::string Xi_massIDprim51("prmrevxiveeNoPixells3hit5_#Xi_" + ss.str());  //prmrefitCL>0.01 && XiCL>0.01 && /\0CL>0.01&& l > 3sigma trkhits > 5
       std::string Xi_massIDprim52("prmrevxiveeNoPixells3Rvtx_#Xi_" + ss.str());  //prmrefitCL>0.01 && XiCL>0.01 && /\0CL>0.01&& l > 3sigma r_cas < r_trk
       std::string Xi_massIDprim53("prmrevxiveeNoPixells3IPTool3_#Xi_" + ss.str());  //prmrefitCL>0.01 && XiCL>0.01 && /\0CL>0.01&& l > 3sigma IP3DTools_Xi trks > 3sigma
//end Pixelless Xi daughter trk
       std::string Xi2ID("chisq_" + ss.str());
       
       std::string Xi_massID1prim2("prmre_#Xi1_" + ss.str());   //   prmrefitCL>0.01
       std::string Xi_massID1prim14("prmrevxi_#Xi1_" + ss.str());         //prmrefitCL>0.01 && XiCL>0.01
       std::string Xi_massID1prim15("prmrevxivee_#Xi1_" + ss.str());      //prmrefitCL>0.01 && XiCL>0.01 && /\0CL>0.01
       std::string Xi_massID1prim17("prmrevxiveels3_#Xi1_" + ss.str());  //prmrefitCL>0.01 && XiCL>0.01 && /\0CL>0.01&& l > 3sigma
       std::string Xi_massID1prim18("prmrevxiveels4_#Xi1_" + ss.str());  //prmrefitCL>0.01 && XiCL>0.01 && /\0CL>0.01&& l > 4sigma
       std::string Xi_massID1prim19("prmrevxiveels3hit3_#Xi1_" + ss.str());  //prmrefitCL>0.01 && XiCL>0.01 && /\0CL>0.01&& l > 3sigma trkhits > 3
       std::string Xi_massID1prim20("prmrevxiveels3hit4_#Xi1_" + ss.str());  //prmrefitCL>0.01 && XiCL>0.01 && /\0CL>0.01&& l > 3sigma trkhits > 4
       std::string Xi_massID1prim21("prmrevxiveels3hit5_#Xi1_" + ss.str());  //prmrefitCL>0.01 && XiCL>0.01 && /\0CL>0.01&& l > 3sigma trkhits > 5
       std::string Xi_massID1prim22("prmrevxiveels3Rvtx_#Xi1_" + ss.str());  //prmrefitCL>0.01 && XiCL>0.01 && /\0CL>0.01&& l > 3sigma r_cas < r_trk
//
//CL histos
//
       std::string Xi_CLIDprim17("prmrevxiveels3_CL#Xi_" + ss.str());  //prmrefitCL>0.01 && XiCL>0.01 && /\0CL>0.01&& l > 3sigma
       std::string Xi_CLIDprim18("prmrevxiveels4_CL#Xi_" + ss.str());  //prmrefitCL>0.01 && XiCL>0.01 && /\0CL>0.01&& l > 4sigma
       std::string Xi_CLIDprim19("prmrevxiveels3hit3_CL#Xi_" + ss.str());  //prmrefitCL>0.01 && XiCL>0.01 && /\0CL>0.01&& l > 3sigma trkhits > 3
       std::string Xi_CLIDprim20("prmrevxiveels3hit4_CL#Xi_" + ss.str());  //prmrefitCL>0.01 && XiCL>0.01 && /\0CL>0.01&& l > 3sigma trkhits > 4
       std::string Xi_CLIDprim21("prmrevxiveels3hit5_CL#Xi_" + ss.str());  //prmrefitCL>0.01 && XiCL>0.01 && /\0CL>0.01&& l > 3sigma trkhits > 5
       std::string Xi_CLIDprim22("prmrevxiveels3Rvtx_CL#Xi_" + ss.str());  //prmrefitCL>0.01 && XiCL>0.01 && /\0CL>0.01&& l > 3sigma r_cas < r_trk

       if (r_cas > i)
       {
	  masshisto.Fill(ximass,Xi_massIDpass0);
	  if (theCascades.size() == 1)masshisto.Fill(ximass,Xi_massID1cas);
	  if (probvxi > 0.01){
	     masshisto.Fill(ximass,Xi_massIDprob0);
          }
	  if (probvee > 0.01){
	     masshisto.Fill(ximass,Xi_massIDprob3);
          }
	  if (trkhits > 5){
             masshisto.Fill(ximass,Xi_massIDhits0);
	  }  
	  if (trkhits > 4){
             masshisto.Fill(ximass,Xi_massIDhits1);
	  }  
	  if (trkhits > 3){
             masshisto.Fill(ximass,Xi_massIDhits2);
	  }  
	  if ( primary.cl() > 0.01){
	     masshisto.Fill(ximass,Xi_massIDprim0);     
	  }
	  if (primary.cl() < NewPrimaryProb && NewPrimaryProb > 0.01){
	     masshisto.Fill(ximass,Xi_massIDprim1);     
	  }
	  if (NewPrimaryProb > 0.01 ){
	     masshisto.Fill(ximass,Xi_massIDprim2);
	     if (theCascades.size() == 1)masshisto.Fill(ximass,Xi_massID1prim2); 
	     if (probvxi > 0.01){
		masshisto.Fill(ximass,Xi_massIDprim14);
		if (theCascades.size() == 1)masshisto.Fill(ximass,Xi_massID1prim14);
		if (probvee > 0.01){
	           masshisto.Fill(ximass,Xi_massIDprim15);
	           if (theCascades.size() == 1)masshisto.Fill(ximass,Xi_massID1prim15);
		   if (Xiprmsep.mag() >  3*XiprmsepE) {
		      masshisto.Fill(ximass,Xi_massIDprim17);
                      if (std::abs(ximass - 1.3216)< 0.020 )clhisto.Fill(probvxi,Xi_CLIDprim17);   
		      if (theCascades.size() == 1)masshisto.Fill(ximass,Xi_massID1prim17);
		      if (IsKs){
			masshisto.Fill(ximass,Xi_massIDprim26);
		      }else{
			masshisto.Fill(ximass,Xi_massIDprim33);
		      }
		      if (probvxi>0.05)masshisto.Fill(ximass,Xi_massIDprim40);
		      if (piontrk->hitPattern().numberOfValidPixelHits()==0)masshisto.Fill(ximass,Xi_massIDprim47);
		      if (trkhits > 3){
			 masshisto.Fill(ximass,Xi_massIDprim19);
                         if (std::abs(ximass - 1.3216)< 0.020 )clhisto.Fill(probvxi,Xi_CLIDprim19);
		         if (theCascades.size() == 1)masshisto.Fill(ximass,Xi_massID1prim19);
		         if (IsKs){
			    masshisto.Fill(ximass,Xi_massIDprim28);
			 }else{
                            masshisto.Fill(ximass,Xi_massIDprim35);
                      	 }
			 if (probvxi>0.05)masshisto.Fill(ximass,Xi_massIDprim42);
		         if (piontrk->hitPattern().numberOfValidPixelHits()==0)masshisto.Fill(ximass,Xi_massIDprim49);
		      }   
		      if (trkhits > 4){
			 masshisto.Fill(ximass,Xi_massIDprim20);
			 if (std::abs(ximass - 1.3216)< 0.020 )clhisto.Fill(probvxi,Xi_CLIDprim20);
		         if (theCascades.size() == 1)masshisto.Fill(ximass,Xi_massID1prim20);
		         if (IsKs){
			    masshisto.Fill(ximass,Xi_massIDprim29);
			 }else{
			    masshisto.Fill(ximass,Xi_massIDprim36);
			 }
			 if (probvxi>0.05)masshisto.Fill(ximass,Xi_massIDprim43);
		         if (piontrk->hitPattern().numberOfValidPixelHits()==0)masshisto.Fill(ximass,Xi_massIDprim50);
		      }   
		      if (trkhits > 5){
			 masshisto.Fill(ximass,Xi_massIDprim21);
			 if (std::abs(ximass - 1.3216)< 0.020 )clhisto.Fill(probvxi,Xi_CLIDprim21);
		         if (theCascades.size() == 1)masshisto.Fill(ximass,Xi_massID1prim21);
		         if (IsKs){
			    masshisto.Fill(ximass,Xi_massIDprim30);
			 }else{
			    masshisto.Fill(ximass,Xi_massIDprim37);
			 }
			 if (probvxi>0.05)masshisto.Fill(ximass,Xi_massIDprim44);
		         if (piontrk->hitPattern().numberOfValidPixelHits()==0)masshisto.Fill(ximass,Xi_massIDprim51);
		      }
		      if (r_cas < r_trk) {
			 masshisto.Fill(ximass,Xi_massIDprim22);
			 if (std::abs(ximass - 1.3216)< 0.020 )clhisto.Fill(probvxi,Xi_CLIDprim22);
		         if (theCascades.size() == 1)masshisto.Fill(ximass,Xi_massID1prim22);
		         if (IsKs){
			    masshisto.Fill(ximass,Xi_massIDprim31);
			 }else{
			    masshisto.Fill(ximass,Xi_massIDprim38);
			 }
			 if (probvxi>0.05)masshisto.Fill(ximass,Xi_massIDprim45);
		         if (piontrk->hitPattern().numberOfValidPixelHits()==0)masshisto.Fill(ximass,Xi_massIDprim52);
		      }   
	              if ((PionVeeD3 > 3*PionVeeD3E)
	               && (ProtonVeeD3 > 3*ProtonVeeD3E)
	               && (PionD3 > 3*PionD3E) ){
	                   masshisto.Fill(ximass,Xi_massIDprim24);
		      }     
	              if ((PionVeeD3Tool > 3)
	               && (ProtonVeeD3Tool > 3)
	               && (PionD3Tool > 3) ){
	                   masshisto.Fill(ximass,Xi_massIDprim25);
			   if (probvxi>0.05)masshisto.Fill(ximass,Xi_massIDprim46);
			   if (piontrk->hitPattern().numberOfValidPixelHits()==0)masshisto.Fill(ximass,Xi_massIDprim53);
		      }     
		   }
		   if (Xiprmsep.mag() >  4*XiprmsepE) {
		      masshisto.Fill(ximass,Xi_massIDprim18);
		      if (std::abs(ximass - 1.3216)< 0.020 )clhisto.Fill(probvxi,Xi_CLIDprim18);
	              if (theCascades.size() == 1)masshisto.Fill(ximass,Xi_massID1prim18);
	              if (IsKs){
			 masshisto.Fill(ximass,Xi_massIDprim27);
		      }else{
			 masshisto.Fill(ximass,Xi_massIDprim34);
		      }
		      if (probvxi>0.05)masshisto.Fill(ximass,Xi_massIDprim41);
		      if (piontrk->hitPattern().numberOfValidPixelHits()==0)masshisto.Fill(ximass,Xi_massIDprim48);
		   }
	           if ((PionVeeD3Tool > 3)
	            && (ProtonVeeD3Tool > 3)
	               && (PionD3Tool > 3) ){
	                 masshisto.Fill(ximass,Xi_massIDprim23);
	                 if (IsKs){
			    masshisto.Fill(ximass,Xi_massIDprim32);
			 }else{
			    masshisto.Fill(ximass,Xi_massIDprim39);
			 }
		   }     
		}
	     }     
  	     if (trkhits > 3){
                masshisto.Fill(ximass,Xi_massIDprim8);
	        if ((PionVeeD3 > 3*PionVeeD3E)
	         && (ProtonVeeD3 > 3*ProtonVeeD3E)
	         && (PionD3 > 3*PionD3E) ){
		    masshisto.Fill(ximass,Xi_massIDprim11);
		}
	     }  
  	     if (trkhits > 4){
                masshisto.Fill(ximass,Xi_massIDprim9);
	        if ((PionVeeD3 > 3*PionVeeD3E)
	         && (ProtonVeeD3 > 3*ProtonVeeD3E)
	         && (PionD3 > 3*PionD3E) ){
		    masshisto.Fill(ximass,Xi_massIDprim12);
		}
	     }  
  	     if (trkhits > 5){
                masshisto.Fill(ximass,Xi_massIDprim10);
	        if ((PionVeeD3 > 3*PionVeeD3E)
	         && (ProtonVeeD3 > 3*ProtonVeeD3E)
	         && (PionD3 > 3*PionD3E) ){
		    masshisto.Fill(ximass,Xi_massIDprim13);
		}
	     }  
	  }
	  if (PionD3 > 3*PionD3E ){
	     masshisto.Fill(ximass,Xi_massIDIP3D1);     
	     if (NewPrimaryProb > 0.01 ){
	        masshisto.Fill(ximass,Xi_massIDprim3);     
	     }
	  }
	  if (PionD0 > 3*PionD0E ){
	     masshisto.Fill(ximass,Xi_massIDIP2D1);     
	     if (NewPrimaryProb > 0.01 ){
	        masshisto.Fill(ximass,Xi_massIDprim4);     
	     }
	  }
	  if ((PionVeeD3 > 3*PionVeeD3E) && (ProtonVeeD3 > 3*ProtonVeeD3E) ){
	     masshisto.Fill(ximass,Xi_massIDIP3D2);     
	     if (NewPrimaryProb > 0.01 ){
	        masshisto.Fill(ximass,Xi_massIDprim5);     
	     }
	  }
	  if ((PionVeeD0 > 3*PionVeeD0E) && (ProtonVeeD0 > 3*ProtonVeeD0E) ){
	     masshisto.Fill(ximass,Xi_massIDIP2D2);     
	     if (NewPrimaryProb > 0.01 ){
	        masshisto.Fill(ximass,Xi_massIDprim6);     
	     }
	  }
	  if ((PionVeeD3 > 3*PionVeeD3E)
	   && (ProtonVeeD3 > 3*ProtonVeeD3E)
	   && (PionD3 > 3*PionD3E) ){
	     masshisto.Fill(ximass,Xi_massIDIP3D);     
	     if (NewPrimaryProb > 0.01 ){
	        masshisto.Fill(ximass,Xi_massIDprim7);
		if (probvxi > 0.01 && probvee > 0.01 ){
		   masshisto.Fill(ximass,Xi_massIDprim16);	
		}     
	     }
	  }
	  
          if (r_cas > r_trk) {
            massrejhisto.Fill(ximass,Xi_massID);
	    masshisto.Fill(ximass,Xi_massIDtype2);
	    if (probvxi > 0.01){
	      massrejhisto.Fill(ximass,Xi_massIDprob1);
              masshisto.Fill(ximass,Xi_massIDprob2);
	      if (probvee > 0.01){
		 if (NewPrimaryProb > 0.01 ){
		     massrejhisto.Fill(ximass,Xi_massIDprim15);
	             if ((PionVeeD3 > 3*PionVeeD3E)
	              && (ProtonVeeD3 > 3*ProtonVeeD3E)
	              && (PionD3 > 3*PionD3E) ){
			 massrejhisto.Fill(ximass,Xi_massIDprim16);    
		     }
		     if (Xiprmsep.mag() >  3*XiprmsepE) {
		         massrejhisto.Fill(ximass,Xi_massIDprim17);   
		     }
		     if (Xiprmsep.mag() >  4*XiprmsepE) {
		         massrejhisto.Fill(ximass,Xi_massIDprim18);  
		     }
		     
		 }
		 
	      }	   
            }
	  }else{
            masshisto.Fill(ximass,Xi_massIDtype1);
	    if (probvxi > 0.01)
              masshisto.Fill(ximass,Xi_massIDprob1);
	  }
          chisqhisto.Fill(chisq,Xi2ID);
	  for (int j=0;j<11;j++){
	    std::stringstream ssj;
	    if ( j == 0 ){
	        ssj << 100;
	    }else{
	        ssj << j;
	    }
	    Xi_massID="#Xi_" + ss.str()+ "#chi^2_"+ssj.str();
	    if (r_cas > r_trk) {
	      if (chisq < j && j>0)
	        massrejhisto.Fill(ximass,Xi_massID);
	      if (chisq < 100 && j==0)
	        massrejhisto.Fill(ximass,Xi_massID);
	    }else{
	      if (chisq < j && j>0)
	        masshisto.Fill(ximass,Xi_massID);
	      if (chisq < 100 && j==0)
	        masshisto.Fill(ximass,Xi_massID);
	    }
	  }
       }
     }
   }//loop cascades

}//end analyzer


// ------------ method called once each job just before starting event loop  ------------
void 
CascadeAnalyzer::beginJob()
{
  edm::Service<TFileService> fs;

  //Directories
  const std::string massdirname("masshistos");
  TFileDirectory veedir     = TFileDirectory( fs->mkdir("veemasshistos") );
  TFileDirectory massdir    = TFileDirectory( fs->mkdir("masshistos") );
  TFileDirectory masscasdir = TFileDirectory( fs->mkdir("mass1cas") );
  TFileDirectory massrejdir = TFileDirectory( fs->mkdir("type2histos") );
  TFileDirectory xi2dir     = TFileDirectory( fs->mkdir("Xi2histos") );
  TFileDirectory cldir     = TFileDirectory( fs->mkdir("clhistos") );

  veemasshisto.Set("#Lambda inv. mass",veedir,"","p#pi mass");
  veemasshisto.Set("#Lambda mass refit",veedir,"","p #pi refit");
  veemasshisto.Set("#Delta#Lambda mass",veedir,"","#Delta p#pi + 1.1");
  veemasshisto.Set("Ks mass",veedir,"","M #pi #pi  -0.497648 + 1.1");
  veemasshisto.Set("#Delta Ks mass",veedir,"","M #pi#pi -0.497648 + 1.1");
  veemasshisto.Set("Ks mass veto",veedir,"","M #pi#pi -0.497648 + 1.1");

   for (int i=0;i<16;i++){
   std::stringstream ss;
   ss << i;
   std::string mytitle0("#Lambda #pi for R > " + ss.str() );
   std::string Xi_massID("#Xi_" + ss.str());
   std::string Xi_CLID("CL#Xi_" + ss.str());
   masshisto.Set(Xi_massID,massdir,"pass0_",mytitle0); //all
   //----Normalizing histos ---- single cuts ------		   
   std::string mytitle(mytitle0 + ", type 1");		   
   masshisto.Set(Xi_massID,massdir,"type1_",mytitle);  //    type 1 
   mytitle = mytitle0 + ", type 2";		   
   masshisto.Set(Xi_massID,massdir,"type2_",mytitle);  //    type 2 
   mytitle = mytitle0 + ",#Xi CL>0.01";
   masshisto.Set(Xi_massID,massdir,"prob0_",mytitle); //    XiCL>0.01
   mytitle = mytitle0 + ",#Lambda CL>0.01";
   masshisto.Set(Xi_massID,massdir,"probv_",mytitle); //    /\0CL>0.01
   mytitle = mytitle0 + ",#pi trkhits>5";
   masshisto.Set(Xi_massID,massdir,"hits6_",mytitle); //   trkhits > 5
   mytitle = mytitle0 + ",#pi trkhits>4";
   masshisto.Set(Xi_massID,massdir,"hits5_",mytitle); //   trkhits > 4
   mytitle = mytitle0 + ",#pi trkhits>3";
   masshisto.Set(Xi_massID,massdir,"hits4_",mytitle); //   trkhits > 3
   mytitle = mytitle0 + ",#pi ip3D>3#sigma";
   masshisto.Set(Xi_massID,massdir,"IP3D(#pi)_",mytitle);//  ip3D pi > 3 sigma
   mytitle = mytitle0 + ",#pi ip2D>3#sigma";
   masshisto.Set(Xi_massID,massdir,"IP2D(#pi)_",mytitle);//  ip2D pi > 3 sigma
   mytitle = mytitle0 + ",#Lambda ip3D>3#sigma";
   masshisto.Set(Xi_massID,massdir,"IP3D(#Lambda)_",mytitle);//  ip3D /\0 > 3 sigma
   mytitle = mytitle0 + ",#Lambda ip2D>3#sigma";
   masshisto.Set(Xi_massID,massdir,"IP2D(#Lambda)_",mytitle);//  ip2D /\0 > 3 sigma
   mytitle = mytitle0 + ",#Xi trks ip3D>3#sigma";
   masshisto.Set(Xi_massID,massdir,"IP3D_",mytitle);   //  ip3D Xi trks > 3 sigma
   mytitle = mytitle0 + ",prmCL>0.01";
   masshisto.Set(Xi_massID,massdir,"prm0_",mytitle);  //    prmCL>0.01
   mytitle = mytitle0 + ",prmrefitCL>0.01";
   masshisto.Set(Xi_massID,massdir,"prmre_",mytitle); //    prmrefitCL>0.01
   //---- histos ---- combined cuts ------		   
   mytitle = mytitle0 + ",#Xi CL>0.01, type 1";
   masshisto.Set(Xi_massID,massdir,"prob1_",mytitle); //    XiCL>0.01 type 1
   mytitle = mytitle0 + ",#Xi CL>0.01, type 2";
   masshisto.Set(Xi_massID,massdir,"prob2_",mytitle); //    XiCL>0.01 type 2
   mytitle = mytitle0 + ",prmCL<prmrefitCL,prmrefitCL>0.01";
   masshisto.Set(Xi_massID,massdir,"prm1_",mytitle);  //    prmCL<prmrefitCL,prmrefitCL>0.01
   mytitle = mytitle0 + ",prmrefitCL>0.01,#pi ip3D>3#sigma";
   masshisto.Set(Xi_massID,massdir,"prmreIP3D(#pi)_",mytitle);    //    prmrefitCL>0.01,ip3D pi > 3 sigma
   mytitle = mytitle0 + ",prmrefitCL>0.01,#pi ip2D>3#sigma";
   masshisto.Set(Xi_massID,massdir,"prmreIP2D(#pi)_",mytitle);    //    prmrefitCL>0.01,ip2D pi > 3 sigma
   mytitle = mytitle0 + ",prmrefitCL>0.01,#Lambda ip3D>3#sigma";
   masshisto.Set(Xi_massID,massdir,"prmreIP3D(#Lambda)_",mytitle);//    prmrefitCL>0.01,ip3D /\0 > 3 sigma
   mytitle = mytitle0 + ",prmrefitCL>0.01,#Lambda ip2D>3#sigma";
   masshisto.Set(Xi_massID,massdir,"prmreIP2D(#Lambda)_",mytitle);//    prmrefitCL>0.01,ip2D /\0 > 3 sigma
   mytitle = mytitle0 + ",prmrefitCL>0.01,#Xi CL>0.01";
   masshisto.Set(Xi_massID,massdir,"prmrevxi_",mytitle);          //    prmrefitCL>0.01,XiCL>0.01
   mytitle = mytitle0 + ",prmrefitCL>0.01,#Xi CL>0.01,#Lambda CL>0.01";
   masshisto.Set(Xi_massID,massdir,"prmrevxivee_",mytitle);       //    prmrefitCL>0.01,XiCL>0.01,/\0CL>0.01
   mytitle = mytitle0 + ",prmrefitCL>0.01,#Xi CL>0.01,#Lambda CL>0.01,#Xi trks ip3D>3#sigma";
   masshisto.Set(Xi_massID,massdir,"prmrevxiveeIP3D_",mytitle);   //    prmrefitCL>0.01,XiCL>0.01,/\0CL>0.01,ip3D Xi trks > 3 sigma
   mytitle = mytitle0 + ",prmrefitCL>0.01,#Xi CL>0.01,#Lambda CL>0.01,#Xi trks ip3DTool>3#sigma";
   masshisto.Set(Xi_massID,massdir,"prmrevxiveeIPTool3_",mytitle);   //    prmrefitCL>0.01,XiCL>0.01,/\0CL>0.01,ip3D Xi trks > 3 sigma
   mytitle = mytitle0 + ",prmrefitCL>0.01,#Xi CL>0.01,#Lambda CL>0.01,#Xi L>3#sigma";
   masshisto.Set(Xi_massID,massdir,"prmrevxiveels3_",mytitle);   //    prmrefitCL>0.01,XiCL>0.01,/\0CL>0.01,ip3D Xi L > 3 sigma
   clhisto.Set(Xi_CLID,cldir,"prmrevxiveels3_",mytitle);
   mytitle = mytitle0 + ",prmrefitCL>0.01,#Xi CL>0.01,#Lambda CL>0.01,#Xi L>3#sigma trkhits>3";
   masshisto.Set(Xi_massID,massdir,"prmrevxiveels3hit3_",mytitle);   //    prmrefitCL>0.01,XiCL>0.01,/\0CL>0.01,ip3D Xi L > 3 sigma trkhits>3
   clhisto.Set(Xi_CLID,cldir,"prmrevxiveels3hit3_",mytitle);
   mytitle = mytitle0 + ",prmrefitCL>0.01,#Xi CL>0.01,#Lambda CL>0.01,#Xi L>3#sigma trkhits>4";
   masshisto.Set(Xi_massID,massdir,"prmrevxiveels3hit4_",mytitle);   //    prmrefitCL>0.01,XiCL>0.01,/\0CL>0.01,ip3D Xi L > 3 sigma trkhits>4
   clhisto.Set(Xi_CLID,cldir,"prmrevxiveels3hit4_",mytitle);
   mytitle = mytitle0 + ",prmrefitCL>0.01,#Xi CL>0.01,#Lambda CL>0.01,#Xi L>3#sigma trkhits>5";
   masshisto.Set(Xi_massID,massdir,"prmrevxiveels3hit5_",mytitle);   //    prmrefitCL>0.01,XiCL>0.01,/\0CL>0.01,ip3D Xi L > 3 sigma trkhits>5
   clhisto.Set(Xi_CLID,cldir,"prmrevxiveels3hit5_",mytitle);
   mytitle = mytitle0 + ",prmrefitCL>0.01,#Xi CL>0.01,#Lambda CL>0.01,#Xi L>3#sigma Rvtx<trkhit";
   masshisto.Set(Xi_massID,massdir,"prmrevxiveels3Rvtx_",mytitle);   //    prmrefitCL>0.01,XiCL>0.01,/\0CL>0.01,ip3D Xi L > 3 sigma trkhits>5		   
   clhisto.Set(Xi_CLID,cldir,"prmrevxiveels3Rvtx_",mytitle);
   
//IsKs
   mytitle = mytitle0 + ",prmrefitCL>0.01,#Xi CL>0.01,#Lambda CL>0.01,IsK_S,#Xi trks ip3DTool>3#sigma";
   masshisto.Set(Xi_massID,massdir,"prmrevxiveeIsKsIPTool3_",mytitle);   //    prmrefitCL>0.01,XiCL>0.01,/\0CL>0.01,Is Ks,ip3D Xi trks > 3 sigma
   mytitle = mytitle0 + ",prmrefitCL>0.01,#Xi CL>0.01,#Lambda CL>0.01,IsK_S,#Xi L>3#sigma";
   masshisto.Set(Xi_massID,massdir,"prmrevxiveeIsKsls3_",mytitle);   //    prmrefitCL>0.01,XiCL>0.01,/\0CL>0.01,Is Ks,Xi L > 3 sigma
   mytitle = mytitle0 + ",prmrefitCL>0.01,#Xi CL>0.01,#Lambda CL>0.01,IsK_S,#Xi L>3#sigma trkhits>3";
   masshisto.Set(Xi_massID,massdir,"prmrevxiveeIsKsls3hit3_",mytitle);   //    prmrefitCL>0.01,XiCL>0.01,/\0CL>0.01,Is Ks,Xi L > 3 sigma trkhits>3
   mytitle = mytitle0 + ",prmrefitCL>0.01,#Xi CL>0.01,#Lambda CL>0.01,IsK_S,#Xi L>3#sigma trkhits>4";
   masshisto.Set(Xi_massID,massdir,"prmrevxiveeIsKsls3hit4_",mytitle);   //    prmrefitCL>0.01,XiCL>0.01,/\0CL>0.01,Is Ks,Xi L > 3 sigma trkhits>4
   mytitle = mytitle0 + ",prmrefitCL>0.01,#Xi CL>0.01,#Lambda CL>0.01,IsK_S,#Xi L>3#sigma trkhits>5";
   masshisto.Set(Xi_massID,massdir,"prmrevxiveeIsKsls3hit5_",mytitle);   //    prmrefitCL>0.01,XiCL>0.01,/\0CL>0.01,Is Ks,Xi L > 3 sigma trkhits>5
   mytitle = mytitle0 + ",prmrefitCL>0.01,#Xi CL>0.01,#Lambda CL>0.01,IsK_S,#Xi L>3#sigma Rvtx<trkhit";
   masshisto.Set(Xi_massID,massdir,"prmrevxiveeIsKsls3Rvtx_",mytitle);   //    prmrefitCL>0.01,XiCL>0.01,/\0CL>0.01,Is Ks,Xi L > 3 sigma r_cas < rtrk		   
   mytitle = mytitle0 + ",prmrefitCL>0.01,#Xi CL>0.01,#Lambda CL>0.01,IsK_S,#Xi L>4#sigma";
   masshisto.Set(Xi_massID,massdir,"prmrevxiveeIsKsls4_",mytitle);   //    prmrefitCL>0.01,XiCL>0.01,/\0CL>0.01,Is Ks,Xi L > 4 sigma
//IsKs  end
//NoKs
   mytitle = mytitle0 + ",prmrefitCL>0.01,#Xi CL>0.01,#Lambda CL>0.01,NoK_S,#Xi trks ip3DTool>3#sigma";
   masshisto.Set(Xi_massID,massdir,"prmrevxiveeNoKsIPTool3_",mytitle);   //    prmrefitCL>0.01,XiCL>0.01,/\0CL>0.01,No Ks,ip3D Xi trks > 3 sigma
   mytitle = mytitle0 + ",prmrefitCL>0.01,#Xi CL>0.01,#Lambda CL>0.01,NoK_S,#Xi L>3#sigma";
   masshisto.Set(Xi_massID,massdir,"prmrevxiveeNoKsls3_",mytitle);   //    prmrefitCL>0.01,XiCL>0.01,/\0CL>0.01,No Ks,Xi L > 3 sigma
   mytitle = mytitle0 + ",prmrefitCL>0.01,#Xi CL>0.01,#Lambda CL>0.01,NoK_S,#Xi L>3#sigma trkhits>3";
   masshisto.Set(Xi_massID,massdir,"prmrevxiveeNoKsls3hit3_",mytitle);   //    prmrefitCL>0.01,XiCL>0.01,/\0CL>0.01,No Ks,Xi L > 3 sigma trkhits>3
   mytitle = mytitle0 + ",prmrefitCL>0.01,#Xi CL>0.01,#Lambda CL>0.01,NoK_S,#Xi L>3#sigma trkhits>4";
   masshisto.Set(Xi_massID,massdir,"prmrevxiveeNoKsls3hit4_",mytitle);   //    prmrefitCL>0.01,XiCL>0.01,/\0CL>0.01,No Ks,Xi L > 3 sigma trkhits>4
   mytitle = mytitle0 + ",prmrefitCL>0.01,#Xi CL>0.01,#Lambda CL>0.01,NoK_S,#Xi L>3#sigma trkhits>5";
   masshisto.Set(Xi_massID,massdir,"prmrevxiveeNoKsls3hit5_",mytitle);   //    prmrefitCL>0.01,XiCL>0.01,/\0CL>0.01,No Ks,Xi L > 3 sigma trkhits>5
   mytitle = mytitle0 + ",prmrefitCL>0.01,#Xi CL>0.01,#Lambda CL>0.01,NoK_S,#Xi L>3#sigma Rvtx<trkhit";
   masshisto.Set(Xi_massID,massdir,"prmrevxiveeNoKsls3Rvtx_",mytitle);   //    prmrefitCL>0.01,XiCL>0.01,/\0CL>0.01,No Ks,Xi L > 3 sigma r_cas < rtrk		   
   mytitle = mytitle0 + ",prmrefitCL>0.01,#Xi CL>0.01,#Lambda CL>0.01,NoK_S,#Xi L>4#sigma";
   masshisto.Set(Xi_massID,massdir,"prmrevxiveeNoKsls4_",mytitle);   //    prmrefitCL>0.01,XiCL>0.01,/\0CL>0.01,No Ks,Xi L > 4 sigma
//NoKs  end
//Xi CL>5%
   mytitle = mytitle0 + ",prmrefitCL>0.01,#Xi CL>0.05,#Lambda CL>0.01,#Xi L>3#sigma";
   masshisto.Set(Xi_massID,massdir,"prmrevxi05veels3_",mytitle);   //    prmrefitCL>0.01,XiCL>0.05,/\0CL>0.01,Xi L > 3 sigma
   mytitle = mytitle0 + ",prmrefitCL>0.01,#Xi CL>0.05,#Lambda CL>0.01,#Xi L>3#sigma trkhits>3";
   masshisto.Set(Xi_massID,massdir,"prmrevxi05veels3hit3_",mytitle);   //    prmrefitCL>0.01,XiCL>0.05,/\0CL>0.01,Xi L > 3 sigma trkhits>3
   mytitle = mytitle0 + ",prmrefitCL>0.01,#Xi CL>0.05,#Lambda CL>0.01,#Xi L>3#sigma trkhits>4";
   masshisto.Set(Xi_massID,massdir,"prmrevxi05veels3hit4_",mytitle);   //    prmrefitCL>0.01,XiCL>0.05,/\0CL>0.01,Xi L > 3 sigma trkhits>4
   mytitle = mytitle0 + ",prmrefitCL>0.01,#Xi CL>0.05,#Lambda CL>0.01,#Xi L>3#sigma trkhits>5";
   masshisto.Set(Xi_massID,massdir,"prmrevxi05veels3hit5_",mytitle);   //    prmrefitCL>0.01,XiCL>0.05,/\0CL>0.01,Xi L > 3 sigma trkhits>5
   mytitle = mytitle0 + ",prmrefitCL>0.01,#Xi CL>0.05,#Lambda CL>0.01,#Xi L>3#sigma Rvtx<trkhit";
   masshisto.Set(Xi_massID,massdir,"prmrevxi05veels3Rvtx_",mytitle);   //    prmrefitCL>0.01,XiCL>0.05,/\0CL>0.01,Xi L > 3 sigma r_cas < rtrk              
   mytitle = mytitle0 + ",prmrefitCL>0.01,#Xi CL>0.05,#Lambda CL>0.01,#Xi L>3#sigma ip3DTool>3#sigma";
   masshisto.Set(Xi_massID,massdir,"prmrevxi05veels3IPTool3_",mytitle);   // prmrefitCL>0.01,XiCL>0.05,/\0CL>0.01,Xi L > 3 sigma IPTool3 > 3
   mytitle = mytitle0 + ",prmrefitCL>0.01,#Xi CL>0.05,#Lambda CL>0.01,#Xi L>4#sigma";
   masshisto.Set(Xi_massID,massdir,"prmrevxi05veels4_",mytitle);   //    prmrefitCL>0.01,XiCL>0.05,/\0CL>0.01,Xi L > 4 sigma
//End Xi CL>%5
//Xi Pixelless trk daughter
   mytitle = mytitle0 + ",prmrefitCL>0.01,#Xi CL>0.01,#Lambda CL>0.01,#Xi L>3#sigma, pixelless";
   masshisto.Set(Xi_massID,massdir,"prmrevxiveeNoPixells3_",mytitle);     // prmrefitCL>0.01,XiCL>0.01,/\0CL>0.01,Xi L > 3 sigma
   mytitle = mytitle0 + ",prmrefitCL>0.01,#Xi CL>0.01,#Lambda CL>0.01,#Xi L>3#sigma trkhits>3 pixelless";
   masshisto.Set(Xi_massID,massdir,"prmrevxiveeNoPixells3hit3_",mytitle); // prmrefitCL>0.01,XiCL>0.01,/\0CL>0.01,Xi L > 3 sigma trkhits>3
   mytitle = mytitle0 + ",prmrefitCL>0.01,#Xi CL>0.01,#Lambda CL>0.01,#Xi L>3#sigma trkhits>4 pixelless";
   masshisto.Set(Xi_massID,massdir,"prmrevxiveeNoPixells3hit4_",mytitle); // prmrefitCL>0.01,XiCL>0.01,/\0CL>0.01,Xi L > 3 sigma trkhits>4
   mytitle = mytitle0 + ",prmrefitCL>0.01,#Xi CL>0.01,#Lambda CL>0.01,#Xi L>3#sigma trkhits>5 pixelless";
   masshisto.Set(Xi_massID,massdir,"prmrevxiveeNoPixells3hit5_",mytitle); // prmrefitCL>0.01,XiCL>0.01,/\0CL>0.01,Xi L > 3 sigma trkhits>5
   mytitle = mytitle0 + ",prmrefitCL>0.01,#Xi CL>0.01,#Lambda CL>0.01,#Xi L>3#sigma Rvtx<trkhit pixelless";
   masshisto.Set(Xi_massID,massdir,"prmrevxiveeNoPixells3Rvtx_",mytitle); // prmrefitCL>0.01,XiCL>0.01,/\0CL>0.01,Xi L > 3 sigma r_cas < rtrk
   mytitle = mytitle0 + ",prmrefitCL>0.01,#Xi CL>0.01,#Lambda CL>0.01,#Xi L>3#sigma ip3DTool>3#sigma pixelless";
   masshisto.Set(Xi_massID,massdir,"prmrevxiveeNoPixells3IPTool3_",mytitle); // prmrefitCL>0.01,XiCL>0.01,/\0CL>0.01,Xi L > 3 sigma IPTool3 > 3
   mytitle = mytitle0 + ",prmrefitCL>0.01,#Xi CL>0.01,#Lambda CL>0.01,#Xi L>4#sigma pixelless";
   masshisto.Set(Xi_massID,massdir,"prmrevxiveeNoPixells4_",mytitle);     // prmrefitCL>0.01,XiCL>0.01,/\0CL>0.01,Xi L > 4 sigma
//End Xi Pixelless trk daughter

   mytitle = mytitle0 + ",prmrefitCL>0.01,#Xi CL>0.01,#Lambda CL>0.01,#Xi L>3#sigma,#Xi trks ip3D>3#sigma";
   masshisto.Set(Xi_massID,massdir,"prmrevxiveels3IP3D_",mytitle);   //    prmrefitCL>0.01,XiCL>0.01,/\0CL>0.01,ip3D Xi trks > 3 sigma
   mytitle = mytitle0 + ",prmrefitCL>0.01,#Xi CL>0.01,#Lambda CL>0.01,#Xi L>3#sigma,#Xi trks ip3DTool>3#sigma";
   masshisto.Set(Xi_massID,massdir,"prmrevxiveels3IPTool3_",mytitle);   //    prmrefitCL>0.01,XiCL>0.01,/\0CL>0.01,ip3D Xi trks > 3 sigma
   mytitle = mytitle0 + ",prmrefitCL>0.01,#Xi CL>0.01,#Lambda CL>0.01,#Xi L>4#sigma";
   masshisto.Set(Xi_massID,massdir,"prmrevxiveels4_",mytitle);   //    prmrefitCL>0.01,XiCL>0.01,/\0CL>0.01,ip3D Xi L > 3 sigma
   clhisto.Set(Xi_CLID,cldir,"prmrevxiveels4_",mytitle);
   mytitle = mytitle0 + ",prmrefitCL>0.01,#Xi trks ip3D>3#sigma";
   masshisto.Set(Xi_massID,massdir,"prmreIP3D_",mytitle);         //    prmrefitCL>0.01,ip3D Xi trks > 3 sigma
   mytitle = mytitle0 + ",prmrefitCL>0.01,#pi trkhits>3";
   masshisto.Set(Xi_massID,massdir,"prmrehits4_",mytitle);        //    prmrefitCL>0.01,trkhits > 3
   mytitle = mytitle0 + ",prmrefitCL>0.01,#pi trkhits>4";
   masshisto.Set(Xi_massID,massdir,"prmrehits5_",mytitle);        //    prmrefitCL>0.01,trkhits > 4
   mytitle = mytitle0 + ",prmrefitCL>0.01,#pi trkhits>5";
   masshisto.Set(Xi_massID,massdir,"prmrehits6_",mytitle);        //    prmrefitCL>0.01,trkhits > 5
   mytitle = mytitle0 + ",prmrefitCL>0.01,#pi trkhits>3,#Xi trks ip3D>3#sigma";
   masshisto.Set(Xi_massID,massdir,"prmrehits4IP3D_",mytitle);    //    prmrefitCL>0.01,trkhits > 3,ip3D Xi trks > 3 sigma
   mytitle = mytitle0 + ",prmrefitCL>0.01,#pi trkhits>4,#Xi trks ip3D>3#sigma";
   masshisto.Set(Xi_massID,massdir,"prmrehits5IP3D_",mytitle);    //    prmrefitCL>0.01,trkhits > 4,ip3D Xi trks > 3 sigma
   mytitle = mytitle0 + ",prmrefitCL>0.01,#pi trkhits>5,#Xi trks ip3D>3#sigma";
   masshisto.Set(Xi_massID,massdir,"prmrehits6IP3D_",mytitle);    //    prmrefitCL>0.01,trkhits > 5,ip3D Xi trks > 3 sigma
   //type 2
   mytitle = mytitle0;
   massrejhisto.Set(Xi_massID,massrejdir,"",mytitle);
   mytitle = mytitle0 + ",#Xi CL>0.01";
   massrejhisto.Set(Xi_massID,massrejdir,"prob1_",mytitle);
   mytitle = mytitle0 + ",prmrefitCL>0.01,#Xi CL>0.01,#Lambda CL>0.01";
   massrejhisto.Set(Xi_massID,massrejdir,"prmrevxivee_",mytitle);       //    prmrefitCL>0.01,XiCL>0.01,/\0CL>0.01
   mytitle = mytitle0 + ",prmrefitCL>0.01,#Xi CL>0.01,#Lambda CL>0.01,#Xi trks ip3D>3#sigma";
   massrejhisto.Set(Xi_massID,massrejdir,"prmrevxiveeIP3D_",mytitle);   //    prmrefitCL>0.01,XiCL>0.01,/\0CL>0.01,ip3D Xi trks > 3 sigma
   mytitle = mytitle0 + ",prmrefitCL>0.01,#Xi CL>0.01,#Lambda CL>0.01,#Xi L>3#sigma";
   massrejhisto.Set(Xi_massID,massrejdir,"prmrevxiveels3_",mytitle);   //    prmrefitCL>0.01,XiCL>0.01,/\0CL>0.01,ip3D Xi L > 3 sigma
   mytitle = mytitle0 + ",prmrefitCL>0.01,#Xi CL>0.01,#Lambda CL>0.01,#Xi L>4#sigma";
   massrejhisto.Set(Xi_massID,massrejdir,"prmrevxiveels4_",mytitle);   //    prmrefitCL>0.01,XiCL>0.01,/\0CL>0.01,ip3D Xi L > 3 sigma
   //1 cascade per event
   std::string Xi_massID1("#Xi1_" + ss.str());
   mytitle = mytitle0;
   masshisto.Set(Xi_massID1,masscasdir,"all1cas_",mytitle); //all
   mytitle = mytitle0 + ",prmrefitCL>0.01";
   masshisto.Set(Xi_massID1,masscasdir,"prmre_",mytitle); //    prmrefitCL>0.01
   mytitle = mytitle0 + ",prmrefitCL>0.01,#Xi CL>0.01";
   masshisto.Set(Xi_massID1,masscasdir,"prmrevxi_",mytitle);          //    prmrefitCL>0.01,XiCL>0.01
   mytitle = mytitle0 + ",prmrefitCL>0.01,#Xi CL>0.01,#Lambda CL>0.01";
   masshisto.Set(Xi_massID1,masscasdir,"prmrevxivee_",mytitle);       //    prmrefitCL>0.01,XiCL>0.01,/\0CL>0.01
   mytitle = mytitle0 + ",prmrefitCL>0.01,#Xi CL>0.01,#Lambda CL>0.01,#Xi L>3#sigma";
   masshisto.Set(Xi_massID1,masscasdir,"prmrevxiveels3_",mytitle);   //    prmrefitCL>0.01,XiCL>0.01,/\0CL>0.01,ip3D Xi L > 3 sigma
   mytitle = mytitle0 + ",prmrefitCL>0.01,#Xi CL>0.01,#Lambda CL>0.01,#Xi L>3#sigma trkhits>3";
   masshisto.Set(Xi_massID1,masscasdir,"prmrevxiveels3hit3_",mytitle);   //    prmrefitCL>0.01,XiCL>0.01,/\0CL>0.01,ip3D Xi L > 3 sigma trkhits>3
   mytitle = mytitle0 + ",prmrefitCL>0.01,#Xi CL>0.01,#Lambda CL>0.01,#Xi L>3#sigma trkhits>4";
   masshisto.Set(Xi_massID1,masscasdir,"prmrevxiveels3hit4_",mytitle);   //    prmrefitCL>0.01,XiCL>0.01,/\0CL>0.01,ip3D Xi L > 3 sigma trkhits>4
   mytitle = mytitle0 + ",prmrefitCL>0.01,#Xi CL>0.01,#Lambda CL>0.01,#Xi L>3#sigma trkhits>5";
   masshisto.Set(Xi_massID1,masscasdir,"prmrevxiveels3hit5_",mytitle);   //    prmrefitCL>0.01,XiCL>0.01,/\0CL>0.01,ip3D Xi L > 3 sigma trkhits>5
   mytitle = mytitle0 + ",prmrefitCL>0.01,#Xi CL>0.01,#Lambda CL>0.01,#Xi L>3#sigma Rvtx<trkhit";
   masshisto.Set(Xi_massID1,masscasdir,"prmrevxiveels3Rvtx_",mytitle);   //    prmrefitCL>0.01,XiCL>0.01,/\0CL>0.01,ip3D Xi L > 3 sigma trkhits>5		   
   mytitle = mytitle0 + ",prmrefitCL>0.01,#Xi CL>0.01,#Lambda CL>0.01,#Xi L>4#sigma";
   masshisto.Set(Xi_massID1,masscasdir,"prmrevxiveels4_",mytitle);   //    prmrefitCL>0.01,XiCL>0.01,/\0CL>0.01,ip3D Xi L > 3 sigma
   
   mytitle ="#chi^2 for R >" + ss.str();
   std::string Xi2ID("chisq_" + ss.str());
   chisqhisto.Set(Xi2ID,xi2dir,"#chi^2_",mytitle);
   for (int j=0;j<11;j++){
      std::stringstream ssj;
      if ( j == 0 ){
        ssj << 100;
      }else{
        ssj << j;
      }
      mytitle ="#Lambda #pi^- for R > " + ss.str() + " #chi^2 <" +ssj.str();
      Xi_massID="#Xi_" + ss.str()+ "#chi^2_"+ssj.str();
      masshisto.Set(Xi_massID,massdir,"",mytitle);
      massrejhisto.Set(Xi_massID,massrejdir,"",mytitle);
   }//loop j
  }//loop i
  clhisto.SetPrefix("");    //prefix inverted this patch it
  // Set to save jet histogram errors
//  masshisto.Sumw2();
//  chisqhisto.Sumw2();
}//begin job

// ------------ method called once each job just after ending the event loop  ------------
void 
CascadeAnalyzer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(CascadeAnalyzer);
