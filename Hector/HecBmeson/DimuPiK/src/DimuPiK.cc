// -*- C++ -*-
//
// Package:    DimuPiK
// Class:      DimuPiK
// 
/**\class DimuPiK DimuPiK.cc HecMeson/DimuPiK/src/DimuPiK.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Hector Mendez
//         Created:  Tue Oct  2 09:00:11 CDT 2012
// $Id: DimuPiK.cc,v 1.1 2013/03/18 19:52:23 mendez Exp $
//
//

// system include files
#include <memory>

// user include files
#include "HecBmeson/DimuPiK/interface/DimuPiK.h"                    //--Hec
int Lb[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
//#include "FWCore/MessageLogger/interface/MessageLogger.h"             //--Hec
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
//#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"             //--Hec

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonTrackLinks.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/CaloMuon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"    //cascade colletcion
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"    //cascade daughter track

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"       //track   collection
//#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
//#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"    //cascade colletcion
//#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"    //cascade daughter track

#include "Analyzers/CascadeProducer/interface/PrimaryInfo.h"             //--Eduardo Cascade & Primary
#include "Analyzers/CascadeProducer/interface/VertexRefit.h"
//--#include "Analyzers/CascadeProducer/interface/KinematicFitDriver.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "Analyzers/CascadeProducer/interface/ClosestApproachOnHelixLine.h"

//--Bfield
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
//#include "DataFormats/VertexReco/interface/VertexFwd.h"
//#include "DataFormats/VertexReco/interface/Vertex.h"


//--kinemactic fitter
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
//--Trigger (from keith)
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"

#include "Analyzers/CascadeProducer/interface/masses.h"

#include "TH1.h"

#include <vector>
#include <utility>
//
// constants, enums and typedefs
  const double kshort_mass_c   = 0.497648;
//
// static data member definitions
//
// constructors and destructor
//
DimuPiK::DimuPiK(const edm::ParameterSet& iConfig)
:
 trackTags_       (iConfig.getParameter<edm::InputTag>("tracks"    )),
 theMuonsLabel_   (iConfig.getParameter<edm::InputTag>("MuonsLabel")),
 VeeAlgo_         (iConfig.getUntrackedParameter<std::string>("VeeAlgo","generalV0Candidates")),
 MyPrint_         (iConfig.getUntrackedParameter<bool>("MyPrint","False")),
 MyPrintMC_        (iConfig.getUntrackedParameter<bool>("MyPrintMC", "false")),
 doMC_             (iConfig.getUntrackedParameter<bool>("doMC", "false") ),
 ChargedKa_       (iConfig.getUntrackedParameter<bool>("ChargedKa","False")),
 hlTriggerResults_(iConfig.getUntrackedParameter<std::string>("HLTriggerResults",std::string("TriggerResults::HLT")) ) 
{
  NeutralKa_ = false;
  if( !ChargedKa_ ){
      NeutralKa_ = true;
      std::cout<<"----doing Neutral Kaon:  B+-->u+u-Pi+K0"<<std::endl;
  } else
      std::cout<<"----doing Charged Kaon:  B0-->u+u-Pi+K-"<<std::endl;
    
  edm::Service<TFileService> fs;
  //--Diagnostic and testing Histograms
  histo_trk    = fs->make<TH1D>("Ntrk"   , "Ntrk"    , 200  , 0 , 200 );
  histo_nmu    = fs->make<TH1D>("NMuons" , "NMuons"  , 200  , 0 , 200 );
  histo_muId   = fs->make<TH1D>("muId"   , "Muon Id" , 40   , 0 , 40  );
  histo_primVtx= fs->make<TH1D>("primVtx", "primVtx" , 100  , 0 , 100 );
  histo_kshort = fs->make<TH1D>("kshort" , "Kshort"  , 100  , 0 , 100 );
    
  int Mnch=6000;
  double Mxi=0., Mxf=150.;
  histo_dimu[0] = fs->make<TH1D>("dimu0", "DiMuom All"     , Mnch , Mxi , Mxf );    //--25 MeV/channel
  histo_dimu[1] = fs->make<TH1D>("dimu1", "DiMuon Closest" , Mnch , Mxi , Mxf );
  histo_dimu[2] = fs->make<TH1D>("dimu2", "DiMuon Vtx"     , Mnch , Mxi , Mxf );
  
  histo_dimuPi[0] = fs->make<TH1D>("dimuPi0", "DiMuon Pi All"  , Mnch , Mxi , Mxf );
  histo_dimuPi[1] = fs->make<TH1D>("dimuPi1", "DiMuon Pi Clo1" , Mnch , Mxi , Mxf );
  histo_dimuPi[2] = fs->make<TH1D>("dimuPi2", "DiMuon Pi Clo2" , Mnch , Mxi , Mxf );
  histo_dimuPi[3] = fs->make<TH1D>("dimuPi3", "DiMuon Pi Clo3" , Mnch , Mxi , Mxf );
  histo_dimuPi[4] = fs->make<TH1D>("dimuPi4", "DiMuon Pi Clo4" , Mnch , Mxi , Mxf );
  histo_dimuPi[5] = fs->make<TH1D>("dimuPi5", "DiMuon Pi Clo5" , Mnch , Mxi , Mxf );
  
  histo_dimuPiK[0] = fs->make<TH1D>("dimuPiK0", "DiMuon Pi K0" , Mnch , Mxi , Mxf );
  histo_dimuPiK[1] = fs->make<TH1D>("dimuPiK1", "DiMuon Pi K1" , Mnch , Mxi , Mxf );
  
  histo_k0s[0] = fs->make<TH1D>("k0s0", "Kshort All 0" , 150, 0.35, 0.65); //--2.0 MeV/channel
  histo_k0s[1] = fs->make<TH1D>("k0s1", "Kshort All 1" , 150, 0.35, 0.65);
  histo_k0s[2] = fs->make<TH1D>("k0s2", "Kshort All 2" , 150, 0.35, 0.65);
  histo_k0s[3] = fs->make<TH1D>("k0s3", "Kshort All 3" , 150, 0.35, 0.65);
  histo_k0s[4] = fs->make<TH1D>("k0s4", "Kshort All 4" , 150, 0.35, 0.65);
  
  histo_lamb = fs->make<TH1D>("lamb", "lamb reflexion" , 200, 1.0, 2.0); //--0.5 MeV/channel
  
  histo_closest  = fs->make<TH1D>("closest","Closest to uu", 300,-3,3);
  histo_MCgen    = fs->make<TH1D>("MCgen","MC generated event", 12,0,12);
  
  histo_NLambda0= fs->make<TH1D>("NLambda0","MC gen Lambda0 daug", 10,0,10);
  histo_L0daug  = fs->make<TH1D>("L0daug","MC gen L0daug id", 4600,-2300,2300);
  
  histo_B0Energy  = fs->make<TH1D>("B0Energy" , "B0Energy" , 300, 0, 10);
  histo_massB0    = fs->make<TH1D>("massB0"   , "massB0"   , 400, 3, 7);//--10 Mev/ch
  histo_massZ4430 = fs->make<TH1D>("massZ4430", "massZ4430", 300, 0, 10);
  histo_massKaon  = fs->make<TH1D>("massKaon" , "massKaon" , 30, 0, 1) ;
  histo_masspsip  = fs->make<TH1D>("masspsip" , "masspsip" , 150, 0, 5);
  histo_masspion  = fs->make<TH1D>("masspion" , "masspion" , 30, 0, 1) ;
  histo_massmuimuj= fs->make<TH1D>("massmuij" , "massmuij" , 150, 0, 5) ;//--33MeV/ch
  histo_massmuijpi= fs->make<TH1D>("massmuijpi" , "massmuijpi" , 600, 3, 6) ;//--5 Mev/ch
  histo_massBch   = fs->make<TH1D>("massBch" , "massBch"  , 400, 3, 7);//--10 Mev/ch
  histo_massKsh   = fs->make<TH1D>("massKsh" , "massKsh"  , 30, 0, 1);//--33 Mev/ch
  histo_combB     = fs->make<TH1D>("gen_combB   " , "gen_combB   "  , 10, 0, 10);
  histo_combpsip  = fs->make<TH1D>("gen_combpsip" , "gen_combpsip"  , 5, 0, 5);
  histo_combpi    = fs->make<TH1D>("gen_combpi  " , "gen_combpi  "  , 5, 0, 5);
  histo_combKch   = fs->make<TH1D>("gen_combKch " , "gen_combKch "  , 5, 0, 5);
  histo_combKsh   = fs->make<TH1D>("gen_combKsh " , "gen_combKsh "  , 5, 0, 5);
  histo_comb3p    = fs->make<TH1D>("gen_combpsippik " , "gen_combpsippik "  , 5, 0, 5);
  histo_primCL    = fs->make<TH1D>(" primCL   " , " primCL   "  , 1000, 0, 100);
  histo_primCLHM  = fs->make<TH1D>(" primCLHM " , " primCLHM "  , 1000, 0, 100);
  histo_primDz    = fs->make<TH1D>(" primDz   " , " primDz   "  , 1000, 0, 100);
  histo_primDxy   = fs->make<TH1D>(" primDxy  " , " primDxy  "  , 1000, 0, 100);
  histo_primDxyz  = fs->make<TH1D>(" primDxyz " , " primDxyz "  , 1000, 0, 100);
}

DimuPiK::~DimuPiK()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
DimuPiK::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
//--using namespace edm;
   
//--Tracks Collection
//edm::Handle<reco::TrackCollection> tracks;
  edm::Handle< std::vector<pat::GenericParticle> >tracks;         //--thePATTrackHandle;
  iEvent.getByLabel(trackTags_,tracks); int Ntk = tracks->size();
 
//--Get Muon Collection
//edm::Handle<reco::MuonCollection> allmuons;
  edm::Handle< std::vector<pat::Muon> >allmuons;
  iEvent.getByLabel(theMuonsLabel_, allmuons); int allMu = allmuons->size();
  
  //--Handles for Std Vees (Lambdas & Kshort) and Load them into a vector of composite candidates  
   std::vector<reco::VertexCompositeCandidate> theVees;
   edm::Handle<reco::VertexCompositeCandidateCollection> theVeeHandle;
 //iEvent.getByLabel(VeeAlgo_, "Lambda", theVeeHandle);
   iEvent.getByLabel(VeeAlgo_, "Kshort", theVeeHandle);

   theVees.insert( theVees.end(), theVeeHandle->begin(), theVeeHandle->end() );
   int allK0s  = theVees.size();
   
//--Handles for Primary Vertex  (from Eduardo)
  PrimaryInfo primary(iEvent,iConfig);
   
//--Handle fron Primary just to get the number of rec Prim vertex (Sep 10, 11)
 edm::Handle<reco::VertexCollection> privtxs;
 //iEvent.getByLabel(thePrimaryVertexLabel, privtxs); int allPrim = privtxs->size();
 //iEvent.getByLabel("PrimaryCollection", privtxs); int allPrim = privtxs->size();
 iEvent.getByLabel("offlinePrimaryVertices", privtxs); int allPrim = privtxs->size();
  
//--Handles for beamspot used for primary refit
  reco::BeamSpot beamSpot;
  edm::Handle<reco::BeamSpot> beamSpotHandle;	     
  iEvent.getByLabel("offlineBeamSpot", beamSpotHandle);
  if( beamSpotHandle.isValid() )beamSpot = *beamSpotHandle;
   
//--Handles for B-field
  edm::ESHandle<MagneticField> bFieldHandle;
  edm::ESHandle<GlobalTrackingGeometry> globTkGeomHandle;
  iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);
 //const MagneticField *BField = bFieldHandle.product();
  iSetup.get<GlobalTrackingGeometryRecord>().get(globTkGeomHandle);
   
//--Handles for Tracker Geometry
  //edm::ESHandle<GlobalTrackingGeometry> globTkGeomHandle;
  //iSetup.get<GlobalTrackingGeometryRecord>().get(globTkGeomHandle);
    
//--Handles for Transient Track Builder 
// edm::ESHandle<TransientTrackBuilder> theTTBHandle;
// iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBHandle);
 
//--HMoreno  //--Load Tracks into a vector of Tracks (from Eduardo CascadeFitter.cc)
//--HMoreno  std::vector<reco::Track> theTracks;
//--HMoreno  theTracks.insert( theTracks.end(), 
//--HMoreno		   tracks->begin(), 
//--HMoreno		tracks->end()  );
 
//--Loop over tracks to convert to transient tracks (from Eduardo CascadeFitter.cc)
  std::vector<reco::TrackRef> theTrackRefs;
  std::vector<reco::TransientTrack> theTransTracks;
  for(std::vector<pat::GenericParticle>::const_iterator itrack  = tracks->begin();itrack != tracks->end();   ++itrack){
     reco::TrackRef tmpRef = itrack->track() ;
     reco::TransientTrack tmpTk2( *tmpRef, &(*bFieldHandle), globTkGeomHandle );
     theTrackRefs.push_back( tmpRef );
     theTransTracks.push_back( tmpTk2 );
  }
  
/*
  edm::LogInfo("HecTauPAT")<<"Hec: Run: "<<evt_.runNb<<" Evt: "<<evt_.eventNb
 		<<" Lumi:"<<evt_.lumiBlock<<" Ntk: "<<evt_.Ntk<<" Mu:"<<evt_.allMu
 		<<" CalMu:"<<evt_.allCalMu
*/
     
//--Looks on Muons  [Oct 1, 2012]
  histo_trk->Fill( Ntk);
  histo_nmu->Fill( allMu );
  histo_primVtx->Fill( allPrim );
  histo_kshort->Fill( allK0s );
  
//std::vector<int> theMuonTrkIndexes;
//theMuonTrkIndexes.push_back(iMuon->innerTrack().index());

  init(0);     //--clean all structures
  HecHltTrig(iEvent);   //--Unpack Trigger

//------------------------------------MC matching Feb 4-2013--------------------------------------------

// Init parameters
 for(int ic1 = 0; ic1<11 ;ic1++)
 {
    for(int ic2 = 0; ic2<11 ;ic2++)
    {
       myP[ic1][ic2] = -9999;
    }
 }
                                 
 double massmuimuj= -1;
 double massmuijpi= -1;
 double Enmuij = -1; double Pmuijx = -1; double Pmuijy = -1; 
 double Pmuijz = -1; double Muij = -1;  

 double Enmuijpi = -1; double Pmuijpix = -1; double Pmuijpiy = -1;
 double Pmuijpiz = -1; double Muijpi   = -1;
 
 double BEnKch = -1; double BPxKch = -1; double BPyKch = -1;
 double BPzKch = -1; double BMassKch = -1;
 
 double BEnKsh = -1; double BPxKsh = -1; double BPyKsh = -1;
 double BPzKsh = -1; double BMassKsh = -1;
 
 double Enpiij = -1; double Ppiijx = -1; double Ppiijy = -1;
 double Ppiijz = -1; double Mpiij = -1;  
 
         
 if( doMC_ ){
     int Tb[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
     //std::string genParticles_="genParticlesPlusSim";  // this need edm::HepMCProduct
     std::string genParticles_="genParticles";
     
     //--get genParticles  
     edm::Handle<reco::GenParticleCollection> genParticles;
      
     iEvent.getByLabel(genParticles_, genParticles);
     if( MyPrintMC_ )
         std::cout << "Size of genParticle collection is " << genParticles->size() << std::endl;
	 
	 
     int comb[6] = {0,0,0,0,0,0}; 
       
         //--check if any of our signals were generated
     for(size_t k = 0; k < genParticles->size(); k++ ){      
         const reco::Candidate  *BCand = &(*genParticles)[ k ];            //--B0 candidate
         //-------------B0 candidate------------------------------------
         if ( abs(BCand->pdgId())==511 && abs(BCand->daughter(0)->pdgId())!=511     //select mother & daughters
	      || abs(BCand->pdgId())==521 && abs(BCand->daughter(0)->pdgId())!=521){// B0: 511, B0bar: -511 
	      
	                                                                          // B+: 521, B-: -521
              bool found_B       = false; //--B										  	      
              bool found_Psip    = false; //--uu
              bool found_Zplus   = false; //--Z+(4430)
              bool found_Kminus  = false; //--K^-
              bool found_Kplus   = false; //--K^+
              bool found_Kshort  = false; //--K^0_short
              bool found_Pionp   = false; //--pi^+       
              bool found_Pionm   = false; //--pi^-	      	      
	      
              if( MyPrintMC_ ){                                                   
                  std::cout <<" Found B " << std::endl;
                  std::cout <<" genparticle " << k << " has pdgid = " << BCand->pdgId() 
                            <<" Ndaughter: "<<BCand->numberOfDaughters()<< std::endl;
			                    
                  for(uint i = 0; i < BCand->numberOfDaughters(); i++ )
                  std::cout<<" B "<<i<<" daughterid: "<<BCand->daughter(i)->pdgId();
                  std::cout <<" "<< std::endl;
              }//--end if( MyPrintMC )
	      std::cout <<" B phys. var " << std::endl;
              
	      calcMyP(0, BCand);//--Quant. for B
	                  
              if( BCand->numberOfDaughters()==3 ){
	          
		  for( uint i = 0; i < BCand->numberOfDaughters(); i++){
		   
		        //----------Psi prime (uu) candidate---------------		  
		       
		       if( abs(BCand->daughter(i)->pdgId())== 100443 ){// found psi2s
		       
		           //comb[1]++;
			   std::cout <<" Psip phys. var " << std::endl;       
			   
			   const reco::Candidate * genDauPsip1 = BCand->daughter(i);      		   		
			   calcMyP(1, genDauPsip1);//--Quant. for Psip
			   
			   std::cout<<" Found Psi_p1 " << std::endl;			      			      
			   std::cout<<" Psi_p1 " << i << " has pdgid = " << genDauPsip1->pdgId() 
				    <<" Ndaughter: "<<genDauPsip1->numberOfDaughters()<< std::endl;
			   
			   
			   bool up0=false, un0=false;		 
			   if(found_Psip)std::cout <<" DOUBLE found : Psip  "<<" "<<k<< std::endl; 
			   	                   	
			   if( genDauPsip1->numberOfDaughters() == 2 ){
			       for( uint j = 0; j < genDauPsip1->numberOfDaughters(); j++){
				    if( MyPrintMC_ )
					 std::cout<<" Psi_p1 "<<j<<" daughterid: "
					          << genDauPsip1->daughter(j)->pdgId()
						  << std::endl;				
				    if( !up0 && genDauPsip1->daughter(j)->pdgId() == 13 ){
				         up0 = true;//--Stands for muon +
			                std::cout <<" iMuon phys. var " << std::endl;
                                        
					
					calcMyP(2, genDauPsip1->daughter(j));//--Quant. for iMuon
					

				    }//-- end if( !up0 && genDauPsip1->daughter(j)->pdgId() == 13 ){
				
				    if( !un0 && genDauPsip1->daughter(j)->pdgId() ==-13 ){
				        un0 = true;
					std::cout <<" jMuon phys. var " << std::endl;	
					
					
					calcMyP(3, genDauPsip1->daughter(j));//--Quant.fot jMuon
					
				    }//-- end if( !un0 && genDauPsip1->daughter(j)->pdgId() == -13 ){
			     
			       }//--end for ( uint j = 0; j < genDauPsip1->numberOfDaughters(); j++){
			   }//--end if( genDauPsip1->numberOfDaughters() == 2){
			   if( up0 && un0 )found_Psip = true;												      
		       }//--end if( abs(BCand.daughter(i)->pdgId())== 100443 ){
		       

		  
		       if( abs(BCand->daughter(i)->pdgId()) == 211) {// found Pion^(+/-)
		           //comb[2]++;
                           if(found_Pionm)std::cout <<" DOUBLE found : Pionm  "<<"  "<<k<< std::endl;
		           if(found_Pionp)std::cout <<" DOUBLE found : Pionp  "<<"  "<<k<< std::endl;
			   std::cout <<" Pion phys. var " << std::endl; 
		           
			   calcMyP(4, BCand->daughter(i));//--Quant. for Pion
			    
			    //if (fabs(BCand.daughter(i)->eta()) < 2.5  &&  
			    //    BCand.daughter(i)->pt() > 0.1) { accept_pion = true;}
			    //-- Changed from 0.4 to 0.1 01-16-13
			    //if( MyPrintMC_ )
			    //	 std::cout <<" Found Pion " << std::endl;     
			   if( /*!found_Pionp  &&*/ BCand->daughter(i)->pdgId()  == 211 )found_Pionp  = true;    
			   if( /*!found_Pionm  &&*/ BCand->daughter(i)->pdgId()  ==-211 )found_Pionm =  true;
		       
		       } //--end if( abs(BCand.daughter(i)->pdgId()) == 211) {
		    
		       if( abs(BCand->pdgId())==511  ){//-- If a B0 or B0bar, found Kplus/minus
		       
			   
		           //--------K^- candidate-----------------
		           // bool up0=false, un0=false;
		           if( abs(BCand->daughter(i)->pdgId())== 321 ){// Found K charged meson
			       //comb[3]++;
                               if( found_Kplus )std::cout <<" DOUBLE found : Kplus  "
			                                  <<"  "<<k<< std::endl;
		               if( found_Kminus )std::cout<<" DOUBLE found : Kminus "
			                                  <<"  "<<k<< std::endl;
			       std::cout <<" K char phys. var " << std::endl;
			       
			       				  			       
			       calcMyP(5, BCand->daughter(i));//--Quant. for K^(+/-)
			       			       
						   			   
			       if(myP[5][5] == 1) found_Kplus   = true;
			       if(myP[5][5] == -1)found_Kminus  = true;
		           }//--end if( abs(BCand.daughter(i)->pdgId())== 321 ){
		    		        
		    
		       }//end  if( abs(BCand.pdgId())==511 for K charged ){
		       
		       if( abs(BCand->pdgId())==521  ){// Found Charged B meson
		       
		           if( abs(BCand->daughter(i)->pdgId())== 310 ){// found K^0_short
			       //comb[4]++;
		               if(found_Kshort)std::cout <<" DOUBLE found : Kshort "
			                                 <<"  "<<k<< std::endl;			    			   
			   			       
			       const reco::Candidate * genDauKsh = BCand->daughter(i);
		 
        		       std::cout <<" Found K_short " << std::endl;			      			       			  			      
			       std::cout <<" K_short " << i << " has pdgid = " << genDauKsh->pdgId() 
				         <<" Ndaughter: "<<genDauKsh->numberOfDaughters()<< std::endl;
			       std::cout <<" K short phys. var " << std::endl;
			       
			       calcMyP(6, genDauKsh);//--Quant. for K^0_s
			       		
			       if( genDauKsh->numberOfDaughters() == 2 ){
			           bool pp0=false, pn0=false;
			           for( uint j = 0; j < genDauKsh->numberOfDaughters(); j++){
				        if( MyPrintMC_ )
					    std::cout<<" K_short "<<j<<" daughterid: "
					             <<genDauKsh->daughter(j)->pdgId()
						     << std::endl;				
				        if( !pp0 && genDauKsh->daughter(j)->pdgId() == 211 ){
				            pp0 = true;//--pion +
			                    std::cout <<" iPion phys. var " << std::endl;
					    
                                            calcMyP(7, genDauKsh->daughter(j));
					    
				        }//end if( !pp0 && genDauKsh->daughter(j)->pdgId() == 211 ){
							     
			                if( !pn0 && genDauKsh->daughter(j)->pdgId() ==-211 ){
				            pn0 = true;//-- pion -	
			                    std::cout <<" jPion phys. var " << std::endl;
					    
                                            calcMyP(8, genDauKsh->daughter(j));

				        }//--end if( !pn0 && genDauKsh->daughter(j)->pdgId() ==-211 ){
			     
			           }//--end for ( uint j = 0; j < genDauKsh->numberOfDaughters(); j++){
				   if( pp0 && pn0 )found_Kshort = true;
			       
			       }//--end if( genDauKsh->numberOfDaughters() == 2){			       												      
		           }//--end if( abs(BCand.daughter(i)->pdgId())== 310 ){//-- Found K short

		       }// end if( abs(BCand.pdgId())==521 B(+/-) )
		  }//--end for( uint i = 0; i < BCand.numberOfDaughters(); i++)
		  
		  if( found_Psip && (found_Pionm || found_Pionp) && (found_Kplus || found_Kminus)){
		  
	              found_B = true;
		      
		      	                            
	          }//end if( found_Psip && (found_Pionm || found_Pionp) && (found_Kplus || found_Kminus)){
		  
              }//--end if(BCand.numberOfDaughters()==3)
	      
	      if( found_B ){
              
	          comb[0]++;
                  //if( found_Psip ){
      
	              Enmuij = myP[2][0] + myP[3][0];
                      Pmuijx = myP[2][1] + myP[3][1];
                      Pmuijy = myP[2][2] + myP[3][2];
	              Pmuijz = myP[2][3] + myP[3][3];
       
                      Muij   = sqrt( Enmuij*Enmuij 
		                   - Pmuijx*Pmuijx - Pmuijy*Pmuijy - Pmuijz*Pmuijz);
		    
                      histo_massmuimuj->Fill(Muij);		     
      
 
	              //if( found_Pionm || found_Pionp){
      
	                  Enmuijpi = Enmuij + myP[4][0];
	                  Pmuijpix = Pmuijx + myP[4][1];
	                  Pmuijpiy = Pmuijy + myP[4][2];
	                  Pmuijpiz = Pmuijz + myP[4][3];
       
	                  Muijpi   = sqrt( Enmuijpi*Enmuijpi 
			                 - Pmuijpix*Pmuijpix - Pmuijpiy*Pmuijpiy - Pmuijpiz*Pmuijpiz  );
			  
		          histo_massmuijpi->Fill(Muijpi);
                      //}//--end if( found_Pionm || found_Pionp){
		  
	          //}//--end if( found_Psip ){
	      }//--end if( found_B ){
	      	      	      	      
         }//--end if ( abs(BCand.pdgId())==511 || abs(BCand.pdgId())==521, B0, B(+/-) )	 
     }//--end for( size_t k = 0; k < genParticles->size(); k++ ) 
     histo_combB->Fill(comb[0]);
//      if( found_Psip ){
//      
//          Enmuij = myP[2][0] + myP[3][0];
// 	 Pmuijx = myP[2][1] + myP[3][1];
// 	 Pmuijy = myP[2][2] + myP[3][2];
//          Pmuijz = myP[2][3] + myP[3][3];
// 	 
// 	 Muij   = sqrt( Enmuij*Enmuij - Pmuijx*Pmuijx
// 	              - Pmuijy*Pmuijy - Pmuijz*Pmuijz);
// 		      
// 	 histo_massmuimuj->Fill(Muij);	               
//      
// 
//          if( found_Pionm || found_Pionp){
//      
//              Enmuijpi = Enmuij + myP[4][0];
// 	     Pmuijpix = Pmuijx + myP[4][1];
// 	     Pmuijpiy = Pmuijy + myP[4][2];
//              Pmuijpiz = Pmuijz + myP[4][3];
// 	 
// 	     Muijpi   = sqrt( Enmuijpi*Enmuijpi - Pmuijpix*Pmuijpix
// 	                - Pmuijpiy*Pmuijpiy - Pmuijpiz*Pmuijpiz  );
//      
//              histo_massmuijpi->Fill(Muijpi);
// 	     
//              if( found_Kplus || found_Kminus ){//--B^0 
//      
//                  BEnKch   = Enmuijpi + myP[5][0]; 
// 	         BPxKch   = Pmuijpix + myP[5][1]; 
// 	         BPyKch   = Pmuijpiy + myP[5][2];
//                  BPzKch   = Pmuijpiz + myP[5][3]; 
// 	         BMassKch = sqrt(BEnKch*BEnKch
// 		            - BPxKch*BPxKch - BPyKch*BPyKch - BPzKch*BPzKch);
// 			    
// 		 histo_massB0    ->Fill(BMassKch);	      
// 			      
// 	     }//--end if( found_Kplus || found_Kminus ){ 
// 	    
//              if( found_Kshort){//--B^(+/-) 
//      
//                  BEnKsh   = Enmuijpi + myP[6][0]; 
// 	         BPxKsh   = Pmuijpix + myP[6][1]; 
// 	         BPyKsh   = Pmuijpiy + myP[6][2];
//                  BPzKsh   = Pmuijpiz + myP[6][3]; 
// 	         BMassKsh = sqrt(BEnKsh*BEnKsh - BPxKsh*BPxKsh 
// 		              - BPyKsh*BPyKsh - BPzKsh*BPzKsh);
// 	    }//--end if( found_Kshort ){
// 	    	    	    
//          }//--end if( found_Pionm || found_Pionp){
//      
//      }//end if( found_Psip ){
//      
//      if( found_Kshort ){
//      
//          Enpiij = myP[7][0] + myP[8][0];
// 	 Ppiijx = myP[7][1] + myP[8][1];
// 	 Ppiijy = myP[7][2] + myP[8][2];
//          Ppiijz = myP[7][3] + myP[8][3];
// 	 
// 	 Mpiij   = sqrt( Enpiij*Enpiij - Ppiijx*Ppiijx
// 	              - Ppiijy*Ppiijy - Ppiijz*Ppiijz);         
//      
//      }
//      
//      
//      
//      histo_massB0    ->Fill(BMassKch);
//      histo_massBch   ->Fill(BMassKsh);
//      histo_massKsh   ->Fill(Mpiij);
     
                               	      
     //if( MyPrintMC_ ){	 
         //if( found_Psip	)std::cout<<found_Psip  <<" Found psi(2S): "<<std::endl;
	     //else	 std::cout<<found_Psip  <<" NO Found psi(2S): "<<std::endl;
	      
	      //if( found_Kminus   )std::cout<<found_Kminus <<" Found K-: "<<std::endl;
	      //else                std::cout<<found_Kminus <<" NO Found K-: "<<std::endl;
	      
	      //if( found_Kplus	)std::cout<<found_Kplus <<" Found K+: "<<std::endl;
	      //else		 std::cout<<found_Kplus <<" NO Found K+: "<<std::endl;
	      
	      //if( found_Pionp )std::cout<<found_Pionp <<" Found Pion+: "<<std::endl;
	      //else	       std::cout<<found_Pionp <<" NO Found Pion+: "<<std::endl;
    
	      //if( found_Pionm )std::cout<<found_Pionm <<" Found Pion-: "<<std::endl;
	      //else	       std::cout<<found_Pionm <<" NO Found Pion-: "<<std::endl;
	  
	      //if( found_Zplus )std::cout<<found_Zplus <<" Found Z4430: "<<std::endl;
	      //else		 std::cout<<found_Zplus <<" NO Found Z4430: "<<std::endl;
 	    
         
     //}//--end if( MyPrintMC )
     
//      histo_combB   ->Fill(comb[0]);
//      histo_combpsip->Fill(comb[1]);
//      histo_combpi  ->Fill(comb[2]);
//      histo_combKch ->Fill(comb[3]);
//      histo_combKsh ->Fill(comb[4]);
//      histo_comb3p  ->Fill(comb[5]);
         
         //-----Counting Generated Bmeson ------        
         
     //if( found_Psip  && found_Kminus && found_Pionp){
     //    Lb[1]++;
     //    Tb[1]++;
     //}
     
     //if( found_Psip  && found_Kplus  && found_Pionm)Lb[2]++,Tb[2]++;             
     //    Lb[3] = Lb[1] + Lb[2], Tb[3] = Tb[1] + Tb[2];
         
         //if( found_Psip1 && found_LambP_0 )Lb[4]++,Tb[4]++;
         //if( found_Psip1 && found_LambA_0 )Lb[5]++,Tb[5]++;
         //Lb[6] = Lb[4] + Lb[5], Tb[6] = Tb[4] + Tb[5];
         
         //if( found_Psip2 && found_LambP_0 )Lb[7]++,Tb[7]++;
         //if( found_Psip2 && found_LambA_0 )Lb[8]++,Tb[8]++;
         //Lb[9] = Lb[7] + Lb[8], Tb[9] = Tb[7] + Tb[8];
         
     for( uint i = 1; i < 10; i++)
          histo_MCgen->Fill( i,Tb[i]); 
  
     if( MyPrintMC_ )
         std::cout <<"||||||||>>>>>>>end of event<<<<<<<||||||||"<< std::endl;
  }//--end if( doMC_)
 //---------------------------------------------------------------------------------------------

    
//--Constant for Vertexing
    float muon_sig   = muon_mass_c*1.e-6;
    float pion_sig   = pion_mass_c*1.e-6;
    float kaon_sig   = kaon_mass_c*1.e-6;
    float chi = 0., ndf = 0.;
    
//--Start Muon Loop     std::cout<<"Hec: Starting diMuon Loop"<<std::endl;
int Comb_uu=0;
for(std::vector<pat::Muon>::const_iterator iMu=allmuons->begin(); iMu!=allmuons->end();++iMu){
   reco::TrackRef iMuTrackRef = iMu->innerTrack(); //--.index();
   HecMuVar(0,iMu,iMuTrackRef);                    //--Variables for Muon i  (0)
   
   histo_muId->Fill( nMuVar[0][1] );
   
   for(std::vector<pat::Muon>::const_iterator jMu=iMu+1; jMu!=allmuons->end();++jMu){
      reco::TrackRef jMuTrackRef = jMu->innerTrack();	//--.index();
      HecMuVar(1,jMu,jMuTrackRef);			//--Variables for Muon j  (1)

      bool goodDimuon = (nMuVar[0][1] == 111 && nMuVar[1][1] >= 11)    //--1 G and 1 (SG or T or ST or G)
                      ||(nMuVar[1][1] == 111 && nMuVar[0][1] >= 11);

      if( goodDimuon && nMuVar[0][0]*nMuVar[1][0]<0 ){         //--Take only good muons and Opp. sign

        //--Save Dimuon P4 in a LorentzVector
        double uuPx = iMu->p4().Px() + jMu->p4().Px();
        double uuPy = iMu->p4().Py() + jMu->p4().Py();
        double uuPz = iMu->p4().Pz() + jMu->p4().Pz();
        double uuEn = sqrt( muon_mass_c*muon_mass_c + iMu->p4().P()*iMu->p4().P() )
	            + sqrt( muon_mass_c*muon_mass_c + jMu->p4().P()*jMu->p4().P() );
        reco::Particle::LorentzVector uuP4(0.0,0.0,0.0,0.0);                          
        uuP4.SetPxPyPzE(uuPx,uuPy,uuPz,uuEn);
      
        double Muu = ( iMu->p4() +  jMu->p4() ).M();           //--Dimuon  mass
        histo_dimu[0]->Fill( Muu );
	HecCutMu(0,iMu);
        HecCutMu(1,jMu);
	
	//bool uuMcut = std::abs( Muu - jpsi_mass_c ) < 0.2;  //--200 MeV
        bool uuMcut = ( Muu>3.45 && Muu<4.0 );                 //--for psi'
                //--||( Muu>0.50 && Muu<1.2 );                //--for eta, rho/omega and Phi
        if( uuMcut ){                                         //--Skip if not Dimuon Selected
	    
	//--do Here dimuon Vertex
	reco::TransientTrack iMuTT( iMuTrackRef, &(*bFieldHandle) );
	reco::TransientTrack jMuTT( jMuTrackRef, &(*bFieldHandle) );
     
	KinematicParticleFactoryFromTransientTrack pFactory;
	 
	std::vector<RefCountedKinematicParticle> uuParticle;
	uuParticle.push_back(pFactory.particle( iMuTT, muon_mass_c, chi, ndf, muon_sig ));    //--Muon i
	uuParticle.push_back(pFactory.particle( jMuTT, muon_mass_c, chi, ndf, muon_sig ));    //--Muon j

	//--Trajectory states to calculate DCA for the 2 tracks
	FreeTrajectoryState iState = pFactory.particle( iMuTT, muon_mass_c, chi, ndf, muon_sig )->currentState().freeTrajectoryState();
	FreeTrajectoryState jState = pFactory.particle( jMuTT, muon_mass_c, chi, ndf, muon_sig )->currentState().freeTrajectoryState();

	if( HecClosest(0, iState, jState) ){   //--iMu close to jMu
            histo_dimu[1]->Fill( Muu );

            KinematicFitDriver uuRec( uuParticle, "MuMu" );//--HM 03/05/13
	    
	    //if( HecDimuVtx(Muu, uuParticle, uuP4) ){      //--Do here DiMuon Vertex //--HM 03/05/13
            if( HecDimuVtx(Muu, uuRec, uuP4) ){      //--Do here DiMuon Vertex
            histo_dimu[2]->Fill( Muu ); 

            //--Trigger matching for i & j muons
	    const pat::Muon *patMuonP = &(*iMu);
	    const pat::Muon *patMuonM = &(*jMu);
            HecHltMuTrig(0, patMuonP);
            HecHltMuTrig(1, patMuonM);
	                
	    Comb_uu++;
	    //--Add here 2 extra Pions to make Psi(2s)-->J/Psi Pi+Pi-
	    //--
	    
	    //--add here bachelor Pion track
	    int Comb_kPi=0;
            for(unsigned int ktrack = 0; ktrack < theTrackRefs.size(); ktrack++){       //--Pion k   
               reco::TrackRef kPiTrk = theTrackRefs[ktrack];
	       bool sameTrk_kPi_ij = kPiTrk==iMuTrackRef || kPiTrk==jMuTrackRef;
	       if( !sameTrk_kPi_ij ){                                                //--Skip same Trk

	       int kPi = 0;
	       HecTrkVar(kPi,pion_mass_c, kPiTrk);        //--Save kPi variables
	       
               reco::Particle::LorentzVector kPiP4(0.0,0.0,0.0,0.0);                     //--Fill the 4-Momentum Vector for kPi 
               kPiP4.SetPxPyPzE( TrkVar[kPi][1], TrkVar[kPi][2], TrkVar[kPi][3], TrkVar[kPi][4] );
	       
               double MuukPi    = (iMu->p4() + jMu->p4() + kPiP4 ).M();  //--Mass: iMu jMu kPi (3 body)(un-fitted muons)
               double MuikPi    = (iMu->p4() + kPiP4 ).M();              //--Mass: iMu kPi (2 body)
               double MujkPi    = (jMu->p4() + kPiP4 ).M();              //--Mass: jMu kPi (2 body)
               double MuukPiFit = (uuP4fit + kPiP4).M();                 //--uuFitted + Pion
	       
	       histo_dimuPi[0]->Fill( MuukPiFit );
	       
	       //--Require Pion close to Dimuon and cut on Pt
               double Ruu_kPi = HecRcone(kPi, uuVtx[25], uuVtx[24]);
	       
               bool good_kPi = TrkVar[kPi][5]>0.25     //--kPi Transverse Momentum 
	                    && Ruu_kPi<1;              //--R dimuon-kPi
               if( good_kPi ){
               histo_dimuPi[1]->Fill( MuukPiFit );

	       //--closest approach between muons and kPion
	       reco::TransientTrack kPiTT( *kPiTrk   , &(*bFieldHandle) );
	       
	       FreeTrajectoryState kPiState = pFactory.particle( kPiTT, pion_mass_c, chi, ndf, pion_sig )->currentState().freeTrajectoryState();

	       if( HecClosest(1, kPiState, iState) ){   //--kPi close to iMu close
	       if( HecClosest(2, kPiState, jState) ){   //--kPi close to jMu close
	       
	       bool goodVtx3Body = false;
	       double MuupFit =-1;
	       double uupV3_CL = 999;
               reco::Particle::LorentzVector  uupV3_P4;
	       reco::Particle::Point	      uupV3_vx;
	       reco::Vertex::CovarianceMatrix uupV3_Cov;
	       if( NeutralKa_ ){    //--do 3 body (mu+mu-pi+) Vertex
		   std::vector<RefCountedKinematicParticle> uupParticle;
		   uupParticle.push_back(pFactory.particle( iMuTT, muon_mass_c, chi, ndf, muon_sig ));    //--Muon i
		   uupParticle.push_back(pFactory.particle( jMuTT, muon_mass_c, chi, ndf, muon_sig ));    //--Muon j
		   uupParticle.push_back(pFactory.particle( kPiTT, pion_mass_c, chi, ndf, pion_sig  ));   //--Pion k
		  
		KinematicFitDriver uupRec( uupParticle, "MuMuMp" );
		if( uupRec.isValid() ){ 		 //--Valid 3 body vertex
		       goodVtx3Body = true;

		    //--Output Vertex fit variables
		       MuupFit  				= uupRec.mass();	  //--3 body (Dimuon Pi) mass
                                                      uupV3_P4  = uupRec.P4();            //--DiMuPi 4-Momentum
		       reco::Particle::LorentzVector  iMuV3_P4  = uupRec.P4FirstChild();  //--Mui
		       reco::Particle::LorentzVector  jMuV3_P4  = uupRec.P4NextChild();   //--Muj 
		       reco::Particle::LorentzVector  kPiV3_P4  = uupRec.P4NextChild();   //--Pion
		                                      uupV3_vx  = uupRec.VertexDecay();   //--Vx,Vy,Vz 
		                                      uupV3_Cov = uupRec.VertexCov();	      
		       double uupV3_Chi2			= uupRec.chi2();
		       double uupV3_Ndof			= uupRec.ndof();
		              uupV3_CL = ChiSquaredProbability( uupV3_Chi2, uupV3_Ndof );  //--Confidence Level uup Vertex	       
		}   //--End if( uupRec.isValid() ){		     //--Valid 3 body vertex
	       }
	       Comb_kPi++;
	       
	       int Comb_kKa=0;
               if( ChargedKa_ ){                   //--Add here Charged Kaon track
               for(unsigned int ntrack = 0; ntrack < theTrackRefs.size(); ntrack++){     //--kaon k  
                  reco::TrackRef kKaTrk = theTrackRefs[ntrack];
	          bool sameTrk_kKa_ijkPi  = kKaTrk==iMuTrackRef || kKaTrk==jMuTrackRef
		                         || kKaTrk==kPiTrk;
                  if( !sameTrk_kKa_ijkPi ){                 //--Skip same tracks
		  
		  int kKa = 1;
	          HecTrkVar(kKa,kaon_mass_c, kKaTrk);
		  
                  reco::Particle::LorentzVector kKaP4(0.0,0.0,0.0,0.0);                     //--Fill the 4-Momentum Vector for kKa
                  kKaP4.SetPxPyPzE( TrkVar[kKa][1], TrkVar[kKa][2], TrkVar[kKa][3], TrkVar[kKa][4] );

                  double Muupk   = (uuP4fit + kPiP4 + kKaP4 ).M();    //--Mass: iMu jMu kPi kKa (4 body)(Fitted muons)
                  double MkPikKa = (kPiP4 + kKaP4 ).M();
		  
                  //--Select good Pion-Kaon combination (too many)
                  double Ruu_kKa = HecRcone(kKa, uuVtx[25], uuVtx[24]);
		    
		  bool goodPair_kPikKa = TrkVar[kPi][0]*TrkVar[kKa][0]<0              //--k-pion opposite charge
		                  //--&& Ruu_kKa<1                                    //--k close to uu
				      && TrkVar[kKa][5]>0.25                          //--kKa->pt>0.25
				      &&(TrkVar[kKa][5] > TrkVar[kPi][5]);            //--momentun cut pt_kaon > pt_pion
                  if( goodPair_kPikKa ){
		  
	          //--closest approach between muons and kPion
                  reco::TransientTrack kKaTT( *kKaTrk    , &(*bFieldHandle) );
		  
                  FreeTrajectoryState kKaState = pFactory.particle( kKaTT, kaon_mass_c, chi, ndf, kaon_sig )->currentState().freeTrajectoryState();
		  
		  if( HecClosest(3, kKaState, iState) ){   //--kKa close to iMu close
	          if( HecClosest(4, kKaState, jState) ){   //--kKa close to jMu close
	          if( HecClosest(5, kKaState, kPiState) ){ //--kKa close to kPi close

	          histo_dimuPi[2]->Fill( MuukPiFit );
		  if( Muupk>5.1 && Muupk<5.5 )    //--B candidate
                      histo_dimuPi[3]->Fill( MuukPiFit );

                  //--Do Vertex  MuMu + Pi + Ka
                  std::vector<RefCountedKinematicParticle> uupkParticle;
                  uupkParticle.push_back(pFactory.particle( iMuTT, muon_mass_c, chi, ndf, muon_sig ));    //--Muon i
                  uupkParticle.push_back(pFactory.particle( jMuTT, muon_mass_c, chi, ndf, muon_sig ));    //--Muon j
                  uupkParticle.push_back(pFactory.particle( kPiTT, pion_mass_c, chi, ndf, pion_sig  ));   //--Pion k
                  uupkParticle.push_back(pFactory.particle( kKaTT, kaon_mass_c, chi, ndf, kaon_sig  ));   //--Kaon k
                  
		  KinematicFitDriver uupkRec( uupkParticle, "MuMuMpMk" );
	          if( uupkRec.isValid() ){                  //--Valid 4 body vertex

	          //--Output Vertex fit variables
                  double MuupkFit                          = uupkRec.mass();           //--4 body (Dimuon Pi ka) mass
                  reco::Particle::LorentzVector  uupkV_P4  = uupkRec.P4();             //--DiMuPiKa 4-Momentum
                  reco::Particle::LorentzVector  iMuV_P4   = uupkRec.P4FirstChild();   //--Mui
                  reco::Particle::LorentzVector  jMuV_P4   = uupkRec.P4NextChild();    //--Muj
                  reco::Particle::LorentzVector  kPiV_P4   = uupkRec.P4NextChild();    //--Pion
                  reco::Particle::LorentzVector  kKaV_P4   = uupkRec.P4NextChild();    //--Kaon
                  reco::Particle::Point	         uupkV_vx  = uupkRec.VertexDecay();    //--Vx,Vy,Vz 
                  reco::Vertex::CovarianceMatrix uupkV_Cov = uupkRec.VertexCov();					  
                  double uupkV_Chi2                        = uupkRec.chi2();
                  double uupkV_Ndof                        = uupkRec.ndof();
                  double uupkV_CL = ChiSquaredProbability( uupkV_Chi2, uupkV_Ndof );  //--Confidence Level uupk Vertex		  
		       
                  //--Do Primary here.  Primary Vertex Re-fitting to exclude the DiMuon, Pion  and Kaon from the primary
                  reco::Vertex refitPrimVertex = *primary.BestVertex();
                  std::vector<reco::TrackRef> ExclusionList;
                  ExclusionList.push_back( iMuTrackRef );
                  ExclusionList.push_back( jMuTrackRef );
                  ExclusionList.push_back( kPiTrk );
                  ExclusionList.push_back( kKaTrk );
                  
		  
                  //VertexRefit RefitPrimary( primary.BestVertex(), ExclusionList, bFieldHandle, beamSpot );
                  //refitPrimVertex = RefitPrimary.Refitted();
		  
		  //--HMoreno 03/04/13-----------------------------------------------------------------------
		  VertexRefit RefitPrimaryHM(primary.BestVertex() ,ExclusionList ,bFieldHandle ,beamSpot );    //--Highest Multiplicity Primary 
                  reco::Vertex refitVertexPrimHM = RefitPrimaryHM.Refitted(); 
                  double primCLHM   = ChiSquaredProbability( refitVertexPrimHM.chi2(), refitVertexPrimHM.ndof());
		  
		  GlobalPoint  secver = uuRec.RKVertex()->position();
                  //-GlobalVector secmom = uuRec.RKParent()->currentState().globalMomentum();       //--only dimuon       
                  GlobalVector secmomTot = GlobalVector(uupkRec.P4().Px()           //--this include dimuon, pion, kaon 
                                                       ,uupkRec.P4().Py()
                                                       ,uupkRec.P4().Pz());
						       
                  if( MyPrint_ ){
                  //std::cout<<"B: vector P="<<secmom<<std::endl;
                  std::cout<<"        : uupkRec   ="<<uupkRec.P4().Px()<<std::endl;
                  std::cout<<"        : ---->     ="<<secmomTot<<std::endl;
                  }
		  
	          //-VertexRefit RefitPrimary(primary.HigherCosAlpha(secver,secmom) ,ExclusionList ,bFieldHandle ,beamSpot ); //--pointing
                  VertexRefit RefitPrimary(primary.HigherCosAlpha(secver,secmomTot) ,ExclusionList ,bFieldHandle ,beamSpot ); //--pointing
                  reco::Vertex refitVertexPrim = RefitPrimary.Refitted();
                  double primCL   = ChiSquaredProbability( refitVertexPrim.chi2(), refitVertexPrim.ndof());
                 
                  //--std::cout<<" CL: High Multiplicity ="<<primCLHM<<" Pointing Angle ="<<primCL<<" ->"<<primCLHM-primCL<<std::endl;
              
                  //--GlobalPoint PosVertex = GlobalPoint(refitVertexPrim.x(), refitVertexPrim.y(), refitVertexPrim.z() );

		  double primDist = sqrt( refitVertexPrim.x()*refitVertexPrim.x()
                                        + refitVertexPrim.y()*refitVertexPrim.y()
                                        + refitVertexPrim.z()*refitVertexPrim.z() );
                  double primR    = sqrt( refitVertexPrim.x()*refitVertexPrim.x()
                                        + refitVertexPrim.y()*refitVertexPrim.y() );
 
                  double primDz   = refitVertexPrim.z() - refitVertexPrimHM.z();
                  double primDxy  = sqrt( (refitVertexPrim.x() - refitVertexPrimHM.x())*(refitVertexPrim.x() - refitVertexPrimHM.x())
                                        + (refitVertexPrim.y() - refitVertexPrimHM.y())*(refitVertexPrim.y() - refitVertexPrimHM.y()) );
                  double primDxyz= sqrt(  (refitVertexPrim.x() - refitVertexPrimHM.x())*(refitVertexPrim.x() - refitVertexPrimHM.x())
                                        + (refitVertexPrim.y() - refitVertexPrimHM.y())*(refitVertexPrim.y() - refitVertexPrimHM.y())
                                        + (refitVertexPrim.z() - refitVertexPrimHM.z())*(refitVertexPrim.z() - refitVertexPrimHM.z()) );
              
                  if( MyPrint_ ){
                      std::cout<<"B0: primCL ="<<primCL<<" High Mult: "<<primCLHM<<std::endl;
                      std::cout<<"B0: primDz ="<<primDz<<" primDxy =  "<<primDxy <<" primDxyz ="<<primDxyz<<std::endl;
                  }
                  histo_primCL  ->Fill(primCL)  ;
		  histo_primCLHM->Fill(primCLHM);
		  histo_primDz  ->Fill(primDz)  ;
		  histo_primDxy ->Fill(primDxy) ;
		  histo_primDxyz->Fill(primDxyz);		  
		  
		  //-- end HMoreno 03/04/13-----------------------------------------------------------------
		  
		  
                  if( refitPrimVertex.isValid() ){                    //--Skip no valid Primary

	          //--Output Prim Vertex fit variables 			       
                  double pV_Chi2 = refitPrimVertex.chi2();
                  double pV_Ndof = refitPrimVertex.ndof();
                  double pV_CL   = ChiSquaredProbability( pV_Chi2, pV_Ndof ); //--CLl Primary

                  //--Check separation from primary (L over Sigma)
                  typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
                  typedef ROOT::Math::SVector<double, 3> SVector3;

                  SMatrixSym3D totalCov = refitPrimVertex.covariance() + uupkV_Cov;
                  SVector3 uupk_PrimarySep3D( uupkV_vx.x() - refitPrimVertex.x(),
                                              uupkV_vx.y() - refitPrimVertex.y(),
                                              uupkV_vx.z() - refitPrimVertex.z() );

                  double uupkL     = ROOT::Math::Mag( uupk_PrimarySep3D );  //--Flag 3dim distance Dimuon-Prim
                  double uupkSigma = sqrt(ROOT::Math::Similarity( totalCov, uupk_PrimarySep3D )) / uupkL;
                  double uupkLoS = uupkL / uupkSigma;  //--Flag 3dim distance/sigma  Dimuon-Primary

                  //--Here Pointing Secondary Vertex  [Sep 9, 2011]
                  GlobalVector uupk_PrimSec( uupkV_vx.x() - refitPrimVertex.x(),
                                             uupkV_vx.y() - refitPrimVertex.y(),
                                             uupkV_vx.z() - refitPrimVertex.z() );	  
                  GlobalVector uupk_P3( uupkRec.P4().x(),
                                        uupkRec.P4().y(),
                                        uupkRec.P4().z() );	    
                  double V_dot_Puupk =uupk_P3.dot(uupk_PrimSec);		      
                  double cos_alpha = V_dot_Puupk /( uupk_PrimSec.mag() * uupk_P3.mag() );
		  
                  histo_dimuPiK[0]->Fill( MuupkFit );     //--B candidate
		  if( cos_alpha>0.9 ){          //--Cosine Sec to Prim cut

		  histo_dimuPiK[1]->Fill( MuupkFit );
                  histo_dimuPi[4]->Fill( MuukPi );
		  if( MuupkFit>5.1 && MuupkFit<5.5 && cos_alpha>0.9)    //--B candidate
                      histo_dimuPi[5]->Fill( MuukPiFit );

	          //--fill n-tuple
		  Comb_kKa++;
                  fill_evt( iEvent, Ntk, allMu, allPrim, allK0s,
		                    Comb_uu, Comb_kPi, Comb_kKa);
                  fill_IuuC();
                  fill_uuC(Muu);
		  
		  fill_kPi(kPi);
		  fill_kKa(kKa);
                  fill_vtxAll(MuukPiFit, MuukPi, MuikPi,  MujkPi,   Ruu_kPi,
                              Muupk, MkPikKa, Ruu_kKa, MuupkFit, uupkV_CL,
                              pV_CL,  uupkLoS, uupkL,    cos_alpha);
                  fill_tightMuons();
                  fill_softMuons();
	          hecmu_tree_->Fill();
	          init(1);   //--don't want to clean the counter Comb
                  }    //--End if( cos_alpha>0.9 ){          //--Cosine Sec to Prim cut 
                  }    //--End if( refitPrimVertex.isValid() ){            //--Skip no valid Primary
                  }    //--End if( uuPiKaRec.isValid() ){               //--Valid 4 body vertex
                  }    //--End if( HecClosest(5, kKaState, kPiState) ){ //--kKa close to kPi close
                  }    //--End if( HecClosest(4, kKaState, jState) ){   //--kKa close to jMu close
                  }    //--End if( HecClosest(3, kKaState, iState) ){   //--kKa close to iMu close
		  }    //--End if( goodPair_kPikKa ){
                  }    //--End if( !sameTrk_kKa_ijkPi ){                //--Skip same tracks
               }       //--End for(unsigned int ntrack = 0; ntrack < theTrackRefs.size(); ntrack++){     //--kaon k
	       }       //--End if( ChargedKa_ )
	       
	       if( NeutralKa_ ){       //--Add here Kshort-->Pi+Pi-
	         for( unsigned int veeindex = 0; veeindex < theVees.size(); veeindex++ ){
                    double Mk0s = theVees[veeindex].mass();
                    histo_k0s[0]->Fill( Mk0s );

                    reco::TrackRef nPiTrk = (dynamic_cast<reco::RecoChargedCandidate *> (theVees[veeindex].daughter(0) ) )->track();
                    reco::TrackRef rPiTrk = (dynamic_cast<reco::RecoChargedCandidate *> (theVees[veeindex].daughter(1) ) )->track();
		    
		    bool sameTrkn = nPiTrk==iMuTrackRef || nPiTrk==jMuTrackRef || nPiTrk==kPiTrk;
		    bool sameTrkr = rPiTrk==iMuTrackRef || rPiTrk==jMuTrackRef || rPiTrk==kPiTrk;		  
                    if( !sameTrkn && !sameTrkr ){                 //--Skip same tracks

                    histo_k0s[1]->Fill( Mk0s );
                    int nPi = 2;
	            HecTrkVar(nPi,pion_mass_c, nPiTrk);        //--Save kPi variables
		    int rPi = 3;
	            HecTrkVar(rPi,pion_mass_c, rPiTrk);

                    //--Reconstruct K0short Vertex
                    reco::TransientTrack nPiTT( nPiTrk, &(*bFieldHandle) );
                    reco::TransientTrack rPiTT( rPiTrk, &(*bFieldHandle) );

                    //--KinematicParticleFactoryFromTransientTrack pFactory;
                    std::vector<RefCountedKinematicParticle> K0sParticle;
                    K0sParticle.push_back(pFactory.particle( nPiTT, pion_mass_c, chi, ndf, pion_sig ));
                    K0sParticle.push_back(pFactory.particle( rPiTT, pion_mass_c, chi, ndf, pion_sig ));

                    KinematicFitDriver K0sRec( K0sParticle, "Kshort" );
                    if( K0sRec.isValid() ){
		    
                    histo_k0s[2]->Fill( Mk0s );
		    
                    //--Output lambda Vtx fit variables
                    double Mk0sVFit   = K0sRec.mass();  //--valid vertex Kshort Mass
                    double k0sV_CL    = ChiSquaredProbability( K0sRec.chi2(), K0sRec.ndof() );
  
                    //--Here do Lambda0 Mass constrains Fit. (proton + pion) to Lambda
                    //--bool k0sMcut = std::abs( Mk0sVFit - kshort_mass_c ) < 0.01;      //--10 MeV  around Nominal kshort
                    K0sRec.AddMassConstraint( kshort_mass_c, kshort_mass_c*1.e-6 );
                    if( K0sRec.isValid() ){
		    
                    //--Output lambda mass constrain fit variables
                    double Mk0sMFit = K0sRec.mass();                                          //--Valid vertex K0s Mass 
                    double k0sVM_CL = ChiSquaredProbability( K0sRec.chi2(), K0sRec.ndof() );  //--Confidence Level K0s Vtx 
                    reco::Particle::LorentzVector k0sVM_P4  = K0sRec.P4();                    //--kshort 4-Momentum with Mass constrains
                    reco::Particle::LorentzVector nPiVM_P4  = K0sRec.P4FirstChild();          //--pion n Momentum with Vtx constrains
                    reco::Particle::LorentzVector rPiVM_P4  = K0sRec.P4NextChild();           //--pion r	 
                    reco::Particle::Point          k0sVM_vx = K0sRec.VertexDecay();           //--Vx,Vy,Vz
                    reco::Vertex::CovarianceMatrix k0sVM_Cov= K0sRec.VertexCov();
                    double k0sVMx = k0sVM_vx.x();
                    double k0sVMy = k0sVM_vx.y();
                    double k0sVMz = k0sVM_vx.z();
                    double k0sDist= sqrt( k0sVMx*k0sVMx + k0sVMy*k0sVMy + k0sVMz*k0sVMz );
                    double k0sR   = sqrt( k0sVMx*k0sVMx + k0sVMy*k0sVMy );
		    
                    histo_k0s[3]->Fill( Mk0sMFit );

                    //--Calculate Lambda_0 invariant mass
		    double lamEpr = 0;
		    double lamEpi = 0;
		    if( nPiVM_P4.Pt() > rPiVM_P4.Pt() ){
		        lamEpr = sqrt( proton_mass_c*proton_mass_c + nPiVM_P4.P()*nPiVM_P4.P() );  //--Assuming Proton for nPi
			lamEpi = sqrt( pion_mass_c*pion_mass_c + rPiVM_P4.P()*rPiVM_P4.P() );      //--rPiVM_P4.E();
		    } else {
		        lamEpr = sqrt( proton_mass_c*proton_mass_c + rPiVM_P4.P()*rPiVM_P4.P() );  //--Assuming Proton for rPi
			lamEpi = sqrt( pion_mass_c*pion_mass_c + nPiVM_P4.P()*nPiVM_P4.P() );      //--rPiVM_P4.E();
		    }
                    double lamE   = lamEpr + lamEpi;      
                    double lamPx  = nPiVM_P4.Px() + rPiVM_P4.Px();         //--Momentum at Vtx
                    double lamPy  = nPiVM_P4.Py() + rPiVM_P4.Py();
                    double lamPz  = nPiVM_P4.Pz() + rPiVM_P4.Pz();

                    reco::Particle::LorentzVector lambda(0.0,0.0,0.0,0.0);  //--Fill Lorentz lambda P4		  
                    lambda.SetPxPyPzE( lamPx, lamPy, lamPz, lamE );
                    double Mlam = lambda.M();
                    histo_lamb->Fill( Mlam );
		   
                    double MuupK0s  = (uuP4fit + kPiP4 + k0sVM_P4 ).M();  //--Mass: iMu jMu kPi kKa (4 body)(Fitted muons)
                    double MkPikK0s = (          kPiP4 + k0sVM_P4 ).M();
                    
		    if( goodVtx3Body ){                    //--Vertex mumupi

                    //--Check if k0short comes from Dimuon-Pion Vertex (calculate cosine in 3D and 2D)
                    GlobalVector Vector_uuV_to_k0sVM( k0sVM_vx.x()- uupV3_vx.x(),
                                                      k0sVM_vx.y()- uupV3_vx.y(),
                                                      k0sVM_vx.z()- uupV3_vx.z() );			   
                    GlobalVector k0s_P3( k0sVM_P4.px(), k0sVM_P4.py(), k0sVM_P4.pz() ); 
                    double cos_uupk0s3D = ( k0s_P3.dot( Vector_uuV_to_k0sVM ) )/( Vector_uuV_to_k0sVM.mag() * k0s_P3.mag() );

		    if( cos_uupk0s3D > 0.90 ){
		    
                    //--Check separation K0short from Dimuon-Pion (L over Sigma) in 3D
                    typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D_uupk0s_L;
                    typedef ROOT::Math::SVector<double, 3> SVector3_uupk0s_L;
 
                    SMatrixSym3D_uupk0s_L totalCov_uupk0s_L = k0sVM_Cov + uupV3_Cov;   //--3D Separation Dimuon - K0short
                    SVector3_uupk0s_L     uupk0sL_PrimarySep3D( uupV3_vx.x() - k0sVM_vx.x(),
                                                                uupV3_vx.y() - k0sVM_vx.y(),
                                                                uupV3_vx.z() - k0sVM_vx.z() ); 
                    double uuk0sL_L3D    = ROOT::Math::Mag( uupk0sL_PrimarySep3D );
                    double uuk0s_Sigma3D = sqrt(ROOT::Math::Similarity( totalCov_uupk0s_L, uupk0sL_PrimarySep3D )) / uuk0sL_L3D;
		    //--double L/S =
		    
		    //--Do Primary here.  Primary Vertex Re-fitting to exclude the DiMuon, Pion  and Kaon from the primary
                    reco::Vertex refitPrimVertex = *primary.BestVertex();
                    std::vector<reco::TrackRef> ExclusionList;
                    ExclusionList.push_back( iMuTrackRef );
                    ExclusionList.push_back( jMuTrackRef );
                    ExclusionList.push_back( kPiTrk );
		    
// 		    //--HMoreno 03/04/13-----------------------------------------------------------------------
// 		    VertexRefit RefitPrimaryHM(primary.BestVertex() ,ExclusionList ,bFieldHandle ,beamSpot );    //--Highest Multiplicity Primary 
//                     reco::Vertex refitVertexPrimHM = RefitPrimaryHM.Refitted(); 
//                     double primCLHM   = ChiSquaredProbability( refitVertexPrimHM.chi2(), refitVertexPrimHM.ndof());
// 		  
// 		    GlobalPoint  secver = uuRec.RKVertex()->position();
//                     //-GlobalVector secmom = uuRec.RKParent()->currentState().globalMomentum();       //--only dimuon       
//                     GlobalVector secmomTot = GlobalVector(uupRec.P4().Px() + K0sRec.P4().Px()           //--this include dimuon, pion, kaon 
//                                                          ,uupRec.P4().Py() + K0sRec.P4().Py()
//                                                          ,uupRec.P4().Pz())+ K0sRec.P4().Pz();
// 						       
//                     if( MyPrint_ ){
//                     //std::cout<<"B: vector P="<<secmom<<std::endl;
//                     std::cout<<"        : uupk0sRec   ="<<uupkRec.P4().Px() + K0sRec.P4().Px()<<std::endl;
//                     std::cout<<"        : ---->     ="<<secmomTot<<std::endl;
//                     }
// 		  
// 	            //-VertexRefit RefitPrimary(primary.HigherCosAlpha(secver,secmom) ,ExclusionList ,bFieldHandle ,beamSpot ); //--pointing
//                     VertexRefit RefitPrimary(primary.HigherCosAlpha(secver,secmomTot) ,ExclusionList ,bFieldHandle ,beamSpot ); //--pointing
//                     reco::Vertex refitVertexPrim = RefitPrimary.Refitted();
//                     double primCL   = ChiSquaredProbability( refitVertexPrim.chi2(), refitVertexPrim.ndof());
//                  
//                     //--std::cout<<" CL: High Multiplicity ="<<primCLHM<<" Pointing Angle ="<<primCL<<" ->"<<primCLHM-primCL<<std::endl;
//               
//                     //--GlobalPoint PosVertex = GlobalPoint(refitVertexPrim.x(), refitVertexPrim.y(), refitVertexPrim.z() );
// 
// 		    double primDist = sqrt( refitVertexPrim.x()*refitVertexPrim.x()
//                                           + refitVertexPrim.y()*refitVertexPrim.y()
//                                           + refitVertexPrim.z()*refitVertexPrim.z() );
//                     double primR    = sqrt( refitVertexPrim.x()*refitVertexPrim.x()
//                                           + refitVertexPrim.y()*refitVertexPrim.y() );
//  
//                     double primDz   = refitVertexPrim.z() - refitVertexPrimHM.z();
//                     double primDxy  = sqrt( (refitVertexPrim.x() - refitVertexPrimHM.x())*(refitVertexPrim.x() - refitVertexPrimHM.x())
//                                           + (refitVertexPrim.y() - refitVertexPrimHM.y())*(refitVertexPrim.y() - refitVertexPrimHM.y()) );
//                     double primDxyz= sqrt(  (refitVertexPrim.x() - refitVertexPrimHM.x())*(refitVertexPrim.x() - refitVertexPrimHM.x())
//                                           + (refitVertexPrim.y() - refitVertexPrimHM.y())*(refitVertexPrim.y() - refitVertexPrimHM.y())
//                                           + (refitVertexPrim.z() - refitVertexPrimHM.z())*(refitVertexPrim.z() - refitVertexPrimHM.z()) );
//               
//                     if( MyPrint_ ){
//                         std::cout<<"B0: primCL ="<<primCL<<" High Mult: "<<primCLHM<<std::endl;
//                         std::cout<<"B0: primDz ="<<primDz<<" primDxy ="<<primDxy<<" primDxyz ="<<primDxyz<<std::endl;
//                     }
// 
// 		  
// 		  
// 		    //-- end HMoreno 03/04/13-----------------------------------------------------------------

                    VertexRefit RefitPrimary( primary.BestVertex(), ExclusionList, bFieldHandle, beamSpot );
                    refitPrimVertex = RefitPrimary.Refitted();
                    if( refitPrimVertex.isValid() ){                    //--Skip no valid Primary

                    //--Output Prim Vertex fit variables 			       
                    double pV_Chi2 = refitPrimVertex.chi2();
                    double pV_Ndof = refitPrimVertex.ndof();
                    double pV_CL   = ChiSquaredProbability( pV_Chi2, pV_Ndof ); //--CLl Primary
 
                    //--Check separation from primary (L over Sigma)
                    typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
                    typedef ROOT::Math::SVector<double, 3> SVector3;

                    SMatrixSym3D totalCov = refitPrimVertex.covariance() + uupV3_Cov;
                    SVector3 uup_PrimarySep3D( uupV3_vx.x() - refitPrimVertex.x(),
                                               uupV3_vx.y() - refitPrimVertex.y(),
                                               uupV3_vx.z() - refitPrimVertex.z() );

                    double uupL     = ROOT::Math::Mag( uup_PrimarySep3D );  //--Flag 3dim distance Dimuon-Prim
                    double uupSigma = sqrt(ROOT::Math::Similarity( totalCov, uup_PrimarySep3D )) / uupL;
                    //double uupLoS   = uupL / uupSigma;    //--Flag 3dim distance/sigma  Dimuon-Primary  //--L/S

                    //--Here Pointing Secondary Vertex  [Sep 9, 2011]
                    GlobalVector uup_PrimSec( uupV3_vx.x() - refitPrimVertex.x(),
                                              uupV3_vx.y() - refitPrimVertex.y(),
                                              uupV3_vx.z() - refitPrimVertex.z() );     
                    GlobalVector uup_P3( uupV3_P4.Px(), uupV3_P4.Py(), uupV3_P4.Pz() );
                    double V_dot_Puup = uup_P3.dot(uup_PrimSec);                   
                    double cos_alpha_uup_prim = V_dot_Puup /( uup_PrimSec.mag() * uup_P3.mag() );  //--cosine (pointing angle)
               
                    if( cos_alpha_uup_prim>0.9 ){    //--Skip not uup point back
                    histo_dimuPiK[0]->Fill( MuupK0s );     //--B candidate

		    //--fill n-Tuple
                    Comb_kKa++;
                    fill_evt( iEvent, Ntk, allMu, allPrim, allK0s,
		                     Comb_uu, Comb_kPi, Comb_kKa);
                    fill_IuuC();
                    fill_uuC(Muu);

		    fill_kPi(kPi);
		    fill_nPi(nPi);
		    fill_rPi(rPi);

		    fill_k0short(Mk0s, Mk0sVFit, k0sV_CL, k0sVM_CL, Mlam, k0sDist, k0sR,
		                 cos_uupk0s3D, uuk0sL_L3D, uuk0s_Sigma3D,
				 MuupFit, uupV3_CL, pV_CL,
                                 MuupK0s, MkPikK0s, uupL, uupSigma, cos_alpha_uup_prim);
		    
	            hecmu_tree_->Fill();
	            init(1);   //--don't want to clean the counter Comb
                    }   //--End if( cos_alpha_uup_prim>0.9 ){    //--Skip not uup point back
                    }   //--End if( refitPrimVertex.isValid() ){      //--Skip no valid Primary
		    }   //--End if( cos_uupk0s3D > 0.9 ){
		    }   //--End if( goodVtx3Body ){
                    }   //--End if( K0sRec.isValid() ){  Kshort Mass constrain
                    }   //--End if( K0sRec.isValid() ){  Vertex
                    }   //--End if( !sameTrkn && !sameTrkr ){  
		 }  //--End for( unsigned int veeindex = 0; veeindex < theVees.size(); veeindex++ )       
	       }//--End if( NeutralKa_ )
	       
	       }   //--End if( HecClosest(2, kPiState, jState) ){   //--kPi close to jMu close
	       }   //--End if( HecClosest(1, kPiState, iState) ){   //--kPi close to iMu close
	       }   //--End if( good_kPi ){
	       }   //--End if( !sameTrk_kPi_ij ){                                                      //--Skip same Trk
            }      //--End for(unsigned int ktrack = 0; ktrack < theTrackRefs.size(); ktrack++){       //--Pion k 
            }  //--End if( HecDimuVtx(std::vector<RefCountedKinematicParticle> uuParticle) ){	   //--Do here DiMuon Vertex
	}      //--End if( HecClosest(0, iState, jState){    //--iMu close to jMu   
	}      //--End if( uuMcut ){                         //--Skip if not Dimuon Selected	
      }	       //--End if( goodDimuon && nMuVar[0][0]*nMuVar[1][0]<0 ) goodMuon and Opp sign	     
   }	       //--End jMu for(std::vector<pat::Muon>::const_iterator jMu=iMu+1;	  jMu!=allmuons->end();++jMu
}              //--End iMu for(std::vector<pat::Muon>::const_iterator iMu=allmuons->begin();iMu!=allmuons->end();++iMu

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}

// ------------ method called once each job just before starting event loop  ------------
void 
DimuPiK::beginJob()
{
  std::cout << " DimuPiK::beginJob():  Exotic Ntuple" << std::endl;
  hecmu_tree_ = new TTree("Exotic","Dettach Ntuple");
  int bufsize = 64000;

  hecmu_tree_->Branch("Dimuon",&uuC_,"iEta/D:iPhi:iPx:iPy:iPz:id0:idz:iCal:iSeg:iIso:jEta:jPhi:jPx:jPy:jPz:jd0:jdz:jCal:jSeg:jIso:Muu:dca:xcpt:ycpt:zcpt:uuDist:uuR:MuuVF:uuVCL:uuVMpx:uuVMpy:uuVMpz:MuuMFit:uuVMCL:uuVx:uuVy:uuVz:eta:phi", bufsize);
  hecmu_tree_->Branch("kPi",&kPiC_,"kPiQ/D:kPiPx:kPiPy:kPiPz:kPiEn:kPiEta:kPiPhi:kPiNhits:kPid0:kPidz", bufsize);
  hecmu_tree_->Branch("TightMu" , &tightMuC_, "iNormChi2/D:iMuMuHits:iMuStations:iMuPxHits:iMuTrkLayer:jNormChi2:jMuMuHits:jMuStations:jMuPxHits:jMuTrkLayer", bufsize);
  hecmu_tree_->Branch("SoftMu"  , &softMuC_,  "iMuInTrkLayer/D:iMuInNormChi2:jMuInTrkLayer:jMuInNormChi2", bufsize);

if( ChargedKa_ ){
  hecmu_tree_->Branch("kKa",&kKaC_,"kKaQ/D:kKaPx:kKaPy:kKaPz:kKaEn:kKaEta:kKaPhi:kKaNhits:kKad0:kKadz", bufsize);
  hecmu_tree_->Branch("vtx",&vtxC_,"MuukPiFit/D:MuukPi:MuikPi:MujkPi:Ruu_kPi:Muupk:MkPikKa:Ruu_kKa:MuupkFit:uupkV_CL:pV_CL:uupkLoS:uupkL:cos_alpha", bufsize);
}
else {
  hecmu_tree_->Branch("nPi",&nPiC_,"nPiQ/D:nPiPx:nPiPy:nPiPz:nPiEn:nPiEta:nPiPhi:nPiNhits:nPid0:nPidz", bufsize);
  hecmu_tree_->Branch("rPi",&rPiC_,"rPiQ/D:rPiPx:rPiPy:rPiPz:rPiEn:rPiEta:rPiPhi:rPiNhits:rPid0:rPidz", bufsize);
  hecmu_tree_->Branch("k0s",&k0sC_,"Mk0s/D:Mk0sVFit:k0sV_CL:k0sVM_CL:Mlam:k0sDist:k0sR:cos_uupk0s3D:uuk0sL_L3D:uuk0s_Sigma3D:MuupFit:uupV3_CL:pV_CL:MuupK0s:MkPikK0s:uupL:uupSigma:cos_uup_prim", bufsize);
}

  //--Integer Branch
  hecmu_tree_->Branch("iHeader",&IevtC_,"run/I:evtnum:lumib:allTrk:allMu:allPrim:allK0s:upTrig:umTrig:upLMtrig:umLMtrig:evtTrig1:evtTrig2:evtTrig3:evtTrig4:Comb_uu:Comb_kPi:Comb_kKa", bufsize);
  hecmu_tree_->Branch("iDimuon",&IuuC_, "iQ/I:imuId:iNtrkHits:jQ:jmuId:jNtrkHits", bufsize);

}

// ------------ method called once each job just after ending the event loop  ------------
void 
DimuPiK::endJob() 
{
std::cout << " DimuPiK::endJob()" << std::endl;
}

// ------------ method called when starting to processes a run  ------------
void 
DimuPiK::beginRun(edm::Run const&, edm::EventSetup const&)
{
std::cout << " DimuPiK::beginRun()" << std::endl;
}

// ------------ method called when ending the processing of a run  ------------
void 
DimuPiK::endRun(edm::Run const&, edm::EventSetup const&)
{
std::cout << " DimuPiK::endRun()" << std::endl;
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
DimuPiK::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
std::cout << " DimuPiK::beginLuminosityBlock()" << std::endl;
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
DimuPiK::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
std::cout << " DimuPiK::endLuminosityBloc()" << std::endl;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DimuPiK::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

 //Specify that only 'tracks' is allowed
 //To use, remove the default given above and uncomment below
 //ParameterSetDescription desc;
 //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
 //descriptions.addDefault(desc);
}
//----------Functions for my N-tuple
void DimuPiK::init(int all)
{
  if( all == 0 )IevtC_.init();    //--clean everything
  IuuC_.init();
   uuC_.init();
   kPiC_.init();
   kKaC_.init();
   vtxC_.init();
   k0sC_.init();
}
void DimuPiK::Ievt::init()
{
  runNb      = -1; eventNb    = -1; lumiBlock  = -1;
  allTrk     = -1; allMu      = -1; allPrim    = -1; allK0s   = -1;
  upTrig     = -1; umTrig     = -1; upLMtrig   = -1; umLMtrig = -1;
  evtTrig1   =  0; evtTrig2   =  0; evtTrig3   =  0; evtTrig4 =  0;
  Comb_uu    =  0; Comb_kPi   =  0; Comb_kKa   =  0;
}
void DimuPiK::fill_evt(const edm::Event& iEvent, int Ntk, int allMu, int allPrim, int allK0s,
                       int Comb_uu, int Comb_kPi, int Comb_kKa)
{
  IevtC_.runNb = iEvent.id().run();  IevtC_.eventNb = iEvent.id().event(); IevtC_.lumiBlock = iEvent.luminosityBlock();
  IevtC_.allTrk     = Ntk;            //--tracks->size();
  IevtC_.allMu      = allMu;          //--->allmuons->size();
  IevtC_.allPrim    = allPrim;
  IevtC_.allK0s     = allK0s;
  IevtC_.upTrig     = nTrig[0][0];    //--upTrig 
  IevtC_.umTrig     = nTrig[1][0];    //--umTrig  
  IevtC_.upLMtrig   = nTrig[0][1];    //--Low Mass trig fro muon Pos 
  IevtC_.umLMtrig   = nTrig[1][1];    //--Low Mass trig fro muon Neg
  IevtC_.evtTrig1   = evtTrig1; IevtC_.evtTrig2   = evtTrig2; IevtC_.evtTrig3   = evtTrig3; IevtC_.evtTrig4   = evtTrig4;
  IevtC_.Comb_uu    = Comb_uu;  IevtC_.Comb_kPi   = Comb_kPi; IevtC_.Comb_kKa   = Comb_kKa;
}
void DimuPiK::IuuC::init()
{
  iQ = 0; imuId = -1; iNtrkHits = -1;
  jQ = 0; jmuId = -1; jNtrkHits = -1;
}
void DimuPiK::fill_IuuC()
{
  IuuC_.iQ = nMuVar[0][0]; IuuC_.imuId = nMuVar[0][1]; IuuC_.iNtrkHits = nMuVar[0][4];
  IuuC_.jQ = nMuVar[1][0]; IuuC_.jmuId = nMuVar[1][1]; IuuC_.jNtrkHits = nMuVar[1][4];
}
void DimuPiK::uuC::init()
{
  iEta =9.9; iPhi =-1.0; iPx =-999; iPy =-999; iPz=-999; id0 =-1.0; idz =-1.0; iCal =-1.0; iSeg =-1.0; iIso =-1.0;
  jEta =9.9; jPhi =-1.0; jPx =-999; jPy =-999; jPz=-999; jd0 =-1.0; jdz =-1.0; jCal =-1.0; jSeg =-1.0; jIso =-1.0;
  Muu = -1.0; dca = -1.0; xcpt = -1.0; ycpt = -1.0; zcpt = -1.0; uuDist = -1.0; uuR = -1.0; 
  MuuVF=-1.0; uuVCL=-1.0; uuVMpx=-1;uuVMpy=-1;uuVMpz=-1;
  MuuMFit=-1; uuVMCL=-1; Vx = 999; Vy = 999; Vz = 999; eta = 9; phi = 9;
}
void DimuPiK::fill_uuC(double Muu)
{
  uuC_.iEta    = MuVar[0][0]; 
  uuC_.iPhi    = MuVar[0][1]; 
  uuC_.iPx     = MuVar[0][4]; 
  uuC_.iPy     = MuVar[0][5];
  uuC_.iPz     = MuVar[0][6]; 
  uuC_.id0     = MuVar[0][11]; 
  uuC_.idz     = MuVar[0][12]; 
  uuC_.iCal    = MuVar[0][7]; 
  uuC_.iSeg    = MuVar[0][8];
  uuC_.iIso    = MuVar[0][9];
  
  uuC_.jEta    = MuVar[1][0]; 
  uuC_.jPhi    = MuVar[1][1]; 
  uuC_.jPx     = MuVar[1][4]; 
  uuC_.jPy     = MuVar[1][5];
  uuC_.jPz     = MuVar[1][6]; 
  uuC_.jd0     = MuVar[1][11];
  uuC_.jdz     = MuVar[1][12];
  uuC_.jCal    = MuVar[1][7]; 
  uuC_.jSeg    = MuVar[1][8];
  uuC_.jIso    = MuVar[1][9]; 
  
  uuC_.Muu     = Muu;
  uuC_.dca     = dcaVar[0][0];
  uuC_.xcpt    = dcaVar[0][1];
  uuC_.ycpt    = dcaVar[0][2];
  uuC_.zcpt    = dcaVar[0][3];
  uuC_.uuDist  = uuVtx[7];
  uuC_.uuR     = uuVtx[8];
  uuC_.MuuVF   = uuVtx[0];     //--MuuFit
  uuC_.uuVCL   = uuVtx[1];     //--uuV_CL
  uuC_.uuVMpx  = uuVtx[13];    //--uuVMpx;
  uuC_.uuVMpy  = uuVtx[14];    //--uuVMpy;
  uuC_.uuVMpz  = uuVtx[15];    //--uuVMpz;
  uuC_.MuuMFit = uuVtx[2];     //--MuuMFit
  uuC_.uuVMCL  = uuVtx[3];     //--uuVM_CL
  uuC_.Vx      = uuVtx[4];
  uuC_.Vy      = uuVtx[5];
  uuC_.Vz      = uuVtx[6];
  uuC_.eta     = uuVtx[24];
  uuC_.phi     = uuVtx[25];
}
void DimuPiK::kPi::init()
{
kPiQ=0; kPiPx=-999; kPiPy=-999; kPiPz=-999; kPiEn=-999; kPiEta=9; kPiPhi=9; kPiNhits=0; kPid0=99; kPidz=99;
}
void DimuPiK::fill_kPi(int n)
{
  kPiC_.kPiQ     = TrkVar[n][0];  
  kPiC_.kPiPx    = TrkVar[n][1];
  kPiC_.kPiPy    = TrkVar[n][2];
  kPiC_.kPiPz    = TrkVar[n][3];   
  kPiC_.kPiEn    = TrkVar[n][4];
  kPiC_.kPiEta   = TrkVar[n][7];
  kPiC_.kPiPhi   = TrkVar[n][8];
  kPiC_.kPiNhits = TrkVar[n][9];
  kPiC_.kPid0    = TrkVar[n][10];
  kPiC_.kPidz    = TrkVar[n][11];
}
void DimuPiK::nPi::init()
{
nPiQ=0; nPiPx=-999; nPiPy=-999; nPiPz=-999; nPiEn=-999; nPiEta=9; nPiPhi=9; nPiNhits=0; nPid0=99; nPidz=99;
}
void DimuPiK::fill_nPi(int n)
{
  nPiC_.nPiQ     = TrkVar[n][0];  
  nPiC_.nPiPx    = TrkVar[n][1];
  nPiC_.nPiPy    = TrkVar[n][2];
  nPiC_.nPiPz    = TrkVar[n][3];   
  nPiC_.nPiEn    = TrkVar[n][4];
  nPiC_.nPiEta   = TrkVar[n][7];
  nPiC_.nPiPhi   = TrkVar[n][8];
  nPiC_.nPiNhits = TrkVar[n][9];
  nPiC_.nPid0    = TrkVar[n][10];
  nPiC_.nPidz    = TrkVar[n][11];
}
void DimuPiK::rPi::init()
{
rPiQ=0; rPiPx=-999; rPiPy=-999; rPiPz=-999; rPiEn=-999; rPiEta=9; rPiPhi=9; rPiNhits=0; rPid0=99; rPidz=99;
}
void DimuPiK::fill_rPi(int n)
{
  rPiC_.rPiQ     = TrkVar[n][0];  
  rPiC_.rPiPx    = TrkVar[n][1];
  rPiC_.rPiPy    = TrkVar[n][2];
  rPiC_.rPiPz    = TrkVar[n][3];   
  rPiC_.rPiEn    = TrkVar[n][4];
  rPiC_.rPiEta   = TrkVar[n][7];
  rPiC_.rPiPhi   = TrkVar[n][8];
  rPiC_.rPiNhits = TrkVar[n][9];
  rPiC_.rPid0    = TrkVar[n][10];
  rPiC_.rPidz    = TrkVar[n][11];
}
void DimuPiK::kKa::init()
{
kKaQ=0; kKaPx=-999; kKaPy=-999; kKaPz=-999; kKaEn=-999; kKaEta=9; kKaPhi=9; kKaNhits=0; kKad0=99; kKadz=99;
}
void DimuPiK::fill_kKa(int n)
{
  kKaC_.kKaQ     = TrkVar[n][0];  
  kKaC_.kKaPx    = TrkVar[n][1];
  kKaC_.kKaPy    = TrkVar[n][2];
  kKaC_.kKaPz    = TrkVar[n][3];   
  kKaC_.kKaEn    = TrkVar[n][4];
  kKaC_.kKaEta   = TrkVar[n][7];
  kKaC_.kKaPhi   = TrkVar[n][8];
  kKaC_.kKaNhits = TrkVar[n][9];
  kKaC_.kKad0    = TrkVar[n][10];
  kKaC_.kKadz    = TrkVar[n][11];
}
void DimuPiK::vtx::init()
{
MuukPiFit=-1; MuukPi=-1;  MuikPi=-1;   MujkPi=-1;   Ruu_kPi=-1;
Muupk=-1;     MkPikKa=-1;  Ruu_kKa=-1; MuupkFit=-1; uupkV_CL=-1;
pV_CL=-1;     uupkLoS=-1; uupkL=-1;    cos_alpha=-1;
}
void DimuPiK::fill_vtxAll(double MuukPiFit, double MuukPi,  double MuikPi,  double MujkPi,   double Ruu_kPi,
                          double Muupk,     double MkPikKa, double Ruu_kKa, double MuupkFit, double uupkV_CL, 
                          double pV_CL,     double uupkLoS, double uupkL,    double cos_alpha)
{
  vtxC_.MuukPiFit  = MuukPiFit;
  vtxC_.MuukPi    = MuukPi;
  vtxC_.MuikPi    = MuikPi;
  vtxC_.MujkPi    = MujkPi;
  vtxC_.Ruu_kPi   = Ruu_kPi;
  vtxC_.Muupk     = Muupk; 
  vtxC_.MkPikKa   = MkPikKa;
  vtxC_.Ruu_kKa   = Ruu_kKa;
  vtxC_.MuupkFit  = MuupkFit; 
  vtxC_.uupkV_CL  = uupkV_CL;
  vtxC_.pV_CL     = pV_CL;
  vtxC_.uupkLoS   = uupkLoS;
  vtxC_.uupkL     = uupkL;
  vtxC_.cos_alpha = cos_alpha;
}
void DimuPiK::k0s::init(){
  Mk0s = -1; Mk0sVFit = -1; k0sV_CL = -1; k0sVM_CL = -1; Mlam = -1;k0sDist=-1; k0sR=-1;
  cos_uupk0s3D=-1; uuk0sL_L3D=-1; uuk0s_Sigma3D=-1;
  MuupFit=-1; uupV3_CL=0; pV_CL=0;
  MuupK0s=-1; MkPikK0s=-1; uupL=-1; uupSigma=-1; cos_alpha_uup_prim=-1;
}
void DimuPiK::fill_k0short(double Mk0s, double Mk0sVFit, double k0sV_CL, double k0sVM_CL, double Mlam, double k0sDist, double k0sR,
                           double cos_uupk0s3D, double uuk0sL_L3D, double uuk0s_Sigma3D,
			   double MuupFit, double uupV3_CL, double pV_CL,
                           double MuupK0s, double MkPikK0s, double uupL, double uupSigma, double cos_alpha_uup_prim){
  k0sC_.Mk0s     = Mk0s;
  k0sC_.Mk0sVFit = Mk0sVFit;
  k0sC_.k0sV_CL  = k0sV_CL;
  k0sC_.k0sVM_CL = k0sVM_CL;
  k0sC_.Mlam     = Mlam;
  k0sC_.k0sDist  =k0sDist;
  k0sC_.k0sR     = k0sR;
  k0sC_.cos_uupk0s3D  = cos_uupk0s3D;
  k0sC_.uuk0sL_L3D    = uuk0sL_L3D;
  k0sC_.uuk0s_Sigma3D = uuk0s_Sigma3D;
  k0sC_.MuupFit       = MuupFit;
  k0sC_.uupV3_CL      = uupV3_CL;
  k0sC_.pV_CL         = pV_CL;
  k0sC_.MuupK0s       = MuupK0s;
  k0sC_.MkPikK0s      = MkPikK0s;
  k0sC_.uupL          = uupL;
  k0sC_.uupSigma      = uupSigma;
  k0sC_.cos_alpha_uup_prim = cos_alpha_uup_prim;
}

//////Tight and soft muons functions///////////////////////

void DimuPiK::tightMuC::init()
{
  iNormChi2 = 0, iMuMuHits = -1, iMuStations = -1, iMuPxHits = -1, iMuTrkLayer = -1;
  jNormChi2 = 0, jMuMuHits = -1, jMuStations = -1, jMuPxHits = -1, jMuTrkLayer = -1;
}            
void DimuPiK::fill_tightMuons( )
{
 tightMuC_.iNormChi2     = CutMu[0][0];
 tightMuC_.iMuMuHits     = CutMu[0][1];
 tightMuC_.iMuStations   = CutMu[0][2];
 tightMuC_.iMuPxHits     = CutMu[0][3];
 tightMuC_.iMuTrkLayer   = CutMu[0][4];
 tightMuC_.jNormChi2     = CutMu[1][0];
 tightMuC_.jMuMuHits     = CutMu[1][1];
 tightMuC_.jMuStations   = CutMu[1][2];
 tightMuC_.jMuPxHits     = CutMu[1][3];
 tightMuC_.jMuTrkLayer   = CutMu[1][4];
}
void DimuPiK::softMuC::init()
{
iMuInTrkLayer = -1, iMuInNormChi2 = 0;
jMuInTrkLayer = -1, jMuInNormChi2 = 0;
}            
void DimuPiK::fill_softMuons( )
{
 softMuC_.iMuInNormChi2   = CutMu[0][5];
 softMuC_.iMuInTrkLayer   = CutMu[0][6];
 softMuC_.jMuInNormChi2   = CutMu[1][5];
 softMuC_.jMuInTrkLayer   = CutMu[1][6];
}
void DimuPiK::tightKMuC::init()
{
 kNormChi2 = 0, kMuMuHits = -1, kMuStations = -1, kMuPxHits = -1, kMuTrkLayer = -1;
}            
void DimuPiK::fill_tightKMuons( )
{
 tightKMuC_.kNormChi2	= CutMu[2][0];
 tightKMuC_.kMuMuHits	= CutMu[2][1];
 tightKMuC_.kMuStations	= CutMu[2][2];
 tightKMuC_.kMuPxHits	= CutMu[2][3];
 tightKMuC_.kMuTrkLayer	= CutMu[2][4];
}
void DimuPiK::softKMuC::init()
{
kMuInTrkLayer = -1, kMuInNormChi2 = 0;
}            
void DimuPiK::fill_softKMuons( )
{
 softKMuC_.kMuInNormChi2   = CutMu[2][5];
 softKMuC_.kMuInTrkLayer   = CutMu[2][6];
}
void DimuPiK::HecCutMu(int m, std::vector<pat::Muon>::const_iterator mMu){
 //--Do here the Tight nMuon Variables  Feb 2013
   CutMu[m][0] = 99;
   CutMu[m][1] = -1;
      if( mMu->isGlobalMuon() ){
         CutMu[m][0] = mMu->globalTrack()->normalizedChi2();
         CutMu[m][1] = mMu->globalTrack()->hitPattern().numberOfValidMuonHits();
      }
   CutMu[m][2] = mMu->numberOfMatchedStations();
   CutMu[m][3] = mMu->innerTrack()->hitPattern().numberOfValidPixelHits();
   CutMu[m][4] = mMu->track()->hitPattern().trackerLayersWithMeasurement();
 
  //--Do here the Soft nMuon Variables Feb 2013
   CutMu[m][5] = 99;
   CutMu[m][6] = -1;
      if( mMu->isTrackerMuon() ){
         CutMu[m][5] = mMu->innerTrack()->normalizedChi2();
         CutMu[m][6] = mMu->innerTrack()->hitPattern().pixelLayersWithMeasurement();
      } 
}
//////end Tight and soft muons functions///////////////////




//--Functions for my Analyzer
void DimuPiK::HecMuVar(int n, std::vector<pat::Muon>::const_iterator nMu, reco::TrackRef nMuTrk){
  //--Double Variables
  MuVar[n][0]  = nMu->eta();
  MuVar[n][1]  = nMu->phi();
  MuVar[n][2]  = nMu->p4().Pt();
  MuVar[n][3]  = nMu->p4().P();
  MuVar[n][4]  = nMu->p4().Px();
  MuVar[n][5]  = nMu->p4().Py();
  MuVar[n][6]  = nMu->p4().Pz();
  MuVar[n][7]  = nMu->caloCompatibility();
  MuVar[n][8]  = muon::segmentCompatibility(*nMu);
  MuVar[n][9]  = nMu->isolationR03().sumPt;
  MuVar[n][10] = nMu->combinedQuality().trkKink;   		//--Chi squared
  MuVar[n][11] = 999;             //--for StandAloneMuon  (for Tracker and Global Muon after the MuonId)
  MuVar[n][12] = 999;	          //--for StandAloneMuon  (for Tracker and Global Muon after the MuonId)
  
  //--Integer Variables
  nMuVar[n][0]  = nMu->charge();
 
  int muId=0, muSA=0, muGL=0, muTK=0; 	      //--Muon Id
  if( nMu->isStandAloneMuon() )muSA = 1;      //--2^0
  if( nMu->isGlobalMuon()     )muGL = 10;     //--2^1
  if( nMu->isTrackerMuon()    )muTK = 100;    //--2^2
  int imuId = muTK + muGL + muSA;                                   //--bynary number
  muId = muTK*pow(2,2)/100 + muGL*pow(2,1)/10 + muSA*pow(2,0)/1;    //--convert to decimal
     
  int iGL = 0;
  if(  nMu->isGlobalMuon()                          )iGL = 1;
  if(!(nMu->isGlobalMuon()) && nMu->isTrackerMuon() )iGL = 2;
  
  nMuVar[n][1]  = imuId;
  nMuVar[n][2]  = muId;       //1:--1  : S
                              //2:--10 : G
                              //3:--11 : G-S
                              //4:--100: T
                              //5:--101: T-S
                              //6:--110: T-G  doesn't exist
                              //7:--111: T-G-S
  nMuVar[n][3]  = iGL;  
  nMuVar[n][4]  = 0;
  
  if( nMuVar[n][1] > 1 ){      //--for Tracker & Global Muons
     MuVar[n][11] = nMuTrk->d0();
     MuVar[n][12] = nMuTrk->dz();
    nMuVar[n][4]  = nMuTrk->numberOfValidHits();
  }
  if( MyPrint_ ){
    std::cout<<"Hec: MuId: iTK:"<<muTK<<" iGL:"<<muGL<<" SA:"<<muSA<<" Dec:"<<muId<<" Bin:"<<imuId<<std::endl;
    std::cout<<"Hec:  Mu Compatibility:"<<MuVar[n][7]<<" MuSeg:"<<MuVar[n][8]<<" MuIso:"<<MuVar[n][9]<<std::endl;
  }
    
}

void DimuPiK::calcMyP(int n,const reco::Candidate *Cand){
   
      myP[n][0] = Cand->energy();
      myP[n][1] = Cand->px();
      myP[n][2] = Cand->py();
      myP[n][3] = Cand->pz();
      myP[n][4] = sqrt(myP[n][0]*myP[n][0] - myP[n][1]*myP[n][1]
   		  - myP[n][2]*myP[n][2] - myP[n][3]*myP[n][3]);
      myP[n][5] = Cand->charge();
      myP[n][6] = sqrt(myP[n][1]*myP[n][1] + myP[n][2]*myP[n][2]);
      myP[n][7] = Cand->eta();
      myP[n][8] = Cand->mass();
      myP[n][9] = Cand->pdgId();
      myP[n][10]= myP[n][8] - myP[n][4];	 

   if( MyPrintMC_){
       std::cout<<" Energy = "     <<myP[n][0]<<" , px = "<<myP[n][1]
   		<<" , py = "       <<myP[n][2]<<" , pz = "<<myP[n][3]<<", Inv. mass = "<<myP[n][4]
     		<<", Charge = "    <<myP[n][5]<<", pT = " <<myP[n][6]<<", eta = "      <<myP[n][7]
     		<<", Gen. mass"    <<myP[n][8]<<", Id :"  <<myP[n][9]
     		<<", Mass Diff = " <<myP[n][10]<<std::endl;				  
       std::cout<<"   "<<std::endl;
   }

}


bool DimuPiK::HecClosest(int Mode, FreeTrajectoryState iState, FreeTrajectoryState jState){
//--Meaning of  dcaVar[Mode][4]:
//--Mode=0:Mui-Muj     0-4:dca, xCpt,yCpt,zCpt
//--     1:kPi-Mui 
//--     2:kPi-Muj 
//--     3:kKa-Mui
//--     4:kKa-Muj
//--     5:kKa-kPi 
//--Measure distance between tracks at their closest approach
    ClosestApproachInRPhi cApp; 
    GlobalPoint cxPt;	   
    double dca= -9999;
     
    cApp.calculate( iState, jState);
    if( !cApp.status() ){
      std::cout<<"dca bad status capp"<<std::endl;
      //--continue;    Need to check where this continue will go ??
    } 
    else{
          dca = fabs( cApp.distance() );
          cxPt = cApp.crossingPoint();
    }
    dcaVar[Mode][0] = dca;	    //--Flag  Closest Approach distance between muons
    dcaVar[Mode][1] = cxPt.x();     //--Crosing point in X
    dcaVar[Mode][2] = cxPt.y();     //--Crosing point in Y
    dcaVar[Mode][3] = cxPt.z();     //--Crosing point in Z
    
    //--Cut on fiducial tracking volume [-120 cm < r < 120 cm  and -260 cm < Z < 260 cm]   r = sqrt( x^2 + y^2 )
    //--Pixel [-20 cm < r < 20 cm  and -60 cm < Z < 60 cm] to run vertex fit
    if( ( dca >= 0. && dca <= 1 ) && 
        ( sqrt( cxPt.x()*cxPt.x() + cxPt.y()*cxPt.y() ) < 120.) &&
          std::abs(cxPt.z()) < 300. ){
        return true;
    }
    return false;
}



void DimuPiK::HecHltMuTrig(int nMu, const pat::Muon* patMu){      
//--unpacks J/Psi trigger bit
      int upTrig2mu3=!patMu->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered3").empty();

      int upTrig2mu3JPsi=!patMu->triggerObjectMatchesByFilter("hltDoubleMu3JpsiL3Filtered").empty();

      int upTrig2mu6p5JPsiDisp=!patMu->triggerObjectMatchesByFilter("hltDimuon6p5JpsiDisplacedL3Filtered").empty();
    		  
      int upTrig2mu7JPsiDisp=!patMu->triggerObjectMatchesByFilter("hltDisplacedmumuFilterJpsi").empty(); //--also the same for HLT_DoubleMu3p5_Jpsi_Displaced_v2

      int upTrig2mu4JPsiDisp=!patMu->triggerObjectMatchesByFilter("hltDisplacedmumuFilterDoubleMu4Jpsi").empty();	       

      int upTrig2mu6p5JPsiPro=!patMu->triggerObjectMatchesByFilter("hltDimuon6p5JpsiL3Filtered").empty();

      int upTrig2mu0JPsiPro=!patMu->triggerObjectMatchesByFilter("hltVertexmumuFilterJpsi").empty();

      int upTrig2mu0JPsiNoVtxPro=!patMu->triggerObjectMatchesByFilter("hltJpsiNoVertexingL3Filtered").empty();

      nTrig[nMu][0] = upTrig2mu3             *1
    	 	    + upTrig2mu3JPsi         *10
    		    + upTrig2mu6p5JPsiDisp   *100
    		    + upTrig2mu7JPsiDisp     *1000
    		    + upTrig2mu4JPsiDisp     *10000
    		    + upTrig2mu6p5JPsiPro    *100000
    		    + upTrig2mu0JPsiPro      *1000000
    		    + upTrig2mu0JPsiNoVtxPro *10000000;
    		
      //--Unpacks Displaced LowMass_Dimuon trigger bit
      int upLMDimuon6p5=!patMu->triggerObjectMatchesByPath("HLT_Dimuon6p5_LowMass_Displaced_v*").empty();

      int upLMDimuon7  =!patMu->triggerObjectMatchesByPath("HLT_Dimuon7_LowMass_Displaced_v*").empty();

      int upLMDoubleMu4=!patMu->triggerObjectMatchesByPath("HLT_DoubleMu4_LowMass_Displaced_v*").empty();

      int upLMDoubleMu4p5=!patMu->triggerObjectMatchesByPath("HLT_DoubleMu4p5_LowMass_Displaced_v*").empty();
     
      int upLMDoubleMu5=!patMu->triggerObjectMatchesByPath("HLT_DoubleMu5_LowMass_Displaced_v*").empty();

      nTrig[nMu][1] = upLMDimuon6p5   *1
    		    + upLMDimuon7     *10
    		    + upLMDoubleMu4   *100
    		    + upLMDoubleMu4p5 *1000
    		    + upLMDoubleMu5   *10000;

      if( MyPrint_ ){
    	  std::cout << "Mu+ Trigger : "<<nTrig[nMu][0]<<" nMu "<<std::endl;
    	  std::cout << "Mu Trigger mu3 "         <<upTrig2mu3	         << std::endl;
    	  std::cout << "Mu TriggerJPsi "         <<upTrig2mu3JPsi        << std::endl;
    	  std::cout << "Mu TriggerJPsiDisp "     <<upTrig2mu6p5JPsiDisp  << std::endl;
    	  std::cout << "Mu Trigger7JPsiDisp "    <<upTrig2mu7JPsiDisp    << std::endl;
    	  std::cout << "Mu Trigger4JPsiDisp "    <<upTrig2mu4JPsiDisp    << std::endl;
    	  std::cout << "Mu Trigger5JPsiDisp "    <<upTrig2mu6p5JPsiPro   << std::endl;
    	  std::cout << "Mu Trigger0JPsiDisp "    <<upTrig2mu0JPsiPro     << std::endl;
    	  std::cout << "Mu TriggerNoVtxJPsiDisp "<<upTrig2mu0JPsiNoVtxPro<< std::endl;
	  std::cout << " "<<std::endl;
          std::cout << " LowMass Triger : "<< nTrig[nMu][1]<<std::endl;
          std::cout << "Dimuon6p5         "<< upLMDimuon6p5   <<std::endl;
          std::cout << "Dimuon7           "<< upLMDimuon7     <<std::endl;
          std::cout << "DoubleMu4         "<< upLMDoubleMu4   <<std::endl;
          std::cout << "DoubleMu4p5       "<< upLMDoubleMu4p5 <<std::endl;
          std::cout << "DoubleMu5         "<< upLMDoubleMu5   <<std::endl;	   
          std::cout << "-----------------------------"<<std::endl;
      }
}
void DimuPiK::HecHltTrig(const edm::Event& iEvent){
//--Trigger
  evtTrig1=0, evtTrig2=0, evtTrig3=0, evtTrig4=0;
 
//--Do trigger here (code from Keith)	
  unsigned int hlt_mu3=0, hlt_mu5=0, hlt_mu8=0;
  unsigned int hlt_2mu0=0,hlt_2mu3=0,hlt_2mu3JPsi_v1=0,hlt_2mu3JPsi_v2=0,hlt_2mu3_quark_v1=0,hlt_2mu3_quark_v2=0; 
  unsigned int hlt_2mu6p5_dis=0,hlt_2mu7_dis=0,hlt_2mu3p5_dis=0,hlt_2mu4_dis=0;
  unsigned int hlt_2mu6p5_pro=0, hlt_2mu0_pro_v1=0, hlt_2mu0_pro_v5=0, hlt_2mu0_pro_v6=0;
  unsigned int hlt_2mu0_pro_NoVtx_v2=0, hlt_2mu0_pro_NoVtx_v3=0;
  unsigned int hlt_mu0trk0=0, hlt_mu3trk0=0, hlt_mu0trkmu0=0, hlt_mu3trkmu0=0, hlt_mu0trkmu0OST=0, hlt_mu3trkmu0OST=0;
  unsigned int hlt_mu0trkmu0OST_tight=0;
  unsigned int hlt_L1muOpen=0, hlt_L12muOpen=0, hlt_L12muOpenTight=0, hlt_2mu0L2=0, hlt_DimuPsi2s=0, hlt_DimuLowMd=0;
    
//--first get HLT results
  edm::Handle<edm::TriggerResults> hltresults;
  try {
   std::string const &trig = std::string("TriggerResults::")+hlTriggerResults_;
   iEvent.getByLabel(edm::InputTag(trig),hltresults);
  }
  catch ( ... ) {
   std::cout << "Couldn't get handle on HLT Trigger!" << std::endl;
  }
  if (!hltresults.isValid()) {
   std::cout << "No Trigger Results!" << std::endl;
  }
  else {
   int ntrigs = hltresults->size();
   if (ntrigs == 0){
     std::cout << "No trigger name given in TriggerResults of the input " << std::endl;
  }	  
//--get hold of trigger names - based on TriggerResults object!
  const edm::TriggerNames & triggerNames_ = iEvent.triggerNames(*hltresults);
  
  //HLTConfigProvider hltConfig_;
  for(int itrig=0; itrig< ntrigs; itrig++){
   TString trigName = triggerNames_.triggerName(itrig);
   int hltflag = (*hltresults)[itrig].accept();     
   int trigPrescale = 0;
   //if( hltflag == 1 ){
       //std::string trigName_cp = triggerNames_.triggerName(itrig);
       //trigPrescale = hltConfig_.prescaleValue(itrig, trigName_cp);
   //}
       
   if( MyPrint_ ) 
       std::cout << "Trigger " <<  trigName << " was passed = " <<  hltflag <<" prs_Scale:"<<trigPrescale<<std::endl;
       
   if(trigName=="HLT_Mu3_v3") hlt_mu3 = hltflag;
   if(trigName=="HLT_Mu5_v5") hlt_mu5 = hltflag;
   if(trigName=="HLT_Mu8_v1") hlt_mu8 = hltflag;
   if(trigName=="HLT_Dimuon0_Jpsi_v1"	     ) hlt_2mu0 = hltflag;
   if(trigName=="HLT_DoubleMu3_v2"	     ) hlt_2mu3 = hltflag;
   if(trigName=="HLT_DoubleMu3_Jpsi_v1"      ) hlt_2mu3JPsi_v1   = hltflag;   // 160329-161176 data
   if(trigName=="HLT_DoubleMu3_Jpsi_v2"      ) hlt_2mu3JPsi_v2   = hltflag;   //  161216-163261 data   also in Brian's MC
   if(trigName=="HLT_DoubleMu3_Quarkonium_v1") hlt_2mu3_quark_v1 = hltflag;   // 160329-161176 data 
   if(trigName=="HLT_DoubleMu3_Quarkonium_v2") hlt_2mu3_quark_v2 = hltflag;   //161216-163261 data    also in Brian's MC
  
   if(trigName=="HLT_Dimuon6p5_Jpsi_Displaced_v1") hlt_2mu6p5_dis = hltflag;  // 163269-163869 data
   if(trigName=="HLT_Dimuon7_Jpsi_Displaced_v1"||
     trigName=="HLT_Dimuon7_Jpsi_Displaced_v3"  ) hlt_2mu7_dis   = hltflag;  // 165088-170248 data
 
   if(trigName=="HLT_DoubleMu3p5_Jpsi_Displaced_v2") hlt_2mu3p5_dis = hltflag; // 170249-173235 data
   if(trigName=="HLT_DoubleMu4_Jpsi_Displaced_v1"  ) hlt_2mu4_dis   = hltflag; // 173236- data
  
   if(trigName=="HLT_Dimuon6p5_Jpsi_v1") hlt_2mu6p5_pro  = hltflag;  // 163269-163869 data  prescaled
   if(trigName=="HLT_Dimuon0_Jpsi_v1"  ) hlt_2mu0_pro_v1 = hltflag;  // 165121-170248 data  prescaled
   if(trigName=="HLT_Dimuon0_Jpsi_v5"  ) hlt_2mu0_pro_v5 = hltflag;  // 170249-173235 data  prescaled
   if(trigName=="HLT_Dimuon0_Jpsi_v6"  ) hlt_2mu0_pro_v6 = hltflag;  // 173263- data  prescaled
 
   if(trigName=="HLT_Dimuon0_Jpsi_NoVertex_v2") hlt_2mu0_pro_NoVtx_v2 = hltflag;  // 170249-173235 data  prescaled
   if(trigName=="HLT_Dimuon0_Jpsi_NoVertex_v3") hlt_2mu0_pro_NoVtx_v3 = hltflag;  // 173263- data  prescaled
   
   if(trigName=="HLT_Mu0_Track0_Jpsi"	) hlt_mu0trk0	   = hltflag;
   if(trigName=="HLT_Mu3_Track0_Jpsi"	) hlt_mu3trk0	   = hltflag;
   if(trigName=="HLT_Mu0_TkMu0_Jpsi"	) hlt_mu0trkmu0    = hltflag;
   if(trigName=="HLT_Mu3_TkMu0_Jpsi"	) hlt_mu3trkmu0    = hltflag;
   if(trigName=="HLT_Mu0_TkMu0_OST_Jpsi") hlt_mu0trkmu0OST = hltflag;
   if(trigName=="HLT_Mu3_TkMu0_OST_Jpsi") hlt_mu3trkmu0OST = hltflag;
   
   if(trigName=="HLT_Mu0_TkMu0_OST_Jpsi_Tight_v1") hlt_mu0trkmu0OST_tight = hltflag;
   
   if(trigName=="HLT_L1MuOpen"  	  ) hlt_L1muOpen       = hltflag;
   if(trigName=="HLT_L1DoubleMuOpen"	  ) hlt_L12muOpen      = hltflag;
   if(trigName=="HLT_L1DoubleMuOpen_Tight") hlt_L12muOpenTight = hltflag;
   if(trigName=="HLT_L2DoubleMu0_v4"	  ) hlt_2mu0L2         = hltflag;
   if(trigName=="HLT_Dimuon7_PsiPrime_v3"	  )hlt_DimuPsi2s  = hltflag;
   if(trigName=="HLT_Dimuon7_LowMass_Displaced_v4")hlt_DimuLowMd  = hltflag;
  } //--for (int itrig=0; itrig< ntrigs; itrig++) {   
  }  //--end else for if( !hltresults.isValid() ){
     
  //--pack event Trigger
  evtTrig1 = hlt_mu3           *1
           + hlt_mu5           *10
           + hlt_mu8           *100
           + hlt_2mu0          *1000
           + hlt_2mu3          *10000
           + hlt_2mu3JPsi_v1   *100000
           + hlt_2mu3JPsi_v2   *1000000
           + hlt_2mu3_quark_v1 *10000000;

  evtTrig2 = hlt_2mu3_quark_v2 *1
           + hlt_2mu6p5_dis    *10
           + hlt_2mu7_dis      *100
           + hlt_2mu3p5_dis    *1000
           + hlt_2mu4_dis      *10000
           + hlt_2mu6p5_pro    *100000
           + hlt_2mu0_pro_v1   *1000000
           + hlt_2mu0_pro_v5   *10000000;

  evtTrig3 = hlt_2mu0_pro_v6       *1
           + hlt_2mu0_pro_NoVtx_v2 *10
           + hlt_2mu0_pro_NoVtx_v3 *100
           + hlt_mu0trk0           *1000
           + hlt_mu3trk0           *10000
           + hlt_mu0trkmu0         *100000
           + hlt_mu3trkmu0         *1000000
           + hlt_mu0trkmu0OST      *10000000;

  evtTrig4 = hlt_mu3trkmu0OST       *1
           + hlt_mu0trkmu0OST_tight *10           
           + hlt_L1muOpen           *100
           + hlt_L12muOpen          *1000
           + hlt_L12muOpenTight     *10000
           + hlt_2mu0L2             *100000
           + hlt_DimuPsi2s          *1000000
           + hlt_DimuLowMd          *10000000;                                     

  if( MyPrint_ ) 
      std::cout << hlt_mu3<<"--------------------------------- "<<std::endl;     
}
//bool DimuPiK::HecDimuVtx(double Muu, std::vector<RefCountedKinematicParticle> uuParticle,reco::Particle::LorentzVector uuP4){ //--Hm 03/05/13
bool DimuPiK::HecDimuVtx(double Muu, KinematicFitDriver uuRec, reco::Particle::LorentzVector uuP4){
    uuP4fit.SetPxPyPzE(0.0,0.0,0.0,0.0);
    //KinematicFitDriver uuRec( uuParticle, "MuMu" );
    if( !uuRec.isValid() )            //--Valid Dimuon vertex
      return false;

//--Output Vertex fit variables 
    uuVtx[0] = uuRec.mass();   //--MuuVFit
    uuVtx[1] = ChiSquaredProbability( uuRec.chi2(), uuRec.ndof() );  //--uuV_CL
    uuVtx[2] = -1;            //--MuuMFit     Mass Fitted Dimuon mass 
    uuVtx[3] = -1;            //--uuVM_CL     Confidence Level uu Vtx & Mass constrains

    //--Here Dimu Mass constraint to massC_value (= jpsi)
    double massC_value = 0;
    bool uu1sMcut = std::abs( uuVtx[0] - jpsi_mass_c )<0.2;    //--200 MeV
    bool uu2sMcut = std::abs( uuVtx[0] - psi2s_mass_c)<0.2;    //--200 MeV
    bool uuPhiMcut= std::abs( uuVtx[0] - 1.020)<0.075;         //--75 MeV  Phi
    if( uu1sMcut )massC_value = jpsi_mass_c;
    if( uu2sMcut )massC_value = psi2s_mass_c;
    if( uuPhiMcut)massC_value = 1.020;
    if( uu1sMcut || uu2sMcut || uuPhiMcut){	  
      uuRec.AddMassConstraint( massC_value, massC_value*1.e-6 );	      
      if( uuRec.isValid() ){ 
        uuVtx[2] = uuRec.mass();
        uuVtx[3] = ChiSquaredProbability( uuRec.chi2(), uuRec.ndof() );
      }  //--End good J/Psi or Psi(2s) mass constrain	     
    }    //--End uu1sMcut || uu2sMcut
       
    //--Dimuon fitted variables
    reco::Particle::LorentzVector  uuV_P4 = uuRec.P4();             //--DiMu 4-Momentum
    reco::Particle::LorentzVector iMuV_P4 = uuRec.P4FirstChild();   //--Mui
    reco::Particle::LorentzVector jMuV_P4 = uuRec.P4NextChild();    //--Muj
    reco::Particle::Point uuV_vx = uuRec.VertexDecay();	            //--Vx,Vy,Vz 
    reco::Vertex::CovarianceMatrix uuV_Cov = uuRec.VertexCov();
    uuVtx[4] = uuV_vx.x();				     //--uuVx	
    uuVtx[5] = uuV_vx.y();				     //--uuVy	
    uuVtx[6] = uuV_vx.z();				     //--uuVz	
    uuVtx[7] = sqrt( uuVtx[4]*uuVtx[4] +  uuVtx[5]*uuVtx[5] +  uuVtx[6]*uuVtx[6] );    //--uuDist	
    uuVtx[8] = sqrt( uuVtx[4]*uuVtx[4] +  uuVtx[5]*uuVtx[5] );		               //--uuR	

    if( MyPrint_ ){
        std::cout<<"Hec: Collection uu M= "<<Muu<<" Mass Vtx ReFit: "<<uuVtx[0]<<std::endl;
        std::cout<<"Hec: CL = "<<uuVtx[1]
	         <<" Vtx (x,y,z): ("<<uuVtx[4]<<","<<uuVtx[5]<<","<<uuVtx[6]<<")"<<std::endl;
    }     
    //--Check dimuon invariant mass
    uuVtx[9]  = iMuV_P4.E()  + jMuV_P4.E();				     //--uuE 
    uuVtx[10] = sqrt( muon_mass_c*muon_mass_c + iMuV_P4.P()*iMuV_P4.P() );   //--iuE 
    uuVtx[11] = sqrt( muon_mass_c*muon_mass_c + jMuV_P4.P()*jMuV_P4.P() );   //--juE
    uuVtx[12] = uuVtx[10] + uuVtx[11];                                       //--uuE = iuE + juE
    uuVtx[13] = iMuV_P4.Px() + jMuV_P4.Px();	     //--Momentum at Vtx   uu Px 
    uuVtx[14] = iMuV_P4.Py() + jMuV_P4.Py();	     //--uuPy
    uuVtx[15] = iMuV_P4.Pz() + jMuV_P4.Pz();	     //--uuPz
    uuVtx[16] = iMuV_P4.Px();         //--iPx  
    uuVtx[17] = iMuV_P4.Py();         //--iPy 
    uuVtx[18] = iMuV_P4.Pz();         //--iPz  
      
    uuVtx[19] = jMuV_P4.Px();	     //--jPx
    uuVtx[20] = jMuV_P4.Py();	     //--jPx
    uuVtx[21] = jMuV_P4.Pz();	     //--jPx
 
    uuVtx[22] = sqrt( uuVtx[9]*uuVtx[9] - uuVtx[13]*uuVtx[13] - uuVtx[14]*uuVtx[14]  - uuVtx[15]*uuVtx[15] );	//--dimuon at Vertex         myuuM 
    uuVtx[23] = sqrt( uuV_P4.E()*uuV_P4.E() - uuV_P4.Px()*uuV_P4.Px()	     //--Mass constrains Momentum  myfuuM
                                            - uuV_P4.Py()*uuV_P4.Py()
                                            - uuV_P4.Pz()*uuV_P4.Pz() );
    uuP4fit.SetPxPyPzE( uuV_P4.Px(), uuV_P4.Py(), uuV_P4.Pz(), uuV_P4.E() );
    
    uuVtx[24] = uuP4.eta();                         //--uuEta
    uuVtx[25] = uuP4.phi();                         //--uuPhi
    uuVtx[26] = sqrt( uuVtx[24]*uuVtx[24] + uuVtx[25]*uuVtx[25]);  //--uuDR 

    if( MyPrint_ ){
        std::cout<<"Muu="<<Muu<<" Vfit[0] M="<<uuVtx[0]<<" [23]M:"<<uuVtx[23]<<" [22]M"<<uuVtx[22]<<" MassC:"<<uuP4fit.M()<<std::endl;
        std::cout<<" "<<uu1sMcut<<" "<<uu2sMcut<<" "<<uuPhiMcut<<" "<<massC_value<<std::endl;
    }
    return true;                   
}
void DimuPiK::HecTrkVar(int n, double MyMass, reco::TrackRef iTrk){
    TrkVar[n][0]  = iTrk->charge();  
    TrkVar[n][1]  = iTrk->px();
    TrkVar[n][2]  = iTrk->py();
    TrkVar[n][3]  = iTrk->pz();   
    TrkVar[n][4]  = sqrt( MyMass*MyMass + TrkVar[n][1]*TrkVar[n][1] + TrkVar[n][2]*TrkVar[n][2] + TrkVar[n][3]*TrkVar[n][3] );    //--Energy
    TrkVar[n][5]  = iTrk->pt();       //--sqrt( Px*Px + Py*Py )
    TrkVar[n][6]  = iTrk->p();
    TrkVar[n][7]  = iTrk->eta();
    TrkVar[n][8]  = iTrk->phi();

    TrkVar[n][9]  = iTrk->numberOfValidHits();  //--piNtrkHits
    TrkVar[n][10] = iTrk->d0();
    TrkVar[n][11] = iTrk->dz();
}
double DimuPiK::HecRcone(int n, double uuPhi, double uuEta){
    double jpiuuR = 999;
    double phi1 = uuPhi;     //--Calculate delta phi between dimuon and the pion
    double phi2 = TrkVar[n][8];
    if( phi1<0 )phi1 = phi1 + 6.2832;
    if( phi2<0 )phi2 = phi2 + 6.2832;
    double jD_Phi = phi1 - phi2;
    if( jD_Phi>3.1416 )
        jD_Phi = jD_Phi - 6.2832;
    else if( jD_Phi<-3.1416 )
	     jD_Phi = jD_Phi + 6.2832;
	     
    jpiuuR = sqrt( (uuEta-TrkVar[n][7])*(uuEta-TrkVar[n][7]) + jD_Phi*jD_Phi );
    return jpiuuR;
}
//define this as a plug-in
DEFINE_FWK_MODULE(DimuPiK);
