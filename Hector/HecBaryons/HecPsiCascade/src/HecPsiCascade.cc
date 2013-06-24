// -*- C++ -*-
//
// Package:    HecPsiCascade
// Class:      HecPsiCascade
// 
/**\class HecPsiCascade HecPsiCascade.cc HecBaryons/HecPsiCascade/src/HecPsiCascade.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Hector Mendez
//         Created:  Tue Jun 14 08:50:59 CDT 2011
// $Id: HecPsiCascade.cc,v 1.18 2011/10/18 16:26:12 mendez Exp $
//
//


// system include files
#include <memory>

// user include files
#include "HecBaryons/HecPsiCascade/interface/HecPsiCascade.h"         //--Hec

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"             //--Hec

#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
//#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"             //--Hec

#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/CaloTopology/interface/EcalBarrelTopology.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonTrackLinks.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/CaloMuon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"       //track   collection
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"    //cascade colletcion
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"    //cascade daughter track

#include "Analyzers/CascadeProducer/interface/PrimaryInfo.h"             //--Eduardo Cascade & Primary
#include "Analyzers/CascadeProducer/interface/VertexRefit.h"
#include "Analyzers/CascadeProducer/interface/KinematicFitDriver.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "Analyzers/CascadeProducer/interface/ClosestApproachOnHelixLine.h"

#include "Analyzers/CascadeProducer/interface/CascadeCandidates.h" 
 

//--kinemactic fitter
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

//--Bfield
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "Analyzers/CascadeProducer/interface/masses.h"
//
// class declaration  <--I Included this in my header file  interface/HecPsiCascade.h   Hec
//
//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
HecPsiCascade::HecPsiCascade(const edm::ParameterSet& iConfig)
:
 trackTags_         (iConfig.getParameter<edm::InputTag>("tracks"        )),
 electronCollection_(iConfig.getParameter<edm::InputTag>("gsfElectrons"  )),
 theMuonsLabel_     (iConfig.getParameter<edm::InputTag>("MuonsLabel"    )),
 theCaloMuonsLabel_ (iConfig.getParameter<edm::InputTag>("CaloMuonsLabel")),
 minTracks_         (iConfig.getUntrackedParameter<unsigned int>("minTracks",0)),
 check_elC_         (iConfig.getUntrackedParameter<bool>("check_elC",true)),
 Pe_match_          (iConfig.getUntrackedParameter<double>("Pe_match",1000000.0)),
 mass_constrain_    (iConfig.getParameter<std::string>("Mass_Constrain")),
 onlyDiMu_          (iConfig.getUntrackedParameter<bool>("onlyDiMu","False")),
 MyPrint_           (iConfig.getUntrackedParameter<bool>("MyPrint","False")),
 CasAlgo_           (iConfig.getUntrackedParameter<std::string>("CasAlgo","MyCascade")),
 CasDecayName_      (iConfig.getUntrackedParameter<std::string>("CasDecayName","Cascade")),
 VeeAlgo_           (iConfig.getUntrackedParameter<std::string>("VeeAlgo","generalV0Candidates"))
{
  
  //--set mass for mass fit constrain to be combine with Xi
  if( mass_constrain_ == "m_jpsi" )
      massC_value = jpsi_mass_c;
  massC_sigma = massC_value * 1.e-6;
  
  std::cout<<"Hec: Hist & Param Initialization: HecPsiCascade::HecPsiCascade: Constrain Mas -->"
           <<massC_value <<std::endl;
   
   //now do what ever initialization is needed
  edm::Service<TFileService> fs;
  //--Diagnostic and testing Histograms
  int Mnch=4000;
  double Mxi=0., Mxf=200.;
  histo_trk  = fs->make<TH1D>("Ntrk"       , "Ntrk"       , 200 , 0 , 200 );
  histo_Pti  = fs->make<TH1D>("Pti"        , "Pt All"     , 200 , 0 , 100 );
  histo_ij   = fs->make<TH1D>("Massij"     , "Mass ij All", Mnch , Mxi , Mxf );
  histo_ijPt = fs->make<TH1D>("MassijPt"   , "Mass ij Pt" , Mnch , Mxi , Mxf );
  
  histo_nmu    = fs->make<TH1D>("NMuons"      , "NMuons"   , 200 , 0 , 200 );
  histo_muId   = fs->make<TH1D>("muId"        , "Muon Id"   , 10 , 0 , 10 );
  histo_dimu   = fs->make<TH1D>("MassDiMu"    , "Mass DiMu All"   , Mnch , Mxi , Mxf );         //--50 MeV/channel
  histo_dimuGG = fs->make<TH1D>("MassDiMuGG"  , "Mass DiMuGG"     , Mnch , Mxi , Mxf );
  histo_dimuTT = fs->make<TH1D>("MassDiMuTT"  , "Mass DiMuTT"     , Mnch , Mxi , Mxf );
  histo_dimuGT = fs->make<TH1D>("MassDiMuGT"  , "Mass DiMuGT"     , Mnch , Mxi , Mxf );
  histo_dimu00 = fs->make<TH1D>("MassDiMu00"  , "Mass DiMu00"     , Mnch , Mxi , Mxf );
  histo_dimu01 = fs->make<TH1D>("MassDiMu01"  , "Mass DiMu01"     , Mnch , Mxi , Mxf );
  histo_dimuCan= fs->make<TH1D>("MassDiMuCan" , "Mass DiMuCan Xi" , Mnch , Mxi , Mxf );
  
  histo_LoS    = fs->make<TH1D>("DiMu LoS"   , "DiMu LoS All", 100 , 0 , 100 );
  
  histo_dimuX0  = fs->make<TH1D>("MassDiMuX0"   , "Mass DiMuX All Pmu", Mnch , Mxi , Mxf );
  histo_dimuX1  = fs->make<TH1D>("MassDiMuX1"   , "Mass DiMuX All Pvx", Mnch , Mxi , Mxf );
  histo_dimuX2  = fs->make<TH1D>("MassDiMuX2"   , "Mass DiMuX All Pmc", Mnch , Mxi , Mxf );
  histo_dimuX3  = fs->make<TH1D>("MassDiMuX3"   , "Mass DiMuX All Pmc", Mnch , Mxi , Mxf );

  histo_cas01 = fs->make<TH1D>("MassCas01"  , "Mass Cas All"            , 160 , 1.25 , 1.41 );  //--1 MeV/channel
  histo_cas02 = fs->make<TH1D>("MassCas02"  , "Mass Cas in DiMu"        , 160 , 1.25 , 1.41 );
  histo_cas03 = fs->make<TH1D>("MassCas03"  , "Mass Cas in DiMu nTrk"   , 160 , 1.25 , 1.41 );
  histo_cas04 = fs->make<TH1D>("MassCas04"  , "Mass Cas in DiMu Sel"    , 160 , 1.25 , 1.41 );
  histo_lam01 = fs->make<TH1D>("MassLam01"  , "Mass Lam All"       , 200 , 1.07 , 1.17 );       //--0.5 MeV/channel
  histo_lam02 = fs->make<TH1D>("MassLam02"  , "Mass Lam Dimu"      , 200 , 1.07 , 1.17 );
  
//histo_ij0  = fs->make<TH1D>("Massij0"   , "Mass ij 0e" , 400 , 0 , 200 );
//histo_ij1  = fs->make<TH1D>("Massij1"   , "Mass ij 1e" , 400 , 0 , 200 );
//histo_ij2  = fs->make<TH1D>("Massij2"   , "Mass ij 2e" , 400 , 0 , 200 );
  histo_elec = fs->make<TH1D>("Nelec"      , "Nelec"      , 100 , 0 , 100 );
//histo_NeTrk= fs->make<TH1D>("NeTrk"      , "NeTrk"      , 100 , 0 , 100 );
//histo_DNel = fs->make<TH1D>("DNel"       , "Delta Nel"  ,  10 ,- 5,   5 );
//histo_Delp = fs->make<TH1D>("Delp"       , "Delta el p" , 100 ,-20,  80 );
//histo_DelpC= fs->make<TH1D>("DelpC"      , "DeltaP CUT" , 100 ,-20,  80 );
//histo_Delq = fs->make<TH1D>("Delq"       , "Delta el q" ,   6 , -3,   3 );
  histo_eiej = fs->make<TH1D>("Masseieja"  , "Mass eiej"  , Mnch , Mxi , Mxf );
//histo_DM   = fs->make<TH1D>("DM"         , "Delta Mass" , 100 ,-1/5000,1/5000 );
  histo_eEn  = fs->make<TH1D>("elecEnergy" , "elec Energy", 800 , 0 , 200 );
  histo_EoP  = fs->make<TH1D>("EoP"        , "elec E/P"   , 100 , 0 , 5 );
//histo_Clu  = fs->make<TH1D>("phCluster" , "Nclus"      , 100 , 0 , 100 );
//histo_Eni  = fs->make<TH1D>("phEnergy"  , "Energy i"   , 800 , 0 , 200 );
  histo_uu_prim_cos = fs->make<TH1D>("uuP_cos", "uuPrim Cosine", 220 , -1.1 , 1.1 );
}


HecPsiCascade::~HecPsiCascade()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

//
// member functions
//

// ------------ method called to for each event  ------------
void
HecPsiCascade::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //using namespace edm;
  //--Get track Collection
  //using reco::TrackCollection;
  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByLabel(trackTags_,tracks); int Ntk = tracks->size();
  //const reco::TrackCollection* trackCollection = tracks.product();

  //--Get electron Collection [from: CMSSW/RecoEgamma/Examples/plugins/GsfElectronDataAnalyzer.cc]
  edm::Handle<reco::GsfElectronCollection> gsfElectrons;
  iEvent.getByLabel(electronCollection_,gsfElectrons);

  //--Get Muon Collection
  edm::Handle<reco::MuonCollection> allmuons;
  iEvent.getByLabel(theMuonsLabel_, allmuons); int allMu = allmuons->size();
  
  //--Get CaloMuon Collection
  edm::Handle<reco::CaloMuonCollection> allcalmuons;
  iEvent.getByLabel(theCaloMuonsLabel_, allcalmuons);int allCalMu = allcalmuons->size();
  
  //--Handles for Primary Vertex  (from Eduardo)
  PrimaryInfo primary(iEvent,iConfig);
  
  //--Handle fron Primary just to get the number of rec Prim vertex (Sep 10, 11)
  edm::Handle<reco::VertexCollection> privtxs;
//iEvent.getByLabel(thePrimaryVertexLabel, privtxs); int allPrim = privtxs->size();
//iEvent.getByLabel("PrimaryCollection", privtxs); int allPrim = privtxs->size();
  iEvent.getByLabel("offlinePrimaryVertices", privtxs); int allPrim = privtxs->size();
  
  //--Cascades and Lambdas from Eduardo                    
  std::vector<reco::VertexCompositeCandidate> theCascades;
  std::vector<reco::VertexCompositeCandidate> theVees;
  CascadeCandidates MyXiTrack( iEvent );
  if( !onlyDiMu_ ){
    //--Handles for Std Vees (Lambdas) and Load them into a vector of composite candidates
    edm::Handle<reco::VertexCompositeCandidateCollection> theVeeHandle;
    iEvent.getByLabel(VeeAlgo_, "Lambda", theVeeHandle);
    theVees.insert( theVees.end(), 
                    theVeeHandle->begin(),
		    theVeeHandle->end() );
          
    //--Handles for Cascades and Load them into a vector of composite candidates
    edm::Handle<reco::VertexCompositeCandidateCollection> theCasHandle;
    iEvent.getByLabel(CasAlgo_, CasDecayName_, theCasHandle);
    theCascades.insert( theCascades.end(), 
                        theCasHandle->begin(),
                        theCasHandle->end() );
    //--Run CascadeCandidates  [Oct 12, 2011]
    MyXiTrack.init(trackTags_,CasAlgo_,CasDecayName_,iEvent,iSetup);  

  }  //--End if( !onlyDiMu_ )
  int allCascade = theCascades.size();
  int allLambda  = theVees.size();
  
  //--Handles for beamspot used for primary refit
  reco::BeamSpot beamSpot;
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByLabel("offlineBeamSpot", beamSpotHandle);
  if( beamSpotHandle.isValid() )beamSpot = *beamSpotHandle;

  //--Handles for B-field
  edm::ESHandle<MagneticField> bFieldHandle;
  iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);
  const MagneticField *BField = bFieldHandle.product();
  
  //--Handles for Tracker Geometry
  edm::ESHandle<GlobalTrackingGeometry> globTkGeomHandle;
  iSetup.get<GlobalTrackingGeometryRecord>().get(globTkGeomHandle);
   
  //--Handles for Transient Track Builder 
  edm::ESHandle<TransientTrackBuilder> theTTBHandle;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBHandle);
   
  //--Fill Ntuple (Branch: evt)
  fill_evt( iEvent, Ntk, allMu, allCalMu, allCascade, allLambda, allPrim );
  /*edm::LogInfo("HecPsiCascade")<<"Hec: Run: "<<evt_.runNb<<" Evt: "<<evt_.eventNb
               <<" Lumi:"<<evt_.lumiBlock<<" Ntk: "<<evt_.Ntk<<" Mu:"<<evt_.allMu
               <<" CalMu:"<<evt_.allCalMu <<" NXi: "<<evt_.allCascade<<" NLambda: "<<evt_.allLambda;*/
                                 
  //--Loop over rec electrons
  histo_elec->Fill( gsfElectrons->size() );
  for(reco::GsfElectronCollection::const_iterator gsfIter=gsfElectrons->begin();
                                                  gsfIter!=gsfElectrons->end(); gsfIter++){     
     double iCharge = gsfIter->charge();
     histo_eEn->Fill( gsfIter->energy() );
     histo_EoP->Fill( gsfIter->eSuperClusterOverP() );
     
     //--gsfIter->overlap();
     for(reco::GsfElectronCollection::const_iterator gsfJter=gsfIter+1; 
                                                     gsfJter!=gsfElectrons->end(); gsfJter++){     
        double jCharge = gsfJter->charge();
        if( iCharge * jCharge >=0 )continue;
        
        //--Di-electron Invariant Mass
        double Masseiej = ( gsfIter->p4() +  gsfJter->p4() ).M();
        histo_eiej->Fill( Masseiej );
        
        //--The next 20 lines are equivalent to the 2 previous lines to calculate di-electron mass
/*      double Pexi = gsfIter->px();
        double Peyi = gsfIter->py();
        double Pezi = gsfIter->pz();
        double Peti = gsfIter->pt();   //--sqrt( Pexi*Pexi + Peyi*Peyi )                 
        double Enei = sqrt( elec_mass_c*elec_mass_c + Pexi*Pexi + Peyi*Peyi + Pezi*Pezi );
        
        double Pexj = gsfJter->px();
        double Peyj = gsfJter->py();
        double Pezj = gsfJter->pz();
        double Enej  = sqrt( elec_mass_c*elec_mass_c + Pexj*Pexj + Peyj*Peyj + Pezj*Pezj );
        
        double Peex = Pexi + Pexj;
        double Peey = Peyi + Peyj;
        double Peez = Pezi + Pezj;
        double Enee = Enei + Enej;
         
        double Massee = Enee*Enee - Peex*Peex - Peey*Peey - Peez*Peez;
        double Massij = -1.0;
        if( Massee >= 0)
            Massij = sqrt( Massee );
        double DeltaM = Masseiej - Massij;
        histo_DM->Fill( DeltaM );
        
        if( DeltaM != 0 )
          std::cout<<"Hec: Wrong Di-electron Mass"
                   <<"Event " << iEvent.id() <<" Delta M"<< DeltaM << std::endl;      
*/
     }    //--electrons gsfJter
  }       //--electrons gsfIter
   
  //--Looks at tracks
  histo_trk->Fill( tracks->size() ); 
  double mass_uno = muon_mass_c; 
  for(reco::TrackCollection::const_iterator itTrack = tracks->begin(); itTrack != tracks->end(); ++itTrack){
     int iCharge = itTrack->charge();
     
     double Pti = itTrack->pt();
   //double Pki  = itTrack->p();
     
     histo_Pti->Fill( Pti );
    
     for(reco::TrackCollection::const_iterator jtTrack = itTrack+1; jtTrack != tracks->end(); ++jtTrack){
        int jCharge = jtTrack->charge();
        
      //double Ptj = jtTrack->pt();
      //double Pkj = jtTrack->p();

        if( iCharge * jCharge >=0 )continue;
                  
        double Pkxi = itTrack->px();
        double Pkyi = itTrack->py();
        double Pkzi = itTrack->pz();
        double Pkti = itTrack->pt();   //--sqrt( Pxi*Pxi + Pyi*Pyi )
        double Enki = sqrt( mass_uno*mass_uno + Pkxi*Pkxi + Pkyi*Pkyi + Pkzi*Pkzi );
         
        double Pkxj = jtTrack->px();
        double Pkyj = jtTrack->py();
        double Pkzj = jtTrack->pz();
        double Pktj = jtTrack->pt();   //--sqrt( Pxj*Pxj + Pyj*Pyj )
        double Enkj = sqrt( muon_mass_c*muon_mass_c + Pkxj*Pkxj + Pkyj*Pkyj + Pkzj*Pkzj );
         
        double Pkkx = Pkxi + Pkxj;
        double Pkky = Pkyi + Pkyj;
        double Pkkz = Pkzi + Pkzj;
        double Enkk = Enki + Enkj;
         
        double Masskk = Enkk*Enkk - Pkkx*Pkkx - Pkky*Pkky - Pkkz*Pkkz;
        double Massij = -1.0;
        if( Masskk >= 0)
            Massij = sqrt( Masskk );
        
        histo_ij->Fill( Massij );
        if( Pkti>=2 && Pktj>=2 )histo_ijPt->Fill( Massij );
               
     }     //--end for(TrackCollection::const_iterator jtTrack = itTrack+1;         Kaon j
  }        //--end for(TrackCollection::const_iterator itTrack = tracks->begin();   Kaon i
 
  //--Looks on Muons  [July 25, 2011]
  histo_nmu->Fill( allmuons->size() );
  //std::vector<int> theMuonTrkIndexes;
  //theMuonTrkIndexes.push_back(iMuon->innerTrack().index());
  for(reco::MuonCollection::const_iterator iMuon = allmuons->begin(); iMuon != allmuons->end(); ++iMuon){
     IuuC_.iCharge = iMuon->charge();
     uuC_.iEta     = iMuon->eta();
     uuC_.iPhi     = iMuon->phi();
     uuC_.iPt      = iMuon->pt();
     uuC_.iP       = iMuon->p();
     uuC_.iCal     = iMuon->caloCompatibility();
     uuC_.iSeg     = muon::segmentCompatibility(*iMuon);
     uuC_.iIso     = iMuon->isolationR03().sumPt;
     reco::TrackRef    iMuTrackRef = iMuon->innerTrack(); //--.index();
     IuuC_.iNtrkHits = -1;
     uuC_.id0        = -1;
     uuC_.idz        = -1;
     
     int muId=0, muSA=0, muGL=0, muTK=0;          //--Muon Id
     if( iMuon->isStandAloneMuon() )muSA = 1;     //--2^0
     if( iMuon->isGlobalMuon()     )muGL = 10;    //--2^1
     if( iMuon->isTrackerMuon()    )muTK = 100;   //--2^2
     IuuC_.imuId = muTK + muGL + muSA;
     
     muId = muTK*pow(2,2)/100 + muGL*pow(2,1)/10 + muSA*pow(2,0)/1;
     histo_muId->Fill( muId );
     
     if( MyPrint_ ){
         std::cout<<"Hec: MuId: TK:"<<muTK<<" GL:"<<muGL<<" SA:"<<muSA<<" Dec:"<<muId<<" Bin:"<<IuuC_.imuId<<std::endl;
         std::cout<<"Hec:  Mu Compatibility:"<<uuC_.iCal<<" MuSeg:"<<uuC_.iSeg<<" MuIso:"<<uuC_.iIso<<std::endl;
     }
     
     IuuC_.iGL = 0;
     //if(  iMuon->isGlobalMuon()  && iMuon->isTrackerMuon() )uuC_.iGL = 1;
     if(  iMuon->isGlobalMuon()  )IuuC_.iGL = 1;
     if(!(iMuon->isGlobalMuon()) && iMuon->isTrackerMuon() )IuuC_.iGL = 2;
     
     for(reco::MuonCollection::const_iterator jMuon = iMuon+1; jMuon != allmuons->end(); ++jMuon){
        IuuC_.jCharge = jMuon->charge();
        uuC_.jEta     = jMuon->eta();
        uuC_.jPhi     = jMuon->phi();
        uuC_.jPt      = jMuon->pt();
        uuC_.jP       = jMuon->p(); 
        uuC_.jCal     = jMuon->caloCompatibility();
        uuC_.jSeg     = muon::segmentCompatibility(*jMuon);
        uuC_.jIso     = jMuon->isolationR03().sumPt;  
        reco::TrackRef    jMuTrackRef = jMuon->innerTrack();  //--.index();
        IuuC_.jNtrkHits = -1;
        uuC_.jd0        = -1;
        uuC_.jdz        = -1;
     
        muId=0, muSA=0, muGL=0, muTK=0;   
        if( jMuon->isStandAloneMuon() )muSA = 1;     //--2^0
        if( jMuon->isGlobalMuon()     )muGL = 10;    //--2^1
        if( jMuon->isTrackerMuon()    )muTK = 100;   //--2^2
        IuuC_.jmuId = muTK + muGL + muSA;
               
        IuuC_.jGL = 0;
        //if(  jMuon->isGlobalMuon()  && jMuon->isTrackerMuon() )uuC_.jGL = 1;
        if(  jMuon->isGlobalMuon()  )IuuC_.jGL = 1;
        if(!(jMuon->isGlobalMuon()) && jMuon->isTrackerMuon() )IuuC_.jGL = 2;
        
        if( IuuC_.iCharge * IuuC_.jCharge >=0 )continue;
        uuC_.Muu = ( iMuon->p4() +  jMuon->p4() ).M();
         
        histo_dimu->Fill( uuC_.Muu );
        if(  IuuC_.iGL==1 && IuuC_.jGL==1 )histo_dimuGG->Fill( uuC_.Muu );  //--Global-Global
                                                                            //--Global no tracker
        if(  IuuC_.iGL==2 && IuuC_.jGL==2 )histo_dimuTT->Fill( uuC_.Muu );  //--Tracker-Tracker [No Global]
        if( (IuuC_.iGL==1 && IuuC_.jGL==2) 
          ||(IuuC_.iGL==2 && IuuC_.jGL==1))histo_dimuGT->Fill( uuC_.Muu );  //--1 Global-1 tracker
        if(  IuuC_.iGL==0 && IuuC_.jGL==0 )histo_dimu00->Fill( uuC_.Muu );  //--No global No tracker (calo - calo ?)
        if( (IuuC_.iGL==0 && IuuC_.jGL!=0) 
          ||(IuuC_.iGL!=0 && IuuC_.jGL==0))histo_dimu01->Fill( uuC_.Muu );  //--1 Global or 1 Tracker whith  1 calo
         
        //--Skip Calo Muon (No Global Muon and No Tracker Muon)
        //--if( uuC_.iGL==0 || uuC_.jGL==0 )continue;
        //--my Flags to read my N-Tuple output are:
        uuC_.dca        = -999;  //--di-muon Closest Aproach
        uuFit_.MuuFit   = -1;    //--re-Fitted Dimuon Mass with Vertex constrain 
        uuFit_.uuL      = -1;    //--L:distance between Sec and prim (Sigma: Error on L on Sec-Primary re-Fit). Primary re-Fit
        uuFit_.uuLoS    = -1;    //--L over Sigma betwen Dimuon Vtx (Sec) and Prim (Prim without Dimuon). Primary re-Fit        
        uuFit_.MuuMFit  = -1;    //--re-Fitted Dimuon Mass with Vertex and Mass constrain
        uuFit_.VuuVpCos = -999;  //--uu Sec Vtx Pointing to the Primary. 
        ILam_.uuPair    = -1;    //--Cascade-dimuon combinations
        
        if( IuuC_.iGL!=0 && IuuC_.jGL!=0 ){    //--Here only good muons
          IuuC_.iNtrkHits = iMuTrackRef->numberOfValidHits();
          uuC_.id0        = iMuTrackRef->d0();
          uuC_.idz        = iMuTrackRef->dz();
          
          IuuC_.jNtrkHits = jMuTrackRef->numberOfValidHits();
          uuC_.jd0        = jMuTrackRef->d0();
          uuC_.jdz        = jMuTrackRef->dz();
          if( MyPrint_ )
              std::cout<<"Hec: ValHits:"<<IuuC_.iNtrkHits<<" d0:"<<uuC_.id0<<" dz:"<<uuC_.idz<<std::endl;

          float muon_sig   = muon_mass_c*1.e-6;
          float pion_sig   = pion_mass_c*1.e-6;
          float proton_sig = proton_mass_c*1.e-6;
          float chi = 0., ndf = 0.;
          
          //--Re-fit dimuons to Vertex
          reco::TransientTrack iMuTT( iMuTrackRef, &(*bFieldHandle) );
          reco::TransientTrack jMuTT( jMuTrackRef, &(*bFieldHandle) );
           
          KinematicParticleFactoryFromTransientTrack pFactory;
           
          std::vector<RefCountedKinematicParticle> uuParticle;
          uuParticle.push_back(pFactory.particle( iMuTT, muon_mass_c, chi, ndf, muon_sig ));
          uuParticle.push_back(pFactory.particle( jMuTT, muon_mass_c, chi, ndf, muon_sig ));
          
          // Trajectory states to calculate DCA for the 2 tracks
          FreeTrajectoryState iState = pFactory.particle( iMuTT, muon_mass_c, chi, ndf, muon_sig )->currentState().freeTrajectoryState();
          FreeTrajectoryState jState = pFactory.particle( jMuTT, muon_mass_c, chi, ndf, muon_sig )->currentState().freeTrajectoryState();

          //--Measure distance between tracks at their closest approach
          ClosestApproachInRPhi cApp;           
          double dca_uu = -9999;
          GlobalPoint cxPt;
          cApp.calculate(iState, jState);
          if( !cApp.status() ) {
              std::cout<<"uu bad status capp"<<std::endl;
              //--continue;    Need to check where this continue will go ??
          } else {
              dca_uu = fabs( cApp.distance() );
              cxPt = cApp.crossingPoint();
          }
          uuC_.dca  = dca_uu;       //--Flag  Closest Approach distance between muons
          uuC_.xcpt = cxPt.x();     //--Crosing point in X
          uuC_.ycpt = cxPt.y();     //--Crosing point in Y
          uuC_.zcpt = cxPt.z();     //--Crosing point in Z
          if( MyPrint_ ){
              std::cout<<"uu dca = "<<dca_uu<<" CroossPoint ="<<cxPt<<std::endl;
              std::cout<<" mui= "<<IuuC_.imuId<<" muj= "<<IuuC_.jmuId<<" mass ="<<uuC_.Muu<<std::endl;
              std::cout<<" -------------------------"  <<std::endl;
          }
          //--Cut on fiducial tracking volume [-120 cm < r < 120 cm  and -260 cm < Z < 260 cm]   r = sqrt( x^2 + y^2 )
          //--Pixel [-20 cm < r < 20 cm  and -60 cm < Z < 60 cm] to run vertex fit
          if( ( dca_uu >= 0. && dca_uu <= 5 ) && 
              ( sqrt( cxPt.x()*cxPt.x() + cxPt.y()*cxPt.y() ) < 120.) &&
                std::abs(cxPt.z()) < 300. ){
	 
          KinematicFitDriver uuRec( uuParticle, "MuMu" );
          if( !uuRec.isValid() ){
              if( MyPrint_ )
                  std::cout<<"Hec: Mass: uu: "<<uuC_.Muu<<" DiMuV BadFit:"<<std::endl;
              //--continue;   //--Goto the end of the if( uuC_.iGL!=0 && uuC_.jGL!=0 ) <--NO TRUE 
          } 
          else {
          //--Output Vertex fit variables
          reco::Particle::LorentzVector uuV_P4;
          reco::Particle::LorentzVector iMuV_P4;
          reco::Particle::LorentzVector jMuV_P4;
          reco::Particle::Point uuV_vx;
          double uuV_Chi2, uuV_Ndof;
          reco::Vertex::CovarianceMatrix uuV_Cov;
          uuFit_.MuuFit = uuRec.mass();          //--Flag Dimuon mass
          uuV_P4        = uuRec.P4();            //--DiMu 4-Momentum
          iMuV_P4       = uuRec.P4FirstChild();  //--Mui
          jMuV_P4       = uuRec.P4NextChild();   //--Muj
          uuV_vx        = uuRec.VertexDecay();   //--Vx,Vy,Vz 
          uuV_Cov       = uuRec.VertexCov();                                   
          uuV_Chi2      = uuRec.chi2();
          uuV_Ndof      = uuRec.ndof();
          uuFit_.uuV_Cl = ChiSquaredProbability( uuV_Chi2, uuV_Ndof );  //--Confidence Level uu Vtx
          uuFit_.uuVx   = uuV_vx.x();
          uuFit_.uuVy   = uuV_vx.y();
          uuFit_.uuVz   = uuV_vx.z();
          if( MyPrint_ )
              std::cout<<"Hec: Collection uu M= "<<uuC_.Muu<<" Vtx ReFit: "<<uuFit_.MuuFit<<std::endl;

          //--Primary Vertex Re-fitting to exclude the DiMuon from the primary
          reco::Vertex refitPrimVertex = *primary.BestVertex();
          std::vector<reco::TrackRef> ExclusionList;
          ExclusionList.push_back( iMuTrackRef );
          ExclusionList.push_back( jMuTrackRef );

          VertexRefit RefitPrimary( primary.BestVertex(), ExclusionList, bFieldHandle, beamSpot );
          refitPrimVertex = RefitPrimary.Refitted();
          if( !refitPrimVertex.isValid() ){
            if( MyPrint_ )
                std::cout<<"Hec: Mass: DiMu: "<<uuC_.Muu<<" Prim BadFit:"<<std::endl;
                //--continue;     //--Goto the end of the if( uuC_.iGL!=0 && uuC_.jGL!=0 ) <--NO TRUE
          } 
          else {
          //--Output Prim Vertex fit variables 
          double pV_Chi2, pV_Ndof;
          uuFit_.xpV   = refitPrimVertex.x();
          uuFit_.ypV   = refitPrimVertex.y();
          uuFit_.zpV   = refitPrimVertex.z();                                    
          pV_Chi2      = refitPrimVertex.chi2();
          pV_Ndof      = refitPrimVertex.ndof();
          uuFit_.pV_Cl = ChiSquaredProbability( pV_Chi2, pV_Ndof );    //--Confidence Level Prim
          
          //--Check separation from primary (L over Sigma)
          typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
          typedef ROOT::Math::SVector<double, 3> SVector3;

          SMatrixSym3D totalCov = refitPrimVertex.covariance() + uuV_Cov;
          SVector3 uu_PrimarySep3D( uuV_vx.x() - refitPrimVertex.x(),
                                    uuV_vx.y() - refitPrimVertex.y(),
                                    uuV_vx.z() - refitPrimVertex.z() );

          uuFit_.uuL     = ROOT::Math::Mag( uu_PrimarySep3D );                            //--Flag 3dim distance Dimuon-Prim
          double uuSigma = sqrt(ROOT::Math::Similarity( totalCov, uu_PrimarySep3D )) / uuFit_.uuL;

          uuFit_.uuLoS = uuFit_.uuL / uuSigma;                                            //--Flag 3dim distance/sigma  Dimuon-Prim
          histo_LoS->Fill( uuFit_.uuLoS );
          
          if( MyPrint_ )
              std::cout<<"Hec: L= "<<uuFit_.uuL<<" LoS= " << uuFit_.uuLoS<<" MassFitted"<<uuFit_.MuuFit<<std::endl;
              
          //--Here Pointing Secondary Vertex  [Sep 9, 2011]
          GlobalVector uu_PrimSec( uuV_vx.x() - refitPrimVertex.x(),
                                   uuV_vx.y() - refitPrimVertex.y(),
                                   uuV_vx.z() - refitPrimVertex.z() );                                   
          GlobalVector uu_P3( uuRec.P4().x(),
                              uuRec.P4().y(),
                              uuRec.P4().z() );                              
          double V_dot_Puu =uu_P3.dot(uu_PrimSec);
                                                                                 
          double cos_alpha = V_dot_Puu /( uu_PrimSec.mag() * uu_P3.mag() );
          uuFit_.VuuVpCos = cos_alpha;
          histo_uu_prim_cos->Fill( uuFit_.VuuVpCos ); 
          if( MyPrint_ )
              std::cout<<"Hec: Pointing cos alpha= "<<uuFit_.VuuVpCos<<std::endl;
               
          //--Here Dimu Mass constraint to massC_value (= jpsi)
          bool uuMcut = std::abs( uuFit_.MuuFit - massC_value )<0.2;    //--200 MeV
          if( uuMcut ){       //continue;  <--goes to the end of the loop
          
          uuRec.AddMassConstraint( massC_value, massC_sigma );          
          if( !uuRec.isValid() ){
            if( MyPrint_ )
                std::cout<<"Hec: DiMuon mass constraint fit failed!"<<std::endl;
                //--continue;   <--goes to the end of the loop. don't want this
          }
          else {
          //--Output mass constrain fit variables 
          double uuVM_Chi2, uuVM_Ndof;        
          uuFit_.MuuMFit = uuRec.mass();                       //--Flag  dimuon mass for mass constrain
          reco::Particle::LorentzVector uuVM_P4 = uuRec.P4();  //--DiMu 4-Momentum with mass constrain                                    
          uuVM_Chi2      = uuRec.chi2();
          uuVM_Ndof      = uuRec.ndof();
          uuFit_.uuVM_Cl = ChiSquaredProbability( uuVM_Chi2, uuVM_Ndof );    //--Confidence Level Prim
          uuFit_.uuVM_px = uuVM_P4.px();
          uuFit_.uuVM_py = uuVM_P4.py();
          uuFit_.uuVM_pz = uuVM_P4.pz();
//--          
//           RefCountedKinematicTree DiMuVertexFitTree = DiMuRec.Tree();
//           DiMuVertexFitTree->movePointerToTheTop();
//           RefCountedKinematicParticle DiMu_Particle = DiMuVertexFitTree->currentParticle();
//           RefCountedKinematicVertex DiMu_vertex_Fit = DiMuVertexFitTree->currentDecayVertex();
// 
//           KinematicParticleFitter massconstraintFitter;
//           KinematicConstraint * DiMu_c = new MassKinematicConstraint(massC_value,massC_sigma);
//           
//            //--add mass constraint to the vee fit
//           DiMuVertexFitTree = massconstraintFitter.fit(DiMu_c,DiMuVertexFitTree); //--use mass constrain
//          if( !DiMuVertexFitTree->isValid() ){ 
//             if( MyPrint_ )std::cout<< "DiMuon mass constraint fit failed!"<<std::endl;
//             continue;
//           }          
//           //--Save Re-fitted values form the mass constraint// 
//           uuFit_.MuuFitM = DiMuVertexFitTree->currentParticle()->currentState().mass();
//           reco::Particle::LorentzVector DiMuVM_P4(
//                  DiMuVertexFitTree->currentParticle()->currentState().globalMomentum().x(),
//                  DiMuVertexFitTree->currentParticle()->currentState().globalMomentum().y(),
//                  DiMuVertexFitTree->currentParticle()->currentState().globalMomentum().z(),
//                  DiMuVertexFitTree->currentParticle()->currentState().kinematicParameters().energy() );       
//------------------------           

          if( MyPrint_ )
              std::cout<<"Hec: DiMuon mass constraint "<<uuFit_.MuuFit<<" --> " <<uuFit_.MuuMFit<<std::endl;
             
          //--Here do Cascade  Xi- --> Lamda(P + Pi-) + Pi-  + DiMuon
          ILam_.uuPair = 0;                           //--Flag  dimuon-cascade combinations
          histo_dimuCan->Fill( uuC_.Muu );
          if( !onlyDiMu_ ){
            /* Find this decay with 3 vertex:
              Vx_P <-----uuLoS----->Vx_S ----------->u+u-
                        [casB]       \    Pi  proton
                               casLoS \     \/
                                       \    /
                                        \  /lamLoS
                                      Vx_T/
                                          \ 
                                           \pi
            */
            for( unsigned int casindex = 0; casindex < theCascades.size(); casindex++ ){
               Lam_.xiM = theCascades[casindex].mass();
               bool xiMcut = std::abs( Lam_.xiM - cascade_mass_c ) < 0.02;   //--20 MeV around Nominal Cascade
               histo_cas02->Fill( Lam_.xiM );

               //--Get TrackRef for the Cascade and skip it if equal to the DiMuon
               //--Pion from Cascade
               reco::TrackRef piontrk = (dynamic_cast<reco::RecoChargedCandidate *>
                                      ( theCascades[casindex].daughter(0) ) )->track();

               //--vee from cascade
               reco::VertexCompositeCandidate* CasVee = dynamic_cast<reco::VertexCompositeCandidate*>
                                   ( theCascades[casindex].daughter(1) );
               reco::TrackRef protonveetrk = (dynamic_cast<reco::RecoChargedCandidate *> 
                                         (CasVee->daughter(0) ) )->track();
               reco::TrackRef pionveetrk   = (dynamic_cast<reco::RecoChargedCandidate *> 
                                         (CasVee->daughter(1) ) )->track();
               //--Flags for Cascades 
               ILam_.xiCan      = -1;   //--No repeated tracks dimuon, pion and Proton
               Lam_.MlamVFit    = -1;   //--valid Lambda vertex 
               Lam_.MlamMFit    = -1;   //--valid Lambda mass constrains
               Cas_.MxiVFit     = -1;   //--valid Cascade vertex
               Cas_.MxiMFit     = -1;   //--valid mass constrain for Cascade
               Cas_.dca         = -999; //--distance closest approach  (cascade_b-->dimuon-cascade)
               Casb_.Mxi_bVFit  = -1;   //--valid Cascade_b vertex
               if( piontrk       == iMuTrackRef || piontrk      == jMuTrackRef )continue; //--goes end of the cascade loop [pion trk = muon trk --> skip]
               if( protonveetrk  == iMuTrackRef || protonveetrk == jMuTrackRef )continue; //--                             [proton trk = muon trk --> skip]
               if( pionveetrk    == iMuTrackRef || pionveetrk   == jMuTrackRef )continue; //--                             [pion trk from lambda = muon trk -->skip]
               ILam_.xiCan = 1;         //--Flag: I have 3 differents tracks

               histo_cas03->Fill( Lam_.xiM );                                 //--Cascade mass
               if( uuMcut )histo_cas04->Fill( Lam_.xiM );

               Lam_.lamM = CasVee->mass();                                    //--lambda mass
               bool lamMcut = std::abs( Lam_.lamM - lambda_mass_c ) < 0.01;   //--10 MeV  around Nominal Lambda
               histo_lam02->Fill( Lam_.lamM );
               
               //--Dimuon + Xi
               Lam_.uuXiM_0 = ( iMuon->p4() +  jMuon->p4() + theCascades[casindex].p4() ).M();    //--un-constrained dimuon mass + un-constrained cascade
               histo_dimuX0->Fill( Lam_.uuXiM_0 );
               
               //--Fitted(Vertex) Dimuon + Xi
               Lam_.uuXiM_1 = ( uuV_P4 + theCascades[casindex].p4() ).M();
               histo_dimuX1->Fill( Lam_.uuXiM_1 );
               
               //--Fitted(Vertex + Mass) Dimuon + Xi
               Lam_.uuXiM_2 = ( uuVM_P4 + theCascades[casindex].p4() ).M();
               histo_dimuX2->Fill( Lam_.uuXiM_2 );
               if( uuMcut && xiMcut && lamMcut )
                   histo_dimuX3->Fill( Lam_.uuXiM_2 );

               //--Fit Dimuon + Cascade to a single vertex
               //--vee (lamdba) from cascade
               reco::TransientTrack pionTveetrk(     pionveetrk, &(*bFieldHandle) );
               reco::TransientTrack protonTveetrk( protonveetrk, &(*bFieldHandle) );

               std::vector<RefCountedKinematicParticle> LamParticle;
               LamParticle.push_back(pFactory.particle( pionTveetrk,   pion_mass_c,   chi, ndf, pion_sig ));
               LamParticle.push_back(pFactory.particle( protonTveetrk, proton_mass_c, chi, ndf, proton_sig ));               
               KinematicFitDriver LamRec( LamParticle, "Lambda" );
               if( !LamRec.isValid() ){
                 if( MyPrint_ )
                     std::cout<< "Hec: Lambda Vertex fit failed!"<<std::endl;
                     continue;   //--goes end of the cascade loop
               } 
               //--Output lambda Vtx fit variables
               reco::Particle::Point lamV_vx;
               reco::Vertex::CovarianceMatrix lamV_Cov;
               Lam_.MlamVFit     = LamRec.mass();          //--Flag:  valid Lambda vertex 
               lamV_vx           = LamRec.VertexDecay();   //--Vx,Vy,Vz 
               lamV_Cov          = LamRec.VertexCov();                                   
               double lamV_Chi2  = LamRec.chi2();
               double lamV_Ndof  = LamRec.ndof();
               Lam_.lamV_Cl      = ChiSquaredProbability( lamV_Chi2, lamV_Ndof );  //--Confidence Level Lambda Vtx
               Lam_.lamVx        = lamV_vx.x();
               Lam_.lamVy        = lamV_vx.y();
               Lam_.lamVz        = lamV_vx.z();
               
               LamRec.AddMassConstraint( lambda_mass_c, lambda_mass_c*1.e-6 );     //--mass constrain (proton + pion) to Lambda mass       
               if( !LamRec.isValid() ){
                 if( MyPrint_ )
                     std::cout<< "Hec: Lambda mass constraint fit failed!"<<std::endl;
                     continue;   //--goes end of the cascade loop
               } 
               //--Output lambda mass constrain fit variables
               Lam_.MlamMFit      = LamRec.mass();                                    //--Flag:  valid Lambda mass constrains                           
               double lamVM_Chi2  = LamRec.chi2();
               double lamVM_Ndof  = LamRec.ndof();
               Lam_.lamVM_Cl      = ChiSquaredProbability( lamVM_Chi2, lamVM_Ndof );  //--Confidence Level
               reco::Particle::LorentzVector lamVM_P4 = LamRec.P4();
               Lam_.lamVM_px = lamVM_P4.px();
               Lam_.lamVM_py = lamVM_P4.py();
               Lam_.lamVM_pz = lamVM_P4.pz();
                
               //--Check separation from Dimuon (L over Sigma)
               typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D_uuL;
               typedef ROOT::Math::SVector<double, 3> SVector3_uuL;
 
               SMatrixSym3D_uuL totalCov_uuL = lamV_Cov + uuV_Cov;                    //--Separation Dimuon - Lambda
               SVector3_uuL uuL_PrimarySep3D( uuV_vx.x() - lamV_vx.x(),
                                              uuV_vx.y() - lamV_vx.y(),
                                              uuV_vx.z() - lamV_vx.z() );
 
               Lam_.uuL_L       = ROOT::Math::Mag( uuL_PrimarySep3D );                //--3 dim distance between Lambda-dimuon
               double uuL_Sigma = sqrt(ROOT::Math::Similarity( totalCov_uuL, uuL_PrimarySep3D )) / Lam_.uuL_L;
 
               Lam_.uuL_LoS = Lam_.uuL_L / uuL_Sigma;
               if( MyPrint_ )
                   std::cout<< "Hec: Lambda mass constraint"<<Lam_.MlamMFit<<std::endl;
 
               //--pi trk from cascade
               reco::TransientTrack pionTtrk( piontrk, &(*bFieldHandle) );
 
               //--Cascade Xi  Vertex  (Lambda + Pion)
               std::vector<RefCountedKinematicParticle> XiParticle;
               XiParticle.push_back(pFactory.particle( pionTtrk, pion_mass_c, chi, ndf, pion_sig ));
               XiParticle.push_back( LamRec.RKParent() );               
               KinematicFitDriver XiRec( XiParticle, "Cascade" );
               if( !XiRec.isValid() ){
                 if( MyPrint_ )
                     std::cout<< "Hec: Cascade Vertex fit failed!"<<std::endl;
                     continue;   //--goes end of the cascade loop
               }
               //--Output cascade Vtx fit variables
               reco::Particle::Point xiV_vx;
               reco::Vertex::CovarianceMatrix xiV_Cov;
               Cas_.MxiVFit     = XiRec.mass();          //--Flag:  valid Cascade vertex
               xiV_vx           = XiRec.VertexDecay();   //--Vx,Vy,Vz 
               xiV_Cov          = XiRec.VertexCov();                                   
               double xiV_Chi2  = XiRec.chi2();
               double xiV_Ndof  = XiRec.ndof();
               Cas_.xiV_Cl      = ChiSquaredProbability( xiV_Chi2, xiV_Ndof );       //--Confidence Level cacade Vtx
               Cas_.xiVx        = xiV_vx.x();
               Cas_.xiVy        = xiV_vx.y();
               Cas_.xiVz        = xiV_vx.z();
               
               //--Need to clean cascades  (used Eduardo Ntuplizer.  CascadeNtuplizer.cc  variable Xiip3d < 10
               
               //--Tag for cascade  (Oct 12, 2011)
               ILam_.NPixelXiTracks = MyXiTrack.NumPixelXiTracks(casindex);
               Cas_.XiCL            = MyXiTrack.XiCL(casindex);
               Cas_.rhoPi = MyXiTrack.RhoTrkDaughter(casindex);  //--rho of pion
               Cas_.rhoXi = MyXiTrack.RhoVertexDecay(casindex);  //--rho of cascade
               
               Cas_.IpPi  = MyXiTrack.BestSignificanceImpactParameterAtPixelTrack(casindex); //--Impact parameter of pion
                                                                                             //--in Xi vtx when NumPixel >0)
               reco::Vertex uuRec_0 = *(uuRec.RKVertex());
               Cas_.IpXiuu = MyXiTrack.SignificanceImpactParameter3DAtVertex(casindex,uuRec_0); //--Impact parametr of Xi
                                                                                                //--at dimuon vertex               
               if( MyPrint_ )                                                                               
                   std::cout<<"Hec: cascade tracks:-->rhoPi:"<<Cas_.rhoPi<<" rhoXi:"<<Cas_.rhoXi
                                                  <<"  IpPi:"<<Cas_.IpPi<<"  Npix:"<<ILam_.NPixelXiTracks 
                                                  <<" IpXiuu:"<<Cas_.IpXiuu<<std::endl;
              
               //--Cascade mass constrains
               XiRec.AddMassConstraint( cascade_mass_c, cascade_mass_c*1.e-6 );       //--Mass constrain Lambda + Pion to cascade   
               if( !XiRec.isValid() ){
                  if( MyPrint_ )
                      std::cout<< "Hec: Cascade mass constraint fit failed!"<<std::endl;
                      continue;   //--goes end of the cascade loop
               }  
               //--Output cascade mass constrain fit variables       
               Cas_.MxiMFit      = XiRec.mass();                                      //--Flag:  valid mass constrain for Cascade                              
               double xiVM_Chi2  = XiRec.chi2();
               double xiVM_Ndof  = XiRec.ndof();
               Cas_.xiVM_Cl      = ChiSquaredProbability( xiVM_Chi2, xiVM_Ndof );     //--Confidence Level
               reco::Particle::LorentzVector xiVM_P4 = XiRec.P4();
               Cas_.xiVM_px = xiVM_P4.px();                                           //--XiRec.P4().x()  is the same ?? yes
               Cas_.xiVM_py = xiVM_P4.py();
               Cas_.xiVM_pz = xiVM_P4.pz();
               if( MyPrint_ )
                   std::cout<< "Hec: Cascade momentum x:"<<xiVM_P4.px()<<" = ? "<<XiRec.P4().x()<<std::endl;
                
               //--Check separation from Dimuon (L over Sigma)
               typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D_uuXi;
               typedef ROOT::Math::SVector<double, 3> SVector3_uuXi;
 
               SMatrixSym3D_uuXi totalCov_uuXi = xiV_Cov + uuV_Cov;                   //--Separacion Dimuon - cascade
               SVector3_uuXi uuXi_PrimarySep3D( uuV_vx.x() - xiV_vx.x(),
                                                uuV_vx.y() - xiV_vx.y(),
                                                uuV_vx.z() - xiV_vx.z() );
 
               Cas_.uuXi_L        = ROOT::Math::Mag( uuXi_PrimarySep3D );
               double uuXi_Sigma  = sqrt(ROOT::Math::Similarity( totalCov_uuXi, uuXi_PrimarySep3D )) / Cas_.uuXi_L;
 
               Cas_.uuXi_LoS = Cas_.uuXi_L / uuXi_Sigma;
               if( MyPrint_ )
                   std::cout<< "Hec: Cascade mass constraint"<<Cas_.MxiMFit<<std::endl;
 
               //--Cascade_b Vertex Fit (this should be in the same point as the Dimuon vertex)
               std::vector<RefCountedKinematicParticle> Xi_bParticle;
               Xi_bParticle.push_back( XiRec.RKParent() );
               Xi_bParticle.push_back( uuRec.RKParent() );
                
               //--Measure distance between tracks at their closest approach               
               ClosestApproachOnHelixLine cApp;
               GlobalPoint crosspoint;
	       std::pair<GlobalPoint, GlobalPoint> points;
               const TrackCharge  q_xi = XiRec.RKParent()->currentState().particleCharge();
               const GlobalVector p_xi = XiRec.RKParent()->currentState().globalMomentum();
               const GlobalPoint  x_xi = XiRec.RKParent()->currentState().globalPosition();
               const GlobalVector p_uu = uuRec.RKParent()->currentState().globalMomentum();
               const GlobalPoint  x_uu = uuRec.RKParent()->currentState().globalPosition();
	       double dca = -999;
               if( cApp.calculate( q_xi, p_xi, x_xi, p_uu, x_uu, 
                                  XiRec.RKParent()->currentState().freeTrajectoryState().parameters().magneticField() ) ){
                  crosspoint=cApp.crossingPoint();
                  points=cApp.points();
                  dca = cApp.distance();
               }
               Cas_.dca  = dca;                         //--Flag:  distance closest approach
               Cas_.xcPt = crosspoint.x();
               Cas_.ycPt = crosspoint.y();
               Cas_.zcPt = crosspoint.z();
               double uu_Cas = uuRec.RKParent()->currentState().globalMomentum().dot(XiRec.RKParent()->currentState().globalMomentum());
               double cos_uu_Cas = uu_Cas /( uuRec.RKParent()->currentState().globalMomentum().mag() *
                                             XiRec.RKParent()->currentState().globalMomentum().mag() );         //--opening angle between dimuon and cascade
               Cas_.VuuVxiCosOp = cos_uu_Cas;                                                                   //--Opening Angle
                                         
               float uuDist = sqrt( uuFit_.uuVx*uuFit_.uuVx + uuFit_.uuVy*uuFit_.uuVy + uuFit_.uuVz*uuFit_.uuVz );
               float xiDist = sqrt( Cas_.xiVx*Cas_.xiVx + Cas_.xiVy*Cas_.xiVy + Cas_.xiVz*Cas_.xiVz );
                
               float uuR = sqrt( uuFit_.uuVx*uuFit_.uuVx + uuFit_.uuVy*uuFit_.uuVy );
               float xiR = sqrt( Cas_.xiVx*Cas_.xiVx + Cas_.xiVy*Cas_.xiVy );
                
               if( MyPrint_ ){
                   std::cout<<" dca = "<<dca<<" CroossPoint ="<<crosspoint<<std::endl;
                   std::cout<<"  r="<<sqrt( crosspoint.x()*crosspoint.x() + crosspoint.y()*crosspoint.y() )<<std::endl;
                   std::cout<<"Xi "<<xiV_Chi2<<" Ndof: "<<xiV_Ndof<<" CL: "<<Cas_.xiV_Cl<<std::endl; 
                   std::cout<<"dimuon: chi**2: "<<uuV_Chi2<<" Ndof: "<<uuV_Ndof<<" CL: "<<uuFit_.uuV_Cl<<std::endl;
                   std::cout<<" mui= "<<IuuC_.imuId<<" muj= "<<IuuC_.jmuId<<std::endl;
                   std::cout<<"Dimuon Vertex:  Vx= "<<uuFit_.uuVx<<" Vy= "<<uuFit_.uuVy<<" Vz= "<<uuFit_.uuVz<<" Dist="<<uuDist<<" R="<<uuR<<std::endl;
                   std::cout<<"Cascade Vertex: Vx= "<<Cas_.xiVx  <<" Vy= "<<Cas_.xiVy  <<" Vz= "<<Cas_.xiVz  <<" Dist="<<xiDist<<" R="<<xiR<<std::endl;
                   std::cout<<"Puu dot  Pxi: "<<uu_Cas<<" "<<cos_uu_Cas<<std::endl;
                   std::cout<<"-----------------------------------"  <<std::endl;
               }

               //--Cut on fiducial tracking volume [-120 cm < r < 120 cm  and -260 cm < Z < 260 cm]   r = sqrt( x^2 + y^2 )
               //--Pixel [-20 cm < r < 20 cm  and -60 cm < Z < 60 cm] to run vertex fit
               if( ( dca >= 0. && dca <= 5 ) && 
                   (sqrt( crosspoint.x()*crosspoint.x() + crosspoint.y()*crosspoint.y() ) < 5.) &&      //--fiducial inside tracker [Aug 24, 2011]
                    std::abs( crosspoint.z() ) < 60. ){
                
               //--Vertex Cascade_b (Cascade + dimuon)   
               KinematicFitDriver Xi_bRec( Xi_bParticle, "Cascade_b" ); 
               
               if( !Xi_bRec.isValid() ){
                 if( MyPrint_ )
                     std::cout<< "Hec: Cascade_b vertex fit failed!"<<std::endl;
                     continue;   //--goes end of the cascade loop
               } 
               //--Output cascade_b Vtx fit variables
               reco::Particle::Point xi_bV_vx;
               reco::Vertex::CovarianceMatrix xi_bV_Cov;
               Casb_.Mxi_bVFit    = Xi_bRec.mass();          //--Flag:  valid Cascade_b vertex
               xi_bV_vx           = Xi_bRec.VertexDecay();   //--Vx,Vy,Vz 
               xi_bV_Cov          = Xi_bRec.VertexCov();                                   
               double xi_bV_Chi2  = Xi_bRec.chi2();
               double xi_bV_Ndof  = Xi_bRec.ndof();
               Casb_.xi_bV_Cl     = ChiSquaredProbability( xi_bV_Chi2, xi_bV_Ndof );  //--Confidence Level cacade Vtx
               Casb_.xi_bVx       = xi_bV_vx.x();
               Casb_.xi_bVy       = xi_bV_vx.y();
               Casb_.xi_bVz       = xi_bV_vx.z();
               
               //--Cascade_b
               reco::Particle::LorentzVector xi_bV_P4 = Xi_bRec.P4();    //--Cascade_b 4-momentum
               Casb_.xi_bV_px = xi_bV_P4.px();                           //--Xi_bRec.P4().x()   Both are the same (sep 9, 11)
               Casb_.xi_bV_py = xi_bV_P4.py();
               Casb_.xi_bV_pz = xi_bV_P4.pz();
               
               
               //--Refit the primary and exclude the Dimuon, Pion and lambda (proton + pion).  5 tracks (exclude Xi and /\0 tracks) 
               std::vector<reco::TrackRef> ExclusionList;
               ExclusionList.push_back(piontrk);	      
               ExclusionList.push_back(pionveetrk);	      
               ExclusionList.push_back(protonveetrk);
               ExclusionList.push_back( iMuTrackRef );
               ExclusionList.push_back( jMuTrackRef ); 
               VertexRefit RefitPrimary(primary.BestVertex() ,ExclusionList ,bFieldHandle ,beamSpot );
               reco::Vertex refitVertexPrim = RefitPrimary.Refitted();
               GlobalPoint PosVertex = GlobalPoint(refitVertexPrim.x(),refitVertexPrim.y(),refitVertexPrim.z());
               // NewPrimaryProb = ChiSquaredProbability( refitVertexPrim.chi2(), refitVertexPrim.ndof());
               
               //--Cascade cleaning  (Eduardo N-Tuplizer)
               const GlobalPoint XiVtxPos( theCascades[casindex].vx(), 
                                           theCascades[casindex].vy(), 
				           theCascades[casindex].vz()  );
               
               GlobalVector XiVectorP( theCascades[casindex].momentum().x(),
                                       theCascades[casindex].momentum().y(),
			               theCascades[casindex].momentum().z() );
                    
               //--Build Xi trajectory state
               FreeTrajectoryState XiFTS =FreeTrajectoryState(XiVtxPos,XiVectorP,piontrk->charge(),BField);
               
               //--use pion track at Xi decay to estimate errors to include in  XiTTrack
               TrajectoryStateClosestToPoint PionTrajAtXiDecay = pionTtrk.trajectoryStateClosestToPoint(XiVtxPos);           
               FreeTrajectoryState piFTS  = PionTrajAtXiDecay.theState();
               XiFTS.setCartesianError(piFTS.cartesianError());

               reco::TransientTrack XiTTrack = (*theTTBHandle).build(XiFTS);
               
               TrajectoryStateClosestToPoint XiTraj = XiTTrack.trajectoryStateClosestToPoint(PosVertex);
               
               double XiD0  = XiTraj.isValid()?XiTraj.perigeeParameters().transverseImpactParameter()
                                   :-1000;
               double XiD0E = XiTraj.hasError()?XiTraj.perigeeError().transverseImpactParameterError()
                                    :-1000;     
               double XiDz  = XiTraj.isValid()?XiTraj.perigeeParameters().longitudinalImpactParameter()
                                   :-1000;
               double XiDzE = XiTraj.hasError()?XiTraj.perigeeError().longitudinalImpactParameterError()
                                    :-1000;     
               double XiD3  = sqrt(XiD0*XiD0 + XiDz*XiDz);
               double XiD3E = sqrt(XiD0E*XiD0E*XiD0*XiD0 + XiDzE*XiDzE*XiDz*XiDz)/XiD3;
               Cas_.Xiip3d = XiD3/XiD3E;
               //-----------------------------------------------------
               //--Check separation Xi_b from primary (L over Sigma)
               typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D_Xi_b;
               typedef ROOT::Math::SVector<double, 3> SVector3_Xi_b;

               SMatrixSym3D_Xi_b totalCov_xi_b = refitPrimVertex.covariance() + xi_bV_Cov;     //--separation cascade_b to primary(without dimuons)
               SVector3_Xi_b xi_bPrimarySep3D( xi_bV_vx.x() - refitPrimVertex.x(),
                                               xi_bV_vx.y() - refitPrimVertex.y(),
                                               xi_bV_vx.z() - refitPrimVertex.z() );

               Casb_.xi_bL       = ROOT::Math::Mag( xi_bPrimarySep3D );
               double xi_bSigma  = sqrt(ROOT::Math::Similarity( totalCov_xi_b, xi_bPrimarySep3D )) / Casb_.xi_bL;

               Casb_.xi_bLoS = Casb_.xi_bL / xi_bSigma;
               
               //--Check separation Xi_b from dimuon (L over Sigma)
               typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D_uuXi_b;
               typedef ROOT::Math::SVector<double, 3> SVector3_uuXi_b;

               SMatrixSym3D_uuXi_b totalCov_uuXi_b = xi_bV_Cov + uuV_Cov;                      //--separation dimuon -cascade_b
               SVector3_uuXi_b uuXi_bPrimarySep3D( uuV_vx.x() - xi_bV_vx.x(),
                                                   uuV_vx.y() - xi_bV_vx.y(),
                                                   uuV_vx.z() - xi_bV_vx.z() );

               Casb_.uuxi_bL       = ROOT::Math::Mag( uuXi_bPrimarySep3D );
               double uuxi_bSigma  = sqrt(ROOT::Math::Similarity( totalCov_uuXi_b, uuXi_bPrimarySep3D )) / Casb_.uuxi_bL;

               Casb_.uuxi_bLoS = Casb_.uuxi_bL / uuxi_bSigma;
               if( MyPrint_ )
                   std::cout<< "Hec: Cascade_b Vtx"<<Casb_.Mxi_bVFit<<std::endl;
            
               //--Here Pointing Secondary Vertex  (dimuon + cascade) to Primary  [Sep 9, 2011]
               GlobalVector xib_PrimSec( xi_bV_vx.x() - refitPrimVertex.x(),
                                         xi_bV_vx.y() - refitPrimVertex.y(),
                                         xi_bV_vx.z() - refitPrimVertex.z() );                           
               GlobalVector xib_P3( xi_bV_P4.px(),
                                    xi_bV_P4.py(),
                                    xi_bV_P4.pz() );                              
               double V_dot_Pxib =xib_P3.dot( xib_PrimSec );
                                                                                 
               double cos_alphab = V_dot_Pxib /( xib_PrimSec.mag() * xib_P3.mag() );
               Casb_.VxibVpCos = cos_alphab;
               if( MyPrint_ ){
                   std::cout<<"Hec: Pointing Xi_b cos alpha= "<<Casb_.VxibVpCos<<std::endl;                                                                                                              
                   std::cout<<"Hec: Last cascade tracks:-->rhoPi:"<<Cas_.rhoPi<<" rhoXi:"<<Cas_.rhoXi
                                                       <<"  IpPi:"<<Cas_.IpPi<<"  Npix:"<<ILam_.NPixelXiTracks 
                                                      <<" IpXiuu:"<<Cas_.IpXiuu<<std::endl;
               }
               //--Fill Ntuple and re-initialize all variables
               ILam_.uuPair++;
               //fill_uuC();
               //fill_uuFit();
               //fill_MyCas;
               hecmu_tree_->Fill();
               
               }  //--Cascade-Dimuon Fiducial vertex volume
            }     //--for( unsigned int casindex
          }       //--End if( !onlyDiMu_ ){
          
          }     //--End:  if( !uuRec.isValid() )  else              //--dimuon mass constrain
          }     //--End:  if( uuMcut )                              //--diMuon Mass Cut for mass constrain
          }     //--End:  if( !refitPrimVertex.isValid() ) else     //--primary vertex without dimuon
          }     //--End:  if ( !uuRec.isValid() ) else              //--Dimuon Vertex
          }     //--End:  if( (dca_uu > 0. && dca < 5 ) && ((sqrt(  //--Geometrical Vertex constrains Fiducial
        }       //--End:  if( uuC_.GLi!=0 && uuC_.GLj!=0 )          //--Skip standalone muon
         
        //--Fill Ntuple and re-initialize all variables
        if( onlyDiMu_ ){
          //fill_uuC();
          //fill_uuFit();
          hecmu_tree_->Fill();
        }
         
     }   //--for(reco::MuonCollection::const_iterator jMuon = iMuon+1;
  }      //--for(reco::MuonCollection::const_iterator iMuon = allmuons->begin();
  if( !onlyDiMu_ ){
    //--Check out here Lambda  (Lambda --> Proton + Pi-)
    for( unsigned int veeindex = 0; veeindex < theVees.size(); veeindex++ ) {
       double recLambda_mass = theVees[veeindex].mass();
       histo_lam01->Fill( recLambda_mass );
    }  //--for( unsigned int casindex
      
    //--Check out here Cascade (Xi- --> Lambda + Pi- )
    for( unsigned int casindex = 0; casindex < theCascades.size(); casindex++ ) {
       double recXi_mass = theCascades[casindex].mass();
       histo_cas01->Fill( recXi_mass );
    }  //--for( unsigned int casindex
  }    //--End if( !onlyDiMu_ )
  
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
HecPsiCascade::beginJob()
{
  std::cout << " HecPsiCascade::beginJob" << std::endl;
  hecmu_tree_ = new TTree("DiMuCascade","Dimuon Cascade Ntuple");
  int bufsize = 64000;
  
  hecmu_tree_->Branch("Dimuon", &uuC_,"iEta/D:iPhi:iPt:iP:id0:idz:iCal:iSeg:iIso:jEta:jPhi:jPt:jP:jd0:jdz:jCal:jSeg:jIso:Muu:dca:xcpt:ycpt:zcpt", bufsize);
  hecmu_tree_->Branch("uuFits", &uuFit_,"MuuFit/D:MuuMFit:uuLoS:uuL:uuV_Cl:pV_Cl:uuVM_Cl:uuVx:uuVy:uuVz:xpV:ypV:zpV:uuVM_px:uuVM_py:uuVM_pz:CosVuuVp", bufsize);

  if( !onlyDiMu_ ){
    hecmu_tree_->Branch("Lambda", &Lam_,"xiM/D:lamM:uuXiM_0:uuXiM_1:uuXiM_2:MlamVFit:lamV_Cl:lamVx:lamVy:lamVz:MlamMFit:lamVM_Cl:lamVM_px:lamVM_py:lamVM_pz:uuL_L:uuL_LoS", bufsize);
    
    hecmu_tree_->Branch("Cascade", &Cas_,"dca/D:xcPt:ycPt:zcPt:MxiVFit:xiV_Cl:xiVx:xiVy:xiVz:MxiMFit:xiVM_Cl:xiVM_px:xiVM_py:xiVM_pz:uuXi_L:uuXi_LoS:Xiip3d:VuuVxiCosOp:XiCL:rhoPi:rhoXi:IpPi:IpXiuu", bufsize);
    
    hecmu_tree_->Branch("Cascade_b", &Casb_,"Mxi_bVFit/D:xi_bV_Cl:xi_bVx:xi_bVy:xi_bVz:xi_bV_px:xi_bV_py:xi_bV_pz:xi_bL:xi_bLoS:uuxi_bL:uuxi_bLoS:VxibVpCos", bufsize); 
  }  
  hecmu_tree_->Branch("iHeader", &Ievt_,  "run/I:evtnum:lumib:Ntk:allMu:allCalMu:allCascade:allLambda:allPrim", bufsize);
  hecmu_tree_->Branch("iDimuon", &IuuC_,"iGL/I:iCharge:imuId:iNtrkHits:jGL:jCharge:jmuId:jNtrkHits", bufsize);
  if( !onlyDiMu_ )
    hecmu_tree_->Branch("iLambda", &ILam_,"uuPair/I:xiCan:NPixelXiTracks", bufsize); 
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HecPsiCascade::endJob() {
  std::cout << " HecPsiCascade::endJob" << std::endl;
}
void HecPsiCascade::init()
{
  Ievt_.init();
  IuuC_.init();
  uuC_.init();
  uuFit_.init();
  ILam_.init();
  Lam_.init();
  Cas_.init();
  Casb_.init();
}
void HecPsiCascade::Ievt::init()
{
  int          dummy_int = -1;
  runNb      = dummy_int;
  eventNb    = dummy_int;
  lumiBlock  = dummy_int;
  Ntk        = dummy_int;
  allMu      = dummy_int;
  allCalMu   = dummy_int;
  allCascade = dummy_int;
  allLambda  = dummy_int;
}
void HecPsiCascade::fill_evt(const edm::Event& iEvent,
     int Ntk, int allMu, int allCalMu, int allCascade, int allLambda, int allPrim ){
     init();
     Ievt_.runNb      = iEvent.id().run(); 
     Ievt_.eventNb    = iEvent.id().event();
     Ievt_.lumiBlock  = iEvent.luminosityBlock();
     Ievt_.Ntk        = Ntk;            //--tracks->size();
     Ievt_.allMu      = allMu;          //--->allmuons->size();
     Ievt_.allCalMu   = allCalMu;       //--allcalmuons->size();
     Ievt_.allCascade = allCascade;     //--theCascades.size();
     Ievt_.allLambda  = allLambda;      //--theVees.size();
     Ievt_.allPrim    = allPrim;
     if( MyPrint_ )
       std::cout<<"Hec: Run: "<<Ievt_.runNb<<" Evt: "<<Ievt_.eventNb<<" Lumi:"<<Ievt_.lumiBlock
                <<" Ntk: "<<Ievt_.Ntk<<" Mu:"<<Ievt_.allMu<<" CalMu:"<<Ievt_.allCalMu
                <<" NXi: "<<Ievt_.allCascade<<" NLambda: "<<Ievt_.allLambda 
                <<" NPrimVtx:"<<Ievt_.allPrim
                <<std::endl;
}
void HecPsiCascade::IuuC::init()
{
  int         dummy_int   = -1;
  iGL       = dummy_int;
  iCharge   = 0;
  imuId     = dummy_int;
  iNtrkHits = dummy_int;
  jGL       = dummy_int;
  jCharge   = 0;
  jmuId     = dummy_int;
  jNtrkHits = dummy_int;
}
void HecPsiCascade::uuC::init()
{
  double    dummy_float = -1.0;
  iEta    = dummy_float;
  iPhi    = dummy_float;
  iPt     = dummy_float;
  iP      = dummy_float;
  id0     = dummy_float;
  idz     = dummy_float;
  iCal    = dummy_float;
  iSeg    = dummy_float;
  iIso    = dummy_float;
  jEta    = dummy_float;
  jPhi    = dummy_float;
  jPt     = dummy_float;
  jP      = dummy_float;
  jd0     = dummy_float;
  jdz     = dummy_float;
  jCal    = dummy_float;
  jSeg    = dummy_float;
  jIso    = dummy_float;
  Muu     = dummy_float;
  dca     = dummy_float;
  xcpt    = dummy_float;
  ycpt    = dummy_float;
  zcpt    = dummy_float;
}
void HecPsiCascade::fill_uuC()
{
//--Already set in the Analyzer   
// iGL,iCharge,imuId;
// jGL,jCharge,jmuId;
// iEta,iPt,iP;
// jEta,jPt,jP;
// double Muu
}  
void HecPsiCascade::uuFit::init()
{
  double    dummy_float = -1.0;
  MuuFit  = dummy_float;
  MuuMFit = dummy_float;
  uuLoS   = dummy_float;
  uuL     = dummy_float;
  uuV_Cl  = dummy_float;
  pV_Cl   = dummy_float;
  uuVM_Cl = dummy_float;
  uuVx    = dummy_float;
  uuVy    = dummy_float;
  uuVz    = dummy_float;
  xpV     = dummy_float;
  ypV     = dummy_float;
  zpV     = dummy_float;
  uuVM_px = dummy_float;
  uuVM_py = dummy_float;
  uuVM_pz = dummy_float;
  VuuVpCos= -999;
}
void HecPsiCascade::fill_uuFit()
{
  //--Already set in the Analyzer 
  //--uuFit_.MuuFit  = refitMass;     //--uu Vtx reFitted mass (=-1 fit failed)
  //--uuFit_.MuuMFit = refitMass_mc;  //--uu Vtx+MassConstrain reFitted mass (=-1 fit failed)
  //--uuFit_.uuLoS   = uuLoS;         //--uu L over Sigma to Primary
  //--uuFit_.uuL     = uuL;           //--L distance between Prim and DiMuon
  //--uuFit_.uuV_Cl  = uu _Cl;        //--Confidence Level Dimuon Vtx
  //--uuFit_.pV_Cl   = pV_Cl;         //--Confidence Level Primary Vtx without Dimuon
  //--uuVM_Cl =                           Confidence Level Vtx & Mass constrain for Dimuon   
  //--uuVx    =                           Vx for Dimuon Vtx
  //--uuVy    =                           Vy
  //--uuVz    =                           Vz
  //--xpV     =                           Primary Vx
  //--ypV     =                           Primary Vy
  //--zpV     =                           Primary Vz
  //--uuVM_px =                           Px for dimuon with Vtx & Mass constrain
  //--uuVM_py =                           Py
  //--uuVM_pz =                           Pz
  //--VuuVpCos=                           pointing back sec to prim
  
}
void HecPsiCascade::ILam::init()
{
  int         dummy_int   = -1;
  uuPair    = dummy_int; 
  xiCan     = dummy_int;
  NPixelXiTracks= dummy_int;
}
void HecPsiCascade::Lam::init()
{
  double      dummy_float = -1.0;
  xiM       = dummy_float;
  lamM      = dummy_float; 
  uuXiM_0   = dummy_float; 
  uuXiM_1   = dummy_float; 
  uuXiM_2   = dummy_float;
  MlamVFit  = dummy_float; 
  lamV_Cl   = dummy_float;  
  lamVx     = dummy_float;    
  lamVy     = dummy_float;    
  lamVz     = dummy_float;
  MlamMFit  = dummy_float; 
  lamVM_Cl  = dummy_float; 
  lamVM_px  = dummy_float; 
  lamVM_py  = dummy_float; 
  lamVM_pz  = dummy_float;
  uuL_L     = dummy_float;
  uuL_LoS   = dummy_float;
}
void HecPsiCascade::Cas::init()
{
  double      dummy_float = -1.0;
  dca       = dummy_float;
  xcPt      = dummy_float;
  ycPt      = dummy_float;
  zcPt      = dummy_float;
  MxiVFit   = dummy_float; 
  xiV_Cl    = dummy_float;  
  xiVx      = dummy_float;    
  xiVy      = dummy_float;    
  xiVz      = dummy_float;
  MxiMFit   = dummy_float; 
  xiVM_Cl   = dummy_float; 
  xiVM_px   = dummy_float; 
  xiVM_py   = dummy_float; 
  xiVM_pz   = dummy_float;
  uuXi_L    = dummy_float; 
  uuXi_LoS  = dummy_float;
  Xiip3d     = 999;
  VuuVxiCosOp = -999;
  XiCL      = -1;
  rhoPi     =-999;
  rhoXi     =-999;
  IpPi      =-999;
  IpXiuu    =-999;
}
void HecPsiCascade::Casb::init()
{
  double      dummy_float = -1.0;
  Mxi_bVFit = dummy_float; 
  xi_bV_Cl  = dummy_float; 
  xi_bVx    = dummy_float; 
  xi_bVy    = dummy_float; 
  xi_bVz    = dummy_float;
  xi_bV_px  = dummy_float; 
  xi_bV_py  = dummy_float; 
  xi_bV_pz  = dummy_float;
  xi_bL     = dummy_float; 
  xi_bLoS   = dummy_float; 
  uuxi_bL   = dummy_float; 
  uuxi_bLoS = dummy_float;
  VxibVpCos = -999;
} 
//define this as a plug-in
DEFINE_FWK_MODULE(HecPsiCascade);
