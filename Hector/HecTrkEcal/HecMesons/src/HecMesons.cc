// -*- C++ -*-
//
// Package:    HecMesons
// Class:      HecMesons
// 
/**\class HecMesons HecMesons.cc HecTrkEcal/HecMesons/src/HecMesons.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Hector Mendez
//         Created:  Wed Jun 30 15:47:55 CDT 2010
// $Id: HecMesons.cc,v 1.1 2010/08/12 18:17:06 mendez Exp $
//
//


// system include files
#include <memory>
#include <string>

// user include files
#include "HecTrkEcal/HecMesons/interface/HecMesons.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

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

//#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//

// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
HecMesons::HecMesons(const edm::ParameterSet& iConfig)
:
 trackTags_(iConfig.getUntrackedParameter<edm::InputTag>("tracks")),
 electronCollection_(iConfig.getUntrackedParameter<edm::InputTag>("gsfElectrons")),
 minTracks_(iConfig.getUntrackedParameter<unsigned int>("minTracks",0)),
 check_elC_(iConfig.getUntrackedParameter<bool>("check_elC",true)),
 Pe_match_(iConfig.getUntrackedParameter<double>("Pe_match",1000000.0)),
 mass_uno_(iConfig.getUntrackedParameter<double>("mass_uno",0))
{
   //now do what ever initialization is needed
  edm::Service<TFileService> fs;
  histo_trk  = fs->make<TH1D>("Ntrk"       , "Ntrk"       , 200 , 0 , 200 );
  histo_Pti  = fs->make<TH1D>("Pti"        , "Pt All"     , 200 , 0 , 100 );
  histo_ij   = fs->make<TH1D>("Mass ij"    , "Mass ij All", 400 , 0 , 200 );
  histo_ijPt = fs->make<TH1D>("Mass ijPt"  , "Mass ij Pt" , 400 , 0 , 200 );
  
  histo_ij0  = fs->make<TH1D>("Mass ij0"   , "Mass ij 0e" , 400 , 0 , 200 );
  histo_ij1  = fs->make<TH1D>("Mass ij1"   , "Mass ij 1e" , 400 , 0 , 200 );
  histo_ij2  = fs->make<TH1D>("Mass ij2"   , "Mass ij 2e" , 400 , 0 , 200 );
  histo_elec = fs->make<TH1D>("Nelec"      , "Nelec"      , 100 , 0 , 100 );
  histo_NeTrk= fs->make<TH1D>("NeTrk"      , "NeTrk"      , 100 , 0 , 100 );
  histo_DNel = fs->make<TH1D>("DNel"       , "Delta Nel"  ,  10 ,- 5,   5 );
  histo_Delp = fs->make<TH1D>("Delp"       , "Delta el p" , 100 ,-20,  80 );
  histo_DelpC= fs->make<TH1D>("DelpC"      , "DeltaP CUT" , 100 ,-20,  80 );
  histo_Delq = fs->make<TH1D>("Delq"       , "Delta el q" ,   6 , -3,   3 );
  histo_eiej = fs->make<TH1D>("Mass eieja" , "Mass eiej"  , 400 , 0 , 200 );
  histo_DM   = fs->make<TH1D>("DM"         , "Delta Mass" , 100 ,-1/5000,1/5000 );
  histo_eEn  = fs->make<TH1D>("elec Energy", "elec Energy", 800 , 0 , 200 );
  histo_EoP  = fs->make<TH1D>("EoP"        , "elec E/P"   , 100 , 0 , 5 );
  histo_Clu  = fs->make<TH1D>("ph Cluster" , "Nclus"      , 100 , 0 , 100 );
  histo_Eni  = fs->make<TH1D>("ph Energy"  , "Energy i"   , 800 , 0 , 200 );

}


HecMesons::~HecMesons()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
HecMesons::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //using namespace edm;
  
  //--Get track Collection
  //using reco::TrackCollection;
  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByLabel(trackTags_,tracks); 
  //const reco::TrackCollection* trackCollection = tracks.product();
  
  //--Get electron Collection [from: CMSSW/RecoEgamma/Examples/plugins/GsfElectronDataAnalyzer.cc]
  edm::Handle<reco::GsfElectronCollection> gsfElectrons;
  iEvent.getByLabel(electronCollection_,gsfElectrons);


  //--Loop over rec electrons
  histo_elec->Fill( gsfElectrons->size() );
  for(reco::GsfElectronCollection::const_iterator gsfIter=gsfElectrons->begin(); gsfIter!=gsfElectrons->end(); gsfIter++){     
     double iCharge = gsfIter->charge();
     
     histo_eEn->Fill( gsfIter->energy() );
     histo_EoP->Fill( gsfIter->eSuperClusterOverP() );
     
     //--gsfIter->overlap();
     for(reco::GsfElectronCollection::const_iterator gsfJter=gsfIter+1; gsfJter!=gsfElectrons->end(); gsfJter++){     
        double jCharge = gsfJter->charge();
        if( iCharge * jCharge >=0 )continue;
        
        //--Di-electron Invariant Mass
        double Masseiej = ( gsfIter->p4() +  gsfJter->p4() ).M();
        histo_eiej->Fill( Masseiej );
        
        //--The next 20 lines are equivalent to the 2 previous lines to calculate di-electron invariant mass
        double mel = 0.00051099891;
        
        double Pexi = gsfIter->px();
        double Peyi = gsfIter->py();
        double Pezi = gsfIter->pz();
      //double Peti = gsfIter->pt();   //--sqrt( Pexi*Pexi + Peyi*Peyi )                 
        double Enei = sqrt( mel*mel + Pexi*Pexi + Peyi*Peyi + Pezi*Pezi );
        
        double Pexj = gsfJter->px();
        double Peyj = gsfJter->py();
        double Pezj = gsfJter->pz();
        double Enej  = sqrt( mel*mel + Pexj*Pexj + Peyj*Peyj + Pezj*Pezj );
        
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
        
        /*
        if( DeltaM != 0 )
          std::cout<<"Hec: Wrong Di-electron Mass"
                   <<"Event " << iEvent.id() <<" Delta M"<< DeltaM << std::endl;      
        */
     }    //--electrons gsfJter
  }       //--electrons gsfIter
  
  //--Match electrons [from: CMSSW/Validation/RecoMET/src/METTester.cc]
  int nEls = 0;
#define TOT_e 500
  unsigned int elecHec[TOT_e]={0};
  double elec_q[TOT_e]={0};
  double elec_p[TOT_e]={0};
//for(int i=0; i<=TOT_e; i++ )elecHec[i] = 0;  
//for(int i=0; i<=TOT_e; i++ )if( elecHec[i] != 0 )std::cout<<" Hec: NO CLEAN ARRAY: "<< i <<" "<< elecHec[i] << std::endl;

//--from:  CMSSW/ RecoEgamma/ EgammaTools/ src/ ConversionFinder.cc
//--   const TrackCollection *ctftracks = track_h.product();
//--   const reco::TrackRef el_ctftrack = gsfElectron.closestCtfTrackRef();
//--   int ctfidx = -999.;
//--   if(el_ctftrack.isNonnull() && gsfElectron.shFracInnerHits() > minFracSharedHits)
//--     ctfidx = static_cast<int>(el_ctftrack.key());
 
//for( edm::View<reco::GsfElectron>::const_iterator eleit = gsfElectrons->begin(); eleit != gsfElectrons->end(); eleit++ ) {
  for( reco::GsfElectronCollection::const_iterator eleit = gsfElectrons->begin(); eleit != gsfElectrons->end(); eleit++ ) {
    //me["helePt"]->Fill( eleit->p4().pt() );  
    //me["heleEta"]->Fill( eleit->p4().eta() );
    //me["heleHoE"]->Fill( eleit->hadronicOverEm() );

      reco::TrackRef el_track = eleit->closestCtfTrackRef();

      unsigned int ele_idx = el_track.isNonnull() ? el_track.key() : 99999;  //ele_idx = 99999 if el_track is null (or bad)      [False]
                                                                             //ele_idx = Trk Number pointing to the generalTrack [True]
    //if( eleit->hadronicOverEm()  < 0.1  && ele_idx < tracks->size() ){
      if( eleit->shFracInnerHits() > 0.45 && ele_idx < tracks->size() ){
    //if( ele_idx < tracks->size() ){
        ++nEls;
        elecHec[ele_idx] = 1;
        elec_q[ele_idx] = eleit->charge();
        elec_p[ele_idx] = eleit->p();
        if( check_elC_ )
          std::cout<<" -Number of elec: "<< nEls
                   <<" trk Number:" << ele_idx 
                   <<" Q:"<< elec_q[ele_idx]
                   <<" P:"<< elec_p[ele_idx]
                   <<" Had/Em Fraction:"   << eleit->hadronicOverEm()
                   <<" HitsInner Fraction:"<< eleit->shFracInnerHits()
                   << std::endl;
      }
  }
  histo_NeTrk->Fill( nEls );  
  for(unsigned int i=0; i<=tracks->size(); i++ )
     if( check_elC_ && elecHec[i] != 0 )
       std::cout<<nEls<<" eTrk: "<< i <<" "<< elecHec[i] << std::endl;
  if( check_elC_ )
    std::cout<<"Delta_elec:"<< gsfElectrons->size() - nEls << std::endl;
  histo_DNel->Fill( gsfElectrons->size() - nEls );
  
  //--Get track Collection
  if( minTracks_ <= tracks->size() ) {
      edm::LogInfo("Demo")<<"Hec: Event " << iEvent.id() << " Number of Tracks "<<tracks->size();
      std::cout<<"Event " << iEvent.id() 
               << " ----Number of Rec Tracks: "<<tracks->size() 
               <<" Input mass:"<<mass_uno_<<" GeV" <<" Pe_match:"<<Pe_match_<<" GeV"
               << std::endl;
  }
  histo_trk->Fill( tracks->size() );
  unsigned int Ntki = 0;
  unsigned int Ntkj = 0;
  unsigned int Ntke = 0;
  for(reco::TrackCollection::const_iterator itTrack = tracks->begin(); itTrack != tracks->end(); ++itTrack){
     int iCharge = itTrack->charge();
     
     double Pti = itTrack->pt();
     histo_Pti->Fill( Pti );
     
     double Pki  = itTrack->p();
     double Delta_Pi = 9999.0;
     unsigned e_it = 0;
     if( elecHec[Ntki] != 0 ){
       e_it = 1;
       ++Ntke;
       Delta_Pi = elec_p[Ntki] - Pki;
       histo_Delp->Fill( Delta_Pi );
       if( abs(Delta_Pi) <= Pe_match_ )histo_DelpC->Fill( Delta_Pi );
       histo_Delq->Fill( elec_q[Ntki] - iCharge );
       if( check_elC_ )
         std::cout<<"->Number of elec: "<< Ntke <<" trk Number:" << Ntki    <<" Q:"<< iCharge <<" P:"<< Pki << std::endl;
     }
     ++Ntki;
     Ntkj = Ntki;
     
     for(reco::TrackCollection::const_iterator jtTrack = itTrack+1; jtTrack != tracks->end(); ++jtTrack){
        int jCharge = jtTrack->charge();
        double Pkj = jtTrack->p();
        unsigned e_jt = 0;
        double Delta_Pj = 9999.0;
        if( elecHec[Ntkj] != 0 ){
          e_jt = 1;
          Delta_Pj = elec_p[Ntkj] - Pkj;
          if( check_elC_ )
            std::cout<<Ntki<<" ->elec trk Number:" << Ntkj    <<" Q:"<< jCharge <<" P:"<< Pkj << std::endl;
        }
        ++Ntkj;
        
        if( iCharge * jCharge >=0 )continue;
         
        double mka = mass_uno_; 
         
        double Pkxi = itTrack->px();
        double Pkyi = itTrack->py();
        double Pkzi = itTrack->pz();
        double Pkti = itTrack->pt();   //--sqrt( Pxi*Pxi + Pyi*Pyi )
        double Enki = sqrt( mka*mka + Pkxi*Pkxi + Pkyi*Pkyi + Pkzi*Pkzi );
         
        double Pkxj = jtTrack->px();
        double Pkyj = jtTrack->py();
        double Pkzj = jtTrack->pz();
        double Pktj = jtTrack->pt();   //--sqrt( Pxj*Pxj + Pyj*Pyj )
        double Enkj = sqrt( mka*mka + Pkxj*Pkxj + Pkyj*Pkyj + Pkzj*Pkzj );
         
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
        
        //--Look here electron Matching (generatTracks vs gsfElectrons). where is the Z peak ?
        if(  e_it == 0 && e_jt == 0 )
          histo_ij0->Fill( Massij );    //--found zero e
          
        if( (e_it == 1 && abs(Delta_Pi)<= Pe_match_ && e_jt == 0)||(e_it == 0 && e_jt == 1 && abs(Delta_Pj) <= Pe_match_) )
          histo_ij1->Fill( Massij );    //--found one  e
        
        if( e_it == 1 && abs(Delta_Pi) <= Pe_match_  && e_jt == 1 && abs(Delta_Pj) <= Pe_match_ )
          histo_ij2->Fill( Massij );    //--found both e
                 
     }     //--end for(TrackCollection::const_iterator jtTrack = itTrack+1; jtTrack != tracks->end(); ++jtTrack)         Kaon j
  }        //--end for(TrackCollection::const_iterator itTrack = tracks->begin(); itTrack != tracks->end(); ++itTrack)   Kaon i


  //--This part on ECal is taken from Carley:
  //--[ /uscms_data/d2/carley/photon/CMSSW_3_5_6/src/RecoEgamma/PhotonIdentification ]
  //--and [ https://twiki.cern.ch/twiki/bin/viewauth/CMS/VGammaFirstPaper ]
  
//   //--Get Ecal and do Pi0 from rec hits
    edm::Handle<EcalRecHitCollection> ecalhitsCollH;
    iEvent.getByLabel("ecalRecHit","EcalRecHitsEB", ecalhitsCollH);
  //iEvent.getByLabel("reducedEcalRecHitsEB","", ecalhitsCollH);
  //const EcalRecHitCollection* rechitsCollection_ = ecalhitsCollH.product();


   //Hm.  What are these?  You'll see in a minute:  some cleverness
   //is required to make sure we know what photon owns what crystal.     
   //--std::map<DetId,int> crysclus;
   
     edm::Handle<reco::BasicClusterCollection> bcHandle;
     iEvent.getByLabel(edm::InputTag("islandBasicClusters:islandBarrelBasicClusters"), bcHandle);
    
     int nClus=0;
     float clusX[1000];
     float clusY[1000];
     float clusZ[1000];
     float clusE[1000];
     float clusEta[1000];
     float clusPhi[1000];
         
     const reco::BasicClusterCollection basicHandle = *bcHandle;
     for(reco::BasicClusterCollection::const_iterator b_iter=basicHandle.begin(); b_iter!=basicHandle.end();++b_iter){
    
         //--std::vector<std::pair<DetId, float> > clusdet = b_iter->hitsAndFractions();
         //--for (int ik=0;ik<int(clusdet.size());++ik){
         //-- crysclus.insert(std::make_pair(clusdet[ik].first, nClus));
         //--}   I don't understand this loop make_pair  Hec  
    
         clusX[nClus]   = b_iter->position().x();
         clusY[nClus]   = b_iter->position().y();
         clusZ[nClus]   = b_iter->position().z();
         clusEta[nClus] = b_iter->position().eta();
         clusPhi[nClus] = b_iter->position().phi();
         clusE[nClus]   = b_iter->energy();
         nClus++; 
         
         histo_Eni->Fill( clusE[nClus] );    
     }
     //std::cout<<" ----Clusters:"<<nClus<<std::endl;
     histo_Clu->Fill( nClus ); 
     if( minTracks_ <= tracks->size() || check_elC_ )
       std::cout<<" ------------------------------------"<< std::endl;
       
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
HecMesons::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HecMesons::endJob() { 

  std::cout<<" -----------Finito----------------"<< std::endl;
}

//define this as a plug-in
DEFINE_FWK_MODULE(HecMesons);
