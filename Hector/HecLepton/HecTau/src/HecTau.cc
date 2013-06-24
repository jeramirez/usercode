	// -*- C++ -*-
//
// Package:    HecTau
// Class:      HecTau
// 
/**\class HecTau HecTau.cc HecLepton/HecTau/src/HecTau.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  carlos malca
//         Created:  Wed Feb  6 15:17:19 CST 2013
// $Id: HecTau.cc,v 1.3 2013/03/15 19:16:02 mendez Exp $
//
//


// system include files
#include <memory>

// user include files
#include "HecLepton/HecTau/interface/HecTau.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

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
#include "Analyzers/CascadeProducer/interface/KinematicFitDriver.h"
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
//

//
// static data member definitions
//

//
// constructors and destructor
//
HecTau::HecTau(const edm::ParameterSet& iConfig)
:
 trackTags_         (iConfig.getParameter<edm::InputTag>("tracks"    )),
 theMuonsLabel_     (iConfig.getParameter<edm::InputTag>("MuonsLabel")),
 MyPrint_           (iConfig.getUntrackedParameter<bool>("MyPrint","False")),
 MyOneKaon_         (iConfig.getUntrackedParameter<bool>("MyOneKaon","False")),
 MyThreeMu_         (iConfig.getUntrackedParameter<bool>("MyThreeMu","False")),
 hlTriggerResults_  (iConfig.getUntrackedParameter<std::string>("HLTriggerResults",std::string("TriggerResults::HLT")) )
{
   //now do what ever initialization is needed
  edm::Service<TFileService> fs;
  //--Diagnostic and testing Histograms
  histo_trk      = fs->make<TH1D>("Ntrk"  , "Ntrk"       , 200 , 0 , 200 );
  histo_nmu      = fs->make<TH1D>("NMuons", "NMuons"     , 200 , 0 , 200 );
  histo_muId     = fs->make<TH1D>("muId"  , "Muon Id"    , 40  , 0 , 40  );

  int Mnch=4000;
  double Mxi=0., Mxf=200.;
  histo_dimu     = fs->make<TH1D>("MassDiMu"    , "Mass DiMu All"   , Mnch , Mxi , Mxf );         //--50 MeV/channel
  histo_dimuk[0] = fs->make<TH1D>("MassDimuK0"  , "Mass DiMu Kaon 0", Mnch , Mxi , Mxf );
  histo_dimuk[1] = fs->make<TH1D>("MassDimuK1"  , "Mass DiMu Kaon 1", Mnch , Mxi , Mxf );
  histo_dimuk[2] = fs->make<TH1D>("MassDimuK2"  , "Mass DiMu Kaon 2", Mnch , Mxi , Mxf );
  histo_dimuk[3] = fs->make<TH1D>("MassDimuK3"  , "Mass DiMu Kaon 3", Mnch , Mxi , Mxf );
  histo_dimuk[4] = fs->make<TH1D>("MassDimuK4"  , "Mass DiMu Kaon 4", Mnch , Mxi , Mxf );

  histo_uuu[0]   = fs->make<TH1D>("uuu0"  , "Mass Three Muons 0", Mnch , Mxi , Mxf );
  histo_uuu[1]   = fs->make<TH1D>("uuu1"  , "Mass Three Muons 1", Mnch , Mxi , Mxf );
  histo_uuu[2]   = fs->make<TH1D>("uuu2"  , "Mass Three Muons 2", Mnch , Mxi , Mxf );

  histo          = fs->make<TH1D>("Ptk"      , "Pt Kaon"             , 100 , 0 , 10 );
  histo_DeltaPku = fs->make<TH1D>("DeltaPku" , "P KaonMuon Matching" , 82 ,-20.5 ,20.5 );
  histo_DeltaQ   = fs->make<TH1D>("DeltaQ"   , "q KaonMuon Matching" , 5 ,-2.5 ,2.5 );

  histo_pV_ClHM  = fs->make<TH1D>(" primCL HM  " , " primCL HM  "  , 1000, 0, 100);
}


HecTau::~HecTau()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
HecTau::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
std::cout<<"-------------------------Carlos Analyzer-------------------- "<<std::endl;
//--Tracks Collection
 edm::Handle< std::vector<pat::GenericParticle> >tracks;    //--thePATTrackHandle;
 iEvent.getByLabel(trackTags_,tracks); int Ntk = tracks->size();


 //--Get Muon Collection
 edm::Handle< std::vector<pat::Muon> >allmuons;
 iEvent.getByLabel(theMuonsLabel_, allmuons); int allMu = allmuons->size();

 //--Handles for Primary Vertex  (from Eduardo)
 PrimaryInfo primary(iEvent,iConfig);

 //--Handle fron Primary just to get the number of rec Prim vertex (Sep 10, 11)
 edm::Handle<reco::VertexCollection> privtxs;
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

 //--Loop over tracks to convert to transient tracks (from Eduardo CascadeFitter.cc)
 std::vector<reco::TrackRef> theTrackRefs;
 std::vector<reco::TransientTrack> theTransTracks;
 for(std::vector<pat::GenericParticle>::const_iterator itrack  = tracks->begin();itrack != tracks->end();   ++itrack){
     reco::TrackRef tmpRef = itrack->track() ;
     reco::TransientTrack tmpTk2( *tmpRef, &(*bFieldHandle), globTkGeomHandle );
     theTrackRefs.push_back( tmpRef );
     theTransTracks.push_back( tmpTk2 );
 }

 //--Looks on Muons  [July 25, 2011]
 histo_trk->Fill( Ntk);
 histo_nmu->Fill( allMu );

 //--HMoreno 05-06-12----------------------------------------------------------------

 init();     //--clean all structures

 HecHltTrig(iEvent);   //-->Unpack Trigger

 for(std::vector<pat::Muon>::const_iterator iMu =allmuons->begin();
                                            iMu!=allmuons->end(); ++iMu){
    reco::TrackRef iMuTrackRef = iMu->innerTrack();   //--index();
    HecMuVar(0,iMu,iMuTrackRef);

    histo_muId->Fill( nMuVar[0][1] );

    for(std::vector<pat::Muon>::const_iterator jMu =iMu+1;
                                               jMu!=allmuons->end(); ++jMu){
       reco::TrackRef jMuTrackRef = jMu->innerTrack();
       HecMuVar(1,jMu,jMuTrackRef);

       bool goodDimuon = (nMuVar[0][1] == 111 && nMuVar[1][1] >= 11)
                       ||(nMuVar[1][1] == 111 && nMuVar[0][1] >= 11);     //--1 G and 1 (SG or T or ST or G)

       bool MySel = true;

       if( MyOneKaon_ ) MySel = nMuVar[0][0]*nMuVar[1][0] < 0;

          if( MySel ){     //--Skip same Sign muons for the kaon loop    [Mar 10, 12]

             if( goodDimuon ){         //--Take only good muons         
                                //--1  : S
                                //--10 : G
                                //--11 : G-S
                                //--100: T
                                //--101: T-S
                                //--110: T-G  doesn't exist
                                //--111: T-G-S      

              double Muu = ( iMu->p4() +  jMu->p4() ).M();
              histo_dimu->Fill( Muu );

              //--Do here the Tight Muon Variables  Feb 2013
              HecCutMu(0,iMu);
              HecCutMu(1,jMu);
	      
              fill_tightMuons();
              fill_softMuons();
            
              //--Constant for Vertexing
              float muon_sig   = muon_mass_c*1.e-6;
              float chi = 0., ndf = 0.;

              //--do Here dimuon Vertex
              reco::TransientTrack iMuTT( iMuTrackRef, &(*bFieldHandle) );
              reco::TransientTrack jMuTT( jMuTrackRef, &(*bFieldHandle) );

              KinematicParticleFactoryFromTransientTrack pFactory;
              /*            
              std::vector<RefCountedKinematicParticle> uuParticle;
              uuParticle.push_back(pFactory.particle( iMuTT, muon_mass_c, chi, ndf, muon_sig ));    //--Muon i
              uuParticle.push_back(pFactory.particle( jMuTT, muon_mass_c, chi, ndf, muon_sig ));    //--Muon j
              */
              //--Trajectory states to calculate DCA for the 2 tracks
              FreeTrajectoryState iState = pFactory.particle( iMuTT, muon_mass_c, chi, ndf, muon_sig )->currentState().freeTrajectoryState();
              FreeTrajectoryState jState = pFactory.particle( jMuTT, muon_mass_c, chi, ndf, muon_sig )->currentState().freeTrajectoryState();

              //--Measure distance between tracks at their closest approach for dimuon (iMu close to jMu)
              if (HecClosest(0, iState, jState)){
                const pat::Muon *patMuonP = &(*iMu);
                const pat::Muon *patMuonM = &(*jMu);

                if( MyOneKaon_ && nMuVar[0][0]<0 ){
                patMuonM = &(*iMu);
                patMuonP = &(*jMu);
                }

                HecHltMuTrig (0, patMuonP);
                HecHltMuTrig (1, patMuonM);

                iComb = 0;     //--to count only one the dimuon (i-j) combination
                //--dimuon mass constrain
                   if( MyOneKaon_ ){     //--Here dimuon with a Kaon (B- -> J/Psi(u+u-) K- )    
                   //bool uuMcut = std::abs( Muu - jpsi_mass_c ) < 0.2;     //--200 MeV
                   bool uuMcut = ( Muu > 2.5 && Muu < 4.0 );     //--200 MeV

                   double MyMass_0 = kaon_mass_c;
                   float  MySig_0  = MyMass_0*1.e-6;

                      if( uuMcut ){
                          //--Skip if no J/Psi Mass candidate
                          //--ITER-unsigned int k = -1;
                          //--for(reco::TrackCollection::const_iterator kTrack = tracks->begin(); kTrack != tracks->end(); ++kTrack) {  //--Kaons
                         for(unsigned int itrack = 0; itrack < theTrackRefs.size(); itrack++){  //--Kaons    
                            reco::TrackRef kTrack = theTrackRefs[itrack];

                            int kCharge = kTrack->charge();
                            double kEta = kTrack->eta();
                            double kPhi = kTrack->phi();
                            double kPt  = kTrack->pt();       //--sqrt( Px*Px + Py*Py )
                            double kPx  = kTrack->px();
                            double kPy  = kTrack->py();
                            double kPz  = kTrack->pz();
                            double kEn  = sqrt( MyMass_0*MyMass_0 + kPx*kPx + kPy*kPy + kPz*kPz);

                            double kCal = 0;
                            double kSeg = 0;
                            double kIso = 0;
                            double kKink= -1;//--Chi squared
                               if( MyPrint_)
                                  std::cout<<"Hec: kKink:"<<kKink<<std::endl;

                            int kGL       = 0;
                            int kmuId     = 0;

                            int kNtrkHits = kTrack->numberOfValidHits();
                            double kd0    = kTrack->d0();
                            double kdz    = kTrack->dz();

                            //--Fill the 4-Momentum Vector for Kaon  -->TLorentzVector a(0.0,0.0,0.0,0.0);
                            reco::Particle::LorentzVector kp4(0.0,0.0,0.0,0.0);
                            kp4.SetPxPyPzE( kTrack->px(), kTrack->py(), kTrack->pz(), kEn );

                            double Muuk = (iMu->p4() + jMu->p4() + kp4 ).M(); //--3 body invariant mass
                            double Muik = (iMu->p4() + kp4 ).M();           //--2 body invariant mass
                            double Mujk = (jMu->p4() + kp4 ).M();           //--2 body invariant mass

                            //Tihgt muons
                            double kNormChi2 = 0; 
                            double kMuMuHits = -1; 
                            double kMuStations = -1; 
                            double kMuPxHits = -1; 
                            double kMuTrkLayer = -1;

                            //Soft muons
                            double kMuInTrkLayer = -1;
                            double kMuInNormChi2 = 0;    

                            //fill_softKMuons(kMuInNormChi2,kMuInTrkLayer);
	      
                            //fill_tightKMuons(kNormChi2,kMuMuHits,kMuStations,kMuPxHits,kMuTrkLayer);

                            //fill_softKMuons(kMuInNormChi2,kMuInTrkLayer);     

                            bool noRepTrk = true;

                               if( kTrack == iMuTrackRef ){
                                   noRepTrk = false;
                                   histo_DeltaPku->Fill( MuVar[0][2] - kPt );
                                   histo_DeltaQ->Fill( nMuVar[0][0] - kCharge );
                                   histo_dimuk[0]->Fill( Muik );
                                   histo_dimuk[1]->Fill( Muuk );
                               }

                               if( kTrack == jMuTrackRef ){
                                   noRepTrk = false;
                                   histo_DeltaPku->Fill( MuVar[1][2] - kPt );
                                   histo_DeltaQ->Fill(  nMuVar[1][0] - kCharge );
                                   histo_dimuk[0]->Fill( Mujk );
                                   histo_dimuk[1]->Fill( Muuk );
                               }

                               if( MyPrint_ ){
                                if( !noRepTrk )std::cout<<"Hec: Kaon = Muon"<<std::endl;
                                double kEn1  = sqrt( MyMass_0*MyMass_0 + kTrack->px()*kTrack->px() + kTrack->py()*kTrack->py() + kTrack->pz()*kTrack->pz() );
                                double Ediff = fabs( kEn-kEn1 );
                                if( Ediff > 1./10. )std::cout<<Ediff<<" Hec:  Energy:"<<kEn<<" "<<kEn1<<std::endl;
                                 std::cout<<"Hec: "<<kp4.energy()<<" = kEn ="<<kEn<<std::endl;
                                 std::cout<<"Hec: Muuk "<<Muuk<<std::endl;
                               }

                               //--Here below combined 3 tracks and do Vertex and fill the ntuple
                               if( noRepTrk ){                             //--Skip repeated tracks 
                                  histo->Fill( kPt );
                                  histo_dimuk[2]->Fill( Muuk );

                                  if( kPt>1 && MuVar[0][2]>2 && MuVar[1][2]>2 ) histo_dimuk[3]->Fill( Muuk );

                                  //--dimuon + Kaon Vertex
                                  reco::TransientTrack kKaTT( *kTrack, &(*bFieldHandle) );

                                  KinematicParticleFactoryFromTransientTrack pFactory;

                                  std::vector<RefCountedKinematicParticle> uukParticle;
                                  uukParticle.push_back(pFactory.particle( iMuTT, muon_mass_c, chi, ndf, muon_sig ));    //--Muon i
                                  uukParticle.push_back(pFactory.particle( jMuTT, muon_mass_c, chi, ndf, muon_sig ));    //--Muon j
                                  uukParticle.push_back(pFactory.particle( kKaTT, MyMass_0   , chi, ndf, MySig_0  ));    //--Kaon k

                                  //--Make sure that all tracks are close together so we can try to find a common vertex
                                  FreeTrajectoryState kState = pFactory.particle( kKaTT, MyMass_0, chi, ndf, MySig_0 )->currentState().freeTrajectoryState();
                                  //--Measure distance between tracks at their closest approach  kmuon - imuon
                                  /*ClosestApproachInRPhi cApp; 
                                  GlobalPoint cxPt;*/
                                  if( HecClosest(1, kState, iState) ){   //--k close to iMu close
                                  if( HecClosest(2, kState, jState) ){   //--k close to jMu close

                                  KinematicFitDriver uukRec( uukParticle, "MuMuMk" );
                                  if( !uukRec.isValid() ){                                         //--Skip not valid mu-mu-ka vertex
                                     if( MyPrint_ )
                                        std::cout<<"Hec: Muu = "<<Muu<<" Muuk ="<<Muuk<<" BadFit uukVertex"<<std::endl;
                                     }
                                     else {
                                     //--Output Vertex fit variables
                                     double MuukFit                               = uukRec.mass();              //--Flag Dimuon kaon mass
                                     reco::Particle::LorentzVector  uukV_P4       = uukRec.P4();                //--DiMuKa 4-Momentum
                                     reco::Particle::LorentzVector  iMuV_P4       = uukRec.P4FirstChild();      //--Mui
                                     reco::Particle::LorentzVector  jMuV_P4       = uukRec.P4NextChild();       //--Muj
                                     reco::Particle::LorentzVector  kKaV_P4       = uukRec.P4NextChild();       //--Kaon
                                     reco::Particle::Point      uukV_vx           = uukRec.VertexDecay();       //--Vx,Vy,Vz 
                                     reco::Vertex::CovarianceMatrix uukV_Cov      = uukRec.VertexCov();
                                     double   uukV_Chi2                        = uukRec.chi2();
                                     double   uukV_Ndof                        = uukRec.ndof();
                                     double   uukV_Cl = ChiSquaredProbability( uukV_Chi2, uukV_Ndof );            //--Confidence Level uuk Vertex
                                        if( MyPrint_ ){
                                           std::cout<<"Hec: Collection uuk M= "<<Muu<<" Vtx ReFit: "<<MuukFit<<std::endl;
                                           std::cout<<"Hec: CL uuk "<<uukV_Cl<<" Vtx ReFit: "<<uukV_vx.x()<<" "<<uukV_vx.y()<<" "<<uukV_vx.z()<<std::endl;
                                        }
                                        //--do dimuon mass constrains

                                        //--do Primary here.  Primary Vertex Re-fitting to exclude the DiMuon and Kaon from the primary
                                        reco::Vertex refitPrimVertexHM = *primary.BestVertex();
                                        std::vector<reco::TrackRef> ExclusionList;
                                        ExclusionList.push_back( iMuTrackRef );
                                        ExclusionList.push_back( jMuTrackRef );
                                        ExclusionList.push_back( kTrack );

                                        VertexRefit RefitPrimaryHM( primary.BestVertex(), ExclusionList, bFieldHandle, beamSpot );   //--add highest multipl prim 11-march
                                        refitPrimVertexHM = RefitPrimaryHM.Refitted();
                                        
                                        //--Output Prim Vertex fit variables			      
                                        double pV_Chi2 = refitPrimVertexHM.chi2();
                                        double pV_Ndof = refitPrimVertexHM.ndof();
                                        double pV_ClHM = ChiSquaredProbability( pV_Chi2, pV_Ndof );//--Confidence Level Prim --mod pV_Cl to pV_ClHM

                                        //--here pointing march 11 2013
                                        GlobalPoint secver     = uukRec.RKVertex()->position(); 
                                        GlobalVector secmomTot = GlobalVector(uukRec.P4().Px(),
                                        				      uukRec.P4().Py(),
                                        				      uukRec.P4().Pz());
                                        
                                        if( MyPrint_ ){
                                        std::cout << "        : uukRec =" << uukRec.P4().Px() << std::endl;
                                        std::cout << "        : -----> =" << secmomTot      << std::endl;
                                        }

                                        VertexRefit RefitPrimary(primary.HigherCosAlpha(secver,secmomTot), ExclusionList, bFieldHandle, beamSpot);
                                        reco::Vertex refitPrimVertex = RefitPrimary.Refitted();
                                        
  
                                           if( !refitPrimVertex.isValid() ){                      //--Skip no valid Primary
                                              if( MyPrint_ )
                                                     std::cout<<"Hec: Mass: DiMu: "<<Muu<<" Prim BadFit:"<<std::endl;
                                                  //--continue; //--Goto the end of the if( uuC_.iGL!=0 && uuC_.jGL!=0 ) <--NO TRUE
                                           }
                                           else{
                                               /*  //--Output Prim Vertex fit variables                          
                                                 double pV_Chi2 = refitPrimVertexHM.chi2();
                                                 double pV_Ndof = refitPrimVertexHM.ndof();
                                                 double pV_ClHM = ChiSquaredProbability( pV_Chi2, pV_Ndof );//--Confidence Level Prim --mod pV_Cl to pV_ClHM */

                                                 //histo_pV_ClHM->Fill(pV_ClHM);

                                                 //--here pointing march 11 2013
                                                 /*GlobalPoint secver     = uukRec.RKVertex()->position(); 
                                                 GlobalVector secmomTot = GlobalVector(uukRec.P4().Px(),
                                                                                       uukRec.P4().Py(),
                                                                                       uukRec.P4().Pz());
                                                 
                                                 if( MyPrint_ ){
                                                 std::cout << "        : uukRec =" << uukRec.P4().Px() << std::endl;
                                                 std::cout << "        : -----> =" << secmomTot      << std::endl;
                                                 }

                                                 VertexRefit RefitPrimary(primary.HigherCosAlpha(secver,secmomTot), ExclusionList, bFieldHandle, beamSpot);
                                                 reco::Vertex refitPrimVertex = RefitPrimary.Refitted();*/
 
                                                 double pV_Cl = ChiSquaredProbability(refitPrimVertex.chi2(), refitPrimVertex.ndof());
                                                
                                                 double primDist = sqrt( refitPrimVertex.x()*refitPrimVertex.x()
                                                                       + refitPrimVertex.y()*refitPrimVertex.y()
                                                                       + refitPrimVertex.z()*refitPrimVertex.z() );
                                                 double primR    = sqrt( refitPrimVertex.x()*refitPrimVertex.x()
                                                                       + refitPrimVertex.y()*refitPrimVertex.y() );
 
                                                 double primDz  = refitPrimVertex.z() - refitPrimVertexHM.z();
                                                 double primDxy = sqrt( (refitPrimVertex.x() - refitPrimVertexHM.x())*(refitPrimVertex.x() - refitPrimVertexHM.x())
                                                                      + (refitPrimVertex.y() - refitPrimVertexHM.y())*(refitPrimVertex.y() - refitPrimVertexHM.y()) );
                                                 double primDxyz= sqrt( (refitPrimVertex.x() - refitPrimVertexHM.x())*(refitPrimVertex.x() - refitPrimVertexHM.x())
                                                                      + (refitPrimVertex.y() - refitPrimVertexHM.y())*(refitPrimVertex.y() - refitPrimVertexHM.y())
                                                                      + (refitPrimVertex.z() - refitPrimVertexHM.z())*(refitPrimVertex.z() - refitPrimVertexHM.z()) );
              
                                                 if( MyPrint_ ){
                                                 std::cout<<"meson_b: primCL ="<<pV_Cl<<" High Mult: "<<pV_ClHM<<std::endl;
                                                 std::cout<<"meson_b: primDz ="<<primDz<<" primDxy ="<<primDxy<<" primDxyz ="<<primDxyz<<std::endl;
                                                 }
                                                 //--here end pointing

                                                 //--Check separation from primary (L over Sigma)
                                                 typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
                                                 typedef ROOT::Math::SVector<double, 3> SVector3;

                                                 SMatrixSym3D totalCov = refitPrimVertex.covariance() + uukV_Cov;
                                                 SVector3 uuk_PrimarySep3D( uukV_vx.x() - refitPrimVertex.x(),
                                                                            uukV_vx.y() - refitPrimVertex.y(),
                                                                            uukV_vx.z() - refitPrimVertex.z() );

                                                 double uukL     = ROOT::Math::Mag( uuk_PrimarySep3D );  //--Flag 3dim distance Dimuon-Prim
                                                 double uukSigma = sqrt(ROOT::Math::Similarity( totalCov, uuk_PrimarySep3D )) / uukL;

                                                 double uukLoS = uukL / uukSigma;  //--Flag 3dim distance/sigma  Dimuon-Prim
                                                    if( MyPrint_ )
                                                       std::cout<<"Hec: L= "<<uukL<<" LoS= " << uukLoS<<" MassFitted"<<MuukFit<<std::endl;

                                                 //--Here Pointing Secondary Vertex  [Sep 9, 2011]
                                                 GlobalVector uuk_PrimSec( uukV_vx.x() - refitPrimVertex.x(),
                                                                           uukV_vx.y() - refitPrimVertex.y(),
                                                                           uukV_vx.z() - refitPrimVertex.z() );
                                                 GlobalVector uuk_P3( uukRec.P4().x(),
                                                                      uukRec.P4().y(),
                                                                      uukRec.P4().z() );

                                                 double V_dot_Puuk = uuk_P3.dot(uuk_PrimSec);
                                                 double cos_alpha = V_dot_Puuk /( uuk_PrimSec.mag() * uuk_P3.mag() );

                                                    if( MyPrint_ )
                                                       std::cout<<"Hec: Pointing cos alpha= "<<cos_alpha<<std::endl;

                                                    if( kPt>1 && MuVar[0][2]>2 && MuVar[1][2]>2 && uukLoS>5 && cos_alpha > 0.9 )
                                                       histo_dimuk[4]->Fill( Muuk );

                                                 //int ukTrig = 0;
                                                 //int mukLMtrig = 0;
   
                                                 histo_pV_ClHM->Fill( pV_ClHM );
						 
                                                 //--fill n-tuple
                                                 iComb++;
                                                 fill_evt( iEvent, Ntk, allMu, allPrim, iComb);          //--Fill Ntuple (Branch: evt)

                                                 fill_IuuC();

                                                 fill_uuC(Muu );

                                                 fill_IuukC(kGL,kCharge,kmuId,kNtrkHits );

                                                 fill_tightKMuons(kNormChi2,kMuMuHits,kMuStations,kMuPxHits,kMuTrkLayer);

                                                 fill_softKMuons(kMuInNormChi2,kMuInTrkLayer);     

                                                 fill_uukC(kEta,kPhi,kPt,kPx,kPy,kPz,
                                                           kd0,kdz,kCal,kSeg,kIso,kKink,
                                                           Muuk,Muik,Mujk,MuukFit,uukV_Cl,pV_Cl,primDist,primR,primDz,primDxy,primDxyz,
                                                           uukL,uukLoS,cos_alpha,iMuV_P4.pt(),jMuV_P4.pt(),kKaV_P4.pt(),
                                                           uukV_vx.x(),uukV_vx.y(),uukV_vx.z(),refitPrimVertex.x(),refitPrimVertex.y(),refitPrimVertex.z() );

                                                 hecmu_tree_->Fill();
                                                 //init();
                                              }  //--End Valid Primary Vertex

                                     }   //--End Valid Dimuon + Kaon Vertex

                                  }     //--End for muon j - kaon Closest Approach

                                  }     //--End for muon i - kaon Closest Approach

                               }        //--End if 3 different tracks if( k != iMuTrackRef.index() && k != jMuTrackRef.index()

                            }           //--End for for(TrackCollection Kaon

                         }              //--End if( uuMcut )         no J/Psi Mass candidate

                      }                 //--End for if( MyOneKaon_ )  B --> J/Psi + K 

                  //===> Here I add the third muon <===
                      if( MyThreeMu_ ){
                         for(std::vector<pat::Muon>::const_iterator kMu=jMu+1;kMu != allmuons->end();++kMu){

                            int kCharge = kMu->charge();
                            double kEta = kMu->eta();
                            double kPhi = kMu->phi();
                            double kPt  = kMu->p4().Pt();
                            double kPx  = kMu->p4().Px();
                            double kPy  = kMu->p4().Py();
                            double kPz  = kMu->p4().Pz();
                            double kCal = kMu->caloCompatibility();
                            double kSeg = muon::segmentCompatibility(*kMu);
                            double kIso = kMu->isolationR03().sumPt;
                            double kKink= kMu->combinedQuality().trkKink;//--Chi squared
                            if( MyPrint_ )
                                std::cout<<"Hec: kKink:"<<kKink<<std::endl;

                            reco::TrackRef kMuTrackRef = kMu->innerTrack();   //--.index();

                            int muId=0, muSA=0, muGL=0, muTK=0;
                            if( kMu->isStandAloneMuon() )muSA = 1;        //--2^0
                            if( kMu->isGlobalMuon()     )muGL = 10;    //--2^1
                            if( kMu->isTrackerMuon()    )muTK = 100;   //--2^2
                            int kmuId = muTK + muGL + muSA;
                            if( MyPrint_ )
                                std::cout<<"Hec: MuId: kTK:"<<muTK<<" kGL:"<<muGL<<" SA:"<<muSA<<" Dec:"<<muId<<" Bin:"<<kmuId<<std::endl;

                            int kGL = 0;
                            if(  kMu->isGlobalMuon()  ) kGL = 1;
                            if(!(kMu->isGlobalMuon()) && kMu->isTrackerMuon() ) kGL = 2;

                            if( kmuId>=11 ){     //--Skip No Global and No Tracker Muons.  Take only good muons 
                                                 //--1  : S
                                                 //--10 : G
                                                 //--11 : G-S
                                                 //--100: T
                                                 //--101: T-S
                                                 //--110: T-G  doesn't exist
                                                 //--111: T-G-S  
                               int kNtrkHits = kMuTrackRef->numberOfValidHits();
                               double kd0    = kMuTrackRef->d0();
                               double kdz    = kMuTrackRef->dz();

                               double Muuk = (iMu->p4() + jMu->p4() + kMu->p4() ).M();  //--3 body invariant mass
                               double Muik = (iMu->p4() + kMu->p4() ).M();                      //--2 body invariant mass
                               double Mujk = (jMu->p4() + kMu->p4() ).M();                      //--2 body invariant mass

                                  if( MyPrint_ ){
                                     std::cout<<"Hec: Third Muon"<<std::endl;
                                     std::cout<<"Hec: Muuk "<<Muuk<<std::endl;
                                  }

                               //--Here below combined 3 tracks and do Vertex and fill the ntuple
                               histo->Fill( kPt );

                               histo_uuu[0]->Fill( Muuk );

                                  if( kPt>1 && MuVar[0][2]>2 && MuVar[1][2]>2 ) histo_uuu[1]->Fill( Muuk );

                               //--Do here the Tight KMuon Variables  Feb 2013
                               //HecCutMu(2,kMu);
                               //Tihgt muons
                               double kNormChi2 = 99; 
                               double kMuMuHits = -1; 
                               
                               if ( kMu->isGlobalMuon() ){
                                  kNormChi2 = kMu->globalTrack()->normalizedChi2();
                                  kMuMuHits = kMu->globalTrack()->hitPattern().numberOfValidMuonHits();
                               }

			       double kMuStations = kMu->numberOfMatchedStations(); 
                               double kMuPxHits   = kMu->innerTrack()->hitPattern().numberOfValidPixelHits(); 
                               double kMuTrkLayer = kMu->track()->hitPattern().trackerLayersWithMeasurement();

                               //Soft muons
                               double kMuInNormChi2 = 99;
                               double kMuInTrkLayer = -1;    

                               if ( kMu->isTrackerMuon() ){
                                  kMuInNormChi2 = kMu->innerTrack()->normalizedChi2();
                                  kMuInTrkLayer = kMu->innerTrack()->hitPattern().pixelLayersWithMeasurement();
                               }

                               fill_tightKMuons(kNormChi2,kMuMuHits,kMuStations,kMuPxHits,kMuTrkLayer);

                               fill_softKMuons(kMuInNormChi2,kMuInTrkLayer);     

                               //--dimuon + Kaon Vertex
                               reco::TransientTrack kKaTT( *kMuTrackRef  , &(*bFieldHandle) );

                               KinematicParticleFactoryFromTransientTrack pFactory;

                               std::vector<RefCountedKinematicParticle> uukParticle;
                               uukParticle.push_back(pFactory.particle( iMuTT, muon_mass_c, chi, ndf, muon_sig ));    //--Muon i
                               uukParticle.push_back(pFactory.particle( jMuTT, muon_mass_c, chi, ndf, muon_sig ));    //--Muon j
                               uukParticle.push_back(pFactory.particle( kKaTT, muon_mass_c, chi, ndf, muon_sig ));    //--Muon k

                               //--Make sure that all tracks are close together so we can try to find a common vertex
                               FreeTrajectoryState kState = pFactory.particle( kKaTT, muon_mass_c, chi, ndf, muon_sig )->currentState().freeTrajectoryState();

                               //--Measure distance between tracks at their closest approach
                               /*ClosestApproachInRPhi cApp; 
                               GlobalPoint cxPt;*/
                                  if( HecClosest(4, kState, iState) ){  //--kMu close to iMu close
                                  if( HecClosest(5, kState, jState) ){  //--kMu close to jMu close

                                     KinematicFitDriver uukRec( uukParticle, "MuMuMk" );
                                     if( !uukRec.isValid() ){                                          //--Skip not valid mu-mu-ka vertex

                                        if( MyPrint_ )
                                           std::cout<<"Hec: Muu = "<<Muu<<" Muuk ="<<Muuk<<" BadFit uukVertex"<<std::endl;
                                     }
                                     else{

                                     //--Output Vertex fit variables
                                     double MuukFit                               = uukRec.mass();                 //--Flag Dimuon kaon mass
                                     reco::Particle::LorentzVector  uukV_P4       = uukRec.P4();                   //--DiMuKa 4-Momentum
                                     reco::Particle::LorentzVector  iMuV_P4       = uukRec.P4FirstChild();         //--Mui
                                     reco::Particle::LorentzVector  jMuV_P4       = uukRec.P4NextChild();          //--Muj
                                     reco::Particle::LorentzVector  kKaV_P4       = uukRec.P4NextChild();          //--Kaon
                                     reco::Particle::Point         uukV_vx        = uukRec.VertexDecay();          //--Vx,Vy,Vz 
                                     reco::Vertex::CovarianceMatrix uukV_Cov      = uukRec.VertexCov();
                                     double   uukV_Chi2                           = uukRec.chi2();
                                     double   uukV_Ndof                           = uukRec.ndof();
                                     double   uukV_Cl = ChiSquaredProbability( uukV_Chi2, uukV_Ndof );               //--Confidence Level uuk Vertex

                                        if( MyPrint_ ){
                                           std::cout<<"Hec: Collection uuk M= "<<Muu<<" Vtx ReFit: "<<MuukFit<<std::endl;
                                           std::cout<<"Hec: CL uuk "<<uukV_Cl<<" Vtx ReFit: "<<uukV_vx.x()<<" "<<uukV_vx.y()<<" "<<uukV_vx.z()<<std::endl;
                                        }
                                     //--do dimuon mass constrains

                                     //--do Primary here.  Primary Vertex Re-fitting to exclude the DiMuon and Kaon from the primary
                                     reco::Vertex refitPrimVertexHM = *primary.BestVertex();
                                     std::vector<reco::TrackRef> ExclusionList;
                                     ExclusionList.push_back( iMuTrackRef );
                                     ExclusionList.push_back( jMuTrackRef );
                                     ExclusionList.push_back( kMuTrackRef );

                                     VertexRefit RefitPrimaryHM( primary.BestVertex(), ExclusionList, bFieldHandle, beamSpot );
                                     refitPrimVertexHM = RefitPrimaryHM.Refitted();
                                     
                                     //--Output Prim Vertex fit variables                                
                                     double pV_Chi2 = refitPrimVertexHM.chi2();
                                     double pV_Ndof = refitPrimVertexHM.ndof();
                                     double pV_ClHM = ChiSquaredProbability( pV_Chi2, pV_Ndof );//--Confidence Level Prim --mod pV_Cl to pV_ClHM

                                     //histo_pV_ClHM->Fill(pV_ClHM);

                                     //--here pointing march 11 2013
                                     GlobalPoint secver     = uukRec.RKVertex()->position(); 
                                     GlobalVector secmomTot = GlobalVector(uukRec.P4().Px(),
                                     					   uukRec.P4().Py(),
                                     					   uukRec.P4().Pz());
                                     
                                     if( MyPrint_ ){
                                     std::cout << "	   : uukRec =" << uukRec.P4().Px() << std::endl;
                                     std::cout << "	   : -----> =" << secmomTot	 << std::endl;
                                     }

                                     VertexRefit RefitPrimary(primary.HigherCosAlpha(secver,secmomTot), ExclusionList, bFieldHandle, beamSpot);
                                     reco::Vertex refitPrimVertex = RefitPrimary.Refitted();
                                        
                                        if( !refitPrimVertex.isValid() ){                        //--Skip no valid Primary
                                           if( MyPrint_ )
                                              std::cout<<"Hec: Mass: DiMu: "<<Muu<<" Prim BadFit:"<<std::endl;
                                              //--continue;     //--Goto the end of the if( uuC_.iGL!=0 && uuC_.jGL!=0 ) <--NO TRUE
                                        }
                                        else {
                                        /*//--Output Prim Vertex fit variables                                
                                        double pV_Chi2 = refitPrimVertexHM.chi2();
                                        double pV_Ndof = refitPrimVertexHM.ndof();
                                        double pV_ClHM = ChiSquaredProbability( pV_Chi2, pV_Ndof );//--Confidence Level Prim --mod pV_Cl to pV_ClHM

                                        //histo_pV_ClHM->Fill(pV_ClHM);

                                        //--here pointing march 11 2013
                                        GlobalPoint secver     = uukRec.RKVertex()->position(); 
                                        GlobalVector secmomTot = GlobalVector(uukRec.P4().Px(),
                                        				      uukRec.P4().Py(),
                                        				      uukRec.P4().Pz());
                                        
                                        if( MyPrint_ ){
                                        std::cout << "        : uukRec =" << uukRec.P4().Px() << std::endl;
                                        std::cout << "        : -----> =" << secmomTot      << std::endl;
                                        }

                                        VertexRefit RefitPrimary(primary.HigherCosAlpha(secver,secmomTot), ExclusionList, bFieldHandle, beamSpot);
                                        reco::Vertex refitPrimVertex = RefitPrimary.Refitted();*/
					
                                        double pV_Cl = ChiSquaredProbability(refitPrimVertex.chi2(), refitPrimVertex.ndof());
                                        
                                        double primDist = sqrt( refitPrimVertex.x()*refitPrimVertex.x()
                                        		      + refitPrimVertex.y()*refitPrimVertex.y()
                                        		      + refitPrimVertex.z()*refitPrimVertex.z() );
                                        double primR	= sqrt( refitPrimVertex.x()*refitPrimVertex.x()
                                        		      + refitPrimVertex.y()*refitPrimVertex.y() );
 
                                        double primDz  = refitPrimVertex.z() - refitPrimVertexHM.z();
                                        double primDxy = sqrt( (refitPrimVertex.x() - refitPrimVertexHM.x())*(refitPrimVertex.x() - refitPrimVertexHM.x())
                                        		     + (refitPrimVertex.y() - refitPrimVertexHM.y())*(refitPrimVertex.y() - refitPrimVertexHM.y()) );
                                        double primDxyz= sqrt( (refitPrimVertex.x() - refitPrimVertexHM.x())*(refitPrimVertex.x() - refitPrimVertexHM.x())
                                        		     + (refitPrimVertex.y() - refitPrimVertexHM.y())*(refitPrimVertex.y() - refitPrimVertexHM.y())
                                        		     + (refitPrimVertex.z() - refitPrimVertexHM.z())*(refitPrimVertex.z() - refitPrimVertexHM.z()) );
              
                                        if( MyPrint_ ){
                                        std::cout<<"meson_b: primCL ="<<pV_Cl<<" High Mult: "<<pV_ClHM<<std::endl;
                                        std::cout<<"meson_b: primDz ="<<primDz<<" primDxy ="<<primDxy<<" primDxyz ="<<primDxyz<<std::endl;
                                        }
                                        //--here end pointing

                                        //--Check separation from primary (L over Sigma)
                                        typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
                                        typedef ROOT::Math::SVector<double, 3> SVector3;

                                        SMatrixSym3D totalCov = refitPrimVertex.covariance() + uukV_Cov;
                                        SVector3 uuk_PrimarySep3D( uukV_vx.x() - refitPrimVertex.x(),
                                                                 uukV_vx.y() - refitPrimVertex.y(),
                                                                 uukV_vx.z() - refitPrimVertex.z() );

                                        double uukL     = ROOT::Math::Mag( uuk_PrimarySep3D );                      //--Flag 3dim distance Dimuon-Prim
                                        double uukSigma = sqrt(ROOT::Math::Similarity( totalCov, uuk_PrimarySep3D )) / uukL;

                                        double uukLoS = uukL / uukSigma;                                              //--Flag 3dim distance/sigma  Dimuon-Prim
                                           if( MyPrint_ )
                                              std::cout<<"Hec: L= "<<uukL<<" LoS= " << uukLoS<<" MassFitted"<<MuukFit<<std::endl;
                                       //--Here Pointing Secondary Vertex  [Sep 9, 2011]
                                        GlobalVector uuk_PrimSec( uukV_vx.x() - refitPrimVertex.x(),
                                                                  uukV_vx.y() - refitPrimVertex.y(),
                                                                  uukV_vx.z() - refitPrimVertex.z() );
                                        GlobalVector uuk_P3( uukRec.P4().x(),
                                                             uukRec.P4().y(),
                                                             uukRec.P4().z() );

                                        double V_dot_Puuk =uuk_P3.dot(uuk_PrimSec);
                                        double cos_alpha = V_dot_Puuk /( uuk_PrimSec.mag() * uuk_P3.mag() );
                                           if( MyPrint_ )
                                              std::cout<<"Hec: Pointing cos alpha= "<<cos_alpha<<std::endl;

                                           if( kPt>1 && MuVar[0][2]>2 && MuVar[1][2]>2 && uukLoS>5 && cos_alpha > 0.9 )
                                              histo_uuu[2]->Fill( Muuk );


                                        //--Trigger from Keith [Apr 2, 12]
                                        const pat::Muon *patMuonK = &(*kMu);

                                        //--unpacks J/Psi trigger bit
                                        HecHltMuTrig(2, patMuonK);

                                        histo_pV_ClHM->Fill( pV_ClHM );

                                        //--fill n-tuple
                                        iComb++;
                                        fill_evt( iEvent, Ntk, allMu, allPrim, iComb);//--Fill Ntuple (Branch: evt)

                                        fill_IuuC();

                                        fill_uuC(Muu);

                                        fill_IuukC(kGL,kCharge,kmuId,kNtrkHits );

                                        fill_uukC(kEta,kPhi,kPt,kPx,kPy,kPz,
                                                  kd0,kdz,kCal,kSeg,kIso,kKink,
                                                  Muuk,Muik,Mujk,MuukFit,uukV_Cl,pV_Cl,primDist,primR,primDz,primDxy,primDxyz,
                                                  uukL,uukLoS,cos_alpha,iMuV_P4.pt(),jMuV_P4.pt(),kKaV_P4.pt(),
                                                  uukV_vx.x(),uukV_vx.y(),uukV_vx.z(),refitPrimVertex.x(),refitPrimVertex.y(),refitPrimVertex.z() );

                                       // fill_softKMuons(kMuInNormChi2,kMuInTrkLayer);
	      
	                               // fill_tightKMuons(kNormChi2,kMuMuHits,kMuStations,kMuPxHits,kMuTrkLayer);
      
					hecmu_tree_->Fill();
                                        //init();

                                        }       //--End Primary Vertex Valid  

                                     }        //--End uuu Vertex Valid   

                                  }         //--End if closest approach muon K - muon j 

                                  }          //--End if closest approach muon K - muon i

                            }             //--End if good Muons k

                         }        //--End for(reco::MuonCollection::const_iterator kMuon

                      }           //--End if( MyThreeMu_ ){  3 muons for tau

        //}       //--End Valid dimuon Vertex
              }              //--End dimuon Closest Approach        
             }               //--End if good Muons
          }               //--End dimuon Opposite charge
    }               //--End for(reco::MuonCollection::const_iterator jMuon
 }                //--End for(reco::MuonCollection::const_iterator iMuon

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
HecTau::beginJob()
{
std::cout << " HecTau::beginJob" << std::endl;
 hecmu_tree_ = new TTree("ThreeBody","Dettach Ntuple");
 int bufsize = 64000;

 hecmu_tree_->Branch("Dimuon"  , &uuC_,      "iEta/D:iPhi:iPt:iPx:iPy:iPz:id0:idz:iCal:iSeg:iIso:iKink:jEta:jPhi:jPt:jPx:jPy:jPz:jd0:jdz:jCal:jSeg:jIso:jKink:Muu:dca:xcpt:ycpt:zcpt",bufsize);
 hecmu_tree_->Branch("TreePart", &uukC_,     "kEta/D:kPhi:kPt:kPx:kPy:kPz:kd0:kdz:kCal:kSeg:kIso:kKink:Muuk:Muik:Mujk:MuukFit:uukV_Cl:pV_Cl:primDist:primR:primDz:primDxy:primDxyz:uukL:uukLoS:cos_alpha:iVpt:jVpt:kVpt:uukVx:uukVy:uukVz:PrimVx:PrimVy:PrimVz", bufsize);
 hecmu_tree_->Branch("TightMu" , &tightMuC_, "iNormChi2/D:iMuMuHits:iMuStations:iMuPxHits:iMuTrkLayer:jNormChi2:jMuMuHits:jMuStations:jMuPxHits:jMuTrkLayer", bufsize);
 hecmu_tree_->Branch("SoftMu"  , &softMuC_,  "iMuInNormChi2/D:iMuInTrkLayer:jMuInNormChi2:jMuInTrkLayer", bufsize);
 hecmu_tree_->Branch("TightKMu", &tightKMuC_,"kNormChi2/D:kMuMuHits:kMuStations:kMuPxHits:kMuTrkLayer", bufsize);
 hecmu_tree_->Branch("SoftKMu" , &softKMuC_, "kMuInNormChi2/D:kMuInTrkLayer", bufsize);
 
 hecmu_tree_->Branch("iHeader", &Ievt_, "run/I:evtnum:lumib:Ntk:allMu:allPrim:iComb:upTrig:umTrig:umLMtrig:upLMtrig:ukTrig:mukLMtrig:tevtTrig1:evtTrig2:evtTrig3:evtTrig4", bufsize);
 hecmu_tree_->Branch("iDimuon", &IuuC_, "iQ/I:imuId:iNtrkHits:jQ:jmuId:jNtrkHits", bufsize);
 hecmu_tree_->Branch("kDimuon", &IuukC_,"kGL/I:kCharge:kmuId:kNtrkHits", bufsize);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HecTau::endJob() 
{
  std::cout << " HecTau::endJob" << std::endl;
}

// ------------ method called when starting to processes a run  ------------
void 
HecTau::beginRun(edm::Run const&, edm::EventSetup const&)
{
  std::cout << " HecTau::beginRun" << std::endl;
}

// ------------ method called when ending the processing of a run  ------------
void 
HecTau::endRun(edm::Run const&, edm::EventSetup const&)
{
  std::cout << " HecTau::endRun" << std::endl;
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
HecTau::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
  std::cout << " HecTau::beginLuminosityBlock" << std::endl;
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
HecTau::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
  std::cout << " HecTau::endLuminosityBlock" << std::endl;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HecTau::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
//--from here below are my functions
void HecTau::init()
{
  Ievt_.init();
  IuuC_.init();
  uuC_.init();
  IuukC_.init();
  uukC_.init();
  tightMuC_.init();
  softMuC_.init();
//  tightKMuC_.init();
//  softKMuC_.init();  
}
void HecTau::Ievt::init()
{
  runNb      = -1;
  eventNb    = -1;
  lumiBlock  = -1;
  Ntk        = -1;
  allMu      = -1;
  allPrim    = -1;
  iComb      =  0;
  upTrig     = -1;
  umTrig     = -1;
  umLMtrig   = -1;
  upLMtrig   = -1;
  ukTrig     = -1;
  mukLMtrig  = -1;
  evtTrig1   =  0;
  evtTrig2   =  0;
  evtTrig3   =  0;
  evtTrig4   =  0;
}
void HecTau::fill_evt(const edm::Event& iEvent, int Ntk, int allMu, int allPrim,int iComb)
{
     Ievt_.runNb      = iEvent.id().run();
     Ievt_.eventNb    = iEvent.id().event();
     Ievt_.lumiBlock  = iEvent.luminosityBlock();
     Ievt_.Ntk        = Ntk;            //--tracks->size();
     Ievt_.allMu      = allMu;          //--->allmuons->size();
     Ievt_.allPrim    = allPrim;
     Ievt_.iComb      = iComb;    //--HMoreno----Ntuple Trig Var
     Ievt_.upTrig     = nTrig[0][0];   //--HMoreno----Ntuple Trig Var
     Ievt_.umTrig     = nTrig[1][0];   //--HMoreno----Ntuple Trig Var
     Ievt_.ukTrig     = nTrig[2][0];   //--HMoreno----Ntuple Trig Var
     Ievt_.upLMtrig   = nTrig[0][1];   //--HMoreno----Ntuple Trig Var
     Ievt_.umLMtrig   = nTrig[1][1];   //--HMoreno----Ntuple Trig Var
     Ievt_.mukLMtrig  = nTrig[2][1];   //--HMoreno----Ntuple Trig Var
     Ievt_.evtTrig1   = evtTrig1; //--HMoreno----Ntuple Trig Var
     Ievt_.evtTrig2   = evtTrig2; //--HMoreno----Ntuple Trig Var
     Ievt_.evtTrig3   = evtTrig3; //--HMoreno----Ntuple Trig Var
     Ievt_.evtTrig4   = evtTrig4; //--HMoreno----Ntuple Trig Var
}
void HecTau::IuuC::init()
{
/*  iGL       = -1;
  iCharge   =  0;
  imuId     = -1;
  iNtrkHits = -1;
  jGL       = -1;
  jCharge   =  0;
  jmuId     = -1;
  jNtrkHits = -1;*/
  iQ = 0; imuId = -1; iNtrkHits = -1;
  jQ = 0; jmuId = -1; jNtrkHits = -1;

}
void HecTau::fill_IuuC()
{
  IuuC_.iQ        = nMuVar[0][0];
  IuuC_.imuId     = nMuVar[0][1];
  IuuC_.iNtrkHits = nMuVar[0][4];
  IuuC_.jQ        = nMuVar[1][0];
  IuuC_.jmuId     = nMuVar[1][1];
  IuuC_.jNtrkHits = nMuVar[1][4];
}
void HecTau::uuC::init()
{
  iEta    = -1.0;
  iPhi    = -1.0;
  iPt     = -1.0;
  iPx     = -1.0;
  iPy     = -1.0;
  iPz     = -1.0;
  id0     = -1.0;
  idz     = -1.0;
  iCal    = -1.0;
  iSeg    = -1.0;
  iIso    = -1.0;
  iKink   = -1.0;
  jEta    = -1.0;
  jPhi    = -1.0;
  jPt     = -1.0;
  jPx     = -1.0;
  jPy     = -1.0;
  jPz     = -1.0;
  jd0     = -1.0;
  jdz     = -1.0;
  jCal    = -1.0;
  jSeg    = -1.0;
  jIso    = -1.0;
  jKink   = -1.0;
  Muu     = -1.0;
  dca     = -1.0;
  xcpt    = -1.0;
  ycpt    = -1.0;
  zcpt    = -1.0;
}
void HecTau::fill_uuC( double Muu)
{
  uuC_.iEta    = MuVar[0][0];
  uuC_.iPhi    = MuVar[0][1];
  uuC_.iPt     = MuVar[0][2];
  uuC_.iPx     = MuVar[0][4];
  uuC_.iPy     = MuVar[0][5];
  uuC_.iPz     = MuVar[0][6];
  uuC_.id0     = MuVar[0][11];
  uuC_.idz     = MuVar[0][12];
  uuC_.iCal    = MuVar[0][7];
  uuC_.iSeg    = MuVar[0][8];
  uuC_.iIso    = MuVar[0][9];
  uuC_.iKink   = MuVar[0][10];//4th part in including variable
  uuC_.jEta    = MuVar[1][0];
  uuC_.jPhi    = MuVar[1][1];
  uuC_.jPt     = MuVar[1][2];
  uuC_.jPx     = MuVar[1][4];
  uuC_.jPy     = MuVar[1][5];
  uuC_.jPz     = MuVar[1][6];
  uuC_.jd0     = MuVar[1][11];
  uuC_.jdz     = MuVar[1][12];
  uuC_.jCal    = MuVar[1][7];
  uuC_.jSeg    = MuVar[1][8];
  uuC_.jIso    = MuVar[1][9];
  uuC_.jKink   = MuVar[1][10];//
  uuC_.Muu     = Muu ;
  uuC_.dca     = dcaVar[0][0];
  uuC_.xcpt    = dcaVar[0][1];
  uuC_.ycpt    = dcaVar[0][2];
  uuC_.zcpt    = dcaVar[0][3];
}
void HecTau::IuukC::init()
{
  kGL       = -1;
  kCharge   =  0;
  kmuId     = -1;
  kNtrkHits = -1;
}
void HecTau::fill_IuukC(int kGL, int kCharge, int kmuId, int kNtrkHits)
{
  IuukC_.kGL       = kGL;
  IuukC_.kCharge   = kCharge;
  IuukC_.kmuId     = kmuId;
  IuukC_.kNtrkHits = kNtrkHits;
}
void HecTau::uukC::init()
{
  kEta      = -1.0;
  kPhi      = -1.0;
  kPt       = -1.0;
  kPx       = -1.0;
  kPy       = -1.0;
  kPz       = -1.0;
  kd0       = -1.0;
  kdz       = -1.0;
  kCal      = -1.0;
  kSeg      = -1.0;
  kIso      = -1.0;
  kKink     = -1.0;
  Muuk      = -1.0;
  Muik      = -1.0;
  Mujk      = -1.0;
  MuukFit   = -1.0;
  uukV_Cl   = -1.0;
  pV_Cl     = -1.0;
  primDist  = -1.0;
  primR     = -1.0;
  primDz    = -1.0;
  primDxy   = -1.0;
  primDxyz  = -1.0;
  uukL      = -1.0;
  uukLoS    = -1.0;
  cos_alpha = -1.5;
  iVpt      = -1.0;
  jVpt      = -1.0;
  kVpt      = -1.0;
  uukVx     = -999.0;
  uukVy     = -999.0;
  uukVz     = -999.0;
  PrimVx    = -999.0;
  PrimVy    = -999.0;
  PrimVz    = -999.0;
}
void HecTau::fill_uukC(double kEta,double kPhi,double kPt,double kPx,double kPy,double kPz,double kd0,double kdz,double kCal,double kSeg,double kIso,
                       double kKink,
                       double Muuk,double Muik,double Mujk,double MuukFit,double uukV_Cl,double pV_Cl,double primDist,double primR,double primDz,double primDxy,double primDxyz,
                       double uukL,double uukLoS,double cos_alpha,double iVpt,double jVpt,double kVpt,
                       double uukVx,double uukVy,double uukVz,double PrimVx,double PrimVy,double PrimVz )
{
  uukC_.kEta      = kEta;
  uukC_.kPhi      = kPhi;
  uukC_.kPx       = kPt;
  uukC_.kPx       = kPx;
  uukC_.kPy       = kPy;
  uukC_.kPz       = kPz;
  uukC_.kd0       = kd0;
  uukC_.kdz       = kdz;
  uukC_.kCal      = kCal;
  uukC_.kSeg      = kSeg;
  uukC_.kIso      = kIso;
  uukC_.kKink     = kKink;
  uukC_.Muuk      = Muuk;
  uukC_.Muik      = Muik;
  uukC_.Mujk      = Mujk;
  uukC_.MuukFit   = MuukFit;
  uukC_.uukV_Cl   = uukV_Cl;
  uukC_.pV_Cl     = pV_Cl;
  uukC_.primDist  = primDist;
  uukC_.primR     = primR;
  uukC_.primDz    = primDz;
  uukC_.primDxy   = primDxy;
  uukC_.primDxyz  = primDxyz;
  uukC_.uukL      = uukL;
  uukC_.uukLoS    = uukLoS;
  uukC_.cos_alpha = cos_alpha;
  uukC_.iVpt      = iVpt;
  uukC_.jVpt      = jVpt;
  uukC_.kVpt      = kVpt;
  uukC_.uukVx     = uukVx;
  uukC_.uukVy     = uukVy;
  uukC_.uukVz     = uukVz;
  uukC_.PrimVx    = PrimVx;
  uukC_.PrimVy    = PrimVy;
  uukC_.PrimVz    = PrimVz;
}
void HecTau::tightMuC::init()
{
  iNormChi2   = 99;
  iMuMuHits   = -1;
  iMuStations = 9696; 
  iMuPxHits   = 9696; 
  iMuTrkLayer = 9696;
  jNormChi2   = 99; 
  jMuMuHits   = -1; 
  jMuStations = 9696; 
  jMuPxHits   = 9696; 
  jMuTrkLayer = 9696;
}            
void HecTau::fill_tightMuons( )
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
void HecTau::softMuC::init()
{
 iMuInNormChi2 = 99;
 iMuInTrkLayer = -1; 
 jMuInNormChi2 = 99;
 jMuInTrkLayer = -1; 
 }            
void HecTau::fill_softMuons( )
{
 softMuC_.iMuInNormChi2   = CutMu[0][5];
 softMuC_.iMuInTrkLayer   = CutMu[0][6];
 softMuC_.jMuInNormChi2   = CutMu[1][5];
 softMuC_.jMuInTrkLayer   = CutMu[1][6];
}
void HecTau::tightKMuC::init()
{
  kNormChi2   = 99; 
  kMuMuHits   = -1; 
  kMuStations = -1; 
  kMuPxHits   = -1; 
  kMuTrkLayer = -1;
}            
void HecTau::fill_tightKMuons(double kNormChi2,double kMuMuHits,double kMuStations,double kMuPxHits,double kMuTrkLayer )
{
 tightKMuC_.kNormChi2	= kNormChi2;
 tightKMuC_.kMuMuHits	= kMuMuHits;
 tightKMuC_.kMuStations	= kMuStations;
 tightKMuC_.kMuPxHits	= kMuPxHits;
 tightKMuC_.kMuTrkLayer	= kMuTrkLayer;
}
void HecTau::softKMuC::init()
{
 kMuInNormChi2 = 99;
 kMuInTrkLayer = -1; 
}            
void HecTau::fill_softKMuons(double kMuInNormChi2,double kMuInTrkLayer )
{
 softKMuC_.kMuInNormChi2   = kMuInNormChi2;
 softKMuC_.kMuInTrkLayer   = kMuInTrkLayer;
}
//function for my Analizer
void HecTau::HecMuVar(int n, std::vector<pat::Muon>::const_iterator nMu, reco::TrackRef nMuTrk){
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
  MuVar[n][10] = nMu->combinedQuality().trkKink;                //--Chi squared
  MuVar[n][11] = 999;             //--for StandAloneMuon  (for Tracker and Global Muon after the MuonId)
  MuVar[n][12] = 999;             //--for StandAloneMuon  (for Tracker and Global Muon after the MuonId)

  //--Integer Variables
  nMuVar[n][0]  = nMu->charge();

  int muId=0, muSA=0, muGL=0, muTK=0;         //--Muon Id
  if( nMu->isStandAloneMuon() )muSA = 1;      //--2^0
  if( nMu->isGlobalMuon()     )muGL = 10;     //--2^1
  if( nMu->isTrackerMuon()    )muTK = 100;    //--2^2
  int imuId = muTK + muGL + muSA;                                   //--bynary number
  muId = muTK*pow(2,2)/100 + muGL*pow(2,1)/10 + muSA*pow(2,0)/1;    //--convert to decimal

  int iGL = 0;
  if(  nMu->isGlobalMuon()                          )iGL = 1;
  if(!(nMu->isGlobalMuon()) && nMu->isTrackerMuon() )iGL = 2;

    nMuVar[n][1]  = imuId;
    nMuVar[n][2]  = muId;     //1:--1  : S
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
bool HecTau::HecClosest(int Mode, FreeTrajectoryState iState, FreeTrajectoryState jState){
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
    dcaVar[Mode][0] = dca;          //--Flag  Closest Approach distance between muons
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
void HecTau::HecHltTrig(const edm::Event& iEvent){
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
   if(trigName=="HLT_Dimuon0_Jpsi_v1"        ) hlt_2mu0 = hltflag;
   if(trigName=="HLT_DoubleMu3_v2"           ) hlt_2mu3 = hltflag;
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

   if(trigName=="HLT_Mu0_Track0_Jpsi"   ) hlt_mu0trk0      = hltflag;
   if(trigName=="HLT_Mu3_Track0_Jpsi"   ) hlt_mu3trk0      = hltflag;
   if(trigName=="HLT_Mu0_TkMu0_Jpsi"    ) hlt_mu0trkmu0    = hltflag;
   if(trigName=="HLT_Mu3_TkMu0_Jpsi"    ) hlt_mu3trkmu0    = hltflag;
   if(trigName=="HLT_Mu0_TkMu0_OST_Jpsi") hlt_mu0trkmu0OST = hltflag;
   if(trigName=="HLT_Mu3_TkMu0_OST_Jpsi") hlt_mu3trkmu0OST = hltflag;

   if(trigName=="HLT_Mu0_TkMu0_OST_Jpsi_Tight_v1") hlt_mu0trkmu0OST_tight = hltflag;

   if(trigName=="HLT_L1MuOpen"            ) hlt_L1muOpen       = hltflag;
   if(trigName=="HLT_L1DoubleMuOpen"      ) hlt_L12muOpen      = hltflag;
   if(trigName=="HLT_L1DoubleMuOpen_Tight") hlt_L12muOpenTight = hltflag;
   if(trigName=="HLT_L2DoubleMu0_v4"      ) hlt_2mu0L2         = hltflag;
   if(trigName=="HLT_Dimuon7_PsiPrime_v3"         )hlt_DimuPsi2s  = hltflag;
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
void HecTau::HecHltMuTrig(int nMu, const pat::Muon* patMu){
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
      int upLMDimuon6p5  =!patMu->triggerObjectMatchesByPath("HLT_Dimuon6p5_LowMass_Displaced_v*").empty();

      int upLMDimuon7    =!patMu->triggerObjectMatchesByPath("HLT_Dimuon7_LowMass_Displaced_v*").empty();

      int upLMDoubleMu4  =!patMu->triggerObjectMatchesByPath("HLT_DoubleMu4_LowMass_Displaced_v*").empty();

      int upLMDoubleMu4p5=!patMu->triggerObjectMatchesByPath("HLT_DoubleMu4p5_LowMass_Displaced_v*").empty();

      int upLMDoubleMu5  =!patMu->triggerObjectMatchesByPath("HLT_DoubleMu5_LowMass_Displaced_v*").empty();

      nTrig[nMu][1] = upLMDimuon6p5   *1
                    + upLMDimuon7     *10
                    + upLMDoubleMu4   *100
                    + upLMDoubleMu4p5 *1000
                    + upLMDoubleMu5   *10000;

      if( MyPrint_ ){
          std::cout << "Mu+ Trigger : "<<nTrig[nMu][0]<<" nMu "<<std::endl;
          std::cout << "Mu Trigger mu3 "         <<upTrig2mu3            << std::endl;
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
void HecTau::HecCutMu(int m, std::vector<pat::Muon>::const_iterator mMu){
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
//define this as a plug-in
DEFINE_FWK_MODULE(HecTau);
