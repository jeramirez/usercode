#include <memory>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <cmath>

#include <TFile.h>
#include <TH1.h>
#include <Math/GenVector/VectorUtil.h>

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/Registry.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "SimDataFormats/JetMatching/interface/JetFlavour.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourMatching.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/JetTracksAssociation.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "HLTriggerOffline/BJet/interface/JetPlots.h"
#include "HLTriggerOffline/BJet/interface/OfflineJetPlots.h"
#include "HLTriggerOffline/BJet/interface/FlavouredJetPlots.h"
#include "HLTriggerOffline/BJet/interface/VertexPlots.h"
#include "HLTriggerOffline/BJet/interface/RatePlots.h"

////////////////////////////////////////////////////////////////////////////////
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidate.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidateCollection.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"


#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "HLTriggerOffline/HLTBeamSpot/interface/HLTBeamSpotAnalyzer.h"
#include "RecoVertex/BeamSpotProducer/interface/BSFitter.h"
#include "CondFormats/BeamSpotObjects/interface/BeamSpotObjects.h"

#include "TMath.h"

////////////////////////////////////////////////////////////////////////////////
// force LogDebug statements to be generated
#define ML_DEBUG
/// http://cmslxr.fnal.gov/lxr/source/HLTrigger/Configuration/python/HLT_2E30_cff.py?v=CMSSW_2_1_0#6534
//hltBLifetimeRegionalCtfWithMaterialTracks


HLTBeamSpotAnalyzer::HLTBeamSpotAnalyzer(const edm::ParameterSet& iConfig) {

	outputfilename_ = iConfig.getUntrackedParameter<std::string>("OutputFileName");
	
  file_ = TFile::Open(outputfilename_.c_str(),"RECREATE");

  ftree_ = new TTree("mytree","mytree");
  ftree_->AutoSave();
  
  ftree_->Branch("pt",&fpt,"fpt/D");
  ftree_->Branch("d0",&fd0,"fd0/D");
  ftree_->Branch("sigmad0",&fsigmad0,"fsigmad0/D");
  ftree_->Branch("phi0",&fphi0,"fphi0/D");
  ftree_->Branch("z0",&fz0,"fz0/D");
  ftree_->Branch("sigmaz0",&fsigmaz0,"fsigmaz0/D");
  ftree_->Branch("theta",&ftheta,"ftheta/D");
  ftree_->Branch("eta",&feta,"feta/D");
  ftree_->Branch("charge",&fcharge,"fcharge/I");
  ftree_->Branch("chi2",&fchi2,"fchi2/D");
  ftree_->Branch("ndof",&fndof,"fndof/D");
  ftree_->Branch("nHit",&fnHit,"fnHit/i");
  ftree_->Branch("nStripHit",&fnStripHit,"fnStripHit/i");
  ftree_->Branch("nPixelHit",&fnPixelHit,"fnPixelHit/i");
  ftree_->Branch("nTIBHit",&fnTIBHit,"fnTIBHit/i");
  ftree_->Branch("nTOBHit",&fnTOBHit,"fnTOBHit/i");
  ftree_->Branch("nTIDHit",&fnTIDHit,"fnTIDHit/i");
  ftree_->Branch("nTECHit",&fnTECHit,"fnTECHit/i");
  ftree_->Branch("nPXBHit",&fnPXBHit,"fnPXBHit/i");
  ftree_->Branch("nPXFHit",&fnPXFHit,"fnPXFHit/i");
  ftree_->Branch("cov",&fcov,"fcov[7][7]/D");

////////////mik//////////// 
  numTracks = new TH1F("numTracks","Number of Tracks",100,0,100);
////////////mik//////////// 
   
  fBSvector.clear();

  //dump to file
  fasciiFileName = outputfilename_.replace(outputfilename_.size()-4,outputfilename_.size(),"txt");
  fasciiFile.open(fasciiFileName.c_str());

  
  // get parameter
 
  //ckfSeedProducerLabel_ = iConfig.getUntrackedParameter<std::string>("ckfSeedProducerLabel");
  //ckfTrackCandidateProducerLabel_ = iConfig.getUntrackedParameter<std::string>("ckfTrackCandidateProducerLabel");
  //ckfTrackProducerLabel_ = iConfig.getUntrackedParameter<std::string>("ckfTrackProducerLabel");

  //sameNumberOfTracks = iConfig.getUntrackedParameter<unsigned int>("sameNumberOfTracks");

  ////////////mik//////////// 
  nSHit = iConfig.getParameter<edm::ParameterSet>("BSAnalyzerParameters").getParameter<int>("minSHit");
  ////////////mik//////////// 

  fptmin = iConfig.getParameter<edm::ParameterSet>("BSAnalyzerParameters").getParameter<double>("MinimumPt");
  fmaxNtracks = iConfig.getParameter<edm::ParameterSet>("BSAnalyzerParameters").getParameter<int>("MaximumNtracks");
  ckfTrackProducerLabel_ = iConfig.getParameter<edm::ParameterSet>("BSAnalyzerParameters").getUntrackedParameter<std::string>("TrackCollection");
   
  write2DB_ = iConfig.getParameter<edm::ParameterSet>("BSAnalyzerParameters").getParameter<bool>("WriteToDB");
  runallfitters_ = iConfig.getParameter<edm::ParameterSet>("BSAnalyzerParameters").getParameter<bool>("RunAllFitters");
  inputBeamWidth_ = iConfig.getParameter<edm::ParameterSet>("BSAnalyzerParameters").getUntrackedParameter<double>("InputBeamWidth",-1.);
  
  ftotal_tracks = 0;
  ftotalevents = 0;
////////////mik//////////// 
  cutTrack =0;
////////////mik//////////// 

}

HLTBeamSpotAnalyzer::~HLTBeamSpotAnalyzer() 
{
}

void HLTBeamSpotAnalyzer::beginJob(const edm::EventSetup & setup) 
{

}


void HLTBeamSpotAnalyzer::analyze(const edm::Event & iEvent, const edm::EventSetup & setup) 
{
	
	ftree_->SetBranchAddress("theta",&ftheta);
	ftree_->SetBranchAddress("pt",&fpt);
	ftree_->SetBranchAddress("eta",&feta);
	ftree_->SetBranchAddress("charge",&fcharge);
	ftree_->SetBranchAddress("chi2",&fchi2);
	ftree_->SetBranchAddress("ndof",&fndof);
	ftree_->SetBranchAddress("d0",&fd0);
	ftree_->SetBranchAddress("sigmad0",&fsigmad0);
	ftree_->SetBranchAddress("phi0",&fphi0);
	ftree_->SetBranchAddress("z0",&fz0);
	ftree_->SetBranchAddress("sigmaz0",&fsigmaz0);
	ftree_->SetBranchAddress("nHit",&fnHit);
	ftree_->SetBranchAddress("nStripHit",&fnStripHit);
	ftree_->SetBranchAddress("nPixelHit",&fnPixelHit);
	ftree_->SetBranchAddress("nTIBHit",&fnTIBHit);
	ftree_->SetBranchAddress("nTOBHit",&fnTOBHit);
	ftree_->SetBranchAddress("nTIDHit",&fnTIDHit);
	ftree_->SetBranchAddress("nTECHit",&fnTECHit);
	ftree_->SetBranchAddress("nPXBHit",&fnPXBHit);
	ftree_->SetBranchAddress("nPXFHit",&fnPXFHit);
	ftree_->SetBranchAddress("cov",&fcov);
  
	// get collections
	
	
	edm::Handle<reco::TrackCollection> ckfTrackCollectionHandle;
	//iEvent.getByLabel(TkTag,tkCollection);
	iEvent.getByLabel(ckfTrackProducerLabel_,ckfTrackCollectionHandle);
	const reco::TrackCollection *ckfTrackCollection = ckfTrackCollectionHandle.product();

        numTracks->Fill(ckfTrackCollection->size());

	
	for ( reco::TrackCollection::const_iterator track = ckfTrackCollection->begin();
		  track != ckfTrackCollection->end();
		  ++track ) {

		
	  fpt = track->pt();
	  //std::cout << "pt= "<< track->pt() << std::endl;
	  //std::cout << "eta= "<< track->eta() << std::endl;
	  
	  feta = track->eta();
	  fphi0 = track->momentum().phi();
	  fcharge = track->charge();
	  fchi2 = track->chi2();
	  fndof = track->ndof();
	  
	  //fsigmaphi0 = track->phi0Error();
	  fd0 = track->d0();
	  fsigmad0 = track->d0Error();
	  fz0 = track->dz();
	  fsigmaz0 = track->dzError();
	  ftheta = track->theta();

	  for (int i=0; i<5; ++i) {
		  for (int j=0; j<5; ++j) {
			  fcov[i][j] = track->covariance(i,j);
		  }
	  }
	  
	  // loop over hits in tracks, count
	  fnHit      = 0;
	  fnStripHit = 0;
	  fnPixelHit = 0;
	  fnTIBHit   = 0;
	  fnTOBHit   = 0;
	  fnTIDHit   = 0;
	  fnTECHit   = 0;
	  fnPXBHit   = 0;
	  fnPXFHit   = 0;

	  for ( trackingRecHit_iterator recHit = track->recHitsBegin();
			recHit != track->recHitsEnd();
			++ recHit ) {

		  
		  ++fnHit;
		  //std::cout << "fnHit="<< fnHit << std::endl;
		  
		  DetId id((*recHit)->geographicalId());

		  if ( (unsigned int)id.subdetId() == StripSubdetector::TIB ) {
			  ++fnStripHit;
			  ++fnTIBHit;
		  } else if ( (unsigned int)id.subdetId() == StripSubdetector::TOB ) {
			  ++fnStripHit;
			  ++fnTOBHit;
		  } else if ( (unsigned int)id.subdetId() == StripSubdetector::TID ) {
			  ++fnStripHit;
			  ++fnTIDHit;
		  } else if ( (unsigned int)id.subdetId() == StripSubdetector::TEC ) {
			  ++fnStripHit;
			  ++fnTECHit;
		  } else if ( (unsigned int)id.subdetId() == PixelSubdetector::PixelBarrel ) {
			  ++fnPixelHit;
			  ++fnPXBHit;
		  } else if ( (unsigned int)id.subdetId() == PixelSubdetector::PixelEndcap ) {
			  ++fnPixelHit;
			  ++fnPXFHit;
		  }
	  }

	  ftree_->Fill();

	  ftotal_tracks++;
////////////mik//////////// 

/*	std::cout << "fnStripHit = "<< fnStripHit << std::endl;
	std::cout << "fnPixelHit = "<< fnPixelHit << std::endl;
	std::cout << "fchi2      = "<< fchi2 << std::endl;
	std::cout << "fndof      = "<< fndof << std::endl;
	std::cout << "chi2/ndof  = "<< fchi2/fndof << std::endl;
	std::cout << "fpt        = "<< fpt << std::endl;
	std::cout << "abs(fd0)   = "<< std::abs(fd0) << std::endl;
*/		 
////////////mik//////////// 
          // track quality
////////////mik//////////// 
	  if (fnStripHit >= nSHit && fnPixelHit >= 2 && fchi2/fndof<5 && fpt>2 && std::abs(fd0)<0.9) {
//	  if (fnStripHit >= 8 && fnPixelHit >= 2 && fchi2/fndof<5 && fpt>2 && std::abs(fd0)<0.9) {
//		  std::cout << "track quality track quality track quality track quality track quality"<< std::endl;
 
 		  //fchi2/fndof<5 && fpt>4 && std::abs(fd0)<0.1 && TMath::Prob(fchi2,((int)fndof))>0.02 ) {
		  //if 
	cutTrack = ((ftotal_tracks/500) - 1);
	for (int k=0; k<15; k++ ){	 
		 if( fBSvector[k].size()<= fmaxNtracks+k*500){
		  fBSvector[k].push_back(BSTrkParameters(fz0,fsigmaz0,fd0,fsigmad0,fphi0,fpt,0.,0.));
		  }
		}  
	  }
////////////mik//////////// 
	  
    
	}
  
	ftotalevents++;
}

void HLTBeamSpotAnalyzer::endJob()
{

////////////mik//////////// 
	file_->cd();
	TTree *newtree = new TTree("mytreebeam","mytreebeam");
	newtree->Branch("ftrack",&ftrack,"ftrack/D");
	newtree->Branch("beam_x0",&beam_x0,"beam_x0/D");
	newtree->Branch("beam_y0",&beam_y0,"beam_y0/D");
	newtree->Branch("beam_z0",&beam_z0,"beam_z0/D");
	newtree->Branch("beam_sigmaZ",&beam_sigmaZ,"beam_sigmaZ/D");
	newtree->Branch("beam_dxdz",&beam_dxdz,"beam_dxdz/D");
	newtree->Branch("beam_dydz",&beam_dydz,"beam_dydz/D");
	newtree->Branch("beam_cov",&beam_cov,"beam_cov[7][7]/D");
	newtree->SetBranchAddress("ftrack",&ftrack);
	newtree->SetBranchAddress("beam_x0",&beam_x0);
	newtree->SetBranchAddress("beam_y0",&beam_y0);
	newtree->SetBranchAddress("beam_z0",&beam_z0);
	newtree->SetBranchAddress("beam_sigmaZ",&beam_sigmaZ);
	newtree->SetBranchAddress("beam_dxdz",&beam_dxdz);
	newtree->SetBranchAddress("beam_dydz",&beam_dydz);
	newtree->SetBranchAddress("beam_cov",&beam_cov);       
////////////mik//////////// 


	for (int k=0; k<15; k++ ){	 
	std::cout << "cutTrack" << cutTrack <<std::endl;
	std::cout << "\n-------------------------------------\n" << std::endl;
	std::cout << "\n Total number of events processed: "<< ftotalevents << std::endl;
	std::cout << "\n-------------------------------------\n\n" << std::endl;
	std::cout << " calculating beam spot..." << std::endl;
	std::cout << " we will use " << fBSvector[k].size() << " good tracks out of " << ftotal_tracks << std::endl;

	// default fit to extract beam spot info
	BSFitter *myalgo1 = new BSFitter( fBSvector[k] );
	reco::BeamSpot beam_default = myalgo1->Fit();
	std::cout << "\n RESULTS OF DEFAULT FIT:" << std::endl;
	std::cout << beam_default << std::endl;

	// dump to file
	fasciiFile << "ntrack " << fBSvector[k].size() << std::endl;
	fasciiFile << "X " << beam_default.x0() << std::endl;
	fasciiFile << "Y " << beam_default.y0() << std::endl;
	fasciiFile << "Z " << beam_default.z0() << std::endl;
	fasciiFile << "sigmaZ " << beam_default.sigmaZ() << std::endl;
	fasciiFile << "dxdz " << beam_default.dxdz() << std::endl;
	fasciiFile << "dydz " << beam_default.dydz() << std::endl;
	if (inputBeamWidth_ > 0 ) {
		fasciiFile << "BeamWidth " << inputBeamWidth_ << std::endl;
	} else {
		fasciiFile << "BeamWidth " << beam_default.BeamWidth() << std::endl;
	}
	

	ftrack = fBSvector[k].size();
	beam_x0 = beam_default.x0();
	beam_y0 = beam_default.y0();
	beam_z0 = beam_default.z0();
	beam_sigmaZ = beam_default.sigmaZ();
	beam_dxdz = beam_default.dxdz();
	beam_dydz = beam_default.dydz();

////////////mik//////////// 
	
	for (int i = 0; i<6; ++i) {
		fasciiFile << "Cov("<<i<<",j) ";
		for (int j=0; j<7; ++j) {
			fasciiFile << beam_default.covariance(i,j) << " ";
	                beam_cov[i][j] = beam_default.covariance(i,j);
		}
		fasciiFile << std::endl;
	}
	// beam width error
	if (inputBeamWidth_ > 0 ) {
		fasciiFile << "Cov(6,j) 0 0 0 0 0 0 " << pow(2.e-4,2) << std::endl;
	} else {
		fasciiFile << "Cov(6,j) 0 0 0 0 0 0 " << beam_default.covariance(6,6) << std::endl;
	}
	newtree->Fill();
	delete myalgo1;

	}
	newtree->Write();
	
/*	BSFitter *myalgo = new BSFitter( fBSvector[0] );
	reco::BeamSpot beam_default = myalgo->Fit();


	if (runallfitters_) {
	
	  
	// add new branches
	std::cout << " add new branches to output file " << std::endl;
	beam_default = myalgo->Fit_d0phi();
	file_->cd();
	TTree *newtree = new TTree("mytreecorr","mytreecorr");
	newtree->Branch("d0phi_chi2",&fd0phi_chi2,"fd0phi_chi2/D");
	newtree->Branch("d0phi_d0",&fd0phi_d0,"fd0phi_d0/D");
	newtree->SetBranchAddress("d0phi_chi2",&fd0phi_chi2);
	newtree->SetBranchAddress("d0phi_d0",&fd0phi_d0);
	std::vector<BSTrkParameters>  tmpvector = myalgo->GetData();
	
	std::vector<BSTrkParameters>::iterator iparam = tmpvector.begin();
	for( iparam = tmpvector.begin() ;
		 iparam != tmpvector.end() ; ++iparam) {
		fd0phi_chi2 = iparam->d0phi_chi2();
		fd0phi_d0   = iparam->d0phi_d0();
		newtree->Fill();
	}
	newtree->Write();

	// iterative
	std::cout << " d0-phi Iterative:" << std::endl;
	BSFitter *myitealgo = new BSFitter( fBSvector[0] );
	myitealgo->Setd0Cut_d0phi(4.0);
	reco::BeamSpot beam_ite = myitealgo->Fit_ited0phi();
	std::cout << beam_ite << std::endl;

	
	std::cout << "\n Now run tests of the different fits\n";
	// from here are just tests
	std::string fit_type = "chi2";
	myalgo->SetFitVariable(std::string("z"));
	myalgo->SetFitType(std::string("chi2"));
	reco::BeamSpot beam_fit_z_chi2 = myalgo->Fit();
	std::cout << " z Chi2 Fit ONLY:" << std::endl;
	std::cout << beam_fit_z_chi2 << std::endl;
	
	
	fit_type = "combined";
	myalgo->SetFitVariable("z");
	myalgo->SetFitType("combined");
	reco::BeamSpot beam_fit_z_lh = myalgo->Fit();
	std::cout << " z Log-Likelihood Fit ONLY:" << std::endl;
	std::cout << beam_fit_z_lh << std::endl;

	
	myalgo->SetFitVariable("d");
	myalgo->SetFitType("d0phi");
	reco::BeamSpot beam_fit_dphi = myalgo->Fit();
	std::cout << " d0-phi0 Fit ONLY:" << std::endl;
	std::cout << beam_fit_dphi << std::endl;

		
	myalgo->SetFitVariable(std::string("d*z"));
	myalgo->SetFitType(std::string("likelihood"));
	reco::BeamSpot beam_fit_dz_lh = myalgo->Fit();
	std::cout << " Log-Likelihood Fit:" << std::endl;
	std::cout << beam_fit_dz_lh << std::endl;

	
	myalgo->SetFitVariable(std::string("d*z"));
	myalgo->SetFitType(std::string("resolution"));
	reco::BeamSpot beam_fit_dresz_lh = myalgo->Fit();
	std::cout << " IP Resolution Fit" << std::endl;
	std::cout << beam_fit_dresz_lh << std::endl;

	std::cout << "c0 = " << myalgo->GetResPar0() << " +- " << myalgo->GetResPar0Err() << std::endl;
	std::cout << "c1 = " << myalgo->GetResPar1() << " +- " << myalgo->GetResPar1Err() << std::endl;

	}
*/
	// let's close everything
    file_->cd();
	ftree_->Write();
	numTracks->Write();
	file_->Close();
    		
	std::cout << "[HLTBeamSpotAnalyzer] endJob done \n" << std::endl;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(HLTBeamSpotAnalyzer);
