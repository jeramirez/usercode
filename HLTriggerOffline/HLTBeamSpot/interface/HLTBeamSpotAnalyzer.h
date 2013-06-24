#ifndef HLTriggerOffline_HLTBeamSpot_HLTBeamSpotAnalyzer_h
#define HLTriggerOffline_HLTBeamSpot_HLTBeamSpotAnalyzer_h

/**_________________________________________________________________
   class:   HLTBeamSpotAnalyzer.h.h
   package: RecoVertex/BeamSpotProducer
   


 author: Francisco Yumiceva, Fermilab (yumiceva@fnal.gov)

 version $Id: HLTBeamSpotAnalyzer.h,v 1.1 2008/09/15 15:50:51 bonnett Exp $

________________________________nedit ________________________________**/


// C++ standard
#include <string>
// CMS
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "RecoVertex/BeamSpotProducer/interface/BSTrkParameters.h"

// ROOT
#include "TFile.h"
#include "TTree.h"

#include<fstream>

class HLTBeamSpotAnalyzer : public edm::EDAnalyzer {
 public:
  explicit HLTBeamSpotAnalyzer(const edm::ParameterSet&);
  ~HLTBeamSpotAnalyzer();

 private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  std::string outputfilename_;
  std::string fasciifilename;
  std::ofstream fasciiFile;
  std::string fasciiFileName;

  TFile* file_;
  TTree* ftree_;
  TH1F * numTracks;
////////////mik//////////// 
  double ftrack;
  double beam_x0;
  double beam_y0;
  double beam_z0;
  double beam_sigmaZ;
  double beam_dxdz;
  double beam_dydz;
  double beam_cov[7][7];
////////////mik///////////  

  int    ftotalevents;
  int    cutTrack;
  double ftheta;
  double fpt;
  double feta;
  int    fcharge;
  double fchi2;
  double fndof;
  double fphi0;
  double fd0;
  double fsigmad0;
  double fz0;
  double fsigmaz0;
  unsigned int fnHit;
  unsigned int fnStripHit;
  unsigned int fnPixelHit;
  unsigned int fnTIBHit;
  unsigned int fnTOBHit;
  unsigned int fnTIDHit;
  unsigned int fnTECHit;
  unsigned int fnPXBHit;
  unsigned int fnPXFHit;
  double fd0phi_chi2;
  double fd0phi_d0;
  double fcov[7][7];
  /////////mik///////////
  std::map<int,std::vector< BSTrkParameters > > fBSvector;
  ////////mik////////////
  std::string ckfSeedProducerLabel_;
  std::string ckfTrackCandidateProducerLabel_;
  std::string ckfTrackProducerLabel_;

  unsigned int sameNumberOfTracks;

  float fptmin;
  int fmaxNtracks;  
  int  nSHit;  

  bool write2DB_;
  bool runallfitters_;
  int ftotal_tracks;
  double inputBeamWidth_;
  
};

#endif
