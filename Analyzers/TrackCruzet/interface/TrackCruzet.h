#ifndef TrackCruzet_h
#define TrackCruzet_h

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/CurrentProcessingContext.h"

#include "TFile.h"
#include "TH1D.h"
#include "TH1I.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TTree.h"

// system include files
#include <iostream>
#include <memory>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include<fstream>
 
/* Collaborating Class Declarations */
class Propagator;
 
//#include "TrackingTools/GeomPropagators/interface/Propagator.h"
//#include "HLTrigger/HLTcore/interface/HLTFilter.h"

class TrackAssociatorBase;

class TrackCruzet : public edm::EDAnalyzer
{
public:
    explicit TrackCruzet(const edm::ParameterSet&);
    ~TrackCruzet();

    /* Operations */
    bool filter(const edm::Event&, const edm::EventSetup& , double radius, double max_z);
//     bool filterloose(const edm::Event&, const edm::EventSetup& );

private:
    virtual void beginJob(const edm::EventSetup&) ;
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;
///////////////////////////////////////////////////////////////////

    std::string theMuCollectionLabel; // label of muons
    std::string thePropagatorName; // name of propagator to be used

//   // radius of cylinder
    double theRadius_SiTrk;
    double theRadius_Pixel;
    double theRadius_IP;
//   // half lenght of cylinder
    double theMaxZ_SiTrk;
    double theMaxZ_Pixel;
    double theMaxZ_IP;

    double posX_;
    double posY_;
    double posZ_;

    int nmuon_;

    int nday_;
    int minRun_;
    int maxRun_;
    
    int yer_;
    int mon_;
    int day_;
    int hor_;
    int min_;
    int sec_;

    mutable Propagator* thePropagator;
    std::ofstream fasciiFile;
    std::string fasciiFileName;
    std::ofstream fileout;
    std::ofstream fileNoPixel;
    std::ofstream fileflagloosmuon;
    std::ofstream fileloosemuonlist;
    std::ofstream filefussy;
    std::ofstream fileshaft;
    std::ofstream fileshaft1;
    std::ofstream date_run_event;
    std::ofstream testmuon;

    edm::InputTag  ctfTrackCollectionTag_;       //combinatorial track finder (cft)
    edm::InputTag  rsTrackCollectionTag_;        //Cosmic Track Finder
    edm::InputTag  cosmicTFTrackCollectionTag_;  //Road Search

//  TFile* rootFile_;
//  TTree* rootTree_;
//  TFile* file_;

    std::map<std::string, TH1D*> Trk_;
    std::map<std::string, TH2D*> Trk2D_;
    std::map<int, std::string> paramt_;
    std::map<int, std::string> cut_;
    std::map<int, std::string> collection_;

    TH1I* Flag;
    TH1I* Flagloose;
    TH1I* Flagloosemuon1;
    TH1I* runNumber;
    TH1I* date_data;
    TH1D* nclusters_total;
    std::map<std::string, TH1I*> MuCollect1I_;
    std::map<std::string, TH1D*> MuCollect1D_;
    std::map<std::string, TH1F*> MuCollect1F_;
    std::map<std::string, TH2D*> MuCollect2D_;
    std::map<int, std::string> cutsta_;
    
    TH1F* Run_Good_muon;
    TH1F* Run_All_muon;
    TH1F* Run_All;
    TH1I* MuCollect_NRhits;
    TH1I* MuCollect_NDThits;
    TH1I* MuCollect_NCSChits;
    TH1I* MuCollect_NRPChits;
  
    TH1D*MuCollect_innerPositionX;
    TH1D*MuCollect_innerPositionY;
    TH1D*MuCollect_innerPositionZ;
    TH1D*MuCollect_outerPositionX;
    TH1D*MuCollect_outerPositionY;
    TH1D*MuCollect_outerPositionZ;
  
    TH1D*MuCollect_innerMomentumX;
    TH1D*MuCollect_innerMomentumY;
    TH1D*MuCollect_innerMomentumZ;
    TH1D*MuCollect_outerMomentumX;
    TH1D*MuCollect_outerMomentumY;
    TH1D*MuCollect_outerMomentumZ;

//    TH1D*MuCollect_delta_p;
//    TH1D*MuCollect_delta_p_muon;

    TH1D* MuCollect_vx;
    TH1D* MuCollect_vy;
    TH1D* MuCollect_vz;
    TH1D* MuCollect_recHitsPositionX;
    TH1D* MuCollect_recHitsPositionY;
    TH1D* MuCollect_recHitsPositionZ;
    TH1D* MuCollect_recHitsGlobalPositionX;
    TH1D* MuCollect_recHitsGlobalPositionY;
    TH1D* MuCollect_recHitsGlobalPositionZ;
    
    //TH2D* eff_muon;

//     TH2D* MuCollect_recHitsGlobalPositionYvsX;
//     TH2D* MuCollect_recHitsGlobalPositionXvsZ;
//     TH2D* MuCollect_recHitsGlobalPositionYvsZ;
    TH3D* MuCollect_recHitsGlobalPositionXYZ;
    
    TH2D* MuCollect_surfaceXZ;
    TH2D* MuCollect_shaftXZ;
    TH2D* MuCollect_shaftXZ1;
    TH1I* MuCollect_nRecHits_shaft;
    TH1D* MuCollect_normalizedchi2_shaft;
       
    TH1I* MuCollect_time_date;	
    TH1I* MuCollect_time_date0;	
    TH1I* MuCollect_time_date1;	
    TH1I* SiTrk_date;	
    TH1I* SiTrk_date1;	
    TH1I* MuCollect_time_hor;	
    TH1I* MuCollect_time_day;	
    	
    TH1D* MuCollect_test;	
    TH1D* MuCollect_Flagtest1;	
     
};

#endif
