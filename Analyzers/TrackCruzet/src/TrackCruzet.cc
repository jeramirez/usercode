#include "Analyzers/TrackCruzet/interface/TrackCruzet.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

//TrackAssociation
#include "DataFormats/TrackCandidate/interface/TrackCandidate.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidateCollection.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonTime.h"

//TrackInfo
#include "AnalysisDataFormats/TrackInfo/interface/TrackInfo.h"
#include "AnalysisDataFormats/TrackInfo/interface/TrackInfoEnum.h"
#include "AnalysisDataFormats/TrackInfo/interface/TrackInfoTrackAssociation.h"

//TFileService headers
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

// system include files
#include <iostream>
#include <memory>
#include <string>
// #include <sstream>
// #include <vector>
// #include <map>
#include "TMath.h"
#include <time.h>
#include "TStyle.h"
#include "TLatex.h"

//FOR CLUSTERINFO
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"
#include "DataFormats/SiStripDigi/interface/SiStripDigi.h"
#include "DataFormats/SiStripCluster/interface/SiStripCluster.h"
//#include "DataFormats/SiStripCluster/interface/SiStripClusterInfo.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
///////////////////////////////////////////////////////
/* Collaborating Class Header */
#include "DataFormats/GeometrySurface/interface/Cylinder.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"

#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"

#include "DataFormats/DTRecHit/interface/DTRecSegment4D.h"
#include "DataFormats/TrackingRecHit/interface/RecSegment.h"

#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"
#include "DataFormats/DTRecHit/interface/DTRecHitCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCRecHit2DCollection.h"

#include "DataFormats/GeometrySurface/interface/GloballyPositioned.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"

/* C++ Headers */
using namespace std;
using namespace edm;


class TrackAssociatorByHits;
class TrackerHitAssociator;

TrackCruzet::TrackCruzet(const edm::ParameterSet& iConfig)
{
    //file_ = new TFile(iConfig.getUntrackedParameter<std::string>("OutputFileName").c_str(),"RECREATE");
    ctfTrackCollectionTag_ = iConfig.getParameter<edm::InputTag>("ctfTrackCollection");
    rsTrackCollectionTag_ = iConfig.getParameter<edm::InputTag>("rsTrackCollection");
    cosmicTFTrackCollectionTag_ = iConfig.getParameter<edm::InputTag>("cosmicTFTrackCollection");

    theMuCollectionLabel = iConfig.getParameter<string>("MuCollectLabel");

    thePropagatorName = iConfig.getParameter<std::string>("PropagatorName");
    thePropagator = 0;

    theRadius_IP = iConfig.getParameter<double>("radius_IP"); // cyl's radius (cm)
    theMaxZ_IP = iConfig.getParameter<double>("maxZ_IP"); // cyl's half lenght (cm)

    theRadius_Pixel = iConfig.getParameter<double>("radius_Pixel"); // cyl's radius (cm)
    theMaxZ_Pixel = iConfig.getParameter<double>("maxZ_Pixel"); // cyl's half lenght (cm)

    theRadius_SiTrk = iConfig.getParameter<double>("radius_SiTrk"); // cyl's radius (cm)
    theMaxZ_SiTrk = iConfig.getParameter<double>("maxZ_SiTrk"); // cyl's half lenght (cm)

    posX_ = iConfig.getParameter<double>("position_x"); // Region X (cm)
    posY_ = iConfig.getParameter<double>("position_y"); // Region Y (cm)
    posZ_ = iConfig.getParameter<double>("position_z"); // Region Z (cm)

    nmuon_ = iConfig.getParameter<unsigned int>("nmuon"); 

    nday_  = iConfig.getParameter<unsigned int>("days");
    
    maxRun_ = iConfig.getParameter<unsigned int>("maxRun");
    minRun_  = iConfig.getParameter<unsigned int>("minRun");

    yer_  = iConfig.getParameter<unsigned int>("date_year");
    mon_  = iConfig.getParameter<unsigned int>("date_month");
    day_  = iConfig.getParameter<unsigned int>("date_day");
    hor_  = iConfig.getParameter<unsigned int>("date_hour");
    min_  = iConfig.getParameter<unsigned int>("date_min");
    sec_  = iConfig.getParameter<unsigned int>("date_sec");

    LogDebug("HLTMuonPointing") << " SALabel : " << theMuCollectionLabel
    << " Radius : " << theRadius_SiTrk
    << " Half lenght : " << theMaxZ_SiTrk;

    fasciiFileName = iConfig.getParameter<string>("ascciFileName");
    fasciiFile.open(fasciiFileName.c_str());
    fileout.open("momentum.txt");
    fileNoPixel.open("NoPixelTrack.txt");
    filefussy.open("fuzzy.txt");
    fileflagloosmuon.open("flagloosmuon.txt");
    fileloosemuonlist.open("muon_filterlist.txt");
    fileshaft.open("shaft.txt");
    fileshaft1.open("flagloose_false.txt");
    date_run_event.open("date_run_event.txt");
    testmuon.open("length1.txt");
//////////////////////////////////////////////////////////////////////////////////////////

    edm::Service<TFileService> fs;

    runNumber  = fs->make<TH1I>("runNumber", "Run Number", 102000, 1000.5, 103000.5);
    date_data  = fs->make<TH1I>("date_data", "date", nday_*86400, 0.0, nday_*86400);        
    Flag       = fs->make<TH1I>("Flag", "Flag per Track", 5, 0, 5);
    Flagloose  = fs->make<TH1I>("Flagloose", "Flag with STAmuon PositionYZ_x RecHits > 20", 5, 0, 5);
    Flagloosemuon1  = fs->make<TH1I>("Flagloosemuon1", "Flag with Flagloose & Numb muon >1", 5, 0, 5);
    nclusters_total  = fs->make<TH1D>("nclusters_total","Total Number of Clusters per Event",10000,0.,10000.);

    collection_[0] = "ctf_";         //combinatorial track finder (CTF)
    collection_[1] = "cosmicTF_";    //Cosmic Track Finder
    collection_[2] = "rs_";          //Road Search
    cut_[0] = "";
    cut_[1] = "SiTrack_";          //flag1
    cut_[2] = "Pixel_";            //flag2
    cut_[3] = "IP_";               //flag3
    cut_[4] = "SiTrack_loose_";    //flagloose1
    cut_[5] = "Pixel_loose_";      //flagloose2
    cut_[6] = "IP_loose_";         //flagloose3

    cut_[7] = "SiTrack_loosemuon1_";    //flagloosemuon11
    cut_[8] = "Pixel_loosemuon1_";      //flagloosemuon12
    cut_[9] = "IP_loosemuon1_";         //flagloosemuon13

    cut_[10] = "No_SiTrack_loose_";    //flagloose1
    cut_[11] = "No_Pixel_loose_";      //flagloose2
    cut_[12] = "No_IP_loose_";         //flagloose3

    for (int j = 0; j < 3; ++j)
    {
        std::string collect = collection_[j];

        TFileDirectory dir = fs->mkdir(collect);

        for (int i = 0; i < 13; ++i)
        {
            std::string cut = cut_[i];
            TFileDirectory subDir = dir.mkdir( cut );
            Trk_[collection_[j] + cut_[i] + "NumTracks"]        = subDir.make<TH1D>(TString(collection_[j] + cut_[i] + "NumTracks"), TString("Number Of Tracks "+collection_[j]), 10, -0.5, 9.5);
            Trk_[collection_[j] + cut_[i] + "NumTracks"]->Sumw2();
            Trk_[collection_[j] + cut_[i] + "date"]        = subDir.make<TH1D>(TString(collection_[j] + cut_[i] + "date"), TString("date "+collection_[j]), nday_*86400, 0.0, nday_*86400);
           // Trk_[collection_[j] + cut_[i] + "date"]->Sumw2();
            Trk_[collection_[j] + cut_[i] + "NumRecHitsPerTrk"] = subDir.make<TH1D>(TString(collection_[j] + cut_[i] + "NumRecHitsPerTrk"), TString("Number Of RecHits Per Track "+collection_[j]), 51, -0.5, 51.5);
            Trk_[collection_[j] + cut_[i] + "NumRecHitsPerTrk"]->Sumw2(); 
            Trk_[collection_[j] + cut_[i] + "Run"] = subDir.make<TH1D>(TString(collection_[j] + cut_[i] + "Run"), TString("Run "+collection_[j]), maxRun_ - minRun_, minRun_, maxRun_);
            Trk_[collection_[j] + cut_[i] + "Run"]->Sumw2(); 

            Trk_[collection_[j] + cut_[i] + "d0"]               = subDir.make<TH1D>(TString(collection_[j] + cut_[i] + "d0"), TString("Impact Parameter (d0) "+collection_[j]), 100, -200, 200);
            Trk_[collection_[j] + cut_[i] + "d0"]->Sumw2(); 
            Trk_[collection_[j] + cut_[i] + "dxy"]              = subDir.make<TH1D>(TString(collection_[j] + cut_[i] + "dxy"),TString("Track dxy "+collection_[j]), 100, -200, 200);
            Trk_[collection_[j] + cut_[i] + "dxy"]->Sumw2();               
            Trk_[collection_[j] + cut_[i] + "dz"]               = subDir.make<TH1D>(TString(collection_[j] + cut_[i] + "dz"), TString("Track dz "+collection_[j]), 100, -500, 500);
            Trk_[collection_[j] + cut_[i] + "dz"]->Sumw2();               
            Trk_[collection_[j] + cut_[i] + "dsz"]              = subDir.make<TH1D>(TString(collection_[j] + cut_[i] + "dsz"),TString("Track dsz "+collection_[j]), 100, -300, 300);
            Trk_[collection_[j] + cut_[i] + "dsz"]->Sumw2();
                           
            Trk_[collection_[j] + cut_[i] + "pt"]               = subDir.make<TH1D>(TString(collection_[j] + cut_[i] + "pt"),      TString("Track p_{T} (GeV) "+collection_[j]), 100, -0.5, 10);
            Trk_[collection_[j] + cut_[i] + "pt"]->Sumw2();               
            Trk_[collection_[j] + cut_[i] + "px"]               = subDir.make<TH1D>(TString(collection_[j] + cut_[i] + "px"),      TString("Track p_{x} (GeV) "+collection_[j]), 100, -10, 10);
            Trk_[collection_[j] + cut_[i] + "px"]->Sumw2();                  
            Trk_[collection_[j] + cut_[i] + "py"]               = subDir.make<TH1D>(TString(collection_[j] + cut_[i] + "py"),      TString("Track p_{y} (GeV) "+collection_[j]), 100, -10, 10);
            Trk_[collection_[j] + cut_[i] + "py"]->Sumw2();               
            Trk_[collection_[j] + cut_[i] + "pz"]               = subDir.make<TH1D>(TString(collection_[j] + cut_[i] + "pz"),      TString("Track p_{z} (GeV) "+collection_[j]), 100, -10, 10);
            Trk_[collection_[j] + cut_[i] + "pz"]->Sumw2();               
            Trk_[collection_[j] + cut_[i] + "p"]                = subDir.make<TH1D>(TString(collection_[j] + cut_[i] + "p"),       TString("Track p (GeV) "+collection_[j]), 100, -1, 10);
            Trk_[collection_[j] + cut_[i] + "p"]->Sumw2();               
            Trk_[collection_[j] + cut_[i] + "eta"]              = subDir.make<TH1D>(TString(collection_[j] + cut_[i] + "eta"),     TString("Track Eta "+collection_[j]), 120, -3.0, 3.0);
            Trk_[collection_[j] + cut_[i] + "eta"]->Sumw2();               
            Trk_[collection_[j] + cut_[i] + "phi"]              = subDir.make<TH1D>(TString(collection_[j] + cut_[i] + "phi"),     TString("Track Phi "+collection_[j]), 100, -TMath::Pi(), TMath::Pi());
            Trk_[collection_[j] + cut_[i] + "phi"]->Sumw2();               
            Trk_[collection_[j] + cut_[i] + "theta"]            = subDir.make<TH1D>(TString(collection_[j] + cut_[i] + "theta"),   TString("Track Theta "+collection_[j]), 100, 0, TMath::Pi());
            Trk_[collection_[j] + cut_[i] + "theta"]->Sumw2();               
            Trk_[collection_[j] + cut_[i] + "Chi2"]             = subDir.make<TH1D>(TString(collection_[j] + cut_[i] + "Chi2"),    TString("Track Chi2 "+collection_[j]), 100, -0.5, 800);
            Trk_[collection_[j] + cut_[i] + "Chi2"]->Sumw2();
            Trk_[collection_[j] + cut_[i] + "Chi2Prob"]         = subDir.make<TH1D>(TString(collection_[j] + cut_[i] + "Chi2Prob"),TString("Chi2 Probability "+collection_[j]), 100, 0, 1);
            Trk_[collection_[j] + cut_[i] + "Chi2Prob"]->Sumw2();
            
            Trk2D_[collection_[j] + cut_[i] + "NumRecHitsPerTrkVsPhi"]   = subDir.make<TH2D>(TString(collection_[j] + cut_[i] + "NumRecHitsPerTrkVsPhi"),   TString("Number Of RecHits Per Track Vs Phi "+collection_[j]),  100, -TMath::Pi(), TMath::Pi(), 35, -0.5, 34.5);
            Trk2D_[collection_[j] + cut_[i] + "NumRecHitsPerTrkVsPhi"]->Sumw2();
            Trk2D_[collection_[j] + cut_[i] + "NumRecHitsPerTrkVsTheta"] = subDir.make<TH2D>(TString(collection_[j] + cut_[i] + "NumRecHitsPerTrkVsTheta"), TString("Number Of RecHits Per Track Vs Theta "+collection_[j]), 100, 0, TMath::Pi(), 35, -0.5, 34.5);
            Trk2D_[collection_[j] + cut_[i] + "NumRecHitsPerTrkVsTheta"]->Sumw2();
            Trk2D_[collection_[j] + cut_[i] + "NumRecHitsPerTrkVsEta"]   = subDir.make<TH2D>(TString(collection_[j] + cut_[i] + "NumRecHitsPerTrkVsEta"),   TString("Number Of RecHits Per Track Vs Eta "+collection_[j]),   52, -2.6, 2.6, 35, -0.5, 34.5);
            Trk2D_[collection_[j] + cut_[i] + "NumRecHitsPerTrkVsEta"]->Sumw2();
            Trk2D_[collection_[j] + cut_[i] + "Chi2overDoFVsTheta"]      = subDir.make<TH2D>(TString(collection_[j] + cut_[i] + "Chi2overDoFVsTheta"), TString("Chi2overDoF Vs Theta "+collection_[j]), 100, 0, TMath::Pi(), 100, 0, 800);
            Trk2D_[collection_[j] + cut_[i] + "Chi2overDoFVsTheta"]->Sumw2();
            Trk2D_[collection_[j] + cut_[i] + "Chi2overDoFVsPhi"]        = subDir.make<TH2D>(TString(collection_[j] + cut_[i] + "Chi2overDoFVsPhi"),   TString("Chi2overDoF Vs Phi "+collection_[j]),  100, -TMath::Pi(), TMath::Pi(), 100, 0, 800);
            Trk2D_[collection_[j] + cut_[i] + "Chi2overDoFVsPhi"]->Sumw2();
            Trk2D_[collection_[j] + cut_[i] + "Chi2overDoFVsEta"]        = subDir.make<TH2D>(TString(collection_[j] + cut_[i] + "Chi2overDoFVsEta"),   TString("Chi2overDoF Vs Eta "+collection_[j]),  100, -3.2, 3.2, 100, -0.5, 800);
            Trk2D_[collection_[j] + cut_[i] + "Chi2overDoFVsEta"]->Sumw2();

            Trk_[collection_[j] + cut_[i] + "Errdxy"]      = subDir.make<TH1D>(TString(collection_[j] + cut_[i] + "Errdxy"),          TString("Track dxy Error "+collection_[j]), 100, 0, 0.45);
            Trk_[collection_[j] + cut_[i] + "Errdxy"]->Sumw2(); 
            Trk_[collection_[j] + cut_[i] + "Errdz"]       = subDir.make<TH1D>(TString(collection_[j] + cut_[i] + "Errdz"),           TString("Track dz Error "+collection_[j]), 100, 0, 25);
            Trk_[collection_[j] + cut_[i] + "Errdz"]->Sumw2();
            Trk_[collection_[j] + cut_[i] + "TrackPtErr"]  = subDir.make<TH1D>(TString(collection_[j] + cut_[i] + "TrackPtErr"),      TString("Track p_{T} Error "+collection_[j]), 100, 0, 3000);
            Trk_[collection_[j] + cut_[i] + "TrackPtErr"]->Sumw2();
            Trk_[collection_[j] + cut_[i] + "TrackPxErr"]  = subDir.make<TH1D>(TString(collection_[j] + cut_[i] + "TrackPxErr"),      TString("Track p_{x} Error "+collection_[j]), 1000, 0, 1000);
            Trk_[collection_[j] + cut_[i] + "TrackPxErr"]->Sumw2();
            Trk_[collection_[j] + cut_[i] + "TrackPyErr"]  = subDir.make<TH1D>(TString(collection_[j] + cut_[i] + "TrackPyErr"),      TString("Track p_{y} Error "+collection_[j]), 1000, 0, 1000);
            Trk_[collection_[j] + cut_[i] + "TrackPyErr"]->Sumw2();
            Trk_[collection_[j] + cut_[i] + "TrackPzErr"]  = subDir.make<TH1D>(TString(collection_[j] + cut_[i] + "TrackPzErr"),      TString("Track p_{z} Error "+collection_[j]), 1000, 0, 1000);
            Trk_[collection_[j] + cut_[i] + "TrackPzErr"]->Sumw2();
            Trk_[collection_[j] + cut_[i] + "TrackPErr"]   = subDir.make<TH1D>(TString(collection_[j] + cut_[i] + "TrackPErr"),       TString("TrackP Error "+collection_[j]),  1000, 0, 1000);
            Trk_[collection_[j] + cut_[i] + "TrackPErr"]->Sumw2();
            Trk_[collection_[j] + cut_[i] + "TrackPhiErr"] = subDir.make<TH1D>(TString(collection_[j] + cut_[i] + "TrackPhiErr"),     TString("TrackPhi Error "+collection_[j]), 100, 0, 0.05);
            Trk_[collection_[j] + cut_[i] + "TrackPhiErr"]->Sumw2();
            Trk_[collection_[j] + cut_[i] + "TrackThetaErr"] = subDir.make<TH1D>(TString(collection_[j] + cut_[i] + "TrackThetaErr"), TString("TrackTheta Error "+collection_[j]), 100, 0, 0.5);
            Trk_[collection_[j] + cut_[i] + "TrackThetaErr"]->Sumw2();
            Trk_[collection_[j] + cut_[i] + "TrackEtaErr"] = subDir.make<TH1D>(TString(collection_[j] + cut_[i] + "TrackEtaErr"),     TString("TrackEta Error "+collection_[j]), 100, 0, 0.5);
            Trk_[collection_[j] + cut_[i] + "TrackEtaErr"]->Sumw2();

        }
    }

    cutsta_[0] = "";
    cutsta_[1] = "looseRechit";    //cut RecHits>20
    cutsta_[2] = "looseXY";        //cut Region X
    cutsta_[3] = "looseY";         //cut Region Y
    cutsta_[4] = "looseZ";         //cut Region Z    
    cutsta_[5] = "looseYZ_x";       //cut Region YZ_x
    cutsta_[6] = "loosetheta";     //cut delta theta <0.25
    cutsta_[7] = "loose_all_no_theta";     //cut delta theta <0.25
    cutsta_[8] = "loose";      
    
    for (int ii = 0; ii < 9; ++ii)
    {
        std::string cutsta = cutsta_[ii]; 
        TFileDirectory dirsta = fs->mkdir("MuCollect");
        TFileDirectory subDirsta = dirsta.mkdir( cutsta );
        MuCollect1I_["MuCollect_Nmuon_" + cutsta_[ii]]          = subDirsta.make<TH1I>(TString("MuCollect_Nmuon_"    +cutsta_[ii]), TString("Number of muon "+cutsta_[ii]), 11, -.5, 10.5);
        MuCollect1I_["MuCollect_Nmuon_" + cutsta_[ii]]->Sumw2();
        MuCollect1I_["MuCollect_Nrechits_" + cutsta_[ii]]       = subDirsta.make<TH1I>(TString("MuCollect_Nrechits_" +cutsta_[ii]), TString("Number of RecHits "+cutsta_[ii]), 51,  -.5,  51.5);
        MuCollect1I_["MuCollect_Nrechits_" + cutsta_[ii]]->Sumw2();
        MuCollect1I_["MuCollect_date_" + cutsta_[ii]]           = subDirsta.make<TH1I>(TString("MuCollect_date_"     +cutsta_[ii]), TString("date " +cutsta_[ii]), nday_*86400, 0.0, nday_*86400);
      //  MuCollect1I_["MuCollect_date_" + cutsta_[ii]]->Sumw2();
        MuCollect1F_["MuCollect_run_" + cutsta_[ii]]            = subDirsta.make<TH1F>(TString("MuCollect_run_"      +cutsta_[ii]), TString("Run " +cutsta_[ii]), maxRun_ - minRun_, minRun_, maxRun_);
        MuCollect1F_["MuCollect_run_" + cutsta_[ii]]->Sumw2();

        MuCollect1D_["MuCollect_d0_" + cutsta_[ii]]           = subDirsta.make<TH1D>(TString("MuCollect_d0_"  +cutsta_[ii]), TString("Impact Parameter (d0) " +cutsta_[ii]), 100, -860, 860);
        MuCollect1D_["MuCollect_d0_" + cutsta_[ii]]->Sumw2();
        MuCollect1D_["MuCollect_dxy_" + cutsta_[ii]]           = subDirsta.make<TH1D>(TString("MuCollect_dxy_" +cutsta_[ii]),TString("Track dxy " +cutsta_[ii]), 100, -860, 860);
        MuCollect1D_["MuCollect_dxy_" + cutsta_[ii]]->Sumw2();               
        MuCollect1D_["MuCollect_dz_"  + cutsta_[ii]]           = subDirsta.make<TH1D>(TString("MuCollect_dz_"  +cutsta_[ii]),TString("Track dz "  +cutsta_[ii]), 100, -1000, 1000);
        MuCollect1D_["MuCollect_dz_"  + cutsta_[ii]]->Sumw2();               
        MuCollect1D_["MuCollect_dsz_" + cutsta_[ii]]           = subDirsta.make<TH1D>(TString("MuCollect_dsz_" +cutsta_[ii]),TString("Track dsz " +cutsta_[ii]), 100, -1300, 1300);
        MuCollect1D_["MuCollect_dsz_" + cutsta_[ii]]->Sumw2();               

        MuCollect1D_["MuCollect_chi2_" + cutsta_[ii]]             = subDirsta.make<TH1D>(TString("MuCollect_chi2_" +cutsta_[ii]), TString("Track Chi2 " +cutsta_[ii]), 100, 0.5, 1100.0);
        MuCollect1D_["MuCollect_chi2_" + cutsta_[ii]]->Sumw2();
        MuCollect1D_["MuCollect_eta_" + cutsta_[ii]]              = subDirsta.make<TH1D>(TString("MuCollect_eta_"  +cutsta_[ii]), TString("Track Eta "  +cutsta_[ii]), 100, -3.0, 3.0);
        MuCollect1D_["MuCollect_eta_" + cutsta_[ii]]->Sumw2();
        MuCollect1D_["MuCollect_theta_" + cutsta_[ii]]            = subDirsta.make<TH1D>(TString("MuCollect_theta_"+cutsta_[ii]), TString("Track Theta "+cutsta_[ii]), 100, 0, 3.2);
        MuCollect1D_["MuCollect_theta_" + cutsta_[ii]]->Sumw2();
        MuCollect1D_["MuCollect_delta_theta_" + cutsta_[ii]]      = subDirsta.make<TH1D>(TString("MuCollect_delta_theta_"+cutsta_[ii]), TString("Delta Theta " +cutsta_[ii]), 100, -3.2, 3.2);
        MuCollect1D_["MuCollect_delta_theta_" + cutsta_[ii]]->Sumw2();
        MuCollect1D_["MuCollect_abs_delta_theta_" + cutsta_[ii]] = subDirsta.make<TH1D>(TString("MuCollect_abs_delta_theta_"+cutsta_[ii]), TString("(Abs) Delta Theta " +cutsta_[ii]), 100, 0, 3.2);
        MuCollect1D_["MuCollect_abs_delta_theta_" + cutsta_[ii]]->Sumw2();
        MuCollect1D_["MuCollect_phi_" + cutsta_[ii]]              = subDirsta.make<TH1D>(TString("MuCollect_phi_"  +cutsta_[ii]), TString("Track Phi " +cutsta_[ii]),  100, -3.2, 3.2);
        MuCollect1D_["MuCollect_phi_" + cutsta_[ii]]->Sumw2(); 
        MuCollect1D_["MuCollect_pt_" + cutsta_[ii]]               = subDirsta.make<TH1D>(TString("MuCollect_pt_"   +cutsta_[ii]), TString("Track pt_{T} (GeV) "+cutsta_[ii]), 100, 0.5, 50.0);
        MuCollect1D_["MuCollect_pt_" + cutsta_[ii]]->Sumw2();
        MuCollect1D_["MuCollect_normalizedchi2_" + cutsta_[ii]]   = subDirsta.make<TH1D>(TString("MuCollect_normalizedchi2_" +cutsta_[ii]), TString("Track normalizedhi2 "+cutsta_[ii]), 100, 0.5, 100.0);
        MuCollect1D_["MuCollect_normalizedchi2_" + cutsta_[ii]]->Sumw2();
        MuCollect1D_["MuCollect_delta_p_" + cutsta_[ii]]          = subDirsta.make<TH1D>(TString("MuCollect_delta_p_" +cutsta_[ii]), TString("Delta p "+cutsta_[ii]), 100, 0.85, 1.05);
        MuCollect1D_["MuCollect_delta_p_" + cutsta_[ii]]->Sumw2();

        MuCollect2D_["MuCollect_PhiVsEta_" + cutsta_[ii]]         = subDirsta.make<TH2D>(TString("MuCollect_PhiVsEta_" +cutsta_[ii]), TString("Phi Vs Eta "+cutsta_[ii]), 100, -2.6, 2.6, 100, -3.2, 3.2);
        MuCollect2D_["MuCollect_PhiVsEta_" + cutsta_[ii]]->Sumw2();

        MuCollect2D_["MuCollect_recHitsGlobPYvsX_" + cutsta_[ii]] = subDirsta.make<TH2D>(TString("MuCollect_recHitsGlobPYvsX_" +cutsta_[ii]),  TString("Globalposition Y vs X of the RecHits "+cutsta_[ii]), 200, -850.0,   850.0, 200, -850.0, 850.0);
        MuCollect2D_["MuCollect_recHitsGlobPXvsZ_" + cutsta_[ii]] = subDirsta.make<TH2D>(TString("MuCollect_recHitsGlobPXvsZ_" +cutsta_[ii]),  TString("Globalposition X vs Z of the RecHits "+cutsta_[ii]), 200, -1200.0, 1200.0, 200, -850.0, 850.0);
        MuCollect2D_["MuCollect_recHitsGlobPYvsZ_" + cutsta_[ii]] = subDirsta.make<TH2D>(TString("MuCollect_recHitsGlobPYvsZ_" +cutsta_[ii]),  TString("Globalposition Y vs Z of the RecHits "+cutsta_[ii]), 200, -1200.0, 1200.0, 200, -850.0, 850.0);

        MuCollect2D_["MuCollect_recHitsGlob_ip_YvsX_" + cutsta_[ii]] = subDirsta.make<TH2D>(TString("MuCollect_recHitsGlob_ip_YvsX_" +cutsta_[ii]),  TString("Globalinnerposition Y vs X of the RecHits "+cutsta_[ii]), 200, -850.0,   850.0, 200, -850.0, 850.0);
        MuCollect2D_["MuCollect_recHitsGlob_ip_XvsZ_" + cutsta_[ii]] = subDirsta.make<TH2D>(TString("MuCollect_recHitsGlob_ip_XvsZ_" +cutsta_[ii]),  TString("Globalinnerposition X vs Z of the RecHits "+cutsta_[ii]), 200, -1200.0, 1200.0, 200, -850.0, 850.0);
        MuCollect2D_["MuCollect_recHitsGlob_ip_YvsZ_" + cutsta_[ii]] = subDirsta.make<TH2D>(TString("MuCollect_recHitsGlob_ip_YvsZ_" +cutsta_[ii]),  TString("Globalinnerposition Y vs Z of the RecHits "+cutsta_[ii]), 200, -1200.0, 1200.0, 200, -850.0, 850.0);
        MuCollect2D_["MuCollect_recHitsGlob_op_YvsX_" + cutsta_[ii]] = subDirsta.make<TH2D>(TString("MuCollect_recHitsGlob_op_YvsX_" +cutsta_[ii]),  TString("Globalouterposition Y vs X of the RecHits "+cutsta_[ii]), 200, -850.0,   850.0, 200, -850.0, 850.0);
        MuCollect2D_["MuCollect_recHitsGlob_op_XvsZ_" + cutsta_[ii]] = subDirsta.make<TH2D>(TString("MuCollect_recHitsGlob_op_XvsZ_" +cutsta_[ii]),  TString("Globalouterposition X vs Z of the RecHits "+cutsta_[ii]), 200, -1200.0, 1200.0, 200, -850.0, 850.0);
        MuCollect2D_["MuCollect_recHitsGlob_op_YvsZ_" + cutsta_[ii]] = subDirsta.make<TH2D>(TString("MuCollect_recHitsGlob_op_YvsZ_" +cutsta_[ii]),  TString("Globalouterposition Y vs Z of the RecHits "+cutsta_[ii]), 200, -1200.0, 1200.0, 200, -850.0, 850.0);

        MuCollect1D_["MuCollect_length2d_"            + cutsta_[ii]]= subDirsta.make<TH1D>(TString("MuCollect_length2d_" +cutsta_[ii]), TString("length x,y "  +cutsta_[ii]), 100, 0, 1200);
        MuCollect1D_["MuCollect_length2d_" + cutsta_[ii]]->Sumw2();
        MuCollect1D_["MuCollect_length3d_"            + cutsta_[ii]]= subDirsta.make<TH1D>(TString("MuCollect_length3d_" +cutsta_[ii]), TString("length x,y,z "+cutsta_[ii]), 100, 0, 1200);
        MuCollect1D_["MuCollect_length3d_" + cutsta_[ii]]->Sumw2();
    }

    TFileDirectory dirsta = fs->mkdir("MuCollect");
    MuCollect_NDThits     = dirsta.make<TH1I>("MuCollect_NDThits", "Number of DT Hits ", 51, -.5, 51.5);
    MuCollect_NCSChits    = dirsta.make<TH1I>("MuCollect_NCSChits", "Number of CSC Hits ", 31, -.5, 31.5);
    MuCollect_NRPChits    = dirsta.make<TH1I>("MuCollect_NRPChits", "Number of RPC Hits ", 11, -.5, 11.5);
    ////MuCollect-----MuCollect
    MuCollect_vx          = dirsta.make<TH1D>("MuCollect_vx", "Track Vertex position x", 200, -650.0, 650.0);
    MuCollect_vx->Sumw2();
    MuCollect_vy          = dirsta.make<TH1D>("MuCollect_vy", "Track Vertex position y", 200, -450.0, 450.0);
    MuCollect_vy->Sumw2();
    MuCollect_vz          = dirsta.make<TH1D>("MuCollect_vz", "Track Vertex position z", 200, -950.0, 950.0);
    MuCollect_vz->Sumw2();
    
    MuCollect_recHitsPositionX = dirsta.make<TH1D>("MuCollect_recHitsPositionX", "X position of the RecHits", 200, -300.0, 300.0);
    MuCollect_recHitsPositionX->Sumw2();
    MuCollect_recHitsPositionY = dirsta.make<TH1D>("MuCollect_recHitsPositionY", "Y position of the RecHits", 200, -300.0, 300.0);
    MuCollect_recHitsPositionY->Sumw2();
    MuCollect_recHitsPositionZ = dirsta.make<TH1D>("MuCollect_recHitsPositionZ", "Z position of the RecHits", 200, -300.0, 300.0);
    MuCollect_recHitsPositionZ->Sumw2();
    
    MuCollect_recHitsGlobalPositionX = dirsta.make<TH1D>("MuCollect_recHitsGlobalPositionX", "X Globalposition of the RecHits", 200, -850.0, 850.0);
    MuCollect_recHitsGlobalPositionX->Sumw2();
    MuCollect_recHitsGlobalPositionY = dirsta.make<TH1D>("MuCollect_recHitsGlobalPositionY", "Y Globalposition of the RecHits", 200, -850.0, 850.0);
    MuCollect_recHitsGlobalPositionY->Sumw2();
    MuCollect_recHitsGlobalPositionZ = dirsta.make<TH1D>("MuCollect_recHitsGlobalPositionZ", "Z Globalposition of the RecHits", 200, -1000.0, 1000.0);
    MuCollect_recHitsGlobalPositionZ->Sumw2();

    MuCollect_innerPositionX = dirsta.make<TH1D>("MuCollect_innerPositionX", "Track innerPositionX", 200, -800.0, 800.0);
    MuCollect_innerPositionX->Sumw2();             
    MuCollect_innerPositionY = dirsta.make<TH1D>("MuCollect_innerPositionY", "Track innerPositionY", 200, -800.0, 800.0);
    MuCollect_innerPositionY->Sumw2();             
    MuCollect_innerPositionZ = dirsta.make<TH1D>("MuCollect_innerPositionZ", "Track innerPositionZ", 200, -800.0, 800.0);
    MuCollect_innerPositionZ->Sumw2();            
    MuCollect_outerPositionX = dirsta.make<TH1D>("MuCollect_outerPositionX", "Track outerPositionX", 200, -800.0, 800.0);
    MuCollect_outerPositionX->Sumw2();            
    MuCollect_outerPositionY = dirsta.make<TH1D>("MuCollect_outerPositionY", "Track outerPositionY", 200, -800.0, 800.0);
    MuCollect_outerPositionY->Sumw2();          
    MuCollect_outerPositionZ = dirsta.make<TH1D>("MuCollect_outerPositionZ", "Track outerPositionZ", 200, -800.0, 800.0);                                            
    MuCollect_outerPositionZ->Sumw2();
    
    MuCollect_innerMomentumX = dirsta.make<TH1D>("MuCollect_innerMomentumX", "Track innerMomentumX", 100, -20.0, 20.0);
    MuCollect_innerMomentumX->Sumw2();             
    MuCollect_innerMomentumY = dirsta.make<TH1D>("MuCollect_innerMomentumY", "Track innerMomentumY", 100, -20.0, 20.0);
    MuCollect_innerMomentumY->Sumw2();             
    MuCollect_innerMomentumZ = dirsta.make<TH1D>("MuCollect_innerMomentumZ", "Track innerMomentumZ", 100, -20.0, 20.0);
    MuCollect_innerMomentumZ->Sumw2();            
    MuCollect_outerMomentumX = dirsta.make<TH1D>("MuCollect_outerMomentumX", "Track outerMomentumX", 100, -20.0, 20.0);
    MuCollect_outerMomentumX->Sumw2();            
    MuCollect_outerMomentumY = dirsta.make<TH1D>("MuCollect_outerMomentumY", "Track outerMomentumY", 100, -20.0, 20.0);
    MuCollect_outerMomentumY->Sumw2();          
    MuCollect_outerMomentumZ = dirsta.make<TH1D>("MuCollect_outerMomentumZ", "Track outerMomentumZ", 100, -20.0, 20.0);                                            
    MuCollect_outerMomentumZ->Sumw2();

//     MuCollect_delta_p = dirsta.make<TH1D>("MuCollect_delta_p", "Delta p", 100, 0.85, 1.05);                                            
//     MuCollect_delta_p->Sumw2();
    
//     MuCollect_recHitsGlobalPositionYvsX = dirsta.make<TH2D>("MuCollect_recHitsGlobalPositionYvsX", "Globalposition Y vs X of the RecHits", 200, -850.0,   850.0, 200, -850.0, 850.0);
//     MuCollect_recHitsGlobalPositionXvsZ = dirsta.make<TH2D>("MuCollect_recHitsGlobalPositionXvsZ", "Globalposition X vs Z of the RecHits", 200, -1200.0, 1200.0, 200, -850.0, 850.0);
//     MuCollect_recHitsGlobalPositionYvsZ = dirsta.make<TH2D>("MuCollect_recHitsGlobalPositionYvsZ", "Globalposition Y vs Z of the RecHits", 200, -1200.0, 1200.0, 200, -850.0, 850.0);
    MuCollect_recHitsGlobalPositionXYZ = dirsta.make<TH3D>("MuCollect_recHitsGlobalPositionXYZ", "Globalposition XYZ of the RecHits", 200, -1200.0, 1200.0, 200, -850.0, 850.0, 200, -850.0, 850.0);

    MuCollect_surfaceXZ = dirsta.make<TH2D>("MuCollect_surfaceXZ", "surface X vs Z (100m)", 300, -24000.0, 20000.0, 300, -22000.0, 22000.0);
    MuCollect_shaftXZ = dirsta.make<TH2D>("MuCollect_shaftXZ", "shaft X vs Z(100m)",  300, -29000.0, 20000.0, 300, -22000.0, 22000.0);
    MuCollect_shaftXZ1 = dirsta.make<TH2D>("MuCollect_shaftXZ1", "shaft X vs Z (100m)",  300, -24000.0, 20000.0, 300, -22000.0, 22000.0);
    
//    eff_muon = dirsta.make<TH2D>("eff_muon", "eff_muon",  600, 50900, 51500, 100, 0, 2);
//    MuCollect_shaftXZ = dirsta.make<TH2D>("MuCollect_shaftXZ", "shaft X vs Z (100m)", 300, -2800.0, 100.0, 300, -600.0, 600.0);
    
    MuCollect_normalizedchi2_shaft = dirsta.make<TH1D>("MuCollect_normalizedchi2_shaft" , "Track normalizedhi2 MuCollectuon shaft" , 100, 0.5, 100.0);
    MuCollect_normalizedchi2_shaft->Sumw2();
    
    MuCollect_nRecHits_shaft = dirsta.make<TH1I>("MuCollect_Nrechits_shaft" , "Number of RecHits MuCollectuon shaft" , 51,  -.5,  51.5);  
    MuCollect_nRecHits_shaft->Sumw2();  
      
     MuCollect_time_date  = fs->make<TH1I>("MuCollect_time_date", "MuCollect_date", nday_*86400, 0.0, nday_*86400);        
     MuCollect_time_date0 = fs->make<TH1I>("MuCollect_time_date0", "MuCollect_date", nday_*86400, 0.0, nday_*86400);        

     MuCollect_time_date1 = fs->make<TH1I>("MuCollect_time_date1", "MuCollect_date", nday_*86400, 0.0, nday_*86400);        
    SiTrk_date = fs->make<TH1I>("SiTrk_date", "SiTrk_date", nday_*86400, 0.0, nday_*86400);        
    SiTrk_date1 = fs->make<TH1I>("SiTrk_date1", "SiTrk_date", nday_*86400, 0.0, nday_*86400);        

    Run_Good_muon = fs->make<TH1F>("Run_Good_muon", "Run_Good_muon",maxRun_ - minRun_, minRun_, maxRun_);
    Run_Good_muon->Sumw2();
    Run_All_muon =  fs->make<TH1F>("Run_All_muon", "Run_All_muon",maxRun_ - minRun_, minRun_, maxRun_);
    Run_All_muon->Sumw2();
    Run_All=  fs->make<TH1F>("Run_All", "Run_All",maxRun_ - minRun_, minRun_, maxRun_);
    Run_All->Sumw2();

}


TrackCruzet::~TrackCruzet()
{
    //if ( file_ != 0 ) {
    //   file_->cd();
    //   file_->Write();
    //   delete file_;
    //}
}
bool TrackCruzet::filter(const edm::Event& iEvent, const edm::EventSetup& iSetup , double radius , double max_z)
{

    int nrun = iEvent.id().run();
    int nevent  = iEvent.id().event();
//    int nSitrk = 0;
//    int nPixel = 0;
//    int nIp = 0;
    std::string ntracker;
    bool accept = false;

    if (!thePropagator)
    {
        ESHandle<Propagator> prop;
        iSetup.get<TrackingComponentsRecord>().get(thePropagatorName, prop);
        thePropagator = prop->clone();
        thePropagator->setPropagationDirection(anyDirection);
    }

    ESHandle<MagneticField> theMGField;
    iSetup.get<IdealMagneticFieldRecord>().get(theMGField);

    ESHandle<GlobalTrackingGeometry> theTrackingGeometry;
    iSetup.get<GlobalTrackingGeometryRecord>().get(theTrackingGeometry);

    // Get the RecTrack collection from the event
    Handle<reco::TrackCollection> staTracks;
    iEvent.getByLabel(theMuCollectionLabel, staTracks);

    if (staTracks->size() == 0)
    {
//   	std::cout<<"    NO theMuCollectuon "<<std::endl;
        return accept;
    }

    reco::TrackCollection::const_iterator staTrack;

    for (staTrack = staTracks->begin(); staTrack != staTracks->end(); ++staTrack)
    {
        reco::TransientTrack track(*staTrack, &*theMGField, theTrackingGeometry);

        TrajectoryStateOnSurface innerTSOS = track.innermostMeasurementState();

        LogDebug("HLTMuonPointing") << " InnerTSOS " << innerTSOS;

        // Get a surface (here a cylinder of radius 1290mm) ECAL
        Cylinder::PositionType pos0;
        Cylinder::RotationType rot0;
        const Cylinder::CylinderPointer cyl = Cylinder::build(pos0, rot0, radius);

        TrajectoryStateOnSurface tsosAtCyl =
            thePropagator->propagate(*innerTSOS.freeState(), *cyl);

        if ( tsosAtCyl.isValid() )
        {
            LogDebug("HLTMuonPointing") << " extrap TSOS " << tsosAtCyl;
            if (fabs(tsosAtCyl.globalPosition().z()) < max_z)
            {

                if (radius == theRadius_SiTrk)
                {
                    ntracker = "SiTrack";
            //        ++nSitrk;
                }
                if (radius == theRadius_Pixel)
                {
                    ntracker = "Pixel  ";
            //        ++nPixel;
                }
                if (radius == theRadius_IP)
                {
                    ntracker = "IP     ";
            //        ++nIp;
                }

             //         std::cout<<ntracker<<"    radio = "<<radius<<" Z= "<<max_z<<"   RUN = "<<nrun<<"   Event = "<<nevent<<std::endl;
                fasciiFile << ntracker << "    radio = " << radius << " Z= " << max_z << "   RUN = " << nrun << "   Event = " << nevent << std::endl;
                accept = true;
                return accept;
            }
            else
            {
                LogDebug("HLTMuonPointing") << " extrap TSOS z too big " << tsosAtCyl.globalPosition().z();
            }
        }
        else
        {
            LogDebug("HLTMuonPointing") << " extrap to cyl failed ";
        }

    }

    return accept;

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 
// ------------ method called to for each event  ------------
void TrackCruzet::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    int flag = -1;
    int flagloose = -1;
    int noflagloose = -1;
    int flagloosemuon1 = -1;

    int run_nr = iEvent.id().run();
    int ev_nr  = iEvent.id().event();
    double time_nr  = iEvent.time().value();
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////   
    
    time_t aTime = iEvent.time().value(); // ctime(&aTime)
    time_t t(iEvent.time().value() >> 32);  // gmtime(&t)
    
    std::string dates(asctime(gmtime(&t) ) );
    std::string dates1(asctime(localtime(&t) ) );

    struct tm* t1 = gmtime(&t);
 //   cout << "hh:mm:ss " << t1->tm_hour << ":" << t1->tm_min << ":"<< t1->tm_sec << endl;
//    cout << "dd-mm-aaaa: " << t1->tm_mday << "-" << t1->tm_mon+1 << "-"<< t1->tm_year+1900 << endl;

//    cout <<" ########################################################"<< endl; 
//    cout <<" fecha utc file                     =   " << dates << endl;
    cout <<" fecha local file                   =   " << dates1 << endl;
//    cout <<" seconds since January 1, 1970 file =   " << time(&aTime) << endl;    
//    cout <<" ########################################################"<< endl;
    
    date_run_event<<"  Run: " << run_nr <<"  Event: " << ev_nr <<"  utc: "<< dates << std::endl;
  
    struct tm timeinfo;
    timeinfo.tm_year = yer_  - 1900;
    timeinfo.tm_mon  = mon_  - 1;
    timeinfo.tm_mday = day_ ;
    timeinfo.tm_hour = hor_ ;
    timeinfo.tm_min  = min_ ;
    timeinfo.tm_sec  = sec_ ;
    timeinfo.tm_isdst = -1; 
    time_t art = mktime(&timeinfo);
    double dif_time = difftime(t, art);
    
    std::string date(asctime(gmtime(&art) ) );
    std::string date2(asctime(localtime(&art) ) );
//    cout <<" ########################################################"<< endl;
//    cout <<" fecha utc in                       =   " << date << endl;
    cout <<" fecha local in                     =   " << date2 << endl;
//    cout <<" seconds since January 1, 1970 in   =   " << time(&art) << endl;
//    cout <<" ########################################################"<< endl;

    cout << "Diferencia local: " << dif_time << " segundos"<< endl;
    
    cout<< "F" << yer_ << "-" << mon_ << "-" << day_ << " " << hor_ << ":" << min_ << ":" << sec_<< endl;
 
    date_data->Fill(dif_time);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////   
    runNumber->Fill(run_nr);

    static const reco::TrackCollection s_empty_ctf;
    static const reco::TrackCollection s_empty_cotf;
    static const reco::TrackCollection s_empty_rs;
    std::map<int, const reco::TrackCollection *>typeTrackCollection;
    typeTrackCollection[0] = &s_empty_ctf;
    typeTrackCollection[1] = &s_empty_cotf;
    typeTrackCollection[2] = &s_empty_rs;
    std::map<int, edm::Handle<reco::TrackCollection> >typeTrackCollectionHandle;
    if ( iEvent.getByLabel(ctfTrackCollectionTag_, typeTrackCollectionHandle[0]))
    {
        typeTrackCollection[0] = typeTrackCollectionHandle[0].product();
    }
    else
    {
        std::cout << "Collection ctfWithMaterialTracks cannot be found -> using empty collection of same type. Run = " << run_nr <<"Event = " << ev_nr << std::endl;
    }
    if ( iEvent.getByLabel(cosmicTFTrackCollectionTag_, typeTrackCollectionHandle[1]))
    {
        typeTrackCollection[1] = typeTrackCollectionHandle[1].product();
    }
    else
    {
        std::cout << "Collection cosmictrackfinder cannot be found -> using empty collection of same type. Run = " << run_nr <<"  Event = " << ev_nr << std::endl;
    }
    if ( iEvent.getByLabel(rsTrackCollectionTag_, typeTrackCollectionHandle[2]))
    {
        typeTrackCollection[2] = typeTrackCollectionHandle[2].product();
    }
    else
    {
        std::cout << "Collection rsWithMaterialTracks cannot be found -> using empty collection of same type. Run = " << run_nr <<"  Event = " << ev_nr << std::endl;
    }
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
    if (typeTrackCollection[0]->size() == 0 && typeTrackCollection[1]->size() == 0 && typeTrackCollection[2]->size() == 0) flag = 0;
    bool pass = filter(iEvent, iSetup, theRadius_SiTrk, theMaxZ_SiTrk);
    if (pass){
        flag = 1;
        SiTrk_date1->Fill(dif_time);  
        SiTrk_date->Fill(dif_time);  
        SiTrk_date->GetXaxis()->SetTimeDisplay(1);
        SiTrk_date->GetXaxis()->SetLabelSize(0.03);
   
        TDatime T0(yer_,mon_,day_,hor_+12,min_,sec_);
        SiTrk_date->GetXaxis()->SetTimeOffset(T0.Convert(),"local");
        SiTrk_date->GetXaxis()->SetTimeFormat("#splitline{%H\:%M\:%S}{%d\/%m\/%y}");
    
    }
    pass = filter(iEvent, iSetup, theRadius_Pixel, theMaxZ_Pixel);
    if (pass) flag = 2;
    pass = filter(iEvent, iSetup, theRadius_IP, theMaxZ_IP);
    if (pass) flag = 3;
    Flag->Fill(flag);

   // Get the RecTrack collection from the event
    ESHandle<MagneticField> theMGField;
    iSetup.get<IdealMagneticFieldRecord>().get(theMGField);

    ESHandle<GlobalTrackingGeometry> theTrackingGeometry;
    iSetup.get<GlobalTrackingGeometryRecord>().get(theTrackingGeometry);
 
    edm::Handle<reco::TrackCollection> staTracks;
    iEvent.getByLabel(theMuCollectionLabel, staTracks);

    reco::TrackCollection::const_iterator staTrack;

    int Nmuonloose1 = 0;
    int Nmuonloose = 0;
    int NmuonlooseY = 0;
    int NmuonlooseZ = 0;
    int NmuonlooseXY = 0;
    int NmuonlooseYZ_x = 0;
    int NmuonlooseRechit = 0;
    int NmuonlooseTheta = 0;
    int Nmuonloosenotheta = 0;
    double recHitsGlobX = 0; 
    double recHitsGlobY = 0;
    double recHitsGlobZ = 0;
    
    int NrecHitYP = 0;
    int NrecHitYN = 0;
    int NrecHitZ = 0;
    int NrecHitXYP = 0;
    int NrecHitXYN = 0;

    double px_in = 0;
    double py_in = 0;
    double pz_in = 0;
    double px_out = 0;
    double py_out = 0;
    double pz_out = 0;
    double Zox = 0;
    double Zoy = 0;
    double Y0 = 0;
    double X0 = 0;
    double X1 = 0;
    double X2 = 0;
    double Z1 = 0;
    double Z2 = 0;
    double Xoy = 0;
    Run_All->Fill(run_nr);
//   edm::LogInfo ("Muons") << "Reconstructed tracks: " << staTracks->size() << endl;
//  std::cout << "MuCollectuons  size  =  " << staTracks->size() << endl;
    MuCollect1I_["MuCollect_Nmuon_" + cutsta_[0]]->Fill(staTracks->size());
    //loop on the Muons//////////////////////////////////////////////////////////////////////////////////////////////////////
    for (staTrack = staTracks->begin(); staTrack != staTracks->end(); ++staTrack)
    {
        reco::TransientTrack track(*staTrack, &*theMGField, theTrackingGeometry);

        int nrechits = 0;
        int nDThits  = 0;
        int nCSChits = 0;
        int nRPChits = 0;
        
        px_in  = staTrack->innerMomentum().x();
        py_in  = staTrack->innerMomentum().y();
        pz_in  = staTrack->innerMomentum().z();
        px_out = staTrack->outerMomentum().x();
        py_out = staTrack->outerMomentum().y();
        pz_out = staTrack->outerMomentum().z();
        
        math::XYZPoint innerPo = staTrack->innerPosition();
        GlobalPoint ip(innerPo.x(), innerPo.y(),innerPo.z());

        math::XYZPoint outerPo = staTrack->outerPosition();
        GlobalPoint op(outerPo.x(), outerPo.y(),outerPo.z());
        
        double z_2d = sqrt((ip.x() - op.x())*(ip.x() - op.x()) + (ip.y() - op.y())*(ip.y() - op.y()));
        double z_3d = sqrt((ip.x() - op.x())*(ip.x() - op.x()) + (ip.y() - op.y())*(ip.y() - op.y())+ (ip.z() - op.z())*(ip.z() - op.z()));
 
        double norm_px_in   = sqrt( px_in*px_in + py_in*py_in + pz_in*pz_in );
        double norm_px_out  = sqrt( px_out*px_out + py_out*py_out + pz_out*pz_out);
        double px_in_px_out = px_in*px_out + py_in*py_out + pz_in*pz_out;
                
        double delta_p = px_in_px_out/(norm_px_in*norm_px_out);
        
        double deltaTheta = staTrack->outerMomentum().Theta() - staTrack->innerMomentum().Theta();

        Run_All_muon->Fill(run_nr);
        
        // loop on the RecHits
        for (trackingRecHit_iterator it = track.recHitsBegin ();  it != track.recHitsEnd (); it++)
        {
            if ((*it)->isValid ())
            {
                edm::LogInfo ("Muons") << "Analyzer:  Aha this looks like a Rechit!" << std::endl;
                if ((*it)->geographicalId().subdetId() == MuonSubdetId::DT)
                {
                    nDThits++;
                }
                else if ((*it)->geographicalId().subdetId() == MuonSubdetId::CSC)
                {
                    nCSChits++;
                }
                else if ((*it)->geographicalId().subdetId() == MuonSubdetId::RPC)
                {
                    nRPChits++;
                }
                else
                {
                    edm::LogInfo ("Muons") <<  "This is an UNKNOWN hit !! " << std::endl;
                }
                
                nrechits++;

                const GeomDet* geomDet = theTrackingGeometry->idToDet((*it)->geographicalId());
                recHitsGlobX = geomDet->toGlobal((*it)->localPosition()).x();
                recHitsGlobY = geomDet->toGlobal((*it)->localPosition()).y();
                recHitsGlobZ = geomDet->toGlobal((*it)->localPosition()).z();
//                cout<<" RecHitsGlobalPosition " <<" x: "<<recHitsGlobX <<" y: "<<recHitsGlobY <<" z: "<<recHitsGlobZ <<endl;

                double X = (*it)->localPosition().x();
                double Y = (*it)->localPosition().y();
                double Z = (*it)->localPosition().z();
//                cout<<" RecHitsLocalPosition  "<<" X: "<<X <<" Y: "<<Y <<" Z: "<<Z <<endl;
                
                if (recHitsGlobY > posY_) 
                {
                    NrecHitYP++;
//                    cout<<" RecHitsGlobY  "<<recHitsGlobY<<" RecHitsY+  "<<NrecHitYP <<endl;
                }
                if (recHitsGlobY < -posY_)
                {
                    NrecHitYN++;
//                    cout<<" RecHitsGlobY  "<<recHitsGlobY<<" RecHitsY-  "<<NrecHitYN <<endl;
                }
                if (fabs(recHitsGlobZ) < posZ_)
                {
                    NrecHitZ++;
//                    cout<<" RecHitsGlobZ  "<<recHitsGlobZ<<" RecHitsZ  "<<NrecHitZ <<endl;
                }
                if (recHitsGlobX > posX_ && fabs(recHitsGlobY) < posY_)
                {
                    NrecHitXYP++;
                }
                if (recHitsGlobX < -posX_ && fabs(recHitsGlobY) < posY_)
                {
                    NrecHitXYN++;
                }
                MuCollect_recHitsPositionX->Fill((*it)->localPosition().x());
                MuCollect_recHitsPositionY->Fill((*it)->localPosition().y());
                MuCollect_recHitsPositionZ->Fill((*it)->localPosition().z());

                MuCollect_recHitsGlobalPositionX->Fill(recHitsGlobX);
                MuCollect_recHitsGlobalPositionY->Fill(recHitsGlobY);
                MuCollect_recHitsGlobalPositionZ->Fill(recHitsGlobZ);

//                 MuCollect_recHitsGlobalPositionYvsX->Fill(recHitsGlobX,recHitsGlobY);
//                 MuCollect_recHitsGlobalPositionXvsZ->Fill(recHitsGlobZ,recHitsGlobX);
//                 MuCollect_recHitsGlobalPositionYvsZ->Fill(recHitsGlobZ,recHitsGlobY);
                MuCollect_recHitsGlobalPositionXYZ->Fill(recHitsGlobZ,recHitsGlobX,recHitsGlobY);
                
                if((recHitsGlobX >= -200) && (recHitsGlobX <= 0) && (recHitsGlobY >= -200) && (recHitsGlobY <= 0))
                {
                    filefussy<<" Fussy  "<< run_nr <<"  Event = " << ev_nr<<"  date GMT = " <<dates << std::endl;
                }

            }                

        }//end loop on the RecHits//////////////////////////////////////////////////////////////////////////////////////////

        if (nrechits > 20)
        {
            NmuonlooseRechit++;
            MuCollect1I_["MuCollect_Nrechits_"         + cutsta_[1]]->Fill(nrechits);      
            MuCollect1I_["MuCollect_date_"             + cutsta_[1]]->Fill(dif_time);
            MuCollect1F_["MuCollect_run_"              + cutsta_[1]]->Fill(run_nr);
            MuCollect1D_["MuCollect_d0_"               + cutsta_[1]]->Fill(staTrack->d0());
            MuCollect1D_["MuCollect_dxy_"              + cutsta_[1]]->Fill(staTrack->dxy());         
            MuCollect1D_["MuCollect_dz_"               + cutsta_[1]]->Fill(staTrack->dz());         
            MuCollect1D_["MuCollect_dsz_"              + cutsta_[1]]->Fill(staTrack->dsz());              
            MuCollect1D_["MuCollect_chi2_"             + cutsta_[1]]->Fill(staTrack->chi2());         
            MuCollect1D_["MuCollect_normalizedchi2_"   + cutsta_[1]]->Fill(staTrack->normalizedChi2());
            MuCollect1D_["MuCollect_eta_"              + cutsta_[1]]->Fill(staTrack->eta());           
            MuCollect1D_["MuCollect_theta_"            + cutsta_[1]]->Fill(staTrack->theta());           
            MuCollect1D_["MuCollect_delta_theta_"      + cutsta_[1]]->Fill(deltaTheta);
            MuCollect1D_["MuCollect_abs_delta_theta_"  + cutsta_[1]]->Fill(abs(deltaTheta));
            MuCollect1D_["MuCollect_phi_"              + cutsta_[1]]->Fill(staTrack->phi());         
            MuCollect1D_["MuCollect_pt_"               + cutsta_[1]]->Fill(staTrack->pt());          
            MuCollect1D_["MuCollect_delta_p_"          + cutsta_[1]]->Fill(delta_p);
            MuCollect2D_["MuCollect_PhiVsEta_"         + cutsta_[1]]->Fill(staTrack->eta(), staTrack->phi());      
            for (trackingRecHit_iterator it = track.recHitsBegin ();  it != track.recHitsEnd (); it++)
            {
                if ((*it)->isValid ())
                {
                    const GeomDet* geomDet = theTrackingGeometry->idToDet((*it)->geographicalId());
                    recHitsGlobX = geomDet->toGlobal((*it)->localPosition()).x();
                    recHitsGlobY = geomDet->toGlobal((*it)->localPosition()).y();
                    recHitsGlobZ = geomDet->toGlobal((*it)->localPosition()).z();
                    MuCollect2D_["MuCollect_recHitsGlobPYvsX_" + cutsta_[1]]->Fill(recHitsGlobX,recHitsGlobY);
                    MuCollect2D_["MuCollect_recHitsGlobPXvsZ_" + cutsta_[1]]->Fill(recHitsGlobZ,recHitsGlobX);
                    MuCollect2D_["MuCollect_recHitsGlobPYvsZ_" + cutsta_[1]]->Fill(recHitsGlobZ,recHitsGlobY);
                }
            }
            MuCollect2D_["MuCollect_recHitsGlob_ip_YvsX_" + cutsta_[1]]->Fill(ip.x(),ip.y());
            MuCollect2D_["MuCollect_recHitsGlob_ip_XvsZ_" + cutsta_[1]]->Fill(ip.z(),ip.x());
            MuCollect2D_["MuCollect_recHitsGlob_ip_YvsZ_" + cutsta_[1]]->Fill(ip.z(),ip.y());
            MuCollect2D_["MuCollect_recHitsGlob_op_YvsX_" + cutsta_[1]]->Fill(op.x(),op.y());
            MuCollect2D_["MuCollect_recHitsGlob_op_XvsZ_" + cutsta_[1]]->Fill(op.z(),op.x());
            MuCollect2D_["MuCollect_recHitsGlob_op_YvsZ_" + cutsta_[1]]->Fill(op.z(),op.y());

            MuCollect1D_["MuCollect_length2d_"            + cutsta_[1]]->Fill(z_2d);
            MuCollect1D_["MuCollect_length3d_"            + cutsta_[1]]->Fill(z_3d);
         }
                
        if (NrecHitXYP == 0 && NrecHitXYN == 0)
        {
	    NmuonlooseXY++;
            MuCollect1I_["MuCollect_Nrechits_"         + cutsta_[2]]->Fill(nrechits);      
            MuCollect1I_["MuCollect_date_"             + cutsta_[2]]->Fill(dif_time);
            MuCollect1F_["MuCollect_run_"              + cutsta_[2]]->Fill(run_nr);
            MuCollect1D_["MuCollect_d0_"               + cutsta_[2]]->Fill(staTrack->d0());
            MuCollect1D_["MuCollect_dxy_"              + cutsta_[2]]->Fill(staTrack->dxy());         
            MuCollect1D_["MuCollect_dz_"               + cutsta_[2]]->Fill(staTrack->dz());         
            MuCollect1D_["MuCollect_dsz_"              + cutsta_[2]]->Fill(staTrack->dsz());         
            MuCollect1D_["MuCollect_chi2_"             + cutsta_[2]]->Fill(staTrack->chi2());         
            MuCollect1D_["MuCollect_normalizedchi2_"   + cutsta_[2]]->Fill(staTrack->normalizedChi2());
            MuCollect1D_["MuCollect_eta_"              + cutsta_[2]]->Fill(staTrack->eta());           
            MuCollect1D_["MuCollect_theta_"            + cutsta_[2]]->Fill(staTrack->theta());           
            MuCollect1D_["MuCollect_delta_theta_"      + cutsta_[2]]->Fill(deltaTheta);
            MuCollect1D_["MuCollect_abs_delta_theta_"  + cutsta_[2]]->Fill(abs(deltaTheta));
            MuCollect1D_["MuCollect_phi_"              + cutsta_[2]]->Fill(staTrack->phi());         
            MuCollect1D_["MuCollect_pt_"               + cutsta_[2]]->Fill(staTrack->pt());          
            MuCollect1D_["MuCollect_delta_p_"          + cutsta_[2]]->Fill(delta_p);
            MuCollect2D_["MuCollect_PhiVsEta_"         + cutsta_[2]]->Fill(staTrack->eta(), staTrack->phi());      
            for (trackingRecHit_iterator it = track.recHitsBegin ();  it != track.recHitsEnd (); it++)
            {
                if ((*it)->isValid ())
                {
                    const GeomDet* geomDet = theTrackingGeometry->idToDet((*it)->geographicalId());
                    recHitsGlobX = geomDet->toGlobal((*it)->localPosition()).x();
                    recHitsGlobY = geomDet->toGlobal((*it)->localPosition()).y();
                    recHitsGlobZ = geomDet->toGlobal((*it)->localPosition()).z();
            MuCollect2D_["MuCollect_recHitsGlobPYvsX_" + cutsta_[2]]->Fill(recHitsGlobX,recHitsGlobY);
            MuCollect2D_["MuCollect_recHitsGlobPXvsZ_" + cutsta_[2]]->Fill(recHitsGlobZ,recHitsGlobX);
            MuCollect2D_["MuCollect_recHitsGlobPYvsZ_" + cutsta_[2]]->Fill(recHitsGlobZ,recHitsGlobY);
                }
            }            
            MuCollect2D_["MuCollect_recHitsGlob_ip_YvsX_" + cutsta_[2]]->Fill(ip.x(),ip.y());
            MuCollect2D_["MuCollect_recHitsGlob_ip_XvsZ_" + cutsta_[2]]->Fill(ip.z(),ip.x());
            MuCollect2D_["MuCollect_recHitsGlob_ip_YvsZ_" + cutsta_[2]]->Fill(ip.z(),ip.y());
            MuCollect2D_["MuCollect_recHitsGlob_op_YvsX_" + cutsta_[2]]->Fill(op.x(),op.y());
            MuCollect2D_["MuCollect_recHitsGlob_op_XvsZ_" + cutsta_[2]]->Fill(op.z(),op.x());
            MuCollect2D_["MuCollect_recHitsGlob_op_YvsZ_" + cutsta_[2]]->Fill(op.z(),op.y());

            MuCollect1D_["MuCollect_length2d_"            + cutsta_[2]]->Fill(z_2d);
            MuCollect1D_["MuCollect_length3d_"            + cutsta_[2]]->Fill(z_3d);
        } 
               
        if (NrecHitYP > 0 && NrecHitYN > 0)
        {
	    NmuonlooseY++;
            MuCollect1I_["MuCollect_Nrechits_"         + cutsta_[3]]->Fill(nrechits);      
            MuCollect1I_["MuCollect_date_"             + cutsta_[3]]->Fill(dif_time);
            MuCollect1F_["MuCollect_run_"              + cutsta_[3]]->Fill(run_nr);
            MuCollect1D_["MuCollect_d0_"               + cutsta_[3]]->Fill(staTrack->d0());
            MuCollect1D_["MuCollect_dxy_"              + cutsta_[3]]->Fill(staTrack->dxy());         
            MuCollect1D_["MuCollect_dz_"               + cutsta_[3]]->Fill(staTrack->dz());         
            MuCollect1D_["MuCollect_dsz_"              + cutsta_[3]]->Fill(staTrack->dsz());         
            MuCollect1D_["MuCollect_chi2_"             + cutsta_[3]]->Fill(staTrack->chi2());         
            MuCollect1D_["MuCollect_normalizedchi2_"   + cutsta_[3]]->Fill(staTrack->normalizedChi2());
            MuCollect1D_["MuCollect_eta_"              + cutsta_[3]]->Fill(staTrack->eta());           
            MuCollect1D_["MuCollect_theta_"            + cutsta_[3]]->Fill(staTrack->theta());           
            MuCollect1D_["MuCollect_delta_theta_"      + cutsta_[3]]->Fill(deltaTheta);
            MuCollect1D_["MuCollect_abs_delta_theta_"  + cutsta_[3]]->Fill(abs(deltaTheta));
            MuCollect1D_["MuCollect_phi_"              + cutsta_[3]]->Fill(staTrack->phi());         
            MuCollect1D_["MuCollect_pt_"               + cutsta_[3]]->Fill(staTrack->pt());          
            MuCollect1D_["MuCollect_delta_p_"          + cutsta_[3]]->Fill(delta_p);
            MuCollect2D_["MuCollect_PhiVsEta_"         + cutsta_[3]]->Fill(staTrack->eta(), staTrack->phi());      
            for (trackingRecHit_iterator it = track.recHitsBegin ();  it != track.recHitsEnd (); it++)
            {
                if ((*it)->isValid ())
                {
                    const GeomDet* geomDet = theTrackingGeometry->idToDet((*it)->geographicalId());
                    recHitsGlobX = geomDet->toGlobal((*it)->localPosition()).x();
                    recHitsGlobY = geomDet->toGlobal((*it)->localPosition()).y();
                    recHitsGlobZ = geomDet->toGlobal((*it)->localPosition()).z();
                    MuCollect2D_["MuCollect_recHitsGlobPYvsX_" + cutsta_[3]]->Fill(recHitsGlobX,recHitsGlobY);
                    MuCollect2D_["MuCollect_recHitsGlobPXvsZ_" + cutsta_[3]]->Fill(recHitsGlobZ,recHitsGlobX);
                    MuCollect2D_["MuCollect_recHitsGlobPYvsZ_" + cutsta_[3]]->Fill(recHitsGlobZ,recHitsGlobY);
                }
            }            
            MuCollect2D_["MuCollect_recHitsGlob_ip_YvsX_" + cutsta_[3]]->Fill(ip.x(),ip.y());
            MuCollect2D_["MuCollect_recHitsGlob_ip_XvsZ_" + cutsta_[3]]->Fill(ip.z(),ip.x());
            MuCollect2D_["MuCollect_recHitsGlob_ip_YvsZ_" + cutsta_[3]]->Fill(ip.z(),ip.y());
            MuCollect2D_["MuCollect_recHitsGlob_op_YvsX_" + cutsta_[3]]->Fill(op.x(),op.y());
            MuCollect2D_["MuCollect_recHitsGlob_op_XvsZ_" + cutsta_[3]]->Fill(op.z(),op.x());
            MuCollect2D_["MuCollect_recHitsGlob_op_YvsZ_" + cutsta_[3]]->Fill(op.z(),op.y());

            MuCollect1D_["MuCollect_length2d_"            + cutsta_[3]]->Fill(z_2d);
            MuCollect1D_["MuCollect_length3d_"            + cutsta_[3]]->Fill(z_3d);
        }
        
        if (NrecHitZ > 0)
        {
	    NmuonlooseZ++;
            MuCollect1I_["MuCollect_Nrechits_"         + cutsta_[4]]->Fill(nrechits);      
            MuCollect1I_["MuCollect_date_"             + cutsta_[4]]->Fill(dif_time);
            MuCollect1F_["MuCollect_run_"              + cutsta_[4]]->Fill(run_nr);
            MuCollect1D_["MuCollect_d0_"               + cutsta_[4]]->Fill(staTrack->d0());
            MuCollect1D_["MuCollect_dxy_"              + cutsta_[4]]->Fill(staTrack->dxy());         
            MuCollect1D_["MuCollect_dz_"               + cutsta_[4]]->Fill(staTrack->dz());         
            MuCollect1D_["MuCollect_dsz_"              + cutsta_[4]]->Fill(staTrack->dsz());         
            MuCollect1D_["MuCollect_chi2_"             + cutsta_[4]]->Fill(staTrack->chi2());         
            MuCollect1D_["MuCollect_normalizedchi2_"   + cutsta_[4]]->Fill(staTrack->normalizedChi2());
            MuCollect1D_["MuCollect_eta_"              + cutsta_[4]]->Fill(staTrack->eta());           
            MuCollect1D_["MuCollect_theta_"            + cutsta_[4]]->Fill(staTrack->theta());           
            MuCollect1D_["MuCollect_delta_theta_"      + cutsta_[4]]->Fill(deltaTheta);
            MuCollect1D_["MuCollect_abs_delta_theta_"  + cutsta_[4]]->Fill(abs(deltaTheta));
            MuCollect1D_["MuCollect_phi_"              + cutsta_[4]]->Fill(staTrack->phi());         
            MuCollect1D_["MuCollect_pt_"               + cutsta_[4]]->Fill(staTrack->pt());          
            MuCollect1D_["MuCollect_delta_p_"          + cutsta_[4]]->Fill(delta_p);
            MuCollect2D_["MuCollect_PhiVsEta_"         + cutsta_[4]]->Fill(staTrack->eta(), staTrack->phi());      
            for (trackingRecHit_iterator it = track.recHitsBegin ();  it != track.recHitsEnd (); it++)
            {
                if ((*it)->isValid ())
                {
                    const GeomDet* geomDet = theTrackingGeometry->idToDet((*it)->geographicalId());
                    recHitsGlobX = geomDet->toGlobal((*it)->localPosition()).x();
                    recHitsGlobY = geomDet->toGlobal((*it)->localPosition()).y();
                    recHitsGlobZ = geomDet->toGlobal((*it)->localPosition()).z();
                    MuCollect2D_["MuCollect_recHitsGlobPYvsX_" + cutsta_[4]]->Fill(recHitsGlobX,recHitsGlobY);
                    MuCollect2D_["MuCollect_recHitsGlobPXvsZ_" + cutsta_[4]]->Fill(recHitsGlobZ,recHitsGlobX);
                    MuCollect2D_["MuCollect_recHitsGlobPYvsZ_" + cutsta_[4]]->Fill(recHitsGlobZ,recHitsGlobY);
                }
            }            
            MuCollect2D_["MuCollect_recHitsGlob_ip_YvsX_" + cutsta_[4]]->Fill(ip.x(),ip.y());
            MuCollect2D_["MuCollect_recHitsGlob_ip_XvsZ_" + cutsta_[4]]->Fill(ip.z(),ip.x());
            MuCollect2D_["MuCollect_recHitsGlob_ip_YvsZ_" + cutsta_[4]]->Fill(ip.z(),ip.y());
            MuCollect2D_["MuCollect_recHitsGlob_op_YvsX_" + cutsta_[4]]->Fill(op.x(),op.y());
            MuCollect2D_["MuCollect_recHitsGlob_op_XvsZ_" + cutsta_[4]]->Fill(op.z(),op.x());
            MuCollect2D_["MuCollect_recHitsGlob_op_YvsZ_" + cutsta_[4]]->Fill(op.z(),op.y());

            MuCollect1D_["MuCollect_length2d_"            + cutsta_[4]]->Fill(z_2d);
            MuCollect1D_["MuCollect_length3d_"            + cutsta_[4]]->Fill(z_3d);
        }
        
        if ((NrecHitXYP == 0 && NrecHitXYN == 0)&& (NrecHitYP > 0 && NrecHitYN > 0) && NrecHitZ > 0)
        {
	    NmuonlooseYZ_x++;
            MuCollect1I_["MuCollect_Nrechits_"         + cutsta_[5]]->Fill(nrechits);      
            MuCollect1I_["MuCollect_date_"             + cutsta_[5]]->Fill(dif_time);
            MuCollect1F_["MuCollect_run_"              + cutsta_[5]]->Fill(run_nr);
            MuCollect1D_["MuCollect_d0_"               + cutsta_[5]]->Fill(staTrack->d0());
            MuCollect1D_["MuCollect_dxy_"              + cutsta_[5]]->Fill(staTrack->dxy());         
            MuCollect1D_["MuCollect_dz_"               + cutsta_[5]]->Fill(staTrack->dz());         
            MuCollect1D_["MuCollect_dsz_"              + cutsta_[5]]->Fill(staTrack->dsz());         
            MuCollect1D_["MuCollect_chi2_"             + cutsta_[5]]->Fill(staTrack->chi2());         
            MuCollect1D_["MuCollect_normalizedchi2_"   + cutsta_[5]]->Fill(staTrack->normalizedChi2());
            MuCollect1D_["MuCollect_eta_"              + cutsta_[5]]->Fill(staTrack->eta());           
            MuCollect1D_["MuCollect_theta_"            + cutsta_[5]]->Fill(staTrack->theta());           
            MuCollect1D_["MuCollect_delta_theta_"      + cutsta_[5]]->Fill(deltaTheta);
            MuCollect1D_["MuCollect_abs_delta_theta_"  + cutsta_[5]]->Fill(abs(deltaTheta));
            MuCollect1D_["MuCollect_phi_"              + cutsta_[5]]->Fill(staTrack->phi());         
            MuCollect1D_["MuCollect_pt_"               + cutsta_[5]]->Fill(staTrack->pt());          
            MuCollect1D_["MuCollect_delta_p_"          + cutsta_[5]]->Fill(delta_p);
            MuCollect2D_["MuCollect_PhiVsEta_"         + cutsta_[5]]->Fill(staTrack->eta(), staTrack->phi());      
            for (trackingRecHit_iterator it = track.recHitsBegin ();  it != track.recHitsEnd (); it++)
            {
                if ((*it)->isValid ())
                {
                    const GeomDet* geomDet = theTrackingGeometry->idToDet((*it)->geographicalId());
                    recHitsGlobX = geomDet->toGlobal((*it)->localPosition()).x();
                    recHitsGlobY = geomDet->toGlobal((*it)->localPosition()).y();
                    recHitsGlobZ = geomDet->toGlobal((*it)->localPosition()).z();
                    MuCollect2D_["MuCollect_recHitsGlobPYvsX_" + cutsta_[5]]->Fill(recHitsGlobX,recHitsGlobY);
                    MuCollect2D_["MuCollect_recHitsGlobPXvsZ_" + cutsta_[5]]->Fill(recHitsGlobZ,recHitsGlobX);
                    MuCollect2D_["MuCollect_recHitsGlobPYvsZ_" + cutsta_[5]]->Fill(recHitsGlobZ,recHitsGlobY);
                }
            }
            
            MuCollect2D_["MuCollect_recHitsGlob_ip_YvsX_" + cutsta_[5]]->Fill(ip.x(),ip.y());
            MuCollect2D_["MuCollect_recHitsGlob_ip_XvsZ_" + cutsta_[5]]->Fill(ip.z(),ip.x());
            MuCollect2D_["MuCollect_recHitsGlob_ip_YvsZ_" + cutsta_[5]]->Fill(ip.z(),ip.y());
            MuCollect2D_["MuCollect_recHitsGlob_op_YvsX_" + cutsta_[5]]->Fill(op.x(),op.y());
            MuCollect2D_["MuCollect_recHitsGlob_op_XvsZ_" + cutsta_[5]]->Fill(op.z(),op.x());
            MuCollect2D_["MuCollect_recHitsGlob_op_YvsZ_" + cutsta_[5]]->Fill(op.z(),op.y());

            MuCollect1D_["MuCollect_length2d_"            + cutsta_[5]]->Fill(z_2d);
            MuCollect1D_["MuCollect_length3d_"            + cutsta_[5]]->Fill(z_3d);
        }
        
        if (abs(deltaTheta) < 0.25)
        {
            NmuonlooseTheta++;
            MuCollect1I_["MuCollect_Nrechits_"         + cutsta_[6]]->Fill(nrechits);      
            MuCollect1I_["MuCollect_date_"             + cutsta_[6]]->Fill(dif_time);
            MuCollect1F_["MuCollect_run_"              + cutsta_[6]]->Fill(run_nr);
            MuCollect1D_["MuCollect_d0_"               + cutsta_[6]]->Fill(staTrack->d0());
            MuCollect1D_["MuCollect_dxy_"              + cutsta_[6]]->Fill(staTrack->dxy());         
            MuCollect1D_["MuCollect_dz_"               + cutsta_[6]]->Fill(staTrack->dz());         
            MuCollect1D_["MuCollect_dsz_"              + cutsta_[6]]->Fill(staTrack->dsz());         
            MuCollect1D_["MuCollect_chi2_"             + cutsta_[6]]->Fill(staTrack->chi2());         
            MuCollect1D_["MuCollect_normalizedchi2_"   + cutsta_[6]]->Fill(staTrack->normalizedChi2());
            MuCollect1D_["MuCollect_eta_"              + cutsta_[6]]->Fill(staTrack->eta());           
            MuCollect1D_["MuCollect_theta_"            + cutsta_[6]]->Fill(staTrack->theta());           
            MuCollect1D_["MuCollect_delta_theta_"      + cutsta_[6]]->Fill(deltaTheta);
            MuCollect1D_["MuCollect_abs_delta_theta_"  + cutsta_[6]]->Fill(abs(deltaTheta));
            MuCollect1D_["MuCollect_phi_"              + cutsta_[6]]->Fill(staTrack->phi());         
            MuCollect1D_["MuCollect_pt_"               + cutsta_[6]]->Fill(staTrack->pt());          
            MuCollect1D_["MuCollect_delta_p_"          + cutsta_[6]]->Fill(delta_p);
            MuCollect2D_["MuCollect_PhiVsEta_"         + cutsta_[6]]->Fill(staTrack->eta(), staTrack->phi());      
            for (trackingRecHit_iterator it = track.recHitsBegin ();  it != track.recHitsEnd (); it++)
            {
                if ((*it)->isValid ())
                {
                    const GeomDet* geomDet = theTrackingGeometry->idToDet((*it)->geographicalId());
                    recHitsGlobX = geomDet->toGlobal((*it)->localPosition()).x();
                    recHitsGlobY = geomDet->toGlobal((*it)->localPosition()).y();
                    recHitsGlobZ = geomDet->toGlobal((*it)->localPosition()).z();
                    MuCollect2D_["MuCollect_recHitsGlobPYvsX_" + cutsta_[6]]->Fill(recHitsGlobX,recHitsGlobY);
                    MuCollect2D_["MuCollect_recHitsGlobPXvsZ_" + cutsta_[6]]->Fill(recHitsGlobZ,recHitsGlobX);
                    MuCollect2D_["MuCollect_recHitsGlobPYvsZ_" + cutsta_[6]]->Fill(recHitsGlobZ,recHitsGlobY);
                }
            }
            
            MuCollect2D_["MuCollect_recHitsGlob_ip_YvsX_" + cutsta_[6]]->Fill(ip.x(),ip.y());
            MuCollect2D_["MuCollect_recHitsGlob_ip_XvsZ_" + cutsta_[6]]->Fill(ip.z(),ip.x());
            MuCollect2D_["MuCollect_recHitsGlob_ip_YvsZ_" + cutsta_[6]]->Fill(ip.z(),ip.y());
            MuCollect2D_["MuCollect_recHitsGlob_op_YvsX_" + cutsta_[6]]->Fill(op.x(),op.y());
            MuCollect2D_["MuCollect_recHitsGlob_op_XvsZ_" + cutsta_[6]]->Fill(op.z(),op.x());
            MuCollect2D_["MuCollect_recHitsGlob_op_YvsZ_" + cutsta_[6]]->Fill(op.z(),op.y());

            MuCollect1D_["MuCollect_length2d_"            + cutsta_[6]]->Fill(z_2d);
            MuCollect1D_["MuCollect_length3d_"            + cutsta_[6]]->Fill(z_3d);
        }
                
        if (nrechits > 20 && (NrecHitXYP == 0 && NrecHitXYN == 0) && (NrecHitYP > 0 && NrecHitYN > 0) && NrecHitZ > 0)
        {
            Nmuonloosenotheta++;
            MuCollect1I_["MuCollect_Nrechits_"            + cutsta_[7]]->Fill(nrechits);      
            MuCollect1I_["MuCollect_date_"                + cutsta_[7]]->Fill(dif_time);
            MuCollect1F_["MuCollect_run_"                 + cutsta_[7]]->Fill(run_nr);
            MuCollect1D_["MuCollect_d0_"                  + cutsta_[7]]->Fill(staTrack->d0());
            MuCollect1D_["MuCollect_dxy_"                 + cutsta_[7]]->Fill(staTrack->dxy());         
            MuCollect1D_["MuCollect_dz_"                  + cutsta_[7]]->Fill(staTrack->dz());         
            MuCollect1D_["MuCollect_dsz_"                 + cutsta_[7]]->Fill(staTrack->dsz());         
            MuCollect1D_["MuCollect_chi2_"                + cutsta_[7]]->Fill(staTrack->chi2());         
            MuCollect1D_["MuCollect_normalizedchi2_"      + cutsta_[7]]->Fill(staTrack->normalizedChi2());
            MuCollect1D_["MuCollect_eta_"                 + cutsta_[7]]->Fill(staTrack->eta());           
            MuCollect1D_["MuCollect_theta_"               + cutsta_[7]]->Fill(staTrack->theta());           
            MuCollect1D_["MuCollect_delta_theta_"         + cutsta_[7]]->Fill(deltaTheta);
            MuCollect1D_["MuCollect_abs_delta_theta_"     + cutsta_[7]]->Fill(abs(deltaTheta));
            MuCollect1D_["MuCollect_phi_"                 + cutsta_[7]]->Fill(staTrack->phi());         
            MuCollect1D_["MuCollect_pt_"                  + cutsta_[7]]->Fill(staTrack->pt());          
            MuCollect1D_["MuCollect_delta_p_"             + cutsta_[7]]->Fill(delta_p);
            MuCollect2D_["MuCollect_PhiVsEta_"            + cutsta_[7]]->Fill(staTrack->eta(), staTrack->phi());      
            for (trackingRecHit_iterator it = track.recHitsBegin ();  it != track.recHitsEnd (); it++)
            {
                if ((*it)->isValid ())
                {
                    const GeomDet* geomDet = theTrackingGeometry->idToDet((*it)->geographicalId());
                    recHitsGlobX = geomDet->toGlobal((*it)->localPosition()).x();
                    recHitsGlobY = geomDet->toGlobal((*it)->localPosition()).y();
                    recHitsGlobZ = geomDet->toGlobal((*it)->localPosition()).z();
                    MuCollect2D_["MuCollect_recHitsGlobPYvsX_"    + cutsta_[7]]->Fill(recHitsGlobX,recHitsGlobY);
                    MuCollect2D_["MuCollect_recHitsGlobPXvsZ_"    + cutsta_[7]]->Fill(recHitsGlobZ,recHitsGlobX);
                    MuCollect2D_["MuCollect_recHitsGlobPYvsZ_"    + cutsta_[7]]->Fill(recHitsGlobZ,recHitsGlobY);
                }
            }            
            MuCollect2D_["MuCollect_recHitsGlob_ip_YvsX_" + cutsta_[7]]->Fill(ip.x(),ip.y());
            MuCollect2D_["MuCollect_recHitsGlob_ip_XvsZ_" + cutsta_[7]]->Fill(ip.z(),ip.x());
            MuCollect2D_["MuCollect_recHitsGlob_ip_YvsZ_" + cutsta_[7]]->Fill(ip.z(),ip.y());
            MuCollect2D_["MuCollect_recHitsGlob_op_YvsX_" + cutsta_[7]]->Fill(op.x(),op.y());
            MuCollect2D_["MuCollect_recHitsGlob_op_XvsZ_" + cutsta_[7]]->Fill(op.z(),op.x());
            MuCollect2D_["MuCollect_recHitsGlob_op_YvsZ_" + cutsta_[7]]->Fill(op.z(),op.y());

            MuCollect1D_["MuCollect_length2d_"            + cutsta_[7]]->Fill(z_2d);
            MuCollect1D_["MuCollect_length3d_"            + cutsta_[7]]->Fill(z_3d);
        }
 //       if (nrechits > 20 && (NrecHitXYP == 0 && NrecHitXYN == 0) && (NrecHitYP > 0 && NrecHitYN > 0) && NrecHitZ > 0)
        if (nrechits > 20 && (abs(deltaTheta) < 0.25) && (NrecHitXYP == 0 && NrecHitXYN == 0) && (NrecHitYP > 0 && NrecHitYN > 0) && NrecHitZ > 0)
        {
            Nmuonloose++;
            MuCollect1I_["MuCollect_Nrechits_"            + cutsta_[8]]->Fill(nrechits);      
            MuCollect1I_["MuCollect_date_"                + cutsta_[8]]->Fill(dif_time);
            MuCollect1F_["MuCollect_run_"                 + cutsta_[8]]->Fill(run_nr);
            MuCollect1D_["MuCollect_d0_"                  + cutsta_[8]]->Fill(staTrack->d0());
            MuCollect1D_["MuCollect_dxy_"                 + cutsta_[8]]->Fill(staTrack->dxy());         
            MuCollect1D_["MuCollect_dz_"                  + cutsta_[8]]->Fill(staTrack->dz());         
            MuCollect1D_["MuCollect_dsz_"                 + cutsta_[8]]->Fill(staTrack->dsz());         
            MuCollect1D_["MuCollect_chi2_"                + cutsta_[8]]->Fill(staTrack->chi2());         
            MuCollect1D_["MuCollect_normalizedchi2_"      + cutsta_[8]]->Fill(staTrack->normalizedChi2());
            MuCollect1D_["MuCollect_eta_"                 + cutsta_[8]]->Fill(staTrack->eta());           
            MuCollect1D_["MuCollect_theta_"               + cutsta_[8]]->Fill(staTrack->theta());           
            MuCollect1D_["MuCollect_delta_theta_"         + cutsta_[8]]->Fill(deltaTheta);
            MuCollect1D_["MuCollect_abs_delta_theta_"     + cutsta_[8]]->Fill(abs(deltaTheta));
            MuCollect1D_["MuCollect_phi_"                 + cutsta_[8]]->Fill(staTrack->phi());         
            MuCollect1D_["MuCollect_pt_"                  + cutsta_[8]]->Fill(staTrack->pt());          
            MuCollect1D_["MuCollect_delta_p_"             + cutsta_[8]]->Fill(delta_p);
            MuCollect2D_["MuCollect_PhiVsEta_"            + cutsta_[8]]->Fill(staTrack->eta(), staTrack->phi());      
            for (trackingRecHit_iterator it = track.recHitsBegin ();  it != track.recHitsEnd (); it++)
            {
                if ((*it)->isValid ())
                {
                    const GeomDet* geomDet = theTrackingGeometry->idToDet((*it)->geographicalId());
                    recHitsGlobX = geomDet->toGlobal((*it)->localPosition()).x();
                    recHitsGlobY = geomDet->toGlobal((*it)->localPosition()).y();
                    recHitsGlobZ = geomDet->toGlobal((*it)->localPosition()).z();
                    MuCollect2D_["MuCollect_recHitsGlobPYvsX_"    + cutsta_[8]]->Fill(recHitsGlobX,recHitsGlobY);
                    MuCollect2D_["MuCollect_recHitsGlobPXvsZ_"    + cutsta_[8]]->Fill(recHitsGlobZ,recHitsGlobX);
                    MuCollect2D_["MuCollect_recHitsGlobPYvsZ_"    + cutsta_[8]]->Fill(recHitsGlobZ,recHitsGlobY);
                }
            }
            
            MuCollect2D_["MuCollect_recHitsGlob_ip_YvsX_" + cutsta_[8]]->Fill(ip.x(),ip.y());
            MuCollect2D_["MuCollect_recHitsGlob_ip_XvsZ_" + cutsta_[8]]->Fill(ip.z(),ip.x());
            MuCollect2D_["MuCollect_recHitsGlob_ip_YvsZ_" + cutsta_[8]]->Fill(ip.z(),ip.y());
            MuCollect2D_["MuCollect_recHitsGlob_op_YvsX_" + cutsta_[8]]->Fill(op.x(),op.y());
            MuCollect2D_["MuCollect_recHitsGlob_op_XvsZ_" + cutsta_[8]]->Fill(op.z(),op.x());
            MuCollect2D_["MuCollect_recHitsGlob_op_YvsZ_" + cutsta_[8]]->Fill(op.z(),op.y());

            MuCollect1D_["MuCollect_length2d_"            + cutsta_[8]]->Fill(z_2d);
            MuCollect1D_["MuCollect_length3d_"            + cutsta_[8]]->Fill(z_3d);
            
            int tmploose = -1;
            int notmploose = -1;
            bool passloose = filter(iEvent, iSetup, theRadius_SiTrk, theMaxZ_SiTrk);
            if (passloose) tmploose = 1;
            else {
                  notmploose = 1;
            }
            passloose = filter(iEvent, iSetup, theRadius_Pixel, theMaxZ_Pixel);
            if (passloose) tmploose = 2;
            else {
                  notmploose = 2;
                  fileNoPixel << "NO MuCollect Track Pixel "<<"   date GMT = " <<dates<<"  Run = " << run_nr <<"  Event = " << ev_nr << std::endl; 
            }
            passloose = filter(iEvent, iSetup, theRadius_IP, theMaxZ_IP);
            if (passloose) tmploose = 3;
            else {
                  notmploose = 3;
            }
            if (flagloose <= tmploose) flagloose = tmploose;
            if (noflagloose <= notmploose) noflagloose = notmploose;
           fileshaft1 <<" flagloose   =  "<< flagloose   <<"  Run = " << run_nr <<"  Event = " << ev_nr  <<std::endl;
           fileshaft1 <<" noflagloose =  "<< noflagloose <<"  Run = " << run_nr <<"  Event = " << ev_nr  <<std::endl;
           cout<<" flagloose   =  "<< flagloose  <<"  Run = " << run_nr <<"  Event = " << ev_nr  <<endl;
           cout<<" noflagloose =  "<< noflagloose<<"  Run = " << run_nr <<"  Event = " << ev_nr  <<endl;

        }
        MuCollect_NDThits->Fill(nDThits);
        MuCollect_NCSChits->Fill(nCSChits);
        MuCollect_NRPChits->Fill(nRPChits);
        MuCollect1I_["MuCollect_Nrechits_"            + cutsta_[0]]->Fill(nrechits);      
        MuCollect1I_["MuCollect_date_"                + cutsta_[0]]->Fill(dif_time);
        MuCollect1F_["MuCollect_run_"                 + cutsta_[0]]->Fill(run_nr);
        MuCollect1D_["MuCollect_d0_"                  + cutsta_[0]]->Fill(staTrack->d0());
        MuCollect1D_["MuCollect_dxy_"                 + cutsta_[0]]->Fill(staTrack->dxy());         
        MuCollect1D_["MuCollect_dz_"                  + cutsta_[0]]->Fill(staTrack->dz());         
        MuCollect1D_["MuCollect_dsz_"                 + cutsta_[0]]->Fill(staTrack->dsz());         
        MuCollect1D_["MuCollect_chi2_"                + cutsta_[0]]->Fill(staTrack->chi2());         
        MuCollect1D_["MuCollect_normalizedchi2_"      + cutsta_[0]]->Fill(staTrack->normalizedChi2());
        MuCollect1D_["MuCollect_eta_"                 + cutsta_[0]]->Fill(staTrack->eta());           
        MuCollect1D_["MuCollect_theta_"               + cutsta_[0]]->Fill(staTrack->theta());           
        MuCollect1D_["MuCollect_delta_theta_"         + cutsta_[0]]->Fill(deltaTheta);
        MuCollect1D_["MuCollect_abs_delta_theta_"     + cutsta_[0]]->Fill(abs(deltaTheta));
        MuCollect1D_["MuCollect_phi_"                 + cutsta_[0]]->Fill(staTrack->phi());         
        MuCollect1D_["MuCollect_pt_"                  + cutsta_[0]]->Fill(staTrack->pt());          
        MuCollect1D_["MuCollect_delta_p_"             + cutsta_[0]]->Fill(delta_p);
        MuCollect2D_["MuCollect_PhiVsEta_"            + cutsta_[0]]->Fill(staTrack->eta(), staTrack->phi());      
        for (trackingRecHit_iterator it = track.recHitsBegin ();  it != track.recHitsEnd (); it++)
        {
            if ((*it)->isValid ())
            {
                const GeomDet* geomDet = theTrackingGeometry->idToDet((*it)->geographicalId());
                recHitsGlobX = geomDet->toGlobal((*it)->localPosition()).x();
                recHitsGlobY = geomDet->toGlobal((*it)->localPosition()).y();
                recHitsGlobZ = geomDet->toGlobal((*it)->localPosition()).z();
                MuCollect2D_["MuCollect_recHitsGlobPYvsX_"    + cutsta_[0]]->Fill(recHitsGlobX,recHitsGlobY);
                MuCollect2D_["MuCollect_recHitsGlobPXvsZ_"    + cutsta_[0]]->Fill(recHitsGlobZ,recHitsGlobX);
                MuCollect2D_["MuCollect_recHitsGlobPYvsZ_"    + cutsta_[0]]->Fill(recHitsGlobZ,recHitsGlobY);
             }
        }
            
        MuCollect2D_["MuCollect_recHitsGlob_ip_YvsX_" + cutsta_[0]]->Fill(ip.x(),ip.y());
        MuCollect2D_["MuCollect_recHitsGlob_ip_XvsZ_" + cutsta_[0]]->Fill(ip.z(),ip.x());
        MuCollect2D_["MuCollect_recHitsGlob_ip_YvsZ_" + cutsta_[0]]->Fill(ip.z(),ip.y());
        MuCollect2D_["MuCollect_recHitsGlob_op_YvsX_" + cutsta_[0]]->Fill(op.x(),op.y());
        MuCollect2D_["MuCollect_recHitsGlob_op_XvsZ_" + cutsta_[0]]->Fill(op.z(),op.x());
        MuCollect2D_["MuCollect_recHitsGlob_op_YvsZ_" + cutsta_[0]]->Fill(op.z(),op.y());

        MuCollect1D_["MuCollect_length2d_"            + cutsta_[0]]->Fill(z_2d);
        MuCollect1D_["MuCollect_length3d_"            + cutsta_[0]]->Fill(z_3d);
 
        MuCollect_innerPositionX->Fill(staTrack->innerPosition().x());
        MuCollect_innerPositionY->Fill(staTrack->innerPosition().y());
        MuCollect_innerPositionZ->Fill(staTrack->innerPosition().z());
        MuCollect_outerPositionX->Fill(staTrack->outerPosition().x());
        MuCollect_outerPositionY->Fill(staTrack->outerPosition().y());
        MuCollect_outerPositionZ->Fill(staTrack->outerPosition().z());

        MuCollect_innerMomentumX->Fill(staTrack->innerMomentum().x());
        MuCollect_innerMomentumY->Fill(staTrack->innerMomentum().y());
        MuCollect_innerMomentumZ->Fill(staTrack->innerMomentum().z());
        MuCollect_outerMomentumX->Fill(staTrack->outerMomentum().x());
        MuCollect_outerMomentumY->Fill(staTrack->outerMomentum().y());
        MuCollect_outerMomentumZ->Fill(staTrack->outerMomentum().z());
        
        MuCollect_time_date0->Fill(dif_time);
        MuCollect_time_date->Fill(dif_time);
        MuCollect_time_date->GetXaxis()->SetTimeDisplay(1);
        MuCollect_time_date->GetXaxis()->SetLabelSize(0.03);
//        MuCollect_test->GetXaxis()->SetTimeFormat("%d\/%m\/%y%F2008-07-12 00:00:00");
        MuCollect_time_date->GetXaxis()->SetTimeFormat("#splitline{%H\:%M\:%S}{%d\/%m\/%y}%F2008-07-01 12:00:00");
   
        MuCollect_time_date1->Fill(dif_time);  
        MuCollect_time_date1->GetXaxis()->SetTimeDisplay(1);
        MuCollect_time_date1->GetXaxis()->SetLabelSize(0.03);
   
        TDatime T0(yer_,mon_,day_,hor_+12,min_,sec_);
        MuCollect_time_date1->GetXaxis()->SetTimeOffset(T0.Convert(),"local");
        MuCollect_time_date1->GetXaxis()->SetTimeFormat("#splitline{%H\:%M\:%S}{%d\/%m\/%y}");
        
        fileout<<" innerMomentum  "<<"px = "<<staTrack->innerMomentum().x()<<" py = "<<staTrack->innerMomentum().y()<<" pz = "<<staTrack->innerMomentum().z() << std::endl;
        fileout<<" outerMomentum  "<<"px = "<<staTrack->outerMomentum().x()<<" py = "<<staTrack->outerMomentum().y()<<" pz = "<<staTrack->outerMomentum().z() << std::endl;
               

//        MuCollect_delta_p->Fill(delta_p);    
        
        MuCollect_vx->Fill(staTrack->vx());
        MuCollect_vy->Fill(staTrack->vy());
        MuCollect_vz->Fill(staTrack->vz());
                      
        Y0 = 8950 -  staTrack->vy();
        X0 = 500 - staTrack->vx();
        
        Z1 = -1800 -  staTrack->vz();
        Z2 = Z1 - 1000;
        X1 = 500 - staTrack->vx();
        X2 = -500 - staTrack->vx();
//        if((staTrack->py() < 0) && (staTrack->pz() > 0))  Zoy = Y0*(staTrack->pz()/staTrack->py());               
//        if((staTrack->pz() > 0))  Zox = X0*(staTrack->pz()/staTrack->px());
        Zoy = Y0*(staTrack->pz()/staTrack->py());               
        Zox = X0*(staTrack->pz()/staTrack->px());
       
//          std::cout << " Y0 = 8950 - vy = " << Y0  << std::endl;
//          std::cout << " X0 = 500 -|vx| = " << X0  << std::endl;
//          std::cout << " Z1 = -1800 - vz = " << Z1  << std::endl;
//          std::cout << " Z2 = Z1 - 1000 = " << Z2  << std::endl;
//          std::cout << " Zox = X0*pz/px = " << Zox << std::endl;
//          std::cout << " Zoy = Y0*pz/py = " << Zoy << std::endl;
 
  //      double x0 = 8950 -  staTrack->vx();
        Xoy = Y0*(staTrack->px()/staTrack->py());
        
//          std::cout << " Xox =  " << Xoy << std::endl;
//          std::cout << " zoy =  " << Zoy << std::endl;
          
        MuCollect_surfaceXZ->Fill(Zoy, Xoy);
         
        if((Zoy <= Z1) && (Zoy >= Z2) && (Xoy <= X1) && (Xoy >= X2))
        {
        //    MuCollect_nRecHits_shaft->Fill(nrechits);
        //    MuCollect_normalizedchi2_shaft->Fill(staTrack->normalizedChi2());
       
            MuCollect_shaftXZ1->Fill(Zoy, Xoy);
            
//             fileshaft1<< " Y0 = 9100 - vy  = " << Y0  << std::endl;
//             fileshaft1<< " X0 = 500 - |vx| = " << X0  << std::endl;
//             fileshaft1<< " Z1 = -1800 - vz = " << Z1  << std::endl;
//             fileshaft1<< " Z2 =  Z1 - 1000 = " << Z2  << std::endl;
//             fileshaft1<< " X1 = 500 - vx = " << X1  << std::endl;
//             fileshaft1<< " X2 =  -500 - vx = " << X2  << std::endl;
//             fileshaft1<< " Zox = X0*pz/px  = " << Zox << std::endl;
//             fileshaft1<< " Zoy = Y0*pz/py  = " << Zoy << std::endl;
//             fileshaft1<< " Xoy = Y0*px/py  = " << Xoy << std::endl;
//             fileshaft1<< " px = " << staTrack->px() << std::endl;
//             fileshaft1<< " py = " << staTrack->py() << std::endl;
//             fileshaft1<< " pz = " << staTrack->pz() << std::endl;
//             fileshaft1<< " pt = " << staTrack->pt() << std::endl;
//             fileshaft1<<"            Muons for shaft  :"<<" Run = "<< run_nr <<"  Event = " << ev_nr <<"   date GMT = " <<dates<< std::endl;
//             fileshaft1<< " " << std::endl;
        }       
        if((Zoy < -1800) && (Zoy > -2800) && (Xoy < 500) && (Xoy > -500))
        {
            MuCollect_nRecHits_shaft->Fill(nrechits);
            MuCollect_normalizedchi2_shaft->Fill(staTrack->normalizedChi2());
       
            MuCollect_shaftXZ->Fill(Zoy, Xoy);
            
            fileshaft<< " Y0 = 9100 - vy  = " << Y0  << std::endl;
            fileshaft<< " X0 = 500 - |vx| = " << X0  << std::endl;
            fileshaft<< " X1 = 500 - vx = " << X1  << std::endl;
            fileshaft<< " X2 =  -500 - vx = " << X2  << std::endl;
            fileshaft<< " Zox = X0*pz/px  = " << Zox << std::endl;
            fileshaft<< " Zoy = Y0*pz/py  = " << Zoy << std::endl;
            fileshaft<< " Xoy = Y0*px/py  = " << Xoy << std::endl;
            fileshaft<< " px = " << staTrack->px() << std::endl;
            fileshaft<< " py = " << staTrack->py() << std::endl;
            fileshaft<< " pz = " << staTrack->pz() << std::endl;
            fileshaft<< " pt = " << staTrack->pt() << std::endl;
            fileshaft<<"            Muons for shaft  :"<<" Run = "<< run_nr <<"  Event = " << ev_nr <<"   date GMT = " <<dates<< std::endl;
            fileshaft<< " " << std::endl;
        }       
        
//        MuCollect_vx->Fill(staTrack->vertex().x());
//        MuCollect_vx->Fill(staTrack->vertex().X());
        
    }//end loop on the Muons/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//    cout<<" loose =  "<<Nmuonloose <<endl;

    if (Nmuonloose > nmuon_)
    {
 //   cout<<" Nmuonloose =  "<< Nmuonloose <<endl;
        int tmploosemuon1 = -1;
        bool passloose1 = filter(iEvent, iSetup, theRadius_SiTrk, theMaxZ_SiTrk);
        if (passloose1) 
        {
        tmploosemuon1 = 1;
        fileflagloosmuon<<" loosmuon  SiTrk  "<< tmploosemuon1 <<"  run = "<< run_nr <<"  Event = " << ev_nr << std::endl;
        fileloosemuonlist<<run_nr <<"  " << ev_nr << std::endl;
        }
        passloose1 = filter(iEvent, iSetup, theRadius_Pixel, theMaxZ_Pixel);
        if (passloose1) 
        {
        tmploosemuon1 = 2;
        fileflagloosmuon<<" loosmuon  Pixel  "<< tmploosemuon1 <<"  run = "<< run_nr <<"  Event = " << ev_nr << std::endl;
        fileloosemuonlist<<run_nr <<"  " << ev_nr << std::endl;
        }
        passloose1 = filter(iEvent, iSetup, theRadius_IP, theMaxZ_IP);
        if (passloose1) 
        {
        tmploosemuon1 = 3;
        fileflagloosmuon<<" loosmuon  IP  "<< tmploosemuon1 <<"  run = "<< run_nr <<"  Event = " << ev_nr << std::endl;
        fileloosemuonlist<<run_nr <<"  " << ev_nr << std::endl;
        }
        if (flagloosemuon1 <= tmploosemuon1) flagloosemuon1 = tmploosemuon1;

        Run_Good_muon->Fill(run_nr);
 //   cout<<" flagloosemuon1 =  "<< flagloosemuon1 <<endl;
//    fileflagloosmuon<<" loosmuon  "<< run_nr <<"  Event = " << ev_nr << std::endl;
    }

    MuCollect1I_["MuCollect_Nmuon_" + cutsta_[1]]->Fill(NmuonlooseRechit);
    MuCollect1I_["MuCollect_Nmuon_" + cutsta_[2]]->Fill(NmuonlooseXY);
    MuCollect1I_["MuCollect_Nmuon_" + cutsta_[3]]->Fill(NmuonlooseY);
    MuCollect1I_["MuCollect_Nmuon_" + cutsta_[4]]->Fill(NmuonlooseZ);
    MuCollect1I_["MuCollect_Nmuon_" + cutsta_[5]]->Fill(NmuonlooseYZ_x);
    MuCollect1I_["MuCollect_Nmuon_" + cutsta_[6]]->Fill(NmuonlooseTheta);
    MuCollect1I_["MuCollect_Nmuon_" + cutsta_[7]]->Fill(Nmuonloosenotheta);
    MuCollect1I_["MuCollect_Nmuon_" + cutsta_[8]]->Fill(Nmuonloose);
    
    Flagloose->Fill(flagloose);
    Flagloosemuon1->Fill(flagloosemuon1);
    
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    double chi2prob, pxerror, pyerror, pzerror, pterror, perror;

    for  (int j = 0; j < 3; ++j)
    {

        Trk_[collection_[j] + cut_[0] + "NumTracks"]->Fill(typeTrackCollection[j]->size());

        if (typeTrackCollection[j]->size() > 1)
        {
            std::cout << " #### " << collection_[j] << " TRACKS " << " >1: " << typeTrackCollection[j]->size() <<"  Run = " << run_nr <<"  Event = " <<
             ev_nr <<"   date GMT = " <<dates << "Diferencia local: " << dif_time <<"s"<< std::endl;
        }
        for ( reco::TrackCollection::const_iterator track = typeTrackCollection[j]->begin(); track != typeTrackCollection[j]->end(); ++track )
        {
            pterror  = (track->pt()) ? track->ptError() / track->pt() : 0.0;
            pxerror  = -1.0;
            pyerror  = -1.0;
            pzerror  = -1.0;
            perror   = -1.0;
            
            Trk_[collection_[j] + cut_[0] + "date"]->Fill(dif_time);
            Trk_[collection_[j] + cut_[0] + "Run"]->Fill(run_nr);            
            Trk_[collection_[j] + cut_[0] + "d0"]->Fill(track->d0());            
            //chi2prob= TMath::Prob(track->chi2(),track->ndof());
            Trk_[collection_[j] + cut_[0] + "Chi2"]->Fill(track->chi2());
            Trk_[collection_[j] + cut_[0] + "Chi2Prob"]->Fill(TMath::Prob(track->chi2(), track->ndof()));

            Trk_[collection_[j] + cut_[0] + "dxy"]->Fill(track->dxy());               
            Trk_[collection_[j] + cut_[0] + "dz"]->Fill(track->dz());               
            Trk_[collection_[j] + cut_[0] + "dsz"]->Fill(track->dsz()); 
            Trk_[collection_[j] + cut_[0] + "pt"]->Fill(track->pt());
            Trk_[collection_[j] + cut_[0] + "px"]->Fill(track->px());
            Trk_[collection_[j] + cut_[0] + "py"]->Fill(track->py());
            Trk_[collection_[j] + cut_[0] + "pz"]->Fill(track->pz());
            Trk_[collection_[j] + cut_[0] + "p"]->Fill(track->p());
            Trk_[collection_[j] + cut_[0] + "eta"]->Fill(track->eta());
            Trk_[collection_[j] + cut_[0] + "phi"]->Fill(track->phi());
            Trk_[collection_[j] + cut_[0] + "theta"]->Fill(track->theta());
            Trk_[collection_[j] + cut_[0] + "NumRecHitsPerTrk"]->Fill(track->recHitsSize());
            
            Trk2D_[collection_[j] + cut_[0] + "NumRecHitsPerTrkVsPhi"]->Fill(track->phi(), track->recHitsSize());
            Trk2D_[collection_[j] + cut_[0] + "NumRecHitsPerTrkVsTheta"]->Fill(track->theta(), track->recHitsSize());
            Trk2D_[collection_[j] + cut_[0] + "NumRecHitsPerTrkVsEta"]->Fill(track->eta(), track->recHitsSize());

            Trk2D_[collection_[j] + cut_[0] + "Chi2overDoFVsTheta"]->Fill(track->theta(), track->normalizedChi2());
            Trk2D_[collection_[j] + cut_[0] + "Chi2overDoFVsPhi"]->Fill(track->phi(), track->normalizedChi2());
            Trk2D_[collection_[j] + cut_[0] + "Chi2overDoFVsEta"]->Fill(track->eta(), track->normalizedChi2());

            Trk_[collection_[j] + cut_[0] + "Errdxy"]->Fill(track->dxyError());
            Trk_[collection_[j] + cut_[0] + "Errdz"]->Fill(track->dzError());
            Trk_[collection_[j] + cut_[0] + "TrackPtErr"]->Fill(pterror);
            Trk_[collection_[j] + cut_[0] + "TrackPxErr"]->Fill(pxerror);
            Trk_[collection_[j] + cut_[0] + "TrackPyErr"]->Fill(pyerror);
            Trk_[collection_[j] + cut_[0] + "TrackPzErr"]->Fill(pzerror);
            Trk_[collection_[j] + cut_[0] + "TrackPErr"]->Fill(perror);
            Trk_[collection_[j] + cut_[0] + "TrackPhiErr"]->Fill(track->phiError());
            Trk_[collection_[j] + cut_[0] + "TrackThetaErr"]->Fill(track->thetaError());
            Trk_[collection_[j] + cut_[0] + "TrackEtaErr"]->Fill(track->etaError());
        }
    }

    for (int j = 0; j < 3; ++j)
    {
        if (typeTrackCollection[j]->size() > 1)
        {
            std::cout << "  " << collection_[j] << " TRACKS " << " >1: " << typeTrackCollection[j]->size() << "  " << run_nr << " / " << ev_nr << std::endl;
        }

//        std::cout << " flag " << flag << " flagloose  " << flagloose << "u" << u << std::endl;

        for (int k = 1; k < 4; ++k)
        {
            //    std::cout<< " flag0 "<< flag<<"   " <<k<<std::endl;
            if (flag == k)
            {
                int collec_cnt_flag = 0;
//                 std::cout << " flagloop " << flag << k << std::endl;

                for ( reco::TrackCollection::const_iterator track = typeTrackCollection[j]->begin(); track != typeTrackCollection[j]->end(); ++track )
                {

                    collec_cnt_flag++;

                    pterror  = (track->pt()) ? track->ptError() / track->pt() : 0.0;
                    pxerror  = -1.0;
                    pyerror  = -1.0;
                    pzerror  = -1.0;
                    perror   = -1.0;

                    Trk_[collection_[j] + cut_[k] + "date"]->Fill(dif_time);
                    Trk_[collection_[j] + cut_[k] + "Run"]->Fill(run_nr);            
                    Trk_[collection_[j] + cut_[k] + "d0"]->Fill(track->d0());            
                    //chi2prob= TMath::Prob(track->chi2(),track->ndof());
                    Trk_[collection_[j] + cut_[k] + "Chi2"]->Fill(track->chi2());
                    Trk_[collection_[j] + cut_[k] + "Chi2Prob"]->Fill(TMath::Prob(track->chi2(), track->ndof()));

                    Trk_[collection_[j] + cut_[k] + "dxy"]->Fill(track->dxy());               
                    Trk_[collection_[j] + cut_[k] + "dz"]->Fill(track->dz());               
                    Trk_[collection_[j] + cut_[k] + "dsz"]->Fill(track->dsz()); 
                    Trk_[collection_[j] + cut_[k] + "pt"]->Fill(track->pt());
                    Trk_[collection_[j] + cut_[k] + "px"]->Fill(track->px());
                    Trk_[collection_[j] + cut_[k] + "py"]->Fill(track->py());
                    Trk_[collection_[j] + cut_[k] + "pz"]->Fill(track->pz());
                    Trk_[collection_[j] + cut_[k] + "p"]->Fill(track->p());
                    Trk_[collection_[j] + cut_[k] + "eta"]->Fill(track->eta());
                    Trk_[collection_[j] + cut_[k] + "phi"]->Fill(track->phi());
                    Trk_[collection_[j] + cut_[k] + "theta"]->Fill(track->theta());
                    Trk_[collection_[j] + cut_[k] + "NumRecHitsPerTrk"]->Fill(track->recHitsSize());

                    Trk2D_[collection_[j] + cut_[k] + "NumRecHitsPerTrkVsPhi"]->Fill(track->phi(), track->recHitsSize());
                    Trk2D_[collection_[j] + cut_[k] + "NumRecHitsPerTrkVsTheta"]->Fill(track->theta(), track->recHitsSize());
                    Trk2D_[collection_[j] + cut_[k] + "NumRecHitsPerTrkVsEta"]->Fill(track->eta(), track->recHitsSize());

                    Trk2D_[collection_[j] + cut_[k] + "Chi2overDoFVsTheta"]->Fill(track->theta(), track->normalizedChi2());
                    Trk2D_[collection_[j] + cut_[k] + "Chi2overDoFVsPhi"]->Fill(track->phi(), track->normalizedChi2());
                    Trk2D_[collection_[j] + cut_[k] + "Chi2overDoFVsEta"]->Fill(track->eta(), track->normalizedChi2());

                    Trk_[collection_[j] + cut_[k] + "Errdxy"]->Fill(track->dxyError());
                    Trk_[collection_[j] + cut_[k] + "Errdz"]->Fill(track->dzError());
                    Trk_[collection_[j] + cut_[k] + "TrackPtErr"]->Fill(pterror);
                    Trk_[collection_[j] + cut_[k] + "TrackPxErr"]->Fill(pxerror);
                    Trk_[collection_[j] + cut_[k] + "TrackPyErr"]->Fill(pyerror);
                    Trk_[collection_[j] + cut_[k] + "TrackPzErr"]->Fill(pzerror);
                    Trk_[collection_[j] + cut_[k] + "TrackPErr"]->Fill(perror);
                    Trk_[collection_[j] + cut_[k] + "TrackPhiErr"]->Fill(track->phiError());
                    Trk_[collection_[j] + cut_[k] + "TrackThetaErr"]->Fill(track->thetaError());
                    Trk_[collection_[j] + cut_[k] + "TrackEtaErr"]->Fill(track->etaError());
                }
//                std::cout << collection_[j] << " collec_cnt_flag" << k << " = " << collec_cnt_flag  << std::endl;
                Trk_[collection_[j] + cut_[k] + "NumTracks"]->Fill(collec_cnt_flag);
            }

            if (flagloose == k)
            {
                int collec_cnt_flagloose = 0;
//                 std::cout << "  flagloose   " <<flagloose<< k << std::endl;

                for ( reco::TrackCollection::const_iterator track = typeTrackCollection[j]->begin(); track != typeTrackCollection[j]->end(); ++track )
                {

                    collec_cnt_flagloose++;

                    pterror  = (track->pt()) ? track->ptError() / track->pt() : 0.0;
                    pxerror  = -1.0;
                    pyerror  = -1.0;
                    pzerror  = -1.0;
                    perror   = -1.0;

                    Trk_[collection_[j] + cut_[k+3] + "date"]->Fill(dif_time);
                    Trk_[collection_[j] + cut_[k+3] + "Run"]->Fill(run_nr);            
                    Trk_[collection_[j] + cut_[k+3] + "d0"]->Fill(track->d0());            
                    //chi2prob= TMath::Prob(track->chi2(),track->ndof());
                    Trk_[collection_[j] + cut_[k+3] + "Chi2"]->Fill(track->chi2());
                    Trk_[collection_[j] + cut_[k+3] + "Chi2Prob"]->Fill(TMath::Prob(track->chi2(), track->ndof()));

                    Trk_[collection_[j] + cut_[k+3] + "dxy"]->Fill(track->dxy());               
                    Trk_[collection_[j] + cut_[k+3] + "dz"]->Fill(track->dz());               
                    Trk_[collection_[j] + cut_[k+3] + "dsz"]->Fill(track->dsz()); 
                    Trk_[collection_[j] + cut_[k+3] + "pt"]->Fill(track->pt());
                    Trk_[collection_[j] + cut_[k+3] + "px"]->Fill(track->px());
                    Trk_[collection_[j] + cut_[k+3] + "py"]->Fill(track->py());
                    Trk_[collection_[j] + cut_[k+3] + "pz"]->Fill(track->pz());
                    Trk_[collection_[j] + cut_[k+3] + "p"]->Fill(track->p());
                    Trk_[collection_[j] + cut_[k+3] + "eta"]->Fill(track->eta());
                    Trk_[collection_[j] + cut_[k+3] + "phi"]->Fill(track->phi());
                    Trk_[collection_[j] + cut_[k+3] + "theta"]->Fill(track->theta());
                    Trk_[collection_[j] + cut_[k+3] + "NumRecHitsPerTrk"]->Fill(track->recHitsSize());

                    Trk2D_[collection_[j] + cut_[k+3] + "NumRecHitsPerTrkVsPhi"]->Fill(track->phi(), track->recHitsSize());
                    Trk2D_[collection_[j] + cut_[k+3] + "NumRecHitsPerTrkVsTheta"]->Fill(track->theta(), track->recHitsSize());
                    Trk2D_[collection_[j] + cut_[k+3] + "NumRecHitsPerTrkVsEta"]->Fill(track->eta(), track->recHitsSize());

                    Trk2D_[collection_[j] + cut_[k+3] + "Chi2overDoFVsTheta"]->Fill(track->theta(), track->normalizedChi2());
                    Trk2D_[collection_[j] + cut_[k+3] + "Chi2overDoFVsPhi"]->Fill(track->phi(), track->normalizedChi2());
                    Trk2D_[collection_[j] + cut_[k+3] + "Chi2overDoFVsEta"]->Fill(track->eta(), track->normalizedChi2());

                    Trk_[collection_[j] + cut_[k+3] + "Errdxy"]->Fill(track->dxyError());
                    Trk_[collection_[j] + cut_[k+3] + "Errdz"]->Fill(track->dzError());
                    Trk_[collection_[j] + cut_[k+3] + "TrackPtErr"]->Fill(pterror);
                    Trk_[collection_[j] + cut_[k+3] + "TrackPxErr"]->Fill(pxerror);
                    Trk_[collection_[j] + cut_[k+3] + "TrackPyErr"]->Fill(pyerror);
                    Trk_[collection_[j] + cut_[k+3] + "TrackPzErr"]->Fill(pzerror);
                    Trk_[collection_[j] + cut_[k+3] + "TrackPErr"]->Fill(perror);
                    Trk_[collection_[j] + cut_[k+3] + "TrackPhiErr"]->Fill(track->phiError());
                    Trk_[collection_[j] + cut_[k+3] + "TrackThetaErr"]->Fill(track->thetaError());
                    Trk_[collection_[j] + cut_[k+3] + "TrackEtaErr"]->Fill(track->etaError());
                }
//                std::cout << collection_[j] << " collec_cnt_flagloose" << k + 3 << " = " << collec_cnt_flagloose  << std::endl;
                Trk_[collection_[j] + cut_[k+3] + "NumTracks"]->Fill(collec_cnt_flagloose);
            }

            if (flagloosemuon1 == k)
            {
                int collec_cnt_flagloosemuon1 = 0;
//                 std::cout << "  flagloose   " <<flagloose<< k << std::endl;

                for ( reco::TrackCollection::const_iterator track = typeTrackCollection[j]->begin(); track != typeTrackCollection[j]->end(); ++track )
                {

                    collec_cnt_flagloosemuon1++;

                    pterror  = (track->pt()) ? track->ptError() / track->pt() : 0.0;
                    pxerror  = -1.0;
                    pyerror  = -1.0;
                    pzerror  = -1.0;
                    perror   = -1.0;

                    Trk_[collection_[j] + cut_[k+6] + "date"]->Fill(dif_time);
                    Trk_[collection_[j] + cut_[k+6] + "Run"]->Fill(run_nr);            
                    Trk_[collection_[j] + cut_[k+6] + "d0"]->Fill(track->d0());            
                     //chi2prob= TMath::Prob(track->chi2(),track->ndof());
                    Trk_[collection_[j] + cut_[k+6] + "Chi2"]->Fill(track->chi2());
                    Trk_[collection_[j] + cut_[k+6] + "Chi2Prob"]->Fill(TMath::Prob(track->chi2(), track->ndof()));

                    Trk_[collection_[j] + cut_[k+6] + "dxy"]->Fill(track->dxy());               
                    Trk_[collection_[j] + cut_[k+6] + "dz"]->Fill(track->dz());               
                    Trk_[collection_[j] + cut_[k+6] + "dsz"]->Fill(track->dsz()); 
                    Trk_[collection_[j] + cut_[k+6] + "pt"]->Fill(track->pt());
                    Trk_[collection_[j] + cut_[k+6] + "px"]->Fill(track->px());
                    Trk_[collection_[j] + cut_[k+6] + "py"]->Fill(track->py());
                    Trk_[collection_[j] + cut_[k+6] + "pz"]->Fill(track->pz());
                    Trk_[collection_[j] + cut_[k+6] + "p"]->Fill(track->p());
                    Trk_[collection_[j] + cut_[k+6] + "eta"]->Fill(track->eta());
                    Trk_[collection_[j] + cut_[k+6] + "phi"]->Fill(track->phi());
                    Trk_[collection_[j] + cut_[k+6] + "theta"]->Fill(track->theta());
                    Trk_[collection_[j] + cut_[k+6] + "NumRecHitsPerTrk"]->Fill(track->recHitsSize());

                    Trk2D_[collection_[j] + cut_[k+6] + "NumRecHitsPerTrkVsPhi"]->Fill(track->phi(), track->recHitsSize());
                    Trk2D_[collection_[j] + cut_[k+6] + "NumRecHitsPerTrkVsTheta"]->Fill(track->theta(), track->recHitsSize());
                    Trk2D_[collection_[j] + cut_[k+6] + "NumRecHitsPerTrkVsEta"]->Fill(track->eta(), track->recHitsSize());

                    Trk2D_[collection_[j] + cut_[k+6] + "Chi2overDoFVsTheta"]->Fill(track->theta(), track->normalizedChi2());
                    Trk2D_[collection_[j] + cut_[k+6] + "Chi2overDoFVsPhi"]->Fill(track->phi(), track->normalizedChi2());
                    Trk2D_[collection_[j] + cut_[k+6] + "Chi2overDoFVsEta"]->Fill(track->eta(), track->normalizedChi2());

                    Trk_[collection_[j] + cut_[k+6] + "Errdxy"]->Fill(track->dxyError());
                    Trk_[collection_[j] + cut_[k+6] + "Errdz"]->Fill(track->dzError());
                    Trk_[collection_[j] + cut_[k+6] + "TrackPtErr"]->Fill(pterror);
                    Trk_[collection_[j] + cut_[k+6] + "TrackPxErr"]->Fill(pxerror);
                    Trk_[collection_[j] + cut_[k+6] + "TrackPyErr"]->Fill(pyerror);
                    Trk_[collection_[j] + cut_[k+6] + "TrackPzErr"]->Fill(pzerror);
                    Trk_[collection_[j] + cut_[k+6] + "TrackPErr"]->Fill(perror);
                    Trk_[collection_[j] + cut_[k+6] + "TrackPhiErr"]->Fill(track->phiError());
                    Trk_[collection_[j] + cut_[k+6] + "TrackThetaErr"]->Fill(track->thetaError());
                    Trk_[collection_[j] + cut_[k+6] + "TrackEtaErr"]->Fill(track->etaError());
                }
//                std::cout << collection_[j] << " collec_cnt_flagloose" << k + 6 << " = " << collec_cnt_flagloose  << std::endl;
                Trk_[collection_[j] + cut_[k+6] + "NumTracks"]->Fill(collec_cnt_flagloosemuon1);
            }

            if (noflagloose == k)
            {
                int collec_cnt_noflagloose = 0;
//                 std::cout << "  flagloose   " <<flagloose<< k << std::endl;

                for ( reco::TrackCollection::const_iterator track = typeTrackCollection[j]->begin(); track != typeTrackCollection[j]->end(); ++track )
                {

                    collec_cnt_noflagloose++;

                    pterror  = (track->pt()) ? track->ptError() / track->pt() : 0.0;
                    pxerror  = -1.0;
                    pyerror  = -1.0;
                    pzerror  = -1.0;
                    perror   = -1.0;

                    Trk_[collection_[j] + cut_[k+9] + "date"]->Fill(dif_time);
                    Trk_[collection_[j] + cut_[k+9] + "Run"]->Fill(run_nr);            
                    Trk_[collection_[j] + cut_[k+9] + "d0"]->Fill(track->d0());            
                     //chi2prob= TMath::Prob(track->chi2(),track->ndof());
                    Trk_[collection_[j] + cut_[k+9] + "Chi2"]->Fill(track->chi2());
                    Trk_[collection_[j] + cut_[k+9] + "Chi2Prob"]->Fill(TMath::Prob(track->chi2(), track->ndof()));

                    Trk_[collection_[j] + cut_[k+9] + "dxy"]->Fill(track->dxy());               
                    Trk_[collection_[j] + cut_[k+9] + "dz"]->Fill(track->dz());               
                    Trk_[collection_[j] + cut_[k+9] + "dsz"]->Fill(track->dsz()); 
                    Trk_[collection_[j] + cut_[k+9] + "pt"]->Fill(track->pt());
                    Trk_[collection_[j] + cut_[k+9] + "px"]->Fill(track->px());
                    Trk_[collection_[j] + cut_[k+9] + "py"]->Fill(track->py());
                    Trk_[collection_[j] + cut_[k+9] + "pz"]->Fill(track->pz());
                    Trk_[collection_[j] + cut_[k+9] + "p"]->Fill(track->p());
                    Trk_[collection_[j] + cut_[k+9] + "eta"]->Fill(track->eta());
                    Trk_[collection_[j] + cut_[k+9] + "phi"]->Fill(track->phi());
                    Trk_[collection_[j] + cut_[k+9] + "theta"]->Fill(track->theta());
                    Trk_[collection_[j] + cut_[k+9] + "NumRecHitsPerTrk"]->Fill(track->recHitsSize());

                    Trk2D_[collection_[j] + cut_[k+9] + "NumRecHitsPerTrkVsPhi"]->Fill(track->phi(), track->recHitsSize());
                    Trk2D_[collection_[j] + cut_[k+9] + "NumRecHitsPerTrkVsTheta"]->Fill(track->theta(), track->recHitsSize());
                    Trk2D_[collection_[j] + cut_[k+9] + "NumRecHitsPerTrkVsEta"]->Fill(track->eta(), track->recHitsSize());

                    Trk2D_[collection_[j] + cut_[k+9] + "Chi2overDoFVsTheta"]->Fill(track->theta(), track->normalizedChi2());
                    Trk2D_[collection_[j] + cut_[k+9] + "Chi2overDoFVsPhi"]->Fill(track->phi(), track->normalizedChi2());
                    Trk2D_[collection_[j] + cut_[k+9] + "Chi2overDoFVsEta"]->Fill(track->eta(), track->normalizedChi2());

                    Trk_[collection_[j] + cut_[k+9] + "Errdxy"]->Fill(track->dxyError());
                    Trk_[collection_[j] + cut_[k+9] + "Errdz"]->Fill(track->dzError());
                    Trk_[collection_[j] + cut_[k+9] + "TrackPtErr"]->Fill(pterror);
                    Trk_[collection_[j] + cut_[k+9] + "TrackPxErr"]->Fill(pxerror);
                    Trk_[collection_[j] + cut_[k+9] + "TrackPyErr"]->Fill(pyerror);
                    Trk_[collection_[j] + cut_[k+9] + "TrackPzErr"]->Fill(pzerror);
                    Trk_[collection_[j] + cut_[k+9] + "TrackPErr"]->Fill(perror);
                    Trk_[collection_[j] + cut_[k+9] + "TrackPhiErr"]->Fill(track->phiError());
                    Trk_[collection_[j] + cut_[k+9] + "TrackThetaErr"]->Fill(track->thetaError());
                    Trk_[collection_[j] + cut_[k+9] + "TrackEtaErr"]->Fill(track->etaError());
                }
//                std::cout << collection_[j] << " collec_cnt_flagloose" << k + 9 << " = " << collec_cnt_flagloose  << std::endl;
                Trk_[collection_[j] + cut_[k+9] + "NumTracks"]->Fill(collec_cnt_noflagloose);
            }
        }
    }
/////////////////////////////////////////////////////////////////////////////////////////////////////
//   rootTree_->Fill();
}


void TrackCruzet::beginJob(const edm::EventSetup& iSetup)
{


}


void TrackCruzet::endJob()
{
    /*   if ( rootTree_ ) {
         rootTree_->GetDirectory()->cd();
         rootTree_->Write();
         delete rootTree_;
       }*/
}


//define this as a plug-in
DEFINE_FWK_MODULE(TrackCruzet);
