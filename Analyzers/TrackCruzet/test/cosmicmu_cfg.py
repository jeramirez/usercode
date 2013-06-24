import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

#process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.load("HCALWork.CosmicMu.CRAFT_V3P_SuperPointing_v3_cff")

process.TFileService = cms.Service("TFileService",
   fileName = cms.string("test.root")
)

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.connect = "frontier://PromptProd/CMS_COND_21X_GLOBALTAG"
process.GlobalTag.globaltag = "CRAFT_V2P::All"
process.prefer("GlobalTag")

process.load("Configuration.StandardSequences.MagneticField_38T_cff")
#process.load("Configuration.StandardSequences.GeometryPilot2_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Geometry.CommonDetUnit.globalTrackingGeometry_cfi")

from TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi import *
from TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAlong_cfi import *
from TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorOpposite_cfi import *

process.load("TrackingTools.TrackAssociator.DetIdAssociatorESProducer_cff")

#process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
#    fileNames = cms.untracked.vstring(
#    'file:../../../HCALWork/METFilter/test/condor/METFLT_66714Cosmics_50_100_1k1M_1.root'
#    'file:/uscmst1b_scratch/lpc1/lpcphys/wclarida/majorana_data/Raw/digireco/majorana200_reco_1_0.root'
#    'file:/uscms/home/wclarida/SIM_PROGS/CMSSW_2_1_10/src/HCALWork/METFilter/test/condor/METFLT_66714Cosmics_50_100_0.root',
#    'file:/uscms/home/wclarida/SIM_PROGS/CMSSW_2_1_10/src/HCALWork/METFilter/test/condor/METFLT_66714Cosmics_50_100_1.root',
#    'file:/uscms/home/wclarida/SIM_PROGS/CMSSW_2_1_10/src/HCALWork/METFilter/test/condor/METFLT_66714Cosmics_50_100_2.root',
#    'file:/uscms/home/wclarida/SIM_PROGS/CMSSW_2_1_10/src/HCALWork/METFilter/test/condor/METFLT_66714Cosmics_50_100_3.root',
#    'file:/uscms/home/wclarida/SIM_PROGS/CMSSW_2_1_10/src/HCALWork/METFilter/test/condor/METFLT_66714Cosmics_50_100_4.root',
#    'file:/uscms/home/wclarida/SIM_PROGS/CMSSW_2_1_10/src/HCALWork/METFilter/test/condor/METFLT_66714Cosmics_50_100_5.root'
#    )
#)
#process.load("TrackingTools.TrackAssociator.default_cfi")

#from RecoMuon.MuonIdentification.muonIdProducerSequence_cff import *

process.demo = cms.EDAnalyzer('CosmicMu',
   Jets = cms.InputTag("iterativeCone5CaloJets"),
   Muons = cms.InputTag("muons"),
   MinMuonPt = cms.untracked.double(0.0),                              
   MinJetPt = cms.untracked.double(1.5),                              
   TrackAssociatorParameters = cms.PSet(
         muonMaxDistanceSigmaX = cms.double(0.0),
         muonMaxDistanceSigmaY = cms.double(0.0),
         CSCSegmentCollectionLabel = cms.InputTag("cscSegments"),
         dRHcal = cms.double(9999.0),
         dREcal = cms.double(9999.0),
         CaloTowerCollectionLabel = cms.InputTag("towerMaker"),
         useEcal = cms.bool(True),
         dREcalPreselection = cms.double(0.05),
         HORecHitCollectionLabel = cms.InputTag("horeco"),
         dRMuon = cms.double(9999.0),
         crossedEnergyType = cms.string('SinglePointAlongTrajectory'),
         propagateAllDirections = cms.bool(True),
         muonMaxDistanceX = cms.double(5.0),
         muonMaxDistanceY = cms.double(5.0),
         useHO = cms.bool(True),
         accountForTrajectoryChangeCalo = cms.bool(False),
         DTRecSegment4DCollectionLabel = cms.InputTag("dt4DSegments"),
         EERecHitCollectionLabel = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
         dRHcalPreselection = cms.double(0.2),
         useMuon = cms.bool(True),
         useCalo = cms.bool(False),
         EBRecHitCollectionLabel = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
         dRMuonPreselection = cms.double(0.2),
         truthMatch = cms.bool(False),
         HBHERecHitCollectionLabel = cms.InputTag("hbhereco"),
         useHcal = cms.bool(True)
     )
)


process.p = cms.Path(process.demo)
