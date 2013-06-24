import FWCore.ParameterSet.Config as cms

process = cms.Process("MyHecTau")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.MagneticField_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.load("Configuration.StandardSequences.MagneticField_cff")
process.source = cms.Source("PoolSource",
    # replace 'myfile.root'' with the source file you want to use
    fileNames = cms.untracked.vstring(
        #'file:myfile.root'   
 '/store/data/Run2011A/BTag//RECO/05Aug2011-v1/0000/FA0870EC-98C0-E011-B5E6-003048D15DB6.root'
,'/store/data/Run2011A/BTag//RECO/05Aug2011-v1/0000/FA104CDE-72C0-E011-B7AC-002618943915.root'
#'/store/data/Run2011A/BTag//RECO/05Aug2011-v1/0000/FA22DF4F-8EC0-E011-B379-00261894392F.root'
#,'/store/data/Run2011A/BTag//RECO/05Aug2011-v1/0000/FA31153F-6EC0-E011-899B-002618943966.root'
#,'/store/data/Run2011A/BTag//RECO/05Aug2011-v1/0000/FA3A966A-8BC0-E011-B055-002618943905.root'
#,'/store/data/Run2011A/BTag//RECO/05Aug2011-v1/0000/FA875B20-82C0-E011-8B18-003048678FE0.root'
#,'/store/data/Run2011A/BTag//RECO/05Aug2011-v1/0000/FAAED355-76C0-E011-9C29-003048678B92.root'
#,'/store/data/Run2011A/BTag//RECO/05Aug2011-v1/0000/FC3F6C46-6EC0-E011-A2B4-00304867BEE4.root'
#,'/store/data/Run2011A/BTag//RECO/05Aug2011-v1/0000/FCB1909F-83C0-E011-AC0C-001A92811732.root'
#,'/store/data/Run2011A/BTag//RECO/05Aug2011-v1/0000/FCF10047-72C0-E011-9918-00261894390C.root'
#,'/store/data/Run2011A/BTag//RECO/05Aug2011-v1/0000/FEA84F8F-77C0-E011-966C-003048679048.root'
#,'/store/data/Run2011A/BTag//RECO/05Aug2011-v1/0000/FEF3522C-82C0-E011-BD60-003048678BB8.root'      
        
    )
)
process.load('Configuration/StandardSequences/Geometry_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'GR_R_42_V19::All' 

process.HecTau = cms.EDAnalyzer('HecTau'
   #, tracks = cms.untracked.InputTag('ctfWithMaterialTracks')
   , tracks         = cms.InputTag('generalTracks')
   , MuonsLabel     = cms.InputTag('muons')
   , MyPrint        = cms.untracked.bool(False)                  #--cout for short run
   , MyOneKaon      = cms.untracked.bool(True)                   #--run B->J/Psi + K
   , MyThreeMu      = cms.untracked.bool(False)                  #--run tau -> 3 muons
)

process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('hec_histo.root')
)

process.p = cms.Path(process.HecTau)
