import FWCore.ParameterSet.Config as cms

process = cms.Process("MyHecPsiLambda")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.MagneticField_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
       #--'file:myfile.root'
 '/store/data/Run2010A/MuOnia/RECO/Apr21ReReco-v1/0000/261CA6C4-526E-E011-B236-E0CB4E29C504.root'
,'/store/data/Run2010A/MuOnia/RECO/Apr21ReReco-v1/0000/3C00ACA7-476E-E011-B98D-E0CB4E19F9A3.root'
,'/store/data/Run2010A/MuOnia/RECO/Apr21ReReco-v1/0000/82439C9B-DF6D-E011-931A-001EC9D8D085.root'
,'/store/data/Run2010A/MuOnia/RECO/Apr21ReReco-v1/0000/54518D04-E26D-E011-8579-485B39800BBB.root'
,'/store/data/Run2010A/MuOnia/RECO/Apr21ReReco-v1/0000/E0FE9A5C-CA6D-E011-AC5A-E0CB4E19F959.root'
,'/store/data/Run2010A/MuOnia/RECO/Apr21ReReco-v1/0000/D67426EE-CD6D-E011-9605-90E6BA19A226.root'
,'/store/data/Run2010A/MuOnia/RECO/Apr21ReReco-v1/0000/848D7942-D36D-E011-82D4-90E6BA0D09B9.root'
,'/store/data/Run2010A/MuOnia/RECO/Apr21ReReco-v1/0000/3E3D73E5-C76D-E011-921E-E0CB4E1A116D.root'
,'/store/data/Run2010A/MuOnia/RECO/Apr21ReReco-v1/0000/7C890594-C56D-E011-9774-0019BB3FF40C.root'
,'/store/data/Run2010A/MuOnia/RECO/Apr21ReReco-v1/0000/1AC00621-C76D-E011-AFD5-90E6BA442EEC.root'
,'/store/data/Run2010A/MuOnia/RECO/Apr21ReReco-v1/0000/808995AE-DF6D-E011-B04F-90E6BA0D09B4.root'
,'/store/data/Run2010A/MuOnia/RECO/Apr21ReReco-v1/0000/6464FBE1-D26D-E011-AEE6-90E6BA442F36.root'
,'/store/data/Run2010A/MuOnia/RECO/Apr21ReReco-v1/0000/C0BC15B5-EF6D-E011-80F3-00221992FF06.root'
,'/store/data/Run2010A/MuOnia/RECO/Apr21ReReco-v1/0000/6A9D842B-C66D-E011-ACA5-E0CB4E1A1182.root'
,'/store/data/Run2010A/MuOnia/RECO/Apr21ReReco-v1/0000/9EB1047D-D46D-E011-A19E-E0CB4E1A119A.root'
,'/store/data/Run2010A/MuOnia/RECO/Apr21ReReco-v1/0000/7C9D3747-E06D-E011-96BB-E0CB4E29C4FD.root'
    )
)

process.load('Configuration/StandardSequences/Geometry_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'GR_R_42_V19::All' 

process.HecJL = cms.EDAnalyzer('HecPsiLambda'
   , tracks     = cms.InputTag('generalTracks')
   , MuonsLabel = cms.InputTag('muons')
   , MyPrint    = cms.untracked.bool(False)                  #--cout for short run
)

process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('histo_jlambda.root')
)

process.p = cms.Path(process.HecJL)
