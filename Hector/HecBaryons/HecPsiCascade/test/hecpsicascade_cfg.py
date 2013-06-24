import FWCore.ParameterSet.Config as cms

process = cms.Process("HecPpsiCascade")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.MagneticField_cff")

#process.load("RecoEgamma.PhotonIdentification.photonId_cff")
#process.load("Geometry.CaloEventSetup.CaloGeometry_cfi")
#process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")

#process.load("Geometry.CaloEventSetup.CaloTopology_cfi")
#process.load("RecoEgamma.EgammaPhotonProducers.photonSequence_cff")
#process.load("RecoEcal.EgammaClusterProducers.islandBasicClusters_cfi")

process.load('Configuration/StandardSequences/Geometry_cff')

#--look for the tag:  https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'GR_R_38X_V15::All' 

#process.load("Configuration/StandardSequences/RawToDigi_Data_cff")
#process.load("Configuration/StandardSequences/Reconstruction_cff")
#process.load('Configuration/EventContent/EventContent_cff')

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.MessageLogger.destinations = cms.untracked.vstring("histo_mm_05.log")

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('histo_mm_05.root')
)

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
   #'file:myfile.root'
   # The file below are at FNAL at:  /pnfs/cms/WAX/11/
    '/store/data/Run2010A/MuOnia/RECO/Nov4ReReco_v1/0000/A6FF30E2-99E9-DF11-9C75-90E6BA0D0990.root'
   ,'/store/data/Run2010A/MuOnia/RECO/Nov4ReReco_v1/0000/DAA02AFF-49E9-DF11-B297-90E6BAE8CC08.root'
   ,'/store/data/Run2010A/MuOnia/RECO/Nov4ReReco_v1/0000/E0359986-36E9-DF11-BC78-E0CB4E29C4CB.root'
   ,'/store/data/Run2010A/MuOnia/RECO/Nov4ReReco_v1/0000/0CA891C2-45E9-DF11-887C-485B39800C10.root'
   ,'/store/data/Run2010A/MuOnia/RECO/Nov4ReReco_v1/0000/20E03B56-39E9-DF11-BE17-485B39800B96.root'
   ,'/store/data/Run2010A/MuOnia/RECO/Nov4ReReco_v1/0000/7ECE5F11-3BE9-DF11-B4C4-E0CB4EA0A8D7.root'
   ,'/store/data/Run2010A/MuOnia/RECO/Nov4ReReco_v1/0000/A23DAC34-97E9-DF11-B416-90E6BA0D09E6.root'
   ,'/store/data/Run2010A/MuOnia/RECO/Nov4ReReco_v1/0000/B0BB0331-ADE9-DF11-8D7B-E0CB4E553636.root'
   ,'/store/data/Run2010A/MuOnia/RECO/Nov4ReReco_v1/0000/FC8EADE9-35E9-DF11-80B3-E0CB4E19F969.root'
   ,'/store/data/Run2010A/MuOnia/RECO/Nov4ReReco_v1/0000/F05C934F-8CE9-DF11-B070-E0CB4E19F97C.root'
    )
)

process.HecMuons = cms.EDAnalyzer('HecPsiCascade'
  #, tracks = cms.untracked.InputTag('ctfWithMaterialTracks')
   , tracks         = cms.InputTag('generalTracks')
   , gsfElectrons   = cms.InputTag('gsfElectrons')
   , MuonsLabel     = cms.InputTag('muons')
   , CaloMuonsLabel = cms.InputTag("calomuons")
   , minTracks      = cms.untracked.uint32(0)
   , check_elC      = cms.untracked.bool(True)                   #--Dump e- matching info
   , Pe_match       = cms.untracked.double(10.0)
   , Mass_Constrain = cms.string("m_jpsi")                       #--Mass constrain fit
   , onlyDiMu       = cms.untracked.bool(True)                   #--run only Dimuons or DiMu + X
   , MyPrint        = cms.untracked.bool(False)                   #--cout for short run
   , CasAlgo        = cms.untracked.string("MyCascade")          #--Label of Cascade collection (My clone)
   , CasDecayName   = cms.untracked.string("Cascade")            #--Collection Label decay (cascade or omega)(Eduardo}
   , VeeAlgo        = cms.untracked.string("looseV0Candidates")  #--Collection Label decay (Lambda) (Eduardo}
)

process.dump=cms.EDAnalyzer('EventContentAnalyzer')

#--Cascade
process.load('Analyzers.CascadeProducer.looseV0Candidates_cff')
process.load('Analyzers.CascadeProducer.cascadeproducer_cff')
process.MyCascade = process.looseCascade.clone()
process.MyCascade.XiDCACut    = cms.double(2.0)
process.MyCascade.VtxFromFit  = cms.bool(True)
process.MyCascade.MassFromFit = cms.bool(True)
#process.MyCascade.v0Algo      = cms.string('generalV0Candidates')   #--cms default V0  L/S > 15
#process.MyCascade.v0Algo      = cms.string('looseV0Candidates')     #--cms V0 with L/S > 5
process.LambdaFilter = cms.EDFilter("VeeCountFilter",
                                                      algo = cms.string('looseV0Candidates'),
                                                      name = cms.string('Lambda'),
                                                      minNumber = cms.uint32(1) 
                                                      )

#process.p = cms.Path(process.HecMuons*process.dump)
#process.p = cms.Path(process.looseV0Candidates*process.LambdaFilter*process.MyCascade*process.HecMuons)
process.p = cms.Path(process.HecMuons)
