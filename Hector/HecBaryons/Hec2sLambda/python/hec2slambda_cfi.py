import FWCore.ParameterSet.Config as cms

demo = cms.EDAnalyzer('Hec2sLambda'
    ,tracks = cms.untracked.InputTag('ctfWithMaterialTracks')
)
