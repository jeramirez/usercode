import FWCore.ParameterSet.Config as cms

demo = cms.EDAnalyzer('HecTau'
    ,tracks = cms.untracked.InputTag('ctfWithMaterialTracks')
)
