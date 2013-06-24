import FWCore.ParameterSet.Config as cms

demo = cms.EDAnalyzer('HecPsiLambda'
    ,tracks = cms.untracked.InputTag('ctfWithMaterialTracks')
)
