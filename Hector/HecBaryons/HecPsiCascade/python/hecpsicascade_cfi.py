import FWCore.ParameterSet.Config as cms

demo = cms.EDAnalyzer('HecPsiCascade'
    ,tracks = cms.untracked.InputTag('ctfWithMaterialTracks')
)
