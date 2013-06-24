import FWCore.ParameterSet.Config as cms

selectedLayer1JetsPt30 = cms.EDFilter("PATJetSelector",
    src = cms.InputTag("selectedLayer1Jets"),
    cut = cms.string('pt > 30. && abs(eta) < 2.4')
)


