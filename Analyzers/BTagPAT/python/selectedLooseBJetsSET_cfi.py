import FWCore.ParameterSet.Config as cms

from Analyzers.BTagPAT.Jetpt30Selector_cfi import *
selectedLooseBJetsSET = cms.EDFilter("BDiscriminatorPatJetSelector",
    src = cms.InputTag("selectedLayer1JetsPt30"),
    disc = cms.string('softElectronBJetTags'),
    discCut = cms.double(0.0)
)


