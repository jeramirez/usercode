import FWCore.ParameterSet.Config as cms

from Analyzers.BTagPAT.Jetpt30Selector_cfi import *
selectedLooseBJetsSSV = cms.EDFilter("BDiscriminatorPatJetSelector",
    src = cms.InputTag("selectedLayer1JetsPt30"),
    disc = cms.string('simpleSecondaryVertexBJetTags'),
    discCut = cms.double(1.2)
)


