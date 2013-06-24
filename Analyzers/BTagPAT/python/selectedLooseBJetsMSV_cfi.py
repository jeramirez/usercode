import FWCore.ParameterSet.Config as cms

from Analyzers.BTagPAT.Jetpt30Selector_cfi import *
selectedLooseBJetsMSV = cms.EDFilter("BDiscriminatorPatJetSelector",
    src = cms.InputTag("selectedLayer1JetsPt30"),
    disc = cms.string('combinedSecondaryVertexMVABJetTags'),
    discCut = cms.double(0.4291)
)


