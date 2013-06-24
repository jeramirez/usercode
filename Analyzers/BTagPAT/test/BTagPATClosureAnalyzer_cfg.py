# Example PAT analysis in CMSSW

import FWCore.ParameterSet.Config as cms

from Analyzers.BTagPAT.BJetOperatingPointsParameters_cfi import *

process = cms.Process("BTagPATClosureAnalyzer")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:PATLayer1_Output.fromAOD_full_ttbar.root')
#    fileNames = cms.untracked.vstring('file:PATLayer1_Output.fromAOD_full_qcd80_120_all.root')
)

process.MessageLogger = cms.Service("MessageLogger")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)
process.BTagPATClosure = cms.EDAnalyzer("BTagPATClosureAnalyzer",
    BJetOperatingPointsParameters,
    jetTag  = cms.untracked.InputTag("selectedLayer1Jets"),
    muonTag = cms.untracked.InputTag("selectedLayer1Muons"),
    BjetTag = cms.PSet(
        verbose = cms.untracked.bool(True),
        tagger = cms.string('trackCountingHighEffBJetTags'),   #away jet tagger
        purity = cms.string('Loose'),                          #away jet tagger
    ),
    jetcuts = cms.PSet(
       MinPt = cms.double(30.0),
       MaxEta = cms.double(2.4),
       MaxDeltaR = cms.double(0.7),
       MinPtRel = cms.double(-1.0) ## not ptrel cut
    ),
    muoncuts = cms.PSet(
       MinNHits = cms.int32(1),         # use default pat gives
       MinMuonPt = cms.double(6.0),
       MaxMuonEta = cms.double(2.4)     # same as jet
    )
)


process.TFileService = cms.Service("TFileService",
    fileName = cms.string('btagpatclosurepy.root')
)

process.p = cms.Path(
	process.BTagPATClosure
)

