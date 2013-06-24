# Example PAT analysis in CMSSW
# Used for BTagPATAnalyzer
#cmsRun BTagPATAnalyzer2009_cfg.py print maxEvents=50 output2=test.root patlayer1=allLayer1Jets

import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

# setup 'standard'  options
options = VarParsing.VarParsing ('standard')
options.register ('output2',
                  'btagpatanalyzerpy.root',                        # default value
                  VarParsing.VarParsing.multiplicity.singleton,  # singleton or list
                  VarParsing.VarParsing.varType.string,          # string, int, or float
                  "Outputfile name I want to")
options.register ('patlayer1',
                  'selectedLayer1Jets',                          # default value
                  VarParsing.VarParsing.multiplicity.singleton,  # singleton or list
                  VarParsing.VarParsing.varType.string,          # string, int, or float
                  "Pat Layer1 Container")

# setup any defaults you want
options.output = 'standard_default_try_3.root'
options.files= 'file:PATLayer1_Output.fromAOD_full_ttbar.root'
options.maxEvents = -1 # -1 means all events
# get and parse the command line arguments
options.parseArguments()

from Analyzers.BTagPAT.BJetOperatingPointsParameters_cfi import *

process = cms.Process("BTagPATAnalyzer")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.files)
)

process.MessageLogger = cms.Service("MessageLogger")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)
process.BTagPATAnalyzerTC2 = cms.EDAnalyzer("BTagPATAnalyzer",
    BJetOperatingPointsParameters,
    jetTag = cms.untracked.InputTag(options.patlayer1),
    BjetTag = cms.PSet(
        verbose = cms.untracked.bool(True),
        tagger = cms.untracked.string('TC2'),
        purity = cms.string('Loose'),
        discriminator = cms.string('trackCountingHighEffBJetTags'),
        maxdiscriminatorcut = cms.untracked.double(30.0),
        mindiscriminatorcut = cms.untracked.double(-10.0)
    )
)


process.BTagPATAnalyzerTC3 = cms.EDAnalyzer("BTagPATAnalyzer",
    BJetOperatingPointsParameters,
    jetTag = cms.untracked.InputTag(options.patlayer1),
    BjetTag = cms.PSet(
        verbose = cms.untracked.bool(False),
        tagger = cms.untracked.string('TC3'),
        purity = cms.string('Loose'),
        discriminator = cms.string('trackCountingHighPurBJetTags'),
        maxdiscriminatorcut = cms.untracked.double(30.0),
        mindiscriminatorcut = cms.untracked.double(-10.0)
    )
)

process.BTagPATAnalyzerTP = cms.EDAnalyzer("BTagPATAnalyzer",
    BJetOperatingPointsParameters,
    jetTag = cms.untracked.InputTag(options.patlayer1),
    BjetTag = cms.PSet(
        verbose = cms.untracked.bool(False),
        tagger = cms.untracked.string('TP'),
        purity = cms.string('Loose'),
        discriminator = cms.string('jetProbabilityBJetTags'),
        maxdiscriminatorcut = cms.untracked.double(2.6),
        mindiscriminatorcut = cms.untracked.double(-0.1)
    )
)

process.BTagPATAnalyzerBTP = cms.EDAnalyzer("BTagPATAnalyzer",
    BJetOperatingPointsParameters,
    jetTag = cms.untracked.InputTag(options.patlayer1),
    BjetTag = cms.PSet(
        verbose = cms.untracked.bool(False),
        tagger = cms.untracked.string('BTP'),
        purity = cms.string('Loose'),
        discriminator = cms.string('jetBProbabilityBJetTags'),
        maxdiscriminatorcut = cms.untracked.double(8.1),
        mindiscriminatorcut = cms.untracked.double(-0.1)
    )
)
process.BTagPATAnalyzerSSV = cms.EDAnalyzer("BTagPATAnalyzer",
    BJetOperatingPointsParameters,
    jetTag = cms.untracked.InputTag(options.patlayer1),
    BjetTag = cms.PSet(
        verbose = cms.untracked.bool(False),
        tagger = cms.untracked.string('SSV'),
        purity = cms.string('Loose'),
        discriminator = cms.string('simpleSecondaryVertexBJetTags'),
        maxdiscriminatorcut = cms.untracked.double(8.0),
        mindiscriminatorcut = cms.untracked.double(0.0)
    )
)

process.BTagPATAnalyzerCSV = cms.EDAnalyzer("BTagPATAnalyzer",
    BJetOperatingPointsParameters,
    jetTag = cms.untracked.InputTag(options.patlayer1),
    BjetTag = cms.PSet(
        verbose = cms.untracked.bool(False),
        tagger = cms.untracked.string('CSV'),
        purity = cms.string('Loose'),
        discriminator = cms.string('combinedSecondaryVertexBJetTags'),
        maxdiscriminatorcut = cms.untracked.double(1.1),
        mindiscriminatorcut = cms.untracked.double(-0.1)
    )
)

process.BTagPATAnalyzerMSV = cms.EDAnalyzer("BTagPATAnalyzer",
    BJetOperatingPointsParameters,
    jetTag = cms.untracked.InputTag(options.patlayer1),
    BjetTag = cms.PSet(
        verbose = cms.untracked.bool(False),
        tagger = cms.untracked.string('MSV'),
        purity = cms.string('Loose'),
        discriminator = cms.string('combinedSecondaryVertexMVABJetTags'),
        maxdiscriminatorcut = cms.untracked.double(1.1),
        mindiscriminatorcut = cms.untracked.double(-0.1)
    )
)

process.BTagPATAnalyzerSETbyIP3d = cms.EDAnalyzer("BTagPATAnalyzer",
    BJetOperatingPointsParameters,
    jetTag = cms.untracked.InputTag(options.patlayer1),
    BjetTag = cms.PSet(
        verbose = cms.untracked.bool(False),
        tagger = cms.untracked.string('SETbyIP3d'),
        purity = cms.string('Loose'),
        discriminator = cms.string('softElectronByIP3dBJetTags'),
        maxdiscriminatorcut = cms.untracked.double(30.0),
        mindiscriminatorcut = cms.untracked.double(-10.0)
    )
)

process.BTagPATAnalyzerSETbyPt = cms.EDAnalyzer("BTagPATAnalyzer",
    BJetOperatingPointsParameters,
    jetTag = cms.untracked.InputTag(options.patlayer1),
    BjetTag = cms.PSet(
        verbose = cms.untracked.bool(False),
        tagger = cms.untracked.string('SETbyPt'),
        purity = cms.string('Loose'),
        discriminator = cms.string('softElectronByPtBJetTags'),
        maxdiscriminatorcut = cms.untracked.double(8.01),
        mindiscriminatorcut = cms.untracked.double(-0.01)
    )
)

process.BTagPATAnalyzerSMT = cms.EDAnalyzer("BTagPATAnalyzer",
    BJetOperatingPointsParameters,
    jetTag = cms.untracked.InputTag(options.patlayer1),
    BjetTag = cms.PSet(
        verbose = cms.untracked.bool(False),
        tagger = cms.untracked.string('SMT'),
        purity = cms.string('Loose'),
        discriminator = cms.string('softMuonBJetTags'),
        maxdiscriminatorcut = cms.untracked.double(1.1),
        mindiscriminatorcut = cms.untracked.double(-0.1)
    )
)

process.BTagPATAnalyzerSMTbyIP3d = cms.EDAnalyzer("BTagPATAnalyzer",
    BJetOperatingPointsParameters,
    jetTag = cms.untracked.InputTag(options.patlayer1),
    BjetTag = cms.PSet(
        verbose = cms.untracked.bool(False),
        tagger = cms.untracked.string('SMTbyIP3d'),
        purity = cms.string('Loose'),
        discriminator = cms.string('softMuonByIP3dBJetTags'),
        maxdiscriminatorcut = cms.untracked.double(30.0),
        mindiscriminatorcut = cms.untracked.double(-10.0)
    )
)

process.BTagPATAnalyzerSMTbyPt = cms.EDAnalyzer("BTagPATAnalyzer",
    BJetOperatingPointsParameters,
    jetTag = cms.untracked.InputTag(options.patlayer1),
    BjetTag = cms.PSet(
        verbose = cms.untracked.bool(False),
        tagger = cms.untracked.string('SMTbyPt'),
        purity = cms.string('Loose'),
        discriminator = cms.string('softMuonByPtBJetTags'),
        maxdiscriminatorcut = cms.untracked.double(8.01),
        mindiscriminatorcut = cms.untracked.double(-0.01)
    )
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(options.output2)
)

process.p = cms.Path(
	process.BTagPATAnalyzerTC2
	*process.BTagPATAnalyzerTC3
	*process.BTagPATAnalyzerBTP
	*process.BTagPATAnalyzerSSV
	*process.BTagPATAnalyzerCSV
	*process.BTagPATAnalyzerMSV
	*process.BTagPATAnalyzerSETbyIP3d
	*process.BTagPATAnalyzerSETbyPt
	*process.BTagPATAnalyzerSMT
	*process.BTagPATAnalyzerSMTbyIP3d
	*process.BTagPATAnalyzerSMTbyPt
	*process.BTagPATAnalyzerTP
)

