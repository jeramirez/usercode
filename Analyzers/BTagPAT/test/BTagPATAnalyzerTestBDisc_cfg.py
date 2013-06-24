# Example PAT analysis in CMSSW
# This example shows how the filter module works (selectedLayer1BJetsPt30)
#  Because the parser does not works for the discriminant
#  we use the filtered data 

import FWCore.ParameterSet.Config as cms
from Analyzers.BTagPAT.BJetOperatingPointsParameters_cfi import *

process = cms.Process("BTagPATAnalyzer")
#uncomment this if want select Loose BJets
#process.load("Analyzers.BTagPAT.selectedLooseBJets_cff")
#
# and add into the path any of the following lines
#process.LooseBJetsTC2     <----
#process.LooseBJetsTC3         |
#process.LooseBJetsTP          |
#process.LooseBJetsSSV         |
#process.LooseBJetsCSV         |
#process.LooseBJetsMSV         |
#process.LooseBJetsSET         |
#process.LooseBJetsSMT         |
#process.selectedLooseBJets <----# If don't know which one select this

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:PATLayer1_Output.fromAOD_full_ttbar.root')
)

#keep the logging output to a nice level
process.MessageLogger = cms.Service("MessageLogger")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

process.BTagPATAnalyzerTC2 = cms.EDAnalyzer("BTagPATAnalyzer",
    BJetOperatingPointsParameters,
# Baseline cut to compute operating points
#    jetTag = cms.untracked.InputTag("selectedLayer1BJetsPt30"),
# This uses the selected b- sample for a given analysis
    jetTag = cms.untracked.InputTag("selectedMyLooseBJetsTC2"),
    BjetTag = cms.PSet(
        verbose = cms.untracked.bool(True),
        tagger = cms.untracked.string('TC2'),
        purity = cms.string('Loose'),
        discriminator = cms.string('trackCountingHighEffBJetTags'),
        maxdiscriminatorcut = cms.untracked.double(15.0),
        mindiscriminatorcut = cms.untracked.double(-1.0)
    )
)

####### The two filters below can be loaded directly when 
####### uncommeting  selectedLooseBJets_cff and adding the selection 
####### into the path and analyzer
#filter to select b-jets 
#see Jetpt30Selector_cfi.py
process.selectedLayer1BJetsPt30 = cms.EDFilter("PATJetSelector",
    src = cms.InputTag("selectedLayer1Jets"),
#  Parser works here
     cut = cms.string('pt > 30. && abs(eta)< 2.4')
#  Parser fail in any of cuts below
#    cut = cms.string('pt > 30. & bDiscriminator( trackCountingHighEffBJetTags ) > 2')
#    cut = cms.string('pt > 30. & jetCharge > 0')
#    cut = cms.string('pt > 30. & bDiscriminator("") > 2')
#    cut = cms.string('pt > 30. & bDiscriminator("trackCountingHighEffBJetTags") > 2')
#    cut = cms.string('pt > 30. & bDiscriminator() > 2')
)

# extra filter to add the discriminant in the selection
# this filter maybe discarted once the parser is fixed.
# see selectedLooseBJetsTC2_cfi.py and comment above.
process.selectedMyLooseBJetsTC2 = cms.EDFilter("BDiscriminatorPatJetSelector",
    src = cms.InputTag("selectedLayer1BJetsPt30"),
    disc = cms.string('trackCountingHighEffBJetTags'),
    discCut = cms.double(1.912)
)

# this is where the output (histogram) file is defined:
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('btagpatanalyzerpy_bdisc.root')
)

process.p = cms.Path(
	process.selectedLayer1BJetsPt30
	*process.selectedMyLooseBJetsTC2
	*process.BTagPATAnalyzerTC2
)

