import os
import sys

import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("Demo")

options   = VarParsing.VarParsing('standard')
#set default
options.maxEvents = -1
options.files     = 'file:cascproducer.root'
options.output    = 'mycasntplelooseskim2010.root'
options.parseArguments()

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Services_cff") 
process.load("Configuration.StandardSequences.Geometry_cff")
#process.load("Geometry.TrackerRecoData.trackerRecoGeometryXML_cfi")
#process.load("Geometry.TrackerNumberingBuilder.trackerNumberingGeometry_cfi")
#process.load("Geometry.TrackerGeometryBuilder.trackerGeometry_cfi")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

readFiles = cms.untracked.vstring(options.files)
outPutFile= cms.string(options.output)
process.source = cms.Source("PoolSource",fileNames = readFiles)

#readFiles.extend((
#        'file:cascproducerlooseveeskimmall2010_0179_09.root',
#        'file:cascproducerlooseveeskimmall2010_0179_10.root',
#        'file:cascproducerlooseveeskimmall2010_0179_11.root',
#        'file:cascproducerlooseveeskimmall2010_0179_12.root',
#        'file:cascproducerlooseveeskimmall2010_0179_13.root',
#        'file:cascproducerlooseveeskimmall2010_0179_14.root',
#        'file:cascproducerlooseveeskimmall2010_0179_15.root',
#        'file:cascproducerlooseveeskimmall2010_0179_16.root',
#))

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = cms.string('GR_R_35X_V8::All')
process.GlobalTag.globaltag = 'GR_R_39X_V5::All'
process.load("Configuration.StandardSequences.MagneticField_cff")

process.dca02 = cms.EDAnalyzer('CascadeNtuplizer',
    VeeAlgo       = cms.string('looseV0Candidates'),
    CasAlgo       = cms.string('CascCandidatesDCA02'),
    CasDecayName  = cms.untracked.string('Cascade'),
    trackingAlgo  = cms.InputTag('generalTracks')
)

process.dca05 = cms.EDAnalyzer('CascadeNtuplizer',
    VeeAlgo       = cms.string('looseV0Candidates'),
    CasAlgo       = cms.string('CascCandidatesDCA05'),
    CasDecayName  = cms.untracked.string('Cascade'),
    trackingAlgo  = cms.InputTag('generalTracks')
)

process.dca1 = cms.EDAnalyzer('CascadeNtuplizer',
    VeeAlgo       = cms.string('looseV0Candidates'),
    CasAlgo       = cms.string('CascCandidatesDCA1'),
    CasDecayName  = cms.untracked.string('Cascade'),
    trackingAlgo  = cms.InputTag('generalTracks')
)

process.dca2 = cms.EDAnalyzer('CascadeNtuplizer',
    VeeAlgo       = cms.string('looseV0Candidates'),
    CasAlgo       = cms.string('CascCandidatesDCA2'),
    CasDecayName  = cms.untracked.string('Cascade'),
    trackingAlgo  = cms.InputTag('generalTracks')
)

process.dca5 = cms.EDAnalyzer('CascadeNtuplizer',
    VeeAlgo       = cms.string('looseV0Candidates'),
    CasAlgo       = cms.string('CascCandidatesDCA5'),
    CasDecayName  = cms.untracked.string('Cascade'),
    trackingAlgo  = cms.InputTag('generalTracks')
)

process.dca1fitvtx = cms.EDAnalyzer('CascadeNtuplizer',
    VeeAlgo       = cms.string('looseV0Candidates'),
    CasAlgo       = cms.string('CascCandidatesDCA1fitvtx'),
    CasDecayName  = cms.untracked.string('Cascade'),
    trackingAlgo  = cms.InputTag('generalTracks')
)
process.dca2fitvtx = cms.EDAnalyzer('CascadeNtuplizer',
    VeeAlgo       = cms.string('looseV0Candidates'),
    CasAlgo       = cms.string('CascCandidatesDCA2fitvtx'),
    CasDecayName  = cms.untracked.string('Cascade'),
    trackingAlgo  = cms.InputTag('generalTracks')
)
process.dca5fitvtx = cms.EDAnalyzer('CascadeNtuplizer',
    VeeAlgo       = cms.string('looseV0Candidates'),
    CasAlgo       = cms.string('CascCandidatesDCA5fitvtx'),
    CasDecayName  = cms.untracked.string('Cascade'),
    trackingAlgo  = cms.InputTag('generalTracks')
)
process.dca02fit2 = cms.EDAnalyzer('CascadeNtuplizer',
    VeeAlgo       = cms.string('looseV0Candidates'),
    CasAlgo       = cms.string('CascCandidatesDCA02fit2'),
    CasDecayName  = cms.untracked.string('Cascade'),
    trackingAlgo  = cms.InputTag('generalTracks')
)
process.dca1fit2 = cms.EDAnalyzer('CascadeNtuplizer',
    VeeAlgo       = cms.string('looseV0Candidates'),
    CasAlgo       = cms.string('CascCandidatesDCA1fit2'),
    CasDecayName  = cms.untracked.string('Cascade'),
    trackingAlgo  = cms.InputTag('generalTracks')
)
process.dca2fit2 = cms.EDAnalyzer('CascadeNtuplizer',
    VeeAlgo       = cms.string('looseV0Candidates'),
    CasAlgo       = cms.string('CascCandidatesDCA2fit2'),
    CasDecayName  = cms.untracked.string('Cascade'),
    trackingAlgo  = cms.InputTag('generalTracks')
)
process.dca5fit2 = cms.EDAnalyzer('CascadeNtuplizer',
    VeeAlgo       = cms.string('looseV0Candidates'),
    CasAlgo       = cms.string('CascCandidatesDCA5fit2'),
    CasDecayName  = cms.untracked.string('Cascade'),
    trackingAlgo  = cms.InputTag('generalTracks')
)
process.dca10fit2 = cms.EDAnalyzer('CascadeNtuplizer',
    VeeAlgo       = cms.string('looseV0Candidates'),
    CasAlgo       = cms.string('CascCandidatesDCA10fit2'),
    CasDecayName  = cms.untracked.string('Cascade'),
    trackingAlgo  = cms.InputTag('generalTracks')
)

process.odca2 = cms.EDAnalyzer('CascadeNtuplizer',
    xmin          = cms.untracked.double(1.62),
    xmax          = cms.untracked.double(1.72),
    VeeAlgo       = cms.string('looseV0Candidates'),
    CasAlgo       = cms.string('CascCandidatesDCA2'),
    CasDecayName  = cms.untracked.string('Omega'),
    trackingAlgo  = cms.InputTag('generalTracks')
)

process.odca5 = cms.EDAnalyzer('CascadeNtuplizer',
    xmin          = cms.untracked.double(1.62),
    xmax          = cms.untracked.double(1.72),
    VeeAlgo       = cms.string('looseV0Candidates'),
    CasAlgo       = cms.string('CascCandidatesDCA5'),
    CasDecayName  = cms.untracked.string('Omega'),
    trackingAlgo  = cms.InputTag('generalTracks')
)

process.odca2fit2 = cms.EDAnalyzer('CascadeNtuplizer',
    xmin          = cms.untracked.double(1.62),
    xmax          = cms.untracked.double(1.72),
    pitrk_pt_cut  = cms.untracked.double(1.0),
    VeeAlgo       = cms.string('looseV0Candidates'),
    CasAlgo       = cms.string('CascCandidatesDCA2fit2'),
    CasDecayName  = cms.untracked.string('Omega'),
    trackingAlgo  = cms.InputTag('generalTracks')
)

process.odca5fit2 = cms.EDAnalyzer('CascadeAnalyzer',
    xmin          = cms.untracked.double(1.62),
    xmax          = cms.untracked.double(1.72),
    pitrk_pt_cut  = cms.untracked.double(1.0), 
    VeeAlgo       = cms.string('looseV0Candidates'),
    CasAlgo       = cms.string('CascCandidatesDCA5fit2'),
    CasDecayName  = cms.untracked.string('Omega'),
    trackingAlgo  = cms.InputTag('generalTracks')
)

##################################bscnobamhalo############################################
process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff')
process.load('HLTrigger/HLTfilters/hltLevel1GTSeed_cfi')

process.L1T1=process.hltLevel1GTSeed.clone()
process.L1T1.L1TechTriggerSeeding = cms.bool(True)
#process.L1T1.L1SeedsLogicalExpression = cms.string('(40 OR 41) AND NOT (36 OR 37 OR 38 OR 39)')
process.L1T1.L1SeedsLogicalExpression = cms.string('0 AND (40 OR 41) AND NOT (36 OR 37 OR 38 OR 39)')

####################apply the scraping event filter here#####################
process.noScraping= cms.EDFilter("FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False), ## Or 'True' to get some per-event info
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.2)
)
process.noscrapingnobeamhalo = cms.Path(process.L1T1*process.noScraping)

process.TFileService = cms.Service("TFileService",
    fileName = outPutFile
#    fileName = cms.string('mycasanalyzerlooseskim2010.root'),
#    SelectEvents = cms.untracked.PSet(
#        SelectEvents = cms.vstring('noscrapingnobeamhalo')
#    )
)

process.p = cms.Path(
process.dca2 + process.dca2fit2 + process.dca2fitvtx
+ process.dca1 + process.dca1fit2 + process.dca1fitvtx
+ process.dca5 +process.dca5fitvtx + process.dca5fit2
+ process.odca2 + process.odca5
+ process.odca2fit2 + process.odca5fit2
#+ process.dca10fit2
+ process.dca02 + process.dca02fit2
#+ process.dca05
)
