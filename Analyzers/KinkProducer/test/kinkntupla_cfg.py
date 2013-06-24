import os
import sys

import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
#from Analyzers.CascadeProducer.Cert_136033_149442_7TeV_Dec22ReReco_Collisions10_JSON_v4_cfi import *

#a copy of T3Ximumuskim7TeVsoft2011_cfg.py
process = cms.Process("DEBUG")

options   = VarParsing.VarParsing('standard')
#set default
options.maxEvents = -1
#options.files     = '/user/jramirez/cascade/424/run2010A_locatest_with_fallback/cascadeskimRun2010A_368_2_8nP.root'
options.output    = 'cascadeskimRun2011A.root'
#--CRAB--don't like this line --#
options.parseArguments()
readFiles = cms.untracked.vstring(options.files)
#process.out is untracked
#outPutFile= cms.untracked.string(options.output)
#process TFileService is tracked
outPutFile= cms.string(options.output)

# import of standard configurations
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Services_cff") 
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.EventContent.EventContent_cff')
process.load('Geometry.TrackerGeometryBuilder.trackerGeometry_cfi')
process.TrackerDigiGeometryESModule.applyAlignment = False
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
# import cascade standard configurations
process.load('Analyzers.CascadeProducer.looseV0Candidates_cff')
process.load('Analyzers.CascadeProducer.cascadeproducer_cff')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )
process.options = cms.untracked.PSet(
     SkipEvent = cms.untracked.vstring('ProductNotFound')
)

# Input source
#  Lumi info
#-------------------------------------------
#Json file:  Cert_136033-149442_7TeV_Dec22ReReco_Collisions10_JSON_v4.txt  written.
#------------------------------------------- 
process.source = cms.Source("PoolSource",fileNames = readFiles)
#process.source.lumisToProcess = lumisToProcess7TeVCollisions10

#skim output
##output module definition
#process.RECOEventContent.outputCommands.append('keep *_looseV0Candidates_*_*')
##process.RECOEventContent.outputCommands.append('keep *_CascCandidates_*_*')
##process.RECOEventContent.outputCommands.append('keep *_looseCascade_*_*')
#process.RECOEventContent.outputCommands.append('keep *_looseXiPixelTracks_*_*')
#process.RECOEventContent.outputCommands.append('keep *_looseOmegaPixelTracks_*_*')
#process.RECOEventContent.outputCommands.append('keep *_XiPixelTracksCand_*_*')
#process.RECOEventContent.outputCommands.append('keep *_OmegaPixelTracksCand_*_*')
#process.out = cms.OutputModule("PoolOutputModule",
#    splitLevel = cms.untracked.int32(0),
#    outputCommands = process.RECOEventContent.outputCommands,
#    fileName = outPutFile,
#    SelectEvents = cms.untracked.PSet(
#        SelectEvents = cms.vstring('casevts', 'casevtskf', 'omeevts', 'omeevtskf')
#    )
#)

############# ntupla output
process.TFileService = cms.Service("TFileService",
    fileName = outPutFile
#    fileName = cms.string('mycasanalyzerlooseskim2010.root'),
#    SelectEvents = cms.untracked.PSet(
#        SelectEvents = cms.vstring('noscrapingnobeamhalo')
#    )
)
##############################

#process.GlobalTag.globaltag = cms.string('GR_R_35X_V8::All')
#process.GlobalTag.globaltag = 'GR_R_39X_V5::All'
process.GlobalTag.globaltag = 'GR_R_42_V19::All'

#Lambda Loose Skim Filter
process.LambdaFilter = cms.EDFilter("VeeCountFilter",
                                                      algo = cms.string('looseV0Candidates'),
                                                      name = cms.string('Lambda'),
                                                      minNumber = cms.uint32(1) 
                                                      )
#Cascade Producer
# Default cascade cuts 
#   DCA = 2cm, 
#   use Vertex from kinematic fitter
#   use Mass from kinematic fitter  
process.CascCandidates = process.looseCascade.clone()
process.CascCandidates.XiDCACut           = cms.double(2.0)
process.CascCandidates.VtxFromFit         = cms.bool(True)
process.CascCandidates.MassFromFit        = cms.bool(True)

#process.looseXiPixelTracks = cms.EDProducer("CascadeProducerPixellessCandidates",
process.looseXiPixelTracks = cms.EDProducer("CascadeCandidatesProducer",
    CasAlgo       = cms.string('looseCascade'),
    CasDecayName  = cms.untracked.string('Cascade'),
    trackingAlgo  = cms.InputTag('generalTracks'),
    XiImpactParameterAtPrimary = cms.double(50.),
    XiCLAtDecayforNonTracks    = cms.double(0.0005),
    XiSignificanceSeparation3D = cms.double(3.),
    WithPixelTracksCandidates  = cms.bool(True)
)
#process.XiPixelTracksCand = cms.EDProducer("CascadeProducerPixellessCandidates",
process.XiPixelTracksCand = cms.EDProducer("CascadeCandidatesProducer",
    CasAlgo       = cms.string('CascCandidates'),
    CasDecayName  = cms.untracked.string('Cascade'),
    trackingAlgo  = cms.InputTag('generalTracks'),
    XiImpactParameterAtPrimary = cms.double(50.),
    XiCLAtDecayforNonTracks    = cms.double(0.0005),
    XiSignificanceSeparation3D = cms.double(3.),
    WithPixelTracksCandidates  = cms.bool(True)
)
#process.looseOmegaPixelTracks = cms.EDProducer("CascadeProducerPixellessCandidates",
process.looseOmegaPixelTracks = cms.EDProducer("CascadeCandidatesProducer",
    CasAlgo       = cms.string('looseCascade'),
    CasDecayName  = cms.untracked.string('Omega'),
    trackingAlgo  = cms.InputTag('generalTracks'),
    XiImpactParameterAtPrimary = cms.double(50.),
    XiCLAtDecayforNonTracks    = cms.double(0.0005),
    XiSignificanceSeparation3D = cms.double(3.),
    WithPixelTracksCandidates  = cms.bool(True)
)
#process.OmegaPixelTracksCand = cms.EDProducer("CascadeProducerPixellessCandidates",
process.OmegaPixelTracksCand = cms.EDProducer("CascadeCandidatesProducer",
    CasAlgo       = cms.string('CascCandidates'),
    CasDecayName  = cms.untracked.string('Omega'),
    trackingAlgo  = cms.InputTag('generalTracks'),
    XiImpactParameterAtPrimary = cms.double(50.),
    XiCLAtDecayforNonTracks    = cms.double(0.0005),
    XiSignificanceSeparation3D = cms.double(3.),
    WithPixelTracksCandidates  = cms.bool(True)
)
process.cascadeproducer = cms.Path(process.looseCascade + process.CascCandidates)

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
    thresh = cms.untracked.double(0.25)
)

####################apply good primary vertex filter here#####################
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                           vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                                           minimumNDOF = cms.uint32(4),
                                           maxAbsZ = cms.double(15), 
                                           maxd0 = cms.double(2) 
)
##################################cas filter #################################
process.CasFilter = process.looseCasFilter.clone()
process.CasFilter.algo = cms.string('CascCandidates')
#process.CasFilter.algo = cms.string('XiPixelTracksCand')
process.looseCasFilter.algo = cms.string('looseXiPixelTracks')

process.OmeFilter = process.looseOmeFilter.clone()
process.OmeFilter.algo = cms.string('CascCandidates')
#process.OmeFilter.algo = cms.string('OmegaPixelTracksCand')
process.looseOmeFilter.algo = cms.string('looseOmegaPixelTracks')

#########analyzers#################
process.dca2loose = cms.EDAnalyzer('TestKinkNtuplizer',
    VeeAlgo       = cms.string('looseV0Candidates'),
#    CasAlgo       = cms.string('looseCascade'),
    CasAlgo       = cms.string('looseXiPixelTracks'),
    CasDecayName  = cms.untracked.string('Cascade'),
    CasMass       = cms.untracked.double(1.32171),
    trackingAlgo  = cms.InputTag('generalTracks')
,  MuonsLabel     = cms.InputTag('muons')
)
process.dca2fit2KF = cms.EDAnalyzer('TestKinkNtuplizer',
    VeeAlgo       = cms.string('looseV0Candidates'),
    CasAlgo       = cms.string('CascCandidates'),
#    CasAlgo       = cms.string('XiPixelTracksCand'),
    CasDecayName  = cms.untracked.string('Cascade'),
    CasMass       = cms.untracked.double(1.32171),
    trackingAlgo  = cms.InputTag('generalTracks')
,  MuonsLabel     = cms.InputTag('muons')
)
process.odca2loose = cms.EDAnalyzer('TestKinkNtuplizer',
    xmin          = cms.untracked.double(1.62),
    xmax          = cms.untracked.double(1.72),
    VeeAlgo       = cms.string('looseV0Candidates'),
#    CasAlgo       = cms.string('looseCascade'),
    CasAlgo       = cms.string('looseOmegaPixelTracks'),
    CasDecayName  = cms.untracked.string('Omega'),
    CasMass       = cms.untracked.double(1.67245),
    trackingAlgo  = cms.InputTag('generalTracks')
,  MuonsLabel     = cms.InputTag('muons')
)
process.odca2fit2KF = cms.EDAnalyzer('TestKinkNtuplizer',
    xmin          = cms.untracked.double(1.62),
    xmax          = cms.untracked.double(1.72),
#    pitrk_pt_cut  = cms.untracked.double(1.0),
    VeeAlgo       = cms.string('looseV0Candidates'),
    CasAlgo       = cms.string('CascCandidates'),
#    CasAlgo       = cms.string('OmegaPixelTracksCand'),
    CasDecayName  = cms.untracked.string('Omega'),
    CasMass       = cms.untracked.double(1.67245),
    trackingAlgo  = cms.InputTag('generalTracks')
,  MuonsLabel     = cms.InputTag('muons')
)
###############################################################
#process.noscrapingnobeamhalo = cms.Path(process.L1T1*process.noScraping)
process.casevts = cms.Path(
#process.L1T1*process.noScraping*
process.primaryVertexFilter
*process.looseV0Candidates*process.LambdaFilter
*process.looseCascade
*process.looseXiPixelTracks
*process.looseCasFilter
*process.dca2loose
)
process.casevtskf = cms.Path( 
#process.L1T1*process.noScraping*
process.primaryVertexFilter
*process.looseV0Candidates*process.LambdaFilter
*process.CascCandidates
#*process.XiPixelTracksCand
*process.CasFilter
*process.dca2fit2KF
)
process.omeevts = cms.Path(
#process.L1T1*process.noScraping*
process.primaryVertexFilter
*process.looseV0Candidates*process.LambdaFilter
*process.looseCascade
#*process.looseOmegaPixelTracks
*process.looseOmeFilter
*process.odca2loose
)
process.omeevtskf = cms.Path(
#process.L1T1*process.noScraping*
process.primaryVertexFilter
*process.looseV0Candidates*process.LambdaFilter
*process.CascCandidates
#*process.OmegaPixelTracksCand
*process.OmeFilter
*process.odca2fit2KF
)

#process.e = cms.EndPath(process.out)
process.schedule= cms.Schedule(
#process.casevts
#,
#process.casevtskf
#,process.omeevts
#,
process.omeevtskf
#,process.e
)

#process.p = cms.Path(
#process.primaryVertexFilter
#*process.looseV0Candidates*process.LambdaFilter
#*(process.looseCascade
#*(process.looseXiPixelTracks*process.looseCasFilter*process.dca2loose
#+process.looseOmegaPixelTracks*process.looseOmeFilter*process.odca2loose)
#+process.CascCandidates
#*(process.XiPixelTracksCand*process.CasFilter*process.dca2fit2KF
#+process.OmegaPixelTracksCand*process.OmeFilter*process.odca2fit2KF)
#)
#)
