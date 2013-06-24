import os
import sys

import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("TRIGGERANDVEESKIM")

options=VarParsing.VarParsing('standard')
#defaults
options.output    = 'triggerandlooseveeskimmed0179.root'
options.maxEvents = 10
options.parseArguments()

process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.3 $'),
    name = cms.untracked.string('$Source: /local/reps/CMSSW/UserCode/PuertoRicoGroup/Analyzers/CascadeProducer/test/TriggerandLooseVeeSkim_cfg.py,v $'),
    annotation = cms.untracked.string('Xi input skim')
)

#
#
# This is for testing purposes.
#
#
# run 123151 lumisection 14

readFiles   = cms.untracked.vstring(options.files)
secFiles    = cms.untracked.vstring()
outPutFile  = cms.untracked.string(options.output)
process.source = cms.Source ("PoolSource",fileNames = readFiles)


process.source.inputCommands = cms.untracked.vstring("keep *", "drop *_MEtoEDMConverter_*_*", "drop L1GlobalTriggerObjectMapRecord_hltL1GtObjectMap__HLT")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)


#------------------------------------------
# Load standard sequences.
#------------------------------------------
process.load('Configuration/StandardSequences/MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration/StandardSequences/Geometry_cff')


process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'GR_R_35X_V6::All' 

process.load("Configuration/StandardSequences/RawToDigi_Data_cff")
process.load("Configuration/StandardSequences/Reconstruction_cff")
process.load('Configuration/EventContent/EventContent_cff')

process.FEVTEventContent.outputCommands.append('drop *_MEtoEDMConverter_*_*')
process.FEVTEventContent.outputCommands.append('keep *_looseV0Candidates_*_*')

##################################bscnobamhalo############################################
process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff')
process.load('HLTrigger/HLTfilters/hltLevel1GTSeed_cfi')

process.L1T1=process.hltLevel1GTSeed.clone()
process.L1T1.L1TechTriggerSeeding = cms.bool(True)
#process.L1T1.L1SeedsLogicalExpression = cms.string('(40 OR 41) AND NOT (36 OR 37 OR 38 OR 39)')
process.L1T1.L1SeedsLogicalExpression = cms.string('0 AND (40 OR 41) AND NOT (36 OR 37 OR 38 OR 39)')

######################################Vee Track#################################################

process.load('Analyzers.CascadeProducer.looseV0Candidates_cff')
process.LambdaFilter = cms.EDFilter("VeeCountFilter",
                                                      algo = cms.string('looseV0Candidates'),
                                                      name = cms.string('Lambda'),
                                                      minNumber = cms.uint32(1) 
                                                      )

#process.VeePath = cms.Path(process.L1T1*process.looseV0Candidates*process.LambdaFilter)
process.VeeLoose = cms.Path(process.looseV0Candidates)
process.VeePath = cms.Path(process.L1T1*process.LambdaFilter)

process.recoveeout = cms.OutputModule("PoolOutputModule",
                               outputCommands = process.FEVTEventContent.outputCommands,
                               SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('VeePath'
										)),
                               fileName = outPutFile
                               )

###########################################################################################

process.outpath = cms.EndPath(process.recoveeout)
process.schedule= cms.Schedule(process.VeeLoose,process.VeePath,process.outpath)


