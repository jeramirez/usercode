import os
import sys

import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("OWNPARTICLES")

#setup options
options = VarParsing.VarParsing('standard')
#setup defaults
options.output    = 'cascproducerlooseveeskimmall2010_09.root'
options.files     = 'file:/hep/alpha01/home/CMSDATA/vee358_2010/RAW-RECO/Vmasscut/triggerandlooseveeskimmed0179_09.root'
options.maxEvents = 1
options.parseArguments()

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

readFiles  = cms.untracked.vstring(options.files)
secFiles   = cms.untracked.vstring()
outPutFile = cms.untracked.string(options.output)
process.source = cms.Source("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
# replace 'myfile.root' with the source file you want to use
#readFiles.extend((
#        'file:/hep/alpha01/home/CMSDATA/vee358_2010/RAW-RECO/Vmasscut/triggerandlooseveeskimmed0179_09.root'
#        'file:/hep/alpha01/home/CMSDATA/vee358_2010/RAW-RECO/Vmasscut/triggerandlooseveeskimmed0179_10.root'
#        'file:/hep/alpha01/home/CMSDATA/vee358_2010/RAW-RECO/Vmasscut/triggerandlooseveeskimmed0179_11.root'
#        'file:/hep/alpha01/home/CMSDATA/vee358_2010/RAW-RECO/Vmasscut/triggerandlooseveeskimmed0179_12.root'
#        'file:/hep/alpha01/home/CMSDATA/vee358_2010/RAW-RECO/Vmasscut/triggerandlooseveeskimmed0179_13.root'
#        'file:/hep/alpha01/home/CMSDATA/vee358_2010/RAW-RECO/Vmasscut/triggerandlooseveeskimmed0179_14.root'
#        'file:/hep/alpha01/home/CMSDATA/vee358_2010/RAW-RECO/Vmasscut/triggerandlooseveeskimmed0179_15.root',
#        'file:/hep/alpha01/home/CMSDATA/vee358_2010/RAW-RECO/Vmasscut/triggerandlooseveeskimmed0179_16.root',
#));

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = cms.string('GR09_R_V5::All')
process.GlobalTag.globaltag = cms.string('GR_R_35X_V8::All')
process.load("Configuration.StandardSequences.MagneticField_cff")
#process.load("RecoVertex.V0Producer.generalV0Candidates_cff")
process.load('Analyzers.CascadeProducer.looseV0Candidates_cff')
process.load('Analyzers.CascadeProducer.cascadeproducer_cff')

process.CascCandidatesDCA02 = process.looseCascade.clone()
process.CascCandidatesDCA02.XiDCACut           = cms.double(0.2)
process.CascCandidatesDCA1  = process.looseCascade.clone()
process.CascCandidatesDCA1.XiDCACut            = cms.double(1.0)
process.CascCandidatesDCA2  = process.looseCascade.clone()
process.CascCandidatesDCA02.XiDCACut           = cms.double(2.0)
process.CascCandidatesDCA5  = process.looseCascade.clone()
process.CascCandidatesDCA5.XiDCACut            = cms.double(5.0)

process.CascCandidatesDCA02fitvtx = process.CascCandidatesDCA02.clone()
process.CascCandidatesDCA02fitvtx.VtxFromFit         = cms.bool(True)

process.CascCandidatesDCA1fitvtx = process.CascCandidatesDCA1.clone()
process.CascCandidatesDCA1fitvtx.VtxFromFit         = cms.bool(True)

process.CascCandidatesDCA2fitvtx = process.CascCandidatesDCA2.clone()
process.CascCandidatesDCA2fitvtx.VtxFromFit         = cms.bool(True)

process.CascCandidatesDCA5fitvtx = process.CascCandidatesDCA5.clone()
process.CascCandidatesDCA5fitvtx.VtxFromFit         = cms.bool(True)

process.CascCandidatesDCA02fit2 = process.CascCandidatesDCA02fitvtx.clone()
process.CascCandidatesDCA02fit2.MassFromFit        = cms.bool(True)

process.CascCandidatesDCA1fit2 = process.CascCandidatesDCA1fitvtx.clone()
process.CascCandidatesDCA1fit2.MassFromFit        = cms.bool(True)

process.CascCandidatesDCA2fit2 = process.CascCandidatesDCA2fitvtx.clone()
process.CascCandidatesDCA2fit2.MassFromFit        = cms.bool(True)

process.CascCandidatesDCA5fit2 = process.CascCandidatesDCA5fitvtx.clone()
process.CascCandidatesDCA5fit2.MassFromFit        = cms.bool(True)

##################################bscnobamhalo############################################
#already done
#process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff')
#process.load('HLTrigger/HLTfilters/hltLevel1GTSeed_cfi')
#
#process.L1T1=process.hltLevel1GTSeed.clone()
#process.L1T1.L1TechTriggerSeeding = cms.bool(True)
##process.L1T1.L1SeedsLogicalExpression = cms.string('(40 OR 41) AND NOT (36 OR 37 OR 38 OR 39)')
#process.L1T1.L1SeedsLogicalExpression = cms.string('0 AND (40 OR 41) AND NOT (36 OR 37 OR 38 OR 39)')
#
#process.bscnobeamhalo = cms.Path(process.L1T1)
#
##################################cas filter ############################################
process.Cas02Filter = process.looseCasFilter.clone()
process.Cas02Filter.algo = cms.string('CascCandidatesDCA02')

process.Cas1Filter = process.looseCasFilter.clone()
process.Cas1Filter.algo = cms.string('CascCandidatesDCA1')

process.Cas2Filter = process.looseCasFilter.clone()
process.Cas2Filter.algo = cms.string('CascCandidatesDCA2')

process.Cas5Filter = process.looseCasFilter.clone()
process.Cas5Filter.algo = cms.string('CascCandidatesDCA5')

process.Ome02Filter = process.looseOmeFilter.clone()
process.Ome02Filter.algo = cms.string('CascCandidatesDCA02')

process.Ome1Filter = process.looseOmeFilter.clone()
process.Ome1Filter.algo = cms.string('CascCandidatesDCA1')

process.Ome2Filter = process.looseOmeFilter.clone()
process.Ome2Filter.algo = cms.string('CascCandidatesDCA2')

process.Ome5Filter = process.looseOmeFilter.clone()
process.Ome5Filter.algo = cms.string('CascCandidatesDCA5')

process.Cas02fitFilter = cms.EDFilter("VeeCountFilter",
                                                      algo = cms.string('CascCandidatesDCA02fit2'),
                                                      name = cms.string('Cascade'),
                                                      minNumber = cms.uint32(1)
                                                      )
process.Cas1fitFilter = cms.EDFilter("VeeCountFilter",
                                                      algo = cms.string('CascCandidatesDCA1fit2'),
                                                      name = cms.string('Cascade'),
                                                      minNumber = cms.uint32(1)
                                                      )
process.Cas2fitFilter = cms.EDFilter("VeeCountFilter",
                                                      algo = cms.string('CascCandidatesDCA2fit2'),
                                                      name = cms.string('Cascade'),
                                                      minNumber = cms.uint32(1)
                                                      )
process.Cas5fitFilter = cms.EDFilter("VeeCountFilter",
                                                      algo = cms.string('CascCandidatesDCA5fit2'),
                                                      name = cms.string('Cascade'),
                                                      minNumber = cms.uint32(1)
                                                      )
process.Ome02fitFilter = cms.EDFilter("VeeCountFilter",
                                                      algo = cms.string('CascCandidatesDCA02fit2'),
                                                      name = cms.string('Omega'),
                                                      minNumber = cms.uint32(1)
                                                      )
process.Ome1fitFilter = cms.EDFilter("VeeCountFilter",
                                                      algo = cms.string('CascCandidatesDCA1fit2'),
                                                      name = cms.string('Omega'),
                                                      minNumber = cms.uint32(1)
                                                      )
process.Ome2fitFilter = cms.EDFilter("VeeCountFilter",
                                                      algo = cms.string('CascCandidatesDCA2fit2'),
                                                      name = cms.string('Omega'),
                                                      minNumber = cms.uint32(1)
                                                      )
process.Ome5fitFilter = cms.EDFilter("VeeCountFilter",
                                                      algo = cms.string('CascCandidatesDCA5fit2'),
                                                      name = cms.string('Omega'),
                                                      minNumber = cms.uint32(1)
                                                      )

process.casevts02 = cms.Path(process.Cas02Filter )
process.casevts1  = cms.Path(process.Cas1Filter )
process.casevts2  = cms.Path(process.Cas2Filter )
process.casevts5  = cms.Path(process.Cas5Filter )
process.omeevts02 = cms.Path(process.Ome02Filter )
process.omeevts1  = cms.Path(process.Ome1Filter )
process.omeevts2  = cms.Path(process.Ome2Filter )
process.omeevts5  = cms.Path(process.Ome5Filter )

process.out = cms.OutputModule("PoolOutputModule",
    fileName = outPutFile,
    SelectEvents = cms.untracked.PSet(
#        SelectEvents = cms.vstring('bscnobeamhalo')
        SelectEvents = cms.vstring('casevts02','casevts1','casevts2','casevts5','omeevts02','omeevts1','omeevts2','omeevts5')
    )
)

  
process.p = cms.Path(process.CascCandidatesDCA02 
                   + process.CascCandidatesDCA02fitvtx
                   + process.CascCandidatesDCA02fit2
                   + process.CascCandidatesDCA1
                   + process.CascCandidatesDCA1fitvtx
                   + process.CascCandidatesDCA1fit2
                   + process.CascCandidatesDCA2
                   + process.CascCandidatesDCA2fitvtx
                   + process.CascCandidatesDCA2fit2
                   + process.CascCandidatesDCA5
                   + process.CascCandidatesDCA5fitvtx
                   + process.CascCandidatesDCA5fit2
)

process.e = cms.EndPath(process.out)
process.schedule= cms.Schedule(process.p
,process.casevts02
,process.casevts1
,process.casevts2
,process.casevts5
,process.omeevts02
,process.omeevts1
,process.omeevts2
,process.omeevts5
,process.e)
