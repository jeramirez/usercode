import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
from HecBmeson.DimuPiK.Cert_160404_180252_7TeV_PromptReco_Collisions11_CMSSWConfig_MuonPhys_Hec_cfi import *

## import skeleton process
#from PhysicsTools.PatAlgos.patTemplate_cfg import *
from PhysicsTools.PatAlgos.tools.coreTools import *

#--The above and below lines comes from Charles twiki page:
#  https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideAboutPythonConfigFile

# setup 'standard'  options
options = VarParsing.VarParsing ('standard')
options.register ('output2',
                  'hola1.root',                                   # default value
                  VarParsing.VarParsing.multiplicity.singleton,  # singleton or list
                  VarParsing.VarParsing.varType.string,          # string, int, or float
                  "Outputfile name I want to")

# setup any defaults you want
options.output = '/uscms/home/cplager/nobackup/outputFiles/try_3.root'
options.files= 'file1.root', 'file2.root'
options.maxEvents = 100              # -1 means all events

# get and parse the command line arguments
options.parseArguments()
#output2 = options.files[0]
#output2 = options.secondaryOutput
myoutput = options.output2
print '    --> this is my output2: ', myoutput
myinput  = options.files[0]
print '    --> this is my input: ', myinput
#---------------------------------------------

process = cms.Process("MyDimuPiK")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

process.load("PhysicsTools.PatAlgos.patSequences_cff")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration/StandardSequences/Geometry_cff')

#--look for the tag:  https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'GR_R_42_V19::All' 

process.MessageLogger.cerr.FwkReport.reportEvery = 10

process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(myoutput+".root")
)
process.source = cms.Source("PoolSource",
    #--skipEvents=cms.untracked.uint32(4000),
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(    
   #'file:myfile.root'
   'file:/eos/uscms/store/user/pivarski/Z4430/RECONSTRUCTED_Z4430_00118.root'
  ,'file:/eos/uscms/store/user/pivarski/Z4430/RECONSTRUCTED_Z4430_00029.root'
   # The file below are at FNAL at:  /pnfs/cms/WAX/11/
   # at charma /data/se/cms/
   #myinput   
    )
)
#--process.source.lumisToProcess = lumisToProcess7TeVCollisions11Nov14   #--JSON File
#------------------PAT
process.load("PhysicsTools.PatAlgos.cleaningLayer1.genericTrackCleaner_cfi")
process.cleanPatTracks.checkOverlaps.muons.requireNoOverlaps = cms.bool(False)
process.cleanPatTracks.checkOverlaps.electrons.requireNoOverlaps = cms.bool(False)

process.load("PhysicsTools.PatAlgos.mcMatchLayer0.muonMatch_cfi")
process.muonMatch.matched = cms.InputTag("genParticlesPlusSim")
process.muonMatch.maxDeltaR = cms.double(0.02)
process.muonMatch.resolveByMatchQuality = cms.bool(True)

## add track candidates
from PhysicsTools.PatAlgos.tools.trackTools import *
makeTrackCandidates(process,
     label        = 'TrackCands',                  
     tracks       = cms.InputTag('generalTracks'), 
     particleType = 'pi+',                         
     preselection = 'pt > 0.1',                     
     selection    = 'pt > 0.1',                     
     isolation    = {},                            
     isoDeposits  = [],                            
     mcAs         = None          
)
removeMCMatching(process, ['All'],outputInProcess = False )

# rerun genparticles
process.load('SimGeneral.HepPDTESSource.pdt_cfi')
process.genParticlesRemake = cms.EDProducer("GenParticleProducer",
   saveBarCodes = cms.untracked.bool(True),
   src = cms.InputTag("generator"),
   abortOnUnknownPDGCode = cms.untracked.bool(False)
)
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.genParticlesPlusSim = cms.EDProducer("GenPlusSimParticleProducer",
	src	      = cms.InputTag("g4SimHits"), # use "famosSimHits" for FAMOS
	setStatus     = cms.int32(8),
	filter        = cms.vstring("pt > 0.0"),  
	genParticles   = cms.InputTag("genParticlesRemake")
)  

#--------------------------------------------
process.load('RecoVertex.V0Producer.generalV0Candidates_cff')
process.generalV0Candidates.innerHitPosCut = cms.double(-1)   #--AOD
process.Myk0s = process.generalV0Candidates.clone()
process.Myk0s.vtxSignificance2DCut  = cms.double(5.0)    #--default is 15
process.Myk0s.impactParameterSigCut = cms.double(0.5)    #--default is 2

process.HecExotic = cms.EDAnalyzer('DimuPiK'
   , tracks     = cms.InputTag('cleanPatTrackCands')
   , MuonsLabel = cms.InputTag('cleanPatMuonsTriggerMatch')
   , VeeAlgo    = cms.untracked.string("Myk0s")
   , MyPrint    = cms.untracked.bool(False)          #--cout for short run
   , ChargedKa  = cms.untracked.bool(False)          #--true=u+u-pi+k-    false:u+u-pi+k0s
   , HLTriggerResults = cms.untracked.string('HLT')
)

#----------------------PAT
process.patDefaultSequence.remove(process.patJetCorrFactors)
process.patDefaultSequence.remove(process.patJetCharge)
process.patDefaultSequence.remove(process.patJetPartonMatch)
process.patDefaultSequence.remove(process.patJetGenJetMatch)
process.patDefaultSequence.remove(process.patJetPartons)
process.patDefaultSequence.remove(process.patJetPartonAssociation)
process.patDefaultSequence.remove(process.patJetFlavourAssociation)
process.patDefaultSequence.remove(process.patJets)
process.patDefaultSequence.remove(process.metJESCorAK5CaloJet)
process.patDefaultSequence.remove(process.metJESCorAK5CaloJetMuons)
process.patDefaultSequence.remove(process.patMETs)
process.patDefaultSequence.remove(process.selectedPatJets)
process.patDefaultSequence.remove(process.cleanPatJets)
process.patDefaultSequence.remove(process.countPatJets)

# do trigger matching between PAT muons and HLT muons
process.cleanMuonTriggerMatchHLT = cms.EDProducer( "PATTriggerMatcherDRLessByR"                  # match by DeltaR only, best match by DeltaR
, src                   = cms.InputTag( "cleanPatMuons" )
, matched               = cms.InputTag( "patTrigger" )    # default producer label as defined in PhysicsTools/PatAlgos/python/triggerLayer1/triggerProducer_cfi.py
, matchedCuts           = cms.string( 'type( "TriggerMuon" )' )
, maxDeltaR             = cms.double( 0.10 )
, resolveAmbiguities    = cms.bool( False )       # only one match per trigger object
, resolveByMatchQuality = cms.bool( False )      # take best match found per reco object: by DeltaR here (s. above)
)

# Embedding in muons
process.cleanPatMuonsTriggerMatch = cms.EDProducer(
  "PATTriggerMatchMuonEmbedder"
, src     = cms.InputTag(  "cleanPatMuons" )
, matches = cms.VInputTag( 'cleanMuonTriggerMatchHLT' )
)
process.trigMatch = cms.Sequence(process.cleanMuonTriggerMatchHLT *process.cleanPatMuonsTriggerMatch )
## let it run
process.pat = cms.Path(
    process.patDefaultSequence
)

# switch on PAT trigger info
from PhysicsTools.PatAlgos.tools.trigTools import *
switchOnTrigger( process, outputModule = '' )

#----------------------------
#process.out = cms.OutputModule("PoolOutputModule",
#                               fileName = cms.untracked.string('patTuple.root'),
                               # save only events passing the full path
                               ##SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
                               # save PAT Layer 1 output; you need a '*' to
                               # unpack the list of commands 'patEventContent'
#                               outputCommands = cms.untracked.vstring('drop *', *patEventContent )
#                               )

#process.outpath = cms.EndPath(process.out)

process.phec = cms.Path(process.trigMatch*process.Myk0s*process.HecExotic)
process.trigMatchPath = cms.Path( process.trigMatch )

#process.schedule = cms.Schedule(process.pat, process.trigMatchPaths)
process.schedule = cms.Schedule(process.pat, process.trigMatchPath, process.phec )
