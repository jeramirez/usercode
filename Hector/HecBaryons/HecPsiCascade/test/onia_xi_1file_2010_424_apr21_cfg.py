import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
from HecBaryons.HecPsiCascade.Cert_136033_149442_7TeV_Dec22ReReco_Collisions10_JSON_v4_cfi import *

#--The above and below lines comes from Charles twiki page:
#  https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideAboutPythonConfigFile

# setup 'standard'  options
options = VarParsing.VarParsing ('standard')
options.register ('output2',
                  'hola.root',                                   # default value
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

process = cms.Process("HecPpsiCascade")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.MagneticField_cff")

#process.load("RecoEgamma.PhotonIdentification.photonId_cff")
#process.load("Geometry.CaloEventSetup.CaloGeometry_cfi")
#process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")

#process.load("Geometry.CaloEventSetup.CaloTopology_cfi")
#process.load("RecoEgamma.EgammaPhotonProducers.photonSequence_cff")
#process.load("RecoEcal.EgammaClusterProducers.islandBasicClusters_cfi")

process.load('Configuration/StandardSequences/Geometry_cff')
process.load('Geometry.TrackerGeometryBuilder.trackerGeometry_cfi')
process.TrackerDigiGeometryESModule.applyAlignment = False

#--look for the tag:  https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'GR_R_42_V19::All' 

#process.load("Configuration/StandardSequences/RawToDigi_Data_cff")
#process.load("Configuration/StandardSequences/Reconstruction_cff")
#process.load('Configuration/EventContent/EventContent_cff')

#process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
#process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 200
#process.MessageLogger.destinations = cms.untracked.vstring(myoutput+".log")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(myoutput+".root")
)

process.source = cms.Source("PoolSource",
    #------skipEvents=cms.untracked.uint32(4000),
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(    
   #'file:myfile.root'
   # The file below are at FNAL at:  /pnfs/cms/WAX/11/
   myinput   
    )
)
process.source.lumisToProcess = lumisToProcess7TeVCollisions10    #--Eduardo Lumi

process.HecMuons = cms.EDAnalyzer('HecPsiCascade'
  #, tracks = cms.untracked.InputTag('ctfWithMaterialTracks')
   , tracks         = cms.InputTag('generalTracks')
   , gsfElectrons   = cms.InputTag('gsfElectrons')
   , MuonsLabel     = cms.InputTag('muons')
   , CaloMuonsLabel = cms.InputTag("calomuons")
   , minTracks      = cms.untracked.uint32(0)
   , check_elC      = cms.untracked.bool(True)                   #--Dump e- matching info
   , Pe_match       = cms.untracked.double(10.0)
   , Mass_Constrain = cms.string("m_jpsi")                       #--Mass constrain fit
   , onlyDiMu       = cms.untracked.bool(False)                  #--run only Dimuons or DiMu + X
   , MyPrint        = cms.untracked.bool(False)                  #--cout for short run
   , CasAlgo        = cms.untracked.string("MyCascade")          #--Label of Cascade collection (My clone)
   , CasDecayName   = cms.untracked.string("Cascade")            #--Collection Label decay (cascade or omega)(Eduardo}
   , VeeAlgo        = cms.untracked.string("looseV0Candidates")  #--Collection Label decay (Lambda) (Eduardo}
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
    thresh = cms.untracked.double(0.25)
)

####################apply good primary vertex filter here#####################
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                           vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                                           minimumNDOF = cms.uint32(4),
                                           maxAbsZ = cms.double(15), 
                                           maxd0 = cms.double(2) 
)

process.dump=cms.EDAnalyzer('EventContentAnalyzer')

#--Cascade
process.load('Analyzers.CascadeProducer.looseV0Candidates_cff')
process.load('Analyzers.CascadeProducer.cascadeproducer_cff')
process.MyCascade = process.looseCascade.clone()
process.MyCascade.XiDCACut    = cms.double(2.0)                      #--closest aproach Lambda-pion < 2 cm
process.MyCascade.VtxFromFit  = cms.bool(True)
process.MyCascade.MassFromFit = cms.bool(True)
process.MyCascade.LambdaMassWidthCut = cms.double(0.030)             #--30 MeV Lambda Window
#process.MyCascade.v0Algo      = cms.string('generalV0Candidates')   #--cms default V0  L/S > 15
#process.MyCascade.v0Algo      = cms.string('looseV0Candidates')     #--cms V0 with L/S > 5
process.LambdaFilter = cms.EDFilter("VeeCountFilter",
                                                      algo = cms.string('looseV0Candidates'),
                                                      name = cms.string('Lambda'),
                                                      minNumber = cms.uint32(1) 
                                                      )

#process.p = cms.Path(process.HecMuons*process.dump)
#process.p = cms.Path(process.L1T1*process.noScraping*process.primaryVertexFilter*process.looseV0Candidates*process.LambdaFilter*process.MyCascade*process.HecMuons)
process.p = cms.Path(process.looseV0Candidates*process.LambdaFilter*process.MyCascade*process.HecMuons)
#process.p = cms.Path(process.HecMuons)
