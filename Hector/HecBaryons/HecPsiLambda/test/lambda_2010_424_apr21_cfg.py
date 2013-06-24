import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
from HecBaryons.HecPsiLambda.Cert_136033_149442_7TeV_Dec22ReReco_Collisions10_JSON_v4_cfi import *

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

process = cms.Process("MyHecPsiLambda")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration/StandardSequences/Geometry_cff')

#--look for the tag:  https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'GR_R_42_V19::All' 

process.MessageLogger.cerr.FwkReport.reportEvery = 200

process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(myoutput+".root")
)
process.source = cms.Source("PoolSource",
    #--skipEvents=cms.untracked.uint32(4000),
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(    
   #'file:myfile.root'
   # The file below are at FNAL at:  /pnfs/cms/WAX/11/
   # at charma /data/se/cms/
   myinput   
    )
)
process.source.lumisToProcess = lumisToProcess7TeVCollisions10    #--Eduardo Lumi

process.HecJL = cms.EDAnalyzer('HecPsiLambda'
   , tracks     = cms.InputTag('generalTracks')
   , MuonsLabel = cms.InputTag('muons')
   , MyPrint    = cms.untracked.bool(False)                  #--cout for short run
)


process.p = cms.Path(process.HecJL)
