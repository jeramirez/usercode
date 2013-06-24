#config file to convert input data into patlayer format
#to be used later as input data for bjatpat analyzer
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

#arguments
options = VarParsing.VarParsing ('standard')
options.register ('mypatsequence',
                  'PhysicsTools.PatAlgos.patSequences_cff',      # default value
                  VarParsing.VarParsing.multiplicity.singleton,  # singleton or list
                  VarParsing.VarParsing.varType.string,          # string, int, or float
                  "PAT Layer 0+1 Sequences")

options.output = 'PATLayer1_Output.fromAOD_ttbar1_full.root'
options.files= '/store/relval/CMSSW_3_1_0_pre9/RelValTTbar/GEN-SIM-RECO/IDEAL_31X_v1/0007/DE732988-5E4F-DE11-82ED-001D09F25208.root',\
       '/store/relval/CMSSW_3_1_0_pre9/RelValTTbar/GEN-SIM-RECO/IDEAL_31X_v1/0007/82FD1C7A-6E4F-DE11-9198-0019B9F72CC2.root',\
       '/store/relval/CMSSW_3_1_0_pre9/RelValTTbar/GEN-SIM-RECO/IDEAL_31X_v1/0007/82B99E83-5E4F-DE11-96FA-001D09F28EA3.root',\
       '/store/relval/CMSSW_3_1_0_pre9/RelValTTbar/GEN-SIM-RECO/IDEAL_31X_v1/0007/345DA7B5-F64E-DE11-8C23-001617DBD556.root',\
       '/store/relval/CMSSW_3_1_0_pre9/RelValTTbar/GEN-SIM-RECO/IDEAL_31X_v1/0007/22620711-524F-DE11-9D47-001617C3B65A.root',\
       '/store/relval/CMSSW_3_1_0_pre9/RelValTTbar/GEN-SIM-RECO/IDEAL_31X_v1/0007/1C645B87-5E4F-DE11-A7BB-000423D985E4.root',\
       '/store/relval/CMSSW_3_1_0_pre9/RelValTTbar/GEN-SIM-RECO/IDEAL_31X_v1/0007/14981583-5E4F-DE11-8A22-001D09F253C0.root'

options.maxEvents = 10
options.mypatsequence = 'Analyzers.BTagPAT.patSequences_cff'
options.parseArguments()

process = cms.Process("PAT")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('PATSummaryTables')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    default          = cms.untracked.PSet( limit = cms.untracked.int32(0)  ),
    PATSummaryTables = cms.untracked.PSet( limit = cms.untracked.int32(-1) )
)
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring( options.files )
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( options.maxEvents ) )

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('IDEAL_31X::All')
process.load("Configuration.StandardSequences.MagneticField_cff")

# PAT Layer 0+1
process.load(options.mypatsequence)
process.content = cms.EDAnalyzer("EventContentAnalyzer")

# replacements currently needed to make the taus work
process.allLayer1Taus.addTauID = False

process.p = cms.Path(
    process.patDefaultSequence  
)

# Output module configuration
#from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
from Analyzers.BTagPAT.patEventContentEd_cff import patEventContent
from Analyzers.BTagPAT.patEventContentEd_cff import patEventContentallLayer1
from Analyzers.BTagPAT.patEventContentEd_cff import patEventContentNoLayer1Cleaning
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string(options.output),
    # save only events passing the full path
    SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
    # save PAT Layer 1 output
    outputCommands = cms.untracked.vstring('drop *', *patEventContentNoLayer1Cleaning + patEventContentallLayer1 ) # you need a '*' to unpack the list of commands 'patEventContent'
)
process.outpath = cms.EndPath(process.out)

