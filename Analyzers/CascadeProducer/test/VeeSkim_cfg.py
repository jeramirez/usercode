import FWCore.ParameterSet.Config as cms

process = cms.Process("VEESKIM")

process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.3 $'),
    name = cms.untracked.string('$Source: /local/reps/CMSSW/UserCode/PuertoRicoGroup/Analyzers/CascadeProducer/test/VeeSkim_cfg.py,v $'),
    annotation = cms.untracked.string('Xi input skim')
)

#
#
# This is for testing purposes.
#
#
# run 123151 lumisection 14
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
#1'/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0103/3C8F9421-AFEE-DE11-8AA7-0024E87687BE.root',
#1'/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0103/36593A1C-AFEE-DE11-A1F3-0024E8768D5B.root',
#1'/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0103/E06A076B-AFEE-DE11-ADE6-001D0967D0FD.root',
#1'/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0103/621A02D0-AFEE-DE11-98CE-0024E876A814.root',
#1'/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0103/A234211B-AFEE-DE11-90D7-0024E86E8DA7.root',
#1'/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0103/72A3CB9C-ADEE-DE11-A3E6-001D0967D0FD.root',
#1'/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0103/C4E0F513-AFEE-DE11-B945-0024E87699E7.root',
#1'/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0103/740F5E1D-AFEE-DE11-AD64-0024E876842C.root',
#1'/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0103/04B672AB-ADEE-DE11-B9A0-0024E8769958.root',
#1'/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0103/66E7B6AB-ADEE-DE11-AFB6-001D0967D698.root'
#2'/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0103/021F9EA5-ADEE-DE11-B8C1-0024E876A814.root',
#2'/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0103/767C15A3-ADEE-DE11-8F50-0024E87680C0.root',
#2'/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0103/EC718E1C-AFEE-DE11-B222-0024E87663E1.root',
#2'/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0103/EE7D0A5B-AFEE-DE11-8461-0024E87687BE.root',
#2'/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0103/D6FC56AE-ADEE-DE11-B2F7-0024E87680E7.root',
#2'/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0103/1A552B1E-AFEE-DE11-858B-0024E8768BD5.root',
#2'/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0103/1C1DE5DA-AFEE-DE11-B5D5-0024E86E8DA7.root',
#2'/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0103/16479F93-AFEE-DE11-AEBC-0024E8768D5B.root',
#2'/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0103/78959DAA-ADEE-DE11-83E2-0024E876A807.root',
#2'/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0103/86869228-AFEE-DE11-9552-001D0967D55D.root'
'/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0103/4018241E-AFEE-DE11-B040-0024E86E8D18.root',
'/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0103/8869E1D0-AFEE-DE11-94DD-001D0967D698.root',
'/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0103/102E591C-AFEE-DE11-AF43-0024E8768C57.root',
'/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0103/B27D7FB0-ADEE-DE11-A452-001D0967D6AC.root',
'/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0103/C64F14AC-FBED-DE11-A823-0015178C6480.root',
'/store/data/BeamCommissioning09/MinimumBias/RECO/Dec19thReReco_336p3_v2/0103/AA39E3C5-2BEE-DE11-B0E8-00151796D4B4.root'

)
)

process.source.inputCommands = cms.untracked.vstring("keep *", "drop *_MEtoEDMConverter_*_*", "drop L1GlobalTriggerObjectMapRecord_hltL1GtObjectMap__HLT")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)


#------------------------------------------
# Load standard sequences.
#------------------------------------------
process.load('Configuration/StandardSequences/MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration/StandardSequences/GeometryIdeal_cff')


process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'GR09_R_V5::All' 

process.load("Configuration/StandardSequences/RawToDigi_Data_cff")
process.load("Configuration/StandardSequences/Reconstruction_cff")
process.load('Configuration/EventContent/EventContent_cff')

process.FEVTEventContent.outputCommands.append('drop *_MEtoEDMConverter_*_*')


######################################Vee Track#################################################

process.LambdaFilter = cms.EDFilter("VeeCountFilter",
                                                      algo = cms.string('generalV0Candidates'),
                                                      name = cms.string('Lambda'),
                                                      minNumber = cms.uint32(1) 
                                                      )

process.VeePath = cms.Path(process.LambdaFilter)

process.recoveeout = cms.OutputModule("PoolOutputModule",
                               outputCommands = process.FEVTEventContent.outputCommands,
                               SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('VeePath'
										)),
                               fileName = cms.untracked.string('veeskimmed.root')
                               )

###########################################################################################

process.outpath = cms.EndPath(process.recoveeout)



