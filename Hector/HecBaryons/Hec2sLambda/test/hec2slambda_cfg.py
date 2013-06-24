import FWCore.ParameterSet.Config as cms

process = cms.Process("MyHec2sLambda")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.MagneticField_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
       #--'file:myfile.root'
       #--'file:/uscms_data/d3/maria22/00390632-5BA1-E011-9427-0017A4770C2C.root'
 #--'/store/data/Run2010A/MuOnia/RECO/Apr21ReReco-v1/0000/261CA6C4-526E-E011-B236-E0CB4E29C504.root'
'/store/user/maria22/store/mc/Summer11/LambdaBToPsiMuMu_2MuPEtaFilter_Tight_7TeV-pythia6-evtgen/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/022C6AB9-57A1-E011-9623-1CC1DE1CDD20.root'
,'/store/user/maria22/store/mc/Summer11/LambdaBToPsiMuMu_2MuPEtaFilter_Tight_7TeV-pythia6-evtgen/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/0239034E-18A1-E011-93E3-00226407B98C.root'
,'/store/user/maria22/store/mc/Summer11/LambdaBToPsiMuMu_2MuPEtaFilter_Tight_7TeV-pythia6-evtgen/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/0258FDB4-15A1-E011-9CC2-001B78E2A8C8.root'
,'/store/user/maria22/store/mc/Summer11/LambdaBToPsiMuMu_2MuPEtaFilter_Tight_7TeV-pythia6-evtgen/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/042B7E01-5FA1-E011-BDDF-1CC1DE0437C8.root'
,'/store/user/maria22/store/mc/Summer11/LambdaBToPsiMuMu_2MuPEtaFilter_Tight_7TeV-pythia6-evtgen/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/08392F7F-32A1-E011-A4F1-00237DA1CD7E.root'
,'/store/user/maria22/store/mc/Summer11/LambdaBToPsiMuMu_2MuPEtaFilter_Tight_7TeV-pythia6-evtgen/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/0462F701-54A1-E011-A808-0022649F01AA.root'
,'/store/user/maria22/store/mc/Summer11/LambdaBToPsiMuMu_2MuPEtaFilter_Tight_7TeV-pythia6-evtgen/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/049EBC3A-42A1-E011-91EC-00237DA0F456.root'
,'/store/user/maria22/store/mc/Summer11/LambdaBToPsiMuMu_2MuPEtaFilter_Tight_7TeV-pythia6-evtgen/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/08496A6E-0CA1-E011-B678-0017A4770418.root'
,'/store/user/maria22/store/mc/Summer11/LambdaBToPsiMuMu_2MuPEtaFilter_Tight_7TeV-pythia6-evtgen/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/08EFE1B0-15A1-E011-8DE9-0017A477141C.root'
,'/store/user/maria22/store/mc/Summer11/LambdaBToPsiMuMu_2MuPEtaFilter_Tight_7TeV-pythia6-evtgen/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/0A858D7D-46A1-E011-9275-00237DA1AC24.root'
,'/store/user/maria22/store/mc/Summer11/LambdaBToPsiMuMu_2MuPEtaFilter_Tight_7TeV-pythia6-evtgen/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/0AC12093-43A1-E011-80D9-1CC1DE1CF1BA.root'
,'/store/user/maria22/store/mc/Summer11/LambdaBToPsiMuMu_2MuPEtaFilter_Tight_7TeV-pythia6-evtgen/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/0CAA8A19-45A1-E011-8728-0017A4770420.root'
,'/store/user/maria22/store/mc/Summer11/LambdaBToPsiMuMu_2MuPEtaFilter_Tight_7TeV-pythia6-evtgen/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/0E803F07-43A1-E011-AF96-0017A477002C.root'
,'/store/user/maria22/store/mc/Summer11/LambdaBToPsiMuMu_2MuPEtaFilter_Tight_7TeV-pythia6-evtgen/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/0E426382-32A1-E011-BB9C-001F296AC6F2.root'
,'/store/user/maria22/store/mc/Summer11/LambdaBToPsiMuMu_2MuPEtaFilter_Tight_7TeV-pythia6-evtgen/GEN-SIM-RECO/PU_S4_START42_V11-v1/0000/0EB93FF1-5BA1-E011-8956-78E7D1E4B4C6.root' 
    )
)

process.load('Configuration/StandardSequences/Geometry_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'GR_R_42_V19::All' 

process.Hec2SL  = cms.EDAnalyzer('Hec2sLambda'
   , tracks     = cms.InputTag('generalTracks')
   , MuonsLabel = cms.InputTag('muons')
   , MyPrint    = cms.untracked.bool(False)                  #--cout for short run
)

process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('histo_2slambda.root')
)

process.p = cms.Path(process.Hec2SL)

