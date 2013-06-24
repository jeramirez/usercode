import FWCore.ParameterSet.Config as cms

process = cms.Process("DemoHec")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("RecoEgamma.PhotonIdentification.photonId_cff")
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi")
process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")
####process.load("MagneticField.Engine.volumeBasedMagneticField_cfi")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Geometry.CaloEventSetup.CaloTopology_cfi")
process.load("RecoEgamma.EgammaPhotonProducers.photonSequence_cff")
process.load("RecoEcal.EgammaClusterProducers.islandBasicClusters_cfi")
process.load('Configuration.StandardSequences.Reconstruction_cff')


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10))
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.MessageLogger.destinations = cms.untracked.vstring("hola_ces.log")#

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
   #'file:/uscms_data/d2/pollack/jobs/cmssw_3_6_3/phase1R3916F/fast/phase1R3916F_fnal_02.root'  #--No Rad
   #'file:/uscms_data/d2/pollack/jobs/cmssw_3_6_3/phase1R3916F/fast/phase1R3916F_fnal_03.root'  #--Rad
   #'file:fasttst_ee_00.root'
   'file:fasttst_ee_02.root'
   #'file:fasttst_ee_01.root'
   #-- 'file:fasttst_mm_00.root' 
   #--no exist  '/store/mc/Winter09/Zee_M20/GEN-SIM-DIGI-RECO/IDEAL_V12_FastSim_v1/0000/FC9369D3-1E16-DE11-864C-001EC9B48D0A.root',
   #--Full sim low E cluster '/store/mc/Spring10/Zee_M20_CTEQ66-powheg/GEN-SIM-RECO/START3X_V26-v2/0022/D0DB5AF9-6B62-DF11-A003-001EC94BA37C.root',
   #'/store/mc/Winter09/ZeeJet_Pt80to120/GEN-SIM-DIGI-RECO/IDEAL_V11_FastSim_v1/0000/F8322527-75ED-DD11-8357-0019B9E4FD9D.root',
#
   #'file:/uscms_data/d2/mendez/../pollack/jobs/cmssw_3_6_3/phase1R3916F/fast/phase1R3916F_fnal_02.root'
   #'file:/uscmst1b_scratch/lpc1/lpceg/carley/JetMET/Photest.root'
   ##'/store/data/Commissioning10/MinimumBias/RECO/v8/000/132/528/867BB35A-2B3E-DF11-AC1D-000423D99996.root',
#
   #-->'file:FEADF2F9-D25D-DF11-91B3-002618943935.root'
   #'file:/afs/cern.ch/cms/Tutorials/TWIKI_DATA/RelValTTbar_RECO_331.root'
#     '/store/relval/CMSSW_3_6_2/RelValMinBias/ALCARECO/START36_V10_ALCARECOEcalCalPhiSym-v1/0003/C2A91967-4371-DF11-83BB-001A92971BDA.root'
#
#   '/store/group/muon/ppMuX/Spring10ReDigi_skimJPsiLoose_v2/3f0fb5c06e0ff4a2dfc570422e29b244/skimJPsiLoose_Spring10ReDigi_9_2.root'
#
   # 'file:/uscms_data/d2/mendez/../pollack/jobs/cmssw_3_6_3/phase1R3916F/fast/phase1R3916F_fnal_02.root'
   # 'file:/uscms_data/d2/mendez/../pollack/jobs/cmssw_3_3_6/phase1R3916F/fast/phase1R3916F_fnal_30.root'
   # 'file:/uscms_data/d2/mendez/../pollack/jobs/cmssw_3_3_6/phase1R3916F/fast/phase1R3916F_fnal_29.root'
# 
#
   ####'/store/data/Commissioning10/MinimumBias/RECO/May6thPDSkim2_SD_EGMonitor-v1/0134/FE0AA62F-975D-DF11-A2F1-002618943831.root',
#   '/store/data/Commissioning10/MinimumBias/RAW-RECO/v8/000/133/532/EC93873A-D74B-DF11-A1B9-00E08179185D.root',
#   '/store/data/Commissioning10/MinimumBias/RAW-RECO/v8/000/133/531/D6E1CE68-ED4B-DF11-A676-003048D45F84.root',
#   '/store/data/Commissioning10/MinimumBias/RAW-RECO/v8/000/133/529/223C34BD-EC4B-DF11-9CA6-003048D476D4.root',
#   '/store/data/Commissioning10/MinimumBias/RAW-RECO/v8/000/133/526/B28AEED6-E94B-DF11-9124-00E08178C103.root',
#   '/store/data/Commissioning10/MinimumBias/RAW-RECO/v8/000/133/521/C0EFDC20-024C-DF11-A82A-00E08178C155.root',
#   '/store/data/Commissioning10/MinimumBias/RAW-RECO/v8/000/133/518/C620697E-B94B-DF11-A125-003048D46090.root'
     #--'/store/relval/2008/6/22/RelVal-RelValSingleGammaPt35-1213986417-IDEAL_V2-2nd/0004/443BCAED-CB40-DD11-AB37-000423D6B48C.root'
#     
    #'/store/data/Commissioning10/MinimumBias/RECO/May6thPDSkim2_SD_JetMETTau-v1/0137/FEADF2F9-D25D-DF11-91B3-002618943935.root',
#     
   #'dcache:/pnfs/cms/WAX/11/store/user/askew/minbias/EJTermFile/Photest003.root',
   #'dcache:/pnfs/cms/WAX/11/store/user/askew/minbias/EJTermFile/Photest001.root',
   # 'dcache:/pnfs/cms/WAX/11/store/user/askew/minbias/EJTermFile/Photest000.root',
   # 'dcache:/pnfs/cms/WAX/11/store/user/askew/minbias/EJTermFile/Photest004.root',
    )
)

process.demo = cms.EDAnalyzer('HecMesons'
  #, tracks = cms.untracked.InputTag('ctfWithMaterialTracks')
   , tracks = cms.untracked.InputTag('generalTracks')
   , gsfElectrons = cms.untracked.InputTag('gsfElectrons')
   , minTracks = cms.untracked.uint32(0)
   , check_elC = cms.untracked.bool(True)                #--Dump e- matching info
   , Pe_match  = cms.untracked.double(10.0)
 # , mass_uno  = cms.untracked.double(0.105)             #--Muons
 # , mass_uno  = cms.untracked.double(0.493667)          #--Kaons
   , mass_uno  = cms.untracked.double(0.00051099891)     #--electrons
)

process.TFileService = cms.Service("TFileService",
    #fileName = cms.string('histo_ee_cesar_phase1_02_norad.root')
    #fileName = cms.string('histo_ee_cesar_phase1_03_rad.root')
    #fileName = cms.string('histo_ee_cesar02.root')
    fileName = cms.string('histo_ee_00.root')
    #fileName = cms.string('histo_ee_01.root')
    #fileName = cms.string('histo_mm_00.root')
)
process.dump=cms.EDAnalyzer('EventContentAnalyzer')
#process.p = cms.Path(process.demo*process.dump)
process.p = cms.Path(process.islandBasicClusters*process.demo)
#process.p = cms.Path(process.demo)
