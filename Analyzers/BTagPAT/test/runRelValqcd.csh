#!/bin/csh
cmsRun patLayer1_fromAOD_full_ttbar_cfg.py print\
 output='PATLayer1_Output.fromAOD_qcd_full.root' \
 files='/store/relval/CMSSW_3_1_0_pre9/RelValQCD_Pt_80_120/GEN-SIM-RECO/IDEAL_31X_v1/0007/589CA925-534F-DE11-AF9E-001D09F2841C.root' \
 files_load='RelVal_qcd_80_120.txt' \
 maxEvents=-1\
 mypatsequence='Analyzers.BTagPAT.patSequences_cff'
#list of input file:
#'/store/relval/CMSSW_3_1_0_pre9/RelValQCD_Pt_80_120/GEN-SIM-RECO/IDEAL_31X_v1/0006/A8CCCCA8-5F4E-DE11-B177-001617C3B70E.root' \
#'/store/relval/CMSSW_3_1_0_pre9/RelValQCD_Pt_80_120/GEN-SIM-RECO/IDEAL_31X_v1/0006/7EE39777-604E-DE11-9605-001D09F2A465.root',\
#'/store/relval/CMSSW_3_1_0_pre9/RelValQCD_Pt_80_120/GEN-SIM-RECO/IDEAL_31X_v1/0006/521E943E-614E-DE11-9E92-001D09F29169.root',\
#'/store/relval/CMSSW_3_1_0_pre9/RelValQCD_Pt_80_120/GEN-SIM-RECO/IDEAL_31X_v1/0006/32847043-624E-DE11-BBBB-001D09F295A1.root',\
#'/store/relval/CMSSW_3_1_0_pre9/RelValQCD_Pt_80_120/GEN-SIM-RECO/IDEAL_31X_v1/0006/0268E2E0-5D4E-DE11-90B3-001D09F253D4.root'\

       
