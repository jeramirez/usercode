#preliminary b-tagging Operating Points
#obtained with  cmssw_3_1_0_pre9
#qcd validation /store/relval/CMSSW_3_1_0_pre9/RelValQCD_Pt_80_120/GEN-SIM-RECO/IDEAL_31X_v1/0007/
#               /store/relval/CMSSW_3_1_0_pre9/RelValQCD_Pt_80_120/GEN-SIM-RECO/IDEAL_31X_v1/0006/
#corrected pt 30 |eta| <2.4 taggability >2
#
import FWCore.ParameterSet.Config as cms

BJetOperatingPointsParameters = cms.PSet(
   BJetOperatingPoints = cms.PSet(
      DefaultBdisc = cms.string('trackCountingHighEffBJetTags'),
      DefaultOp = cms.string('Loose'),
        discCutTight = cms.vdouble(
            14.96,  3.747,         #TCHE, TCHP,
            0.7191, 3.471,         #JTP,  JBTP, 
            3.513,  0.9278,        #SSV,  CSV, 
            0.9349,                #MSV,  
            15.63,  3.187,         #SETByIP3d,SETByPt
            0.3138, 22.83,  1.993  #SMT,  SMTByIP3d,SMTByPt  
        ),
        discCutMedium = cms.vdouble(
            4.575,  2.451,
            0.5042, 2.296,
            2.148,  0.8239,
            0.7785,
            2.403,  1.055,
            0.1373, 0.3761, 0.5866
        ),
        discCutLoose = cms.vdouble(
            1.992,  1.673,
            0.2466, 1.178,
            1.2,    0.3578,
            0.3769, 
           -9.0,    0.0,
            0.0,   -9.0,    0.0
        ),
        bdiscriminators = cms.vstring(
            'trackCountingHighEffBJetTags','trackCountingHighPurBJetTags',
            'jetProbabilityBJetTags','jetBProbabilityBJetTags',
            'simpleSecondaryVertexBJetTags','combinedSecondaryVertexBJetTags',
            'combinedSecondaryVertexMVABJetTags',
            'softElectronByIP3dBJetTags','softElectronByPtBJetTags',
            'softMuonBJetTags','softMuonByIP3dBJetTags','softMuonByPtBJetTags'
        )
   )
)
