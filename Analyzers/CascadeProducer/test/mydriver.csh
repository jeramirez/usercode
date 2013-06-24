#!/bin/csh
#cmsDriver.py SingleElectronFlatPt5To100.cfi -s GEN:ProductionFilterSequence,SIM,DIGI,L1,DIGI2RAW,RAW2DIGI,RECO -n 10 --conditions FrontierConditions_GlobalTag,MC_31X_V2::All --datatier 'GEN-SIM-RAW-RECO' --eventcontent FEVT --no_exec
#cmsDriver.py PythiaH190ZZ4mu_cfi.py -s GEN:ProductionFilterSequence --conditions FrontierConditions_GlobalTag,IDEAL_31X::All --datatier 'GEN-SIM-RAW' --eventcontent RAWSIM -n 1000 --no_exec 
#cmsDriver.py Configuration/Generator/python/JpsiMM_Pt_20_inf_cfi.py -s GEN:ProductionFilterSequence,FASTSIM \
#cmsDriver.py Configuration/Generator/python/bJpsiX_cfi.py -s GEN:ProductionFilterSequence,FASTSIM \
#cmsDriver.py Analyzers/CascadeProducer/test/myparticlegun_cff.py -s GEN:ProductionFilterSequence,FASTSIM \
#cmsDriver.py Configuration/Generator/python/bJpsiX_EXTRAS_cff.py -s GEN:ProductionFilterSequence,FASTSIM \
#cmsDriver.py SingleElectronFlatPt5To100.cfi -s GEN:ProductionFilterSequence,SIM,DIGI,L1,DIGI2RAW,RAW2DIGI,RECO \
cmsDriver.py SingleElectronFlatPt5To100.cfi -s GEN:ProductionFilterSequence,SIM,DIGI,L1,DIGI2RAW,RAW2DIGI,RECO \
--conditions FrontierConditions_GlobalTag,STARTUP_31X::All \
--datatier GEN-SIM-DIGI-RECO \
--eventcontent FEVT \
-n 10 --no_exec
#--pileup=NoPileUp \
#--beamspot=Early10TeVCollision \
#--eventcontent AODSIM \

