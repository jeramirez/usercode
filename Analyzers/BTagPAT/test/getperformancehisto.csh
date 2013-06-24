#!/bin/csh
#
cmsRun BTagPATAnalyzer_cfg.py print maxEvents=-1 \
output2=selectedLayer1_qcd.root \
patlayer1=selectedLayer1Jets \
files=file:PATLayer1_Output.fromAOD_qcd_full.root 
