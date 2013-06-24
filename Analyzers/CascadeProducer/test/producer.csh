#! /bin/csh

setenv SCRAM_ARCH slc5_ia32_gcc434
#++condor
if ($?PATH)then
  echo "PATH is $PATH"
else
  setenv HOME /usr/users/$USER
  setenv PATH "/bin:/usr/bin:/usr/local/bin:/usr/X11R6/bin"
endif
#--condor
source ~/cms_sl5.csh
#cd /uscms_data/d2/jervar/CMSSW/CMSSW_3_5_8/src 
cd /hep/alpha01/home/eduardo/CMSSW/CMSSW_3_5_8/src 
cmsenv
cd Analyzers/CascadeProducer/test 

setenv MYINPUTFILE file:/hep/alpha01/home/CMSDATA/vee358_2010/RAW-RECO/Vmasscut/triggerandlooseveeskimmed0${1}_${2}.root
setenv MYOUTPUT cascproducerlooseveeskimmall2010_0${1}_${2}.root
cmsRun cascadeproducer_cfg.py print files=$MYINPUTFILE output=$MYOUTPUT maxEvents=-1
