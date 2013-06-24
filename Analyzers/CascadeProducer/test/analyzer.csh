#! /bin/csh

setenv SCRAM_ARCH slc5_ia32_gcc434
setenv MYCMSPYTHON cascadeanalyzer_cfg.py
source /sharesoft/osg/app/cmssoft/cms/cmsset_default.csh

cd /home/eduardo/CMSSW/CMSSW_3_9_7/src 
cmsenv
cd Analyzers/CascadeProducer/test 

setenv MYINPUTPATH /user/eduardo/reprocess/cascade/397
setenv MYOUTPUTPATH /mnt/hadoop/user/eduardo/cascade/analysis/397
setenv MYINPUTFILE $MYINPUTPATH/T3cascproducerlooseveeskimmall2010_0${1}_${2}.root
setenv MYOUTPUT T3casanalyzerlooseskim2010_0${1}_${2}.root

if ($?TMPDIR) then
   echo "environment ready, now mv to $TMPDIR"
   rm -f $TMPDIR/$MYCMSPYTHON
   cp $MYCMSPYTHON $TMPDIR/$MYCMSPYTHON
   cd $TMPDIR
endif

cmsRun $MYCMSPYTHON print files=$MYINPUTFILE output=$MYOUTPUT maxEvents=-1

echo "cmsRun finished ..."
echo "now ... copying output to $MYOUTPUTPATH" 
rm -f $MYOUTPUTPATH/$MYOUTPUT 
cp -f $MYOUTPUT $MYOUTPUTPATH/$MYOUTPUT 
 
