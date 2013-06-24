#! /bin/csh

setenv SCRAM_ARCH slc5_amd64_gcc434
setenv MYCMSPYTHON kinkntupla_cfg.py
setenv MYRUN 2011A
setenv MYCMSHOME /home/eduardo/CMSSW/CMSSW_4_2_4/src
setenv MYCMSTESTDIR Analyzers/KinkProducer/test
#input dir
setenv MYHADOOPPATH /mnt/hadoop/data
setenv MYINPUTPATH /store/data/Run${MYRUN}/MuOnia/RECO/PromptReco-v1/000/176/886
#file output name
setenv MYPREFIXOUTPUT T3kink${MYRUN}
setenv MYOUTPUTPATH /home/eduardo/CMSSW/CMSSW_4_2_4/src/Analyzers/KinkProducer/test

source /sharesoft/osg/app/cmssoft/cms/cmsset_default.csh

cd $MYCMSHOME 
cmsenv
cd $MYCMSTESTDIR 

#input parameters
#first file number sorted by creation time
#last file number sorted by creation time
#input path will replace test
if ($#argv < 1 ) then
  echo "arguments missing example: $0 1 1"
endif
set myfirst  = `expr $1`
set mylast   = `expr $2`
if ($#argv > 2 ) then
  setenv MYINPUTPATH $3
endif

#BIG LOOP script
set jobid  = 0;
foreach file ( `ls -tr $MYHADOOPPATH/$MYINPUTPATH`)
 setenv MYINPUTFILE $MYINPUTPATH/$file
 set jobid  = `expr $jobid + 1`
 if ( $jobid >= $myfirst && $jobid <= $mylast ) then
    if ( -e $MYHADOOPPATH/$MYINPUTFILE ) then
      echo "file ${MYHADOOPPATH}$MYINPUTFILE"
      setenv MYOUTPUT ${MYPREFIXOUTPUT}_${jobid}_1.root
      echo "starting cmsRun ..."
      if ($?TMPDIR) then
         echo "environment ready, now mv to $TMPDIR"
         rm -f $TMPDIR/$MYCMSPYTHON
         cp $MYCMSPYTHON $TMPDIR/$MYCMSPYTHON
         cd $TMPDIR
      endif
#      cmsRun $MYCMSPYTHON print files=file:${MYHADOOPPATH}$MYINPUTFILE output=$MYOUTPUT maxEvents=-1 >& _condor_stdout

      cmsRun $MYCMSPYTHON print files=file:${MYHADOOPPATH}$MYINPUTFILE output=$MYOUTPUT maxEvents=-1
       if ($?TMPDIR) then
         echo "now ... copying output to $MYOUTPUTPATH"
         rm -f $MYOUTPUTPATH/$MYOUTPUT
         cp -f $MYOUTPUT $MYOUTPUTPATH/$MYOUTPUT
         echo "output size for $MYOUTPUT"
         ls -l $MYOUTPUT
         rm -f $MYOUTPUT
         echo "now ... removing tmp copy of  $MYCMSPYTHON"
         rm -f $TMPDIR/$MYCMSPYTHON
       endif

    endif
 endif
end

echo "exiting job..."
