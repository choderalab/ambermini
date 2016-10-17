#!/bin/bash
## Setup information for the GLYCAM Run_tests.bash script
##
## Please see 00_README in this directory for more documentation.
##
##########
##
# These file types will be instrumental in the tests
##
# This file is particularly important.  
# Tests will be performed for any files having this
#	string as a suffix.	
LEAPINSUFF=".leapin"
# The rest of these files will be written using the
#	prefix found for each file matching LEAPINSUFF
TOPSUFF=".parm7"
CRDSUFF=".rst7"
LEAPOUTSUFF=".tleap_out"
LEAPERRSUFF=".tleap_error"
DIFFOUTSUFF=".dacdif_out"
TEMPLOGSUFF=".localLog"
# This suffix is not expected to change as it is the standard
#	comparison suffix for AMBER tests
SAVESUFF=".save"
# The temp files DACDIF writes will have this prefix
DACTEMPPREF="ddtmp."
# Overall info will go here
TESTOUT="Test_Results"
# This file sometimes gets made
DUMMYLEAPLOG="leap.log"
##########
# These files can be removed after successful tests
CLEANUPSUFFIXES="${TESTSUFF} ${LEAPERRSUFF} ${LEAPOUTSUFF} ${DIFFOUTSUFF} ${TEMPLOGSUFF} ${TOPSUFF} ${CRDSUFF}"
# These files should be removed even after unsuccessful tests
ErrorCLEANUPSUFFIXES="${TEMPLOGSUFF}"
##


##
## Check if we have a particular task to do 
##
CLEAN=0
EVALUATE=0
GENERATE=0
if [ "$1" = "clean" ] ; then
  CLEAN=1
elif [ "$1" = "evaluate" ] ; then
  EVALUATE=1
elif [ "$1" = "generate" ] ; then
  GENERATE=1
elif [ $1 ] ; then
  echo "Unknown run option \"${1}\".  Exiting."
  exit 1
fi
if [ $2 ] ; then
  echo "Too many arguments on command line.  Exiting."
  exit 1
fi

##


##
## If this run is not just for cleaning, then perform some basic setup
##
if [ $CLEAN -eq  0 ] ; then  # Do some other setup for the tests
  if [ -z $AMBERHOME ] ; then
    echo "Error: The GLYCAM tleap tests require AMBERHOME to be defined."
    exit 1
  fi
  DACDIF="$AMBERHOME/AmberTools/test/dacdif"
  TLEAP="$AMBERHOME/bin/tleap"
  if [ ! -x $DACDIF ] ; then
    echo "$DACDIF not found." 
    exit 1
  fi
  if [ ! -x $TLEAP ] ; then
    echo "$TLEAP not found."
    exit 1
  fi
fi
