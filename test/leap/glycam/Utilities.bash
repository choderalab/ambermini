#!/bin/bash
## Utility functions for the GLYCAM Run_tests.bash script
##
## Please see 00_README in this directory for more documentation.
##


##
## A little function to run tleap
##
RunTLEAP() {
  $TLEAP -f $1 >> $LEAPOUT 2>>$LEAPERR
}
##


##
## Set up the current filenames to use
##
##	Usage:  SetFILENAMES(filename.leapinsuffix)
##
SetFILENAMES() {
  FILEBASE=${1%${LEAPINSUFF}}
  PRMTOP=${FILEBASE}${TOPSUFF}
  INPCRD=${FILEBASE}${CRDSUFF}
  LEAPOUT=${FILEBASE}${LEAPOUTSUFF}
  LEAPERR=${FILEBASE}${LEAPERRSUFF}
  DIFFOUT=${FILEBASE}${DIFFOUTSUFF}
  DIFFERR=${FILEBASE}${DIFFERRSUFF}
  TEMPLOG=${FILEBASE}${TEMPLOGSUFF} 
  SAVEPRMTOP=${FILEBASE}${TOPSUFF}${SAVESUFF}
  SAVEINPCRD=${FILEBASE}${CRDSUFF}${SAVESUFF}
}
##


##
## Functions to facilitate error reporting
##
ErrorReport() {
 echo "During the test of files ${1} and ${2}, there was an
Error:  ${3}
==============================================================" | tee -a $TESTOUT
}
AbortReport() {
 echo "The ${1} test must be aborted because:
Error:  ${2}
==============================================================" | tee -a $TESTOUT
}
##


##
## Functions for removing things
##
CleanMe() { # Cleans a parameter set directory
  CleanSuffixes "${CLEANUPSUFFIXES}"
  if [ -e $DUMMYLEAPLOG ] ; then 
	rm $DUMMYLEAPLOG
  fi
  if [ -e $TESTOUT ] ; then 
  	rm $TESTOUT
  fi
}
CleanError() { # Cleanup if there was an error
  CleanSuffixes "${ErrorCLEANUPSUFFIXES}"
}
CleanSuffixes() { # remove all files with the approved suffixes
  for i in ${1} ; do 
	files="$(ls *${i} 2> /dev/null | wc -l)"
	if [ "${files}" != "0" ] ; then
		for j in $(ls *${i} 2> /dev/null) ; do
       			rm $j
  		done
  	fi
  done
}
CleanDACDIFTempFiles() {
  files="$(ls ${DACTEMPPREF}* 2> /dev/null | wc -l)"
  if [ "${files}" != "0" ] ; then
	for j in $(ls ${DACTEMPPREF}* 2>/dev/null) ; do
		rm $j
	done
  	fi
}
##


##
## Check test directories and contents. 
##
CheckAndSetupTestDirectory() {
  #
  # Make sure that the directory contains at least one .leapin file
  #
  numLEAPIN=0
  if [ "${CLEAN}" = 0 ] ; then
	files=$(ls *$LEAPINSUFF 2> /dev/null | wc -l)
	if [ "${files}" != "0" ] ; then
    		for chk in $(ls *$LEAPINSUFF 2>/dev/null) ; do
			numLEAPIN=$((numLEAPIN+1))
    		done
	fi
    if [ "${numLEAPIN}" -eq 0 ] ; then
      AbortReport ${IAM} "There are no leap input files to test."
      AbortTest="y"
    fi
    touch dum
    if [ ! -f dum ] ; then
	AbortReport ${IAM} "Cannot write to test directory."
	AbortTest="y"
    fi
    rm dum
  fi
}

