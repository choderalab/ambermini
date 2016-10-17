#!/bin/bash
## Functions for the GLYCAM Run_tests.bash script that do the main testing
##
## Please see 00_README in this directory for more documentation.
##


##
## The central test function
##
RunParmset() { # This is run for each parmset being tested
  #
  # Make sure the testing environment is good
  #
  CheckAndSetupTestDirectory 
  #
  # If we can't do the tests, bail
  #
  if [ "${AbortTest}" = "y" ] ; then
    NUMTEST=$((NUMTEST+1))
    ERR=$((ERR+1))
    exit 1
  fi

  #
  # Make sure we start clean
  #
  CleanMe
  #
  # If that's all you want, go to the next parameter set
  #
  if  [ "$CLEAN" -eq 1 ]  ; then
    exit 0

  #
  # Announce our intentions
  #
  elif [ "$GENERATE" -eq 1 ] ; then 
  echo "**
Generating new test prmtop and inpcrd built by tleap with GLYCAM_${IAM}.
**" 
  else
  echo "**
Testing prmtop and inpcrd built by tleap with GLYCAM_${IAM}.
**" | tee $TESTOUT

  fi

  #
  # Run each of the tleap input files 
  #
  for i in $(ls *${LEAPINSUFF} 2>/dev/null) ; do
	#
	# Get the filenames for this run
	#
	SetFILENAMES ${i}
	#
	# Run leap on the leap input file
	#
	RunTLEAP ${i}
	#
	# If we aren't generating new .save files, run tests
	#
	if [ "$GENERATE" -eq 0 ] ; then
		#
		# Run dacdif to see if the files are good
		#
		DoTest ${SAVEPRMTOP} ${PRMTOP} 
		DoTest ${SAVEINPCRD} ${INPCRD} 
	#
	# If we are generating new .save files, move the new standards
	#
	else
		mv ${PRMTOP} ${SAVEPRMTOP}
		mv ${INPCRD} ${SAVEINPCRD}
	fi
  done
}

##
## The function that runs dacdif to see if files were generated properly
##
DoTest() {
  LOCALERR=0
  NUMTEST=$((NUMTEST+1))
  #
  # See if dacdif is defined -- if not, bail.
  #
  if [ -z $DACDIF ] ; then
    AbortReport "GLYCAM" "dacdif ($DACDIF) not found."
    LOCALERR=$((LOCALERR+1))
    exit
  fi
  #
  # Check to see if the files to be diff'd are present.
  #
  if [ ! -e $1 ] ; then
    ErrorReport $1 $2 "Standard file $1 not found."
    LOCALERR=$((LOCALERR+1))
    elif [ ! -e $2 ] ; then
      ErrorReport $1 $2 "Test output file $2 not found."
      LOCALERR=$((LOCALERR+1))
  fi
  ERR=$((ERR+LOCALERR))
  #
  # If there are no errors in this test, go ahead with dacdif
  #
  if [ $LOCALERR -eq 0 ] ; then
    #
    # Run dacdif, copying output to a local temporary file
    #
    $DACDIF -k -a 1.6e-1 $1 $2 | tee ${TEMPLOG}
    # Copy results into the main log
    cat $TEMPLOG >> $DIFFOUT
    #
    # See if the temporary file reports failure & act accordingly
    #
    OK=$(grep -c "FAILURE" ${TEMPLOG})
    if [ $OK -ne 0 ] ; then
      FAIL=$((FAIL+1))
    fi
  fi
}
##


##
## This function counts/records how many tests passed, failed or had errors
##
EndTest() {
  if [ $FAIL -gt 0 ] ; then
    echo "  $FAIL out of $NUMTEST comparisons failed." >> $TESTOUT
  fi
  if [ $ERR -gt 0 ]  ; then
    echo "  $ERR out of $NUMTEST tests had errors." >> $TESTOUT
  fi
  PASS=$((NUMTEST-$ERR-$FAIL))
  echo "  $PASS of $NUMTEST comparisons passed." >> $TESTOUT
  echo ""
  #
  # Clean up after, deleting output files (e.g., PRMTOP and INPCRD) if the tests passed
  #
  if [ "$EVALUATE" -eq 0 ] ; then # if the user hasn't asked us to save files
	if [ $PASS -eq $NUMTEST ] ; then
    		CleanMe
	else # delete files that will only confuse things if there are errors
    		CleanError 
	fi
  fi
  #
  # In any case, clean up the dacdif temp files
  #
  CleanDACDIFTempFiles
  # 
  # If a "leap.log" got generated, remove it
  #
  if [ -e $DUMMYLEAPLOG ] ; then
	rm $DUMMYLEAPLOG 
  fi
}
