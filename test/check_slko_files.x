#!/bin/csh -f
#checks if cc.spl file exists - if it does then it
#assumes that all slko files are installed.

if( -r $AMBERHOME/dat/slko/C-C.skf ) then
    #Exists
    exit(0)
else
   echo "DFTB SLKO files not found - Skipping Test..."
   echo ""
   exit(1)
endif


