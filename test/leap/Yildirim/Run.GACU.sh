#!/bin/bash

# Test Yildirim et al Chi parameters for RNA.
# Point energy tests are disabled as AmberMini does not build sander

. ../TestCommon.sh

#CleanFiles leap.in leap.out pp.in rGACU.out restrt mdinfo rGACU.parm7 rGACU.rst7
CleanFiles leap.in leap.out rGACU.parm7 rGACU.rst7

# NOTE: The following logic is included so this test can be used with
#       AmberTools 15. It can (should?) eventually be removed.
if [ -f "$AMBERHOME/dat/leap/cmd/oldff/leaprc.ff14SB" ] ; then
  echo "AmberTools 16"
  LEAPRC="leaprc.RNA.YIL"
else
  echo "AmberTools <= 15"
  LEAPRC="oldff/leaprc.parmCHI_YIL.bsc"
fi

printf "\nTest Yildirim et al Chi parameters for RNA r(GACU).\n\n"
cat > leap.in <<EOF
source $LEAPRC
m = sequence { G A C U } 
saveamberparm m rGACU.parm7 rGACU.rst7
quit
EOF
RunTleap rGACU.parm7
# Single point energy
#cat > pp.in <<EOF
#single minimization step
#&cntrl
#   imin = 1, ntx = 1, irest = 0, ntwx = 0,
#   ntc = 1, ntf = 1, ntb = 0, cut = 9999.0,
#   igb = 1, ntxo = 1, ioutfm = 0, ntwr = 500,
#&end
#EOF
#$TESTsander -i pp.in -p rGACU.parm7 -c rGACU.rst7 -o rGACU.out
$DACDIF rGACU.parm7.save rGACU.parm7
$DACDIF rGACU.rst7.save rGACU.rst7
#$DACDIF rGACU.out.save rGACU.out

#CleanFiles leap.in leap.out pp.in restrt mdinfo
CleanFiles leap.in leap.out

exit 0
