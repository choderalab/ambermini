#!/bin/bash

# Test BSC0 parameters for DNA
# Point energy tests are disabled as AmberMini does not build sander
. ../TestCommon.sh

#CleanFiles leap.in leap.out pp.in mdout new.mdout restrt mdinfo DDD.bsc0.parm7 DDD.bsc0.rst7
CleanFiles leap.in leap.out DDD.bsc0.parm7 DDD.bsc0.rst7

printf "\nTest BSC0 parameters.\n\n"
# NOTE: The following logic is included so this test can be used with
#       AmberTools 15. It can (should?) eventually be removed.
if [ -f "$AMBERHOME/dat/leap/cmd/oldff/leaprc.ff14SB" ] ; then
  echo "AmberTools 16"
  LEAPRC="oldff/leaprc.ff14SB"
else
  echo "AmberTools <= 15"
  LEAPRC="leaprc.ff14SB"
fi
cat > leap.in <<EOF
source $LEAPRC
m = loadpdb DDD.pdb
saveamberparm m DDD.bsc0.parm7 DDD.bsc0.rst7
quit
EOF
RunTleap DDD.bsc0.parm7
## Single point energy
#cat > pp.in <<EOF
#single minimization step
#&cntrl
#   imin = 1, ntx = 1, irest = 0, ntwx = 0,
#   ntc = 1, ntf = 1, ntb = 0, cut = 9999.0,
#   igb = 1, ioutfm = 0, ntxo = 1, ntwr = 500,
#&end
#EOF
#$TESTsander -i pp.in -p DDD.bsc0.parm7 -c DDD.bsc0.rst7
$DACDIF DDD.bsc0.parm7.save DDD.bsc0.parm7
$DACDIF DDD.bsc0.rst7.save DDD.bsc0.rst7
#$DACDIF mdout.save mdout

#CleanFiles leap.in leap.out pp.in restrt mdinfo
CleanFiles leap.in leap.out

exit 0
