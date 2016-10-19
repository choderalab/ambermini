#!/bin/bash

# Test different GB radii sets
# Point energy tests are disabled as AmberMini does not build sander

. ../TestCommon.sh

CleanFiles leap.in leap.out *.parm7 1.rst7 *.out 

# NOTE: The following logic is included so this test can be used with
#       AmberTools 15. It can (should?) eventually be removed.
if [ -f "$AMBERHOME/dat/leap/cmd/oldff/leaprc.ff14SB" ] ; then
  echo "AmberTools 16"
  LEAPRC="leaprc.protein.ff14SB"
else
  echo "AmberTools <= 15"
  LEAPRC="leaprc.ff14SB"
fi

printf "\nTest GB radii assignments.\n\n"
for PBRADII in bondi amber6 mbondi mbondi2 parse mbondi3 ; do
  echo "Radii set: $PBRADII"
  cat > leap.in <<EOF
source $LEAPRC
set default pbradii $PBRADII
m = sequence {ACE ALA ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP CTYR}
saveamberparm m $PBRADII.parm7 1.rst7
quit
EOF
  RunTleap $PBRADII.parm7
#  # Single point energy
#  cat > pp.in <<EOF
#single minimization step
#&cntrl
#   imin = 1, ntx = 1, irest = 0, ntwx = 0,
#   ntc = 1, ntf = 1, ntb = 0, cut = 9999.0,
#   igb = 1, ntxo = 1, ioutfm = 0, ntwr = 500,
#&end
#EOF
#  CleanFiles mdinfo restrt
#  $TESTsander -i pp.in -p $PBRADII.parm7 -c 1.rst7 -o $PBRADII.out
  $DACDIF $PBRADII.parm7.save $PBRADII.parm7
#  $DACDIF $PBRADII.out.save $PBRADII.out
done
#CleanFiles leap.in leap.out pp.in restrt mdinfo 1.rst7
CleanFiles leap.in leap.out 1.rst7

exit 0
