#!/bin/bash

# Test parameters for Protein/DNA/RNA system
# Point energy tests are disabled as AmberMini does not build sander
. ../TestCommon.sh

CleanFiles leap.in leap.out 4py5.parm7 new.4py5.parm7 4py5.rst7 pp.in \
           mdout new.mdout restrt mdinfo bsc1.mdout bsc1.4py5.parm7

printf "\nTest parameters for protein/DNA/RNA system.\n\n"
# NOTE: The following logic is included so this test can be used with
#       AmberTools 15. It can (should?) eventually be removed.
if [ -f "$AMBERHOME/dat/leap/cmd/oldff/leaprc.ff14SB" ] ; then
  echo "AmberTools 16"
  LEAPRC="oldff/leaprc.ff14SB"
else
  echo "AmberTools <= 15"
  LEAPRC="leaprc.ff14SB"
fi

# Input for single point energy
#cat > pp.in <<EOF
#single minimization step
#&cntrl
#   imin = 1, ntx = 1, irest = 0, ntwx = 0,
#   ntc = 2, ntf = 2, ntb = 0, cut = 9999.0,
#   igb = 5, ioutfm = 0, ntxo = 1, ntwr = 500,
#&end
#EOF

# Use the old ff14SB containing BSC0/OL3/14SB
TestOldFF14SB() {
  cat > leap.in <<EOF
source $LEAPRC
set default pbradii mbondi2
m = loadpdb Amber.4py5.mod.pdb
saveamberparm m 4py5.parm7 4py5.rst7
quit
EOF
  RunTleap 4py5.parm7

#  $TESTsander -i pp.in -p 4py5.parm7 -c 4py5.rst7
  $DACDIF 4py5.parm7.save 4py5.parm7
  $DACDIF 4py5.rst7.save 4py5.rst7
#  $DACDIF mdout.save mdout
}

# Use the new BSC0/OL3/14SB separate leaprc files
TestOldBSC0() {
if [ -f "$AMBERHOME/dat/leap/cmd/oldff/leaprc.DNA.bsc0" ] ; then
  cat > leap.in <<EOF
source leaprc.protein.ff14SB
source oldff/leaprc.DNA.bsc0
source leaprc.RNA.OL3
set default pbradii mbondi2
m = loadpdb Amber.4py5.mod.pdb
saveamberparm m new.4py5.parm7 4py5.rst7
quit
EOF
  RunTleap new.4py5.parm7
  CleanFiles mdinfo restrt
#  $TESTsander -i pp.in -p new.4py5.parm7 -c 4py5.rst7 -o new.mdout
  $DACDIF 4py5.parm7.save new.4py5.parm7
  $DACDIF 4py5.rst7.save 4py5.rst7
#  $DACDIF mdout.save new.mdout
fi
}

#TestOldFF14SB
#TestOldBSC0

# Use the new BSC1/OL3/14SB separate leaprc files
cat > leap.in <<EOF
source leaprc.protein.ff14SB
source leaprc.DNA.bsc1
source leaprc.RNA.OL3
m = loadpdb Amber.4py5.mod.pdb
saveamberparm m bsc1.4py5.parm7 4py5.rst7
quit
EOF
RunTleap bsc1.4py5.parm7
CleanFiles mdinfo restrt
#$TESTsander -i pp.in -p bsc1.4py5.parm7 -c 4py5.rst7 -o bsc1.mdout
$DACDIF bsc1.4py5.parm7.save bsc1.4py5.parm7
$DACDIF 4py5.rst7.save 4py5.rst7
#$DACDIF bsc1.mdout.save bsc1.mdout

#CleanFiles leap.in leap.out pp.in restrt mdinfo
CleanFiles leap.in leap.out

exit 0
