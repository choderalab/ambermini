#!/bin/bash

# Test OL3 Chi parameters for RNA.
# Point energy tests are disabled as AmberMini does not build sander

. ../TestCommon.sh

#CleanFiles leap.in leap.out rGACC.tip3p.parm7 rGACC.nomin.rst7 pp.in mdout restrt mdinfo \
#           tp3.mdout
CleanFiles leap.in leap.out rGACC.tip3p.parm7 rGACC.nomin.rst7

# Test loading leaprc.RNA.OL3 with solvent and ions
# NOTE: The following logic is included so this test can be used with
#       AmberTools 15. It can (should?) eventually be removed.
if [ -f "$AMBERHOME/dat/leap/cmd/oldff/leaprc.ff14SB" ] ; then
  echo "AmberTools 16"
  LEAPRC="leaprc.RNA.OL3"
else
  echo "AmberTools <= 15"
  LEAPRC="leaprc.ff14SB"
fi

printf "\nTest OL3 Chi parameters for RNA.\n\n"
cat > leap.in <<EOF
source $LEAPRC
source leaprc.water.tip3p
# TIP3P ions
loadamberparams frcmod.ionsjc_tip3p
m = loadpdb ../Yildirim/rGACC.pdb
solvateoct m TIP3PBOX 15.54988860294556361621 0.9
addions m Na+ 1
addions m Na+ 1
addions m Na+ 1
saveamberparm m rGACC.tip3p.parm7 rGACC.nomin.rst7
quit
EOF
RunTleap rGACC.tip3p.parm7
# Single point energy
#cat > pp.in <<EOF
#single minimization step
#&cntrl
#   imin = 1, ntx = 1, irest = 0, ntwx = 0,
#   ntc = 1, ntf = 1, ntb = 1, cut = 8.0,
#   ntxo = 1, ioutfm = 0, ntwr = 500,
#&end
#EOF
#$TESTsander -i pp.in -p rGACC.tip3p.parm7 -c rGACC.nomin.rst7 -o tp3.mdout
$DACDIF rGACC.tip3p.parm7.save rGACC.tip3p.parm7
$DACDIF rGACC.nomin.rst7.save rGACC.nomin.rst7
#$DACDIF tp3.mdout.save tp3.mdout

# Test that loading a protein leaprc does not affect RNA params
cat > leap.in <<EOF
source $LEAPRC
source leaprc.protein.ff14SB
m = loadpdb ../Yildirim/rGACC.pdb
saveamberparm m rGACC.OL3.parm7 rGACC.OL3.rst7
quit
EOF
RunTleap rGACC.OL3.parm7
# Single point energy
#cat > pp.in <<EOF
#single minimization step
#&cntrl
#   imin = 1, ntx = 1, irest = 0, ntwx = 0,
#   ntc = 2, ntf = 2, ntb = 0, cut = 9999.0,
#   igb = 1, ntxo = 1, ioutfm = 0, ntwr = 500,
#&end
#EOF
#rm mdinfo restrt
#$TESTsander -i pp.in -p rGACC.OL3.parm7 -c rGACC.OL3.rst7
$DACDIF rGACC.OL3.parm7.save rGACC.OL3.parm7
$DACDIF rGACC.OL3.rst7.save rGACC.OL3.rst7
#$DACDIF mdout.save mdout

#CleanFiles pp.in restrt mdinfo leap.in leap.out
CleanFiles leap.in leap.out

exit 0
