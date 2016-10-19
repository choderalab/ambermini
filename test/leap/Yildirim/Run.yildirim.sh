#!/bin/bash

# Test Yildirim et al Chi parameters for RNA.
# Point energy tests are disabled as AmberMini does not build sander

. ../TestCommon.sh

#CleanFiles leap.in leap.out rGACC.yildirim.tip3p.parm7 rGACC.nomin.rst7 \
#           rGACC.Yil.parm7 rGACC.Yil.rst7 pp.in mdout restrt mdinfo tp3.mdout
CleanFiles leap.in leap.out rGACC.yildirim.tip3p.parm7 rGACC.nomin.rst7 \
           rGACC.Yil.parm7 rGACC.Yil.rst7

# Test loading leaprc.RNA.YIL with solvent and ions
# NOTE: The following logic is included so this test can be used with
#       AmberTools 15. It can (should?) eventually be removed.
if [ -f "$AMBERHOME/dat/leap/cmd/oldff/leaprc.ff14SB" ] ; then
  echo "AmberTools 16"
  LEAPRC="leaprc.RNA.YIL"
else
  echo "AmberTools <= 15"
  LEAPRC="oldff/leaprc.parmCHI_YIL.bsc"
fi

printf "\nTest Yildirim et al Chi parameters for RNA.\n\n"
cat > leap.in <<EOF
source $LEAPRC 
set default pbradii mbondi2
# TIP3P water and ions
source leaprc.water.tip3p
m = loadpdb rGACC.pdb
solvateoct m TIP3PBOX 15.54988860294556361621 0.9
addions m Na+ 1
addions m Na+ 1
addions m Na+ 1
saveamberparm m rGACC.yildirim.tip3p.parm7 rGACC.nomin.rst7
quit
EOF
RunTleap rGACC.yildirim.tip3p.parm7
## Single point energy
#cat > pp.in <<EOF
#single minimization step
#&cntrl
#   imin = 1, ntx = 1, irest = 0, ntwx = 0,
#   ntc = 1, ntf = 1, ntb = 1, cut = 8.0,
#   ntxo = 1, ioutfm = 0, ntwr = 500,
#&end
#EOF
#$TESTsander -i pp.in -p rGACC.yildirim.tip3p.parm7 -c rGACC.nomin.rst7 -o tp3.mdout
$DACDIF rGACC.yildirim.tip3p.parm7.save rGACC.yildirim.tip3p.parm7
$DACDIF rGACC.nomin.rst7.save rGACC.nomin.rst7
#$DACDIF tp3.mdout.save tp3.mdout

# Test that loading a protein leaprc does not affect Yildirim params
cat > leap.in <<EOF
source $LEAPRC
source leaprc.protein.ff14SB
m = loadpdb rGACC.pdb
saveamberparm m rGACC.Yil.parm7 rGACC.Yil.rst7
quit
EOF
RunTleap rGACC.Yil.parm7
## Single point energy
#cat > pp.in <<EOF
#single minimization step
#&cntrl
#   imin = 1, ntx = 1, irest = 0, ntwx = 0,
#   ntc = 2, ntf = 2, ntb = 0, cut = 9999.0,
#   igb = 1, ntxo = 1, ioutfm = 0, ntwr = 500,
#&end
#EOF
#rm mdinfo restrt
#$TESTsander -i pp.in -p rGACC.Yil.parm7 -c rGACC.Yil.rst7
$DACDIF rGACC.Yil.parm7.save rGACC.Yil.parm7
$DACDIF rGACC.Yil.rst7.save rGACC.Yil.rst7
#$DACDIF mdout.save mdout

#CleanFiles leap.in leap.out pp.in restrt mdinfo
CleanFiles leap.in leap.out

exit 0
