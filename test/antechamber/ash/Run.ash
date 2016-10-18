#!/bin/csh -f


../../../bin/antechamber -i ash.pdb -fi pdb -o ash.mol2 -fo mol2 -c bcc \
     >& antechamber.out || goto error
../../../bin/parmchk -i ash.mol2 -f mol2 -o frcmod || goto error
#  ../../../bin/tleap -s -f leap.in > leap.out || goto error

../../dacdif -a 1.5e-3 ash.mol2.save ash.mol2
../../dacdif frcmod.save frcmod

/bin/rm -f ANTE* ATOMTYPE.INF BCCTYPE.INF FOR* NEWPDB.PDB PREP.INF \
      leap.log prmcrd divcon.dmx divcon.rst divcon.in

exit(0)

error:
echo "  ${0}:  Program error"
exit(1)
