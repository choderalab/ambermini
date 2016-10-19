#!/bin/csh -f

echo "*** NMA Genetic Algorithm Perfect Fit ***"

../../../bin/paramfit -i Job_Control.in -p NMA.prmtop -c mdcrd_creation/mdcrd -q mdcrd_creation/amber_energy.dat --random-seed 5000 > prog_out.txt

awk 'BEGIN{FS="\t"}{print $2}' < energy.dat > energy.out
../../dacdif -r 1.e-5 saved_output/energy.out.saved energy.out
../../dacdif saved_output/frcmod.saved frcmod

exit(0)

error:
echo "  ${0}: Program error"
exit(1)

