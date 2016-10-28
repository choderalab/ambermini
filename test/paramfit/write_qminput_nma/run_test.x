#!/bin/csh -f

echo "*** Create Gaussian Input Files (new prmtop) ***"

../../../bin/paramfit -i Job_Control.in -p mdcrd_calc/NMA.prmtop -c mdcrd_calc/mdcrd > prog_out.txt || goto error

#../../dacdif saved_output/prog_out.saved prog_out.txt
../../dacdif saved_output/Job.0.gjf.saved Job.0.gjf
../../dacdif saved_output/Job.1.gjf.saved Job.1.gjf
../../dacdif saved_output/Job.2.gjf.saved Job.2.gjf
../../dacdif saved_output/Job.3.gjf.saved Job.3.gjf
../../dacdif saved_output/Job.4.gjf.saved Job.4.gjf

exit(0)

error:
echo "  ${0}: Program error"
exit(1)

