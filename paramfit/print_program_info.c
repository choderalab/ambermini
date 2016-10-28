/*****************************************************
 * AMBER Bond Angle and Dihedral Parameter Optimiser *
 *                                                   *
 *           Written by: Robin Betz  (2011)          *
 *                       Ross Walker (2004)          *
 *                   UC San Diego                    *
 *           San Diego Supercomputer Center          *
 *            La Jolla, California, 92092            *
 *                       USA                         *
 *****************************************************/

/*print_program_info.c*/

/*

This module prints information about the program - either
the program header info or the development history

void print_program_info()
void print_program_history()
  
*/

#include <stdio.h>

#include "function_def.h"

void print_program_info(void)
{
printf("\n");
printf("                *****************************************************\n");
printf("                * AMBER Bond Angle and Dihedral Parameter Optimiser *\n");
printf("                *                                                   *\n");
printf("                *                      v3.0.0                       *\n");
printf("                *                                                   *\n");
printf("                *                    Written by:                    *\n");
printf("                *                 Robin Betz (2011)                 *\n");
printf("                *                 Ross Walker (2004)                *\n");
printf("                *          The Walker Molecular Dynamics Lab        *\n");
printf("                *         University of California, San Diego       *\n");
printf("                *            La Jolla, California, 92092            *\n");
printf("                *                       USA                         *\n");
printf("                *****************************************************\n");
printf("\n");
fflush(stdout); /*Flush the printf buffer*/
}

void print_program_history(void)
{
printf("\n");
printf("        -------------------------- PROGRAM HISTORY ---------------------------\n");
printf("        |                            date order                              |\n");
printf("        |                                                                    |\n");
printf("        |                     CURRENT VERSION = v3.0.0                       |\n");
printf("        |                          28 / Feb / 2012                           |\n");
printf("        |                                                                    |\n");
printf("        | Version 0.1.0 Alpha 1                                              |\n");
printf("        | ---------------------                                              |\n");
printf("        | 1) Initial development of program modules.                         |\n");
printf("        | 2) Initial creation of routines for reading job control and data   |\n");
printf("        |    files.                                                          |\n");
printf("        | 3) Initial implementation of the Simplex algorithm for             |\n");
printf("        |    minimisation.                                                   |\n");
printf("        |                                                                    |\n");
printf("        | Version 0.1.1 Alpha 1                                              |\n");
printf("        | ---------------------                                              |\n");
printf("        | 1) Completed implementation of simplex minimiser.                  |\n");
printf("        |    Minimisation is available for BONDS, ANGLES and DIHEDRALS.      |\n");
printf("        | 2) Completed implementation of sum of squares minimisation for     |\n");
printf("        |    standard amber force field equation.                            |\n");
printf("        | 3) Added extra options to job control file:                        |\n");
printf("        |    CONV_LIMIT, BONDFC_dx, BONDEQ_dx, ANGLEFC_dx, ANGLEEQ_dx,       |\n");
printf("        |    DIHEDRALBH_dx, DIHEDRALN_dx and DIHEDRALG_dx                    |\n");
printf("        | 4) Added routines to calculate the R squared value of the fit.     |\n");
printf("        |                                                                    |\n");
printf("        | Version 0.1.2 Beta 1                                               |\n");
printf("        | --------------------                                               |\n");
printf("        | 1) Created 3 test cases to test bond, angle and dihedral fits.     |\n");
printf("        |                                                                    |\n");
printf("        | Version 0.2.0 Alpha 1                                              |\n");
printf("        | ---------------------                                              |\n");
printf("        | 1) Updated entire minimisation system to handle 3 different options|\n");
printf("        |    for K - K can now be set to a float value in the job control    |\n");
printf("        |    file which means K will be kept constant throughout. It can be  |\n");
printf("        |    set to READ in which case K is expected after the QM energy.    |\n");
printf("        |    In this situation K can be different for each structure but will|\n");
printf("        |    not be changed during the fitting procedure. Finally K can also |\n");
printf("        |    be set to FIT in which case a starting guess should be given as |\n");
printf("        |    the first float in the start point file. K will then be adjusted|\n");
printf("        |    during the fitting.                                             |\n");
printf("        |                                                                    |\n");
printf("        | Version 0.2.1 Alpha 1                                              |\n");
printf("        | ---------------------                                              |\n");
printf("        | 1) Modified convergence checking to now test how many cycles it has|\n");
printf("        |    been since the minimum value found changed. If this is more than|\n");
printf("        |    100 outer cycles (restarts) convergence is assumed.             |\n");
printf("        | 2) Minor changes to dx variable defaults                           |\n");
printf("        |                                                                    |\n");
printf("        | Version 0.2.2 Alpha 1                                              |\n");
printf("        | ---------------------                                              |\n");
printf("        | 1) Formalised the testing procedure and introduced several new test|\n");
printf("        |    cases that fully test bond, angle and dihedral fitting.         |\n");
printf("        | 2) Created utility program called proc_g98_output that will take a |\n");
printf("        |    series of gaussian files and extract the specified energy, bond |\n");
printf("        |    lenghts, angles sizes and dihedral values etc. This program will|\n");
printf("        |    then print to standard out the data required for the -q QM_data |\n");
printf("        |    file used by the parameter fitter.                              |\n");
printf("        |                                                                    |\n");
printf("        | Version 0.2.2 Alpha 2                                              |\n");
printf("        | ---------------------                                              |\n");
printf("        | 1) Fixed bug in angle fitting that was leading to the angle force  |\n");
printf("        |    constant being treated as being in KCal/(mol deg)^2 instead of  |\n");
printf("        |    KCal/(mol rad)^2.                                               |\n");
printf("        |                                                                    |\n");
printf("        | Version 0.2.2 Alpha 3                                              |\n");
printf("        | ---------------------                                              |\n");
printf("        | 1) Modifications to proc_g98_output to improve numerical stability |\n");
printf("        |    in dihedral calculations. Minor tweaks to try and avoid         |\n");
printf("        |    situations that result in divide by zero errors.                |\n");
printf("        |                                                                    |\n");
printf("        | Version 0.2.2 Alpha 4                                              |\n");
printf("        | ---------------------                                              |\n");
printf("        | 1) Fixed bug in mishandling of dihedralBH_dx, dihedralN_dx and     |\n");
printf("        |    dihedralG_dx parameters by simplex routine.                     |\n");
printf("        | 2) Fixed minor bugs concerning the use of K=FIT in simplex routine.|\n");
printf("        |                                                                    |\n");
printf("        | Version 2.0.0 Alpha 1                                              |\n");
printf("        | ---------------------                                              |\n");
printf("        | 1) Work commenced on a completely new version of the code that     |\n");
printf("        |    will be much more closely coupled to amber. Reading prmtop and  |\n");
printf("        |    mdcrd files. In this way the concept of an atom type should     |\n");
printf("        |    be much better preserved.                                       |\n");
printf("        |                                                                    |\n");
printf("        | Version 2.0.1 Beta 1                                               |\n");
printf("        | --------------------                                               |\n");
printf("        | 1) Commenced testing of parameter fitting with the simplex routine |\n");
printf("        |    and the amber standard force field. Started creating a test     |\n");
printf("        |    suite of example runs for checking that the program is behaving |\n");
printf("        |    as expected.                                                    |\n");
printf("        | 2) Added routine to print the amber energy generated from the final|\n");
printf("        |    parameters.                                                     |\n");
printf("        | 3) Changed the default parameters fitting to include: For Bonds Kr |\n");
printf("        |    and Req are fit. For Angles Kt and THeq are fit. For Dihedrals  |\n");
printf("        |    only Kp is fit by default.                                      |\n");
printf("        |                                                                    |\n");
printf("        | Version 2.0.2 Beta 1                                               |\n");
printf("        | --------------------                                               |\n");
printf("        | 1) Completed the first set of test routines and confirmed that     |\n");
printf("        |    basic fitting of bond, angle and dihedral fitting works.        |\n");
printf("        | 2) Introduced a system whereby the use can manually specify what   |\n");
printf("        |    what parameters they want fitted and which they want to keep    |\n");
printf("        |    constant.                                                       |\n");
printf("        |                                                                    |\n");
printf("        | Version 2.0.3 Alpha 1                                              |\n");
printf("        | ---------------------                                              |\n");
printf("        | 1) Started development of MPI version of code. Currently only      |\n");
printf("        |    the evaluation of the eval_sum_squares_amber_std routine is     |\n");
printf("        |    parallelised.                                                   |\n");
printf("        | 2) Minor code clean up and portability improvement.                |\n");
printf("        |                                                                    |\n");
printf("        | Version 2.0.3 Alpha 2                                              |\n");
printf("        | ---------------------                                              |\n");
printf("        | 1) Minor bug fixes and continued testing of parallel code.         |\n");
printf("        |                                                                    |\n");
printf("        | Version 2.0.3                                                      |\n");
printf("        | -------------                                                      |\n");
printf("        | 1) MPI Code now passes all tests and appears to work well.         |\n");
printf("        |                                                                    |\n");
printf("        | Version 2.0.4                                                      |\n");
printf("        | -------------                                                      |\n");
printf("        | 1) Fixed overflow bug in create input. Program will now no longer  |\n");
printf("        |    spuriously quit with an error claiming the filename to be       |\n");
printf("        |    written would exceed 1023 characters.                           |\n");
printf("        | 2) Fixed incorrect error message in read_qm_energy.                |\n");
printf("        | 3) Changed the calc_r_squared routine to now use the Pearson       |\n");
printf("        |    correlation formula to find R and then square it.               |\n");
printf("        | 4) Fixed minor bug in the calculation of the average function value|\n");
printf("        |    in the simplex routine.                                         |\n");
printf("        | 5) Fixed bug involving the testing of simplex convergence. In some |\n");
printf("        |    cases func_hi could actually contain parameters lower than      |\n");
printf("        |    func_low. This has been fixed by adding a sort routine prior to |\n");
printf("        |    checking for convergence.                                       |\n");
printf("        |                                                                    |\n");
printf("        | Version 2.0.5                                                      |\n");
printf("        | -------------                                                      |\n");
printf("        | 1) Migrated code to Amber CVS tree.                                |\n");
printf("        | 2) Modified output slightly to make it compatible with amber's     |\n");
printf("        |    dacdif.                                                         |\n");
printf("        | 3) Added 1 to the printing of the number of simplex outer cycles   |\n");
printf("        |    since the program counts from zero internally but it makes more |\n");
printf("        |    sense to count from 1 when printing info to the output file.    |\n");
printf("        |                                                                    |\n");
printf("        | Version 2.0.6                                                      |\n");
printf("        | -------------                                                      |\n");
printf("        | 1) Reduced the precision to which some numbers are printed in the  |\n");
printf("        |    output file to make automated testing on different platforms    |\n");
printf("        |    easier.                                                         |\n");
printf("        |                                                                    |\n");
printf("        | Version 2.1.0                                                      |\n");
printf("        | -------------                                                      |\n");
printf("        | 1) Added the option to load or save the dimensions to be fit to a  |\n");
printf("        |    file so that they do not have to be typed in for each run.      |\n");
printf("        | 2) Added the genetic algorithm fitting function as an alternative  |\n");
printf("        |    to simplex for fitting.                                         |\n");
printf("        | 3) Set the random number generator to seed on each run.            |\n");
printf("        |                                                                    |\n");
printf("        | Version 3.0.0                                                      |\n");
printf("        | -------------                                                      |\n");
printf("        | 1) Added bounds checking on input structures that can warn when    |\n");
printf("        |    input structures do not properly span the parameters to be      |\n");
printf("        |    optimised.                                                      |\n");
printf("        | 2) Numerous improvements to the genetic algorithm                  |\n");
printf("        | 3) Implemented parallel version for genetic algorithm and fixed    |\n");
printf("        |    problems with the simplex one.                                  |\n");
printf("        | 4) Rewrote test cases                                              |\n");
printf("        ----------------------------------------------------------------------\n");
printf("\n\n");
}

