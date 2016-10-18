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

/*

  This code will attempt least squares fitting of a series of bond,
  angle and dihedral parameters to a set of QM energies.
  
  It attempts to minimise the following function using the default
  amber bond, angle and dihedral energy equations
  
  f = Sum[1->N]((Sum[1->nb]Bonds + Sum[1->na]Angles + Sum[1->nd]Dihedrals + K - E[QMn])^2)

  Ebond = BOND_FC (BOND - BOND_EQL)^2
  Eangle = ANGLE_FC (ANGLE - ANGLE_EQL)^2
  Edihedral = (DIHEDRAL_BH/2)(1+Cos(DIHEDRAL_MIN_NUM*DIHEDRAL-DIHEDRAL_PHASE))

  Currently it calculates the non-bond terms as part of the energy evaluation but
  does not attempt to optimise these parameters. Note, there is currently no cutoff
  available for the non-bonds.


  FILE FORMATS:

  1) QM_ENERGY_DATA_FILE - this is only needed for fitting and should contain NSTRUCTURES worth of floats with
                           one energy per line in the units specified by ENERGY_UNITS, default is hartrees.

  2) JOB_CONTROL_FILE - This consists of lines beginning with a # which are comments and ignored vs lines beginning
                        with KEYWORD=VALUE.

  3) PRMTOP - Old or New format prmtop is acceptable - see amber manual for format details.

  4) MDCRD FILE - should contain at least NSTRUCTURES worth of structures in the same format as an
                  amber gas phase simulation - NO BOX INFO IS ALLOWED.

  
  
*/

/*Standard Header Files*/
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <signal.h>

#ifdef OPENMP
/* Header file for parallel */
#include <omp.h>
#endif


/*************** CODE NOTES ***************/
/*
   All function prototypes and structure
   definitions are included in the
   function_def.h header file
*/
/******************************************/

#include "function_def.h"

/********** START **********/

int main(int argc, char *argv[])
{
  int int_retval;                              // Return value from subroutines
  global_options_struct global_options;        // Structure containing all of the global options
  parm_struct* parm_datas=NULL;                // Array of parmtop data structures
  coord_set* coords_datas=NULL;                // Array containing the coordinate info from each mdcrd
  bounds_struct bounds_data;

  // Variables for time measurements
  time_t starttime = time(NULL);
  time_t endtime;
  
  // Register signal handling functions
  signal(SIGINT, handle_sigint);
  signal(SIGSEGV, print_backtrace);

  // Initialise memory counters
  global_options.mem_allocated = 0;
  
  /*Setup default options - may be overidden later when we read the job control file,
    the command line options etc.*/
  int_retval = set_default_options(&global_options);
  process_retval(int_retval,global_options.VERBOSITY);

  /*The first job is to process the command line and work out what our job control,
    QM data and start point files are, determine what is on the command line and assign variables
    as necessary, print program history, print program help, etc.
  */
  int_retval = process_command_line(argc, argv, &global_options, &parm_datas, &coords_datas);
  process_retval(int_retval,global_options.VERBOSITY);
  if (global_options.VERBOSITY>=MEDIUM)
  {
    print_program_info();
    printf("*************************************************************************************\n");
    printf("                                Execution started at:\n");
    printf("|                             %s\n",asctime(localtime(&starttime)));
  }
  printf("|\n");
#ifdef OPENMP
#pragma omp parallel
{
  #pragma omp master
  printf("|                 Running OpenMP version of code using %2d processors\n\n", omp_get_num_threads());
}
#endif
  /* Seed the random number generator- allow manual seeding to replicate runs for debugging*/
  if (global_options.RANDOM_SEED == 0) {
    printf("| Random seed = %d\n", (int)starttime);
    srand(starttime);
  } else {
    printf("| Specified random seed = %d\n", global_options.RANDOM_SEED);
    srand(global_options.RANDOM_SEED);
  }

  /* Read in the options in the Job Control File  if specified, otherwise call the wizard*/
  if (global_options.job_control_filename)
    int_retval=read_job_control_file(&global_options,parm_datas, coords_datas);
  else
    int_retval=job_control_wizard(&global_options, parm_datas); 
  process_retval(int_retval,global_options.VERBOSITY);
  
  // Read the prmtop files and allocate memory for the parm structures 
  int_retval=read_prmtops(&global_options, &parm_datas);
  process_retval(int_retval,global_options.VERBOSITY);
  int_retval=process_prmtops(&global_options, parm_datas);
  process_retval(int_retval,global_options.VERBOSITY);
  
  /* If fitting, calculate the dimensions of fit and sanity checks */
  if (global_options.RUNTYPE==FIT || global_options.RUNTYPE==SET_PARAMS) {
//     /* Add more dihedral terms if desired - CURRENTLY IN TESTING*/
//     if (global_options.NDIHEDRALS > 0)
//      not_enough_dihedrals(&parm_data, global_options.NDIHEDRALS);
    
    calc_fit_dimensions(&global_options, parm_datas);
    if (parm_datas[0].ndimensions<1) {
      printf("   ERROR IN MAIN - NDIMENSIONS VALUE OF %d IS LESS THAN 1.\n",parm_datas[0].ndimensions);
      int_retval=ABORT;
      process_retval(int_retval,global_options.VERBOSITY);
    }
  }
  
  // If it's just a parameter-defining run, print out a message and quit- no need to read in anything else
  if (global_options.RUNTYPE==SET_PARAMS) {
    if (global_options.num_prmtops>1)
      print_multiprmtop_summary(&global_options,parm_datas);
    else
      print_parameter_summary(&global_options, &(parm_datas[0]));
    printf("*   Successfully saved desired parameters.\n");
    printf("*   Set runtype to FIT to do a fit with these parameters.\n");
    return SUCCESS;
  }

  // We're doing either the creation of quantum input files or a fit
  // so read in all input structures (will allocate memory inside)
  int_retval=read_mdcrds(&global_options, parm_datas, &coords_datas); 
  process_retval(int_retval,global_options.VERBOSITY);  
  
  // Read in QM data
  if (global_options.RUNTYPE == FIT) {
    int_retval=read_qm(&global_options,parm_datas,coords_datas);
    process_retval(int_retval,global_options.VERBOSITY);
  }
  if (global_options.VERBOSITY>=MEDIUM)
    print_job_control_summary(&global_options, parm_datas, coords_datas); /*summarises the current job control options*/
      
  // end initialization code
    
  // Bounds check if creating input files or doing a fit
  if (global_options.RUNTYPE != SET_PARAMS && global_options.CHECK_BOUNDS!=NO) {
      int_retval = calculate_structure_diversity(&global_options, &bounds_data, &parm_datas[0], &coords_datas[0]);
      process_retval(int_retval,global_options.VERBOSITY);
      int_retval = check_dihedrals(&global_options, &parm_datas[0], &coords_datas[0], &bounds_data);
      process_retval(int_retval,global_options.VERBOSITY);
  }
  
  // If we're creating input then make the input files, and finish
  if (global_options.RUNTYPE == CREATE_INPUT) {
    int_retval=create_qm_input(&global_options, parm_datas, coords_datas);
    process_retval(int_retval,global_options.VERBOSITY);
    clean_up_bounds(&bounds_data, &parm_datas[0]);
    global_unlock(&global_options, &parm_datas, &coords_datas);
    return SUCCESS;
  }
  // Otherwise, do the fit
  int_retval=do_fit(&global_options, parm_datas, coords_datas);
  process_retval(int_retval,global_options.VERBOSITY);
  
  /* Check that the results are reasonable based on the given data in the input structures.       *
  * Dihedrals are already checked at the beginning, but do bonds and angles here to make sure    *
  * the result is in a well-represented area in the input data.                                  */
  if (global_options.CHECK_BOUNDS) {
    int_retval = check_angles(&global_options, &parm_datas[0], &coords_datas[0], &bounds_data);
    process_retval(int_retval,global_options.VERBOSITY);
    int_retval = check_bonds(&global_options, &parm_datas[0], &coords_datas[0], &bounds_data);
    process_retval(int_retval,global_options.VERBOSITY);
    clean_up_bounds(&bounds_data, &parm_datas[0]);
  }

  if (global_options.VERBOSITY>=MEDIUM)
  {
    endtime = time(NULL);
    printf("\n|             Program Execution Completed at: %s",asctime(localtime(&endtime)));
    printf("|                            Elapsed Time = %.2f seconds\n",difftime(endtime, starttime));
    printf("*************************************************************************************\n");
  }
  global_unlock(&global_options, &parm_datas, &coords_datas);
  
  exit(SUCCESS);
}

