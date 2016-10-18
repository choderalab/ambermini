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

/** @file defaults.c
 * Contains the default option for everything.
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "function_def.h"
#include "constants.h"

/**
 * Initializes variables representing program options to their defaults.
 * @param[in,out] global_options The global options structure where defaults will be set.
 * @return Integer indicating success or failure.
 */
int set_default_options(global_options_struct *global_options)
{
  /*SET 1 - COMMAND LINE SWITCHES*/                                        
  global_options->VERBOSITY = MEDIUM;
  global_options->RANDOM_SEED = 0; // no seed specifed by default

  // No default job control filename- call the wizard if this happens 
  global_options->job_control_filename = NULL;

  // Prmtop and mdcrd list are null
  global_options->prmtop_list=NULL;
  global_options->mdcrd_list=NULL;

  /*SET 2 - JOB CONTROL OPTIONS*/
  /*Default RUNTYPE is to fit*/
  global_options->RUNTYPE=FIT;
  global_options->function_calls=0;
  
  /*Default QM File Format is Gaussian*/
  global_options->QMFILEFORMAT=GAUSSIAN;
  global_options->COORD_FORMAT=TRAJECTORY;
  
  global_options->QMFILEOUTSTART=NULL;
  global_options->QMFILEOUTEND=NULL;
  
  global_options->TOTAL_STRUCTURES=0;      /*Number of structures to fit to - default = zero so we return an error later*/
  global_options->SORT_ENERGY=NO;        // By default don't sort
  
  global_options->ALGORITHM=SIMPLEX;                    /*fitting routine to be used*/
  global_options->FUNC_TO_FIT=SUM_SQUARES_AMBER_STANDARD;       /*The function to be used to fit to the energy surface*/

  global_options->SCEE=1.2;
  global_options->SCNB=2.0;

  global_options->PARAMETERS_TO_FIT=DEFAULT;
  global_options->PARAMETER_FILE_NAME=NULL;
  
  global_options->K_FIT=NO;
  
  /*options for simplex routine*/
  global_options->CONV_LIMIT=1.0E-15;           /*Convergence limit for fitting routine*/
  
  /*options for bounds checking*/
  global_options->CHECK_BOUNDS=WARN;
  global_options->ANGLE_LIMIT=PI/20.0;
  global_options->DIHEDRAL_SPAN=10;
  global_options->BOND_LIMIT=0.1;
  
  /*Note, the smaller the dx values are the longer convergence will take*/
  global_options->BONDFC_dx=5.0;
  global_options->BONDEQ_dx=0.02;
  global_options->ANGLEFC_dx=1.0;
  global_options->ANGLEEQ_dx=0.05;
  global_options->DIHEDRALBH_dx=0.2;
  global_options->DIHEDRALN_dx=0.01;
  global_options->DIHEDRALG_dx=0.05;
  global_options->K_dx=10.0;

  /* Options for the genetic algorithm */
  global_options->NOPTIMIZATIONS=50;
  global_options->MAX_GENERATIONS=10000;
  global_options->GENERATIONS_TO_CONVERGE=20;
  global_options->SEARCH_SPACE=-1.0; // search all of parameter spaces
  global_options->MUTATION_RATE=0.10;
  global_options->PARENT_PERCENT=0.25;
  global_options->GENERATIONS_TO_SIMPLEX=5;
  global_options->GENERATIONS_WITHOUT_SIMPLEX=10;
  
  /* Options for the dihedral least squares fit */
  global_options->FIT_PHASE=NO;
  
  global_options->SCATTERPLOTS=NO;
  
  global_options->QMHEADER=NULL;
  global_options->QM_SYSTEM_CHARGE=0;
  global_options->QM_SYSTEM_MULTIPLICITY=1;
  global_options->QM_ENERGY_UNITS=HARTREE;
  global_options->QM_FORCE_UNITS=KCALMOL_ANGSTROM;
  global_options->WRITE_FRCMOD=NULL;
  global_options->WRITE_ENERGY=NULL;
  global_options->NDIHEDRALS=0;
  return SUCCESS;
}
