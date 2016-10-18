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

/** @file fitting_control.c
 * Contains the main routine that controls the fitting process,
 * which calls the relevant algorithm after doing setup appropriate
 * to the options the user has selected.
 */

#include <stdio.h>
#include <stdlib.h>

#include "function_def.h"

/**
 * Conducts the fit.
 * Prints a parameter summary, does bounds checking, calls appropriate algorithm,
 * verifies and prints final parameters, and writes output file formats according to the
 * options the user has chosen.
 * 
 * @param[in] global_options The global options structure containing user defined options
 * @param[in,out] parm_data Array of initial parameters and what to fit for each molecule, will have final parameters inserted.
 * @param[in] coords_data Array of coordinate sets containing the atomic coordinates and QM data for input structures.
 * @return Integer indicating success or failure
 */ 
int do_fit(global_options_struct *global_options, parm_struct *parm_datas, coord_set *coords_datas)
{	
  double initial_function_value;
  double initial_r_squared_value;
  double initial_energy_struc_1;
  double final_r_squared_value;
  double final_function_value;
  double final_energy_struc_1;
  int retval, i;

  /*The initial parameters and details on the fitting etc will have been listed by the options_summary
    routine. So the first thing we do here is do an initial evaluation of our function*/

  if (global_options->VERBOSITY >= MEDIUM || global_options->RANDOM_SEED != 0) // don't print for test cases
  {
    printf("\n   ------------------------------- INITIAL PARAMETERS --------------------------------\n");
    if (global_options->num_prmtops >1 ) {
      printf("   --------------------- ( printing only parameters to be fit ) ----------------------\n");
      print_multiprmtop_summary(global_options, parm_datas);
    } else
      print_parameter_summary(global_options, parm_datas);
    printf("   -----------------------------------------------------------------------------------\n");
  }
  if (global_options->FUNC_TO_FIT==SUM_SQUARES_AMBER_STANDARD)
  {   
      /* Calculate how many bonds, angles, and dihedrals will be fit. */
      // using one parameter set rather than the whole array of them in these functions.
      global_options->BOND_PARAMS = calculate_no_fit_params(&parm_datas[0], BONDS);
      global_options->ANGLE_PARAMS = calculate_no_fit_params(&parm_datas[0], ANGLES);
      global_options->DIHEDRAL_PARAMS = calculate_no_fit_params(&parm_datas[0], DIHEDRALS);
      
      // Calculate initial function and R^2 value if desired
      if (global_options->VERBOSITY >= MEDIUM) {
        initial_function_value = eval_sum_squares_amber_std_multiprmtop(global_options, parm_datas, coords_datas, NULL);
        initial_r_squared_value = calc_r_squared_multiprmtop(global_options, parm_datas, coords_datas);
        initial_energy_struc_1=eval_amber_std_for_single_struct(global_options, &(parm_datas[0]), coords_datas[0].struc);
        printf("   Sum of squares for initial parameters = %15.10f kcal^2/mol^2\n",initial_function_value);
        printf("   R^2 value for initial parameters      = %10.6f\n",initial_r_squared_value);
        printf("   Calculated energy with initial parameters for structure 1 = %10.6f KCal/mol\n",initial_energy_struc_1+parm_datas[0].K);
        printf("   Actual energy for structure 1 should be                   = %10.6f KCal/mol\n\n",coords_datas[0].struc[0].energy);
      }
   }
  else if (global_options->FUNC_TO_FIT==AMBER_FORCES) {
     /* Calculate how many bonds, angles, and dihedrals will be fit */
    global_options->BOND_PARAMS = calculate_no_fit_params(&parm_datas[0], BONDS);
    global_options->ANGLE_PARAMS = calculate_no_fit_params(&parm_datas[0], ANGLES);
    global_options->DIHEDRAL_PARAMS = calculate_no_fit_params(&parm_datas[0], DIHEDRALS);
      
    // Mark the atoms that are to be fit
    for (i=0; i<global_options->num_prmtops; ++i) {
      parm_datas[i].fit_atom = (int *)malloc(coords_datas[0].natoms*sizeof(int));
      if (!parm_datas[i].fit_atom)
        malloc_failure_int("do_fit", "parm_data->fit_atom", coords_datas[0].natoms);
      mark_relevant_atoms(global_options, &parm_datas[i]);
    }
     
    // Calculate initial function and R^2 value if desired
    if (global_options->VERBOSITY >= MEDIUM) {
      initial_function_value = eval_sum_amber_forces_multiprmtop(global_options, parm_datas, coords_datas) 
        + eval_sum_squares_amber_std_multiprmtop(global_options, parm_datas, coords_datas, NULL);
      initial_r_squared_value = calc_r_squared_multiprmtop(global_options, parm_datas, coords_datas);
      initial_energy_struc_1=eval_amber_std_for_single_struct(global_options, &parm_datas[0], coords_datas[0].struc);
      
      force_struct *eval = (force_struct*)malloc(coords_datas[0].natoms*sizeof(force_struct));
      if (eval==NULL) return ALLOC_FAIL;
      eval_amber_forces_single_struct(global_options, &parm_datas[0], &coords_datas[0], eval, 0);
    
      printf("   Sum of force difference magnitudes for initial parameters = %15.10f kcal/mol-A\n", initial_function_value);
      printf("   R^2 value for initial parameters                          = %10.6f\n", initial_r_squared_value);
      printf("   Calculated energy with initial parameters for structure 1 = %10.6f KCal/mol\n",initial_energy_struc_1);
      printf("   Calculated forces with initial parameters atom 1, structure 1 = %3.4f %3.4f %3.4f kcal/mol-A\n\n", eval[0].x, eval[0].y, eval[0].z);
      free(eval);
    }
  }
  else if (global_options->FUNC_TO_FIT==DIHEDRAL_LEAST_SQUARES) {
    /* Calculate how many bonds, angles, and dihedrals will be fit. */
    printf("!  Warning- setting all dihedral terms to fit pk\n");
    if (global_options->FIT_PHASE) printf("!  Warning- setting all dihedral terms to fit phase\n");
    else printf("!  Warning- fitting no dihedral phases\n");
    int j,k;
    for (i=0; i<global_options->num_prmtops; ++i) {
      for (k=0; k<parm_datas[i].unique_dihedrals_found; ++k) {
        for (j=0; j<parm_datas[i].dihedral_data[k].num_terms; ++j) {
          parm_datas[i].dihedral_data[k].term[j].DO_FIT_KP=YES;
          if(global_options->FIT_PHASE)
            parm_datas[i].dihedral_data[k].term[j].DO_FIT_PHASE=YES;
          else
            parm_datas[i].dihedral_data[k].term[j].DO_FIT_PHASE=NO;
        }
      }
    }
    
    global_options->BOND_PARAMS = calculate_no_fit_params(&parm_datas[0], BONDS);
    global_options->ANGLE_PARAMS = calculate_no_fit_params(&parm_datas[0], ANGLES);
    global_options->DIHEDRAL_PARAMS = calculate_no_fit_params(&parm_datas[0], DIHEDRALS);
    printf("   Using dihedral minimization algorithm described by Chad Hopkins and Adrian Roitberg.\n");
  }
  else {
    printf("   ERROR IN DO_FIT - FUNCTION %d IS NOT YET IMPLEMENTED\n",global_options->FUNC_TO_FIT);
    return NOT_IMPLEMENTED;
  }
   
  // If K is being fit along with other parameters, give an error
  if (global_options->K_FIT==YES && global_options->PARAMETERS_TO_FIT!=K_ONLY) {
    printf("   ERROR: K is being fit along with other parameters. This will produce inaccurate\n");
    printf("         results and is disallowed. Specify a value for K or set PARAMETERS_TO_FIT=K_ONLY\n");
    return INVALID_DATA;
  }
  
  if (global_options->FUNC_TO_FIT==DIHEDRAL_LEAST_SQUARES) { // fitting function and algorithm combined
    retval=dihedral_least_squares(global_options, parm_datas, coords_datas, global_options->FIT_PHASE);
  } 
  else if (global_options->ALGORITHM==SIMPLEX) {
     retval=minimise_function_simplex(global_options, parm_datas, coords_datas);
  } 
  else if (global_options->ALGORITHM==GENETIC) {
    if (global_options->VERBOSITY>=HIGH)
      print_parameter_summary(global_options,&parm_datas[0]);
    retval = minimise_function_genetic(global_options, parm_datas, coords_datas);
  }
  else if (global_options->ALGORITHM==BOTH)
  {
    retval = minimise_function_genetic(global_options, parm_datas, coords_datas);
    retval = minimise_function_simplex(global_options, parm_datas, coords_datas);
  }
  else if (global_options->ALGORITHM==NONE)
  {
    retval=SUCCESS;
    // Filler for printing information or something
  }
  else
  {
    printf("   ERROR IN do_fit() - UNKNOWN FITTING FUNCTION: %d\n",global_options->ALGORITHM);
    return UNKNOWN_OPT;
  }
  
  if (retval!=SUCCESS ) {
      /*We can trap one error here - that is if we quit due to exceeding maximum iterations*/
      /*Other errors are fatal*/
      if (retval==EXCEEDEDMAXITERATIONS || retval==MINSTATIC)
        printf("!  Warning - final parameters do NOT represent a converged fit.\n");
      else
        return(retval);
   }

    // Print out a final summary
    if (global_options->VERBOSITY>=MEDIUM)
    {
      printf("   ------------------------------- FINAL PARAMETERS ---------------------------------\n");
      if (global_options->num_prmtops>1) {
        printf("   --------------------- ( printing only parameters to be fit ) ---------------------\n");
        print_multiprmtop_summary(global_options,parm_datas);
      } else
        print_parameter_summary(global_options,&parm_datas[0]);
        printf("   ----------------------------------------------------------------------------------\n");
      
      if (global_options->VERBOSITY>=MEDIUM)
        printf("|   Called the fitness function %d times.\n", global_options->function_calls);
      
      if (global_options->FUNC_TO_FIT==SUM_SQUARES_AMBER_STANDARD) {
        final_function_value = eval_sum_squares_amber_std_multiprmtop(global_options, parm_datas, coords_datas, NULL);
        final_r_squared_value=calc_r_squared_multiprmtop(global_options, parm_datas, coords_datas);
        final_energy_struc_1=eval_amber_std_for_single_struct(global_options, &parm_datas[0], coords_datas[0].struc);
        
        printf("   Function value with fitted parameters  =  %12.4f, R^2 = %12.4f\n",final_function_value,final_r_squared_value);
        printf("   Calculated energy with fitted parameters for structure 1 = %11.4f KCal/mol\n",final_energy_struc_1);
        printf("\n");
      }
      else if (global_options->FUNC_TO_FIT==AMBER_FORCES) {
        double finale = eval_sum_squares_amber_std_multiprmtop(global_options, parm_datas, coords_datas, NULL);
        final_function_value = eval_sum_amber_forces_multiprmtop(global_options, parm_datas, coords_datas) 
                             + finale;
        final_r_squared_value=calc_r_squared_multiprmtop(global_options, parm_datas, coords_datas);
        final_energy_struc_1=eval_amber_std_for_single_struct(global_options, &parm_datas[0], coords_datas[0].struc);
        
        force_struct *eval = (force_struct*)malloc(coords_datas[0].natoms*sizeof(force_struct));
        if (eval==NULL) return ALLOC_FAIL;
        eval_amber_forces_single_struct(global_options, &parm_datas[0], &coords_datas[0], eval, 0);
        
        printf("   Final force value = %12.4f, final energy function = %12.4f\n", final_function_value-finale, finale);
        printf("   Function value with fitted parameters  =  %12.4f, R^2 = %12.4f\n",final_function_value,final_r_squared_value);
        printf("   Calculated energy with fitted parameters for structure 1 = %11.4f KCal/mol\n",final_energy_struc_1);
        printf("   Calculated forces with fitted parameters atom 1, structure 1 = %3.4f %3.4f %3.4f kcal/mol-A\n\n", eval[0].x, eval[0].y, eval[0].z);
        printf("\n");
        free(eval);
      
      printf("   Experimental: Writing theta.txt and mag.txt with difference for atoms that are fit\n");
      print_forces(global_options, &parm_datas[0], &coords_datas[0]);
      }
    }
    
    // if desired, save the frcmod file
    if (global_options->WRITE_FRCMOD)
      write_frcmod(global_options, &parm_datas[0]);
    
    // if desired, write the list of final energies for comparison
    if (global_options->WRITE_ENERGY) {
        for (i=0; i<global_options->num_prmtops; ++i)
          write_energy(global_options,&parm_datas[i],&coords_datas[i]);
    }
    
   return SUCCESS;
}


