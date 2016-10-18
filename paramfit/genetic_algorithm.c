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

/** @file genetic_algorithm.c
 * Contains functions used by the genetic algorithm fitting.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "function_def.h"
#include "constants.h"

/**
 * Conducts the majority of the genetic algorith minimization.
 * 
 * @param[in] global_options The global options structure, with algorithm options inside
 * @param[in,out] parm_datas Array containing the parameters for each molecule, will be updated each generation 
 * @param[in] coords_data Array containing coordinates and QM data for each input structure for each molecule
 * @return Integer indicating success or failure.
 */
int minimise_function_genetic(global_options_struct *global_options, parm_struct *parm_datas, coord_set *coords_data)
{
  /* variables to use go here */
  double **data_matrix;           // will hold all of the possible optimizations
  double *function_values;       // will hold the amber function results for each optimization
  double previous_best;           // used for calculating the convergence ratio
  double initial_result;          // for comparison with the converged results
  double mean=0.0;                // useful statistics
  double mutation_rate = global_options->MUTATION_RATE;
                                  // percentage of values randomly changed each generation
  double parent_percent = global_options->PARENT_PERCENT;
                                  // this constitutes the top percent of solutions
  double weighting;               // used to change values by a certain amount
  

  int generation;                 // how many generations have passed
  int converged;                  // how generations in a row have passed with a low first derivative
  int gen_since_simplex=0;        // how many generations have passed since a simplex iteration
  int col, row;
  int retval;
  int k_offset = (global_options->K_FIT==YES) ? 1 : 0;
  int i;                          // counting variable for loops
  
  /* Set up the global options structure for the simplex calls. This structure has simplex set up so that
   * only a few iterations will be run */
  double old_conv_limit;
  double old_verbosity = global_options->VERBOSITY;
  if (global_options->ALGORITHM==BOTH) old_conv_limit = global_options->CONV_LIMIT;
  global_options->CONV_LIMIT=10.0;
  
  /* print out a header with some information */
  if (global_options->VERBOSITY>=MEDIUM)
  {
    printf("   --------------------------GENETIC ALGORITHM MINIMISATION ---------------------------\n");
    if (global_options->FUNC_TO_FIT==SUM_SQUARES_AMBER_STANDARD)
      printf("   Minimising function SUM_SQUARES_AMBER_STANDARD, using the GENETIC ALGORITHM\n");
    else if (global_options->FUNC_TO_FIT==AMBER_FORCES)
      printf("   Minimising function AMBER_FORCES using the GENETIC ALGORITHM\n");
    else
      printf("   Minimising function UNKNOWN, using the GENETIC ALGORITHM\n");      
    if (global_options->GENERATIONS_TO_SIMPLEX>0)
      printf("   Running SIMPLEX REFINEMENT every %d converged gen, then break for %d\n",
             global_options->GENERATIONS_TO_SIMPLEX,global_options->GENERATIONS_WITHOUT_SIMPLEX);
    else
      printf("   Running pure genetic algorithm with NO SIMPLEX REFINEMENT\n");
    printf("   -------------------------------------------------------------------------------------\n");
    // Print out algorithm parameters
    printf("   GENERATIONS_TO_CONVERGE = %d    MAX_GENERATIONS = %d    CONV_LIMIT = %4.2E\n",
           global_options->GENERATIONS_TO_CONVERGE, global_options->MAX_GENERATIONS, global_options->CONV_LIMIT);
    printf("   OPTIMIZATIONS = %d    PARENT_PERCENT = %4.2f\n", global_options->NOPTIMIZATIONS, global_options->PARENT_PERCENT);
    printf("   SEARCH_SPACE = %.2f    MAX_GENERATIONS = %d\n",
           global_options->SEARCH_SPACE, global_options->MAX_GENERATIONS);
    if (global_options->VERBOSITY>=MEDIUM)
      printf("   ------------------------------------ CONVERGENCE ------------------------------------\n");
    else
      printf("   -------------------------------------------------------------------------------------\n");
    fflush(stdout); /*Flush the printf buffer*/
  }
 
 /* Allocate the 1D matrix to hold the function results for each optimization */
 if (global_options->VERBOSITY>=HIGH)
   printf("Allocating %d bytes for *function_values.\n",(int)((global_options->NOPTIMIZATIONS)*sizeof(double)));
 function_values=(double *)calloc(global_options->NOPTIMIZATIONS, sizeof(double));
 if (function_values==NULL)
 {
   malloc_failure_double("minimise_function_genetic", "function_values", (global_options->NOPTIMIZATIONS));
   return ALLOC_FAIL;
 }
 
  /* Allocate space for the optimizations */
  /* The data matrix is 2 dimensional, with num_optimizations rows and NDIMENSIONS cols */
  data_matrix=alloc_data_matrix(global_options->NOPTIMIZATIONS, parm_datas[0].ndimensions);
  if (data_matrix==NULL)
  {
    printf("   ERROR - MALLOC FAILED FOR data_matrix\n");
    printf("   COMMAND WAS data_matrix=alloc_2D_double(global_options->NOPTIMIZATIONS,parm_datas[0].ndimensions)\n");
    return ALLOC_FAIL;
  }
  if ( (int)(global_options->NOPTIMIZATIONS*parent_percent) <= 2)
  {
    printf("   ERROR: NOPTIMIZATIONS=%d is too small to conduct the algorithm.\n", global_options->NOPTIMIZATIONS);
    free_data_matrix(data_matrix, global_options->NOPTIMIZATIONS);
    free(function_values);
    return INVALID_DATA;
  }
  
// end memory allocation
  
  /* Initially fill the data matrix with the initial parameters. The first *
    * potential solution will be the starting parameters and will be        *
    * changed only at the end of the initial generation creation.           */
  // Parameters to fit should be identical across all input structures, so just use the first
  retval = modify_params_scratch_data(global_options, &(parm_datas[0]), data_matrix[0], READ);
  if (retval != parm_datas[0].ndimensions)
  {
    printf("   ERROR IN minimise_function_genetic()\n");
    printf("             modify_params_scratch_data() DOES NOT MATCH\n");
    printf("             param dimensions OF %d. ABORTING RUN.\n",parm_datas[0].ndimensions);
    fflush(stdout); /*Flush the printf buffer*/
    free_data_matrix(data_matrix, global_options->NOPTIMIZATIONS);
    free(function_values);
    /*This is most probably a bug so return an unknown error*/
    return(FAILURE);
  }
  for (row=1; row<global_options->NOPTIMIZATIONS; ++row)
    for (col=0; col<parm_datas[0].ndimensions; ++col)
      data_matrix[row][col] = data_matrix[0][col];
 
  /* Evaluate the function for the initial parameters - works with only 1 prmtop also*/
  initial_result = eval_sum_squares_amber_std_multiprmtop(global_options, parm_datas, coords_data, NULL);
  int params_found = 0;
  int count = 0;
  double delta = 0;

  if (global_options->VERBOSITY>=HIGH) {
    printf("ORIGINAL:\n");
    for (col=0; col<parm_datas[0].ndimensions; ++col)
      printf("%8.3f", data_matrix[0][col]);
    printf("\tFUNCTION: %f\n---\n", initial_result);
  }
  
  /* Create the initial generation randomly. Some values, like K, can be any number to start.      *
  * Others will be random within an allowed domain, such as angles and dihedral phases, using the *
  * initial parameters as a starting point.                                                       */
  for (row=0; row < global_options->NOPTIMIZATIONS; ++row)
  {
    params_found = 0;
    for(col=0; col<parm_datas[0].ndimensions; ++col)
    {
      // search space value determines how far to look or to look everywhere
      if (global_options->SEARCH_SPACE <= 0.0) // search everything
      {
        delta = 0.0;
        data_matrix[row][col] = 0.0;
        if (global_options->K_FIT==YES)
          parm_datas[0].K = (rand() & 1) ? rand() : -1.0*rand();
      }
      else
      {
        // deviation from initial parameter is random, can be positive or negative
        delta = global_options->SEARCH_SPACE * ( (double)rand()/(double)RAND_MAX );
        if (rand() & 1) delta*= -1.0;
      }

      // K can be any value, positive or negative, is not based on initial parameters unless
      // an initial guess is provided
      if (col==0 && global_options->K_FIT==YES && global_options->SEARCH_SPACE > 0.0)
      {
          data_matrix[row][col] += data_matrix[row][col] * delta;
      }
      else if (col-k_offset<global_options->BOND_PARAMS)
      {
        /*Current column represents a bonding parameter. Loop through the bonds until we find it*/
        params_found=0;
        for (count=0;count<parm_datas[0].unique_bonds_found;++count)
        {
          // Base Kr off the initial value or make it between 100 and 1000
          if (parm_datas[0].bond_data[count].DO_FIT_KR==YES)
          {
            /*At this point we check to see if our params_found to date match the column we are on*/
            if (params_found==col-k_offset)
            {
              if (data_matrix[row][col] == 0.0)
          data_matrix[row][col] = 900*((double)rand()/(double)RAND_MAX) + 100.0;
              else
          data_matrix[row][col] += data_matrix[row][col] * delta;
              /*no point continuing for this column*/
              break;
            }
            ++params_found;
          }
          // Req is based on the initial value or between 0 and 3
          if (parm_datas[0].bond_data[count].DO_FIT_REQ==YES)
          {
            if (params_found==col-k_offset)
            {
              if (data_matrix[row][col] == 0.0)
                data_matrix[row][col] = 3.0*((double)rand()/(double)RAND_MAX);
              else
                data_matrix[row][col] += data_matrix[row][col] * delta;
              break;
            }
            ++params_found;
          }
        }
      }
      else if (col-k_offset<global_options->ANGLE_PARAMS+global_options->BOND_PARAMS)
      {
        /*Current column represents an angle parameter. Loop through the angles until we find it*/
        params_found=0;
        for (count=0;count<parm_datas[0].unique_angles_found;++count)
        {
          // Kt based off initial parameters or between 0 and 170
          if (parm_datas[0].angle_data[count].DO_FIT_KT==YES)
          {
            if (params_found+global_options->BOND_PARAMS==col-k_offset)
            {
              if (data_matrix[row][col] == 0.0)
                data_matrix[row][col] = 170.0*((double)rand()/(double)RAND_MAX);
              else
                data_matrix[row][col] += data_matrix[row][col] * delta;
              break;
            }
            ++params_found;
          }
          // Theta needs to stay between zero and PI
          if (parm_datas[0].angle_data[count].DO_FIT_THEQ==YES)
          {
            if (params_found+global_options->BOND_PARAMS==col-k_offset)
            {
              data_matrix[row][col] += data_matrix[row][col]*delta;
              if (data_matrix[row][col] <= 0.0 || data_matrix[row][col] > PI)
                data_matrix[row][col] = ((double)rand()/(double)RAND_MAX) * PI;
              break;
            }
            ++params_found;
          }
        }
      }
      else if (col-k_offset<global_options->DIHEDRAL_PARAMS+global_options->ANGLE_PARAMS+global_options->BOND_PARAMS)
      {
        /*Current column represents a dihedral parameter. Loop through the dihedrals until we find it*/
        params_found=0;
        for (count=0;count<parm_datas[0].unique_dihedrals_found;++count) {
          for (i=0; i<parm_datas[0].dihedral_data[count].num_terms; ++i) {
            // Kp based on initial value or between -30 and 30
            if (parm_datas[0].dihedral_data[count].term[i].DO_FIT_KP==YES)
            {
              if (params_found+global_options->BOND_PARAMS+global_options->ANGLE_PARAMS==col-k_offset)
              {
                if (data_matrix[row][col] <= 0.0) // 0 to 10
                  data_matrix[row][col] = -30.0 + 60.0*((double)rand()/(double)RAND_MAX);
                else
                  data_matrix[row][col] += data_matrix[row][col]*delta;
                break;
              }
              ++params_found;
            }
            // Np based off initial value, is an integer, or between 1 and 5
            if (parm_datas[0].dihedral_data[count].term[i].DO_FIT_NP==YES) {
              if (params_found+global_options->BOND_PARAMS+global_options->ANGLE_PARAMS==col-k_offset) {
                if (fabs(data_matrix[row][col])>5.0 || data_matrix[row][col]==0.0)
                {
                  data_matrix[row][col] = (double)(rand() % 6) + 1.0;
                  if (rand()&1) data_matrix[row][col] *= -1.0;
                }
                else
                  data_matrix[row][col] += (double)( (data_matrix[row][col]*delta) );
                break;
              }
              ++params_found;
            }
            // Phase between 0 and PI, based off initial value or random
            if (parm_datas[0].dihedral_data[count].term[i].DO_FIT_PHASE==YES) {
              if (params_found+global_options->BOND_PARAMS+global_options->ANGLE_PARAMS==col-k_offset) {
                data_matrix[row][col] += data_matrix[row][col]*delta;
                if (data_matrix[row][col] <= 0.0 || data_matrix[row][col] > PI)
                  data_matrix[row][col] =  PI*((double)rand()/(double)RAND_MAX);
                break;
              }
              ++params_found;
            }
          }
        }
      }
    }
  }
  
/* Evaluate the function for each of the initial optimisations and sort them by this fitness so    *
 * that the top row contains the best optimisation so far.                                         */
for(row=0; row<global_options->NOPTIMIZATIONS; ++row)
{
  // put the current parameters into parm datas for the function evaluation
  retval = update_prmtop_data(global_options, parm_datas, data_matrix[row]);
  process_retval(retval, global_options->VERBOSITY);
  
  // evaluate function for the row- evaluation done in parallel within the function
  if (global_options->FUNC_TO_FIT==SUM_SQUARES_AMBER_STANDARD) {
    function_values[row] = eval_sum_squares_amber_std_multiprmtop(global_options, parm_datas, coords_data, NULL);
  } else if (global_options->FUNC_TO_FIT==AMBER_FORCES) {
    function_values[row] = eval_sum_amber_forces_multiprmtop(global_options, parm_datas, coords_data);
    function_values[row] += eval_sum_squares_amber_std_multiprmtop(global_options, parm_datas, coords_data, NULL); // TEST
  } else {
    printf("   ERROR: FUNCTION %d IS NOT YET IMPLEMENTED\n",global_options->FUNC_TO_FIT);
    free_data_matrix(data_matrix, global_options->NOPTIMIZATIONS);
    free(function_values);
    return(NOT_IMPLEMENTED);
  }
} // end fitness function

  // sort the results by fitness (insertion sort)
  // data[0] will point to the row with smallest function evaluation
  for (row=1; row<global_options->NOPTIMIZATIONS; ++row) {
    double value = function_values[row];
    double *temp = data_matrix[row];
    
    col = row-1;
    while(col >=0 && function_values[col] > value) {
      data_matrix[col+1] = data_matrix[col];
      function_values[col+1] = function_values[col];
      col--;
    }
    data_matrix[col+1] = temp;
    function_values[col+1] = value;
  }
  
  generation = 0;
  converged  = 0;

/* This is the main loop, which continues through the process of recombination, function          *
 * evaluation, and bounds checking until there is neglible improvement for                        *
 * global_options->GENERATIONS_TO_CONVERGE or the maximum number of generations is hit.           */
//    while (generation < global_options->MAX_GENERATIONS && (converged < global_options->GENERATIONS_TO_CONVERGE || generation < 50))
   while (generation < global_options->MAX_GENERATIONS && converged < global_options->GENERATIONS_TO_CONVERGE)
  {
    // hold on to stats from the previous generation
    previous_best = function_values[0];
      
  /* Do recombination and mutation row-by row in this large loop */
    int parent1, parent2;
    for (row=(int)(global_options->NOPTIMIZATIONS*parent_percent)+1; row<global_options->NOPTIMIZATIONS; ++row) {
      // Do recombination, check for clones and mutate on a column-by-column basis- treat variables as independent of one another.
      int recomb_method = rand() & 1;
      
      if (recomb_method == 1) weighting = rand() % parm_datas[0].ndimensions; // pick split point
      
      for (col=0; col<parm_datas[0].ndimensions; ++col)
      {
        parent2 = rand() % global_options->NOPTIMIZATIONS*parent_percent;
        /* Recombination */
        do {
          parent1 = rand() % (int)(global_options->NOPTIMIZATIONS*parent_percent);
        } while (parent1==parent2);

        /* There are three ways of doing recombination. Which one this child gets is determined randomly. *
         * For RANGE, the child values are chosen to be random values between that of the two parents.    *
         * For SHUFFLE, the child values are taken from either parent at random.                          *
         * For CHOP, a split point is chosen. All values before come from one parent, all after from the  *
         * other.
         * These methods allow more diversity to be introduced, because there is no one way of combining  *
         * parent values for this problem due to coupling, etc.                                           *
         * 
         * Testing has shown that using the RANGE and CHOP methods produces the best results so these are *
         * the ones used by the code.                                                                     */
        // RANGE METHOD
        // Child value is somewhere in the range between the two parents
        if (recomb_method == 0 ) {
          weighting = ( (double)rand() / (double)RAND_MAX );
          if (data_matrix[parent1][col] > data_matrix[parent2][col])
            data_matrix[row][col] = data_matrix[parent2][col] + weighting*(data_matrix[parent1][col]-data_matrix[parent2][col]);
          else
            data_matrix[row][col] = data_matrix[parent1][col] + weighting*(data_matrix[parent2][col]-data_matrix[parent1][col]);
        }
        // CHOP METHOD
        // Split point was chosen before the column loop
        else {
          if (col < weighting) data_matrix[row][col] = data_matrix[parent1][col];
          else                 data_matrix[row][col] = data_matrix[parent2][col];
        }
          
        /*
        // SHUFFLE METHOD - not used
        // Child values taken from either parent at random
        else if (recomb_method == 1) {
          if (rand() & 1) data_matrix[row][col] = data_matrix[parent1][col];
          else            data_matrix[row][col] = data_matrix[parent2][col];
        }
        */
       
        /* Mutation */
        if (rand() % (int)(100*mutation_rate) == 0)
          do_mutation(global_options, parm_datas, data_matrix[row], col, YES);
        
      } // end col loop

  } // end recombination counter
  
  /* Force the whole bottom 20% to be completely mutated- this can pop it out of local minima. *
  * The fitness of these won't be evaluated until after the next child is created, but that's  *
  * okay since they're too low in the hierarchy to get to be parents most of the time.         */
  for (row=0.8*global_options->NOPTIMIZATIONS; row<global_options->NOPTIMIZATIONS; ++row) {
    for (col=0; col<parm_datas[0].ndimensions; ++col) {
      do_mutation(global_options, parm_datas, data_matrix[row], col, YES);
  } }

  /* Run simplex on a random 5% for a very short number of iterations if more than 5 generations have  *
   * passed. This should accelerate the algorithm, especially near the end. Does not work with         *
   * multiprmtop fits. Will not work with forces since force evaluations are so expensive.             */
    if (global_options->GENERATIONS_TO_SIMPLEX>0 && global_options->FUNC_TO_FIT==SUM_SQUARES_AMBER_STANDARD &&
        converged >= global_options->GENERATIONS_TO_SIMPLEX && gen_since_simplex > global_options->GENERATIONS_WITHOUT_SIMPLEX) {
      global_options->VERBOSITY=OFF; // dont want a ton of simplex output
      for (count=0; count<(int)(0.05*global_options->NOPTIMIZATIONS)+1; ++count) {
        do {
          row = rand() % global_options->NOPTIMIZATIONS;
        } while (row==0);
        update_prmtop_data(global_options, parm_datas, data_matrix[row]);
        minimise_function_simplex(global_options , &(parm_datas[0]), coords_data);
        update_data_matrix(global_options, parm_datas, data_matrix[row]);
      } 
      global_options->VERBOSITY=old_verbosity;
      gen_since_simplex=0;
    } else {
      gen_since_simplex++;
    }
   
    // Do a fitness evaluation for each solution
    for (row=0; row<global_options->NOPTIMIZATIONS; ++row) {
      retval = update_prmtop_data(global_options, parm_datas, data_matrix[row]);
      process_retval(retval, global_options->VERBOSITY);
        
      if (global_options->FUNC_TO_FIT==SUM_SQUARES_AMBER_STANDARD) {
        int num_altered = 0;
        function_values[row] = eval_sum_squares_amber_std_multiprmtop(global_options, parm_datas, coords_data, &num_altered);
        if (num_altered) {
      printf("UPDATING row %f\n", row);
          retval = update_data_matrix(global_options, parm_datas, data_matrix[row]);
          process_retval(retval, global_options->VERBOSITY);
        }
      } else if (global_options->FUNC_TO_FIT==AMBER_FORCES) {
        function_values[row] = eval_sum_amber_forces_multiprmtop(global_options, parm_datas, coords_data);
        function_values[row] += eval_sum_squares_amber_std_multiprmtop(global_options, parm_datas, coords_data, NULL); // TEST
      } else {
        printf("   ERROR: FUNCTION %d IS NOT YET IMPLEMENTED\n",global_options->FUNC_TO_FIT);
        free_data_matrix(data_matrix, global_options->NOPTIMIZATIONS);
        free(function_values);
        return(NOT_IMPLEMENTED);
      }
    }
    
    // sort the array by fitness
    for (row=1; row<global_options->NOPTIMIZATIONS; ++row)
    {
      double value = function_values[row];
      double *temp = data_matrix[row];
      col = row-1;
      while(col >=0 && function_values[col] > value)
      {
        data_matrix[col+1] = data_matrix[col];
        function_values[col+1] = function_values[col];
        col--;
      }
      data_matrix[col+1] = temp;
      function_values[col+1] = value;
    }
  
    /* Print the top five -- for debugging */
//     for (col=0; col<parm_datas[0].ndimensions; ++col)
//       printf("%8.3f", data_matrix[0][col]);
//     printf("\tFUNCTION: %f\n---\n", function_values[0]);

  /* Check for convergence- after a while the algorithm will begin to stagnate, improving the best  *
   * optimization less and less often. When there is no improvement in the best solution for many   *
   * generations in a row (when the first derivative is zero), the algorithm is considered to have  *
   * converged. The number of generations that must happen with this condition is held in           *
   * global_options->GENERATIONS_TO_CONVERGE.                                                       */

  if ( generation >= 3 && fabs(previous_best - function_values[0]) < 1.0e-5 )
      ++converged;
    else
      converged = 0;
    // calculate some useful statistics
    mean = 0.0;
    for (row=0; row<global_options->NOPTIMIZATIONS; ++row)
      mean += function_values[row];
    mean /= (double)global_options->NOPTIMIZATIONS;
    

    if (generation < 3)
    printf("   Gen %4i:  Best= %10.5f  \tMean= %10.5g Elapsed=%10d/%d\n",
           generation, function_values[0], mean, generation, 3);
    else
      printf("   Gen %4i:  Best= %10.5f  \tMean= %10.5g Conv=   %10d/%d\n",
      generation, function_values[0], mean, converged, global_options->GENERATIONS_TO_CONVERGE);
    fflush(stdout);
    
    ++generation;
  } // end main generation loop
  
  // Once convergence has been reached, copy the answer into the parm_data struct
  if (generation < global_options->MAX_GENERATIONS)
    printf("| Took    %d generations to converge.\n", generation);
  retval = update_prmtop_data(global_options, parm_datas, data_matrix[0]);
  process_retval(retval, global_options->VERBOSITY);

  // Do some cleanup
  free_data_matrix(data_matrix, global_options->NOPTIMIZATIONS);
  free(function_values);
  function_values = NULL;
  
  // Reset the global options simplex setting if we changed it
  if (global_options->ALGORITHM==BOTH) global_options->CONV_LIMIT = old_conv_limit;
    
  // Warn if it's exited because of having passed the max number of generations
  if (generation >= global_options->MAX_GENERATIONS)
  {
    printf("\n");
    printf("   WARNING: Generation count of %d exceeded maximum of %d.\n", generation, global_options->MAX_GENERATIONS);
    printf("            Either 1) Fit fewer dimensions or             \n");
    printf("                   2) Increase MAX_GENERATIONS.         \n\n");
    fflush(stdout);
    return EXCEEDEDMAXITERATIONS;
  }

  return SUCCESS;
 }

/**
 * Conducts mutation or a validity check on the given element in the given row
 * of the data matrix. This is called in the genetic algorithm on an element with
 * a set probability, most elements will never see this function.
 * 
 * @param[in] global_options The global options structure
 * @param[in] parm_data The parameter structure, used to look up what kind of parameter element to mutate is.
 * @param[in,out] row Array representing a row in the data matrix where mutation will happen
 * @param[in] col Index in the row to mutate
 * @param[in] do_mutate True if value is to be changed, false to just conduct a bounds check
 * @return Integer indicating success or failure
 */
int do_mutation(global_options_struct *global_options, parm_struct *parm_data, double *row, int col, bool_t do_mutate)
{
  double mutation_amount = 0.10;      // how much mutation can change each value             
  double delta = 0.0;                 // how much value will change by (Check bounds = 0)
  int k_offset = (global_options->K_FIT==YES) ? 1 : 0;
  int params_found = 0;
  int count = 0;
  int i;

  if (do_mutate==YES)
  {
    delta = (double)rand() / ((double)RAND_MAX+1.) * mutation_amount; // a value > 0, < rec_deviation
    if (rand() & 1) delta*=-1.0;
    // If the value is zero, change it a bit so mutation will work and it will stay nonzero
    if (row[col] == 0.0) row[col] = (rand() & 1) ? 1.0 : -1.0;
  }

  /* This part either mutates or does a bounds check. If it is set to mutate, delta is nonzero
  * and so the data point will be changed. Otherwise, delta is zero so nothing is added to the
  * data and a bounds check happens.
  */
  if (col==0 && global_options->K_FIT==YES)
  {
    // Don't mutate K too much since it is probably a big number and will just kill your fitness
    // since K doesn't have a local minimum and is fit pretty much instantly. Instead let it
    // change by +- 10.0 for small adjustments when other parameters change.
    row[col] += ( (double)rand() / ((double)RAND_MAX+1.) )*10.0;
  }
  else if (col-k_offset<global_options->BOND_PARAMS)
  {
    params_found=0;
    for (count=0; count<parm_data->unique_bonds_found; ++count)
    {
      if (parm_data->bond_data[count].DO_FIT_KR==YES)
      {
        if (params_found==col-k_offset)
        {
          row[col] += row[col] * delta;
          // Bond force constant between 100 and 1000
          if (row[col] < 100.0 || row[col] > 1000.0)
          row[col] = 900.0*( (double)rand()/(double)RAND_MAX ) + 100.0;
          break;
        }
        ++params_found;
      }
      if (parm_data->bond_data[count].DO_FIT_REQ==YES)
      {
        if (params_found==col-k_offset)
        {
	  // Bond length between 0 and 3 Angstroms
	  if (row[col] <= 0.0 || row[col] > 3.0)
	    row[col] = 3.0*( (double)rand()/(double)RAND_MAX );
          row[col] += row[col]*delta;
          break;
        }
        ++params_found;
      }
    }
  }
  else if (col-k_offset<global_options->ANGLE_PARAMS+global_options->BOND_PARAMS)
  {
    params_found=0;
    for (count=0;count<parm_data->unique_angles_found;++count)
    {
      if (parm_data->angle_data[count].DO_FIT_KT==YES)
      {
        if (params_found+global_options->BOND_PARAMS==col-k_offset)
        {
	  // Angle force constant between 0 and 170
	  if (row[col] <= 0.0 || row[col] > 170.0)
	    row[col] = 170.0*( (double)rand() / (double)RAND_MAX );
          row[col] += row[col]*delta;
          break;
        }
        ++params_found;
      }
      // Theta needs to stay between zero and PI
      if (parm_data->angle_data[count].DO_FIT_THEQ==YES)
      {
        if (params_found+global_options->BOND_PARAMS==col-k_offset)
        {
          row[col] += row[col]*delta;
          
          if (row[col] < 0.0 || row[col] > PI)
          {
            row[col] = ((double)rand()/(double)RAND_MAX) * PI;
            if (do_mutate==NO) printf("Fixed angle %i theta!\n", col);
          }
          break;
        }
        ++params_found;
      }
    }
  }
  else if (col-k_offset<global_options->DIHEDRAL_PARAMS+global_options->ANGLE_PARAMS+global_options->BOND_PARAMS)
  {
    params_found=0;
    for (count=0;count<parm_data->unique_dihedrals_found;++count) {
      for (i=0; i<parm_data->dihedral_data[i].num_terms; ++i) {
        // Kr based on initial value
        if (parm_data->dihedral_data[count].term[i].DO_FIT_KP==YES) {
          if (params_found+global_options->BOND_PARAMS+global_options->ANGLE_PARAMS==col-k_offset) {
            row[col] += row[col]*delta;
            // Dihedral force constant between 0 and 30
            if (row[col] < -30.0 || row[col] > 30.0)
            row[col] = -30.0 + 60.0*( (double)rand() / (double)RAND_MAX );
            break;
          }
          ++params_found;
        }
        // Np positive or negative, integer in range 1-5
        if (parm_data->dihedral_data[count].term[i].DO_FIT_NP==YES)
        {
          if (params_found+global_options->BOND_PARAMS+global_options->ANGLE_PARAMS==col-k_offset)
          {
            row[col] += row[col]*delta;
            
            if (row[col] < -5.0 || row[col] > 5.0)
            {
              row[col] = (double)((rand() % 6) + 1);
              if (rand()&1) row[col]*=-1.0;
            }
      
      // ensure that it is an integer
        row[col] = (double)((int)(row[col]+0.5));
            break;
          }
          ++params_found;
        }
        // Phase between 0 and PI, based off initial value
        if (parm_data->dihedral_data[count].term[i].DO_FIT_PHASE==YES)
        {
          if (params_found+global_options->BOND_PARAMS+global_options->ANGLE_PARAMS==col-k_offset)
          {
            row[col] += row[col]*delta;
            
            if (row[col] < 0.0 || row[col] > PI)
            {
              if (do_mutate==NO) printf("Fixed dihedral %i phase!\n", col);
              row[col] = ((double)rand()/(double)RAND_MAX) * PI;
            }
            break;
          }
          ++params_found;
        }
      }
    }
  }
  return SUCCESS;
}


/**
 * Allocates a more traditional, non-contiguous 2D array for the genetic algorithm.
 * We don't use alloc_2d_double because that gets a contiguous block of memory and
 * we specifically need non-contiguous since we will be swapping the rows around
 * in sorting each generation, and free-ing the contiguous matrix later is a pain
 * since it is now unclear which was the original first row.
 * @see free_data_matrix
 * 
 * @param[in] rows The number of rows to be in the array
 * @param[in] cols The number of columns to be in the array
 * @return Pointer to a 2D double array of the specified size
 */
 double **alloc_data_matrix(int rows, int cols) {
   double **dm = (double **)malloc(rows*sizeof(double *));
   if (dm==NULL) return NULL;
   
   int i;
   for (i=0; i<rows; ++i) {
     dm[i] = (double *)calloc(cols, sizeof(double));
     if (dm[i]==NULL) return NULL;
   }
   return dm;
 }
 
 /**
  * Frees the data matrix used for the genetic algorithm.
  * @see alloc_data_matrix
  * 
  * @param[in,out] dm Pointer to the 2D array to be freed
  * @param[in] rows The number of rows in the array
  */
 void free_data_matrix(double **dm, int rows) {
   int i;
   for (i=0; i<rows; ++i)
     free(dm[i]);
   free(dm);
 }
     













