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

/** @file simplex.c 
 * Contains the simplex minimization algorithm
 * Requires (N+1)*N (doubles) Storage for the simplex array as well as (3*NDIMENSIONS)+10 doubles as scratch.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "function_def.h"

/**
 * The simplex minimization function.
 * Written mostly by Ross Walker
 * 
 * @param[in] global_options The global options structure
 * @param[in,out] parm_data Pointer to array of parameter structures, updated at each iteration
 * @param[in] coords_data Array of coordinate sets containing all the input structures with QM data
 * @return Integer indicating success or failure
 */
int minimise_function_simplex(global_options_struct *global_options, parm_struct *parm_datas, coord_set *coords_data)	
{
    double **vertex_matrix;
    double *function_results;
    double *parameters;   /*used for scratch storage of the variable parameters - NDIMENSIONS LONG*/
    int inner_itr;  /*Number of inner iterations to perform for each outer iteration*/
    int outer_itr_max;  /*maximum number of outer iterations to conduct before quitting*/
    int outer_itr_count;
    int inner_itr_count;
    int col,row;
    short int loop_exit_condition;
    int function_calls;

    int func_low, func_hi, func_next_hi; /*Holds the function_results location for the lowest func eval, highes func eval and
                                            next highest func eval*/

    double conv_ratio;  /*used for checking convergence - stores ratio of [Y(highest) - Y(lowest)]/[(Y(Highest)+Y(lowest)]*/

    double function_results_sum; /*Stores the sum of a complete set of N+1 function evaluations - used for printing the average Y value*/
    
    double *avg_point;                /*Average of the simplex points except the highest point - needs to be NDIM columns long*/
    double *reflect_point;            /*The point that represents the highest point through the centre of the simplex, NDIM columns long*/
    double reflect_point_func;        /*Value of function at reflect_point*/
    double *extended_reflect_point;   /*used as an extension of our reflect point to extrude further and see if the function still drops*/
    double extended_reflect_point_func; /*Value of function at extended reflect point*/

    double pointa,pointb,pointc,pointd;  /*extra points used instead of default simplex contraction*/
    double expMin;  /*Used in the polynomial expansion*/
      
    double reflect_ratio;
    double extension_length;
    double reflection_reduction_ratio;

    int min_unchanged_counter;
    double previous_min; /*These two variables are used to check if our minimum value has changed in the last 100 steps and quit if they haven't*/

    int retval;

    short int k_offset;  /*set to 1 if K_FIT = true - used to calculate the offset for the beginning of the bonds, angles etc in the parameter scratch array*/
    int params_found;
    int count;

    /*various parameters for tweaking the method - probably best to leave at the defaults*/
    reflect_ratio = 1.0;
    extension_length = 2.0;
    reflection_reduction_ratio=0.5;  /*ratio by which to reduce the reflection point when extended reflection point did not get us to a better minimum*/
    /*end*/
    
    min_unchanged_counter=0;
    previous_min=0.0;

    function_calls=0;
    
    inner_itr=NSIMPLEX_INNER_PER_DIM*parm_datas[0].ndimensions;
    outer_itr_max=NSIMPLEX_OUTER_MAX;
    inner_itr_count=0;
    outer_itr_count=0;
    func_low=0;
    func_hi=0;
    function_results_sum=0.0;
    
    loop_exit_condition=EXCEEDEDMAXITERATIONS;
    if (global_options->VERBOSITY>=MEDIUM) {
      printf("   --------------------------------- SIMPLEX MINIMISATION ----------------------------\n");
      if (global_options->FUNC_TO_FIT==SUM_SQUARES_AMBER_STANDARD)
        printf("   Minimising function SUM_SQUARES_AMBER_STANDARD, using the SIMPLEX METHOD\n");
      else if (global_options->FUNC_TO_FIT==AMBER_FORCES)
        printf("   Minimising function AMBER_FORCES, using the SIMPLEX METHOD\n");
      else
        printf("   Minimising function UNKNOWN, using the SIMPLEX METHOD\n");
      if (global_options->VERBOSITY >= MEDIUM)
        printf("   -------------------------------------- CONVERGENCE --------------------------------\n");
      else
        printf("   -----------------------------------------------------------------------------------\n");
      fflush(stdout); /*Flush the printf buffer*/
    }
    /*The parameters are stored in the form of unique bond terms, angle terms and dihedral terms*/

    /*Stage 1 - we need a matrix of simplex vertices - This is a 2D matrix with NDIM+1 rows each of NDIM columns*/
    vertex_matrix=alloc_2D_double(parm_datas[0].ndimensions+1,parm_datas[0].ndimensions);
    if (vertex_matrix==NULL)
    {
      printf("*** ERROR - MALLOC FAILED FOR vertex_matrix\n");
      printf("*** COMMAND WAS vertex_matrix=alloc_2D_double(parm_datas[0].ndimensions+1,parm_datas[0].ndimensions)\n");
       return ALLOC_FAIL;
    }
    
    /*We also need some scratch space for other results etc.*/
    /*First of all we need N+1 doubles for a function_results array, N doubles for our avg point and N doubles for our reflect point
      and N doubles for our extended_reflect_point*/
    if (global_options->VERBOSITY>=HIGH )
      printf("Allocating %d bytes for *function_results.\n",(int)((parm_datas[0].ndimensions+1)*sizeof(double)));
    function_results=(double *)calloc((parm_datas[0].ndimensions+1),sizeof(double));
    if (function_results==NULL)
    {
      /*malloc failure*/
      malloc_failure_double("minimise_function_simplex", "function_results", (parm_datas[0].ndimensions+1));
      return ALLOC_FAIL;
    }
    /*Allocate memory for our avg_point*/
    if (global_options->VERBOSITY>=HIGH)
      printf("Allocating %d bytes for *avg_point.\n",(int)(parm_datas[0].ndimensions*sizeof(double)));
    avg_point=(double *)calloc(parm_datas[0].ndimensions,sizeof(double));
    if (avg_point==NULL)
    {
      /*malloc failure*/
      malloc_failure_double("minimise_function_simplex", "avg_point", parm_datas[0].ndimensions);
      return ALLOC_FAIL;
    }
    /*Allocate memory for our reflect_point*/
    if (global_options->VERBOSITY>=HIGH)
        printf("Allocating %d bytes for *reflect_point.\n",(int)(parm_datas[0].ndimensions*sizeof(double)));
    reflect_point=(double *)calloc(parm_datas[0].ndimensions,sizeof(double));
    if (reflect_point==NULL)
    {
      /*malloc failure*/
      malloc_failure_double("minimise_function_simplex", "reflect_point", parm_datas[0].ndimensions);
      return ALLOC_FAIL;
    }
    /*Allocate memory for our extended_reflect_point*/
    if (global_options->VERBOSITY>=HIGH)
      printf("Allocating %d bytes for *extended_reflect_point.\n",(int)(parm_datas[0].ndimensions*sizeof(double)));
    extended_reflect_point=(double *)calloc(parm_datas[0].ndimensions,sizeof(double));
    if (extended_reflect_point==NULL)
    {
      /*malloc failure*/
      malloc_failure_double("minimise_function_simplex", "extended_reflect_point", parm_datas[0].ndimensions);
      return ALLOC_FAIL;
    }
    /*allocate memory for our parameter scratch array*/
    if (global_options->VERBOSITY>=HIGH)
      printf("Allocating %d bytes for *parameters.\n",(int)(parm_datas[0].ndimensions*sizeof(double)));
    parameters=(double *)calloc(parm_datas[0].ndimensions,sizeof(double));
    if (parameters==NULL)
    {
      /*malloc failure*/
      malloc_failure_double("minimise_function_simplex", "parameters", parm_datas[0].ndimensions);
      return ALLOC_FAIL;
    }
    
    if (global_options->K_FIT==YES)
      k_offset=1;
    else
      k_offset=0;

    /*Sanity check here - sum of above should equal the number of dimensions (per prmtop)*/
    if (parm_datas[0].ndimensions != k_offset+global_options->BOND_PARAMS+global_options->ANGLE_PARAMS+global_options->DIHEDRAL_PARAMS )
    {
      printf("*** ERROR IN minimise_function_simplex() - SUM OF k_offset(%d), bond_params(%d),\n",k_offset, global_options->BOND_PARAMS);
      printf("***          angle_params(%d) and dihedral_params(%d) DOES NOT EQUAL\n",global_options->ANGLE_PARAMS,global_options->DIHEDRAL_PARAMS);
      printf("***          NDIMENSIONS OF %d. ABORTING RUN.\n",parm_datas[0].ndimensions);
       fflush(stdout); /*Flush the printf buffer*/
       /*This is most probably a bug so return an unknown error*/
       return(FAILURE);
    }
    
    /*Start main outer loop*/
    /*This loop essentially restarts the SIMPLEX minimisation after every inner_itr loops using the lambda values:
    K_dx - if K_FIT = YES
    BONDFC_dx
    BONDEQ_dx
    ANGLEFC_dx
    ANGLEEQ_dx
    DIHEDRALBH_dx
    DIHEDRALN_dx
    DIHEDRALG_dx
    These are adjusted by (1-(RAND_RATIO/2)) + RAND_RATIO*rand()/RAND_MAX in order to remove any symmetry in the matrix
    that might lead to convergence problems
    */

/*BEGIN SIMPLEX OUTER LOOP*/
    while (outer_itr_count<outer_itr_max)
    {
      /*Before we fill the simplex array we ideally need a linear array of parameters - this should be
        refilled from the parm_struct each time we loop here. We are best having an external routine to
        do this.
      */
      // Take the initial parameters from the first prmtop 
      retval=modify_params_scratch_data(global_options, &(parm_datas[0]), parameters, READ);
      /*Check we get back the number of dimensions we expect*/
      if (retval!=parm_datas[0].ndimensions)
      {
        printf("*** ERROR IN minimise_function_simplex() - RETURN VALUE OF %d FROM\n",retval);
        printf("***          modify_params_scratch_data() DOES NOT MATCH\n");
        printf("***          parm dimensions OF %d. ABORTING RUN.\n",parm_datas[0].ndimensions);
         fflush(stdout); /*Flush the printf buffer*/
         /*This is most probably a bug so return an unknown error*/
         return(FAILURE);
      }
      for(col=0;col<parm_datas[0].ndimensions;++col) /*Loop for filling the simplex matrix*/
      {
        for(row=0;row<=parm_datas[0].ndimensions;++row) /*loop over our N+1 rows*/
           vertex_matrix[row][col]=parameters[col];     /*Fill this row of the simplex matrix with the [col]th variable parameter*/
           
        /*we now add the gamma value to each of our points in turn
          to get our N+1 starting points. i.e
             Pi = P0 + Gamma*ei
          where ei are the unit vectors

          since we can currently have different values of Gamma for bonds, angles and dihedrals and the
          parameters they contain we have to make sure we add the correct correction to each
        */
        if (col==0 && global_options->K_FIT==TRUE)
        {
          /*The first parameter of our parameters array is K*/
          vertex_matrix[col][col]+=global_options->K_dx*( (1-(RAND_RATIO/2)) + RAND_RATIO*rand()/(double)RAND_MAX );
        }
        else if (col-k_offset<global_options->BOND_PARAMS)
        {
          /*Current column represents a bonding parameter. Loop through the bonds until we find it*/
          params_found=0;
          for (count=0;count<parm_datas[0].unique_bonds_found;++count)
          {
             if (parm_datas[0].bond_data[count].DO_FIT_KR==YES)
             {
               /*At this point we check to see if our params_found to date match the column we are on*/
               if (params_found==col-k_offset)
               {
                  /*found it and it is a bond force constant*/
                  vertex_matrix[col][col]+=global_options->BONDFC_dx*( (1-(RAND_RATIO/2)) + RAND_RATIO*rand()/(double)RAND_MAX );
                  /*no point continuing for this column*/
                  break;
               }
               ++params_found;
             }
             if (parm_datas[0].bond_data[count].DO_FIT_REQ==YES)
             {
               /*At this point we check to see if our params_found to date match the column we are on*/
               if (params_found==col-k_offset)
               {
                  /*found it and it is a bond eq constant*/
                  vertex_matrix[col][col]+=global_options->BONDEQ_dx*( (1-(RAND_RATIO/2)) + RAND_RATIO*rand()/(double)RAND_MAX );
                  /*no point continuing for this column*/
                  break;
                }
                ++params_found;
             }
          }
        }
        else if (col-k_offset<global_options->ANGLE_PARAMS+global_options->BOND_PARAMS )
        {
          /*Current column represents an angle parameter. Loop through the angles until we find it*/
          params_found=0;
          for (count=0;count<parm_datas[0].unique_angles_found;++count)
          {
             if (parm_datas[0].angle_data[count].DO_FIT_KT==YES)
             {
               /*At this point we check to see if our params_found to date match the column we are on*/
               if (params_found+global_options->BOND_PARAMS==col-k_offset)
               {
                  /*found it and it is an angle force constant*/
                  vertex_matrix[col][col]+=global_options->ANGLEFC_dx*( (1-(RAND_RATIO/2)) + RAND_RATIO*rand()/(double)RAND_MAX );
                  /*no point continuing for this column*/
                  break;
               }
               ++params_found;
             }
             if (parm_datas[0].angle_data[count].DO_FIT_THEQ==YES)
             {
               /*At this point we check to see if our params_found to date match the column we are on*/
               if (params_found+global_options->BOND_PARAMS==col-k_offset)
               {
                  /*found it and it is an angle eq constant*/
                  vertex_matrix[col][col]+=global_options->ANGLEEQ_dx*( (1-(RAND_RATIO/2)) + RAND_RATIO*rand()/(double)RAND_MAX );
                  /*no point continuing for this column*/
                  break;
                }
                ++params_found;
             }
          }
        }
        else if (col-k_offset<global_options->DIHEDRAL_PARAMS+global_options->ANGLE_PARAMS+global_options->BOND_PARAMS )
        {
          /*Current column represents a dihedral parameter. Loop through the dihedrals until we find it*/
          params_found=0; int tm=0;
          for (count=0;count<parm_datas[0].unique_dihedrals_found;++count) {
            for (tm=0; tm<parm_datas[0].dihedral_data[count].num_terms; ++tm) {
              if (parm_datas[0].dihedral_data[count].term[tm].DO_FIT_KP==YES)
              {
                /*At this point we check to see if our params_found to date match the column we are on*/
                if (params_found+global_options->BOND_PARAMS+global_options->ANGLE_PARAMS==col-k_offset)
                {
                    /*found it and it is a dihedral barrier height*/
                    vertex_matrix[col][col]+=global_options->DIHEDRALBH_dx*( (1-(RAND_RATIO/2)) + RAND_RATIO*rand()/(double)RAND_MAX );
                    /*no point continuing for this column*/
                    break;
                }
                ++params_found;
              }
              if (parm_datas[0].dihedral_data[count].term[tm].DO_FIT_NP==YES)
              {
                /*At this point we check to see if our params_found to date match the column we are on*/
                if (params_found+global_options->BOND_PARAMS+global_options->ANGLE_PARAMS==col-k_offset)
                {
                    /*found it and it is an dihedral N constant*/
                    vertex_matrix[col][col]+=global_options->DIHEDRALN_dx*( (1-(RAND_RATIO/2)) + RAND_RATIO*rand()/(double)RAND_MAX );
                    /*no point continuing for this column*/
                    break;
                  }
                  ++params_found;
              }
              if (parm_datas[0].dihedral_data[count].term[tm].DO_FIT_PHASE==YES)
              {
                /*At this point we check to see if our params_found to date match the column we are on*/
                if (params_found+global_options->BOND_PARAMS+global_options->ANGLE_PARAMS==col-k_offset)
                {
                    /*found it and it is an dihedral PHASE constant*/
                    vertex_matrix[col][col]+=global_options->DIHEDRALG_dx*( (1-(RAND_RATIO/2)) + RAND_RATIO*rand()/(double)RAND_MAX );
                    /*no point continuing for this column*/
                    break;
                  }
                  ++params_found;
              }
            }
          }
        }
      } /*End of for (col=0;col<parm_datas[0].ndimensions;++col) - loop for filling simplex matrix*/
      /*
        What we now have at this point is a matrix that consists of NDIM+1 parameter sets made up
        by our existing parameters and each parameter multiplied by a value dx
      */
      if (global_options->VERBOSITY>=HIGH)
      {
        /*OUTPUT THE COMPLETE MATRIX TO THE SCREEN*/
        printf("vertex_matrix for outer step %d:\n",outer_itr_count);
        for(row=0;row<=parm_datas[0].ndimensions;row++)
        {
          for(col=0;col<parm_datas[0].ndimensions;col++)
          {
            printf("row(%d), col(%d): %8.6f\n",row,col,vertex_matrix[row][col]);
          }
          printf("\n");
        }
        printf("\n");
      } 

      /*The next stage is to evaluate our function for each set of parameters - for this we will adjust our
        temporary scratch array of parameters, copy the data back into the parm_data and then evaluate the function.

        Note, the function evaluation is done in parallel,using the data in parm_data.
        
        We work by taking each row of the simplex array, which represents a set of parameters, in turn. The evaluated
        function we store in our *function_results_array. This should already have had enough memory allocated such that
        it can store NDIMENSIONS+1 worth of function evaluations. Every time we call our function evaluator we will
        also update the counter.

        Since the last row of our simplex matrix contains our parameters without the dx additions so we
        should get *parameters and consequently parm_data filled with the parameters we started this loop with.

      */
      for (row=0;row<=parm_datas[0].ndimensions;++row) {
        for (col=0;col<parm_datas[0].ndimensions;++col) {
           /*fill the *parameters scratch array with this set of parameters*/
           parameters[col]=vertex_matrix[row][col];
           /*remember to copy this across to parm_data since this is what the function evaluator works from*/
           retval=update_prmtop_data(global_options,parm_datas,parameters);
           process_retval(retval, global_options->VERBOSITY);
        }
        /*Evaluate the function for this parameter row*/
        if (global_options->FUNC_TO_FIT==SUM_SQUARES_AMBER_STANDARD)
        {
           function_results[row] = eval_sum_squares_amber_std_multiprmtop(global_options, parm_datas, coords_data, NULL);
           ++function_calls;
        }
        else if (global_options->FUNC_TO_FIT==AMBER_FORCES) {
          function_results[row] = eval_sum_amber_forces_multiprmtop(global_options, parm_datas, coords_data);
          ++function_calls;
        }
        else
        {
          printf("*** FUNCTION %d IS NOT YET IMPLEMENTED\n",global_options->FUNC_TO_FIT);
          return(NOT_IMPLEMENTED);
        }
      } /*End of loop for filling our function_results array for our matrix of params for(row=0;row<parm_datas[0].ndimensions;++row)*/

      /*
        Now we have our first set of starting points we run our inner loop - which is the actual
        Simplex downhill method - we run this for inner_itr steps after which we check for convergence and if
        we are converged we break out of this outer loop - otherwise we iterate the outer loop again...
      */
/*INNER LOOP*/
      for (inner_itr_count=1;inner_itr_count<=inner_itr;++inner_itr_count)
      {
        function_results_sum=0.0;
        func_low=0;
        if (function_results[0]>function_results[1])
        {
          /*Our first set of parameters give a higher result than the second set*/
          func_hi=0;
          func_next_hi=1;
        }
        else
        {
          func_hi=1;
          func_next_hi=0;
        }
        for (row=0;row<=parm_datas[0].ndimensions;++row) /*Note we loop over a total of NDIMENSIONS+1 here*/
        {
           /*loop over each row = each parameter set*/
           if ( function_results[row] < function_results[func_low] ) /*Test if function with current row params is less than the lowest found so far*/
           {
             func_low=row; /*If it is lower change our low locator to this row*/
           }
           if ( function_results[row] > function_results[func_hi] ) /*Repeat test to see if it is higher than our current highest*/
           {
             /*If it is higher then update the highest locator and the next highest locator*/
             func_next_hi=func_hi;
             func_hi=row;
           }
           else if ( function_results[row] > function_results[func_next_hi] )
           {
             /*We are between highest and second highest*/
             /*only update if we aren't currently at the highest value*/
             if (row!=func_hi)
             {
               func_next_hi=row;
             }
           }
        } /*End of loop over rows (parameter sets) for (row=0;row<=parm_datas[0].ndimensions;++row)*/

        /*If this is our very first time through the outer loop and the inner loop print the very first step*/
        if (global_options->VERBOSITY>=MEDIUM && inner_itr_count==1 && outer_itr_count==0)
        {
          conv_ratio=2.0*fabs((function_results[func_hi]-function_results[func_low]))/(fabs(function_results[func_hi])+fabs(function_results[func_low]));
          /*Calculate average*/
          function_results_sum=0.0;
          for (row=0;row<=parm_datas[0].ndimensions;++row) /*Note we loop over a total of NDIMENSIONS+1 here*/
          {
             function_results_sum+=function_results[row]; /*used for printing the average over all the functions with the different params, should converge*/
          }                                                                              
          printf("   Step %5d: Conv=%12.4E min=%12.4f,max=%12.4f avg%12.4f\n",outer_itr_count,
          conv_ratio,function_results[func_low],function_results[func_hi],function_results_sum/(parm_datas[0].ndimensions+1));
          fflush(stdout); /*Flush the printf buffer*/
        }
        /*If diagnostics are debug we will print a status line showing our convergence here for all inner loops*/
        if (global_options->VERBOSITY>=HIGH)
        {
          conv_ratio=2.0*fabs(function_results[func_hi]-function_results[func_low])
                    /(fabs(function_results[func_hi])+fabs(function_results[func_low]));
          function_results_sum=0.0;
          for (row=0;row<=parm_datas[0].ndimensions;++row) /*Note we loop over a total of NDIMENSIONS+1 here*/
            function_results_sum+=function_results[row]; /*used for printing the average over all the functions with the different params, should converge*/
         
          printf("RST: %d - INNER: %d > Conv=%10.8E Ymin=%10.8f, Ymax=%10.8f, Yavg=%10.8f\n",outer_itr_count,inner_itr_count,
          conv_ratio,function_results[func_low],function_results[func_hi],function_results_sum/(parm_datas[0].ndimensions+1));
          fflush(stdout);
        }
        /*Make sure our avg point is cleared before use*/
        for(col=0;col<parm_datas[0].ndimensions;++col)
          avg_point[col]=0.0;
        
        /*Our avg_point should be the average of all points except the highest one*/
        for(row=0;row<=parm_datas[0].ndimensions;++row)
          if(row!=func_hi)
            for(col=0;col<parm_datas[0].ndimensions;++col)
              avg_point[col]+=vertex_matrix[row][col];
          
        /*now convert avg_point from the sum to the average and at the same time make our reflection point*/
        for(col=0;col<parm_datas[0].ndimensions;++col)
        {
          avg_point[col]/= (double)parm_datas[0].ndimensions;
          reflect_point[col]=(1.0+reflect_ratio)*avg_point[col]-reflect_ratio*vertex_matrix[func_hi][col]; /*Each parameter for our row of params that gives the higest function result*/
          parameters[col]=reflect_point[col];  /*Copy the reflection point to our parameters array*/
        }
        /*Need to update the parm_data arrays*/
        retval=update_prmtop_data(global_options, parm_datas, parameters);
        process_retval(retval, global_options->VERBOSITY);
                  
        /*Now we need to evaluate our function at our reflection point and see where it lies amongst all our other points*/
        if (global_options->FUNC_TO_FIT==SUM_SQUARES_AMBER_STANDARD)
        {
          reflect_point_func = eval_sum_squares_amber_std_multiprmtop(global_options, parm_datas, coords_data, NULL);
          ++function_calls;
        }
        else if (global_options->FUNC_TO_FIT==AMBER_FORCES) {
          reflect_point_func = eval_sum_amber_forces_multiprmtop(global_options, parm_datas, coords_data);
          ++function_calls;
        }
        else
        {
          printf("*** FUNCTION %d IS NOT YET IMPLEMENTED\n", global_options->FUNC_TO_FIT);
          return(NOT_IMPLEMENTED);
        }

        /* *** NOW BEGINS THE TESTING - WE NEED TO TEST WHERE OUR NEW REFLECT POINT FITS IN WITH OUR CURRENT POINTS ****/
        if(reflect_point_func<=function_results[func_low]) /*Is our reflect_point lower than all the other points?*/
        {
          /*it is a better estimation of the minima - extrude point further and see if it drops again*/
          for(col=0;col<parm_datas[0].ndimensions;++col)
          {
            extended_reflect_point[col]=extension_length*reflect_point[col]+(1.0-extension_length)*(avg_point[col]);
            /*And update our parameters with this new extension of the reflected point*/
            parameters[col]=extended_reflect_point[col];
          }
          
          retval=update_prmtop_data(global_options, parm_datas, parameters);
          process_retval(retval, global_options->VERBOSITY);
          
            /*Now evaluate our function with the extended point to see if it is lower than the original reflection*/
          if (global_options->FUNC_TO_FIT==SUM_SQUARES_AMBER_STANDARD)
          {
            extended_reflect_point_func = eval_sum_squares_amber_std_multiprmtop(global_options, parm_datas, coords_data, NULL);
            ++function_calls;
          }
          else if (global_options->FUNC_TO_FIT==AMBER_FORCES) {
            extended_reflect_point_func = eval_sum_amber_forces_multiprmtop(global_options, parm_datas, coords_data);
            ++function_calls;
          }
          else
          {
            printf("*** FUNCTION %d IS NOT YET IMPLEMENTED\n",global_options->FUNC_TO_FIT);
            return(NOT_IMPLEMENTED);
          }

          /*Test the return of the function and see if it is lower than the reflect_point -> I.e are we heading in the corect
            direction?*/
          if(reflect_point_func>extended_reflect_point_func)
          {
            /*We have an even better estimation of our minimum*/
            /*Replace the highest set of parameters in our vertex matrix with the new extended point = now the lowest point*/
            for (col=0;col<parm_datas[0].ndimensions;++col)
            {
              vertex_matrix[func_hi][col]=extended_reflect_point[col];
            }
            /*Also update the function value to save us having to do yet another function evaluation*/
            function_results[func_hi]=extended_reflect_point_func;
          } /*End if(reflect_point_func>extended_reflect_point_func)*/
          else
          { /*the extended point took us away from the minimum*/
            for(col=0;col<parm_datas[0].ndimensions;++col)
            {
              /*Stick with our reflected point, replace the highest set of parameters in our vertex matrix with the new reflected
                point parameters = The lowest found point to date*/
              vertex_matrix[func_hi][col]=reflect_point[col];
            }
            function_results[func_hi]=reflect_point_func;
          } /*End else of if(reflect_point_func>extended_reflect_point_func)*/
        } /* End of if(reflect_point_func<=(*(function_results+func_low))) Is our reflect_point lower than all the other points?*/
        else  /*This is from the if command above after we calculated our reflected point - here this means our reflected point*/
        {    /*is not lower than the lowest point we have found so far*/
          /*It could be higher than all the points we have, including the original, or it could be somewhere in
            between the existing points. Hence we test each condition and react as appropriate*/
          if( reflect_point_func>=function_results[func_next_hi] ) /*Test if it is between our current highest point (the one that was reflected) and the next highest point*/
          {
            if( reflect_point_func<function_results[func_hi] )
            { /*it is between the highest point and the next highest*/
              /*Record it as being our new highest point*/
              for(col=0;col<parm_datas[0].ndimensions;++col)
              {
                vertex_matrix[func_hi][col]=reflect_point[col];
              }
              /*Also update our function results to avoid an extra function evaluation*/
              function_results[func_hi]=reflect_point_func;
            }
            for(col=0;col<parm_datas[0].ndimensions;++col) /*Test the function evaluation point of the a reduced value of the reflection point and*/
            {                                /*see where this lies*/
              extended_reflect_point[col]=reflection_reduction_ratio*vertex_matrix[func_hi][col]+(1.0-reflection_reduction_ratio)*(avg_point[col]);
              parameters[col]=extended_reflect_point[col];
            }
            /*Need to update the parm_data arrays*/
            retval=update_prmtop_data(global_options, parm_datas, parameters);
            process_retval(retval, global_options->VERBOSITY);
            /*Now evaluate our function with the extended point to see if it is lower than the original reflection*/
            if (global_options->FUNC_TO_FIT==SUM_SQUARES_AMBER_STANDARD)
            {
              extended_reflect_point_func = eval_sum_squares_amber_std_multiprmtop(global_options, parm_datas, coords_data, NULL);
              ++function_calls;
            }
            else if (global_options->FUNC_TO_FIT==AMBER_FORCES) {
              extended_reflect_point_func = eval_sum_amber_forces_multiprmtop(global_options, parm_datas, coords_data);
              ++function_calls;
            }
            else
            {
	      printf("*** FUNCTION %d IS NOT YET IMPLEMENTED\n",global_options->FUNC_TO_FIT);
              return(NOT_IMPLEMENTED);
            }
            if(extended_reflect_point_func<function_results[func_hi]) /*it now becomes our highest point since we have already*/
            {                                                            /*tested that it is higher than all the known points*/
              for(col=0;col<parm_datas[0].ndimensions;++col)
              {
                vertex_matrix[func_hi][col]=extended_reflect_point[col];
              }
              /*update the function results with the value for what is now our highest point*/
              function_results[func_hi]=extended_reflect_point_func;
            }
            else /*else of if(extended_reflect_point_func<function_results[func_hi])*/
            {  /*extended_reflect_point_func > the the highest function point we already have - therefore we moved away from the minimum*/
              for(col=0;col<parm_datas[0].ndimensions;++col)
              { /*pic parameters that are half way between the highest point and the lowest*/
                reflect_point[col]=0.5*(vertex_matrix[func_hi][col]+vertex_matrix[func_low][col]);
                parameters[col]=reflect_point[col];
              }
              /*Need to update the parm_data arrays*/
              retval=update_prmtop_data(global_options, parm_datas, parameters);
              process_retval(retval, global_options->VERBOSITY);
              /*Now evaluate our function with the extended point to see if it is lower than the original reflection*/
              if (global_options->FUNC_TO_FIT==SUM_SQUARES_AMBER_STANDARD)
              {
                reflect_point_func = eval_sum_squares_amber_std_multiprmtop(global_options, parm_datas, coords_data, NULL);
                ++function_calls;
              }
              else if (global_options->FUNC_TO_FIT==AMBER_FORCES) {
                reflect_point_func = eval_sum_amber_forces_multiprmtop(global_options, parm_datas, coords_data);
                ++function_calls;
              }
              else
              {
		printf("*** FUNCTION %d IS NOT YET IMPLEMENTED\n",global_options->FUNC_TO_FIT);
                return(NOT_IMPLEMENTED);
              }
              /*Now test the results and see if we are now less than our highest function to date*/
              if(reflect_point_func<function_results[func_hi])
              {
                for(col=0;col<parm_datas[0].ndimensions;++col)
                {
                  vertex_matrix[func_hi][col]=reflect_point[col];
                }
                /*Update the results and avoid another function evaluation*/
                function_results[func_hi]=reflect_point_func;
              }
              else /*else of if(reflect_point_func<function_results[func_hi])*/
              { /*the point was higher than all our known points to date*/
                for(col=0;col<parm_datas[0].ndimensions;++col)
                {
                  extended_reflect_point[col]=-vertex_matrix[func_hi][col]+2.0*vertex_matrix[func_low][col];
                  parameters[col]=extended_reflect_point[col];
                }
                /*Need to update the parm_data arrays*/
                retval=update_prmtop_data(global_options, parm_datas, parameters);
                process_retval(retval, global_options->VERBOSITY);
                /*Now evaluate our function with the extended point to see if it is lower than the original reflection*/
                if (global_options->FUNC_TO_FIT==SUM_SQUARES_AMBER_STANDARD)
                {
                  reflect_point_func = eval_sum_squares_amber_std_multiprmtop(global_options, parm_datas, coords_data, NULL);
                  ++function_calls;
                }
                else if (global_options->FUNC_TO_FIT==AMBER_FORCES) {
                  reflect_point_func = eval_sum_amber_forces_multiprmtop(global_options, parm_datas, coords_data);
                  ++function_calls;
                }
                else
                {
		  printf("*** FUNCTION %d IS NOT YET IMPLEMENTED\n",global_options->FUNC_TO_FIT);
                  return(NOT_IMPLEMENTED);
                }
                if(extended_reflect_point_func<function_results[func_hi])
                {
                  for(col=0;col<parm_datas[0].ndimensions;++col)
                  {
                    vertex_matrix[func_hi][col]=extended_reflect_point[col];
                  }
                  function_results[func_hi]=extended_reflect_point_func;
                }
                else /*Else of if(extended_reflect_point_func<function_results[func_hi])*/
                {
                  /*Original simplex method would at this point reduce itself in all dimensions - but this is
                    not very good for NDimensions>10 so instead we will test a series of other points. We check the function
                    in the point point1/point1' and if this yields a function lower than the highest point / the reflect_point the
                    highest point us changed to point to point1/point1'. If it is not lower the function is checked at point2/point2'
                    and if that is lower than the highest point/reflected point we change the highest point to point2/point2'.
                    Else - try something else... Expansion based on 3rd order polynomial
                  */
                  pointa=3.0*function_results[func_hi]-8.0*reflect_point_func+6.0*function_results[func_low]-extended_reflect_point_func;
                  pointb=function_results[func_hi]-2.0*function_results[func_low]+extended_reflect_point_func;
                  pointc=-0.5*function_results[func_hi]+8.0*(reflect_point_func/3.0)-2*function_results[func_low]+extended_reflect_point_func/6.0;
                  pointd=(pointb*pointb)-(4*pointa*pointc);

                  if(pointd>0.0) {
                    expMin=0.5*(-pointb-sqrt(pointd))/pointa;
                    for(col=0;col<parm_datas[0].ndimensions;++col) {
                      reflect_point[col]=expMin*vertex_matrix[func_hi][col]+(1-expMin)*vertex_matrix[func_low][col];
                      parameters[col]=reflect_point[col];
                    }
                    
                    retval=update_prmtop_data(global_options, parm_datas, parameters);
                    process_retval(retval, global_options->VERBOSITY);
                    
                    /*Now evaluate our function with the extended point to see if it is lower than the original reflection*/
                    if (global_options->FUNC_TO_FIT==SUM_SQUARES_AMBER_STANDARD) {
                      reflect_point_func = eval_sum_squares_amber_std_multiprmtop(global_options, parm_datas, coords_data, NULL);
                      ++function_calls;
                    } else if (global_options->FUNC_TO_FIT==AMBER_FORCES) {
                      reflect_point_func = eval_sum_amber_forces_multiprmtop(global_options, parm_datas, coords_data);
                      ++function_calls;
                    } else {
                      printf("*** FUNCTION %d IS NOT YET IMPLEMENTED\n",global_options->FUNC_TO_FIT);
                      return(NOT_IMPLEMENTED);
                    }
                  }
                  if(reflect_point_func<function_results[func_hi]) /*Is better than the highest point*/ {
                    for(col=0;col<parm_datas[0].ndimensions;++col)
                      vertex_matrix[func_hi][col]=reflect_point[col];
                    function_results[func_hi]=reflect_point_func;
                  } else {
                    for(col=0;col<parm_datas[0].ndimensions;++col)
                      vertex_matrix[func_hi][col]=vertex_matrix[func_low][col];
                    function_results[func_hi]=function_results[func_low];
                  }
                } /*End of else of if(extended_reflect_point_func<function_results[func_hi])*/
              } /*End of else of if(reflect_point_func<function_results[func_hi])*/
            } /*End of else of if(extended_reflect_point_func<function_results[func_hi])*/
          } /*End of if( reflect_point_func>=function_results[func_next_hi] )*/
          else /*else of if( reflect_point_func>=function_results[func_next_hi] ) Test if it is between our current highest point (the one that was reflected) and the next highest point*/
          {
            for(col=0;col<parm_datas[0].ndimensions;++col)
            {
              vertex_matrix[func_hi][col]=reflect_point[col];
            }
            function_results[func_hi]=reflect_point_func;
          } /*End of else of if(reflect_point_func>=(*(function_results+func_next_hi))) Test if it is between our current highest point (the one that was reflected) and the next highest point*/
        } /*End of else of if(reflect_point_func<=function_results[func_low])*/
      } /*End of inner loop -> inner_itr_count=1;inner_itr_count<=inner_itr*/
/*END INNER LOOP*/      

      /*Check to see if our minimum value has changed since the last loop*/
      if (previous_min==function_results[func_low])
        ++min_unchanged_counter;
      else
        min_unchanged_counter=0;

      previous_min=function_results[func_low];

      /*now we have finished our inner simplex loop we need to test for convergence
        and stop if necessary or continue for another outer loop by restarting with
        modified input parameters (adjusted by dx)
      */

      /*Note, this convergence routine will break down if the function is very close to zero - i.e.
        the sum of the func_hi and func_low values are close to zero.
        Therefore we also need to test if the func_hi value - func_low value is close zero*/

      /*At this point it is possible that our func_hi and func_low data may be scrambled - for example func_hi may actually contain
        the lowest value - Thus we do a resort here to get things correct. This is the same sort as is done at the beginning of
        every inner loop*/
      /*BEGIN SORT*/
      func_low=0;
      if (function_results[0]>function_results[1])
      {
        /*Our first set of parameters give a higher result than the second set*/
        func_hi=0;
        func_next_hi=1;
      }
      else
      {
        func_hi=1;
        func_next_hi=0;
      }
      for (row=0;row<=parm_datas[0].ndimensions;++row) /*Note we loop over a total of NDIMENSIONS+1 here*/
      {
         /*loop over each row = each parameter set*/
         if ( function_results[row] < function_results[func_low] ) /*Test if function with current row params is less than the lowest found so far*/
         {
           func_low=row; /*If it is lower change our low locator to this row*/
         }
         if ( function_results[row] > function_results[func_hi] ) /*Repeat test to see if it is higher than our current highest*/
         {
           /*If it is higher then update the highest locator and the next highest locator*/
           func_next_hi=func_hi;
           func_hi=row;
         }
         else if ( function_results[row] > function_results[func_next_hi] )
         {
            /*We are between highest and second highest*/
            /*only update if we aren't currently at the highest value*/
           if (row!=func_hi)
           {
             func_next_hi=row;
           }
         }
      } /*End of loop over rows (parameter sets) for (row=0;row<=parm_datas[0].ndimensions;++row)*/
      /*We should now have func_low containing the lowest value and func_hi the highest*/
      conv_ratio=2.0*fabs((function_results[func_hi]-function_results[func_low]))/(fabs(function_results[func_hi])+fabs(function_results[func_low]));
      function_results_sum=0.0;
      for (row=0;row<=parm_datas[0].ndimensions;++row) /*Note we loop over a total of NDIMENSIONS+1 here*/
      {
        function_results_sum+=function_results[row]; /*used for printing the average over all the functions with the different params, should converge*/
      }
      if (global_options->VERBOSITY>=MEDIUM)
      {
        if (global_options->VERBOSITY>=HIGH) /*If debug is on print to higher precision */
          {
            printf("   Step %5d: Conv=%16.14E min=%16.14f max=%16.14f avg=%16.14f\n",outer_itr_count+1,
            conv_ratio,function_results[func_low],function_results[func_hi],function_results_sum/(parm_datas[0].ndimensions+1));
          }
          else
          {
            printf("   Step %5d: Conv=%12.4E min=%12.4f max=%12.4f avg=%12.4f\n",outer_itr_count+1,
            conv_ratio,function_results[func_low],function_results[func_hi],function_results_sum/(parm_datas[0].ndimensions+1));
          }
          fflush(stdout); /*Flush the printf buffer*/
      }
      /*Now we test to see if we have converged and if we have we break the outer loop*/
      
  if (conv_ratio<=global_options->CONV_LIMIT) {
	
    if(global_options->VERBOSITY>=MEDIUM) {
      printf("   -----------------------------------------------------------------------------------\n");
      printf("\n   Convergence ratio of %12.4E is better than\n   convergence criteria of %12.4E.\n",conv_ratio,global_options->CONV_LIMIT);
      printf("   Function Converged - Total function evaluations = %d\n\n",function_calls);
      fflush(stdout); /*Flush the printf buffer*/
    }
    loop_exit_condition=SUCCESS;
    break;
  }   /*We also need to check if our function is very close
        to zero - if it is the convergence equation above breaks down - i.e our fit is exceptionally good*/
  else if (fabs(function_results[func_low])<=global_options->CONV_LIMIT) {
    if(global_options->VERBOSITY>=MEDIUM) {
      printf("   -----------------------------------------------------------------------------------\n");
      printf("\n   Convergence of %12.4E in function value is better than\n   convergence criteria of %12.4E.\n",fabs(function_results[func_low]),global_options->CONV_LIMIT);
      printf("   Function Converged - Total function evaluations = %d\n\n",function_calls);
      fflush(stdout); /*Flush the printf buffer*/
    }
    loop_exit_condition=SUCCESS;
    break;
  }
  else if (min_unchanged_counter>=100) { /*Our minimum hasn't changed in over 100 steps so abort*/
    if(global_options->VERBOSITY>=MEDIUM) {
      printf("   -----------------------------------------------------------------------------------\n");
      printf("\n   Minimum function value of %14.6f has not changed in %d cycles\n",previous_min,min_unchanged_counter);
      printf("   Assuming Function Has Converged - Total function evaluations = %d\n\n",function_calls);
      fflush(stdout); /*Flush the printf buffer*/
    }
    loop_exit_condition=MINSTATIC;
    break;
  }
      /*If we got to here we haven't met our convergence criteria so we run through it all again*/
      func_low=0; /*Find out what our new lowest value now is*/
      for(row=0;row<=parm_datas[0].ndimensions;++row)
      {
        if(function_results[row]<function_results[func_low])
        {
          func_low=row;
        }
      }
      for(col=0;col<parm_datas[0].ndimensions;col++)
      {
        parameters[col]=vertex_matrix[func_low][col];
      }
      /*Write the new "better" parameters to the parm_data array*/
      retval=update_prmtop_data(global_options, parm_datas, parameters);
      process_retval(retval, global_options->VERBOSITY);
      ++outer_itr_count;      /*Increment the restart loop counter*/
    } /* End of outer simplex loop while (outer_itr_count<outer_itr_max)
         Note: We will have either exited this loop via a break statement or because we exceeded our max iterations.
         Either way we should examine the loop exit code
      */
/*END SIMPLEX OUTER LOOP*/

    /*Parameters scratch array now contains the optimised parameters*/
    /*Copy it to the parm data*/
    /*Write the new "better" parameters to the parm_data array*/
   retval=update_prmtop_data(global_options, parm_datas, parameters);
   process_retval(retval, global_options->VERBOSITY);
    if (loop_exit_condition==EXCEEDEDMAXITERATIONS)
    {
       /*We quit because we did not converge in the maximum number of iterations defined by outer_itr_max*/
       /*print an informational message for the user*/
	printf("*** Aborting run due to convergence failure.\n");
	printf("*** Exceeded outer_itr_max of %d\n",outer_itr_max);
	printf("*** Try rerunning the minimisation with a different start point and/or increase the\n");
	printf("*** the maximum number of iterations allowed and/or relax the convergence criteria.\n");
	printf("*** PARAMETERS FOR FINAL STEP FOLLOW:\n");
	fflush(stdout); /*Flush the printf buffer*/
    }
    else if (loop_exit_condition==MINSTATIC)
    {
      /*We quit because our minimum function was no longer changing*/
      printf("!  Warning - Convergence criteria of %10.4E was not met.\n",global_options->CONV_LIMIT);
      printf("!            You should check the converged parameters carefully.\n\n");
      printf("!  PARAMETERS FOR FINAL STEP FOLLOW:\n");
      fflush(stdout); /*Flush the printf buffer*/
    }
    else if (loop_exit_condition==SUCCESS)
    {
      /*We successfully converged*/
      if (global_options->VERBOSITY>=MEDIUM)
      {
        printf("   Convergence to %10.4E in Simplex routine achieved after %d cycles.\n",global_options->CONV_LIMIT,(outer_itr_count+1)*inner_itr);
        printf("   (%d INNER x %d OUTER CYCLES)\n\n",inner_itr,outer_itr_count+1); /*+1 since it counts from zero*/
        fflush(stdout); /*Flush the printf buffer*/
      }
    }
    else
    {
       /*Unknown exit condition*/
       printf("*** UNKNOWN EXIT CONDITION %d ENCOUNTERED IN minimise_function_simplex().\n",loop_exit_condition);
       printf("*** CONDITION - FATAL\n");
       fflush(stdout); /*Flush the printf buffer*/
       return(FAILURE);
    }
        
    free(parameters);
    free(extended_reflect_point);
    free(reflect_point);
    free(avg_point);
    free(function_results);
    double_2D_array_free(vertex_matrix);
    
    parameters = NULL;
    extended_reflect_point = NULL;
    reflect_point = NULL;
    avg_point = NULL;
    function_results = NULL;

    return(loop_exit_condition);

}

