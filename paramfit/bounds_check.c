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

/** @file bounds_check.c
 * Contains functions that ensure the minimization algorithms have
 * not wandered off into areas where there are too little data in the
 * input structures.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constants.h"

#include "function_def.h"

/**
 * Collects up information about bonds, angles, and dihedrals in the input
 * structures in an easy to use data structure.
 * @param[in] global_options The global options structure
 * @param[out] bounds_data The data structure with the data about input structure distributions
 * @param[in] parm_data Pointer to a single parameter data structure
 * @param[in] coords_data Pointer to a single input coordinates data structure
 * @return Integer indicating success or failure
 */
int calculate_structure_diversity(global_options_struct *global_options, bounds_struct *bounds_data, parm_struct *parm_data, coord_set *coords_data)
{
  int structure=0, count=0;
  int i, j;
  double tempx1, tempy1, tempz1, tempx2, tempy2, tempz2;
  double tempx3, tempy3, tempz3, tempx4, tempy4, tempz4;  
  
  /* Allocate necessary memory for the arrays */
  /* Each array has one entry for a unique bond, angle, or dihedral type. *
   * This entry is an array of every value for this type over all of the  *
   * input structures. This prevents averaging the values for each        *
   * structure and instead storing all of them.                           */
  bounds_data->bond_lengths=(double**)malloc(parm_data->unique_bonds_found*sizeof(double*));
  if (bounds_data->bond_lengths == NULL) {
    malloc_failure_double("calculate_structure_diversity", "bond_lengths", parm_data->unique_bonds_found);
    return ALLOC_FAIL;
  }
  
  bounds_data->angle_thetas=(double**)malloc(parm_data->unique_angles_found*sizeof(double*));
  if (bounds_data->angle_thetas==NULL) {
    malloc_failure_double("calculate_structure_diversity", "angle_thetas", parm_data->unique_angles_found);
    return ALLOC_FAIL;
  }
  
  bounds_data->dihedral_thetas=(double**)malloc(parm_data->unique_dihedrals_found*sizeof(double*));
  if (bounds_data->dihedral_thetas==NULL) {
    malloc_failure_double("calculate_structure_diversity", "dihedral_thetas", parm_data->unique_dihedrals_found);
    return ALLOC_FAIL;
  }
  bounds_data->mem_allocated += (parm_data->unique_bonds_found+parm_data->unique_angles_found+parm_data->unique_dihedrals_found)*sizeof(double*);
 
  // Calculate all available bond lengths and put them in the array
  for (i=0; i<parm_data->unique_bonds_found; ++i) {
    // Allocate memory for all of this bond- each structure has the same number of this type of bond so we can allocate the whole array now.
    bounds_data->bond_lengths[i] = (double*)malloc(coords_data->num_coords*parm_data->bond_data[i].number*sizeof(double));
    if (bounds_data->bond_lengths[i] == NULL) {
      malloc_failure_double("calculate_structure_diversity", "bounds_data->bond_lengths[i]", coords_data->num_coords*parm_data->bond_data[i].number);
      return ALLOC_FAIL;
    }
    bounds_data->mem_allocated += coords_data->num_coords*parm_data->bond_data[i].number*sizeof(double);
   
    // Fill in the array, one structure at a time. Use count to keep track of array index.
    count=0;
    for (structure=0; structure<coords_data->num_coords; ++structure) {
      for (j=0; j<parm_data->bond_data[i].number; ++j) {
         tempx1=coords_data->struc[structure].x_coord[parm_data->bond_data[i].atom1[j]-1];
        tempy1=coords_data->struc[structure].y_coord[parm_data->bond_data[i].atom1[j]-1];
        tempz1=coords_data->struc[structure].z_coord[parm_data->bond_data[i].atom1[j]-1];
        tempx2=coords_data->struc[structure].x_coord[parm_data->bond_data[i].atom2[j]-1];
        tempy2=coords_data->struc[structure].y_coord[parm_data->bond_data[i].atom2[j]-1];
        tempz2=coords_data->struc[structure].z_coord[parm_data->bond_data[i].atom2[j]-1];
        
        bounds_data->bond_lengths[i][count] = calc_bond_length(tempx1,tempy1,tempz1,tempx2,tempy2,tempz2);
        ++count;
    } }
  }
  
  // Do the same for the angles
  for (i=0; i<parm_data->unique_angles_found; ++i) {
    // Allocate memory
    bounds_data->angle_thetas[i] = (double*)malloc(coords_data->num_coords*parm_data->angle_data[i].number*sizeof(double));
    if (bounds_data->angle_thetas[i]==NULL) {
      malloc_failure_double("calculate_structure_diversity", "bounds_data->angle_thetas[i]", coords_data->num_coords*parm_data->angle_data[i].number);
      return ALLOC_FAIL;
    }
    bounds_data->mem_allocated += coords_data->num_coords*parm_data->angle_data[i].number*sizeof(double);
    
    // Fill in array
    count=0;
    for (structure=0; structure < coords_data->num_coords; ++structure) {
      for (j=0;j<parm_data->angle_data[i].number;++j) {
        tempx1=coords_data->struc[structure].x_coord[parm_data->angle_data[i].atom1[j]-1];
        tempy1=coords_data->struc[structure].y_coord[parm_data->angle_data[i].atom1[j]-1];
        tempz1=coords_data->struc[structure].z_coord[parm_data->angle_data[i].atom1[j]-1];
        tempx2=coords_data->struc[structure].x_coord[parm_data->angle_data[i].atom2[j]-1];
        tempy2=coords_data->struc[structure].y_coord[parm_data->angle_data[i].atom2[j]-1];
        tempz2=coords_data->struc[structure].z_coord[parm_data->angle_data[i].atom2[j]-1];
        tempx3=coords_data->struc[structure].x_coord[parm_data->angle_data[i].atom3[j]-1];
        tempy3=coords_data->struc[structure].y_coord[parm_data->angle_data[i].atom3[j]-1];
        tempz3=coords_data->struc[structure].z_coord[parm_data->angle_data[i].atom3[j]-1];

        /*find angle*/
        double rad = calc_angle_radians(tempx1,tempy1,tempz1,tempx2,tempy2,tempz2,tempx3,tempy3,tempz3);
        if (rad != rad) { // check for NaN error
          printf("ERROR: NaN error:\n");
          printf("Unique angle %i, structure %i, mult %i\n", i, structure, j);
          printf("(%f,%f,%f) (%f,%f,%f) (%f,%f,%f) = %f\n", tempx1,tempy1,tempz1,tempx2,tempy2,tempz2,tempx3,tempy3,tempz3, rad);
          return ABORT;
        }
        bounds_data->angle_thetas[i][count] = rad;
        
        ++count;
    } }
  }
  
  // Finally, the dihedrals
  for (i=0; i<parm_data->unique_dihedrals_found; ++i) {
    // We ignore the terms, since that refers to parameter usage. We just want all the dihedral angles in each structure.
    bounds_data->dihedral_thetas[i] = (double*)malloc(coords_data->num_coords*parm_data->dihedral_data[i].number*sizeof(double));
    if (bounds_data->dihedral_thetas[i] == NULL) {
      malloc_failure_double("calculate_structure_diversity", "bounds_data->dihedral_thetas[i]", 
                            parm_data->dihedral_data[i].number*coords_data->num_coords);
      return ALLOC_FAIL;
    }
    bounds_data->mem_allocated += coords_data->num_coords*parm_data->dihedral_data[i].number*sizeof(double);
     
    // Fill in array
    count=0;
    for (structure=0; structure < coords_data->num_coords; ++structure) {
      for (j=0; j<parm_data->dihedral_data[i].number; ++j) {
        /*First off we need to find the dihedral between the atoms in the angle*/
        /*Get the coordinates*/
        tempx1=coords_data->struc[structure].x_coord[parm_data->dihedral_data[i].atom1[j]-1];
        tempy1=coords_data->struc[structure].y_coord[parm_data->dihedral_data[i].atom1[j]-1];
        tempz1=coords_data->struc[structure].z_coord[parm_data->dihedral_data[i].atom1[j]-1];
        tempx2=coords_data->struc[structure].x_coord[parm_data->dihedral_data[i].atom2[j]-1];
        tempy2=coords_data->struc[structure].y_coord[parm_data->dihedral_data[i].atom2[j]-1];
        tempz2=coords_data->struc[structure].z_coord[parm_data->dihedral_data[i].atom2[j]-1];
        tempx3=coords_data->struc[structure].x_coord[parm_data->dihedral_data[i].atom3[j]-1];
        tempy3=coords_data->struc[structure].y_coord[parm_data->dihedral_data[i].atom3[j]-1];
        tempz3=coords_data->struc[structure].z_coord[parm_data->dihedral_data[i].atom3[j]-1];
        tempx4=coords_data->struc[structure].x_coord[parm_data->dihedral_data[i].atom4[j]-1];
        tempy4=coords_data->struc[structure].y_coord[parm_data->dihedral_data[i].atom4[j]-1];
        tempz4=coords_data->struc[structure].z_coord[parm_data->dihedral_data[i].atom4[j]-1];
          
        /*find dihedral*/
        bounds_data->dihedral_thetas[i][count] = PI-calc_dihedral_radians(tempx1,tempy1,tempz1,tempx2,tempy2,tempz2,tempx3,tempy3,tempz3,tempx4,tempy4,tempz4);

        ++count;
    } }
  }

  // Now sort each of the arrays in ascending order
  // Sort each type of bond 
  for (i=0; i<parm_data->unique_bonds_found; ++i) {
    for (count=1; count<coords_data->num_coords*parm_data->bond_data[i].number; ++count) {
      double temp = bounds_data->bond_lengths[i][count];
      
      j = count-1;
      while(j >=0 && bounds_data->bond_lengths[i][j] > temp) {
        bounds_data->bond_lengths[i][j+1] = bounds_data->bond_lengths[i][j];
        --j;
      }
      bounds_data->bond_lengths[i][j+1] = temp;
    }
  }
  // Angles
  for (i=0; i<parm_data->unique_angles_found; ++i) {
    for (count=1; count<coords_data->num_coords*parm_data->angle_data[i].number; ++count) {
      double temp = bounds_data->angle_thetas[i][count];
      
      j = count-1;
      while(j >=0 && bounds_data->angle_thetas[i][j] > temp) {
        bounds_data->angle_thetas[i][j+1] = bounds_data->angle_thetas[i][j];
        --j;
      }
      bounds_data->angle_thetas[i][j+1] = temp;
    }
  }
  // Dihedrals
  for (i=0; i<parm_data->unique_dihedrals_found; ++i) {
    for (count=1; count<coords_data->num_coords*parm_data->dihedral_data[i].number; ++count) {
      double temp = bounds_data->dihedral_thetas[i][count];
      
      j = count-1;
      while(j >=0 && bounds_data->dihedral_thetas[i][j] > temp) {
        bounds_data->dihedral_thetas[i][j+1] = bounds_data->dihedral_thetas[i][j];
        --j;
      }
      bounds_data->dihedral_thetas[i][j+1] = temp;
    }
  }
  
 return SUCCESS;
}

/** 
 * Deletes the table of bond, angle, and dihedral structure data for clean up.
 * @param[in,out] bounds_data The structure from which to delete the tables
 * @param[in] parm_data The parameter data structure, will be used to determine array dimensions to delete.
 */
void clean_up_bounds(bounds_struct *bounds_data, parm_struct *parm_data)
{
  int i;
  for (i=0; i<parm_data->unique_bonds_found; ++i)
    free(bounds_data->bond_lengths[i]);
  for (i=0; i<parm_data->unique_angles_found; ++i)
    free(bounds_data->angle_thetas[i]);
  for (i=0; i<parm_data->unique_dihedrals_found; ++i)
    free(bounds_data->dihedral_thetas[i]);
  free(bounds_data->bond_lengths);
  free(bounds_data->angle_thetas);
  free(bounds_data->dihedral_thetas);
}

/**
 * Checks the dihedral equilibrium phases are well represented in the input structures
 * @see calculate_structure_diversity
 * @param[in] global_options The global options structure
 * @param[in] parm_data Pointer to single parameter struture containing the dihedral parameters
 * @param[in] coords_data Pointer to single coordinate set
 * @param[in] bounds_data The table of bond, angle, and dihedral values from the structures
 * @return SUCCESS and a warning if out of bounds for error or ignore checking options, else FAILURE
 */
int check_dihedrals(global_options_struct *global_options, parm_struct *parm_data, coord_set *coords_data, bounds_struct *bounds_data)
{
  double increment = PI/global_options->DIHEDRAL_SPAN;     // there should be a dihedral value every this many rads
  double rad=0;
  double invalid_min=-1.0, invalid_max=-1.0;               // to give information about where data is missing.
  int dihedral=0; int term=0;
  int num_missing=0;
  short valid=YES; short go=NO;
  int theta = 0;        // keeps track of unique dihedrals
  
  for (theta=0; theta<parm_data->unique_dihedrals_found; ++theta) // check each dihedral is well represented
  {
    // Check if this dihedral is being altered
    go=NO;
    for(term=0; term<parm_data->dihedral_data[theta].num_terms; ++term) {
      if (parm_data->dihedral_data[theta].term[term].DO_FIT_PHASE==YES ||
          parm_data->dihedral_data[theta].term[term].DO_FIT_KP==YES ){ // only check ones that are being fitted
        go=YES;
        break;
    } }
    if (go==YES) {  
      // Create scatter plot file if desired
      if (global_options->SCATTERPLOTS==TRUE)
        {
        FILE *plot;
        char filename[10];
        sprintf(filename, "%d.diheq", theta);
        printf("  Writing dihedral scatter plot %s\n", filename);
        plot = fopen(filename, "w");
        fprintf(plot, "#%s-%s-%s-%s\n", parm_data->dihedral_data[theta].atom_type1, parm_data->dihedral_data[theta].atom_type2,
                parm_data->dihedral_data[theta].atom_type3, parm_data->dihedral_data[theta].atom_type4);
        for (dihedral=0; dihedral<coords_data->num_coords; ++dihedral) {
          fprintf(plot, "%f  1\n", bounds_data->dihedral_thetas[theta][dihedral]);
        }
        fclose(plot);
      }
      rad = 0;
      dihedral = 0;
      valid = YES;
      num_missing = 0;
      
      while (rad+increment <= PI) {
        dihedral = 0;
        valid = YES;
        int done = NO;
        // attempt to find a dihedral with value in this range
        while ( dihedral < coords_data->num_coords && done==NO ) {
          if ( (bounds_data->dihedral_thetas[theta][dihedral]> rad) && (bounds_data->dihedral_thetas[theta][dihedral] < rad+increment) )
            done = YES;
          else
            ++dihedral;
        }
        
        // check if a structure could be found
        if (dihedral==coords_data->num_coords) {
          if (invalid_min==-1.0)
            invalid_min = rad;
          invalid_max = rad+increment;
          ++num_missing;
          valid = NO;
          
          if (global_options->VERBOSITY>=HIGH)
            printf("%s-%s-%s-%s dihedral FAILED span check.  No data in range %.4f - %.4f rads.\n",  parm_data->dihedral_data[theta].atom_type1,
                                                                                                        parm_data->dihedral_data[theta].atom_type2,
                                                                                                        parm_data->dihedral_data[theta].atom_type3,
                                                                                                        parm_data->dihedral_data[theta].atom_type4,
                                                                                                        rad, rad+increment);
        }
        rad+=increment;
      }
      
      // print out results if it passed
      if (valid==YES && global_options->VERBOSITY>=HIGH)
        printf("%s-%s-%s-%s dihedral PASSED span check for range %.4f - %.4f rads.\n",  parm_data->dihedral_data[theta].atom_type1,
               parm_data->dihedral_data[theta].atom_type2,
               parm_data->dihedral_data[theta].atom_type3,
               parm_data->dihedral_data[theta].atom_type4,
               rad, rad+increment);
               
      if (valid==NO && global_options->VERBOSITY>=MEDIUM)
      {
        if (global_options->CHECK_BOUNDS==YES)
          printf("   ERROR: ");
        else
          printf("   WARNING: ");
        printf("%s-%s-%s-%s dihedral is missing %.i data points in the range %.4f to %.4f radians.\n",  parm_data->dihedral_data[theta].atom_type1,
                                                                                                           parm_data->dihedral_data[theta].atom_type2,
                                                                                                           parm_data->dihedral_data[theta].atom_type3,
                                                                                                           parm_data->dihedral_data[theta].atom_type4,
                                                                                                           num_missing, invalid_min, invalid_max);
      }
    }
  }
  
  if (valid==YES)
  {
    if (global_options->VERBOSITY >= MEDIUM)
      printf(" * Input structures passed dihedral span check.\n");
    return SUCCESS;
  }
  else
  {
    printf("\n\n");
    if (global_options->CHECK_BOUNDS==YES)
      printf("   ERROR: ");
    else
      printf("   WARNING: ");
    printf("Insufficient dihedral information in sample structures.\n");
    printf("            Your settings require at least %i samples with data    \n", global_options->DIHEDRAL_SPAN);
    printf("            at least every %.3f radians (%.2f degrees).            \n", increment, increment*RADIAN_TO_DEGREE);
    printf("            Either 1) Add the missing input data or                \n");
    printf("                   2) Set DIHEDRAL_SPAN to a smaller value or      \n");
    printf("                   3) Set BOUNDS_CHECK to warn (not recommended).  \n");
    printf("            Please read the help and/or documentation.             \n");
    printf("\n");
    fflush(stdout);
    if (global_options->CHECK_BOUNDS==YES)
      return FAILURE;
    else
      return SUCCESS;
  }
}

/**
 * Checks the angle equilibrium value is well defined in the input structures.
 * This means that the value must be within ANGLE_LIMIT of an input structure angle value. This function
 * is intended to be run with a converged set of parameters.
 * @see calculate_structure_diversity
 * @param[in] global_options The global options structure
 * @param[in] parm_data Pointer to the parameter data structure, including the angle value to check
 * @param[in] coords_data Pointer to the coordinate set for these parameters
 * @param[in] bounds_data The table of bonds, angles, and dihedrals in the input structures
 * @return SUCCESS and a warning if out of bounds with check set to warn or ignore, else FAILURE
 */
int check_angles(global_options_struct *global_options, parm_struct *parm_data, coord_set *coords_data, bounds_struct *bounds_data)
{
  // ensure that there are angles around the existing minimum
  int param, pos;
  double diff_next, diff_prev;
  short passed=YES;
  
  for (param=0; param<parm_data->unique_angles_found; ++param)
  {
    // only check angles that are being optimised
    if (parm_data->angle_data[param].DO_FIT_THEQ==YES)
    {
      // Create scatter plot file if desired
      FILE *plot;
      if (global_options->SCATTERPLOTS==TRUE) {
        char filename[10];
        sprintf(filename, "%d.angleq", param);
        printf("  Writing angle scatter plot %s\n", filename);
        plot = fopen(filename, "w");
        fprintf(plot, "#%s-%s-%s\n", parm_data->angle_data[param].atom_type1, parm_data->angle_data[param].atom_type2, parm_data->angle_data[param].atom_type3);
        for (pos=0; pos<coords_data->num_coords; ++pos) {
          fprintf(plot, "%f  1\n", bounds_data->angle_thetas[param][pos]);
        }
        fclose(plot);
      }
      
      // adjust phase to be between 0 and PI if necessary
      if (parm_data->angle_data[param].teq > PI)
        parm_data->angle_data[param].teq = fmod(parm_data->angle_data[param].teq, PI);
      else if (parm_data->angle_data[param].teq < 0.0)
        parm_data->angle_data[param].teq = PI - fmod(parm_data->angle_data[param].teq, PI);
      
      // find position in parameter array
      pos=0;
      while (pos < coords_data->num_coords && parm_data->angle_data[param].teq > bounds_data->angle_thetas[param][pos])
        ++pos;
      
      // find distance to neighbors in array
      if (pos==0) {
        diff_next = bounds_data->angle_thetas[param][pos] - parm_data->angle_data[param].teq;
        diff_prev = -1.0;
      } else if (pos==coords_data->num_coords) {
        diff_next = -1.0;
        diff_prev = parm_data->angle_data[param].teq - bounds_data->angle_thetas[param][pos-1];
      } else {
        diff_next = bounds_data->angle_thetas[param][pos] - parm_data->angle_data[param].teq;
        diff_prev = parm_data->angle_data[param].teq - bounds_data->angle_thetas[param][pos-1];
      }
      
      if (diff_next > global_options->ANGLE_LIMIT || diff_prev > global_options->ANGLE_LIMIT || diff_next < 0.0 || diff_prev < 0.0) {
        passed = NO;
        if (global_options->VERBOSITY>=MEDIUM)
        {
          if (global_options->CHECK_BOUNDS==YES)
            printf("   ERROR: ");
          else
            printf("   WARNING: ");
          if (diff_next == -1.0)
            printf("%s-%s-%s angle has no sample structures after it, nearest sample is %.4f radians before.\n",  
                   parm_data->angle_data[param].atom_type1,
                   parm_data->angle_data[param].atom_type2,
                   parm_data->angle_data[param].atom_type3,
                   diff_prev);
          else if (diff_prev == -1.0)
            printf("%s-%s-%s angle has no sample structures before it, nearest sample is %.4f radians after.\n", 
                   parm_data->angle_data[param].atom_type1,
                   parm_data->angle_data[param].atom_type2,
                   parm_data->angle_data[param].atom_type3,
                   diff_next);
          else
            printf("%s-%s-%s angle has no sample structures within %.4f radians, nearest are %.4f and %.4f radians away.\n",  
                   parm_data->angle_data[param].atom_type1,
                   parm_data->angle_data[param].atom_type2,
                   parm_data->angle_data[param].atom_type3,
                   global_options->ANGLE_LIMIT, diff_prev, diff_next);
        }
      }
      else
        if (global_options->VERBOSITY>=HIGH)
          printf("%s-%s-%s angle PASSED: data exist within %.4f and %.4f radians, passes limit of %.4f.\n",  parm_data->angle_data[param].atom_type1,
                 parm_data->angle_data[param].atom_type2,
                 parm_data->angle_data[param].atom_type3,
                 diff_prev, diff_next, global_options->ANGLE_LIMIT);
    }
  }
  
  if (passed==YES)
  {
    if (global_options->VERBOSITY >= MEDIUM)
      printf(" * Result passed angle validity check.\n");
    return SUCCESS;
  }
  else
  {
    printf("\n\n");
    if (global_options->CHECK_BOUNDS==YES)
      printf("   ERROR: ");
    else
      printf("   WARNING: ");
    printf("Insufficient angle information in sample structures.   \n");
    printf("            Your settings require that sample data exist within    \n");
    printf("            %.3f radians of the algorithm's result.                \n", global_options->ANGLE_LIMIT);
    printf("            Either 1) Add the missing input data                 or\n");
    printf("                   2) Set ANGLE_LIMIT to a larger value          or\n");
    printf("                   3) Set BOUNDS_CHECK to warn (not recommended).  \n");
    printf("            Please read the help and/or documentation.             \n");
    printf("\n");
    fflush(stdout);
    if (global_options->CHECK_BOUNDS==YES)
      return FAILURE;
    else
      return SUCCESS;
  }
  
}

/**
 * Checks if the bond parameters are adequately represented in the input structures.
 * The generated Keq should be within global_options->BOND_LIMIT of an input structure
 * bond equilibrium distance. Returns SUCCESS or FAILURE.
 * 
 * @see calculate_structure_diversity
 * @param[in] global_options The global options structure
 * @param[in] parm_data Pointer to the parameter structure with the bond parameters to check inside
 * @param[in] coords_data Pointer to coordinate structure that bounds data was gathered from
 * @param[in] bounds_data Table of bond, angle, and dihedral data in input structures
 * @return SUCCESS and prints a warning if bounds checking is set to ignore or warn, else FAILURE
 */
int check_bonds(global_options_struct *global_options, parm_struct *parm_data, coord_set *coords_data, bounds_struct *bounds_data)
{
  // ensure that there are bonds around the existing minimum
  int param, pos;
  double diff_next, diff_prev;
  short passed=YES;
  FILE *plot;
  
  for (param=0; param<parm_data->unique_bonds_found; ++param)
  {
    // only check bonds that are being optimised
    if (parm_data->bond_data[param].DO_FIT_REQ==YES)
    {
      // Create scatter plot file if desired
      if (global_options->SCATTERPLOTS==TRUE) {
        char filename[10];
        sprintf(filename, "%d.bondeq\n", param);
        plot = fopen(filename, "w");
        fprintf(plot, filename, "%s-%s", parm_data->bond_data[param].atom_type1, parm_data->bond_data[param].atom_type2);
        for (pos=0; pos<coords_data->num_coords; ++pos) {
          fprintf(plot, "%f\n", bounds_data->bond_lengths[param][pos]);
        }
        fclose(plot);
      }
      
      // find position in parameter array
      pos=0;
      while (pos < coords_data->num_coords && parm_data->bond_data[param].req > bounds_data->bond_lengths[param][pos])
        ++pos;
      // find distance to neighbors in array
      if (pos==0) {
        diff_next = bounds_data->bond_lengths[param][pos] - parm_data->bond_data[param].req;
        diff_prev = -1.0;
      } else if (pos==coords_data->num_coords) {
        diff_next = -1.0;
        diff_prev = parm_data->bond_data[param].req - bounds_data->bond_lengths[param][pos-1];
      } else {
        diff_next = bounds_data->bond_lengths[param][pos] - parm_data->bond_data[param].req;
        diff_prev = parm_data->bond_data[param].req - bounds_data->bond_lengths[param][pos-1];
      }
      
      if (diff_next > global_options->BOND_LIMIT || diff_prev > global_options->BOND_LIMIT || diff_next < 0.0 || diff_prev < 0.0) {
        passed = NO;
        if (global_options->VERBOSITY>=MEDIUM)
        {
          if (global_options->CHECK_BOUNDS==YES)
            printf("   ERROR: ");
          else
            printf("   WARNING: ");
          
          if (diff_next == -1.0)
            printf("%s-%s bond has no sample structures larger than it, nearest sample is %.4f A smaller.\n", parm_data->bond_data[param].atom_type1,
                                                                                                                          parm_data->bond_data[param].atom_type2,
                                                                                                                          diff_prev);
          else if (diff_prev == -1.0)
            printf("%s-%s bond has no sample structures smaller than it, nearest sample is %.4f A larger.\n", parm_data->bond_data[param].atom_type1,
                                                                                                                          parm_data->bond_data[param].atom_type2, 
                                                                                                                          diff_next);
          else
            printf("%s-%s bond has no sample structures within %f A, nearest are %f and %f A\n", parm_data->bond_data[param].atom_type1,
                                                                                                              parm_data->bond_data[param].atom_type2,
                                                                                                              global_options->BOND_LIMIT, diff_prev, diff_next);
        }
        
      }
      else
      {
        if (global_options->VERBOSITY>=HIGH)
          printf("PASSED %s-%s bond: Data exist within %f and %f A, passes limit of %f.\n", parm_data->bond_data[param].atom_type1,
                                                                                            parm_data->bond_data[param].atom_type2,
                                                                                            diff_prev, diff_next, global_options->BOND_LIMIT);
      }
    }
  }
  
  if (passed==YES)
  {
    if (global_options->VERBOSITY >= MEDIUM)
      printf(" * Result passed bond validity check.\n");
    return SUCCESS;
  }
  else
  {
    printf("\n\n");
    if (global_options->CHECK_BOUNDS==YES)
      printf("   mERROR: ");
    else
      printf("   WARNING: ");
    printf("Insufficient bond information in sample structures.    \n");
    printf("            Your settings require that sample data exist within    \n");
    printf("            %.3f A of the algorithm's result.                      \n", global_options->BOND_LIMIT);
    printf("            Either 1) Add the missing input data                 or\n");
    printf("                   2) Set BOND_LIMIT to a larger value           or\n");
    printf("                   3) Set BOUNDS_CHECK to warn (not recommended).  \n");
    printf("            Please read the help and/or documentation.             \n");
    printf("\n");
    fflush(stdout);
    if (global_options->CHECK_BOUNDS==YES)
      return FAILURE;
    else
      return SUCCESS;
  }
} 

/** Checks that the bonds, angles, and dihedrals are in a valid range before
 * conducting a function evaluation.
 *
 * If they are not, it will change any invalid
 * parameter to a random value within the valid range. This prevents the algorithms
 * (especially the simplex algorithm) from crawling into corners of the valid solution
 * space and getting stuck there because there is gradient that points into an 
 * invalid area.
 * 
 * @param[in] global_options The global options structure
 * @param[in,out] parm_data The parameter structure
 * @return The number of parameters that were changed
 */
int check_range(global_options_struct *global_options, parm_struct *parm_data)
{
  // Do not bounds check the simplex, or it will mess up. :(
  if (global_options->ALGORITHM==SIMPLEX) return 0;
  
  int count=0;
  int i,j;
  for (i=0; i<parm_data->unique_bonds_found; ++i) {
    // Bond force constant between 100 and 1000
    if (parm_data->bond_data[i].rk < 100.0 || parm_data->bond_data[i].rk > 1000.0) {
      parm_data->bond_data[i].rk = 900.0*((double)rand()/(double)RAND_MAX) + 100.0;
      ++count;
    }
    
    // Bond length between 0 and 3 Angstroms
    if (parm_data->bond_data[i].req <= 0.0 || parm_data->bond_data[i].req > 3.0) {
      parm_data->bond_data[i].req = 3.0*((double)rand()/((double)RAND_MAX+1));
      ++count;
    }
  }
  for (i=0; i<parm_data->unique_angles_found; ++i) {
    // Angle force constant between 0 and 200
    if (parm_data->angle_data[i].tk < 0.0 || parm_data->angle_data[i].tk > 200.0){
      parm_data->angle_data[i].tk = 200.0*((double)rand()/(double)RAND_MAX);
      ++count;
    }
    
    // Angle should stay between 0 and PI, wrap it around if it's wrong
    if (parm_data->angle_data[i].teq > PI) {
      parm_data->angle_data[i].teq = fmod(parm_data->angle_data[i].teq, PI);
      ++count;
    } else if (parm_data->angle_data[i].teq < 0.0) {
      parm_data->angle_data[i].teq = PI - fmod(parm_data->angle_data[i].teq, PI);
      ++count;
    }
  }
  for (i=0; i<parm_data->unique_dihedrals_found; ++i) {
    for (j=0; j<parm_data->dihedral_data[i].num_terms; ++j) {
      // Dihedral force constant between -30 and 30 (Glycam has some negative Kp for example)
      if (parm_data->dihedral_data[i].term[j].pk < -30.0 || parm_data->dihedral_data[i].term[j].pk > 30.0) {
        parm_data->dihedral_data[i].term[j].pk = 60*((double)rand()/(double)RAND_MAX) + 30.0; 
        ++count;
      }
      
      // Dihedral periodicity between 0 and 5
      if (parm_data->dihedral_data[i].term[j].pn < 0.0 || parm_data->dihedral_data[i].term[j].pn > 5.0) {
        parm_data->dihedral_data[i].term[j].pn = 10.0*((double)rand()/(double)RAND_MAX) + 5.0;
        ++count;
      }
      
      // Dihedral phase between 0 and 2*PI, wrap around if wrong
      if (parm_data->dihedral_data[i].term[j].phase >= 2.0*PI) {
        parm_data->dihedral_data[i].term[j].phase = fmod(parm_data->dihedral_data[i].term[j].phase, 2.0*PI);
        ++count;
      } else if (parm_data->dihedral_data[i].term[j].phase < 0.0) {
        parm_data->dihedral_data[i].term[j].phase = 2.0*PI - fmod(parm_data->dihedral_data[i].term[j].phase, 2.0*PI);
        ++count;
      }
    }
  }
  return count;
}

