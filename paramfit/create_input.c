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

/** @file create_input.c
 * This file contains the routines for making the job files from
 * a prmtop and mdcrd file that can then be used to obtain the energy
 * data to be fitted against
 */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <stdlib.h>
// #include unistd.h
#include "function_def.h"

/**
 * Main input creating function for multiple prmtops
 * Creates input directories if non existent, creates the files to write in the
 * dircctory, and calls the appropriate function to write the correct format. If there
 * is only one molecule, falls through to old fashioned create_input for single prmtops.
 * 
 * @param[in] global_options The global options structure
 * @param[in] parm_datas Pointer to the array of parm structures
 * @param[in] coords_datas Pointer to array of coordinate structures
 * @return Integer indicating success or failure
 */
int create_qm_input(global_options_struct *global_options, parm_struct *parm_datas, coord_set *coords_datas)
{
  int retval;
  char *filename;
  FILE *fptr;
  int i,j;
  
  // Sanity check that relevant job control variables are defined
  if (!global_options->QMFILEOUTSTART) {
    printf("ERROR: QMFILEOUTSTART not defined in job control file.\n");
    return INVALID_DATA;
  } else if (!global_options->QMFILEOUTEND) {
    printf("ERROR: QMFILEOUTEND not defined in job control file.\n");
    return INVALID_DATA;
  } else if (!global_options->QMHEADER) {
    printf("ERROR: QM_HEADER not defined in job control file.\n");
    return INVALID_DATA;
  }
  
  // Print out a message if fitting forces
  if (global_options->FUNC_TO_FIT==AMBER_FORCES) {
    if (global_options->QMFILEFORMAT==GAUSSIAN)
      printf("!  Will fit forces- ensure the route line in your header has the force keyword!\n");
    else {
      printf("*** ERROR: Force calculation input file creation is currently only supported for Gaussian.\n");
      return INVALID_DATA;
  } }
  
  // If single prmtop, just do old fashioned create input
  if (global_options->num_prmtops==1)
    return create_input_single_prmtop(global_options, parm_datas, coords_datas);
   
  // Loop through each prmtop
  for (i=0; i<global_options->num_prmtops; ++i) {
    // Check if directory already exists for these files. If not, create it
    struct stat s;
    if (!stat(coords_datas[i].energy_filename, &s)) { // returns 0 on success
#if WIN32
      if (!S_ISDIR(s.st_mode)) {
#else
      if (!S_ISDIR(s.st_mode) && !S_ISLNK(s.st_mode)) {
#endif
        printf("ERROR! Energy directory for coordinate set '%s' is '%s'\n", coords_datas[i].filename, coords_datas[i].energy_filename);
        printf("       Not a directory or symlink\n");
        return INVALID_DATA;
    } }
    else { // create the directory
      printf("   Creating directory %s\n", coords_datas[i].energy_filename);
#if WIN32
      if (mkdir(coords_datas[i].energy_filename)) { // returns 0 on success
#else
      if (mkdir(coords_datas[i].energy_filename, 0700)) { // returns 0 on success
#endif
        printf("ERROR! Could not create directory '%s' for coordinate set '%s'\n", coords_datas[i].energy_filename, coords_datas[i].filename);
        return FAILURE;
    } }
    
    // Loop over each structure for this prmtop
    for (j=0; j<coords_datas[i].num_coords; ++j) {
      // Create the filename and open the file
      filename = (char*)malloc(strlen(coords_datas[i].energy_filename)+strlen(global_options->QMFILEOUTEND)+strlen(global_options->QMFILEOUTSTART)+10);
      if (!filename) {
        malloc_failure_char("create_qm_input","filename", strlen(coords_datas[i].energy_filename)+strlen(global_options->QMFILEOUTEND)+strlen(global_options->QMFILEOUTSTART)+10);
        return ALLOC_FAIL;
      }
      sprintf(filename,"%s/%s%d%s",coords_datas[i].energy_filename,global_options->QMFILEOUTSTART,j,global_options->QMFILEOUTEND);
      if((fptr=fopen(filename,"w"))==NULL) {
        file_open_failure("create_input", filename);
        return FILE_OPEN_FAIL;
      }
      
    // Write the file with the correct format
     if (global_options->QMFILEFORMAT==GAUSSIAN)
        retval=write_input_gaussian(global_options, &parm_datas[i], &coords_datas[i], j, fptr);
     else if (global_options->QMFILEFORMAT==ADF)
       retval = write_input_adf(global_options, &parm_datas[i], &coords_datas[i], j, fptr);
     else if (global_options->QMFILEFORMAT==GAMESS)
       retval = write_input_gamess(global_options, &parm_datas[i], &coords_datas[i], j, fptr);
     else {
       printf("   ERROR IN create_input SUBROUTINE\n");
       printf("   UNKNOWN FILE FORMAT TO WRITE: %d.\n",global_options->QMFILEFORMAT);
      return UNKNOWN_OPT;
     }
     if (retval!=SUCCESS) return retval;
     fclose(fptr);
     free(filename);
     
    }
    printf("    Successfully wrote %i input files for coordinate set '%s' in directory '%s'\n", coords_datas[i].num_coords, coords_datas[i].filename,
           coords_datas[i].energy_filename);
  }
  return SUCCESS;
}

/**
 * Main input creating function for a single prmtop
 * Sets up input files to write and calls the appropriate function to write it in the
 * correct format.
 * @param[in] global_options The global options structure
 * @param[in] parm_data Pointer to a single set of parameters for this molecule
 * @param[in] coords_data Pointer to a single coordinate set for this molecule
 * @return Integer indicating success or failure
 */
int create_input_single_prmtop(global_options_struct *global_options, parm_struct *parm_data, coord_set *coords_data)
{
  FILE *fptr;
  int current_struct=0;
  int retval;
  char filename_to_write[1024];

  if (global_options->VERBOSITY>=HIGH) {
     printf(" Creating input files. Will write a total of %d sequentially numbered files.\n",coords_data->num_coords);
     printf("           %6d > .",0);
  }

  /* Read the mdcrd if necessary */
  if (coords_data == NULL) {
    retval=read_single_mdcrd(global_options,coords_data);
    if (retval!=SUCCESS) {
      printf("*** ERROR IN create_input SUBROUTINE\n");
      printf("*** FAILED TO READ STRUCTURE %d FROM MDCRD FILE.\n",current_struct);
      return FILE_READ_FAIL;
    }
  }
  
  /* Print out a warning if fitting forces */
  if (global_options->FUNC_TO_FIT==AMBER_FORCES) {
    if (global_options->QMFILEFORMAT==GAUSSIAN)
    printf("!  Will fit forces- ensure the route line in your header has the force keyword!\n");
    else {
      printf("*** ERROR: Force calculation input file creation is currently only supported for Gaussian.\n");
      return INVALID_DATA;
    }
  }
  
  /*loop over all the structures*/
  for (current_struct=0; current_struct<coords_data->num_coords;++current_struct)
  {
     /*Open the output file with the relevant name*/
     /*filename to write is made up of QMFILEOUTSTART, current_struct and QMFILEOUTEND
     the char variable filename_to_write has been allocated as 1024 bytes long but we should
     still really have a sanity check here to ensure we don't exceed it - not that many people
     would ever want a filename over 1023 characters long*/
     if (strlen(global_options->QMFILEOUTSTART)+strlen(global_options->QMFILEOUTEND)+5 >= 1024) {
       printf("ERROR! Probable overflow in quantum input filenames.\n");
       printf("       Reduce the number of characters in QMFILEOUTSTART and QMFILEOUTEND\n");
       return INVALID_DATA;
     }

     /*Note, a potential overflow could occur here - currently NOT checked for*/
     
     /*Now fill the filename_to_write array*/
     sprintf(filename_to_write,"%s%d%s",global_options->QMFILEOUTSTART,current_struct,global_options->QMFILEOUTEND);
     if (global_options->VERBOSITY>=HIGH)
           printf("   Filename to be written is: %s\n",filename_to_write);
     
     if((fptr=fopen(filename_to_write,"w"))==NULL) {
        file_open_failure("create_input", filename_to_write);
        return FILE_OPEN_FAIL;
     }

     /*print info for each file we write*/
     if (global_options->VERBOSITY>=HIGH) {
       /*check if current_struct is divisible by 50, if it is print a new line*/
       if(current_struct%50==0)
         printf("\n           %6d > .",current_struct);
       else
         printf(".");
       fflush(stdout); /*Flush the printf buffer*/
     }

     /*Now we have the file open call the relevant routines for writing the data*/
     if (global_options->QMFILEFORMAT==GAUSSIAN)
        retval=write_input_gaussian(global_options, parm_data, coords_data, current_struct, fptr);
     else if (global_options->QMFILEFORMAT==ADF)
       retval = write_input_adf(global_options, parm_data, coords_data, current_struct, fptr);
     else if (global_options->QMFILEFORMAT==GAMESS)
       retval = write_input_gamess(global_options, parm_data, coords_data, current_struct, fptr);
     else {
       /*Unknown option*/
       printf("   ERROR IN create_input SUBROUTINE\n");
       printf("   UNKNOWN FILE FORMAT TO WRITE: %d.\n",global_options->QMFILEFORMAT);
      return UNKNOWN_OPT;
     }
     if (retval!=SUCCESS)
       return retval;
     fclose(fptr);
  }

  printf("*   Successfully wrote %i files\n",coords_data->num_coords); 
  return SUCCESS;
}

