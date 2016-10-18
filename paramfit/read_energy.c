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
#include <stdlib.h>
#include <sys/stat.h>
#include <string.h>
#include "constants.h"
#include "function_def.h"

/** @file read_energy.c
 * Contains routines for reading in energy data into the coordinate structures.
 * This includes support for reading in each format of QM file, AMBER files (TODO) and
 * straight up lists of energies. This is all handled by the unified function read_qm
 * which calls the appropriate function to the qm data type.
 */

/**
 * The master function for all reading of QM input data. Will call functions to read in 
 * either forces or energies from a list or from raw output files from the appropriate program.
 * 
 * @param[in] global_options The global options structure, which specifies the data format and energy or forces
 * @param[in] parm_datas Array of parameter sets corresponding to coordinates to read in
 * @param[in,out] coords_datas Array of coordinate sets that each have the filename/folder to read from initialized.
 * These will be updated with the QM energies or forces
 * @return Integer indicating success or failure
 */
int read_qm(global_options_struct *global_options, parm_struct *parm_datas, coord_set *coords_datas)
{
  int i, retval;
  int nread=0;
  // Process each coordinate set individually
  for (i=0; i<global_options->num_prmtops; ++i) {
    if (global_options->FUNC_TO_FIT==AMBER_FORCES) {
      if (global_options->QMFILEFORMAT==GAUSSIAN) {
        retval=read_gaussian_forces(global_options, &parm_datas[i], &coords_datas[i]);
        if (retval != SUCCESS) printf("ERROR reading Gaussian forces for coord set %s\n", coords_datas[i].filename);
      }
      else {
        printf("ERROR: only Gaussian is supported for Amber force fitting currently\n");
        return INVALID_DATA;
    } } 
    retval=read_qm_energy_list(global_options, &parm_datas[i], &coords_datas[i]);
    if (retval==SUCCESS) ++nread;
    else {
      printf("ERROR reading QM data for coord set %s\n", coords_datas[i].filename);
      return INVALID_DATA;
    }
  }
  return SUCCESS;
}

/**
 * Reads in the quantum energies and attaches one to each coordinate structure.
 * Performs unit conversion if necessary. This is used when evaluating energies, NOT forces.
 * The QM data file should consist of NSTRUCTURES worth of doubles, with one value per line.
 * 
 * @param[in] global_options The global options structure, containing filename with energies
 * @param[in] parm_data Pointer to single parameter structure corresponding to this coordinate set
 * @param[in,out] coords_data Pointer to single coordinate set that energies will be attached to.
 * @return Integer indicating success or failure.
 */
int read_qm_energy_list(global_options_struct *global_options, parm_struct *parm_data, coord_set *coords_data)
{
  int i;
  FILE *fptr;
  int retval;
  double conversion_factor;

  if (global_options->QM_ENERGY_UNITS==HARTREE)
    conversion_factor = HARTREE_TO_KCALMOL;
  else if (global_options->QM_ENERGY_UNITS==KJMOL)
    conversion_factor = KJMOL_TO_KCALMOL;
  else
    conversion_factor = 1.0;
  
  if (global_options->VERBOSITY >= MEDIUM)
    printf("  Reading energy file or directory  : %s\n",coords_data->energy_filename);
  /*Initially perform a sanity check to see if we are reading at least 1 structure*/
  if ( coords_data->num_coords < 1 ) {
     printf("*** ERROR IN read_qm_energy SUBROUTINE\n");
     printf("*** INVALID NUMBER OF STRUCTURES FOR MDCRD %s.\n",coords_data->filename);
     printf("*** NUMBER OF STRUCTURES SPECIFIED WAS %d WHICH IS NOT > 0\n",coords_data->num_coords);
     return INVALID_DATA;
  }

  /* Check if this is a directory and if so process appropriately */
  struct stat s;
  if (stat(coords_data->energy_filename,&s)!=-1 && S_ISDIR(s.st_mode) ) {
    retval = read_qm_directory(global_options, parm_data, coords_data);
    process_retval(retval, global_options->VERBOSITY);
    return SUCCESS;
  }
  /*Check we can open the file*/
  if((fptr=fopen(coords_data->energy_filename,"r"))==NULL) {
     file_open_failure("read_qm_energy", coords_data->energy_filename);
     return FILE_OPEN_FAIL;
  }
  
  // Allocate memory for one line
  char* line = (char*)malloc(BUFFER_SIZE);
  if (line==NULL) {
    malloc_failure_char("read_prmtops","line",BUFFER_SIZE);
    return ALLOC_FAIL;
  }
  
  /*Now we read in NSTRUCTURES worth of QM energies and check we don't hit EOF while running*/
  for (i=0;i<coords_data->num_coords;++i)
  {
    // Get a line from the dataset
    retval=s_getline( line, BUFFER_SIZE, fptr );  /*blank line*/
    if ( retval==FILE_READ_FAIL ) {
      printf("ERROR: Cannot read line %d in energy filename '%s'\n",i,coords_data->energy_filename);
      free(line);
      return FILE_READ_FAIL;
    }
    
    // Check if there are weights
    if (sscanf(line, "%lf %lf",&coords_data->struc[i].energy, &coords_data->struc[i].weight)!=2) {
      if (sscanf(line, "%lf", &coords_data->struc[i].energy )==1)
        coords_data->struc[i].weight = 1.0;
      else {
        printf("ERROR in read_qm_energy - failed to read energy for mdcrd '%s' structure #%d\n",
               coords_data->energy_filename, i);
        free(line);
        return FILE_READ_FAIL;
      }
     } else if (coords_data->struc[i].weight != 1.0) {
       printf("Structure %i weight = %f\n", i, coords_data->struc[i].weight);
     }
     coords_data->struc[i].energy *= conversion_factor;
     coords_data->struc[i].init_energy=0.0;
  }
  fclose(fptr);
  free(line);
  
  /* Sort the structures in order of increasing energy */
  if (global_options->SORT_ENERGY==YES) {
    printf("  Sorting structures in order of increasing energy\n");
    qsort( coords_data->struc, coords_data->num_coords, sizeof(coords_struct), compare_energy);
  }
  return SUCCESS;
}

/** 
 * Reads all qm output files in a given directory. The names of the files are according to the
 * convention paramfit writes them in with CREATE_INPUT mode- that is prmtop.QMFILEOUTEND.## where prmtop
 * is prmtop->filename. This will support a variety of file formats, which one it is is set in 
 * global_options->QMFORMAT or similar.
 * 
 * @param[in] global_options The global options structure
 * @param[in,out] coords_data Pointer to single coordinate set. Each structure in the set will be populated with 
 * the energy / forces from the matching QM output file.
 * @return Integer indicating success or failure
 */
int read_qm_directory(global_options_struct *global_options, parm_struct *parm_data, coord_set *coords_data)
{
  int i;
  int retval;
  
  // Sanity check that filename things are defined 
  if (!global_options->QMFILEOUTSTART) {
    printf("ERROR: QMFILEOUTSTART not defined in job control file.\n");
    return INVALID_DATA;
  } else if (!global_options->QMFILEOUTEND) {
    printf("ERROR: QMFILEOUTEND not defined in job control file.\n");
    return INVALID_DATA;
  }
  
  // Read in coordinates for this directory depending on type
  if (global_options->FUNC_TO_FIT==AMBER_FORCES) {
    if (global_options->QMFILEFORMAT==GAUSSIAN) {
      retval = read_gaussian_forces(global_options, parm_data, coords_data);
      process_retval(retval, global_options->VERBOSITY);
    } else {
      printf("ERROR: Forces currently only supported with Gaussian output files.\n");
      return INVALID_FORMAT;
    }
  } else {
    if (global_options->QMFILEFORMAT==GAUSSIAN) {
      retval = read_gaussian_energy(global_options, parm_data, coords_data);
      process_retval(retval, global_options->VERBOSITY);
    } else {
      printf("ERROR: Reading QM output files only supported for Gaussian at this time.\n");
      return INVALID_FORMAT;
    }
  }
  return SUCCESS;
}

/**
 * Reads in Gaussian output files with forces for all structures.
 * The files need to have been generated in the CREATE_INPUT mode so that the atom numbering is
 * consistent. This also lets us assume that the files are named QMFILEOUTSTART###QMFILEOUTEND.out
 *
 * @param[in] global_options The global options structure
 * @param[in] parm_data Pointer to single parameter file
 * @param[in,out] coords_data Pointer to single coordinate set, will have forces attached to it.
 * @return Integer indicating success or failure.
 */
int read_gaussian_forces(global_options_struct *global_options, parm_struct *parm_data, coord_set *coords_data)
{ 
  // Sanity check- single prmtop fit only for now
  if (global_options->num_prmtops != 1) {
    printf("*** ERROR: Only a single prmtop can be used with forces at this time.\n");
    return INVALID_DATA;
  }
  
  // Determine if a unit conversion is necessary- internal force units are kcal/mol-A
  double units; char *filename;
  if (global_options->QM_FORCE_UNITS==HARTREE_BOHR)
    units = HARTREE_TO_KCALMOL*BOHR_TO_ANGSTROM;
  else if (global_options->QM_FORCE_UNITS==KCALMOL_ANGSTROM)
    units = 1.0;
  else {
    printf("*** ERROR: Unsupported force units for Gaussian input files.\n");
    return INVALID_DATA;
  }
    
  int i;
  for(i=0; i<coords_data->num_coords; ++i) {
    // Create the filename and open the file
    filename = (char*)malloc(strlen(global_options->QMFILEOUTSTART)+strlen(global_options->QMFILEOUTEND)+5);
    if (!filename) {
      malloc_failure_char("read_gaussian_forces","filename", strlen(coords_data->energy_filename)+strlen(global_options->QMFILEOUTEND)+strlen(global_options->QMFILEOUTSTART)+10);
      return ALLOC_FAIL;
    }
    sprintf(filename,"%s%d%s",global_options->QMFILEOUTSTART,i,global_options->QMFILEOUTEND);
    
    // Open the file
    FILE *fptr = fopen(filename, "r");
    if(!fptr) {
      file_open_failure("read_gaussian_forces", filename);
      return FILE_OPEN_FAIL;
    }
    
    /* Read and discard lines until axes line is found, discard 4 lines, begin reading.
     * Output will look like this:
     *  ***** Axes restored to original set *****
     * -------------------------------------------------------------------
     * Center     Atomic                   Forces (Hartrees/Bohr)
     * Number     Number              X              Y              Z
     * -------------------------------------------------------------------
     *      1        1           0.011759647    0.003870626    0.010325240
     * ...
     */
    char line[200];
    int check;
    do {
      check = s_getline(line, 200, fptr);
      if ( check==FILE_READ_FAIL ) {
        printf("*** ERROR IN read_gaussian_forces SUBROUTINE\n");
        printf("*** Failed to find forces section in output file: %s\n\n", filename);
        printf("*** Check that this is the correct filename and that\n");
        printf("    the Gaussian job completed successfully!\n");
        return FILE_READ_FAIL;
      }
    } while( !strstr(line, "Axes restored to original set") ); 

    // Discard 3 lines of table header and start on 4th
    int atom;
    for (atom=0; atom<4; ++atom) {
      s_getline(line, 200, fptr);
    }
    
    // Read in the data for each atom
    for (atom=0; atom<coords_data->natoms; ++atom) {
      int element = find_atomic_number_from_parm(parm_data, atom+1);
      if ( fscanf(fptr, "%*d %d %15lf %15lf %15lf",  &check, &coords_data->struc[i].force[atom].x,
                  &coords_data->struc[i].force[atom].y, &coords_data->struc[i].force[atom].z) != 4 ) {
        printf("*** ERROR IN read_gaussian_forces SUBROUTINE\n");
        printf("*** Error in reading forces for atom %d\n", atom);
        printf("*** Filename: %s\n", filename);
        return FILE_READ_FAIL;
      }
      if ( check != element ) {
        printf("*** ERROR IN read_gaussian_forces SUBROUTINE\n");
        printf("*** Unexpected element in file: %s\n", filename);
        printf("*** Expected: %d \t Received: %d \n\n", element, check);
        printf("*** Check that the input files used for these output files were\n");
        printf("    generated with paramfit for this same prmtop.\n");
        return FILE_READ_FAIL;
      }
      // Unit conversion
      coords_data->struc[i].force[atom].x *= units;
      coords_data->struc[i].force[atom].y *= units;
      coords_data->struc[i].force[atom].z *= units;
    }
    free(filename);
    fclose(fptr);
  }
  return SUCCESS;
}

/** 
 * Reads in the gaussian energy from a certain set of files into a coordinate structure
 * This lets paramfit write the file, run gaussian, read the file without any need for 
 * the user to go through the gaussian output.
 * @param[in] global_options The global options structure
 * @param[in] parm_data Pointer to single parm struct 
 * @param[in,out] coords_data Pointer to single coord struct to attach energy to
 * @return Integer indicating success or failure
 */
int read_gaussian_energy(global_options_struct *global_options, parm_struct *parm_data, coord_set *coords_data)
{
  int i, retval;
  short normal=NO;
  char *filename, *line;
  FILE *fptr;
  
  for (i=0; i<coords_data->num_coords; ++i) {
    // Create the filename
    filename = (char*)malloc(strlen(coords_data->energy_filename)+strlen(global_options->QMFILEOUTSTART)+strlen(global_options->QMFILEOUTEND)+10);
    if (!filename) {
      malloc_failure_char("read_gaussian_energy","filename",strlen(coords_data->energy_filename)+strlen(global_options->QMFILEOUTSTART)+
      strlen(global_options->QMFILEOUTEND)+10);
      return ALLOC_FAIL;
    }
    sprintf(filename, "%s/%s%d%s", coords_data->energy_filename, global_options->QMFILEOUTSTART, i, global_options->QMFILEOUTEND);
    // Open the file
    if ( (fptr=fopen(filename, "r"))==NULL) {
      file_open_failure("read_gaussian_energy", filename);
      return FILE_OPEN_FAIL;
    }
    line = (char*)malloc(BUFFER_SIZE);
    if (!line) {
      malloc_failure_char("read_gaussian_energy", "line", BUFFER_SIZE);
      return ALLOC_FAIL;
    }
    // Go through line by line
    while (s_getline( line, BUFFER_SIZE, fptr) != FILE_READ_FAIL) {
      if (strstr(line, "SCF Done:")) 
        sscanf(line, "%*s%*s%*s%*s%lf", &(coords_data->struc[i].energy));
      else if (strstr(line, "Normal termination of Gaussian"))
        normal=YES;
    }
    if (normal==NO) {
      printf("ERROR: File '%s' is not a normal Gaussian output file\n",filename);
      return INVALID_DATA;
    }
    
    free(line);
    free(filename);
    fclose(fptr);
  }
  return SUCCESS;
}
