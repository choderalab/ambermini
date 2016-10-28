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

/** @file read_mdcrd.c
 * Contains functions relating to reading and writing coordinate
 * structures.
 */
#include "function_def.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * Reads in each coordinate file from the list of coordinate file.
 * @param[in] global_options The global options structure
 * @param[in] parm_datas The array of parameter sets to use
 * @param[out] struc_datas Pointer to array of structure sets to be initialized. The array is 
 * size 1 already if -c command line option for single mdcrd has been specified, otherwise *s_datas=NULL
 * @return Integer indicating success or failure
 */
int read_mdcrds(global_options_struct *global_options, parm_struct *parm_datas, coord_set **s_datas)
{ 
  // Check that a coordinate file of some sort is defined 
  if (!global_options->mdcrd_list && !*s_datas) {
    printf("*** ERROR - must specify coordinate file(s) with -c or -cf! ***\n");
    command_line_help("paramfit");
    return UNKNOWN_OPT;
  }
  
  int retval;
  FILE *fptr;
  char *filename;
  // Open the file list or allocate and read in single mdcrd
  if (global_options->prmtop_list) {
    if (global_options->mdcrd_list) {
      printf(" Reading mdcrd file list: %s\n",global_options->mdcrd_list);
    } else {
      printf("ERROR: No mdcrd file list specified!\n");
      return INVALID_DATA;
    }
  } else {
    printf(" Reading mdcrd file: %s\n",(*s_datas)[0].filename);
  }
  
  if (global_options->mdcrd_list != NULL) {
    if((fptr=fopen(global_options->mdcrd_list,"r"))==NULL) {
      file_open_failure("read_mdcrds", global_options->mdcrd_list);
      return FILE_OPEN_FAIL;
    }
  } else {
    (*s_datas)[0].natoms = parm_datas[0].NTOTAT;
    retval = alloc_coords(global_options, &(*s_datas)[0]);
    process_retval(retval, global_options->VERBOSITY);
    retval = read_single_mdcrd(global_options, &(*s_datas)[0]);
    process_retval(retval,global_options->VERBOSITY);
    return SUCCESS;
  }
  
  // Allocate memory for one line
  char* line = (char*)malloc(BUFFER_SIZE);
  if (line==NULL) {
    malloc_failure_char("read_mdcrds","line",BUFFER_SIZE);
    return ALLOC_FAIL;
  }
  
  // Read in list of files, one line at a time
  int files_read=0;
  int ns;
  while (s_getline( line, BUFFER_SIZE, fptr) != FILE_READ_FAIL) {
    ++files_read;
    // Allocate memory for new structure set or inital set (*s_datas initally=null so realloc acts like malloc)
    *s_datas=(coord_set*)realloc(*s_datas,files_read*sizeof(coord_set));
    if (*s_datas==NULL) {
      printf("ERROR allocating %d bytes in read_mdcrds for s_datas\n", files_read*sizeof(coord_set));
      return ALLOC_FAIL;
    }
   
    // Allocate temporary filenames string to same length as line so we can be sure it will fit the filename
    char *fn = (char*)malloc(strlen(line)+1);
    if (!fn) {
      malloc_failure_char("read_mdcrds","fn",strlen(line)+1);
      return ALLOC_FAIL;
    }
    char *en = (char*)malloc(strlen(line)+1);
    if (!en) {
      malloc_failure_char("read_mdcrds","en",strlen(line)+1);
      return ALLOC_FAIL;
    }
    // Set the mdcrd filename and number of structures
    if (sscanf(line, "%s %d %s", fn, &ns, en)!=3) {
      printf("ERROR: malformed line in mdcrd list %s line %d\nLine was: %s", global_options->mdcrd_list,
             files_read, line);
      free(en); free(fn);
      return INVALID_DATA;
    }
    // Allocate exactly strlen(fn)+1 characters since add null terminator to it
    (*s_datas)[files_read-1].filename = (char*)malloc(strlen(fn)+1);
    if (!(*s_datas)[files_read-1].filename) {
      malloc_failure_char("read_mdcrd","coord_set->filename",strlen(fn)+1);
      free(en); free(fn);
      return ALLOC_FAIL;
    }
    (*s_datas)[files_read-1].energy_filename=(char*)malloc(strlen(en)+1);
    if (!(*s_datas)[files_read-1].filename) {
      malloc_failure_char("read_mdcrd","coord_set->energy_filename",strlen(en)+1);
      free(en); free(fn);
      return ALLOC_FAIL;
    }
    // Copy everything, and null terminate the strings
    strncpy((*s_datas)[files_read-1].filename, fn, strlen(fn));
    strncpy((*s_datas)[files_read-1].energy_filename, en, strlen(en));
    (*s_datas)[files_read-1].filename[strlen(fn)]='\0';
    (*s_datas)[files_read-1].energy_filename[strlen(en)]='\0';
    // Fillin necessary data
    (*s_datas)[files_read-1].num_coords = ns;
    global_options->TOTAL_STRUCTURES+=ns;
    (*s_datas)[files_read-1].natoms = parm_datas[files_read-1].NTOTAT;
    (*s_datas)[files_read-1].mem_allocated = sizeof(coord_set) + strlen(fn) + strlen(en) + 2;
    free(fn); free(en);
    
    // Allocate memory in the mdcrd  and read in the data
    retval = alloc_coords(global_options, &((*s_datas)[files_read-1]));
    process_retval(retval,global_options->VERBOSITY);
    retval = read_single_mdcrd(global_options, &((*s_datas)[files_read-1]));
    process_retval(retval,global_options->VERBOSITY);
    
  }
  free(line);
  fclose(fptr);
  // Sanity check- #prmtops = #mdcrds
  if (files_read != global_options->num_prmtops) {
    printf("ERROR! Read in %d prmtops but %d mdcrds!\n", global_options->num_prmtops, files_read);
    return INVALID_DATA;
  }
  return SUCCESS;
}
 
/**
 * Reads the coordinate file with all the input structures.
 * @param[in] global_options Contains information on the expected format of the coordinate set
 * @param[in,out] coords_data Pointer to coordinate set structure that will be populated. Filename and number
 * of structures needs to be already initialized. The coordinate data set should also have been allocated.
 * @see alloc_coords
 * @see read_mdcrds
 * @return Integer indicating success or failure.
 */
int read_single_mdcrd(global_options_struct *global_options, coord_set *coords_data)
{
  int i;
  FILE *fptr;
  char mdcrd_temp[BUFFER_SIZE];
  int retval;
  int structure;

  // Check we can open the file
  printf("  Reading mdcrd file    : %s\n", coords_data->filename);
  if((fptr=fopen(coords_data->filename,"r"))==NULL) {
     file_open_failure("read_single_mdcrd", coords_data->filename);
     return FILE_OPEN_FAIL;
  }

  retval=s_getline( mdcrd_temp, BUFFER_SIZE, fptr );  /*blank line*/
  if ( retval==FILE_READ_FAIL ) {
      printf("*** ERROR IN read_mdcrd SUBROUTINE\n");
      printf("*** HIT EOF WHILE READING INITIAL BLANK LINE OF MDCRD FILE: %s\n",coords_data->filename);
      return FILE_READ_FAIL;
  } 
 
  // Check the coordinate file format is as expected
  int int_temp;
  if (global_options->COORD_FORMAT==TRAJECTORY) {
    // Save current file pointer position so we can peek ahead at next line
    fpos_t pos;
    fgetpos(fptr, &pos);
   
    // Peek at next line-- should be first 3 coordinates, not #atoms in system
    s_getline(mdcrd_temp, BUFFER_SIZE, fptr);
    retval=sscanf(mdcrd_temp, "%d", &int_temp);
    if (retval == 1 && int_temp==coords_data->natoms) {
      printf("\nERROR! Coordinate file seems to be RESTART format but you specified TRAJECTORY format!\n");
      printf("       Please check the file and your setting and re-run.\n");
      printf("       Filename \"%s\"\n",coords_data->filename);
      return INVALID_FORMAT;
    }
    // Return file pointer to correct place, as file format check passed
    else {
      fsetpos(fptr, &pos);
      if (global_options->VERBOSITY>=MEDIUM)
        printf("    Coordinate file passed format check\n");
    }
  } else if (global_options->COORD_FORMAT==RESTART) {
    // discard line with number of atoms
    s_getline(mdcrd_temp, BUFFER_SIZE, fptr);
    retval=sscanf(mdcrd_temp, "%d", &int_temp);
    if (retval != 1 || int_temp != coords_data->natoms) {
      printf("\nERROR! Restart file seems to be formatted incorrectly!\n");
      printf("       Check your setting of COORDINATE_FORMAT and the file.\n");
      printf("       Filename \"%s\". Line was:\n", coords_data->filename);
      printf("%s\n", mdcrd_temp);
      return INVALID_FORMAT;
    } else {
      if (global_options->VERBOSITY>=MEDIUM)
        printf("    Coordinate file passed format check\n");
    }
  }
  
  // Read in all of the structures at once
  for (structure=0; structure<coords_data->num_coords; ++structure) {
    for (i=0;i<coords_data->natoms; ++i) {
      /*We should keep reading 3 floats at a time in the order x,y,z
        note, currently there is no allowance for box info in the mdcrd file */
      
      retval=fscanf(fptr, "%lf %lf %lf", &coords_data->struc[structure].x_coord[i], &coords_data->struc[structure].y_coord[i],
                    &coords_data->struc[structure].z_coord[i]);
      if ( retval==EOF || retval != 3 ) {
        printf("ERROR! Cannot read input mdcrd %f!\n", mdcrd_temp);
        printf(" Could not read structure %d atom %d\n", structure, i);
        return FILE_READ_FAIL;
      }
      /*
      retval=fscanf(fptr,"%lf ",&coords_data->struc[structure].x_coord[i]);
      if (retval==EOF) {
        printf("*** ERROR IN READ_MDCRD - HIT EOF WHILE READING ELEMENT: %d\n",i);
        return FILE_READ_FAIL;
      }
      else if (retval!=1) {
        printf("*** ERROR IN READ_MDCRD - FAILED TO READ X_COORD FOR ELEMENT: %d\n",i);
        return INVALID_FORMAT;
      }
    printf("atom 1 = %f\n", coords_data->struc[structure].x_coord[i]);
      retval=fscanf(fptr,"%lf",&coords_data->struc[structure].y_coord[i]);
      if (retval==EOF) {
        printf("*** ERROR IN READ_MDCRD - HIT EOF WHILE READING ELEMENT: %d\n",i);
        return FILE_READ_FAIL;
      }
      else if (retval!=1) {
        printf("*** ERROR IN READ_MDCRD - FAILED TO READ Y_COORD FOR ELEMENT: %d\n",i);
        return INVALID_FORMAT;
      }
      retval=fscanf(fptr,"%lf",&coords_data->struc[structure].z_coord[i]);
      if (retval==EOF) {
        printf("*** ERROR IN READ_MDCRD - HIT EOF WHILE READING ELEMENT: %d\n",i);
        return FILE_READ_FAIL;
      }
      else if (retval!=1) {
        printf("*** ERROR IN READ_MDCRD - FAILED TO READ Z_COORD FOR ELEMENT: %d\n",i);
        return INVALID_FORMAT;
      }
      */
    }
  }
 fclose(fptr);
 return SUCCESS;
}
