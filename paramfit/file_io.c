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

/** @file file_io.c
 * Responsible for all file input with the exception of processing the
 * prmtop file.
 * @see read_prmtop.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "function_def.h"
#include "constants.h"

/**
 * Reads the job control file and puts the options into the options structure.
 * @param[in,out] global_options The global options structure, which contains the location of the job control file.
 * @param[in,out] parm_datas Pointer to parameter struct
 * @param[in,out] coords_data Pointer to array of coords struct already allocated to size 1 that may have its number of
 * @return Integer indicating success or failure
 */
int read_job_control_file(global_options_struct *global_options, parm_struct* parm_data, coord_set *coords_data)
{
   int  int_retval;
   int num_lines; // number of lines in job control file
   int longest_line;
   int number_settings;
   int line_length;
   register int temp_count;
   char *line=NULL;
   char temp_char;
   FILE *fptr=NULL;
   
   num_lines=0;
   longest_line=0;
   temp_count=0;
   number_settings=0;
   /*We open the file in *ptr_job_control_filename and then set up our settings*/

   /*We read a line at a time from the control file*/
   if (global_options->VERBOSITY>=MEDIUM)
      printf(" Reading job control file: %s\n",global_options->job_control_filename);

   /*Read the entire settings file into memory and then work from there*/
   // Open the file
   if((fptr=fopen(global_options->job_control_filename,"r"))==NULL)
   {
      file_open_failure("read_job_control_file", global_options->job_control_filename);
      return FILE_OPEN_FAIL;
   }

  // Count the number of lines in the file and the length of the longest one
  while (fscanf(fptr,"%c",&temp_char)!=EOF)
  {
  if (temp_char == '\n')
  {
    ++num_lines;
    if (temp_count > longest_line)
    {
      longest_line = temp_count;
      temp_count = 0;
    }
  }
  else
    ++temp_count;
  }
  rewind(fptr);
  
  // Allocate memory - we will read in one line at a time so array should be
  //                   as big as the longest line plus a null character
    line = (char*)calloc(longest_line+1,sizeof(char));
    if (line==NULL)
    {
      malloc_failure_char("read_job_control_file", "line", longest_line+1);
      return ALLOC_FAIL;
    }
    
  // Go through the file and grab and process one line at a time
  for (temp_count=0; temp_count < num_lines; ++temp_count)
  {
    fgets(line, longest_line+1, fptr);
    
    // Get the length of the line
    int_retval=0;
    line_length=0;
    while (int_retval < longest_line+1)
    {
      if (line[int_retval] != '\n' && line[int_retval] != '\0')
	++line_length;
      
      ++int_retval;
    }
    
    // Process the line if it is a valid setting
    if (line[0] != '#' && line_length > 2)
    {
	int_retval=process_job_control_setting(line, line_length, &number_settings, global_options, parm_data, coords_data);
      if (int_retval!=SUCCESS)
      {
	if (int_retval==INVALID_FORMAT)
	  printf("!  Error near line %d of job control file, Unknown variable name - ignored.\n",temp_count);
	if (int_retval==INVALID_DATA)
	  printf("!  Error near line %d of job control file, Invalid data for selected variable - ignored.\n",temp_count);
	if (int_retval==INVALID_LINE)
	  printf("!  Error near line %d of settings file, no '=' found - ignored.\n",temp_count);
	if (int_retval==ALLOC_FAIL)
	{
	  printf("!  MALLOC FAILURE in process_setting(), attempting to ignore...\n");
	  printf("!  WARNING - CHECK PARAMETERS - RESULTS MAY BE INVALID\n");
	}
      }
//       else
// 	++number_settings;
    }
    
    // Clean out the line buffer
    for (int_retval=0; int_retval < longest_line+1; ++int_retval)
      line[int_retval] = '\0';
  }
  
  fclose(fptr);

  if (global_options->VERBOSITY>=MEDIUM)
    printf(" Job Control: Read a total of %d lines from job_control file. %d options set.\n\n",temp_count+1, number_settings);
  
  free(line);
  line = NULL;
  
  return(SUCCESS);
}

/**
 * Reads in a previously saved list of parameters to fit.
 * Parameters are listed by name, so can be used between and across prmtops.
 * @param[in] global_options The global options structure, containing the file name
 * @param[in,out] parm_data The parameter data file, will be updated with parameters to fit.
 * @return Integer indicating success or failure.
 */
int read_parameter_file_v2(global_options_struct *global_options, parm_struct *parm_datas)
{
  FILE *fptr;
  char *line = calloc(128,sizeof(char));
  char a1[NAME_SIZE],a2[NAME_SIZE],a3[NAME_SIZE],a4[NAME_SIZE];
  short int temp1, temp2, temp3;
  int i,term;
  
   // open the file
  if((fptr=fopen(global_options->PARAMETER_FILE_NAME,"r"))==NULL)
  {
    file_open_failure("read_parameter_file", global_options->PARAMETER_FILE_NAME);
    free(line);
    return FILE_OPEN_FAIL;
  }
  // Go to the section with bonds
  do {
    fgets(line, 127, fptr);
    if (feof(fptr)) {
      printf("ERROR: Premature EOF looking for BOND section in parameter file!\n");
      return INVALID_DATA;
    }
  } while (strncmp(line, "#### BOND", 9) !=0);
  
  // Read in the bonds
  fgets(line,127,fptr);
  while (strncmp(line, "#",1)!=0) {
    sscanf(line,"%s %s %hd %hd", a1, a2, &temp1, &temp2); 
    if (strlen(a1)==1) sprintf(a1,"%s ",a1);
    if (strlen(a2)==1) sprintf(a2,"%s ",a2);
    for (i=0; i<parm_datas[0].unique_bonds_found; ++i) {
      // append spaces if necessary
      if (strncmp(parm_datas[0].bond_data[i].atom_type1,a1,NAME_SIZE)==0 &&
          strncmp(parm_datas[0].bond_data[i].atom_type2,a2,NAME_SIZE)==0 )
        break;
    }
    if (i<parm_datas[0].unique_bonds_found) {
      parm_datas[0].bond_data[i].DO_FIT_REQ=temp1;
      parm_datas[0].bond_data[i].DO_FIT_KR=temp2;
    } else {
      printf("No match for bond %s-%s in prmtop!\n",a1,a2);
      return INVALID_DATA;
    }
    fgets(line,127,fptr);
  }
  
  // Now angles
  do {
    fgets(line, 127, fptr);
    if (feof(fptr)) {
      printf("ERROR: Premature EOF looking for ANGLE section in parameter file!\n");
      return INVALID_DATA;
    }
  } while (strncmp(line, "#### ANGL", 9)!=0);
  
  // Read in the angles
  fgets(line, 127, fptr);
  while (strncmp(line,"#",1)!=0) {
    sscanf(line, "%s %s %s %hd %hd", a1, a2, a3, &temp1, &temp2);
    if (strlen(a1)==1) sprintf(a1,"%s ",a1);
    if (strlen(a2)==1) sprintf(a2,"%s ",a2);
    if (strlen(a3)==1) sprintf(a3,"%s ",a3);
    for (i=0; i<parm_datas[0].unique_angles_found; ++i) {
      if (strncmp(parm_datas[0].angle_data[i].atom_type1,a1,NAME_SIZE)==0 && 
          strncmp(parm_datas[0].angle_data[i].atom_type2,a2,NAME_SIZE)==0 && 
          strncmp(parm_datas[0].angle_data[i].atom_type3,a3,NAME_SIZE)==0 )
        break;
    }
    if (i<parm_datas[0].unique_angles_found) {
      parm_datas[0].angle_data[i].DO_FIT_KT=temp1;
      parm_datas[0].angle_data[i].DO_FIT_THEQ=temp2;
    } else {
      printf("No match for angle %s-%s-%s in prmtop!\n",a1,a2,a3);
      return INVALID_DATA;
    }
    fgets(line,127,fptr);
  }
  
  // Go to dihedrals
  do {
    fgets(line, 127, fptr);
    if (feof(fptr)) {
      printf("ERROR: Premature EOF looking for DIHEDRAL section in parameter file!\n");
      return INVALID_DATA;
    }
  } while (strncmp(line, "#### DIHE", 9)!=0);
  
  // Read in the dihedrals
  fgets(line,127,fptr);
  while (!feof(fptr)) {
    sscanf(line, "%s %s %s %s %d %hd %hd %hd", a1, a2, a3, a4, &term, &temp1, &temp2, &temp3);
    if (strlen(a1)==1) sprintf(a1,"%s ",a1); // if single letter, add a space following it
    if (strlen(a2)==1) sprintf(a2,"%s ",a2);
    if (strlen(a3)==1) sprintf(a3,"%s ",a3);
    if (strlen(a4)==1) sprintf(a4,"%s ",a4);
    for (i=0; i<parm_datas[0].unique_dihedrals_found; ++i) {
      if (strncmp(parm_datas[0].dihedral_data[i].atom_type1,a1,NAME_SIZE)==0 &&
          strncmp(parm_datas[0].dihedral_data[i].atom_type2,a2,NAME_SIZE)==0 &&
          strncmp(parm_datas[0].dihedral_data[i].atom_type3,a3,NAME_SIZE)==0 &&
          strncmp(parm_datas[0].dihedral_data[i].atom_type4,a4,NAME_SIZE)==0)
        break;
    }
    if (i<parm_datas[0].unique_dihedrals_found && term<parm_datas[0].dihedral_data[i].num_terms) {
      parm_datas[0].dihedral_data[i].term[term].DO_FIT_KP=temp1;
      parm_datas[0].dihedral_data[i].term[term].DO_FIT_NP=temp2;
      parm_datas[0].dihedral_data[i].term[term].DO_FIT_PHASE=temp3;
    } else {
      printf("No match for dihedral %s-%s-%s-%s term %d in prmtop!\n",a1,a2,a3,a4,term);
      return INVALID_DATA;
    }
    fgets(line,127,fptr);
  }
  free(line);
  fclose(fptr);
  if (global_options->VERBOSITY>=MEDIUM)
    printf("   Prmtop     (info): Successfully read in saved parameter information\n");
  return SUCCESS;
}
/**
 * Reads in a previously saved list of parameters to fit.
 * The parameters are in the exact order as the prmtop, making this non-transferable between
 * prmtops. This is now included for backwards-compatibility. It will check if you have a nwe
 * format parameter file and call the v2 function to read it if found.
 * @see read_parameter_file_v2
 * 
 * @param[in] global_options The global options structure, containing the file name
 * @param[in,out] parm_data The parameter data file, will be updated with parameters to fit.
 * @return Integer indicating success or failure.
 */
int read_parameter_file(global_options_struct *global_options, parm_struct *parm_data)
{

  FILE *fptr;
  char *line = (char*)malloc(128);
  int line_number = 0;
  short int verify = 0;
  int i = 0; int j=0;
  short int temp1, temp2, temp3;
  short new=NO;
  
  // open the file
  if((fptr=fopen(global_options->PARAMETER_FILE_NAME,"r"))==NULL) {
    file_open_failure("read_parameter_file", global_options->PARAMETER_FILE_NAME);
    free(line);
    return FILE_OPEN_FAIL;
  }
  // Check if this is a new or old format parameter file 
  fgets(line,127,fptr);
  if (strncmp(line, "#V2", 3)==0) {
    free(line);
    fclose(fptr);
    i = read_parameter_file_v2(global_options, parm_data);
    process_retval(i, global_options->VERBOSITY);
    return SUCCESS;
  } else if (global_options->num_prmtops>1) { // require new format file for multiprmtops
    printf("\nERROR! Must have V2 format parameter file for multiprmtop fits.\n");
    printf("       Re-run with RUNTYPE=SET_PARAMS to create this file.\n\n");
    return INVALID_FORMAT;
  }
  
  // locate the section with bond information
  // Skip two lines to get to the section with bonds
  for (i=0; i<3; i++) {
    fgets(line,127,fptr);
    ++line_number;
  }
  
  for (i=0; i<parm_data->unique_bonds_found; ++i) {
    if (fscanf(fptr, "%hd %hd %hd", &verify, &temp1, &temp2) != 3) {
      printf("ERROR reading data from parameter file for bond #%d.\n", i);
      free(line);
      return INVALID_DATA;
    }
    parm_data->bond_data[i].DO_FIT_REQ = temp1;
    parm_data->bond_data[i].DO_FIT_KR = temp2;
  }
  
  // Skip 4 lines to get to the section with angles
  for (i=0; i<4; i++) {
    fgets(line,127,fptr);
    ++line_number;
  }
  
  // now put in the angle information
    for (i=0; i<parm_data->unique_angles_found; ++i)
    {
      if(fscanf(fptr, "%hd %hd %hd\n", &verify, &temp1, &temp2) != 3)
      {
        printf("ERROR reading data from parameter file for angle #%d.\n", i);
        free(line);
        return INVALID_DATA;
      }
      parm_data->angle_data[i].DO_FIT_KT = temp1;
      parm_data->angle_data[i].DO_FIT_THEQ = temp2;
    }
    
    // Skip 4 lines to get to the section with angles
    for (i=0; i<3; i++)
    {
      fgets(line,127,fptr);
      ++line_number;
    }
    
    // do the dihedrals
    for (i=0; i<parm_data->unique_dihedrals_found; ++i)
    {
      for (j=0; j<parm_data->dihedral_data[i].num_terms; ++j) {
        if (global_options->VERBOSITY >= HIGH) {
          printf("Reading in data for dihedral %s-%s-%s-%s term %d\n",
                parm_data->dihedral_data[i].atom_type1,
                parm_data->dihedral_data[i].atom_type2,
                parm_data->dihedral_data[i].atom_type3,
                parm_data->dihedral_data[i].atom_type4, j);
        }
        if (fscanf(fptr, "%hd %hd %hd %hd", &verify, &temp1, &temp2, &temp3) != 4) // number of arguments read
        {
          printf("ERROR reading data from parameter file for dihedral #%d term %d.\n", i,j);
          free(line);
          return INVALID_DATA;
        }
        parm_data->dihedral_data[i].term[j].DO_FIT_KP = temp1;
        parm_data->dihedral_data[i].term[j].DO_FIT_NP = temp2;
        if (parm_data->dihedral_data[i].improper == YES) // don't fit improper dihedral phase
          parm_data->dihedral_data[i].term[j].DO_FIT_PHASE = NO;
        else
          parm_data->dihedral_data[i].term[j].DO_FIT_PHASE = temp3;
      }
    }
  
  free(line);
  fclose(fptr);
  if (global_options->VERBOSITY>=MEDIUM)
    printf("   Prmtop     (info): Successfully read in old format saved parameter information\n");
  return SUCCESS;
}
