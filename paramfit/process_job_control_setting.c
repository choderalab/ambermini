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

/** @file process_job_control_setting.c
 * Includes big function to process one line of a job control file.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "function_def.h"

/**
 * Processes one line of the job control file and sets corresponding options.
 * This function is called repeatedly for each line of the job control file until
 * all settings are processed.
 * @see read_job_control_file
 * 
 * @param[in] setting_line The raw line from the job control file as a character array
 * @param[in] length The length of this line
 * @param[in,out] number_settings The number of settings that have been set. Will be updated if successful.
 * @param[in,out] global_options The global options structure to input settings to
 * @param[in,out] parm_data Pointer to array of parm structs allocated to size 1, may have K set if specified in job control file
 * @param[in,out] coords_data Pointer to array of coords struct already allocated to size 1 that may have its number of
 * structures set if NSTRUCTURES is specified in the job control file, for single prmtop fits
 * @return Integer indicating success or various types of failure (invalid setting/data, alloc fail, etc)
 */
int process_job_control_setting(char *setting_line, int length, int *number_settings, global_options_struct *global_options, parm_struct *parm_data,
                                coord_set *coords_data)
{
  int temp_count;
  int var_count;
  char *test_string;
  char *data_string;
  int int_temp;
  double double_temp1;

  var_count=0;
  /*This actually processes the line in the settings line*/
  /*we read until we get an =*/
  while (*(setting_line+var_count)!='=') {
    /*check we haven't exceeded the line length, if we had then the line is invalid*/
    ++var_count;
    if (var_count>length) {
      return INVALID_LINE;
    }
  }
  /*var_count is currently 1 too big due to counting the = sign*/
  --var_count;
  /*we should now have a value of var_count telling us how many characters
    of *setting_line represent our variable name*/

  /*allocate memory for test string*/
  /*+1 since counter starts from zero and +1 again for null character*/
  test_string=(char *)malloc(var_count+2);
  if (test_string==NULL) {
    malloc_failure_char("process_setting", "test_string", (var_count+2));
    return ALLOC_FAIL;
  }
  /*allocate memory for data string*/
  /*(var_count+1 to allow for = sign which we don't want to store*/
  data_string=(char *)malloc(length-(var_count+1));
  if (data_string==NULL) {
    malloc_failure_char("process_setting", "data_string", (var_count+2));
    return ALLOC_FAIL;
  }
  /*copy in the variable name*/
  for (temp_count=0;temp_count<=var_count;temp_count++) {
    *(test_string+temp_count)=*(setting_line+temp_count);
  }
  *(test_string+var_count+1)=0;  /*add a null character to allow string comparisons*/

  /*copy in the data*/
  /*Note: it is possible that the last character of our data can actually
    contain a carriage return '\r' since we may have a dos formated text file that has
    \r\n at the end of each line rather than just \n or \n\r.

    So, we need to test for this and remove it if necessary
  */
  for (temp_count=(var_count+2);temp_count<length;temp_count++) {
      *(data_string+(temp_count-(var_count+2)))=*(setting_line+temp_count);
  }
  *(data_string+(temp_count-(var_count+2)))=0;  /*add a null character to allow string comparissons*/

  /*test for rogue carriage return*/
  if (*(data_string+(temp_count-(var_count+3)))=='\r') {
    if (global_options->VERBOSITY>=MEDIUM) {
      printf("!  Rogue carriage return found at end of test data array, char pos=%d\n",(temp_count-(var_count+1)));
      printf("!  Correcting error\n");
    }
    /*we have a rogue carriage return, remove it*/
    *(data_string+(temp_count-(var_count+3)))=0;
  }
  /*now we go through our list of variables looking for a match*/
  
  // RUNTYPE
  if (!(strcmp(test_string,"RUNTYPE")))
  { /*Valid settings for RUNTYPE are FIT or CREATE_INPUT*/
    /*Remember strcmp returns 0 on success*/
    if (!(strcmp(data_string,"FIT")))
      global_options->RUNTYPE=FIT;
    else if (!(strcmp(data_string,"CREATE_INPUT")))
      global_options->RUNTYPE=CREATE_INPUT;
    else if (!(strcmp(data_string,"SET_PARAMS")))
      global_options->RUNTYPE=SET_PARAMS;
    else {
      /*Invalid setting*/
      printf("!   Invalid setting for RUNTYPE\n");
      exit(UNKNOWN_OPT);
    }
  }
  // QMFILEFORMAT
  else if (!(strcmp(test_string,"QMFILEFORMAT")))
  { /*Valid settings for QMFILEFORMAT are GAUSSIAN, GAMESS, ADF,*/
    /*Remember strcmp returns 0 on success*/
    if (!(strcmp(data_string,"GAUSSIAN"))) {
      global_options->QMFILEFORMAT=GAUSSIAN;
      global_options->QM_FORCE_UNITS=HARTREE_BOHR;
    }
    else if (!(strcmp(data_string, "ADF"))) {
      global_options->QMFILEFORMAT=ADF;
    }
    else if (!(strcmp(data_string,"GAMESS"))) {
      global_options->QMFILEFORMAT=GAMESS;
    }
    else {
      /*Invalid setting*/
      printf("!   Invalid setting for QMFILEFORMAT\n");
      exit(UNKNOWN_OPT);
    }
  }
  // NSTRUCTURES
  else if (!(strcmp(test_string,"NSTRUCTURES"))) {
    if (sscanf(data_string, "%d", &int_temp) != 1) {
      printf("!   Invalid setting for NSTRUCTURES\n");
      exit(UNKNOWN_OPT);
    }
    /*&int_temp now contains the value for NSTRUCTURES*/
    /*Check it is within the limits of available options for NSTRUCTURES*/
    if (global_options->mdcrd_list!=NULL) {
      printf("\n  ERROR: Cannot specify mdcrd list %s and NSTRUCTURES!\n",global_options->mdcrd_list);
      free(test_string); free(data_string);
      exit(INVALID_DATA);
    } else if (!(coords_data)) {
      printf("\n  ERROR: NSTRUCTURES specified in job control, but no input coordinate set specified.\n\n");
      exit(INVALID_DATA);
    } else if (int_temp<1) {
      printf("\n  ERROR: Invalid setting for NSTRUCTURES- value must be >=1\n");
      exit(UNKNOWN_OPT);
    }
    else if (int_temp<20) { /*Print a warning about insufficient data*/
      printf("!  WARNING - %d STRUCTURES IS PROBABLY TOO SMALL FOR AN ACCURATE FIT\n",int_temp);
    }
    
    global_options->TOTAL_STRUCTURES=int_temp;
    coords_data[0].num_coords=int_temp;
  }
  // ALGORITHM (or deprecated FITTING_FUNCTION)
  else if (!(strcmp(test_string,"ALGORITHM")) || !(strcmp(test_string,"FITTING_FUNCTION")) )
  { 
    if (!(strcmp(data_string,"SIMPLEX"))) {
      global_options->ALGORITHM=SIMPLEX;
    } else if (!(strcmp(data_string, "GENETIC"))) {
      global_options->ALGORITHM=GENETIC;
    } else if (!(strcmp(data_string, "BOTH"))) { 
      global_options->ALGORITHM=BOTH;
    } else if (!(strcmp(data_string, "NONE"))) {
      global_options->ALGORITHM=NONE;
    } else {
      printf("!   Invalid setting for ALGORITHM\n");
      exit(UNKNOWN_OPT);
    }
  }
  // COORDINATE_FORMAT
  else if (!(strcmp(test_string,"COORDINATE_FORMAT"))) {
    if (!(strcmp(data_string,"TRAJECTORY"))) {
      global_options->COORD_FORMAT=TRAJECTORY;
    } else if (!(strcmp(data_string,"RESTART"))) {
      global_options->COORD_FORMAT=RESTART;
    } else {
      printf("!   Invalid setting for COORDINATE_FORMAT: \"%s\"\n",data_string);
      exit(UNKNOWN_OPT);
    }
  }
  // SCATTERPLOTS
  else if (!(strcmp(test_string,"SCATTERPLOTS"))) {
    printf("  Will print scatterplots!\n");
    global_options->SCATTERPLOTS=TRUE;
  }
  // FUNC_TO_FIT
  else if (!(strcmp(test_string,"FUNC_TO_FIT"))) { 
    if (!(strcmp(data_string,"SUM_SQUARES_AMBER_STANDARD"))) {
      global_options->FUNC_TO_FIT=SUM_SQUARES_AMBER_STANDARD;
    } else if (!(strcmp(data_string,"AMBER_FORCES"))) {
      printf("!  WARNING: Forces fitting is under development!\n");
      global_options->FUNC_TO_FIT=AMBER_FORCES;
    } else if (!(strcmp(data_string,"DIHEDRAL_LEAST_SQUARES"))) {
      global_options->FUNC_TO_FIT=DIHEDRAL_LEAST_SQUARES;
    } else {
      printf("!   Invalid setting for FUNC_TO_FIT: \"%s\"\n", data_string);
      exit(UNKNOWN_OPT);
    }
  }
  // FIT_PHASE
  else if (!(strcmp(test_string,"FIT_PHASE"))) {
    if (!(strcmp(data_string,"YES") && strcmp(data_string,"TRUE"))) {
      global_options->FIT_PHASE=YES;
    } else if (!(strcmp(data_string,"NO") && strcmp(data_string,"FALSE"))) {
      global_options->FIT_PHASE=NO;
    } else {
      printf("!   Invalid setting for FIT_PHASE: \"%s\"\n", data_string);
      exit(UNKNOWN_OPT);
    }
  }
  // WRITE_FRCMOD
  else if (!(strcmp(test_string,"WRITE_FRCMOD"))) {
    int_temp = strlen(data_string);
    if (int_temp > 0 && check_for_valid_filename(data_string,int_temp)==SUCCESS) {
      global_options->WRITE_FRCMOD = calloc( int_temp+1, sizeof(char));
      if (global_options->WRITE_FRCMOD == NULL) {
        exit(ALLOC_FAIL);
      }
      global_options->mem_allocated+=(int_temp+1);
      strcpy(global_options->WRITE_FRCMOD,data_string);
    } else {
      printf("!   Invalid setting for WRITE_FRCMOD\n");
      exit(UNKNOWN_OPT);
    }
  }
  // WRITE_ENERGY
  else if (!(strcmp(test_string, "WRITE_ENERGY"))) {
    int_temp = strlen(data_string);
    if (int_temp > 0 && check_for_valid_filename(data_string,int_temp)==SUCCESS) {
      global_options->WRITE_ENERGY = calloc( int_temp+1, sizeof(char));
      if (global_options->WRITE_ENERGY==NULL) {
        exit(ALLOC_FAIL);
      }
      global_options->mem_allocated+=(int_temp+1);
      strcpy(global_options->WRITE_ENERGY, data_string);
    } else {
      printf("!   Invalid setting for WRITE_ENERGY\n");
      exit(UNKNOWN_OPT);
    }
  }
  // SCNB
  else if (!(strcmp(test_string,"SCNB"))) {
    if (sscanf(data_string, "%lf", &double_temp1) != 1) {
      printf("!   Invalid setting for SCNB\n");
      exit(UNKNOWN_OPT);
    }
    global_options->SCNB=double_temp1;
  }
  // SCEE
  else if (!(strcmp(test_string,"SCEE"))) {
    if (sscanf(data_string, "%lf", &double_temp1) != 1) {
      printf("!   Invalid setting for SCEE\n");
      exit(UNKNOWN_OPT);
    }
    global_options->SCEE=double_temp1;
  }
  // QM_ENERGY_UNITS
  else if (!(strcmp(test_string,"QM_ENERGY_UNITS")) || 
           !(strcmp(test_string,"ENERGY_UNITS")))
  { 
    if (!(strcmp(data_string,"HARTREE"))) {
      global_options->QM_ENERGY_UNITS=HARTREE;
    } else if (!(strcmp(data_string,"KCALMOL"))) {
      global_options->QM_ENERGY_UNITS=KCALMOL;
    } else if (!(strcmp(data_string,"KJMOL"))) {
      global_options->QM_ENERGY_UNITS=KJMOL;
    } else {
      printf("!   Invalid setting for QM_ENERGY_UNITS\n");
      exit(UNKNOWN_OPT);
    }
  }
  // QM_FORCE_UNITS
  else if (!(strcmp(test_string,"QM_FORCE_UNITS")) || 
           !(strcmp(test_string,"FORCE_UNITS")))
    {
    if (!(strcmp(data_string,"HARTREE_BOHR"))) {
      global_options->QM_FORCE_UNITS=HARTREE_BOHR;
    } else if (!(strcmp(data_string,"KCALMOL_ANGSTROM"))) {
      global_options->QM_FORCE_UNITS=KCALMOL_ANGSTROM;
    } else {
      printf("!   Invalid setting for QM_FORCE_UNITS\n");
      exit(UNKNOWN_OPT);
    }
  }
  // PARAMETERS_TO_FIT 
  else if (!(strcmp(test_string,"PARAMETERS_TO_FIT"))) {
    if (!(strcmp(data_string,"DEFAULT"))) {
      global_options->PARAMETERS_TO_FIT=DEFAULT;
    } else if (!(strcmp(data_string,"LOAD"))) {
      global_options->PARAMETERS_TO_FIT=LOAD;
    } else if (!(strcmp(data_string, "K_ONLY"))) {
      global_options->PARAMETERS_TO_FIT=K_ONLY;
      global_options->K_FIT=YES;
    } else {
      printf("!   Invalid setting for PARAMETERS_TO_FIT\n");
      exit(UNKNOWN_OPT);
    }
  }
  // PARAMETER_FILE_NAME
    // specifies the file name to save the file to if there is a parameter file
    // to be loaded or saved.
  else if (!(strcmp(test_string,"PARAMETER_FILE_NAME"))) {
    int_temp = strlen(data_string);
    if (int_temp > 0 && check_for_valid_filename(data_string,int_temp)==SUCCESS) {
      global_options->PARAMETER_FILE_NAME = calloc( int_temp+1, sizeof(char));
      if (global_options->PARAMETER_FILE_NAME==NULL)
        return ALLOC_FAIL;
      global_options->mem_allocated+=(int_temp+1);
      strcpy(global_options->PARAMETER_FILE_NAME, data_string);
      if (global_options->VERBOSITY>=HIGH)
        printf (" Job control: Saved PARAMETER_FILE_NAME is: %s\n", global_options->PARAMETER_FILE_NAME);
    } else { // no filename specified
      printf("!   Invalid setting for PARAMETER_FILE_NAME\n");
      exit(UNKNOWN_OPT);
    }
  }
  // MUTATION_RATE
  else if (!(strcmp(test_string,"MUTATION_RATE"))) {
    if (sscanf(data_string, "%lf", &double_temp1) != 1 || double_temp1>1 || double_temp1<0) {
      printf("!   Invalid setting for MUTATION_RATE\n");
      exit(UNKNOWN_OPT);
    }
    if (double_temp1 == 0.0) {
      printf("!   MUTATION_RATE must be nonzero!\n");
      exit(INVALID_DATA);
    }
    global_options->MUTATION_RATE=double_temp1;
  }
  // PARENT_PERCENT
  else if (!(strcmp(test_string,"PARENT_PERCENT"))) {
    if (sscanf(data_string, "%lf", &double_temp1)!=1 || double_temp1>1 || double_temp1<0) {
      printf("!   Invalid setting for PARENT_PERCENT\n");
      exit(UNKNOWN_OPT);
    }
    global_options->PARENT_PERCENT=double_temp1;
  }
  // K
  else if (!(strcmp(test_string,"K"))) {
    // Cannot specify a global K value if fitting multiple prmtops
    if (global_options->prmtop_list!=NULL) {
      printf("\nERROR: cannot set K in Job Control file when fitting to multiple prmtops.\n");
      printf("       Run paramfit in single-prmtop mode and fit K for each prmtop individually.\n");
      printf("       Write the value of K for each prmtop following the filename in the prmtop list.\n");
      printf("       Example: test.prmtop 9000.0\n");
      free(test_string); free(data_string);
      exit(INVALID_DATA);
    }
    /*Valid settings for K are FIT or a float*
    /*Test first of all to see if it is FIT*/
    /*Remember strcmp returns 0 on success*/
    if (!(strcmp(data_string,"FIT"))) {
      global_options->K_FIT=YES;
    } else {
      if (sscanf(data_string, "%lf", &double_temp1) != 1) {
        printf("!   Invalid setting for K\n");
        exit(UNKNOWN_OPT);
      }
      parm_data[0].K=double_temp1;
      global_options->K_FIT=NO;
    }
  }
  // BONDFC_dx
  else if (!(strcmp(test_string,"BONDFC_dx"))) {
    if (sscanf(data_string, "%lf", &double_temp1) != 1) {
      printf("!   Invalid setting for BONDFC_dx\n");
      exit(UNKNOWN_OPT);
    }
    global_options->BONDFC_dx=double_temp1;
  }
  // BONDEQ_dx
  else if (!(strcmp(test_string,"BONDEQ_dx"))) {
    if (sscanf(data_string, "%lf", &double_temp1) != 1) {
      printf("!   Invalid setting for BONDEQ_dx\n");
      exit(UNKNOWN_OPT);
    }
    global_options->BONDEQ_dx=double_temp1;
  }
  // ANGLEFC_dx
  else if (!(strcmp(test_string,"ANGLEFC_dx"))) {
    if (sscanf(data_string, "%lf", &double_temp1) != 1) {
      printf("!   Invalid setting for ANGLEFC_dx\n");
      exit(UNKNOWN_OPT);
    }
    global_options->ANGLEFC_dx=double_temp1;
  }
  // ANGLEEQ_dx
  else if (!(strcmp(test_string,"ANGLEEQ_dx"))) {
    if (sscanf(data_string, "%lf", &double_temp1) != 1) {
      printf("!   Invalid setting for ANGLEEQ_dx\n");
      exit(UNKNOWN_OPT);
    }
    global_options->ANGLEEQ_dx=double_temp1;
  }
  // DIHEDRALBH_dx
  else if (!(strcmp(test_string,"DIHEDRALBH_dx"))) {
    if (sscanf(data_string, "%lf", &double_temp1) != 1) {
      printf("!   Invalid setting for DIHEDRALBH_dx\n");
      exit(UNKNOWN_OPT);
    }
    global_options->DIHEDRALBH_dx=double_temp1;
  }
  // DIHEDRALN_dx
  else if (!(strcmp(test_string,"DIHEDRALN_dx"))) {
    if (sscanf(data_string, "%lf", &double_temp1) != 1) {
      printf("!   Invalid setting for DIHEDRALN_dx\n");
      exit(UNKNOWN_OPT);
    }
    global_options->DIHEDRALN_dx=double_temp1;
  }
  // DIHEDRALG_dx
  else if (!(strcmp(test_string,"DIHEDRALG_dx"))) {
    if (sscanf(data_string, "%lf", &double_temp1) != 1) {
      printf("!   Invalid setting for DIHEDRALG_dx\n");
      exit(UNKNOWN_OPT);
    }
    global_options->DIHEDRALG_dx=double_temp1;
  }
  // K_dx
  else if (!(strcmp(test_string,"K_dx"))) {
    if (sscanf(data_string, "%lf", &double_temp1) != 1) {
      printf("!   Invalid setting for K_dx\n"); 
      exit(UNKNOWN_OPT);
    }
    global_options->K_dx=double_temp1;
  }
  // OPTIMIZATIONS
  else if (!(strcmp(test_string,"OPTIMIZATIONS"))) {
    if (sscanf(data_string, "%d", &int_temp) != 1) {
      printf("!   Invalid setting for OPTIMIZATIONS\n");
      exit(UNKNOWN_OPT);
    }
    if (int_temp<=0) {
      printf("!   Invalid setting for OPTIMIZATIONS\n");
      exit(INVALID_DATA);
    }
    global_options->NOPTIMIZATIONS=int_temp;
  }
  // NDIHEDRALS currently not supported
  /*
  else if (!(strcmp(test_string, "NDIHEDRALS"))) {
    if (sscanf(data_string, "%d", &int_temp) != 1) {
      printf("!   Invalid setting for NDIHEDRALS\n");
      exit(UNKNOWN_OPT);
    }
    if (int_temp<=0) {
      printf("!   Invalid setting for NDIHEDRALS\n");
      exit(INVALID_DATA);
    }
    global_options->NDIHEDRALS=int_temp;
  }
  */
  // SEARCH SPACE
  else if (!(strcmp(test_string,"SEARCH_SPACE"))) {
    if (sscanf(data_string, "%lf", &double_temp1) != 1) {
      printf("!   Invalid setting for SEARCH_SPACE\n");
      exit(UNKNOWN_OPT);
    }
    global_options->SEARCH_SPACE = double_temp1;
  }
  // MAX GENERATIONS
  else if (!(strcmp(test_string,"MAX_GENERATIONS"))) {
    if(sscanf(data_string, "%d", &int_temp) != 1) {
      printf("!   Invalid setting for MAX_GENERATIONS\n");
      exit(UNKNOWN_OPT);
    }
    if (int_temp < 0) {
      printf("!   MAX_GENERATIONS cannot be negative!\n");
      exit(INVALID_DATA);
    }
    global_options->MAX_GENERATIONS=int_temp;
  }
  // GENERATIONS TO CONVERGE
  else if (!(strcmp(test_string,"GENERATIONS_TO_CONV")) || !(strcmp(test_string,"GENERATIONS_TO_CONVERGE")) ) {
    if (sscanf(data_string, "%d", &int_temp) != 1) {
      printf("!   Invalid setting for GENERATIONS_TO_CONV\n");
      exit(UNKNOWN_OPT);
    }
    if (int_temp < 0) {
      printf("!   MAX_GENERATIONS cannot be negative!\n");
      exit(INVALID_DATA);
    }
    global_options->GENERATIONS_TO_CONVERGE=int_temp;
  }
  // CONVERGENCE LIMIT
  else if (!(strcmp(test_string,"CONV_LIMIT"))) {
    if (sscanf(data_string, "%lf", &double_temp1) != 1) {
      /*Conversion to double did not work = invalid data for setting*/
      printf("!   Invalid setting for CONV_LIMIT\n");
      exit(UNKNOWN_OPT);
    }
    if (double_temp1 < 0.) {
      printf("!   CONV_LIMIT cannot be negative!\n");
      exit(INVALID_DATA);
    }
    global_options->CONV_LIMIT=double_temp1;    
  }
  // CHECK BOUNDS
  else if (!(strcmp(test_string, "CHECK_BOUNDS"))) {
    if (!(strcmp(data_string,"ON"))) {
      global_options->CHECK_BOUNDS=YES;
    } else if (!(strcmp(data_string, "WARN"))) {
      global_options->CHECK_BOUNDS=WARN;
    } else {
      printf("!   Invalid setting for CHECK_BOUNDS\n");
      exit(UNKNOWN_OPT);
    }
  }
  // SORT MDCRDS
  else if (!(strcmp(test_string,"SORT_MDCRDS"))) {
    if (!(strcmp(data_string,"ON")) || !(strcmp(data_string,"YES"))) {
      global_options->SORT_ENERGY=YES;
    } else if (!(strcmp(data_string,"OFF")) || !(strcmp(data_string,"NO"))) {
      global_options->SORT_ENERGY=NO;
    } else {
      printf("!   Invalid setting for SORT_MDCRDS\n");
      free(test_string);
      free(data_string);
      exit(UNKNOWN_OPT);
    }
  }
  // ANGLE LIMIT
  else if (!(strcmp(test_string, "ANGLE_LIMIT"))) {
    if (sscanf(data_string, "%lf", &double_temp1) != 1) {
      printf("!   Invalid setting for ANGLE_LIMIT\n");
      exit(UNKNOWN_OPT);
    }
    global_options->ANGLE_LIMIT=double_temp1;
  }
  // BOND_LIMIT
  else if (!(strcmp(test_string, "BOND_LIMIT"))) {
    if (sscanf(data_string, "%lf", &double_temp1) != 1) {
      printf("!   Invalid setting for BOND_LIMIT\n");
      exit(UNKNOWN_OPT);
    }
    global_options->BOND_LIMIT=double_temp1;
  }
  // DIHEDRAL_SPAN
  else if (!(strcmp(test_string, "DIHEDRAL_SPAN"))) {
    if (sscanf(data_string, "%lf", &double_temp1) != 1) {
      printf("!   Invalid setting for DIHEDRAL_SPAN\n");
      exit(UNKNOWN_OPT);
    }
    global_options->DIHEDRAL_SPAN=(int)double_temp1;
  }
  // QMHEADER
  else if (!(strcmp(test_string, "QMHEADER"))) {
    int_temp=strlen(data_string);
    if (int_temp > 0 && check_for_valid_filename(data_string,int_temp)==SUCCESS) {
      global_options->QMHEADER = calloc( int_temp+1, sizeof(char));
      if (global_options->QMHEADER==NULL)
        return ALLOC_FAIL;
      global_options->mem_allocated+=(int_temp+1);
      strcpy(global_options->QMHEADER, data_string);
    }
    else {
      printf("!   Invalid setting for QMHEADER\n");
      exit(UNKNOWN_OPT);
    }
  }
  // GENERATIONS_TO_SIMPLEX
  else if (!(strcmp(test_string,"GENERATIONS_TO_SIMPLEX"))) {
    if (sscanf(data_string, "%d", &int_temp) != 1) {
      printf("! Invalid setting for GENERATIONS_TO_SIMPLEX\n");
      exit(UNKNOWN_OPT);
    }
    global_options->GENERATIONS_TO_SIMPLEX=int_temp;
    if (global_options->GENERATIONS_TO_SIMPLEX<0) {
      printf("! GENERATIONS_TO_SIMPLEX cannot be negative!\n");
      exit(INVALID_DATA);
    }
  }
  // GENERATIONS_WITHOUT_SIMPLEX
  else if (!(strcmp(test_string,"GENERATIONS_WITHOUT_SIMPLEX"))) {
    if (sscanf(data_string, "%d", &int_temp) != 1) {
      printf("! Invalid setting for GENERATIONS_WITHOUT_SIMPLEX\n");
      exit(UNKNOWN_OPT);
    }
    global_options->GENERATIONS_WITHOUT_SIMPLEX=int_temp;
    if (global_options->GENERATIONS_WITHOUT_SIMPLEX<0) {
      printf("! GENERATIONS_WITHOUT_SIMPLEX cannot be negative!\n");
      exit(INVALID_DATA);
    }
  }
  // QM_SYSTEM_CHARGE
  else if (!(strcmp(test_string,"QM_SYSTEM_CHARGE"))) {
    if (sscanf(data_string, "%d", &int_temp) != 1) {
      /*Conversion to int did not work = invalid data for setting*/
      printf("!   Invalid setting for QM_SYSTEM_CHARGE\n");
      exit(UNKNOWN_OPT);
    }
    /*Check it is sensible*/
    if (abs(int_temp)>10) {
      printf("!  WARNING - A QM_SYSTEM_CHARGE OF %d SEEMS STRANGE. CONTINUING ANYWAY.\n",int_temp);
    }
    global_options->QM_SYSTEM_CHARGE=int_temp;
  }
  // QM_SYSTEM_MULTIPLICITY
  else if (!(strcmp(test_string,"QM_SYSTEM_MULTIPLICITY"))) {
    if (sscanf(data_string, "%d", &int_temp) != 1) {
      printf("!   Invalid setting for QM_SYSTEM_MULTIPLICITY\n");
      exit(UNKNOWN_OPT);
    }
    if (int_temp>5) {
      printf("!  WARNING - A QM_SYSTEM_MULTIPLICITY OF %d SEEMS STRANGE. CONTINUING ANYWAY.\n",int_temp);
    }
    else if (int_temp<1) {
      printf("!   QM_SYSTEM_MULTIPLICITY cannot be less than 1!\n");
      exit(INVALID_DATA);
    }
    global_options->QM_SYSTEM_MULTIPLICITY=int_temp;
  }
  // QMFILEOUTSTART
  else if (!(strcmp(test_string,"QMFILEOUTSTART"))) {
     /*check for invalid filename characters*/
     int_temp = strlen(data_string);
     if(int_temp > 0 && check_for_valid_filename(data_string,int_temp)==SUCCESS) {
       global_options->QMFILEOUTSTART = (char *)malloc(int_temp+1);
       if (global_options->QMFILEOUTSTART == NULL) {
         exit(ALLOC_FAIL);
       }
       global_options->mem_allocated+=(int_temp+1);
       strcpy(global_options->QMFILEOUTSTART,data_string);
     } else {
      printf("!   Invalid setting for QMFILEOUTSTART\n");
        exit(UNKNOWN_OPT);
     }
  }
  // QMFILEOUTEND
  else if (!(strcmp(test_string,"QMFILEOUTEND"))) {
     /*check for invalid filename characters*/
     int_temp = strlen(data_string);
     if(int_temp > 0 && check_for_valid_filename(data_string, int_temp)==SUCCESS) {
       global_options->QMFILEOUTEND = (char *)malloc(int_temp+1);
       if (global_options->QMFILEOUTEND == NULL)
         return ALLOC_FAIL;
       global_options->mem_allocated+=(int_temp+1);
       strcpy(global_options->QMFILEOUTEND,data_string);
     } else {
        printf("!   Invalid setting for QMFILEOUTEND\n");
        free(test_string);
        free(data_string);    
        test_string = NULL;
        data_string = NULL;
        exit(UNKNOWN_OPT);
     }
  }
  else {
      /*Unknown variable name*/
      printf("!   Unknown variable \"%s\". Check spelling and paramfit version.\n", test_string);
      free(test_string);
      free(data_string);
      exit(UNKNOWN_OPT);
  }
  /*If we succesfully set a value then we must increase number_settings by 1*/
  *number_settings=*number_settings+1;
  free(test_string);
  free(data_string);
  return SUCCESS;
}

