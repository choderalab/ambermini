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

/** @file process_command_line.c
 * Utilities for processing command line options
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "function_def.h"

#define MAX_CMDLINE_OPTIONS 13  /**< Prog name + 12 options*/
/**
 * Processes all command line options.
 * Possible options combinations are:
 * paramfit -i [Job control file] -p [prmtop] -c [mdcrd] -q [quantum data] -d [ON/OFF/DEBUG] --random-seed [seed] 
 * paramfit -i [Job control file] -pf [prmtop list] -cf [mdcrd list] -qf [quantum list] -d [ON/OFF/DEBUG] --random-seed [seed] 
 * @param[in] argv Pointer to array with all options, argv[0] is the name of the program (paramfit)
 * @param[in,out] global_options Global options structure that will be updated
 * @param[in,out] parm_datas Unallocated array of parm structs, will be updated with initial blank one if -p specified
 * @param[in,out] coords_datas Unallocated of coords structs, will be updated with initial blank one if -c specified
 * @return Integer indicating success, failure, options problems, etc.
 */ 
int process_command_line(int argc, char **argv, global_options_struct *global_options, parm_struct **parm_datas, coord_set **coords_datas)
{
  int i;
  /*Step 1, check to see if too many options have been specified */
  if (argc > MAX_CMDLINE_OPTIONS)
  {
     printf("*** ERROR - Too many command line options specified.\n");
     command_line_help();
     return TOO_MANY_OPT;
  }

/*If the user specified no command line options I will switch on diagnostics
and inform the user that we are using default options*/
  if (argc<2)
  {
      printf("\n!  No command line options given - LOADING DEFAULTS\n");
      printf("!  Setting verbosity to medium.\n\n");
      global_options->VERBOSITY=MEDIUM;
      return SUCCESS;
  }

  if (!(strcmp(argv[1],"/?")) || !(strcmp(argv[1],"-?")) || !(strcmp(argv[1],"/help")) || !(strcmp(argv[1],"-help")) || !(strcmp(argv[1],"--help")))
  {
     /*User has requested help */
     command_line_help();
     return CMD_HELP_REQ;
  }
  if (!(strcmp(argv[1],"/history")) || !(strcmp(argv[1],"/HISTORY")) || !(strcmp(argv[1],"/History")))
  {
     /*User has requested program history */
     print_program_history();
     return HIST_REQ;
  }
  /*step 2, check each command line option to see if it is valid*/
  for (i=1; i<argc; i++)
  {
     /*Command line option 1 = -i Job_Control.in*/
     if(!(strcmp(argv[i],"-i")))
     {
        /*option is for job_control_filename*/
        ++i; /*Increment i to get next command line option */

        if (i >= argc) /*User didn't specify anything after the -i */
        {
          printf("*** ERROR - Unknown command line option.\n");
          command_line_help();
          return UNKNOWN_OPT;
        }
        /* global_options->job_control_filename has not been allocated by default so we alloc it.
           so we have to use realloc */
        global_options->job_control_filename = (char *)malloc( ( strlen(argv[i])+1 ) * sizeof(char) );
        if (global_options->job_control_filename == NULL)
        {
           malloc_failure_char("process_command_line", "global_options->job_control_filename", (strlen(argv[i])+1));
           return ALLOC_FAIL;
        }
        global_options->mem_allocated+=(((strlen(argv[i])+1)*sizeof(char)));
        strcpy(global_options->job_control_filename,argv[i]);
     }
     // Single prmtop filename- allocate the single prmtop here
     else if(!(strcmp(argv[i],"-p"))) {
       ++i; /*Increment i to get next command line option */

       if (i >= argc) { // User didn't specify anything after the -p
          printf("*** ERROR - Unknown command line option.\n");
          command_line_help();
          return UNKNOWN_OPT;
       }
       // Check -pf and -cf wasn't either wasn't also specified first
       if (global_options->prmtop_list || global_options->mdcrd_list) {
         printf("ERROR! Cannot specify both -p and -pf or -cf\n");
         command_line_help();
       }
       *parm_datas=(parm_struct*)malloc(sizeof(parm_struct));
       (*parm_datas)[0].filename = (char*)malloc(strlen(argv[i])+1);
      if (!(*parm_datas)[0].filename) {
        malloc_failure_char("process_command_line", "(*parm_datas)[0].filename", (strlen(argv[i])+1));
        return ALLOC_FAIL;
      }
      strcpy((*parm_datas)[0].filename, argv[i]);
      global_options->num_prmtops=1;
      global_options->mem_allocated+=(strlen(argv[i])+1);
     }
     // Command line option for file containing list of prmtops
     else if (!(strcmp(argv[i],"-pf"))) {
       ++i;
       if (i>=argc) {
         printf("*** ERROR - Unknown command line option.\n");
         command_line_help();
         return UNKNOWN_OPT;
       }
       // check -p or -c wasn't also specified first
       if (*parm_datas || *coords_datas) {
         printf("ERROR! Cannot specify both -p or -c and -pf\n");
         command_line_help();
         return INVALID_DATA;
       }
       global_options->prmtop_list = (char*)malloc(strlen(argv[i])+1);
       if (global_options->prmtop_list==NULL) {
         malloc_failure_char("process_command_line", "global_options->prmtop_list", (strlen(argv[i])+1));
         return ALLOC_FAIL;
       }
       global_options->mem_allocated+=(strlen(argv[i])+1);
       strcpy(global_options->prmtop_list,argv[i]);
     }
     // Command line options for file containing list of mdcrds
     else if (!(strcmp(argv[i],"-cf"))) {
       ++i;
       if (i>=argc) {
         printf("*** ERROR - Unknown command line option.\n");
         command_line_help();
         return UNKNOWN_OPT;
       }
       // Check -p and -c haven't been set
       if (*parm_datas || *coords_datas) {
         printf("ERROR! Cannot specify -cf and -p or -c!\n");
         command_line_help();
         return INVALID_DATA;
       }
       global_options->mdcrd_list = (char*)malloc(strlen(argv[i])+1);
       if (!global_options->mdcrd_list) {
         malloc_failure_char("process_command_line", "global_options->mdcrd_list", strlen(argv[i])+1);
         return ALLOC_FAIL;
       }
       global_options->mem_allocated+=strlen(argv[i])+1;
       strcpy(global_options->mdcrd_list,argv[i]);
     }
     // Single coordinate set- allocate the set and set its name
     else if(!(strcmp(argv[i],"-c"))) {
       ++i;
       if (i >= argc) {
          printf("*** ERROR - Unknown command line option.\n");
          command_line_help();
          return UNKNOWN_OPT;
       }
       // Check -pf and -cf unset
       if (global_options->prmtop_list || global_options->mdcrd_list) {
         printf("ERROR! Cannot specify -c and -pf or -cf\n");
         command_line_help();
         return INVALID_DATA;
       }
       // Allocate single coord set if we haven't done so yet (another option could have done this already)
       if (!*coords_datas) *coords_datas=(coord_set*)malloc(sizeof(coord_set));
       if (!*coords_datas) {
         printf("ERROR! Failed to allocate %lu bytes for *coords_datas\n", sizeof(coord_set));
         return ALLOC_FAIL;
       }
       (*coords_datas)[0].filename = (char*)malloc(strlen(argv[i])+1);
       if (!(*coords_datas)[0].filename) {
         malloc_failure_char("process_command_line", "(*coords_datas)[0].filename", strlen(argv[i])+1);
         return ALLOC_FAIL;
       }
       global_options->mem_allocated+=(strlen(argv[i])+1); 
       strcpy((*coords_datas)[0].filename, argv[i]);
       (*coords_datas)[0].energy_filename=NULL;
     }
     else if(!(strcmp(argv[i],"-q")))
     {   /* Option 4 = -q QM_data.dat */
       /*option is for energy filename*/
       ++i; /*Increment i to get next command line option */

       if (i >= argc) /*User didn't specify anything after the -q */ {
         printf("*** ERROR - Unknown command line option.\n");
         command_line_help();
         return UNKNOWN_OPT;
       }
       // Check -pf and -cf unset
       if (global_options->prmtop_list || global_options->mdcrd_list) {
         printf("ERROR! Cannot specify -q and -pf or -cf\n");
         command_line_help();
         return INVALID_DATA;
       }
       // Allocate single coord set if we haven't done so yet (another option could have done this already)
       if (!*coords_datas) *coords_datas=(coord_set*)malloc(sizeof(coord_set));
       if (!*coords_datas) {
         printf("ERROR! Failed to allocate %lu bytes for *coords_datas\n", sizeof(coord_set));
         return ALLOC_FAIL;
       }
       (*coords_datas)[0].energy_filename=(char*)malloc(strlen(argv[i])+1);
       strcpy((*coords_datas)[0].energy_filename, argv[i]);
       global_options->mem_allocated+=strlen(argv[i]+1);
     }
     else if(!(strcmp(argv[i],"-v")))
     { /* Option 5 = -d Verbosity setting */
        ++i;
        if (i >= argc) /*User didn't specify anything after the -v */
        {
          printf("*** ERROR - Unknown command line option.\n");
          command_line_help();
          return UNKNOWN_OPT;
        }             
        /*Valid options are HIGH, MEDIUM, LOW, or lower case versions */
        if ((strcmp(argv[i],"HIGH")) && (strcmp(argv[i],"MEDIUM")) && (strcmp(argv[i],"LOW")) && 
            (strcmp(argv[i],"high")) && (strcmp(argv[i],"medium")) && (strcmp(argv[i],"low")))
        {
          printf("*** ERROR - Unknown command line option.\n");
          command_line_help();
          return UNKNOWN_OPT; /*User specified unknown debug level*/
        }
        /*determine what level the user has specified*/
        if (!(strcmp(argv[i],"LOW"))||!(strcmp(argv[i],"low")))
        {
          global_options->VERBOSITY=LOW;
        }
        else if (!(strcmp(argv[i],"MEDIUM"))||!(strcmp(argv[i],"medium")))
        {
          global_options->VERBOSITY=MEDIUM;
        }
        else if (!(strcmp(argv[i],"HIGH"))||!(strcmp(argv[i],"high")))
        {
          global_options->VERBOSITY=HIGH;
        }
        else
        {
          /*Unknown option*/
          printf("*** ERROR - Unknown command line option.\n");
          command_line_help();
          return UNKNOWN_OPT;
        }
      }
      else if (!(strcmp(argv[i], "--random-seed")))
      {
        ++i; 
        if (i>=argc)
        {
          printf("*** ERROR - Unknown command line option.\n");
          command_line_help();
          return UNKNOWN_OPT;
        }
        if (sscanf(argv[i], "%d", &global_options->RANDOM_SEED) != 1) 
        {
          // Invalid data for the random seet
          printf("*** ERROR - invalid random seed.\n");
          command_line_help();
          return UNKNOWN_OPT;
        }
      }
      else if (!(strcmp(argv[i], "-d")))
      {
        // support for deprecated -d option
        ++i;
        if (i >= argc)
        {
          printf("*** ERROR - Unknown command line option.\n");
          command_line_help();
          return UNKNOWN_OPT;
        }
        if(!(strcmp(argv[i],"ON")))
          global_options->VERBOSITY=MEDIUM;
        else if (!(strcmp(argv[i],"OFF")))
          global_options->VERBOSITY=LOW;
        else if (!(strcmp(argv[i],"DEBUG")))
          global_options->VERBOSITY=HIGH;
        else
        {
          printf("*** ERROR - Unknown command line option.\n");
          command_line_help();
          return UNKNOWN_OPT;
        }
        printf("!  WARNING - The -d command line parameter is deprecated. Use the -v option instead.\n");
      }
      else
      {
        printf("*** ERROR - Unknown command line option.\n");
        command_line_help();
        return UNKNOWN_OPT; /*Unknown command line option*/
      }
  }
  
  // Do a few sanity checks to ensure there is a valid combination of options:
  // If job control file specified, we need either a prmtop or a prmtop list
  if (global_options->job_control_filename && 
      (!global_options->prmtop_list && !(*parm_datas)) ) {
    printf("*** ERROR - Cannot specify a job control file without a topology file to operate on.\n");
    printf("            Requires -p or -pf arguments\n");
    command_line_help();
    return INVALID_DATA;
  }
    

  return SUCCESS;
}




