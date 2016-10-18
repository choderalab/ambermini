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

/** @file error_messages.c
 * Contains a number of functions for printing out error messages.
 * These functions are usually called when a failure occurs, and they
 * print out more information about the failure and then exit.
 */

#include <stdio.h>
#include <stdlib.h>
#include "function_def.h"

/**
 * Processes the return value of a function and exits if an error
 * @param[in] err_code The return value from the function
 * @param[in] VERBOSITY How verbose the program is
 */
void process_retval(int err_code, verbosity_t VERBOSITY)
{
  if (err_code != SUCCESS)
  {
    /*Failure or abort flag due to printing of program history / help etc..*/
    /*Lower (less -ve) than ABORT = not a failure, an abort signal*/
    if (VERBOSITY>=HIGH && err_code <= ABORT)
      printf("*** PROGRAM ABORT CODE: %d\n",err_code);
    else if (VERBOSITY>=HIGH)
      printf("*** PROGRAM ERROR CODE: %d\n",err_code);
    
    exit(err_code);
  }  
}

/**
 * Prints a standard malloc failure message for char data types then exits.
 * Includes the routine name, variable name, and number of bytes requested
 * for char data types
 * @param[in] routine The function where the failure occured
 * @param[in] var_name The variable that failed to be allocated
 * @param[in] chars_requested The number of characters requested to be allocated
 */ 
void malloc_failure_char(char *routine, char *var_name, int chars_requested)
{
   printf("\n*** FATAL ERROR IN %s FUNCTION\n",routine);
   printf("m*** COMMAND WAS: %s = malloc(%d * sizeof(char));\n",var_name,chars_requested);
   printf("*** STATUS - FAILED TO ALLOCATE %d bytes\n",chars_requested);
  exit(ALLOC_FAIL);
}

/**
 * Prints out a standard malloc failure message for int data types then exits.
 * Includes the routine name, variable name, and number of bytes requested
 * for int data types.
 * @param[in] routine The function where the failure occured
 * @param[in] var_name The variable that failed to be allocated
 * @param[in] ints_requested The number of integers that were to be allocated
 */
void malloc_failure_int(char *routine, char *var_name, int ints_requested)
{
   printf("\n*** FATAL ERROR IN %s FUNCTION\n",routine);
   printf("*** COMMAND WAS: %s = malloc(%d * sizeof(int));\n",var_name,ints_requested);
   printf("*** STATUS - FAILED TO ALLOCATE %d bytes\n",(int)(ints_requested*sizeof(int)));
   exit(ALLOC_FAIL);
}

/**
 * Prints out a standard malloc failure message for short int data types then exits.
 * Includes the routine name, variable name, and number of bytes requested for short
 * int data types.
 * @param[in] routine The function where the failure occured
 * @param[in] var_name The variable that failed to be allocated
 * @param[in] short_ints_requested The number of short ints that were to be allocated
 */
void malloc_failure_short_int(char *routine, char *var_name, int short_ints_requested)
{
   printf("\n*** FATAL ERROR IN %s FUNCTION\n",routine);
   printf("*** COMMAND WAS: %s = malloc(%d * sizeof(short int));\n",var_name,short_ints_requested);
   printf("*** STATUS - FAILED TO ALLOCATE %d bytes\n",(int)(short_ints_requested*sizeof(short int)));
   exit(ALLOC_FAIL);
}

/**
 * Prints out a standard malloc failure message for double data types then exits.
 * Includes the routine name, variable name, and number of bytes requested for double
 * data types.
 * @param[in] routine The function where the failure occured
 * @param[in] var_name The variable that failed to be allocated
 * @param[in] doubles_requested The number of doubles that were to be allocated
 */
void malloc_failure_double(char *routine, char *var_name, int doubles_requested)
{
   printf("\n*** FATAL ERROR IN %s FUNCTION\n",routine);
   printf("*** COMMAND WAS: %s = malloc(%d * sizeof(double));\n",var_name,doubles_requested);
   printf("*** STATUS - FAILED TO ALLOCATE %d bytes\n",(int)(doubles_requested*sizeof(double)));
   exit(ALLOC_FAIL);
}

/**
 * Prints a standard file open failure message then exits.
 * Includes the function name and the file that was attempted to be opened.
 * @param[in] routine The function where the failure occured
 * @param[in] var_name The filename that could not be opened
 */
void file_open_failure(char *routine, char *var_name)
{
   printf("\n*** FATAL ERROR IN %s FUNCTION\n",routine);
   printf("*** STATUS - FAILED TO OPEN FILE: '%s'\n",var_name);
   exit(FILE_OPEN_FAIL);
}
