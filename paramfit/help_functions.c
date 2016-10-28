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

/*process_command_line.c*/
/*This module contains a number of Help routines for printing out help messages*/

#include <stdio.h>

#include "function_def.h"

/*This module prints out the command line help when the user either requests
  help with /? -? /help or -help or if the user provides unknown command line
  arguments */

void command_line_help()
{
  printf("Usage is:\n");
  printf("  paramfit -i Job_Control.in -p prmtop -c mdcrd -q QM_data.dat -v [LOW/MEDIUM/HIGH]\n");
  printf("                                       --- OR --- \n");
  printf("  paramfit -i Job_Control.in -pf prmtop_list -cf mdcrd_list -v [LOW/MEDIUM/HIGH]\n"); 
  printf("\n");
  printf("Valid switches include:\n");
  printf("     -i Job_Control.in\tJob control file location (mandatory)\n");
  printf("     -p prmtop\tParameter file location --OR--\n");
  printf("    -pf prmtop list\tList of multiple parameter files to use, their K values\n");
  printf("     -c mdcrd\tCoordinate file location --OR--\n");
  printf("    -cf mdcrd list\tList of mdcrd files, number of structres, qm files location\n");
  printf("     -q QM_data.dat\tList of quantum energies (for single fits only)\n");
  printf("     -v MEDIUM\n");
  printf("     --random-seed 0 (for debugging only, no default value)\n\n");
  printf("     /history prints program development history\n\n");
  printf("     For HELP please see the documentation\n");
  printf("\n");
}
