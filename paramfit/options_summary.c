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

/** @file options_summary.c
 * Prints out routines for printing options summaries at program start
 */

#include <stdio.h>
#include <string.h>
#include "function_def.h"

/**
 * Prints a summary of the options in the job control file
 * @param[in] global_options The global options structure
 * @param[in] parm_data Array of parameter files
 * @param[in] coords_datas Array of coodinate sets of input structures
 */
void print_job_control_summary(global_options_struct *global_options, parm_struct *parm_datas, coord_set *coords_datas)
{
  int i; /*Stores number returned by printf = no. chars actually printed*/
  int j; int p; // For memory totals
  int line_length;
  int mem_total;
  i=printf("        ---------------------------------------------------------------------\n");
  printf("        |                          OPTIONS SUMMARY                         |\n");
  printf("        |                          ---------------                          |\n");

  line_length=i-2;

  /*printf returns the number of characters actually printed so I will use that for
    my line completion command, storing the return value in i*/

  /*We cycle through for every available option */
  /*Fork here depending on whether we are doing a fit or creating QM input files*/
  print_open_line_box(&i);
  i+=printf("Summary of Run Type Options:");
  print_close_line_box(line_length-i);
  if (global_options->RUNTYPE==CREATE_INPUT)
  { /*We are creating input files from an mdcrd file*/
     print_open_line_box(&i);
     i+=printf("             Mode = CREATE_INPUT, ");
     print_close_line_box(line_length-i);
     print_open_line_box(&i);
     i+=printf("      File format = ");
     switch (global_options->QMFILEFORMAT)
     {
        case GAUSSIAN:
          i+=printf("GAUSSIAN");
          break;
        case GAMESS:
          i+=printf("GAMESS");
          break;
        case ADF:
          i+=printf("ADF");
          break;
        default:
          i+=printf("UNKNOWN FORMAT");
          break;
     }
     i+=printf(", Structures to write = %d",coords_datas[0].num_coords);
     print_close_line_box(line_length-i);
     print_open_line_box(&i);
     i+=printf("  Filename format = %s##%s",global_options->QMFILEOUTSTART,global_options->QMFILEOUTEND);
     print_close_line_box(line_length-i);
     print_open_line_box(&i);
     print_close_line_box(line_length-i);
     print_open_line_box(&i);
     i+=printf("Summary of structural info:");
     print_close_line_box(line_length-i);
     print_open_line_box(&i);
     i+=printf("   QMCHARGE = %-4d MULTIPLICITY = %-4d",global_options->QM_SYSTEM_CHARGE,global_options->QM_SYSTEM_MULTIPLICITY);
     print_close_line_box(line_length-i);
     print_open_line_box(&i);
     i+=printf("     NATOMS = %-4d   NATOMTYPES = %-4d    AMBERCHARGE = %-.4f",coords_datas[0].natoms, parm_datas[0].NTYPES,parm_datas[0].AMBER_SYSTEM_CHARGE);
     print_close_line_box(line_length-i);
     print_open_line_box(&i);
     i+=printf("     NBONDS = %-4d      NANGLES = %-4d     NDIHEDRALS = %-4d",parm_datas[0].NBONH+parm_datas[0].NBONA,
                                                                              parm_datas[0].NTHETH+parm_datas[0].NTHETA,
                                                                              parm_datas[0].NPHIH+parm_datas[0].NPHIA);
     print_close_line_box(line_length-i);
     print_open_line_box(&i);
     i+=printf(" NBONDTYPES = %-4d  NANGLETYPES = %-4d NDIHEDRALTYPES = %-4d",parm_datas[0].MUMBND,parm_datas[0].MUMANG,parm_datas[0].MPTRA);
     print_close_line_box(line_length-i);
     print_open_line_box(&i);
     print_close_line_box(line_length-i);
     print_open_line_box(&i);
     i+=printf("Estimate Memory Usage (per cpu):");
     print_close_line_box(line_length-i);
     print_open_line_box(&i);
     
     i+=printf(" Coordinate info will be read from disk as required.");
     print_close_line_box(line_length-i);
     mem_total=0;
     print_open_line_box(&i);
     i+=printf("               OPTION STORAGE = %-d bytes",global_options->mem_allocated);
     mem_total+=global_options->mem_allocated;
     print_close_line_box(line_length-i);
     print_open_line_box(&i);
     i+=printf("               PRMTOP STORAGE = %-d bytes",parm_datas[0].mem_allocated);
     mem_total+=parm_datas[0].mem_allocated;
     print_close_line_box(line_length-i);
     print_open_line_box(&i);
     int j; int p=0;
     for (j=0; j<global_options->num_prmtops; ++j)
       p+=coords_datas[j].mem_allocated;
     i+=printf("           COORDINATE STORAGE = %-d bytes",p);
     mem_total+=p;
     print_close_line_box(line_length-i);
     print_open_line_box(&i);
     i+=printf(" TOTAL ESTIMATED MEMORY USAGE = %-d bytes",mem_total);
     print_close_line_box(line_length-i);
  }
  else
  { /*we are doing a fit*/
    print_open_line_box(&i);
    if (global_options->ALGORITHM == SIMPLEX)
      i+=printf("  Run Mode = FIT, Minimiser = SIMPLEX");
    else if (global_options->ALGORITHM == GENETIC)
      i+=printf("  Run Mode = FIT, Minimiser = GENETIC");
    else if (global_options->ALGORITHM == BOTH)
      i+=printf("  Run Mode = FIT, Minimiser = GENETIC THEN SIMPLEX");
    else
      i+=printf("  Run Mode = FIT, Minimiser = UNKNOWN");     
    print_close_line_box(line_length-i);
    print_open_line_box(&i);
    if (global_options->FUNC_TO_FIT==SUM_SQUARES_AMBER_STANDARD)
      i+=printf("  Function to be Fit = SUM_SQUARES_AMBER_STANDARD");
    else if (global_options->FUNC_TO_FIT==AMBER_FORCES)
      i+=printf("  Function to be Fit = AMBER_FORCES");
    else
       i+=printf("  Function to be Fit = UNKNOWN");
    print_close_line_box(line_length-i);
    print_open_line_box(&i);
    print_close_line_box(line_length-i);
    print_open_line_box(&i);
    i+=printf("Terms to be fit:");
    print_close_line_box(line_length-i);
    print_open_line_box(&i);
    if (global_options->K_FIT==YES) {
      i+=printf("  K = 1, UNIQUE_BONDS = %d, UNIQUE ANGLES = %d",parm_datas[0].unique_bonds_found,parm_datas[0].unique_angles_found);
      print_close_line_box(line_length-i);
      print_open_line_box(&i);
      i+=printf("  UNIQUE DIHEDRALS = %d",parm_datas[0].unique_dihedrals_found);       
      i+=printf("  TOTAL DIHEDRAL TERMS = %d", parm_datas[0].unique_dihedral_terms);
    }
    else {
      i+=printf("  UNIQUE_BONDS = %d, UNIQUE ANGLES = %d, UNIQUE DIHEDRALS = %d",parm_datas[0].unique_bonds_found,parm_datas[0].unique_angles_found,parm_datas[0].unique_dihedrals_found);     
		}
    print_close_line_box(line_length-i);
    print_open_line_box(&i);
    i+=printf("  NBONDS = %d, NANGLES = %d, NDIHEDRALS = %d",parm_datas[0].NBONH+parm_datas[0].NBONA,
                                                             parm_datas[0].NTHETH+parm_datas[0].NTHETA,
                                                             parm_datas[0].NPHIH+parm_datas[0].NPHIA);     
    print_close_line_box(line_length-i);
    print_open_line_box(&i);
    i+=printf("                                Total dimensions of fit = %d",parm_datas[0].ndimensions);
    print_close_line_box(line_length-i);
    print_open_line_box(&i);
    print_close_line_box(line_length-i);
    print_open_line_box(&i);
    if (global_options->FUNC_TO_FIT==SUM_SQUARES_AMBER_STANDARD)
      i+=printf("Sample structures for least squares fit = %d",global_options->TOTAL_STRUCTURES);
    else if (global_options->FUNC_TO_FIT==AMBER_FORCES)
      i+=printf("Sample structures for forces fit = %d", global_options->TOTAL_STRUCTURES);
    else
      i+=printf("Sample structures for UNKNOWN fit = %d",global_options->TOTAL_STRUCTURES);
    print_close_line_box(line_length-i);
    
    print_open_line_box(&i);
    if (global_options->SORT_ENERGY==YES)
      i+=printf("Sorting input structures by ascending energy");
    print_close_line_box(line_length-i);
    
    if (parm_datas[0].ndimensions>3*global_options->TOTAL_STRUCTURES) {
      print_open_line_box(&i);
      i+=printf("!  WARNING - INSUFFICIENT DATA FOR ACCURATE FIT");
      print_close_line_box(line_length-i);
    }
    print_open_line_box(&i);
    print_close_line_box(line_length-i);
    print_open_line_box(&i);
    if (global_options->K_FIT==YES)
      i+=printf("Energy Correction Term (K) = TO BE FIT");
    else
      i+=printf("Energy Correction Term (K) = %10.8f",parm_datas[0].K);
    print_close_line_box(line_length-i);
    print_open_line_box(&i);
    print_close_line_box(line_length-i);
    if (global_options->ALGORITHM==SIMPLEX)
    {
      print_open_line_box(&i);
      i+=printf("          BONDFC_dx = %10.4f",global_options->BONDFC_dx);
      print_close_line_box(line_length-i);
      print_open_line_box(&i);
      i+=printf("          BONDEQ_dx = %10.4f",global_options->BONDEQ_dx);
      print_close_line_box(line_length-i);
      print_open_line_box(&i);
      i+=printf("         ANGLEFC_dx = %10.4f",global_options->ANGLEFC_dx);
      print_close_line_box(line_length-i);
      print_open_line_box(&i);
      i+=printf("         ANGLEEQ_dx = %10.4f",global_options->ANGLEEQ_dx);
      print_close_line_box(line_length-i);
      print_open_line_box(&i);
      i+=printf("      DIHEDRALBH_dx = %10.4f",global_options->DIHEDRALBH_dx);
      print_close_line_box(line_length-i);
      print_open_line_box(&i);
      i+=printf("       DIHEDRALN_dx = %10.4f",global_options->DIHEDRALN_dx);
      print_close_line_box(line_length-i);
      print_open_line_box(&i);
      i+=printf("       DIHEDRALG_dx = %10.4f",global_options->DIHEDRALG_dx);
      print_close_line_box(line_length-i);
      if (global_options->K_FIT==YES)
      {
        print_open_line_box(&i);
        i+=printf("             K_dx = %10.4f", global_options->K_dx);
        print_close_line_box(line_length-i);
      }
      print_open_line_box(&i);
      print_close_line_box(line_length-i);
      print_open_line_box(&i);
      i+=printf("Convergence requested to within %10.4E",global_options->CONV_LIMIT);
      print_close_line_box(line_length-i);
    }
    if (global_options->ALGORITHM==GENETIC)
    {
      print_open_line_box(&i);
      i+=printf("    OPTIMIZATIONS = %10d", global_options->NOPTIMIZATIONS);
      print_close_line_box(line_length-i);
      print_open_line_box(&i);
      i+=printf("  MAX GENERATIONS = %10d", global_options->MAX_GENERATIONS);
      print_close_line_box(line_length-i);
      print_open_line_box(&i);
      i+=printf("    MUTATION RATE = %10.2f", global_options->MUTATION_RATE);
      print_close_line_box(line_length-i);
    }
    print_open_line_box(&i);
    print_close_line_box(line_length-i);
    print_open_line_box(&i);
    i+=printf("Estimate Memory Usage (per cpu):");
    print_close_line_box(line_length-i);
    print_open_line_box(&i);
    i+=printf(" Coordinate info will be read from disk as required.");
    print_close_line_box(line_length-i);
    mem_total=0;
    print_open_line_box(&i);
    i+=printf("               OPTION STORAGE = %-d bytes",global_options->mem_allocated);
    mem_total+=global_options->mem_allocated;
    print_close_line_box(line_length-i);
    
    print_open_line_box(&i);
    p=0;
    for (j=0; j<global_options->num_prmtops; ++j) 
      p+= parm_datas[j].mem_allocated;
    i+=printf("               PRMTOP STORAGE = %-d bytes",p);
    mem_total += p;
    print_close_line_box(line_length-i);
    
    print_open_line_box(&i);
    p=0;
    for (j=0; j<global_options->num_prmtops; ++j)
      p+=coords_datas[j].mem_allocated;
    i+=printf("           COORDINATE STORAGE = %-d bytes",p);
    mem_total+=p;
    print_close_line_box(line_length-i);
    if (global_options->ALGORITHM==SIMPLEX)
    {
      print_open_line_box(&i);
      i+=printf("        SIMPLEX ARRAY STORAGE = %-d bytes",
      (int) (((parm_datas[0].ndimensions+1)*parm_datas[0].ndimensions+ (5*parm_datas[0].ndimensions) + 14)*sizeof(double)));
      /*SIMPLEX ROUTINE NEEDS (N+1)*N (doubles) Storage or the simplex array as well as (5*NDIMENSIONS)+14 doubles as scratch*/
      mem_total+=((parm_datas[0].ndimensions+1)*parm_datas[0].ndimensions+ (5*parm_datas[0].ndimensions) + 14)*sizeof(double);
      print_close_line_box(line_length-i);
    }
    if (global_options->ALGORITHM==GENETIC)
    {
      print_open_line_box(&i);
      i+=printf("             GA ARRAY STORAGE = %-d bytes",parm_datas[0].ndimensions*global_options->NOPTIMIZATIONS*sizeof(double));
      mem_total+=parm_datas[0].ndimensions*global_options->NOPTIMIZATIONS*sizeof(double);
      print_close_line_box(line_length-i);
    }
    if (global_options->CHECK_BOUNDS != NO)
    {
      print_open_line_box(&i);
      //i+=printf("      C
      print_close_line_box(line_length-i);
    }
    print_open_line_box(&i);
    i+=printf(" TOTAL ESTIMATED MEMORY USAGE = %-d bytes",mem_total);
    print_close_line_box(line_length-i);
    
  }

  print_open_line_box(&i);
  print_close_line_box(line_length-i);
  printf("        ---------------------------------------------------------------------\n\n");
  fflush(stdout); /*Flush the printf buffer*/
}
