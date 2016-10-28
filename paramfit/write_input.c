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
/** @file write_input.c
 * Contains routines for writing output files.
 * Contains routines for writing all quantum input files
 * as well as various other outputs from fitting
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "function_def.h"
#include "constants.h"

/**
 * Writes a gaussian input file for one structure.
 * Should not be called directly. 
 * Uses a header file containing all of the quantum options that
 * the atom coordinates are appended to, so this can be used to make input
 * files for force or energy calculations.
 * 
 * @see write_input
 * @param[in] global_options The global options structure
 * @param[in] parm_data The parameter file
 * @param[in] current_struct Pointer to coordinate set to use
 * @param[in] num The number of the structure to write from the coordinate set
 * @param[in,out] fptr Pointer to the file to write to
 * @return Integer indicating success or failure
 */
int write_input_gaussian(global_options_struct *global_options, parm_struct *parm_data, coord_set *current_struct, int num, FILE *fptr)
{

  int i,j,k;
  int atomic_number;
  /*The output file should already be open here so we just start writing the data*/
  /* First open the header file and copy it to the top of the input file */
  if (global_options->QMHEADER)
  {
    FILE *header = fopen(global_options->QMHEADER, "r");
    if (!header) {
      printf("! ERROR opening QM header file at %s\n", global_options->QMHEADER);
      return INVALID_DATA;
    }

    // Get the size of the header
    fseek(header, 0, SEEK_END);
    long size = ftell(header);
    rewind(header);
    
    // allocate buffer, read in data, and write it to the file
    unsigned char *buffer=(unsigned char *)malloc((size_t)size);
    fread(buffer, size, 1, header);
    fwrite(buffer, size, 1, fptr);
    
    // clean up
    free(buffer);
    fclose(header);
  }
  
  /*Now in gaussian we need to write a blank line, a title and then a blank line, use the current_struct for the title*/
  fprintf(fptr,"\n");
  fprintf(fptr,"Structure: %d\n",num);
  fprintf(fptr,"\n");
  /*Now the charge and multiplicity*/
  fprintf(fptr,"  %d  %d\n",global_options->QM_SYSTEM_CHARGE,global_options->QM_SYSTEM_MULTIPLICITY);

  /*Loop over all atoms*/
  for (i=0;i<current_struct->natoms;++i)
  {
    /*Stage 1 - get the element*/
    /*we really only have 2 ways to do this, either from the mass or from the first letter of the
      atom name - both are unreliable*/
    /*Lets just use the first letter of the atom name - note, this can give misreads with for example
      N and Na, we can try and correct for these by checking the mass.
    */
    atomic_number=find_atomic_number_from_parm(parm_data, i+1);
    if (atomic_number==UNKNOWN_ELEMENT)
    {
      /*Element not recognised*/
      printf("\n*** ERROR IN WRITE_INPUT FROM CALL TO FIND_ATOMIC_NUMBER_FROM_PARM\n");
      printf("*** UNKNOWN ELEMENT FOR ATOM NAME: %s WITH MASS: %f\n",parm_data->atom[i].igraph, parm_data->atom[i].amass);
      printf("*** CHECK elements.c\n");
      return (UNKNOWN_ELEMENT);
    }
    /*Now write the element as its symbol*/
    fprintf(fptr,"  ");
    print_atomic_number_as_symbol(fptr, atomic_number);
    fprintf(fptr,"   ");
    /*Now the coordinates*/
    fprintf(fptr,"%10.6f %10.6f %10.6f\n",
	    current_struct->struc[num].x_coord[i], current_struct->struc[num].y_coord[i], current_struct->struc[num].z_coord[i]);
  }
  fprintf(fptr,"\n");
  
  // Now write connectivity data
  for (i=0; i<current_struct->natoms; ++i) {
    if (i>0) fprintf(fptr,"%d ", i);
    for (j=0; j<parm_data->unique_bonds_found; ++j) {
      for (k=0; k<parm_data->bond_data[j].number; ++k) {
        if (parm_data->bond_data[j].atom1[k] == i) {
            fprintf(fptr, "%d 1.0 ",  parm_data->bond_data[j].atom2[k]);
        }
        /*
        if (parm_data->bond_data[j].atom2[k] == i) {
            fprintf(fptr, "%d 1.0 ",  parm_data->bond_data[j].atom1[k]+1);
        }
        if (parm_data->bond_data[j].atom2[k] == i) {
          fprintf(fptr, "%d %d 1.0 ", i, parm_data->bond_data[j].atom1[k]);
          write=TRUE;
        }
        */
    } }
      if (i>0) fprintf(fptr, "\n");
  }
  fprintf(fptr, "%d\n", current_struct->natoms);
    
  
  /*Print a final carriage return*/
  fprintf(fptr,"\n");
  return (SUCCESS);
}

/**
 * Writes a file of which parameters are to be fit.
 * This can be read in later to fit multiple runs with one set of parameters to fit.
 * Assumes if multiple prmtops, the parameters are the same in each one
 * 
 * @param[in] global_options Global options structure, containing filename to save ase
 * @param[in] parm_data The parameter file containing which parameters are to be fit
 * @return Integer indicating success or failure
 */
int write_input_parameters(global_options_struct *global_options, parm_struct *parm_data)
{
  FILE *fptr;
  int i,j;

  /*Now fill the filename_to_write array*/
  printf("   Filename to be written is: %s\n",global_options->PARAMETER_FILE_NAME);
  
  if((fptr=fopen(global_options->PARAMETER_FILE_NAME,"w"))==NULL)
  {
    file_open_failure("create_input", global_options->PARAMETER_FILE_NAME);
    return FILE_OPEN_FAIL;
  }
  // print a header line indicating the data that go with this particular file
  fprintf(fptr,"#V2 AUTO-GENERATED PARAMETER FILE - DO NOT MODIFY\n");
  // put information about the desired bond parameters into the file
  fprintf(fptr, "#\n# BOND INFORMATION:\n");
  fprintf(fptr, "#### BOND\tREQ\tKR ####\n");
  for (i=0; i<parm_data->unique_bonds_found; i++)
    if (parm_data->bond_data[i].DO_FIT_REQ==YES || parm_data->bond_data[i].DO_FIT_KR==YES)
      fprintf(fptr, "\t%s %s\t%d\t%d\n", parm_data->bond_data[i].atom_type1, parm_data->bond_data[i].atom_type2,
              parm_data->bond_data[i].DO_FIT_REQ, parm_data->bond_data[i].DO_FIT_KR);
  
  // put information about the desired angle parameters into the file
  fprintf(fptr, "#\n# ANGLE PARAMETERS:\n");
  fprintf(fptr, "#### ANGLE\tKT\tTHEQ ####\n");
  for (i=0; i<parm_data->unique_angles_found; i++)
    if (parm_data->angle_data[i].DO_FIT_KT==YES || parm_data->angle_data[i].DO_FIT_THEQ==YES)
      fprintf(fptr, "\t%s %s %s\t%d\t%d\n", parm_data->angle_data[i].atom_type1, parm_data->angle_data[i].atom_type2,
              parm_data->angle_data[i].atom_type3, parm_data->angle_data[i].DO_FIT_KT, parm_data->angle_data[i].DO_FIT_THEQ);
  
  // put information about the desired dihedral parameters into the file
  fprintf(fptr, "#\n# DIHEDRAL PARAMETERS:\n");
  fprintf(fptr, "#### DIHEDRAL\tTERM\tKP\tNP\tPHASE ####\n");
  for (i=0; i<parm_data->unique_dihedrals_found; i++) {
    for (j=0; j<parm_data->dihedral_data[i].num_terms; ++j) {
      if (parm_data->dihedral_data[i].term[j].DO_FIT_KP==YES || parm_data->dihedral_data[i].term[j].DO_FIT_NP==YES ||
          parm_data->dihedral_data[i].term[j].DO_FIT_PHASE==YES)
        fprintf(fptr, "\t%s %s %s %s\t%d\t%d\t%d\t%d\n", parm_data->dihedral_data[i].atom_type1, parm_data->dihedral_data[i].atom_type2,
                parm_data->dihedral_data[i].atom_type3, parm_data->dihedral_data[i].atom_type4, j, parm_data->dihedral_data[i].term[j].DO_FIT_KP, 
                parm_data->dihedral_data[i].term[j].DO_FIT_NP, parm_data->dihedral_data[i].term[j].DO_FIT_PHASE);
  } }
	    
  return SUCCESS;
}

/**
 * Write a force field modification file.
 * To be used to save the results of a fit for reading into Leap
 * 
 * @param[in] global_options The global options structure, containing filename to write
 * @param[in] parm_data The parameter file containing fitted results
 * @return Integer indicating success or failure
 */
int write_frcmod(global_options_struct *global_options, parm_struct *parm_data)
{
  FILE *fptr;
  int i,j;
  
  /* Make sure the settings are correct */
  if((fptr=fopen(global_options->WRITE_FRCMOD,"w"))==NULL)
  {
    file_open_failure("write_ffrcmod", global_options->WRITE_FRCMOD);
    return FILE_OPEN_FAIL;
  }
  if (global_options->VERBOSITY>=MEDIUM)
    printf(" * Saving ffrcmod file to %s\n", global_options->WRITE_FRCMOD);
  
  // print a header line indicating the data that go with this particular file
  fprintf(fptr, "Generated frcmod with paramfit for %s\n", parm_data->filename);
  
  // print the bond information, if any
  if (global_options->BOND_PARAMS > 0)
    fprintf(fptr,"BOND\n");
  for (i=0; i<parm_data->unique_bonds_found; ++i)
  {
    if (parm_data->bond_data[i].DO_FIT_KR==YES || parm_data->bond_data[i].DO_FIT_REQ==YES)
    {
      fprintf(fptr, "%s-%s %14.4f %14.4f\n", parm_data->bond_data[i].atom_type1, parm_data->bond_data[i].atom_type2,
                                               parm_data->bond_data[i].rk, parm_data->bond_data[i].req);
    }
  }
  fprintf(fptr, "\n");
  
  // now deal with the angles
  if (global_options->ANGLE_PARAMS > 0)
    fprintf(fptr, "ANGL\n");
  for (i=0; i<parm_data->unique_angles_found; ++i)
  {
    if (parm_data->angle_data[i].DO_FIT_KT==YES || parm_data->angle_data[i].DO_FIT_THEQ==YES)
    {
      fprintf(fptr, "%s-%s-%s %14.4f %14.4f\n", parm_data->angle_data[i].atom_type1, parm_data->angle_data[i].atom_type2,
                                                      parm_data->angle_data[i].atom_type3, parm_data->angle_data[i].tk,
                                                      parm_data->angle_data[i].teq*RADIAN_TO_DEGREE);
    }
  }
  fprintf(fptr, "\n");
  
 // Now the dihedrals
  int num_fit;
  if (global_options->DIHEDRAL_PARAMS > 0)
   fprintf(fptr, "DIHE\n");
  int impropers = 0;
  for (i=0; i<parm_data->unique_dihedrals_found; ++i) {
    // Need to count # of fitted dihedrals for this term so negative sign on periodicity is correct
    num_fit=0;
    for (j=0; j<parm_data->dihedral_data[i].num_terms; ++j) {
      if (parm_data->dihedral_data[i].term[j].DO_FIT_KP==YES || parm_data->dihedral_data[i].term[j].DO_FIT_NP==YES ||
          parm_data->dihedral_data[i].term[j].DO_FIT_PHASE==YES)
        ++num_fit;
    }
    for (j=0; j<parm_data->dihedral_data[i].num_terms; ++j) {
      if (parm_data->dihedral_data[i].term[j].DO_FIT_KP==YES || parm_data->dihedral_data[i].term[j].DO_FIT_NP==YES ||
          parm_data->dihedral_data[i].term[j].DO_FIT_PHASE==YES) {
        if (parm_data->dihedral_data[i].improper==YES) {
          ++impropers;
          break;
        } else {
          // pn negative to indicate more terms follow for this dihedral
          double pn = parm_data->dihedral_data[i].term[j].pn;
          if (num_fit!=1) pn *= -1;
          fprintf(fptr, "%s-%s-%s-%s    1 %14.4f %14.4f %14.4f\n", parm_data->dihedral_data[i].atom_type1,
              parm_data->dihedral_data[i].atom_type2, parm_data->dihedral_data[i].atom_type3,
              parm_data->dihedral_data[i].atom_type4, parm_data->dihedral_data[i].term[j].pk,
              parm_data->dihedral_data[i].term[j].phase*RADIAN_TO_DEGREE, pn);
          --num_fit;
        }
      } 
    }
  }
  fprintf(fptr, "\n");

  if (impropers > 0)
    fprintf(fptr, "IMPR\n");
  for (i=0; i<parm_data->unique_dihedrals_found; ++i) {
    // Skip if not improper dihedral
    if (parm_data->dihedral_data[i].improper==NO) continue;
    
    // Need to count # of fitted dihedrals for this term so negative sign on periodicity is correct
    num_fit=0;
    for (j=0; j<parm_data->dihedral_data[i].num_terms; ++j) {
      if (parm_data->dihedral_data[i].term[j].DO_FIT_KP==YES || parm_data->dihedral_data[i].term[j].DO_FIT_NP==YES ||
          parm_data->dihedral_data[i].term[j].DO_FIT_PHASE==YES)
        ++num_fit;
    }
    for (j=0; j<parm_data->dihedral_data[i].num_terms; ++j) {
      if (parm_data->dihedral_data[i].term[j].DO_FIT_KP==YES || parm_data->dihedral_data[i].term[j].DO_FIT_NP==YES ||
          parm_data->dihedral_data[i].term[j].DO_FIT_PHASE==YES)
      {
        // pn has - sign to indicate more terms follow for this dihedral 
        double pn = parm_data->dihedral_data[i].term[j].pn;
        if (num_fit!=1) pn *= -1;  
        fprintf(fptr, "%s-%s-%s-%s     %14.4f %14.4f %14.4f\n", parm_data->dihedral_data[i].atom_type1,
                parm_data->dihedral_data[i].atom_type2, parm_data->dihedral_data[i].atom_type3,
                parm_data->dihedral_data[i].atom_type4, parm_data->dihedral_data[i].term[j].pk,
                parm_data->dihedral_data[i].term[j].phase*RADIAN_TO_DEGREE, pn);
        --num_fit;
      }
    }
  }
  fprintf(fptr, "\n");

 return SUCCESS;
}

/**
 * Writes an input file for ADF for a structure.
 * The file is initialized in write quantum.
 * Note that this file is only for energy calculations, not forces.
 * 
 * @param[in] global_options The global options structure
 * @param[in] parm_data Pointer to the parameter structure that describes these coordinates
 * @param[in] current_struct Pointer to the coordinate set that will be used
 * @param[in] num Number of the structure in teh coordinate set to write
 * @param[in,out] fptr Pointer to the file to write to, initialized elsewhere
 * @return Integer indicating success or failure
 */
int write_input_adf(global_options_struct *global_options, parm_struct *parm_data, coord_set *current_struct, int num, FILE *fptr)
{
  // write some basic header information
  fprintf(fptr, "TITLE structure %d generated with paramfit\n", num);
  if (global_options->QMHEADER)
    fprintf(fptr, "INLINE %s\n!\n", global_options->QMHEADER);
  fprintf(fptr, "UNITS\n length Angstrom\n angle Radian\nEnd\n!\n");
  
  // write the charge and multiplicity
  fprintf(fptr, "CHARGE %d %d\n!\n", global_options->QM_SYSTEM_CHARGE, global_options->QM_SYSTEM_MULTIPLICITY);
  fprintf(fptr, "TOTALENERGY\n!\n");
  
  // write the cartesian matrix
  fprintf(fptr, "ATOMS Cartesian\n");
  int i, atomic_number;
  for (i=0;i<current_struct->natoms;++i)
  {
    atomic_number=find_atomic_number_from_parm(parm_data, i+1);
    if (atomic_number==UNKNOWN_ELEMENT)
    {
      /*Element not recognised*/
      printf("\n*** ERROR IN WRITE_INPUT FROM CALL TO FIND_ATOMIC_NUMBER_FROM_PARM\n");
      printf("*** UNKNOWN ELEMENT FOR ATOM NAME: %s WITH MASS: %f\n",parm_data->atom[i].igraph, parm_data->atom[i].amass);
      printf("*** CHECK elements.c\n");
      return (UNKNOWN_ELEMENT);
    }
    fprintf(fptr, " %i ", i);
    print_atomic_number_as_symbol(fptr, atomic_number);
    fprintf(fptr, " %.6f %.6f %.6f\n",
	    current_struct->struc[num].x_coord[i], current_struct->struc[num].y_coord[i], current_struct->struc[num].z_coord[i]);
  }
  fprintf(fptr, "End\n");
  return SUCCESS;
}

/**
 * Writes a quantum input file for GAMESS.
 * This is for energy calculation only.
 * 
 * @param[in] global_options The global options structure
 * @param[in] parm_data The parameter structure
 * @param[in] current_structure Pointer to coordinate set to write from
 * @param[in] num The number of the structure in the coordinate set that will be written
 * @param[in,out] fptr File to write to, should be already initialized
 * @return Integer indicating success or failure.
 */
int write_input_gamess(global_options_struct *global_options, parm_struct *parm_data, coord_set *current_struct, int num, FILE *fptr)
{
  // include the user's header, which will probably be the $SYSTEM section and the $BASIS section
  time_t now = time(NULL);
  fprintf(fptr, "! Gamess input file for structure %d automatically generated with paramfit\n", num);
  fprintf(fptr, "! Generated on %s for prmtop at %s", asctime(localtime(&now)), parm_data->filename);
  if (global_options->QMHEADER)
  {
    FILE *header = fopen(global_options->QMHEADER, "r");
    if (!header)
    {
      printf("! ERROR opening QM header file at %s\n", global_options->QMHEADER);
      return INVALID_DATA;
    }
    // Get the size of the header
    fseek(header, 0, SEEK_END);
    long size = ftell(header);
    rewind(header);
    
    // allocate buffer, read in data, and write it to the file
    unsigned char *buffer=(unsigned char *)malloc((size_t)size);
    fread(buffer, size, 1, header);
    fwrite(buffer, size, 1, fptr);
    
    // clean up
    free(buffer);
    fclose(header);
  }
  
  // write the basic $CONTRL section
  fprintf(fptr, " $CONTRL\n");
  fprintf(fptr, "  RUNTYP=ENERGY\n");
  fprintf(fptr, "  COORD=UNIQUE\n");
  fprintf(fptr, "  ICHARG=%d\n  MULT=%d\n", global_options->QM_SYSTEM_CHARGE, global_options->QM_SYSTEM_MULTIPLICITY);
  fprintf(fptr, " $END\n!\n");
  
  // write the $DATA section
  fprintf(fptr, " $DATA\n");
  fprintf(fptr, "  Structure %d generated by paramfit\n", num);
  fprintf(fptr, "  C1\n"); // don't specify symmetry since molecule could be anything
  
  // write all the atoms
  int i, atomic_number;
  for (i=0;i<current_struct->natoms;++i)
  {
    atomic_number=find_atomic_number_from_parm(parm_data, i+1);
    if (atomic_number==UNKNOWN_ELEMENT)
    {
      /*Element not recognised*/
      printf("\n*** ERROR IN WRITE_INPUT FROM CALL TO FIND_ATOMIC_NUMBER_FROM_PARM\n");
      printf("*** UNKNOWN ELEMENT FOR ATOM NAME: %s WITH MASS: %f\n",parm_data->atom[i].igraph, parm_data->atom[i].amass);
      printf("*** CHECK elements.c\n");
      return (UNKNOWN_ELEMENT);
    }
    fprintf(fptr, "  ");
    print_atomic_number_as_symbol(fptr, atomic_number);
    fprintf(fptr, " %d %10.6f %10.6f %10.6f\n", atomic_number,
	    current_struct->struc[num].x_coord[i], current_struct->struc[num].y_coord[i], current_struct->struc[num].z_coord[i]);
  }
  fprintf(fptr, " $END\n");
  return SUCCESS;
}

/**
 * Writes a data file with amber and quantum energies for each structure with given parameters.
 * Useful to plot the results of a calculation to see the quality of fit.
 * 
 * @param[in] global_options The global options structure, with filename to save to
 * @param[in] parm_data  Pointer to single parameter set to use when calculating amber energies
 * @param[in] coords_data Pointer to single coordinate set with attached quantum energy
 * @param[in] generation The generation of the genetic algorithm, so can be called iteratively
 * @return Integer indicating success or failure
 */
int write_energy(global_options_struct *global_options, parm_struct *parm_data, coord_set *coords_data)
{
  FILE *fptr;
  int i;
  double temp_energy;
  char filename[1000];
  
  // create the file
  if (global_options->num_prmtops > 1)
    sprintf(filename, "%s_%s", parm_data->filename, global_options->WRITE_ENERGY);
  else
    sprintf(filename, "%s", global_options->WRITE_ENERGY);
  
  if((fptr=fopen(filename,"w"))==NULL)
  {
    file_open_failure("write_energies", global_options->WRITE_ENERGY);
    return FILE_OPEN_FAIL;
  }
  if (global_options->VERBOSITY>=MEDIUM)
    printf(" * Saving energy file with %d structures to %s\n", coords_data->num_coords, filename);

  // print a header line indicating the data that go with this particular file
  fprintf(fptr, "# Final energies generated by paramfit for %s\n", parm_data->filename);
  fprintf(fptr, "# Num\tAmber+K\t\tQuantum\t\tInitial Amber+K\n");
  
  // evaluate the energy for all of the input structures and print it to the file
  // here we put K back in since it represents the difference between quantum and classical so everything lines up
  for (i=0; i<coords_data->num_coords; ++i) {
    temp_energy = eval_amber_std_for_single_struct(global_options, parm_data, &coords_data->struc[i]);
    fprintf(fptr, "%d\t%10.5f\t%10.5f\t%10.5f\n", i, (temp_energy+parm_data->K), coords_data->struc[i].energy, (coords_data->struc[i].init_energy+parm_data->K));
  }
  
  return SUCCESS;
}
