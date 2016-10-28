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
/** @file param_summary.c
 * Prints out a summary of all the parameters
 */

#include <stdio.h>
#include "function_def.h"
#include "constants.h"

/**
 * Prints out a summary of only the parameters to be fit over a set of multiple prmtops.
 * Assumes that each prmtop contains every single parameter to be fit, because only the data
 * from the first one is used.
 * @see verify_prmtops
 * @param[in] global_options The global options structure
 * @param[in] parm_datas Array of parameter structures
 */
void print_multiprmtop_summary(global_options_struct *global_options, parm_struct *parm_datas)
{
  int param,term;
  // Print out K for each prmtop
  printf("       Prmtop                       K \n");
  printf("      ------------------          --------\n");
  for (param=0; param<global_options->num_prmtops; ++param) {
    printf("      %18s %18.6f kcal/mol\n", parm_datas[param].filename, parm_datas[param].K);
  }
  printf("\n");
  // Print out bond terms
  for (param=0; param<parm_datas[0].unique_bonds_found; ++param) {
    if (parm_datas[0].bond_data[param].DO_FIT_KR==YES || parm_datas[0].bond_data[param].DO_FIT_REQ==YES) {
      printf("            *(%s-%s) ",parm_datas[0].bond_data[param].atom_type1,parm_datas[0].bond_data[param].atom_type2);
      printf("Kr = %8.4f kcal/(mol A)^2,",parm_datas[0].bond_data[param].rk);
      printf("r_eq = %8.4f A \n",parm_datas[0].bond_data[param].req);
  } }
  // Angle terms
  for (param=0; param<parm_datas[0].unique_angles_found; ++param) {
    if (parm_datas[0].angle_data[param].DO_FIT_KT==YES || parm_datas[0].angle_data[param].DO_FIT_THEQ==YES) {
      printf("       *(%s-%s-%s) ",parm_datas[0].angle_data[param].atom_type1,parm_datas[0].angle_data[param].atom_type2,
              parm_datas[0].angle_data[param].atom_type3);
      printf("Kt = %8.4f kcal/(mol rad)^2, ",parm_datas[0].angle_data[param].tk);
      printf("th_eq = %8.4f deg \n", parm_datas[0].angle_data[param].teq*RADIAN_TO_DEGREE);
  } }
  // Dihedral terms
  for (param=0; param<parm_datas[0].unique_dihedrals_found; ++param) {
    for (term=0; term<parm_datas[0].dihedral_data[param].num_terms;++term) {
      if (parm_datas[0].dihedral_data[param].term[term].DO_FIT_KP==YES ||
          parm_datas[0].dihedral_data[param].term[term].DO_FIT_NP==YES ||
          parm_datas[0].dihedral_data[param].term[term].DO_FIT_PHASE==YES) {
        if (parm_datas[0].dihedral_data[param].improper==YES) printf("  IMP ");
        else printf("      ");
        printf("*(%s-%s-%s-%s) ",parm_datas[0].dihedral_data[param].atom_type1,parm_datas[0].dihedral_data[param].atom_type2,
               parm_datas[0].dihedral_data[param].atom_type3,parm_datas[0].dihedral_data[param].atom_type4);
        printf("Kp = %8.4f kcal/mol, ",parm_datas[0].dihedral_data[param].term[term].pk);
        printf("Np = %6.4f, ",parm_datas[0].dihedral_data[param].term[term].pn);
        printf("Phase = %8.4f Deg \n",parm_datas[0].dihedral_data[param].term[term].phase*RADIAN_TO_DEGREE);    
  } } }
}

/**
 * Prints out a summary of all parameters currently stored in one parm_data.
 * A * will mark which parameters are being fit.
 * @param global_options The global options structure
 * @param parm_data The parameter structure to print
 */
void print_parameter_summary(global_options_struct *global_options, parm_struct *parm_data)
{
  int param, term;
  
  if (global_options->FUNC_TO_FIT==SUM_SQUARES_AMBER_STANDARD || 
      global_options->FUNC_TO_FIT==AMBER_FORCES ||
      global_options->FUNC_TO_FIT==DIHEDRAL_LEAST_SQUARES )
  {
    printf("   Parameters for force field equation: AMBER_STANDARD:\n");
    if (global_options->VERBOSITY >= HIGH)
      printf("   (* means parameter is NOT constant during fit)\n");
    else
      printf("   (* means parameter is NOT constant during fit)\n");

    /*First off, if K was fitted print the value of K*/
    if (global_options->K_FIT==YES)
       printf("                         *K = %8.6f kcal/mol\n",parm_data->K);
    else
      printf("                         K = %8.6f kcal/mol\n", parm_data->K);

    /*next, print the bond parameters*/
    for (param=0;param<parm_data->unique_bonds_found;++param) {
      printf("             (%s-%s)",parm_data->bond_data[param].atom_type1,parm_data->bond_data[param].atom_type2);
      if (parm_data->bond_data[param].DO_FIT_KR==YES) printf("*");
      else printf(" ");
      printf("Kr = %8.4f kcal/(mol A)^2,",parm_data->bond_data[param].rk);
      
      if (parm_data->bond_data[param].DO_FIT_REQ==YES) printf("*");
      else printf(" ");
      printf("r_eq = %8.4f A \n",parm_data->bond_data[param].req);
    }
     /*next, print the angle parameters*/
     for (param=0;param<parm_data->unique_angles_found;++param)
     {
        printf("        (%s-%s-%s)",parm_data->angle_data[param].atom_type1,parm_data->angle_data[param].atom_type2,
               parm_data->angle_data[param].atom_type3);
        if (parm_data->angle_data[param].DO_FIT_KT==YES) printf("*");
        else printf(" ");
        printf("Kt = %8.4f kcal/(mol rad)^2, ",parm_data->angle_data[param].tk);
        
        if (parm_data->angle_data[param].DO_FIT_THEQ==YES) printf("*");
        else printf(" ");
        printf("th_eq = %8.4f deg \n", parm_data->angle_data[param].teq*RADIAN_TO_DEGREE);
      }
     /*finally, print the dihedral parameters*/
     for (param=0;param<parm_data->unique_dihedrals_found;++param)
     {
       for (term=0; term<parm_data->dihedral_data[param].num_terms; ++term) {
          if (parm_data->dihedral_data[param].improper==YES) printf("   IMP ");
          else printf("       ");
          printf("(%s-%s-%s-%s)",parm_data->dihedral_data[param].atom_type1,parm_data->dihedral_data[param].atom_type2,
                 parm_data->dihedral_data[param].atom_type3,parm_data->dihedral_data[param].atom_type4);
          
          if (parm_data->dihedral_data[param].term[term].DO_FIT_KP==YES) printf("*");
          else printf(" ");
          printf("Kp = %8.4f kcal/mol, ",parm_data->dihedral_data[param].term[term].pk);
          
          if (parm_data->dihedral_data[param].term[term].DO_FIT_NP==YES) printf("*");
          else printf(" ");
          printf("Np = %6.4f, ",parm_data->dihedral_data[param].term[term].pn);
          
          if (parm_data->dihedral_data[param].term[term].DO_FIT_PHASE==YES) printf("*");
          else printf(" ");
          printf("Phase = %8.4f Deg \n",parm_data->dihedral_data[param].term[term].phase*RADIAN_TO_DEGREE);    
       }
     }
   }
  else
  {
     printf("  FORCE FIELD EQUATION IS NOT CURRENTLY IMPLEMENTED.\n");
  }
}

