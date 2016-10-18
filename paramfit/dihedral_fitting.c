/*****************************************************
 * AMBER Bond Angle and Dihedral Parameter Optimiser *
 *                                                   *
 *           Written by: Robin Betz  (2013)          *
 *                       Ross Walker (2004)          *
 *                   UC San Diego                    *
 *           San Diego Supercomputer Center          *
 *            La Jolla, California, 92092            *
 *                       USA                         *
 *****************************************************/

/** @file dihedral_fitting.c
 * Implements dihedral fitting as a linear minimization
 * Uses the method described by Chad Hopkins and Adrian Roitberg to represent the dihedral
 * fitting problem as linear minimization. Sets up the problem, constructs the linear problem,
 * solves it with Cholesky decomposition and back-substitution, and returns new dihedral data
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "function_def.h"
#include "constants.h"

/**
 * Wrapper function for all dihedral fitting functions.
 *
 * Calls all the methods necessary to spit out new parameters so you don't have to remember which ones
 * to use.
 *
 * @see construct_initial_dihedral_params, conduct_dihedral_least_squares
 *
 * @param[in] global_options The global options structure
 * @param[in,out] parm_data Array containing the dihedral and parameter information at the beginning, and will be updated with the new parameters at the end.
 * @param[in] coords_data Pointer to single structure containing coordinates and QM energy of each input structure
 * @param[in] fit_phase Whether or not phases will be fit or just dihedral force constants
 * @return Integer indicating success or failure. Parm_data is updated with the new parameters if successful
 */
int dihedral_least_squares(global_options_struct *global_options, parm_struct *parm_datas, coord_set *coords_datas, bool_t fit_phase)
{
  conduct_dihedral_least_squares(global_options, parm_datas, coords_datas, fit_phase);
  return SUCCESS;
}


/**
 * Constructs and solves the linear system for the dihedral fitting algorithm.
 *
 * Uses the transformations described by Hopkins and Roitberg to present dihedral fitting
 * as a linear system and solves it using a Cholesky decomposition and back-substitution,
 * then undoes the transformations to get back to normal dihedral space.
 *
 * @see dihedral_least_squares
 *
 * @param[in] global_options The global options structure
 * @param[in,out] parm_data Single parameter structure, will be update
 * @param[in] coords_data Single input coordinate set with associated QM energies
 * @param[in] fit_phase Whether to fit phases as well or just dihedral force constants
 * @return Integer indicating success or failure
 */
int conduct_dihedral_least_squares(global_options_struct *global_options, parm_struct *parm_data, coord_set *coords_data,
                                       bool_t fit_phase)
{
  if (global_options->num_prmtops != 1) {
    printf("!  ERROR: Multiple molecule fits not supported for dihedral least squares fit!");
    return INVALID_DATA;
  }

  double *k;         // for matrix construction
  double *T;         // for mm energies
  double **C;        // matrix for the linear system to solve
  int dihe, structure;
  int i,j,n;

  int phaseoffset=parm_data->unique_dihedral_terms;
  int size = parm_data->unique_dihedral_terms;
  if (fit_phase == YES) size *= 2;

  C = alloc_2D_double(size, size);
  if (C==NULL) malloc_failure_double("dihedral_least_squares", "C", size*size);

  // Calculate the difference between qm and mm energy of each structure
  T = (double *)malloc(coords_data->num_coords*sizeof(double));
  if (T==NULL) malloc_failure_double("dihedral_least_squares", "mm", coords_data->num_coords);
  for (j=0; j<coords_data->num_coords; ++j) {
    T[j] = coords_data->struc[j].energy - eval_amber_std_for_single_struct( global_options, parm_data, &(coords_data->struc[j]) );
  }

  // Construct F
  double tempx1, tempy1, tempz1, tempx2, tempy2, tempz2;
  double tempx3, tempy3, tempz3, tempx4, tempy4, tempz4;
  double **F = alloc_2D_double(size, coords_data->num_coords);

  for (structure=0; structure<coords_data->num_coords; ++structure) {
    dihe=0;
    for (i=0; i<parm_data->unique_dihedrals_found; ++i) {
      F[dihe][structure]=0.;
      if (fit_phase==YES) F[dihe+i+phaseoffset][structure]=0;
        for (n=0; n<parm_data->dihedral_data[i].number; ++n) {
          tempx1=coords_data->struc[structure].x_coord[parm_data->dihedral_data[i].atom1[n]-1];
          tempy1=coords_data->struc[structure].y_coord[parm_data->dihedral_data[i].atom1[n]-1];
          tempz1=coords_data->struc[structure].z_coord[parm_data->dihedral_data[i].atom1[n]-1];
          tempx2=coords_data->struc[structure].x_coord[parm_data->dihedral_data[i].atom2[n]-1];
          tempy2=coords_data->struc[structure].y_coord[parm_data->dihedral_data[i].atom2[n]-1];
          tempz2=coords_data->struc[structure].z_coord[parm_data->dihedral_data[i].atom2[n]-1];
          tempx3=coords_data->struc[structure].x_coord[parm_data->dihedral_data[i].atom3[n]-1];
          tempy3=coords_data->struc[structure].y_coord[parm_data->dihedral_data[i].atom3[n]-1];
          tempz3=coords_data->struc[structure].z_coord[parm_data->dihedral_data[i].atom3[n]-1];
          tempx4=coords_data->struc[structure].x_coord[parm_data->dihedral_data[i].atom4[n]-1];
          tempy4=coords_data->struc[structure].y_coord[parm_data->dihedral_data[i].atom4[n]-1];
          tempz4=coords_data->struc[structure].z_coord[parm_data->dihedral_data[i].atom4[n]-1];
          double angle =  PI-calc_dihedral_radians(tempx1,tempy1,tempz1,tempx2,tempy2,tempz2,tempx3,tempy3,tempz3,tempx4,tempy4,tempz4);

          for (j=0; j<parm_data->dihedral_data[i].num_terms; ++j) {
            F[dihe+j][structure] += cos(parm_data->dihedral_data[i].term[j].pn * angle);
            if (fit_phase==YES) F[dihe+phaseoffset][structure] += sin(parm_data->dihedral_data[i].term[j].pn * angle );
          }
        }
      dihe+=parm_data->dihedral_data[i].num_terms;
    }
  }

  // Construct K
  k = (double *)malloc(size*sizeof(double));
  if (k==NULL) malloc_failure_double("dihedral_least_squares", "k", size);
  for (i=0; i<size; ++i) {
    k[i] = 0.;
    n=0;
    for (j=0; j<coords_data->num_coords; ++j) {
        k[i] += F[i][j]*T[j];
    }
  }
  // We no longer need mm array so free it
  free(T); T = NULL;

  // construct the linear matrix C = F*F
  for (i=0; i<size; ++i) {
    for (j=0; j<size; ++j) {
      C[i][j] = 0.;
      for (n=0; n<global_options->TOTAL_STRUCTURES; ++n) {
        C[i][j] += F[i][n]*F[j][n];
      }
    }
  }
  // We no longer need F so free it
  double_2D_array_free(F);

  // Decompose the symmetric matrix C into C=LL* using the Cholesky-Banachiewicz algorithm
  // DEBUG: this is correct
  double **L = alloc_2D_double(size, size);
  if (L==NULL) malloc_failure_double("dihedral_least_squares", "L", size*size);
  for (i=0; i<size; ++i) {
    for (j=0; j<size; ++j) {
      if (i<j) {
        L[i][j] = 0; // lower triangular so upper triangle is 0
      }
      else if (i==j) { // on the diagonal
        L[i][j] = C[i][j];
        for (n=0; n<j; ++n) {
          L[i][j] -= (L[j][n]*L[j][n]);
        }
        if (L[i][j] <= 0.) {
          printf("\n\n! ERROR: Singular matrix in dihedral_least_squares\n");
          printf("  Add more input structures\n" );
          exit(INVALID_DATA);
        }
        L[i][j] = sqrt(L[i][j]);
      }
      else { // below the diagonal
        L[i][j] = C[i][j];
        for (n=0; n<j; ++n)
          L[i][j] -= L[i][n]*L[j][n];
        L[i][j] /= L[j][j];
      }
    }
  }

  // Solve Ly=k by forward substitution
  // DEBUG: this is correct
  double *y = (double*)malloc(size*sizeof(double));
  if (y==NULL) malloc_failure_double("dihedral_least_squares","y", size);
  for (i=0; i<size; ++i) {
    y[i] = k[i];
    for (j=0; j<i; ++j) {
      y[i] -= L[i][j]*y[j];
    }
    y[i] /= L[i][i];
  }

  // Now solve L*a=y by backward substitution
  // DEBUG: this is correct
  double *delta = (double*)calloc(2*parm_data->unique_dihedral_terms,sizeof(double));
  if (delta==NULL) malloc_failure_double("dihedral_least_squares","delta", 2*parm_data->unique_dihedral_terms);
  for (i=size-1; i>=0; --i) {
    delta[i] = y[i];
    for (j=i+1; j<size; ++j) {
      delta[i] -= L[j][i]*delta[j];
    }
    delta[i] /= L[i][i];
  }

  // Calculate the final result
  // Do it exactly as pseudocode for debugging purposes even though it is a little less
  // time and space efficient this way
  // DEBUG: this is correct
  double *ini_alt=(double*)malloc(2*parm_data->unique_dihedral_terms*sizeof(double));
  double *a_alt=(double*)malloc(2*parm_data->unique_dihedral_terms*sizeof(double));
  dihe=0;
  for (i=0; i<parm_data->unique_dihedrals_found; ++i) {
    for (j=0; j<parm_data->dihedral_data[i].num_terms; ++j) {
      ini_alt[dihe] = parm_data->dihedral_data[i].term[j].pk*cos(parm_data->dihedral_data[i].term[j].phase);
      ini_alt[dihe+phaseoffset] = parm_data->dihedral_data[i].term[j].pk*sin(parm_data->dihedral_data[i].term[j].phase);
      ++dihe;
    }
  }

  for (i=0; i<parm_data->unique_dihedral_terms*2; ++i) {
    a_alt[i] = ini_alt[i] + delta[i];
  }
  double *A = (double*)calloc(2*parm_data->unique_dihedral_terms, sizeof(double));
  dihe=0;
  for (i=0; i<parm_data->unique_dihedrals_found; ++i) {
    for (j=0; j<parm_data->dihedral_data[i].num_terms; ++j) {
      A[dihe] = sqrt(a_alt[dihe]*a_alt[dihe] + a_alt[dihe+phaseoffset]*a_alt[dihe+phaseoffset]);

      if (parm_data->dihedral_data[i].term[j].pk == 0.) // atan will return NaN but this is defined as 0. in this context for the phase
        A[dihe+phaseoffset] = 0.0;
      else
        A[dihe+phaseoffset] = atan2(a_alt[dihe+phaseoffset], a_alt[dihe]);

      parm_data->dihedral_data[i].term[j].pk = A[dihe];
      parm_data->dihedral_data[i].term[j].phase = A[dihe+phaseoffset];
      parm_data->dihedral_data[i].term[j].phase = fmod(parm_data->dihedral_data[i].term[j].phase, 2.0*PI);
      ++dihe;
    }
  }

  // clean up and exit
  free(delta);
  free(ini_alt); free(a_alt);
  free(A);

  free(k);
  free(y);
  double_2D_array_free(C);
  double_2D_array_free(L);

  return SUCCESS;
}
