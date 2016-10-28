/*****************************************************
 * AMBER Bond Angle and Dihedral Parameter Optimiser *
 *                                                   *
 *           Written by: Robin Betz  (2012)          *
 *                   UC San Diego                    *
 *           San Diego Supercomputer Center          *
 *            La Jolla, California, 92092            *
 *                       USA                         *
 *****************************************************/

/* forces.h */
#ifndef pFORCES_H
#define pFORCES_H

/* Force structure is just a vector */
typedef struct _force_struct {
  double x, y, z;
} force_struct;

#endif

