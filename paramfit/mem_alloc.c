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

/*mem_alloc.c*/
/*Memory allocation routines*/

#include <stdlib.h>
#include <stdio.h>

#include "function_def.h"

double **alloc_2D_double(int nrows, int ncolumns)
{
  /*Allocates a 2d_double_array consisting of a series of pointers pointing to each
    row that are then allocated to be ncolumns long each.*/

  /*Uses calloc - slower but all locations will be zeroed*/
  
  /*Tries to keep contents contiguous - thus reallocation is difficult!*/

  /*Returns the pointer **array. Returns NULL on error*/
  int i;
  
  double **array = (double **)calloc(nrows,sizeof(double *));
  if (array==NULL)
    return NULL;
  array[0] = (double *)calloc(nrows * ncolumns,sizeof(double));
  if (array[0]==NULL)
     return NULL;
  
  for (i = 1; i < nrows; ++i)
    array[i] = array[0] + i * ncolumns;

  return array;
}

void double_2D_array_free(double **array)
{
  /*Frees the memory previously allocated by alloc_2D_double*/
  free(array[0]);
  free(array);
}

void global_unlock(global_options_struct *global_options, parm_struct **parm_datas, coord_set **coords_data)
{
  int i;
  free(global_options->mdcrd_list);    
  free(global_options->prmtop_list);
  
  if (global_options->job_control_filename) {
    free(global_options->job_control_filename); 
    global_options->job_control_filename = NULL;
  }
  if (global_options->WRITE_FRCMOD)
  {
    free(global_options->WRITE_FRCMOD);
    global_options->WRITE_FRCMOD=NULL;
  }
  if (global_options->QMHEADER)
  {
    free(global_options->QMHEADER);
    global_options->QMHEADER=NULL;
  }
  if (global_options->WRITE_ENERGY)
  {
    free(global_options->WRITE_ENERGY);
    global_options->WRITE_ENERGY=NULL;
  }
  if (global_options->PARAMETER_FILE_NAME)
  {
    free(global_options->PARAMETER_FILE_NAME);
    global_options->PARAMETER_FILE_NAME = NULL;
  }
  
  if (global_options->QMFILEOUTSTART) {
    free(global_options->QMFILEOUTSTART);
    global_options->QMFILEOUTSTART = NULL;
  }
  if (global_options->QMFILEOUTEND) {
    free(global_options->QMFILEOUTEND);
    global_options->QMFILEOUTEND = NULL;
  }

  for (i=0; i<global_options->num_prmtops; ++i) {
    free_prmtop(global_options, &(*parm_datas)[i]);
    if (coords_data)
      free_coords(global_options, &(*coords_data)[i]);
  }
  free(*parm_datas);
  free(*coords_data); 
}

void free_prmtop(global_options_struct *global_options, parm_struct *parm_data) 
{
  int i;
  free(parm_data->filename);
  parm_data->filename=NULL;
  if (global_options->RUNTYPE==FIT) {
    free(parm_data->bond_data);
    free(parm_data->angle_data);
    for (i=0; i<parm_data->unique_dihedrals_found; ++i)
      free(parm_data->dihedral_data[i].term);
    free(parm_data->dihedral_data);
  
    if (global_options->FUNC_TO_FIT==AMBER_FORCES)
      free(parm_data->fit_atom);
  }
  if ( parm_data->NHB > 0 ) {
      free(parm_data->bg);  
      free(parm_data->ag); 
      
      parm_data->bg = NULL;
      parm_data->ag = NULL;
  }

  if ( parm_data->NEXT > 0 ) {
    free(parm_data->natex); 
    parm_data->natex = NULL;
  }
  if (parm_data->MPHIA > 0) {
    free(parm_data->pdihedral); 
    parm_data->pdihedral = NULL;
  }
  if (parm_data->NPHIH > 0) {
    free(parm_data->pdihedralH); 
    parm_data->pdihedralH = NULL;
  }

  if (parm_data->MTHETS > 0) {
    free(parm_data->pangle);
    parm_data->pangle = NULL;
  }
  if (parm_data->NTHETH > 0 ) {
    free(parm_data->pangleH);
    parm_data->pangleH = NULL;
  }
  if (parm_data->NBONA > 0 ) {
    free(parm_data->pbond);
    parm_data->pbond = NULL;
  }
  if (parm_data->NBONH > 0) {
    free(parm_data->pbondH);
    parm_data->pbondH = NULL;
  }

  if (parm_data->NTYPES > 0)
  {
    free(parm_data->cn2);
    free(parm_data->cn1);
    
    parm_data->cn2 = NULL;
    parm_data->cn1 = NULL;
  }

  if (parm_data->NATYP > 0)
  {
    free(parm_data->solty);   
    parm_data->solty = NULL;
  }
  if (parm_data->MPTRA > 0)
  {
    free(parm_data->phase);  
    free(parm_data->pn);    
    free(parm_data->pk);
  
    parm_data->phase = NULL;
    parm_data->pn = NULL;
    parm_data->pk = NULL;
  }
  if (parm_data->MUMANG)
  {
    free(parm_data->teq);
    free(parm_data->tk);
  
    parm_data->teq = NULL;
    parm_data->tk = NULL;
  }
  if (parm_data->MUMBND)
  {
    free(parm_data->req);
    free(parm_data->rk);
  
    parm_data->req = NULL;
    parm_data->rk = NULL;
  }
  if (parm_data->NTOTRS > 0)
  {
    free(parm_data->residue);
    parm_data->residue = NULL;
  }  
  if (parm_data->NTYPES > 0)
  {
      free(parm_data->nno);
      parm_data->nno = NULL;
  }
  if ( parm_data->NTOTAT > 0)
  {
      free(parm_data->atom);   
      parm_data->atom = NULL;
  }
  
  free(parm_data->title);
  parm_data->title = NULL;
    
}

/**
 * Allocates memory for one set of coordinates.
 * Allocates a structure for each input conformation, and for each of those allocates
 * space for the x,y,z coordinates of each atom.
 * @param[in] global_options The global options structure
 * @param[in,out] c Pointer to coordinate set to allocate for, with filename,natoms and nstructures specified
 * @return Integer indicating success or failure
 */
int alloc_coords(global_options_struct *global_options, coord_set *c)
{ 
  int i,j;
  // Allocate set of structures
  c->struc = (coords_struct*)malloc(c->num_coords*sizeof(coords_struct));
  if (!c->struc) {
    printf("ERROR! Failed to allocate %d bytes for coords_data\n",c->num_coords*sizeof(coords_struct));
    return ALLOC_FAIL;
  }
  c->mem_allocated = sizeof(coord_set);
    
  for (i=0; i<c->num_coords; ++i)
  {
    c->struc[i].mem_allocated=sizeof(coords_struct);
    c->struc[i].x_coord = (double *) malloc(c->natoms*sizeof(double));
    if (!c->struc[i].x_coord) {
      malloc_failure_double("alloc_coords", "coords_data->x_coord",(c->natoms*sizeof(double)));
      return ALLOC_FAIL;
    }
    c->struc[i].y_coord = (double *)malloc(c->natoms*sizeof(double));
    if (!c->struc[i].y_coord) {
      malloc_failure_char("alloc_coords", "coords_data->y_coord",(c->natoms*sizeof(double)));
      return ALLOC_FAIL;
    }
    c->struc[i].z_coord = (double *)malloc(c->natoms*sizeof(double));
    if (!c->struc[i].z_coord){
    malloc_failure_double("alloc_coords", "coords_data->z_coord",(c->natoms*sizeof(double)));
    return ALLOC_FAIL;
    }
    c->struc[i].mem_allocated += 3*sizeof(double)*c->natoms;
    
    // If fitting forces, allocate space for the force on each atom
    if (global_options->FUNC_TO_FIT==AMBER_FORCES) {
      c->struc[i].force = (force_struct *) malloc(c->natoms*sizeof(force_struct));
      if (!c->struc[i].force) {
        malloc_failure_char("alloc_coords", "coords_data->force",c->natoms*sizeof(force_struct));
        return ALLOC_FAIL;
      }
      c->struc[i].mem_allocated += sizeof(force_struct)*c->natoms;
    }
    else c->struc[i].force = NULL;
    
    c->mem_allocated+=c->struc[i].mem_allocated;
  }
  return SUCCESS;
}

/**
 * Frees the memory used by one set of coordinates
 * @param[in] global_options The global options structure
 * @param[in,out] c Pointer to the coordinate set whose contents will be freed
 * @return Integer indicating success or failure
 */
int free_coords(global_options_struct *global_options, coord_set *c)
{
  int i;
  for (i=0; i<c->num_coords; ++i)
  {
    free(c->struc[i].x_coord);
    free(c->struc[i].y_coord);
    free(c->struc[i].z_coord);
    if (global_options->FUNC_TO_FIT==AMBER_FORCES)
      free(c->struc[i].force);
  }
  free(c->struc);
  free(c->filename);
  free(c->energy_filename);
  return SUCCESS;
}
