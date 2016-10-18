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

/** @file misc_utils.c
 * Contains a number of small utilities used in the program that don't
 * really have anywhere else to go.
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "constants.h"

#include "function_def.h"

/**
 * Prints the closing of a line in a box styled like
 * --------------
 * |            |
 * --------------
 * @see print_open_line_box
 * @param[in] no_spaces The number of spaces to go before the line
 */
void print_close_line_box(int no_spaces)
{
  register int i;
  for (i=0;i<no_spaces;i++)
  {
    printf(" ");
  }
  printf("|\n");
}

/**
 * Prints the beginning of an ASCII box.
 * @see print_close_line_box
 * @param[in,out] i Contains number of characters printed
 */
void print_open_line_box(int *i)
{
  *i=printf("        | ");
}

/**
 * Checks for invalid characters in a filename string.
 * @param[in] data_string The string to check
 * @param[in] length The length of the string
 * @return Integer representing SUCCESS or INVALID_FORMAT
 */
int check_for_valid_filename(const char *data_string, const int length)
{
  /*This won't catch everything but is useful*/
  int temp_count;
  for (temp_count=0;temp_count<length;temp_count++)
  {
    if (
       (*(data_string+temp_count)<32)
       ||
       (((*(data_string+temp_count))>33)&&((*(data_string+temp_count))<43))
       ||
       (((*(data_string+temp_count))>57)&&((*(data_string+temp_count))<65))
       ||
       (*(data_string+temp_count)==92)
       ||
       (*(data_string+temp_count)==124)
       ||
       (*(data_string+temp_count)>126)
       )
       return INVALID_FORMAT; /*error invalid character for filename*/
  }
  return SUCCESS;
}

/**
 * "Safe" version of getline, exits with an error if EOF is found.
 * @param[in,out] line String that will contain the read line
 * @param[in] max The maximum number of characters to read
 * @param[in] fp File pointer to the file to read a line from
 * @return The number of characters read, or FILE_READ_FAIL if EOF found
 */
int s_getline( char *line, int max, FILE *fp )
{
  if( fgets( line, max, fp ) == NULL )
  {
    return FILE_READ_FAIL;
  }
  else
  {
    return strlen( line );    
  }
}


 /**
  * Used to locate a title in a prmtop.
  * Leaves the file pointer positioned at the line below the
  * selectected label. Rewinds the file pointer before searching.
  * 
  * @param[in,out] fptr The file pointer to use
  * @param[in] label The title to find
  * @return Integer SUCCESS, or INVALID_LINE if label not found.
  */
int find_flag( FILE *fptr, char *label )
{

  int retval;
  int i;
  char *buffer;

  buffer = (char *)calloc(BUFFER_SIZE,sizeof(char));
  if (buffer==NULL)
  {
  malloc_failure_char("find_flag", "buffer", BUFFER_SIZE);
  return ALLOC_FAIL;
  }
  
  rewind(fptr);
  while (TRUE) /*loop forever since this will be broken by an internal return statement for either success or an EOF*/
  {
      while (strncmp(buffer, "%FLAG", 5) != 0) /*loop through until we hit a flag*/
      {
        retval = s_getline( buffer, BUFFER_SIZE, fptr );
        /*s_getline returns FILE_READ_FAIL on EOF*/
        if ( retval==FILE_READ_FAIL )
        {
          free(buffer);
          buffer = NULL;
          return INVALID_LINE;
        }
      }

      /*found the flag, now check if it is the correct label*/
      if ( !strncmp( buffer+6, label, strlen(label) ) )
      {
        /*We found the label we were looking for*/
        free(buffer);
        buffer = NULL;
        return SUCCESS;
      }
      /*if we got to here then it wasn't the correct flag
        clear the beginning of the buffer so the next strncmp doesn't trigger on the current contents
      */
      for (i=0;i<=11;++i)
        buffer[i]=' ';
  }   
}

/**
 * Desigend to read 4 characters from a file stream and null terminate the string.
 * @see prmtop_params.h for definition of NAME_SIZE
 * @param[in,out] fptr File pointer to read from
 * @param[in,out] stringp The string to read into
 * @return Integer indicating SUCCESS or INVALID_LINE
 */
int name_copy( FILE *fptr, char *stringp)
{
  int i;
  int ierror;
  for (i = 0; i < NAME_SIZE - 1; ++i)
  {
    do
    {
      ierror = fscanf(fptr, "%c", &stringp[i]);
      if ( ierror == EOF ) /*Hit end of file*/
        return INVALID_LINE;
      if (stringp[i] == ' ' && i!=1) stringp[i] = (char) 0; // no spaces allowed in atom name after 1 char

    } while ( stringp[i] == '\n');

    if ( ierror == 0 ) /*Got a carriage return as part of the atom name*/
        return INVALID_LINE;
  }
  stringp[NAME_SIZE-1] = (char) 0;
  return SUCCESS;
}

/**
 * Translates atom numbers from prmtop format into understandable numbers.
 * PARM for some reason obfuscates each atom number by (when positive),
 * at = (at+3) /3. This function returns the true atom number. For example,
 * unObfuscateAtom(24) = 9.
 * @see ObfuscateAtom
 * @param[in] at The obfuscated atom to translate
 * @return The true atom number
 */
int unObfuscateAtom( int at )
{
  if ( at < 0 )
    at = ( -at + 3 ) / 3;
  else
    at = ( at + 3 ) / 3;
  return( at );
}

/**
 * Re-obfuscates an atom from atomic number into prmtop format.
 * This can be used in writing a prmtop.
 * @see unObfuscateAtom
 * @param[in] at The atom to obfuscate
 * @return The prmtop obfuscated atom number.
 */
int ObfuscateAtom( int at )
{
  at = (3*at) - 3;
  return at;
}

/**
 * Calculates the distance between 2 3-dimenstional points.
 * 
 * \f$ d=\sqrt{ (x_1-x_2)^2 + (y_1-y_2)^2 + (z_1-z_2)^2 } \f$
 */
double calc_bond_length(double bond1x, double bond1y, double bond1z, double bond2x, double bond2y, double bond2z)
{

  return(sqrt((bond1x-bond2x)*(bond1x-bond2x)+(bond1y-bond2y)*(bond1y-bond2y)+(bond1z-bond2z)*(bond1z-bond2z)));

}


/**
 * Calculates the angle between 3 3-dimenstional points
 * Angle between \f$ \vec{a_1},\vec{a_2},\vec{a_3} \f$ is defined as:
 * \f$ \vec{v_1} = \vec{a_1} - \vec{a_2} \f$
 * \f$ \vec{v_2} = \vec{a_3} - \vec{a_2} \f$
 * 
 * \f$ \theta = \cos^{-1}\frac{\vec{v_1}\cdot\vec{v_2}}{|\vec{v_1}||\vec{v_2}|} \f$
 * 
 * @return The angle between the three atoms 
 */
double calc_angle_radians(double atom1x, double atom1y, double atom1z, double atom2x, double atom2y, double atom2z,
                          double atom3x, double atom3y, double atom3z)
{
  double temp_angle;
  double vector1[3];
  double vector2[3];
  double modvector1;
  double modvector2;
  double vectordot;

  vector1[0]=atom1x-atom2x;
  vector1[1]=atom1y-atom2y;
  vector1[2]=atom1z-atom2z;
  vector2[0]=atom3x-atom2x;
  vector2[1]=atom3y-atom2y;
  vector2[2]=atom3z-atom2z;

  modvector1=sqrt((vector1[0]*vector1[0])+(vector1[1]*vector1[1])+(vector1[2]*vector1[2]));
  modvector2=sqrt((vector2[0]*vector2[0])+(vector2[1]*vector2[1])+(vector2[2]*vector2[2]));

  vectordot=(vector1[0]*vector2[0])+(vector1[1]*vector2[1])+(vector1[2]*vector2[2]);

  temp_angle=acos(vectordot/(modvector1*modvector2));
  return temp_angle;

}

/** Calculates the dihedral between 4 points in 3 dimensions
 * Given four atoms A,B,C & D:
 * @verbatim
   A      D
    \    /
     B--C
 * @endverbatim
 * The torsion angle along the torsion axis B--C is the angle between the planes ABC and BCD.
 * The best way to calculate this is in terms of 3 vectors:
 *
 * a = A->B
 * b = B->C
 * c = C->D
 *
 * In terms of these vectors the dihedral angle is then:
 * \f$ \theta = \cos^{-1}\frac{(\vec{a}\times\vec{b})(\vec{b}\times\vec{c})}{|\vec{a}\times\vec{b}||\vec{b}\times\vec{c}|} \f$
 *
 * The sign of the dihedral is then found from the triple scalar product calculated by evaluating the
 * determinant of the matrix:
 * @verbatim
   | a(x) a(y) a(z) |    | a(x) a(y) a(z) |
   | b(x) b(y) b(z) | -> | b(x) b(y) b(z) |
   | c(x) c(y) c(z) |    | c(x) c(y) c(z) |
 * @endverbatim
 * which is:
 *  a[x] * (b[y]*c[z]-c[y]*b[z]) - a[y]*(b[x]*c[z]-c[x]*b[z]) + a[z]*(b[x]*c[y]-c[x]*b[y])
 * @return The angle between the input atoms, in radians
 */
double calc_dihedral_radians(double atom1x, double atom1y, double atom1z, double atom2x, double atom2y, double atom2z,
                          double atom3x, double atom3y, double atom3z, double atom4x, double atom4y, double atom4z)
{
 double vector1x, vector1y, vector1z;
 double vector2x, vector2y, vector2z;
 double vector3x, vector3y, vector3z;
 double v1Xv2xv2Xv3;
 double v1Xv2;
 double v2Xv3;
 double cosdihe;
 double dihedral;
 double sign;
 
 vector1x = atom2x - atom1x;
 vector1y = atom2y - atom1y;
 vector1z = atom2z - atom1z;

 vector2x = atom2x - atom3x;
 vector2y = atom2y - atom3y;
 vector2z = atom2z - atom3z;

 vector3x = atom3x - atom4x;
 vector3y = atom3y - atom4y;
 vector3z = atom3z - atom4z;

 v1Xv2xv2Xv3 = (vector1y * vector2z - vector1z * vector2y) * (vector2y * vector3z - vector2z * vector3y) +
               (vector2x * vector1z - vector1x * vector2z) * (vector3x * vector2z - vector2x * vector3z) +
               (vector1x * vector2y - vector2x * vector1y) * (vector2x * vector3y - vector3x * vector2y);
               


 v1Xv2 = (vector1y * vector2z - vector1z * vector2y) * (vector1y * vector2z - vector1z * vector2y) +
         (vector2x * vector1z - vector1x * vector2z) * (vector2x * vector1z - vector1x * vector2z) +
         (vector1x * vector2y - vector2x * vector1y) * (vector1x * vector2y - vector2x * vector1y);

 v2Xv3 = (vector2y * vector3z - vector2z * vector3y) * (vector2y * vector3z - vector2z * vector3y ) +
         (vector3x * vector2z - vector2x * vector3z) * (vector3x * vector2z - vector2x * vector3z ) +
         (vector2x * vector3y - vector3x * vector2y) * (vector2x * vector3y - vector3x * vector2y );
         
 cosdihe = v1Xv2xv2Xv3 / sqrt (v1Xv2 * v2Xv3);
 if (cosdihe > 1.0) // take care of some floating point issues that will sometimes cause an out of bounds to the arccos
   cosdihe = 1.0;
 if (cosdihe < -1.0)
   cosdihe = -1.0;
 dihedral = acos(cosdihe);

 sign = vector1x * (vector2y * vector3z - vector3y * vector2z) - vector1y * (vector2x * vector3z - vector3x * vector2z)
        + vector1z * (vector2x * vector3y - vector3x * vector2y);
        
 if (sign>0)
 {
    if (dihedral>=0) return dihedral;
    if (dihedral<0)  return -dihedral;
 }
 else if (sign<0)
 {
    if (dihedral>0) return -dihedral;
    if (dihedral<=0) return dihedral;
 }

 return dihedral; /*Should never actually gets to this line due to above but this suppresses compiler warnings about
                    hitting the end of a non-void function*/
}

/**
 * Calculates the total number of parameters to be fit- the dimensionality of the problem.
 * TODO - update with multiple prmtops as right now it is a simple loop
 * @param[in,out] global_options The global options structure, where NDIMENSIONS is to be updated
 * @param[in] parm_datas Array of parameter data structures
 */
void calc_fit_dimensions(global_options_struct *global_options, parm_struct *parm_datas)
{
  int i,j,p;
  /*Calculates the number of dimensions of the fit*/
  for (p=0; p<global_options->num_prmtops; ++p) {
    parm_datas[p].ndimensions=0;
  if (global_options->K_FIT==YES)
    ++parm_datas[p].ndimensions;
    for (i=0;i<parm_datas[p].unique_bonds_found;++i) {
      if (parm_datas[p].bond_data[i].DO_FIT_KR==YES)
        ++parm_datas[p].ndimensions;
      if (parm_datas[p].bond_data[i].DO_FIT_REQ==YES)    
        ++parm_datas[p].ndimensions;
    }
    for (i=0;i<parm_datas[p].unique_angles_found;++i) {
      if (parm_datas[p].angle_data[i].DO_FIT_KT==YES)
        ++parm_datas[p].ndimensions;
      if (parm_datas[p].angle_data[i].DO_FIT_THEQ==YES)
        ++parm_datas[p].ndimensions;
    }
    for (i=0;i<parm_datas[p].unique_dihedrals_found;++i) {
      for (j=0; j<parm_datas[p].dihedral_data[i].num_terms; ++j) {
        if (parm_datas[p].dihedral_data[i].term[j].DO_FIT_KP==YES)
          ++parm_datas[p].ndimensions;
        if (parm_datas[p].dihedral_data[i].term[j].DO_FIT_NP==YES)
          ++parm_datas[p].ndimensions;
        if (parm_datas[p].dihedral_data[i].term[j].DO_FIT_PHASE==YES)
          ++parm_datas[p].ndimensions;
      }
    }
  }
}

/**
 * Puts the parameters in a data array into each prmtop. Used in multi-prmtop
 * fits to update each of them before conducting an energy evaluation.
 * @param[in] global_options The global options structure
 * @param[in,out] parm_data Array of parameter data structures that will be updated.
 * @param[in] parameters Pre-allocated array size NDIMENSIONS to read from
 * @return Integer indicating success or failure.
 */
int update_prmtop_data(global_options_struct *global_options, parm_struct *parm_datas, double *parameters)
{
  int i;
  int num;
  int pnum=0;
  for (i=0; i<global_options->num_prmtops; ++i) {
    num = modify_params_scratch_data(global_options, &(parm_datas[i]), parameters, WRITE);
    if (pnum && num!=pnum) {
      printf("ERROR! Differing number of parameters to fit in prmtops %d and %d\n", i,i-1);
      return INVALID_DATA;
    }
    pnum=num;
  }
  return SUCCESS;
}

/**
 * Takes the parameters from the prmtop and puts them in a data array. Used in multi-prmtop fits to
 * update the data matrix following changes to any prmtop (such as bounds check or initial settings)
 * 
 * Uses only the parameters from the first prmtop because they are all kept consistent 
 * 
 * @param[in] global_options The global options structure
 * @param[in,out] parm_data Array of parameter data structures that will be updated.
 * @param[in] parameters Pre-allocated array size NDIMENSIONS to read from
 * @return Integer indicating success or failure.
 */
int update_data_matrix(global_options_struct *global_options, parm_struct *parm_datas, double *parameters) 
{
  if (parameters == NULL) {
    printf("ERROR: Attempting to update uninitialized data matrix in update_data_matrix(%p, %p, %p)",
           global_options, parm_datas, parameters);
    return FAILURE;
  }
 
  int retval = modify_params_scratch_data(global_options, &(parm_datas[0]), parameters, READ);
  if (retval==parm_datas[0].ndimensions) return SUCCESS;
  else return FAILURE;
}

/**
 * Reads and writes parameters in the parm_struct data structure.
 * Extracts parameters marked as variable from the prmtop and puts them in the linear array, or takes
 * parameters from the array and puts them in the data structure.
 * The order in which the extraction done is as follows:
 *  K
 *  BOND x (KR, KEQ)
 *  BOND y (KR, KEQ)
 *  .
 *  ANGLE x (KT, THEQ)
 *  ANGLE y (KT, THEQ)
 *  .
 *  DIHEDRAL x (KP, PN, PHASE)
 *  DIHEDRAL y (KP, PN, PHASE)
 *  .
 * @param[in] global_options The global options structure
 * @param[in,out] parm_data The parameters data structure to read or write to
 * @param[in,out] parameters Pre-allocated array size NDIMENSIONS to read or write to
 * @param[in] MODE if READ, copy the data from parm_data to parameters. if WRITE, copy from parameter to parm_data.
 * @return The number of parameters successfully extracted.
 * 
 *  I am aware that this whole procedure here is clunky and slow but it makes it significantly easier to understand
 *  and debug. At some point I will replace this whole system with a much more efficient method
 */
int modify_params_scratch_data(global_options_struct *global_options, parm_struct *parm_data, double *parameters, readwrite_t MODE)
{
  int number_extracted;
  int i,j;

  number_extracted=0;
// 
  if (MODE==READ)
  {
    if (global_options->K_FIT==YES)
    {
      parameters[number_extracted]=parm_data->K;
      ++number_extracted;
    }
    /*Now do all the bonds*/
    for (i=0;i<parm_data->unique_bonds_found;++i)
    {
      /*For each bond in turn check if KR and/or REQ are variable*/
      if (parm_data->bond_data[i].DO_FIT_KR==YES)
      {
         parameters[number_extracted]=parm_data->bond_data[i].rk;
         ++number_extracted;
      }
      if (parm_data->bond_data[i].DO_FIT_REQ==YES)
      {
         parameters[number_extracted]=parm_data->bond_data[i].req;
         ++number_extracted;
      }
    }
    /*Now do the angles*/
    for (i=0;i<parm_data->unique_angles_found;++i)
    {
      /*For each angle in turn check if KT and/or THEQ are variable*/
      if (parm_data->angle_data[i].DO_FIT_KT==YES)
      {
         parameters[number_extracted]=parm_data->angle_data[i].tk;
         ++number_extracted;
      }
      if (parm_data->angle_data[i].DO_FIT_THEQ==YES)
      {
         parameters[number_extracted]=parm_data->angle_data[i].teq;
         ++number_extracted;
      }
    }
    /*now do the dihedrals*/
    for (i=0;i<parm_data->unique_dihedrals_found;++i) {
      for (j=0; j<parm_data->dihedral_data[i].num_terms; ++j) {
        /*For each diheedral in turn check if KP and/or NP and/or PHASE are variable*/
        if (parm_data->dihedral_data[i].term[j].DO_FIT_KP==YES)
        {
          parameters[number_extracted]=parm_data->dihedral_data[i].term[j].pk;
          ++number_extracted;
        }
        if (parm_data->dihedral_data[i].term[j].DO_FIT_NP==YES)
        {
          parameters[number_extracted]=parm_data->dihedral_data[i].term[j].pn;
          ++number_extracted;
        }
        if (parm_data->dihedral_data[i].term[j].DO_FIT_PHASE==YES)
        {
          parameters[number_extracted]=parm_data->dihedral_data[i].term[j].phase;
          ++number_extracted;
        }
      }
    }
  }
  else if (MODE==WRITE)
  {
    /* Write the variables into the parm structure */
    
    if (global_options->K_FIT==YES)
    {
      parm_data->K=parameters[number_extracted];
      ++number_extracted;
    }
    
    /*Now do all the bonds*/
    for (i=0;i<parm_data->unique_bonds_found;++i)
    {
      /* KR stays between 100 and 1000 */
      if (parm_data->bond_data[i].DO_FIT_KR==YES)
      {
        parm_data->bond_data[i].rk=parameters[number_extracted];
        ++number_extracted;
      }
      /* Length between 0 and 3 Angstroms */
      if (parm_data->bond_data[i].DO_FIT_REQ==YES)
      {
        parm_data->bond_data[i].req=parameters[number_extracted];
        ++number_extracted;
      }
    }
    /*Now do the angles*/
    for (i=0;i<parm_data->unique_angles_found;++i)
    {
      if (parm_data->angle_data[i].DO_FIT_KT==YES) {
        parm_data->angle_data[i].tk=parameters[number_extracted];
        ++number_extracted;
      }
      if (parm_data->angle_data[i].DO_FIT_THEQ==YES) {
        parm_data->angle_data[i].teq=parameters[number_extracted];
        ++number_extracted;
      }
    }
    /*now do the dihedrals*/
    for (i=0;i<parm_data->unique_dihedrals_found;++i) {
      for (j=0; j<parm_data->dihedral_data[i].num_terms; ++j) {
        /*For each diheedral in turn check if KP and/or NP and/or PHASE are variable*/
        // PK stays between -30 and 30 (Glycam ff lists some negative kp, for example)
        if (parm_data->dihedral_data[i].term[j].DO_FIT_KP==YES) {
          parm_data->dihedral_data[i].term[j].pk=parameters[number_extracted];
          ++number_extracted;
        }
        if (parm_data->dihedral_data[i].term[j].DO_FIT_NP==YES) {
          parm_data->dihedral_data[i].term[j].pn = parameters[number_extracted];
          ++number_extracted;
        }
        if (parm_data->dihedral_data[i].term[j].DO_FIT_PHASE==YES) {
          parm_data->dihedral_data[i].term[j].phase=parameters[number_extracted];
          ++number_extracted;
        }
      }
    }
  }
  else
  {
    printf("   ERROR IN modify_params_scratch_data() - MODE = %d - UNKNOWN MODE\n",MODE);
    printf("          NON-FATAL, ATTEMPTING TO CONTINUE.\n");
    fflush(stdout); /*Flush the printf buffer*/
  }

  return(number_extracted);
}

/**
 * Prints a backtrace for debugging.
 * This does not work in Cygwin and so is ifdef'd out in that case.
 * @param[in] signal The signal that is caught
 */ 
void print_backtrace(int signal)
{
#if defined(CYGWIN) || defined(WIN32)
  fprintf(stderr, "Error: signal %d\n",signal);
#else
  void *array[10];
  size_t size;
  size = backtrace(array, 10);
  
  fprintf(stderr, "Error: signal %d:\n", signal);
  backtrace_symbols_fd(array,size,2);
#endif
  exit(1);
}

/**
 * Prints an error message on a signal.
 * This may someday print something more useful.
 * @param[in] signal The signal that is caught
 */
void handle_sigint(int param)
{
  printf("\n!   ERROR: Program terminated.\n");
  exit(ABORT);
}

/**
 * Calculates the number of parameters to be fit according to MODE
 * Used to find number of bond, angle, and dihedral params separately.
 * 
 * @param[in] parm_data The parameter structure to examine
 * @param[in] MODE BONDS, ANGLES, or DIHEDRALS indicating parameters to calculate
 * @return The number of parameters of type MODE that are to be fit
 */
int calculate_no_fit_params(parm_struct *parm_data, short int MODE)
{
  int number_params;
  int i,j;
  number_params=0;

  if (MODE==BONDS)
  {
    for (i=0;i<parm_data->unique_bonds_found;++i)
    {
      /*For each bond in turn check if KR and/or REQ are variable*/
      if (parm_data->bond_data[i].DO_FIT_KR==YES)
         ++number_params;
      if (parm_data->bond_data[i].DO_FIT_REQ==YES)
         ++number_params;
    }
  }
  else if (MODE==ANGLES)
  {
    for (i=0;i<parm_data->unique_angles_found;++i)
    {
      if (parm_data->angle_data[i].DO_FIT_KT==YES)
         ++number_params;
      if (parm_data->angle_data[i].DO_FIT_THEQ==YES)
         ++number_params;
    }
  }
  else if (MODE==DIHEDRALS)
  {
    for (i=0;i<parm_data->unique_dihedrals_found;++i)
    {
      for (j=0; j<parm_data->dihedral_data[i].num_terms; ++j) {
        if (parm_data->dihedral_data[i].term[j].DO_FIT_KP==YES)
          ++number_params;
        if (parm_data->dihedral_data[i].term[j].DO_FIT_NP==YES)
          ++number_params;
        if (parm_data->dihedral_data[i].term[j].DO_FIT_PHASE==YES)
          ++number_params;
      }
    }
  }
  else
  {
    printf("   ERROR IN calculate_no_fit_params() - RW = %d - UNKNOWN MODE\n",MODE);
    printf("             NON-FATAL, ATTEMPTING TO CONTINUE.\n");
    fflush(stdout); /*Flush the printf buffer*/
  }

  return(number_params);
}

/**
 * Checks if two dihedral types are the same or different.
 * Checks with name only, not by the value of the parameters, and only
 * checks in one direction.
 * 
 * @param[in] first Pointer to first dihedral_type_struct to examine
 * @param[in] second Pointer to dihedral_type_struct to compare to
 * @return YES if they are the same, NO if not
 */
int dihedral_types_equal(dihedral_type_struct *first, dihedral_type_struct *second)
{
  if (  !strcmp(first->atom_type1, second->atom_type1) &&
        !strcmp(first->atom_type2, second->atom_type2) &&
        !strcmp(first->atom_type3, second->atom_type3) &&
        !strcmp(first->atom_type4, second->atom_type4) )
    return YES;
  else
    return NO;
}

/**
 * Pretty prints one dihedral value and parameters
 * @param[in] type Pointer to the dihedral_type_struct to print
 * @param[in] term Which term of this dihedral the parameters should be printed for
 */
void print_dihedral(dihedral_type_struct *type, int term)
{
  printf("%s-%s-%s-%s: KP: %.3f NP: %.3f Phase: %.3f\n", type->atom_type1, type->atom_type2, type->atom_type3, type->atom_type4, 
         type->term[term].pk, type->term[term].pn, type->term[term].phase*RADIAN_TO_DEGREE);
}

/**
 * Compare function for qsort in the coordinate structures.
 * Used to compare two coords structures by energy.
 * Please don't call this directly, just pass a function pointer to qsort.
 */
int compare_energy(const void *a, const void *b)
{
  coords_struct *ia = (coords_struct*)a;
  coords_struct *ib = (coords_struct*)b;
  
  if (ia->energy > ib->energy)
    return 1;
  else if (ia->energy < ib->energy)
    return -1;
  else
    return 0;
}

/**
 * Potentially adds more dihedral terms to existing ones so you can get a better fit.
 * This is completely deprecated by dihedral data refactoring and so is 
 * commented out for now.
 */
int not_enough_dihedrals(parm_struct *parm_data, int n)
{
  printf("Not_enough_dihedrals currently deprecated!\n");
  return FAILURE;
/*
  // re-alloc the dihedral data
  dihedral_data_struct *more = (dihedral_data_struct*)malloc(sizeof(dihedral_data_struct)*parm_data->unique_dihedrals_found);
  if (more==NULL)
  {
    printf("   ERROR: Failed to reallocate dihedral structures.\n");
    exit(ALLOC_FAIL);
  }
  
  // copy the data into the larger structure
  int i, j, k;
  short int new_copy;
  int counter=0;
  int n_dihedrals = parm_data->unique_dihedrals_found; // since parm_data's will change as more dihedral terms are added
  
  for (i=0; i<n_dihedrals; ++i)
  {
    // only add to dihedrals that will be fit that don't have enough terms already
    if (parm_data->dihedral_data[i].num_terms < n && 
      (parm_data->dihedral_data[i].DO_FIT_KP || parm_data->dihedral_data[i].DO_FIT_NP || parm_data->dihedral_data[i].DO_FIT_PHASE) )
    {
      // Realloc the temporary dihedral data structure for room for n-num more
      parm_data->unique_dihedrals_found += (n-parm_data->dihedral_data[i].num_terms);
      more = (dihedral_data_struct*)realloc(more, sizeof(dihedral_data_struct)*parm_data->unique_dihedrals_found);
      new_copy = YES;
    }
    else // just copy this one over
      new_copy = NO;
    
    for(j=0; j<( new_copy==YES ? n-parm_data->dihedral_data[i].num_terms+1 : 1); ++j)
    {
      // copy basic info about this dihedral
      strcpy(more[counter].atom_type1, parm_data->dihedral_data[i].atom_type1);
      strcpy(more[counter].atom_type2, parm_data->dihedral_data[i].atom_type2);
      strcpy(more[counter].atom_type3, parm_data->dihedral_data[i].atom_type3);
      strcpy(more[counter].atom_type4, parm_data->dihedral_data[i].atom_type4);
      more[counter].improper = parm_data->dihedral_data[i].improper;
      
      // copy info about which parameters will be fit
      more[counter].phase = parm_data->dihedral_data[i].phase;
      more[counter].pk = parm_data->dihedral_data[i].pk;
      more[counter].DO_FIT_PHASE = parm_data->dihedral_data[i].DO_FIT_PHASE;
      more[counter].DO_FIT_KP = parm_data->dihedral_data[i].DO_FIT_KP;
      if (new_copy==YES)
      {
        // give it a random periodicity
        more[counter].pn = (rand()&1) ? (double)(rand()%6)+1.0 : -1.0*((double)(rand()%6)+1.0);
        more[counter].DO_FIT_NP=YES;
      }
      else
      {
        more[counter].pn = parm_data->dihedral_data[i].pn;
        more[counter].DO_FIT_NP = parm_data->dihedral_data[i].DO_FIT_NP;
      }
      
      // copy info about atoms in this type of dihedral
      more[counter].number = parm_data->dihedral_data[i].number;
      for (k=0; k<MAX_DIHEDRALS_PER_TYPE; ++k)
      {
        more[counter].atom1[k] = parm_data->dihedral_data[i].atom1[k];
        more[counter].atom2[k] = parm_data->dihedral_data[i].atom2[k];
        more[counter].atom3[k] = parm_data->dihedral_data[i].atom3[k];
        more[counter].atom4[k] = parm_data->dihedral_data[i].atom4[k];
      }
      ++counter;
    }
  }
  // remove the old array and replace it with our new one - this needs to be done so  
  // all our dihedrals stay in order.                                                  
  free(parm_data->dihedral_data);
  parm_data->dihedral_data = more;
  
  printf("*  Successfully set %d dihedral terms.\n", n);
  return SUCCESS;
*/
}


