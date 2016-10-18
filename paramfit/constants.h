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
 
/*Constants.h*/

/*Contains conversion parameters and constants*/

/*1 Hartree is 4.3597482x10^-18J x 6.02214x10^23 mol-1
  1 Cal = 4.184 J (Exact)
  Therefore 1 Hartree = 627509.8954 calories/mol = 627.5098954 kcal/mol
*/
#define HARTREE_TO_KCALMOL        627.5098954
#define KJMOL_TO_KCALMOL          1.0/4.184   /*Exact*/
#define BOHR_TO_ANGSTROM          1.0/0.529

#define PI                        3.141592653589793238462643383279 /*OK it's only in fun:-)*/

/*360 degree = 2PI radians
  Therefore 1 degree = 2PI/360 radians
*/
#define DEGREE_TO_RADIAN          (2.0*PI)/360.0
#define RADIAN_TO_DEGREE          360.0/(2.0*PI)


