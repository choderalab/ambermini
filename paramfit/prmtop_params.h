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

/** @file prmtop_params.h
 * Contains parameters pertaining to prmtop files, including information necessary
 * to reconstruct them, and the data structures into which the prmtop bonds, angles,
 * and dihedrals are inserted in order to make them ready for fitting
 * 
 * Information is duplicated in two forms-- the parm for (to allow re-creation of the
 * parm file if required at some point), and the more obvious form, where each "bond"
 * has associated parameters and optional scale factors.
 */

#define BUFFER_SIZE 1024              /**< Length of buffer for reading in prmtop data */

#define PRMTOP_TITLE_LENGTH 81        /**< Title length of old format prmtop file */

#define MAX_BONDS_PER_TYPE 256        /**< Maximum number of bonds of a given type allowed */
#define MAX_ANGLES_PER_TYPE 512       /**< Maximum number of angles of a given type allowed */
#define MAX_DIHEDRALS_PER_TYPE 1024   /**< Maximum number of dihedrals of a given type allowed */

/**The parameter file assumes the atom, residue, symbol, etc. names to be         *
 * four characters (we will store them as strings, requiring a null terminator,   *
 * hence the size is 5).                                                          */
#define NAME_SIZE 5
#define NAME_DEFAULT "    "
typedef char Name[NAME_SIZE];

/** Defines an atom with data from the prmtop */
typedef struct _atom_struct
{
  double chrg;                        /**<  the charge on the atom            */
  double amass;                       /**<  the mass of the atom              */
  int join;                           /**<  the tree joining info             */
  int irotat;                         /**<  last at to move if cur at moved   */
  int iac;                            /**<  atom types involved in L-J        */
  int numex;                          /**<  index into excluded atom list     */
  int res;                            /**<  the residue number                */
  Name igraph;                        /**<  the true atom name                */
  Name isymbl;                        /**<  the atom type                     */
  Name itree;                         /**<  the atom tree symbol              */
} atom_struct;

/** Defines the residues in the prmtop.                                      */
typedef struct _residue
{
  Name labres;                        /**<  the residue name                 */
  int ipres;                          /**<  the pntr list or all residues    */
} residue;

/** Unprocessed bond information as it appears in the prmtop
 * Sorted for uniquness and turned into a bond_data_struct
 * @see process_prmtop
 */
typedef struct _parmbond_struct {
  int ib;                             /**< first atom in bond               */
  int jb;                             /**< first  second atom in bond              */
  int icb;                            /**< first  pointer to parameters            */
} parmbond_struct;

/** Unprocessed angle information as it appears in the prmtop
 * Sorted for uniqueness and turned into a angle_data_struct
 * @see process_prmtop */
typedef struct _parmangle_struct {
  int it;                             /**< first  first atom in angle              */
  int jt;                             /**< first  second atom in angle             */
  int kt;                             /**< first  third atom in angle              */
  int ict;                            /**< first  pointer to parameters            */
} parmangle_struct;

/** Unprocessed dihedral information as it appears in the prmtop 
 * Sorted for uniqueness and turned into a dihedral_data_struct
 * @see process_prmtop */
typedef struct _parmdihedral_struct {
  int ip;                             /**< first  first atom in dihedral              */
  int jp;                             /**< first  second atom in dihedral             */
  int kp;                             /**< first  third atom in dihedral              */
  int lp;                             /**< first  fourth atom in dihedral             */
  int icp;                            /**< first  pointer to parameters            */
} parmdihedral_struct;

/** Coordinates and parameters of all instances of one unique bond */
typedef struct _bond_data_struct {
  int  number;                        /**< first number of bonds of this type*/
  short int DO_FIT_KR;                /**< first Whether or not to vary this bond term in the fitting - default is yes*/
  short int DO_FIT_REQ;
  char atom_type1[NAME_SIZE];         /**< first The first atom type for this bond type*/
  char atom_type2[NAME_SIZE];         /**< first The second atom type for this bond type*/
  double rk;                          /**< first Force constant for this bond type*/
  double req;                         /**< first Equilibrium bond length for this bond type*/
  int atom1[MAX_BONDS_PER_TYPE];      /**< first The atoms involved in this bond type, - should really be a pointer here but this is quicker*/
  int atom2[MAX_BONDS_PER_TYPE];
} bond_data_struct;

/** Coordinates and parameters of all instances of one unique angle */
typedef struct _angle_data_struct {
  int  number;                        /**< first number of angles of this type*/
  short int DO_FIT_KT;                /**< first Whether or not to vary this angle term in the fitting - default is yes*/
  short int DO_FIT_THEQ;
  char atom_type1[NAME_SIZE];         /**< first The first atom type for this angle type*/
  char atom_type2[NAME_SIZE];         /**< first The second atom type for this angle type*/
  char atom_type3[NAME_SIZE];         /**< first The third atom type for this angle type*/  
  double tk;                          /**< first Force constant for this angle type*/
  double teq;                         /**< first Equilibrium angle for this angle type*/
  int atom1[MAX_ANGLES_PER_TYPE];     /**< first The atoms involved in this angle type, - should really be a pointer here but this is quicker*/
  int atom2[MAX_ANGLES_PER_TYPE];
  int atom3[MAX_ANGLES_PER_TYPE];
} angle_data_struct;

/** Coordinates and parameters of all instances of one unique dihedral. */
typedef struct _dihedral_data_struct {
  short int DO_FIT_KP;                /**< Whether or not to vary this dihedral term in the fitting - default is yes for KP, no for NP and PHASE*/
  short int DO_FIT_NP;
  short int DO_FIT_PHASE; 
  double pk;                          /**< Force constant for this dihedral type, note this is actually Vn/2 */
  double pn;                          /**< periodicity of this dihedral type */
  double phase;                       /**< phase of this dihedral type */  
} dihedral_data_struct;

/** All information about instances of one unique dihedral type
 * Each term of this dihedral is a dihedral_data_struct array in term array
 * This allows for easy addition of more terms to a dihedral and the identification
 * of how many terms there are, etc.
 */
typedef struct _dihedral_type_struct {
  int atom1[MAX_DIHEDRALS_PER_TYPE];  /**< The atoms involved in this dihedral type, - should really be a pointer here but this is quicker */
  int atom2[MAX_DIHEDRALS_PER_TYPE];
  int atom3[MAX_DIHEDRALS_PER_TYPE];
  int atom4[MAX_DIHEDRALS_PER_TYPE];  
  
  char atom_type1[NAME_SIZE];         /**< The first atom type for this dihedral type */
  char atom_type2[NAME_SIZE];         /**< The second atom type for this dihedral type */
  char atom_type3[NAME_SIZE];         /**< The third atom type for this dihedral type */
  char atom_type4[NAME_SIZE];         /**< The fourth atom type for this dihedral type */  
  int num_terms;                      /**< The number of terms of this dihedral (# dihedral_data_struct arrays) */
  short int improper;                 /**<  whether or not this is an improper dihedral */
  dihedral_data_struct *term;         /**< Array of dihedral data structures with data for each term */
  int  number;                        /**< number of dihedrals of this type */
} dihedral_type_struct;

/** Stores all of the parmtop data */
typedef struct _parm_struct {
  int mem_allocated;                  /**< Updated with number of bytes allocated for pointer inside this structure*/
  short int newparm;                  /**< YES if the prmtop format is >=v7.0 else NO*/
  char *title;                        /**< The title in the prmtop file*/
  char *filename;                     /**< The path to this file                                  */
  /* The integer control variables                         */
  int  NTOTAT;                        /**<  total number of atoms in the system                   */
  int  NTYPES;                        /**<  number of AMBER atom types used, max is 60            */
  int  NBONH;                         /**<  number of bonds containing hydrogen                   */
  int  NBONA;                         /**<  number of bonds without hydrogen                      */
  int  NTHETH;                        /**<  number of angles containing hydrogen                  */
  int  NTHETA;                        /**<  number of angles not containing hydrogen              */
  int  NPHIH;                         /**<  number of dihedrals containing hydrogen               */
  int  NPHIA;                         /**<  number of dihedrals not containing hydrogen           */
  int  JHPARM;                        /**<  NOT USED                                              */
  int  JPARM;                         /**<  NOT USED                                              */
  int  NEXT;                          /**<  total number of excluded atoms                        */
  int  NTOTRS;                        /**<  total number of residues                              */
  int  MBONA;                         /**<  NBONA + number of constraint bonds                    */
  int  MTHETS;                        /**<  NTHETS (sic) + number of constraint angles            */
  int  MPHIA;                         /**<  NPHIA + number of constraint dihedral angles          */
  int  MUMBND;                        /**<  total number of unique bond types                     */
  int  MUMANG;                        /**<  total number of unique angle types                    */
  int  MPTRA;                         /**<  total number of unique dihedral types                 */
  int  NATYP;                         /**<  number of "atoms" defined in parameter file           */
  int  NHB;                           /**<  number of types of hydrogen bonded pair interactions  */
  int  IFPERT;                        /**<  =1 if perturbation info is to be read =0 otherwise    */
  int  NBPER;                         /**<  number of bonds to be perturbed                       */
  int  NGPER;                         /**<  number of angles to be perturbed                      */
  int  NDPER;                         /**<  number of dihedrals to be perturbed                   */
  int  MBPER;                         /**<  num of pert bonds across boundary to non-pert groups  */
  int  MGPER;                         /**<  num of pert angles across boundary to non-pert groups */
  int  MDPER;                         /**<  num of pert dihedrals across bndry to non-pert groups */
  int  IFBOX;                         /**<  =1 if periodic box info to be read =0 otherwise       */
  int  NMXRS;                         /**<  number of atoms in the largest residue                */
  int  IFCAP;                         /**<  =1 if CAP option was used in edit, =0 otherwise       */
  int  NUMEXTRA;                      /**<  number of extra points (aka lone pairs)               */
  atom_struct *atom;                  /**< Structure containing the info for each atom*/
  double AMBER_SYSTEM_CHARGE;         /**< Total charge of the system - calculated by summing the charges of the atoms */
  int *nno;                           /**<  index for non-bond of each type                       */
  residue *residue;                   /**< Stores the residue info                                */
  double *rk;                         /**< Bond force constants                                   */
  double *req;                        /**< Bond equilibrium constants                             */
  double *tk;                         /**< Angle force constants                                  */
  double *teq;                        /**< Angle equilibrium constants                            */
  double *pk;                         /**< dihedral force constants - note these are actually vn/2*/
  double *pn;                         /**< Dihedral periodicity constants                         */
  double *phase;                      /**< Dihedral phase constants                               */
  double *solty;                      /**< NOT USED BUT READ FROM PRMTOP ANYWAY                   */
  double *cn1;                        /**<  L-J r**12 and r**6 for all pos...                     */
  double *cn2;                        /**<  ...atom type interactions                             */
  parmbond_struct *pbondH;            /**<  bonds with hydrogen                                   */
  parmbond_struct *pbond;             /**<  bonds without hydrogen                                */
  parmangle_struct *pangleH;          /**<  angles with hydrogen                                  */
  parmangle_struct *pangle;           /**<  angles without hydrogen                               */
  parmdihedral_struct *pdihedralH;    /**<  dihedrals with hydrogen                               */
  parmdihedral_struct *pdihedral;     /**<  dihedrals without hydrogen                            */
  int *natex;                         /**<  excluded atom list                                    */
  double *ag;                         /**<  H-bond r**12 and r**10...                             */
  double *bg;                         /**<  ...                                                   */
  double *hbcut;                      /**<  NO LONGER USED                                        */
  double K;                           /**< Constant to adjust axis to be same for quantum and MM  */
  bond_data_struct *bond_data;        /**< Contains the bond data split into specific bond types  */
  angle_data_struct *angle_data;      /**< Contains the angle data split into specific angle types*/
  dihedral_type_struct *dihedral_data;/**< Contains the dihedral data split into specific dihedral types*/  
  int unique_bonds_found;   
  int unique_angles_found;   
  int unique_dihedrals_found;         /**<  Total number of unique dihedral types. Each can have any number of terms */
  int unique_dihedral_terms;          /**<  Total number of unique dihedral terms                 */
  int *fit_atom;                      /**<  Marks atoms as being relevant if forces are to be fit */
  int ndimensions;                    /**<  The number of dimensions of fit in this structure     */
} parm_struct;

