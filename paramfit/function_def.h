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

/*function_def.h*/

#include "prmtop_params.h"
#include "forces.h"

#include <stdio.h>

/*Contains constant declarations, function definitions etc.*/

/*Return codes*/
/*
      Note, greater than ABORT = not a failure, an exit code
      Less than ABORT is a real failure if not zero
      0 = SUCCESS - success
      -1 = FAILURE - Generic failure
      -2 = NOT_IMPLEMENTED not implemented
      -101 = UNKNOWN_OPT failure unknown option
      -201 = TOO_MANY_OPT failure too many options
      -301 = ALLOC_FAIL allocation failure
      -401 = FILE_OPEN_FAIL File open failure - Read
      -1000 = ABORT
      -1001 = CMD_HELP_REQ User requested command line help
      -1002 = HIST_REQ User requested program history
      -1003 = HELP_REQ User requested a help topic, print it and return 6
      -2001 = INVALID_FORMAT
      -2002 = INVALID_DATA
      -2003 = INVALID_LINE
      -3001 = EXCEEDEDMAXITERATIONS - EXCEEDED MAXIMUM ITERATIONS
      -3002 = MINSTATIC Function no longer varying
      -3003 = DATA_OVERFLOW
      -3004 = UNKNOWN_ELEMENT
      
Note lower than -1000 are abort codes, not error codes
*/
#define _GNU_SOURCE

#define SUCCESS  -0
#define FAILURE  -1
#define NOT_IMPLEMENTED -2
#define UNKNOWN_OPT    -101
#define TOO_MANY_OPT   -201
#define ALLOC_FAIL     -301
#define FILE_OPEN_FAIL -401
#define FILE_READ_FAIL -501
#define ABORT        -1000
#define CMD_HELP_REQ -1001
#define HIST_REQ     -1002
#define HELP_REQ     -1003
#define INVALID_FORMAT -2001
#define INVALID_DATA   -2002
#define INVALID_LINE   -2003
#define EXCEEDEDMAXITERATIONS -3001
#define MINSTATIC -3002
#define DATA_OVERFLOW -3003
#define UNKNOWN_ELEMENT -3004

#define OFF 0
#define ON  1
#define DEBUG 2
#define WARN 3

#define BONDS 1
#define ANGLES 2
#define DIHEDRALS 3

/*MAXIMUM VALUES FOR ITERATIONS ETC*/
#define NSIMPLEX_INNER_PER_DIM 25 /*This is multiplied by the number of dimensions in order to work out how many inner loops to run*/
#define NSIMPLEX_OUTER_MAX 100000
#define RAND_RATIO    0.2  /*Ratio by which to multiply the random number received from rand() when
adjusting simplex vertices by gamma*/                                                                

/** Available minimising functions */
typedef enum
{
  SUM_SQUARES_AMBER_STANDARD,
  AMBER_FORCES,
  DIHEDRAL_LEAST_SQUARES
} function_t;

/** Available minimising algorithms */
typedef enum
{
  SIMPLEX,
  GENETIC,
  BOTH,
  NONE
} algorithm_t;

/** Input QM file formats enum */
typedef enum
{
  GAUSSIAN,
  ADF,
  GAMESS
} qm_format_t ;

/** Input coordinate file formats enum */
typedef enum
{
  TRAJECTORY,
  RESTART
} coord_format_t ;

/** Create a boolean type since this is in C */
typedef enum
{
  YES  = 1,
  TRUE = 1,
  NO   = 0,
  FALSE = 0
} bool_t;

/** Read or write, used for file i/o */
typedef enum
{
  READ,
  WRITE
} readwrite_t;

/** Settings for modes for parameter setting or getting */
typedef enum
{
  DEFAULT,                               // fit a default set of parameters
  LOAD,                                  // read in parameters from a previously created file
  SAVE,                                  // save parameters to an output file
  K_ONLY                                 // only fit K
} parameter_mode_t;

/** Enum for paramfit's three run types */
typedef enum
{
  CREATE_INPUT,                          // write input files for quantum package
  FIT,                                   // do the actual fitting
  SET_PARAMS                             // set which parameters should be fit
} runtype_t;

/** Units for energy available */
typedef enum
{
  HARTREE,
  KCALMOL,
  KJMOL
} energy_t;

/** Units for force available */
typedef enum
{
  HARTREE_BOHR,
  KCALMOL_ANGSTROM
} force_t;

/** Level of verbosity available */
typedef enum
{
  LOW,
  MEDIUM,
  HIGH
} verbosity_t;

/*DEFINITIONS OF STRUCTURES*/

/** Contains all of the global options used by the program. 
 * Only one instance of this structure exists, and it is passed to most of the other functions and such */
typedef struct _global_options_struct
{
   int mem_allocated;                 /**< Updated with number of bytes allocated for pointer inside this structure */
   
   /* Command line switches */
   verbosity_t VERBOSITY;             /**< How verbose to be - low, medium, or high */
   char *job_control_filename;        /**< Contains the path and filename of the control file */
   char *prmtop_list;                 /**< Path and filename of file containing list of prmtop files */
   char *mdcrd_list;                  /**< Path and filename of file containing list of mdcrd files */
   int num_prmtops;                   /**< How many prmtops there are */

   /* General options */
   char *PARAMETER_FILE_NAME;         /**< Filename with parameters to be fit */
   char *WRITE_ENERGY;                /**< Filename, if any, to save final qm and md energies of structures to */
   short SORT_ENERGY;                 /**< Whether or not mdcrd should be sorted in order of energy before writing */
   char *WRITE_FRCMOD;                /**< Filename, if any, to save ffrcmod to */
   int RANDOM_SEED;                   /**< for duplicating runs if necessary, for debugging usually */
   runtype_t RUNTYPE;
   qm_format_t QMFILEFORMAT;
   coord_format_t COORD_FORMAT;
   energy_t QM_ENERGY_UNITS;          /**< Unit of energy to expect from QM data file */
   force_t QM_FORCE_UNITS;            /**< Unit of force from QM data file*/
   parameter_mode_t PARAMETERS_TO_FIT;
   algorithm_t ALGORITHM;             /**< Fitting routine to be used */
   function_t FUNC_TO_FIT;            /**< The function to be used to fit to the energy surface */
   bool_t K_FIT;                      /**< whether or not to fit the K parameter*/
   int function_calls;                /**< Number of times the fitness function has been called

   /* Things to be fit */
   int BOND_PARAMS;                   /**< stores number of parameters of each type to be fit*/
   int ANGLE_PARAMS;
   int DIHEDRAL_PARAMS;
   int NDIHEDRALS;                    /**< number of terms to give fitted dihedrals, at minimum*/
   
   /* Genetic algorithm options */
   int NOPTIMIZATIONS;                /**< number of parameter sets per generation */
   int MAX_GENERATIONS;               /**< maximum number of generations to run */
   int GENERATIONS_TO_CONVERGE;       /**< number generations in a row without a change to end algorithm */
   int GENERATIONS_TO_SIMPLEX;        /**< number generations in a row without a change to trigger simplex refinement */
   int GENERATIONS_WITHOUT_SIMPLEX;   /**< minimum number of generations between simplex refinements, allows GA time to incorporate new results */
   double SEARCH_SPACE;               /**< distance away from initial parameter set to search */
   double MUTATION_RATE;              /**< percentage of values that may be mutated in a given generation */
   double PARENT_PERCENT;             /**< percentage of values that are allowed to enter recombination */
   
   /* Simplex algorithm options */
   double CONV_LIMIT;                 /**< convergence limit */
   double BONDFC_dx;                  /**< simplex step sizes */
   double BONDEQ_dx;                  /**< simplex step sizes */
   double ANGLEFC_dx;                 /**< simplex step sizes */
   double ANGLEEQ_dx;                 /**< simplex step sizes */
   double DIHEDRALBH_dx;              /**< simplex step sizes. Dihedral BH in prmtop is actually Vn/2 */
   double DIHEDRALN_dx;               /**< simplex step sizes */
   double DIHEDRALG_dx;               /**< simplex step sizes */
   double K_dx;                       /**< simplex step sizes */
   
   /* Bounds checking options */
   bool_t CHECK_BOUNDS;               /**< whether or not to check bounds */
   double ANGLE_LIMIT;                /**< converged angle theta must be this close to values in input structures */
   double BOND_LIMIT;                 /**< converged bond length must be this close to values in input structures */
   int DIHEDRAL_SPAN;                 /**< each dihedral must be spanned by this many input structures */
   bool_t SCATTERPLOTS;               /**< whether or not to write scatter plots with input and output equilibrium parameters */
   
   /* Molecule information options */
   int TOTAL_STRUCTURES;              /**< number of input structures over all molecules */
   double SCNB;                       /**< 1-4 scaling factors for use with standard amber force field equation */
   double SCEE; 

   /* Quantum input file creation options */
   char *QMFILEOUTSTART;              /**< filename for QM input files- format is startNNNend where NNN is structure number */
   char *QMFILEOUTEND;
   char *QMHEADER;                    /**< stores the location of a file to go at the beginning of qm input */
   int QM_SYSTEM_CHARGE;              /**< integral charge of the system */
   int QM_SYSTEM_MULTIPLICITY;        /**< integral multiplicity of the system */
   
   /* Dihedral least squares fit options */
   bool_t FIT_PHASE;                 /**< Whether or not to fit dihedral phases in dihedral least squares function */
   
} global_options_struct;

/** Contains data relating to one input conformation.
 * One of these structures is created for each input conformation, and the array of these structures
 * is collected in coords_struct *coords_data passed to most functions. */
typedef struct _coords_struct
{
  int mem_allocated;                 /**< Updated for amount of memory allocated per structure */
  double *x_coord;                   /**< Array of x-coordinates for each atom in the structure */
  double *y_coord;                   /**< Array of y-coordinates for each atom in the structure */
  double *z_coord;                   /**< Array of z-coordinates for each atom in the structure */
  double energy;                     /**< QM energy for this coordinate set in Kcal/mol */
  double init_energy;                /**< Amber energy calculated with initial parameters */
  double weight;                     /**< The relative importance of this structure in the data set */
  force_struct *force;               /**< Contains forces on each atom if fitting forces */
   
} coords_struct;

/** Contains all structures for one input mdcrd.
 * This is useful in multiprmtop fitting as having a 2D array of coord_structs introduces a lot
 * more opportunities for bugs. This also helps keep track of the filename.
 */
typedef struct _coord_set
{
  coords_struct *struc;              /**< Array containing one coordinate structure for each input conformation */
  char *filename;                    /**< Name of the mdcrd this was read from */
  char *energy_filename;             /**< Name of the directory or file where qm output files/energies are */
  int num_coords;                    /**< The number of structures in the array */
  int natoms;                        /**< The number of atoms in each structure */
  int mem_allocated;                 /**< Tracks amount of memory allocated in this structure */
} coord_set;

/** Contains data used by the bounds checking functions.
 * Only one of these is created. */
typedef struct _bounds_struct
{
  // these hold the values at which bond length, angles, and dihedrals are defined
  double **bond_lengths;             /**< 2D array of bonds x all bond lengths found in the input conformations of the given bond */
  double **angle_thetas;             /**< 2D array of angles x all angle thetas found in the input conformations of the given angle */
  double **dihedral_thetas;          /**< 2D array of dihedrals x all dihedral phases found in the input conformations of a given dihedral */
  int mem_allocated;                 /**< How much memory is being used */
} bounds_struct;

/*Function defs*/
void print_program_info(void);
void print_program_history(void);
int set_default_options(global_options_struct *global_options);
void process_retval(int err_code, verbosity_t VERBOSITY);
void malloc_failure_char(char *routine, char *var_name, int chars_requested);
void malloc_failure_int(char *routine, char *var_name, int ints_requested);
void malloc_failure_short_int(char *routine, char *var_name, int short_ints_requested);
void malloc_failure_double(char *routine, char *var_name, int doubles_requested);
void file_open_failure(char *routine, char *var_name);
double **alloc_2D_double(int nrows, int ncolumns);
void global_unlock(global_options_struct *global_options, parm_struct **parm_data, coord_set **coords_data);
void print_close_line_box(int no_spaces);
void print_open_line_box(int *i);
int check_for_valid_filename(const char *data_string, const int length);
int s_getline(char *line, int max, FILE *fp);
int find_flag( FILE *fptr, char *label );
int name_copy( FILE *fptr, char *stringp);
int find_atomic_number_from_parm(parm_struct *parm_data, int atom);
void print_atomic_number_as_symbol(FILE *fptr, int atomic_number);
double calc_r_squared_multiprmtop(global_options_struct *global_options, parm_struct *parm_datas, coord_set *coords_data);
int unObfuscateAtom( int at );
int ObfuscateAtom( int at ); 
void print_parameter_summary(global_options_struct *global_options, parm_struct *parm_data);
double calc_bond_length(double bond1x, double bond1y, double bond1z, double bond2x, double bond2y, double bond2z);
double calc_angle_radians(double atom1x, double atom1y, double atom1z, double atom2x, double atom2y, double atom2z,
                          double atom3x, double atom3y, double atom3z);
double calc_dihedral_radians(double atom1x, double atom1y, double atom1z, double atom2x, double atom2y, double atom2z,
                          double atom3x, double atom3y, double atom3z, double atom4x, double atom4y, double atom4z);
void calc_fit_dimensions(global_options_struct *global_options, parm_struct *parm_data);
int modify_params_scratch_data(global_options_struct *global_options, parm_struct *parm_data, double *parameters, readwrite_t MODE);
int calculate_no_fit_params(parm_struct *parm_data, short int MODE);
void double_2D_array_free(double **array);

int write_input_parameters(global_options_struct *global_options, parm_struct *parm_data);
int read_input_parameters(global_options_struct *global_options, parm_struct *parm_data);
int read_parameter_file(global_options_struct *global_options, parm_struct *parm_data);
int write_frcmod(global_options_struct *global_options, parm_struct *parm_data);
void handle_sigint(int param);
void print_backtrace(int signal);
int dihedral_types_equal(dihedral_type_struct *first, dihedral_type_struct *second);
void print_dihedral(dihedral_type_struct *type, int term);
int compare_energy(const void *a, const void *b); // this is for qsort.
int not_enough_dihedrals(parm_struct *parm_data, int n);
int do_fit(global_options_struct *global_options, parm_struct *parm_datas, coord_set *coords_datas);
int minimise_function_simplex(global_options_struct *global_options, parm_struct *parm_datas, coord_set *coords_data);  
int check_range(global_options_struct *global_options, parm_struct *parm_data);

/* Bounds checking functions */
int calculate_structure_diversity(global_options_struct *global_options, bounds_struct *bounds_data, parm_struct *parm_data, coord_set *coords_data);
int check_bonds(global_options_struct *global_options, parm_struct *parm_data, coord_set *coords_data, bounds_struct *bounds_data);
int check_angles(global_options_struct *global_options, parm_struct *parm_data, coord_set *coords_data, bounds_struct *bounds_data);
int check_dihedrals(global_options_struct *global_options, parm_struct *parm_data, coord_set *coords_data, bounds_struct *bounds_data);
void clean_up_bounds(bounds_struct *bounds_data, parm_struct *parm_data);

/* Job control and command line parameter functions */
int process_job_control_setting(char *setting_line, int length, int *number_settings, global_options_struct *global_options, parm_struct *parm_datas, coord_set *coords_data);
int read_job_control_file(global_options_struct *global_options, parm_struct *parm_data, coord_set *coords_data);
void print_job_control_summary(global_options_struct *global_options, parm_struct *parm_datas, coord_set *coords_datas);
int process_command_line(int argc, char **argv, global_options_struct *global_options, parm_struct **parm_datas, coord_set **coords_datas);
void command_line_help();

/* QM i/o functions */
int read_qm(global_options_struct *global_options, parm_struct *parm_datas, coord_set *coords_datas);
int read_qm_directory(global_options_struct *global_options, parm_struct *parm_data, coord_set *coords_data);
int read_qm_energy_list(global_options_struct *global_options, parm_struct *parm_data, coord_set *coords_data);
int read_gaussian_forces(global_options_struct *global_options, parm_struct *parm_data, coord_set *coords_data);
int read_gaussian_energy(global_options_struct *global_options, parm_struct *parm_data, coord_set *coords_data);

int write_input_adf(global_options_struct *global_options, parm_struct *parm_data, coord_set *current_struct, int num, FILE *fptr);
int write_input_gamess(global_options_struct *global_options, parm_struct *parm_data, coord_set *current_struct, int num, FILE *fptr);
int write_input_gaussian(global_options_struct *global_options, parm_struct *parm_data, coord_set *current_struct, int num, FILE *fptr);

/* Output producing functions */
int create_qm_input(global_options_struct *global_options, parm_struct *parm_datas, coord_set *coords_datas);
int create_input_single_prmtop(global_options_struct *global_options, parm_struct *parm_data, coord_set *coords_data);
int write_energy(global_options_struct *global_options, parm_struct *parm_data, coord_set *coords_data);

/* Prmtop functions */
void set_dihedral_fit(global_options_struct *global_options, dihedral_type_struct *s, int t,
                     bool_t d_kp, bool_t d_np, bool_t d_phase);
int read_prmtops(global_options_struct *global_options, parm_struct **parm_datas);
int read_single_prmtop(global_options_struct *global_options, parm_struct *parm_data);
int process_prmtops(global_options_struct *global_options, parm_struct *parm_datas);
int process_single_prmtop(global_options_struct *global_options, parm_struct *parm_data);
int verify_prmtops(global_options_struct *global_options, parm_struct *parm_datas);
int bondcomparator(const void *a, const void *b);
int anglecomparator(const void *a, const void *b);
int dihedralcomparator(const void *a, const void *b);
void print_multiprmtop_summary(global_options_struct *global_options, parm_struct *parm_datas);
void free_prmtop(global_options_struct *global_options, parm_struct *parm_data);
int update_prmtop_data(global_options_struct *global_options, parm_struct *parm_datas, double *parameters);
int update_data_matrix(global_options_struct *global_options, parm_struct *parm_datas, double *parameters);

/* Mdcrd functions */
int read_mdcrds(global_options_struct *global_options, parm_struct *parm_datas, coord_set **s_datas);
int read_single_mdcrd(global_options_struct *global_options, coord_set *coords_data);
int alloc_coords(global_options_struct *global_options, coord_set *c);
int free_coords(global_options_struct *global_options, coord_set *c);

/* Energy evaluation functions */
double eval_amber_std_for_single_struct(global_options_struct *global_options, parm_struct *parm_data, coords_struct *coords_datas);
double eval_sum_squares_amber_std(global_options_struct *global_options, parm_struct *parm_data, coord_set *coords_data, int *range_altered);
double eval_sum_squares_amber_std_multiprmtop(global_options_struct *global_options, parm_struct *parm_datas, coord_set *coords_datas, int *range_altered);

/* Forces evaluation functions*/
int eval_amber_forces_single_struct(global_options_struct *global_options, parm_struct *parm_data, coord_set *coords_datas, force_struct *forces, int structure);
double eval_sum_amber_forces(global_options_struct *global_options, parm_struct *parm_data, coord_set *coords_data); 
void print_forces(global_options_struct *global_options, parm_struct *parm_data, coord_set *coords_data);
double eval_sum_amber_forces_multiprmtop(global_options_struct *global_options, parm_struct *parm_datas, coord_set *coords_datas);

/* Genetic algorithm functions */
int minimise_function_genetic(global_options_struct *global_options, parm_struct *parm_data, coord_set *coords_data);
int do_mutation(global_options_struct *global_options, parm_struct *parm_data, double *row, int col, bool_t do_mutate);
double **alloc_data_matrix(int rows, int cols);
void free_data_matrix(double **dm, int rows);

/* Wizard functions */
int job_control_wizard(global_options_struct *global_options, parm_struct *parm_datas) ;
int get_option(int min, int max);
double get_float();
void genetic_wizard(global_options_struct *global_options, FILE *file);
void simplex_wizard(global_options_struct *global_options, FILE *file);

/* New dihedral fittings */
int dihedral_least_squares(global_options_struct *global_options, parm_struct *parm_data, coord_set *coords_data, bool_t fit_phase);
int conduct_dihedral_least_squares(global_options_struct *global_options, parm_struct *parm_data, coord_set *coords_data, 
                                       bool_t fit_phase);
