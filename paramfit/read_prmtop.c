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

/** @file read_prmtop.c
  The routines here are responsible for opening, reading
  and processing the prmtop file.

  Note, we allocate memory and read all of the options in the
  prmtop file even though a number, such as hydrogen bonding
  are not actually used in the fitting. The reason for this
  is that it means other parts of the program can be updated
  at a later date without worrying if the data is actually read
  from the prmtop file.
  */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "function_def.h"

/**
 * Reads in all prmtop data.
 * Starts by getting filenames one line at a time from the prmtop list and then
 * reading in data from each one of those.
 * @see read_single_prmtop
 * @param[in] global_options The global options structure, with the parm list filename
 * @param[in,out] parm_data Pointer to an array of parm_structs, will be reallocd.
 * @return Integer indicating success or failure
 */
int read_prmtops(global_options_struct *global_options, parm_struct **parm_datas)
{
  int retval;
  FILE *fptr;

  // Check a prmtop is provided in command line. If not, prompt for one
  if (!global_options->prmtop_list && !*parm_datas) {

    printf("*** ERROR - must specify prmtop file with -p or prmtop list with -pf! *** \n");
    printf("            Re-run with appropriate arguments.\n");
    return UNKNOWN_OPT;
  }

  // Open the file list or read in single prmtop
  if (global_options->prmtop_list) {
    if (global_options->VERBOSITY>=MEDIUM)
      printf(" Reading prmtop file list: %s\n",global_options->prmtop_list);
    if((fptr=fopen(global_options->prmtop_list,"r"))==NULL) {
      file_open_failure("read_prmtops", global_options->prmtop_list);
      return FILE_OPEN_FAIL;
    }
  } else {
    global_options->num_prmtops=1;
    retval = read_single_prmtop(global_options, &(*parm_datas)[0]);
    process_retval(retval,global_options->VERBOSITY);
    return SUCCESS;
  }

  // Allocate memory for one line and initial parm data
  char* line = (char*)malloc(BUFFER_SIZE);
  if (line==NULL) {
    malloc_failure_char("read_prmtops","line",BUFFER_SIZE);
    return ALLOC_FAIL;
  }

  // Read in list of files, one line at a time
  global_options->num_prmtops=0;
  int files_read=0;
  while (s_getline( line, BUFFER_SIZE, fptr) != FILE_READ_FAIL) {
    // Allocate memory for new parm struct (initially parm_datas=NULL so realloc acts like malloc)
    ++global_options->num_prmtops;
    *parm_datas=(parm_struct*)realloc(*parm_datas,global_options->num_prmtops*sizeof(parm_struct));
    if (parm_datas==NULL) {
      printf("ERROR allocating %lu bytes in read_prmtops for parm_datas\n", global_options->num_prmtops*sizeof(parm_struct));
      return ALLOC_FAIL;
    }
    // Set the filename of this structure
    // No +1 because we will remove trailing newline character
    (*parm_datas)[global_options->num_prmtops-1].filename= (char*)malloc(strlen(line));
    if ((*parm_datas)[global_options->num_prmtops-1].filename==NULL) {
      malloc_failure_char("read_prmtops","parm_data->filename",strlen(line));
      return ALLOC_FAIL;
    }
    // If K is not being fit, values for K for each prmtop should be provided in the prmtop list.
    (*parm_datas)[global_options->num_prmtops-1].K = 0.0;
    if (global_options->K_FIT==NO) {
      if ( sscanf(line,"%s %lf", (*parm_datas)[global_options->num_prmtops-1].filename, &((*parm_datas)[global_options->num_prmtops-1].K) ) !=2) {
        printf("ERROR: malformed line in prmtop list %s line %d\nLine was: %s", global_options->prmtop_list,
            files_read, line);
        free((*parm_datas)[global_options->num_prmtops-1].filename);
        return INVALID_FORMAT;
      }
    } else if (global_options->K_FIT==YES) {
      printf("ERROR: cannot fit K with multiple prmtops.\n");
      printf("       Run paramfit in single-prmtop mode and fit K for each prmtop individually.\n");
      printf("       Write the value of K for each prmtop following the filename in the prmtop list.\n");
      printf("       Example: %s 9000.0\n", (*parm_datas)[global_options->num_prmtops-1].filename);
      return INVALID_FORMAT;
    }
    //strncpy((*parm_datas)[global_options->num_prmtops-1].filename, line, strlen(line)-1);
    (*parm_datas)[global_options->num_prmtops-1].filename[strlen(line)-1]='\0';

    // Read in data into new prmtop
    retval = read_single_prmtop(global_options, &(*parm_datas)[global_options->num_prmtops-1]);
    process_retval(retval,global_options->VERBOSITY);
    ++files_read;
  }

  // Clean up
  fclose(fptr);
  free(line);
  return SUCCESS;
}

/**
 * Reads the raw data from the prmtop file.
 * Sorts out how many atoms, atom types, and basic info. Does not do any processing
 * into the fitting data structure.
 * @see process_prmtop.c
 * 
 * @param[in] global_options The global options structure
 * @param[in,out] parm_data The parameter structure that will have info put in
 * @param[in] parm_data->filename The path to the prmtop file to read
 */
int read_single_prmtop(global_options_struct *global_options, parm_struct *parm_data)
{
  FILE *fptr;
  fpos_t pos;
  char *prmtop_data;  /*buffer to store the prmtop in memory for quicker reading*/
  int linelen;
  int retval; /*Process return values*/
  int skip_section;
  int i,j;
  int len;

  parm_data->AMBER_SYSTEM_CHARGE=0.0;
  parm_data->mem_allocated = sizeof(parm_struct);

  /*I will allocate BUFFER_SIZE for our prmtop buffer - no prmtop line should
    be longer than this*/
  if (global_options->VERBOSITY>=HIGH)
    printf("   Allocating %d bytes for *prmtop_data\n",BUFFER_SIZE);
  prmtop_data = (char *) malloc(BUFFER_SIZE);
  if (prmtop_data == NULL)
  {
    malloc_failure_char("read_prmtop","prmtop_data",BUFFER_SIZE);
    return ALLOC_FAIL;
  }

  /*Step 1 - check if we can open the file*/
  /*We read a line at a time from the prmtop file*/
  if (global_options->VERBOSITY>=MEDIUM)
    printf("  Reading prmtop file    : %s\n",parm_data->filename);

  if((fptr=fopen(parm_data->filename,"r"))==NULL)
  {
    file_open_failure("read_prmtop", parm_data->filename);
    return FILE_OPEN_FAIL;
  }

  retval=s_getline( prmtop_data, BUFFER_SIZE, fptr);  /*%VERSION???*/
  if ( retval==FILE_READ_FAIL )
  {
    printf("*** ERROR IN readparm SUBROUTINE\n");
    printf("*** HIT EOF WHILE DETERMINING VERSION OF PRMTOP FILE: %s\n",parm_data->filename);
    return FILE_READ_FAIL;
  }


  /*First things first, see if this is a "new" >v7.0 format prmtop file
    and read the title*/
  if ( !strncmp( prmtop_data, "%VERSION", 8 ) )
  {
    parm_data->newparm=YES;
    retval=s_getline( prmtop_data, BUFFER_SIZE, fptr);  /*%FLAG TITLE*/
    if ( retval==FILE_READ_FAIL )
    {
      printf("*** ERROR IN readparm SUBROUTINE\n");
      printf("*** HIT EOF WHILE TRYING TO READ TITLE CARD IN PRMTOP FILE: %s\n",parm_data->filename);
      return FILE_READ_FAIL;
    }
    if ( strncmp( prmtop_data+6, "TITLE", 5 ) && strncmp( prmtop_data+6, "CTITLE", 6) )
    {
      printf("*** ERROR IN readparm SUBROUTINE\n");
      printf("*** FAILED TO FIND TITLE CARD IN PRMTOP FILE: %s\n",parm_data->filename);
      return INVALID_FORMAT;
    }
    // Skip format or other comment line, get the title
    prmtop_data[0]='%';
    while (strncmp(prmtop_data,"%",1)) {
      if (fgets (prmtop_data,BUFFER_SIZE,fptr) != NULL) {
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** HIT EOF WHILE TRYING TO READ TITLE IN PRMTOP FILE: %s\n",parm_data->filename);
        return FILE_READ_FAIL;
      }
    }
    linelen=strlen(prmtop_data);
    /*allocate the memory for the parm title card*/
    if (global_options->VERBOSITY>=HIGH)
      printf("   Allocating %d bytes for parm_data->*title\n",(linelen+1));
    parm_data->title = (char *) malloc(linelen+1);
    if (parm_data->title == NULL)
    {
      malloc_failure_char("read_prmtop", "parm_data->title", (linelen+1));
      return ALLOC_FAIL;
    }
    parm_data->mem_allocated+=(linelen+1);
    strncpy(parm_data->title, prmtop_data, linelen);
    /*add null character so we can handle it as a string*/
    parm_data->title[linelen]='\0';

    if (global_options->VERBOSITY>=HIGH)
    {
      printf("   Prmtop   (format): This is a NEW format prmtop file.\n");
      printf("   Prmtop    (title): %s",parm_data->title);
    }
  }
  else
  {
    parm_data->newparm=NO; /*Old format so no %FORMAT(xxxx) cards*/
    /*allocate the memory for the parm title card*/
    if (global_options->VERBOSITY>=HIGH)
      printf("   Allocating %d bytes for parm_data->*title\n",(PRMTOP_TITLE_LENGTH+1));
    parm_data->title = (char *) malloc(PRMTOP_TITLE_LENGTH+1);
    if (parm_data->title == NULL)
    {
      malloc_failure_char("read_prmtop", "parm_data->title", (PRMTOP_TITLE_LENGTH+1));
      return ALLOC_FAIL;
    }
    parm_data->mem_allocated+=(PRMTOP_TITLE_LENGTH+1);
    strncpy(parm_data->title, prmtop_data, PRMTOP_TITLE_LENGTH);
    /*add null character so we can handle it as a string*/
    parm_data->title[PRMTOP_TITLE_LENGTH]='\0';

    if (global_options->VERBOSITY>=MEDIUM)
    {
      printf("   Prmtop   (format): This is an OLD format prmtop file.\n");
      printf("   Prmtop    (title): %s",parm_data->title);
    }
  }

  /*Next step is to read in the integer control variables*/
  if ( parm_data->newparm )
  {
    /*Test if we can find the POINTERS flag*/
    retval=find_flag( fptr, "POINTERS" );

    if ( retval!=SUCCESS ) {
      printf("*** ERROR IN readparm SUBROUTINE\n");
      printf("*** FAILED TO FIND POINTERS CARD IN PRMTOP FILE: %s\n",parm_data->filename);
      return INVALID_FORMAT;
    }
    // Skip comment lines
    prmtop_data[0]='%';
    while (prmtop_data[0]=='%') {
      fgetpos(fptr, &pos); // Get position of beginning of line in case this is a good line
      if (fgets (prmtop_data,BUFFER_SIZE,fptr) == NULL) {
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** HIT EOF WHILE TRYING TO READ FORMAT CARD IN PRMTOP FILE: %s\n",parm_data->filename);
        return FILE_READ_FAIL;
      }
    }
    fsetpos(fptr, &pos); // Go to beginning of good line
  }
  /*New prmtop file we should read 31 integers, old prmtop = 30 integers*/
  fscanf(fptr,"%d",&parm_data->NTOTAT);
  fscanf(fptr,"%d",&parm_data->NTYPES);
  // NATOMTYPES = parm_data->NTYPES;

  fscanf(fptr,"%d",&parm_data->NBONH);  /*Number of bonds involving Hydrogen*/
  fscanf(fptr,"%d",&parm_data->NBONA);  /*Number of bonds not involving Hydrogen*/
  // NBONDS = (parm_data->NBONH + parm_data->NBONA);

  fscanf(fptr,"%d",&parm_data->NTHETH);  /*Number of angles involving Hydrogen*/
  fscanf(fptr,"%d",&parm_data->NTHETA);  /*Number of angles not involving Hydrogen*/
  // NANGLES = (parm_data->NTHETH + parm_data->NTHETA);

  fscanf(fptr,"%d",&parm_data->NPHIH);  /*Number of dihedrals involving Hydrogen*/
  fscanf(fptr,"%d",&parm_data->NPHIA);  /*Number of dihedrals not involving Hydrogen*/
  // NDIHEDRALS = (parm_data->NPHIH + parm_data->NPHIA); 

  fscanf(fptr,"%d",&parm_data->JHPARM);
  fscanf(fptr,"%d",&parm_data->JPARM);
  fscanf(fptr,"%d",&parm_data->NEXT);
  fscanf(fptr,"%d",&parm_data->NTOTRS);
  fscanf(fptr,"%d",&parm_data->MBONA);
  fscanf(fptr,"%d",&parm_data->MTHETS);
  fscanf(fptr,"%d",&parm_data->MPHIA);

  fscanf(fptr,"%d",&parm_data->MUMBND);  /*Number of unique bond types*/
  // NBONDTYPES = (parm_data->MUMBND);

  fscanf(fptr,"%d",&parm_data->MUMANG);  /*Number of unique angle types*/
  // NANGLETYPES = (parm_data->MUMANG);

  fscanf(fptr,"%d",&parm_data->MPTRA);  /*Number of unique dihedral types*/
  // NDIHEDRALTYPES = (parm_data->MPTRA);

  fscanf(fptr,"%d",&parm_data->NATYP);
  fscanf(fptr,"%d",&parm_data->NHB);
  fscanf(fptr,"%d",&parm_data->IFPERT);
  fscanf(fptr,"%d",&parm_data->NBPER);
  fscanf(fptr,"%d",&parm_data->NGPER);
  fscanf(fptr,"%d",&parm_data->NDPER);
  fscanf(fptr,"%d",&parm_data->MBPER);
  fscanf(fptr,"%d",&parm_data->MGPER);
  fscanf(fptr,"%d",&parm_data->MDPER);
  fscanf(fptr,"%d",&parm_data->IFBOX);
  fscanf(fptr,"%d",&parm_data->NMXRS);
  fscanf(fptr,"%d",&parm_data->IFCAP);
  fscanf(fptr,"%d",&parm_data->NUMEXTRA);

  if (global_options->VERBOSITY>=HIGH) {
    printf("   Prmtop (pointers): Read in pointer control variables.\n");
  }
  if (global_options->VERBOSITY>=HIGH) {
    printf("   Prmtop (pointers): NTOTAT NTYPES  NBONH  NBONA NTHETH NTHETA  NPHIH  NPHIA\n");
    printf("   Prmtop (pointers): %6d %6d %6d %6d %6d %6d %6d %6d\n"
        ,parm_data->NTOTAT,parm_data->NTYPES,parm_data->NBONH,parm_data->NBONA,parm_data->NTHETH
        ,parm_data->NTHETA,parm_data->NPHIH,parm_data->NPHIA);
    printf("   Prmtop (pointers): JHPARM  JPARM   NEXT NTOTRS  MBONA MTHETS  MPHIA MUMBND\n");
    printf("   Prmtop (pointers): %6d %6d %6d %6d %6d %6d %6d %6d\n"
        ,parm_data->JHPARM,parm_data->JPARM,parm_data->NEXT,parm_data->NTOTRS,parm_data->NBONA
        ,parm_data->MTHETS,parm_data->MPHIA,parm_data->MUMBND);
    printf("   Prmtop (pointers): MUMANG  MPTRA  NATYP    NHB IFPERT  NBPER  NGPER  NDPER\n");
    printf("   Prmtop (pointers): %6d %6d %6d %6d %6d %6d %6d %6d\n"
        ,parm_data->MUMANG,parm_data->MPTRA,parm_data->NATYP,parm_data->NHB,parm_data->IFPERT
        ,parm_data->NBPER,parm_data->NGPER,parm_data->NDPER);
    printf("   Prmtop (pointers):  MBPER  MGPER  MDPER  IFBOX  NMXRS  IFCAP NUMEXTRA\n");
    printf("   Prmtop (pointers): %6d %6d %6d %6d %6d %6d %6d\n"
        ,parm_data->MBPER,parm_data->MGPER,parm_data->MDPER,parm_data->IFBOX,parm_data->NMXRS
        ,parm_data->IFCAP,parm_data->NUMEXTRA);
  }

  /*Next step is to read in the atom names*/
  if ( parm_data->NTOTAT > 0 ) {
    if ( parm_data->newparm ) {
      /*Test if we can find the ATOM_NAME flag*/
      retval=find_flag( fptr, "ATOM_NAME" );
      if ( retval!=SUCCESS ) {
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** FAILED TO FIND ATOM_NAME CARD IN PRMTOP FILE: %s\n",parm_data->filename);
        return INVALID_FORMAT;
      }
      // Skip comment lines
      prmtop_data[0]='%';
      while (prmtop_data[0]=='%') {
        fgetpos(fptr, &pos); // Get position of beginning of line in case this is line with data
        if (fgets (prmtop_data,BUFFER_SIZE,fptr) == NULL) {
          printf("*** ERROR IN readparm SUBROUTINE\n");
          printf("*** HIT EOF WHILE TRYING TO READ ATOM_NAME COMMENTS IN PRMTOP FILE: %s\n",parm_data->filename);
          return FILE_READ_FAIL;
        }
      }
      fsetpos(fptr, &pos); // Go to beginning of the line

    }
    /*Allocate the memory to store the atom info*/
    if (global_options->VERBOSITY>=HIGH)
      printf("   Allocating %d bytes for parm_data->*atom\n",(int)(parm_data->NTOTAT*sizeof( atom_struct )));
    parm_data->atom = (atom_struct *) malloc(parm_data->NTOTAT*sizeof( atom_struct ));
    if (parm_data->atom == NULL) {
      malloc_failure_char("read_prmtop", "parm_data->atom", (parm_data->NTOTAT*sizeof( atom_struct )));
      return ALLOC_FAIL;
    }
    parm_data->mem_allocated+=(parm_data->NTOTAT*sizeof( atom_struct ));
    for (i=0; i < parm_data->NTOTAT; ++i) {
      retval=name_copy(fptr, parm_data->atom[i].igraph);
      if ( retval != SUCCESS) {
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** READ INVALID ATOM NAME FROM PRMTOP FILE: %s\n",parm_data->filename);
        return INVALID_FORMAT;
      }
    }
    if (global_options->VERBOSITY>=HIGH) {
      printf("   Prmtop    (atoms): Read in atom names.\n");
    }
  }
  /*Next step, read in the charges*/
  if ( parm_data->NTOTAT > 0 ) {
    if ( parm_data->newparm ) {
      /*Test if we can find the CHARGE flag*/
      retval=find_flag( fptr, "CHARGE" );

      if ( retval!=SUCCESS ) {
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** FAILED TO FIND CHARGE CARD IN PRMTOP FILE: %s\n",parm_data->filename);
        return INVALID_FORMAT;
      }
      // Skip comment lines
      prmtop_data[0]='%';
      while (prmtop_data[0]=='%') {
        fgetpos(fptr, &pos); // Get position of beginning of line in case this is a good line
        if (fgets (prmtop_data,BUFFER_SIZE,fptr) == NULL) {
          printf("*** ERROR IN readparm SUBROUTINE\n");
          printf("*** HIT EOF WHILE TRYING TO READ FORMAT CARD IN PRMTOP FILE: %s\n",parm_data->filename);
          return FILE_READ_FAIL;
        }
      }
      fsetpos(fptr, &pos); // Go to beginning of good line

    }

    for (i=0;i < parm_data->NTOTAT; ++i) {
      /*fscanf returns number of arguments read so should return 1*/
      if ( fscanf(fptr,"%lf",&parm_data->atom[i].chrg)!=1 )
      {
        /*Error failed to find atom charge - Issue a warning and set to zero*/
        printf("!  WARNING - Failed to find charge for atom %d - assuming charge is zero.\n",i);
        parm_data->atom[i].chrg=0.0;
      }
      /*update the amber total charge*/
      parm_data->AMBER_SYSTEM_CHARGE+=parm_data->atom[i].chrg;
    }

    if (global_options->VERBOSITY>=HIGH) {
      printf("   Prmtop    (atoms): Read in atom charges.\n");
    }
  }
  /*Next step, read in the masses*/
  if ( parm_data->NTOTAT > 0 ) {
    if ( parm_data->newparm ) {
      /*Test if we can find the MASS flag*/
      retval=find_flag( fptr, "MASS" );

      if ( retval!=SUCCESS ) {
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** FAILED TO FIND MASS CARD IN PRMTOP FILE: %s\n",parm_data->filename);
        return INVALID_FORMAT;
      }
      // Skip comment lines
      prmtop_data[0]='%';
      while (prmtop_data[0]=='%') {
        fgetpos(fptr, &pos); // Get position of beginning of line in case this is a good line
        if (fgets (prmtop_data,BUFFER_SIZE,fptr) == NULL) {
          printf("*** ERROR IN readparm SUBROUTINE\n");
          printf("*** HIT EOF WHILE TRYING TO READ FORMAT CARD IN PRMTOP FILE: %s\n",parm_data->filename);
          return FILE_READ_FAIL;
        }
      }
      fsetpos(fptr, &pos); // Go to beginning of good line

    }

    for (i=0;i < parm_data->NTOTAT; ++i) {
      /*fscanf returns number of arguments read so should return 1*/
      if ( fscanf(fptr,"%lf",&parm_data->atom[i].amass)!=1 ) {
        /*Error failed to find atom mass - Issue a warning and set to 1.0*/
        printf("!  WARNING - Failed to find mass for atom %d - assuming mass is 1.0.\n",i);
        parm_data->atom[i].amass=1.0;
      }
    }

    if (global_options->VERBOSITY>=HIGH) {
      printf("   Prmtop    (atoms): Read in atom masses.\n");
    }
  }
  /*Next step, read in the atom type array*/
  if ( parm_data->NTOTAT )
  {
    if ( parm_data->newparm )
    {
      /*Test if we can find the ATOM_TYPE_INDEX flag*/
      retval=find_flag( fptr, "ATOM_TYPE_INDEX" );

      if ( retval!=SUCCESS ) {
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** FAILED TO FIND ATOM_TYPE_INDEX CARD IN PRMTOP FILE: %s\n",parm_data->filename);
        return INVALID_FORMAT;
      }
      // Skip comment lines
      prmtop_data[0]='%';
      while (prmtop_data[0]=='%') {
        fgetpos(fptr, &pos); // Get position of beginning of line in case this is a good line
        if (fgets (prmtop_data,BUFFER_SIZE,fptr) == NULL) {
          printf("*** ERROR IN readparm SUBROUTINE\n");
          printf("*** HIT EOF WHILE TRYING TO READ FORMAT CARD IN PRMTOP FILE: %s\n",parm_data->filename);
          return FILE_READ_FAIL;
        }
      }
      fsetpos(fptr, &pos); // Go to beginning of good line
    }

    for (i=0;i < parm_data->NTOTAT; ++i) {
      /*fscanf returns number of arguments read so should return 1*/
      if ( fscanf(fptr,"%d",&parm_data->atom[i].iac)!=1 )
      {
        /*Fatal Error failed to find atom type*/
        printf("*** ERROR - FAILED TO FIND ATOM TYPE FOR ATOM %d\n",i);
        return INVALID_FORMAT;
      }
    }

    if (global_options->VERBOSITY>=HIGH) {
      printf("   Prmtop    (types): Read in atom type index.\n");
    }
  }
  /*Next step, read in the excluded atoms array*/
  if ( parm_data->NTOTAT ) {
    if ( parm_data->newparm ) {
      /*Test if we can find the NUMBER_EXCLUDED_ATOMS flag*/
      retval=find_flag( fptr, "NUMBER_EXCLUDED_ATOMS" );

      if ( retval!=SUCCESS )
      {
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** FAILED TO FIND NUMBER_EXCLUDED_ATOMS CARD IN PRMTOP FILE: %s\n",parm_data->filename);
        return INVALID_FORMAT;
      }
      // Skip comment lines
      prmtop_data[0]='%';
      while (prmtop_data[0]=='%') {
        fgetpos(fptr, &pos); // Get position of beginning of line in case this is a good line
        if (fgets (prmtop_data,BUFFER_SIZE,fptr) == NULL) {
          printf("*** ERROR IN readparm SUBROUTINE\n");
          printf("*** HIT EOF WHILE TRYING TO READ FORMAT CARD IN PRMTOP FILE: %s\n",parm_data->filename);
          return FILE_READ_FAIL;
        }
      }
      fsetpos(fptr, &pos); // Go to beginning of good line
    }

    for (i=0;i < parm_data->NTOTAT; ++i) {
      if ( fscanf(fptr,"%d",&parm_data->atom[i].numex)!=1 ) {
        printf("*** ERROR - FAILED TO FIND NUMBER EXCLUDED FOR ATOM %d\n",i);
        return INVALID_FORMAT;
      }
    }

    if (global_options->VERBOSITY>=HIGH) {
      printf("   Prmtop (excluded): Read in NUMEX (index to excl atom list).\n");
    }

  }
  /*Next step, read in the nonbonded parm index array*/
  if ( parm_data->NTYPES > 0 ) {
    if ( parm_data->newparm ) {
      /*Test if we can find the NONBONDED_PARM_INDEX flag*/
      retval=find_flag( fptr, "NONBONDED_PARM_INDEX" );

      if ( retval!=SUCCESS ) {
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** FAILED TO FIND NONBONDED_PARM_INDEX CARD IN PRMTOP FILE: %s\n",parm_data->filename);
        return INVALID_FORMAT;
      }
      // Skip comment lines
      prmtop_data[0]='%';
      while (prmtop_data[0]=='%') {
        fgetpos(fptr, &pos); // Get position of beginning of line in case this is a good line
        if (fgets (prmtop_data,BUFFER_SIZE,fptr) == NULL) {
          printf("*** ERROR IN readparm SUBROUTINE\n");
          printf("*** HIT EOF WHILE TRYING TO READ FORMAT CARD IN PRMTOP FILE: %s\n",parm_data->filename);
          return FILE_READ_FAIL;
        }
      }
      fsetpos(fptr, &pos); // Go to beginning of good line
    }

    /*Non bonded parm index consists of NTYPES^2*/
    len = (parm_data->NTYPES * parm_data->NTYPES);
    /*allocate memory first of all*/
    if (global_options->VERBOSITY>=HIGH)
      printf("   Allocating %d bytes for parm_data->*nno\n",(int)(len*sizeof(int)));
    parm_data->nno = (int *) malloc(len*sizeof(int));
    if (parm_data->nno == NULL) {
      malloc_failure_char("read_prmtop", "parm_data->nno", (len*sizeof(int)));
      return ALLOC_FAIL;
    }
    parm_data->mem_allocated+=(len*sizeof(int));

    for (i=0;i < len; ++i) {
      if ( fscanf(fptr,"%d",&parm_data->nno[i])!=1 ) {
        /*Fatal Error failed to find atom type*/
        printf("*** ERROR - FAILED TO FIND NONBONDED INDEX FOR ELEMENT %d\n",i);
        return INVALID_FORMAT;
      }
    }

    if (global_options->VERBOSITY>=HIGH) {
      printf("   Prmtop  (nonbond): Read in NNO (index for nonbond of @type).\n");
    }
  }

  /*Next step, read in the residue labels*/
  if ( parm_data->NTOTRS > 0 )
  {
    if ( parm_data->newparm ) {
      /*Test if we can find the RESIDUE_LABEL flag*/
      retval=find_flag( fptr, "RESIDUE_LABEL" );

      if ( retval!=SUCCESS )
      {
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** FAILED TO FIND RESIDUE_LABEL CARD IN PRMTOP FILE: %s\n",parm_data->filename);
        return INVALID_FORMAT;
      }
      // Skip comment lines
      prmtop_data[0]='%';
      while (prmtop_data[0]=='%') {
        fgetpos(fptr, &pos); // Get position of beginning of line in case this is a good line
        if (fgets (prmtop_data,BUFFER_SIZE,fptr) == NULL) {
          printf("*** ERROR IN readparm SUBROUTINE\n");
          printf("*** HIT EOF WHILE TRYING TO READ FORMAT CARD IN PRMTOP FILE: %s\n",parm_data->filename);
          return FILE_READ_FAIL;
        }
      }
      fsetpos(fptr, &pos); // Go to beginning of good line

    }

    /*allocate memory first of all +1 for later*/
    if (global_options->VERBOSITY>=HIGH)
      printf("   Allocating %d bytes for parm_data->*residue\n",(int) ((parm_data->NTOTRS+1)*sizeof(residue)));
    parm_data->residue = (residue *) malloc((parm_data->NTOTRS+1)*sizeof(residue));
    if (parm_data->residue == NULL)
    {
      malloc_failure_char("read_prmtop", "parm_data->residue", (parm_data->NTOTRS+1)*sizeof(residue));
      return ALLOC_FAIL;
    }
    parm_data->mem_allocated+=(parm_data->NTOTRS+1)*sizeof(residue);

    for (i=0; i < parm_data->NTOTRS; ++i)
    {
      retval=name_copy(fptr, parm_data->residue[i].labres);
      if ( retval != SUCCESS) {
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** READ INVALID RESIDUE NAME FROM PRMTOP FILE: %s\n",parm_data->filename);
        return INVALID_FORMAT;
      }
    }
    if (global_options->VERBOSITY>=HIGH)
      printf("   Prmtop (residues): Read in residue labels.\n");
  }

  if ( parm_data->NTOTRS > 0)
  {
    /*Now read in ipres*/
    if ( parm_data->newparm )
    {
      /*Test if we can find the RESIDUE_POINTER flag*/
      retval=find_flag( fptr, "RESIDUE_POINTER" );

      if ( retval!=SUCCESS )
      {
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** FAILED TO FIND RESIDUE_POINTER CARD IN PRMTOP FILE: %s\n",parm_data->filename);
        return INVALID_FORMAT;
      }
      // Skip comment lines
      prmtop_data[0]='%';
      while (prmtop_data[0]=='%') {
        fgetpos(fptr, &pos); // Get position of beginning of line in case this is a good line
        if (fgets (prmtop_data,BUFFER_SIZE,fptr) == NULL) {
          printf("*** ERROR IN readparm SUBROUTINE\n");
          printf("*** HIT EOF WHILE TRYING TO READ FORMAT CARD IN PRMTOP FILE: %s\n",parm_data->filename);
          return FILE_READ_FAIL;
        }
      }
      fsetpos(fptr, &pos); // Go to beginning of good line
    }

    for (i=0;i < parm_data->NTOTRS; ++i)
    {
      /*fscanf returns number of arguments read so should return 1*/
      if ( fscanf(fptr,"%8d",&parm_data->residue[i].ipres)!=1 )
      {
        /*Fatal Error failed to find atom type*/
        printf("*** ERROR - FAILED TO FIND RESIDUE POINTER FOR RESIDUE %d\n",i);
        return INVALID_FORMAT;
      }
    }
    parm_data->residue[parm_data->NTOTRS].ipres = parm_data->NTOTAT + 1;
    /*Now fill in the residues for each atom*/
    for (i=0; i < parm_data->NTOTRS; ++i)
      for (j=parm_data->residue[i].ipres-1; j < parm_data->residue[i+1].ipres - 1; ++j)
        parm_data->atom[j].res = i;
    parm_data->atom[parm_data->NTOTAT-1].res = parm_data->NTOTRS-1;

    if (global_options->VERBOSITY>=HIGH)
      printf("   Prmtop (residues): Read in residue to atom pointer list.\n");
  }

  /*Next step is to read in the bond parameters and allocate memory*/
  if ( parm_data->MUMBND > 0 )
  {
    if ( parm_data->newparm )
    {
      /*Test if we can find the BOND_FORCE_CONSTANT flag*/
      retval=find_flag( fptr, "BOND_FORCE_CONSTANT" );

      if ( retval!=SUCCESS )
      {
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** FAILED TO FIND BOND_FORCE_CONSTANT IN PRMTOP FILE: %s\n",parm_data->filename);
        return INVALID_FORMAT;
      }
      // Skip comment lines
      prmtop_data[0]='%';
      while (prmtop_data[0]=='%') {
        fgetpos(fptr, &pos); // Get position of beginning of line in case this is a good line
        if (fgets (prmtop_data,BUFFER_SIZE,fptr) == NULL) {
          printf("*** ERROR IN readparm SUBROUTINE\n");
          printf("*** HIT EOF WHILE TRYING TO READ FORMAT CARD IN PRMTOP FILE: %s\n",parm_data->filename);
          return FILE_READ_FAIL;
        }
      }
      fsetpos(fptr, &pos); // Go to beginning of good line

    }

    /*Allocate memory for the force constants*/
    if (global_options->VERBOSITY>=HIGH)
      printf("   Allocating %d bytes for parm_data->*rk\n",(int) ((parm_data->MUMBND)*sizeof(double)));
    parm_data->rk = (double *) malloc((parm_data->MUMBND)*sizeof(double));
    if (parm_data->rk== NULL)
    {
      malloc_failure_char("read_prmtop", "parm_data->rk", (parm_data->MUMBND)*sizeof(double));
      return ALLOC_FAIL;
    }
    parm_data->mem_allocated+=(parm_data->MUMBND)*sizeof(double);

    for (i=0; i < parm_data->MUMBND; ++i)
    {
      /*fscanf returns number of arguments read so should return 1*/
      if ( fscanf(fptr,"%lf",&parm_data->rk[i])!=1 )
      {
        /*Error failed to find force constant - Issue a warning and set to 0.0*/
        printf("!  WARNING - Failed to find force constant for bond %d - assuming it is 0.0.\n",i);
        parm_data->rk[i]=0.0;
      }
    }
    if (global_options->VERBOSITY>=HIGH)
      printf("   Prmtop    (bonds): Read in bond force constants RK.\n");

    if ( parm_data->newparm )
    {
      /*Test if we can find the BOND_EQUIL_VALUE flag*/
      retval=find_flag( fptr, "BOND_EQUIL_VALUE" );

      if ( retval!=SUCCESS )
      {
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** FAILED TO FIND BOND_EQUIL_VALUE IN PRMTOP FILE: %s\n",parm_data->filename);
        return INVALID_FORMAT;
      }
      // Skip comment lines
      prmtop_data[0]='%';
      while (prmtop_data[0]=='%') {
        fgetpos(fptr, &pos); // Get position of beginning of line in case this is a good line
        if (fgets (prmtop_data,BUFFER_SIZE,fptr) == NULL) {
          printf("*** ERROR IN readparm SUBROUTINE\n");
          printf("*** HIT EOF WHILE TRYING TO READ FORMAT CARD IN PRMTOP FILE: %s\n",parm_data->filename);
          return FILE_READ_FAIL;
        }
      }
      fsetpos(fptr, &pos); // Go to beginning of good line

    }

    /*Allocate memory for the bond eq constants*/
    if (global_options->VERBOSITY>=HIGH)
      printf("   Allocating %d bytes for parm_data->*req\n",(int)((parm_data->MUMBND)*sizeof(double)));
    parm_data->req = (double *) malloc((parm_data->MUMBND)*sizeof(double));
    if (parm_data->req== NULL)
    {
      malloc_failure_char("read_prmtop", "parm_data->req", (parm_data->MUMBND)*sizeof(double));
      return ALLOC_FAIL;
    }
    parm_data->mem_allocated+=(parm_data->MUMBND)*sizeof(double);
    for (i=0; i < parm_data->MUMBND; ++i)
    {
      /*fscanf returns number of arguments read so should return 1*/
      if ( fscanf(fptr,"%lf",&parm_data->req[i])!=1 )
      {
        /*Error failed to find bond eq constant - Issue a warning and set to 1.0*/
        printf("!  WARNING - Failed to find bond eq constant for bond %d - assuming it is 1.0.\n",i);
        parm_data->req[i]=1.0;
      }
    }
    if (global_options->VERBOSITY>=HIGH)
      printf("   Prmtop    (bonds): Read in bond equilibrium constants REQ.\n");
  }

  /*NOW THE ANGLES*/
  if ( parm_data->MUMANG > 0 )
  {
    if ( parm_data->newparm )
    {
      /*Test if we can find the ANGLE_FORCE_CONSTANT flag*/
      retval=find_flag( fptr, "ANGLE_FORCE_CONSTANT" );

      if ( retval!=SUCCESS )
      {
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** FAILED TO FIND ANGLE_FORCE_CONSTANT IN PRMTOP FILE: %s\n",parm_data->filename);
        return INVALID_FORMAT;
      }
      // Skip comment lines
      prmtop_data[0]='%';
      while (prmtop_data[0]=='%') {
        fgetpos(fptr, &pos); // Get position of beginning of line in case this is a good line
        if (fgets (prmtop_data,BUFFER_SIZE,fptr) == NULL) {
          printf("*** ERROR IN readparm SUBROUTINE\n");
          printf("*** HIT EOF WHILE TRYING TO READ FORMAT CARD IN PRMTOP FILE: %s\n",parm_data->filename);
          return FILE_READ_FAIL;
        }
      }
      fsetpos(fptr, &pos); // Go to beginning of good line
    }

    /*Allocate memory for the force constants*/
    if (global_options->VERBOSITY>=HIGH)
      printf("   Allocating %d bytes for parm_data->*tk\n",(int)((parm_data->MUMANG)*sizeof(double)));
    parm_data->tk = (double *) malloc((parm_data->MUMANG)*sizeof(double));
    if (parm_data->tk== NULL)
    {
      malloc_failure_char("read_prmtop", "parm_data->tk", (parm_data->MUMANG)*sizeof(double));
      return ALLOC_FAIL;
    }
    parm_data->mem_allocated+=(parm_data->MUMANG)*sizeof(double);
    for (i=0; i < parm_data->MUMANG; ++i)
    {
      /*fscanf returns number of arguments read so should return 1*/
      if ( fscanf(fptr,"%lf",&parm_data->tk[i])!=1 )
      {
        /*Error failed to find force constant - Issue a warning and set to 0.0*/
        printf("!  WARNING - Failed to find force constant for angle %d - assuming it is 0.0.\n",i);
        parm_data->tk[i]=0.0;
      }
    }
    if (global_options->VERBOSITY>=HIGH)
      printf("   Prmtop   (angles): Read in angle force constants TK.\n");

    if ( parm_data->newparm )
    {
      /*Test if we can find the ANGLE_EQUIL_VALUE flag*/
      retval=find_flag( fptr, "ANGLE_EQUIL_VALUE" );

      if ( retval!=SUCCESS )
      {
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** FAILED TO FIND ANGLE_EQUIL_VALUE IN PRMTOP FILE: %s\n",parm_data->filename);
        return INVALID_FORMAT;
      }
      // Skip comment lines
      prmtop_data[0]='%';
      while (prmtop_data[0]=='%') {
        fgetpos(fptr, &pos); // Get position of beginning of line in case this is a good line
        if (fgets (prmtop_data,BUFFER_SIZE,fptr) == NULL) {
          printf("*** ERROR IN readparm SUBROUTINE\n");
          printf("*** HIT EOF WHILE TRYING TO READ FORMAT CARD IN PRMTOP FILE: %s\n",parm_data->filename);
          return FILE_READ_FAIL;
        }
      }
      fsetpos(fptr, &pos); // Go to beginning of good line

    }

    /*Allocate memory for the angle eq constants*/
    if (global_options->VERBOSITY>=HIGH)
      printf("   Allocating %d bytes for parm_data->*teq\n",(int)((parm_data->MUMANG)*sizeof(double)));
    parm_data->teq = (double *) malloc((parm_data->MUMANG)*sizeof(double));
    if (parm_data->teq== NULL)
    {
      malloc_failure_char("read_prmtop", "parm_data->teq", (parm_data->MUMANG)*sizeof(double));
      return ALLOC_FAIL;
    }
    parm_data->mem_allocated+=(parm_data->MUMANG)*sizeof(double);
    for (i=0; i < parm_data->MUMANG; ++i)
    {
      /*fscanf returns number of arguments read so should return 1*/
      if ( fscanf(fptr,"%lf",&parm_data->teq[i])!=1 )
      {
        /*Error failed to find bond eq constant - Issue a warning and set to 1.0*/
        printf("!  WARNING - Failed to find angle eq constant for angle %d - assuming it is 0.0.\n",i);
        parm_data->teq[i]=0.0;
      }
    }
    if (global_options->VERBOSITY>=HIGH)
      printf("   Prmtop   (angles): Read in angle equilibrium constants TEQ.\n");
  }

  /*NOW THE DIHEDRAL PARAMS*/
  if ( parm_data->MPTRA > 0 )
  {
    if ( parm_data->newparm )
    {
      /*Test if we can find the DIHEDRAL_FORCE_CONSTANT flag*/
      retval=find_flag( fptr, "DIHEDRAL_FORCE_CONSTANT" );

      if ( retval!=SUCCESS )
      {
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** FAILED TO FIND DIHEDRAL_FORCE_CONSTANT IN PRMTOP FILE: %s\n",parm_data->filename);
        return INVALID_FORMAT;
      }
      // Skip comment lines
      prmtop_data[0]='%';
      while (prmtop_data[0]=='%') {
        fgetpos(fptr, &pos); // Get position of beginning of line in case this is a good line
        if (fgets (prmtop_data,BUFFER_SIZE,fptr) == NULL) {
          printf("*** ERROR IN readparm SUBROUTINE\n");
          printf("*** HIT EOF WHILE TRYING TO READ FORMAT CARD IN PRMTOP FILE: %s\n",parm_data->filename);
          return FILE_READ_FAIL;
        }
      }
      fsetpos(fptr, &pos); // Go to beginning of good line
    }

    /*Allocate memory for the force constants*/
    if (global_options->VERBOSITY>=HIGH)
      printf("   Allocating %d bytes for parm_data->*pk\n",(int)((parm_data->MPTRA)*sizeof(double)));
    parm_data->pk = (double *) malloc((parm_data->MPTRA)*sizeof(double));
    if (parm_data->pk== NULL)
    {
      malloc_failure_char("read_prmtop", "parm_data->pk", (parm_data->MPTRA)*sizeof(double));
      return ALLOC_FAIL;
    }
    parm_data->mem_allocated+=(parm_data->MPTRA)*sizeof(double);
    for (i=0; i < parm_data->MPTRA; ++i)
    {
      /*fscanf returns number of arguments read so should return 1*/
      if ( fscanf(fptr,"%lf",&parm_data->pk[i])!=1 )
      {
        /*Error failed to find force constant - Issue a warning and set to 0.0*/
        printf("WARNING - Failed to find force constant for dihedral %d - assuming it is 0.0.\n",i);
        parm_data->pk[i]=0.0;
      }
    }
    if (global_options->VERBOSITY>=HIGH)
      printf("   Prmtop(dihedrals): Read in dihedral force constants PK.\n");

    if ( parm_data->newparm )
    {
      /*Test if we can find the DIHEDRAL_PERIODICITY flag*/
      retval=find_flag( fptr, "DIHEDRAL_PERIODICITY" );

      if ( retval!=SUCCESS )
      {
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** FAILED TO FIND DIHEDRAL_PERIODICITY IN PRMTOP FILE: %s\n",parm_data->filename);
        return INVALID_FORMAT;
      }
      retval=s_getline( prmtop_data, BUFFER_SIZE, fptr ); /*%FORMAT(5E16.8)*/
      if ( retval==FILE_READ_FAIL )
      {
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** HIT EOF WHILE TRYING TO READ FORMAT CARD IN PRMTOP FILE: %s\n",parm_data->filename);
        return FILE_READ_FAIL;
      }
    }

    /*Allocate memory for the dihedral periodicity constants*/
    if (global_options->VERBOSITY>=HIGH)
      printf("   Allocating %d bytes for parm_data->*pn\n",(int)((parm_data->MPTRA)*sizeof(double)));
    parm_data->pn = (double *) malloc((parm_data->MPTRA)*sizeof(double));
    if (parm_data->pn== NULL)
    {
      malloc_failure_char("read_prmtop", "parm_data->pn", (parm_data->MPTRA)*sizeof(double));
      return ALLOC_FAIL;
    }
    parm_data->mem_allocated+=(parm_data->MPTRA)*sizeof(double);
    for (i=0; i < parm_data->MPTRA; ++i)
    {
      /*fscanf returns number of arguments read so should return 1*/
      if ( fscanf(fptr,"%lf",&parm_data->pn[i])!=1 )
      {
        /*Error failed to find periodicity - Issue a warning and set to 1.0*/
        printf("!  WARNING - Failed to find periodicity for dihedral %d - assuming it is 1.0.\n",i);
        parm_data->pn[i]=1.0;
      }
    }
    if (global_options->VERBOSITY>=HIGH)
      printf("   Prmtop(dihedrals): Read in dihedral periodicity constants PN.\n");

    if ( parm_data->newparm )
    {
      /*Test if we can find the DIHEDRAL_PHASE flag*/
      retval=find_flag( fptr, "DIHEDRAL_PHASE" );

      if ( retval!=SUCCESS )
      {
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** FAILED TO FIND DIHEDRAL_PHASE IN PRMTOP FILE: %s\n",parm_data->filename);
        return INVALID_FORMAT;
      }
      // Skip comment lines
      prmtop_data[0]='%';
      while (prmtop_data[0]=='%') {
        fgetpos(fptr, &pos); // Get position of beginning of line in case this is a good line
        if (fgets (prmtop_data,BUFFER_SIZE,fptr) == NULL) {
          printf("*** ERROR IN readparm SUBROUTINE\n");
          printf("*** HIT EOF WHILE TRYING TO READ FORMAT CARD IN PRMTOP FILE: %s\n",parm_data->filename);
          return FILE_READ_FAIL;
        }
      }
      fsetpos(fptr, &pos); // Go to beginning of good line
    }

    /*Allocate memory for the dihedral phase constants*/
    if (global_options->VERBOSITY>=HIGH)
      printf("   Allocating %d bytes for parm_data->*phase\n",(int)((parm_data->MPTRA)*sizeof(double)));
    parm_data->phase = (double *) malloc((parm_data->MPTRA)*sizeof(double));
    if (parm_data->phase== NULL)
    {
      malloc_failure_char("read_prmtop", "parm_data->phase", (parm_data->MPTRA)*sizeof(double));
      return ALLOC_FAIL;
    }
    parm_data->mem_allocated+=(parm_data->MPTRA)*sizeof(double);
    for (i=0; i < parm_data->MPTRA; ++i)
    {
      /*fscanf returns number of arguments read so should return 1*/
      if ( fscanf(fptr,"%lf",&parm_data->phase[i])!=1 )
      {
        /*Error failed to find periodicity - Issue a warning and set to 1.0*/
        printf("!  WARNING - Failed to find phase for dihedral %d - assuming it is 0.0.\n",i);
        parm_data->phase[i]=0.0;
      }
    }
    if (global_options->VERBOSITY>=HIGH)
      printf("   Prmtop(dihedrals): Read in dihedral phase constants PHASE.\n");
  }

  /*NEXT READ IN SOLTY*/
  if (parm_data->NATYP > 0)
  {
    skip_section=NO;
    if ( parm_data->newparm )
    {
      /*Test if we can find the SOLTY flag*/
      retval=find_flag( fptr, "SOLTY" );

      if ( retval!=SUCCESS )
      {
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** FAILED TO FIND SOLTY CARD IN PRMTOP FILE: %s\n",parm_data->filename);
        printf("*** NOT FATAL, attempting to continue.\n");
        skip_section = YES;
      }
      // Skip comment lines
      prmtop_data[0]='%';
      while (prmtop_data[0]=='%') {
        fgetpos(fptr, &pos); // Get position of beginning of line in case this is a good line
        if (fgets (prmtop_data,BUFFER_SIZE,fptr) == NULL) {
          printf("*** ERROR IN readparm SUBROUTINE\n");
          printf("*** HIT EOF WHILE TRYING TO READ FORMAT CARD IN PRMTOP FILE: %s\n",parm_data->filename);
          return FILE_READ_FAIL;
        }
      }
      fsetpos(fptr, &pos); // Go to beginning of good line
    }

    /*Allocate memory for the solty*/
    if (global_options->VERBOSITY>=HIGH)
      printf("   Allocating %d bytes for parm_data->*solty\n",(int)((parm_data->NATYP)*sizeof(double)));
    parm_data->solty = (double *) malloc((parm_data->NATYP)*sizeof(double));
    if (parm_data->solty== NULL)
    {
      malloc_failure_char("read_prmtop", "parm_data->solty", (parm_data->NATYP)*sizeof(double));
      return ALLOC_FAIL;
    }
    parm_data->mem_allocated+=(parm_data->NATYP)*sizeof(double);
    for (i=0; i < parm_data->NATYP; ++i)
    {
      /*fscanf returns number of arguments read so should return 1*/
      if (skip_section!=YES) {
        if ( fscanf(fptr,"%lf",&parm_data->solty[i])!=1 ) {
          /*Error failed to solty coefficient - Issue a warning and set to 0.0*/
          printf("!  WARNING - Failed to find solty coefficient for atom type %d - assuming it is 0.0.\n",i);
          parm_data->solty[i]=0.0;
        }
      } else {
        /*No solty flag - set all to zero*/
        printf("!  WARNING - Failed to find solty coefficient for atom type %d - assuming it is 0.0.\n",i);
        parm_data->solty[i]=0.0;
      }
    }
    if (global_options->VERBOSITY>=HIGH) {
      printf("   Prmtop    (solty): Read in solty coefficients.\n");
    }
  }

  /*NEXT STEP IS TO READ IN LENNARD JONES COEFFICIENTS FOR CN1 (a) and CN2 (b)*/
  if (parm_data->NTYPES>0)
  {
    if ( parm_data->newparm )
    {
      /*Test if we can find the LENNARD_JONES_ACOEF flag*/
      retval=find_flag( fptr, "LENNARD_JONES_ACOEF" );

      if ( retval!=SUCCESS ) {
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** FAILED TO FIND LENNARD_JONES_ACOEF CARD IN PRMTOP FILE: %s\n",parm_data->filename);
        return INVALID_FORMAT;
      }
      // Skip comment lines
      prmtop_data[0]='%';
      while (prmtop_data[0]=='%') {
        fgetpos(fptr, &pos); // Get position of beginning of line in case this is a good line
        if (fgets (prmtop_data,BUFFER_SIZE,fptr) == NULL) {
          printf("*** ERROR IN readparm SUBROUTINE\n");
          printf("*** HIT EOF WHILE TRYING TO READ FORMAT CARD IN PRMTOP FILE: %s\n",parm_data->filename);
          return FILE_READ_FAIL;
        }
      }
      fsetpos(fptr, &pos); // Go to beginning of good line
    }

    /*Allocate memory for the CN1*/
    /*Requires (NTYPES^2 + NTYPES)/2*/
    len = parm_data->NTYPES;
    len = len * (len + 1) / 2;
    if (global_options->VERBOSITY>=HIGH)
      printf("   Allocating %d bytes for parm_data->*cn1\n",(int)(len*sizeof(double)));
    parm_data->cn1 = (double *) malloc(len*sizeof(double));
    if (parm_data->cn1== NULL)
    {
      malloc_failure_double("read_prmtop", "parm_data->cn1", len*sizeof(double));
      return ALLOC_FAIL;
    }
    parm_data->mem_allocated+=len*sizeof(double);
    for (i=0; i < len; ++i)
    {
      /*fscanf returns number of arguments read so should return 1*/
      if ( fscanf(fptr,"%lf",&parm_data->cn1[i])!=1 )
      {
        /*Error failed to cn1 coefficient - Issue a warning and set to 0.0*/
        printf("!  WARNING - Failed to find cn1 coefficient for element %d - assuming it is 0.0.\n",i);
        parm_data->cn1[i]=0.0;
      }
    }
    if (global_options->VERBOSITY>=HIGH)
      printf("   Prmtop  (L-J CN1): Read in J-N CN1 coefficients.\n");

    if ( parm_data->newparm )
    {
      /*Test if we can find the LENNARD_JONES_BCOEF flag*/
      retval=find_flag( fptr, "LENNARD_JONES_BCOEF" );

      if ( retval!=SUCCESS ) {
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** FAILED TO FIND LENNARD_JONES_BCOEF CARD IN PRMTOP FILE: %s\n",parm_data->filename);
        return INVALID_FORMAT;
      }
      // Skip comment lines
      prmtop_data[0]='%';
      while (prmtop_data[0]=='%') {
        fgetpos(fptr, &pos); // Get position of beginning of line in case this is a good line
        if (fgets (prmtop_data,BUFFER_SIZE,fptr) == NULL) {
          printf("*** ERROR IN readparm SUBROUTINE\n");
          printf("*** HIT EOF WHILE TRYING TO READ FORMAT CARD IN PRMTOP FILE: %s\n",parm_data->filename);
          return FILE_READ_FAIL;
        }
      }
      fsetpos(fptr, &pos); // Go to beginning of good line
    }

    /*Allocate memory for CN2*/
    if (global_options->VERBOSITY>=HIGH)
      printf("   Allocating %d bytes for parm_data->*cn2\n",(int)(len*sizeof(double)));
    parm_data->cn2 = (double *) malloc(len*sizeof(double));
    if (parm_data->cn2== NULL)
    {
      malloc_failure_double("read_prmtop", "parm_data->cn2", len*sizeof(double));
      return ALLOC_FAIL;
    }
    parm_data->mem_allocated+=len*sizeof(double);
    for (i=0; i < len; ++i)
    {
      /*fscanf returns number of arguments read so should return 1*/
      if ( fscanf(fptr,"%lf",&parm_data->cn2[i])!=1 )
      {
        /*Error failed to cn2 coefficient - Issue a warning and set to 0.0*/
        printf("!  WARNING - Failed to find cn2 coefficient for element %d - assuming it is 0.0.\n",i);
        parm_data->cn2[i]=0.0;
      }
    }
    if (global_options->VERBOSITY>=HIGH)
      printf("   Prmtop  (L-J CN2): Read in J-N CN2 coefficients.\n");
  }

  /*NOW WE READ IN INFO FOR BONDS WITH HYDROGEN*/
  if ( parm_data->NBONH > 0)
  {
    /*step 1 allocate memory*/
    if (global_options->VERBOSITY>=HIGH)
      printf("   Allocating %d bytes for parm_data->*pbondH\n",(int)(parm_data->NBONH*sizeof(parmbond_struct)));
    parm_data->pbondH = (parmbond_struct *) malloc(parm_data->NBONH*sizeof(parmbond_struct));
    if (parm_data->pbondH== NULL)
    {
      malloc_failure_char("read_prmtop", "parm_data->pbondH", parm_data->NBONH*sizeof(parmbond_struct));
      return ALLOC_FAIL;
    }
    parm_data->mem_allocated+=parm_data->NBONH*sizeof(parmbond_struct);
    if ( parm_data->newparm )
    {
      /*Test if we can find the BONDS_INC_HYDROGEN flag*/
      retval=find_flag( fptr, "BONDS_INC_HYDROGEN" );

      if ( retval!=SUCCESS )
      {
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** FAILED TO FIND BONDS_INC_HYDROGEN CARD IN PRMTOP FILE: %s\n",parm_data->filename);
        return INVALID_FORMAT;
      }
      // Skip comment lines
      prmtop_data[0]='%';
      while (prmtop_data[0]=='%') {
        fgetpos(fptr, &pos); // Get position of beginning of line in case this is a good line
        if (fgets (prmtop_data,BUFFER_SIZE,fptr) == NULL) {
          printf("*** ERROR IN readparm SUBROUTINE\n");
          printf("*** HIT EOF WHILE TRYING TO READ FORMAT CARD IN PRMTOP FILE: %s\n",parm_data->filename);
          return FILE_READ_FAIL;
        }
      }
      fsetpos(fptr, &pos); // Go to beginning of good line
    }

    for (i=0; i < parm_data->NBONH; ++i) {
      /*3 items per bond, 1st = 1st atom, 2nd = 2nd atom, 3rd = point to params*/
      /*fscanf returns number of arguments read so should return 1*/
      if ( fscanf(fptr,"%d",&parm_data->pbondH[i].ib)!=1 )
      {
        /*Error failed to find first atom in bond - Issue an error and quit*/
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** FAILED TO FIND FIRST ATOM IN BOND %d\n",i);
        return INVALID_FORMAT;
      }
      if ( fscanf(fptr,"%d",&parm_data->pbondH[i].jb)!=1 )
      {
        /*Error failed to find second atom in bond - Issue an error and quit*/
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** FAILED TO FIND SECOND ATOM IN BOND %d\n",i);
        return INVALID_FORMAT;
      }
      if ( fscanf(fptr,"%d",&parm_data->pbondH[i].icb)!=1 )
      {
        /*Error failed to find bond param pointer for bond - Issue an error and quit*/
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** FAILED TO FIND BOND PARAMETER POINTER FOR BOND %d\n",i);
        return INVALID_FORMAT;
      }
    }
    if (global_options->VERBOSITY>=HIGH) {
      printf("   Prmtop    (bonds): Read in info for bonds with hydrogen.\n");
    }
  }

  /*NOW WE READ IN INFO FOR BONDS WITHOUT HYDROGEN*/
  if (parm_data->NBONA > 0)
  {
    /*step 1 allocate memory*/
    if (global_options->VERBOSITY>=HIGH)
      printf("   Allocating %d bytes for parm_data->*pbond\n",(int)(parm_data->NBONA*sizeof(parmbond_struct)));
    parm_data->pbond = (parmbond_struct *) malloc(parm_data->NBONA*sizeof(parmbond_struct));
    if (parm_data->pbond== NULL)
    {
      malloc_failure_char("read_prmtop", "parm_data->pbond", parm_data->NBONA*sizeof(parmbond_struct));
      return ALLOC_FAIL;
    }
    parm_data->mem_allocated+=parm_data->NBONA*sizeof(parmbond_struct);
    if ( parm_data->newparm )
    {
      /*Test if we can find the BONDS_WITHOUT_HYDROGEN flag*/
      retval=find_flag( fptr, "BONDS_WITHOUT_HYDROGEN" );

      if ( retval!=SUCCESS )
      {
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** FAILED TO FIND BONDS_WITHOUT_HYDROGEN CARD IN PRMTOP FILE: %s\n",parm_data->filename);
        return INVALID_FORMAT;
      }
      // Skip comment lines
      prmtop_data[0]='%';
      while (prmtop_data[0]=='%') {
        fgetpos(fptr, &pos); // Get position of beginning of line in case this is a good line
        if (fgets (prmtop_data,BUFFER_SIZE,fptr) == NULL) {
          printf("*** ERROR IN readparm SUBROUTINE\n");
          printf("*** HIT EOF WHILE TRYING TO READ FORMAT CARD IN PRMTOP FILE: %s\n",parm_data->filename);
          return FILE_READ_FAIL;
        }
      }
      fsetpos(fptr, &pos); // Go to beginning of good line
    }

    for (i=0; i < parm_data->NBONA; ++i) {
      /*3 items per bond, 1st = 1st atom, 2nd = 2nd atom, 3rd = point to params*/
      /*fscanf returns number of arguments read so should return 1*/
      if ( fscanf(fptr,"%d",&parm_data->pbond[i].ib)!=1 )
      {
        /*Error failed to find first atom in bond - Issue an error and quit*/
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** FAILED TO FIND FIRST ATOM IN BOND %d\n",i);
        return INVALID_FORMAT;
      }
      if ( fscanf(fptr,"%d",&parm_data->pbond[i].jb)!=1 )
      {
        /*Error failed to find second atom in bond - Issue an error and quit*/
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** FAILED TO FIND SECOND ATOM IN BOND %d\n",i);
        return INVALID_FORMAT;
      }
      if ( fscanf(fptr,"%d",&parm_data->pbond[i].icb)!=1 )
      {
        /*Error failed to find bond param pointer for bond - Issue an error and quit*/
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** FAILED TO FIND BOND PARAMETER POINTER FOR BOND %d\n",i);
        return INVALID_FORMAT;
      }
    }
    if (global_options->VERBOSITY>=HIGH)
      printf("   Prmtop    (bonds): Read in info for bonds without hydrogen.\n");
  }

  /*NOW WE READ IN INFO FOR ANGLES WITH HYDROGEN*/
  if (parm_data->NTHETH > 0 )
  {
    /*step 1 allocate memory*/
    if (global_options->VERBOSITY>=HIGH)
      printf("   Allocating %d bytes for parm_data->*pangleH\n",(int)(parm_data->NTHETH*sizeof(parmangle_struct)));
    parm_data->pangleH = (parmangle_struct *) malloc(parm_data->NTHETH*sizeof(parmangle_struct));
    if (parm_data->pangleH== NULL)
    {
      malloc_failure_char("read_prmtop", "parm_data->pangleH", parm_data->NTHETH*sizeof(parmangle_struct));
      return ALLOC_FAIL;
    }
    parm_data->mem_allocated+=parm_data->NTHETH*sizeof(parmangle_struct);
    if ( parm_data->newparm )
    {
      /*Test if we can find the ANGLES_INC_HYDROGEN flag*/
      retval=find_flag( fptr, "ANGLES_INC_HYDROGEN" );

      if ( retval!=SUCCESS )
      {
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** FAILED TO FIND ANGLES_INC_HYDROGEN CARD IN PRMTOP FILE: %s\n",parm_data->filename);
        return INVALID_FORMAT;
      }
      // Skip comment lines
      prmtop_data[0]='%';
      while (prmtop_data[0]=='%') {
        fgetpos(fptr, &pos); // Get position of beginning of line in case this is a good line
        if (fgets (prmtop_data,BUFFER_SIZE,fptr) == NULL) {
          printf("*** ERROR IN readparm SUBROUTINE\n");
          printf("*** HIT EOF WHILE TRYING TO READ FORMAT CARD IN PRMTOP FILE: %s\n",parm_data->filename);
          return FILE_READ_FAIL;
        }
      }
      fsetpos(fptr, &pos); // Go to beginning of good line
    }

    for (i=0; i < parm_data->NTHETH; ++i) {
      /*4 items per angle, 1st = 1st atom, 2nd = 2nd atom, 3rd = 3rd atom, 4th = pointer to param*/
      /*fscanf returns number of arguments read so should return 1*/
      if ( fscanf(fptr,"%d",&parm_data->pangleH[i].it)!=1 )
      {
        /*Error failed to find first atom in angle - Issue an error and quit*/
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** FAILED TO FIND FIRST ATOM IN ANGLE %d\n",i);
        return INVALID_FORMAT;
      }
      if ( fscanf(fptr,"%d",&parm_data->pangleH[i].jt)!=1 )
      {
        /*Error failed to find second atom in angle - Issue an error and quit*/
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** FAILED TO FIND SECOND ATOM IN ANGLE %d\n",i);
        return INVALID_FORMAT;
      }
      if ( fscanf(fptr,"%d",&parm_data->pangleH[i].kt)!=1 )
      {
        /*Error failed to find third atom in angle - Issue an error and quit*/
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** FAILED TO FIND THIRD ATOM IN ANGLE %d\n",i);
        return INVALID_FORMAT;
      }
      if ( fscanf(fptr,"%d",&parm_data->pangleH[i].ict)!=1 )
      {
        /*Error failed to find bond param pointer for angle - Issue an error and quit*/
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** FAILED TO FIND ANGLE PARAMETER POINTER FOR ANGLE %d\n",i);
        return INVALID_FORMAT;
      }
    }
    if (global_options->VERBOSITY>=HIGH)
      printf("   Prmtop   (angles): Read in info for angles with hydrogen.\n");
  }
  /*NOW WE READ IN INFO FOR ANGLES WITHOUT HYDROGEN*/
  if (parm_data->MTHETS > 0)
  {
    /*step 1 allocate memory*/
    if (global_options->VERBOSITY>=HIGH)
      printf("   Allocating %d bytes for parm_data->*pangle\n",(int)(parm_data->MTHETS*sizeof(parmangle_struct)));
    parm_data->pangle = (parmangle_struct *) malloc(parm_data->MTHETS*sizeof(parmangle_struct));
    if (parm_data->pangle== NULL)
    {
      malloc_failure_char("read_prmtop", "parm_data->pangle", parm_data->MTHETS*sizeof(parmangle_struct));
      return ALLOC_FAIL;
    }
    parm_data->mem_allocated+=parm_data->MTHETS*sizeof(parmangle_struct);
    if ( parm_data->newparm )
    {
      /*Test if we can find the ANGLES_WITHOUT_HYDROGEN flag*/
      retval=find_flag( fptr, "ANGLES_WITHOUT_HYDROGEN" );

      if ( retval!=SUCCESS )
      {
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** FAILED TO FIND ANGLES_WITHOUT_HYDROGEN CARD IN PRMTOP FILE: %s\n",parm_data->filename);
        return INVALID_FORMAT;
      }
      // Skip comment lines
      prmtop_data[0]='%';
      while (prmtop_data[0]=='%') {
        fgetpos(fptr, &pos); // Get position of beginning of line in case this is a good line
        if (fgets (prmtop_data,BUFFER_SIZE,fptr) == NULL) {
          printf("*** ERROR IN readparm SUBROUTINE\n");
          printf("*** HIT EOF WHILE TRYING TO READ FORMAT CARD IN PRMTOP FILE: %s\n",parm_data->filename);
          return FILE_READ_FAIL;
        }
      }
      fsetpos(fptr, &pos); // Go to beginning of good line
    }

    for (i=0; i < parm_data->MTHETS; ++i) {
      /*4 items per angle, 1st = 1st atom, 2nd = 2nd atom, 3rd = 3rd atom, 4th = pointer to param*/
      /*fscanf returns number of arguments read so should return 1*/
      if ( fscanf(fptr,"%d",&parm_data->pangle[i].it)!=1 )
      {
        /*Error failed to find first atom in angle - Issue an error and quit*/
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** FAILED TO FIND FIRST ATOM IN ANGLE %d\n",i);
        return INVALID_FORMAT;
      }
      if ( fscanf(fptr,"%d",&parm_data->pangle[i].jt)!=1 )
      {
        /*Error failed to find second atom in angle - Issue an error and quit*/
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** FAILED TO FIND SECOND ATOM IN ANGLE %d\n",i);
        return INVALID_FORMAT;
      }
      if ( fscanf(fptr,"%d",&parm_data->pangle[i].kt)!=1 )
      {
        /*Error failed to find third atom in angle - Issue an error and quit*/
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** FAILED TO FIND THIRD ATOM IN ANGLE %d\n",i);
        return INVALID_FORMAT;
      }
      if ( fscanf(fptr,"%d",&parm_data->pangle[i].ict)!=1 )
      {
        /*Error failed to find bond param pointer for angle - Issue an error and quit*/
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** FAILED TO FIND ANGLE PARAMETER POINTER FOR ANGLE %d\n",i);
        return INVALID_FORMAT;
      }
    }
    if (global_options->VERBOSITY>=HIGH)
      printf("   Prmtop   (angles): Read in info for angles without hydrogen.\n");
  }

  /*NOW WE READ IN INFO FOR DIHEDRALS INC HYDROGEN*/

  parm_data->pdihedralH = NULL; //default
  if (parm_data->NPHIH > 0)
  {
    /*step 1 allocate memory*/
    if (global_options->VERBOSITY>=HIGH)
      printf("   Allocating %d bytes for parm_data->*pdihedralH\n",(int)(parm_data->NPHIH*sizeof(parmdihedral_struct)));
    parm_data->pdihedralH = (parmdihedral_struct *) malloc(parm_data->NPHIH*sizeof(parmdihedral_struct));

    if (parm_data->pdihedralH== NULL)
    {
      malloc_failure_char("read_prmtop", "parm_data->pdihedralH", parm_data->NPHIH*sizeof(parmdihedral_struct));
      return ALLOC_FAIL;
    }
    parm_data->mem_allocated+=parm_data->NPHIH*sizeof(parmdihedral_struct);
    if ( parm_data->newparm )
    {
      /*Test if we can find the DIHEDRALS_INC_HYDROGEN flag*/
      retval=find_flag( fptr, "DIHEDRALS_INC_HYDROGEN" );

      if ( retval!=SUCCESS )
      {
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** FAILED TO FIND DIHEDRALS_INC_HYDROGEN CARD IN PRMTOP FILE: %s\n",parm_data->filename);
        return INVALID_FORMAT;
      }
      // Skip comment lines
      prmtop_data[0]='%';
      while (prmtop_data[0]=='%') {
        fgetpos(fptr, &pos); // Get position of beginning of line in case this is a good line
        if (fgets (prmtop_data,BUFFER_SIZE,fptr) == NULL) {
          printf("*** ERROR IN readparm SUBROUTINE\n");
          printf("*** HIT EOF WHILE TRYING TO READ FORMAT CARD IN PRMTOP FILE: %s\n",parm_data->filename);
          return FILE_READ_FAIL;
        }
      }
      fsetpos(fptr, &pos); // Go to beginning of good line
    }
    for (i=0; i < parm_data->NPHIH; ++i)
    {
      /*5 items per dihedral, 1st = 1st atom, 2nd = 2nd atom, 3rd = 3rd atom, 4th = 4th atom 5th = pointer to param*/
      /*fscanf returns number of arguments read so should return 1*/
      if ( fscanf(fptr,"%d",&parm_data->pdihedralH[i].ip)!=1 )
      {
        /*Error failed to find first atom in angle - Issue an error and quit*/
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** FAILED TO FIND FIRST ATOM IN DIHEDRAL %d\n",i);
        return INVALID_FORMAT;
      }
      if ( fscanf(fptr,"%d",&parm_data->pdihedralH[i].jp)!=1 )
      {
        /*Error failed to find second atom in angle - Issue an error and quit*/
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** FAILED TO FIND SECOND ATOM IN DIHEDRAL %d\n",i);
        return INVALID_FORMAT;
      }
      if ( fscanf(fptr,"%d",&parm_data->pdihedralH[i].kp)!=1 )
      {
        /*Error failed to find third atom in angle - Issue an error and quit*/
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** FAILED TO FIND THIRD ATOM IN DIHEDRAL %d\n",i);
        return INVALID_FORMAT;
      }
      if ( fscanf(fptr,"%d",&parm_data->pdihedralH[i].lp)!=1 )
      {
        /*Error failed to find fourth atom in angle - Issue an error and quit*/
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** FAILED TO FIND FOURTH ATOM IN DIHEDRAL %d\n",i);
        return INVALID_FORMAT;
      }
      if ( fscanf(fptr,"%d",&parm_data->pdihedralH[i].icp)!=1 )
      {
        /*Error failed to find  param pointer for dihedral - Issue an error and quit*/
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** FAILED TO FIND PARAMETER POINTER FOR DIHEDRAL %d\n",i);
        return INVALID_FORMAT;
      }
    }
    if (global_options->VERBOSITY>=HIGH) {
      printf("   Prmtop(dihedrals): Read in info for dihedrals with hydrogen.\n");
    }
  }

  /*NOW WE READ IN INFO FOR DIHEDRALS WITHOUT HYDROGEN*/
  parm_data->pdihedral = NULL; // default
  if (parm_data->MPHIA > 0)
  {
    /*step 1 allocate memory*/
    if (global_options->VERBOSITY>=HIGH)
      printf("   Allocating %d bytes for parm_data->*pdihedral\n",(int)(parm_data->MPHIA*sizeof(parmdihedral_struct)));
    parm_data->pdihedral = (parmdihedral_struct *) malloc(parm_data->MPHIA*sizeof(parmdihedral_struct));
    if (parm_data->pdihedral== NULL)
    {
      malloc_failure_char("read_prmtop", "parm_data->pdihedral", parm_data->MPHIA*sizeof(parmdihedral_struct));
      return ALLOC_FAIL;
    }
    parm_data->mem_allocated+=parm_data->MPHIA*sizeof(parmdihedral_struct);
    if ( parm_data->newparm )
    {
      /*Test if we can find the DIHEDRALS_WITHOUT_HYDROGEN flag*/
      retval=find_flag( fptr, "DIHEDRALS_WITHOUT_HYDROGEN" );

      if ( retval!=SUCCESS )
      {
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** FAILED TO FIND DIHEDRALS_WITHOUT_HYDROGEN CARD IN PRMTOP FILE: %s\n",parm_data->filename);
        return INVALID_FORMAT;
      }
      // Skip comment lines
      prmtop_data[0]='%';
      while (prmtop_data[0]=='%') {
        fgetpos(fptr, &pos); // Get position of beginning of line in case this is a good line
        if (fgets (prmtop_data,BUFFER_SIZE,fptr) == NULL) {
          printf("*** ERROR IN readparm SUBROUTINE\n");
          printf("*** HIT EOF WHILE TRYING TO READ FORMAT CARD IN PRMTOP FILE: %s\n",parm_data->filename);
          return FILE_READ_FAIL;
        }
      }
      fsetpos(fptr, &pos); // Go to beginning of good line
    }
    for (i=0; i < parm_data->MPHIA; ++i)
    {
      /*5 items per dihedral, 1st = 1st atom, 2nd = 2nd atom, 3rd = 3rd atom, 4th = 4th atom 5th = pointer to param*/
      /*fscanf returns number of arguments read so should return 1*/
      if ( fscanf(fptr,"%d",&parm_data->pdihedral[i].ip)!=1 )
      {
        /*Error failed to find first atom in angle - Issue an error and quit*/
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** FAILED TO FIND FIRST ATOM IN DIHEDRAL %d\n",i);
        return INVALID_FORMAT;
      }
      if ( fscanf(fptr,"%d",&parm_data->pdihedral[i].jp)!=1 )
      {
        /*Error failed to find second atom in angle - Issue an error and quit*/
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** FAILED TO FIND SECOND ATOM IN DIHEDRAL %d\n",i);
        return INVALID_FORMAT;
      }
      if ( fscanf(fptr,"%d",&parm_data->pdihedral[i].kp)!=1 )
      {
        /*Error failed to find third atom in angle - Issue an error and quit*/
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** FAILED TO FIND THIRD ATOM IN DIHEDRAL %d\n",i);
        return INVALID_FORMAT;
      }
      if ( fscanf(fptr,"%d",&parm_data->pdihedral[i].lp)!=1 )
      {
        /*Error failed to find fourth atom in angle - Issue an error and quit*/
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** FAILED TO FIND FOURTH ATOM IN DIHEDRAL %d\n",i);
        return INVALID_FORMAT;
      }
      if ( fscanf(fptr,"%d",&parm_data->pdihedral[i].icp)!=1 )
      {
        /*Error failed to find  param pointer for dihedral - Issue an error and quit*/
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** FAILED TO FIND PARAMETER POINTER FOR DIHEDRAL %d\n",i);
        return INVALID_FORMAT;
      }
    }
    if (global_options->VERBOSITY>=HIGH)
      printf("   Prmtop(dihedrals): Read in info for dihedrals without hydrogen.\n");
  }

  /*NOW WE READ IN THE EXCLUDED ATOM LIST*/
  if ( parm_data->NEXT > 0 )
  {
    if ( parm_data->newparm )
    {
      /*Test if we can find the EXCLUDED_ATOMS_LIST flag*/
      retval=find_flag( fptr, "EXCLUDED_ATOMS_LIST" );

      if ( retval!=SUCCESS )
      {
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** FAILED TO FIND EXCLUDED_ATOMS_LIST CARD IN PRMTOP FILE: %s\n",parm_data->filename);
        return INVALID_FORMAT;
      }
      // Skip comment lines
      prmtop_data[0]='%';
      while (prmtop_data[0]=='%') {
        fgetpos(fptr, &pos); // Get position of beginning of line in case this is a good line
        if (fgets (prmtop_data,BUFFER_SIZE,fptr) == NULL) {
          printf("*** ERROR IN readparm SUBROUTINE\n");
          printf("*** HIT EOF WHILE TRYING TO READ FORMAT CARD IN PRMTOP FILE: %s\n",parm_data->filename);
          return FILE_READ_FAIL;
        }
      }
      fsetpos(fptr, &pos); // Go to beginning of good line
    }

    /*Allocate memory for the natex*/
    if (global_options->VERBOSITY>=HIGH)
      printf("   Allocating %d bytes for parm_data->*natex\n",(int)((parm_data->NEXT)*sizeof(int)));
    parm_data->natex = (int *) malloc((parm_data->NEXT)*sizeof(int));
    if (parm_data->natex== NULL)
    {
      malloc_failure_char("read_prmtop", "parm_data->natex", (parm_data->NEXT)*sizeof(int));
      return ALLOC_FAIL;
    }
    parm_data->mem_allocated+=parm_data->NEXT*sizeof(int);
    for (i=0; i < parm_data->NEXT; ++i)
    {
      /*fscanf returns number of arguments read so should return 1*/
      if ( fscanf(fptr,"%d",&parm_data->natex[i])!=1 )
      {
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** FAILED TO FIND EXCLUDED_ATOMS_LIST ENTRY FOR ELEMENT: %d\n",i);
        return INVALID_FORMAT;
      }
    }
    if (global_options->VERBOSITY>=HIGH)
      printf("   Prmtop (excluded): Read in excluded atom list.\n");
  }

  /*NEXT WE READ IN THE HBOND PARAMETERS*/
  /*Note these are not currently used in the fitting so we don't actually quit if we
    don't find the flags*/
  if ( parm_data->NHB > 0)
  {
    /*Start with HBOND_ACOEF*/
    skip_section=NO;
    if ( parm_data->newparm )
    {
      /*Test if we can find the HBOND_ACOEF flag*/
      retval=find_flag( fptr, "HBOND_ACOEF" );

      if ( retval!=SUCCESS )
      {
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** FAILED TO FIND HBOND_ACOEF CARD IN PRMTOP FILE: %s\n",parm_data->filename);
        printf("*** NOT FATAL, attempting to continue.\n");
        skip_section = YES;
      }
      // Skip comment lines
      prmtop_data[0]='%';
      while (prmtop_data[0]=='%') {
        fgetpos(fptr, &pos); // Get position of beginning of line in case this is a good line
        if (fgets (prmtop_data,BUFFER_SIZE,fptr) == NULL) {
          printf("*** ERROR IN readparm SUBROUTINE\n");
          printf("*** HIT EOF WHILE TRYING TO READ FORMAT CARD IN PRMTOP FILE: %s\n",parm_data->filename);
          return FILE_READ_FAIL;
        }
      }
      fsetpos(fptr, &pos); // Go to beginning of good line
    }

    /*Allocate memory for the ag*/
    if (global_options->VERBOSITY>=HIGH)
      printf("   Allocating %d bytes for parm_data->*ag\n",(int)((parm_data->NHB)*sizeof(double)));
    parm_data->ag = (double *) malloc((parm_data->NHB)*sizeof(double));
    if (parm_data->ag == NULL)
    {
      malloc_failure_char("read_prmtop", "parm_data->ag", (parm_data->NHB)*sizeof(double));
      return ALLOC_FAIL;
    }
    parm_data->mem_allocated+=parm_data->NHB*sizeof(double);
    for (i=0; i < parm_data->NHB; ++i)
    {
      /*fscanf returns number of arguments read so should return 1*/
      if (skip_section!=YES)
      {
        if ( fscanf(fptr,"%lf",&parm_data->ag[i])!=1 )
        {
          /*Error failed to HBOND_ACOEF coefficient - Issue a warning and set to 0.0*/
          printf("!  WARNING - Failed to find HBOND_ACOEF coefficient for bond %d - assuming it is 0.0.\n",i);
          parm_data->ag[i]=0.0;
        }
      }
      else
      {
        /*No HBOND_ACOEF flag - set all to zero*/
        printf("!  WARNING - Failed to find HBOND_ACOEF coefficient for bond %d - assuming it is 0.0.\n",i);
        parm_data->ag[i]=0.0;
      }
    }
    if (global_options->VERBOSITY>=HIGH)
      printf("   Prmtop   (hbonds): Read in h-bond A coefficients (AG).\n");

    /*Start with HBOND_BCOEF*/
    skip_section=NO;
    if ( parm_data->newparm )
    {
      /*Test if we can find the HBOND_BCOEF flag*/
      retval=find_flag( fptr, "HBOND_BCOEF" );

      if ( retval!=SUCCESS )
      {
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** FAILED TO FIND HBOND_BCOEF CARD IN PRMTOP FILE: %s\n",parm_data->filename);
        printf("*** NOT FATAL, attempting to continue.\n");
        skip_section = YES;
      }
      // Skip comment lines
      prmtop_data[0]='%';
      while (prmtop_data[0]=='%') {
        fgetpos(fptr, &pos); // Get position of beginning of line in case this is a good line
        if (fgets (prmtop_data,BUFFER_SIZE,fptr) == NULL) {
          printf("*** ERROR IN readparm SUBROUTINE\n");
          printf("*** HIT EOF WHILE TRYING TO READ FORMAT CARD IN PRMTOP FILE: %s\n",parm_data->filename);
          return FILE_READ_FAIL;
        }
      }
      fsetpos(fptr, &pos); // Go to beginning of good line
    }

    /*Allocate memory for the bg*/
    if (global_options->VERBOSITY>=HIGH)
      printf("   Allocating %d bytes for parm_data->*bg\n",(int)((parm_data->NHB)*sizeof(double)));
    parm_data->bg = (double *) malloc((parm_data->NHB)*sizeof(double));
    if (parm_data->bg == NULL)
    {
      malloc_failure_char("read_prmtop", "parm_data->bg", (parm_data->NHB)*sizeof(double));
      return ALLOC_FAIL;
    }
    parm_data->mem_allocated+=parm_data->NHB*sizeof(double);
    for (i=0; i < parm_data->NHB; ++i)
    {
      /*fscanf returns number of arguments read so should return 1*/
      if (skip_section!=YES)
      {
        if ( fscanf(fptr,"%lf",&parm_data->bg[i])!=1 )
        {
          /*Error failed to HBOND_ACOEF coefficient - Issue a warning and set to 0.0*/
          printf("!  WARNING - Failed to find HBOND_BCOEF coefficient for bond %d - assuming it is 0.0.\n",i);
          parm_data->bg[i]=0.0;
        }
      }
      else
      {
        /*No HBOND_BCOEF flag - set all to zero*/
        printf("!  WARNING - Failed to find HBOND_BCOEF coefficient for bond %d - assuming it is 0.0.\n",i);
        parm_data->bg[i]=0.0;
      }
    }
    if (global_options->VERBOSITY>=HIGH)
      printf("   Prmtop   (hbonds): Read in h-bond B coefficients (BG).\n");

    /*Now do HBCUT*/
    skip_section=NO;
    if ( parm_data->newparm )
    {
      /*Test if we can find the HBCUT flag*/
      retval=find_flag( fptr, "HBCUT" );

      if ( retval!=SUCCESS )
      {
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** FAILED TO FIND HBCUT CARD IN PRMTOP FILE: %s\n",parm_data->filename);
        printf("*** NOT FATAL, attempting to continue.\n");
        skip_section = YES;
      }
      // Skip comment lines
      prmtop_data[0]='%';
      while (prmtop_data[0]=='%') {
        fgetpos(fptr, &pos); // Get position of beginning of line in case this is a good line
        if (fgets (prmtop_data,BUFFER_SIZE,fptr) == NULL) {
          printf("*** ERROR IN readparm SUBROUTINE\n");
          printf("*** HIT EOF WHILE TRYING TO READ FORMAT CARD IN PRMTOP FILE: %s\n",parm_data->filename);
          return FILE_READ_FAIL;
        }
      }
      fsetpos(fptr, &pos); // Go to beginning of good line
    }

    /*Allocate memory for the hbcut*/
    if (global_options->VERBOSITY>=HIGH)
      printf("   Allocating %d bytes for parm_data->*hbcut\n",(int)((parm_data->NHB)*sizeof(double)));
    parm_data->hbcut = (double *) malloc((parm_data->NHB)*sizeof(double));
    if (parm_data->hbcut == NULL)
    {
      malloc_failure_char("read_prmtop", "parm_data->hbcut", (parm_data->NHB)*sizeof(double));
      return ALLOC_FAIL;
    }
    parm_data->mem_allocated+=parm_data->NHB*sizeof(double);
    for (i=0; i < parm_data->NHB; ++i)
    {
      /*fscanf returns number of arguments read so should return 1*/
      if (skip_section!=YES) {
        if ( fscanf(fptr,"%lf",&parm_data->hbcut[i])!=1 ) {
          /*Error failed to HBCUT coefficient - Issue a warning and set to 0.0*/
          printf("!  WARNING - Failed to find HBCUT coefficient for bond %d - assuming it is 0.0.\n",i);
          parm_data->hbcut[i]=0.0;
        }
      } else {
        /*No HBCUT flag - set all to zero*/
        printf("!  WARNING - Failed to find HBCUT coefficient for bond %d - assuming it is 0.0.\n",i);
        parm_data->hbcut[i]=0.0;
      }
    }
    if (global_options->VERBOSITY>=HIGH)
      printf("   Prmtop   (hbonds): Read in h-bond cut coefficients (HBCUT).\n");
  }

  /*NOW WE READ THE AMBER ATOM TYPES*/
  if ( parm_data->NTOTAT > 0 )
  {
    if ( parm_data->newparm )
    {
      /*Test if we can find the AMBER_ATOM_TYPE flag*/
      retval=find_flag( fptr, "AMBER_ATOM_TYPE" );

      if ( retval!=SUCCESS )
      {
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** FAILED TO FIND AMBER_ATOM_TYPE CARD IN PRMTOP FILE: %s\n",parm_data->filename);
        return INVALID_FORMAT;
      }
      // Skip comment lines
      prmtop_data[0]='%';
      while (prmtop_data[0]=='%') {
        fgetpos(fptr, &pos); // Get position of beginning of line in case this is a good line
        if (fgets (prmtop_data,BUFFER_SIZE,fptr) == NULL) {
          printf("*** ERROR IN readparm SUBROUTINE\n");
          printf("*** HIT EOF WHILE TRYING TO READ FORMAT CARD IN PRMTOP FILE: %s\n",parm_data->filename);
          return FILE_READ_FAIL;
        }
      }
      fsetpos(fptr, &pos); // Go to beginning of good line
    }

    for (i=0; i < parm_data->NTOTAT; ++i) {
      retval=name_copy(fptr, parm_data->atom[i].isymbl);
      if ( retval != SUCCESS)
      {
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** READ INVALID ATOM TYPE SYMBOL FROM PRMTOP FILE: %s\n",parm_data->filename);
        return INVALID_FORMAT;
      }
    }
    if (global_options->VERBOSITY>=HIGH)
      printf("   Prmtop    (atoms): Read in atom symbols (types).\n");
  }

  /*NEXT, WE READ THE TREE INFO*/
  if ( parm_data->NTOTAT > 0 )
  {
    if ( parm_data->newparm )
    {
      /*Test if we can find the TREE_CHAIN_CLASSIFICATION flag*/
      retval=find_flag( fptr, "TREE_CHAIN_CLASSIFICATION" );

      if ( retval!=SUCCESS )
      {
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** FAILED TO FIND TREE_CHAIN_CLASSIFICATION CARD IN PRMTOP FILE: %s\n",parm_data->filename);
        return INVALID_FORMAT;
      }
      // Skip comment lines
      prmtop_data[0]='%';
      while (prmtop_data[0]=='%') {
        fgetpos(fptr, &pos); // Get position of beginning of line in case this is a good line
        if (fgets (prmtop_data,BUFFER_SIZE,fptr) == NULL) {
          printf("*** ERROR IN readparm SUBROUTINE\n");
          printf("*** HIT EOF WHILE TRYING TO READ FORMAT CARD IN PRMTOP FILE: %s\n",parm_data->filename);
          return FILE_READ_FAIL;
        }
      }
      fsetpos(fptr, &pos); // Go to beginning of good line
    }

    for (i=0; i < parm_data->NTOTAT; ++i) {
      retval=name_copy(fptr, parm_data->atom[i].itree);
      if ( retval != SUCCESS)
      {
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** READ INVALID ATOM TREE CHAIN ID FROM PRMTOP FILE: %s\n",parm_data->filename);
        return INVALID_FORMAT;
      }
    }
    if (global_options->VERBOSITY>=HIGH) {
      printf("   Prmtop    (atoms): Read in tree information.\n");
    }
  }

  /*NOW READ JOIN ARRAY*/
  if (parm_data->NTOTAT > 0)
  {
    if ( parm_data->newparm )
    {
      /*Test if we can find the JOIN_ARRAY flag*/
      retval=find_flag( fptr, "JOIN_ARRAY" );

      if ( retval!=SUCCESS )
      {
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** FAILED TO FIND JOIN_ARRAY CARD IN PRMTOP FILE: %s\n",parm_data->filename);
        return INVALID_FORMAT;
      }
      // Skip comment lines
      prmtop_data[0]='%';
      while (prmtop_data[0]=='%') {
        fgetpos(fptr, &pos); // Get position of beginning of line in case this is a good line
        if (fgets (prmtop_data,BUFFER_SIZE,fptr) == NULL) {
          printf("*** ERROR IN readparm SUBROUTINE\n");
          printf("*** HIT EOF WHILE TRYING TO READ FORMAT CARD IN PRMTOP FILE: %s\n",parm_data->filename);
          return FILE_READ_FAIL;
        }
      }
      fsetpos(fptr, &pos); // Go to beginning of good line
    }

    for (i=0; i < parm_data->NTOTAT; ++i) {
      /*fscanf returns number of arguments read so should return 1*/
      if ( fscanf(fptr,"%d",&parm_data->atom[i].join)!=1 )
      {
        /*Fatal Error failed to find atom type*/
        printf("*** ERROR - FAILED TO FIND JOIN INFO FOR ATOM %d\n",i);
        return INVALID_FORMAT;
      }
    }

    if (global_options->VERBOSITY>=HIGH) {
      printf("   Prmtop    (atoms): Read in atom join info.\n");
    }
  }

  /*NOW DO IROTAT*/
  if (parm_data->NTOTAT > 0)
  {
    if ( parm_data->newparm )
    {
      /*Test if we can find the IROTAT flag*/
      retval=find_flag( fptr, "IROTAT" );

      if ( retval!=SUCCESS ) {
        printf("*** ERROR IN readparm SUBROUTINE\n");
        printf("*** FAILED TO FIND IROTAT CARD IN PRMTOP FILE: %s\n",parm_data->filename);
        return INVALID_FORMAT;
      }
      // Skip comment lines
      prmtop_data[0]='%';
      while (prmtop_data[0]=='%') {
        fgetpos(fptr, &pos); // Get position of beginning of line in case this is a good line
        if (fgets (prmtop_data,BUFFER_SIZE,fptr) == NULL) {
          printf("*** ERROR IN readparm SUBROUTINE\n");
          printf("*** HIT EOF WHILE TRYING TO READ FORMAT CARD IN PRMTOP FILE: %s\n",parm_data->filename);
          return FILE_READ_FAIL;
        }
      }
      fsetpos(fptr, &pos); // Go to beginning of good line
    }

    for (i=0; i < parm_data->NTOTAT; ++i) {
      /*fscanf returns number of arguments read so should return 1*/
      if ( fscanf(fptr,"%d",&parm_data->atom[i].irotat)!=1 )
      {
        /*Fatal Error failed to find irotat*/
        printf("*** ERROR - FAILED TO FIND IROTAT INFO FOR ATOM %d\n",i);
        return INVALID_FORMAT;
      }
    }

    if (global_options->VERBOSITY>=HIGH) {
      printf("   Prmtop    (atoms): Read in atom irotat info.\n");
    }
  }

  /*We now have loads of other bits and pieces like solvent caps, box info, perturbation info, LES, etc etc.*/
  /*None of this is used by this program though so don't bother reading it*/
  fclose(fptr);
  free(prmtop_data);
  prmtop_data = NULL;

  return SUCCESS;
}
