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

/** @file process_prmtop.c
 * This routine is responsible for processing the prmtop file into seperate parameters for
 * each bond, angle and dihedral type.
 *
 * This is necessary because the prmtop file and thus prmtop storage arrays do not keep seperate
 * parameters for different bond / angle or dihedral setups that share the same values for these parameters
 *
 * e.g. O-C-N-H and O-C-N-CH3 share the same parameters and so this is stored in the same spot in memory
 * We need to split this in order to optimise them individually
 * @see read_prmtop.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "function_def.h"
#include "constants.h"

/**
 * Processes all the parmtop data.
 * Processes each individual parm_struct.
 * @param[in] global_options The global options structure, specifies number of prmtops
 * @param[in,out] parm_datas Array of parm structures to process
 * @return Integer indicating success or failure
 */
int process_prmtops(global_options_struct *global_options, parm_struct *parm_datas)
{
  int i; int retval;
  for (i=0; i<global_options->num_prmtops; ++i) {
    retval = process_single_prmtop(global_options, &parm_datas[i]);
    process_retval(retval, global_options->VERBOSITY);
  }
  if (verify_prmtops(global_options,parm_datas)!=SUCCESS) {
    printf("ERROR: prmtops do not contain same parameters to fit!\n");
    return INVALID_DATA;
  }
  return SUCCESS;
}

/**
 * Checks that each prmtop contains all of the parameters that are to be fit.
 * This is necessary for consistency in the fitting, as there is no point in doing a multi-molecule
 * fit if a parameter is only present in one of the molecules.
 * @param[in] global_options The global options structure
 * @param[in] parm_datas Array of parameter structure
 * @return Integer indicating success (prmtops okay) or failure (inconsistency detected or other error)
 */
int verify_prmtops(global_options_struct *global_options, parm_struct *parm_datas)
{
  int i,j,k,l;
  int nbondt,nanglet,ndihet;
  int nbondp,nanglep,ndihep;
  nbondp=nanglep=ndihep=0;
  nbondt=nanglet=ndihet=0;
  
  // Allocate memory for parameter comparison arrays- 2D array of pointers to parameters in prmtop
  bond_data_struct*** bonds = (bond_data_struct***)malloc(global_options->num_prmtops*sizeof(bond_data_struct**));
  angle_data_struct*** angles = (angle_data_struct***)malloc(global_options->num_prmtops*sizeof(angle_data_struct**));
  dihedral_type_struct*** dihes = (dihedral_type_struct***)malloc(global_options->num_prmtops*sizeof(dihedral_type_struct**));
  if ( !bonds || !angles || !dihes ) {
    printf("ERROR: Unable to allocate memory for parameter comparison arrays!\n");
    return ALLOC_FAIL;
  }
  
  // Pull out all parameters to be fit in each, expanding array as we go
  for (i=0; i<global_options->num_prmtops; ++i) {
    // Initial allocation so realloc will work okay
    bonds[i] = (bond_data_struct**)malloc(sizeof(bond_data_struct*));
    angles[i] = (angle_data_struct**)malloc(sizeof(angle_data_struct*));
    dihes[i] = (dihedral_type_struct**)malloc(sizeof(dihedral_type_struct));
    if (!bonds[i] || !angles[i] || !dihes[i]) {
      printf("ERROR: Unable to allocate memory for parameter comparison arrays!\n");
      return ALLOC_FAIL;
    }
    // Put in the bonds
    nbondt=0;
    for (j=0; j<parm_datas[i].unique_bonds_found; ++j)
      if (parm_datas[i].bond_data[j].DO_FIT_KR==YES || parm_datas[i].bond_data[j].DO_FIT_REQ==YES) {
        ++nbondt;
        if (nbondt>1) bonds[i] = (bond_data_struct**)realloc(bonds[i],nbondt*sizeof(bond_data_struct*));
        bonds[i][nbondt-1] = &parm_datas[i].bond_data[j];
      }
    // Put in the angles
    nanglet=0;
    for (j=0; j<parm_datas[i].unique_angles_found; ++j)
      if (parm_datas[i].angle_data[j].DO_FIT_KT==YES || parm_datas[i].angle_data[j].DO_FIT_THEQ==YES) {
        ++nanglet;
        if (nanglet>1) angles[i] = (angle_data_struct**)realloc(angles[i],nanglet*sizeof(angle_data_struct*));
        angles[i][nanglet-1] = &parm_datas[i].angle_data[j];
      }
    // Put in the dihedral types (terms will be examined later)
    ndihet=0;
    for (j=0; j<parm_datas[i].unique_dihedrals_found; ++j)
      for (k=0; k<parm_datas[i].dihedral_data[j].num_terms; ++k)
        if (parm_datas[i].dihedral_data[j].term[k].DO_FIT_KP==YES || parm_datas[i].dihedral_data[j].term[k].DO_FIT_NP==YES ||
            parm_datas[i].dihedral_data[j].term[k].DO_FIT_PHASE==YES) {
          ++ndihet;
          if (ndihet>1) dihes[i] = (dihedral_type_struct**)realloc(dihes[i],ndihet*sizeof(dihedral_type_struct*));
          dihes[i][ndihet-1] = &parm_datas[i].dihedral_data[j];
        }
    // Sanity check - number of terms to be fit is same in this structure as previous
    if ( (nbondp && nanglep && ndihep) ) {
      if (nbondt != nbondp) {
        printf("  ERROR: Different number of bond terms to be fit across prmtops!\n");
        printf("   Prmtops were %s and %s\n", parm_datas[i].filename, parm_datas[i-1].filename);
        printf("   Number of bond terms was %d and %d\n", nbondt, nbondp);
        return FAILURE;
      } else if (nanglet != nanglep) {
        printf("  ERROR: Different number of angle terms to be fit across prmtops!\n");
        printf("   Prmtops were %s and %s\n", parm_datas[i].filename, parm_datas[i-1].filename);
        printf("   Number of angle terms was %d and %d\n", nanglet, nanglep);
        return FAILURE;
      } else if (ndihet != ndihep) {
        printf("  ERROR: Different number of dihedral terms to be fit across prmtops!\n");
        printf("   Prmtops were %s and %s\n", parm_datas[i].filename, parm_datas[i-1].filename);
        printf("   Number of dihedral terms was %d and %d\n", nanglet, nanglep);
        return FAILURE;
      }
    } else {
      nbondp = nbondt;
      nanglep = nanglet;
      ndihep = ndihet;
    }
    // Sort these arrays by name using custom comparators so they ideally will all be identical in order
    if (nbondt>1) qsort(bonds[i], nbondt, sizeof(bond_data_struct*),bondcomparator);
    if (nanglet>1) qsort(angles[i], nanglet, sizeof(angle_data_struct*),anglecomparator);
    if (ndihet>1) qsort(dihes[i], ndihet, sizeof(dihedral_type_struct*),dihedralcomparator);
  }
  // Compare fit parameter collections
  short fail = NO;
  for (i=0; i<global_options->num_prmtops; ++i) {
    for (j=i+1; j<global_options->num_prmtops; ++j) {
      fail = NO;
      // Compare bonds
      for (k=0; k<nbondt; ++k) {
        if ( strcmp(bonds[i][k]->atom_type1,bonds[j][k]->atom_type1) || strcmp(bonds[i][k]->atom_type2,bonds[j][k]->atom_type2) ||
             bonds[i][k]->DO_FIT_KR != bonds[j][k]->DO_FIT_KR || bonds[i][k]->DO_FIT_REQ != bonds[j][k]->DO_FIT_REQ ) {
          fail=YES;
          break;
      } }
      // Compare angles
      if (fail==NO) {
        for (k=0; k<nanglet; ++k) {
          if ( strcmp(angles[i][k]->atom_type1,angles[j][k]->atom_type1) || strcmp(angles[i][k]->atom_type2,angles[j][k]->atom_type2) ||
               strcmp(angles[i][k]->atom_type3,angles[j][k]->atom_type3) || angles[i][k]->DO_FIT_KT != angles[j][k]->DO_FIT_KT ||
               angles[i][k]->DO_FIT_THEQ != angles[j][k]->DO_FIT_THEQ ) {
            fail=YES;
            break;
      } } }
      // Compare dihedrals- check the types are equal then check the fitting options for each term match
      if (fail==NO) {
        for (k=0; k<ndihet; ++k) {
          if (dihedral_types_equal(dihes[i][k], dihes[j][k])==NO || dihes[i][k]->num_terms != dihes[j][k]->num_terms) {
            fail=YES;
            break;
          } else {
            for (l=0; l<dihes[i][k]->num_terms; ++l) {
              if ( dihes[i][k]->term[l].DO_FIT_KP != dihes[j][k]->term[l].DO_FIT_KP || dihes[i][k]->term[l].DO_FIT_NP != dihes[j][k]->term[l].DO_FIT_NP ||
                   dihes[i][k]->term[l].DO_FIT_PHASE != dihes[j][k]->term[l].DO_FIT_PHASE ) {
                fail=YES;
                break;
      } } } } }
      // Print out parameters of failing prmtops
      if (fail==YES) {
        printf(" ERROR: Prmtop parameters to be fit do not match!\n");
        printf("  PRMTOP 1 = %s\n", parm_datas[i].filename);
        print_parameter_summary(global_options, &(parm_datas[i]));
        printf("\n  PRMTOP 2 = %s\n", parm_datas[i].filename);
        print_parameter_summary(global_options, &(parm_datas[j]));
        // Clean up memory
        for (i=0; i<parm_datas[0].unique_bonds_found; ++i) free(bonds[i]);
        for (i=0; i<parm_datas[0].unique_angles_found; ++i) free(angles[i]);
        for (i=0; i<parm_datas[0].unique_dihedrals_found; ++i) free(dihes[i]);
        free(bonds);
        free(angles);
        free(dihes);
        return FAILURE;
      }
    }
  }
  // Clean up memory
  for (i=0; i<global_options->num_prmtops; ++i) free(bonds[i]);
  for (i=0; i<global_options->num_prmtops; ++i) free(angles[i]);
  for (i=0; i<global_options->num_prmtops; ++i) free(dihes[i]);
  free(bonds);
  free(angles);
  free(dihes);
  return SUCCESS;
}

/**
 * Compares names of two bonds, for qsort.
 */
int bondcomparator(const void *a, const void *b)
{
  bond_data_struct** b1 = (bond_data_struct**)a;
  bond_data_struct** b2 = (bond_data_struct**)b;
  
  char name1[NAME_SIZE*2+1+1]; // XXXX-XXXX\0
  char name2[NAME_SIZE*2+1+1]; // = 2 atoms, 1 dash, null terminator
  sprintf(name1,"%s-%s",(*b1)->atom_type1,(*b1)->atom_type2);
  sprintf(name2,"%s-%s",(*b2)->atom_type1,(*b2)->atom_type2);
  return strcmp(name1,name2);
}

/**
 * Compares names of two angles, for qsort
 */
int anglecomparator(const void *a, const void *b)
{
  angle_data_struct** a1 = (angle_data_struct**)a;
  angle_data_struct** a2 = (angle_data_struct**)b;
  
  char name1[NAME_SIZE*3+2+1]; // 3 atoms, 2 dashes, null terminator
  char name2[NAME_SIZE*3+2+1];
  sprintf(name1,"%s-%s-%s",(*a1)->atom_type1,(*a1)->atom_type2,(*a1)->atom_type3);
  sprintf(name2,"%s-%s-%s",(*a2)->atom_type1,(*a2)->atom_type2,(*a2)->atom_type3);
  return strcmp(name1, name2);
}

/**
 * Compares names of two dihedrals, for qsort
 */
int dihedralcomparator(const void *a, const void *b)
{
  dihedral_type_struct **d1 = (dihedral_type_struct**)a;
  dihedral_type_struct **d2 = (dihedral_type_struct**)b;
  
  char name1[NAME_SIZE*4+3+1]; // 4 atoms, 3 dashes, null terminator
  char name2[NAME_SIZE*4+3+1];
  sprintf(name1,"%s-%s-%s-%s",(*d1)->atom_type1,(*d1)->atom_type2,(*d1)->atom_type3,(*d1)->atom_type4);
  sprintf(name2,"%s-%s-%s-%s",(*d2)->atom_type1,(*d2)->atom_type2,(*d2)->atom_type3,(*d2)->atom_type4);
  return strcmp(name1,name2);
}

/**
 * Processes raw prmtop data into arrays of unique bond, angle, and dihedral parameters that are
 * in structures that are much easier to optimize. Done on one prmtop at a time
 * @see read_prmtops
 * @see process_prmtops
 * 
 * @param[in,out] global_options The global options structure. unique_[bonds,angles,dihedrals]_found will be updated
 * @param[in,out] parm_data Contains the raw prmtop data, will also be populated with processed prmtop data
 * @return Integer indicating success or failure
 */
int process_single_prmtop(global_options_struct *global_options, parm_struct *parm_data)
{
  int i,j,k;
  int atom1, atom2, atom3, atom4;
  double rk, req, tk, teq, pk, pn, phase;
  char type1[NAME_SIZE];
  char type2[NAME_SIZE];
  char type3[NAME_SIZE];
  char type4[NAME_SIZE];    
  short int is_unique;
  short int is_improper;
  short int is_new_term;
  char bond_make_up1[(NAME_SIZE*2)-1];
  char bond_make_up2[(NAME_SIZE*2)-1];
  char angle_make_up1[(NAME_SIZE*3)-2];
  char angle_make_up2[(NAME_SIZE*3)-2];
  char dihedral_make_up1[(NAME_SIZE*4)-3];
  char dihedral_make_up2[(NAME_SIZE*4)-3];    
  int check1,check2;
  char response = '\0';

  /* Ask the user some general categories so as to overwhelm them with less output */
  short int prompt_bonds=YES, prompt_angles=YES, prompt_dihedrals=YES;
  short int b_req=NO, b_kr=NO;
  short int a_kt=NO, a_theq=NO;
  short int d_kp=NO, d_np=NO, d_phase=NO;
  
  /* Make sure a file name has been specified to save in if necessary */
  if (global_options->RUNTYPE==SET_PARAMS && global_options->PARAMETER_FILE_NAME==NULL)
  {
    printf("\n\n!  ERROR: No location specified to save parameter file. Please\n");
    printf("          read the documentation on how to use this runtype and   \n");
    printf("          specify a PARAMETER_FILE_NAME.\n");
    return FAILURE;
  }
  /* If a filename has been specified and it's not save, warn the user of ambiguity */
  if ( global_options->PARAMETERS_TO_FIT==DEFAULT && global_options->PARAMETER_FILE_NAME!=NULL && 
       global_options->RUNTYPE != SET_PARAMS && global_options->VERBOSITY >=MEDIUM )
  {
      printf("! Ambiguous parameters to fit options- parameter file and DEFAULT parameters specified.\n");
      return FAILURE;
  }
  
  /* If a filename has been specified and it's K only, warn the user */
  if (global_options->PARAMETERS_TO_FIT==K_ONLY && global_options->PARAMETER_FILE_NAME!=NULL)
  {
    printf("! Ambiguous parameters to fit options- parameter file and K_ONLY parameters specified.\n");
    return FAILURE;
  }
  
  if (global_options->RUNTYPE==SET_PARAMS)
  {
    global_options->PARAMETERS_TO_FIT = SAVE;
    
    printf("  Would you like to optimise any bond parameters? (y/n):\n");
    while( response!='y' && response!='n' )
      scanf("%c",&response);
    prompt_bonds = (response == 'y') ? YES : NO;
    response = '\0';
    fflush(stdout);

    if (prompt_bonds == YES)
    {
      printf("    - Any bond REQ? (y/n): \n");
      while(response != 'y' && response != 'n')
        scanf("%c", &response);
      b_req = (response == 'y') ? YES : NO;
      response = '\0';
      fflush(stdout);
      
      printf("    - Any bond KR? (y/n): \n");
      while(response != 'y' && response != 'n')
        scanf("%c", &response);
      b_kr = (response == 'y') ? YES : NO;
      response = '\0';
      fflush(stdout);
    }
    
    printf("  Would you like to optimise any angle parameters? (y/n): \n");
    fflush(stdout);
    while( response!='y' && response!='n' )
      scanf("%c",&response);
    prompt_angles = (response == 'y') ? YES : NO;
    response = '\0';
    if (prompt_angles == YES)
    {
      printf("    - Any angle KT? (y/n): \n");
      while(response != 'y' && response != 'n')
        scanf("%c", &response);
      a_kt = (response == 'y') ? YES : NO;
      response = '\0';
      fflush(stdout);
      
      printf("    - Any angle THEQ? (y/n): \n");
      while(response != 'y' && response != 'n')
        scanf("%c", &response);
      a_theq = (response == 'y') ? YES : NO;
      response = '\0';
      fflush(stdout);
    }
   
   printf("  Would you like to optimise any dihedral parameters? (y/n): \n");
    while( response!='y' && response!='n' )
      scanf("%c",&response);
    prompt_dihedrals = (response == 'y') ? YES : NO;    
    fflush(stdout);
    
    response = '\0';
    if (prompt_dihedrals == YES)
    {
      printf("    - Any dihedral KP? (y/n): \n");
      while(response != 'y' && response != 'n')
        scanf("%c", &response);
      d_kp = (response == 'y') ? YES : NO;
      fflush(stdout);
      
      response = '\0';
      printf("    - Any dihedral NP? (y/n): \n" );
      while(response != 'y' && response != 'n')
        scanf("%c", &response);
      d_np = (response == 'y') ? YES : NO;
      fflush(stdout);
      
      response = '\0';
      printf("    - Any dihedral PHASE? (y/n): \n");
      while(response != 'y' && response != 'n')
        scanf("%c", &response);
      d_phase = (response == 'y') ? YES : NO;
      response = '\0';
      fflush(stdout);
    }
  }

  // Set bond, angle, dihedral data structs initially to NULL so realloc acts like malloc the first time
  parm_data->bond_data=NULL;
  parm_data->angle_data=NULL;
  parm_data->dihedral_data=NULL;

  if (global_options->VERBOSITY>=HIGH)
    printf("   Prmtop   (unique): Processing prmtop to find unique bonds, angles & dihedrals.\n");
  
  parm_data->unique_bonds_found=0;
  parm_data->unique_angles_found=0;
  parm_data->unique_dihedrals_found=0;    
  parm_data->unique_dihedral_terms=0;    
  /*First of all we will loop over all the bonds and split them into their unique types*/
  /*We will go through each bond in turn and assign it to a location in bond_data. The problem
    we have is that in theory two bond terms could exist that have the same types but different
    parameters, this is not an accute problem for bonds and angles but is for dihedrals*/
  
  for (i=0;i<parm_data->NBONH;++i)
  {
     /*Loop over all bonds involving hydrogen*/
     /*Get the atoms involved in the bond we are looking at*/
     atom1=unObfuscateAtom(parm_data->pbondH[i].ib);
     atom2=unObfuscateAtom(parm_data->pbondH[i].jb);
     /*Get the parameters involved*/
     rk=parm_data->rk[parm_data->pbondH[i].icb-1];
     req=parm_data->req[parm_data->pbondH[i].icb-1];     
     /*Get the atom types*/
     strcpy(type1,parm_data->atom[atom1-1].isymbl);
     strcpy(type2,parm_data->atom[atom2-1].isymbl);
     
     /*Ok, now we have the data for this bond*/
     /*Now we need to determine if it is unique*/
     /*loop over all unique bonds to date and see if we get a match*/
     is_unique=YES;
     for (j=0;j<parm_data->unique_bonds_found;++j)
     {
       /*A bond is unique if both atoms types don't match, and both parameters don't match*/
       /*Note, we can do the parameters in both directions*/
       /*First try the bond one way, then the reverse since directionality is irrelevant*/
       strcpy(bond_make_up1,type1);
       strcat(bond_make_up1,type2);
       strcpy(bond_make_up2,parm_data->bond_data[j].atom_type1);
       strcat(bond_make_up2,parm_data->bond_data[j].atom_type2);
       check1=!strcmp(bond_make_up1,bond_make_up2);             /*strcmp returns zero on a match - so check1 will be 1 if matched*/
       strcpy(bond_make_up1,type2);
       strcat(bond_make_up1,type1);
       strcpy(bond_make_up2,parm_data->bond_data[j].atom_type1);
       strcat(bond_make_up2,parm_data->bond_data[j].atom_type2);
       check2=!strcmp(bond_make_up1,bond_make_up2);             /*strcmp returns zero on a match - so check2 will be 1 if reverse matched*/
       /*if the atom type is not unique and matches here see if the parameters also match
         if they do then we have found the matching bond type for this bond so add the data and break
         else we keep looping*/
       if (check1 || check2)
       {
         if (rk == parm_data->bond_data[j].rk && req == parm_data->bond_data[j].req)
         {
           /*bond type is not unique*/
           is_unique=NO;
           /*fill in the data*/
           ++(parm_data->bond_data[j].number);
           /*WE NEED TO CHECK HERE THAT WE HAVEN'T OVERFLOWED OUR ARRAY*/
           if (parm_data->bond_data[j].number>=MAX_BONDS_PER_TYPE)
           {
             printf("*** ERROR IN PROCESS_PRMTOP - MAX_BONDS_PER_TYPE OF %d\n",MAX_BONDS_PER_TYPE);
             printf("*** EXCEEDED FOR parm_data->bond_data[j], j = %d\n",j);
             printf("*** EDIT MAX_BONDS_PER_TYPE IN prmtop_params.h AND RECOMPILE.\n");
             return DATA_OVERFLOW;
           }
           parm_data->bond_data[j].atom1[(parm_data->bond_data[j].number)-1]=atom1;
           parm_data->bond_data[j].atom2[(parm_data->bond_data[j].number)-1]=atom2;
           
           break; /*Quit the loop, no point doing the rest if we have matched at this point*/
         }
       }
     } /*essentially if we finish this loop without a hit it must be unique*/
     if (is_unique==YES)
     {
       /*Bond is unique, have to extend number of bond_data structures by 1*/
       ++(parm_data->unique_bonds_found);
       
       /*reallocate memory for the bond_data*/
        if (global_options->VERBOSITY>=HIGH)
           printf("   Reallocating %d bytes for parm_data->*bond_data\n",(int) sizeof(struct _bond_data_struct)*parm_data->unique_bonds_found);
        parm_data->bond_data = (struct _bond_data_struct *) realloc(parm_data->bond_data,sizeof(struct _bond_data_struct)*parm_data->unique_bonds_found);
        if (parm_data->bond_data == NULL) {
          malloc_failure_char("process_prmtop", "parm_data->bond_data", sizeof(struct _bond_data_struct)*parm_data->unique_bonds_found);
          return ALLOC_FAIL;
        }

       parm_data->mem_allocated+=sizeof(struct _bond_data_struct);
       /*now put the data in the bond_structure*/
       parm_data->bond_data[parm_data->unique_bonds_found-1].number=1;
       strcpy(parm_data->bond_data[parm_data->unique_bonds_found-1].atom_type1,type1);
       strcpy(parm_data->bond_data[parm_data->unique_bonds_found-1].atom_type2,type2);
       parm_data->bond_data[parm_data->unique_bonds_found-1].rk=rk;
       parm_data->bond_data[parm_data->unique_bonds_found-1].req=req;
       parm_data->bond_data[parm_data->unique_bonds_found-1].atom1[0]=atom1;
       parm_data->bond_data[parm_data->unique_bonds_found-1].atom2[0]=atom2;
       
       /*Need to check if we are to fit this parameter or not*/
       parm_data->bond_data[parm_data->unique_bonds_found-1].DO_FIT_KR=NO;
       parm_data->bond_data[parm_data->unique_bonds_found-1].DO_FIT_REQ=NO;
       if ( global_options->RUNTYPE==SET_PARAMS && prompt_bonds ) {
        /*prompt the user whether to fit or not*/
        if (b_kr==YES) {
          response='\0';
          printf("Fit Parameter: BOND (%s-%s) KR? (y/n): ",
          parm_data->bond_data[parm_data->unique_bonds_found-1].atom_type1,
          parm_data->bond_data[parm_data->unique_bonds_found-1].atom_type2);
          fflush(stdout);
          while( response!='y' && response!='n' )
            scanf("%c",&response);
          if (response=='y')
              parm_data->bond_data[parm_data->unique_bonds_found-1].DO_FIT_KR=YES;
        }
        if (b_req==YES) {
          response='\0';
          printf("Fit Parameter: BOND (%s-%s) REQ? (y/n): ",
                parm_data->bond_data[parm_data->unique_bonds_found-1].atom_type1,
                parm_data->bond_data[parm_data->unique_bonds_found-1].atom_type2);
          fflush(stdout);
          while( response!='y' && response!='n' )
            scanf("%c",&response);
          if (response=='y')
            parm_data->bond_data[parm_data->unique_bonds_found-1].DO_FIT_REQ=YES;
        }
      }
      else if (global_options->PARAMETERS_TO_FIT==DEFAULT) {
        parm_data->bond_data[parm_data->unique_bonds_found-1].DO_FIT_KR=YES;  /*DEFAULT IS TO FIT EVERYTHING, this will be adjusted as necessary by other routines*/
        parm_data->bond_data[parm_data->unique_bonds_found-1].DO_FIT_REQ=YES;
      } else if (global_options->PARAMETERS_TO_FIT==K_ONLY) {
        parm_data->bond_data[parm_data->unique_bonds_found-1].DO_FIT_KR=NO;
        parm_data->bond_data[parm_data->unique_bonds_found-1].DO_FIT_REQ=NO;
      }
    }
  } /*End of bonds with Hydrogen*/

  /*Now do bonds without hydrogen*/
  for (i=0;i<parm_data->NBONA;++i)
  {
    /*Get the atoms involved in the bond we are looking at*/
    atom1=unObfuscateAtom(parm_data->pbond[i].ib);
    atom2=unObfuscateAtom(parm_data->pbond[i].jb);
    /*Get the parameters involved*/
    rk=parm_data->rk[parm_data->pbond[i].icb-1];
    req=parm_data->req[parm_data->pbond[i].icb-1];
    /*Get the atom types*/
    strcpy(type1,parm_data->atom[atom1-1].isymbl);
    strcpy(type2,parm_data->atom[atom2-1].isymbl);
    
    /*Ok, now we have the data for this bond*/
    /*Now we need to determine if it is unique*/
    /*loop over all unique bonds to date and see if we get a match*/
    is_unique=YES;
    for (j=0;j<parm_data->unique_bonds_found;++j) {
      /*A bond is unique if both atoms types don't match, and both parameters don't match*/
      /*Note, we can do the parameters in both directions*/
      /*First try the bond one way, then the reverse since directionality is irrelevant*/
      strcpy(bond_make_up1,type1);
      strcat(bond_make_up1,type2);
      strcpy(bond_make_up2,parm_data->bond_data[j].atom_type1);
      strcat(bond_make_up2,parm_data->bond_data[j].atom_type2);
      check1=!strcmp(bond_make_up1,bond_make_up2);             /*strcmp returns zero on a match - so check1 will be 1 if matched*/
      strcpy(bond_make_up1,type2);
      strcat(bond_make_up1,type1);
      strcpy(bond_make_up2,parm_data->bond_data[j].atom_type1);
      strcat(bond_make_up2,parm_data->bond_data[j].atom_type2);
      check2=!strcmp(bond_make_up1,bond_make_up2);             /*strcmp returns zero on a match - so check2 will be 1 if reverse matched*/
      /*if the atom type is not unique and matches here see if the parameters also match
        if they do then we have found the matching bond type for this bond so add the data and break
        else we keep looping*/
      if (check1 || check2) {
        if (rk == parm_data->bond_data[j].rk && req == parm_data->bond_data[j].req) {
          /*bond type is not unique*/
          is_unique=NO;
          /*fill in the data*/
          ++(parm_data->bond_data[j].number);
          /*WE NEED TO CHECK HERE THAT WE HAVEN'T OVERFLOWED OUR ARRAY*/
          if (parm_data->bond_data[j].number>=MAX_BONDS_PER_TYPE) {
            printf("*** ERROR IN PROCESS_PRMTOP - MAX_BONDS_PER_TYPE OF %d\n",MAX_BONDS_PER_TYPE);
            printf("*** EXCEEDED FOR parm_data->bond_data[j], j = %d\n",j);
            printf("*** EDIT MAX_BONDS_PER_TYPE IN prmtop_params.h AND RECOMPILE.\n");
            return DATA_OVERFLOW;
          }
          parm_data->bond_data[j].atom1[(parm_data->bond_data[j].number)-1]=atom1;
          parm_data->bond_data[j].atom2[(parm_data->bond_data[j].number)-1]=atom2;
           
          break; /*Quit the loop, no point doing the rest if we have matched at this point*/
        }
      }
    } /*essentially if we finish this loop without a hit it must be unique*/
    if (is_unique==YES) {
      /*Bond is unique, have to extend number of bond_data structures by 1*/
      ++(parm_data->unique_bonds_found);
      /*reallocate memory for the bond_data*/
      if (global_options->VERBOSITY>=HIGH)
        printf("   Reallocating %d bytes for parm_data->*bond_data\n",(int) (sizeof(struct _bond_data_struct)*parm_data->unique_bonds_found));
      parm_data->bond_data = (struct _bond_data_struct *) realloc(parm_data->bond_data,sizeof(struct _bond_data_struct)*parm_data->unique_bonds_found);
      if (parm_data->bond_data == NULL) {
        malloc_failure_char("process_prmtop", "parm_data->bond_data", sizeof(struct _bond_data_struct)*parm_data->unique_bonds_found);
        return ALLOC_FAIL;
      }
      parm_data->mem_allocated+=sizeof(struct _bond_data_struct);
      /*now put the data in the bond_structure*/
      parm_data->bond_data[parm_data->unique_bonds_found-1].number=1;
      strcpy(parm_data->bond_data[parm_data->unique_bonds_found-1].atom_type1,type1);
      strcpy(parm_data->bond_data[parm_data->unique_bonds_found-1].atom_type2,type2);
      parm_data->bond_data[parm_data->unique_bonds_found-1].rk=rk;
      parm_data->bond_data[parm_data->unique_bonds_found-1].req=req;
      parm_data->bond_data[parm_data->unique_bonds_found-1].atom1[0]=atom1;
      parm_data->bond_data[parm_data->unique_bonds_found-1].atom2[0]=atom2;
      
      /*Need to check if we are to fit this parameter or not*/
      parm_data->bond_data[parm_data->unique_bonds_found-1].DO_FIT_REQ=NO;
      parm_data->bond_data[parm_data->unique_bonds_found-1].DO_FIT_KR=NO;
      if (global_options->RUNTYPE==SET_PARAMS) {
        /*prompt the user whether to fit or not*/
        if (b_kr == YES) {
          response='\0';
          printf("Fit Parameter: BOND (%s-%s) KR? (y/n): ",
          parm_data->bond_data[parm_data->unique_bonds_found-1].atom_type1,
          parm_data->bond_data[parm_data->unique_bonds_found-1].atom_type2);
          fflush(stdout);
          while( response!='y' && response!='n' )
            scanf("%c",&response);
          if (response=='y')
            parm_data->bond_data[parm_data->unique_bonds_found-1].DO_FIT_KR=YES;
        }
        if (b_req==YES) {
          response='\0';
          printf("Fit Parameter: BOND (%s-%s) REQ? (y/n): ",
                parm_data->bond_data[parm_data->unique_bonds_found-1].atom_type1,
                parm_data->bond_data[parm_data->unique_bonds_found-1].atom_type2);
          fflush(stdout);
          while( response!='y' && response!='n' )
            scanf("%c",&response);
          if (response=='y')
            parm_data->bond_data[parm_data->unique_bonds_found-1].DO_FIT_REQ=YES;
        }
      } else if (global_options->PARAMETERS_TO_FIT==DEFAULT) {
        parm_data->bond_data[parm_data->unique_bonds_found-1].DO_FIT_KR=YES;  /*DEFAULT IS TO FIT EVERYTHING, this will be adjusted as necessary by other routines*/
        parm_data->bond_data[parm_data->unique_bonds_found-1].DO_FIT_REQ=YES;
      } else if (global_options->PARAMETERS_TO_FIT==K_ONLY) {
        parm_data->bond_data[parm_data->unique_bonds_found-1].DO_FIT_KR=NO;
        parm_data->bond_data[parm_data->unique_bonds_found-1].DO_FIT_REQ=NO;
      }
    }
  } /*End of bonds withOUT Hydrogen*/
  
  /*Now we do angles with hydrogen*/
  for (i=0;i<parm_data->NTHETH;++i)
  {
    /*Loop over all angles involving hydrogen*/
    /*Get the atoms involved in the angle we are looking at*/
    atom1=unObfuscateAtom(parm_data->pangleH[i].it);
    atom2=unObfuscateAtom(parm_data->pangleH[i].jt);
    atom3=unObfuscateAtom(parm_data->pangleH[i].kt);
    /*Get the parameters involved*/
    tk=parm_data->tk[parm_data->pangleH[i].ict-1];
    teq=parm_data->teq[parm_data->pangleH[i].ict-1];
    /*Get the atom types*/
    strcpy(type1,parm_data->atom[atom1-1].isymbl);
    strcpy(type2,parm_data->atom[atom2-1].isymbl);
    strcpy(type3,parm_data->atom[atom3-1].isymbl);
    
    /*Ok, now we have the data for this angle*/
    /*Now we need to determine if it is unique*/
    /*loop over all unique angles to date and see if we get a match*/
    is_unique=YES;
    for (j=0;j<parm_data->unique_angles_found;++j) {
      /*An angle is unique if all atoms types don't match, and both parameters don't match*/
      /*Note, we can do the atom types in both directions*/
      /*First try the angle one way, then the reverse since directionality is irrelevant*/
      strcpy(angle_make_up1,type1);
      strcat(angle_make_up1,type2);
      strcat(angle_make_up1,type3);       
      strcpy(angle_make_up2,parm_data->angle_data[j].atom_type1);
      strcat(angle_make_up2,parm_data->angle_data[j].atom_type2);
      strcat(angle_make_up2,parm_data->angle_data[j].atom_type3);       
      check1=!strcmp(angle_make_up1,angle_make_up2);             /*strcmp returns zero on a match - so check1 will be 1 if matched*/
      strcpy(angle_make_up1,type3);
      strcat(angle_make_up1,type2);
      strcat(angle_make_up1,type1);
      strcpy(angle_make_up2,parm_data->angle_data[j].atom_type1);
      strcat(angle_make_up2,parm_data->angle_data[j].atom_type2);
      strcat(angle_make_up2,parm_data->angle_data[j].atom_type3);
      check2=!strcmp(angle_make_up1,angle_make_up2);             /*strcmp returns zero on a match - so check2 will be 1 if reverse matched*/
      /*if the atom type is not unique and matches here see if the parameters also match
        if they do then we have found the matching angle type for this angle so add the data and break
        else we keep looping*/
      if (check1 || check2) {
        if (tk == parm_data->angle_data[j].tk && teq == parm_data->angle_data[j].teq) {
          /*angle type is not unique*/
          is_unique=NO;
          /*fill in the data*/
          ++(parm_data->angle_data[j].number);
          /*WE NEED TO CHECK HERE THAT WE HAVEN'T OVERFLOWED OUR ARRAY*/
          if (parm_data->angle_data[j].number>=MAX_ANGLES_PER_TYPE) {
            printf("*** ERROR IN PROCESS_PRMTOP - MAX_ANGLES_PER_TYPE OF %d\n",MAX_ANGLES_PER_TYPE);
            printf("*** EXCEEDED FOR parm_data->angle_data[j], j = %d\n",j);
            printf("*** EDIT MAX_ANGLES_PER_TYPE IN prmtop_params.h AND RECOMPILE.\n");
            return DATA_OVERFLOW;
          }
          parm_data->angle_data[j].atom1[(parm_data->angle_data[j].number)-1]=atom1;
          parm_data->angle_data[j].atom2[(parm_data->angle_data[j].number)-1]=atom2;
          parm_data->angle_data[j].atom3[(parm_data->angle_data[j].number)-1]=atom3;
          
          break; /*Quit the loop, no point doing the rest if we have matched at this point*/
        }
      }
    } /*essentially if we finish this loop without a hit it must be unique*/
    if (is_unique==YES) {
      /*Angle is unique, have to extend number of angle_data structures by 1*/
      ++(parm_data->unique_angles_found);
      /*reallocate memory for the angle_data*/
      if (global_options->VERBOSITY>=HIGH)
        printf("   Reallocating %d bytes for parm_data->*angle_data\n",(int) (sizeof(struct _angle_data_struct)*parm_data->unique_angles_found));
      parm_data->angle_data = (struct _angle_data_struct *) realloc(parm_data->angle_data,sizeof(struct _angle_data_struct)*parm_data->unique_angles_found);
      if (parm_data->angle_data == NULL) {
        malloc_failure_char("process_prmtop", "parm_data->angle_data", sizeof(struct _angle_data_struct)*parm_data->unique_angles_found);
        return ALLOC_FAIL;
      }
      parm_data->mem_allocated+=sizeof(struct _angle_data_struct);
      /*now put the data in the angle_structure*/
      parm_data->angle_data[parm_data->unique_angles_found-1].number=1;
      strcpy(parm_data->angle_data[parm_data->unique_angles_found-1].atom_type1,type1);
      strcpy(parm_data->angle_data[parm_data->unique_angles_found-1].atom_type2,type2);
      strcpy(parm_data->angle_data[parm_data->unique_angles_found-1].atom_type3,type3);
      parm_data->angle_data[parm_data->unique_angles_found-1].tk=tk;
      parm_data->angle_data[parm_data->unique_angles_found-1].teq=teq;
      parm_data->angle_data[parm_data->unique_angles_found-1].atom1[0]=atom1;
      parm_data->angle_data[parm_data->unique_angles_found-1].atom2[0]=atom2;
      parm_data->angle_data[parm_data->unique_angles_found-1].atom3[0]=atom3;
      
      /*Need to check if we are to fit this parameter or not*/
      parm_data->angle_data[parm_data->unique_angles_found-1].DO_FIT_THEQ=NO;
      parm_data->angle_data[parm_data->unique_angles_found-1].DO_FIT_KT=NO;
      if (global_options->RUNTYPE==SET_PARAMS) {
        /*prompt the user whether to fit or not*/
        if (a_kt == YES) {
          response='\0';
          printf("Fit Parameter: ANGLE (%s-%s-%s) KT? (y/n): ",
          parm_data->angle_data[parm_data->unique_angles_found-1].atom_type1,
          parm_data->angle_data[parm_data->unique_angles_found-1].atom_type2,
          parm_data->angle_data[parm_data->unique_angles_found-1].atom_type3);
          fflush(stdout);
          while( response!='y' && response!='n' )
            scanf("%c",&response);
          if (response=='y')
            parm_data->angle_data[parm_data->unique_angles_found-1].DO_FIT_KT=YES;
        }
        if (a_theq == YES) {
          response='\0';
          printf("Fit Parameter: ANGLE (%s-%s-%s) THEQ? (y/n): ",
                 parm_data->angle_data[parm_data->unique_angles_found-1].atom_type1,
                 parm_data->angle_data[parm_data->unique_angles_found-1].atom_type2,
                 parm_data->angle_data[parm_data->unique_angles_found-1].atom_type3);
          fflush(stdout);
          while( response!='y' && response!='n' )
          {
            scanf("%c",&response);
          }
          if (response=='y')
            parm_data->angle_data[parm_data->unique_angles_found-1].DO_FIT_THEQ=YES;
        }
      } else if (global_options->PARAMETERS_TO_FIT==DEFAULT) {
        parm_data->angle_data[parm_data->unique_angles_found-1].DO_FIT_KT=YES; 
        parm_data->angle_data[parm_data->unique_angles_found-1].DO_FIT_THEQ=YES;
      } else if (global_options->PARAMETERS_TO_FIT==K_ONLY) {
        parm_data->angle_data[parm_data->unique_angles_found-1].DO_FIT_KT=NO; 
        parm_data->angle_data[parm_data->unique_angles_found-1].DO_FIT_THEQ=NO;
      }
    }
  } /*End of angles with Hydrogen*/

  /*Now we do angles withOUT hydrogen*/
  for (i=0;i<parm_data->NTHETA;++i) {
    /*Loop over all angles not involving hydrogen*/
    /*Get the atoms involved in the angle we are looking at*/
    atom1=unObfuscateAtom(parm_data->pangle[i].it);
    atom2=unObfuscateAtom(parm_data->pangle[i].jt);
    atom3=unObfuscateAtom(parm_data->pangle[i].kt);
    /*Get the parameters involved*/
    tk=parm_data->tk[parm_data->pangle[i].ict-1];
    teq=parm_data->teq[parm_data->pangle[i].ict-1];
    /*Get the atom types*/
    strcpy(type1,parm_data->atom[atom1-1].isymbl);
    strcpy(type2,parm_data->atom[atom2-1].isymbl);
    strcpy(type3,parm_data->atom[atom3-1].isymbl);

    /*Ok, now we have the data for this angle*/
    /*Now we need to determine if it is unique*/
    /*loop over all unique angles to date and see if we get a match*/
    is_unique=YES;
    for (j=0;j<parm_data->unique_angles_found;++j) {
      /*An angle is unique if all atoms types don't match, and both parameters don't match*/
      /*Note, we can do the atom types in both directions*/
      /*First try the angle one way, then the reverse since directionality is irrelevant*/
      strcpy(angle_make_up1,type1);
      strcat(angle_make_up1,type2);
      strcat(angle_make_up1,type3);
      strcpy(angle_make_up2,parm_data->angle_data[j].atom_type1);
      strcat(angle_make_up2,parm_data->angle_data[j].atom_type2);
      strcat(angle_make_up2,parm_data->angle_data[j].atom_type3);
      check1=!strcmp(angle_make_up1,angle_make_up2);             /*strcmp returns zero on a match - so check1 will be 1 if matched*/
      strcpy(angle_make_up1,type3);
      strcat(angle_make_up1,type2);
      strcat(angle_make_up1,type1);
      strcpy(angle_make_up2,parm_data->angle_data[j].atom_type1);
      strcat(angle_make_up2,parm_data->angle_data[j].atom_type2);
      strcat(angle_make_up2,parm_data->angle_data[j].atom_type3);
      check2=!strcmp(angle_make_up1,angle_make_up2);             /*strcmp returns zero on a match - so check2 will be 1 if reverse matched*/
      /*if the atom type is not unique and matches here see if the parameters also match
        if they do then we have found the matching angle type for this angle so add the data and break
        else we keep looping*/
      if (check1 || check2) {
        if (tk == parm_data->angle_data[j].tk && teq == parm_data->angle_data[j].teq) {
          /*angle type is not unique*/
          is_unique=NO;
          /*fill in the data*/
          ++(parm_data->angle_data[j].number);
          /*WE NEED TO CHECK HERE THAT WE HAVEN'T OVERFLOWED OUR ARRAY*/
          if (parm_data->angle_data[j].number>=MAX_ANGLES_PER_TYPE) {
            printf("*** ERROR IN PROCESS_PRMTOP - MAX_ANGLES_PER_TYPE OF %d\n",MAX_ANGLES_PER_TYPE);
            printf("*** EXCEEDED FOR parm_data->angle_data[j], j = %d\n",j);
            printf("*** EDIT MAX_ANGLES_PER_TYPE IN prmtop_params.h AND RECOMPILE.\n");
            return DATA_OVERFLOW;
          }
          parm_data->angle_data[j].atom1[(parm_data->angle_data[j].number)-1]=atom1;
          parm_data->angle_data[j].atom2[(parm_data->angle_data[j].number)-1]=atom2;
          parm_data->angle_data[j].atom3[(parm_data->angle_data[j].number)-1]=atom3;
          
          break; /*Quit the loop, no point doing the rest if we have matched at this point*/
        }
      }
    } /*essentially if we finish this loop without a hit it must be unique*/
    if (is_unique==YES) {
      /*Angle is unique, have to extend number of angle_data structures by 1*/
      ++(parm_data->unique_angles_found);
      /*reallocate memory for the angle_data*/
      if (global_options->VERBOSITY>=HIGH)
        printf("   Reallocating %d bytes for parm_data->*angle_data\n",(int) (sizeof(struct _angle_data_struct)*parm_data->unique_angles_found));
      parm_data->angle_data = (struct _angle_data_struct *) realloc(parm_data->angle_data,sizeof(struct _angle_data_struct)*parm_data->unique_angles_found);
      if (parm_data->angle_data == NULL) {
        malloc_failure_char("process_prmtop", "parm_data->angle_data", sizeof(struct _angle_data_struct)*parm_data->unique_angles_found);
        return ALLOC_FAIL;
      }
      parm_data->mem_allocated+=sizeof(struct _angle_data_struct);
      /*now put the data in the angle_structure*/
      parm_data->angle_data[parm_data->unique_angles_found-1].number=1;
      strcpy(parm_data->angle_data[parm_data->unique_angles_found-1].atom_type1,type1);
      strcpy(parm_data->angle_data[parm_data->unique_angles_found-1].atom_type2,type2);
      strcpy(parm_data->angle_data[parm_data->unique_angles_found-1].atom_type3,type3);
      parm_data->angle_data[parm_data->unique_angles_found-1].tk=tk;
      parm_data->angle_data[parm_data->unique_angles_found-1].teq=teq;
      parm_data->angle_data[parm_data->unique_angles_found-1].atom1[0]=atom1;
      parm_data->angle_data[parm_data->unique_angles_found-1].atom2[0]=atom2;
      parm_data->angle_data[parm_data->unique_angles_found-1].atom3[0]=atom3;
      
      /*Need to check if we are to fit this parameter or not*/
      parm_data->angle_data[parm_data->unique_angles_found-1].DO_FIT_THEQ=NO;
      parm_data->angle_data[parm_data->unique_angles_found-1].DO_FIT_KT=NO;
      if (global_options->RUNTYPE==SET_PARAMS) {
        /*prompt the user whether to fit or not*/
        if (a_kt==YES) {
          response='\0';
          printf("Fit Parameter: ANGLE (%s-%s-%s) KT? (y/n): ",
          parm_data->angle_data[parm_data->unique_angles_found-1].atom_type1,
          parm_data->angle_data[parm_data->unique_angles_found-1].atom_type2,
          parm_data->angle_data[parm_data->unique_angles_found-1].atom_type3);
          fflush(stdout);
          while( response!='y' && response!='n' )
            scanf("%c",&response);
          if (response=='y')
            parm_data->angle_data[parm_data->unique_angles_found-1].DO_FIT_KT=YES;
        }
        if (a_theq==YES) {
          response='\0';
          printf("Fit Parameter: ANGLE (%s-%s-%s) THEQ? (y/n): ",
                 parm_data->angle_data[parm_data->unique_angles_found-1].atom_type1,
                 parm_data->angle_data[parm_data->unique_angles_found-1].atom_type2,
                 parm_data->angle_data[parm_data->unique_angles_found-1].atom_type3);
          fflush(stdout);
          while( response!='y' && response!='n' )
            scanf("%c",&response);
          if (response=='y')
            parm_data->angle_data[parm_data->unique_angles_found-1].DO_FIT_THEQ=YES;
        }
      } else if (global_options->PARAMETERS_TO_FIT == DEFAULT) {
        parm_data->angle_data[parm_data->unique_angles_found-1].DO_FIT_KT=YES;
        parm_data->angle_data[parm_data->unique_angles_found-1].DO_FIT_THEQ=YES;
      } else if (global_options->PARAMETERS_TO_FIT == K_ONLY) {
        parm_data->angle_data[parm_data->unique_angles_found-1].DO_FIT_KT=NO;
        parm_data->angle_data[parm_data->unique_angles_found-1].DO_FIT_THEQ=NO;
      }
    }
  } /*End of angles withOUT Hydrogen*/
  
  /*Now we do dihedrals with hydrogen*/
  for (i=0;i<parm_data->NPHIH;++i)
  {
    /*Loop over all dihedrals involving hydrogen*/
    // Check if it's an improper dihedral
    if (parm_data->pdihedralH[i].lp < 0)
      is_improper=YES;
    else
      is_improper=NO;
    
    /*Get the atoms involved in the dihedral we are looking at*/
    atom1=unObfuscateAtom(parm_data->pdihedralH[i].ip);
    atom2=unObfuscateAtom(parm_data->pdihedralH[i].jp);
    atom3=unObfuscateAtom(parm_data->pdihedralH[i].kp);
    atom4=unObfuscateAtom(parm_data->pdihedralH[i].lp);     
    /*Get the parameters involved*/
    pk=parm_data->pk[parm_data->pdihedralH[i].icp-1];
    pn=parm_data->pn[parm_data->pdihedralH[i].icp-1];
    phase=parm_data->phase[parm_data->pdihedralH[i].icp-1];
    
    /*Get the atom types*/
    strcpy(type1,parm_data->atom[atom1-1].isymbl);
    strcpy(type2,parm_data->atom[atom2-1].isymbl);
    strcpy(type3,parm_data->atom[atom3-1].isymbl);
    strcpy(type4,parm_data->atom[atom4-1].isymbl);

    /*Ok, now we have the data for this dihedral*/
    /*Now we need to determine if it is unique*/
    /*loop over all unique dihedrals to date and see if we get a match*/
    is_unique=YES;
    for (j=0;j<parm_data->unique_dihedrals_found;++j) {
      /*A dihedral is unique if all atom types don't match*/
      /*Note, we can do the atom types in both directions*/
      /*First try the dihedral one way, then the reverse since directionality is irrelevant*/
      strcpy(dihedral_make_up1,type1);
      strcat(dihedral_make_up1,type2);
      strcat(dihedral_make_up1,type3);
      strcat(dihedral_make_up1,type4);       
      strcpy(dihedral_make_up2,parm_data->dihedral_data[j].atom_type1);
      strcat(dihedral_make_up2,parm_data->dihedral_data[j].atom_type2);
      strcat(dihedral_make_up2,parm_data->dihedral_data[j].atom_type3);
      strcat(dihedral_make_up2,parm_data->dihedral_data[j].atom_type4);       
      check1=!strcmp(dihedral_make_up1,dihedral_make_up2);             /*strcmp returns zero on a match - so check1 will be 1 if matched*/
      strcpy(dihedral_make_up1,type4);
      strcat(dihedral_make_up1,type3);
      strcat(dihedral_make_up1,type2);
      strcat(dihedral_make_up1,type1);
      strcpy(dihedral_make_up2,parm_data->dihedral_data[j].atom_type1);
      strcat(dihedral_make_up2,parm_data->dihedral_data[j].atom_type2);
      strcat(dihedral_make_up2,parm_data->dihedral_data[j].atom_type3);
      strcat(dihedral_make_up2,parm_data->dihedral_data[j].atom_type4);
      check2=!strcmp(dihedral_make_up1,dihedral_make_up2);             /*strcmp returns zero on a match - so check1 will be 1 if matched*/
      
      // Nonunique dihedral if the names match and the improper state matches an existing dihedral
      if ( (check1 || check2) && (is_improper==parm_data->dihedral_data[j].improper) ) {
        is_unique=NO;
        
        // Check if this term has already been found
        for (k=0; k<parm_data->dihedral_data[j].num_terms; ++k) {
          if (pk==parm_data->dihedral_data[j].term[k].pk && pn==parm_data->dihedral_data[j].term[k].pn &&
              phase==parm_data->dihedral_data[j].term[k].phase) break;
        }
        if (k<parm_data->dihedral_data[j].num_terms){ // broke out of loop <=> dihedral term+type nonunique
          /* The prmtop lists each dihedral found in the molecule once for each term of that dihedral. We need
           * to keep track of the atoms involved in the dihedral only for the first term of it, because if we
           * did it for every term we would over count the dihedral by num_terms. So, check if the periodicity
           * of the dihedral we found matches that of the first one. If not, just don't add the atom data to the
           * array. If this block is entered we already know about this dihedral term so we don't have to worry
           * about any of the data it contains, it's already in a dihedral_term_struct somewhere.
           */
          if (parm_data->dihedral_data[j].term[0].pn == pn) {
            ++parm_data->dihedral_data[j].number;
            if (parm_data->dihedral_data[j].number>=MAX_DIHEDRALS_PER_TYPE) { // check array isn't overflowed
              printf("*** ERROR IN PROCESS_PRMTOP - MAX_DIHEDRALS_PER_TYPE OF %d\n",MAX_DIHEDRALS_PER_TYPE);
              printf("*** EXCEEDED FOR parm_data->dihedral_data[j], j = %d\n",j);
              printf("*** EDIT MAX_DIHEDRALS_PER_TYPE IN prmtop_params.h AND RECOMPILE.\n");
              return DATA_OVERFLOW;
            }
            parm_data->dihedral_data[j].atom1[(parm_data->dihedral_data[j].number)-1]=atom1;
            parm_data->dihedral_data[j].atom2[(parm_data->dihedral_data[j].number)-1]=atom2;
            parm_data->dihedral_data[j].atom3[(parm_data->dihedral_data[j].number)-1]=atom3;
            parm_data->dihedral_data[j].atom4[(parm_data->dihedral_data[j].number)-1]=atom4;
          }
        } else { // dihedral already exists but params unique, so add a new term to this dihedral
          ++parm_data->dihedral_data[j].num_terms;
          ++parm_data->unique_dihedral_terms;
          parm_data->dihedral_data[j].term = (dihedral_data_struct *)realloc(parm_data->dihedral_data[j].term,
                                                                       sizeof(dihedral_data_struct)*parm_data->dihedral_data[j].num_terms);
          if (parm_data->dihedral_data[j].term == NULL) {
            malloc_failure_char("process_prmtop", "parm_data->dihedral_data[j].term", sizeof(dihedral_data_struct)*parm_data->dihedral_data[j].num_terms);
            return ALLOC_FAIL;
          }
          parm_data->mem_allocated+=sizeof(dihedral_data_struct);
          k = parm_data->dihedral_data[j].num_terms-1;
          // Initial values for this term
          parm_data->dihedral_data[j].term[k].pk=pk;
          parm_data->dihedral_data[j].term[k].pn=pn;
          parm_data->dihedral_data[j].term[k].phase=phase;       
          // Set terms to be fit
          set_dihedral_fit(global_options, &(parm_data->dihedral_data[j]), k, d_kp, d_np, d_phase);
        }
        break; /*Quit the loop, no point doing the rest if we have matched at this point*/
      }
    } /*essentially if we finish this loop without a hit it must be unique*/
    if (is_unique==YES) {
      /*dihedral is unique, have to extend number of dihedral_data structures by 1*/
      ++parm_data->unique_dihedrals_found;
      ++parm_data->unique_dihedral_terms;
      k=parm_data->unique_dihedrals_found-1;
      parm_data->dihedral_data = (dihedral_type_struct *)realloc(parm_data->dihedral_data,
                                                                 sizeof(dihedral_type_struct)*parm_data->unique_dihedrals_found);
      if (parm_data->dihedral_data == NULL) {
        malloc_failure_char("process_prmtop", "parm_data->dihedral_data", sizeof(dihedral_type_struct)*parm_data->unique_dihedrals_found);
        return ALLOC_FAIL;
      }
      parm_data->mem_allocated+=sizeof(dihedral_type_struct);
      // allocate first term
      parm_data->dihedral_data[k].term=malloc(sizeof(dihedral_data_struct));
      if (parm_data->dihedral_data[k].term==NULL) {
        malloc_failure_char("process_prmtop", "parm_data->dihedral_data[k].term", sizeof(dihedral_data_struct));
        return ALLOC_FAIL;
      }
      parm_data->dihedral_data[j].num_terms = 1;
      parm_data->mem_allocated+=sizeof(dihedral_data_struct);
      /*now put the data in the dihedral_structure*/
      strcpy(parm_data->dihedral_data[k].atom_type1,type1);
      strcpy(parm_data->dihedral_data[k].atom_type2,type2);
      strcpy(parm_data->dihedral_data[k].atom_type3,type3);
      strcpy(parm_data->dihedral_data[k].atom_type4,type4);       
      parm_data->dihedral_data[k].improper=is_improper;
      parm_data->dihedral_data[k].number=1;
      
      parm_data->dihedral_data[k].term[0].pk=pk;
      parm_data->dihedral_data[k].term[0].pn=pn;
      parm_data->dihedral_data[k].term[0].phase=phase;       
      parm_data->dihedral_data[k].atom1[0]=atom1;
      parm_data->dihedral_data[k].atom2[0]=atom2;
      parm_data->dihedral_data[k].atom3[0]=atom3;
      parm_data->dihedral_data[k].atom4[0]=atom4;   
      // Set terms to be fit
      set_dihedral_fit(global_options, &(parm_data->dihedral_data[k]), 0, d_kp, d_np, d_phase);
    }
  } /*End of dihedrals with Hydrogen*/
  
  /*Now we do dihedrals withOUT hydrogen*/
  for (i=0;i<parm_data->NPHIA;++i) {
    /*Loop over all dihedrals not involving hydrogen*/
    /*Check if it's an improper dihedral*/
    if (parm_data->pdihedral[i].lp < 0)
      is_improper=YES;
    else
      is_improper=NO;
    
    /*Get the atoms involved in the dihedral we are looking at*/
    atom1=unObfuscateAtom(parm_data->pdihedral[i].ip);
    atom2=unObfuscateAtom(parm_data->pdihedral[i].jp);
    atom3=unObfuscateAtom(parm_data->pdihedral[i].kp);
    atom4=unObfuscateAtom(parm_data->pdihedral[i].lp);
    /*Get the parameters involved*/
    pk=parm_data->pk[parm_data->pdihedral[i].icp-1];
    pn=parm_data->pn[parm_data->pdihedral[i].icp-1];
    phase=parm_data->phase[parm_data->pdihedral[i].icp-1];
    
    /*Get the atom types*/
    strcpy(type1,parm_data->atom[atom1-1].isymbl);
    strcpy(type2,parm_data->atom[atom2-1].isymbl);
    strcpy(type3,parm_data->atom[atom3-1].isymbl);
    strcpy(type4,parm_data->atom[atom4-1].isymbl);

    /*Ok, now we have the data for this dihedral*/
    /*Now we need to determine if it is unique*/
    /*loop over all unique dihedrals to date and see if we get a match*/
    is_unique=YES;
    for (j=0;j<parm_data->unique_dihedrals_found;++j) {
      /*A dihedral is unique if all atom types don't match*/
      /*Note, we can do the atom types in both directions*/
      /*First try the dihedral one way, then the reverse since directionality is irrelevant*/
      strcpy(dihedral_make_up1,type1);
      strcat(dihedral_make_up1,type2);
      strcat(dihedral_make_up1,type3);
      strcat(dihedral_make_up1,type4);
      strcpy(dihedral_make_up2,parm_data->dihedral_data[j].atom_type1);
      strcat(dihedral_make_up2,parm_data->dihedral_data[j].atom_type2);
      strcat(dihedral_make_up2,parm_data->dihedral_data[j].atom_type3);
      strcat(dihedral_make_up2,parm_data->dihedral_data[j].atom_type4);
      check1=!strcmp(dihedral_make_up1,dihedral_make_up2);             /*strcmp returns zero on a match - so check1 will be 1 if matched*/
      strcpy(dihedral_make_up1,type4);
      strcat(dihedral_make_up1,type3);
      strcat(dihedral_make_up1,type2);
      strcat(dihedral_make_up1,type1);
      strcpy(dihedral_make_up2,parm_data->dihedral_data[j].atom_type1);
      strcat(dihedral_make_up2,parm_data->dihedral_data[j].atom_type2);
      strcat(dihedral_make_up2,parm_data->dihedral_data[j].atom_type3);
      strcat(dihedral_make_up2,parm_data->dihedral_data[j].atom_type4);
      check2=!strcmp(dihedral_make_up1,dihedral_make_up2);             /*strcmp returns zero on a match - so check1 will be 1 if matched*/
      
      // We have a nonunique dihedral if the name matches, and the improper state matches
      if ( (check1 || check2) && (is_improper==parm_data->dihedral_data[j].improper) ) {
        is_unique=NO;
        // Check if this term has already been found
        for (k=0; k<parm_data->dihedral_data[j].num_terms; ++k) {
          if (pk==parm_data->dihedral_data[j].term[k].pk && pn==parm_data->dihedral_data[j].term[k].pn &&
              phase==parm_data->dihedral_data[j].term[k].phase) break;
        }
        if (k<parm_data->dihedral_data[j].num_terms){ // broke out of loop <=> dihedral term+type nonunique
          
          // Please see the long comment in dihedrals with H above that explains the logic here
          if (parm_data->dihedral_data[j].term[0].pn == pn) {
            ++parm_data->dihedral_data[j].number;
            if (parm_data->dihedral_data[j].number>=MAX_DIHEDRALS_PER_TYPE) {
              printf("*** ERROR IN PROCESS_PRMTOP - MAX_DIHEDRALS_PER_TYPE OF %d\n",MAX_DIHEDRALS_PER_TYPE);
              printf("*** EXCEEDED FOR parm_data->dihedral_data[j].term[k], j = %d k=%d\n",j,k);
              printf("*** EDIT MAX_DIHEDRALS_PER_TYPE IN prmtop_params.h AND RECOMPILE.\n");
              return DATA_OVERFLOW;
            }
            parm_data->dihedral_data[j].atom1[(parm_data->dihedral_data[j].number)-1]=atom1;
            parm_data->dihedral_data[j].atom2[(parm_data->dihedral_data[j].number)-1]=atom2;
            parm_data->dihedral_data[j].atom3[(parm_data->dihedral_data[j].number)-1]=atom3;
            parm_data->dihedral_data[j].atom4[(parm_data->dihedral_data[j].number)-1]=atom4;
          }
        }
        else { // dihedral already exists but params unique, so add a new term to this dihedral
          ++parm_data->unique_dihedral_terms;
          ++parm_data->dihedral_data[j].num_terms;
          parm_data->dihedral_data[j].term = (dihedral_data_struct *)realloc(parm_data->dihedral_data[j].term,
                                                                       sizeof(dihedral_data_struct)*parm_data->dihedral_data[j].num_terms);
          if (parm_data->dihedral_data[j].term == NULL) {
            malloc_failure_char("process_prmtop", "parm_data->dihedral_data[j].term", sizeof(dihedral_data_struct)*parm_data->dihedral_data[j].num_terms);
            return ALLOC_FAIL;
          }
          parm_data->mem_allocated+=sizeof(dihedral_data_struct);
          k = parm_data->dihedral_data[j].num_terms-1;
          // Initial values for this term
          parm_data->dihedral_data[j].term[k].pk=pk;
          parm_data->dihedral_data[j].term[k].pn=pn;
          parm_data->dihedral_data[j].term[k].phase=phase;       
          // Set terms to be fit
          set_dihedral_fit(global_options, &(parm_data->dihedral_data[j]), k, d_kp, d_np, d_phase);
        }
        break; /*Quit the loop, no point doing the rest if we have matched at this point*/
      }
    } /*essentially if we finish this loop without a hit it must be unique*/
    if (is_unique==YES) {
      /*dihedral is unique, have to extend number of dihedral_data structures by 1*/
      ++(parm_data->unique_dihedrals_found);
      ++(parm_data->unique_dihedral_terms);
      k=parm_data->unique_dihedrals_found-1;
      parm_data->dihedral_data = (dihedral_type_struct*)realloc(parm_data->dihedral_data,sizeof(dihedral_type_struct)*parm_data->unique_dihedrals_found);
      if (parm_data->dihedral_data == NULL) {
        malloc_failure_char("process_prmtop", "parm_data->dihedral_data", sizeof(dihedral_type_struct)*parm_data->unique_dihedrals_found);
        return ALLOC_FAIL;
      }
      parm_data->mem_allocated+=sizeof(dihedral_type_struct);
      // allocate initial term
      parm_data->dihedral_data[j].term = (dihedral_data_struct*)malloc(sizeof(dihedral_data_struct));
      if (parm_data->dihedral_data[j].term==NULL) {
        malloc_failure_char("process_prmtop","parm_data->dihedral_data[k]",sizeof(dihedral_data_struct));
        return ALLOC_FAIL;
      }
      parm_data->dihedral_data[j].num_terms = 1;
      parm_data->mem_allocated+=sizeof(dihedral_data_struct);
      // Put data in the dihedral structure
      strcpy(parm_data->dihedral_data[k].atom_type1,type1);
      strcpy(parm_data->dihedral_data[k].atom_type2,type2);
      strcpy(parm_data->dihedral_data[k].atom_type3,type3);
      strcpy(parm_data->dihedral_data[k].atom_type4,type4);
      parm_data->dihedral_data[k].improper=is_improper;
      parm_data->dihedral_data[k].number=1;
      
      parm_data->dihedral_data[k].term[0].pk=pk;
      parm_data->dihedral_data[k].term[0].pn=pn;
      parm_data->dihedral_data[k].term[0].phase=phase;
      parm_data->dihedral_data[k].atom1[0]=atom1;
      parm_data->dihedral_data[k].atom2[0]=atom2;
      parm_data->dihedral_data[k].atom3[0]=atom3;
      parm_data->dihedral_data[k].atom4[0]=atom4;
      // Set fitting options for this dihedral
      set_dihedral_fit(global_options, &(parm_data->dihedral_data[k]), 0, d_kp, d_np, d_phase);
    }
  } /*End of dihedrals withOUT Hydrogen*/
 
  // Save the inputted options to a file if desired
  if (global_options->RUNTYPE==SET_PARAMS) { 
    if (write_input_parameters(global_options, parm_data) != SUCCESS)
      return FAILURE;
  }
  
  // Load options if indicated
  if (global_options->PARAMETERS_TO_FIT == LOAD) {
   if (read_parameter_file(global_options, parm_data) != SUCCESS)
     return FAILURE;
  }
  
  if (global_options->VERBOSITY>=HIGH) {
    printf("Contents of parm_data->bond_data array:\n");
    for(i=0;i<parm_data->unique_bonds_found;++i) {
      printf("%d: number=%d, %s-%s, rk=%.4f, req=%.4f\n",i,parm_data->bond_data[i].number,parm_data->bond_data[i].atom_type1
      ,parm_data->bond_data[i].atom_type2,parm_data->bond_data[i].rk,parm_data->bond_data[i].req);
      printf("   ATOMS: ");
      for (j=0;j<parm_data->bond_data[i].number;++j)
        printf("%d-%d ",parm_data->bond_data[i].atom1[j],parm_data->bond_data[i].atom2[j]);
      printf("\n");
    }
    
    printf("Contents of parm_data->angle_data array:\n");
    for(i=0;i<parm_data->unique_angles_found;++i) {
      printf("%d: number=%d, %s-%s-%s, tk=%.4f, teq=%.4f\n",i,parm_data->angle_data[i].number,parm_data->angle_data[i].atom_type1
      ,parm_data->angle_data[i].atom_type2,parm_data->angle_data[i].atom_type3,parm_data->angle_data[i].tk,parm_data->angle_data[i].teq*RADIAN_TO_DEGREE);
      printf("   ATOMS: ");
      for (j=0;j<parm_data->angle_data[i].number;++j)
        printf("%d-%d-%d ",parm_data->angle_data[i].atom1[j],parm_data->angle_data[i].atom2[j],parm_data->angle_data[i].atom3[j]);
      printf("\n");
    }
    
    printf("Contents of parm_data->dihedral_data array:\n");
    for(i=0;i<parm_data->unique_dihedrals_found;++i) {
      for (j=0; j<parm_data->dihedral_data[i].num_terms; ++j) {
        printf("%d: number=%d, %s-%s-%s-%s, pk=%.4f, phase=%.4f, pn=%.4f\n",i,parm_data->dihedral_data[i].number,parm_data->dihedral_data[i].atom_type1
        ,parm_data->dihedral_data[i].atom_type2,parm_data->dihedral_data[i].atom_type3,parm_data->dihedral_data[i].atom_type4
        ,parm_data->dihedral_data[i].term[j].pk,parm_data->dihedral_data[i].term[j].phase*RADIAN_TO_DEGREE,parm_data->dihedral_data[i].term[j].pn);
        printf("   ATOMS: ");
        for (k=0;k<parm_data->dihedral_data[i].number;++k) {
          printf("%d-%d-%d-%d ",parm_data->dihedral_data[i].atom1[k],parm_data->dihedral_data[i].atom2[k],
                 parm_data->dihedral_data[i].atom3[k],parm_data->dihedral_data[i].atom4[k]);
        }
        printf("\n");
      }
    }
  }
  
  if (global_options->VERBOSITY>=MEDIUM) {
    printf("   Prmtop   (unique): Found %d unique bonds.\n",parm_data->unique_bonds_found);
    printf("   Prmtop   (unique): Found %d unique angles.\n",parm_data->unique_angles_found);
    printf("   Prmtop   (unique): Found %d unique dihedrals.\n",parm_data->unique_dihedrals_found);
    printf("   Prmtop   (unique): Found %d unique dihedral terms.\n",parm_data->unique_dihedral_terms);
  }
    
  return SUCCESS;
}

/**
 * Sets dihedral fitting options for a given dihedral.
 * All parameters for this term should have already been initialized. 
 * @param[in] global_options The global options structure with prompting options
 * @param[in,out] s The dihedral type to set options for
 * @param[in] t Integer indicating the term to fit
 * @param[in] d_kp Whether to prompt to fit dihedral KP
 * @param[in] d_np Whether to prompt to fit dihedral NP
 * @param[in] d_phase Whether to prompt to fit dihedral PHASE
 */
void set_dihedral_fit(global_options_struct *global_options, dihedral_type_struct *s, int t,
                     bool_t d_kp, bool_t d_np, bool_t d_phase)
{
  /*Need to check if we are to fit this parameter or not*/
  char response = '\0';
  s->term[t].DO_FIT_NP=NO;
  s->term[t].DO_FIT_KP=NO;
  s->term[t].DO_FIT_PHASE=NO;
  
  if (global_options->RUNTYPE==SET_PARAMS)
  {
    /*prompt the user whether to fit or not*/
    if (d_kp == YES)
    {
      response='\0';
      printf("Fit Parameter: ");
      if (s->improper==YES)
        printf("IMPROPER ");
      printf("DIHEDRAL (%s-%s-%s-%s) term %d KP? (y/n): ",
             s->atom_type1,
             s->atom_type2,
             s->atom_type3,
             s->atom_type4,
             t+1);
      fflush(stdout);
      while( response!='y' && response!='n' )
        scanf("%c",&response);
      if (response=='y')
        s->term[t].DO_FIT_KP=YES;
    }
    
    if (d_np == YES)
    {
      response='\0';
        printf("Fit Parameter: ");
        if (s->improper==YES)
          printf("IMPROPER ");
        printf("DIHEDRAL (%s-%s-%s-%s) term %d NP? (y/n): ",
             s->atom_type1,
             s->atom_type2,
             s->atom_type3,
             s->atom_type4,
             t+1);
        fflush(stdout);
        while( response!='y' && response!='n' )
          scanf("%c",&response);
        if (response=='y')
          s->term[t].DO_FIT_NP=YES;
    }
    
    if (d_phase==YES)
    {
      response='\0';
      printf("Fit Parameter: ");
      if (s->improper==YES) {
        printf("IMPROPER -- skipping phase\n");
        s->term[t].DO_FIT_PHASE=NO;
      } else {
        printf("DIHEDRAL (%s-%s-%s-%s) term %d PHASE? (y/n): ",
               s->atom_type1,
               s->atom_type2,
               s->atom_type3,
               s->atom_type4,
               t+1);
        fflush(stdout);
        while( response!='y' && response!='n' )
          scanf("%c",&response);
        if (response=='y')
          s->term[t].DO_FIT_PHASE=YES;
      }
    }
  }
  else if (global_options->PARAMETERS_TO_FIT == DEFAULT)
  {
    /*DEFAULT FOR DIHEDRALS IS TO FIT ONLY KP, this will be adjusted as necessary by other routines*/
    s->term[t].DO_FIT_KP=YES;
    s->term[t].DO_FIT_NP=NO;
    s->term[t].DO_FIT_PHASE=NO;
  } 
  else if (global_options->PARAMETERS_TO_FIT == K_ONLY) {
    s->term[t].DO_FIT_KP=NO;
    s->term[t].DO_FIT_NP=NO;
    s->term[t].DO_FIT_PHASE=NO;
  }
}

