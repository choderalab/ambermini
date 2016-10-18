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

/** @file eval_amber_std.c
 * Function evaluation code for energies with the AMBER equation.
 * Evaluates the sum of the squares values using the standard AMBER force field,
 * using the values in the parameter file for each structure in the coordinate data.
 */

#include <stdio.h>
#include <math.h>
#include "function_def.h"
#include "constants.h"

/**
 * Conducts the energy evaluation over every prmtop.
 * @param[in] global_options The global options structure
 * @param[in] parm_datas The array of parm structures
 * @param[in] coords_datas The array of coordinate sets
 * @param[out] range_altered Pointer to int that will be set to the number of parameters altered due to their being out of bounds.
 *                           Functions should check this flag to see if they need to copy parm_data back into their internal data structures.
 * @return Sum of squares of energy difference over all structures.
 */
double eval_sum_squares_amber_std_multiprmtop(global_options_struct *global_options, parm_struct *parm_datas, coord_set *coords_datas, int *range_altered)
{
  double score = 0.0;
  int i;
  if (range_altered != NULL) range_altered = 0;
  for (i=0; i<global_options->num_prmtops; ++i) {
    score += eval_sum_squares_amber_std(global_options, &(parm_datas[i]), &(coords_datas[i]), range_altered);
  }
  ++global_options->function_calls;
  
  return score;
}

/**
 * Evaluates the sum squares difference between quantum and calculated energy for each structure.
 * The energy calculation for each structure is done in parallel with openmp.
 * 
 * @param[in] global_options The global options structure
 * @param[in] parm_data The parameter file containing parameter set to be scored
 * @param[in] coords_data Contains coordinates of atoms in all input structures, and qm energies (pointer to ONE coordinate set)
 * @param[out] range_altered Pointer to int that will be set to the number of parameters altered due to their being out of bounds.
 *                           Functions should check this flag to see if they need to copy parm_data back into their internal data structures.
 *                           If set to NULL, no bounds checking will be performed
 * @return The sum of the squares of the difference in energy between AMBER function evaluation and quantum value.
 */
double eval_sum_squares_amber_std(global_options_struct *global_options, parm_struct *parm_data, coord_set *coords_data, int *range_altered)
{
  // Before anything else is done, check that the parameters to be evaluated
  // are valid, and if not, correct them.
  if (range_altered != NULL)
    *range_altered += check_range(global_options, parm_data);
  
  // Parallel portion of code - each thread evaluates a subset of the structures

  /* The force field equation is:

     Ebond     = Kr(r-req)^2     - Kr is in KCal/(mol A)^2 and req is in angstroms
     Eangle    = Ktheta(theta - thetaeq)^2   - Ktheta is in KCal/mol rad)^2, thetaeq is in radians
     Edihedral = 0.5Vn(1+cos(n*phi-gamma))  Vn is in KCal/mol, phi and gamma are in radians
     Eelec     = Sum(i<j)(qiqj/eRij)
     Evdw      = Sum(i<j)([Aij/Rij^12]-[Bij/Rij^6])

     Etot = sum over all bonds, angles, dihedrals, vdw and Eelec

     In terms of minimisation what we shall do is evaluate the sum of the squares for a range of different
     structures, the more the better.

     Ideally for any structure

     EQM = EAMBER (but the origins are not the same so this is strictly dEQM = dEAMBER)

     Therefore (Ebonds + Eangles + Edihedrals + Eelec + Evdw) - EQM = 0

     Here we aren't fitting van der waals or electrostatics but we will be calculating them. Note there are
     no periodics here and no cut offs.

     Since QM energy minimum does not necessarily match Amber energy minimum so we have a value of K. Thus our
     equation is:

     f = Sum[1->N]((Sum[1->nb]Bonds + Sum[1-na]Angles + Sum[1->nd]Dihedrals + Sum[i<j]Elec + Sum[i<j]VDW + K - E[QMn])^2)

  */

  double result;
  int structure;
  
  result=0.0;
#pragma omp parallel for reduction(+:result)
  for (structure=0; structure<coords_data->num_coords; ++structure)
  {
  double individual_sum;
  int i,j,k;
  double tempx1, tempy1, tempz1, tempx2, tempy2, tempz2;
  double tempx3, tempy3, tempz3, tempx4, tempy4, tempz4;
  double bond_length_temp;
  double angle_temp;
  double dihedral_temp;
  double elec_temp;
  double vdw_temp;
  short int excluded;
  int excluded_temp;
  int excluded_offset;
  int vdw_offset;
  double qi, qj, Aij, Bij, Rij;
  
    /*Add K, subtract the QM energy, add all the bonds, angles, dihedrals nonbonds etc and then square*/
    individual_sum = parm_data->K;
    /*This should subtract energy[structure] from the sum*/
    
    individual_sum-=coords_data->struc[structure].energy;

  /*Now do the electrostatics and van der waals - do both together for speed*/
  /*note we avoid double counting in the electrostatics since we do sum(i<j)[qiqj]
                                                                             ----
                                                                             Rij*/
  excluded_offset=0;
  for (i=0;i<parm_data->NTOTAT;++i) /*loop over all atoms*/
  {
     excluded_temp=0; /*number that have been exlcuded for this atom to date*/
     for (j=i+1;j<parm_data->NTOTAT;++j)
     {
       /*This loops as follows: atom 1 ( 2 to NTOTAT), atom 2 ( 3 to NTOTAT), atom3 ( 4 to NTOTAT) etc.*/
       /*Only calculate if this pair is not excluded in the prmtop file - since it is a 1-2, 1-3 or 1-4.*/
       excluded=NO;
       if (excluded_temp<parm_data->atom[i].numex)
       {
         /*check if this pair is excluded*/
         if ( j == (parm_data->natex[excluded_offset+excluded_temp]-1) )
         {
           ++excluded_temp;
           excluded=YES;
         }
       }
       if (excluded==NO)
       {
           /*not excluded, calculate elec and vdw for this atom pair*/
           qi = parm_data->atom[i].chrg;
           qj = parm_data->atom[j].chrg;
	   
           /*Note, the way Aij and Bij are looked up is quite complex. Essentially we have a diagonal matrix of
             Ntypes^2 in *nno

             For example for NMA there are 7 types, so nno is 28 ints long

             1  2  4  7 11 16 22
             2  3  5  8 12 17 23
             4  5  6  9 13 18 24
             7  8  9 10 14 19 25
            11 12 13 14 15 10 26
            16 17 18 19 20 21 27
            22 23 24 25 26 27 28

             We then look this up based on the atom types of the two atoms in the pair. So, we use the first atom type in the
             pair to get the offset into this array in the form ntypes*atom_type_number. Then we move the number of types represented
             by the integer type of atom 2 further into the array. This gives us the index for the Aij and Bij coefficient arrays.

              E.g if we were looking at atoms 5 and 10 we first get the atom type of 5 which is 3, so our offset is 7 * 2 = 14 which positions
              us at the first element of row 3. Then we check atom 10's type which is 7 so we move to column 7, offset = 7*(type1-1)+(type2 - 1)
              = 20 - Not 21 since we are counting from 0. This returns the number 24 which is used to find the 24th element of the A coefficient
              and B coefficient arrays*/
	   
           vdw_offset=(parm_data->NTYPES*(parm_data->atom[i].iac - 1))+(parm_data->atom[j].iac - 1);
           Aij= parm_data->cn1[(parm_data->nno[vdw_offset]-1)];
           Bij= parm_data->cn2[(parm_data->nno[vdw_offset]-1)];
	   
           tempx1=coords_data->struc[structure].x_coord[i];
           tempy1=coords_data->struc[structure].y_coord[i];
           tempz1=coords_data->struc[structure].z_coord[i];
           tempx2=coords_data->struc[structure].x_coord[j];
           tempy2=coords_data->struc[structure].y_coord[j];
           tempz2=coords_data->struc[structure].z_coord[j];
           /*find distance between*/
           Rij = calc_bond_length(tempx1,tempy1,tempz1,tempx2,tempy2,tempz2);
	   
           /*now find the electronic energy for this pair*/
           individual_sum+=(qi*qj)/(Rij);
           /*and now the vdw energy*/
           individual_sum+=(Aij / pow(Rij,12)) - (Bij / pow(Rij,6));
       }
     }
     excluded_offset+=parm_data->atom[i].numex;
    }
    
    /*Note, above does not do 1-4 interactions, we do this by considering each dihedral term for which ALL the atom types are +ve*/
    /*dihedrals with H first*/
    for (i=0;i<parm_data->NPHIH;++i)
    {
       /*Check all atoms are +ve*/
       if (parm_data->pdihedralH[i].ip>=0 && parm_data->pdihedralH[i].jp >=0 && parm_data->pdihedralH[i].kp >=0 && parm_data->pdihedralH[i].lp >=0)
       {
           /*all are +ve so calculate 1-4 EE and 1-4 VDW, don't forget to scale*/
           qi = parm_data->atom[unObfuscateAtom(parm_data->pdihedralH[i].ip)-1].chrg;
           qj = parm_data->atom[unObfuscateAtom(parm_data->pdihedralH[i].lp)-1].chrg;
           /*Note, the way Aij and Bij are looked up is quite complex. Essentially we have a diagonal matrix of
             Ntypes^2 in *nno

             For example for NMA there are 7 types, so nno is 28 ints long

             1  2  4  7 11 16 22
             2  3  5  8 12 17 23
             4  5  6  9 13 18 24
             7  8  9 10 14 19 25
            11 12 13 14 15 10 26
            16 17 18 19 20 21 27
            22 23 24 25 26 27 28

             We then look this up based on the atom types of the two atoms in the pair. So, we use the first atom type in the
             pair to get the offset into this array in the form ntypes*atom_type_number. Then we move the number of types represented
             by the integer type of atom 2 further into the array. This gives us the index for the Aij and Bij coefficient arrays.

              E.g if we were looking at atoms 5 and 10 we first get the atom type of 5 which is 3, so our offset is 7 * 2 = 14 which positions
              us at the first element of row 3. Then we check atom 10's type which is 7 so we move to column 7, offset = 7*(type1-1)+(type2 - 1)
              = 20 - Not 21 since we are counting from 0. This returns the number 24 which is used to find the 24th element of the A coefficient
              and B coefficient arrays*/
           vdw_offset=(parm_data->NTYPES*(parm_data->atom[unObfuscateAtom(parm_data->pdihedralH[i].ip)-1].iac - 1))
                      +(parm_data->atom[unObfuscateAtom(parm_data->pdihedralH[i].lp)-1].iac - 1);
           Aij= parm_data->cn1[(parm_data->nno[vdw_offset]-1)];
           Bij= parm_data->cn2[(parm_data->nno[vdw_offset]-1)];
           tempx1=coords_data->struc[structure].x_coord[unObfuscateAtom(parm_data->pdihedralH[i].ip)-1];
           tempy1=coords_data->struc[structure].y_coord[unObfuscateAtom(parm_data->pdihedralH[i].ip)-1];
           tempz1=coords_data->struc[structure].z_coord[unObfuscateAtom(parm_data->pdihedralH[i].ip)-1];
           tempx2=coords_data->struc[structure].x_coord[unObfuscateAtom(parm_data->pdihedralH[i].lp)-1];
           tempy2=coords_data->struc[structure].y_coord[unObfuscateAtom(parm_data->pdihedralH[i].lp)-1];
           tempz2=coords_data->struc[structure].z_coord[unObfuscateAtom(parm_data->pdihedralH[i].lp)-1];
           /*find distance between*/
           Rij = calc_bond_length(tempx1,tempy1,tempz1,tempx2,tempy2,tempz2);

           /*now find the electronic energy for this pair*/
           elec_temp=(qi*qj)/(Rij);
           elec_temp/=global_options->SCEE;
           individual_sum+=elec_temp;
           /*and now the vdw energy*/
           vdw_temp=(Aij / pow(Rij,12)) - (Bij / pow(Rij,6));
           vdw_temp/=global_options->SCNB;
           individual_sum+=vdw_temp;
       }
    }
    /*dihedrals without H*/
    for (i=0;i<parm_data->NPHIA;++i)
    {
       /*Check all atoms are +ve*/
       if (parm_data->pdihedral[i].ip>=0 && parm_data->pdihedral[i].jp >=0 && parm_data->pdihedral[i].kp >=0 && parm_data->pdihedral[i].lp >=0)
       {
           /*all are +ve so calculate 1-4 EE and 1-4 VDW, don't forget to scale*/
           /*all are +ve so calculate 1-4 EE and 1-4 VDW, don't forget to scale*/
           qi = parm_data->atom[unObfuscateAtom(parm_data->pdihedral[i].ip)-1].chrg;
           qj = parm_data->atom[unObfuscateAtom(parm_data->pdihedral[i].lp)-1].chrg;
           /*Note, the way Aij and Bij are looked up is quite complex. Essentially we have a diagonal matrix of
             Ntypes^2 in *nno

             For example for NMA there are 7 types, so nno is 28 ints long

             1  2  4  7 11 16 22
             2  3  5  8 12 17 23
             4  5  6  9 13 18 24
             7  8  9 10 14 19 25
            11 12 13 14 15 10 26
            16 17 18 19 20 21 27
            22 23 24 25 26 27 28

             We then look this up based on the atom types of the two atoms in the pair. So, we use the first atom type in the
             pair to get the offset into this array in the form ntypes*atom_type_number. Then we move the number of types represented
             by the integer type of atom 2 further into the array. This gives us the index for the Aij and Bij coefficient arrays.

              E.g if we were looking at atoms 5 and 10 we first get the atom type of 5 which is 3, so our offset is 7 * 2 = 14 which positions
              us at the first element of row 3. Then we check atom 10's type which is 7 so we move to column 7, offset = 7*(type1-1)+(type2 - 1)
              = 20 - Not 21 since we are counting from 0. This returns the number 24 which is used to find the 24th element of the A coefficient
              and B coefficient arrays*/
           vdw_offset=(parm_data->NTYPES*(parm_data->atom[unObfuscateAtom(parm_data->pdihedral[i].ip)-1].iac - 1))
                      +(parm_data->atom[unObfuscateAtom(parm_data->pdihedral[i].lp)-1].iac - 1);
           Aij= parm_data->cn1[(parm_data->nno[vdw_offset]-1)];
           Bij= parm_data->cn2[(parm_data->nno[vdw_offset]-1)];
           tempx1=coords_data->struc[structure].x_coord[unObfuscateAtom(parm_data->pdihedral[i].ip)-1];
           tempy1=coords_data->struc[structure].y_coord[unObfuscateAtom(parm_data->pdihedral[i].ip)-1];
           tempz1=coords_data->struc[structure].z_coord[unObfuscateAtom(parm_data->pdihedral[i].ip)-1];
           tempx2=coords_data->struc[structure].x_coord[unObfuscateAtom(parm_data->pdihedral[i].lp)-1];
           tempy2=coords_data->struc[structure].y_coord[unObfuscateAtom(parm_data->pdihedral[i].lp)-1];
           tempz2=coords_data->struc[structure].z_coord[unObfuscateAtom(parm_data->pdihedral[i].lp)-1];
           /*find distance between*/
           Rij = calc_bond_length(tempx1,tempy1,tempz1,tempx2,tempy2,tempz2);

           /*now find the electronic energy for this pair*/
           elec_temp=(qi*qj)/(Rij);
           elec_temp/=global_options->SCEE;
           individual_sum+=elec_temp;
           /*and now the vdw energy*/
           vdw_temp=(Aij / pow(Rij,12)) - (Bij / pow(Rij,6));
           vdw_temp/=global_options->SCNB;
           individual_sum+=vdw_temp;
       }
    }
    
    /*Now do the bonds*/
    for (i=0;i<parm_data->unique_bonds_found;++i)
    {
       for (j=0;j<parm_data->bond_data[i].number;++j)
       {
         /*First off we need to find the distance between the atoms in the bond*/
         /*Get the coordinates*/
         tempx1=coords_data->struc[structure].x_coord[parm_data->bond_data[i].atom1[j]-1];
         tempy1=coords_data->struc[structure].y_coord[parm_data->bond_data[i].atom1[j]-1];
         tempz1=coords_data->struc[structure].z_coord[parm_data->bond_data[i].atom1[j]-1];
         tempx2=coords_data->struc[structure].x_coord[parm_data->bond_data[i].atom2[j]-1];
         tempy2=coords_data->struc[structure].y_coord[parm_data->bond_data[i].atom2[j]-1];
         tempz2=coords_data->struc[structure].z_coord[parm_data->bond_data[i].atom2[j]-1];
         /*find bond length*/
         bond_length_temp=calc_bond_length(tempx1,tempy1,tempz1,tempx2,tempy2,tempz2);
	 
         /*Now calculate the energy contribution for this bond and add it to the individual sum*/
         individual_sum+=(parm_data->bond_data[i].rk)*pow(bond_length_temp-parm_data->bond_data[i].req,2);
	                           /*BONDjFC*/                    /*R*/                         /*Req*/
       }
    }
    
    
    /*Angles*/
    for (i=0;i<parm_data->unique_angles_found;++i)
    {
       for (j=0;j<parm_data->angle_data[i].number;++j)
       {
         /*First off we need to find the angle between the atoms in the angle*/
         /*Get the coordinates*/
         tempx1=coords_data->struc[structure].x_coord[parm_data->angle_data[i].atom1[j]-1];
         tempy1=coords_data->struc[structure].y_coord[parm_data->angle_data[i].atom1[j]-1];
         tempz1=coords_data->struc[structure].z_coord[parm_data->angle_data[i].atom1[j]-1];
         tempx2=coords_data->struc[structure].x_coord[parm_data->angle_data[i].atom2[j]-1];
         tempy2=coords_data->struc[structure].y_coord[parm_data->angle_data[i].atom2[j]-1];
         tempz2=coords_data->struc[structure].z_coord[parm_data->angle_data[i].atom2[j]-1];
         tempx3=coords_data->struc[structure].x_coord[parm_data->angle_data[i].atom3[j]-1];
         tempy3=coords_data->struc[structure].y_coord[parm_data->angle_data[i].atom3[j]-1];
         tempz3=coords_data->struc[structure].z_coord[parm_data->angle_data[i].atom3[j]-1];

         /*find angle*/
         angle_temp=calc_angle_radians(tempx1,tempy1,tempz1,tempx2,tempy2,tempz2,tempx3,tempy3,tempz3);
	 
         /*Now calculate the energy contribution for this angle and add it to the individual sum*/
         individual_sum+=(parm_data->angle_data[i].tk)*pow(angle_temp-parm_data->angle_data[i].teq,2);
                              /*ANGLEjFC*/                    /*T*/                         /*Teq*/
       }
    }


    /*Dihedrals*/
    for (i=0;i<parm_data->unique_dihedrals_found;++i) // each dihedral type
    {
       for (j=0;j<parm_data->dihedral_data[i].num_terms;++j) // each term of that dihedral
       {
          for (k=0; k<parm_data->dihedral_data[i].number; ++k) { // each instance of that term
          /*First off we need to find the dihedral between the atoms in the angle*/
          /*Get the coordinates*/
          tempx1=coords_data->struc[structure].x_coord[parm_data->dihedral_data[i].atom1[k]-1];
          tempy1=coords_data->struc[structure].y_coord[parm_data->dihedral_data[i].atom1[k]-1];
          tempz1=coords_data->struc[structure].z_coord[parm_data->dihedral_data[i].atom1[k]-1];
          tempx2=coords_data->struc[structure].x_coord[parm_data->dihedral_data[i].atom2[k]-1];
          tempy2=coords_data->struc[structure].y_coord[parm_data->dihedral_data[i].atom2[k]-1];
          tempz2=coords_data->struc[structure].z_coord[parm_data->dihedral_data[i].atom2[k]-1];
          tempx3=coords_data->struc[structure].x_coord[parm_data->dihedral_data[i].atom3[k]-1];
          tempy3=coords_data->struc[structure].y_coord[parm_data->dihedral_data[i].atom3[k]-1];
          tempz3=coords_data->struc[structure].z_coord[parm_data->dihedral_data[i].atom3[k]-1];
          tempx4=coords_data->struc[structure].x_coord[parm_data->dihedral_data[i].atom4[k]-1];
          tempy4=coords_data->struc[structure].y_coord[parm_data->dihedral_data[i].atom4[k]-1];
          tempz4=coords_data->struc[structure].z_coord[parm_data->dihedral_data[i].atom4[k]-1];

          /*find dihedral*/
          dihedral_temp=calc_dihedral_radians(tempx1,tempy1,tempz1,tempx2,tempy2,tempz2,tempx3,tempy3,tempz3,tempx4,tempy4,tempz4);
          dihedral_temp=PI-dihedral_temp; /*Due to amber calculating the angle in the opposite direction to me*/
          
          /*Note, pk is actually pk/2*/
          /*Now calculate the energy contribution for this dihedral and add it to the individual sum*/
          individual_sum+=(parm_data->dihedral_data[i].term[j].pk)*
          (1+cos(dihedral_temp*parm_data->dihedral_data[i].term[j].pn-parm_data->dihedral_data[i].term[j].phase));
         }
       }
    }

    /*All done, now square the result and add it to our running total*/
//     printf("Energy for %d = %f - %f\n", structure, individual_sum+coords_data[structure].energy-global_options->K, coords_data[structure].energy);
  
    result+=(individual_sum*individual_sum*coords_data->struc[structure].weight);
    // Save first energy evaluation for each structure
    if (coords_data->struc[structure].init_energy==0.0) coords_data->struc[structure].init_energy=individual_sum;
  }  
  /*All threads are returning COMPLETE result here as openMP reduction occurs at end of for block*/

  return result ;
}

/**
 * Calculates the AMBER energy for a single structure.
 * @param[in] global_options The global options structure
 * @param[in] parm_data Pointer to a single set of parameters to use in the calculation
 * @param[in] coords_data Pointer to a single coordinate structure to evaluate
 * @return The energy of the structure
 */
double eval_amber_std_for_single_struct(global_options_struct *global_options, parm_struct *parm_data, coords_struct *coords_data)
{
  /*This function evaluates the standard amber energy for a single structure*/
  /*Note structures are numbered from 1 to NSTRUCTURES*/

  /*Note there is no adjustment for K in this routine*/
  double individual_sum;
  int i,j,k;
  double tempx1, tempy1, tempz1, tempx2, tempy2, tempz2;
  double tempx3, tempy3, tempz3, tempx4, tempy4, tempz4;
  double bond_length_temp;
  double angle_temp;
  double dihedral_temp;
  double vdw_temp;
  double elec_temp;
  short int excluded;
  int excluded_temp;
  int excluded_offset;
  double bond_energy;
  double angle_energy;
  double dihedral_energy;
  double elec_energy;
  double vdw_energy;
  double vdw14_energy;
  double elec14_energy;
  int vdw_offset;
  double qi, qj, Aij, Bij, Rij;
  int no_pairs;

  individual_sum=0.0;
  excluded_offset=0;
  bond_energy=0.0;
  angle_energy=0.0;
  dihedral_energy=0.0;
  elec_energy=0.0;
  vdw_energy=0.0;
  vdw14_energy=0.0;
  elec14_energy=0.0;
  no_pairs=0;
    
  /*Now do the electrostatics and van der waals - do both together for speed*/
  /*note we avoid double counting in the electrostatics since we do sum(i<j)[qiqj]
                                                                            ----
                                                                            Rij*/
  for (i=0;i<parm_data->NTOTAT;++i) /*loop over all atoms*/
  {
    excluded_temp=0; /*number that have been exlcuded for this atom to date*/
    for (j=i+1;j<parm_data->NTOTAT;++j)
    {
      /*This loops as follows: atom 1 ( 2 to NTOTAT), atom 2 ( 3 to NTOTAT), atom3 ( 4 to NTOTAT) etc.*/
      /*Only calculate if this pair is not excluded in the prmtop file - since it is a 1-2, 1-3 or 1-4.*/
      excluded=NO;
      if (excluded_temp<parm_data->atom[i].numex)
      {
        /*check if this pair is excluded*/
        if ( j == (parm_data->natex[excluded_offset+excluded_temp]-1) )
        {
          ++excluded_temp;
          excluded=YES;
        }
      }
      if (excluded==NO)
      {
        /*not excluded, calculate elec and vdw for this atom pair*/
        ++no_pairs;
        qi = parm_data->atom[i].chrg;
        qj = parm_data->atom[j].chrg;
  
        /*Note, the way Aij and Bij are looked up is quite complex. Essentially we have a diagonal matrix of
          Ntypes^2 in *nno

          For example for NMA there are 7 types, so nno is 28 ints long

          1  2  4  7 11 16 22
          2  3  5  8 12 17 23
          4  5  6  9 13 18 24
          7  8  9 10 14 19 25
          11 12 13 14 15 10 26
          16 17 18 19 20 21 27
          22 23 24 25 26 27 28

          We then look this up based on the atom types of the two atoms in the pair. So, we use the first atom type in the
          pair to get the offset into this array in the form ntypes*atom_type_number. Then we move the number of types represented
          by the integer type of atom 2 further into the array. This gives us the index for the Aij and Bij coefficient arrays.

            E.g if we were looking at atoms 5 and 10 we first get the atom type of 5 which is 3, so our offset is 7 * 2 = 14 which positions
            us at the first element of row 3. Then we check atom 10's type which is 7 so we move to column 7, offset = 7*(type1-1)+(type2 - 1)
            = 20 - Not 21 since we are counting from 0. This returns the number 24 which is used to find the 24th element of the A coefficient
            and B coefficient arrays*/
        vdw_offset=(parm_data->NTYPES*(parm_data->atom[i].iac - 1))+(parm_data->atom[j].iac - 1);
        Aij= parm_data->cn1[(parm_data->nno[vdw_offset]-1)];
        Bij= parm_data->cn2[(parm_data->nno[vdw_offset]-1)];
  
        tempx1=coords_data->x_coord[i];
        tempy1=coords_data->y_coord[i];
        tempz1=coords_data->z_coord[i];
        tempx2=coords_data->x_coord[j];
        tempy2=coords_data->y_coord[j];
        tempz2=coords_data->z_coord[j];
        /*find distance between*/
        Rij = calc_bond_length(tempx1,tempy1,tempz1,tempx2,tempy2,tempz2);

        /*now find the electronic energy for this pair*/
        elec_temp=(qi*qj)/(Rij);
        elec_energy+=elec_temp;
        individual_sum+=elec_temp;
        /*and now the vdw energy*/
        vdw_temp=(Aij / pow(Rij,12)) - (Bij / pow(Rij,6));
        vdw_energy+=vdw_temp;
        individual_sum+=vdw_temp;
      }
    }
    excluded_offset+=parm_data->atom[i].numex;
  }
  
  /*Note, above does not do 1-4 interactions, we do this by considering each dihedral term for which ALL the atom types are +ve*/
  /*dihedrals with H first*/
  for (i=0;i<parm_data->NPHIH;++i)
  {
    /*Check all atoms are +ve*/
    if (parm_data->pdihedralH[i].ip>=0 && parm_data->pdihedralH[i].jp >=0 && parm_data->pdihedralH[i].kp >=0 && parm_data->pdihedralH[i].lp >=0)
    {
        /*all are +ve so calculate 1-4 EE and 1-4 VDW, don't forget to scale*/
        qi = parm_data->atom[unObfuscateAtom(parm_data->pdihedralH[i].ip)-1].chrg;
        qj = parm_data->atom[unObfuscateAtom(parm_data->pdihedralH[i].lp)-1].chrg;
        /*Note, the way Aij and Bij are looked up is quite complex. Essentially we have a diagonal matrix of
          Ntypes^2 in *nno

          For example for NMA there are 7 types, so nno is 28 ints long

          1  2  4  7 11 16 22
          2  3  5  8 12 17 23
          4  5  6  9 13 18 24
          7  8  9 10 14 19 25
          11 12 13 14 15 10 26
          16 17 18 19 20 21 27
          22 23 24 25 26 27 28

          We then look this up based on the atom types of the two atoms in the pair. So, we use the first atom type in the
          pair to get the offset into this array in the form ntypes*atom_type_number. Then we move the number of types represented
          by the integer type of atom 2 further into the array. This gives us the index for the Aij and Bij coefficient arrays.

            E.g if we were looking at atoms 5 and 10 we first get the atom type of 5 which is 3, so our offset is 7 * 2 = 14 which positions
            us at the first element of row 3. Then we check atom 10's type which is 7 so we move to column 7, offset = 7*(type1-1)+(type2 - 1)
            = 20 - Not 21 since we are counting from 0. This returns the number 24 which is used to find the 24th element of the A coefficient
            and B coefficient arrays*/
        vdw_offset=(parm_data->NTYPES*(parm_data->atom[unObfuscateAtom(parm_data->pdihedralH[i].ip)-1].iac - 1))
                    +(parm_data->atom[unObfuscateAtom(parm_data->pdihedralH[i].lp)-1].iac - 1);
        Aij= parm_data->cn1[(parm_data->nno[vdw_offset]-1)];
        Bij= parm_data->cn2[(parm_data->nno[vdw_offset]-1)];
        tempx1=coords_data->x_coord[unObfuscateAtom(parm_data->pdihedralH[i].ip)-1];
        tempy1=coords_data->y_coord[unObfuscateAtom(parm_data->pdihedralH[i].ip)-1];
        tempz1=coords_data->z_coord[unObfuscateAtom(parm_data->pdihedralH[i].ip)-1];
        tempx2=coords_data->x_coord[unObfuscateAtom(parm_data->pdihedralH[i].lp)-1];
        tempy2=coords_data->y_coord[unObfuscateAtom(parm_data->pdihedralH[i].lp)-1];
        tempz2=coords_data->z_coord[unObfuscateAtom(parm_data->pdihedralH[i].lp)-1];
        /*find distance between*/
        Rij = calc_bond_length(tempx1,tempy1,tempz1,tempx2,tempy2,tempz2);

        /*now find the electronic energy for this pair*/
        elec_temp=(qi*qj)/(Rij);
        elec_temp/=global_options->SCEE;         
        elec14_energy+=elec_temp;
        individual_sum+=elec_temp;
        /*and now the vdw energy*/
        vdw_temp=(Aij / pow(Rij,12)) - (Bij / pow(Rij,6));
        vdw_temp/=global_options->SCNB;
        vdw14_energy+=vdw_temp;
  individual_sum+=vdw_temp;
    }
  }
  /*dihedrals without H*/
  for (i=0;i<parm_data->NPHIA;++i)
  {
    /*Check all atoms are +ve*/
    if (parm_data->pdihedral[i].ip>=0 && parm_data->pdihedral[i].jp >=0 && parm_data->pdihedral[i].kp >=0 && parm_data->pdihedral[i].lp >=0)
    {
        /*all are +ve so calculate 1-4 EE and 1-4 VDW, don't forget to scale*/
        /*all are +ve so calculate 1-4 EE and 1-4 VDW, don't forget to scale*/
        qi = parm_data->atom[unObfuscateAtom(parm_data->pdihedral[i].ip)-1].chrg;
        qj = parm_data->atom[unObfuscateAtom(parm_data->pdihedral[i].lp)-1].chrg;
        /*Note, the way Aij and Bij are looked up is quite complex. Essentially we have a diagonal matrix of
          Ntypes^2 in *nno

          For example for NMA there are 7 types, so nno is 28 ints long

          1  2  4  7 11 16 22
          2  3  5  8 12 17 23
          4  5  6  9 13 18 24
          7  8  9 10 14 19 25
          11 12 13 14 15 10 26
          16 17 18 19 20 21 27
          22 23 24 25 26 27 28

          We then look this up based on the atom types of the two atoms in the pair. So, we use the first atom type in the
          pair to get the offset into this array in the form ntypes*atom_type_number. Then we move the number of types represented
          by the integer type of atom 2 further into the array. This gives us the index for the Aij and Bij coefficient arrays.

            E.g if we were looking at atoms 5 and 10 we first get the atom type of 5 which is 3, so our offset is 7 * 2 = 14 which positions
            us at the first element of row 3. Then we check atom 10's type which is 7 so we move to column 7, offset = 7*(type1-1)+(type2 - 1)
            = 20 - Not 21 since we are counting from 0. This returns the number 24 which is used to find the 24th element of the A coefficient
            and B coefficient arrays*/
        vdw_offset=(parm_data->NTYPES*(parm_data->atom[unObfuscateAtom(parm_data->pdihedral[i].ip)-1].iac - 1))
                    +(parm_data->atom[unObfuscateAtom(parm_data->pdihedral[i].lp)-1].iac - 1);
        Aij= parm_data->cn1[(parm_data->nno[vdw_offset]-1)];
        Bij= parm_data->cn2[(parm_data->nno[vdw_offset]-1)];
        tempx1=coords_data->x_coord[unObfuscateAtom(parm_data->pdihedral[i].ip)-1];
        tempy1=coords_data->y_coord[unObfuscateAtom(parm_data->pdihedral[i].ip)-1];
        tempz1=coords_data->z_coord[unObfuscateAtom(parm_data->pdihedral[i].ip)-1];
        tempx2=coords_data->x_coord[unObfuscateAtom(parm_data->pdihedral[i].lp)-1];
        tempy2=coords_data->y_coord[unObfuscateAtom(parm_data->pdihedral[i].lp)-1];
        tempz2=coords_data->z_coord[unObfuscateAtom(parm_data->pdihedral[i].lp)-1];
        /*find distance between*/
        Rij = calc_bond_length(tempx1,tempy1,tempz1,tempx2,tempy2,tempz2);

        /*now find the electronic energy for this pair*/
        elec_temp=(qi*qj)/(Rij);
        elec_temp/=global_options->SCEE;
        elec14_energy+=elec_temp;
        individual_sum+=elec_temp;
        /*and now the vdw energy*/
        vdw_temp=(Aij / pow(Rij,12)) - (Bij / pow(Rij,6));
        vdw_temp/=global_options->SCNB;
        vdw14_energy+=vdw_temp;
  individual_sum+=vdw_temp;
    }     
  }
  
  /*Now do the bonds*/
  for (i=0;i<parm_data->unique_bonds_found;++i)
  {
    for (j=0;j<parm_data->bond_data[i].number;++j)
    {
      /*First off we need to find the distance between the atoms in the bond*/
      /*Get the coordinates*/
      tempx1=coords_data->x_coord[parm_data->bond_data[i].atom1[j]-1];
      tempy1=coords_data->y_coord[parm_data->bond_data[i].atom1[j]-1];
      tempz1=coords_data->z_coord[parm_data->bond_data[i].atom1[j]-1];
      tempx2=coords_data->x_coord[parm_data->bond_data[i].atom2[j]-1];
      tempy2=coords_data->y_coord[parm_data->bond_data[i].atom2[j]-1];
      tempz2=coords_data->z_coord[parm_data->bond_data[i].atom2[j]-1];

      /*find bond length*/
      bond_length_temp=calc_bond_length(tempx1,tempy1,tempz1,tempx2,tempy2,tempz2);;
      /*Now calculate the energy contribution for this bond and add it to the individual sum*/
      individual_sum+=(parm_data->bond_data[i].rk)*pow(bond_length_temp-parm_data->bond_data[i].req,2);
                            /*BONDjFC*/                    /*R*/                         /*Req*/
      if (global_options->VERBOSITY>=HIGH)
            bond_energy+=(parm_data->bond_data[i].rk)*pow(bond_length_temp-parm_data->bond_data[i].req,2);
      
    }
  }

  /*Angles*/
  for (i=0;i<parm_data->unique_angles_found;++i)
  {
    for (j=0;j<parm_data->angle_data[i].number;++j)
    {
      /*First off we need to find the angle between the atoms in the angle*/
      /*Get the coordinates*/
      tempx1=coords_data->x_coord[parm_data->angle_data[i].atom1[j]-1];
      tempy1=coords_data->y_coord[parm_data->angle_data[i].atom1[j]-1];
      tempz1=coords_data->z_coord[parm_data->angle_data[i].atom1[j]-1];
      tempx2=coords_data->x_coord[parm_data->angle_data[i].atom2[j]-1];
      tempy2=coords_data->y_coord[parm_data->angle_data[i].atom2[j]-1];
      tempz2=coords_data->z_coord[parm_data->angle_data[i].atom2[j]-1];
      tempx3=coords_data->x_coord[parm_data->angle_data[i].atom3[j]-1];
      tempy3=coords_data->y_coord[parm_data->angle_data[i].atom3[j]-1];
      tempz3=coords_data->z_coord[parm_data->angle_data[i].atom3[j]-1];

      /*find angle*/
      angle_temp=calc_angle_radians(tempx1,tempy1,tempz1,tempx2,tempy2,tempz2,tempx3,tempy3,tempz3);
      /*Now calculate the energy contribution for this angle and add it to the individual sum*/
      individual_sum+=(parm_data->angle_data[i].tk)*pow(angle_temp-parm_data->angle_data[i].teq,2);
                            /*ANGLEjFC*/                    /*T*/                         /*Teq*/
            if (global_options->VERBOSITY>=HIGH)
        angle_energy+=(parm_data->angle_data[i].tk)*pow(angle_temp-parm_data->angle_data[i].teq,2);    
    }
  }
  
  /*Dihedrals*/
  for (i=0;i<parm_data->unique_dihedrals_found;++i)
  {
    for (k=0; k<parm_data->dihedral_data[i].number; ++k) {
        /*First off we need to find the dihedral between the atoms in the angle*/
        /*Get the coordinates*/
        tempx1=coords_data->x_coord[parm_data->dihedral_data[i].atom1[k]-1];
        tempy1=coords_data->y_coord[parm_data->dihedral_data[i].atom1[k]-1];
        tempz1=coords_data->z_coord[parm_data->dihedral_data[i].atom1[k]-1];
        tempx2=coords_data->x_coord[parm_data->dihedral_data[i].atom2[k]-1];
        tempy2=coords_data->y_coord[parm_data->dihedral_data[i].atom2[k]-1];
        tempz2=coords_data->z_coord[parm_data->dihedral_data[i].atom2[k]-1];
        tempx3=coords_data->x_coord[parm_data->dihedral_data[i].atom3[k]-1];
        tempy3=coords_data->y_coord[parm_data->dihedral_data[i].atom3[k]-1];
        tempz3=coords_data->z_coord[parm_data->dihedral_data[i].atom3[k]-1];
        tempx4=coords_data->x_coord[parm_data->dihedral_data[i].atom4[k]-1];
        tempy4=coords_data->y_coord[parm_data->dihedral_data[i].atom4[k]-1];
        tempz4=coords_data->z_coord[parm_data->dihedral_data[i].atom4[k]-1];

        /*find dihedral*/
        dihedral_temp=calc_dihedral_radians(tempx1,tempy1,tempz1,tempx2,tempy2,tempz2,tempx3,tempy3,tempz3,tempx4,tempy4,tempz4);
        dihedral_temp=PI-dihedral_temp; /*Due to amber calculating the angle in the opposite direction to me*/
        
      for (j=0;j<parm_data->dihedral_data[i].num_terms;++j) {
        /*Now calculate the energy contribution for this dihedral and add it to the individual sum*/
        /*Note, pk is actually pk/2*/
        individual_sum+=(parm_data->dihedral_data[i].term[j].pk)*
                        (1+cos(dihedral_temp*parm_data->dihedral_data[i].term[j].pn-parm_data->dihedral_data[i].term[j].phase));
        if (global_options->VERBOSITY>=HIGH)
          dihedral_energy+=(parm_data->dihedral_data[i].term[j].pk)*
                          (1+cos(dihedral_temp*parm_data->dihedral_data[i].term[j].pn-parm_data->dihedral_data[i].term[j].phase));   
      }
    }
  }
  if(global_options->VERBOSITY>=HIGH) {
    printf("AMBER STD indiv = %.4f, (BND=%.4f, ANG=%.4f, DIHE=%.4f, EEL=%.4f,\n                         1-4EEL=%.4f, VDW=%.4f, 1-4VDW=%.4f, NPAIRS=%d)\n"
    ,individual_sum,bond_energy,angle_energy,dihedral_energy, elec_energy,elec14_energy, vdw_energy, vdw14_energy, no_pairs);
    fflush(stdout); /*Flush the printf buffer*/
  }
  return individual_sum;

}


