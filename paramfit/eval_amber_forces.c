/*****************************************************
 * AMBER Bond Angle and Dihedral Parameter Optimiser *
 *                                                   *
 *           Written by: Robin Betz  (2012)          *
 *                   UC San Diego                    *
 *           San Diego Supercomputer Center          *
 *            La Jolla, California, 92092            *
 *                       USA                         *
 *****************************************************/
/** @file eval_amber_forces.c
 * Functions for evaluating forces on each atom using the AMBER equation.
 * This has been taken pretty much directly from the pmemd code.
 */

#include "function_def.h"
#include "constants.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/**
 * Calculates forces on a single input structure.
 * Uses the parameters in the parm struct and returns forces in the force struct for each atom.
 * 
 * @param[in] global_options The global options structure
 * @param[in] parm_struct The parameter structure with force constants and equilibria to use
 * @param[in] coords_data Pointer to a single coordinate data set
 * @param[in,out] forces Pre-allocated to array[NATOMS], will be populated with forces on each atom
 * @param[in] structure Index of the structure in coords_data to use for atom positions
 * @return Integer representing success or failure
 */
int eval_amber_forces_single_struct(global_options_struct *global_options, parm_struct *parm_data, coord_set *coords_data, force_struct *forces, int structure)
{
  int i, j, k;
  double df;
  // Zero the forces to begin with
  for (i=0; i<coords_data->natoms; ++i) {
    forces[i].x=0.0;
    forces[i].y=0.0;
    forces[i].z=0.0;
  }
  // Calculate the non-bonded contributions
  int excluded_offset=0, excluded=NO, excluded_temp=0;
  
  for (i=0;i<coords_data->natoms;++i)
  {
    excluded_temp=0;
    for (j=i+1;j<coords_data->natoms;++j)
    {
      //This loops as follows: atom 1 ( 2 to NTOTAT), atom 2 ( 3 to NTOTAT), atom3 ( 4 to NTOTAT) etc.
      //Only calculate if this pair is not excluded in the prmtop file - since it is a 1-2, 1-3 or 1-4.
      excluded=NO;
      if (excluded_temp<parm_data->atom[i].numex)
      {
        //check if this pair is excluded
        if ( j == (parm_data->natex[excluded_offset+excluded_temp]-1) )
        {
          ++excluded_temp;
          excluded=YES;
        }
      }
      if (excluded==NO)
      {
        //not excluded, calculate elec and vdw for this atom pair
        double qi = parm_data->atom[i].chrg;
        double qj = parm_data->atom[j].chrg;
        int vdw_offset=(parm_data->NTYPES*(parm_data->atom[i].iac - 1))+(parm_data->atom[j].iac - 1);
        int Aij= parm_data->cn1[(parm_data->nno[vdw_offset]-1)];
        int Bij= parm_data->cn2[(parm_data->nno[vdw_offset]-1)];
        double tempx1=coords_data->struc[structure].x_coord[i];
        double tempy1=coords_data->struc[structure].y_coord[i];
        double tempz1=coords_data->struc[structure].z_coord[i];
        double tempx2=coords_data->struc[structure].x_coord[j];
        double tempy2=coords_data->struc[structure].y_coord[j];
        double tempz2=coords_data->struc[structure].z_coord[j];
        double Rij = calc_bond_length(tempx1,tempy1,tempz1,tempx2,tempy2,tempz2);
        double Rx = (tempx2-tempx1)/(Rij*Rij); // components of normalized vector between the two atoms
        double Ry = (tempy2-tempy1)/(Rij*Rij);
        double Rz = (tempz2-tempz1)/(Rij*Rij);
        
        double magnitude = (qi*qj)/(Rij*Rij);
        double fx = magnitude * Rx;
        double fy = magnitude * Ry;
        double fz = magnitude * Rz;
        
        // Force is relative to atom1: if the charges are opposite, positive force on 1,
        // negative on atom 2. The reverse is true if the force is repulsive due to same charges.
        // So if fx is negative, move atom i towards atom j
        forces[i].x -= fx;
        forces[i].y -= fy;
        forces[i].z -= fz;
        forces[j].x += fx;
        forces[j].y += fy;
        forces[j].z += fz;
        
        // Now calculate Lennard-Jones forces.
        magnitude = -6*( (2*Aij / pow(Rij,13)) - (Bij / pow(Rij,7)) );
        forces[i].x += magnitude * Rx;
        forces[i].y += magnitude * Ry;
        forces[i].z += magnitude * Rz;
        forces[j].x -= magnitude * Rx;
        forces[j].y -= magnitude * Ry;
        forces[j].z -= magnitude * Rz;
      }
    }
    excluded_offset+=parm_data->atom[i].numex;
  }
  // Do 1-4 dihedral interactions now, starting with dihedrals with H
  for (i=0; i<parm_data->NPHIH; ++i)
  {
    //Check all atoms are +ve
    if (parm_data->pdihedralH[i].ip>=0 && parm_data->pdihedralH[i].jp >=0 && parm_data->pdihedralH[i].kp >=0 && parm_data->pdihedralH[i].lp >=0)
    {
      //all are +ve so calculate 1-4 EE and 1-4 VDW, don't forget to scale
      int a = unObfuscateAtom(parm_data->pdihedralH[i].ip)-1;
      int b = unObfuscateAtom(parm_data->pdihedralH[i].lp)-1;
      double qi = parm_data->atom[a].chrg;
      double qj = parm_data->atom[b].chrg;
      int vdw_offset=(parm_data->NTYPES*(parm_data->atom[a].iac - 1)) +
                     (parm_data->atom[b].iac - 1);
      double Aij= parm_data->cn1[(parm_data->nno[vdw_offset]-1)];
      double Bij= parm_data->cn2[(parm_data->nno[vdw_offset]-1)];
      double tempx1=coords_data->struc[structure].x_coord[a];
      double tempy1=coords_data->struc[structure].y_coord[a];
      double tempz1=coords_data->struc[structure].z_coord[a];
      double tempx2=coords_data->struc[structure].x_coord[b];
      double tempy2=coords_data->struc[structure].y_coord[b];
      double tempz2=coords_data->struc[structure].z_coord[b];
      double Rij = calc_bond_length(tempx1,tempy1,tempz1,tempx2,tempy2,tempz2); // distance
      double Rx = (tempx2-tempx1)/(Rij*Rij); // components of normalized vector between the two atoms
      double Ry = (tempy2-tempy1)/(Rij*Rij);
      double Rz = (tempz2-tempz1)/(Rij*Rij);
      
      double magnitude = (qi*qj)/(Rij*Rij) / global_options->SCEE; // <-- note scaling factor
      double fx = magnitude * Rx;
      double fy = magnitude * Ry;
      double fz = magnitude * Rz;
      
      forces[a].x -= fx;
      forces[a].y -= fy;
      forces[a].z -= fz;
      forces[b].x += fx;
      forces[b].y += fy;
      forces[b].z += fz;
      
      // Now for the VDW
      magnitude = -6*( (2*Aij / pow(Rij,13)) - (Bij / pow(Rij,7)) ) / global_options->SCNB;
      forces[a].x += magnitude * Rx;
      forces[a].y += magnitude * Ry;
      forces[a].z += magnitude * Rz;
      forces[b].x -= magnitude * Rx;
      forces[b].y -= magnitude * Ry;
      forces[b].z -= magnitude * Rz;
    }
  }
  // Now dihedrals without H
  for (i=0;i<parm_data->NPHIA;++i)
  {
    //Check all atoms are +ve
    if (parm_data->pdihedral[i].ip>=0 && parm_data->pdihedral[i].jp >=0 && parm_data->pdihedral[i].kp >=0 && parm_data->pdihedral[i].lp >=0)
    {
      //all are +ve so calculate 1-4 EE and 1-4 VDW, don't forget to scale
      int a = unObfuscateAtom(parm_data->pdihedral[i].ip)-1;
      int b = unObfuscateAtom(parm_data->pdihedral[i].lp)-1;
      double qi = parm_data->atom[a].chrg;
      double qj = parm_data->atom[b].chrg;
      int vdw_offset=(parm_data->NTYPES*(parm_data->atom[a].iac - 1)) +
                     (parm_data->atom[b].iac - 1);
      double Aij= parm_data->cn1[(parm_data->nno[vdw_offset]-1)];
      double Bij= parm_data->cn2[(parm_data->nno[vdw_offset]-1)];
      double tempx1=coords_data->struc[structure].x_coord[a];
      double tempy1=coords_data->struc[structure].y_coord[a];
      double tempz1=coords_data->struc[structure].z_coord[a];
      double tempx2=coords_data->struc[structure].x_coord[b];
      double tempy2=coords_data->struc[structure].y_coord[b];
      double tempz2=coords_data->struc[structure].z_coord[b];
      double Rij = calc_bond_length(tempx1,tempy1,tempz1,tempx2,tempy2,tempz2); // distance
      double Rx = (tempx2-tempx1)/(Rij*Rij); // components of normalized vector between the two atoms
      double Ry = (tempy2-tempy1)/(Rij*Rij);
      double Rz = (tempz2-tempz1)/(Rij*Rij);
      
      double magnitude = (qi*qj)/(Rij*Rij) / global_options->SCEE; // <-- note scaling factor
      double fx = magnitude * Rx;
      double fy = magnitude * Ry;
      double fz = magnitude * Rz;
      
      forces[a].x -= fx;
      forces[a].y -= fy;
      forces[a].z -= fz;
      forces[b].x += fx;
      forces[b].y += fy;
      forces[b].z += fz;
      
      // Now for the VDW
      magnitude = -6*( (2*Aij / pow(Rij,13)) - (Bij / pow(Rij,7)) ) / global_options->SCNB;
      forces[a].x += magnitude * Rx;
      forces[a].y += magnitude * Ry;
      forces[a].z += magnitude * Rz;
      forces[b].x -= magnitude * Rx;
      forces[b].y -= magnitude * Ry;
      forces[b].z -= magnitude * Rz;
    }
  }
//   printf("electrostatics print:\n");
//   for (i=0; i<global_options->NATOMS; ++i)
//     printf("Atom %i: ( %f, %f, %f ) Force %.4f\t%.4f\t%.4f\n", i, coords_data[structure].x_coord[i], coords_data[structure].y_coord[i], coords_data[structure].x_coord[i],
//            forces[i].x, forces[i].y, forces[i].z);
//          printf("\n");
  // Bonds, from bonds.F90
  double rij;
  double tempx1, tempy1, tempz1, tempx2, tempy2, tempz2;
  for (i=0;i<parm_data->unique_bonds_found;++i) {
    for (j=0;j<parm_data->bond_data[i].number;++j) {
      // Get all of the atoms in the bond
      int I = parm_data->bond_data[i].atom1[j]-1;
      int J = parm_data->bond_data[i].atom2[j]-1;
      tempx1=coords_data->struc[structure].x_coord[I];
      tempy1=coords_data->struc[structure].y_coord[I];
      tempz1=coords_data->struc[structure].z_coord[I];
      tempx2=coords_data->struc[structure].x_coord[J];
      tempy2=coords_data->struc[structure].y_coord[J];
      tempz2=coords_data->struc[structure].z_coord[J];
      // Calculate the bond length and force
      rij=calc_bond_length(tempx1,tempy1,tempz1,tempx2,tempy2,tempz2);
      double df = (rij-parm_data->bond_data[i].req)*parm_data->bond_data[i].rk;
      // Update forces
      forces[I].x -= (df+df)/rij * (tempx1-tempx2);
      forces[I].y -= (df+df)/rij * (tempy1-tempy2);
      forces[I].z -= (df+df)/rij * (tempz1-tempz2);
      forces[J].x += (df+df)/rij * (tempx1-tempx2);
      forces[J].y += (df+df)/rij * (tempy1-tempy2);
      forces[J].z += (df+df)/rij * (tempz1-tempz2);
    }
  }
  
  // Angles, from angles.F90
  int a, b, c;
  double rXab, rYab, rZab, rXcb, rYcb, rZcb;
  double cac, sth, caa, ccc;
  for (i=0;i<parm_data->unique_angles_found;++i)
  {
    for (j=0;j<parm_data->angle_data[i].number;++j)
    {
      // Get distances between atoms
      a = parm_data->angle_data[i].atom1[j]-1;
      b = parm_data->angle_data[i].atom2[j]-1;
      c = parm_data->angle_data[i].atom3[j]-1;
      rXab = coords_data->struc[structure].x_coord[a] - coords_data->struc[structure].x_coord[b];
      rYab = coords_data->struc[structure].y_coord[a] - coords_data->struc[structure].y_coord[b];
      rZab = coords_data->struc[structure].z_coord[a] - coords_data->struc[structure].z_coord[b];
      rXcb = coords_data->struc[structure].x_coord[c] - coords_data->struc[structure].x_coord[b];
      rYcb = coords_data->struc[structure].y_coord[c] - coords_data->struc[structure].y_coord[b];
      rZcb = coords_data->struc[structure].z_coord[c] - coords_data->struc[structure].z_coord[b];
      double rab = rXab*rXab + rYab*rYab + rZab*rZab;
      double rcb = rXcb*rXcb + rYcb*rYcb + rZcb*rZcb;
      double rac = sqrt(rab*rcb);
      
      // Get angle force
      double cst = fmin(0.999, fmax(-0.999, (rXab*rXcb + rYab*rYcb + rZab*rZcb)/rac));      
      double ant = acos(cst);
      df = -2*( (ant - parm_data->angle_data[i].teq) * parm_data->angle_data[i].tk ) / sin(ant);
      cac = df / rac;
      sth = df * cst;
      caa = sth / rab;
      ccc = sth / rcb;
      
      // Calculate overall forces
      double dt1 = cac * rXcb - caa * rXab;
      double dt2 = cac * rYcb - caa * rYab;
      double dt3 = cac * rZcb - caa * rZab;
      double dt4 = cac * rXab - ccc * rXcb;
      double dt5 = cac * rYab - ccc * rYcb;
      double dt6 = cac * rZab - ccc * rZcb;
      
      // Update forces
      forces[a].x -= dt1;
      forces[a].y -= dt2;
      forces[a].z -= dt3;
      forces[b].x += (dt1 + dt4);
      forces[b].y += (dt2 + dt5);
      forces[b].z += (dt3 + dt6);
      forces[c].x -= dt4;
      forces[c].y -= dt5;
      forces[c].z -= dt6;
    }
  }
 
//  printf("angles print:\n");
//    printf("Atom %i: ( %f, %f, %f ) Force %.4f\t%.4f\t%.4f\n", i, coords_data[structure].x_coord[i], coords_data[structure].y_coord[i], coords_data[structure].x_coord[i],
//           forces[i].x, forces[i].y, forces[i].z);
// printf("\n");

  // Dihedrals, from dihedrals.F90
  for (i=0;i<parm_data->unique_dihedrals_found;++i) {
    for(k=0;k<parm_data->dihedral_data[i].number; ++k) {
        // Get the indices for atoms in this dihedral
        int I = parm_data->dihedral_data[i].atom1[k]-1;
        int J = parm_data->dihedral_data[i].atom2[k]-1;
        int K = parm_data->dihedral_data[i].atom3[k]-1;
        int L = parm_data->dihedral_data[i].atom4[k]-1;
        
        // Get relevant distances between atoms
        double xij = coords_data->struc[structure].x_coord[I] - coords_data->struc[structure].x_coord[J];
        double yij = coords_data->struc[structure].y_coord[I] - coords_data->struc[structure].y_coord[J];
        double zij = coords_data->struc[structure].z_coord[I] - coords_data->struc[structure].z_coord[J];
        double xkj = coords_data->struc[structure].x_coord[K] - coords_data->struc[structure].x_coord[J];
        double ykj = coords_data->struc[structure].y_coord[K] - coords_data->struc[structure].y_coord[J];
        double zkj = coords_data->struc[structure].z_coord[K] - coords_data->struc[structure].z_coord[J];
        double xkl = coords_data->struc[structure].x_coord[K] - coords_data->struc[structure].x_coord[L];
        double ykl = coords_data->struc[structure].y_coord[K] - coords_data->struc[structure].y_coord[L];
        double zkl = coords_data->struc[structure].z_coord[K] - coords_data->struc[structure].z_coord[L];
        
        // Get the normal vector
        double dx = yij * zkj - zij * ykj;
        double dy = zij * xkj - xij * zkj;
        double dz = xij * ykj - yij * xkj;
        double gx = zkj * ykl - ykj * zkl;
        double gy = xkj * zkl - zkj * xkl;
        double gz = ykj * xkl - xkj * ykl;
        double fxi =  sqrt(dx * dx + dy * dy + dz * dz + 1e-18);
        double fyi = sqrt(gx * gx + gy * gy + gz * gz + 1e-18);
        double ct = dx * gx + dy * gy + dz * gz;
        
        // Handle linear dihedral
        double z1, z2, fzi;
        if (fxi > 1e-3) z1 = 1.0 / fxi;
        else z1 = 0.0;
        if (fyi > 1e-3) z2 = 1.0/fyi;
        else z2 = 0.0;
        double z12 = z1*z2;
        if (z12 != 0) fzi = 1.0;
        else fzi = 0;
        double s = xkj * (dz * gy - dy * gz) + ykj * (dx * gz - dz * gx) +
                  zkj * (dy * gx - dx * gy);
        double ap;
        if (s<0) ap = PI + fabs( acos(fmax(-1.0, fmin(1.0, ct * z12))) ); // Pi -> same sign as s
        else ap = PI - fabs( acos(fmax(-1.0, fmin(1.0, ct * z12))) );
        double cphi = cos(ap);
        double sphi = sin(ap);
        double phase = calc_dihedral_radians(coords_data->struc[structure].x_coord[I], coords_data->struc[structure].y_coord[I], coords_data->struc[structure].z_coord[I], 
                                            coords_data->struc[structure].x_coord[J], coords_data->struc[structure].y_coord[J], coords_data->struc[structure].z_coord[J],
                                            coords_data->struc[structure].x_coord[K], coords_data->struc[structure].y_coord[K], coords_data->struc[structure].z_coord[K],
                                            coords_data->struc[structure].x_coord[L], coords_data->struc[structure].y_coord[L], coords_data->struc[structure].z_coord[L]);
    for (j=0;j<parm_data->dihedral_data[i].num_terms;++j) {
        double gamc = parm_data->dihedral_data[i].term[j].pk*cos(phase);
        double gams = parm_data->dihedral_data[i].term[j].pk*sin(phase);
        
        // Calculate the first derivatives with respect to cos(phi)
        int inc = (int)(parm_data->dihedral_data[i].term[j].pn+1.0e-3);
        double cosnp = cos(parm_data->dihedral_data[i].term[j].pn *ap);
        double sinnp = sin(parm_data->dihedral_data[i].term[j].pn*ap);
        double epw = (parm_data->dihedral_data[i].term[j].pk + cosnp*gamc + sinnp*gams)*fzi;
        double dums = sphi >= 0 ? sphi + 1.0e-18 : sphi - 1.0e-18;
        double df;
        double gmul = ((int)parm_data->dihedral_data[i].term[j].pn)&1 ? 0.0 : parm_data->dihedral_data[i].term[j].pn;
        if (fabs(dums) < 1.0e-6)
          df = fzi*gamc*(parm_data->dihedral_data[i].term[j].pn - gmul + gmul*cphi);
        else
          df = fzi*parm_data->dihedral_data[i].term[j].pn*(gamc*sinnp - gams*cosnp)/dums;
        
        // Set up the "first derivatives of cos(phi) with respect to cartesian differences"
        double z11 = z1 * z1;
        z12 = z1 * z2;
        double z22 = z2 * z2;
        double dc1 = -gx * z12 - cphi * dx * z11;
        double dc2 = -gy * z12 - cphi * dy * z11;
        double dc3 = -gz * z12 - cphi * dz * z11;
        double dc4 =  dx * z12 + cphi * gx * z22;
        double dc5 =  dy * z12 + cphi * gy * z22;
        double dc6 =  dz * z12 + cphi * gz * z22;
        
        // Update the first derivative array
        double dr1 = df * ( dc3 * ykj - dc2 * zkj);
        double dr2 = df * ( dc1 * zkj - dc3 * xkj);
        double dr3 = df * ( dc2 * xkj - dc1 * ykj);
        double dr4 = df * ( dc6 * ykj - dc5 * zkj);
        double dr5 = df * ( dc4 * zkj - dc6 * xkj);
        double dr6 = df * ( dc5 * xkj - dc4 * ykj);
        double drx = df * (-dc2 * zij + dc3 * yij + dc5 * zkl - dc6 * ykl);
        double dry = df * ( dc1 * zij - dc3 * xij - dc4 * zkl + dc6 * xkl);
        double drz = df * (-dc1 * yij + dc2 * xij + dc4 * ykl - dc5 * xkl);
        
        
        double fxl = -dr4;
        double fyl = -dr5;
        double fzl = -dr6;
        
        // Sum up the forces
        forces[I].x -= dr1;
        forces[I].y -= dr2;
        forces[I].z -= dr3;
        forces[J].x += (-drx + dr1);
        forces[J].y += (-dry + dr2);
        forces[J].z += (-drz + dr3);
        forces[K].x += (drx + dr4);
        forces[K].y += (dry + dr5);
        forces[K].z += (drz + dr6);
        forces[L].x -= dr4;
        forces[L].y -= dr5;
        forces[L].z -= dr6;
      }
    }
  }
//   printf("Dihedrals print:\n");
//   for (i=0; i<global_options->NATOMS; ++i)
//     printf("Atom %i: ( %f, %f, %f ) Force %.4f\t%.4f\t%.4f\n", i, coords_data[structure].x_coord[i], coords_data[structure].y_coord[i], coords_data[structure].x_coord[i],
//            forces[i].x, forces[i].y, forces[i].z);
//          printf("\n");
//          
  // DEBUG
//   for (i=0; i<coords_data->natoms; ++i)
//     printf("%8.4f %8.4f %8.4f\n",forces[i].x, forces[i].y, forces[i].z);
  return SUCCESS;
}

/**
 * Scores a set of parameters cross multiple molecules according to force comparison.
 * Essentially just sums the force evaluation value over all prmtops. This is safe to use if there
 * is only one molecule so is recommended for all force evaluation calls.
 * @param[in] global_options The global options structure
 * @param[in] parm_datas Array of parameter sets, one for each molecule
 * @param[in] coords_datas Array of coordinate data sets, one for each molecule
 * @return Double representing average difference in force magnitude over all molecules and structures
 */
double eval_sum_amber_forces_multiprmtop(global_options_struct *global_options, parm_struct *parm_datas, coord_set *coords_datas)
{
  int i;
  double score = 0.0;
  for (i=0; i<global_options->num_prmtops; ++i) {
    score += eval_sum_amber_forces(global_options, &parm_datas[i], &coords_datas[i]);
  }
  ++global_options->function_calls;
  
  return score;
}

/**
 * Scores a set of parameters according to force comparison.
 * Allocates all necessary force structures, runs the force calculation (in parallel with openmp)
 * compares the result. The comparison is the average difference in force magnitude per structure.
 * 
 * @param[in] global_options The global options structure
 * @param[in] parm_data The parameter set to evaluate
 * @param[in] coords_data Coordinate set containing atom coordinates for all input structures, and QM forces
 * @return Scalar representing average difference in force magnitude per structure.
 */
double eval_sum_amber_forces(global_options_struct *global_options, parm_struct *parm_data, coord_set *coords_data) 
{
  // Before anything else is done, check that the parameters to be evaluated
  // are valid, and if not, correct them.
  check_range(global_options, parm_data);

  int i, j;
  double score= 0.0;
  
  // For every input conformation, do the force evaluation
  force_struct **eval = (force_struct **) malloc(coords_data->num_coords*sizeof(force_struct*));
  if (eval==NULL) exit(ALLOC_FAIL);
  for (i=0; i<coords_data->num_coords; ++i) {
    eval[i] = (force_struct *) malloc(coords_data->natoms*sizeof(force_struct));
    if (eval[i]==NULL) {
      printf("*** ERROR allocating memory for temporary forces for structure %i\n", i);
      exit(ALLOC_FAIL);
    }
  }
  #pragma omp parallel for
  for (i=0; i<coords_data->num_coords; ++i) 
  {
    if (eval_amber_forces_single_struct(global_options, parm_data, coords_data, eval[i], i) != SUCCESS) {
      printf("*** ERROR evaluating forces for structure %i\n", i);
      exit(FAILURE);
    }
  }
 // do this in serial 
  for (i=0; i<coords_data->num_coords; ++i) 
  {
    // Subtract the quantum value for the forces from the amber value to get the difference in force for each atom for this conformation
    // and add the absolute value of that to our vector sum for each atom
    /* DO NOT move this before the force evaluation because eval_amber_forces will zero the input force struct! */
    double structure_result= 0.0;
    for (j=0; j<coords_data->natoms; ++j) {
      if (parm_data->fit_atom[j] == YES) {
        // Calculate the angle between the two vectors
        double xd = ( eval[i][j].x * coords_data->struc[i].force[j].x );
        double yd = ( eval[i][j].y * coords_data->struc[i].force[j].y );
        double zd = ( eval[i][j].z * coords_data->struc[i].force[j].z );
        double theta = (xd+yd+zd) / sqrt(eval[i][j].x*eval[i][j].x + eval[i][j].y*eval[i][j].y + eval[i][j].z*eval[i][j].z);  // NaNs prolly from here
        // Now calculate the difference in magnitude
        xd = ( eval[i][j].x - coords_data->struc[i].force[j].x );
        yd = ( eval[i][j].y - coords_data->struc[i].force[j].y );
        zd = ( eval[i][j].z - coords_data->struc[i].force[j].z );
    //    double mag = sqrt(xd*xd + yd*yd + zd*zd);
        // Resulting fitness is the angle (scaled so 0 is best), times the magnitude
//        structure_result += (theta + mag);
        structure_result += (xd*xd + yd*yd + zd*zd); // squared magnitude of the difference bettween the two
  
      }
    }
    score += structure_result; // add thread result to energy
    
    // Debug output
  if (global_options->VERBOSITY>=DEBUG) {
     printf("After evaluating struct %i:\n", i);
     int r;
     for (r=0; r<coords_data->natoms; ++r)
       printf("Atom %i: ( %f, %f, %f ) Force %6.4f %6.4f %6.4f\n", r, coords_data->struc[i].x_coord[r], coords_data->struc[i].y_coord[r], coords_data->struc[i].x_coord[r],
            coords_data->struc[i].force[r].x, coords_data->struc[i].force[r].y, coords_data->struc[i].force[r].z);
  }
  }
  
  for (i=0; i<coords_data->num_coords; ++i) free(eval[i]);
  free(eval);
  
  // Average the score to return average score per structure since this is more meaningful to the user than some huge number
  score /= (double)coords_data->num_coords; 
  // Return the magnitude of the force difference for the atoms because current fit algorithms would like scalars
  return score;
}

/** Marks atoms that are involved in parameters to be fit.
 * This is done when fitting forces so that only the forces on atoms involved in bonds,
 * angles, or dihedrals to be fit are considered in the function evaluation. This marks those
 * atoms so they may be easily accessed later.
 *
 * @param[in] global_options The global options structure
 * @param[in,out] parm_data The parameter data, including parameters to fit. Atoms will be marked in here.
 * @return Integer indicating success or failure
 */
int mark_relevant_atoms(global_options_struct *global_options, parm_struct *parm_data)
{
  // Initialize to no atoms marked
  memset(parm_data->fit_atom, 0, parm_data->NTOTAT*sizeof(int));
  
  // Mark all atoms involved in fitted bonds
  int i, j, k;
  for (i=0; i<parm_data->unique_bonds_found; ++i) {
    if (parm_data->bond_data[i].DO_FIT_KR==YES || parm_data->bond_data[i].DO_FIT_REQ==YES) {
      for (j=0; j<parm_data->bond_data[i].number; ++j) {
        parm_data->fit_atom[ parm_data->bond_data[i].atom1[j] ] = YES;
        parm_data->fit_atom[ parm_data->bond_data[i].atom2[j] ] = YES;
      }
    }
  }
  
  // Mark atoms involved in fitted angles
  for (i=0; i<parm_data->unique_angles_found; ++i) {
    if (parm_data->angle_data[i].DO_FIT_KT==YES || parm_data->angle_data[i].DO_FIT_THEQ==YES) {
      for (j=0; j<parm_data->angle_data[i].number; ++j) {
        parm_data->fit_atom[ parm_data->angle_data[i].atom1[j] ] = YES;
        parm_data->fit_atom[ parm_data->angle_data[i].atom2[j] ] = YES;
        parm_data->fit_atom[ parm_data->angle_data[i].atom3[j] ] = YES;
      }
    }
  }
  
  // Mark atoms involved in fitted dihedrals
  for (i=0; i<parm_data->unique_dihedrals_found; ++i) {
    for (j=0; j<parm_data->dihedral_data[i].num_terms; ++j) {
      if (parm_data->dihedral_data[i].term[j].DO_FIT_KP==YES || parm_data->dihedral_data[i].term[j].DO_FIT_NP==YES ||
          parm_data->dihedral_data[i].term[j].DO_FIT_PHASE==YES) {
        for (k=0; k<parm_data->dihedral_data[i].number; ++k) {
          parm_data->fit_atom[ parm_data->dihedral_data[i].atom1[k] ] = YES;
          parm_data->fit_atom[ parm_data->dihedral_data[i].atom2[k] ] = YES;
          parm_data->fit_atom[ parm_data->dihedral_data[i].atom3[k] ] = YES;
          parm_data->fit_atom[ parm_data->dihedral_data[i].atom4[k] ] = YES;
        }
      }
    }
  }
 return SUCCESS;
}

/**
 * Prints out a nice table with forces in amber and quantum for one atom over each structure.
 * It will conduct the forces evaluation itself, including allocation of force structures.
 * @param[in] global_options The global options structure
 * @param[in] parm_data The parameter set to use in the function
 * @param[in] coords_data Pointer to single set of atom coordinates in input structure
 * @param[in] atom Index of the atom to print out. 0 is usually a good choice.
 */
void print_forces(global_options_struct *global_options, parm_struct *parm_data, coord_set *coords_data)
{
  // Create the output file
  FILE *mag, *theta;
  mag=fopen("mag.dat", "w");
  theta=fopen("theta.dat", "w");
  fprintf(theta, "STRUCT\tDIFF\n");
  
  // Create the temporary force structure
  force_struct *single_eval= (force_struct *)malloc(coords_data->natoms*sizeof(force_struct));
  if (single_eval==NULL) {
    printf("*** ERROR allocating memory for temporary forces for structure\n");
    exit(ALLOC_FAIL);
  }
  
  // Run and print to file for each candidate
  int i, atom;
  for (atom=0; atom<coords_data->natoms; ++atom) {
    if (parm_data->fit_atom[atom] == YES) {
      fprintf(mag, "%s ",parm_data->atom[atom].isymbl);
      fprintf(theta, "%s ",parm_data->atom[atom].isymbl);
  } }
  fprintf(mag,"\n");
  fprintf(theta,"\n");
      
      
  for (i=0; i<coords_data->num_coords; ++i) {
      eval_amber_forces_single_struct(global_options, parm_data, coords_data, single_eval, i);
      for (atom=0; atom<coords_data->natoms; ++atom) {
        if (parm_data->fit_atom[atom] == YES) {
        // Calculate the magnitude of the force on our given atom
        double amag = sqrt(single_eval[atom].x*single_eval[atom].x + single_eval[atom].y*single_eval[atom].y + single_eval[atom].z*single_eval[atom].z);
        double qmag = sqrt(coords_data->struc[i].force[atom].x*coords_data->struc[i].force[atom].x + 
                        coords_data->struc[i].force[atom].y*coords_data->struc[i].force[atom].y +
                        coords_data->struc[i].force[atom].z*coords_data->struc[i].force[atom].z );
        fprintf(mag, " %.3f", fabs((qmag-amag)/qmag)*100);
     
        // Calculate the angle between the two vectors
        double xd = ( single_eval[atom].x * coords_data->struc[i].force[atom].x );
        double yd = ( single_eval[atom].y * coords_data->struc[i].force[atom].y );
        double zd = ( single_eval[atom].z * coords_data->struc[i].force[atom].z );
        double dtheta = (xd+yd+zd) / amag; 
        fprintf(theta, " %.3f", dtheta);
      }
    }
    fprintf(mag,"\n");
    fprintf(theta,"\n");
  }
  fclose(mag);
  fclose(theta);
  free(single_eval);
}
