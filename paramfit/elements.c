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
 
/*elements.c*/

/*Contains routines for trying to work out which element an atom is*/

#include <stdio.h>
#include <string.h>
#include "function_def.h"

int find_atomic_number_from_parm(parm_struct *parm_data, int atom)
{
   /*
     This routine tries to work out an elements atomic number from a combination
     of the first letter of it's name and it's mass. Note, this is still unreliable
     but currently the only option I can see.
   */

   /*NOTE, atom numbering starts from 1, but the arrays begin at zero so -1 is subtracted
   from the atom number specified to get the array location*/

   /*First off get the letter from the atom name and then resolve conflicts based on mass*/

   if (!strncmp(parm_data->atom[atom-1].igraph,"A",1))
   {
     /*Could be Al, Ar, As, Ag, At, Ac or Am*/
     if (parm_data->atom[atom-1].amass>20.0 && parm_data->atom[atom-1].amass<30.0)
       return 13; /*Al*/
     else if (parm_data->atom[atom-1].amass>36.0 && parm_data->atom[atom-1].amass<44.0)
       return 18; /*Ar*/
     else if (parm_data->atom[atom-1].amass>70.0 && parm_data->atom[atom-1].amass<78.0)
       return 33; /*As*/
     else if (parm_data->atom[atom-1].amass>100.0 && parm_data->atom[atom-1].amass<115.0)
       return 47; /*Ag*/
     else if (parm_data->atom[atom-1].amass>200.0 && parm_data->atom[atom-1].amass<220.0)
       return 85; /*At*/
     else if (parm_data->atom[atom-1].amass>220.0 && parm_data->atom[atom-1].amass<235.0)
       return 89; /*Ac*/
     else if (parm_data->atom[atom-1].amass>235.0 && parm_data->atom[atom-1].amass<250.0)
       return 95; /*Am*/
     else
       return UNKNOWN_ELEMENT;
   }
   else if (!strncmp(parm_data->atom[atom-1].igraph,"B",1))
   {
     /*Could be Be, B, Br, Ba, Bi, Bk*/
     if (parm_data->atom[atom-1].amass>8.0 && parm_data->atom[atom-1].amass<10.0)
       return 4; /*Be*/
     else if (parm_data->atom[atom-1].amass>10.0 && parm_data->atom[atom-1].amass<12.0)
       return 5; /*B*/
     else if (parm_data->atom[atom-1].amass>75.0 && parm_data->atom[atom-1].amass<84.0)
       return 35; /*Br*/
     else if (parm_data->atom[atom-1].amass>135.0 && parm_data->atom[atom-1].amass<140.0)
       return 56; /*Ba*/
     else if (parm_data->atom[atom-1].amass>200.0 && parm_data->atom[atom-1].amass<215.0)
       return 83; /*Bi*/
     else if (parm_data->atom[atom-1].amass>240.0 && parm_data->atom[atom-1].amass<255.0)
       return 97; /*Bk*/
     else
       return UNKNOWN_ELEMENT;
   }
   else if (!strncmp(parm_data->atom[atom-1].igraph,"C",1))
   {
     /*Could be C, Cl, Ca, Cr, Co, Cu, Cd, Cs, Ce, Cm, Cf*/
     if (parm_data->atom[atom-1].amass>11.0 && parm_data->atom[atom-1].amass<15.0)
       return 6; /*C*/
     else if (parm_data->atom[atom-1].amass>32.0 && parm_data->atom[atom-1].amass<38.0)
       return 17; /*Cl*/
     else if (parm_data->atom[atom-1].amass>38.0 && parm_data->atom[atom-1].amass<44.0)
       return 20; /*Ca*/
     else if (parm_data->atom[atom-1].amass>48.0 && parm_data->atom[atom-1].amass<55.0)
       return 24; /*Cr*/
     else if (parm_data->atom[atom-1].amass>56.0 && parm_data->atom[atom-1].amass<60.0)
       return 27; /*Co*/
     else if (parm_data->atom[atom-1].amass>62.0 && parm_data->atom[atom-1].amass<66.0)
       return 29; /*Cu*/
     else if (parm_data->atom[atom-1].amass>110.0 && parm_data->atom[atom-1].amass<118.0)
       return 48; /*Cd*/
     else if (parm_data->atom[atom-1].amass>128.0 && parm_data->atom[atom-1].amass<136.0)
       return 55; /*Cs*/
     else if (parm_data->atom[atom-1].amass>138.0 && parm_data->atom[atom-1].amass<148.0)
       return 58; /*Ce*/
     else if (parm_data->atom[atom-1].amass>240.0 && parm_data->atom[atom-1].amass<250.0)
       return 96; /*Cm*/
     else if (parm_data->atom[atom-1].amass>250.0 && parm_data->atom[atom-1].amass<256.0)
       return 98; /*Cf*/
     else
       return UNKNOWN_ELEMENT;
   }
   else if (!strncmp(parm_data->atom[atom-1].igraph,"D",1))
   {
     /*Could be Dy*/
     if (parm_data->atom[atom-1].amass>150.0 && parm_data->atom[atom-1].amass<170.0)
       return 66; /*Dy*/
     else
       return UNKNOWN_ELEMENT;
   }
   else if (!strncmp(parm_data->atom[atom-1].igraph,"E",1))
   {
     /*Could be Eu, Er, Es*/
     if (parm_data->atom[atom-1].amass>146.0 && parm_data->atom[atom-1].amass<158.0)
       return 63; /*Eu*/
     else if (parm_data->atom[atom-1].amass>160.0 && parm_data->atom[atom-1].amass<174.0)
       return 68; /*Er*/
     else if (parm_data->atom[atom-1].amass>240.0 && parm_data->atom[atom-1].amass<260.0)
       return 99; /*Es*/
     else
       return UNKNOWN_ELEMENT;
   }
   else if (!strncmp(parm_data->atom[atom-1].igraph,"F",1))
   {
     /*Could be F, Fe, Fr, Fm*/
     if (parm_data->atom[atom-1].amass>16.0 && parm_data->atom[atom-1].amass<20.0)
       return 9; /*F*/
     else if (parm_data->atom[atom-1].amass>50.0 && parm_data->atom[atom-1].amass<60.0)
       return 26; /*Fe*/
     else if (parm_data->atom[atom-1].amass>210.0 && parm_data->atom[atom-1].amass<230.0)
       return 87; /*Fr*/
     else if (parm_data->atom[atom-1].amass>250.0 && parm_data->atom[atom-1].amass<270.0)
       return 100; /*Fm*/
     else
       return UNKNOWN_ELEMENT;
   }
   else if (!strncmp(parm_data->atom[atom-1].igraph,"G",1))
   {
     /*Could be Ga, Ge, Gd */
     if (parm_data->atom[atom-1].amass>67.0 && parm_data->atom[atom-1].amass<71.0)
       return 31; /*Ga*/
     else if (parm_data->atom[atom-1].amass>71.0 && parm_data->atom[atom-1].amass<76.0)
       return 32; /*Ge*/
     else if (parm_data->atom[atom-1].amass>150.0 && parm_data->atom[atom-1].amass<165.0)
       return 64; /*Gd*/
     else
       return UNKNOWN_ELEMENT;
   }
   else if (!strncmp(parm_data->atom[atom-1].igraph,"H",1))
   {
     /*Could be H, He, Hf, Hg, or Ho*/
     if (parm_data->atom[atom-1].amass>0.0 && parm_data->atom[atom-1].amass<3.5)
       return 1; /*No distinguishing between Hydrogen and Deuterium and Tritium*/
     else if (parm_data->atom[atom-1].amass>3.5 && parm_data->atom[atom-1].amass<6.0)
       return 2;
     else if (parm_data->atom[atom-1].amass>150.0 && parm_data->atom[atom-1].amass<170.0)
       return 67;
     else if (parm_data->atom[atom-1].amass>170.0 && parm_data->atom[atom-1].amass<190.0)
       return 72;
     else if (parm_data->atom[atom-1].amass>190.0 && parm_data->atom[atom-1].amass<210.0)
       return 80;
     else
       return UNKNOWN_ELEMENT;
   }
   else if (!strncmp(parm_data->atom[atom-1].igraph,"I",1))
   {
     /*Could be In, I, Ir*/
     if (parm_data->atom[atom-1].amass>110.0 && parm_data->atom[atom-1].amass<118.0)
       return 49; /*In*/
     else if (parm_data->atom[atom-1].amass>122.0 && parm_data->atom[atom-1].amass<130.0)
       return 53; /*I*/
     else if (parm_data->atom[atom-1].amass>186.0 && parm_data->atom[atom-1].amass<196.0)
       return 77; /*Ir*/
     else
       return UNKNOWN_ELEMENT;
   }
   else if (!strncmp(parm_data->atom[atom-1].igraph,"J",1))
   {
     /*Could be NONE*/
     return UNKNOWN_ELEMENT;
   }
   else if (!strncmp(parm_data->atom[atom-1].igraph,"K",1))
   {
     /*Could be K, Kr*/
     if (parm_data->atom[atom-1].amass>36.0 && parm_data->atom[atom-1].amass<44.0)
       return 19; /*K*/
     else if (parm_data->atom[atom-1].amass>80.0 && parm_data->atom[atom-1].amass<88.0)
       return 36; /*Kr*/
     else
       return UNKNOWN_ELEMENT;
   }
   else if (!strncmp(parm_data->atom[atom-1].igraph,"L",1))
   {
     /*Could be Li, La, Lu, Lr*/
     if (parm_data->atom[atom-1].amass>5.0 && parm_data->atom[atom-1].amass<8.0)
       return 3; /*Li*/
     else if (parm_data->atom[atom-1].amass>130.0 && parm_data->atom[atom-1].amass<142.0)
       return 57; /*La*/
     else if (parm_data->atom[atom-1].amass>168.0 && parm_data->atom[atom-1].amass<180.0)
       return 71; /*Lu*/
     else if (parm_data->atom[atom-1].amass>240.0 && parm_data->atom[atom-1].amass<280.0)
       return 103; /*Lr*/
     else
       return UNKNOWN_ELEMENT;
   }
   else if (!strncmp(parm_data->atom[atom-1].igraph,"M",1))
   {
     /*Could be Mg, Mn, Mo, Md*/
     if (parm_data->atom[atom-1].amass>20.0 && parm_data->atom[atom-1].amass<30.0)
       return 12; /*Mg*/
     else if (parm_data->atom[atom-1].amass>50.0 && parm_data->atom[atom-1].amass<58.0)
       return 25; /*Mn*/
     else if (parm_data->atom[atom-1].amass>90.0 && parm_data->atom[atom-1].amass<100.0)
       return 42; /*Mo*/
     else if (parm_data->atom[atom-1].amass>240.0 && parm_data->atom[atom-1].amass<260.0)
       return 101; /*Md*/
     else
       return UNKNOWN_ELEMENT;
   }
   else if (!strncmp(parm_data->atom[atom-1].igraph,"N",1))
   {
     /*Could be N, Ne, Na, Ni, Nb, Nd, Np, No*/
     if (parm_data->atom[atom-1].amass>12.0 && parm_data->atom[atom-1].amass<18.0)
       return 7; /*N*/
     else if (parm_data->atom[atom-1].amass>18.0 && parm_data->atom[atom-1].amass<21.0)
       return 10; /*Ne*/
     else if (parm_data->atom[atom-1].amass>21.0 && parm_data->atom[atom-1].amass<26.0)
       return 11; /*Na*/
     else if (parm_data->atom[atom-1].amass>54.0 && parm_data->atom[atom-1].amass<62.0)
       return 28; /*Ni*/
     else if (parm_data->atom[atom-1].amass>86.0 && parm_data->atom[atom-1].amass<98.0)
       return 41; /*Nb*/
     else if (parm_data->atom[atom-1].amass>140.0 && parm_data->atom[atom-1].amass<150.0)
       return 60; /*Nd*/
     else if (parm_data->atom[atom-1].amass>220.0 && parm_data->atom[atom-1].amass<248.0)
       return 93; /*Np*/
     else if (parm_data->atom[atom-1].amass>248.0 && parm_data->atom[atom-1].amass<260.0)
       return 102; /*No*/
     else
       return UNKNOWN_ELEMENT;
   }
   else if (!strncmp(parm_data->atom[atom-1].igraph,"O",1))
   {
     /*Could be O, Os*/
     if (parm_data->atom[atom-1].amass>12.0 && parm_data->atom[atom-1].amass<18.0)
       return 8; /*O*/
     else if (parm_data->atom[atom-1].amass>180.0 && parm_data->atom[atom-1].amass<200.0)
       return 76; /*Os*/
     else
       return UNKNOWN_ELEMENT;
   }
   else if (!strncmp(parm_data->atom[atom-1].igraph,"P",1))
   {
     /*Could be P, Pd, Pr, Pm, Pt, Po, Pa, Pu*/
     if (parm_data->atom[atom-1].amass>26.0 && parm_data->atom[atom-1].amass<34.0)
       return 15; /*P*/
     else if (parm_data->atom[atom-1].amass>100.0 && parm_data->atom[atom-1].amass<110.0)
       return 46; /*Pd*/
     else if (parm_data->atom[atom-1].amass>136.0 && parm_data->atom[atom-1].amass<142.0)
       return 59; /*Pr*/
     else if (parm_data->atom[atom-1].amass>142.0 && parm_data->atom[atom-1].amass<150.0)
       return 61; /*Pm*/
     else if (parm_data->atom[atom-1].amass>190.0 && parm_data->atom[atom-1].amass<200.0)
       return 78; /*Pt*/
     else if (parm_data->atom[atom-1].amass>205.0 && parm_data->atom[atom-1].amass<215.0)
       return 84; /*Po*/
     else if (parm_data->atom[atom-1].amass>225.0 && parm_data->atom[atom-1].amass<238.0)
       return 91; /*Pa*/
     else if (parm_data->atom[atom-1].amass>238.0 && parm_data->atom[atom-1].amass<254.0)
       return 94; /*Pu*/
     else
       return UNKNOWN_ELEMENT;
   }
   else if (!strncmp(parm_data->atom[atom-1].igraph,"Q",1))
   {
     /*Could be NONE*/
       return UNKNOWN_ELEMENT;
   }
   else if (!strncmp(parm_data->atom[atom-1].igraph,"R",1))
   {
     /*Could be Rb, Ru, Rh, Re, Rn, Ra*/
     if (parm_data->atom[atom-1].amass>80.0 && parm_data->atom[atom-1].amass<90.0)
       return 37; /*Rb*/
     else if (parm_data->atom[atom-1].amass>96.0 && parm_data->atom[atom-1].amass<102.0)
       return 44; /*Ru*/
     else if (parm_data->atom[atom-1].amass>102.0 && parm_data->atom[atom-1].amass<106.0)
       return 45; /*Rh*/
     else if (parm_data->atom[atom-1].amass>180.0 && parm_data->atom[atom-1].amass<190.0)
       return 75; /*Re*/
     else if (parm_data->atom[atom-1].amass>210.0 && parm_data->atom[atom-1].amass<225.0)
       return 86; /*Rn*/
     else if (parm_data->atom[atom-1].amass>225.0 && parm_data->atom[atom-1].amass<230.0)
       return 88; /*Ra*/
     else
       return UNKNOWN_ELEMENT;
   }
   else if (!strncmp(parm_data->atom[atom-1].igraph,"S",1))
   {
     /*Could be Si, S, Sc, Se, Sr, Sn, Sb, Sm*/
     if (parm_data->atom[atom-1].amass>24.0 && parm_data->atom[atom-1].amass<30.0)
       return 14; /*Si*/
     else if (parm_data->atom[atom-1].amass>30.0 && parm_data->atom[atom-1].amass<36.0)
       return 16; /*S*/
     else if (parm_data->atom[atom-1].amass>40.0 && parm_data->atom[atom-1].amass<48.0)
       return 21; /*Sc*/
     else if (parm_data->atom[atom-1].amass>70.0 && parm_data->atom[atom-1].amass<84.0)
       return 34; /*Se*/
     else if (parm_data->atom[atom-1].amass>84.0 && parm_data->atom[atom-1].amass<91.0)
       return 38; /*Sr*/
     else if (parm_data->atom[atom-1].amass>112.0 && parm_data->atom[atom-1].amass<120.0)
       return 50; /*Sn*/
     else if (parm_data->atom[atom-1].amass>120.0 && parm_data->atom[atom-1].amass<126.0)
       return 51; /*Sb*/
     else if (parm_data->atom[atom-1].amass>140.0 && parm_data->atom[atom-1].amass<160.0)
       return 62; /*Sm*/
     else
       return UNKNOWN_ELEMENT;
   }
   else if (!strncmp(parm_data->atom[atom-1].igraph,"T",1))
   {
     /*Could be Ti, Tc, Te, Tb, Tm, Ta, Tl Th, */
     if (parm_data->atom[atom-1].amass>46.0 && parm_data->atom[atom-1].amass<50.0)
       return 22; /*Ti*/
     else if (parm_data->atom[atom-1].amass>94.0 && parm_data->atom[atom-1].amass<102.0)
       return 43; /*Tc*/
     else if (parm_data->atom[atom-1].amass>120.0 && parm_data->atom[atom-1].amass<134.0)
       return 52; /*Te*/
     else if (parm_data->atom[atom-1].amass>150.0 && parm_data->atom[atom-1].amass<166.0)
       return 65; /*Tb*/
     else if (parm_data->atom[atom-1].amass>164.0 && parm_data->atom[atom-1].amass<174.0)
       return 69; /*Tm*/
     else if (parm_data->atom[atom-1].amass>176.0 && parm_data->atom[atom-1].amass<185.0)
       return 73; /*Ta*/
     else if (parm_data->atom[atom-1].amass>200.0 && parm_data->atom[atom-1].amass<208.0)
       return 81; /*Tl*/
     else if (parm_data->atom[atom-1].amass>224.0 && parm_data->atom[atom-1].amass<240.0)
       return 90; /*Th*/
     else
       return UNKNOWN_ELEMENT;
   }
   else if (!strncmp(parm_data->atom[atom-1].igraph,"U",1))
   {
     /*Could be U*/
     if (parm_data->atom[atom-1].amass>230.0 && parm_data->atom[atom-1].amass<250.0)
       return 92; /*U*/
     else
       return UNKNOWN_ELEMENT;
   }
   else if (!strncmp(parm_data->atom[atom-1].igraph,"V",1))
   {
     /*Could be V*/
     if (parm_data->atom[atom-1].amass>46.0 && parm_data->atom[atom-1].amass<54.0)
       return 23; /*V*/
     else
       return UNKNOWN_ELEMENT;
   }
   else if (!strncmp(parm_data->atom[atom-1].igraph,"W",1))
   {
     /*Could be W*/
     if (parm_data->atom[atom-1].amass>176.0 && parm_data->atom[atom-1].amass<186.0)
       return 74; /*W*/
     else
       return UNKNOWN_ELEMENT;
   }
   else if (!strncmp(parm_data->atom[atom-1].igraph,"X",1))
   {
     /*Could be Xe */
     if (parm_data->atom[atom-1].amass>126.0 && parm_data->atom[atom-1].amass<136.0)
       return 54; /*Xe*/
     else
       return UNKNOWN_ELEMENT;
   }
   else if (!strncmp(parm_data->atom[atom-1].igraph,"Y",1))
   {
     /*Could be Y, Yb */
     if (parm_data->atom[atom-1].amass>84.0 && parm_data->atom[atom-1].amass<92.0)
       return 39; /*Y*/
     else if (parm_data->atom[atom-1].amass>168.0 && parm_data->atom[atom-1].amass<180.0)
       return 70; /*Yb*/
     else
       return UNKNOWN_ELEMENT;
   }
   else if (!strncmp(parm_data->atom[atom-1].igraph,"Z",1))
   {
     /*Could be Zn, Zr */
     if (parm_data->atom[atom-1].amass>60.0 && parm_data->atom[atom-1].amass<70.0)
       return 30; /*Zn*/
     else if (parm_data->atom[atom-1].amass>88.0 && parm_data->atom[atom-1].amass<94.0)
       return 40; /*Zr*/
     else
       return UNKNOWN_ELEMENT;
   }

   /*If we get here it must be an unknown element*/
   return UNKNOWN_ELEMENT;

}

void print_atomic_number_as_symbol(FILE *fptr, int atomic_number)
{
  /*Prints the atomic symbol represented by the atomic_number to the file pointed to by fptr*/

  switch (atomic_number)
  {
     case 1:  fprintf(fptr,"H");  break;
     case 2:  fprintf(fptr,"He"); break;
     case 3:  fprintf(fptr,"Li"); break;
     case 4:  fprintf(fptr,"Be"); break;
     case 5:  fprintf(fptr,"B");  break;
     case 6:  fprintf(fptr,"C");  break;
     case 7:  fprintf(fptr,"N");  break;
     case 8:  fprintf(fptr,"O");  break;
     case 9:  fprintf(fptr,"F");  break;
     case 10: fprintf(fptr,"Ne"); break;
     case 11: fprintf(fptr,"Na"); break;
     case 12: fprintf(fptr,"Mg"); break;
     case 13: fprintf(fptr,"Al"); break;
     case 14: fprintf(fptr,"Si"); break;
     case 15: fprintf(fptr,"P");  break;
     case 16: fprintf(fptr,"S");  break;
     case 17: fprintf(fptr,"Cl");  break;
     case 18: fprintf(fptr,"Ar");  break;          
     case 19: fprintf(fptr,"K");  break;
     case 20: fprintf(fptr,"Ca");  break;
     case 21: fprintf(fptr,"Sc");  break;
     case 22: fprintf(fptr,"Ti");  break;
     case 23: fprintf(fptr,"V");  break;
     case 24: fprintf(fptr,"Cr");  break;
     case 25: fprintf(fptr,"Mn");  break;
     case 26: fprintf(fptr,"Fe");  break;
     case 27: fprintf(fptr,"Co");  break;
     case 28: fprintf(fptr,"Ni");  break;
     case 29: fprintf(fptr,"Cu");  break;
     case 30: fprintf(fptr,"Zn");  break;
     case 31: fprintf(fptr,"Ga");  break;
     case 32: fprintf(fptr,"Ge");  break;
     case 33: fprintf(fptr,"As");  break;
     case 34: fprintf(fptr,"Se");  break;
     case 35: fprintf(fptr,"Br");  break;
     case 36: fprintf(fptr,"Kr");  break;
     case 37: fprintf(fptr,"Rb");  break;
     case 38: fprintf(fptr,"Sr");  break;
     case 39: fprintf(fptr,"Y");   break;
     case 40: fprintf(fptr,"Zr");  break;
     case 41: fprintf(fptr,"Nb");  break;
     case 42: fprintf(fptr,"Mo");  break;
     case 43: fprintf(fptr,"Tc");  break;
     case 44: fprintf(fptr,"Ru");  break;
     case 45: fprintf(fptr,"Rh");  break;
     case 46: fprintf(fptr,"Pd");  break;
     case 47: fprintf(fptr,"Ag");  break;
     case 48: fprintf(fptr,"Cd");  break;
     case 49: fprintf(fptr,"In");  break;
     case 50: fprintf(fptr,"Sn");  break;
     case 51: fprintf(fptr,"Sb");  break;
     case 52: fprintf(fptr,"Te");  break;
     case 53: fprintf(fptr,"I");   break;
     case 54: fprintf(fptr,"Xe");  break;
     case 55: fprintf(fptr,"Cs");  break;
     case 56: fprintf(fptr,"Ba");  break;
     case 57: fprintf(fptr,"La");  break;
     case 58: fprintf(fptr,"Ce");  break;
     case 59: fprintf(fptr,"Pr");  break;
     case 60: fprintf(fptr,"Nd");  break;
     case 61: fprintf(fptr,"Pm");  break;
     case 62: fprintf(fptr,"Sm");  break;
     case 63: fprintf(fptr,"Eu");  break;
     case 64: fprintf(fptr,"Gd");  break;
     case 65: fprintf(fptr,"Tb");  break;
     case 66: fprintf(fptr,"Dy");  break;
     case 67: fprintf(fptr,"Ho");  break;
     case 68: fprintf(fptr,"Er");  break;
     case 69: fprintf(fptr,"Tm");  break;
     case 70: fprintf(fptr,"Yb");  break;
     case 71: fprintf(fptr,"Lu");  break;
     case 72: fprintf(fptr,"Hf");  break;
     case 73: fprintf(fptr,"Ta");  break;
     case 74: fprintf(fptr,"W");   break;
     case 75: fprintf(fptr,"Re");  break;
     case 76: fprintf(fptr,"Os");  break;
     case 77: fprintf(fptr,"Ir");  break;
     case 78: fprintf(fptr,"Pt");  break;
     case 79: fprintf(fptr,"Au");  break;
     case 80: fprintf(fptr,"Hg");  break;
     case 81: fprintf(fptr,"Tl");  break;
     case 82: fprintf(fptr,"Pb");  break;
     case 83: fprintf(fptr,"Bi");  break;
     case 84: fprintf(fptr,"Po");  break;
     case 85: fprintf(fptr,"At");  break;
     case 86: fprintf(fptr,"Rn");  break;
     case 87: fprintf(fptr,"Fr");  break;
     case 88: fprintf(fptr,"Ra");  break;
     case 89: fprintf(fptr,"Ac");  break;
     case 90: fprintf(fptr,"Th");  break;
     case 91: fprintf(fptr,"Pa");  break;
     case 92: fprintf(fptr,"U");   break;
     case 93: fprintf(fptr,"Np");  break;
     case 94: fprintf(fptr,"Pu");  break;
     case 95: fprintf(fptr,"Am");  break;
     case 96: fprintf(fptr,"Cm");  break;
     case 97: fprintf(fptr,"Bk");  break;
     case 98: fprintf(fptr,"Cf");  break;
     case 99: fprintf(fptr,"Es");  break;
    case 100: fprintf(fptr,"Md");  break;
    case 101: fprintf(fptr,"No");  break;
    case 102: fprintf(fptr,"Lr");  break;
    case 103: fprintf(fptr,"Rb");  break;
    default:
       printf("\n!  WARNING: UNKNOWN ATOMIC NUMBER IN PRINT_ATOMIC_NUMBER_AS_SYMBOL - %d\n",atomic_number);
       fprintf(fptr,"XXX");
       break;       
  }

  

}


