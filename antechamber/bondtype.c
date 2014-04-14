/*
************************************************************************
*           All Copyright Reserved!                                    *
*                                                                      *
*  Prog:    bondtype                                                   *
*  Version: version 1.0                                                *
*  Author:  Junmei Wang                                                *
*                                                                      *
*  Department of Pharmaceutical Chemistry                              *
*  School of Pharmacy                                                  *
*  University of California                                            *
*  San Francisco   CA 94143                                            *
*  Octomber, 2001                                                      *
************************************************************************
*/

#include <assert.h>

#include "ncsu-penalties.h" /* before define.h */

char *amberhome;
# include "common.h"
# include "define.h"
# include "atom.h"
# include "utility.c"
# include "common.c"
# include "ring.c"
# include "ac.c"
# include "mol2.c"
# define debug 0

static ATOM *atom;
static AROM *arom;
static BOND *bond;
static RING *ring;
static int ringnum = 0;

static MOLINFO minfo;
static CONTROLINFO cinfo;

static char ifilename[MAXCHAR];
static char ofilename[MAXCHAR];
static int atomnum = 0;
static int bondnum = 0;
static int j, k;
static int iformat = 1;
static AV *av;
static AV avH;
static int *valence;
static int *va_best;
static int *conjatom;
static int maxaps = 0;
static int judge_flag = 0;

typedef struct {
	int resno;
	char name[10];
} RES;

static int va_allocated = 0;

static int vastatenum = 0;

static int *va_num = NULL;
static int *va_penalty = NULL;
static int **va_valence = NULL;
static int **va_addr = NULL;

static int *con_num;
static int *ind_bondnum;
static int *ind_valence;
static int *bondindex;
static int *org_bondtype;
static int *bondtype_bak;
static int *bondindex_bak;
static int *apstype; 
static int *seq;

static int resnum;
static ATOM *atom2;
static BOND *bond2;
static AV *av2;
static int *valence2;
static int *va_best2;
static int *apstype2; 
static int atomnum2;
static int bondnum2;
static int maxatomnum2;
static int maxbondnum2;
/* 
i  iaps = 0, frozen unknown aps
   iaps = 1, stop running
*/
static int iaps = 0; 
static RES *res;

static void assignav()
{
	FILE *fp;
	int i;
	int tmpint;
	int index;
	int index1, index2, index3, index4, index5, index6, index7, index8, index9, index10;
	int con1, con2, con3;
	char line[MAXCHAR];
	char filename[MAXCHAR];
	char tmpc1[5], tmpc2[5], tmpc3[5], tmpc4[5], tmpc5[5];
	char tmpc6[5], tmpc7[5], tmpc8[5], tmpc9[5], tmpc10[5];

	build_dat_path(filename, "APS.DAT", sizeof filename, 0);

	if ((fp = fopen(filename, "r")) == NULL) {
		fprintf(stdout, "Cannot open the APS.DAT file: %s, exit\n", filename);
		exit(1);
	}
	for (i = 0; i < atomnum; i++) {
		index1 = 0;
		index2 = 0;
		index3 = 0;
		index4 = 0;
		index5 = 0;
		index6 = 0;
		index7 = 0;
		index8 = 0;
		index9 = 0;
		index10 = 0;
		av[i].aps[0] = 9999;
		av[i].aps[1] = 9999;
		av[i].aps[2] = 9999;
		av[i].aps[3] = 9999;
		av[i].aps[4] = 9999;
		av[i].aps[5] = 9999;
		av[i].aps[6] = 9999;
		av[i].aps[7] = 9999;
		if (atom[i].atomicnum == 6 && atom[i].connum == 3) {
			if (atom[atom[i].con[0]].atomicnum == 8
				&& atom[atom[i].con[0]].connum == 1)
				index1++;
			if (atom[atom[i].con[1]].atomicnum == 8
				&& atom[atom[i].con[1]].connum == 1)
				index1++;
			if (atom[atom[i].con[2]].atomicnum == 8
				&& atom[atom[i].con[2]].connum == 1)
				index1++;
			if (atom[atom[i].con[0]].atomicnum == 16
				&& atom[atom[i].con[0]].connum == 1)
				index1++;
			if (atom[atom[i].con[1]].atomicnum == 16
				&& atom[atom[i].con[1]].connum == 1)
				index1++;
			if (atom[atom[i].con[2]].atomicnum == 16
				&& atom[atom[i].con[2]].connum == 1)
				index1++;
		}
		if (atom[i].atomicnum == 15 && atom[i].connum == 4) {
			if (atom[atom[i].con[0]].atomicnum == 8
				&& atom[atom[i].con[0]].connum == 1)
				index2++;
			if (atom[atom[i].con[1]].atomicnum == 8
				&& atom[atom[i].con[1]].connum == 1)
				index2++;
			if (atom[atom[i].con[2]].atomicnum == 8
				&& atom[atom[i].con[2]].connum == 1)
				index2++;
			if (atom[atom[i].con[3]].atomicnum == 8
				&& atom[atom[i].con[3]].connum == 1)
				index2++;
			if (atom[atom[i].con[0]].atomicnum == 16
				&& atom[atom[i].con[0]].connum == 1)
				index2++;
			if (atom[atom[i].con[1]].atomicnum == 16
				&& atom[atom[i].con[1]].connum == 1)
				index2++;
			if (atom[atom[i].con[2]].atomicnum == 16
				&& atom[atom[i].con[2]].connum == 1)
				index2++;
			if (atom[atom[i].con[3]].atomicnum == 16
				&& atom[atom[i].con[3]].connum == 1)
				index2++;
		}
		if (atom[i].atomicnum == 16 && atom[i].connum == 4) {
			if (atom[atom[i].con[0]].atomicnum == 8
				&& atom[atom[i].con[0]].connum == 1)
				index3++;
			if (atom[atom[i].con[1]].atomicnum == 8
				&& atom[atom[i].con[1]].connum == 1)
				index3++;
			if (atom[atom[i].con[2]].atomicnum == 8
				&& atom[atom[i].con[2]].connum == 1)
				index3++;
			if (atom[atom[i].con[3]].atomicnum == 8
				&& atom[atom[i].con[3]].connum == 1)
				index3++;
			if (atom[atom[i].con[0]].atomicnum == 16
				&& atom[atom[i].con[0]].connum == 1)
				index3++;
			if (atom[atom[i].con[1]].atomicnum == 16
				&& atom[atom[i].con[1]].connum == 1)
				index3++;
			if (atom[atom[i].con[2]].atomicnum == 16
				&& atom[atom[i].con[2]].connum == 1)
				index3++;
			if (atom[atom[i].con[3]].atomicnum == 16
				&& atom[atom[i].con[3]].connum == 1)
				index3++;
		}
/* for N1- */
		if (atom[i].atomicnum == 7 && atom[i].connum == 1) {
			if (atom[atom[i].con[0]].atomicnum == 7 && 
				atom[atom[i].con[0]].connum == 2) 
				index4 = 1;
		}
/* for N2+ */
		if (atom[i].atomicnum == 7 && atom[i].connum == 2) {
			if ((atom[atom[i].con[0]].atomicnum == 7 || atom[atom[i].con[0]].atomicnum == 6)&& 
				atom[atom[i].con[0]].connum == 1) 
				index5 = 1 ;
			if ((atom[atom[i].con[1]].atomicnum == 7 || atom[atom[i].con[1]].atomicnum == 6)&& 
				atom[atom[i].con[1]].connum == 1) 
				index5 = 1 ;
		}
/* for NO2 */
		if (atom[i].atomicnum == 7 && atom[i].connum == 3) {
			if (atom[atom[i].con[0]].atomicnum == 8
				&& atom[atom[i].con[0]].connum == 1)
				index6++;
			if (atom[atom[i].con[1]].atomicnum == 8
				&& atom[atom[i].con[1]].connum == 1)
				index6++;
			if (atom[atom[i].con[2]].atomicnum == 8
				&& atom[atom[i].con[2]].connum == 1)
				index6++;
			if (atom[atom[i].con[0]].atomicnum == 16
				&& atom[atom[i].con[0]].connum == 1)
				index6++;
			if (atom[atom[i].con[1]].atomicnum == 16
				&& atom[atom[i].con[1]].connum == 1)
				index6++;
			if (atom[atom[i].con[2]].atomicnum == 16
				&& atom[atom[i].con[2]].connum == 1)
				index6++;
		}
/* for N2+ */
		if (atom[i].atomicnum == 7 && atom[i].connum == 3 && index6 < 2) {
			if (atom[atom[i].con[0]].atomicnum == 8 && 
				atom[atom[i].con[0]].connum == 1) 
				index7=1 ;
			if (atom[atom[i].con[1]].atomicnum == 8 && 
				atom[atom[i].con[1]].connum == 1) 
				index7=1 ;
			if (atom[atom[i].con[2]].atomicnum == 8 && 
				atom[atom[i].con[2]].connum == 1) 
				index7=1 ;
			if (atom[atom[i].con[0]].atomicnum == 16 && 
				atom[atom[i].con[0]].connum == 1) 
				index7=1 ;
			if (atom[atom[i].con[1]].atomicnum == 16 && 
				atom[atom[i].con[1]].connum == 1) 
				index7=1 ;
			if (atom[atom[i].con[2]].atomicnum == 16 && 
				atom[atom[i].con[2]].connum == 1) 
				index7=1 ;
		}
/* for O1- */
		if (atom[i].atomicnum == 8 && atom[i].connum == 1)  
			if (atom[atom[i].con[0]].atomicnum == 7 && 
				atom[atom[i].con[0]].connum == 3)  {
				con1 = atom[atom[i].con[0]].con[0];
				con2 = atom[atom[i].con[0]].con[1];
				con3 = atom[atom[i].con[0]].con[2];
				if (atom[con1].atomicnum == 8 && atom[con1].connum == 1) index8++;	
				if (atom[con2].atomicnum == 8 && atom[con2].connum == 1) index8++;	
				if (atom[con3].atomicnum == 8 && atom[con3].connum == 1) index8++;	
				if (atom[con1].atomicnum == 16 && atom[con1].connum == 1) index8++;	
				if (atom[con2].atomicnum == 16 && atom[con2].connum == 1) index8++;	
				if (atom[con3].atomicnum == 16 && atom[con3].connum == 1) index8++;	
				if(index8 <= 1)
					index8 =1;
				else
					index8 = 0;
			}
/* for S1- */
		if (atom[i].atomicnum == 16 && atom[i].connum == 1)  
			if (atom[atom[i].con[0]].atomicnum == 7 && 
				atom[atom[i].con[0]].connum == 3)  {
				con1 = atom[atom[i].con[0]].con[0];
				con2 = atom[atom[i].con[0]].con[1];
				con3 = atom[atom[i].con[0]].con[2];
				if (atom[con1].atomicnum == 8 && atom[con1].connum == 1) index9++;	
				if (atom[con2].atomicnum == 8 && atom[con2].connum == 1) index9++;	
				if (atom[con3].atomicnum == 8 && atom[con3].connum == 1) index9++;	
				if (atom[con1].atomicnum == 16 && atom[con1].connum == 1) index9++;	
				if (atom[con2].atomicnum == 16 && atom[con2].connum == 1) index9++;	
				if (atom[con3].atomicnum == 16 && atom[con3].connum == 1) index9++;	
				if(index9 <= 1)
					index9 =1;
				else
					index9 = 0;
			}
/* for C#N- */
		if (atom[i].atomicnum == 6 && atom[i].connum == 1) 
			if (atom[atom[i].con[0]].atomicnum == 7
				&& atom[atom[i].con[0]].connum == 2)
				index10 = 1;

		conjatom[i] = 0;
		if (index1 >= 2 || index2 >= 2)
			conjatom[i] = 1;
/*		if (index3 >= 2 || index4 >= 3)   */
/*	        index3 >=3 eliminates the possibility of S=O bonds in SO2 to be delocalized bonds  
                index4 cannot be larger than 1 according to the code above, it is a bug !
		index6 >=2 makes NO2 bonds are delocalized bonds
*/

   		if (index3 >= 3 || index6 >= 2)
			conjatom[i] = 1;
		for (;;) {
			if (fgets(line, MAXCHAR, fp) == NULL)
				break;
			if ((strncmp("APS", line, 3) == 0 && index1 < 2 && index2 < 2 && index3 < 2 && index4 < 1 && index5 < 1 
							  && index6 < 2 && index7 < 1 && index8 <1 && index9 < 1 && index10 <1)
				|| (strncmp("APSCO2", line, 6) == 0 && index1 >= 2)
				|| (strncmp("APSPO2", line, 6) == 0 && index2 == 2)
				|| (strncmp("APSPO3", line, 6) == 0 && index2 > 2)
                                || (strncmp("APSSO2", line, 6) == 0 && index3 == 2)   
                                || (strncmp("APSSO3", line, 6) == 0 && index3 == 3) 
				|| (strncmp("APSSO4", line, 6) == 0 && index3 == 4)   
				|| (strncmp("APSN1-", line, 6) == 0 && index4 == 1)
				|| (strncmp("APSN2+", line, 6) == 0 && index5 == 1)
				|| (strncmp("APSNO2", line, 6) == 0 && index6 >= 2)
				|| (strncmp("APSN3+", line, 6) == 0 && index7 == 1)
				|| (strncmp("APSO1-", line, 6) == 0 && index8 == 1)   
				|| (strncmp("APSS1-", line, 6) == 0 && index9 == 1)   
				|| (strncmp("APSC1+", line, 6) == 0 && index10 == 1)) {
				sscanf(&line[7], "%s%d%s%s%s%s%s%s%s%s%s", tmpc1, &tmpint,
					   tmpc2, tmpc3, tmpc4, tmpc5, tmpc6, tmpc7, tmpc8,
					   tmpc9, tmpc10);
				if (tmpint == atom[i].atomicnum
					&& (tmpc2[0] == '*'
						|| atoi(tmpc2) == atom[i].connum)) {
					if (tmpc3[0] == '*')
						av[i].aps[0] = 9999;
					else
						av[i].aps[0] = atoi(tmpc3);

					if (tmpc4[0] == '*')
						av[i].aps[1] = 9999;
					else
						av[i].aps[1] = atoi(tmpc4);

					if (tmpc5[0] == '*')
						av[i].aps[2] = 9999;
					else
						av[i].aps[2] = atoi(tmpc5);

					if (tmpc6[0] == '*')
						av[i].aps[3] = 9999;
					else
						av[i].aps[3] = atoi(tmpc6);

					if (tmpc7[0] == '*')
						av[i].aps[4] = 9999;
					else
						av[i].aps[4] = atoi(tmpc7);

					if (tmpc8[0] == '*')
						av[i].aps[5] = 9999;
					else
						av[i].aps[5] = atoi(tmpc8);

					if (tmpc9[0] == '*')
						av[i].aps[6] = 9999;
					else
						av[i].aps[6] = atoi(tmpc9);
					if (tmpc10[0] == '*')
						av[i].aps[7] = 9999;
					else
						av[i].aps[7] = atoi(tmpc10);
					break;
				}
			}
		}
		rewind(fp);
	}
	fclose(fp);

	if (debug) {
		printf("\nList of atomic valence information\n");
		for (i = 0; i < atomnum; i++)
			printf("\n%4d%5s%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d", i + 1,
				   atom[i].name, atom[i].atomicnum, atom[i].connum,
				   av[i].aps[0], av[i].aps[1], av[i].aps[2], av[i].aps[3],
				   av[i].aps[4], av[i].aps[5], av[i].aps[6], av[i].aps[7]);
	}
	maxaps = 0;
	for (i = 0; i < atomnum; i++) {
		index = 0;
		for (j = 0; j < 8; j++) {
			if (av[i].aps[j] == 0) {
				index = 1;
				apstype[i] = 1;
				valence[i] = j;
				va_best[i] = j;
			}
			if (av[i].aps[j] > PSCUTOFF)
				continue;
			if (maxaps < av[i].aps[j])
				maxaps = av[i].aps[j];
		}
		if (index == 0) {
			if(iaps == 0)
				printf("\nFor atom[%d]:%s, the best APS is not zero, bonds involved by this atom are frozen\n", 
				   i + 1, atom[i].name);
			else if(iaps == 1) {
				printf("\nFor atom[%d]:%s, the best APS is not zero, exit\n",
				   i + 1, atom[i].name);
				exit(1);
			}
		}
	}
	if (debug)
		printf("\n\nThe Maximum APS is %d (PSCUTOFF = %d)\n", maxaps, PSCUTOFF);

}

static void memory(int atomsize, int bondsize) {
atom2 = (ATOM *) malloc(sizeof(ATOM) * atomsize);
if (atom2 == NULL) {       
        fprintf(stdout, "memory allocation error for *atom2\n");
        exit(1);
}

bond2 = (BOND *) malloc(sizeof(BOND) * bondsize);
if (bond2 == NULL) { 
        fprintf(stdout, "memory allocation error for *bond2\n");
        exit(1);
}   
int i;
for (i = 0; i < bondsize; ++i) {
        bond2[i].jflag = -1; /* bond type has not been assigned */
}

av2 = (AV *) malloc(sizeof(AV) * atomsize);
if (av2 == NULL) {
        fprintf(stdout, "memory allocation error for *av2\n");
        exit(1);
}

valence2 = (int *) malloc(sizeof(int) * atomsize);
if (valence2 == NULL) {
        fprintf(stdout, "memory allocation error for *valence2\n");
        exit(1);
}

va_best2 = (int *) malloc(sizeof(int) * atomsize);
if (va_best2 == NULL) {
        fprintf(stdout, "memory allocation error for *va_best2\n");
        exit(1);
}

apstype2 = (int *) malloc(sizeof(int) * atomsize);
if (apstype2 == NULL) {
        fprintf(stdout, "memory allocation error for *apstype2\n");
        exit(1);
}

con_num = (int *) malloc(sizeof(int) * atomsize);
if (con_num == NULL) {
        fprintf(stdout, "memory allocation error for *con_num\n");
        exit(1);
}
ind_bondnum = (int *) malloc(sizeof(int) * atomsize);
if (ind_bondnum == NULL) {
        fprintf(stdout, "memory allocation error for *ind_bondnum\n");
        exit(1);
}
ind_valence = (int *) malloc(sizeof(int) * atomsize);
if (ind_valence == NULL) {
        fprintf(stdout, "memory allocation error for *ind_valence\n");
        exit(1);
}
bondindex = (int *) malloc(sizeof(int) * bondsize);
if (bondindex == NULL) {
        fprintf(stdout, "memory allocation error for *bondindex\n");
        exit(1);
}
bondindex_bak = (int *) malloc(sizeof(int) * bondsize);
if (bondindex_bak == NULL) {
        fprintf(stdout, "memory allocation error for *bondindex_bak\n");
        exit(1);
}
bondtype_bak = (int *) malloc(sizeof(int) * bondsize);
if (bondtype_bak == NULL) {
        fprintf(stdout, "memory allocation error for *bondtype_bak\n");
        exit(1);
}
}

static void current_va(int id)
{
	int i;
	if (debug) {
		printf
			("\n\nThe current valence state (%d) is listed as the following",
			 id);
		for (i = 0; i < atomnum2; i++)
			printf("\n%5d%5s%5d%5d", i + 1, atom2[i].name, valence2[i],
				   av2[i].aps[valence2[i]]);
		printf("\n\n");
	}
}

static int jbo_score()
{
	int i;
	int violation;

	violation = 0;

	for (i = 0; i < atomnum2; i++) {
		if(apstype2[i] == -1) continue;
		if (ind_bondnum[i] == 0 && ind_valence[i] != 0)
			violation++;
		if (ind_bondnum[i] > ind_valence[i])
			violation++;
	}
/*we have determined all of the bond types involved in this atom, 
however, the ind_valence[i] is not zero */
	return violation;
}


static void jbo_induce(void)
{
/*find the induced bond types, called by jbo_iteration*/
	int index;
	int i, j, k;
	int tmpint;
	int at1, at2;
	int num;
	index = 1;
	while (index > 0) {
		index = 0;
		num = 0;
		for (i = 0; i < atomnum2; i++)
			if ((ind_bondnum[i] == 1 && ind_valence[i] > 0)
				|| ((ind_bondnum[i] == ind_valence[i]) &&
					ind_valence[i] != 0))
				for (j = 0; j < atom2[i].connum; j++) {
					at1 = i;
					at2 = atom2[i].con[j];
					for (k = 0; k < bondnum2; k++) {
						if (bondindex[k] == 1)
							continue;
						if ((bond2[k].bondi == at1 && bond2[k].bondj == at2)
							|| (bond2[k].bondi == at2
								&& bond2[k].bondj == at1)) {
							tmpint = ind_valence[i] / ind_bondnum[i];
							bond2[k].type = tmpint;
							bondindex[k] = 1;
							ind_valence[at1] -= tmpint;
							ind_valence[at2] -= tmpint;
							ind_bondnum[at1]--;
							ind_bondnum[at2]--;
							num++;
							break;
						}
					}
				}
		if (num != 0)
			index = 1;			/*prevent the program falling into dead cycle */
	}
}


static int jbo_iteration()
{
	int at1, at2;
	int i, k;
	int flag;
	int while_flag = 1;
	int score = 0;
/*
for a certain valence state, we need to find one possible arrangement of bond types
for any undetermined bond, it can be single, double or triple bond, we first assume it
is a single bond, then double then triple
*/
	while (while_flag == 1) {
/*first assume it is a single bond */

		for (k = 0; k < bondnum2; k++) {
			bondtype_bak[k] = bond2[k].type;
			bondindex_bak[k] = bondindex[k];
		}
		for (k = 0; k < atomnum2; k++) {
			valence2[k] = ind_valence[k];
			con_num[k] = ind_bondnum[k];
		}
		for (i = 0; i < bondnum2; i++)
			if (bondindex[i] == -1) {
				at1 = bond2[i].bondi;
				at2 = bond2[i].bondj;
				ind_valence[at1]--;
				ind_valence[at2]--;
				ind_bondnum[at1]--;
				ind_bondnum[at2]--;
				bond2[i].type = 1;
				bondindex[i] = 1;
				jbo_induce();
				score = jbo_score();
				if (debug == 1) {
					for (k = 0; k < bondnum2; k++)
						printf("\nsingle_b%5s%5s%5d",
							   atom2[bond2[k].bondi].name,
							   atom2[bond2[k].bondj].name, bond2[k].type);
					for (k = 0; k < atomnum2; k++)
						printf("\nsingle_a%5s%5d%5d%5d%5d", atom2[k].name,
							   ind_bondnum[k], con_num[k], ind_valence[k],
							   valence2[k]);
					printf("\nsingle %5d %5s %5s %5d", i + 1,
						   atom2[at1].name, atom2[at2].name, score);
				}
				flag = 0;
				for (k = 0; k < bondnum2; k++)
					if (bond2[k].type == -1)
						flag++;
				if (score == 0 && flag == 0)
					return score;	/*all the bond types are set */

				if (score > 0) {	/*error happens, bond type is set to 2 */
					for (k = 0; k < atomnum2; k++) {
						ind_valence[k] = valence2[k];
						ind_bondnum[k] = con_num[k];
					}
					for (k = 0; k < bondnum2; k++) {
						bond2[k].type = bondtype_bak[k];
						bondindex[k] = bondindex_bak[k];
					}
					ind_valence[at1] -= 2;
					ind_valence[at2] -= 2;
					ind_bondnum[at1]--;
					ind_bondnum[at2]--;
					bond2[i].type = 2;
					bondindex[i] = 1;
					jbo_induce();
					score = jbo_score();
					if (debug == 1) {
						for (k = 0; k < bondnum2; k++)
							printf("\ndouble_b%5s%5s%5d",
								   atom2[bond2[k].bondi].name,
								   atom2[bond2[k].bondj].name, bond2[k].type);
						for (k = 0; k < atomnum2; k++)
							printf("\ndouble_a%5s%5d%5d%5d%5d",
								   atom2[k].name, ind_bondnum[k],
								   con_num[k], ind_valence[k], valence2[k]);
						printf("\ndouble %5d %5s %5s %5d", i + 1,
							   atom2[at1].name, atom2[at2].name, score);
					}
					flag = 0;
					for (k = 0; k < bondnum2; k++)
						if (bond2[k].type == -1)
							flag++;
					if (score == 0 && flag == 0)
						return score;
					if (score > 0) {	/*error happens, bond type is set to 3 */
						for (k = 0; k < atomnum2; k++) {
							ind_valence[k] = valence2[k];
							ind_bondnum[k] = con_num[k];
						}
						for (k = 0; k < bondnum2; k++) {
							bond2[k].type = bondtype_bak[k];
							bondindex[k] = bondindex_bak[k];
						}
						ind_valence[at1] -= 3;
						ind_valence[at2] -= 3;
						ind_bondnum[at1]--;
						ind_bondnum[at2]--;
						bond2[i].type = 3;
						bondindex[i] = 1;
						jbo_induce();
						score = jbo_score();
						if (debug == 1) {
							for (k = 0; k < bondnum2; k++)
								printf("\ndouble_b%5s%5s%5d",
									   atom2[bond2[k].bondi].name,
									   atom2[bond2[k].bondj].name,
									   bond2[k].type);
							for (k = 0; k < atomnum2; k++)
								printf("\ndouble_a%5s%5d%5d%5d%5d",
									   atom2[k].name, ind_bondnum[k],
									   con_num[k], ind_valence[k],
									   valence2[k]);
							printf("\ndouble %5d %5s %5s %5d", i + 1,
								   atom2[at1].name, atom2[at2].name, score);
						}
						flag = 0;
						for (k = 0; k < bondnum2; k++)
							if (bond2[k].type == -1)
								flag++;
						if (score == 0 && flag == 0)
							return score;

						if (score > 0) {
							if (debug == 1)
								printf
									("\nCannot assign bond types for the current valence state");
							for (k = 0; k < bondnum2; k++) {
								bond2[k].type = bondtype_bak[k];
								bondindex[k] = bondindex_bak[k];
							}
							return score;
						}
					}			/* end of bond =3 */
				}				/* end of bond =2 */
			}					/* end of bond =1 */
		while_flag = 0;
		for (k = 0; k < bondnum2; k++)
			if (bond2[k].type == -1) {
				while_flag = 1;
				break;
			}
	}
	return score;
}

static int judgebt(void)
{
	int i, j;
	int fail_flag;
	int flag;
/*initalization*/
	for (i = -1; i < vastatenum; i++) {
		if (i != -1) {
			for (j = 0; j < atomnum2; j++)
				valence2[j] = va_best2[j];
			for (j = 0; j < va_num[i]; j++)
				valence2[va_addr[i][j] - 1] = va_valence[i][j];
		}
		current_va(i + 1);
		for (j = 0; j < atomnum2; j++) {
			con_num[j] = atom2[j].connum;
			ind_valence[j] = valence2[j];
			ind_bondnum[j] = con_num[j];
		}
		for (j = 0; j < bondnum2; j++) 
			if(bond2[j].jflag == 1) {
				bondindex[j] = 1;
				con_num[bond2[j].bondi]	--;
				con_num[bond2[j].bondj]	--;
				ind_valence[bond2[j].bondi] -= bond2[j].type;
				ind_valence[bond2[j].bondj] -= bond2[j].type;
				ind_bondnum[bond2[j].bondi] = con_num[bond2[j].bondi];
				ind_bondnum[bond2[j].bondj] = con_num[bond2[j].bondj];
			}
			else {
				bond2[j].type = -1;
				bondindex[j] = -1;
			}
		jbo_induce();
		fail_flag = 0;
		for (j = 0; j < bondnum2; j++)
			if (bondindex[j] == -1) {
				fail_flag = 1;
				break;
			}
		if (fail_flag == 1)
			fail_flag = jbo_iteration();
		flag = 0;
		for (j = 0; j < atomnum2; j++)
			if (ind_valence[j] != 0) {
				flag = 1;
				break;
			}
		if (flag == 1)
			fail_flag = 1;
		flag = 0;
		for (j = 0; j < atomnum2; j++)
			if (ind_bondnum[j] != 0) {
				flag = 1;
				break;
			}
		if (flag == 1)
			fail_flag = 1;

		if (fail_flag == 0) {
			if (debug) {
				printf("\n Assigned bond types for valence state No %d", i + 1);
				for (j = 0; j < bondnum2; j++) {
					if(bond2[j].type2 < 0)
						printf("\n%5d%5d%5d%5d", j + 1, bondindex[j],
						   	bond2[j].type, -1);
					else
						printf("\n%5d%5d%5d%5d", j + 1, bondindex[j],
						   	bond2[j].type, org_bondtype[bond2[j].type2]);
				}
			}
			for (j = 0; j < bondnum2; j++)  
				if(bond2[j].type2 != -1)
					bond[bond2[j].type2].type = bond2[j].type;
			if (i != -1)
				printf
					("\nInfo: Bond types are assigned for valence state %d with penalty of %d\n",
					 i + 1, va_penalty[i]);
			return 1;
		}
	}
	return 0;
}
static void finalize(void)
{
	int i, j, k;
	int bondi, bondj;
	int flag, flag0, flag1, flag2 ;
	int num;
	int index[6];
	if (judge_flag == 0) {
		for (i = 0; i < bondnum; i++) {
			if(bond[i].type != 10) 
				continue;
			bondi = bond[i].bondi;
			bondj = bond[i].bondj;
			flag0 = 0;
			for (j = 0; j < bondnum; j++) {
				if(i==j) continue;
				if(bond[j].bondi == bondi || bond[j].bondj == bondi ||
				   bond[j].bondi == bondj || bond[j].bondj == bondj) 
					if(bond[j].type == 2 || bond[j].type == 8) {
						flag0 = 1;
						break;
					}
			}
			if(flag0 == 0)
				bond[i].type = 8;
			else
				bond[i].type = 7;
		}
	}
	for (i = 0; i < bondnum; i++) {
		bondi = bond[i].bondi;
		bondj = bond[i].bondj;
/*part1*/
		if ((arom[bondi].ar1 > 0 && arom[bondj].ar1 > 0) ||
			(arom[bondi].ar2 > 0 && arom[bondj].ar2 > 0))
			for (j = 0; j < ringnum; j++) {
				if (ring[j].num >= 7)
					continue;
				if (ring[j].num <= 4)
					continue;
				flag = 0;
				for (k = 0; k < ring[j].num; k++)
					if (arom[ring[j].atomno[k]].ar1 <= 0
						&& arom[ring[j].atomno[k]].ar2 <= 0) {
						flag = 1;
						break;
					}
				if (flag == 1)
					continue;
				flag1 = 0;
				flag2 = 0;
				for (k = 0; k < ring[j].num; k++) {
					if (ring[j].atomno[k] == bondi)
						flag1 = 1;
					if (ring[j].atomno[k] == bondj)
						flag2 = 1;
				}
				if (flag1 == 1 && flag2 == 1) {
					if (bond[i].type == 1)
						bond[i].type = 7;
					if (bond[i].type == 2)
						bond[i].type = 8;
					break;
				}
			}
/*part2*/
		if (bond[i].type == 2) {
			if (conjatom[bondi] == 1 ) {
				if (atom[bondj].connum == 1 && atom[bondj].atomicnum == 8) {
					bond[i].type = 9;
					continue;
				}
				if (atom[bondj].connum == 1 && atom[bondj].atomicnum == 16) {
					bond[i].type = 9;
					continue;
				}
			}
			if (conjatom[bondj] == 1) {
				if (atom[bondi].connum == 1 && atom[bondi].atomicnum == 8) {
					bond[i].type = 9;
					continue;
				}
				if (atom[bondi].connum == 1 && atom[bondi].atomicnum == 16) {
					bond[i].type = 9;
					continue;
				}
			}
		}
/*part3*/
		if ((atom[bondi].connum == 2 || atom[bondi].connum == 3)
			&& atom[bondi].atomicnum == 7
			&& ((atom[bondj].connum == 1 && atom[bondj].atomicnum == 8)
				|| (atom[bondj].connum == 1
					&& atom[bondj].atomicnum == 16))) {
			for(j = 0; j <= 5; j++) {
				if(atom[bondi].con[j] < 0) 
					break;
				if(atom[bondi].con[j] == bondj) 
					continue;
				if(atom[atom[bondi].con[j]].connum == 1 && 
					(atom[atom[bondi].con[j]].atomicnum == 8 || atom[atom[bondi].con[j]].atomicnum == 16))
					bond[i].type = 6;
					break;
			}
			if(bond[i].type == 6)
				continue;
		}
		if ((atom[bondj].connum == 2 || atom[bondj].connum == 3)
			&& atom[bondj].atomicnum == 7
			&& ((atom[bondi].connum == 1 && atom[bondi].atomicnum == 8)
				|| (atom[bondi].connum == 1
					&& atom[bondi].atomicnum == 16))) {
                        for(j = 0; j <= 5; j++) {
                                if(atom[bondj].con[j] < 0)
                                        break;
                                if(atom[bondj].con[j] == bondi)
                                        continue;
                                if(atom[atom[bondj].con[j]].connum == 1 &&
                                        (atom[atom[bondj].con[j]].atomicnum == 8 || atom[atom[bondj].con[j]].atomicnum == 16))
                                        bond[i].type = 6;
                                        break;
                        }
			bond[i].type = 6;
			continue;
		}
/*part4*/
		if (bond[i].type == 1)
			if ((atom[bondi].atomicnum == 8 && atom[bondi].connum == 1) ||
				(atom[bondj].atomicnum == 8 && atom[bondj].connum == 1) ||
				(atom[bondi].atomicnum == 16 && atom[bondi].connum == 1) ||
				(atom[bondj].atomicnum == 16 && atom[bondj].connum == 1)) {
				bond[i].type = 9;
				continue;
			}
	}
/*part5*/
	for (i = 0; i < atomnum; i++) {
		if (atom[i].connum == 3 && atom[i].atomicnum == 16) {
			num = 0;
			for (j = 0; j <= 5; j++)
				index[j] = 0;
			for (j = 0; j <= 2; j++) {
				bondi = atom[i].con[j];
				if (atom[bondi].connum == 1
					&& (atom[bondi].atomicnum == 8
						|| atom[bondi].atomicnum == 16)) {
					num++;
					index[j] = 1;
				}
			}
			if (num == 2)
				for (j = 0; j < bondnum; j++) {
					bondi = bond[j].bondi;
					bondj = bond[j].bondj;
					for (k = 0; k <= 2; k++)
						if (index[k] == 1
							&& ((bondi == i && bondj == atom[i].con[k])
								|| (bondj == i
									&& bondi == atom[i].con[k])))
							bond[j].type = 9;
				}
		}
		if (atom[i].connum == 4 && atom[i].atomicnum == 15) {
			num = 0;
			for (j = 0; j <= 5; j++)
				index[j] = 0;
			for (j = 0; j <= 3; j++) {
				bondi = atom[i].con[j];
				if (atom[bondi].connum == 1
					&& (atom[bondi].atomicnum == 8
						|| atom[bondi].atomicnum == 16)) {
					num++;
					index[j] = 1;
				}
			}
			if (num >= 2)
				for (j = 0; j < bondnum; j++) {
					bondi = bond[j].bondi;
					bondj = bond[j].bondj;
					for (k = 0; k <= 3; k++)
						if (index[k] == 1
							&& ((bondi == i && bondj == atom[i].con[k])
								|| (bondj == i
									&& bondi == atom[i].con[k])))
							bond[j].type = 9;
				}
		}
	}
}

/* the valence is <= 8, hence */
static size_t point_to_size_t(int addr, int valence)
{
    assert(addr > 0);
    assert(valence >= 0 && valence < 8);

    size_t sa = (size_t) addr;
    size_t sv = (size_t) valence;

    return (sa << 3) | (0x7 & sv);
}

static void size_t_to_point(size_t x, int *addr, int *valence)
{
    assert(x > 7);
    assert(addr != NULL);
    assert(valence != NULL);

    *addr = (int) ((x & ~0x7) >> 3);
    *valence = (int) (0x7 & x);
}

static void *ncsu_malloc(size_t n)
{
    void *ptr;

    assert(n > 0);

    ptr = malloc(n);
    if (ptr == NULL) {
        fputs("** Error ** : malloc() failed\n", stdout);
        exit(EXIT_FAILURE);
    }

    return ptr;
}

static void *ncsu_realloc(void* old, size_t n)
{
    void *ptr;

    assert(n > 0);

    ptr = (old == NULL ? malloc(n) : realloc(old, n));

    if (ptr == NULL) {
        fputs("** Error ** : realloc() failed\n", stdout);
        exit(EXIT_FAILURE);
    }

    return ptr;
}

static void va_arrays_filler
    (size_t npairs, const ncsu_penalty_pair_t *pairs, void *data)
{
    size_t n;
    size_t **points;

    assert(npairs != 0);
    assert(pairs != NULL);
    assert(data != NULL);

    points = (size_t **) data;

    if (vastatenum == va_allocated) {
        va_allocated = (vastatenum == 0 ? 256 : va_allocated << 1);

        va_num = (int*) ncsu_realloc(va_num, va_allocated*sizeof(int));

        va_addr = (int**) ncsu_realloc(va_addr, va_allocated*sizeof(int*));
        va_valence =
            (int**) ncsu_realloc(va_valence, va_allocated*sizeof(int*));

        va_penalty = (int*) ncsu_realloc(va_penalty, va_allocated*sizeof(int));
    }

    va_num[vastatenum] = (int) npairs;

    va_addr[vastatenum] = (int *) ncsu_malloc(sizeof(int)*npairs);
    va_valence[vastatenum] = (int *) ncsu_malloc(sizeof(int)*npairs);

    va_penalty[vastatenum] = 0;

    for (n = 0; n < npairs; ++n) {
        size_t p = pairs[n].set_id + 1;
        size_t o = pairs[n].element_id;

        size_t_to_point(points[p][o],
            va_addr[vastatenum] + n, va_valence[vastatenum] + n);
        va_penalty[vastatenum] += (int) p;

        if (debug)
            printf(" (%d, %d, %d)", va_addr[vastatenum][n],
                   (int) p, va_valence[vastatenum][n]);
    }

    if (debug)
        putchar('\n');

    ++vastatenum;
}

static void process(void)
{
    size_t m;

    size_t **points;
    size_t *npoints, *ii;

    if (maxaps <= 0)
        goto judge;

    npoints = (size_t *) ncsu_malloc(sizeof(size_t)*(maxaps + 1));
    for (m = 0; m < maxaps + 1; ++m)
        npoints[m] = 0;

    { /* count the number of points with 0 < penalty <= maxaps */
        int i, j;

        for (i = 0; i < atomnum2; ++i)
            for (j = 0; j < 8; ++j) {
                int p = av2[i].aps[j];
                if (p > 0 && p <= maxaps)
                   ++(npoints[p]);
            }
    }

    ii = (size_t *) ncsu_malloc(sizeof(size_t)*(maxaps + 1));
    points = (size_t **) ncsu_malloc(sizeof(size_t*)*(maxaps + 1));
    points[0] = NULL;
    for (m = 1; m < maxaps + 1; ++m) {
        ii[m] = 0;
        points[m] = (npoints[m] == 0 ? NULL :
           (size_t *) ncsu_malloc(sizeof(size_t)*npoints[m]));
    }

    { /* fill in the "points" */
        int i, j;

        for (i = 0; i < atomnum2; ++i)
            for (j = 0; j < 8; ++j) {
                int p = av2[i].aps[j];
                if (p > 0 && p <= maxaps) {
                    int addr = i + 1;
                    int valence = j;

                    assert(ii[p] < npoints[p]);
                    points[p][ii[p]] = point_to_size_t(addr, valence);
                    ++(ii[p]);
                }
            }
    }

    /* fill the va_* arrays */
    for (m = 1; m <= maxaps; ++m) {
        if (debug)
            printf("\nProcessing Penalty of %lu ...\n", (unsigned long) m);

        ncsu_penalties(maxaps, npoints + 1, m, va_arrays_filler, points);

        if (debug)
            printf("\nThe number of unduplicated VA states are %d so far",
                   vastatenum);
    } /* for (m = 1; m <= maxaps; ++m) */

    if (debug) {
        int i, j;

        printf("\n\nList of valence states with penalty score from "
               "the lowest to highest\n");
        for (i = 0; i < vastatenum; ++i)
            for (j = 0; j < va_num[i]; ++j)
                printf("\n%5d%5d%5d%5d%5d", i + 1, j + 1,
                       va_penalty[i], va_addr[i][j], va_valence[i][j]);
    }

    for (m = 0; m < maxaps + 1; ++m) {
        if (points[m] != NULL) {
            free(points[m]);
        }
    }

    free(points);
    free(ii);
    free(npoints);

judge:
    if (judgebt() == 0) {
        puts("\nWarning: the assigned bond types may be wrong, please : "
             "\n(1) double check the structure (the connectivity) and/or "
             "\n(2) adjust atom valence penalty parameters in APS.DAT, and/or "
             "\n(3) increase PSCUTOFF in define.h and recompile bondtype.c"
             "\n    Be cautious, use a large value of PSCUTOFF (>100) will "
             "significantly increase the computation time");
        wac(ofilename, atomnum, atom, bondnum, bond, cinfo, minfo);
//	this is just a warning message, the bond types/atom types are likely right
//	should not exit from here, otherwise we get no output at all.
//       exit(EXIT_FAILURE);
    }
}

int main(int argc, char *argv[])
{
	int i, j;
	int tmpint1, tmpint2;
	int diff;
	int connum;
	int resno;
	int tvalence;
	int overflow_flag = 0;			/*if overflow_flag ==1, reallocate memory */

    amberhome = (char *) getenv("AMBERHOME");
    if( amberhome == NULL ){
       fprintf( stdout, "AMBERHOME is not set!\n" );
       exit(1);
    }
	if (strcmp(COLORTEXT, "YES") == 0 || strcmp(COLORTEXT, "yes") == 0) {
		if (argc == 2
			&& (strcmp(argv[1], "-h") == 0
				|| strcmp(argv[1], "-H") == 0)) {
			printf("[31mUsage: bondtype -i[0m input file name \n"
				   "[31m                -o[0m output file name \n"
				   "[31m                -f[0m file format (ac or mol2)\n"
				   "[31m                -j[0m judge bond type level option, default is part\n"
				   "[32m                   full [0m full judgement\n"
				   "[32m                   part [0m partial judgement, only do reassignment according\n"
				   "		         to known bond type information in the input file\n"
				   "[31m                -s[0m stop running if APS (atomic penality score) is not available\n"
				   "[32m                   0 [0m- no, all the related bonds are frozen, the default\n"
				   "[32m                   1 [0m- yes\n");
			exit(1);
		}
		if (argc != 7 && argc != 9 && argc != 11) {
			printf("[31mUsage: bondtype -i[0m input file name \n"
				   "[31m                -o[0m output file name \n"
				   "[31m                -f[0m file format (ac or mol2)\n"
				   "[31m                -j[0m judge bond type level option, default is part\n"
				   "[32m                   full [0m full judgement\n"
				   "[32m                   part [0m partial judgement, only do reassignment according\n"
				   "		         to known bond type information in the input file\n"
				   "[31m                -s[0m stop running if APS (atomic penality score) is not available\n"
				   "[32m                   0 [0m- no, all the related bonds are frozen, the default\n"
				   "[32m                   1 [0m- yes\n");
			exit(1);
		}
	} else {
		if (argc == 2)
			if (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "-H") == 0) {
				printf("Usage: bondtype -i  input file name \n");
				printf("                -o  output file name\n");
				printf("                -f  file format (ac or mol2) \n");
				printf("                -j  judge bond type level option, default is part\n");
				printf("                    full: full judgement\n");
				printf("                    part: partial judgement, only do reassignment according\n");
				printf(" 		          to known bond type information in the input file\n");
				printf("                -s  stop running if APS (atomic penality score) is not available\n");
				printf(" 	            0 - no, all the related bonds are frozen, the default; 1 - yes\n");
				exit(1);
			}
		if (argc != 7 && argc != 9) {
			printf("Usage: bondtype -i  input file name \n");
			printf("                -o  output file name\n");
			printf("                -f  file format (ac or mol2) \n");
			printf("                -j  judge bond type level option, default is part\n");
			printf("                    full: full judgement\n"); 
			printf("                    part: partial judgement, only do reassignment according\n");
			printf(" 		          to known bond type information in the input file\n");
			printf("                -s  stop running if APS (atomic penality score) is not available\n");
			printf(" 	            0 - no, all the related bonds are frozen, the default; 1 - yes\n");
			exit(1);
		}
	}
	for (i = 1; i < argc; i += 2) {
		if (strcmp(argv[i], "-i") == 0) {
			strcpy(ifilename, argv[i + 1]);
			continue;
		}
		if (strcmp(argv[i], "-o") == 0) {
			strcpy(ofilename, argv[i + 1]);
			continue;
		}
		if (strcmp(argv[i], "-f") == 0) {
			if (strcmp(argv[i + 1], "ac") == 0
				|| strcmp(argv[i + 1], "AC") == 0)
				iformat = 1;
			if (strcmp(argv[i + 1], "mol2") == 0
				|| strcmp(argv[i + 1], "MOL2") == 0)
				iformat = 2;
			continue;
		}
		if (strcmp(argv[i], "-j") == 0) {
			if (strcmp(argv[i + 1], "full") == 0
				|| strcmp(argv[i + 1], "Full") == 0
				|| strcmp(argv[i + 1], "FULL") == 0)
				judge_flag = 1;
			if (strcmp(argv[i + 1], "part") == 0
				|| strcmp(argv[i + 1], "Part") == 0
				|| strcmp(argv[i + 1], "PART") == 0)
				judge_flag = 0;
			continue;
		}
		if (strcmp(argv[i], "-s") == 0) {
			iaps = atoi(argv[i + 1]);
			if(iaps != 0 && iaps != 1) {
				fprintf(stdout, "-s flag must be 0 or 1\n");
				break;
			}
			continue;
		}
	}

	avH.aps[0] = 64;
	avH.aps[1] = 0;
	avH.aps[2] = 64;
	avH.aps[3] = 9999;
	avH.aps[4] = 9999;
	avH.aps[5] = 9999;
	avH.aps[6] = 9999;
	avH.aps[7] = 9999;

	default_cinfo(&cinfo);
	default_minfo(&minfo);
	atom = (ATOM *) malloc(sizeof(ATOM) * cinfo.maxatom);
	if (atom == NULL) {
		fprintf(stdout, "memory allocation error for *atom\n");
		exit(1);
	}
	arom = (AROM *) malloc(sizeof(AROM) * cinfo.maxatom);
	if (arom == NULL) {
		fprintf(stdout, "memory allocation error for *arom\n");
		exit(1);
	}
	bond = (BOND *) malloc(sizeof(BOND) * cinfo.maxbond);
	if (bond == NULL) {
		fprintf(stdout, "memory allocation error for *bond\n");
		exit(1);
	}
	for (i = 0; i < cinfo.maxbond; ++i) {
		bond[i].jflag = -1; /* bond type has not been assigned */
	}
	ring = (RING *) malloc(sizeof(RING) * cinfo.maxring);
	if (ring == NULL) {
		fprintf(stdout, "memory allocation error for *ring\n");
		exit(1);
	}
	if (iformat == 1)
		overflow_flag =
			rac(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo);
	if (iformat == 2)
		overflow_flag =
			rmol2(ifilename, &atomnum, atom, &bondnum, bond, &cinfo,
				  &minfo, 0);
	if (overflow_flag) {
		cinfo.maxatom = atomnum + 10;
		cinfo.maxbond = bondnum + 10;
		free(atom);
		free(arom);
		free(bond);
		atom = (ATOM *) malloc(sizeof(ATOM) * cinfo.maxatom);
		if (atom == NULL) {
			fprintf(stdout, "memory allocation error for *atom\n");
			exit(1);
		}
		arom = (AROM *) malloc(sizeof(AROM) * cinfo.maxatom);
		if (arom == NULL) {
			fprintf(stdout, "memory allocation error for *arom\n");
			exit(1);
		}
		bond = (BOND *) malloc(sizeof(BOND) * cinfo.maxbond);
		if (bond == NULL) {
			fprintf(stdout, "memory allocation error for *bond\n");
			exit(1);
		}
		int i;
		for (i = 0; i < cinfo.maxbond; ++i) {
			bond[i].jflag = -1; /* bond type has not been assigned */
		}
		if (iformat == 1)
			rac(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo);
		if (iformat == 2)
			rmol2(ifilename, &atomnum, atom, &bondnum, bond, &cinfo,
				  &minfo, 0);
	}
	atomicnum(atomnum, atom);
	adjustatomname(atomnum, atom, 1); 
	overflow_flag = ringdetect(atomnum, atom, bondnum, bond, &ringnum, ring, arom,
			   cinfo.maxatom, cinfo.maxring, minfo.inf_filename, 1);
	if(overflow_flag == 1) {
		cinfo.maxring = ringnum + 10;
		free(ring);
		ring = (RING *) malloc(sizeof(RING) * cinfo.maxring);
		if (ring == NULL) {
			fprintf(stdout, "memory allocation error for *ring\n");
			exit(1);
		}
		overflow_flag = ringdetect(atomnum, atom, bondnum, bond, &ringnum, ring, arom,
			   	cinfo.maxatom, cinfo.maxring, minfo.inf_filename, 1);
	}
	if (debug)
		info(atomnum, atom, bondnum, bond, arom, cinfo, minfo);

/* allocate memory for org_bondtype, av, valence, va_best, conjatom, seq and res*/
	org_bondtype = (int *) malloc(sizeof(int) * bondnum);
	if (org_bondtype == NULL) {
		fprintf(stdout, "memory allocation error for *org_bondtype\n");
		exit(1);
	}

	av = (AV *) malloc(sizeof(AV) * atomnum);
	if (av == NULL) {
		fprintf(stdout, "memory allocation error for *av\n");
		exit(1);
	}

	valence = (int *) malloc(sizeof(int) * atomnum);
	if (valence == NULL) {
		fprintf(stdout, "memory allocation error for *valence\n");
		exit(1);
	}

	va_best = (int *) malloc(sizeof(int) * atomnum);
	if (va_best == NULL) {
		fprintf(stdout, "memory allocation error for *va_best\n");
		exit(1);
	}

	conjatom = (int *) malloc(sizeof(int) * atomnum);
	if (conjatom == NULL) {
		fprintf(stdout, "memory allocation error for *conjatom\n");
		exit(1);
	}

	seq = (int *) malloc(sizeof(int) * atomnum);
	if (seq == NULL) {
		fprintf(stdout, "memory allocation error for *seq\n");
		exit(1);
	}

	apstype = (int *) malloc(sizeof(int) * atomnum);
	if (apstype == NULL) {
		fprintf(stdout, "memory allocation error for *apstype\n");
		exit(1);
	}

	res = (RES *) malloc(sizeof(res) * MAXRES);
	if (res == NULL) {
		fprintf(stdout, "memory allocation error for *res\n");
		exit(1);
	}

	for (i = 0; i < bondnum; i++) 
        	org_bondtype[i] = bond[i].type;
	for (i = 0; i < atomnum; i++) {
		atom[i].type = -1;
		apstype[i] = -1;
		seq[i] = -1;
	}
/* assign atomic valences */
	assignav();
	if (judge_flag == 0) {
		finalize();
		wac(ofilename, atomnum, atom, bondnum, bond, cinfo, minfo);
		printf("\n");
		return (0);
	}
/* if judge_flag == 1 */
/*find maximum atomnum and bondnum in a residue */
	resno = -9999999;
	resnum = 0;
	for(i=0;i<atomnum;i++) {
		if(atom[i].resno != resno) {
			res[resnum].resno = atom[i].resno;
			strcpy(res[resnum].name, atom[i].aa);
			resno = atom[i].resno;
			resnum++;
                        if(resnum == MAXRES) {
                                res = (RES *) realloc(res, sizeof(RES) * (MAXRES + resnum));
                                if (res == NULL) {
                                        fprintf(stdout, "Error: memory reallocation with realloc for *res\n");
                                        exit(1);
                                }
                        }
		}
	}
	if(debug) 
		for(i=0;i<resnum;i++) 
			printf("\nRES %5d %5d %10s",  i+1, res[i].resno, res[i].name);

	maxatomnum2 = 0;
	maxbondnum2 = 0;
	for(i=0;i<resnum;i++) {
		tmpint1 = 0;
		tmpint2 = 0;
		for(j=0;j<atomnum;j++) 
			if(atom[j].resno == res[i].resno) tmpint1++;
		for(j=0;j<bondnum;j++) 
			if(atom[bond[j].bondi].resno == res[i].resno && 
			   atom[bond[j].bondj].resno == res[i].resno) 
				tmpint2++;
		if(tmpint1 > maxatomnum2) 
			maxatomnum2 = tmpint1;
		if(tmpint2 > maxbondnum2) 
			maxbondnum2 = tmpint2;
	}

/* allocate memory, automatically increased 10 in case of adding cap Hydrogens */	
	memory(maxatomnum2 + 10, maxbondnum2 + 10); 

/* for those bonds linking two residues, atom type is automatically set to 1 */
	for(i=0;i<bondnum;i++) 
		if(bond[i].jflag != 1 && atom[bond[i].bondi].resno != atom[bond[i].bondj].resno) {
			bond[i].type = 1;
			atom[bond[i].bondi].type = 1;	
			atom[bond[i].bondj].type = 1;	
		}
/* judge bond types for each residue */
        for(i=0;i<resnum;i++) {
/*first prepare atom2 and bond2 */
		if(resnum > 1) 
			printf("\n---Judge bond type for Residue %d with ID of %d and Name of %s ---\n", i+1, res[i].resno, res[i].name);
		tmpint1 = 0;
		diff = 0;
        	for(j=0;j<atomnum;j++) 
			if(atom[j].resno == res[i].resno) {
				atom2[tmpint1] = atom[j];
				av2[tmpint1] = av[j];
				valence2[tmpint1] = valence[j];
				va_best2[tmpint1] = va_best[j];
				apstype2[tmpint1] = apstype[j];
				seq[j] = diff;
				tmpint1++;
			}
			else diff++;

		atomnum2 = tmpint1;

/* adjust the bonded atoms atom[].con[], delete neighbours that belong to other residues */
        	for (j = 0; j < atomnum2; j++)  {
			connum = 0;
                	for (k = 0; k < 6; k++) {
				if(atom2[j].con[k] == -1) continue;
				if(atom[atom2[j].con[k]].resno != res[i].resno) continue; 
				atom2[j].con[connum++] = atom2[j].con[k] - seq[atom2[j].con[k]];
			}
			atom2[j].connum = connum;
		}
/* find maxaps */
	        maxaps = 0;
        	for (j = 0; j < atomnum2; j++) 
                	for (k = 0; k < 8; k++) {
                        	if (av2[j].aps[k] > PSCUTOFF) continue;
                        	if (maxaps < av2[j].aps[k])
                                	maxaps = av2[j].aps[k];
			}	
		if (debug) printf("\n\nThe Maximum APS is %5d\n", maxaps);
/* handle bond now */
		tmpint2 = 0;

        	for(j=0;j<bondnum;j++) 
			if(atom[bond[j].bondi].resno == res[i].resno && atom[bond[j].bondj].resno == res[i].resno) {
				bond2[tmpint2] = bond[j];
				bond2[tmpint2].bondi -= seq[bond[j].bondi];
				bond2[tmpint2].bondj -= seq[bond[j].bondj];
/*bond2[].type2 referring to the original bond id in bond[] */
				bond2[tmpint2].type2 = j;
				tmpint2++;
			}
		bondnum2 = tmpint2; 

/* handle those linking atoms of two residues, such as C and N in amino acid residues, replaceing those atom with H*/
        	for(j=0;j<bondnum;j++) {
			if(atom[bond[j].bondi].resno == res[i].resno && atom[bond[j].bondj].resno != res[i].resno) {
				bond2[bondnum2] = bond[j];
				bond2[bondnum2].bondi -= seq[bond[j].bondi];
				bond2[bondnum2].bondj = atomnum2 ;
				bond2[bondnum2].type = 1;
				bond2[bondnum2].type2 = -1; /* not refer to original bond */
/* create a new atom H*/
				atom2[atomnum2] = atom[bond[j].bondj];	
				strcpy(atom2[atomnum2].name, "H");
				strcpy(atom2[atomnum2].element, "H");
				strcpy(atom2[atomnum2].aa, res[i].name);
				strcpy(atom2[atomnum2].ambername, "H");
				atom2[atomnum2].resno = res[i].resno;
				atom2[atomnum2].atomicnum = 1;
				atom2[atomnum2].connum = 1;
				atom2[atomnum2].con[0] = bond2[bondnum2].bondi;
				atom2[atomnum2].con[1] = -1;
				atom2[atomnum2].con[2] = -1;
				atom2[atomnum2].con[3] = -1;
				atom2[atomnum2].con[4] = -1;
				atom2[atomnum2].con[5] = -1;
				atom2[bond2[bondnum2].bondi].con[atom2[bond2[bondnum2].bondi].connum] = atomnum2;
				atom2[bond2[bondnum2].bondi].connum ++;
				av2[atomnum2] = avH;
				valence2[atomnum2] = 1;
				va_best2[atomnum2] = 1;
				apstype2[atomnum2] = 1;
				atomnum2++;
				bondnum2++;
			}
			if(atom[bond[j].bondi].resno != res[i].resno && atom[bond[j].bondj].resno == res[i].resno) {
				bond2[bondnum2] = bond[j];
				bond2[bondnum2].bondi = atomnum2 ;
				bond2[bondnum2].bondj -= seq[bond[j].bondj];
				bond2[bondnum2].type = 1;
				bond2[bondnum2].type2 = -1;
/* create a new atom H*/
				atom2[atomnum2] = atom[bond[j].bondi];	
				strcpy(atom2[atomnum2].name, "H");
				strcpy(atom2[atomnum2].element, "H");
				strcpy(atom2[atomnum2].aa, res[i].name);
				strcpy(atom2[atomnum2].ambername, "H");
				atom2[atomnum2].resno = res[i].resno;
				atom2[atomnum2].atomicnum = 1;
				atom2[atomnum2].connum = 1;
				atom2[atomnum2].con[0] = bond2[bondnum2].bondj;
				atom2[atomnum2].con[1] = -1;
				atom2[atomnum2].con[2] = -1;
				atom2[atomnum2].con[3] = -1;
				atom2[atomnum2].con[4] = -1;
				atom2[atomnum2].con[5] = -1;
				atom2[bond2[bondnum2].bondj].con[atom2[bond2[bondnum2].bondj].connum] = atomnum2;
				atom2[bond2[bondnum2].bondj].connum ++;
				av2[atomnum2] = avH;
				valence2[atomnum2] = 1;
				va_best2[atomnum2] = 1;
				apstype2[atomnum2] = 1;
				atomnum2++;
				bondnum2++;
			}

		}

/* handle pre-defined bond types */
		for(j=0; j<bondnum2; j++) {
			if(apstype2[bond2[j].bondi] == -1 || apstype2[bond2[j].bondj] == -1)
				bond2[j].jflag = 1;
/* there are other cases that jflag reads in from the input file*/
			if(bond2[j].jflag == 1) {
				if(bond2[j].type == 7) bond2[j].type = 1;
				if(bond2[j].type == 8) bond2[j].type = 2;
				if (bond2[j].type != 1 && bond2[j].type != 2 && bond2[j].type != 3 &&
				    bond2[j].type != 7 && bond2[j].type != 8) {
					printf("\nThe frozen atom type can only be 1, 2, 3, 7 (aromatic single), 8 (aromatic double)"); 
					exit(1);
				}
			}
		}
		for(j=0; j<atomnum2; j++) {
			tvalence = 0;
			connum = 0;
			for(k=0; k<bondnum2; k++) 
				if(bond2[k].bondi == j || bond2[k].bondj == j) {
					tvalence += bond2[k].type;
					if(bond2[k].jflag == 1) connum ++;
				}
			if(apstype2[j] == -1) {
/* no av is assigned */
				av2[j].aps[tvalence] = 0;
				valence2[j] = tvalence;
				va_best2[j] = tvalence;
			}
			else if(atom2[j].connum == connum) {
/* all the related bonds are frozen, then reassign av2*/
				av2[j].aps[0] = 9999;
				av2[j].aps[1] = 9999;
				av2[j].aps[2] = 9999;
				av2[j].aps[3] = 9999;
				av2[j].aps[4] = 9999;
				av2[j].aps[5] = 9999;
				av2[j].aps[6] = 9999;
				av2[j].aps[7] = 9999;
				av2[j].aps[tvalence] = 0;	
				valence2[j] = tvalence;
				va_best2[j] = tvalence;
			}
		}
		process();   
}	
	for(i=0;i<bondnum;i++)
		if(org_bondtype[i] == 0) {
		        overflow_flag = ringdetect(atomnum, atom, bondnum, bond, &ringnum, ring, arom,
                        	   cinfo.maxatom, cinfo.maxring, minfo.inf_filename, 1);
			if(overflow_flag == 1) {
				cinfo.maxring = ringnum + 10;
				free(ring);
				ring = (RING *) malloc(sizeof(RING) * cinfo.maxring);
				if (ring == NULL) {
					fprintf(stdout, "memory allocation error for *ring\n");
					exit(1);
				}
		        	overflow_flag = ringdetect(atomnum, atom, bondnum, bond, &ringnum, ring, arom,
                        	   		cinfo.maxatom, cinfo.maxring, minfo.inf_filename, 1);
			}
			break;
		}

	finalize();
	wac(ofilename, atomnum, atom, bondnum, bond, cinfo, minfo);
//	putchar('\n');
//	putchar('\n');
	return (0);
}
