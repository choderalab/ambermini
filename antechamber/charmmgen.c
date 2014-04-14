/*
************************************************************************
*           All Copyright Reserved!                                    *
*                                                                      *
*  Prog:    charmmgen                                                  *
*  Version: version 2.0                                                *
*  Author:  Victor E. Bazterra                                         *
*                                                                      *
*  Center for high performance computing                               *
*  University of Utah  		                                       *
*  Salt Lake City   UT 84114                                           *
*  June, 2004                                                          *
************************************************************************
*/

char *amberhome;
# include <math.h>
# include "common.h"
# include "define.h"
# include "atom.h"
# include "utility.c"
# include "common.c"
# include "rotate.c"
# include "ac.c"
# include "mol2.c"
# include "prep.c"

#define MAX_FF_CORR 250
#define MAX_FF_ATOMTYPE 250
#define MAX_FF_VDW 250
#define MAX_FF_BOND 2000
#define MAX_FF_ANGLE 5000
#define MAX_FF_TORSION 1500
#define MAX_FF_IMPROPER 500
#define MAX_EQU_VDW 20
#define MAX_EQU_VDW_NUM 50

typedef struct {
	int atid1;
	int atid2;
	int atid3;
	int atid4;
} IMPROPERID;

typedef struct {
	char name1[5];
	char name2[5];
	char name3[5];
	char name4[5];
	int num;
} CORR;

typedef struct {
	char name[5];
	double mass;
	double pol;
        char flag;
} ATOMTYPE;

typedef struct {
	char name1[5];
	char name2[5];
	double length;
	double force;
        char flag;
} BOND_PARAM;

typedef struct {
	char name1[5];
	char name2[5];
	char name3[5];
	double angle;
	double force;
        char flag;
} ANGLE;

typedef struct {
	char name1[5];
	char name2[5];
	char name3[5];
	char name4[5];
	int mul;
	double fterm;
	double phase;
	double force;
        char flag;
} TORSION;

typedef struct {
	char name1[5];
	char name2[5];
	char name3[5];
	char name4[5];
	double fterm;
	double phase;
	double force;
	int mul;
        char flag;
} IMPROPER;

typedef struct {
	char name[5];
	double radius;
	double pot;
        char flag;
} VDW;

typedef struct {
	char name[MAX_EQU_VDW_NUM][5];
	int num;
} EQU_VDW;

MOLINFO minfo;
CONTROLINFO cinfo;
ATOM *atom;
BOND *bond_array;
int atomnum = 0;
int bondnum = 0;
NAME *corrname;
NAME *similarname;
NAME *similarname2;
int *corrindex;
int *similarindex;
int *similarindex2;
IMPROPERID *improper;
char ifilename[MAXCHAR];
char ofilename[MAXCHAR];
char ffilename[MAXCHAR];
char hfilename[MAXCHAR];
char pfilename[MAXCHAR];
char cfilename[MAXCHAR];
int cindex;
char line[MAXCHAR];
double charge;
int impropernum = 0;
int i, j, k, l;
int maxcorr;
int maxatomtype;
int maxvdwparm;
int maxbondparm;
int maxangleparm;
int maxtorsionparm;
int maximproperparm;
int corrnum;
int atomtypenum;
int vdwparmnum;
int bondparmnum;
int angleparmnum;
int torsionparmnum;
int improperparmnum;
/*H-1, C-2, N-3, O-4, F-5, Cl-6, Br-7, I-8, S-9, P-10*/

int *improperindex;

CORR *corr;
ATOMTYPE   *atomtype;
BOND_PARAM *bondparm;
ANGLE      *angleparm;
TORSION    *torsionparm;
IMPROPER   *improperparm;
VDW        *vdwparm;

/* for equivalent vdw types */
EQU_VDW equ_vdw[MAX_EQU_VDW];
int equ_vdw_num = 0;

FILE *fp, *fpout, *ffout, *fhout;

char resname[MAXCHAR];

void improper_prediction(void)
{
	int i, j;
	for (i = 0; i < atomnum; i++)
		for (j = 0; j < corrnum; j++)
			if (strcmp(atom[i].ambername, corr[j].name1) == 0
				&& improperindex[j] == 1)
				atom[i].improper = 1;
}

void improper_id1(char *filename)
{
	FILE *fpin;
	int i;
	int readindex;
	char tmpchar1[MAXCHAR];
	char tmpchar2[MAXCHAR];
	char tmpchar3[MAXCHAR];
	char tmpchar4[MAXCHAR];

	if ((fpin = fopen(filename, "r")) == NULL) {
		fprintf(stdout, "Cannot open file %s in improper_id1(), exit\n", filename);
		exit(1);
	}
	readindex = 0;
	for (;;) {
		if (fgets(line, MAXCHAR, fpin) == NULL)
			break;
		sscanf(line, "%s", tmpchar1);
		if (strcmp("IMPROPER", tmpchar1) == 0) {
			readindex = 1;
			continue;
		}
		if (spaceline(line) == 1 && readindex == 1)
			readindex = 0;
		if (strcmp(tmpchar1, "DONE") == 0 && readindex == 1)
			readindex = 0;
		if (strcmp(tmpchar1, "STOP") == 0 && readindex == 1)
			readindex = 0;
		if (strcmp(tmpchar1, "CHARGE") == 0 && readindex == 1)
			readindex = 0;
		if (strcmp(tmpchar1, "LOOP") == 0 && readindex == 1)
			readindex = 0;
		if (strcmp(tmpchar1, "IMPROPER") == 0 && readindex == 1)
			readindex = 0;

		if (readindex == 1) {
			sscanf(line, "%s%s%s%s", tmpchar1, tmpchar2, tmpchar3,
				   tmpchar4);
			if (strcmp(tmpchar1, "-M") == 0)
				continue;
			if (strcmp(tmpchar2, "-M") == 0)
				continue;
			if (strcmp(tmpchar3, "-M") == 0)
				continue;
			if (strcmp(tmpchar4, "-M") == 0)
				continue;
			if (strcmp(tmpchar1, "+M") == 0)
				continue;
			if (strcmp(tmpchar2, "+M") == 0)
				continue;
			if (strcmp(tmpchar3, "+M") == 0)
				continue;
			if (strcmp(tmpchar4, "+M") == 0)
				continue;
			for (i = 0; i < atomnum; i++) {
				if (strcmp(tmpchar1, atom[i].name) == 0) {
					strcpy(tmpchar1, atom[i].ambername);
					improper[impropernum].atid1 = i;
					continue;
				}
				if (strcmp(tmpchar2, atom[i].name) == 0) {
					strcpy(tmpchar2, atom[i].ambername);
					improper[impropernum].atid2 = i;
					continue;
				}
				if (strcmp(tmpchar3, atom[i].name) == 0) {
					strcpy(tmpchar3, atom[i].ambername);
					improper[impropernum].atid3 = i;
					continue;
				}
				if (strcmp(tmpchar4, atom[i].name) == 0) {
					strcpy(tmpchar4, atom[i].ambername);
					improper[impropernum].atid4 = i;
					continue;
				}
			}
			impropernum++;
		}
	}
	fclose(fpin);
}


void improper_id2()
{
	int i;
	improper_prediction();
	impropernum = 0;
	for (i = 0; i < atomnum; i++)
		if (atom[i].improper == 1) {
			improper[impropernum].atid3 = i;
			improper[impropernum].atid1 = atom[i].con[0];
			improper[impropernum].atid2 = atom[i].con[1];
			improper[impropernum].atid4 = atom[i].con[2];
                        if(atom[i].con[0] <0 || atom[i].con[1] <0 || atom[i].con[2] <0)
                                continue;
			impropernum++;
		}
}


void readparm(char *filename)
{
	int mindex = -1;
	int bindex = 0;
	int aindex = 0;
	int tindex = 0;
	int iindex = 0;
	int vindex = 0;
	int num = 0;
	int tmpnum;
	int i, j, k;
	int flag;
	int vdwparmnum_old;
	FILE *fp;
	char line[MAXCHAR];
	char tmpc[MAXCHAR];
	char tmpc1[MAXCHAR];
	char tmpc2[MAXCHAR];
	char tmpc3[MAXCHAR];
	char tmpc4[MAXCHAR];
	if ((fp = fopen(filename, "r")) == NULL) {
		fprintf(stdout, "Cannot open file %s in readparm(), exit\n", filename);
		exit(1);
	}
	for (;;) {
		if (fgets(line, MAXCHAR, fp) == NULL)
			break;
		num++;
		if (mindex == -1 && num == 1) {
			mindex = 1;
			continue;
		}
		if (mindex == 1 && spaceline(line) == 1) {
			mindex = 0;
			num = 0;
			bindex = -1;
			continue;
		}
		if (bindex == -1 && num == 1) {
			bindex = 1;
			continue;
		}
		if (bindex == 1 && spaceline(line) == 1) {
			bindex = 0;
			aindex = 1;
			continue;
		}
		if (aindex == 1 && spaceline(line) == 1) {
			aindex = 0;
			tindex = 1;
			continue;
		}
		if (tindex == 1 && spaceline(line) == 1) {
			tindex = 0;
			iindex = 1;
			continue;
		}
		if (iindex == 1 && spaceline(line) == 1) {
			iindex = 0;
			vindex = -1;
			num = 0;
			continue;
		}
		if (vindex == -1 && num == 2) {
			vindex = 1;
			continue;
		}
		if (vindex == 2 && spaceline(line) == 1) {
			vindex = 0;
			continue;
		}
		if (vindex == 1)
			if (strncmp(line, "MOD4", 4) == 0) {
				vindex = 2;
				continue;
			}
		if (mindex == 1) {
			sscanf(line, "%s%lf%lf", atomtype[atomtypenum].name,
				   &atomtype[atomtypenum].mass,
				   &atomtype[atomtypenum].pol);
			atomtypenum++;
			if (atomtypenum >= maxatomtype) {
				maxatomtype += MAX_FF_ATOMTYPE;
				atomtype =
					(ATOMTYPE *) realloc(atomtype,
										 sizeof(ATOMTYPE) * maxatomtype);
				if (atomtype == NULL) {
					fprintf(stdout,
							"memory allocation error for *atomtype\n");
					exit(1);
				}
			}
		}
		if (bindex == 1) {
			bondparm[bondparmnum].name1[0] = line[0];
			bondparm[bondparmnum].name1[1] = line[1];
			bondparm[bondparmnum].name2[0] = line[3];
			bondparm[bondparmnum].name2[1] = line[4];
			sscanf(&line[5], "%lf%lf", &bondparm[bondparmnum].force,
				   &bondparm[bondparmnum].length);
			bondparmnum++;
			if (bondparmnum >= maxbondparm) {
				maxbondparm += MAX_FF_BOND;
				bondparm =
					(BOND_PARAM *) realloc(bondparm,
										sizeof(BOND_PARAM) * maxbondparm);
				if (bondparm == NULL) {
					fprintf(stdout,
							"memory allocation error for *bondparm\n");
					exit(1);
				}
			}
		}
		if (aindex == 1) {
			angleparm[angleparmnum].name1[0] = line[0];
			angleparm[angleparmnum].name1[1] = line[1];
			angleparm[angleparmnum].name2[0] = line[3];
			angleparm[angleparmnum].name2[1] = line[4];
			angleparm[angleparmnum].name3[0] = line[6];
			angleparm[angleparmnum].name3[1] = line[7];
			sscanf(&line[8], "%lf%lf", &angleparm[angleparmnum].force,
				   &angleparm[angleparmnum].angle);
			angleparmnum++;
			if (angleparmnum >= maxangleparm) {
				maxangleparm += MAX_FF_ANGLE;
				angleparm =
					(ANGLE *) realloc(angleparm,
									  sizeof(ANGLE) * maxangleparm);
				if (angleparm == NULL) {
					fprintf(stdout,
							"memory allocation error for *angleparm\n");
					exit(1);
				}
			}
		}
		if (tindex == 1) {
			torsionparm[torsionparmnum].name1[0] = line[0];
			torsionparm[torsionparmnum].name1[1] = line[1];
			torsionparm[torsionparmnum].name2[0] = line[3];
			torsionparm[torsionparmnum].name2[1] = line[4];
			torsionparm[torsionparmnum].name3[0] = line[6];
			torsionparm[torsionparmnum].name3[1] = line[7];
			torsionparm[torsionparmnum].name4[0] = line[9];
			torsionparm[torsionparmnum].name4[1] = line[10];
			sscanf(&line[11], "%d%lf%lf%lf",
				   &torsionparm[torsionparmnum].mul,
				   &torsionparm[torsionparmnum].force,
				   &torsionparm[torsionparmnum].phase,
				   &torsionparm[torsionparmnum].fterm);
			torsionparmnum++;
			if (torsionparmnum >= maxtorsionparm) {
				maxtorsionparm += MAX_FF_TORSION;
				torsionparm =
					(TORSION *) realloc(torsionparm,
										sizeof(TORSION) * maxtorsionparm);
				if (torsionparm == NULL) {
					fprintf(stdout,
							"memory allocation error for *torsionparm\n");
					exit(1);
				}
			}
		}
		if (iindex == 1) {
			tmpc1[0] = line[0];
			tmpc1[1] = line[1];
                        tmpc1[2] = '\0';
			tmpc2[0] = line[3];
			tmpc2[1] = line[4];
                        tmpc2[2] = '\0';
			tmpc3[0] = line[6];
			tmpc3[1] = line[7];
                        tmpc3[2] = '\0';
			tmpc4[0] = line[9];
			tmpc4[1] = line[10];
                        tmpc4[2] = '\0';
	                if(strcmp(tmpc1, tmpc2) > 0)  {
       	                	strcpy(tmpc, tmpc2);
       	                	strcpy(tmpc2, tmpc1);
                        	strcpy(tmpc1, tmpc);
                	}
                	if(strcmp(tmpc1, tmpc4) > 0)  {
                        	strcpy(tmpc, tmpc4);
                        	strcpy(tmpc4, tmpc1);
                        	strcpy(tmpc1, tmpc);
                	}
                	if(strcmp(tmpc2, tmpc4) > 0)  {
                        	strcpy(tmpc, tmpc4);
                        	strcpy(tmpc4, tmpc2);
                        	strcpy(tmpc2, tmpc);
                	}
                	strcpy(improperparm[improperparmnum].name1, tmpc1);
                	strcpy(improperparm[improperparmnum].name2, tmpc2);
                	strcpy(improperparm[improperparmnum].name3, tmpc3);
                	strcpy(improperparm[improperparmnum].name4, tmpc4);
			sscanf(&line[11], "%lf%lf%lf",
				   &improperparm[improperparmnum].force,
				   &improperparm[improperparmnum].phase,
				   &improperparm[improperparmnum].fterm);
			improperparm[improperparmnum].mul = 1; 
			improperparmnum++;
			if (improperparmnum >= maximproperparm) {
				maximproperparm += MAX_FF_IMPROPER;
				improperparm =
					(IMPROPER *) realloc(improperparm,
										 sizeof(IMPROPER) *
										 maximproperparm);
				if (improperparm == NULL) {
					fprintf(stdout,
							"memory allocation error for *improperparm\n");
					exit(1);
				}
			}
		}
		if (vindex == 1) {
			if (spaceline(line) == 1)
				continue;
			equ_vdw[equ_vdw_num].num = 0;
			tmpnum = 0;
			flag = 1;
			while (flag) {
				flag = 0;
				sscanf(&line[tmpnum], "%s",
					   equ_vdw[equ_vdw_num].name[equ_vdw[equ_vdw_num].
												 num++]);
				if (equ_vdw[equ_vdw_num].num >= MAX_EQU_VDW_NUM) {
					printf
						("\nError: number of equivalent vdw atoms exceeds MAX_EQU_VDW_NUM, exit\n");
					exit(1);
				}
				for (i = tmpnum; i < strlen(line) - 3; i++) {
					if (line[i - 1] != ' ' && i >= 1 && line[i] != ' ')
						continue;
					if (line[i] != ' ' && line[i + 1] != ' ') {
						tmpnum = i + 2;
						flag = 1;
						break;
					}
					if (line[i] != ' ' && line[i + 1] == ' ') {
						tmpnum = i + 1;
						flag = 1;
						break;
					}
				}
			}
			if (equ_vdw[equ_vdw_num].num >= 2)
				equ_vdw_num++;
			if (equ_vdw_num >= MAX_EQU_VDW_NUM) {
				printf
					("\nError: number of equivalent vdw parameters exceeds MAX_EQU_VDW, exit\n");
				exit(1);
			}
		}
		if (vindex == 2) {
			sscanf(line, "%s%lf%lf", vdwparm[vdwparmnum].name,
				   &vdwparm[vdwparmnum].radius, &vdwparm[vdwparmnum].pot);
			vdwparmnum++;
			if (vdwparmnum >= maxvdwparm) {
				maxvdwparm += MAX_FF_VDW;
				vdwparm =
					(VDW *) realloc(vdwparm, sizeof(VDW) * maxvdwparm);
				if (vdwparm == NULL) {
					fprintf(stdout,
							"memory allocation error for *vdwparm\n");
					exit(1);
				}
			}
		}
	}
	fclose(fp);

	if (equ_vdw_num > 0) {
		vdwparmnum_old = vdwparmnum;
		for (i = 0; i < equ_vdw_num; i++)
			for (j = 1; j < equ_vdw[i].num; j++) {
				if (strlen(equ_vdw[i].name[j]) < 1)
					continue;
				for (k = 0; k < vdwparmnum_old; k++)
					if (strcmp(vdwparm[k].name, equ_vdw[i].name[0]) == 0) {
						strcpy(vdwparm[vdwparmnum].name,
							   equ_vdw[i].name[j]);
						vdwparm[vdwparmnum].radius = vdwparm[k].radius;
						vdwparm[vdwparmnum].pot = vdwparm[k].pot;
						vdwparmnum++;
						break;
					}
			}
	}

/*
for(i=0;i<atomtypenum;i++)
 printf("\n%s %9.4lf %9.4lf", atomtype[i].name, atomtype[i].mass,atomtype[i].pol);

for(i=0;i<bondparmnum;i++)
 printf("\n%s %s %9.4lf %9.4lf", bondparm[i].name1, bondparm[i].name2, 
        bondparm[i].force, bondparm[i].length);

for(i=0;i<angleparmnum;i++)
 printf("\n%s %s %s %9.4lf %9.4lf", angleparm[i].name1, angleparm[i].name2, 
        angleparm[i].name3, angleparm[i].force, angleparm[i].angle);

for(i=0;i<torsionparmnum;i++)
 printf("\n%s %s %s %s %9.4lf %9.4lf %9.4lf", torsionparm[i].name1, 
        torsionparm[i].name2, torsionparm[i].name3,torsionparm[i].name4,
        torsionparm[i].phase, torsionparm[i].force, torsionparm[i].fterm);

for(i=0;i<improperparmnum;i++)
 printf("\n%s %s %s %s %9.4lf %9.4lf %9.4lf", improperparm[i].name1, 
        improperparm[i].name2, improperparm[i].name3,improperparm[i].name4,
        improperparm[i].phase, improperparm[i].force, improperparm[i].fterm);
for(i=0;i<vdwparmnum;i++)
 printf("\n%s %9.4lf %9.4lf", vdwparm[i].name, vdwparm[i].radius, vdwparm[i].pot);
*/
}

void readcorr(char *filename)
{
	FILE *fp;
	char line[MAXCHAR];
	char tmpchar[6][MAXCHAR];
	int num;

	if ((fp = fopen(filename, "r")) == NULL) {
		fprintf(stdout, "Cannot open file %s in readcorr(), exit\n", filename);
		exit(1);
	}

	corrnum = 0;

	for (;;) {
		if (fgets(line, MAXCHAR, fp) == NULL)
			break;
		if (strncmp("CORR", &line[0], 4) == 0) {
			for (i = 0; i < 6; i++)
				tmpchar[i][0] = '\0';
			num = 0;
			for (i = 3; i < strlen(line) - 2; i++) {
				if (line[i + 1] == '\0')
					break;
				if (line[i] == ' ' && line[i + 1] != ' ') {
					sscanf(&line[i], "%s", tmpchar[num]);
					num++;
				}
			}
			strcpy(corr[corrnum].name1, tmpchar[0]);
			improperindex[corrnum] = atoi(tmpchar[1]);
			if (num == 3) {
				strcpy(corr[corrnum].name2, tmpchar[2]);
				corr[corrnum].num = 1;
			}
			if (num == 4) {
				strcpy(corr[corrnum].name2, tmpchar[2]);
				strcpy(corr[corrnum].name3, tmpchar[3]);
				corr[corrnum].num = 2;
			}
			if (num >= 5) {
				strcpy(corr[corrnum].name2, tmpchar[2]);
				strcpy(corr[corrnum].name3, tmpchar[3]);
				strcpy(corr[corrnum].name4, tmpchar[4]);
				corr[corrnum].num = 3;
			}
			corrnum++;

			if (corrnum >= maxcorr) {
				maxcorr += MAX_FF_CORR;
				corr = (CORR *) realloc(corr, sizeof(CORR) * maxcorr);
				if (corr == NULL) {
					fprintf(stdout, "memory allocation error for *corr\n");
					exit(1);
				}
				improperindex =
					(int *) realloc(improperindex, sizeof(int) * maxcorr);
				if (improperindex == NULL) {
					fprintf(stdout,
							"memory allocation error for *improperindex\n");
					exit(1);
				}
			}
		}
	}
	fclose(fp);
/*
for(i=0;i<corrnum;i++)
printf("\n%5d %5d %9s %9s %9s %9s %9d", i+1, corr[i].num, corr[i].name1, corr[i].name2, corr[i].name3,corr[i].name4, improperindex[i]);
exit(1);
*/
}

void prepare(void)
{
	int i, j;
	for (i = 0; i < atomnum; i++) {
		strcpy(corrname[i].name, "");
		strcpy(similarname[i].name, "");
		strcpy(similarname2[i].name, "");
		corrindex[i] = 0;
		similarindex[i] = 0;
		similarindex2[i] = 0;
		for (j = 0; j < corrnum; j++)
			if (strcmp(corr[j].name1, atom[i].ambername) == 0) {
				corrindex[i] = 1;
				strcpy(corrname[i].name, corr[j].name2);
				if (corr[j].num == 2) {
					similarindex[i] = 1;
					strcpy(similarname[i].name, corr[j].name3);
				}
				if (corr[j].num == 3) {
					similarindex[i] = 1;
					strcpy(similarname[i].name, corr[j].name3);
					similarindex2[i] = 1;
					strcpy(similarname2[i].name, corr[j].name4);
				}
				break;
			}
	}
	/*
	   for(i =0; i<atomnum; i++)
	   printf("%5d%5s%5s%5s%5s%5d%5d%5d\n", i+1, atom[i].ambername, corrname[i].name, similarname[i].name,similarname2[i].name,
	   corrindex[i], similarindex[i], similarindex2[i]);
	   exit(1);
	 */
}

int empangle(char *tmpc1, char *tmpc2, char *tmpc3, char *name1,
			 char *name2, char *name3, int id1, int id2, int id3)
{
	int num1 = -1, num2 = -1;
	double bondlength1 = 0.0;
	double bondlength2 = 0.0;
	double cparm, dparm, zparm1, zparm2;
	double angle;
	double force;

	if (tmpc1[1] == '\0') {
		tmpc1[1] = ' ';
		tmpc1[2] = '\0';
	}
	if (tmpc2[1] == '\0') {
		tmpc2[1] = ' ';
		tmpc2[2] = '\0';
	}
	if (tmpc3[1] == '\0') {
		tmpc3[1] = ' ';
		tmpc3[2] = '\0';
	}
	for (i = 0; i < angleparmnum; i++)
		if (angleparm[i].name1[0] == tmpc1[0]
			&& angleparm[i].name1[1] == tmpc1[1]
			&& angleparm[i].name2[0] == tmpc2[0]
			&& angleparm[i].name2[1] == tmpc2[1]
			&& angleparm[i].name3[0] == tmpc1[0]
			&& angleparm[i].name3[1] == tmpc1[1]) {
			num1 = i;
			break;
		}
	if (num1 == -1)
		return 0;
	for (i = 0; i < angleparmnum; i++)
		if (angleparm[i].name1[0] == tmpc3[0]
			&& angleparm[i].name1[1] == tmpc3[1]
			&& angleparm[i].name2[0] == tmpc2[0]
			&& angleparm[i].name2[1] == tmpc2[1]
			&& angleparm[i].name3[0] == tmpc3[0]
			&& angleparm[i].name3[1] == tmpc3[1]) {
			num2 = i;
			break;
		}
	if (num2 == -1)
		return 0;

	angle = 0.5 * (angleparm[num1].angle + angleparm[num2].angle);
	for (i = 0; i < bondparmnum; i++)
		if ((bondparm[i].name1[0] == tmpc1[0]
			 && bondparm[i].name1[1] == tmpc1[1]
			 && bondparm[i].name2[0] == tmpc2[0]
			 && bondparm[i].name2[1] == tmpc2[1])
			|| (bondparm[i].name2[0] == tmpc1[0]
				&& bondparm[i].name2[1] == tmpc1[1]
				&& bondparm[i].name1[0] == tmpc2[0]
				&& bondparm[i].name1[1] == tmpc2[1])) {
			bondlength1 = bondparm[i].length;
			break;
		}
	if (bondlength1 == 0.0)
		return 0;

	for (i = 0; i < bondparmnum; i++)
		if ((bondparm[i].name1[0] == tmpc2[0]
			 && bondparm[i].name1[1] == tmpc2[1]
			 && bondparm[i].name2[0] == tmpc3[0]
			 && bondparm[i].name2[1] == tmpc3[1])
			|| (bondparm[i].name2[0] == tmpc2[0]
				&& bondparm[i].name2[1] == tmpc2[1]
				&& bondparm[i].name1[0] == tmpc3[0]
				&& bondparm[i].name1[1] == tmpc3[1])) {
			bondlength2 = bondparm[i].length;
			break;
		}
	if (bondlength2 == 0.0)
		return 0;

	/* calculate the bond angle force  */
	zparm1 = 0.0;
	if (id1 == 1)
		zparm1 = 0.784;
	if (id1 == 6)
		zparm1 = 1.183;
	if (id1 == 7)
		zparm1 = 1.212;
	if (id1 == 8)
		zparm1 = 1.219;
	if (id1 == 9)
		zparm1 = 1.166;
	if (id1 == 17)
		zparm1 = 1.272;
	if (id1 == 35)
		zparm1 = 1.378;
	if (id1 == 53)
		zparm1 = 1.398;
	if (id1 == 15)
		zparm1 = 1.620;
	if (id1 == 16)
		zparm1 = 1.280;

	zparm2 = 0.0;
	if (id3 == 1)
		zparm2 = 0.784;
	if (id3 == 6)
		zparm2 = 1.183;
	if (id3 == 7)
		zparm2 = 1.212;
	if (id3 == 8)
		zparm2 = 1.219;
	if (id3 == 9)
		zparm2 = 1.166;
	if (id3 == 17)
		zparm2 = 1.272;
	if (id3 == 35)
		zparm2 = 1.378;
	if (id3 == 53)
		zparm2 = 1.398;
	if (id3 == 15)
		zparm2 = 1.620;
	if (id3 == 16)
		zparm2 = 1.280;

	cparm = 0.0;
	if (id2 == 1)
		cparm = 0.0;
	if (id2 == 6)
		cparm = 1.339;
	if (id2 == 7)
		cparm = 1.300;
	if (id2 == 8)
		cparm = 1.249;
	if (id2 == 9)
		cparm = 0.0;
	if (id2 == 17)
		cparm = 0.0;
	if (id2 == 35)
		cparm = 0.0;
	if (id2 == 53)
		cparm = 0.0;
	if (id2 == 15)
		cparm = 1.448;
	if (id2 == 16)
		cparm = 0.906;

	dparm = (bondlength1 - bondlength2) * (bondlength1 - bondlength2);
	dparm =
		dparm / (bondlength1 + bondlength2) * (bondlength1 + bondlength2);
	force =
		143.9 * zparm1 * cparm * zparm2 * exp(-2 * dparm) / (bondlength1 +
															 bondlength2);
	force /= sqrt(angle * 3.1415926 / 180.0);

	fprintf(ffout, "%2s %2s %2s%9.3lf%12.3lf",
		name1, name2, name3, force, angle);
        fprintf(ffout, "   ! Calculated with empirical approach\n");

        angleparm[angleparmnum].flag = 1;
	strcpy(angleparm[angleparmnum].name1, name1);
	strcpy(angleparm[angleparmnum].name2, name2);
	strcpy(angleparm[angleparmnum].name3, name3);
	angleparm[angleparmnum].angle = angle;
	angleparm[angleparmnum].force = force;
	angleparmnum++;
	if (angleparmnum >= maxangleparm) {
		maxangleparm += MAX_FF_ANGLE;
		angleparm =
			(ANGLE *) realloc(angleparm, sizeof(ANGLE) * maxangleparm);
		if (angleparm == NULL) {
			fprintf(stdout, "memory allocation error for *angleparm\n");
			exit(1);
		}
	}
	return 1;
}

void wrt_atomtype(void)
{
	int i, j;
	int suc;
        int count = 1;
	char tmpc[5];
	for (i = 0; i < atomnum; i++) {
		suc = 0;
		for (j = 0; j < atomtypenum; j++)
			if (strcmp(atom[i].ambername, atomtype[j].name) == 0) {
				if (atomtype[j].flag == 0) {
					fprintf(fpout, "MASS %5d %-5s %10.6lf\n",
						count++, atomtype[j].name, atomtype[j].mass);
					atomtype[j].flag = 1;
				}
				suc = 1;
				break;
			}
		if (suc == 0 && corrindex[i] == 1) {
			strcpy(tmpc, corrname[i].name);
			/*  if(tmpc[1]=='\0') tmpc[1]=' '; */
			for (j = 0; j < atomtypenum; j++)
				if (strcmp(tmpc, atomtype[j].name) == 0) {
					suc = 1;
				        fprintf(fpout, "MASS %5d %-5s %10.6lf\n",
				              count++, atomtype[j].name, atomtype[j].mass);
					fprintf(fpout, "                          %s %-3s\n", "same as",
							atomtype[j].name);
					atomtype[atomtypenum] = atomtype[j];
					strcpy(atomtype[atomtypenum].name, atom[i].ambername);
					atomtypenum++;
					if (atomtypenum >= maxatomtype) {
						maxatomtype += MAX_FF_ATOMTYPE;
						atomtype =
							(ATOMTYPE *) realloc(atomtype,
												 sizeof(ATOMTYPE) *
												 maxatomtype);
						if (atomtype == NULL) {
							fprintf(stdout,
									"memory allocation error for *atomtype\n");
							exit(1);
						}
					}
					break;
				}
		}
		if (suc == 0) {
		        fprintf(fpout, "MASS %5d %-5s %4.5f\n",
		              count++, atomtype[j].name, atomtype[j].mass);
			fprintf(fpout, "                          %5s\n", "ATTN, need revision");
			strcpy(atomtype[atomtypenum].name, atom[i].ambername);
			atomtypenum++;
			if (atomtypenum >= maxatomtype) {
				maxatomtype += MAX_FF_ATOMTYPE;
				atomtype =
					(ATOMTYPE *) realloc(atomtype,
										 sizeof(ATOMTYPE) * maxatomtype);
				if (atomtype == NULL) {
					fprintf(stdout,
							"memory allocation error for *atomtype\n");
					exit(1);
				}
			}
		}
	}
}

int vdw(char *at_name, char *corr_name)
{
	int suc = 0;
	int j;
	for (j = 0; j < vdwparmnum; j++)
		if (strcmp(vdwparm[j].name, corr_name) == 0) {
			suc = 1;
                        if (vdwparm[j].flag == 0) {
			  fprintf(ffout, "%2s%10.2lf%10.4lf%10.4lf%10.2lf%10.4lf%10.4lf\n",
			   	  at_name, 0.0, -vdwparm[j].pot  , vdwparm[j].radius,
				  0.0, -vdwparm[j].pot/2, vdwparm[j].radius);
			  fprintf(ffout, "      ! %s %-3s\n", "same as",
				  vdwparm[j].name);
			  vdwparm[j].flag = 1;
			}
			vdwparm[vdwparmnum] = vdwparm[j];
			strcpy(vdwparm[vdwparmnum].name, at_name);
			vdwparmnum++;
			if (vdwparmnum >= maxvdwparm) {
				maxvdwparm += MAX_FF_VDW;
				vdwparm =
					(VDW *) realloc(vdwparm, sizeof(VDW) * maxvdwparm);
				if (vdwparm == NULL) {
					fprintf(stdout,
							"memory allocation error for *vdwparm\n");
					exit(1);
				}
			}
			break;
		}
	return suc;
}

void wrt_vdw(void)
{
	int i, j;
	int suc;
	fprintf(ffout, "\n%s\n", "NONBONDED  NBXMOD 5  GROUP SWITCH CDIEL -");
	fprintf(ffout, "CUTNB 14.0  CTOFNB 12.0  CTONNB 10.0  EPS 1.0  E14FAC 0.83333333  WMIN 1.4\n");
	fprintf(ffout, "!                Emin     Rmin/2              Emin/2     Rmin  (for 1-4's)\n");
	fprintf(ffout, "!             (kcal/mol)    (A)\n");
	for (i = 0; i < atomnum; i++) {
		suc = 0;
		for (j = 0; j < vdwparmnum; j++)
			if (strcmp(vdwparm[j].name, atom[i].ambername) == 0) {
				if (vdwparm[j].flag == 0) {
					fprintf(ffout, "%2s%10.2lf%10.4lf%10.4lf%10.2lf%10.4lf%10.4lf\n",
						atom[i].ambername, 0.0, -vdwparm[j].pot  , vdwparm[j].radius,
								      0.0, -vdwparm[j].pot/2, vdwparm[j].radius);
					vdwparm[j].flag = 1;
				}
                                suc = 1;
				break;
			}
		if (suc == 0 && similarindex[i] == 1)
			suc = vdw(atom[i].ambername, similarname[i].name);

		if (suc == 0 && similarindex2[i] == 1)
			suc = vdw(atom[i].ambername, similarname2[i].name);

		if (suc == 0 && corrindex[i] == 1)
			suc = vdw(atom[i].ambername, corrname[i].name);
		if (suc == 0) {
			fprintf(ffout, "%2s%10.2lf%10.4lf%10.4lf%10.2lf%10.4lf%10.4lf\n",
				atom[i].ambername, 0.0, 0.0, 0.0,
				  	           0.0, 0.0, 0.0);
			fprintf(fpout, "     ! %s\n", "ATTN, need revision");
			strcpy(vdwparm[vdwparmnum].name, atom[i].ambername);
			vdwparmnum++;
			if (vdwparmnum >= maxvdwparm) {
				maxvdwparm += MAX_FF_VDW;
				vdwparm =
					(VDW *) realloc(vdwparm, sizeof(VDW) * maxvdwparm);
				if (vdwparm == NULL) {
					fprintf(stdout,
							"memory allocation error for *vdwparm\n");
					exit(1);
				}
			}
		}
	}
	fprintf(fpout, "\n\n\n");
}

int bond(char *at_name1, char *at_name2, char *corr_name1,
         char *corr_name2, int index)
{
	int suc = 0;
	int k;
	for (k = 0; k < bondparmnum; k++)
		if ((bondparm[k].name1[0] == corr_name1[0]
			 && bondparm[k].name1[1] == corr_name1[1]
			 && bondparm[k].name2[0] == corr_name2[0]
			 && bondparm[k].name2[1] == corr_name2[1])
			|| (bondparm[k].name2[0] == corr_name1[0]
				&& bondparm[k].name2[1] == corr_name1[1]
				&& bondparm[k].name1[0] == corr_name2[0]
				&& bondparm[k].name1[1] == corr_name2[1])) {
			suc = 1;
                        if (index == 1) {
                                fprintf(ffout, "%2s %2s%8.2lf%8.3lf", at_name1,
                                        at_name2, bondparm[k].force, bondparm[k].length);
                                fprintf(ffout, "	! %s  %2s %2s\n",
                                        "same as", bondparm[k].name1, bondparm[k].name2);
                                bondparm[bondparmnum].flag = 1;
 				bondparm[bondparmnum] = bondparm[k];
				strcpy(bondparm[bondparmnum].name1, at_name1);
				strcpy(bondparm[bondparmnum].name2, at_name2);
				bondparmnum++;
				if (bondparmnum >= maxbondparm) {
					maxbondparm += MAX_FF_BOND;
					bondparm =
						(BOND_PARAM *) realloc(bondparm,
											sizeof(BOND_PARAM) * maxbondparm);
					if (bondparm == NULL) {
						fprintf(stdout,
								"memory allocation error for *bondparm\n");
						exit(1);
					}
				}
				break;
			} else if (bondparm[k].flag == 0) {
  		          fprintf(ffout, "%2s %2s%8.2lf%8.3lf\n", bondparm[k].name1,
			          bondparm[k].name2, bondparm[k].force, bondparm[k].length);
			  bondparm[k].flag = 1;
			}
		}
	return suc;
}

void wrt_bond(void)
{
	int i, j;
	int suc;
        float bindv;
	char tmpc1[5], tmpc2[5];
	char name1[5], name2[5];

        fprintf(fpout, "\n");
        fprintf(ffout, "\n%s\n", "BOND");

	for (i = 0; i < atomnum; i++)
		for (j = i + 1; j < atomnum; j++)
			if (atom[i].con[0] == j || atom[i].con[1] == j
				|| atom[i].con[2] == j || atom[i].con[3] == j
				|| atom[i].con[4] == j || atom[i].con[5] == j) {
				suc = 0;
				strcpy(tmpc1, atom[i].ambername);
				strcpy(tmpc2, atom[j].ambername);
				if (tmpc1[1] == '\0') {
					tmpc1[1] = ' ';
					tmpc1[2] = '\0';
				}
				if (tmpc2[1] == '\0') {
					tmpc2[1] = ' ';
					tmpc2[2] = '\0';
				}
				strcpy(name1, tmpc1);
				strcpy(name2, tmpc2);

                                bindv =
  		                (atom[i].x - atom[j].x) *
			        (atom[i].x - atom[j].x);
			        bindv +=
                                (atom[i].y - atom[j].y) *
			        (atom[i].y - atom[j].y);
			        bindv +=
			        (atom[i].z - atom[j].z) *
			        (atom[i].z - atom[j].z);
                                bindv = sqrt(bindv);

                                fprintf(fpout, "BOND %-5s %-5s", atom[i].name, atom[j].name);
                                fprintf(fpout, "        	! dist	%10.4lf\n", bindv);

				suc = bond(name1, name2, tmpc1, tmpc2, 0);

				if (suc == 0 && similarindex[i] == 1) {
					strcpy(tmpc1, similarname[i].name);
					if (tmpc1[1] == '\0') {
						tmpc1[1] = ' ';
						tmpc1[2] = '\0';
					}
					if (tmpc2[1] == '\0') {
						tmpc2[1] = ' ';
						tmpc2[2] = '\0';
					}
					suc = bond(name1, name2, tmpc1, tmpc2, 1);
				}
				if (suc == 0 && similarindex[j] == 1) {
					strcpy(tmpc1, atom[i].ambername);
					strcpy(tmpc2, similarname[j].name);
					if (tmpc1[1] == '\0') {
						tmpc1[1] = ' ';
						tmpc1[2] = '\0';
					}
					if (tmpc2[1] == '\0') {
						tmpc2[1] = ' ';
						tmpc2[2] = '\0';
					}
					suc = bond(name1, name2, tmpc1, tmpc2, 1);
				}
				if (suc == 0 && similarindex[i] == 1
					&& similarindex[j] == 1) {
					strcpy(tmpc1, similarname[i].name);
					strcpy(tmpc2, similarname[j].name);
					if (tmpc1[1] == '\0') {
						tmpc1[1] = ' ';
						tmpc1[2] = '\0';
					}
					if (tmpc2[1] == '\0') {
						tmpc2[1] = ' ';
						tmpc2[2] = '\0';
					}
					suc = bond(name1, name2, tmpc1, tmpc2, 1);
				}

				if (suc == 0) {
					strcpy(tmpc1, atom[i].ambername);
					strcpy(tmpc2, atom[j].ambername);
					if (tmpc1[1] == '\0') {
						tmpc1[1] = ' ';
						tmpc1[2] = '\0';
					}
					if (tmpc2[1] == '\0') {
						tmpc2[1] = ' ';
						tmpc2[2] = '\0';
					}
					strcpy(name1, tmpc1);
					strcpy(name2, tmpc2);
				}
				if (suc == 0 && similarindex2[i] == 1) {
					strcpy(tmpc1, similarname2[i].name);
					if (tmpc1[1] == '\0') {
						tmpc1[1] = ' ';
						tmpc1[2] = '\0';
					}
					if (tmpc2[1] == '\0') {
						tmpc2[1] = ' ';
						tmpc2[2] = '\0';
					}
					suc = bond(name1, name2, tmpc1, tmpc2, 1);
				}
				if (suc == 0 && similarindex2[j] == 1) {
					strcpy(tmpc1, atom[i].ambername);
					strcpy(tmpc2, similarname2[j].name);
					if (tmpc1[1] == '\0') {
						tmpc1[1] = ' ';
						tmpc1[2] = '\0';
					}
					if (tmpc2[1] == '\0') {
						tmpc2[1] = ' ';
						tmpc2[2] = '\0';
					}
					suc = bond(name1, name2, tmpc1, tmpc2, 1);
				}
				/* combination of similar and similar2 */
				if (suc == 0 && similarindex[i] == 1
					&& similarindex2[j] == 1) {
					strcpy(tmpc1, similarname[i].name);
					strcpy(tmpc2, similarname2[j].name);
					if (tmpc1[1] == '\0') {
						tmpc1[1] = ' ';
						tmpc1[2] = '\0';
					}
					if (tmpc2[1] == '\0') {
						tmpc2[1] = ' ';
						tmpc2[2] = '\0';
					}
					suc = bond(name1, name2, tmpc1, tmpc2, 1);
				}
				if (suc == 0 && similarindex2[i] == 1
					&& similarindex[j] == 1) {
					strcpy(tmpc1, similarname2[i].name);
					strcpy(tmpc2, similarname[j].name);
					if (tmpc1[1] == '\0') {
						tmpc1[1] = ' ';
						tmpc1[2] = '\0';
					}
					if (tmpc2[1] == '\0') {
						tmpc2[1] = ' ';
						tmpc2[2] = '\0';
					}
					suc = bond(name1, name2, tmpc1, tmpc2, 1);
				}
				if (suc == 0 && similarindex2[i] == 1
					&& similarindex2[j] == 1) {
					strcpy(tmpc1, similarname2[i].name);
					strcpy(tmpc2, similarname2[j].name);
					if (tmpc1[1] == '\0') {
						tmpc1[1] = ' ';
						tmpc1[2] = '\0';
					}
					if (tmpc2[1] == '\0') {
						tmpc2[1] = ' ';
						tmpc2[2] = '\0';
					}
					suc = bond(name1, name2, tmpc1, tmpc2, 1);
				}
				/* end */
				if (suc == 0) {
					strcpy(tmpc1, atom[i].ambername);
					strcpy(tmpc2, atom[j].ambername);
					if (tmpc1[1] == '\0') {
						tmpc1[1] = ' ';
						tmpc1[2] = '\0';
					}
					if (tmpc2[1] == '\0') {
						tmpc2[1] = ' ';
						tmpc2[2] = '\0';
					}
					strcpy(name1, tmpc1);
					strcpy(name2, tmpc2);
				}

				if (suc == 0 && corrindex[i] == 1) {
					strcpy(tmpc1, corrname[i].name);
					strcpy(tmpc2, atom[j].ambername);
					if (tmpc1[1] == '\0') {
						tmpc1[1] = ' ';
						tmpc1[2] = '\0';
					}
					if (tmpc2[1] == '\0') {
						tmpc2[1] = ' ';
						tmpc2[2] = '\0';
					}
					suc = bond(name1, name2, tmpc1, tmpc2, 1);
				}
				if (suc == 0 && corrindex[j] == 1) {
					strcpy(tmpc1, atom[i].ambername);
					strcpy(tmpc2, corrname[j].name);
					if (tmpc1[1] == '\0') {
						tmpc1[1] = ' ';
						tmpc1[2] = '\0';
					}
					if (tmpc2[1] == '\0') {
						tmpc2[1] = ' ';
						tmpc2[2] = '\0';
					}
					suc = bond(name1, name2, tmpc1, tmpc2, 1);
				}
				if (suc == 0 && corrindex[i] == 1 && corrindex[j] == 1) {
					strcpy(tmpc1, corrname[i].name);
					strcpy(tmpc2, corrname[j].name);
					if (tmpc1[1] == '\0') {
						tmpc1[1] = ' ';
						tmpc1[2] = '\0';
					}
					if (tmpc2[1] == '\0') {
						tmpc2[1] = ' ';
						tmpc2[2] = '\0';
					}
					suc = bond(name1, name2, tmpc1, tmpc2, 1);
				}
				if (suc == 0) {
					if (bondparm[k].flag == 0) {
						fprintf(ffout, "%2s %2s%8.2lf%8.3lf", name1, name2, 0.0, 0.0);
						fprintf(ffout, "	!%s\n", "ATTN, need revision");
						bondparm[k].flag = 1;
					}
					strcpy(bondparm[bondparmnum].name1, name1);
					strcpy(bondparm[bondparmnum].name2, name2);
					bondparmnum++;
					if (bondparmnum >= maxbondparm) {
						maxbondparm += MAX_FF_BOND;
						bondparm =
							(BOND_PARAM *) realloc(bondparm,
												sizeof(BOND_PARAM) *
												maxbondparm);
						if (bondparm == NULL) {
							fprintf(stdout,
									"memory allocation error for *bondparm\n");
							exit(1);
						}
					}
				}
			}
}

int angle(char *at_name1, char *at_name2, char *at_name3, char *corr_name1,
          char *corr_name2, char *corr_name3, int index)
{
	int suc = 0;
	int l;

	for (l = 0; l < angleparmnum; l++)
		if ((angleparm[l].name1[0] == corr_name1[0]
			 && angleparm[l].name1[1] == corr_name1[1]
			 && angleparm[l].name2[0] == corr_name2[0]
			 && angleparm[l].name2[1] == corr_name2[1]
			 && angleparm[l].name3[0] == corr_name3[0]
			 && angleparm[l].name3[1] == corr_name3[1])
			|| (angleparm[l].name3[0] == corr_name1[0]
				&& angleparm[l].name3[1] == corr_name1[1]
				&& angleparm[l].name2[0] == corr_name2[0]
				&& angleparm[l].name2[1] == corr_name2[1]
				&& angleparm[l].name1[0] == corr_name3[0]
				&& angleparm[l].name1[1] == corr_name3[1])) {
			suc = 1;
			if (index == 1) {
				fprintf(ffout, "%2s %2s %2s%9.3lf%12.3lf",
					at_name1, at_name2, at_name3,
     					angleparm[l].force,
					angleparm[l].angle);
                		fprintf(ffout, "	! %s %2s %2s %2s\n",
					"same as", angleparm[l].name1,
					angleparm[l].name2,
					angleparm[l].name3);
				angleparm[angleparmnum] = angleparm[l];
				strcpy(angleparm[angleparmnum].name1, at_name1);
				strcpy(angleparm[angleparmnum].name2, at_name2);
				strcpy(angleparm[angleparmnum].name3, at_name3);
				angleparmnum++;
				if (angleparmnum >= maxangleparm) {
					maxangleparm += MAX_FF_ANGLE;
					angleparm =
						(ANGLE *) realloc(angleparm,
										  sizeof(ANGLE) * maxangleparm);
					if (angleparm == NULL) {
						fprintf(stdout,
								"memory allocation error for *angleparm\n");
						exit(1);
					}
				}
				break;
			} else if (angleparm[l].flag == 0) {
                 	  fprintf(ffout, "%2s %2s %2s%9.3lf%12.3lf\n",
				  at_name1, at_name2, at_name3,
				  angleparm[l].force,
				  angleparm[l].angle);
			  angleparm[l].flag = 1;
			}
		}
	return suc;
}


void wrt_angle(void)
{

	int i, j, k;
	int suc;
	float anglev;
	char tmpc1[5], tmpc2[5], tmpc3[5];
	char name1[5], name2[5], name3[5];

        fprintf(fpout, "\n");
        fprintf(ffout, "\n%s\n", "ANGLE");

	/* NB: non-standard indentation in next four lines; for readability  */
	for (i = 0; i < atomnum; i++) {
		for (j = 0; j < atomnum; j++) {
			for (k = i+1; k < atomnum; k++) {
					if (atom[i].con[0] == j || atom[i].con[1] == j
						|| atom[i].con[2] == j || atom[i].con[3] == j
						|| atom[i].con[4] == j || atom[i].con[5] == j) {
						if (atom[j].con[0] == k || atom[j].con[1] == k
							|| atom[j].con[2] == k || atom[j].con[3] == k
							|| atom[j].con[4] == k
							|| atom[j].con[5] == k) {
							suc = 0;
							strcpy(tmpc1, atom[i].ambername);
							strcpy(tmpc2, atom[j].ambername);
							strcpy(tmpc3, atom[k].ambername);
							if (tmpc1[1] == '\0') {
								tmpc1[1] = ' ';
								tmpc1[2] = '\0';
							}
							if (tmpc2[1] == '\0') {
								tmpc2[1] = ' ';
								tmpc2[2] = '\0';
							}
							if (tmpc3[1] == '\0') {
								tmpc3[1] = ' ';
								tmpc3[2] = '\0';
							}
							strcpy(name1, tmpc1);
							strcpy(name2, tmpc2);
							strcpy(name3, tmpc3);
                                                        anglev = anglecal(atom[i], atom[j], atom[k]);
                                    	                fprintf(fpout, "ANGL %-5s %-5s %-5s ",
                                                                atom[i].name, atom[j].name, atom[k].name);
   		                                        fprintf(fpout, "        	! angle	%10.4lf\n", anglev);
							suc =
								angle(name1, name2, name3, tmpc1, tmpc2,
									  tmpc3, 0);
							if (suc == 0 && similarindex[i] == 1) {
								strcpy(tmpc1, similarname[i].name);
								strcpy(tmpc2, atom[j].ambername);
								strcpy(tmpc3, atom[k].ambername);
								if (tmpc1[1] == '\0') {
									tmpc1[1] = ' ';
									tmpc1[2] = '\0';
								}

								if (tmpc2[1] == '\0') {
									tmpc2[1] = ' ';
									tmpc2[2] = '\0';
								}
								if (tmpc3[1] == '\0') {
									tmpc3[1] = ' ';
									tmpc3[2] = '\0';
								}
								suc =
									angle(name1, name2, name3, tmpc1,
										  tmpc2, tmpc3, 1);
							}

							if (suc == 0 && similarindex[j] == 1) {
								strcpy(tmpc1, atom[i].ambername);
								strcpy(tmpc2, similarname[j].name);
								strcpy(tmpc3, atom[k].ambername);
								if (tmpc1[1] == '\0') {
									tmpc1[1] = ' ';
									tmpc1[2] = '\0';
								}
								if (tmpc2[1] == '\0') {
									tmpc2[1] = ' ';
									tmpc2[2] = '\0';
								}
								if (tmpc3[1] == '\0') {
									tmpc3[1] = ' ';
									tmpc3[2] = '\0';
								}
								suc =
									angle(name1, name2, name3, tmpc1,
										  tmpc2, tmpc3, 1);
							}

							if (suc == 0 && similarindex[k] == 1) {
								strcpy(tmpc1, atom[i].ambername);
								strcpy(tmpc2, atom[j].ambername);
								strcpy(tmpc3, similarname[k].name);
								if (tmpc1[1] == '\0') {
									tmpc1[1] = ' ';
									tmpc1[2] = '\0';
								}
								if (tmpc2[1] == '\0') {
									tmpc2[1] = ' ';
									tmpc2[2] = '\0';
								}
								if (tmpc3[1] == '\0') {
									tmpc3[1] = ' ';
									tmpc3[2] = '\0';
								}
								suc =
									angle(name1, name2, name3, tmpc1,
										  tmpc2, tmpc3, 1);
							}


							/* from here estimate the bond angle parameter using empirical method */
							if (suc == 0) {
								strcpy(tmpc1, atom[i].ambername);
								strcpy(tmpc2, atom[j].ambername);
								strcpy(tmpc3, atom[k].ambername);
								if (tmpc1[1] == '\0') {
									tmpc1[1] = ' ';
									tmpc1[2] = '\0';
								}
								if (tmpc2[1] == '\0') {
									tmpc2[1] = ' ';
									tmpc2[2] = '\0';
								}
								if (tmpc3[1] == '\0') {
									tmpc3[1] = ' ';
									tmpc3[2] = '\0';
								}
								suc =
									empangle(tmpc1, tmpc2, tmpc3, name1,
											 name2, name3,
											 atom[i].atomicnum,
											 atom[j].atomicnum,
											 atom[k].atomicnum);
							}
							if (suc == 0) {
								strcpy(tmpc1, similarname[i].name);
								strcpy(tmpc2, atom[j].ambername);
								strcpy(tmpc3, atom[k].ambername);
								if (tmpc1[1] == '\0') {
									tmpc1[1] = ' ';
									tmpc1[2] = '\0';
								}
								if (tmpc2[1] == '\0') {
									tmpc2[1] = ' ';
									tmpc2[2] = '\0';
								}
								if (tmpc3[1] == '\0') {
									tmpc3[1] = ' ';
									tmpc3[2] = '\0';
								}
								suc =
									empangle(tmpc1, tmpc2, tmpc3, name1,
											 name2, name3,
											 atom[i].atomicnum,
											 atom[j].atomicnum,
											 atom[k].atomicnum);
							}
							if (suc == 0) {
								strcpy(tmpc1, atom[i].ambername);
								strcpy(tmpc2, similarname[j].name);
								strcpy(tmpc3, atom[k].ambername);
								if (tmpc1[1] == '\0') {
									tmpc1[1] = ' ';
									tmpc1[2] = '\0';
								}
								if (tmpc2[1] == '\0') {
									tmpc2[1] = ' ';
									tmpc2[2] = '\0';
								}
								if (tmpc3[1] == '\0') {
									tmpc3[1] = ' ';
									tmpc3[2] = '\0';
								}
								suc =
									empangle(tmpc1, tmpc2, tmpc3, name1,
											 name2, name3,
											 atom[i].atomicnum,
											 atom[j].atomicnum,
											 atom[k].atomicnum);
							}
							if (suc == 0) {
								strcpy(tmpc1, atom[i].ambername);
								strcpy(tmpc2, atom[j].ambername);
								strcpy(tmpc3, similarname[k].name);
								if (tmpc1[1] == '\0') {
									tmpc1[1] = ' ';
									tmpc1[2] = '\0';
								}
								if (tmpc2[1] == '\0') {
									tmpc2[1] = ' ';
									tmpc2[2] = '\0';
								}
								if (tmpc3[1] == '\0') {
									tmpc3[1] = ' ';
									tmpc3[2] = '\0';
								}
								suc =
									empangle(tmpc1, tmpc2, tmpc3, name1,
											 name2, name3,
											 atom[i].atomicnum,
											 atom[j].atomicnum,
											 atom[k].atomicnum);
							}
							/* only use similarname */


							if (suc == 0 && similarindex[i] == 1
								&& similarindex[j] == 1
								&& similarindex[k] == 1) {
								strcpy(tmpc1, similarname[i].name);
								strcpy(tmpc2, similarname[j].name);
								strcpy(tmpc3, similarname[k].name);
								if (tmpc1[1] == '\0') {
									tmpc1[1] = ' ';
									tmpc1[2] = '\0';
								}
								if (tmpc2[1] == '\0') {
									tmpc2[1] = ' ';
									tmpc2[2] = '\0';
								}
								if (tmpc3[1] == '\0') {
									tmpc3[1] = ' ';
									tmpc3[2] = '\0';
								}
								suc =
									angle(name1, name2, name3, tmpc1,
										  tmpc2, tmpc3, 1);
							}

							/* from here estimate the bond angle parameters with empirical method for similaresponding names */
							if (suc == 0) {
								strcpy(tmpc1, similarname[i].name);
								strcpy(tmpc2, similarname[j].name);
								strcpy(tmpc3, similarname[k].name);
								if (tmpc1[1] == '\0') {
									tmpc1[1] = ' ';
									tmpc1[2] = '\0';
								}
								if (tmpc2[1] == '\0') {
									tmpc2[1] = ' ';
									tmpc2[2] = '\0';
								}
								if (tmpc3[1] == '\0') {
									tmpc3[1] = ' ';
									tmpc3[2] = '\0';
								}
								suc =
									empangle(tmpc1, tmpc2, tmpc3, name1,
											 name2, name3,
											 atom[i].atomicnum,
											 atom[j].atomicnum,
											 atom[k].atomicnum);
							}


							if (suc == 0) {
								strcpy(tmpc1, atom[i].ambername);
								strcpy(tmpc2, atom[j].ambername);
								strcpy(tmpc3, atom[k].ambername);
								if (tmpc1[1] == '\0') {
									tmpc1[1] = ' ';
									tmpc1[2] = '\0';
								}
								if (tmpc2[1] == '\0') {
									tmpc2[1] = ' ';
									tmpc2[2] = '\0';
								}
								if (tmpc3[1] == '\0') {
									tmpc3[1] = ' ';
									tmpc3[2] = '\0';
								}
								strcpy(name1, tmpc1);
								strcpy(name2, tmpc2);
								strcpy(name3, tmpc3);
							}


							suc =
								angle(name1, name2, name3, tmpc1, tmpc2,
									  tmpc3, 0);
							if (suc == 0 && similarindex2[i] == 1) {
								strcpy(tmpc1, similarname2[i].name);
								strcpy(tmpc2, atom[j].ambername);
								strcpy(tmpc3, atom[k].ambername);
								if (tmpc1[1] == '\0') {
									tmpc1[1] = ' ';
									tmpc1[2] = '\0';
								}
								if (tmpc2[1] == '\0') {
									tmpc2[1] = ' ';
									tmpc2[2] = '\0';
								}
								if (tmpc3[1] == '\0') {
									tmpc3[1] = ' ';
									tmpc3[2] = '\0';
								}
								suc =
									angle(name1, name2, name3, tmpc1,
										  tmpc2, tmpc3, 1);
							}

							if (suc == 0 && similarindex2[j] == 1) {
								strcpy(tmpc1, atom[i].ambername);
								strcpy(tmpc2, similarname2[j].name);
								strcpy(tmpc3, atom[k].ambername);
								if (tmpc1[1] == '\0') {
									tmpc1[1] = ' ';
									tmpc1[2] = '\0';
								}
								if (tmpc2[1] == '\0') {
									tmpc2[1] = ' ';
									tmpc2[2] = '\0';
								}
								if (tmpc3[1] == '\0') {
									tmpc3[1] = ' ';
									tmpc3[2] = '\0';
								}
								suc =
									angle(name1, name2, name3, tmpc1,
										  tmpc2, tmpc3, 1);
							}

							if (suc == 0 && similarindex2[k] == 1) {
								strcpy(tmpc1, atom[i].ambername);
								strcpy(tmpc2, atom[j].ambername);
								strcpy(tmpc3, similarname2[k].name);
								if (tmpc1[1] == '\0') {
									tmpc1[1] = ' ';
									tmpc1[2] = '\0';
								}
								if (tmpc2[1] == '\0') {
									tmpc2[1] = ' ';
									tmpc2[2] = '\0';
								}
								if (tmpc3[1] == '\0') {
									tmpc3[1] = ' ';
									tmpc3[2] = '\0';
								}
								suc =
									angle(name1, name2, name3, tmpc1,
										  tmpc2, tmpc3, 1);
							}


							/* from here estimate the bond angle parameter using empirical method */
							if (suc == 0) {
								strcpy(tmpc1, similarname2[i].name);
								strcpy(tmpc2, atom[j].ambername);
								strcpy(tmpc3, atom[k].ambername);
								if (tmpc1[1] == '\0') {
									tmpc1[1] = ' ';
									tmpc1[2] = '\0';
								}
								if (tmpc2[1] == '\0') {
									tmpc2[1] = ' ';
									tmpc2[2] = '\0';
								}
								if (tmpc3[1] == '\0') {
									tmpc3[1] = ' ';
									tmpc3[2] = '\0';
								}
								suc =
									empangle(tmpc1, tmpc2, tmpc3, name1,
											 name2, name3,
											 atom[i].atomicnum,
											 atom[j].atomicnum,
											 atom[k].atomicnum);
							}
							if (suc == 0) {
								strcpy(tmpc1, atom[i].ambername);
								strcpy(tmpc2, similarname2[j].name);
								strcpy(tmpc3, atom[k].ambername);
								if (tmpc1[1] == '\0') {
									tmpc1[1] = ' ';
									tmpc1[2] = '\0';
								}
								if (tmpc2[1] == '\0') {
									tmpc2[1] = ' ';
									tmpc2[2] = '\0';
								}
								if (tmpc3[1] == '\0') {
									tmpc3[1] = ' ';
									tmpc3[2] = '\0';
								}
								suc =
									empangle(tmpc1, tmpc2, tmpc3, name1,
											 name2, name3,
											 atom[i].atomicnum,
											 atom[j].atomicnum,
											 atom[k].atomicnum);
							}
							if (suc == 0) {
								strcpy(tmpc1, atom[i].ambername);
								strcpy(tmpc2, atom[j].ambername);
								strcpy(tmpc3, similarname2[k].name);
								if (tmpc1[1] == '\0') {
									tmpc1[1] = ' ';
									tmpc1[2] = '\0';
								}
								if (tmpc2[1] == '\0') {
									tmpc2[1] = ' ';
									tmpc2[2] = '\0';
								}
								if (tmpc3[1] == '\0') {
									tmpc3[1] = ' ';
									tmpc3[2] = '\0';
								}
								suc =
									empangle(tmpc1, tmpc2, tmpc3, name1,
											 name2, name3,
											 atom[i].atomicnum,
											 atom[j].atomicnum,
											 atom[k].atomicnum);
							}
							/* only use similarname2 */

							/* combination of similar and similar2 */
							if (suc == 0 && similarindex[i] == 1
								&& similarindex2[j] == 1
								&& similarindex2[k] == 1) {
								strcpy(tmpc1, similarname[i].name);
								strcpy(tmpc2, similarname2[j].name);
								strcpy(tmpc3, similarname2[k].name);
								if (tmpc1[1] == '\0') {
									tmpc1[1] = ' ';
									tmpc1[2] = '\0';
								}
								if (tmpc2[1] == '\0') {
									tmpc2[1] = ' ';
									tmpc2[2] = '\0';
								}
								if (tmpc3[1] == '\0') {
									tmpc3[1] = ' ';
									tmpc3[2] = '\0';
								}
								suc =
									angle(name1, name2, name3, tmpc1,
										  tmpc2, tmpc3, 1);
							}
							if (suc == 0 && similarindex2[i] == 1
								&& similarindex[j] == 1
								&& similarindex2[k] == 1) {
								strcpy(tmpc1, similarname2[i].name);
								strcpy(tmpc2, similarname[j].name);
								strcpy(tmpc3, similarname2[k].name);
								if (tmpc1[1] == '\0') {
									tmpc1[1] = ' ';
									tmpc1[2] = '\0';
								}
								if (tmpc2[1] == '\0') {
									tmpc2[1] = ' ';
									tmpc2[2] = '\0';
								}
								if (tmpc3[1] == '\0') {
									tmpc3[1] = ' ';
									tmpc3[2] = '\0';
								}
								suc =
									angle(name1, name2, name3, tmpc1,
										  tmpc2, tmpc3, 1);
							}
							if (suc == 0 && similarindex2[i] == 1
								&& similarindex2[j] == 1
								&& similarindex[k] == 1) {
								strcpy(tmpc1, similarname2[i].name);
								strcpy(tmpc2, similarname2[j].name);
								strcpy(tmpc3, similarname[k].name);
								if (tmpc1[1] == '\0') {
									tmpc1[1] = ' ';
									tmpc1[2] = '\0';
								}
								if (tmpc2[1] == '\0') {
									tmpc2[1] = ' ';
									tmpc2[2] = '\0';
								}
								if (tmpc3[1] == '\0') {
									tmpc3[1] = ' ';
									tmpc3[2] = '\0';
								}
								suc =
									angle(name1, name2, name3, tmpc1,
										  tmpc2, tmpc3, 1);
							}
							if (suc == 0 && similarindex[i] == 1
								&& similarindex[j] == 1
								&& similarindex2[k] == 1) {
								strcpy(tmpc1, similarname[i].name);
								strcpy(tmpc2, similarname[j].name);
								strcpy(tmpc3, similarname2[k].name);
								if (tmpc1[1] == '\0') {
									tmpc1[1] = ' ';
									tmpc1[2] = '\0';
								}
								if (tmpc2[1] == '\0') {
									tmpc2[1] = ' ';
									tmpc2[2] = '\0';
								}
								if (tmpc3[1] == '\0') {
									tmpc3[1] = ' ';
									tmpc3[2] = '\0';
								}
								suc =
									angle(name1, name2, name3, tmpc1,
										  tmpc2, tmpc3, 1);
							}
							if (suc == 0 && similarindex[i] == 1
								&& similarindex2[j] == 1
								&& similarindex[k] == 1) {
								strcpy(tmpc1, similarname[i].name);
								strcpy(tmpc2, similarname2[j].name);
								strcpy(tmpc3, similarname[k].name);
								if (tmpc1[1] == '\0') {
									tmpc1[1] = ' ';
									tmpc1[2] = '\0';
								}
								if (tmpc2[1] == '\0') {
									tmpc2[1] = ' ';
									tmpc2[2] = '\0';
								}
								if (tmpc3[1] == '\0') {
									tmpc3[1] = ' ';
									tmpc3[2] = '\0';
								}
								suc =
									angle(name1, name2, name3, tmpc1,
										  tmpc2, tmpc3, 1);
							}
							if (suc == 0 && similarindex2[i] == 1
								&& similarindex[j] == 1
								&& similarindex[k] == 1) {
								strcpy(tmpc1, similarname2[i].name);
								strcpy(tmpc2, similarname[j].name);
								strcpy(tmpc3, similarname[k].name);
								if (tmpc1[1] == '\0') {
									tmpc1[1] = ' ';
									tmpc1[2] = '\0';
								}
								if (tmpc2[1] == '\0') {
									tmpc2[1] = ' ';
									tmpc2[2] = '\0';
								}
								if (tmpc3[1] == '\0') {
									tmpc3[1] = ' ';
									tmpc3[2] = '\0';
								}
								suc =
									angle(name1, name2, name3, tmpc1,
										  tmpc2, tmpc3, 1);
							}
							if (suc == 0 && similarindex2[i] == 1
								&& similarindex2[j] == 1
								&& similarindex2[k] == 1) {
								strcpy(tmpc1, similarname2[i].name);
								strcpy(tmpc2, similarname2[j].name);
								strcpy(tmpc3, similarname2[k].name);
								if (tmpc1[1] == '\0') {
									tmpc1[1] = ' ';
									tmpc1[2] = '\0';
								}
								if (tmpc2[1] == '\0') {
									tmpc2[1] = ' ';
									tmpc2[2] = '\0';
								}
								if (tmpc3[1] == '\0') {
									tmpc3[1] = ' ';
									tmpc3[2] = '\0';
								}
								suc =
									angle(name1, name2, name3, tmpc1,
										  tmpc2, tmpc3, 1);
							}
							/* end */
							/* from here estimate the bond angle parameters with empirical method for similaresponding names */
							/* combination of similar and similar2 */
							if (suc == 0) {
								strcpy(tmpc1, similarname[i].name);
								strcpy(tmpc2, similarname2[j].name);
								strcpy(tmpc3, similarname2[k].name);
								if (tmpc1[1] == '\0') {
									tmpc1[1] = ' ';
									tmpc1[2] = '\0';
								}
								if (tmpc2[1] == '\0') {
									tmpc2[1] = ' ';
									tmpc2[2] = '\0';
								}
								if (tmpc3[1] == '\0') {
									tmpc3[1] = ' ';
									tmpc3[2] = '\0';
								}
								suc =
									empangle(tmpc1, tmpc2, tmpc3, name1,
											 name2, name3,
											 atom[i].atomicnum,
											 atom[j].atomicnum,
											 atom[k].atomicnum);
							}
							if (suc == 0) {
								strcpy(tmpc1, similarname2[i].name);
								strcpy(tmpc2, similarname[j].name);
								strcpy(tmpc3, similarname2[k].name);
								if (tmpc1[1] == '\0') {
									tmpc1[1] = ' ';
									tmpc1[2] = '\0';
								}
								if (tmpc2[1] == '\0') {
									tmpc2[1] = ' ';
									tmpc2[2] = '\0';
								}
								if (tmpc3[1] == '\0') {
									tmpc3[1] = ' ';
									tmpc3[2] = '\0';
								}
								suc =
									empangle(tmpc1, tmpc2, tmpc3, name1,
											 name2, name3,
											 atom[i].atomicnum,
											 atom[j].atomicnum,
											 atom[k].atomicnum);
							}
							if (suc == 0) {
								strcpy(tmpc1, similarname2[i].name);
								strcpy(tmpc2, similarname2[j].name);
								strcpy(tmpc3, similarname[k].name);
								if (tmpc1[1] == '\0') {
									tmpc1[1] = ' ';
									tmpc1[2] = '\0';
								}
								if (tmpc2[1] == '\0') {
									tmpc2[1] = ' ';
									tmpc2[2] = '\0';
								}
								if (tmpc3[1] == '\0') {
									tmpc3[1] = ' ';
									tmpc3[2] = '\0';
								}
								suc =
									empangle(tmpc1, tmpc2, tmpc3, name1,
											 name2, name3,
											 atom[i].atomicnum,
											 atom[j].atomicnum,
											 atom[k].atomicnum);
							}
							if (suc == 0) {
								strcpy(tmpc1, similarname[i].name);
								strcpy(tmpc2, similarname[j].name);
								strcpy(tmpc3, similarname2[k].name);
								if (tmpc1[1] == '\0') {
									tmpc1[1] = ' ';
									tmpc1[2] = '\0';
								}
								if (tmpc2[1] == '\0') {
									tmpc2[1] = ' ';
									tmpc2[2] = '\0';
								}
								if (tmpc3[1] == '\0') {
									tmpc3[1] = ' ';
									tmpc3[2] = '\0';
								}
								suc =
									empangle(tmpc1, tmpc2, tmpc3, name1,
											 name2, name3,
											 atom[i].atomicnum,
											 atom[j].atomicnum,
											 atom[k].atomicnum);
							}
							if (suc == 0) {
								strcpy(tmpc1, similarname[i].name);
								strcpy(tmpc2, similarname2[j].name);
								strcpy(tmpc3, similarname[k].name);
								if (tmpc1[1] == '\0') {
									tmpc1[1] = ' ';
									tmpc1[2] = '\0';
								}
								if (tmpc2[1] == '\0') {
									tmpc2[1] = ' ';
									tmpc2[2] = '\0';
								}
								if (tmpc3[1] == '\0') {
									tmpc3[1] = ' ';
									tmpc3[2] = '\0';
								}
								suc =
									empangle(tmpc1, tmpc2, tmpc3, name1,
											 name2, name3,
											 atom[i].atomicnum,
											 atom[j].atomicnum,
											 atom[k].atomicnum);
							}
							if (suc == 0) {
								strcpy(tmpc1, similarname2[i].name);
								strcpy(tmpc2, similarname[j].name);
								strcpy(tmpc3, similarname[k].name);
								if (tmpc1[1] == '\0') {
									tmpc1[1] = ' ';
									tmpc1[2] = '\0';
								}
								if (tmpc2[1] == '\0') {
									tmpc2[1] = ' ';
									tmpc2[2] = '\0';
								}
								if (tmpc3[1] == '\0') {
									tmpc3[1] = ' ';
									tmpc3[2] = '\0';
								}
								suc =
									empangle(tmpc1, tmpc2, tmpc3, name1,
											 name2, name3,
											 atom[i].atomicnum,
											 atom[j].atomicnum,
											 atom[k].atomicnum);
							}
							if (suc == 0) {
								strcpy(tmpc1, similarname2[i].name);
								strcpy(tmpc2, similarname2[j].name);
								strcpy(tmpc3, similarname2[k].name);
								if (tmpc1[1] == '\0') {
									tmpc1[1] = ' ';
									tmpc1[2] = '\0';
								}
								if (tmpc2[1] == '\0') {
									tmpc2[1] = ' ';
									tmpc2[2] = '\0';
								}
								if (tmpc3[1] == '\0') {
									tmpc3[1] = ' ';
									tmpc3[2] = '\0';
								}
								suc =
									empangle(tmpc1, tmpc2, tmpc3, name1,
											 name2, name3,
											 atom[i].atomicnum,
											 atom[j].atomicnum,
											 atom[k].atomicnum);
							}

							/*end */
							if (suc == 0) {
								strcpy(tmpc1, atom[i].ambername);
								strcpy(tmpc2, atom[j].ambername);
								strcpy(tmpc3, atom[k].ambername);
								if (tmpc1[1] == '\0') {
									tmpc1[1] = ' ';
									tmpc1[2] = '\0';
								}
								if (tmpc2[1] == '\0') {
									tmpc2[1] = ' ';
									tmpc2[2] = '\0';
								}
								if (tmpc3[1] == '\0') {
									tmpc3[1] = ' ';
									tmpc3[2] = '\0';
								}
								strcpy(name1, tmpc1);
								strcpy(name2, tmpc2);
								strcpy(name3, tmpc3);
							}
							if (suc == 0 && corrindex[i] == 1) {
								strcpy(tmpc1, corrname[i].name);
								strcpy(tmpc2, atom[j].ambername);
								strcpy(tmpc3, atom[k].ambername);
								if (tmpc1[1] == '\0') {
									tmpc1[1] = ' ';
									tmpc1[2] = '\0';
								}
								if (tmpc2[1] == '\0') {
									tmpc2[1] = ' ';
									tmpc2[2] = '\0';
								}
								if (tmpc3[1] == '\0') {
									tmpc3[1] = ' ';
									tmpc3[2] = '\0';
								}
								suc =
									angle(name1, name2, name3, tmpc1,
										  tmpc2, tmpc3, 1);
							}

							if (suc == 0 && corrindex[j] == 1) {
								strcpy(tmpc1, atom[i].ambername);
								strcpy(tmpc2, corrname[j].name);
								strcpy(tmpc3, atom[k].ambername);
								if (tmpc1[1] == '\0') {
									tmpc1[1] = ' ';
									tmpc1[2] = '\0';
								}
								if (tmpc2[1] == '\0') {
									tmpc2[1] = ' ';
									tmpc2[2] = '\0';
								}
								if (tmpc3[1] == '\0') {
									tmpc3[1] = ' ';
									tmpc3[2] = '\0';
								}
								suc =
									angle(name1, name2, name3, tmpc1,
										  tmpc2, tmpc3, 1);
							}

							if (suc == 0 && corrindex[k] == 1) {
								strcpy(tmpc1, atom[i].ambername);
								strcpy(tmpc2, atom[j].ambername);
								strcpy(tmpc3, corrname[k].name);
								if (tmpc1[1] == '\0') {
									tmpc1[1] = ' ';
									tmpc1[2] = '\0';
								}
								if (tmpc2[1] == '\0') {
									tmpc2[1] = ' ';
									tmpc2[2] = '\0';
								}
								if (tmpc3[1] == '\0') {
									tmpc3[1] = ' ';
									tmpc3[2] = '\0';
								}
								suc =
									angle(name1, name2, name3, tmpc1,
										  tmpc2, tmpc3, 1);
							}


							/* from here estimate the bond angle parameter using empirical method */
							if (suc == 0) {
								strcpy(tmpc1, corrname[i].name);
								strcpy(tmpc2, atom[j].ambername);
								strcpy(tmpc3, atom[k].ambername);
								if (tmpc1[1] == '\0') {
									tmpc1[1] = ' ';
									tmpc1[2] = '\0';
								}
								if (tmpc2[1] == '\0') {
									tmpc2[1] = ' ';
									tmpc2[2] = '\0';
								}
								if (tmpc3[1] == '\0') {
									tmpc3[1] = ' ';
									tmpc3[2] = '\0';
								}
								suc =
									empangle(tmpc1, tmpc2, tmpc3, name1,
											 name2, name3,
											 atom[i].atomicnum,
											 atom[j].atomicnum,
											 atom[k].atomicnum);
							}
							if (suc == 0) {
								strcpy(tmpc1, atom[i].ambername);
								strcpy(tmpc2, corrname[j].name);
								strcpy(tmpc3, atom[k].ambername);
								if (tmpc1[1] == '\0') {
									tmpc1[1] = ' ';
									tmpc1[2] = '\0';
								}
								if (tmpc2[1] == '\0') {
									tmpc2[1] = ' ';
									tmpc2[2] = '\0';
								}
								if (tmpc3[1] == '\0') {
									tmpc3[1] = ' ';
									tmpc3[2] = '\0';
								}
								suc =
									empangle(tmpc1, tmpc2, tmpc3, name1,
											 name2, name3,
											 atom[i].atomicnum,
											 atom[j].atomicnum,
											 atom[k].atomicnum);
							}
							if (suc == 0) {
								strcpy(tmpc1, atom[i].ambername);
								strcpy(tmpc2, atom[j].ambername);
								strcpy(tmpc3, corrname[k].name);
								if (tmpc1[1] == '\0') {
									tmpc1[1] = ' ';
									tmpc1[2] = '\0';
								}
								if (tmpc2[1] == '\0') {
									tmpc2[1] = ' ';
									tmpc2[2] = '\0';
								}
								if (tmpc3[1] == '\0') {
									tmpc3[1] = ' ';
									tmpc3[2] = '\0';
								}
								suc =
									empangle(tmpc1, tmpc2, tmpc3, name1,
											 name2, name3,
											 atom[i].atomicnum,
											 atom[j].atomicnum,
											 atom[k].atomicnum);
							}
							/* only use corrname */


							if (suc == 0 && corrindex[i] == 1
								&& corrindex[j] == 1
								&& corrindex[k] == 1) {
								strcpy(tmpc1, corrname[i].name);
								strcpy(tmpc2, corrname[j].name);
								strcpy(tmpc3, corrname[k].name);
								if (tmpc1[1] == '\0') {
									tmpc1[1] = ' ';
									tmpc1[2] = '\0';
								}
								if (tmpc2[1] == '\0') {
									tmpc2[1] = ' ';
									tmpc2[2] = '\0';
								}
								if (tmpc3[1] == '\0') {
									tmpc3[1] = ' ';
									tmpc3[2] = '\0';
								}
								suc =
									angle(name1, name2, name3, tmpc1,
										  tmpc2, tmpc3, 1);
							}

							/* from here estimate the bond angle parameters with empirical method for corresponding names */
							if (suc == 0) {
								strcpy(tmpc1, corrname[i].name);
								strcpy(tmpc2, corrname[j].name);
								strcpy(tmpc3, corrname[k].name);
								suc =
									empangle(tmpc1, tmpc2, tmpc3, name1,
											 name2, name3,
											 atom[i].atomicnum,
											 atom[j].atomicnum,
											 atom[k].atomicnum);
							}
							if (suc == 0) {

                                                                if (angleparm[l].flag == 0) {
						                  fprintf(fpout, "%2s %2s %2s%9.3lf  %10.3lf",
							                  name1, name2, name3, 0.0, 0.0);
						                  fprintf(fpout, "      !%s\n",
							                  "ATTN, need revision");
						                  angleparm[l].flag = 1;
					                        }
								strcpy(angleparm[angleparmnum].name1,
									   name1);
								strcpy(angleparm[angleparmnum].name2,
									   name2);
								strcpy(angleparm[angleparmnum].name3,
									   name3);
								angleparmnum++;
								if (angleparmnum >= maxangleparm) {
									maxangleparm += MAX_FF_ANGLE;
									angleparm =
										(ANGLE *) realloc(angleparm,
														  sizeof(ANGLE) *
														  maxangleparm);
									if (angleparm == NULL) {
										fprintf(stdout,
												"memory allocation error for *angleparm\n");
										exit(1);
									}
								}
                                                        }
				                }
                                        }
			}
		}
	}
}



int torsion(char *at_name1, char *at_name2, char *at_name3, char *at_name4,
			char *corr_name1, char *corr_name2, char *corr_name3,
			char *corr_name4, int index)
{
	int suc = 0;
	int m, n;
	for (m = 0; m < torsionparmnum; m++)
		if ((torsionparm[m].name1[0] == corr_name1[0]
			 && torsionparm[m].name1[1] == corr_name1[1]
			 && torsionparm[m].name2[0] == corr_name2[0]
			 && torsionparm[m].name2[1] == corr_name2[1]
			 && torsionparm[m].name3[0] == corr_name3[0]
			 && torsionparm[m].name3[1] == corr_name3[1]
			 && torsionparm[m].name4[0] == corr_name4[0]
			 && torsionparm[m].name4[1] == corr_name4[1])
			|| (torsionparm[m].name4[0] == corr_name1[0]
				&& torsionparm[m].name4[1] == corr_name1[1]
				&& torsionparm[m].name3[0] == corr_name2[0]
				&& torsionparm[m].name3[1] == corr_name2[1]
				&& torsionparm[m].name2[0] == corr_name3[0]
				&& torsionparm[m].name2[1] == corr_name3[1]
				&& torsionparm[m].name1[0] == corr_name4[0]
				&& torsionparm[m].name1[1] == corr_name4[1])
			|| (torsionparm[m].name1[0] == 'X'
				&& torsionparm[m].name1[1] == ' '
				&& torsionparm[m].name4[0] == 'X'
				&& torsionparm[m].name4[1] == ' '
				&& torsionparm[m].name2[0] == corr_name2[0]
				&& torsionparm[m].name2[1] == corr_name2[1]
				&& torsionparm[m].name3[0] == corr_name3[0]
				&& torsionparm[m].name3[1] == corr_name3[1])
			|| (torsionparm[m].name1[0] == 'X'
				&& torsionparm[m].name1[1] == ' '
				&& torsionparm[m].name4[0] == 'X'
				&& torsionparm[m].name4[1] == ' '
				&& torsionparm[m].name2[0] == corr_name3[0]
				&& torsionparm[m].name2[1] == corr_name3[1]
				&& torsionparm[m].name3[0] == corr_name2[0]
				&& torsionparm[m].name3[1] == corr_name2[1])) {
			suc = 1;
			n = m;
			if (index == 1) {
				while (torsionparm[n].fterm < 0) {
                       	               fprintf(ffout, "%2s %2s %2s %2s%10.3lf%10.0lf%10.1lf",
  					      at_name1, at_name2, at_name3, at_name4,
			                      torsionparm[m].force / torsionparm[m].mul,
			                      fabs(torsionparm[m].fterm), torsionparm[m].phase);
			 	       fprintf(ffout,
				  		     "      !%s %2s %2s %2s %2s\n",
						     "same as", torsionparm[n].name1,
						     torsionparm[n].name2,
						     torsionparm[n].name3,
						     torsionparm[n].name4);
                               		torsionparm[torsionparmnum] = torsionparm[n];
                                 	torsionparm[torsionparmnum].flag = 1;
					strcpy(torsionparm[torsionparmnum].name1, at_name1);
					strcpy(torsionparm[torsionparmnum].name2, at_name2);
					strcpy(torsionparm[torsionparmnum].name3, at_name3);
					strcpy(torsionparm[torsionparmnum].name4, at_name4);
					torsionparmnum++;
					if (torsionparmnum >= maxtorsionparm) {
						maxtorsionparm += MAX_FF_TORSION;
						torsionparm =
							(TORSION *) realloc(torsionparm,
												sizeof(TORSION) *
												maxtorsionparm);
						if (torsionparm == NULL) {
							fprintf(stdout,
									"memory allocation error for *torsionparm\n");
							exit(1);
						}
					}
					n++;
				}

                    	        fprintf(ffout, "%2s %2s %2s %2s%10.3lf%10.0lf%10.1lf",
  					      at_name1, at_name2, at_name3, at_name4,
			                      torsionparm[m].force / torsionparm[m].mul,
			                      fabs(torsionparm[m].fterm), torsionparm[m].phase);
				fprintf(ffout,
						"      ! %s %2s %2s %2s %2s\n",
						"same as",
						torsionparm[n].name1,
				  	        torsionparm[n].name2,
						torsionparm[n].name3,
						torsionparm[n].name4);
				torsionparm[torsionparmnum] = torsionparm[n];
				torsionparm[torsionparmnum].flag = 1;
                                strcpy(torsionparm[torsionparmnum].name1, at_name1);
				strcpy(torsionparm[torsionparmnum].name2, at_name2);
				strcpy(torsionparm[torsionparmnum].name3, at_name3);
				strcpy(torsionparm[torsionparmnum].name4, at_name4);
				torsionparmnum++;
				if (torsionparmnum >= maxtorsionparm) {
					maxtorsionparm += MAX_FF_TORSION;
					torsionparm =
						(TORSION *) realloc(torsionparm,
											sizeof(TORSION) *
											maxtorsionparm);
					if (torsionparm == NULL) {
						fprintf(stdout,
								"memory allocation error for *torsionparm\n");
						exit(1);
					}
				}

				break;

			} else if (torsionparm[m].flag == 0) {
		           fprintf(ffout, "%2s %2s %2s %2s%10.3lf%10.0lf%10.1lf\n",
		           torsionparm[m].name1,
			   torsionparm[m].name2,
			   torsionparm[m].name3,
			   torsionparm[m].name4,
			   torsionparm[m].force / torsionparm[m].mul,
			   fabs(torsionparm[m].fterm), torsionparm[m].phase);
		           torsionparm[m].flag = 1;
		        }
		}

	return suc;
}

void wrt_torsion(void)
{

	int i, j, k, l;
	int suc;
        float twistv;
	char tmpc1[5], tmpc2[5], tmpc3[5], tmpc4[5];
	char name1[5], name2[5], name3[5], name4[5];

        fprintf(fpout, "\n");
        fprintf(ffout, "\n%s\n", "DIHEDRAL");

	/* NB: non-standard indentation in next four lines; for readability  */
        for (j = 0; j < atomnum; j++) {
        	for (k = j+1; k < atomnum; k++) {
	 	       for (l = 0; l < atomnum; l++) {
	                      for (i = 0; i < atomnum; i++) {
					if (i != k && l != j)   {
						if (atom[i].con[0] == j || atom[i].con[1] == j
							|| atom[i].con[2] == j || atom[i].con[3] == j
							|| atom[i].con[4] == j
							|| atom[i].con[5] == j) {
							if (atom[j].con[0] == k || atom[j].con[1] == k
								|| atom[j].con[2] == k
								|| atom[j].con[3] == k
								|| atom[j].con[4] == k
								|| atom[j].con[5] == k)
								if (atom[l].con[0] == k
									|| atom[l].con[1] == k
									|| atom[l].con[2] == k
									|| atom[l].con[3] == k
									|| atom[l].con[4] == k
									|| atom[l].con[5] == k) {
									
                                                                        /* New rule to prevent 3 atoms
                                                                           dihedrals VEB CHPC(2004)*/
                                                                        if (i == l) continue;

                                                                        suc = 0;
									strcpy(tmpc1, atom[i].ambername);
									strcpy(tmpc2, atom[j].ambername);
									strcpy(tmpc3, atom[k].ambername);
									strcpy(tmpc4, atom[l].ambername);

									if (tmpc1[1] == '\0') {
										tmpc1[1] = ' ';
										tmpc1[2] = '\0';
									}
									if (tmpc2[1] == '\0') {
										tmpc2[1] = ' ';
										tmpc2[2] = '\0';
									}
									if (tmpc3[1] == '\0') {
										tmpc3[1] = ' ';
										tmpc3[2] = '\0';
									}
									if (tmpc4[1] == '\0') {
										tmpc4[1] = ' ';
										tmpc4[2] = '\0';
									}

									strcpy(name1, tmpc1);
									strcpy(name2, tmpc2);
									strcpy(name3, tmpc3);
									strcpy(name4, tmpc4);

                                                                        twistv = rotate(atom[i], atom[j], atom[k], &atom[l]);
                                                                        fprintf(fpout, "DIHE %-5s %-5s %-5s %-5s"
                                                                                , atom[i].name, atom[j].name, atom[k].name, atom[l].name);
                                                                                fprintf(fpout, "	! dihe	%10.4lf\n", twistv);

									suc =
										torsion(name1, name2, name3, name4,
												tmpc1, tmpc2, tmpc3, tmpc4,
												0);

									if (suc == 0 && similarindex[i] == 1) {
										strcpy(tmpc1, similarname[i].name);
										if (tmpc1[1] == '\0') {
											tmpc1[1] = ' ';
											tmpc1[2] = '\0';
										}
										suc =
											torsion(name1, name2, name3,
													name4, tmpc1, tmpc2,
													tmpc3, tmpc4, 1);
									}

									if (suc == 0 && similarindex[j] == 1) {
										strcpy(tmpc1, atom[i].ambername);
										strcpy(tmpc2, similarname[j].name);
										if (tmpc1[1] == '\0') {
											tmpc1[1] = ' ';
											tmpc1[2] = '\0';
										}
										if (tmpc2[1] == '\0') {
											tmpc2[1] = ' ';
											tmpc2[2] = '\0';
										}
										suc =
											torsion(name1, name2, name3,
													name4, tmpc1, tmpc2,
													tmpc3, tmpc4, 1);
									}
									if (suc == 0 && similarindex[k] == 1) {
										strcpy(tmpc2, atom[j].ambername);
										strcpy(tmpc3, similarname[k].name);
										if (tmpc2[1] == '\0') {
											tmpc2[1] = ' ';
											tmpc2[2] = '\0';
										}
										if (tmpc3[1] == '\0') {
											tmpc3[1] = ' ';
											tmpc3[2] = '\0';
										}
										suc =
											torsion(name1, name2, name3,
													name4, tmpc1, tmpc2,
													tmpc3, tmpc4, 1);
									}

									if (suc == 0 && similarindex[l] == 1) {
										strcpy(tmpc3, atom[k].ambername);
										strcpy(tmpc4, similarname[l].name);
										if (tmpc3[1] == '\0') {
											tmpc3[1] = ' ';
											tmpc3[2] = '\0';
										}
										if (tmpc4[1] == '\0') {
											tmpc4[1] = ' ';
											tmpc4[2] = '\0';
										}
										suc =
											torsion(name1, name2, name3,
													name4, tmpc1, tmpc2,
													tmpc3, tmpc4, 1);
									}

									if (suc == 0 && similarindex[i] == 1
										&& similarindex[j] == 1
										&& similarindex[k] == 1
										&& similarindex[l] == 1) {
										strcpy(tmpc1, similarname[i].name);
										strcpy(tmpc2, similarname[j].name);
										strcpy(tmpc3, similarname[k].name);
										strcpy(tmpc4, similarname[l].name);
										if (tmpc1[1] == '\0') {
											tmpc1[1] = ' ';
											tmpc1[2] = '\0';
										}
										if (tmpc2[1] == '\0') {
											tmpc2[1] = ' ';
											tmpc2[2] = '\0';
										}
										if (tmpc3[1] == '\0') {
											tmpc3[1] = ' ';
											tmpc3[2] = '\0';
										}
										if (tmpc4[1] == '\0') {
											tmpc4[1] = ' ';
											tmpc4[2] = '\0';
										}
										suc =
											torsion(name1, name2, name3,
													name4, tmpc1, tmpc2,
													tmpc3, tmpc4, 1);
									}
									if (suc == 0) {
										strcpy(tmpc1, atom[i].ambername);
										strcpy(tmpc2, atom[j].ambername);
										strcpy(tmpc3, atom[k].ambername);
										strcpy(tmpc4, atom[l].ambername);

										if (tmpc1[1] == '\0') {
											tmpc1[1] = ' ';
											tmpc1[2] = '\0';
										}
										if (tmpc2[1] == '\0') {
											tmpc2[1] = ' ';
											tmpc2[2] = '\0';
										}
										if (tmpc3[1] == '\0') {
											tmpc3[1] = ' ';
											tmpc3[2] = '\0';
										}
										if (tmpc4[1] == '\0') {
											tmpc4[1] = ' ';
											tmpc4[2] = '\0';
										}

										strcpy(name1, tmpc1);
										strcpy(name2, tmpc2);
										strcpy(name3, tmpc3);
										strcpy(name4, tmpc4);
									}

									suc =
										torsion(name1, name2, name3, name4,
												tmpc1, tmpc2, tmpc3, tmpc4,
												0);

									if (suc == 0 && similarindex2[i] == 1) {
										strcpy(tmpc1,
											   similarname2[i].name);
										if (tmpc1[1] == '\0') {
											tmpc1[1] = ' ';
											tmpc1[2] = '\0';
										}
										suc =
											torsion(name1, name2, name3,
													name4, tmpc1, tmpc2,
													tmpc3, tmpc4, 1);
									}

									if (suc == 0 && similarindex2[j] == 1) {
										strcpy(tmpc1, atom[i].ambername);
										strcpy(tmpc2,
											   similarname2[j].name);
										if (tmpc1[1] == '\0') {
											tmpc1[1] = ' ';
											tmpc1[2] = '\0';
										}
										if (tmpc2[1] == '\0') {
											tmpc2[1] = ' ';
											tmpc2[2] = '\0';
										}
										suc =
											torsion(name1, name2, name3,
													name4, tmpc1, tmpc2,
													tmpc3, tmpc4, 1);
									}
									if (suc == 0 && similarindex2[k] == 1) {
										strcpy(tmpc2, atom[j].ambername);
										strcpy(tmpc3,
											   similarname2[k].name);
										if (tmpc2[1] == '\0') {
											tmpc2[1] = ' ';
											tmpc2[2] = '\0';
										}
										if (tmpc3[1] == '\0') {
											tmpc3[1] = ' ';
											tmpc3[2] = '\0';
										}
										suc =
											torsion(name1, name2, name3,
													name4, tmpc1, tmpc2,
													tmpc3, tmpc4, 1);
									}

									if (suc == 0 && similarindex2[l] == 1) {
										strcpy(tmpc3, atom[k].ambername);
										strcpy(tmpc4,
											   similarname2[l].name);
										if (tmpc3[1] == '\0') {
											tmpc3[1] = ' ';
											tmpc3[2] = '\0';
										}
										if (tmpc4[1] == '\0') {
											tmpc4[1] = ' ';
											tmpc4[2] = '\0';
										}
										suc =
											torsion(name1, name2, name3,
													name4, tmpc1, tmpc2,
													tmpc3, tmpc4, 1);
									}

									/* different combinations for similar and similar2 */
									if (suc == 0 && similarindex[i] == 1
										&& similarindex2[j] == 1
										&& similarindex2[k] == 1
										&& similarindex2[l] == 1) {
										strcpy(tmpc1, similarname[i].name);
										strcpy(tmpc2,
											   similarname2[j].name);
										strcpy(tmpc3,
											   similarname2[k].name);
										strcpy(tmpc4,
											   similarname2[l].name);
										if (tmpc1[1] == '\0') {
											tmpc1[1] = ' ';
											tmpc1[2] = '\0';
										}
										if (tmpc2[1] == '\0') {
											tmpc2[1] = ' ';
											tmpc2[2] = '\0';
										}
										if (tmpc3[1] == '\0') {
											tmpc3[1] = ' ';
											tmpc3[2] = '\0';
										}
										if (tmpc4[1] == '\0') {
											tmpc4[1] = ' ';
											tmpc4[2] = '\0';
										}
										suc =
											torsion(name1, name2, name3,
													name4, tmpc1, tmpc2,
													tmpc3, tmpc4, 1);
									}
									if (suc == 0 && similarindex2[i] == 1
										&& similarindex[j] == 1
										&& similarindex2[k] == 1
										&& similarindex2[l] == 1) {
										strcpy(tmpc1,
											   similarname2[i].name);
										strcpy(tmpc2, similarname[j].name);
										strcpy(tmpc3,
											   similarname2[k].name);
										strcpy(tmpc4,
											   similarname2[l].name);
										if (tmpc1[1] == '\0') {
											tmpc1[1] = ' ';
											tmpc1[2] = '\0';
										}
										if (tmpc2[1] == '\0') {
											tmpc2[1] = ' ';
											tmpc2[2] = '\0';
										}
										if (tmpc3[1] == '\0') {
											tmpc3[1] = ' ';
											tmpc3[2] = '\0';
										}
										if (tmpc4[1] == '\0') {
											tmpc4[1] = ' ';
											tmpc4[2] = '\0';
										}
										suc =
											torsion(name1, name2, name3,
													name4, tmpc1, tmpc2,
													tmpc3, tmpc4, 1);
									}
									if (suc == 0 && similarindex2[i] == 1
										&& similarindex2[j] == 1
										&& similarindex[k] == 1
										&& similarindex2[l] == 1) {
										strcpy(tmpc1,
											   similarname2[i].name);
										strcpy(tmpc2,
											   similarname2[j].name);
										strcpy(tmpc3, similarname[k].name);
										strcpy(tmpc4,
											   similarname2[l].name);
										if (tmpc1[1] == '\0') {
											tmpc1[1] = ' ';
											tmpc1[2] = '\0';
										}
										if (tmpc2[1] == '\0') {
											tmpc2[1] = ' ';
											tmpc2[2] = '\0';
										}
										if (tmpc3[1] == '\0') {
											tmpc3[1] = ' ';
											tmpc3[2] = '\0';
										}
										if (tmpc4[1] == '\0') {
											tmpc4[1] = ' ';
											tmpc4[2] = '\0';
										}
										suc =
											torsion(name1, name2, name3,
													name4, tmpc1, tmpc2,
													tmpc3, tmpc4, 1);
									}
									if (suc == 0 && similarindex2[i] == 1
										&& similarindex2[j] == 1
										&& similarindex2[k] == 1
										&& similarindex[l] == 1) {
										strcpy(tmpc1,
											   similarname2[i].name);
										strcpy(tmpc2,
											   similarname2[j].name);
										strcpy(tmpc3,
											   similarname2[k].name);
										strcpy(tmpc4, similarname[l].name);
										if (tmpc1[1] == '\0') {
											tmpc1[1] = ' ';
											tmpc1[2] = '\0';
										}
										if (tmpc2[1] == '\0') {
											tmpc2[1] = ' ';
											tmpc2[2] = '\0';
										}
										if (tmpc3[1] == '\0') {
											tmpc3[1] = ' ';
											tmpc3[2] = '\0';
										}
										if (tmpc4[1] == '\0') {
											tmpc4[1] = ' ';
											tmpc4[2] = '\0';
										}
										suc =
											torsion(name1, name2, name3,
													name4, tmpc1, tmpc2,
													tmpc3, tmpc4, 1);
									}
									if (suc == 0 && similarindex[i] == 1
										&& similarindex[j] == 1
										&& similarindex2[k] == 1
										&& similarindex2[l] == 1) {
										strcpy(tmpc1, similarname[i].name);
										strcpy(tmpc2, similarname[j].name);
										strcpy(tmpc3,
											   similarname2[k].name);
										strcpy(tmpc4,
											   similarname2[l].name);
										if (tmpc1[1] == '\0') {
											tmpc1[1] = ' ';
											tmpc1[2] = '\0';
										}
										if (tmpc2[1] == '\0') {
											tmpc2[1] = ' ';
											tmpc2[2] = '\0';
										}
										if (tmpc3[1] == '\0') {
											tmpc3[1] = ' ';
											tmpc3[2] = '\0';
										}
										if (tmpc4[1] == '\0') {
											tmpc4[1] = ' ';
											tmpc4[2] = '\0';
										}
										suc =
											torsion(name1, name2, name3,
													name4, tmpc1, tmpc2,
													tmpc3, tmpc4, 1);
									}
									if (suc == 0 && similarindex[i] == 1
										&& similarindex2[j] == 1
										&& similarindex[k] == 1
										&& similarindex2[l] == 1) {
										strcpy(tmpc1, similarname[i].name);
										strcpy(tmpc2,
											   similarname2[j].name);
										strcpy(tmpc3, similarname[k].name);
										strcpy(tmpc4,
											   similarname2[l].name);
										if (tmpc1[1] == '\0') {
											tmpc1[1] = ' ';
											tmpc1[2] = '\0';
										}
										if (tmpc2[1] == '\0') {
											tmpc2[1] = ' ';
											tmpc2[2] = '\0';
										}
										if (tmpc3[1] == '\0') {
											tmpc3[1] = ' ';
											tmpc3[2] = '\0';
										}
										if (tmpc4[1] == '\0') {
											tmpc4[1] = ' ';
											tmpc4[2] = '\0';
										}
										suc =
											torsion(name1, name2, name3,
													name4, tmpc1, tmpc2,
													tmpc3, tmpc4, 1);
									}
									if (suc == 0 && similarindex[i] == 1
										&& similarindex2[j] == 1
										&& similarindex2[k] == 1
										&& similarindex[l] == 1) {
										strcpy(tmpc1, similarname[i].name);
										strcpy(tmpc2,
											   similarname2[j].name);
										strcpy(tmpc3,
											   similarname2[k].name);
										strcpy(tmpc4, similarname[l].name);
										if (tmpc1[1] == '\0') {
											tmpc1[1] = ' ';
											tmpc1[2] = '\0';
										}
										if (tmpc2[1] == '\0') {
											tmpc2[1] = ' ';
											tmpc2[2] = '\0';
										}
										if (tmpc3[1] == '\0') {
											tmpc3[1] = ' ';
											tmpc3[2] = '\0';
										}
										if (tmpc4[1] == '\0') {
											tmpc4[1] = ' ';
											tmpc4[2] = '\0';
										}
										suc =
											torsion(name1, name2, name3,
													name4, tmpc1, tmpc2,
													tmpc3, tmpc4, 1);
									}
									if (suc == 0 && similarindex2[i] == 1
										&& similarindex[j] == 1
										&& similarindex[k] == 1
										&& similarindex2[l] == 1) {
										strcpy(tmpc1,
											   similarname2[i].name);
										strcpy(tmpc2, similarname[j].name);
										strcpy(tmpc3, similarname[k].name);
										strcpy(tmpc4,
											   similarname2[l].name);
										if (tmpc1[1] == '\0') {
											tmpc1[1] = ' ';
											tmpc1[2] = '\0';
										}
										if (tmpc2[1] == '\0') {
											tmpc2[1] = ' ';
											tmpc2[2] = '\0';
										}
										if (tmpc3[1] == '\0') {
											tmpc3[1] = ' ';
											tmpc3[2] = '\0';
										}
										if (tmpc4[1] == '\0') {
											tmpc4[1] = ' ';
											tmpc4[2] = '\0';
										}
										suc =
											torsion(name1, name2, name3,
													name4, tmpc1, tmpc2,
													tmpc3, tmpc4, 1);
									}
									if (suc == 0 && similarindex2[i] == 1
										&& similarindex[j] == 1
										&& similarindex2[k] == 1
										&& similarindex[l] == 1) {
										strcpy(tmpc1,
											   similarname2[i].name);
										strcpy(tmpc2, similarname[j].name);
										strcpy(tmpc3,
											   similarname2[k].name);
										strcpy(tmpc4, similarname[l].name);
										if (tmpc1[1] == '\0') {
											tmpc1[1] = ' ';
											tmpc1[2] = '\0';
										}
										if (tmpc2[1] == '\0') {
											tmpc2[1] = ' ';
											tmpc2[2] = '\0';
										}
										if (tmpc3[1] == '\0') {
											tmpc3[1] = ' ';
											tmpc3[2] = '\0';
										}
										if (tmpc4[1] == '\0') {
											tmpc4[1] = ' ';
											tmpc4[2] = '\0';
										}
										suc =
											torsion(name1, name2, name3,
													name4, tmpc1, tmpc2,
													tmpc3, tmpc4, 1);
									}
									if (suc == 0 && similarindex2[i] == 1
										&& similarindex2[j] == 1
										&& similarindex[k] == 1
										&& similarindex[l] == 1) {
										strcpy(tmpc1,
											   similarname2[i].name);
										strcpy(tmpc2,
											   similarname2[j].name);
										strcpy(tmpc3, similarname[k].name);
										strcpy(tmpc4, similarname[l].name);
										if (tmpc1[1] == '\0') {
											tmpc1[1] = ' ';
											tmpc1[2] = '\0';
										}
										if (tmpc2[1] == '\0') {
											tmpc2[1] = ' ';
											tmpc2[2] = '\0';
										}
										if (tmpc3[1] == '\0') {
											tmpc3[1] = ' ';
											tmpc3[2] = '\0';
										}
										if (tmpc4[1] == '\0') {
											tmpc4[1] = ' ';
											tmpc4[2] = '\0';
										}
										suc =
											torsion(name1, name2, name3,
													name4, tmpc1, tmpc2,
													tmpc3, tmpc4, 1);
									}
									if (suc == 0 && similarindex[i] == 1
										&& similarindex[j] == 1
										&& similarindex[k] == 1
										&& similarindex2[l] == 1) {
										strcpy(tmpc1, similarname[i].name);
										strcpy(tmpc2, similarname[j].name);
										strcpy(tmpc3, similarname[k].name);
										strcpy(tmpc4,
											   similarname2[l].name);
										if (tmpc1[1] == '\0') {
											tmpc1[1] = ' ';
											tmpc1[2] = '\0';
										}
										if (tmpc2[1] == '\0') {
											tmpc2[1] = ' ';
											tmpc2[2] = '\0';
										}
										if (tmpc3[1] == '\0') {
											tmpc3[1] = ' ';
											tmpc3[2] = '\0';
										}
										if (tmpc4[1] == '\0') {
											tmpc4[1] = ' ';
											tmpc4[2] = '\0';
										}
										suc =
											torsion(name1, name2, name3,
													name4, tmpc1, tmpc2,
													tmpc3, tmpc4, 1);
									}
									if (suc == 0 && similarindex[i] == 1
										&& similarindex[j] == 1
										&& similarindex2[k] == 1
										&& similarindex[l] == 1) {
										strcpy(tmpc1, similarname[i].name);
										strcpy(tmpc2, similarname[j].name);
										strcpy(tmpc3,
											   similarname2[k].name);
										strcpy(tmpc4, similarname[l].name);
										if (tmpc1[1] == '\0') {
											tmpc1[1] = ' ';
											tmpc1[2] = '\0';
										}
										if (tmpc2[1] == '\0') {
											tmpc2[1] = ' ';
											tmpc2[2] = '\0';
										}
										if (tmpc3[1] == '\0') {
											tmpc3[1] = ' ';
											tmpc3[2] = '\0';
										}
										if (tmpc4[1] == '\0') {
											tmpc4[1] = ' ';
											tmpc4[2] = '\0';
										}
										suc =
											torsion(name1, name2, name3,
													name4, tmpc1, tmpc2,
													tmpc3, tmpc4, 1);
									}
									if (suc == 0 && similarindex[i] == 1
										&& similarindex2[j] == 1
										&& similarindex[k] == 1
										&& similarindex[l] == 1) {
										strcpy(tmpc1, similarname[i].name);
										strcpy(tmpc2,
											   similarname2[j].name);
										strcpy(tmpc3, similarname[k].name);
										strcpy(tmpc4, similarname[l].name);
										if (tmpc1[1] == '\0') {
											tmpc1[1] = ' ';
											tmpc1[2] = '\0';
										}
										if (tmpc2[1] == '\0') {
											tmpc2[1] = ' ';
											tmpc2[2] = '\0';
										}
										if (tmpc3[1] == '\0') {
											tmpc3[1] = ' ';
											tmpc3[2] = '\0';
										}
										if (tmpc4[1] == '\0') {
											tmpc4[1] = ' ';
											tmpc4[2] = '\0';
										}
										suc =
											torsion(name1, name2, name3,
													name4, tmpc1, tmpc2,
													tmpc3, tmpc4, 1);
									}
									if (suc == 0 && similarindex2[i] == 1
										&& similarindex[j] == 1
										&& similarindex[k] == 1
										&& similarindex[l] == 1) {
										strcpy(tmpc1,
											   similarname2[i].name);
										strcpy(tmpc2, similarname[j].name);
										strcpy(tmpc3, similarname[k].name);
										strcpy(tmpc4, similarname[l].name);
										if (tmpc1[1] == '\0') {
											tmpc1[1] = ' ';
											tmpc1[2] = '\0';
										}
										if (tmpc2[1] == '\0') {
											tmpc2[1] = ' ';
											tmpc2[2] = '\0';
										}
										if (tmpc3[1] == '\0') {
											tmpc3[1] = ' ';
											tmpc3[2] = '\0';
										}
										if (tmpc4[1] == '\0') {
											tmpc4[1] = ' ';
											tmpc4[2] = '\0';
										}
										suc =
											torsion(name1, name2, name3,
													name4, tmpc1, tmpc2,
													tmpc3, tmpc4, 1);
									}
									if (suc == 0 && similarindex2[i] == 1
										&& similarindex2[j] == 1
										&& similarindex2[k] == 1
										&& similarindex2[l] == 1) {
										strcpy(tmpc1,
											   similarname2[i].name);
										strcpy(tmpc2,
											   similarname2[j].name);
										strcpy(tmpc3,
											   similarname2[k].name);
										strcpy(tmpc4,
											   similarname2[l].name);
										if (tmpc1[1] == '\0') {
											tmpc1[1] = ' ';
											tmpc1[2] = '\0';
										}
										if (tmpc2[1] == '\0') {
											tmpc2[1] = ' ';
											tmpc2[2] = '\0';
										}
										if (tmpc3[1] == '\0') {
											tmpc3[1] = ' ';
											tmpc3[2] = '\0';
										}
										if (tmpc4[1] == '\0') {
											tmpc4[1] = ' ';
											tmpc4[2] = '\0';
										}
										suc =
											torsion(name1, name2, name3,
													name4, tmpc1, tmpc2,
													tmpc3, tmpc4, 1);
									}
									/* end */

									if (suc == 0) {
										strcpy(tmpc1, atom[i].ambername);
										strcpy(tmpc2, atom[j].ambername);
										strcpy(tmpc3, atom[k].ambername);
										strcpy(tmpc4, atom[l].ambername);

										if (tmpc1[1] == '\0') {
											tmpc1[1] = ' ';
											tmpc1[2] = '\0';
										}
										if (tmpc2[1] == '\0') {
											tmpc2[1] = ' ';
											tmpc2[2] = '\0';
										}
										if (tmpc3[1] == '\0') {
											tmpc3[1] = ' ';
											tmpc3[2] = '\0';
										}
										if (tmpc4[1] == '\0') {
											tmpc4[1] = ' ';
											tmpc4[2] = '\0';
										}

										strcpy(name1, tmpc1);
										strcpy(name2, tmpc2);
										strcpy(name3, tmpc3);
										strcpy(name4, tmpc4);
									}
									if (suc == 0 && corrindex[i] == 1) {
										strcpy(tmpc1, corrname[i].name);
										if (tmpc1[1] == '\0') {
											tmpc1[1] = ' ';
											tmpc1[2] = '\0';
										}
										suc =
											torsion(name1, name2, name3,
													name4, tmpc1, tmpc2,
													tmpc3, tmpc4, 1);
									}

									if (suc == 0 && corrindex[j] == 1) {
										strcpy(tmpc1, atom[i].ambername);
										strcpy(tmpc2, corrname[j].name);
										if (tmpc1[1] == '\0') {
											tmpc1[1] = ' ';
											tmpc1[2] = '\0';
										}
										if (tmpc2[1] == '\0') {
											tmpc2[1] = ' ';
											tmpc2[2] = '\0';
										}
										suc =
											torsion(name1, name2, name3,
													name4, tmpc1, tmpc2,
													tmpc3, tmpc4, 1);
									}
									if (suc == 0 && corrindex[k] == 1) {
										strcpy(tmpc2, atom[j].ambername);
										strcpy(tmpc3, corrname[k].name);
										if (tmpc2[1] == '\0') {
											tmpc2[1] = ' ';
											tmpc2[2] = '\0';
										}
										if (tmpc3[1] == '\0') {
											tmpc3[1] = ' ';
											tmpc3[2] = '\0';
										}
										suc =
											torsion(name1, name2, name3,
													name4, tmpc1, tmpc2,
													tmpc3, tmpc4, 1);
									}

									if (suc == 0 && corrindex[l] == 1) {
										strcpy(tmpc3, atom[k].ambername);
										strcpy(tmpc4, corrname[l].name);
										if (tmpc3[1] == '\0') {
											tmpc3[1] = ' ';
											tmpc3[2] = '\0';
										}
										if (tmpc4[1] == '\0') {
											tmpc4[1] = ' ';
											tmpc4[2] = '\0';
										}
										suc =
											torsion(name1, name2, name3,
													name4, tmpc1, tmpc2,
													tmpc3, tmpc4, 1);
									}

									if (suc == 0 && corrindex[i] == 1
										&& corrindex[j] == 1
										&& corrindex[k] == 1
										&& corrindex[l] == 1) {
										strcpy(tmpc1, corrname[i].name);
										strcpy(tmpc2, corrname[j].name);
										strcpy(tmpc3, corrname[k].name);
										strcpy(tmpc4, corrname[l].name);
										if (tmpc1[1] == '\0') {
											tmpc1[1] = ' ';
											tmpc1[2] = '\0';
										}
										if (tmpc2[1] == '\0') {
											tmpc2[1] = ' ';
											tmpc2[2] = '\0';
										}
										if (tmpc3[1] == '\0') {
											tmpc3[1] = ' ';
											tmpc3[2] = '\0';
										}
										if (tmpc4[1] == '\0') {
											tmpc4[1] = ' ';
											tmpc4[2] = '\0';
										}
										suc =
											torsion(name1, name2, name3,
													name4, tmpc1, tmpc2,
													tmpc3, tmpc4, 1);
									}

									if (suc == 0) {

									   	fprintf(fpout,
							                                "%2s %2s %2s %2s%4d%9.3lf%14.3lf%16.3lf",
							                                name1, name2, name3, name4, 1,
							                                0.0, 0.0, 0.0);
					                                        fprintf(fpout, "      ! %s\n", "ATTN, need revision");

										strcpy(torsionparm[torsionparmnum].
											   name1, name1);
										strcpy(torsionparm[torsionparmnum].
											   name2, name2);
										strcpy(torsionparm[torsionparmnum].
											   name3, name3);
										strcpy(torsionparm[torsionparmnum].
											   name4, name4);
										torsionparmnum++;
										if (torsionparmnum >=
											maxtorsionparm) {
											maxtorsionparm +=
												MAX_FF_TORSION;
											torsionparm = (TORSION *)
												realloc(torsionparm,
														sizeof(TORSION) *
														maxtorsionparm);
											if (torsionparm == NULL) {
												fprintf(stdout,
														"memory allocation error for *torsionparm\n");
												exit(1);
											}
										}
									}
								}
						}
					}
                                }
			}
		}
	}
}


void wrt_improper(void)
{
	int i, j;
	int suc;
	int tmpflag = 0;
        char tmpc[5], tmpn[5];
	char tmpc1[5], tmpc2[5], tmpc3[5], tmpc4[5];
	char tmpn1[5], tmpn2[5], tmpn3[5], tmpn4[5];
	char name1[5], name2[5], name3[5], name4[5];

	for (i = 0; i < impropernum; i++) {
		suc = 0;
                if (!tmpflag) {
	                fprintf(fpout, "\n");
        		fprintf(ffout, "\n%s\n", "IMPHI");
			tmpflag = 1;
    	        }
		strcpy(tmpc1, atom[improper[i].atid1].ambername);
		strcpy(tmpc2, atom[improper[i].atid2].ambername);
		strcpy(tmpc3, atom[improper[i].atid3].ambername);
		strcpy(tmpc4, atom[improper[i].atid4].ambername);

           	strcpy(tmpn1, atom[improper[i].atid1].name);
		strcpy(tmpn2, atom[improper[i].atid2].name);
		strcpy(tmpn3, atom[improper[i].atid3].name);
		strcpy(tmpn4, atom[improper[i].atid4].name);

		if (tmpc1[1] == '\0') {
			tmpc1[1] = ' ';
			tmpc1[2] = '\0';
		}
		if (tmpc2[1] == '\0') {
			tmpc2[1] = ' ';
			tmpc2[2] = '\0';
		}
		if (tmpc3[1] == '\0') {
			tmpc3[1] = ' ';
			tmpc3[2] = '\0';
		}
		if (tmpc4[1] == '\0') {
			tmpc4[1] = ' ';
			tmpc4[2] = '\0';
		}
		
                if(strcmp(tmpc1, tmpc2) > 0)  {
       	           strcpy(tmpc, tmpc2);
       	           strcpy(tmpn, tmpn2);
       	           strcpy(tmpc2, tmpc1);
       	           strcpy(tmpn2, tmpn1);
                   strcpy(tmpc1, tmpc);
                   strcpy(tmpn1, tmpn);
                }
                if(strcmp(tmpc1, tmpc4) > 0)  {
                   strcpy(tmpc, tmpc4);
                   strcpy(tmpn, tmpn4);
                   strcpy(tmpc4, tmpc1);
                   strcpy(tmpn4, tmpn1);
                   strcpy(tmpc1, tmpc);
                   strcpy(tmpn1, tmpn);		   
                }
                if(strcmp(tmpc2, tmpc4) > 0)  {
                   strcpy(tmpc, tmpc4);
		   strcpy(tmpn, tmpn4);
                   strcpy(tmpc4, tmpc2);
		   strcpy(tmpn4, tmpn2);
                   strcpy(tmpc2, tmpc);
		   strcpy(tmpn2, tmpn);
                }
		
		strcpy(name1, tmpc1);
		strcpy(name2, tmpc2);
		strcpy(name3, tmpc3);
		strcpy(name4, tmpc4);


		for (j = 0; j < improperparmnum; j++)
                    if (((improperparm[j].name1[0] == tmpc1[0] && improperparm[j].name1[1] == tmpc1[1]) ||
                         (improperparm[j].name1[0] == 'X'      && improperparm[j].name1[1] == ' '))          &&
                        ((improperparm[j].name2[0] == tmpc2[0] && improperparm[j].name2[1] == tmpc2[1]) ||
                         (improperparm[j].name2[0] == 'X'      && improperparm[j].name2[1] == ' '))          && 
                        ((improperparm[j].name3[0] == tmpc3[0] && improperparm[j].name3[1] == tmpc3[1]) ||
                         (improperparm[j].name3[0] == 'X' && improperparm[j].name3[1] == ' '))               && 
                        ((improperparm[j].name4[0] == tmpc4[0] && improperparm[j].name4[1] == tmpc4[1]) || 
                         (improperparm[j].name4[0] == 'X' && improperparm[j].name4[1] == ' '))) {
                                    fprintf(fpout, "IMPH %-5s %-5s %-5s %-5s\n",
 	                                    tmpn1, tmpn2, tmpn3, tmpn4);
                                    if (improperparm[j].flag == 0) {
                                                fprintf(ffout, "%2s %2s %2s %2s%10.3lf%10.0lf%10.1lf\n",
                                                	improperparm[j].name1,
						        improperparm[j].name2,
					  	        improperparm[j].name3,
						        improperparm[j].name4,
						        improperparm[j].force,
						        improperparm[j].fterm,
						        improperparm[j].phase);
					        improperparm[j].flag = 1;
				    }
				    suc = 1;
				    break;
                    }
		if (suc == 0) {
                        fprintf(fpout, "IMPH %-5s %-5s %-5s %-5s\n",
		                tmpn1, tmpn2, tmpn3, tmpn4);
			fprintf(ffout, "%2s %2s %2s %2s%10.3lf%10.0lf%10.1lf",
				name1, name2, name3, name4, 1.1, 2.0, 180.0);
			fprintf(ffout, "	! Using default value\n");
			improperparm[improperparmnum].flag = 1;
   		        strcpy(improperparm[improperparmnum].name1, name1);
			strcpy(improperparm[improperparmnum].name2, name2);
			strcpy(improperparm[improperparmnum].name3, name3);
			strcpy(improperparm[improperparmnum].name4, name4);
			improperparmnum++;
			if (improperparmnum >= maximproperparm) {
				maximproperparm += MAX_FF_IMPROPER;
				improperparm =
					(IMPROPER *) realloc(improperparm,
										 sizeof(IMPROPER) *
										 maximproperparm);
				if (improperparm == NULL) {
					fprintf(stdout,
							"memory allocation error for *improperparm\n");
					exit(1);
				}
			}
		}
	}
}

void wrt_charmm()
{
	int i;
	fprintf(fhout, "* %s main file generated by charmmgen\n",  resname);
	fprintf(fhout, "* author: Victor Bazterra 2004 CHPC\n");
	fprintf(fhout, "* \n\n");
	fprintf(fhout, "Read rtf card name %s\n\n", ofilename);
	fprintf(fhout, "Read param card name %s\n\n", ffilename);
	fprintf(fhout, "Read sequence card\n");
	fprintf(fhout, "* Reading the coordinates of residue\n");
	fprintf(fhout, "* \n");
	fprintf(fhout, "    1\n");
	fprintf(fhout, "%3s\n", resname);
	fprintf(fhout, "Generate %3s setup\n", resname);
	fprintf(fhout, "Read coor card free\n");
	fprintf(fhout,"* Residues coordinate\n");
	fprintf(fhout,"* \n");

	fprintf(fhout,"%5d\n", atomnum);
	for (i = 0; i < atomnum; i++)
		fprintf(fhout,"%5d%5d %-4s %-4s%10.5lf%10.5lf%10.5lf\n",
			i+1, 1, resname, atom[i].name, atom[i].x, atom[i].y, atom[i].z);

	fprintf(fhout,"Energy ihbfrq 0 inbfrq -1\n");
 	fprintf(fhout,"Open write card unit 14 name tmp.out\n");
	fprintf(fhout,"Write energy card unit 14\n");
	fprintf(fhout,"* %s energy before the optimization\n", resname);
	fprintf(fhout,"* \n");
	fprintf(fhout,"Minimise Abnr Nstep 5000 Tolgrd 0.00001\n");
	fprintf(fhout,"Write energy card unit 14\n");
	fprintf(fhout,"* %s energy after optimization\n", resname);
	fprintf(fhout,"* \n");
	fprintf(fhout,"Write coor card unit 14\n");
	fprintf(fhout,"* %s energy after optimization\n", resname);
	fprintf(fhout,"* \n");
}

void wrt_residue()
{
	double tcharge = 0;

	for (i = 0; i < atomnum; i++)
		tcharge += atom[i].charge;

      	fprintf(fpout, "\nRESI %-4s %2.3f\n", resname, tcharge);
	fprintf(fpout, "GROUP\n");
	for (i = 0; i < atomnum; i++)
		fprintf(fpout,"ATOM %-5s %-5s %10.6lf\n",
			atom[i].name, atom[i].ambername, atom[i].charge);
}

int main(int argc, char *argv[])
{
	int i;
	int format;
	int overflow_flag = 0;			/*if overflow_flag ==1, reallocate memory */

	default_cinfo(&cinfo);
	default_minfo(&minfo);
    amberhome = (char *) getenv("AMBERHOME");
    if( amberhome == NULL ){
       fprintf( stdout, "AMBERHOME is not set!\n" );
       exit(1);
    }
	build_dat_path(minfo.connect_file, "CONNECT.TPL",
		sizeof minfo.connect_file, 0);

	build_path(pfilename, "/dat/leap/parm/", "gaff.dat",
		sizeof pfilename, 0);

	if (strcmp(COLORTEXT, "YES") == 0 || strcmp(COLORTEXT, "yes") == 0) {
		if (argc == 2
			&& (strcmp(argv[1], "-h") == 0
				|| strcmp(argv[1], "-H") == 0)) {
			printf("[31mUsage: charmmgen -i[0m input file name\n"
				   "[31m                 -o[0m frcmod file name\n"
				   "[31m                 -f[0m input file format (prepi, prepc, ac ,mol2) \n"
				   "[31m                 -p[0m ff parmfile\n"
				   "[31m                 -r[0m residue name\n"
				   "[31m                 -c[0m atom type correspondening file, default is ATCOR.DAT\n");
			exit(1);
		}
		if (argc != 7 && argc != 9 && argc != 11 && argc != 13) {
			printf("[31mUsage: charmmgen -i[0m input file name\n"
				   "[31m                 -o[0m frcmod file name\n"
				   "[31m                 -f[0m input file format (prepi, prepc, ac ,mol2) \n"
				   "[31m                 -p[0m ff parmfile\n"
				   "[31m                 -r[0m residue name\n"
				   "[31m                 -c[0m atom type correspondening file, default is ATCOR.DAT\n");
			exit(1);
		}
	} else {
		if (argc == 2
			&& (strcmp(argv[1], "-h") == 0
				|| strcmp(argv[1], "-H") == 0)) {
			printf("Usage: charmmgen -i   input file name\n");
			printf("                 -o   frcmod file name\n");
			printf("                 -f   input file format (prepi, prepc, ac) \n");
			printf("                 -p   ff parmfile\n");
			printf("                 -r   residue name \n");
			printf("                 -c   atom type correspondence file \n");
			printf("                      (default is ATCOR.DAT)\n");
			exit(1);
		}
		if (argc != 7 && argc != 9 && argc != 11 && argc != 13) {
			printf("Usage: charmmgen -i   input file name\n");
			printf("                 -o   frcmod file name\n");
			printf("                 -f   input file format (prepi, prepc, ac) \n");
			printf("                 -p   ff parmfile\n");
			printf("                 -r   residue name \n");
			printf("                 -c   atom type correspondence file \n");
			printf("                      (default is ATCOR.DAT)\n");
			exit(1);
		}
	}
	format = 0;
	cindex = 0;
	for (i = 1; i < argc; i += 2) {
		if (strcmp(argv[i], "-i") == 0)
			strcpy(ifilename, argv[i + 1]);
		if (strcmp(argv[i], "-p") == 0)
			strcpy(pfilename, argv[i + 1]);
		if (strcmp(argv[i], "-o") == 0) {
       			strcpy(ofilename, argv[i + 1]);
			strcpy(ffilename, argv[i + 1]);
			strcpy(hfilename, argv[i + 1]);
                        strcat(ofilename, ".rtf");
                        strcat(ffilename, ".prm");
                        strcat(hfilename, ".inp");
		}
                if (strcmp(argv[i], "-c") == 0) {
			strcpy(cfilename, argv[i + 1]);
			cindex = 1;
		}
                if (strcmp(argv[i], "-r") == 0) 
			strcpy(resname, argv[i + 1]);
		if (strcmp(argv[i], "-f") == 0) {
			if (strcmp(argv[i + 1], "prepi") == 0)
				format = 0;
			if (strcmp(argv[i + 1], "prepc") == 0)
				format = 1;
			if (strcmp(argv[i + 1], "ac") == 0)
				format = 2;
			if (strcmp(argv[i + 1], "mol2") == 0)
				format = 3;
		}
	}

	if (cindex == 0) {
		build_dat_path(cfilename, "ATCOR.DAT", sizeof cfilename, 0);
		cindex = 1;
	}
	if (cindex == 0)
		if ((fpout = fopen("ATCOR.DAT", "r")) != NULL) {
			strcpy(cfilename, "ATCOR.DAT");
			cindex = 1;
		}

	/*  memory allocation */
	/* initialize */
	maxcorr = MAX_FF_CORR;
	maxatomtype = MAX_FF_ATOMTYPE;
	maxvdwparm = MAX_FF_VDW;
	maxbondparm = MAX_FF_BOND;
	maxangleparm = MAX_FF_ANGLE;
	maxtorsionparm = MAX_FF_TORSION;
	maximproperparm = MAX_FF_IMPROPER;

	atomtypenum = 0;
	vdwparmnum = 0;
	bondparmnum = 0;
	angleparmnum = 0;
	torsionparmnum = 0;
	improperparmnum = 0;

	if(strlen(resname) < 1) strcpy(resname, "MOL");

	/*read in prep or ac file */
	atom = (ATOM *) malloc(sizeof(ATOM) * cinfo.maxatom);
	if (atom == NULL) {
		fprintf(stdout, "memory allocation error for *atom\n");
		exit(1);
	}
	bond_array = (BOND *) malloc(sizeof(BOND) * cinfo.maxbond);
	if (bond_array == NULL) {
		fprintf(stdout, "memory allocation error for *bond_array\n");
		exit(1);
	}
	for (i = 0; i < cinfo.maxbond; ++i) {
		bond_array[i].jflag = -1; /* bond type has not been assigned */
	}

	if (format == 0)
		overflow_flag = rprepi(ifilename, &atomnum, atom, &bondnum, bond_array, &cinfo, &minfo);
	if (format == 1)
		overflow_flag = rprepc(ifilename, &atomnum, atom, &cinfo, &minfo);
	if (format == 2)
		overflow_flag =
			rac(ifilename, &atomnum, atom, &bondnum, bond_array, &cinfo,
				&minfo);
	if (format == 3)
		overflow_flag =
			rmol2(ifilename, &atomnum, atom, &bondnum, bond_array, &cinfo,
				  &minfo, 1);
	if (overflow_flag) {
		cinfo.maxatom = atomnum + 10;
		cinfo.maxbond = bondnum + 10;
		free(atom);
		free(bond_array);
		atom = (ATOM *) malloc(sizeof(ATOM) * cinfo.maxatom);
		if (atom == NULL) {
			fprintf(stdout, "memory allocation error for *atom\n");
			exit(1);
		}
		bond_array = (BOND *) malloc(sizeof(BOND) * cinfo.maxbond);
		if (bond_array == NULL) {
			fprintf(stdout, "memory allocation error for *bond_array\n");
			exit(1);
		}
		int i;
		for (i = 0; i < cinfo.maxbond; ++i) {
			bond_array[i].jflag = -1; /* bond type has not been assigned */
		}
		if (format == 0)
			overflow_flag =
				rprepi(ifilename, &atomnum, atom, &bondnum, bond_array, &cinfo, &minfo);
		if (format == 1)
			overflow_flag =
				rprepc(ifilename, &atomnum, atom, &cinfo, &minfo);
		if (format == 2)
			overflow_flag =
				rac(ifilename, &atomnum, atom, &bondnum, bond_array, &cinfo,
					&minfo);
		if (format == 3)
			overflow_flag =
				rmol2(ifilename, &atomnum, atom, &bondnum, bond_array,
					  &cinfo, &minfo, 1);

	}

	if (format == 0) {
		atomicnum(atomnum, atom);
		adjustatomname(atomnum, atom, 1);
	}
	if (format == 1) {
		atomicnum(atomnum, atom);
		overflow_flag =
			connect(minfo.connect_file, atomnum, atom, &bondnum,
					bond_array, cinfo.maxbond);
		adjustatomname(atomnum, atom, 1);
	}
	if (format == 2) {
		atomicnum(atomnum, atom);
/*
		overflow_flag =
			connect(minfo.connect_file, atomnum, atom, &bondnum,
					bond_array, cinfo.maxbond);
*/
	}
	if (format == 3) {
		atomicnum(atomnum, atom);
		default_inf(atomnum, atom, 0);
	}
	corrname = (NAME *) malloc(sizeof(NAME) * atomnum);
	if (corrname == NULL) {
		fprintf(stdout, "memory allocation error for *corrname\n");
		exit(1);
	}
	similarname = (NAME *) malloc(sizeof(NAME) * atomnum);
	if (similarname == NULL) {
		fprintf(stdout, "memory allocation error for *similarname\n");
		exit(1);
	}
	similarname2 = (NAME *) malloc(sizeof(NAME) * atomnum);
	if (similarname2 == NULL) {
		fprintf(stdout, "memory allocation error for *similarname2\n");
		exit(1);
	}
	corrindex = (int *) malloc(sizeof(int) * atomnum);
	if (corrindex == NULL) {
		fprintf(stdout, "memory allocation error for *corrindex\n");
		exit(1);
	}
	similarindex = (int *) malloc(sizeof(int) * atomnum);
	if (similarindex == NULL) {
		fprintf(stdout, "memory allocation error for *similarindex\n");
		exit(1);
	}
	similarindex2 = (int *) malloc(sizeof(int) * atomnum);
	if (similarindex2 == NULL) {
		fprintf(stdout, "memory allocation error for *similarindex2\n");
		exit(1);
	}

	corr = (CORR *) calloc(maxcorr, sizeof(CORR));
	if (corr == NULL) {
		fprintf(stdout, "memory allocation error for *corr\n");
		exit(1);
	}
	improperindex = (int *) calloc(maxcorr, sizeof(int));
	if (improperindex == NULL) {
		fprintf(stdout, "memory allocation error for *improperindex\n");
		exit(1);
	}

	atomtype = (ATOMTYPE *) calloc(maxatomtype, sizeof(ATOMTYPE));
	if (atomtype == NULL) {
		fprintf(stdout, "memory allocation error for *atomtype\n");
		exit(1);
	}
	bondparm = (BOND_PARAM *) calloc(maxbondparm, sizeof(BOND_PARAM));
	if (bondparm == NULL) {
		fprintf(stdout, "memory allocation error for *bondparm\n");
		exit(1);
	}

	angleparm = (ANGLE *) calloc(maxangleparm, sizeof(ANGLE));
	if (angleparm == NULL) {
		fprintf(stdout, "memory allocation error for *angleparm\n");
		exit(1);
	}
	torsionparm = (TORSION *) calloc(maxtorsionparm, sizeof(TORSION));
	if (torsionparm == NULL) {
		fprintf(stdout, "memory allocation error for *torsionparm\n");
		exit(1);
	}
	improperparm = (IMPROPER *) calloc(maximproperparm, sizeof(IMPROPER));
	if (improperparm == NULL) {
		fprintf(stdout, "memory allocation error for *improperparm\n");
		exit(1);
	}

	vdwparm = (VDW *) calloc(maxvdwparm, sizeof(VDW));
	if (vdwparm == NULL) {
		fprintf(stdout, "memory allocation error for *vdwparm\n");
		exit(1);
	}

	improper = (IMPROPERID *) malloc(sizeof(IMPROPERID) * maxcorr);
	if (improper == NULL) {
		fprintf(stdout, "memory allocation error for *similarname2\n");
		exit(1);
	}


	if ((fpout = fopen(ofilename, "w")) == NULL) {
		fprintf(stdout, "Cannot open the ofilename file %s to write, exit\n", ofilename);
		exit(1);
	}

      	if ((ffout = fopen(ffilename, "w")) == NULL) {
		fprintf(stdout, "Cannot open the ffilename file %s to write, exit\n", ofilename);
		exit(1);
	}

      	if ((fhout = fopen(hfilename, "w")) == NULL) {
		fprintf(stdout, "Cannot open the hfilename file %s to write, exit\n", ofilename);
		exit(1);
	}


      	fprintf(fpout, "* Topology File.\n");
	fprintf(fpout, "* \n");
	fprintf(fpout, "   99   1\n");
	fprintf(ffout, "* Force Field Parameter File.\n");
	fprintf(ffout, "* \n");


	/* read in parameters */

	if (cindex == 1)
		readcorr(cfilename);	/*atom type *corresponding file */
	readparm(pfilename);		/*principle parameter file */

	if (format == 0 || format == 1)
		improper_id1(ifilename);
	if (format == 2 || format == 3)
		improper_id2();

	prepare();
        wrt_atomtype();
        wrt_residue();
	wrt_bond();
	wrt_angle();
	wrt_torsion();
	wrt_improper();
	wrt_vdw();
        wrt_charmm();

	fclose(fpout);
	fclose(ffout);
	fclose(fhout);

	free(atom);
	free(bond_array);
	free(corrname);
	free(similarname);
	free(similarname2);
	free(corrindex);
	free(similarindex);
	free(similarindex2);
	free(improper);
	free(improperindex);
	free(corr);
	free(atomtype);
	free(bondparm);
	free(angleparm);
	free(torsionparm);
	free(improperparm);
	free(vdwparm);

	return (0);
}
