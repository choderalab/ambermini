/*
************************************************************************
*           All Copyright Reserved!                                    *
*                                                                      *
*  Prog:    parmchk                                                    *
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
char *amberhome;
# include "common.h"
# include "define.h"
# include "atom.h"
# include "utility.c"
# include "common.c"
# include "rotate.c"
# include "ac.c"
# include "mol2.c"
# include "prep.c"

#define MAXPARM 250   /*maximum atom types */
#define MAXEQUA 25    /*maximum equal atom types for each atom type */
#define MAXCORR 50    /*maximum corresponding atom types for each atom type */
#define MAX_FF_ATOMTYPE 250
#define MAX_FF_VDW 250
#define MAX_FF_BOND 2000
#define MAX_FF_ANGLE 5000
#define MAX_FF_TORSION 1500
#define MAX_FF_IMPROPER 500
#define MAX_EQU_VDW 20
#define MAX_EQU_VDW_NUM 50
#define INITSCORE       999999
#define debug 0

typedef struct {
	int atid1;
	int atid2;
	int atid3;
	int atid4;
} IMPROPERID;

typedef struct {
	char atomtype[10];
	int  pid;
} EQUA;

typedef struct {
	char atomtype[10];
	int  pid;
	int  type;
	double bl;  	/*penalty of bond length */
	double blf; 	/*penalty of bond stretching force constant */
	double ba;  	/*penalty of bond angle */
	double baf; 	/*penalty of bond bending force constant */
	double cba; 	/*penalty of bond angle when replace angle central atom */
	double cbaf;	/*penalty of bond bending force constant of replacing angle central atom */
	double ctor;	/*penalty of torsional angle when replace central atoms */
	double tor; 	/*penalty of torsional angle when replace torsional atoms */
	double ps;      /*overall penalty score*/
	double improper;/* penalty of improper angle when replace improper atoms */
} CORR;

typedef struct {
	char atomtype[10];
	int  group;
	int  equtype;  /*equivalent atom type, cc/ce/cg/nc/ne/pc/pe belong to TYPE 1, cd/cf/ch/nd/nf/pd/pf belong to TYPE 2 and others are 0*/
	int  improper;
	int  ncorr;
	int  nequa;
	double mass;
	CORR corr[MAXCORR];
	EQUA equa[MAXEQUA];
} PARM;

typedef struct {
	char name[10];
	double mass;
	double pol;
} ATOMTYPE;

typedef struct {
	char name1[10];
	char name2[10];
	char name3[10];
	double angle;
	double force;
} ANGLE;

typedef struct {
	char name1[10];
	char name2[10];
	char name3[10];
	char name4[10];
	int mul;
	double fterm;
	double phase;
	double force;
} TORSION;

typedef struct {
	char name1[10];
	char name2[10];
	char name3[10];
	char name4[10];
	int mul;
	int numX;
	double fterm;
	double phase;
	double force;
} IMPROPER;

typedef struct {
	char name[10];
	double radius;
	double pot;
} VDW;

typedef struct {
	char name[MAX_EQU_VDW_NUM][10];
	int num;
} EQU_VDW;

typedef struct {
	double BL;
	double BLF;
	double BA;
	double BAF;
	double X;
	double X3;
	double BA_CTR;
	double TOR_CTR;
	double IMPROPER;
	double GROUP;
	double EQUTYPE;
} WT;

typedef struct {
	double BL;
	double BLF;
	double BA;
	double BAF;
	double BA_CTR;
	double BAF_CTR;
	double TOR;
	double TOR_CTR;
	double FRACT1;
	double FRACT2;
} DEFAULT_VALUE;

MOLINFO minfo;
CONTROLINFO cinfo;
ATOM *atom;
BOND *bond_array;
WT   wt;
DEFAULT_VALUE dv;

int atomnum = 0;
int bondnum = 0;
IMPROPERID *improper;
char ifilename[MAXCHAR];
char ofilename[MAXCHAR];
char pfilename[MAXCHAR];
char cfilename[MAXCHAR];
char line[MAXCHAR];
int  cindex = 0;
int impropernum = 0;
int output_improper_flag = 1;
int i, j, k, l;
int maxatomtype = 0;
int maxvdwparm = 0;
int maxbondparm = 0;
int maxangleparm = 0;
int maxtorsionparm = 0;
int maximproperparm = 0;
int atomtypenum = 0;
int parmnum = 0;
int maxparmnum;
int vdwparmnum = 0;
int bondparmnum = 0;
int angleparmnum = 0;
int torsionparmnum = 0;
int improperparmnum = 0;

double THRESHOLD_BA; 
double equtype_penalty_score;
/* 
  atomtypenum2 etc are numbers of parameters from force field files
  not including those newly found in the file
*/
int atomtypenum2 = 0; 
int vdwparmnum2 = 0;
int bondparmnum2 = 0;
int angleparmnum2 = 0;
int torsionparmnum2 = 0;
int improperparmnum2 = 0;
/*H-1, C-2, N-3, O-4, F-5, Cl-6, Br-7, I-8, S-9, P-10*/

int *parmid; 
PARM *parm;
ATOMTYPE *atomtype;
BOND_FF *bondparm;
ANGLE *angleparm;
TORSION *torsionparm;
IMPROPER *improperparm;
VDW *vdwparm;
/* for equivalent vdw types */
EQU_VDW equ_vdw[MAX_EQU_VDW];
int equ_vdw_num = 0;

FILE *fp, *fpout;

int allparm_flag = 0;
int pformat = 1;

ANGLE   bestba;
int	bestblid = -1;
int 	bestbaid = -1;
int 	besttorid= -1;
int 	bestimproperid=-1;
double bestscore;

/* set up improper atoms for prepi/prepc file */
void improper_id1(char *filename) {
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

/* set up improper atoms for ac or mol2 file */
void improper_id2() {
	int i;
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

void readfrcmod(char *filename) {
int mindex = 0;
int bindex = 0;
int aindex = 0;
int tindex = 0;
int iindex = 0;
int vindex = 0;
int pos_tor = 0;
int tmpnum  = 0;
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
		if (strncmp(line, "MASS", 4) == 0) {
			mindex = 1;
			continue;
		}
		if (strncmp(line, "BOND", 4) == 0) { 
			bindex = 1;
			continue;
		}
		if (strncmp(line, "ANGLE", 5) == 0){ 
			aindex = 1;
			continue;
		}
		if (strncmp(line, "DIHE", 4) == 0) { 
			pos_tor= ftell(fp);
			tindex = 1;
			continue;
		}
		if (strncmp(line, "IMPROPER", 8) == 0) { 
			iindex = 1;
			continue;
		}
		if (strncmp(line, "NONBON", 6) == 0) { 
			vindex = 1;
			continue;
		}
		if (strlen(line) <= 2) {
			if(mindex == 1) mindex = 0;
			if(bindex == 1) bindex = 0;
			if(aindex == 1) aindex = 0;
			if(tindex == 1) {
				fseek(fp,pos_tor,0);
				tindex = 2;
				continue;
			}
			if(tindex == 2) tindex = 0;
			if(iindex == 1) iindex = 0;
			if(vindex == 1) vindex = 0;
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
			if(line[1]== ' ') 
                        	bondparm[bondparmnum].name1[1] = '\0';
			else {
				bondparm[bondparmnum].name1[1] = line[1];
				bondparm[bondparmnum].name1[2] = '\0';
			}
                        bondparm[bondparmnum].name2[0] = line[3];
			if(line[4]== ' ') 
                        	bondparm[bondparmnum].name2[1] = '\0';
			else {
				bondparm[bondparmnum].name2[1] = line[4];
				bondparm[bondparmnum].name2[2] = '\0';
			}
				
                        sscanf(&line[5], "%lf%lf", &bondparm[bondparmnum].force,
                                   &bondparm[bondparmnum].length);
                        bondparmnum++;
                        if (bondparmnum >= maxbondparm) {
                                maxbondparm += MAX_FF_BOND;
                                bondparm =
                                        (BOND_FF *) realloc(bondparm,
                                                                                sizeof(BOND_FF) * maxbondparm);
                                if (bondparm == NULL) {
                                        fprintf(stdout,
                                                        "memory allocation error for *bondparm\n");
                                        exit(1);
                                }
                        }
                }
                if (aindex == 1) {
                        angleparm[angleparmnum].name1[0] = line[0];
			if(line[1]== ' ') 
                        	angleparm[angleparmnum].name1[1] = '\0';
			else {
				angleparm[angleparmnum].name1[1] = line[1];
				angleparm[angleparmnum].name1[2] = '\0';
			}
                        angleparm[angleparmnum].name2[0] = line[3];
			if(line[4]== ' ') 
                        	angleparm[angleparmnum].name2[1] = '\0';
			else {
				angleparm[angleparmnum].name2[1] = line[4];
				angleparm[angleparmnum].name2[2] = '\0';
			}
                        angleparm[angleparmnum].name3[0] = line[6];
			if(line[7]== ' ') 
                        	angleparm[angleparmnum].name3[1] = '\0';
			else {
				angleparm[angleparmnum].name3[1] = line[7];
				angleparm[angleparmnum].name3[2] = '\0';
			}
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
                if (tindex == 1) { /*we first only read special torsional angle parameters*/
                        if(line[0] == 'X' || line[3] == 'X' || line[6] == 'X' || line[9] == 'X')
                                continue;
                        torsionparm[torsionparmnum].name1[0] = line[0];
			if(line[1]== ' ') 
                        	torsionparm[torsionparmnum].name1[1] = '\0';
			else {
				torsionparm[torsionparmnum].name1[1] = line[1];
				torsionparm[torsionparmnum].name1[2] = '\0';
			}
                        torsionparm[torsionparmnum].name2[0] = line[3];
			if(line[4]== ' ') 
                        	torsionparm[torsionparmnum].name2[1] = '\0';
			else {
				torsionparm[torsionparmnum].name2[1] = line[4];
				torsionparm[torsionparmnum].name2[2] = '\0';
			}
                        torsionparm[torsionparmnum].name3[0] = line[6];
			if(line[7]== ' ') 
                        	torsionparm[torsionparmnum].name3[1] = '\0';
			else {
				torsionparm[torsionparmnum].name3[1] = line[7];
				torsionparm[torsionparmnum].name3[2] = '\0';
			}
                        torsionparm[torsionparmnum].name4[0] = line[9];
			if(line[10]== ' ') 
                        	torsionparm[torsionparmnum].name4[1] = '\0';
			else {
				torsionparm[torsionparmnum].name4[1] = line[10];
				torsionparm[torsionparmnum].name4[2] = '\0';
			}
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
                if (tindex == 2) {
                        if(line[0] != 'X' && line[3] != 'X' && line[6] != 'X' && line[9] != 'X')
                                continue;
                        torsionparm[torsionparmnum].name1[0] = line[0];
			if(line[1]== ' ') 
                        	torsionparm[torsionparmnum].name1[1] = '\0';
			else {
				torsionparm[torsionparmnum].name1[1] = line[1];
				torsionparm[torsionparmnum].name1[2] = '\0';
			}
                        torsionparm[torsionparmnum].name2[0] = line[3];
			if(line[4]== ' ') 
                        	torsionparm[torsionparmnum].name2[1] = '\0';
			else {
				torsionparm[torsionparmnum].name2[1] = line[4];
				torsionparm[torsionparmnum].name2[2] = '\0';
			}
                        torsionparm[torsionparmnum].name3[0] = line[6];
			if(line[7]== ' ') 
                        	torsionparm[torsionparmnum].name3[1] = '\0';
			else {
				torsionparm[torsionparmnum].name3[1] = line[7];
				torsionparm[torsionparmnum].name3[2] = '\0';
			}
                        torsionparm[torsionparmnum].name4[0] = line[9];
			if(line[10]== ' ') 
                        	torsionparm[torsionparmnum].name4[1] = '\0';
			else {
				torsionparm[torsionparmnum].name4[1] = line[10];
				torsionparm[torsionparmnum].name4[2] = '\0';
			}
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
                        tmpnum = 0;
                        tmpc1[0] = line[0];
			if(line[1] == ' ')
                        	tmpc1[1] = '\0';
			else {
                        	tmpc1[1] = line[1];
                        	tmpc1[2] = '\0';
			}
                        tmpc2[0] = line[3];
			if(line[4] == ' ')
                        	tmpc2[1] = '\0';
			else {
                        	tmpc2[1] = line[4];
                        	tmpc2[2] = '\0';
			}
                        tmpc3[0] = line[6];
			if(line[7] == ' ')
                        	tmpc3[1] = '\0';
			else {
                        	tmpc3[1] = line[7];
                        	tmpc3[2] = '\0';
			}
                        tmpc4[0] = line[9];
			if(line[10] == ' ')
                        	tmpc4[1] = '\0';
			else {
                        	tmpc4[1] = line[10];
                        	tmpc4[2] = '\0';
			}
                        if(line[0] == 'X') tmpnum++;
                        if(line[3] == 'X') tmpnum++;
                        if(line[6] == 'X') tmpnum++;
                        if(line[9] == 'X') tmpnum++;
/*	we should not change the positions of 'X' */
			if(tmpnum == 0) {
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
			}
			if(tmpnum == 1 || (tmpnum == 2 && line[6] == 'X')) {
				if(line[0] == 'X') 
                        		if(strcmp(tmpc2, tmpc4) > 0)  {
                                		strcpy(tmpc, tmpc4);
                                		strcpy(tmpc4, tmpc2);
                                		strcpy(tmpc2, tmpc);
                        		}
				if(line[3] == 'X') 
                        		if(strcmp(tmpc1, tmpc4) > 0)  {
                                		strcpy(tmpc, tmpc4);
                                		strcpy(tmpc4, tmpc1);
                                		strcpy(tmpc1, tmpc);
                        		}
				if(line[9] == 'X') 
                        		if(strcmp(tmpc1, tmpc2) > 0)  {
                                		strcpy(tmpc, tmpc1);
                                		strcpy(tmpc1, tmpc2);
                                		strcpy(tmpc2, tmpc);
                        		}
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
                        improperparm[improperparmnum].numX = tmpnum; 
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
if(debug == 1) {
        printf("\nMASS\n");
        for(i=0;i<atomtypenum;i++)
                printf("\n%-2s %9.4lf %9.4lf", atomtype[i].name, atomtype[i].mass,atomtype[i].pol);

        printf("\n\nBOND\n");
        for(i=0;i<bondparmnum;i++)
                printf("\n%-2s %-2s %9.4lf %9.4lf", bondparm[i].name1, bondparm[i].name2,
                bondparm[i].force, bondparm[i].length);

        printf("\n\nANGLE\n");
        for(i=0;i<angleparmnum;i++)
                printf("\n%-2s %-2s %-2s %9.4lf %9.4lf", angleparm[i].name1, angleparm[i].name2,
                angleparm[i].name3, angleparm[i].force, angleparm[i].angle);

        printf("\n\nTORSION\n");
        for(i=0;i<torsionparmnum;i++)
                printf("\n%-2s %-2s %-2s %-2s %9.4lf %9.4lf %9.4lf", torsionparm[i].name1,
                torsionparm[i].name2, torsionparm[i].name3,torsionparm[i].name4,
                torsionparm[i].phase, torsionparm[i].force, torsionparm[i].fterm);

        printf("\n\nIMPROPER\n");
        for(i=0;i<improperparmnum;i++)
                printf("\n%-2s %-2s %-2s %-2s %9.4lf %9.4lf %9.4lf", improperparm[i].name1,
                improperparm[i].name2, improperparm[i].name3,improperparm[i].name4,
                improperparm[i].phase, improperparm[i].force, improperparm[i].fterm);
        printf("\n\nVDW\n");
        for(i=0;i<vdwparmnum;i++)
                printf("\n%-2s %9.4lf %9.4lf", vdwparm[i].name, vdwparm[i].radius, vdwparm[i].pot);

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
	int pos_tor = 0;
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
			pos_tor= ftell(fp);
			continue;
		}
		if (tindex == 1 && spaceline(line) == 1) {
			tindex = 2;
			fseek(fp,pos_tor,0);
			continue;
		}
		if (tindex == 2 && spaceline(line) == 1) {
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
                        if(line[1]== ' ')
                                bondparm[bondparmnum].name1[1] = '\0';
                        else {
                                bondparm[bondparmnum].name1[1] = line[1];
                                bondparm[bondparmnum].name1[2] = '\0';
                        }
			bondparm[bondparmnum].name2[0] = line[3];
                        if(line[4]== ' ')
                                bondparm[bondparmnum].name2[1] = '\0';
                        else {
                                bondparm[bondparmnum].name2[1] = line[4];
                                bondparm[bondparmnum].name2[2] = '\0';
                        }
			sscanf(&line[5], "%lf%lf", &bondparm[bondparmnum].force,
				   &bondparm[bondparmnum].length);
			bondparmnum++;
			if (bondparmnum >= maxbondparm) {
				maxbondparm += MAX_FF_BOND;
				bondparm =
					(BOND_FF *) realloc(bondparm,
										sizeof(BOND_FF) * maxbondparm);
				if (bondparm == NULL) {
					fprintf(stdout,
							"memory allocation error for *bondparm\n");
					exit(1);
				}
			}
		}
		if (aindex == 1) {
			angleparm[angleparmnum].name1[0] = line[0];
                        if(line[1]== ' ')
                                angleparm[angleparmnum].name1[1] = '\0';
                        else {
                                angleparm[angleparmnum].name1[1] = line[1];
                                angleparm[angleparmnum].name1[2] = '\0';
                        }
			angleparm[angleparmnum].name2[0] = line[3];
                        if(line[4]== ' ')
                                angleparm[angleparmnum].name2[1] = '\0';
                        else {
                                angleparm[angleparmnum].name2[1] = line[4];
                                angleparm[angleparmnum].name2[2] = '\0';
                        }
			angleparm[angleparmnum].name3[0] = line[6];
                        if(line[7]== ' ')
                                angleparm[angleparmnum].name3[1] = '\0';
                        else {
                                angleparm[angleparmnum].name3[1] = line[7];
                                angleparm[angleparmnum].name3[2] = '\0';
                        }
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
		if (tindex == 1) { /*we first only read special torsional angle parameters*/
			if(line[0] == 'X' || line[3] == 'X' || line[6] == 'X' || line[9] == 'X')
				continue;
			torsionparm[torsionparmnum].name1[0] = line[0];
                        if(line[1]== ' ')
                                torsionparm[torsionparmnum].name1[1] = '\0';
                        else {
                                torsionparm[torsionparmnum].name1[1] = line[1];
                                torsionparm[torsionparmnum].name1[2] = '\0';
                        }
			torsionparm[torsionparmnum].name2[0] = line[3];
                        if(line[4]== ' ')
                                torsionparm[torsionparmnum].name2[1] = '\0';
                        else {
                                torsionparm[torsionparmnum].name2[1] = line[4];
                                torsionparm[torsionparmnum].name2[2] = '\0';
                        }
			torsionparm[torsionparmnum].name3[0] = line[6];
                        if(line[7]== ' ')
                                torsionparm[torsionparmnum].name3[1] = '\0';
                        else {
                                torsionparm[torsionparmnum].name3[1] = line[7];
                                torsionparm[torsionparmnum].name3[2] = '\0';
                        }
			torsionparm[torsionparmnum].name4[0] = line[9];
                        if(line[10]== ' ')
                                torsionparm[torsionparmnum].name4[1] = '\0';
                        else {
                                torsionparm[torsionparmnum].name4[1] = line[10];
                                torsionparm[torsionparmnum].name4[2] = '\0';
                        }
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
                if (tindex == 2) {
			if(line[0] != 'X' && line[3] != 'X' && line[6] != 'X' && line[9] != 'X')
				continue;
                        torsionparm[torsionparmnum].name1[0] = line[0];
                        if(line[1]== ' ')
                                torsionparm[torsionparmnum].name1[1] = '\0';
                        else {
                                torsionparm[torsionparmnum].name1[1] = line[1];
                                torsionparm[torsionparmnum].name1[2] = '\0';
                        }
                        torsionparm[torsionparmnum].name2[0] = line[3];
                        if(line[4]== ' ')
                                torsionparm[torsionparmnum].name2[1] = '\0';
                        else {
                                torsionparm[torsionparmnum].name2[1] = line[4];
                                torsionparm[torsionparmnum].name2[2] = '\0';
                        }
                        torsionparm[torsionparmnum].name3[0] = line[6];
                        if(line[7]== ' ') 
                                torsionparm[torsionparmnum].name3[1] = '\0';
                        else {
                                torsionparm[torsionparmnum].name3[1] = line[7];
                                torsionparm[torsionparmnum].name3[2] = '\0';
                        }
                        torsionparm[torsionparmnum].name4[0] = line[9];
                        if(line[10]== ' ') 
                                torsionparm[torsionparmnum].name4[1] = '\0';
                        else {
                                torsionparm[torsionparmnum].name4[1] = line[10];
                                torsionparm[torsionparmnum].name4[2] = '\0';
                        }
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
			tmpnum = 0;
			tmpc1[0] = line[0];
                        if(line[1] == ' ') 
                                tmpc1[1] = '\0';
                        else {
                                tmpc1[1] = line[1];
                                tmpc1[2] = '\0';
                        }
			tmpc2[0] = line[3];
                        if(line[4] == ' ') 
                                tmpc2[1] = '\0';
                        else {
                                tmpc2[1] = line[4];
                                tmpc2[2] = '\0';
                        }
			tmpc3[0] = line[6];
                        if(line[7] == ' ') 
                                tmpc3[1] = '\0';
                        else {
                                tmpc3[1] = line[7];
                                tmpc3[2] = '\0';
                        }
			tmpc4[0] = line[9];
                        if(line[10] == ' ')
                                tmpc4[1] = '\0';
                        else {
                                tmpc4[1] = line[10];
                                tmpc4[2] = '\0';
                        }
			if(line[0] == 'X') tmpnum++;
			if(line[3] == 'X') tmpnum++;
			if(line[6] == 'X') tmpnum++;
			if(line[9] == 'X') tmpnum++;

/*      we should not change the position of 'X'*/
                        if(tmpnum == 0) {
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
                        }
                        if(tmpnum == 1 || (tmpnum == 2 && line[6] == 'X')) {
                                if(line[0] == 'X')
                                        if(strcmp(tmpc2, tmpc4) > 0)  {
                                                strcpy(tmpc, tmpc4);
                                                strcpy(tmpc4, tmpc2);
                                                strcpy(tmpc2, tmpc);
                                        }
                                if(line[3] == 'X') 
                                        if(strcmp(tmpc1, tmpc4) > 0)  {
                                                strcpy(tmpc, tmpc4);
                                                strcpy(tmpc4, tmpc1);
                                                strcpy(tmpc1, tmpc);
                                        }
                                if(line[9] == 'X')
                                        if(strcmp(tmpc1, tmpc2) > 0)  {
                                                strcpy(tmpc, tmpc1);
                                                strcpy(tmpc1, tmpc2);
                                                strcpy(tmpc2, tmpc);
                                        }
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
			improperparm[improperparmnum].numX = tmpnum; 
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
	if(debug == 1) {
                printf("\nMASS\n");
                for(i=0;i<atomtypenum;i++)
                        printf("\n%-2s %9.4lf %9.4lf", atomtype[i].name, atomtype[i].mass,atomtype[i].pol);

                printf("\n\nBOND\n");
                for(i=0;i<bondparmnum;i++)
                        printf("\n%-2s %-2s %9.4lf %9.4lf", bondparm[i].name1, bondparm[i].name2,
                        bondparm[i].force, bondparm[i].length);

                printf("\n\nANGLE\n");
                for(i=0;i<angleparmnum;i++)
                        printf("\n%-2s %-2s %-2s %9.4lf %9.4lf", angleparm[i].name1, angleparm[i].name2,
                        angleparm[i].name3, angleparm[i].force, angleparm[i].angle);

                printf("\n\nTORSION\n");
                for(i=0;i<torsionparmnum;i++)
                        printf("\n%-2s %-2s %-2s %-2s %9.4lf %9.4lf %9.4lf", torsionparm[i].name1,
                        torsionparm[i].name2, torsionparm[i].name3,torsionparm[i].name4,
                        torsionparm[i].phase, torsionparm[i].force, torsionparm[i].fterm);

                printf("\n\nIMPROPER\n");
                for(i=0;i<improperparmnum;i++)
                        printf("\n%-2s %-2s %-2s %-2s %9.4lf %9.4lf %9.4lf", improperparm[i].name1,
                        improperparm[i].name2, improperparm[i].name3,improperparm[i].name4,
                        improperparm[i].phase, improperparm[i].force, improperparm[i].fterm);
                printf("\n\nVDW\n");
                for(i=0;i<vdwparmnum;i++)
                        printf("\n%-2s %9.4lf %9.4lf", vdwparm[i].name, vdwparm[i].radius, vdwparm[i].pot);
	}
}

void read_parmchk_parm(char *filename)
{
	FILE *fp;
	char line[MAXCHAR];
	char atomtype[10];
	char equaname[10];
	char corrname[10];
	int  i,j, k;
	int  group;
	int  improper;
	int  ncorr;
	int  nequa;
	int  equtype;
	double bl,blf,ba,baf, cba, cbaf, tor, ctor, itor, ps;
	double mass;
	if ((fp = fopen(filename, "r")) == NULL) {
		fprintf(stdout, "Cannot open file %s in read_parmchk_parm(), exit\n", filename);
		exit(1);
	}

	parmnum = 0;

	for (;;) {
		if (fgets(line, MAXCHAR, fp) == NULL)
			break;
		if (strncmp("PARM", &line[0], 4) == 0) {
			equtype = 0;
			sscanf(&line[4], "%s%d%d%lf%d", atomtype, &improper, &group, &mass, &equtype) ;
			strcpy(parm[parmnum].atomtype, atomtype);
			parm[parmnum].improper = improper;
			parm[parmnum].group    = group;
			parm[parmnum].mass     = mass;
			parm[parmnum].equtype  = equtype;
/*	Initalization */
/*	for the sake of simplicity in the coding, an atom type equals to and corresponds to itself*/
/*      Equal atom types must alos be corresponding atom types */
                        strcpy(parm[parmnum].equa[0].atomtype, atomtype);
                        strcpy(parm[parmnum].corr[0].atomtype, atomtype);
                        parm[parmnum].corr[0].bl   = 0;
                        parm[parmnum].corr[0].blf  = 0;
                        parm[parmnum].corr[0].ba   = 0;
                        parm[parmnum].corr[0].baf  = 0;
                        parm[parmnum].corr[0].cba  = 0;
                        parm[parmnum].corr[0].cbaf = 0;
                        parm[parmnum].corr[0].tor  = 0;
                        parm[parmnum].corr[0].ctor = 0;
                        parm[parmnum].corr[0].improper = 0;
                        parm[parmnum].corr[0].ps    = 0;
                        parm[parmnum].corr[0].type = 0;
                        parm[parmnum].corr[0].pid  = parmnum;
			nequa = 1;
			ncorr = 1;
			parm[parmnum].nequa = 1;
			parm[parmnum].ncorr = 1;
			parmnum++;
			if (parmnum >= maxparmnum) {
				maxparmnum += MAXPARM;
				parm = (PARM *) realloc(parm, sizeof(PARM) * maxparmnum);
				if (parm == NULL) {
					fprintf(stdout, "memory allocation error for *parm\n");
					exit(1);
				}
			}
		}
		if (strncmp("EQUA", &line[0], 4) == 0) {
			sscanf(&line[4], "%s",  equaname); 
			strcpy(parm[parmnum-1].equa[nequa].atomtype, equaname);
			parm[parmnum-1].nequa ++;
			nequa ++;
			if(parm[parmnum-1].nequa > MAXEQUA) {
				fprintf(stdout, "Too many equal atom types for Atom Type %d (%s)\n", parmnum, parm[parmnum-1].atomtype);
				exit(1);
			}
			strcpy(parm[parmnum-1].corr[ncorr].atomtype, equaname);
			parm[parmnum-1].corr[ncorr].bl  =0;
			parm[parmnum-1].corr[ncorr].blf =0;
			parm[parmnum-1].corr[ncorr].ba  =0;
			parm[parmnum-1].corr[ncorr].baf =0;
			parm[parmnum-1].corr[ncorr].cba =0;
			parm[parmnum-1].corr[ncorr].cbaf=0;
			parm[parmnum-1].corr[ncorr].tor =0;
			parm[parmnum-1].corr[ncorr].ctor=0;
			parm[parmnum-1].corr[ncorr].ps  =0;
			parm[parmnum-1].corr[ncorr].improper=0;
			parm[parmnum-1].corr[ncorr].pid = -1;
			parm[parmnum-1].corr[ncorr].type = 1;
			parm[parmnum-1].ncorr ++;
			ncorr ++;
			if(parm[parmnum-1].ncorr > MAXCORR) {
				fprintf(stdout, "Too many corresponding atom types for Atom Type %d (%s)\n", parmnum, parm[parmnum-1].atomtype);
				exit(1);
			}
		}
		if (strncmp("CORR", &line[0], 4) == 0) {
			bl  = 0;
			blf = 0;
			ba  = 0;
			baf = 0;
			cba = 0;
			cbaf= 0;
			tor = 0;
			ctor= 0;
			itor= 0;
			ps=0;
			sscanf(&line[4], "%s%lf%lf%lf%lf%lf%lf%lf%lf%lf",  
				corrname, &bl, &blf, &cba, &cbaf, &ba, &baf, &ctor, &tor, &ps); 
			strcpy(parm[parmnum-1].corr[ncorr].atomtype, corrname);
			parm[parmnum-1].corr[ncorr].bl  =bl;
			parm[parmnum-1].corr[ncorr].blf =blf;
			parm[parmnum-1].corr[ncorr].ba  =ba;
			parm[parmnum-1].corr[ncorr].baf =baf;
			parm[parmnum-1].corr[ncorr].cba =cba;
			parm[parmnum-1].corr[ncorr].cbaf=cbaf;
			parm[parmnum-1].corr[ncorr].tor =tor;
			parm[parmnum-1].corr[ncorr].ctor=ctor;
			parm[parmnum-1].corr[ncorr].improper=ps; /*assign genearl score to improper parameter*/
			parm[parmnum-1].corr[ncorr].ps=ps;
			parm[parmnum-1].corr[ncorr].pid = -1;
			parm[parmnum-1].corr[ncorr].type = 2;
			parm[parmnum-1].ncorr ++;
			ncorr ++;
			if(parm[parmnum-1].ncorr > MAXCORR) {
				fprintf(stdout, "Too many corresponding atom types for Atom Type %d (%s)\n", parmnum, parm[parmnum-1].atomtype);
				exit(1);
			}
		}
		if (strncmp("WEIGHT_BL", &line[0], 9) == 0)   	sscanf(&line[9], "%lf", &wt.BL);
		if (strncmp("WEIGHT_BLF", &line[0], 10) == 0) 	sscanf(&line[10], "%lf", &wt.BLF);
		if (strncmp("WEIGHT_BA", &line[0], 9) == 0)   	sscanf(&line[9], "%lf", &wt.BA);
		if (strncmp("WEIGHT_BAF", &line[0], 10) == 0) 	sscanf(&line[10], "%lf", &wt.BAF);
		if (strncmp("WEIGHT_X", &line[0], 8) == 0)   	sscanf(&line[8], "%lf", &wt.X);
		if (strncmp("WEIGHT_X3", &line[0], 9) == 0)  	sscanf(&line[9], "%lf", &wt.X3);
		if (strncmp("WEIGHT_BA_CTR",  &line[0], 13) == 0)  sscanf(&line[13], "%lf", &wt.BA_CTR);
		if (strncmp("WEIGHT_TOR_CTR", &line[0], 14) == 0)  sscanf(&line[14], "%lf", &wt.TOR_CTR);
		if (strncmp("WEIGHT_IMPROPER", &line[0], 15) == 0) sscanf(&line[15], "%lf", &wt.IMPROPER);
		if (strncmp("WEIGHT_GROUP", &line[0], 12) == 0)	   sscanf(&line[12], "%lf", &wt.GROUP);
		if (strncmp("WEIGHT_EQUTYPE", &line[0], 14) == 0)  sscanf(&line[14], "%lf", &wt.EQUTYPE);
		if (strncmp("THRESHOLD_BA", &line[0], 12) == 0)    sscanf(&line[12], "%lf", &THRESHOLD_BA);
		
		if (strncmp("DEFAULT_BL", &line[0], 10) == 0)   sscanf(&line[10], "%lf", &dv.BL);
		if (strncmp("DEFAULT_BLF", &line[0], 11) == 0) 	sscanf(&line[11], "%lf", &dv.BLF);
		if (strncmp("DEFAULT_BA", &line[0], 10) == 0)   sscanf(&line[10], "%lf", &dv.BA);
		if (strncmp("DEFAULT_BAF", &line[0], 11) == 0) 	sscanf(&line[11], "%lf", &dv.BAF);
		if (strncmp("DEFAULT_BA_CTR", &line[0], 14) == 0)  	sscanf(&line[14], "%lf", &dv.BA_CTR);
		if (strncmp("DEFAULT_BAF_CTR", &line[0], 15) == 0)  	sscanf(&line[15], "%lf", &dv.BAF_CTR);
		if (strncmp("DEFAULT_TOR", &line[0], 11) == 0) 		sscanf(&line[11], "%lf", &dv.TOR);
		if (strncmp("DEFAULT_TOR_CTR", &line[0], 15) == 0) 	sscanf(&line[15], "%lf", &dv.TOR_CTR);
		if (strncmp("DEFAULT_FRACT1", &line[0], 14) == 0) 	sscanf(&line[14], "%lf", &dv.FRACT1);
		if (strncmp("DEFAULT_FRACT2", &line[0], 14) == 0) 	sscanf(&line[14], "%lf", &dv.FRACT2);
	}
	fclose(fp);
	for(i=0;i<parmnum;i++) 
		for(j=0;j<parm[i].ncorr;j++) {
			if(parm[i].corr[j].bl   < 0) parm[i].corr[j].bl  = dv.BL; 
			if(parm[i].corr[j].blf  < 0) parm[i].corr[j].blf = dv.BLF; 
			if(parm[i].corr[j].ba   < 0) parm[i].corr[j].ba  = dv.BA ;
			if(parm[i].corr[j].baf  < 0) parm[i].corr[j].baf = dv.BAF;
			if(parm[i].corr[j].cba  < 0) parm[i].corr[j].cba = dv.BA_CTR;
			if(parm[i].corr[j].cbaf < 0) parm[i].corr[j].cbaf= dv.BAF_CTR;
			if(parm[i].corr[j].tor  < 0) parm[i].corr[j].tor = dv.TOR;
			if(parm[i].corr[j].ctor < 0) parm[i].corr[j].ctor= dv.TOR_CTR;
			parm[i].corr[j].ctor = parm[i].corr[j].ctor * dv.FRACT1  + parm[i].corr[j].ps * dv.FRACT2; 
		}
/*	now finding the pid of equal/corresponding atom types */
	for(i=0;i<parmnum;i++) {
		for(j=0;j<parm[i].nequa;j++) 
			for(k=0;k<parmnum;k++) 
				if(strcmp(parm[i].equa[j].atomtype, parm[k].atomtype) == 0) {
					parm[i].equa[j].pid = k;
					break;
				}
		for(j=0;j<parm[i].ncorr;j++) 
			for(k=0;k<parmnum;k++) 
				if(strcmp(parm[i].corr[j].atomtype, parm[k].atomtype) == 0) {
					parm[i].corr[j].pid = k;
					break;
				}
	}
	if(debug == 1) {
		for(i=0;i<parmnum;i++) {
			printf("PARM %5d %5s %5d %5d %9.4lf %5d\n", i+1, parm[i].atomtype, parm[i].improper, parm[i].group, parm[i].mass, parm[i].equtype);
			for(j=0;j<parm[i].nequa;j++)
				printf("EQUA %5d %5s %5d \n", 
				j+1, parm[i].equa[j].atomtype, parm[i].equa[j].pid + 1); 
			for(j=0;j<parm[i].ncorr;j++)
				printf("CORR %5d %5s %5d %5.1lf %5.1lf %5.1lf %5.1lf %5.1lf %5.1lf %5.1lf %5.1lf %5.1lf %5.1lf %5d\n", 
				j+1, parm[i].corr[j].atomtype, parm[i].corr[j].pid + 1, 
                                parm[i].corr[j].bl, parm[i].corr[j].blf, parm[i].corr[j].cba, parm[i].corr[j].cbaf,
				parm[i].corr[j].ba, parm[i].corr[j].baf, parm[i].corr[j].ctor, parm[i].corr[j].tor, 
				parm[i].corr[j].improper, parm[i].corr[j].ps, parm[i].corr[j].type); 
			printf("\n");
		}
		printf("WEIGHTS\n");
		printf("WT for BL: %9.2lf\n", wt.BL);
		printf("WT for BLF: %9.2lf\n", wt.BLF);
		printf("WT for BA: %9.2lf\n", wt.BA);
		printf("WT for BAF: %9.2lf\n", wt.BAF);
		printf("WT for X: %9.2lf\n", wt.X);
		printf("WT for X3: %9.2lf\n", wt.X3);
		printf("WT for BA_CTR: %9.2lf\n", wt.BA_CTR);
		printf("WT for TOR_CTR: %9.2lf\n", wt.TOR_CTR);
		printf("WT for IMPROPER: %9.2lf\n", wt.IMPROPER);
		printf("WT for GROUP: %9.2lf\n", wt.GROUP);
		printf("WT for EQUTYPE: %9.2lf\n", wt.EQUTYPE);
		printf("\nDEFAULT VALUES\n");
		printf("DEFAULT for BL: %9.2lf\n", dv.BL);
		printf("DEFAULT for BLF: %9.2lf\n", dv.BLF);
		printf("DEFAULT for BA: %9.2lf\n", dv.BA);
		printf("DEFAULT for BAF: %9.2lf\n", dv.BAF);
		printf("DEFAULT for BA_CTR: %9.2lf\n", dv.BA_CTR);
		printf("DEFAULT for TOR: %9.2lf\n", dv.TOR);
		printf("DEFAULT for TOR_CTR: %9.2lf\n", dv.TOR_CTR);
		printf("DEFAULT for FRACT1: %9.2lf\n", dv.FRACT1);
		printf("DEFAULT for FRACT2: %9.2lf\n", dv.FRACT2);
		printf("THRESHOLD for BA: %9.2lf\n", THRESHOLD_BA);
	}
}

void assign_parmid(void)
{
	int i, j;
	for (i = 0; i < atomnum; i++) {
		parmid[i] = -1;
		for (j = 0; j < parmnum; j++)
			if (strcmp(parm[j].atomtype, atom[i].ambername) == 0) {
				parmid[i] = j;
				if(parm[j].improper == 1) 
					atom[i].improper = 1;
				break;
			}
	}
}

int empangle(char *tmpc1, char *tmpc2, char *tmpc3, char *name1,
			 char *name2, char *name3, int id1, int id2, int id3, int flag)
{
	int num1 = -1, num2 = -1;
	double bondlength1 = 0.0;
	double bondlength2 = 0.0;
	double cparm, dparm, zparm1, zparm2;
	double angle;
	double force;

	for (i = 0; i < angleparmnum; i++)
		if (strcmp(angleparm[i].name1, tmpc1) == 0 && 
		    strcmp(angleparm[i].name2, tmpc2) == 0 && 
		    strcmp(angleparm[i].name3, tmpc1) == 0) {
			num1 = i;
			break;
		}
	if (num1 == -1)
		return 0;
	for (i = 0; i < angleparmnum; i++)
		if (strcmp(angleparm[i].name1, tmpc3) == 0 && 
		    strcmp(angleparm[i].name2, tmpc2) == 0 && 
		    strcmp(angleparm[i].name3, tmpc3) == 0) {
			num2 = i;
			break;
		}
	if (num2 == -1)
		return 0;

	angle = 0.5 * (angleparm[num1].angle + angleparm[num2].angle);
	for (i = 0; i < bondparmnum; i++)
		if ((strcmp(bondparm[i].name1, tmpc1) == 0 &&
                     strcmp(bondparm[i].name2, tmpc2) == 0) ||
                    (strcmp(bondparm[i].name2, tmpc1) == 0 && 
                     strcmp(bondparm[i].name1, tmpc2) == 0)) {
			bondlength1 = bondparm[i].length;
			break;
		}
	if (bondlength1 == 0.0)
		return 0;

	for (i = 0; i < bondparmnum; i++)
		if ((strcmp(bondparm[i].name1, tmpc2) == 0 &&
                     strcmp(bondparm[i].name2, tmpc3) == 0) ||
                    (strcmp(bondparm[i].name2, tmpc2) == 0 && 
                     strcmp(bondparm[i].name1, tmpc3) == 0)) {
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
		dparm / ((bondlength1 + bondlength2) * (bondlength1 + bondlength2));
	force =
		143.9 * zparm1 * cparm * zparm2 * exp(-2 * dparm) / (bondlength1 +
															 bondlength2);
	force /= sqrt(angle * 3.1415926 / 180.0);
	if(flag == 0) {
		fprintf(fpout, "%-2s-%-2s-%-2s%9.3lf%12.3lf", name1, name2, name3, force, angle);
		fprintf(fpout, "   Calculated with empirical approach for %s-%s-%s\n", tmpc1, tmpc2, tmpc3);
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
	}
	if(flag == 1) {
		bestba.angle = angle;
		bestba.force = force;
		strcpy(bestba.name1, tmpc1);
		strcpy(bestba.name2, tmpc2);
		strcpy(bestba.name3, tmpc3);
		bestbaid = 999999;
	}
	return 1;
}

double equtype_penalty(int id1, int id2, int id3, int id4){
int tn1; 
int tn2;
	if(parm[id1].equtype == 0 && parm[id2].equtype == 0) 
		return 0.0;

	tn1 = parm[id1].equtype + parm[id2].equtype;
	tn2 = parm[id3].equtype + parm[id4].equtype;
	if(tn2 == 3) {
		if(tn1 == 3) 
			return 0;
		else
			return wt.EQUTYPE;
	}	
	if(tn1 == 3 && tn2 != 3) return wt.EQUTYPE;
	return 0.0;
}

void chk_atomtype(void)
{
/*	in this version, only polarizability parameter is replaced, as mass must be present in the PARMCHK.PARM file */
	int i, j;
	int suc;
	int pid;
	int pid2;
	char tmpc[5];
	fprintf(fpout, "%s\n", "Remark line goes here");
	fprintf(fpout, "%s\n", "MASS");
	for (i = 0; i < atomnum; i++) {
		suc = 0;
		pid=parmid[i];
		for (j = 0; j < atomtypenum; j++)
			if (strcmp(atom[i].ambername, atomtype[j].name) == 0) {
				suc = 1;
				if(allparm_flag == 1) {
					fprintf(fpout, "%-2s %-8.3lf   %8.3lf\n",
						atom[i].ambername, atomtype[j].mass,
						atomtype[j].pol);
				}
				break;
			}
/*	check equal atom types	*/	
/*	equal atom types	*/
		if (suc == 0 && parm[pid].nequa > 1) 
                        for(j=1;j<parm[pid].nequa;j++) {
                                pid2 = parm[pid].equa[j].pid;
                                strcpy(tmpc, parm[pid2].atomtype);
                                for (k = 0; k < atomtypenum2; k++)
                                        if (strcmp(tmpc, atomtype[k].name) == 0) {
                                                suc = 1;
                                                fprintf(fpout, "%-2s %-8.3lf   %8.3lf",
                                                                atom[i].ambername, parm[pid].mass,
                                                                atomtype[k].pol);
                                                fprintf(fpout, "               %s %-3s\n", "same as", atomtype[k].name);
                                                atomtype[atomtypenum] = atomtype[k];
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
                                if(suc == 1) break;
                        }
/*	corresponding atom types	*/
		if (suc == 0 && parm[pid].ncorr > 1) 
			for(j=1;j<parm[pid].ncorr;j++) {
				if(parm[pid].corr[j].type <= 1) continue;
				pid2 = parm[pid].corr[j].pid;
				strcpy(tmpc, parm[pid2].atomtype);
				for (k = 0; k < atomtypenum2; k++)
					if (strcmp(tmpc, atomtype[k].name) == 0) {
						suc = 1;
						fprintf(fpout, "%-2s %-8.3lf   %8.3lf",
								atom[i].ambername, parm[pid].mass,
								atomtype[k].pol);
						fprintf(fpout, "               %s %-3s\n", "same as", atomtype[k].name);
						atomtype[atomtypenum] = atomtype[k];
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
				if(suc == 1) break;
			}
		if (suc == 0) {
			fprintf(fpout, "%-2s %-8.3lf   %8.3lf", atom[i].ambername, parm[pid].mass, 0.0);
			fprintf(fpout, "               %5s\n", "ATTN, no polarizability parameter");
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

int vdw(char *at_name, char *corr_name, int nparm)
{
	int suc = 0;
	int i;
	for (i = 0; i < nparm; i++)
		if (strcmp(vdwparm[i].name, corr_name) == 0) {
			suc = 1;
			fprintf(fpout, "  %-2s%16.4lf%8.4lf",
					at_name, vdwparm[i].radius, vdwparm[i].pot);
			fprintf(fpout, "             %s %-3s\n", "same as",
					vdwparm[i].name);
			vdwparm[vdwparmnum] = vdwparm[i];
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

void chk_vdw(void)
{
	int i, j;
	int pid, pid2;
	int suc;
	fprintf(fpout, "\n%s\n", "NONBON");
	for (i = 0; i < atomnum; i++) {
		suc = 0;
		pid = parmid[i];
		for (j = 0; j < vdwparmnum; j++)
			if (strcmp(vdwparm[j].name, atom[i].ambername) == 0) {
				suc = 1;
				if(allparm_flag == 1) {
					fprintf(fpout, "  %-2s%16.4lf%8.4lf\n",
						vdwparm[j].name, vdwparm[j].radius, vdwparm[j].pot);
				}
				break;
			}
		if (suc == 0 && parm[pid].nequa > 1) 
			for(j=1; j< parm[pid].nequa; j++) {
				pid2 = parm[pid].equa[j].pid;	
				suc = vdw(atom[i].ambername, parm[pid2].atomtype, vdwparmnum2);
				if(suc == 1) break;
			}
		if (suc == 0 && parm[pid].ncorr > 1) 
			for(j=1; j< parm[pid].ncorr; j++) {
				if(parm[pid].corr[j].type <= 1) continue;
				pid2 = parm[pid].corr[j].pid;	
				suc = vdw(atom[i].ambername, parm[pid2].atomtype, vdwparmnum2);
				if(suc == 1) break;
			}
		if (suc == 0) {
			fprintf(fpout, "  %-2s%16.4lf%8.4lf", atom[i].ambername, 0.0,
					0.0);
			fprintf(fpout, "             %s\n", "ATTN, need revision");
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
		 char *corr_name2, int nparm, int index)
{
	int suc = 0;
	int k;
	for (k = 0; k < nparm; k++)
		if ((strcmp(bondparm[k].name1, corr_name1) == 0 &&
		     strcmp(bondparm[k].name2, corr_name2) == 0) ||
                    (strcmp(bondparm[k].name1, corr_name2) == 0 &&
		     strcmp(bondparm[k].name2, corr_name1) == 0)) {
			suc = 1;
			if (index == 0 && allparm_flag == 1) {
				fprintf(fpout, "%-2s-%-2s%8.2lf%8.3lf\n", 
						at_name1, at_name2, 
						bondparm[k].force, bondparm[k].length);
			}
			if (index == 1) {
				bestblid = k;
				break;
			}
		}
	return suc;
}

void chk_bond(void)
{
	int i, j, m,n;
	int suc, suc2;
	int pid1, pid2;
	int pid3, pid4;
	char tmpc[5];
	char tmpc1[5], tmpc2[5];
	char name1[5], name2[5];
	double score, score1, score2;
	fprintf(fpout, "\n%s\n", "BOND");

	for (i = 0; i < atomnum; i++)
		for (j = i + 1; j < atomnum; j++)
			if (atom[i].con[0] == j || atom[i].con[1] == j
				|| atom[i].con[2] == j || atom[i].con[3] == j
				|| atom[i].con[4] == j || atom[i].con[5] == j) {
				suc = 0;
				strcpy(tmpc1, atom[i].ambername);
				strcpy(tmpc2, atom[j].ambername);
                                if(strcmp(tmpc1, tmpc2) > 0)  {
                                        strcpy(tmpc, tmpc2);
                                        strcpy(tmpc2, tmpc1);
                                        strcpy(tmpc1, tmpc);
                                }
				strcpy(name1, tmpc1);
				strcpy(name2, tmpc2);

				suc = bond(name1, name2, tmpc1, tmpc2, bondparmnum, 0);
				if(suc == 1) continue;
/*	for equal atom types	*/
                                if(suc == 0) {
					bestblid = -1;
                                        pid1 = parmid[i];
                                        pid2 = parmid[j];
                                        for(m=0;m<parm[pid1].nequa; m++) {
						if(suc == 1) break;
                                                pid3 = parm[pid1].equa[m].pid;
                                                strcpy(tmpc1, parm[pid3].atomtype);
                                                for(n=0;n<parm[pid2].nequa; n++) {
							if(suc == 1) break;
                                                        if(m==0 && n== 0) continue;  
                                                        pid4 = parm[pid2].equa[n].pid;
                                                        strcpy(tmpc2, parm[pid4].atomtype);
                                                        suc2 = bond(name1, name2, tmpc1, tmpc2, bondparmnum2, 1);
							if(suc2 == 1) {
								bestscore = 0;
								suc = 1;
							}
                                        	}
					}
				}
/*	for corresponding atom types	*/
				if(suc == 0) {
					pid1 = parmid[i];
					pid2 = parmid[j];
					bestblid = -1;
					bestscore = INITSCORE;	
					for(m=0;m<parm[pid1].ncorr; m++) {
						pid3 = parm[pid1].corr[m].pid;
						strcpy(tmpc1, parm[pid3].atomtype);	
						score1= parm[pid1].corr[m].bl * wt.BL +
						        parm[pid1].corr[m].blf* wt.BLF;
						for(n=0;n<parm[pid2].ncorr; n++) {
							if(parm[pid1].corr[m].type <= 1 && 
                                                           parm[pid2].corr[n].type <= 1) 
								continue;	
							pid4 = parm[pid2].corr[n].pid;
							strcpy(tmpc2, parm[pid4].atomtype);	
							score2= parm[pid2].corr[n].bl * wt.BL +
						        	parm[pid2].corr[n].blf* wt.BLF;
							score = score1 + score2;
							if(parm[pid3].group != parm[pid4].group) score += wt.GROUP;
							equtype_penalty_score = equtype_penalty(pid1, pid2, pid3, pid4);	
							score += equtype_penalty_score;
							if(score < bestscore) {
								suc2= bond(name1, name2, tmpc1, tmpc2, bondparmnum2, 1);
								if(suc2== 1) {
									bestscore = score;
									suc = 1;
								}
							}
						}
					}
				}
				if(suc==1 && bestblid >= 0) {
                               		fprintf(fpout, "%-2s-%-2s%8.2lf%8.3lf", name1,
                                               		name2, bondparm[bestblid].force, bondparm[bestblid].length);
                               		fprintf(fpout, "       same as %2s-%2s, penalty score=%5.1lf\n",
							bondparm[bestblid].name1, bondparm[bestblid].name2, bestscore);
                               		bondparm[bondparmnum] = bondparm[bestblid]; 
                               		strcpy(bondparm[bondparmnum].name1, name1);
                               		strcpy(bondparm[bondparmnum].name2, name2);
                               		bondparmnum++;
                               		if (bondparmnum >= maxbondparm) {
                                       		maxbondparm += MAX_FF_BOND;
                                       		bondparm =
                                               		(BOND_FF *) realloc(bondparm,
                                                                                       sizeof(BOND_FF) * maxbondparm);
                                       		if (bondparm == NULL) {
                                               		fprintf(stdout,
                                                      	         	"memory allocation error for *bondparm\n");
                                               		exit(1);
                                       		}
                               		}
					bestbaid = -1;
					continue;
                        	}
				else {
					fprintf(fpout, "%-2s-%-2s%8.2lf%8.3lf", name1, name2, 0.0, 0.0);
					fprintf(fpout, "       %s\n", "ATTN, need revision");
					strcpy(bondparm[bondparmnum].name1, name1);
					strcpy(bondparm[bondparmnum].name2, name2);
					bondparmnum++;
					if (bondparmnum >= maxbondparm) {
						maxbondparm += MAX_FF_BOND;
						bondparm =
							(BOND_FF *) realloc(bondparm,
												sizeof(BOND_FF) *
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
		  char *corr_name2, char *corr_name3, int nparm, int index)
{
	int suc = 0;
	int l=0;
	for (l = 0; l < nparm; l++)
		if ((strcmp(angleparm[l].name1, corr_name1) == 0 && 
                     strcmp(angleparm[l].name2, corr_name2) == 0 &&
                     strcmp(angleparm[l].name3, corr_name3) == 0) ||
		    (strcmp(angleparm[l].name3, corr_name1) == 0 && 
                     strcmp(angleparm[l].name2, corr_name2) == 0 &&
                     strcmp(angleparm[l].name1, corr_name3) == 0)) {
			suc = 1;
			if (allparm_flag == 1 && index == 0) {
				fprintf(fpout,
						"%-2s-%-2s-%-2s%9.3lf%12.3lf\n",
						at_name1, at_name2, at_name3,
						angleparm[l].force, angleparm[l].angle);
			}
			if (index == 1) {
				bestba = angleparm[l];
				bestbaid = l;
				break;
			}
		}
	return suc;
}



void chk_angle(void)
{

	int i, j, k, m, n, o;
	int suc, suc2;
	int pid1, pid2, pid3;
	int pid4, pid5, pid6;
	char tmpc[5];
	char tmpc1[5], tmpc2[5], tmpc3[5];
	char tmpc4[5], tmpc5[5], tmpc6[5];
	char name1[5], name2[5], name3[5];
	double score, score1, score2, score3;
	fprintf(fpout, "\n%s\n", "ANGLE");

	/* NB: non-standard indentation in next four lines; for readability  */
	for (i = 0; i < atomnum; i++) {
		for (j = 0; j < atomnum; j++) {
			for (k = 0; k < atomnum; k++) {
				if (i != k) {
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
                                			if(strcmp(tmpc1, tmpc3) > 0)  {
                                        			strcpy(tmpc, tmpc3);
                                        			strcpy(tmpc3, tmpc1);
                                        			strcpy(tmpc1, tmpc);
                                			}
							strcpy(name1, tmpc1);
							strcpy(name2, tmpc2);
							strcpy(name3, tmpc3);
							suc = angle(name1, name2, name3, tmpc1, tmpc2, tmpc3, angleparmnum, 0);
							if(suc == 1) continue;

/* for equal atom types */
							if(suc == 0) {
								bestbaid = -1;
                                        			pid1 = parmid[j]; /* it is j not i, as we do not want to replace atom type for central atom */
                                        			pid2 = parmid[i]; /* it is i not j */
                                        			pid3 = parmid[k];
                                                                for(m=0;m<parm[pid1].nequa; m++) {
									if(suc == 1) break;
                                                                        pid4 = parm[pid1].equa[m].pid;
                                                                        strcpy(tmpc4, parm[pid4].atomtype);
                                                                        for(n=0;n<parm[pid2].nequa; n++) {
										if(suc == 1) break;
                                                                                pid5 = parm[pid2].equa[n].pid;
                                                                                strcpy(tmpc5, parm[pid5].atomtype);
                                                                                for(o=0;o<parm[pid3].nequa; o++) {
											if(suc == 1) break;
                                                                                        if(m== 0 && n==0 && o==0) continue;
                                                                                        pid6 = parm[pid3].corr[o].pid;
                                                                                        strcpy(tmpc6, parm[pid6].atomtype);
                                                                                        suc2 = angle(name1, name2, name3, tmpc5, tmpc4, tmpc6, angleparmnum2, 1);
											if(suc2 == 1) {
												bestscore = 0;
												suc = 1;
											}
                                                                                }
                                                                        }
                                                                }
							}
/*	for corresponding atom types */
                                			if(suc == 0) {
                                                                pid1 = parmid[i];
                                                                pid2 = parmid[j];
                                                                pid3 = parmid[k];
								bestbaid = -1;
                                        			bestscore = INITSCORE;
                                        			for(m=0;m<parm[pid1].ncorr; m++) {
                                                			pid4 = parm[pid1].corr[m].pid;
                                                			strcpy(tmpc4, parm[pid4].atomtype);
                                                			score1= parm[pid1].corr[m].ba * wt.BA +
                                                        			parm[pid1].corr[m].baf* wt.BAF;
                                                			for(n=0;n<parm[pid2].ncorr; n++) {
                                                        			pid5 = parm[pid2].corr[n].pid;
                                                        			strcpy(tmpc5, parm[pid5].atomtype);
                                                        			score2= parm[pid2].corr[n].cba * wt.BA +
                                                                			parm[pid2].corr[n].cbaf* wt.BAF;
										score2 *= wt.BA_CTR;
                                                				for(o=0;o<parm[pid3].ncorr; o++) {
											if(parm[pid1].corr[m].type <= 1 && 
                                                           				   parm[pid2].corr[n].type <= 1 &&
                                                           				   parm[pid3].corr[o].type <= 1) 
												continue;	
                                                        				pid6 = parm[pid3].corr[o].pid;
                                                        				strcpy(tmpc6, parm[pid6].atomtype);
                                                        				score3= parm[pid3].corr[o].ba * wt.BA +
                                                                				parm[pid3].corr[o].baf* wt.BAF;
											score = score1 + score2 + score3 + wt.GROUP;
											if(parm[pid4].group == parm[pid5].group && 
											   parm[pid4].group == parm[pid6].group)
                                                                                        	score -= wt.GROUP;
                                                        				if(score < bestscore && score <=THRESHOLD_BA) {
                                                                				suc2= angle(name1, name2, name3, tmpc4, tmpc5, tmpc6, angleparmnum2, 1);
                                                                				if(suc2== 1) {
													bestscore = score;
													suc = 1;
												}
                                                        				}
										}
                                                			}
								}
							}
							if(suc == 1 && bestbaid >= 0) {
                                				fprintf(fpout, "%-2s-%-2s-%-2s%9.3lf%12.3lf",
                                                				name1, name2, name3,
										bestba.force, bestba.angle);
                                				fprintf(fpout, "   same as %-2s-%-2s-%-2s, penalty score=%5.1lf\n",
                                                				bestba.name1, bestba.name2, bestba.name3, bestscore); 
                                				angleparm[angleparmnum] = bestba;
                                				strcpy(angleparm[angleparmnum].name1, name1);
                                				strcpy(angleparm[angleparmnum].name2, name2);
                                				strcpy(angleparm[angleparmnum].name3, name3);
                                				angleparmnum++;
                                				if (angleparmnum >= maxangleparm) {
                                        				maxangleparm += MAX_FF_ANGLE;
                                        				angleparm = (ANGLE *) realloc(angleparm,
                                                                                 	sizeof(ANGLE) * maxangleparm);
                                        				if (angleparm == NULL) {
                                                				fprintf(stdout,
                                                               				"memory allocation error for *angleparm\n");
                                                				exit(1);
                                        				}
								}
								bestbaid = -1;
								continue;
                                        		}
/* from here estimate the bond angle parameters with empirical method for corresponding names */
							if (suc == 0) {
								suc = empangle(name1, name2, name3, name1, name2, name3,
											 atom[i].atomicnum,
											 atom[j].atomicnum,
											 atom[k].atomicnum, 0);
								if(suc == 1) continue;
							}

/* for equal atom types*/
                                                        if (suc == 0) {
								bestbaid = -1;
                                                                pid1 = parmid[j]; /* it is j, not i*/
                                                                pid2 = parmid[i]; /* it is i, not j*/
                                                                pid3 = parmid[k];
                                                                for(m=0;m<parm[pid1].nequa; m++) {
									if(suc == 1) break;
                                                                        pid4 = parm[pid1].equa[m].pid;
                                                                        strcpy(tmpc4, parm[pid4].atomtype);
                                                                        for(n=0;n<parm[pid2].nequa; n++) {
										if(suc == 1) break;
                                                                                pid5 = parm[pid2].equa[n].pid;
                                                                                strcpy(tmpc5, parm[pid5].atomtype);
                                                                                for(o=0;o<parm[pid3].nequa; o++) {
											if(suc == 1) break;
                                                                                        if(m==0 && n== 0 && o==0) continue;
                                                                                        pid6 = parm[pid3].equa[o].pid;
                                                                                        strcpy(tmpc6, parm[pid6].atomtype);
                                                                                        suc2= empangle(tmpc5, tmpc4, tmpc6, name1, name2, name3,
                                                                                                                atom[i].atomicnum,
                                                                                                                atom[j].atomicnum,
                                                                                                                atom[k].atomicnum, 1);
											if(suc2== 1) {
												bestscore = 0;
												suc = 1;
											}
                                                                                }
                                                                        }
                                                                }
                                                        }
/*	for corresponding atom types*/
							if (suc == 0) {
                                                                pid1 = parmid[i];
                                                                pid2 = parmid[j];
                                                                pid3 = parmid[k];
                                                                bestscore = INITSCORE;
								bestbaid  = -1;
                                                                for(m=0;m<parm[pid1].ncorr; m++) {
                                                                        pid4 = parm[pid1].corr[m].pid;
                                                                        strcpy(tmpc4, parm[pid4].atomtype);
                                                                        score1= parm[pid1].corr[m].ba * wt.BA +
                                                                                parm[pid1].corr[m].baf* wt.BAF;
                                                                        for(n=0;n<parm[pid2].ncorr; n++) {
                                                                                pid5 = parm[pid2].corr[n].pid;
                                                                                strcpy(tmpc5, parm[pid5].atomtype);
                                                                                score2= parm[pid2].corr[n].cba * wt.BA +
                                                                                        parm[pid2].corr[n].cbaf* wt.BAF;
                                                                                score2 *= wt.BA_CTR; 
                                                                                for(o=0;o<parm[pid3].ncorr; o++) {
											if(parm[pid1].corr[m].type <= 1 && 
                                                           				   parm[pid2].corr[n].type <= 1 &&
                                                           				   parm[pid3].corr[o].type <= 1) 
												continue;	
                                                                                        pid6 = parm[pid3].corr[o].pid;
                                                                                        strcpy(tmpc6, parm[pid6].atomtype);
                                                                                        score3= parm[pid3].corr[o].ba * wt.BA +
                                                                                                parm[pid3].corr[o].baf* wt.BAF;
                                                                                        score = score1 + score2 + score3 + wt.GROUP;
											if(parm[pid4].group == parm[pid5].group &&
                                                                                           parm[pid4].group == parm[pid6].group)
                                                                                        	score -= wt.GROUP;
                                                                                        if(score < bestscore) {
												suc2= empangle(tmpc4, tmpc5, tmpc6, name1, name2, name3,
											 			atom[i].atomicnum,
											 			atom[j].atomicnum,
											 			atom[k].atomicnum, 1);
                                                                                                if(suc2== 1) {
													bestscore = score;
													suc = 1;
												}
                                                                                        }       
                                                                                }       
                                                                        }       
								}
							}
							if(suc==1 && bestbaid >= 0) {
                                				fprintf(fpout, "%-2s-%-2s-%-2s%9.3lf%12.3lf",
                                                				name1, name2, name3,
										bestba.force, bestba.angle);
                                				fprintf(fpout, "   Calculated using %2s-%2s-%2s, penalty score=%5.1lf\n",
                                                				bestba.name1, bestba.name2, bestba.name3, bestscore); 
                                				angleparm[angleparmnum] = bestba;
                                				strcpy(angleparm[angleparmnum].name1, name1);
                                				strcpy(angleparm[angleparmnum].name2, name2);
                                				strcpy(angleparm[angleparmnum].name3, name3);
                                				angleparmnum++;
                                				if (angleparmnum >= maxangleparm) {
                                        				maxangleparm += MAX_FF_ANGLE;
                                        				angleparm = (ANGLE *) realloc(angleparm,
                                                                                 	sizeof(ANGLE) * maxangleparm);
                                        				if (angleparm == NULL) {
                                                				fprintf(stdout,
                                                               				"memory allocation error for *angleparm\n");
                                                				exit(1);
                                        				}
								}
								bestbaid = -1;
								continue;
							}
							else {
								fprintf(fpout,
										"%-2s-%-2s-%-2s%9.3lf  %10.3lf",
										name1, name2, name3, 0.0, 0.0);
								fprintf(fpout, "   %s\n",
										"ATTN, need revision");
								strcpy(angleparm[angleparmnum].name1,
									   name1);
								strcpy(angleparm[angleparmnum].name2,
									   name2);
								strcpy(angleparm[angleparmnum].name3,
									   name3);
								angleparmnum++;
								if (angleparmnum >= maxangleparm) {
									maxangleparm += MAX_FF_ANGLE;
									angleparm = (ANGLE *) realloc(angleparm, sizeof(ANGLE) * maxangleparm);
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
}

int torsion(char *at_name1, char *at_name2, char *at_name3, char *at_name4,
			char *corr_name1, char *corr_name2, char *corr_name3,
			char *corr_name4, int nparm, int index)
{
	int suc = 0;
	int m, n;
	for (m = 0; m < nparm; m++)
		if ((strcmp(torsionparm[m].name1, corr_name1) == 0 &&
                     strcmp(torsionparm[m].name2, corr_name2) == 0 &&
                     strcmp(torsionparm[m].name3, corr_name3) == 0 &&
                     strcmp(torsionparm[m].name4, corr_name4) == 0) ||
		    (strcmp(torsionparm[m].name4, corr_name1) == 0 &&
                     strcmp(torsionparm[m].name3, corr_name2) == 0 &&
                     strcmp(torsionparm[m].name2, corr_name3) == 0 &&
                     strcmp(torsionparm[m].name1, corr_name4) == 0)) {
			suc = 1;
			n = m;
			if (allparm_flag == 1 && index == 0) {
				while (torsionparm[n].fterm < 0) {
                        		fprintf(fpout,
                                		"%-2s-%-2s-%-2s-%-2s%4d%9.3lf%14.3lf%16.3lf\n",
                                        	at_name1, at_name2, at_name3, at_name4, 
						torsionparm[n].mul, torsionparm[n].force,
                                        	torsionparm[n].phase, torsionparm[n].fterm);
					n++;
				}
				fprintf(fpout,
						"%-2s-%-2s-%-2s-%-2s%4d%9.3lf%14.3lf%16.3lf\n",
						at_name1, at_name2, at_name3, at_name4, 
						torsionparm[n].mul, torsionparm[n].force,
						torsionparm[n].phase, torsionparm[n].fterm);
			}
			if (index == 1) {
				besttorid = m;
				break;
			}
		}
	return suc;
}

void print_torsion(int besttorid, char *name1, char *name2, char *name3, char *name4) {

	while (torsionparm[besttorid].fterm < 0) {
		fprintf(fpout, "%-2s-%-2s-%-2s-%-2s%4d%9.3lf%14.3lf%16.3lf",
                       	name1, name2, name3, name4, 
			torsionparm[besttorid].mul,   torsionparm[besttorid].force,
                       	torsionparm[besttorid].phase, torsionparm[besttorid].fterm);
        	fprintf(fpout, "      same as %-2s-%-2s-%-2s-%-2s\n",
                       	torsionparm[besttorid].name1,
                       	torsionparm[besttorid].name2,
                       	torsionparm[besttorid].name3,
                       	torsionparm[besttorid].name4);
        	torsionparm[torsionparmnum] = torsionparm[besttorid];
        	strcpy(torsionparm[torsionparmnum].name1, name1);
        	strcpy(torsionparm[torsionparmnum].name2, name2);
        	strcpy(torsionparm[torsionparmnum].name3, name3);
        	strcpy(torsionparm[torsionparmnum].name4, name4);
        	torsionparmnum++;
        	if (torsionparmnum >= maxtorsionparm) {
        		maxtorsionparm += MAX_FF_TORSION;
                	torsionparm = (TORSION *) realloc(torsionparm, sizeof(TORSION) * maxtorsionparm);
                	if (torsionparm == NULL) {
                		fprintf(stdout, "memory allocation error for *torsionparm\n");
                        	exit(1);
                	}
        	}
        	besttorid ++;
	}
	fprintf(fpout, "%-2s-%-2s-%-2s-%-2s%4d%9.3lf%14.3lf%16.3lf",
               	name1, name2, name3, name4, 
		torsionparm[besttorid].mul,   torsionparm[besttorid].force,
                torsionparm[besttorid].phase, torsionparm[besttorid].fterm);
	fprintf(fpout, "      same as %-2s-%-2s-%-2s-%-2s, penalty score=%5.1lf\n",
               	torsionparm[besttorid].name1,
               	torsionparm[besttorid].name2,
               	torsionparm[besttorid].name3,
               	torsionparm[besttorid].name4, bestscore);
	torsionparm[torsionparmnum] = torsionparm[besttorid];
	strcpy(torsionparm[torsionparmnum].name1, name1);
	strcpy(torsionparm[torsionparmnum].name2, name2);
	strcpy(torsionparm[torsionparmnum].name3, name3);
	strcpy(torsionparm[torsionparmnum].name4, name4);
	torsionparmnum++;
	if (torsionparmnum >= maxtorsionparm) {
        	maxtorsionparm += MAX_FF_TORSION;
                torsionparm = (TORSION *) realloc(torsionparm, sizeof(TORSION) * maxtorsionparm);
                if (torsionparm == NULL) {
                	fprintf(stdout, "memory allocation error for *torsionparm\n");
                         exit(1);
                }
         }
}

void chk_torsion(void)
{

	int i, j, k, l;
	int m, n, p, q;
	int pid1, pid2, pid3, pid4;
	int pid5, pid6, pid7, pid8;
	int suc, suc2;
	char tmpc[5];
	char tmpc1[5], tmpc2[5], tmpc3[5], tmpc4[5];
	char tmpc5[5], tmpc6[5], tmpc7[5], tmpc8[5];
	char name1[5], name2[5], name3[5], name4[5];
	double score, score1, score2, score3, score4;
	fprintf(fpout, "\n%s\n", "DIHE");

	/* NB: non-standard indentation in next four lines; for readability  */
	for (i = 0; i < atomnum; i++) {
		for (j = 0; j < atomnum; j++) {
			for (k = 0; k < atomnum; k++) {
				for (l = 0; l < atomnum; l++) {
					if (i != k && l != j) {
						if (atom[i].con[0] == j || atom[i].con[1] == j
							|| atom[i].con[2] == j || atom[i].con[3] == j
							|| atom[i].con[4] == j
							|| atom[i].con[5] == j) {
							if (atom[j].con[0] == k || atom[j].con[1] == k
								|| atom[j].con[2] == k
								|| atom[j].con[3] == k
								|| atom[j].con[4] == k
								|| atom[j].con[5] == k) {
								if (atom[l].con[0] == k
									|| atom[l].con[1] == k
									|| atom[l].con[2] == k
									|| atom[l].con[3] == k
									|| atom[l].con[4] == k
									|| atom[l].con[5] == k) {
									suc = 0;
									strcpy(tmpc1, atom[i].ambername);
									strcpy(tmpc2, atom[j].ambername);
									strcpy(tmpc3, atom[k].ambername);
									strcpy(tmpc4, atom[l].ambername);

                                					if(strcmp(tmpc2, tmpc3) > 0)  {
                                        					strcpy(tmpc, tmpc3);
                                        					strcpy(tmpc3, tmpc2);
                                        					strcpy(tmpc2, tmpc);
                                        					strcpy(tmpc, tmpc4);
                                        					strcpy(tmpc4, tmpc1);
                                        					strcpy(tmpc1, tmpc);
                                					} else if (strcmp(tmpc2, tmpc3) == 0) {
										if(strcmp(tmpc1, tmpc4) > 0) {
                                        						strcpy(tmpc, tmpc4);
                                        						strcpy(tmpc4, tmpc1);
                                        						strcpy(tmpc1, tmpc);
										}
									}
									strcpy(name1, tmpc1);
									strcpy(name2, tmpc2);
									strcpy(name3, tmpc3);
									strcpy(name4, tmpc4);
/* Step 1 check if the special torsional parameter exists or not */
									suc = torsion(name1, name2, name3, name4,
											tmpc1, tmpc2, tmpc3, tmpc4, torsionparmnum, 0);
									if(suc == 1) continue;
/* Step 2 check special torsional parameters using equal atom types */
                                                                        if (suc == 0) {
                                                                        	besttorid = -1;
                                                                                pid1 = parmid[j]; /*central atoms in the first two layers of loops */
                                                                                pid2 = parmid[k];
                                                                                pid3 = parmid[i];
                                                                                pid4 = parmid[l];
                                                                                for(m=0; m<parm[pid1].nequa; m++) {
											if(suc == 1) break;
                                                                                        pid5 = parm[pid1].equa[m].pid;
                                                                                        strcpy(tmpc5, parm[pid5].atomtype);
                                                                                        for(n=0;n<parm[pid2].nequa; n++) {
												if(suc == 1) break;
                                                                                                pid6 = parm[pid2].equa[n].pid;
                                                                                                strcpy(tmpc6, parm[pid6].atomtype);
                                                                                                for(p=0;p<parm[pid3].nequa; p++) {
													if(suc == 1) break;
                                                                                                        pid7 = parm[pid3].equa[p].pid;
                                                                                                        strcpy(tmpc7, parm[pid7].atomtype);
                                                                                                        for(q=0;q<parm[pid4].nequa; q++) {
														if(suc == 1) break;
                                                                                                                if(m==0 && n==0 && p== 0 && q==0) continue;
                                                                                                                pid8 = parm[pid4].equa[q].pid;
                                                                                                                strcpy(tmpc8, parm[pid8].atomtype);
                                                                                                                suc2= torsion(name1, name2, name3, name4, tmpc7, tmpc5, tmpc6, tmpc8, torsionparmnum2, 1);
														if(suc2== 1) {
															bestscore = 0;
															suc = 1;
														}
                                                                                                        }
                                                                                                }
                                                                                        }
                                                                                }
										if(suc == 1 && besttorid >=0) {
											print_torsion(besttorid, name1, name2, name3, name4);
											besttorid = -1;
											continue;
										}
                                                                        }
									
/* Step 3. check general torsional angle terms*/
									if(suc == 0) {
										suc = torsion(name1, name2, name3, name4,
											"X", tmpc2, tmpc3, "X", torsionparmnum2, 0);
										if(suc == 1) continue;
									}

/* Step 4 check general torsional angle terms using equal atom types*/
                                                                        if (suc == 0) {
                                                                        	besttorid = -1;
                                                                                pid2 = parmid[j];
                                                                                pid3 = parmid[k];
                                                                                for(n=0;n<parm[pid2].nequa; n++) {
                                                                                        pid6 = parm[pid2].equa[n].pid;
                                                                                        strcpy(tmpc6, parm[pid6].atomtype);
                                                                                        for(p=0;p<parm[pid3].nequa; p++) {
												if(n==0 && p==0) continue;
                                                                                                pid7 = parm[pid3].equa[p].pid;
                                                                                                strcpy(tmpc7, parm[pid7].atomtype);
                                                                                                suc2= torsion(name1, name2, name3, name4, "X", tmpc6, tmpc7, "X", torsionparmnum2, 1);
												if(suc2== 1) {
													bestscore = 0;
													suc = 1;
												}
                                                                                        }
                                                                                }
										if(suc == 1 && besttorid >=0) {
											print_torsion(besttorid, name1, name2, name3, name4);
											besttorid = -1;
											continue;
										}
                                                                        }
/* Step 5. check special torsional parameters using corresponding atom types */
									if (suc == 0) {
										pid1 = parmid[i];	
										pid2 = parmid[j];	
										pid3 = parmid[k];	
										pid4 = parmid[l];	
                                                                		bestscore = INITSCORE;
										besttorid = -1;
										for(m=0; m<parm[pid1].ncorr; m++) {
											pid5 = parm[pid1].corr[m].pid;
                                                                        		strcpy(tmpc5, parm[pid5].atomtype);
                                                                        		score1= parm[pid1].corr[m].tor ; 
                                                                        		for(n=0;n<parm[pid2].ncorr; n++) {
												pid6 = parm[pid2].corr[n].pid;
                                                                        			strcpy(tmpc6, parm[pid6].atomtype);
                                                                        			score2= parm[pid2].corr[n].ctor * wt.TOR_CTR; 
                                                                        			for(p=0;p<parm[pid3].ncorr; p++) {
													pid7 = parm[pid3].corr[p].pid;
                                                                        				strcpy(tmpc7, parm[pid7].atomtype);
                                                                        				score3= parm[pid3].corr[p].ctor * wt.TOR_CTR; 
                                                                        				for(q=0;q<parm[pid4].ncorr; q++) {
														if(parm[pid1].corr[m].type <= 1 && 
                                                           				   			   parm[pid2].corr[n].type <= 1 &&
                                                           				   			   parm[pid3].corr[p].type <= 1 && 
                                                           				   			   parm[pid4].corr[q].type <= 1) 
															continue;	
														pid8 = parm[pid4].corr[q].pid;
                                                                        					strcpy(tmpc8, parm[pid8].atomtype);
                                                                        					score4= parm[pid4].corr[q].tor; 
														score = score1 + score2 + score3 + score4;
														score += wt.GROUP;
														equtype_penalty_score = equtype_penalty(pid2, pid3, pid6, pid7);	
														score += equtype_penalty_score;
														if(parm[pid5].group == parm[pid6].group && 
                                                                                                                   parm[pid5].group == parm[pid7].group &&
                                                                                                                   parm[pid5].group == parm[pid8].group)
															score -= wt.GROUP;
														if(score < bestscore) {
															suc2= torsion(name1, name2, name3, name4, tmpc5, tmpc6, tmpc7, tmpc8, torsionparmnum2, 1);
															if(suc2== 1) {
																bestscore = score;
																suc = 1;
															}
														}
													}
												}
											}
										}
										if(suc == 1 && besttorid >=0) {
											print_torsion(besttorid, name1, name2, name3, name4);
											besttorid = -1;
											continue;
										}
									}
/* Step 6. check general torsional parameters using corresponding atom types */
                                                                        if (suc == 0) {
                                                                                pid2 = parmid[j];
                                                                                pid3 = parmid[k];
                                                                                bestscore = INITSCORE;
										besttorid = -1;
                                                                                for(n=0;n<parm[pid2].ncorr; n++) {
                                                                                        pid6 = parm[pid2].corr[n].pid;
                                                                                        strcpy(tmpc6, parm[pid6].atomtype);
                                                                                        score2= parm[pid2].corr[n].ctor * wt.TOR_CTR;
                                                                                        for(p=0;p<parm[pid3].ncorr; p++) {
                                                                                                pid7 = parm[pid3].corr[p].pid;
												if(parm[pid2].corr[n].type <= 1 && 
                                                           				   	   parm[pid3].corr[p].type <= 1) 
													continue;	
                                                                                                strcpy(tmpc7, parm[pid7].atomtype);
                                                                                                score3= parm[pid3].corr[p].ctor * wt.TOR_CTR;
                                                                                                score = score2 + score3;
												equtype_penalty_score = equtype_penalty(pid2, pid3, pid6, pid7);	
												score += equtype_penalty_score;
												if(parm[pid6].group != parm[pid7].group) 
													score += wt.GROUP;
                                                                                                if(score < bestscore) {
                                                                                                        suc2= torsion(name1, name2, name3, name4, "X", tmpc6, tmpc7, "X", torsionparmnum2, 1);
													if(suc2== 1) {
														bestscore = score;
														suc = 1;
													}
                                                                                                }
                                                                                        }
                                                                                }
										if(suc == 1 && besttorid >=0) {
											print_torsion(besttorid, name1, name2, name3, name4);
											besttorid = -1;
											continue;
										}
                                                                        }

									if (suc == 0) {
										fprintf(fpout,
												"%-2s-%-2s-%-2s-%-2s%4d%9.3lf%14.3lf%16.3lf",
												name1, name2, name3, name4,
												1, 0.0, 0.0, 0.0);
										fprintf(fpout, "      %s\n",
												"ATTN, need revision");
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
}

void chk_improper(void)
{
	int i,j,k;
	int m,n,p,q;
	int suc;
	int tmpnum;
	int index1, index2, index4;
        char tmpc[5];
        char tmpc1[5], tmpc2[5], tmpc3[5], tmpc4[5];
        char name1[5], name2[5], name3[5], name4[5];
        int pid1, pid2, pid3, pid4;
        int pid5, pid6, pid7, pid8;
	double score, score1, score2, score3, score4;
	fprintf(fpout, "\n%s\n", "IMPROPER");

	for (i = 0; i < impropernum; i++) {
		suc = 0;
		strcpy(tmpc1, atom[improper[i].atid1].ambername);
		strcpy(tmpc2, atom[improper[i].atid2].ambername);
		strcpy(tmpc3, atom[improper[i].atid3].ambername);
		strcpy(tmpc4, atom[improper[i].atid4].ambername);

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
		strcpy(name1, tmpc1);
		strcpy(name2, tmpc2);
		strcpy(name3, tmpc3);
		strcpy(name4, tmpc4);
/*	step 1 check directly	*/
		for (j = 0; j < improperparmnum; j++)
			if(strcmp(improperparm[j].name1, name1) == 0 && 
			   strcmp(improperparm[j].name2, name2) == 0 && 
			   strcmp(improperparm[j].name3, name3) == 0 && 
			   strcmp(improperparm[j].name4, name4) == 0) {
				suc = 1;
				if(allparm_flag == 1) {
                                	fprintf(fpout, "%-2s-%-2s-%-2s-%-2s %11.1lf%15.1lf%12.1lf\n", name1,
                                        	name2, name3, name4, improperparm[j].force, improperparm[j].phase,
                                                improperparm[j].fterm);
				}
				break;
			}
		if(suc == 1) continue;
/* 	Step 2 check special improper torsional parameters using equal atom types */
                if(suc == 0) {                                         
                        pid1 = parmid[improper[i].atid1]; 
                        pid2 = parmid[improper[i].atid2];
                        pid3 = parmid[improper[i].atid3];
                        pid4 = parmid[improper[i].atid4];
                	bestimproperid = -1;
                        for(m=0; m<parm[pid1].nequa; m++) {
				if(suc == 1) break;
                                pid5 = parm[pid1].equa[m].pid;
                                strcpy(tmpc1, parm[pid5].atomtype);
                                for(n=0;n<parm[pid2].nequa; n++) {
					if(suc == 1) break;
                                        pid6 = parm[pid2].equa[n].pid;
                                        strcpy(tmpc2, parm[pid6].atomtype);
                                        for(p=0;p<parm[pid3].nequa; p++) {
						if(suc == 1) break;
                                                pid7 = parm[pid3].equa[p].pid;
                                                strcpy(tmpc3, parm[pid7].atomtype);
                                                for(q=0;q<parm[pid4].nequa; q++) {
							if(suc == 1) break;
                                                        if(m==0 && n==0 && p== 0 && q==0) continue;
                                                        pid8 = parm[pid4].equa[q].pid;
                                                        strcpy(tmpc4, parm[pid8].atomtype);
                                                        for (k = 0; k < improperparmnum2; k++) {
                                                                if (improperparm[k].numX > 0) continue;
                                                                if(strcmp(improperparm[k].name1, tmpc1) == 0 &&
                                                                   strcmp(improperparm[k].name2, tmpc2) == 0 &&
                                                                   strcmp(improperparm[k].name3, tmpc3) == 0 &&
                                                                   strcmp(improperparm[k].name4, tmpc4) == 0) {
                                                                        bestimproperid = k;
                                                                        bestscore = 0;
                                                                        suc = 1;
                                                                        break;
                                                        	}
                                                	}
                                        	}
                                	}
                        	}
                	}
                	if(suc == 1 && bestimproperid >= 0) {
                        	if(output_improper_flag == 1 || allparm_flag == 1) {
                                	fprintf(fpout, "%-2s-%-2s-%-2s-%-2s %11.1lf%15.1lf%12.1lf", name1,
                                                	name2, name3, name4, improperparm[bestimproperid].force,
                                                	improperparm[bestimproperid].phase,
                                                	improperparm[bestimproperid].fterm);
                                	fprintf(fpout, "          Same as %-2s-%-2s-%-2s-%-2s, penalty score=%5.1lf)\n",
                                                	improperparm[bestimproperid].name1,
                                                	improperparm[bestimproperid].name2,
                                                	improperparm[bestimproperid].name3,
                                                	improperparm[bestimproperid].name4, bestscore);
                        	}               
                        	strcpy(improperparm[improperparmnum].name1, name1);
                        	strcpy(improperparm[improperparmnum].name2, name2);    
                        	strcpy(improperparm[improperparmnum].name3, name3);    
                        	strcpy(improperparm[improperparmnum].name4, name4);
                        	improperparm[improperparmnum].phase = improperparm[bestimproperid].phase;
                        	improperparm[improperparmnum].fterm = improperparm[bestimproperid].fterm;
                        	improperparm[improperparmnum].force = improperparm[bestimproperid].force;
                        	improperparmnum++;
                        	if (improperparmnum >= maximproperparm) {
                                	maximproperparm += MAX_FF_IMPROPER;
                                	improperparm = (IMPROPER *) realloc(improperparm, sizeof(IMPROPER) * maximproperparm);
                                	if (improperparm == NULL) {
                                        	fprintf(stdout, "memory allocation error for *improperparm\n");
                                        	exit(1);
                                	} 
                        	}
				bestimproperid = -1;
				continue;
			}
                }    
/*	step 3 considering general atom types */
		if (suc == 0) {
			bestimproperid = -1;
                        bestscore = INITSCORE;
			for (j = 0; j < improperparmnum2; j++) {
				if(improperparm[j].numX == 0) continue;
				if(strcmp(improperparm[j].name1, name1) != 0 && strcmp(improperparm[j].name1, "X") != 0) 
					continue;	
				if(strcmp(improperparm[j].name2, name2) != 0 && strcmp(improperparm[j].name2, "X") != 0) 
					continue;	
				if(strcmp(improperparm[j].name3, name3) != 0 && strcmp(improperparm[j].name3, "X") != 0) 
					continue;	
				if(strcmp(improperparm[j].name4, name4) != 0 && strcmp(improperparm[j].name4, "X") != 0) 
					continue;	
				score1 = 0;
				score2 = 0;
				score3 = 0;
				score4 = 0;
				if(strcmp(improperparm[j].name1, "X") == 0) score1 = wt.X;
				if(strcmp(improperparm[j].name2, "X") == 0) score2 = wt.X;
				if(strcmp(improperparm[j].name3, "X") == 0) score3 = wt.X3;
				if(strcmp(improperparm[j].name4, "X") == 0) score4 = wt.X;
				score = score1 + score2 + score3 + score4;
				if(score < bestscore) {
					bestimproperid = j;
					bestscore = score;
					suc = 1;
				}
			}
                	if(suc == 1 && bestimproperid >= 0) {
                        	if(output_improper_flag == 1 || allparm_flag == 1) {
                               		fprintf(fpout, "%-2s-%-2s-%-2s-%-2s %11.1lf%15.1lf%12.1lf", name1,
                                        	      	name2, name3, name4, improperparm[bestimproperid].force,
                                              		improperparm[bestimproperid].phase,
                                              		improperparm[bestimproperid].fterm);
                               		fprintf(fpout, "          Using general improper torsional angle %2s-%2s-%2s-%2s, penalty score=%5.1lf)\n",
                                       	        	improperparm[bestimproperid].name1,
                                       		       	improperparm[bestimproperid].name2,
                                               		improperparm[bestimproperid].name3,
                                               		improperparm[bestimproperid].name4, bestscore);
                        	}
                       	 	strcpy(improperparm[improperparmnum].name1, name1);
                        	strcpy(improperparm[improperparmnum].name2, name2);
                        	strcpy(improperparm[improperparmnum].name3, name3);
                        	strcpy(improperparm[improperparmnum].name4, name4);
                        	improperparm[improperparmnum].phase = improperparm[bestimproperid].phase;
                        	improperparm[improperparmnum].fterm = improperparm[bestimproperid].fterm;
                        	improperparm[improperparmnum].force = improperparm[bestimproperid].force;
                        	improperparmnum++;
                        	if (improperparmnum >= maximproperparm) {
                                	maximproperparm += MAX_FF_IMPROPER;
                                	improperparm = (IMPROPER *) realloc(improperparm, sizeof(IMPROPER) * maximproperparm);
                                	if (improperparm == NULL) {
                                        	fprintf(stdout, "memory allocation error for *improperparm\n");
                                        	exit(1);
                                	}
                        	}
				bestimproperid = -1;
				continue;
               		}
		}
/*	Step 4 considering both equal atom types and general terms	*/
                if(suc == 0) {          
                        pid1 = parmid[improper[i].atid1];
                        pid2 = parmid[improper[i].atid2];
                        pid3 = parmid[improper[i].atid3];
                        pid4 = parmid[improper[i].atid4];
                        bestscore = INITSCORE;  
                        bestimproperid = -1;            
                        for(m=0; m<parm[pid1].nequa; m++) {
                                pid5 = parm[pid1].equa[m].pid;
                                strcpy(tmpc1, parm[pid5].atomtype);
                                for(n=0;n<parm[pid2].nequa; n++) {
                                        pid6 = parm[pid2].equa[n].pid;
                                        strcpy(tmpc2, parm[pid6].atomtype);
                                        for(p=0;p<parm[pid3].nequa; p++) {
                                                pid7 = parm[pid3].equa[p].pid; 
                                                strcpy(tmpc3, parm[pid7].atomtype);
                                                for(q=0;q<parm[pid4].nequa; q++) {
							if(m==0 && n==0 && p==0 && q==0) continue;  
                                                        pid8 = parm[pid4].equa[q].pid;
                                                        strcpy(tmpc4, parm[pid8].atomtype);
                                                        for (k = 0; k < improperparmnum2; k++){
                                                                if (improperparm[k].numX <= 0) continue;
                                                                if(strcmp(improperparm[k].name1, tmpc1) != 0 && strcmp(improperparm[k].name1, "X") != 0)
                                                                        continue;
                                                                if(strcmp(improperparm[k].name2, tmpc2) != 0 && strcmp(improperparm[k].name2, "X") != 0)
                                                                        continue;
                                                                if(strcmp(improperparm[k].name3, tmpc3) != 0 && strcmp(improperparm[k].name3, "X") != 0)
                                                                        continue;       
                                                                if(strcmp(improperparm[k].name4, tmpc4) != 0 && strcmp(improperparm[k].name4, "X") != 0)
                                                                        continue;  
								score1 = 0;
								score2 = 0;
								score3 = 0;
								score4 = 0;
                                                                if(strcmp(improperparm[k].name1, "X") == 0) score1= wt.X;
                                                                if(strcmp(improperparm[k].name2, "X") == 0) score2= wt.X;
                                                                if(strcmp(improperparm[k].name3, "X") == 0) score3= wt.X3;
                                                                if(strcmp(improperparm[k].name4, "X") == 0) score4= wt.X;
								score = score1 + score2 + score3 + score4;
                                                                if(score < bestscore) {
                                                                        bestimproperid = k;
                                                                        bestscore = score;
                                                                        suc = 1;
                                                                }
                                                                break;
                                                        }
                                                }
                                        }
                                }
                        }
                	if(suc == 1 && bestimproperid >= 0) {
                        	if(output_improper_flag == 1 || allparm_flag == 1) { 
                                	fprintf(fpout, "%-2s-%-2s-%-2s-%-2s %11.1lf%15.1lf%12.1lf", name1,
                                                	name2, name3, name4, improperparm[bestimproperid].force,
                                                	improperparm[bestimproperid].phase,
                                                	improperparm[bestimproperid].fterm);
                                	fprintf(fpout, "          Same as %-2s-%-2s-%-2s-%-2s, penalty score=%5.1lf (use general term))\n",
                                                	improperparm[bestimproperid].name1,
                                                	improperparm[bestimproperid].name2,
                                                	improperparm[bestimproperid].name3,
                                                	improperparm[bestimproperid].name4, bestscore);
                        	}       
                        	strcpy(improperparm[improperparmnum].name1, name1); 
                        	strcpy(improperparm[improperparmnum].name2, name2); 
                        	strcpy(improperparm[improperparmnum].name3, name3); 
                        	strcpy(improperparm[improperparmnum].name4, name4); 
                        	improperparm[improperparmnum].phase = improperparm[bestimproperid].phase;
                        	improperparm[improperparmnum].fterm = improperparm[bestimproperid].fterm;
                        	improperparm[improperparmnum].force = improperparm[bestimproperid].force;
                        	improperparmnum++;
                        	if (improperparmnum >= maximproperparm) {
                                	maximproperparm += MAX_FF_IMPROPER;
                                	improperparm = (IMPROPER *) realloc(improperparm, sizeof(IMPROPER) * maximproperparm);
                                	if (improperparm == NULL) { 
                                        	fprintf(stdout, "memory allocation error for *improperparm\n");
                                        	exit(1);
                                	}       
                        	}       
				bestimproperid = -1;
				continue;
                	}     
		}

/*	Step 5 considering corresponding atom types for specific improper parameters	*/
		if(suc == 0) {
                	pid1 = parmid[improper[i].atid1];
                	pid2 = parmid[improper[i].atid2];
                	pid3 = parmid[improper[i].atid3];
                	pid4 = parmid[improper[i].atid4];
                        bestscore = INITSCORE;
			bestimproperid = -1;
                        for(m=0; m<parm[pid1].ncorr; m++) {
                                pid5 = parm[pid1].corr[m].pid;
                                strcpy(tmpc1, parm[pid5].atomtype);
                                score1= parm[pid1].corr[m].improper;
                                for(n=0;n<parm[pid2].ncorr; n++) {
                                        pid6 = parm[pid2].corr[n].pid;
                                        strcpy(tmpc2, parm[pid6].atomtype);
                                	score2= parm[pid2].corr[n].improper;
                                        for(p=0;p<parm[pid3].ncorr; p++) {
                                                pid7 = parm[pid3].corr[p].pid;
                                                strcpy(tmpc3, parm[pid7].atomtype);
                                		score3= parm[pid3].corr[p].improper * wt.IMPROPER;
                                                for(q=0;q<parm[pid4].ncorr; q++) {
                                                	if(parm[pid1].corr[m].type <= 1 &&
                                                           parm[pid2].corr[n].type <= 1 &&
                                                           parm[pid3].corr[p].type <= 1 &&
                                                           parm[pid4].corr[q].type <= 1)
                                                        	continue;
                                                        pid8 = parm[pid4].corr[q].pid;
                                                        strcpy(tmpc4, parm[pid8].atomtype);
                                			score4= parm[pid4].corr[q].improper;
                                                        score = score1 + score2 + score3 + score4;
							score += wt.GROUP;
							if(parm[pid5].group == parm[pid6].group && 
                                   			   parm[pid5].group == parm[pid7].group &&
                                   			   parm[pid5].group == parm[pid8].group)
			   	       			   score -= wt.GROUP;
                                                        if(score < bestscore) { 
								for (k = 0; k < improperparmnum2; k++)
									if (improperparm[k].numX > 0) continue; 
									if(strcmp(improperparm[k].name1, tmpc1) == 0 && 
			   						   strcmp(improperparm[k].name2, tmpc2) == 0 && 
			   						   strcmp(improperparm[k].name3, tmpc3) == 0 && 
			   						   strcmp(improperparm[k].name4, tmpc4) == 0) {
										bestimproperid = k;
                                                                		bestscore = score;
										suc = 1;
										break;
								}
							}
                                                }
                                        }
                                }
                        }
                	if(suc == 1 && bestimproperid >= 0) {
                        	if(output_improper_flag == 1 || allparm_flag == 1) {
                               		fprintf(fpout, "%-2s-%-2s-%-2s-%-2s %11.1lf%15.1lf%12.1lf", name1,
                                       		       	name2, name3, name4, improperparm[bestimproperid].force,
                                              		improperparm[bestimproperid].phase,
                                              		improperparm[bestimproperid].fterm);
                               		fprintf(fpout, "          Same as %-2s-%-2s-%-2s-%-2s, penalty score=%5.1lf)\n",
                                               		improperparm[bestimproperid].name1,
                                               		improperparm[bestimproperid].name2,
                                               		improperparm[bestimproperid].name3,
                                               		improperparm[bestimproperid].name4, bestscore);
                        	}
                        	strcpy(improperparm[improperparmnum].name1, name1);
                        	strcpy(improperparm[improperparmnum].name2, name2);
                        	strcpy(improperparm[improperparmnum].name3, name3);
                        	strcpy(improperparm[improperparmnum].name4, name4);
                        	improperparm[improperparmnum].phase = improperparm[bestimproperid].phase;
                        	improperparm[improperparmnum].fterm = improperparm[bestimproperid].fterm;
                        	improperparm[improperparmnum].force = improperparm[bestimproperid].force;
                        	improperparmnum++;
                        	if (improperparmnum >= maximproperparm) {
                               		maximproperparm += MAX_FF_IMPROPER;
                                	improperparm = (IMPROPER *) realloc(improperparm, sizeof(IMPROPER) * maximproperparm);
                                	if (improperparm == NULL) {
                                        	fprintf(stdout, "memory allocation error for *improperparm\n");
                                        	exit(1);
                                	}
                        	}
				bestimproperid = -1;
				continue;
                	}
		}
/*	Step6  considering corresponding atom types for general improper parameters	*/
		if(suc == 0) {
                	pid1 = parmid[improper[i].atid1];
                	pid2 = parmid[improper[i].atid2];
                	pid3 = parmid[improper[i].atid3];
                	pid4 = parmid[improper[i].atid4];
                        bestscore = INITSCORE;
			bestimproperid = -1;
                        for(m=0; m<parm[pid1].ncorr; m++) {
                                pid5 = parm[pid1].corr[m].pid;
                                strcpy(tmpc1, parm[pid5].atomtype);
                                score1= parm[pid1].corr[m].improper ; 
                                for(n=0;n<parm[pid2].ncorr; n++) {
                                        pid6 = parm[pid2].corr[n].pid;
                                        strcpy(tmpc2, parm[pid6].atomtype);
                                        score2= parm[pid2].corr[n].improper;
                                        for(p=0;p<parm[pid3].ncorr; p++) {
                                                pid7 = parm[pid3].corr[p].pid;
                                                strcpy(tmpc3, parm[pid7].atomtype);
                                                score3= parm[pid3].corr[p].improper + wt.IMPROPER;
                                                for(q=0;q<parm[pid4].ncorr; q++) {
                                                	if(parm[pid1].corr[m].type <= 1 &&
                                                           parm[pid2].corr[n].type <= 1 &&
                                                           parm[pid3].corr[p].type <= 1 &&
                                                           parm[pid4].corr[q].type <= 1)
                                                        	continue;
                                                        pid8 = parm[pid4].corr[q].pid;
                                                        strcpy(tmpc4, parm[pid8].atomtype);
                                                        score4= parm[pid4].corr[q].improper; 
                                                        score = score1 + score2 + score3 + score4;
							score += wt.GROUP;
							if(parm[pid5].group == parm[pid6].group && 
                                   			   parm[pid5].group == parm[pid7].group &&
                                   			   parm[pid5].group == parm[pid8].group)
			   	       			   score -= wt.GROUP;
							for (k = 0; k < improperparmnum2; k++){
								if (improperparm[k].numX <= 0) continue; 
								if(strcmp(improperparm[k].name1, tmpc1) != 0 && strcmp(improperparm[k].name1, "X") != 0) 
									continue;	
								if(strcmp(improperparm[k].name2, tmpc2) != 0 && strcmp(improperparm[k].name2, "X") != 0) 
									continue;	
								if(strcmp(improperparm[k].name3, tmpc3) != 0 && strcmp(improperparm[k].name3, "X") != 0) 
									continue;	
								if(strcmp(improperparm[k].name4, tmpc4) != 0 && strcmp(improperparm[k].name4, "X") != 0) 
									continue;	
								if(strcmp(improperparm[k].name1, "X") == 0) score += wt.X;
								if(strcmp(improperparm[k].name2, "X") == 0) score += wt.X;
								if(strcmp(improperparm[k].name3, "X") == 0) score += wt.X3;
								if(strcmp(improperparm[k].name4, "X") == 0) score += wt.X;
								if(score < bestscore) {
									bestimproperid = k;
                                                               		bestscore = score;
									suc = 1;
								}
								break;
							}
                                                }
                                        }
                                }
                        }
                	if(suc == 1 && bestimproperid >= 0) {
                        	if(output_improper_flag == 1 || allparm_flag == 1) {
                               		fprintf(fpout, "%-2s-%-2s-%-2s-%-2s %11.1lf%15.1lf%12.1lf", name1,
                                       	        name2, name3, name4, improperparm[bestimproperid].force,
                                       	        improperparm[bestimproperid].phase,
                                      	        improperparm[bestimproperid].fterm);
                               		fprintf(fpout, "          Same as %-2s-%-2s-%-2s-%-2s, penalty score=%5.1lf (use general term))\n",
                                       	        improperparm[bestimproperid].name1,
                                       	        improperparm[bestimproperid].name2,
                                       	        improperparm[bestimproperid].name3,
                                       	        improperparm[bestimproperid].name4, bestscore);
                        	}
                        	strcpy(improperparm[improperparmnum].name1, name1);
                        	strcpy(improperparm[improperparmnum].name2, name2);
                        	strcpy(improperparm[improperparmnum].name3, name3);
                        	strcpy(improperparm[improperparmnum].name4, name4);
                        	improperparm[improperparmnum].phase = improperparm[bestimproperid].phase;
                        	improperparm[improperparmnum].fterm = improperparm[bestimproperid].fterm;
                        	improperparm[improperparmnum].force = improperparm[bestimproperid].force;
                        	improperparmnum++;
                        	if (improperparmnum >= maximproperparm) {
                                	maximproperparm += MAX_FF_IMPROPER;
                                	improperparm = (IMPROPER *) realloc(improperparm, sizeof(IMPROPER) * maximproperparm);
                                	if (improperparm == NULL) {
                                        	fprintf(stdout, "memory allocation error for *improperparm\n");
                                        	exit(1);
                                	}
                        	}
				bestimproperid = -1;
				continue;
                	}
		}
		if (suc == 0) {
			fprintf(fpout, "%-2s-%-2s-%-2s-%-2s %11.1lf%15.1lf%12.1lf", name1,
					name2, name3, name4, 1.1, 180.0, 2.0);
			fprintf(fpout, "          Using the default value\n");
			strcpy(improperparm[improperparmnum].name1, name1);
			strcpy(improperparm[improperparmnum].name2, name2);
			strcpy(improperparm[improperparmnum].name3, name3);
			strcpy(improperparm[improperparmnum].name4, name4);
			improperparm[improperparmnum].phase = 180.0;
			improperparm[improperparmnum].fterm = 2.0;
			improperparm[improperparmnum].force = 1.1;
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

void cleanup_frcmod(char * filename) {
	FILE *fp1, *fp2;
	char command[MAXCHAR] = "cp -f ";
	int status;
	typedef struct {
        	char name[10];
	} ATOMTYPE_CL;
	typedef struct {
        	char name1[10];
        	char name2[10];
	} BOND_CL;
	typedef struct {
        	char name1[10];
        	char name2[10];
        	char name3[10];
	} ANGLE_CL;
	typedef struct {
        	char name1[10];
        	char name2[10];
        	char name3[10];
        	char name4[10];
		double fterm;
	} TORSION_CL;
	ATOMTYPE_CL at[MAX_FF_ATOMTYPE];
	int atnum = 0;
	BOND_CL bond[MAX_FF_BOND];
	int bondnum = 0;
	ANGLE_CL angle[MAX_FF_ANGLE];
	int anglenum = 0;
	TORSION_CL tor[MAX_FF_TORSION];
	int tornum = 0;
	TORSION_CL improp[MAX_FF_IMPROPER];
	int impropnum = 0;
	ATOMTYPE_CL vdw[MAX_FF_VDW];	
	int vdwnum = 0;
	int writeflag = 1;
	int typeflag = 0;
	int tmpint;
	char name1[10];
	char name2[10];
	char name3[10];
	char name4[10];
	char line[MAXCHAR];
	char line2[MAXCHAR];
	char tmpchar[MAXCHAR];
	double tmpf1, tmpf2, tmpf3;

	strcat(command, filename);
	strcat(command, " ANTECHAMBER.FRCMOD");
	status = system(command);
        if(status != 0) {
                fprintf(stdout, "Error: cannot run \"%s\" in cleanup_frcmod function properly, exit\n", command);
                exit(1);
        }
	
        if ((fp1 = fopen("ANTECHAMBER.FRCMOD", "r")) == NULL) {
		fprintf(stdout, "Error: Cannot open the redundant frcmod file - ANTECHAMBER.FRCMOD to do clean up\n"); 
		exit(1);
        }
        if ((fp2 = fopen(filename, "w")) == NULL) {
		fprintf(stdout, "Error: Cannot open the frcmod file %s to write out\n", filename); 
		exit(1);
        }
	for(;;) {
                if (fgets(line, MAXCHAR, fp1) == NULL) break;
                sscanf(line, "%s", tmpchar);
                if (strcmp("MASS", tmpchar) == 0) 
			typeflag = 1;
                if (strcmp("BOND", tmpchar) == 0) 
			typeflag = 2;
                if (strcmp("ANGLE", tmpchar) == 0) 
			typeflag = 3;
                if (strcmp("DIHE", tmpchar) == 0) 
			typeflag = 4;
                if (strcmp("IMPROPER", tmpchar) == 0) 
			typeflag = 5;
                if (strcmp("NONBON", tmpchar) == 0) 
			typeflag = 6;
		if(strlen(line) <= 1) {
			typeflag = 0;
			writeflag = 1;
		}
		if(typeflag == 1) {
			writeflag = 1;
			for(i=0;i<atnum;i++)  {
				if(strcmp(tmpchar, at[i].name) == 0) {
					writeflag = 0;
					break;
				}
			}
			if(writeflag == 1) {
				strcpy(at[atnum].name, tmpchar);
				atnum++;
			}
		}

		if(typeflag == 2) {
			writeflag = 1;
			strcpy(line2, line);
			for(i=0;i<6;i++)  /* we do not want to delete other '-' after the 6th columns */
				if(line2[i] == '-') line2[i] = ' ';
			sscanf(line2, "%s%s", name1, name2);
			for(i=0;i<bondnum;i++)  {
				if(strcmp(name1, bond[i].name1) == 0 && strcmp(name2, bond[i].name2) == 0) {
					writeflag = 0;
					break;
				}
				if(strcmp(name1, bond[i].name2) == 0 && strcmp(name2, bond[i].name1) == 0) {
					writeflag = 0;
					break;
				}
			}
			if(writeflag == 1) {
				strcpy(bond[bondnum].name1, name1);
				strcpy(bond[bondnum].name2, name2);
				bondnum++;
			}
		}

                if(typeflag == 3) {
                        writeflag = 1;
			strcpy(line2, line);
			for(i=0;i<9;i++)  /* we do not want to delete other '-' after the 9th columns */
				if(line2[i] == '-') line2[i] = ' ';
                        sscanf(line2, "%s%s%s", name1, name2,name3);
                        for(i=0;i<anglenum;i++)  {
                                if(strcmp(name1, angle[i].name1) == 0 && strcmp(name2, angle[i].name2) == 0 && strcmp(name3, angle[i].name3) == 0 ) {
                                        writeflag = 0;
                                        break;
                                }
                                if(strcmp(name1, angle[i].name3) == 0 && strcmp(name2, angle[i].name2) == 0 && strcmp(name3, angle[i].name1) == 0 ) {
                                        writeflag = 0;
                                        break;
                                }
                        }
                        if(writeflag == 1) {
                                strcpy(angle[anglenum].name1, name1);
                                strcpy(angle[anglenum].name2, name2);
                                strcpy(angle[anglenum].name3, name3);
                                anglenum++;
                        }
                }

                if(typeflag == 4) {
                        writeflag = 1;
			strcpy(line2, line);
			for(i=0;i<12;i++)  /* we do not want to delete '-' after 12th column, esp., '-' before fterm */
				if(line2[i] == '-') line2[i] = ' ';
                        sscanf(line2, "%s%s%s%s%d%lf%lf%lf", name1, name2,name3,name4, &tmpint, &tmpf1, &tmpf2, &tmpf3);
                        for(i=0;i<tornum;i++)  {
                                if(strcmp(name1, tor[i].name1) == 0 && strcmp(name2, tor[i].name2) == 0 && 
				   strcmp(name3, tor[i].name3) == 0 && strcmp(name4, tor[i].name4) == 0 && 
				   tor[i].fterm == tmpf3) {
                                        writeflag = 0;
                                        break;
                                }
                                if(strcmp(name1, tor[i].name4) == 0 && strcmp(name2, tor[i].name3) == 0 && 
				   strcmp(name3, tor[i].name2) == 0 && strcmp(name4, tor[i].name1) == 0 &&
                                   tor[i].fterm == tmpf3) {
                                        writeflag = 0;
                                        break;
                                }
                        }
                        if(writeflag == 1) {
                                strcpy(tor[tornum].name1, name1);
                                strcpy(tor[tornum].name2, name2);
                                strcpy(tor[tornum].name3, name3);
                                strcpy(tor[tornum].name4, name4);
				tor[tornum].fterm = tmpf3;
                                tornum++;
                        }
                }

                if(typeflag == 5) {
                        writeflag = 1;
			strcpy(line2, line);
			for(i=0;i<12;i++)  /* we do not want to delete '-' after 12th column*/
				if(line2[i] == '-') line2[i] = ' ';
                        sscanf(line2, "%s%s%s%s", name1, name2,name3,name4);
                        for(i=0;i<impropnum;i++)  {
                                if(strcmp(name1, improp[i].name1) == 0 && strcmp(name2, improp[i].name2) == 0 &&
                                   strcmp(name3, improp[i].name3) == 0 && strcmp(name4, improp[i].name4) == 0) {
                                        writeflag = 0;
                                        break;
                                }
                                if(strcmp(name1, improp[i].name4) == 0 && strcmp(name2, improp[i].name3) == 0 &&
                                   strcmp(name3, improp[i].name2) == 0 && strcmp(name4, improp[i].name1) == 0) {
                                        writeflag = 0;
                                        break;
                                }
                        }
                        if(writeflag == 1) {
                                strcpy(improp[impropnum].name1, name1);
                                strcpy(improp[impropnum].name2, name2);
                                strcpy(improp[impropnum].name3, name3);
                                strcpy(improp[impropnum].name4, name4);
                                impropnum++;
                        }
                }

		if(typeflag == 6) {
			writeflag = 1;
			for(i=0;i<vdwnum;i++) 
				if(strcmp(tmpchar, vdw[i].name) == 0) {
					writeflag = 0;
					break;
				}
			if(writeflag == 1) {
				strcpy(vdw[vdwnum].name, tmpchar);
				vdwnum++;
			}
		}

		if(writeflag == 1) 
			fprintf(fp2, "%s", line);
	}
	fclose(fp1);
	fclose(fp2);
}

int main(int argc, char *argv[])
{
	int i;
	int format;
	FILE *fptmp;
	int overflow_flag = 0;			/*if overflow_flag ==1, reallocate memory */

	default_cinfo(&cinfo);
	default_minfo(&minfo);
    amberhome = (char *) getenv("AMBERHOME");
    if( amberhome == NULL ){
       fprintf( stdout, "AMBERHOME is not set!\n" );
       exit(1);
    }
    minfo.connect_file[0] = '\0';
    build_dat_path(minfo.connect_file, "CONNECT.TPL",
    	sizeof minfo.connect_file, 0);
	build_path(pfilename, "/dat/leap/parm/", "gaff.dat",
		sizeof pfilename, 0);

	if (strcmp(COLORTEXT, "YES") == 0 || strcmp(COLORTEXT, "yes") == 0) {
		if (argc == 2
			&& (strcmp(argv[1], "-h") == 0
				|| strcmp(argv[1], "-H") == 0)) {
			printf("[31mUsage: parmchk -i [0m input file name\n"
				   "[31m               -o [0m frcmod file name\n"
				   "[31m               -f [0m input file format (prepi, prepc, ac, mol2) \n"
				   "[31m               -p [0m ff parmfile\n"
				   "[31m               -pf[0m parmfile format, 1 for amber FF data file (the default) and 2 for additional force field parameter file\n"
				   "[31m               -c [0m atom type corresponding score file, default is PARMCHK.DAT\n" 
				   "[31m               -a [0m print out all force field parameters including those in the parmfile\n"
				   "[31m                  [0m can be 'Y' (yes) or 'N' (no) default is 'N' \n"
				   "[31m               -w [0m print out parameters that matching improper dihedral parameters\n"
				   "[31m                  [0m that contain 'X' in the force field parameter file, can be 'Y' (yes)\n"
				   "[31m                  [0m or 'N' (no), default is 'Y'\n");
			exit(1);
		}
		if (argc != 7 && argc != 9 && argc != 11 && argc != 13 && argc != 15 && argc != 17) {
			printf("[31mUsage: parmchk -i[0m input file name\n"
				   "[31m               -o[0m frcmod file name\n"
				   "[31m               -f[0m input file format (prepi, prepc, ac, mol2) \n"
				   "[31m               -p[0m ff parmfile\n"
				   "[31m               -pf[0m parmfile format, 1 for amber FF data file (the default) and 2 for additional force field parameter file\n"
				   "[31m               -c [0m atom type corresponding score file, default is PARMCHK.DAT\n" 
				   "[31m               -a[0m print out all force field parameters including those in the parmfile\n"
				   "[31m                 [0m can be 'Y' (yes) or 'N' (no) default is 'N' \n"
				   "[31m               -w[0m print out parameters that matching improper dihedral parameters\n"
				   "[31m                 [0m that contain 'X' in the force field parameter file, can be 'Y' (yes)\n"
				   "[31m                 [0m or 'N' (no), default is 'Y'\n");
			exit(1);
		}
	} else {
		if (argc == 2
			&& (strcmp(argv[1], "-h") == 0
				|| strcmp(argv[1], "-H") == 0)) {
			printf("Usage: parmchk -i   input file name\n");
			printf("               -o   frcmod file name\n");
			printf("               -f   input file format (prepi, prepc, ac, mol2) \n");
			printf("               -p   ff parmfile\n");
		        printf("	       -pf  parmfile format, 1 for amber FF data file (the default) and 2 for additional force field parameter file\n");
			printf("               -c   atom type corresponding score file \n");
			printf("                    (default is PARMCHK.DAT)\n");
		        printf("               -a   print out all force field parameters including those in the parmfile \n");
		        printf("                    can be 'Y' (yes) or 'N' (no) default is 'N' \n");
		        printf("               -w   print out parameters that matching improper dihedral parameters\n");
		        printf("		    that contain 'X' in the force field parameter file, can be 'Y' (yes)\n");
		        printf("		    or 'N' (no), default is 'Y'\n");
			exit(1);
		}
		if (argc != 7 && argc != 9 && argc != 11 && argc != 13) {
			printf("Usage: parmchk -i   input file name\n");
			printf("               -o   frcmod file name\n");
			printf("               -f   input file format (prepi, prepc, ac, mol2) \n");
			printf("               -p   ff parmfile\n");
		        printf("	       -pf  parmfile format, 1 for amber FF data file (the default) and 2 for additional force field parameter file\n");
			printf("               -c   atom type corresponding score file \n");
			printf("                    (default is PARMCHK.DAT)\n");
		        printf("               -a   print out all force field parameters including those in the parmfile \n");
		        printf("                    can be 'Y' (yes) or 'N' (no) default is 'N' \n");
		        printf("               -w   print out parameters that matching improper dihedral parameters\n");
		        printf("		    that contain 'X' in the force field parameter file, can be 'Y' (yes)\n");
		        printf("		    or 'N' (no), default is 'Y'\n");
			exit(1);
		}
	}
	format = -1;  /* input file format; -1 is an invalid format */
	cindex = 0;
	for (i = 1; i < argc; i += 2) {
		if (strcmp(argv[i], "-i") == 0)
			strcpy(ifilename, argv[i + 1]);
		if (strcmp(argv[i], "-p") == 0)
			strcpy(pfilename, argv[i + 1]);
		if (strcmp(argv[i], "-pf") == 0)
			pformat = atoi(argv[i+1]);
		if (strcmp(argv[i], "-o") == 0)
			strcpy(ofilename, argv[i + 1]);
		if (strcmp(argv[i], "-c") == 0) {
			strcpy(cfilename, argv[i + 1]);
			cindex = 1;
		}
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
		if (strcmp(argv[i], "-a") == 0) {
			if(argv[i + 1][0] == 'Y' ||argv[i + 1][0] == 'y') 
				allparm_flag = 1;
			if(argv[i + 1][0] == 'N' ||argv[i + 1][0] == 'n') 
				allparm_flag = 0;
		}
		if (strcmp(argv[i], "-w") == 0) {
			if(argv[i + 1][0] == 'Y' ||argv[i + 1][0] == 'y') 
				output_improper_flag = 1;
			if(argv[i + 1][0] == 'N' ||argv[i + 1][0] == 'n') 
				output_improper_flag = 0;
		}
	}

	if(pformat != 1 && pformat !=2) pformat = 1;
	if (cindex == 0) {
		build_dat_path(cfilename, "PARMCHK.DAT", sizeof cfilename, 0);
		cindex = 1;
	}
	if (cindex == 0){
		if ((fptmp = fopen("PARMCHK.DAT", "r")) != NULL) {
			strcpy(cfilename, "PARMCHK.DAT");
			cindex = 1;
			fclose(fptmp);
		} else {
			fprintf(stderr, "No parmchk corresponding score file exists, exit\n");
			exit(1);
		}
	}

	/*  memory allocation */
	/* initialize */
	maxparmnum = MAXPARM;
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

	switch (format) {
	case 0:
		overflow_flag = rprepi(ifilename, &atomnum, atom, &bondnum,
			bond_array, &cinfo, &minfo);
		break;
	case 1:
		overflow_flag = rprepc(ifilename, &atomnum, atom, &cinfo, &minfo);
		break;
	case 2:
		overflow_flag = rac(ifilename, &atomnum, atom, &bondnum,
			bond_array, &cinfo, &minfo);
		break;
	case 3:
		overflow_flag = rmol2(ifilename, &atomnum, atom, &bondnum,
			bond_array, &cinfo, &minfo, 1);
		break;
	default:
		printf("Error: invalid input file format"
			"(valid: prepi, prepc, ac, mol2)!\n");
		exit(1);
		break;
}

	if (overflow_flag) {
		cinfo.maxatom = atomnum + 10;
                /* warning bondnum is not correct if overflow_flag != 0 */
                /* so add to maxbond a fudge based on the number of atoms */
                cinfo.maxbond = bondnum + 4*atomnum + 10;
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
/*
		overflow_flag =
			connect(minfo.connect_file, atomnum, atom, &bondnum,
					bond_array, cinfo.maxbond);
*/
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

        parm = (PARM *) calloc(maxparmnum, sizeof(PARM));
        if (parm == NULL) {
                fprintf(stdout, "memory allocation error for *parm\n");
                exit(1);
        }
        parmid = (int *) malloc(sizeof(int) * atomnum);
        if (parmid == NULL) {
                fprintf(stdout, "memory allocation error for *parmid\n");
                exit(1);
        }
	improper = (IMPROPERID *) malloc(sizeof(IMPROPERID) * atomnum);
	if (improper == NULL) {
		fprintf(stdout, "memory allocation error for *improper\n");
		exit(1);
	}

	atomtype = (ATOMTYPE *) calloc(maxatomtype, sizeof(ATOMTYPE));
	if (atomtype == NULL) {
		fprintf(stdout, "memory allocation error for *atomtype\n");
		exit(1);
	}

	bondparm = (BOND_FF *) calloc(maxbondparm, sizeof(BOND_FF));
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

	if ((fpout = fopen(ofilename, "w")) == NULL) {
		fprintf(stdout, "Cannot open a file to write: %s, exit\n", ofilename);
		exit(1);
	}

	/* read in parameters */

	if (cindex == 1)
		read_parmchk_parm(cfilename);	/*atom type *corresponding file */
	assign_parmid();
	if(pformat == 1) readparm(pfilename);		/*principle parameter file */
	if(pformat == 2) readfrcmod(pfilename);	        /*principle parmaeter file in frcmod format*/

	if (format == 0 || format == 1)
		improper_id1(ifilename);
	if (format == 2 || format == 3)
		improper_id2();

	atomtypenum2 = atomtypenum; 
	vdwparmnum2 = vdwparmnum;
	bondparmnum2 = bondparmnum;
	angleparmnum2 = angleparmnum;
	torsionparmnum2 = torsionparmnum;
	improperparmnum2 = impropernum;
	
	chk_atomtype();
	chk_bond();
	chk_angle();
	chk_torsion();
	chk_improper();
	chk_vdw();
	fclose(fpout);
	if(allparm_flag == 1) 
		cleanup_frcmod(ofilename);
/*
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
*/
	return (0);
}
