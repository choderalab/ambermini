# include "common.h"
# include "define.h"
# include "utility.c"
char pfilename[MAXCHAR];
char mfilename[MAXCHAR];
char ofilename[MAXCHAR];
char line[MAXCHAR];
int i, j, k, l;
FILE *fp;

#define MAX_FF_ATOMTYPE 250
#define MAX_FF_VDW 250
#define MAX_FF_BOND 2000
#define MAX_FF_ANGLE 5000
#define MAX_FF_TORSION 1500
#define MAX_FF_IMPROPER 500

typedef struct {
	char name[5];
	char parm[90];
} ATOMTYPE;

typedef struct {
	char name1[5];
	char name2[5];
	char parm[90];
} BOND_FF;

typedef struct {
	char name1[5];
	char name2[5];
	char name3[5];
	char parm[90];
} ANGLE;

typedef struct {
	char name1[5];
	char name2[5];
	char name3[5];
	char name4[5];
	char parm[90];
} TORSION;

typedef struct {
	char name1[5];
	char name2[5];
	char name3[5];
	char name4[5];
	char parm[90];
} IMPROPER;

typedef struct {
	char name[5];
	char parm[90];
} VDW;

int atomtypenum;
int vdwparmnum;
int bondparmnum;
int angleparmnum;
int torsionparmnum;
int improperparmnum;

int maxatomtype;
int maxvdwparm;
int maxbondparm;
int maxangleparm;
int maxtorsionparm;
int maximproperparm;

ATOMTYPE *atomtype;
BOND_FF *bondparm;
ANGLE *angleparm;
TORSION *torsionparm;
IMPROPER *improperparm;
VDW *vdwparm;

int *atomtypeindex;
int *bondparmindex;
int *angleparmindex;
int *torsionparmindex;
int *improperparmindex;
int *vdwparmindex;

void readfrcmod(char *filename)
{
	int mindex = 0;
	int bindex = 0;
	int aindex = 0;
	int tindex = 0;
	int iindex = 0;
	int vindex = 0;
	FILE *fp;
	char line[MAXCHAR];
	if ((fp = fopen(filename, "r")) == NULL) {
		printf("\n Cannot open teh frcmod file %s in readfrcmod(), exit", filename);
		return;
	}
	for (;;) {
		if (fgets(line, MAXCHAR, fp) == NULL)
			break;
		if (strncmp("MASS", &line[0], 4) == 0) {
			mindex = 1;
			continue;
		}
		if (strncmp("BOND", &line[0], 4) == 0) {
			bindex = 1;
			continue;
		}
		if (strncmp("ANGL", &line[0], 4) == 0) {
			aindex = 1;
			continue;
		}
		if (strncmp("DIHE", &line[0], 4) == 0) {
			tindex = 1;
			continue;
		}
		if (strncmp("IMPR", &line[0], 4) == 0) {
			iindex = 1;
			continue;
		}
		if (strncmp("NONB", &line[0], 4) == 0) {
			vindex = 1;
			continue;
		}
		if (mindex == 1 && spaceline(line) == 1)
			mindex = 0;
		if (bindex == 1 && spaceline(line) == 1)
			bindex = 0;
		if (aindex == 1 && spaceline(line) == 1)
			aindex = 0;
		if (tindex == 1 && spaceline(line) == 1)
			tindex = 0;
		if (iindex == 1 && spaceline(line) == 1)
			iindex = 0;
		if (vindex == 1 && spaceline(line) == 1)
			vindex = 0;
		if (mindex == 1) {
			sscanf(line, "%s", atomtype[atomtypenum].name);
			strcpy(atomtype[atomtypenum].parm, line);
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
			strcpy(bondparm[bondparmnum].parm, line);
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
			angleparm[angleparmnum].name1[1] = line[1];
			angleparm[angleparmnum].name2[0] = line[3];
			angleparm[angleparmnum].name2[1] = line[4];
			angleparm[angleparmnum].name3[0] = line[6];
			angleparm[angleparmnum].name3[1] = line[7];
			strcpy(angleparm[angleparmnum].parm, line);
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
			strcpy(torsionparm[torsionparmnum].parm, line);
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
			improperparm[improperparmnum].name1[0] = line[0];
			improperparm[improperparmnum].name1[1] = line[1];
			improperparm[improperparmnum].name2[0] = line[3];
			improperparm[improperparmnum].name2[1] = line[4];
			improperparm[improperparmnum].name3[0] = line[6];
			improperparm[improperparmnum].name3[1] = line[7];
			improperparm[improperparmnum].name4[0] = line[9];
			improperparm[improperparmnum].name4[1] = line[10];
			strcpy(improperparm[improperparmnum].parm, line);
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
			sscanf(line, "%s", vdwparm[vdwparmnum].name);
			strcpy(vdwparm[vdwparmnum].parm, line);
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
}


void check(char *filename)
{
	int mindex = -1;
	int bindex = 0;
	int aindex = 0;
	int tindex = 0;
	int iindex = 0;
	int vindex = 0;
	int num = 0;
	FILE *fp;
	char line[MAXCHAR];
	char tmpchar[5];

	if ((fp = fopen(filename, "r")) == NULL) {
		printf("\n Cannot open file %s to read in check(), exit", filename);
		return;
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
			} else
				vindex = 3;
		if (mindex == 1) {
			tmpchar[0] = '\0';
			sscanf(line, "%s", tmpchar);
			for (i = 0; i < atomtypenum; i++)
				if (strcmp(tmpchar, atomtype[i].name) == 0)
					atomtypeindex[i] = 1;
		}
		if (bindex == 1)
			for (i = 0; i < bondparmnum; i++)
				if ((bondparm[i].name1[0] == line[0]
					 && bondparm[i].name1[1] == line[1]
					 && bondparm[i].name2[0] == line[3]
					 && bondparm[i].name2[1] == line[4])
					|| (bondparm[i].name2[0] == line[0]
						&& bondparm[i].name2[1] == line[1]
						&& bondparm[i].name1[0] == line[3]
						&& bondparm[i].name1[1] == line[4]))
					bondparmindex[i] = 1;
		if (aindex == 1)
			for (i = 0; i < angleparmnum; i++)
				if ((angleparm[i].name1[0] == line[0]
					 && angleparm[i].name1[1] == line[1]
					 && angleparm[i].name2[0] == line[3]
					 && angleparm[i].name2[1] == line[4]
					 && angleparm[i].name3[0] == line[6]
					 && angleparm[i].name3[1] == line[7])
					|| (angleparm[i].name3[0] == line[0]
						&& angleparm[i].name3[1] == line[1]
						&& angleparm[i].name2[0] == line[3]
						&& angleparm[i].name2[1] == line[4]
						&& angleparm[i].name1[0] == line[6]
						&& angleparm[i].name1[1] == line[7]))
					angleparmindex[i] = 1;
		if (tindex == 1)
			for (i = 0; i < torsionparmnum; i++)
				if ((torsionparm[i].name1[0] == line[0]
					 && torsionparm[i].name1[1] == line[1]
					 && torsionparm[i].name2[0] == line[3]
					 && torsionparm[i].name2[1] == line[4]
					 && torsionparm[i].name3[0] == line[6]
					 && torsionparm[i].name3[1] == line[7]
					 && torsionparm[i].name4[0] == line[9]
					 && torsionparm[i].name4[1] == line[10])
					|| (torsionparm[i].name4[0] == line[0]
						&& torsionparm[i].name4[1] == line[1]
						&& torsionparm[i].name3[0] == line[3]
						&& torsionparm[i].name3[1] == line[4]
						&& torsionparm[i].name2[0] == line[6]
						&& torsionparm[i].name2[1] == line[7]
						&& torsionparm[i].name1[0] == line[9]
						&& torsionparm[i].name1[1] == line[10]))
					torsionparmindex[i] = 1;
		if (iindex == 1)
			if ((improperparm[i].name1[0] == line[0]
				 && improperparm[i].name1[1] == line[1]
				 && improperparm[i].name2[0] == line[3]
				 && improperparm[i].name2[1] == line[4]
				 && improperparm[i].name3[0] == line[6]
				 && improperparm[i].name3[1] == line[7]
				 && improperparm[i].name4[0] == line[9]
				 && improperparm[i].name4[1] == line[10])
				|| (improperparm[i].name4[0] == line[0]
					&& improperparm[i].name4[1] == line[1]
					&& improperparm[i].name3[0] == line[3]
					&& improperparm[i].name3[1] == line[4]
					&& improperparm[i].name2[0] == line[6]
					&& improperparm[i].name2[1] == line[7]
					&& improperparm[i].name1[0] == line[9]
					&& improperparm[i].name1[1] == line[10]))
				improperparmindex[i] = 1;

		if (vindex == 3) {
			printf
				("\nThis program does not support vdw equivalent describations");
			printf
				("\nPlease list the equivalented vdw parameters explicitly");
			exit(1);
		}
		if (vindex == 2) {
			tmpchar[0] = '\0';
			sscanf(line, "%s", tmpchar);
			if (strcmp(vdwparm[i].name, tmpchar) == 0)
				vdwparmindex[i] = 1;
		}
	}
	fclose(fp);
}


void join(char *filename1, char *filename2)
{
	FILE *fp1, *fp2;
	char line[MAXCHAR];
	int mindex = 0;
	int bindex = 0;
	int aindex = 0;
	int tindex = 0;
	int iindex = 0;
	int vindex = 0;

	if ((fp1 = fopen(filename1, "r")) == NULL) {
		printf("\n Cannot open file %s to read in join(), exit", filename1);
		return;
	}
	if ((fp2 = fopen(filename2, "w")) == NULL) {
		printf("\n Cannot open file %s to write in join(), exit", filename2);
		return;
	}

	for (;;) {
		if (fgets(line, MAXCHAR, fp1) == NULL)
			break;
		if (vindex == 1 && strncmp(line, "MOD", 3) == 0)
			vindex = 2;
		if (spaceline(line) == 1) {
			if (mindex == 0)
				mindex = 1;
			if (mindex == 2 && bindex == 0)
				bindex = 1;
			if (bindex == 2 && aindex == 0)
				aindex = 1;
			if (aindex == 2 && tindex == 0)
				tindex = 1;
			if (tindex == 2 && iindex == 0)
				iindex = 1;
			if (iindex == 2 && vindex == 0)
				vindex = 1;
			if (vindex == 2)
				vindex = 3;
		}
		if (vindex == 2 && strncmp(line, "END", 3) == 0)
			vindex = 3;
		if (mindex == 1) {
			for (i = 0; i < atomtypenum; i++)
				if (atomtypeindex[i] == 0)
					fprintf(fp2, "%s", atomtype[i].parm);
			fprintf(fp2, "\n");
			mindex = 2;
			continue;
		}
		if (bindex == 1) {
			for (i = 0; i < bondparmnum; i++)
				if (bondparmindex[i] == 0)
					fprintf(fp2, "%s", bondparm[i].parm);
			fprintf(fp2, "\n");
			bindex = 2;
			continue;
		}

		if (aindex == 1) {
			for (i = 0; i < angleparmnum; i++)
				if (angleparmindex[i] == 0)
					fprintf(fp2, "%s", angleparm[i].parm);
			fprintf(fp2, "\n");
			aindex = 2;
			continue;
		}

		if (tindex == 1) {
			for (i = 0; i < torsionparmnum; i++)
				if (torsionparmindex[i] == 0)
					fprintf(fp2, "%s", torsionparm[i].parm);
			fprintf(fp2, "\n");
			tindex = 2;
			continue;
		}
		if (iindex == 1) {
			for (i = 0; i < improperparmnum; i++)
				if (improperparmindex[i] == 0)
					fprintf(fp2, "%s", improperparm[i].parm);
			fprintf(fp2, "\n");
			iindex = 2;
			continue;
		}
		if (vindex == 3) {
			for (i = 0; i < vdwparmnum; i++)
				if (vdwparmindex[i] == 0)
					fprintf(fp2, "%s", vdwparm[i].parm);
			fprintf(fp2, "\n");
			vindex = -1;
			continue;
		}
		fprintf(fp2, "%s", line);
	}
	fclose(fp1);
	fclose(fp2);
}
int main(int argc, char *argv[])
{
	int i;

	if (strcmp(COLORTEXT, "YES") == 0 || strcmp(COLORTEXT, "yes") == 0) {
		if (argc == 2
			&& (strcmp(argv[1], "-h") == 0
				|| strcmp(argv[1], "-H") == 0)) {
			printf("[31mUsage: parmjoin -p[0m parm file -input \n"
				   "[31m                -m[0m frcmod file \n"
				   "[31m                -o[0m parm file - output\n");
			exit(1);
		}
		if (argc != 7 && argc != 5) {
			printf("[31mUsage: parmjoin -p[0m parm file -input \n"
				   "[31m                -m[0m frcmod file \n"
				   "[31m                -o[0m parm file - output\n");
			exit(1);
		}
	}

	else {
		if (argc == 2
			&& (strcmp(argv[1], "-h") == 0
				|| strcmp(argv[1], "-H") == 0)) {
			printf("Usage: parmjoin -p   parm file - input \n");
			printf("                -m   frcmod file\n");
			printf("                -o   parm file - output\n");
			exit(1);
		}
		if (argc != 7) {
			printf("Usage: parmjoin -p   parm file - input \n");
			printf("                -m   frcmod file\n");
			printf("                -o   parm file - output\n");
			exit(1);
		}
	}
	for (i = 1; i < argc; i += 2) {
		if (strcmp(argv[i], "-p") == 0)
			strcpy(pfilename, argv[i + 1]);
		if (strcmp(argv[i], "-m") == 0)
			strcpy(mfilename, argv[i + 1]);
		if (strcmp(argv[i], "-o") == 0)
			strcpy(ofilename, argv[i + 1]);
	}

/*	allocate memory using calloc*/
	maxatomtype = MAX_FF_ATOMTYPE;
	maxvdwparm = MAX_FF_VDW;
	maxbondparm = MAX_FF_BOND;
	maxangleparm = MAX_FF_ANGLE;
	maxtorsionparm = MAX_FF_TORSION;
	maximproperparm = MAX_FF_IMPROPER;

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

	readfrcmod(mfilename);		/*frcmod */

	atomtypeindex = (int *) calloc(atomtypenum, sizeof(int));
	if (atomtypeindex == NULL) {
		fprintf(stdout, "memory allocation error for *atomtypeindex\n");
		exit(1);
	}
	bondparmindex = (int *) calloc(bondparmnum, sizeof(int));
	if (bondparmindex == NULL) {
		fprintf(stdout, "memory allocation error for *bondparmindex\n");
		exit(1);
	}

	angleparmindex = (int *) calloc(angleparmnum, sizeof(int));
	if (angleparmindex == NULL) {
		fprintf(stdout, "memory allocation error for *angleparmindex\n");
		exit(1);
	}
	torsionparmindex = (int *) calloc(torsionparmnum, sizeof(int));
	if (torsionparmindex == NULL) {
		fprintf(stdout, "memory allocation error for *torsionparmindex\n");
		exit(1);
	}
	improperparmindex = (int *) calloc(improperparmnum, sizeof(int));
	if (improperparmindex == NULL) {
		fprintf(stdout,
				"memory allocation error for *improperparmindex\n");
		exit(1);
	}

	vdwparmindex = (int *) calloc(vdwparmnum, sizeof(int));
	if (vdwparmindex == NULL) {
		fprintf(stdout, "memory allocation error for *vdwparmindex\n");
		exit(1);
	}

	for (i = 0; i < atomtypenum; i++)
		atomtypeindex[i] = 0;
	for (i = 0; i < vdwparmnum; i++)
		vdwparmindex[i] = 0;
	for (i = 0; i < torsionparmnum; i++)
		torsionparmindex[i] = 0;
	for (i = 0; i < improperparmnum; i++)
		improperparmindex[i] = 0;
	for (i = 0; i < bondparmnum; i++)
		bondparmindex[i] = 0;
	for (i = 0; i < angleparmnum; i++)
		angleparmindex[i] = 0;

	check(pfilename);			/*check if the parameters from mfilename already exist in pfilename */
	join(pfilename, ofilename);	/*append the unduplicated parameters in mfilename to pfilename */
/*
	free(atomtypeindex);
	free(bondparmindex);
	free(angleparmindex);
	free(torsionparmindex);
	free(improperparmindex);
	free(vdwparmindex);
	free(atomtype);
	free(bondparm);
	free(angleparm);
	free(torsionparm);
	free(improperparm);
	free(vdwparm);
*/
	return (0);
}
