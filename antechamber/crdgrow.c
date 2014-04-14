/*
************************************************************************
*           All Copyright Reserved!                                    *
*                                                                      *
*  Prog:    crdgrow                                                    *
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
# include "rotate.c"
# include "common.c"
# include "pdb.c"
# include "ac.c"
# include "prep.c"
# define NEGTOR -999999999.0
ATOM *atom1;
ATOM *atom2;
ATOM *atom3;
BOND *bond_tmp;
int atomnum = 0;
int pdbnum = 0;
int bondnum_tmp = 0;
CONTROLINFO cinfo;
MOLINFO minfo;
int resid = 0;
double dist = 0.0;

char line[MAXCHAR];
char ifilename[MAXCHAR];
char ofilename[MAXCHAR];
char pfilename[MAXCHAR];

FILE *fpin;
FILE *fpout;
FILE *fp;

char connect_file[MAXCHAR] = "CONNECT.TPL";



double tor(double torsion)
{
	if (torsion >= NEGTOR - 1 && torsion <= NEGTOR + 1)
		return NEGTOR;
	if (torsion >= 360.)
		torsion = fmod(torsion, 360);
	if (torsion <= -360.)
		torsion = fmod(torsion, 360);
	if (torsion > 180.)
		torsion = -360.0 + torsion;
	if (torsion <= -180.0)
		torsion = 360.0 + torsion;
	return torsion;
}

int main(int argc, char *argv[])
{
	int i, j;
	int *index;
	int ref1, ref2, ref3, tmpint1, tmpint2, tmpint3, tmpint4, tmpint5;
	int num;
	int format = 0;
	char tmpchar[20];
	double bond, angle, torsion, torsion0;
	double torsion1, torsion2, torsion3, torsion4, torsion5;
	double torsion1a, torsion2a, torsion3a, torsion4a, torsion5a;
	double diff;
	double tmpfloat;

    amberhome = (char *) getenv("AMBERHOME");
    if( amberhome == NULL ){
       fprintf( stdout, "AMBERHOME is not set!\n" );
       exit(1);
    }
	if (strcmp(COLORTEXT, "YES") == 0 || strcmp(COLORTEXT, "yes") == 0) {
		if (argc == 2
			&& (strcmp(argv[1], "-h") == 0
				|| strcmp(argv[1], "-H") == 0)) {
			printf("[31mUsage: crdgrow -i[0m input file name \n"
				   "[31m               -o[0m output file name \n"
				   "[31m               -p[0m prepin file name\n"
				   "[31m               -f[0m prepin file format: prepi (the default) or prepc\n\n");
			exit(1);
		}
		if (argc != 9) {
			printf("[31mUsage: crdgrow -i[0m input file name \n"
				   "[31m               -o[0m output file name \n"
				   "[31m               -p[0m prepin file name\n"
				   "[31m               -f[0m prepin file format: prepi (the default) or prepc\n\n");
			exit(1);
		}
	} else {
		if (argc == 2
			&& (strcmp(argv[1], "-h") == 0
				|| strcmp(argv[1], "-H") == 0)) {
			printf("Usage: crdgrow -i input file name (pdb)\n");
			printf("               -o output file name (pdb)\n");
			printf("               -p prepin file name\n");
			printf
				("               -f prepin file format: prepi (default) or prepc\n\n");
			exit(1);
		}
		if (argc != 9) {
			printf("Usage: crdgrow -i input file name (pdb)\n");
			printf("               -o output file name (pdb)\n");
			printf("               -p prepin file name\n");
			printf
				("               -f prepin file format (prepi or prepc)\n\n");
			exit(1);
		}
	}
	tmpchar[0] = '\0';
	for (i = 1; i < argc; i += 2) {
		if (strcmp(argv[i], "-i") == 0)
			strcpy(ifilename, argv[i + 1]);
		if (strcmp(argv[i], "-o") == 0)
			strcpy(ofilename, argv[i + 1]);
		if (strcmp(argv[i], "-p") == 0)
			strcpy(pfilename, argv[i + 1]);
		if (strcmp(argv[i], "-f") == 0) {
			for (j = 0; j <= strlen(argv[i + 1]); j++)
				tmpchar[j] = toupper(argv[i + 1][j]);
			if (strncmp("PREPC", tmpchar, 5) == 0)
				format = 1;
		}
	}

/* for connect.tpl file */

	build_dat_path(connect_file, "CONNECT.TPL", sizeof connect_file, 0);

/* memory allocation*/
	default_cinfo(&cinfo);
	default_minfo(&minfo);
	atom1 = (ATOM *) malloc(sizeof(ATOM) * cinfo.maxatom);
	if (atom1 == NULL) {
		fprintf(stdout, "memory allocation error for *atom\n");
		exit(1);
	}
	atom2 = (ATOM *) malloc(sizeof(ATOM) * cinfo.maxatom);
	if (atom2 == NULL) {
		fprintf(stdout, "memory allocation error for *atom\n");
		exit(1);
	}
	atom3 = (ATOM *) malloc(sizeof(ATOM) * cinfo.maxatom);
	if (atom3 == NULL) {
		fprintf(stdout, "memory allocation error for *atom\n");
		exit(1);
	}
	bond_tmp = (BOND *) malloc(sizeof(BOND) * cinfo.maxbond);
	if (bond_tmp == NULL) {
		fprintf(stdout, "memory allocation error for *bond_tmp\n");
		exit(1);
	}
	for (i = 0; i < cinfo.maxbond; ++i) {
		bond_tmp[i].jflag = -1; /* bond type has not been assigned */
	}
	rpdb(ifilename, &pdbnum, atom2, cinfo, minfo, 0);
	if (format == 0)
		rprepi(pfilename, &atomnum, atom1, &bondnum_tmp, bond_tmp, &cinfo, &minfo);
	else
		rprepc(pfilename, &atomnum, atom1, &cinfo, &minfo);
	atomicnum(atomnum, atom1);
	if (format != 0)
		connect(connect_file, atomnum, atom1, &bondnum_tmp, bond_tmp,
			cinfo.maxbond);

	index = (int *) malloc(sizeof(int) * (atomnum +10));
	if (index == NULL) {
		fprintf(stdout, "memory allocation error for *index\n");
		exit(1);
	}

	for (i = 0; i < atomnum; i++) {
		atom3[i].charge = 0.0;
		index[i] = 0;
	}
	for (i = 0; i < pdbnum; i++)
		for (j = 0; j < atomnum; j++)
			if (strcmp(atom2[i].name, atom1[j].name) == 0) {
				index[j] = 1;
				atom3[j] = atom2[i];
				strcpy(atom3[j].aa, atom1[j].aa);
			}

	for (i = 0; i < atomnum; i++)
		if (index[i] == 0) {
			ref1 = -1;
			ref2 = -1;
			ref3 = -1;
			num = 0;
			for (j = 0; j < atom1[i].connum; j++)
				if (index[atom1[i].con[j]] == 1) {
					ref1 = atom1[i].con[j];
					num++;
					break;
				}
			for (j = 0; j < atom1[ref1].connum; j++)
				if (index[atom1[ref1].con[j]] == 1) {
					ref2 = atom1[ref1].con[j];
					num++;
					break;
				}
			for (j = 0; j < atom1[ref2].connum; j++)
				if (index[atom1[ref2].con[j]] == 1
					&& atom1[ref2].con[j] != ref1) {
					ref3 = atom1[ref2].con[j];
					num++;
					break;
				}
			if (num != 3)
				for (j = 0; j < atomnum; j++)
					if (index[j] == 1) {
						if (ref1 == -1 && j != ref2 && j != ref3) {
							ref1 = j;
							num++;
							if (num == 3)
								break;
						}
						if (ref2 == -1 && j != ref1 && j != ref3) {
							ref2 = j;
							num++;
							if (num == 3)
								break;
						}
						if (ref3 == -1 && j != ref1 && j != ref2) {
							ref3 = j;
							num++;
							if (num == 3)
								break;
						}
					}
			if (num != 3) {
				printf
					("\nAt least three atoms are need in the refer pdb file, exit");
				exit(1);
			}
			atom3[i] = atom1[i];
			bond = distance(atom1[i], atom1[ref1]);
			angle = anglecal(atom1[i], atom1[ref1], atom1[ref2]);
			torsion =
				rotate(atom1[i], atom1[ref1], atom1[ref2], &atom1[ref3]);
			torsion1 = NEGTOR;
			torsion2 = NEGTOR;
			torsion3 = NEGTOR;
			torsion4 = NEGTOR;
			torsion5 = NEGTOR;
			torsion1a = NEGTOR;
			torsion2a = NEGTOR;
			torsion3a = NEGTOR;
			torsion4a = NEGTOR;
			torsion5a = NEGTOR;
			tmpint1 = atom1[ref1].con[0];
			if (i != tmpint1 && ref2 != tmpint1 && tmpint1 != -1) {
				if (index[tmpint1] == 1)
					torsion1 =
						rotate(atom3[tmpint1], atom3[ref1], atom3[ref2],
							   &atom3[ref3]);
				torsion1a =
					rotate(atom1[tmpint1], atom1[ref1], atom1[ref2],
						   &atom1[ref3]);
			}
			tmpint2 = atom1[ref1].con[1];
			if (i != tmpint2 && ref2 != tmpint2 && tmpint2 != -1) {
				if (index[tmpint2] == 1)
					torsion2 =
						rotate(atom3[tmpint2], atom3[ref1], atom3[ref2],
							   &atom3[ref3]);
				torsion2a =
					rotate(atom1[tmpint2], atom1[ref1], atom1[ref2],
						   &atom1[ref3]);
			}
			tmpint3 = atom1[ref1].con[2];
			if (i != tmpint3 && ref2 != tmpint3 && tmpint3 != -1) {
				if (index[tmpint3] == 1)
					torsion3 =
						rotate(atom3[tmpint3], atom3[ref1], atom3[ref2],
							   &atom3[ref3]);
				torsion3a =
					rotate(atom1[tmpint3], atom1[ref1], atom1[ref2],
						   &atom1[ref3]);
			}
			tmpint4 = atom1[ref1].con[3];
			if (i != tmpint4 && ref2 != tmpint4 && tmpint4 != -1) {
				if (index[tmpint4] == 1)
					torsion4 =
						rotate(atom3[tmpint4], atom3[ref1], atom3[ref2],
							   &atom3[ref3]);
				torsion4a =
					rotate(atom1[tmpint4], atom1[ref1], atom1[ref2],
						   &atom1[ref3]);
			}
			tmpint5 = atom1[ref1].con[4];
			if (i != tmpint5 && ref2 != tmpint5 && tmpint5 != -1) {
				if (index[tmpint5] == 1)
					torsion5 =
						rotate(atom3[tmpint5], atom3[ref1], atom3[ref2],
							   &atom3[ref3]);
				torsion5a =
					rotate(atom1[tmpint5], atom1[ref1], atom1[ref2],
						   &atom1[ref3]);
			}
			torsion1 = tor(torsion1);
			torsion2 = tor(torsion2);
			torsion3 = tor(torsion3);
			torsion4 = tor(torsion4);
			torsion5 = tor(torsion5);
			torsion1a = tor(torsion1a);
			torsion2a = tor(torsion2a);
			torsion3a = tor(torsion3a);
			torsion4a = tor(torsion4a);
			torsion5a = tor(torsion5a);

			num = 0;
			diff = 0;
			if (torsion1 > NEGTOR + 1 && torsion1a > NEGTOR + 1) {
				tmpfloat = torsion1 - torsion1a;
				tmpfloat = tor(tmpfloat);
				diff += tmpfloat;
				num++;
			}
			if (torsion2 > NEGTOR + 1 && torsion2a > NEGTOR + 1) {
				tmpfloat = torsion2 - torsion2a;
				tmpfloat = tor(tmpfloat);
				diff += tmpfloat;
				num++;
			}
			if (torsion3 > NEGTOR + 1 && torsion3a > NEGTOR + 1) {
				tmpfloat = torsion3 - torsion3a;
				tmpfloat = tor(tmpfloat);
				diff += tmpfloat;
				num++;
			}
			if (torsion4 > NEGTOR + 1 && torsion4a > NEGTOR + 1) {
				tmpfloat = torsion4 - torsion4a;
				tmpfloat = tor(tmpfloat);
				diff += tmpfloat;
				num++;
			}
			if (torsion5 > NEGTOR + 1 && torsion5a > NEGTOR + 1) {
				tmpfloat = torsion5 - torsion5a;
				tmpfloat = tor(tmpfloat);
				diff += tmpfloat;
				num++;
			}

			if (num >= 1) {
				diff /= num;
				torsion0 = torsion + diff;
			} else
				torsion0 = torsion;
			rotate0(atom3[ref3], atom3[ref2], atom3[ref1], &atom3[i], bond,
					angle, torsion0);
			index[i] = 1;
		}
	wpdb(ofilename, atomnum, atom3);
/*
	free(atom1);
	free(atom2);
	free(atom3);
	free(bond_tmp);
*/
	return (0);
/* The end */
}
