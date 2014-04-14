/*
************************************************************************
*           All Copyright Reserved!                                    *
*                                                                      *
*  Prog:    prepgen                                                    *
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
# include "ring.c"
# include "ac.c"

char ifilename[MAXCHAR];
char ofilename[MAXCHAR];
char inf_filename[MAXCHAR] = "PREP.INF";
char mc_filename[MAXCHAR];
char pdbfilename[MAXCHAR] = "NEWPDB.PDB";
char mcfilename[30] = "MAINCHAIN.DAT";
FILE *fpin;
FILE *fpout;
FILE *fprep;
FILE *fptest;
FILE *fppdb;

char line[MAXCHAR];

int i, j, k, l;
int mainatomnum = 0;
int selectnum = 0;
int connum = 0;
int countatom = 4;
int cartindex = 0; /*0 for internal format and 1 for cartesian format*/
char head[6];
char tail[6];
char preheadtype[6];
char posttailtype[6];
int read_mc_index = 0;
int headno = 0;
int tailno = 0;
int omitnum = 0;
double charge = 9999.;

int atomnum = 0;
int bondnum = 0;
int ringnum = 0;

ATOM *atom;
AROM *arom;
BOND *bond;
RING *ring;
ATOM *newatom;
ATOM *atoma;
NAME *mc_name;
NAME *omit;
int *selectchain;
int *mainatom;
int *selectindex;
int *seq1;
int *seq2;
int *connumbak;
int *omitno;
int *posindex;
MOLINFO minfo;
CONTROLINFO cinfo;

int adjust_atom_seq_order = 1; 
/* sometimes no valid prepc or prepi can be generated without adjusting the sequence order of atoms,
   should verify the topolgy file using xleap */

int jbond (int at1, int at2) {
int i;
        for(i=0;i<bondnum;i++)
                if((bond[i].bondi == at1 && bond[i].bondj == at2) || 
                   (bond[i].bondi == at2 && bond[i].bondj == at1)) {
                        return 1;
                }
        return 0;
}

void readmc(void)
{
	FILE *fpin;
	char line[MAXCHAR];
	int number1 = 0;
	int number2 = 0;
	int i, j;
	int num_ter_atom = 0;
	int num_omitted_atom = 0;
	int head_flag = 0;
	int tail_flag = 0;
	if ((fpin = fopen(mcfilename, "r")) == NULL) {
		printf("\n Cannot open the main chain file: %s, exit", mcfilename);
		return;
	}
	read_mc_index = 0;
	for (;;) {
		if (fgets(line, 120, fpin) == NULL) {
			/*   printf("\nFinished reading %s file.", mc_filename); */
			break;
		}
		if (strncmp("PRE_HEAD_TYPE", line, 13) == 0) {
			sscanf(&line[14], "%s", preheadtype);
			printf("\nPRE_HEAD_TYPE is %5s", preheadtype);
		}
		if (strncmp("HEAD_NAME", line, 9) == 0)
			sscanf(&line[10], "%s", head);
	
		if (strncmp("POST_TAIL_TYPE", line, 14) == 0) {
			sscanf(&line[15], "%s", posttailtype);
			printf("\nPOST_TAIL_TYPE is %5s", posttailtype);
		}
		if (strncmp("TAIL_NAME", line, 9) == 0)
			sscanf(&line[10], "%s", tail);

		if (strncmp("MAIN_CHAIN", line, 10) == 0) {
			sscanf(&line[11], "%s", mc_name[number1].name);
			number1++;
		}
		if (strncmp("OMIT_NAME", line, 9) == 0) {
			sscanf(&line[10], "%s", omit[number2].name);
			number2++;
		}
		if (strncmp("CHARGE", line, 6) == 0) {
			sscanf(&line[7], "%lf", &charge);
			printf("\nNet charge of truncated molecule is %8.2lf", charge);
		}
	}
	for (i = 0; i < atomnum; i++) {
		if (strcmp(head, atom[i].name) == 0) {
			mainatom[0] = i;
			headno = i;
			head_flag = 1;
			num_ter_atom ++;
			printf("\nHEAD_ATOM  %5d%5s", i + 1, atom[i].name);
		}
		if (strcmp(tail, atom[i].name) == 0) {
			tailno = i;
			tail_flag = 1;
			num_ter_atom ++;
			printf("\nTAIL_ATOM  %5d%5s", i + 1, atom[i].name);
		}
	}
	for (j = 0; j < number1; j++)
		for (i = 0; i < atomnum; i++) {
			if (strcmp(mc_name[j].name, atom[i].name) == 0) {
				if(head_flag == 1) 
					mainatom[j + 1] = i;
				else 
					mainatom[j] = i;
			}
		}
	if(tail_flag == 1) 
		mainatom[number1 + num_ter_atom -1] = tailno;
		
	if (number1 > 0)
		read_mc_index = 2;
	else
		read_mc_index = 1;
	for (i = 0; i < number1 + num_ter_atom; i++)
		printf("\nMAIN_CHAIN %5d%5d%5s", i + 1, mainatom[i] + 1,
			   atom[mainatom[i]].name);

	for (j = 0; j < number2; j++)
		for (i = 0; i < atomnum; i++)
			if (strcmp(omit[j].name, atom[i].name) == 0) {
				omitno[j] = i;
				num_omitted_atom ++;
				printf("\nOMIT_ATOM  %5d%5d%5s", num_omitted_atom, omitno[j] + 1, atom[i].name);
				break;
			}
	mainatomnum = number1 + num_ter_atom;
	omitnum = number2;
	fclose(fpin);
	printf
		("\nNumber of mainchain atoms (including head and tail atom): %5d",
		 mainatomnum);
	printf("\nNumber of omited atoms: %5d", omitnum);
}

void adjustid(void)
{
	int i, j;
	int num;
	num = 0;
	for (i = 0; i < atomnum; i++)
		if (strcmp(head, atom[i].name) == 0) {
			headno = i;
			mainatom[num] = headno;
			num++;
		}
	for (j = 0; j < mainatomnum; j++)
		for (i = 0; i < atomnum; i++)
			if (strcmp(mc_name[j].name, atom[i].name) == 0) {
				mainatom[num] = i;
				num++;
				break;
			}
	for (i = 0; i < atomnum; i++)
		if (strcmp(tail, atom[i].name) == 0) {
			tailno = i;
			mainatom[num] = tailno;
			num++;
		}
}


void improper(void)
{
	int i, j;
	int index;
	int tmpint;
	for (i = 0; i < atomnum; i++) {

/*for sp2 carbon */
		if (atom[i].atomicnum == 6 && atom[i].connum == 3)
			atom[i].improper = 1;
/*for planar nitrogen */
		if (atom[i].atomicnum == 7 && atom[i].connum == 3
			&& (arom[i].ar1 >= 1 || arom[i].ar2 >= 1 || arom[i].ar3 >= 1))
			atom[i].improper = 1;

/* amide nitrogen */
		if ((atom[i].ambername[0] == 'c' && atom[i].ambername[1] == ' '   ) ||
		    (atom[i].ambername[0] == 'c' && strlen(atom[i].ambername) == 1) || 
		    (atom[i].ambername[0] == 'C' && atom[i].ambername[1] == ' '  )  ||
		    (atom[i].ambername[0] == 'C' && strlen(atom[i].ambername) == 1)) {
			index = 0;
			for (j = 0; j < 3; j++) {
				tmpint = atom[i].con[j];
				if (tmpint == -1)
					break;
				if (atom[tmpint].atomicnum == 7) {
					index = 1;
					break;
				}
			}
			if (index == 1) {
				atom[i].improper = 1;
				atom[tmpint].improper = 1;
			}
		}
/*special treatment for those prep files saved from leap*/
		if (atom[i].atomicnum == 6 && atom[i].connum == 2 && atom[i].ambername[0] == 'C'
			&& (atom[i].ambername[1] == ' '
				|| atom[i].ambername[1] == '\0'))
			atom[i].improper = 1;

		if (atom[i].atomicnum == 7 && atom[i].connum == 2 && atom[i].ambername[0] == 'N'
			&& (atom[i].ambername[1] == ' '
				|| atom[i].ambername[1] == '\0'))
			atom[i].improper = 1;

/* check for NH2 group attached to aromatic ring:  */
		if (atom[i].atomicnum == 7 && arom[i].ar1 == 0 && arom[i].ar2 == 0
			&& arom[i].ar3 == 0 &&atom[i].connum == 3) {
			index = 0;
			for (j = 0; j < 3; j++) {
				tmpint = atom[i].con[j];
				if (tmpint == -1)
					break;
				if (arom[tmpint].ar1 >= 1 || arom[tmpint].ar2 >= 1
					|| arom[tmpint].ar3 >= 1) {
					index = 1;
					break;
				}
			}
			if (index == 1)
				atom[i].improper = 1;
		}
/*
		fprintf( stdout, "%5d%5d%5d%5d%5d %s\n",
			i,atom[i].atomicnum,atom[i].connum,atom[i].arom,
			atom[i].improper, atom[i].ambername );
*/
	}
}


void postprep()
{
	int i, j, k;
	int initindex;
	int tmpint1, tmpint2, tmpint3;
	int bond, angle, twist;
	int flag;
	double bondv, anglev, twistv;
	double tmpx, tmpy, tmpz;
	char tmpchar[MAXCHAR];
	char array[3][MAXCHAR];
	char array0[3][MAXCHAR];
	int *connum;
	NAME *aatype;


	aatype = (NAME *) malloc(sizeof(NAME) * (atomnum + 10));
	if (aatype == NULL) {
		fprintf(stdout, "memory allocation error for *aatype\n");
		exit(1);
	}
	connum = (int *) malloc(sizeof(int) * (atomnum + 10));
	if (connum == NULL) {
		fprintf(stdout, "memory allocation error for *connum\n");
		exit(1);
	}

	if ((fprep = fopen(ofilename, "w")) == NULL) {
		printf("\n Cannot open a file %s to write in postprep(), exit", ofilename);
		return;
	}

	for (i = 0; i < atomnum + 3; i++) {
		connum[i] = -1;
		aatype[i].name[0] = 'X';
	}
	for (i = 0; i < atomnum + 3; i++) 
		connumbak[i] = newatom[i].connum;
	for (i = 3; i < atomnum + 3; i++) {
		for (j = 0; j < 6; j++)
			if (newatom[i].con[j] != -1)
				newatom[i].con[j] = seq2[newatom[i].con[j]];
	}
	for (i = 3; i < atomnum + 3; i++) {
		connum[i] = 0;
		for (j = 3; j < atomnum + 3; j++) {
			if (seq2[newatom[i].bondatom - 3] == j)
				connum[i]++;
			if (seq2[newatom[j].bondatom - 3] == i)
				connum[i]++;
		}
	}
	for (i = 3; i < atomnum + 3; i++)
		for (j = i + 1; j < atomnum + 3; j++) {
			if (seq2[newatom[j].bondatom - 3] == i)
				continue;
			if (connum[i] < connumbak[i] && connum[j] < connumbak[j])
				for (k = 0; k < 6; k++)
					if (newatom[i].con[k] == j) {
						newatom[i].connum--;
						newatom[j].connum--;
					}
		}
	for (j = 3; j < atomnum + 3; j++) {
		if (newatom[j].connum == 2)
			aatype[j].name[0] = 'S';
		if (newatom[j].connum == 1)
			aatype[j].name[0] = 'E';
		if (newatom[j].connum == 3)
			aatype[j].name[0] = 'B';
		if (newatom[j].connum == 4)
			aatype[j].name[0] = '3';
		if (newatom[j].connum == 5)
			aatype[j].name[0] = '4';
		if (newatom[j].connum == 6)
			aatype[j].name[0] = '5';
		if (newatom[j].mainatom == 1)
			aatype[j].name[0] = 'M';
	}

	fprintf(fprep, "%5d%5d%5d\n\n", 0, 0, 2);
	fprintf(fprep, "%s\n", "This is a remark line");
	fprintf(fprep, "%s\n", minfo.resfilename);
	if (cartindex == 1) {
		fprintf(fprep, "%s   XYZ  0\n", minfo.resname);
		fprintf(fprep, "%s\n", "CHANGE     OMIT DU   BEG");
	} else {
		fprintf(fprep, "%s   INT  0\n", minfo.resname);
		fprintf(fprep, "%s\n", "CORRECT     OMIT DU   BEG");
	}
	fprintf(fprep, "%s\n", "  0.0000");
	if (cartindex == 0) {
		fprintf(fprep, "%s\n",
				"   1  DUMM  DU    M    0  -1  -2     0.000      .0        .0      .00000");
		fprintf(fprep, "%s\n",
				"   2  DUMM  DU    M    1   0  -1     1.449      .0        .0      .00000");
		fprintf(fprep, "%s\n",
				"   3  DUMM  DU    M    2   1   0     1.523   111.21       .0      .00000");
	}
	if (cartindex == 1) {
		fprintf(fprep, "%s\n",
				"   1  DUMM  DU    M        999.000     999.0      -999.0           .00000");
		fprintf(fprep, "%s\n",
				"   2  DUMM  DU    M        999.000    -999.0       999.0           .00000");
		fprintf(fprep, "%s\n",
				"   3  DUMM  DU    M       -999.000     999.0       999.0           .00000");
	}
	tmpx = newatom[3].x;
	tmpy = newatom[3].y;
	tmpz = newatom[3].z;

	for (i = 3; i < atomnum + 3; i++) {
		newatom[i].x = newatom[i].x + 3.54 - tmpx;
		newatom[i].y = newatom[i].y + 1.42 - tmpy;
/*    newatom[i].z=newatom[i].z+0.0000001-tmpz; */
		newatom[i].z = newatom[i].z - tmpz;
	}
	initindex = 0;
	for (j = 3; j < atomnum + 3; j++) {
		if (initindex == 0) {
			bond = 2;
			angle = 1;
			twist = 0;
			aatype[j].name[0] = 'M';
			initindex = 1;
		} else {
			bond = newatom[j].bondatom - 3;
			bond = seq2[bond];
			angle = newatom[bond].bondatom;
			if (angle >= 3) {
				angle = angle - 3;
				angle = seq2[angle];
			}
			twist = newatom[angle].bondatom;
			if (twist >= 3) {
				twist = twist - 3;
				twist = seq2[twist];
			}
		}
		bondv =
			(newatom[bond].x - newatom[j].x) * (newatom[bond].x -
												newatom[j].x);
		bondv +=
			(newatom[bond].y - newatom[j].y) * (newatom[bond].y -
												newatom[j].y);
		bondv +=
			(newatom[bond].z - newatom[j].z) * (newatom[bond].z -
												newatom[j].z);
		bondv = sqrt(bondv);
		anglev = anglecal(newatom[j], newatom[bond], newatom[angle]);
		twistv =
			rotate(newatom[j], newatom[bond], newatom[angle],
				   &newatom[twist]);
		/* new code for building residue */
		if (cartindex == 0) {
			fprintf(fprep, "%4d  %-5s %-5s %c", j + 1, newatom[j].name,
					newatom[j].ambername, aatype[j].name[0]);
			fprintf(fprep, "%5d%4d%4d%10.3f%10.3f%10.3f%10.6f\n", bond + 1,
					angle + 1, twist + 1, bondv, anglev, twistv,
					newatom[j].charge);
		}
		if (cartindex == 1) {
			fprintf(fprep, "%4d  %-5s %-5s %c", j + 1, newatom[j].name,
					newatom[j].ambername, aatype[j].name[0]);
			fprintf(fprep, "        %10.6f%12.6f%12.6f%12.6f\n",
					atom[seq1[j]].x, atom[seq1[j]].y, atom[seq1[j]].z,
					newatom[j].charge);
		}
	}

	fprintf(fprep, "\n\n%s\n", "LOOP");
	for (i = 3; i < atomnum + 3; i++)
		for (j = i + 1; j < atomnum + 3; j++) {
			if (seq2[newatom[j].bondatom - 3] == i)
				continue;
			if (connum[i] < connumbak[i] && connum[j] < connumbak[j])
				for (k = 0; k < 6; k++)
					if (newatom[i].con[k] == j)
						fprintf(fprep, "%5s%5s\n", newatom[j].name,
								newatom[i].name);
		}
	fprintf(fprep, "\n%s\n", "IMPROPER");
	for (i = 3; i < atomnum + 3; i++)
		if (newatom[i].improper != 0) {
			tmpint1 = newatom[i].con[0];
			tmpint2 = newatom[i].con[1];
			tmpint3 = newatom[i].con[2];

			flag = 1;

			if (tmpint1 != -1) {
				strcpy(array[0], newatom[tmpint1].ambername);
				strcpy(array0[0], newatom[tmpint1].name);
				/* newitoa(tmpint1, tmpchar); */
				sprintf(tmpchar, "%d", tmpint1);
				strcat(array[0], tmpchar);
			} else if (strcmp(head, newatom[i].name) == 0) {
				strcpy(array[0], preheadtype);
				strcpy(array0[0], "-M");
				/* newitoa(tmpint1, tmpchar); */
				sprintf(tmpchar, "%d", tmpint1);
				strcat(array[0], tmpchar);
			} else if (strcmp(tail, newatom[i].name) == 0) {
				strcpy(array[0], posttailtype);
				strcpy(array0[0], "+M");
				/* newitoa(tmpint1, tmpchar); */
				sprintf(tmpchar, "%d", tmpint1);
				strcat(array[0], tmpchar);
			} else if (newatom[i].atomicnum == 6 && newatom[i].connum == 2 && newatom[i].ambername[0] == 'C' && (newatom[i].ambername[1] == ' '
                                || newatom[i].ambername[1] == '\0')) {
				strcpy(array0[0], "+M");
				sprintf(tmpchar, "%d", tmpint1);
				strcat(array[0], tmpchar);
			} else if (newatom[i].atomicnum == 7 && newatom[i].connum == 2 && newatom[i].ambername[0] == 'N' && (newatom[i].ambername[1] == ' '
                                || newatom[i].ambername[1] == '\0')) {
				strcpy(array0[0], "-M");
				sprintf(tmpchar, "%d", tmpint1);
				strcat(array[0], tmpchar);
			} else {
				strcpy(array0[0], "M");
				printf("\nWarning: ATOM %s has unfilled valence, assuming the linked atom name (in other residue) is \"M\"", newatom[i].name);
				printf("\n         change \"M\" to \"-M\" if ATOM %s is linked to the immediate previous residue", newatom[i].name);
				printf("\n         change \"M\" to \"+M\" if ATOM %s is linked to the immediate after residue", newatom[i].name);
				sprintf(tmpchar, "%d", tmpint1);
				strcat(array[0], tmpchar);
				flag = 0;
			}
			if (tmpint2 != -1) {
				strcpy(array[1], newatom[tmpint2].ambername);
				strcpy(array0[1], newatom[tmpint2].name);
				sprintf(tmpchar, "%d", tmpint2);
				/* newitoa(tmpint2, tmpchar); */
				strcat(array[1], tmpchar);
			} else if (strcmp(head, newatom[i].name) == 0) {
				strcpy(array[1], preheadtype);
				strcpy(array0[1], "-M");
				/* newitoa(tmpint2, tmpchar); */
				sprintf(tmpchar, "%d", tmpint2);
				strcat(array[1], tmpchar);
			} else if (strcmp(tail, newatom[i].name) == 0) {
				strcpy(array[1], posttailtype);
				strcpy(array0[1], "+M");
				/* newitoa(tmpint2, tmpchar); */
				sprintf(tmpchar, "%d", tmpint2);
				strcat(array[1], tmpchar);
			} else if (newatom[i].atomicnum == 6 && newatom[i].connum == 2 && newatom[i].ambername[0] == 'C' && (newatom[i].ambername[1] == ' '
                                || newatom[i].ambername[1] == '\0')) {
				strcpy(array0[1], "+M");
				sprintf(tmpchar, "%d", tmpint1);
				strcat(array[1], tmpchar);
			} else if (newatom[i].atomicnum == 7 && newatom[i].connum == 2 && newatom[i].ambername[0] == 'N' && (newatom[i].ambername[1] == ' '
                                || newatom[i].ambername[1] == '\0')) {
				strcpy(array0[1], "-M");
				sprintf(tmpchar, "%d", tmpint1);
				strcat(array[1], tmpchar);
			} else {
				strcpy(array0[1], "M");
				printf("\nWarning: ATOM %s has unfilled valence, assuming the linked atom name (in other residue) is \"M\"", newatom[i].name);
				printf("\n         change \"M\" to \"-M\" if ATOM %s is linked to the immediate previous residue", newatom[i].name);
				printf("\n         change \"M\" to \"+M\" if ATOM %s is linked to the immediate after residue", newatom[i].name);
				sprintf(tmpchar, "%d", tmpint2);
				strcat(array[1], tmpchar);
				flag = 0;
			}

			if (tmpint3 != -1) {
				strcpy(array[2], newatom[tmpint3].ambername);
				strcpy(array0[2], newatom[tmpint3].name);
				/* newitoa(tmpint3, tmpchar); */
				sprintf(tmpchar, "%d", tmpint3);
				strcat(array[2], tmpchar);
			} else if (strcmp(head, newatom[i].name) == 0) {
				strcpy(array[2], preheadtype);
				strcpy(array0[2], "-M");
				/* newitoa(tmpint3, tmpchar); */
				sprintf(tmpchar, "%d", tmpint3);
				strcat(array[2], tmpchar);
			} else if (strcmp(tail, newatom[i].name) == 0) {
				strcpy(array[2], posttailtype);
				strcpy(array0[2], "+M");
				/* newitoa(tmpint3, tmpchar); */
				sprintf(tmpchar, "%d", tmpint3);
				strcat(array[2], tmpchar);
			} else if (newatom[i].atomicnum == 6 && newatom[i].connum == 2 && newatom[i].ambername[0] == 'C' && (newatom[i].ambername[1] == ' '
                                || newatom[i].ambername[1] == '\0')) {
				strcpy(array0[2], "+M");
				sprintf(tmpchar, "%d", tmpint1);
				strcat(array[2], tmpchar);
			} else if (newatom[i].atomicnum == 7 && newatom[i].connum == 2 && newatom[i].ambername[0] == 'N' && (newatom[i].ambername[1] == ' '
                                || newatom[i].ambername[1] == '\0')) {
				strcpy(array0[2], "-M");
				sprintf(tmpchar, "%d", tmpint1);
				strcat(array[2], tmpchar);
			} else {
				strcpy(array0[2], "M");
				printf("\nWarning: ATOM %s has unfilled valence, assuming the linked atom name (in other residue) is \"M\"", newatom[i].name);
				printf("\n         change \"M\" to \"-M\" if ATOM %s is linked to the immediate previous residue", newatom[i].name);
				printf("\n         change \"M\" to \"+M\" if ATOM %s is linked to the immediate after residue", newatom[i].name);
				sprintf(tmpchar, "%d", tmpint3);
				strcat(array[2], tmpchar);
				flag = 0;
			}

			for (j = 0; j < 3; j++)
				for (k = j + 1; k < 3; k++) {
					if (strcmp(array[k], array[j]) < 0) {
						strcpy(tmpchar, array[j]);
						strcpy(array[j], array[k]);
						strcpy(array[k], tmpchar);
						strcpy(tmpchar, array0[j]);
						strcpy(array0[j], array0[k]);
						strcpy(array0[k], tmpchar);
					}
				}
			if(flag == 1) {
				fprintf(fprep, "%5s%5s%5s%5s\n", array0[0], array0[1],
						newatom[i].name, array0[2]);
			}
		}
	fprintf(fprep, "\n%s\n%s\n", "DONE", "STOP");
	free(connum);
	free(aatype);
	fclose(fprep);
}

void adjust(ATOM atom[], int startmainnum)
{
	int i;
	int selectflag;
	for (i = 0; i < 6; i++) {
		if (atom[startmainnum].con[i] == -1)
			break;
		if (atom[atom[startmainnum].con[i]].bondatom == -1 &&
			atom[atom[startmainnum].con[i]].mainatom == -1) {
			selectflag = atom[startmainnum].con[i];
			atom[selectflag].bondatom = startmainnum + 3;
			newatom[countatom++] = atom[selectflag];
			seq1[countatom - 1] = selectflag;
			/*printf("\n i %5d", i); */
			adjust(atom, selectflag);
		}
	}
}

void getmainatom(ATOM atm[], int selectnum, int startnum)
{
/*determain the main atoms of the molecules*/
	int i, j, k, ii, jj;
	int visitindex = 0;
	int startnum1;
	selectchain[selectnum++] = startnum;
	selectindex[startnum] = selectnum - 1;
	if (cartindex == 1) {
		for (ii = 0; ii < selectnum; ii++)
			for (jj = ii + 1; jj < selectnum; jj++)
				if (posindex[selectchain[ii]] > posindex[selectchain[jj]])
					return;
	}
	if (selectnum > mainatomnum) {
		for (j = 0; j < selectnum; j++)
			mainatom[j] = selectchain[j];
		mainatomnum = selectnum;
	}
	if (selectnum > atomnum)
		return;
	for (i = 0; i < 6; i++) {
		if (atm[startnum].con[i] == -1)
			return;
		if (atm[atm[startnum].con[i]].atomicnum == 0)
			continue;
		if (atm[atm[startnum].con[i]].atomicnum == 1)
			continue;
		startnum1 = atm[startnum].con[i];
		if (selectindex[startnum1] != -1)
			continue;
		for (k = 0; k < selectnum; k++)
			if (startnum1 == selectchain[k]) {
				visitindex = 1;
				break;
			}
		if (visitindex == 1) {
			visitindex = 0;
			selectnum = selectindex[selectchain[k]] + 2;
			for (l = selectnum; l < atomnum; l++)
				selectchain[l] = -1;
			continue;
		}
		getmainatom(atm, selectnum, startnum1);
	}
}

void getmainatomspec(ATOM atm[], int selectnum, int startnum)
{
/* 
find a path linking head and tail atoms
determain the aromatic atom through an iterative way 
*/
	int i, j, k, ii, jj;
	int start;
	selectchain[selectnum++] = startnum;
	selectindex[startnum] = 1;
	if (cartindex == 1) {
		for (ii = 0; ii < selectnum; ii++)
			for (jj = ii + 1; jj < selectnum; jj++)
				if (posindex[selectchain[ii]] > posindex[selectchain[jj]])
					return;
	}
	for (i = 0; i < 6; i++) {
		if (atm[startnum].con[i] == -1)
			return;
		start = atm[startnum].con[i];
		if (atm[start].atomicnum == 1)
			continue;
		for (k = 0; k < selectnum - 2; k++)
			if (start == selectchain[k])
				return;
		/* we have already visited this atom */
		if (start == tailno
			&& (mainatomnum == 0 || mainatomnum > selectnum + 1)) {
			for (j = 0; j < selectnum; j++)
				mainatom[j] = selectchain[j];
			mainatomnum = selectnum + 1;
			mainatom[selectnum] = tailno;
			return;
		}
		getmainatomspec(atm, selectnum, start);
	}
}

void getmainatom2(ATOM atm[]) 
{
/*determain the main atoms of the molecules without changing the sequence order of atoms*/
	int i, j;
	int suc;
	int current_id;
	mainatom[0] = 0; 
	current_id = 0;
/*first atom is always a main atom, otherwise, an atom becomes a main atom if it can satisfy the
following condition:
(1) bonded to the latest main atom
(2) the rest atoms does not bond to other main atoms except the latest one*/
	for (i = 1; i < atomnum; i++) {
		if (atom[i].name[0] == 'H')
			continue;
		if (atom[i].atomicnum == 0)
			continue;
		if (jbond(current_id, i) == 0)
			continue;	
		suc = 1;
		for(j=i+1; j<atomnum;j++) 
			if (jbond(current_id, j) == 1) {
				suc = 0;
				break;	
			}	
		if(suc == 1) {
			mainatom[mainatomnum]= i;
			mainatomnum++;
		}
	}	
}

void prep(void)
{
int i,j,k;
int suc ;
/*The following code determine the main chain of the molecule */
	if(adjust_atom_seq_order == 0) 
		getmainatom2(atom);
	else { 	
	for (i = 0; i < atomnum; i++) {
		if (atom[i].name[0] == 'H')
			continue;
		if (atom[i].atomicnum == 0)
			continue;
		selectnum = 0;
		for (j = 0; j < atomnum; j++) {
			selectchain[j] = -1;
			selectindex[j] = -1;
		}
		if (read_mc_index == 0 || (read_mc_index == 1 && headno == -1 && tailno == -1))
			getmainatom(atom, selectnum, i);
	}

/*	if head is known, but tail is not known  */
        if (read_mc_index == 1 && headno != -1 && tailno == -1) {
                selectnum = 0;
                mainatomnum = 0;
                for (j = 0; j < atomnum; j++) {
                        selectchain[j] = -1;
                        selectindex[j] = -1;
                }
                getmainatom(atom, 0, headno);
        }

/*	if head is not known, but tail is known */
        if (read_mc_index == 1 && headno == -1 && tailno != -1) {
                selectnum = 0;
                mainatomnum = 0;
                for (j = 0; j < atomnum; j++) {
                        selectchain[j] = -1;
                        selectindex[j] = -1;
                }
                getmainatom(atom, 0, tailno);
        }

/*	both head and tail are not known */
	if (read_mc_index == 1 && headno != -1 && tailno != -1) {
		selectnum = 0;
		mainatomnum = 0;
		for (j = 0; j < atomnum; j++) {
			selectchain[j] = -1;
			selectindex[j] = -1;
		}
		getmainatomspec(atom, 0, headno);
	}
	}
	fprintf(fptest, "\n\n%s\n", "  -- The Main Chain Of The Molecule --");
	for (i = 0; i < mainatomnum; i++)
		fprintf(fptest, "\n     %6s%5d%5d", "MAIN ATOM", i + 1,
				mainatom[i] + 1);
/*
   now, the main chain of the molecule has determined
   The following code calculate the inter-coordinate and output the file
*/
	for (i = 0; i < atomnum; i++) {
		atom[i].mainatom = -1;
		for (j = 0; j < mainatomnum; j++)
			if (i == mainatom[j]) {
				atom[i].mainatom = 1;
				break;
			}
		atom[i].bondatom = -1;
		atom[i].angleatom = -1;
		atom[i].twistatom = -1;
 /*printf("\n mainatom %5d", atom[i].mainatom);*/
	}
	strcpy(newatom[0].name, "DU");
	newatom[0].x = 0.0;
	newatom[0].y = 0.0;
	newatom[0].z = 0.0;
	strcpy(newatom[1].name, "DU");
	newatom[1].x = 1.4490;
	newatom[1].y = 0.00;
	newatom[1].z = 0.000001;
	strcpy(newatom[2].name, "DU");
	newatom[2].x = 2.00;
	newatom[2].y = 1.42;
	newatom[2].z =-0.000001;
	atom[mainatom[0]].bondatom = 2;
	atom[mainatom[0]].angleatom = 1;
	atom[mainatom[0]].twistatom = 0;
	newatom[0].bondatom = -1;
	newatom[0].angleatom = -1;
	newatom[0].twistatom = -1;
	newatom[0].mainatom = 0;
	newatom[1].bondatom = 0;
	newatom[1].angleatom = -1;
	newatom[1].twistatom = -1;
	newatom[1].mainatom = 0;
	newatom[2].bondatom = 1;
	newatom[2].angleatom = 0;
	newatom[2].twistatom = -1;
	newatom[2].mainatom = 0;
	for (i = 0; i < 3; i++) {
		newatom[i].connum = 0;
		newatom[i].resno = 1;
		newatom[i].con[0] = -1;
		newatom[i].con[1] = -1;
		newatom[i].con[2] = -1;
		newatom[i].con[3] = -1;
		newatom[i].con[4] = -1;
		newatom[i].con[5] = -1;
		newatom[i].charge = 0.0;
		newatom[i].improper = 0;
		strcpy(newatom[i].aa, minfo.resname);
		strcpy(newatom[i].ambername, "DU");
	}

	newatom[3] = atom[mainatom[0]];
	for (i = 0; i < atomnum + 3; i++) {
		seq1[i] = -1;
		seq2[i] = -1;
	}
	seq2[2] = 2;
	seq2[1] = 1;
	seq2[0] = 0;
	seq1[3] = mainatom[0];

	if(adjust_atom_seq_order == 0) {
		seq1[3] = 0;
		newatom[3] = atom[0];
		newatom[3].bondatom = 2;
		for (i= 1; i < atomnum; i++) {
			suc = 0;
			for (j= 0; j < i; j++) {
				if(suc == 1) break;
				for (k= 0; k < 6; k++) {
					if(atom[j].con[k] < 0) continue;
					if(i == atom[j].con[k]) {	
						atom[i].bondatom = j;
						suc = 1;
						break;
					}
				}
			}
			if(suc == 0) {
				fprintf(stderr, "It is impossible to generate a prepi or prepc file without adjusting the sequence order of atoms!");
				exit(1);
			}
			newatom[i+3] = atom[i];	
			newatom[i+3].bondatom += 3;	
			seq1[i+3] = i;
		}
	}
	else {
		adjust(atom, mainatom[0]);
		for (i = 1; i < mainatomnum; i++) {
			seq1[countatom] = mainatom[i];
			atom[mainatom[i]].bondatom = mainatom[i - 1] + 3;
			newatom[countatom++] = atom[mainatom[i]];
			adjust(atom, mainatom[i]);
		}
	}
	for (i = 0; i < atomnum; i++)
		for (j = 3; j < atomnum + 3; j++)
			if (seq1[j] == i)
				seq2[i] = j;
	for (i = 0; i < atomnum; i++)
		fprintf(fppdb, "%4s%7d  %-4s%-4s%5d%12.3f%8.3f%8.3f%10.6f\n",
				"ATOM", i + 1, newatom[i + 3].name, minfo.resname, 1,
				newatom[i + 3].x, newatom[i + 3].y, newatom[i + 3].z,
				newatom[i + 3].charge);
	postprep();
/* The following code write the prep input file*/
	fprintf(fptest, "\n\n%s\n", "  -- The New Coordinates  --");
	for (i = 0; i < atomnum + 3; i++)
		fprintf(fptest,
				"\n %6d %6s %4d %4d %4d %4d %4d %4d %4d %4d %10.6f %10.6f %10.6f %9.6f %5d",
				i, newatom[i].name, newatom[i].con[0], newatom[i].con[1],
				newatom[i].con[2], newatom[i].con[3], newatom[i].mainatom,
				newatom[i].connum, seq1[i], seq2[i], newatom[i].x,
				newatom[i].y, newatom[i].z, newatom[i].charge, newatom[i].bondatom);
}

int main(int argc, char *argv[])
{
	int i, j, k;
	int index;
	int num;
	int tmpint;
	double pcharge = 0.0;
	double ncharge = 0.0;
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
			printf("[31mUsage: prepgen -i[0m  input file name (ac) \n"
				   "[31m               -o[0m  output file name\n"
				   "[31m               -f[0m  output file format (car or int, default: int)\n"
				   "[31m               -m[0m  mainchain file name\n"
				   "[31m               -rn[0m residue name (default: MOL)\n"
				   "[31m               -rf[0m residue file name (default: molecule.res)\n"
				   "[31m               -s [0m sequence order of atoms changeable? ('yes' or 'no', default is 'yes')\n"
				   "		   -f -m -rn -rf -s are optional\n"
				   "		   -s only works when no mainchain file is read in\n");
			exit(1);
		}
		if (argc != 19 && argc != 17 && argc != 15 && argc != 13 && argc != 11
			&& argc != 9 && argc != 7 && argc != 5) {
			printf("[31mUsage: prepgen -i[0m  input file name(ac) \n"
				   "[31m               -o[0m  output file name \n"
				   "[31m               -f[0m  output file format (car or int, default: int)\n"
				   "[31m               -m[0m  mainchain file name\n"
				   "[31m               -rn[0m residue name (default: MOL)\n"
				   "[31m               -rf[0m residue file name (default: molecule.res)\n"
				   "[31m               -s [0m sequence order of atoms changeable? ('yes' or 'no', default is 'yes')\n"
				   "		   -f -m -rn -rf -s are optional\n"
				   "		   -s only works when no mainchain file is read in\n");
			exit(1);
		}
	} else {
		if (argc == 2)
			if (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "-H") == 0) {
				printf("Usage: prepgen -i  input file name(ac) \n");
				printf("               -o  output file name\n");
				printf("               -f  output file format (car or int, default is int)\n");
				printf("               -m  mainchain file name\n");
				printf("               -rn residue name (default MOL)\n");
				printf("               -rf residue file_name (default molecule.res)\n");
				printf("               -s  sequence order of atoms changeable? ('yes' or 'no', default is 'yes')\n");
				printf("               -f -m -rn -rf -s are optional\n");
				printf("               -s only works when no mainchain file is read in\n");
				exit(1);
			}
		if (argc != 19 && argc != 17 && argc != 15 && argc != 13 && argc != 11
			&& argc != 9 && argc != 7 && argc != 5) {
			printf("Usage: prepgen -i  input file name(ac) \n");
			printf("               -o  output file name\n");
			printf("               -f  output file format (car or int, default is int)\n");
			printf("               -m  mainchain file name\n");
			printf("               -rn residue name (default MOL)\n");
			printf("               -rf residue file_name (default molecule.res)\n");
			printf("               -s  sequence order of atoms changeable? ('yes' or 'no', default is 'yes')\n");
			printf("               -f -m -rn -rf are optional\n");
			printf("               -s only works when no mainchain file is read in\n");
			exit(1);
		}
	}
	index = 0;
	read_mc_index = 0;			/* normal case */

	default_cinfo(&cinfo);
	default_minfo(&minfo);

	for (i = 1; i < argc; i += 2) {
		if (strcmp(argv[i], "-i") == 0)
			strcpy(ifilename, argv[i + 1]);
		if (strcmp(argv[i], "-o") == 0)
			strcpy(ofilename, argv[i + 1]);
		if (strcmp(argv[i], "-f") == 0){
			if (strcmp(argv[i + 1], "car") == 0){
				cartindex = 1;
			} else {
				cartindex = 0;
            		}
        	}
		if (strcmp(argv[i], "-p") == 0)
			strcpy(pdbfilename, argv[i + 1]);
		if (strcmp(argv[i], "-m") == 0) {
			strcpy(mcfilename, argv[i + 1]);
			index = 1;
		}
		if (strcmp(argv[i], "-rn") == 0)
			strcpy(minfo.resname, argv[i + 1]);
		if (strcmp(argv[i], "-rf") == 0)
			strcpy(minfo.resfilename, argv[i + 1]);
		if (strcmp(argv[i], "-inf") == 0)
			strcpy(inf_filename, argv[i + 1]);
		if (strcmp(argv[i], "-s") == 0){
			if(strcmp(argv[i + 1], "NO") == 0 || strcmp(argv[i + 1], "no") == 0 ||
			   strcmp(argv[i + 1], "No") == 0 || strcmp(argv[i + 1], "N") == 0  ||
			   strcmp(argv[i + 1], "n") == 0)
			adjust_atom_seq_order = 0;	
		}

	}
	if ((fppdb = fopen(pdbfilename, "w")) == NULL) {
		printf("\n Cannot open the pdb file %s in the main(), exit", pdbfilename);
		exit(1);
	}
	if ((fptest = fopen(inf_filename, "w")) == NULL) {
		printf("\n Cannot open the information file %s in the main(), exit", inf_filename);
		exit(1);
	}

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
	overflow_flag =
		rac(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo);
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
		overflow_flag =
			rac(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo);
	}

	atomicnum(atomnum, atom);
	adjustatomname(atomnum, atom, 1);
/*	wac("test.ac", atomnum, atom, bondnum, bond, cinfo, minfo); */
/*      detect ring property and aromatic property, check file inf_filename for details */
	overflow_flag = ringdetect(atomnum, atom, bondnum, bond, &ringnum, ring, arom,
			   cinfo.maxatom, cinfo.maxring, minfo.inf_filename, 0);
	if(overflow_flag) {
		cinfo.maxring = ringnum + 10;
                free(ring);
        	ring = (RING *) malloc(sizeof(RING) * cinfo.maxring);
        	if (ring == NULL) {
                	fprintf(stdout, "memory allocation error for *ring\n");
                	exit(1);
        	}
		overflow_flag = ringdetect(atomnum, atom, bondnum, bond, &ringnum, ring, arom,
			   	cinfo.maxatom, cinfo.maxring, minfo.inf_filename, 0);
	}
	improper();
	newatom = (ATOM *) malloc(sizeof(ATOM) * (atomnum + 10));
	if (newatom == NULL) {
		fprintf(stdout, "memory allocation error for *newatom\n");
		exit(1);
	}
	atoma = (ATOM *) malloc(sizeof(ATOM) * (atomnum + 10));
	if (atoma == NULL) {
		fprintf(stdout, "memory allocation error for *atoma\n");
		exit(1);
	}
	mc_name = (NAME *) malloc(sizeof(NAME) * (atomnum + 10));
	if (mc_name == NULL) {
		fprintf(stdout, "memory allocation error for *mc_name\n");
		exit(1);
	}
	omit = (NAME *) malloc(sizeof(NAME) * (atomnum + 10));
	if (omit == NULL) {
		fprintf(stdout, "memory allocation error for *omit\n");
		exit(1);
	}
	selectchain = (int *) malloc(sizeof(int) * (atomnum + 10));
	if (selectchain == NULL) {
		fprintf(stdout, "memory allocation error for *selectchain\n");
		exit(1);
	}
	mainatom = (int *) malloc(sizeof(int) * (atomnum + 10));
	if (mainatom == NULL) {
		fprintf(stdout, "memory allocation error for *mainatom\n");
		exit(1);
	}
	selectindex = (int *) malloc(sizeof(int) * (atomnum + 10));
	if (selectindex == NULL) {
		fprintf(stdout, "memory allocation error for *selectindex\n");
		exit(1);
	}
	seq1 = (int *) malloc(sizeof(int) * (atomnum + 10));
	if (seq1 == NULL) {
		fprintf(stdout, "memory allocation error for *seq1\n");
		exit(1);
	}
	seq2 = (int *) malloc(sizeof(int) * (atomnum + 10));
	if (seq2 == NULL) {
		fprintf(stdout, "memory allocation error for *seq2\n");
		exit(1);
	}
	connumbak = (int *) malloc(sizeof(int) * (atomnum + 10));
	if (connumbak == NULL) {
		fprintf(stdout, "memory allocation error for *connumbak\n");
		exit(1);
	}
	omitno = (int *) malloc(sizeof(int) * (atomnum + 10));
	if (omitno == NULL) {
		fprintf(stdout, "memory allocation error for *omitno\n");
		exit(1);
	}
	posindex = (int *) malloc(sizeof(int) * (atomnum + 10));
	if (posindex == NULL) {
		fprintf(stdout, "memory allocation error for *posindex\n");
		exit(1);
	}

	for (i = 0; i < atomnum; i++)
		omitno[i] = -1;

	for (i = 0; i < atomnum; i++) {
		connum = 0;
		for (j = 0; j < 6; j++)
			if (atom[i].con[j] != -1)
				connum++;
		atom[i].connum = connum;
	}
	if (index == 1) {
		readmc();
		adjust_atom_seq_order = 1;
	}
	else {
		for (i = 0; i < atomnum; i++)
			mainatom[i] = -1;
		mainatomnum = 0;
	}
	for (i = 0; i < atomnum; i++)
		posindex[i] = i + 1;
	if (omitnum > 0)
		for (i = 0; i < atomnum; i++) {
/*	ignor omitted atoms */
			index = 0;
			for (j = 0; j < omitnum; j++)
				if (omitno[j] == i) {
					index = 1;
					break;
				}
			if (index == 1)
				continue;
			for (j = 0; j < 6; j++) {
				tmpint = atom[i].con[j];
				if (atom[i].con[j] < 0)
					continue;
/*	if atom i links to headno and tailno, break the bond*/
				for (k = 0; k < omitnum; k++) {
					if ((i == headno || i == tailno)
						&& omitno[k] == tmpint) {
						atom[i].con[j] = -1;
						continue;
					}
					if (omitno[k] < tmpint && atom[i].con[j] >= 0)
						atom[i].con[j]--;
					if (omitno[k] == tmpint) {
						atom[i].connum --;
						atom[i].con[j] = -1;
						printf("\nInfo: there is a bond linking a non-head and non-tail residue atom (%s) and an omitted atom (%s).", atom[i].name, atom[omitno[k]].name);
						printf("\n      You need to specifically add this bond in leap using the command 'bond <atom1> <atom2> [order]'");
						printf("\n      to link %s to an atom in another residue (similar to disulfide bonds)!\n", atom[i].name);
					}
				}
			}
		}

	for (i = 0; i < atomnum; i++) {
		num = 0;
		for (j = 0; j < 6; j++) {
			if (atom[i].con[j] == -1)
				continue;
			atom[i].con[num] = atom[i].con[j];
			num++;
		}
		atom[i].connum = num;
		for (j = atom[i].connum; j < 6; j++)
			atom[i].con[j] = -1;
	}

	if (omitnum > 0) {
		num = 0;
		for (i = 0; i < atomnum; i++) {
			index = 0;
			for (j = 0; j < omitnum; j++)
				if (omitno[j] == i) {
					index = 1;
					break;
				}
			if (index == 1)
				continue;
			atoma[num] = atom[i];
			num++;
		}
		atomnum = num;
		for (i = 0; i < atomnum; i++)
			atom[i] = atoma[i];
		adjustid();				/*adjust headno, tailno, mainatom etc */
	}

	if (omitnum > 0 && charge <= 9990.) {
		pcharge = 0.0;
		ncharge = 0.0;
		for (i = 0; i < atomnum; i++) {
			if (atom[i].charge > 0.0)
				pcharge += atom[i].charge;
			if (atom[i].charge < 0.0)
				ncharge += atom[i].charge;
		}
		if (pcharge == 0.0)
			printf("\nTotal postive charge is 0.0");
		if (ncharge == 0.0)
			printf("\nTotal negtive charge is 0.0");
		if (pcharge != 0.0 && ncharge != 0.0)
			for (i = 0; i < atomnum; i++) {
				if (atom[i].charge >= 0)
					atom[i].charge += atom[i].charge *
						(charge - pcharge - ncharge) / (pcharge - ncharge);
				else
					atom[i].charge -= atom[i].charge *
						(charge - pcharge - ncharge) / (pcharge - ncharge);
				atoma[i].charge = atom[i].charge;
			}
	}
	prep();
	fclose(fptest);
	fclose(fppdb);
/*
	free(atom);
	free(arom);
	free(bond);
        free(newatom);
        free(atoma);
        free(mc_name);
        free(omit);
        free(selectchain);
        free(mainatom);
        free(selectindex);
        free(seq1);
        free(seq2); 
        free(connumbak);
        free(omitno);
        free(posindex);
*/
	return (0);
}
