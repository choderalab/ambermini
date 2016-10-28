# include "rotate.c"
# define MAXPATHATOMNUM 1000
# define MAX_EQU_ATOM 1000
# define MAXGEOMPARM  100
# define MAXGEOM      1000

typedef struct {
        int id1;
        int id2;
        int id3;
        int id4;
        int type; /*0 : torsional angle id1-id2-id3-id4 is 0, 1: the torsional angle is 180 */
} GEOM;

typedef struct {
        char at1[20];
        char at2[20];
} GEOMPARM;
GEOMPARM geomparm[MAXGEOMPARM];
int ngeomparm ;

char *amberhome;

int *input_selectindex;
int *input_pathnum;
int *input_selectelement;
int *input_pathatomnum;
double **input_pathscore;
int *ref_selectindex;
int *ref_pathnum;
int *ref_selectelement;
int *ref_pathatomnum;
double **ref_pathscore;

int selectnum = 0;
int pathnumindex = 0;
int atomindex = 0;
int i;

void read_geom_parm() {
FILE *fp;
char filename[MAXCHAR];
char line[MAXCHAR];

amberhome = (char *) getenv("AMBERHOME");
if( amberhome == NULL ) {
        fprintf( stdout, "AMBERHOME is not set!\n" );
        strcpy(amberhome, "");
}
else {
        strcpy(filename, amberhome);
        strcat(filename, "/dat/antechamber/");
}
strcat(filename, "ATOM_EQU.TYPE");

if ((fp = fopen(filename, "r")) == NULL) {
        printf("\n Cannot open file %s to read in geometry(), exit", filename);
        return;
}
ngeomparm = 0;
for(;;) {
        if (fgets(line, MAXCHAR, fp) == NULL) break;
        if (strncmp("PARM", line, 4) == 0) {
        	sscanf(&line[4], "%s%s", geomparm[ngeomparm].at1, geomparm[ngeomparm].at2);
                ngeomparm++;
		if(ngeomparm >= MAXGEOMPARM) {
			fprintf(stderr,"Too many parameter defined, extending MAXGEOMPARM in equatom.c and recompile the program\n");	
			exit(1);
		}
        }
}
fclose (fp);
/*
for(i=0;i<ngeomparm;i++)
        printf("parm %5d %5s %5s\n", i+1, geomparm[i].at1, geomparm[i].at2);
*/
}

void geometry(ATOM *atom, int atomnum, BOND *bond, int bondnum, GEOM *geom, int *ngeom) {

int i,j,m,n;
int at1, at2, at3, at4;
int suc;
int count;
double tor;
/* find how many E/Z configurations*/
count = 0;
for(i=0;i<bondnum;i++) {
        at2 = bond[i].bondi;
        at3 = bond[i].bondj;
        suc = 0;
        for(j=0;j<ngeomparm;j++) {
                if(strcmp(geomparm[j].at1, atom[at2].ambername) == 0 &&
                   strcmp(geomparm[j].at2, atom[at3].ambername) == 0) {
                        suc = 1;
                        break;
                }
                if(strcmp(geomparm[j].at2, atom[at2].ambername) == 0 &&
                   strcmp(geomparm[j].at1, atom[at3].ambername) == 0) {
                        suc = 1;
                        break;
                }
        }

        if(suc == 1) {
                for(m=0;m<6;m++) {
                        if(atom[at2].con[m] < 0) break;
                        if(atom[at2].con[m] == at3) continue;
                        at1 = atom[at2].con[m];
                        for(n=0;n<6;n++) {
                                if(atom[at3].con[n] < 0) break;
                                if(atom[at3].con[n] == at2) continue;
                                at4 = atom[at3].con[n];
                                if(at1 < at4) {
                                        geom[count].id1 = at1;
                                        geom[count].id2 = at2;
                                        geom[count].id3 = at3;
                                        geom[count].id4 = at4;
                                }
                                else {
                                        geom[count].id4 = at1;
                                        geom[count].id3 = at2;
                                        geom[count].id2 = at3;
                                        geom[count].id1 = at4;
                                }
                                tor = rotate(atom[at1], atom[at2], atom[at3], &atom[at4]);
                                if(tor <= 90 && tor  >= -90)
                                        geom[count].type = 0;
                                else
                                        geom[count].type = 1;
                                count ++;
				if(count >= MAXGEOM) {
					fprintf(stderr,"Too many E/Z configuration found, extending MAXGEOM in equatom.c and recompile the program\n");	
					exit(1);
				}
                        }
                }
        }
}
*ngeom = count;
/*
for(i=0;i<ngeom;i++)
        printf("geom %5d %5s %5s %5s %5s %5d\n", i+1, atom[geom[i].id1].ambername, atom[geom[i].id2].ambername,
                                                      atom[geom[i].id3].ambername, atom[geom[i].id4].ambername, geom[i].type);
*/
}

void sort(double array[], int elemnum)
{
	int i, j;
	double tmp;
	for (i = 0; i < elemnum; i++)
		for (j = i + 1; j < elemnum; j++) {
			if (array[j] < array[i]) {
				tmp = array[i];
				array[i] = array[j];
				array[j] = tmp;
			}
		}
}

void scorepath_input(ATOM atom[], int selectnum, int startnum, int max_path_length, int ngeom, GEOM geom[])
{
	int i, j, k, m;
	int start;
	double score;
	double coef;
	int resetindex;
	start = -1;
	resetindex = -1;
	input_selectindex[startnum] = selectnum;
	input_selectelement[selectnum++] = startnum;
	if(max_path_length != -1 && selectnum > max_path_length)
		return;
	for (i = 0; i < 6; i++) {
		if (atom[startnum].con[i] == -1) {
			score = 0.0;
			for (j = 0; j < selectnum; j++) {
				coef = 1.0;
				if(ngeom > 0) {
					for(m=0; m< ngeom; m++) {
						if(j>=3 && input_selectelement[j - 3] == geom[m].id1 && 
                                                           input_selectelement[j - 2] == geom[m].id2 &&
                                                           input_selectelement[j - 1] == geom[m].id3 &&
                                                           input_selectelement[j]     == geom[m].id4) {
						 	if(geom[m].type == 1) 
								coef = 1.01;		
							else 
								coef = 0.99;
							break;
						}
						if(j>=3 && input_selectelement[j - 3] == geom[m].id4 && 
                                                           input_selectelement[j - 2] == geom[m].id3 &&
                                                           input_selectelement[j - 1] == geom[m].id2 &&
                                                           input_selectelement[j]     == geom[m].id1) {
						 	if(geom[m].type == 1) 
								coef = 1.01;		
							else 
								coef = 0.99;
							break;
						}
					}

					for(m=0; m< ngeom; m++) {
						if((selectnum - j)>3 && input_selectelement[j + 3] == geom[m].id4 && 
                                                                        input_selectelement[j + 2] == geom[m].id3 &&
                                                                        input_selectelement[j + 1] == geom[m].id2 &&
                                                                        input_selectelement[j]     == geom[m].id1) {
						 	if(geom[m].type == 1) 
								coef = 1.01;		
							else 
								coef = 0.99;
							break;
						}
						if((selectnum - j)>3 && input_selectelement[j + 3] == geom[m].id1 && 
                                                                        input_selectelement[j + 2] == geom[m].id2 &&
                                                                        input_selectelement[j + 1] == geom[m].id3 &&
                                                                        input_selectelement[j]     == geom[m].id4) {
						 	if(geom[m].type == 1) 
								coef = 1.01;		
							else 
								coef = 0.99;
							break;
						}
					}
				}
				score += coef*((j + 1) * 0.11 +
					       atom[input_selectelement[j]].atomicnum * 0.08);
			}
			input_pathscore[atomindex][pathnumindex++] = score;
			if (pathnumindex >= input_pathatomnum[atomindex]) {
				input_pathatomnum[atomindex] += MAXPATHATOMNUM;
				fprintf
					(stdout, "\nInfo: the number of the path atoms exceeds MAXPATHATOMNUM(%d) for atom[%d],extend the size and reallocate the memory automatically",
					 input_pathatomnum[atomindex], atomindex);
				input_pathscore[atomindex] =
					(double *) realloc(input_pathscore[atomindex],
									   input_pathatomnum[atomindex] *
									   sizeof(double));
				if (input_pathscore[atomindex] == NULL) {
					fprintf(stdout,
							" reallocate memory for pathscore[%d] failed\n",
							atomindex);
					exit(1);
				}
			}
			return;
		}
		start = atom[startnum].con[i];
		for (k = 0; k < selectnum; k++)
			if (start == input_selectelement[k]) {
				resetindex = 1;
				break;
			}
		if (resetindex == 1) {
			resetindex = -1;
			continue; 
		}
		if (start == -1)
			return;
		scorepath_input(atom, selectnum, start, max_path_length, ngeom, geom);
		/* we have already visited this atom */
	}
}
void scorepath_ref(ATOM atom[], int selectnum, int startnum, int max_path_length, int ngeom, GEOM geom[])
{
	int i, j, k, m;
	int start;
	double score;
	double coef;
	int resetindex;
	start = -1;
	resetindex = -1;
	ref_selectindex[startnum] = selectnum;
	ref_selectelement[selectnum++] = startnum;
	if(max_path_length != -1 && selectnum > max_path_length)
		return;
	for (i = 0; i < 6; i++) {
		if (atom[startnum].con[i] == -1) {
			score = 0.0;
			for (j = 0; j < selectnum; j++) {
				coef = 1.0;
				if(ngeom > 0) {
					for(m=0; m< ngeom; m++) {
						if(j>=3 && ref_selectelement[j - 3] == geom[m].id1 && 
                                                           ref_selectelement[j - 2] == geom[m].id2 &&
                                                           ref_selectelement[j - 1] == geom[m].id3 &&
                                                           ref_selectelement[j]     == geom[m].id4) {
						 	if(geom[m].type == 1) 
								coef = 1.01;		
							else 
								coef = 0.99;
							break;
						}
						if(j>=3 && ref_selectelement[j - 3] == geom[m].id4 && 
                                                           ref_selectelement[j - 2] == geom[m].id3 &&
                                                           ref_selectelement[j - 1] == geom[m].id2 &&
                                                           ref_selectelement[j]     == geom[m].id1) {
						 	if(geom[m].type == 1) 
								coef = 1.01;		
							else 
								coef = 0.99;
							break;
						}
					}

					for(m=0; m< ngeom; m++) {
						if((selectnum - j)>3 && ref_selectelement[j + 3] == geom[m].id4 && 
                                                                        ref_selectelement[j + 2] == geom[m].id3 &&
                                                                        ref_selectelement[j + 1] == geom[m].id2 &&
                                                                        ref_selectelement[j]     == geom[m].id1) {
						 	if(geom[m].type == 1) 
								coef = 1.01;		
							else 
								coef = 0.99;
							break;
						}
						if((selectnum - j)>3 && ref_selectelement[j + 3] == geom[m].id1 && 
                                                                        ref_selectelement[j + 2] == geom[m].id2 &&
                                                                        ref_selectelement[j + 1] == geom[m].id3 &&
                                                                        ref_selectelement[j]     == geom[m].id4) {
						 	if(geom[m].type == 1) 
								coef = 1.01;		
							else 
								coef = 0.99;
							break;
						}
					}
				}
				score += coef*((j + 1) * 0.11 +
					       atom[ref_selectelement[j]].atomicnum * 0.08);
			}
			ref_pathscore[atomindex][pathnumindex++] = score;
			if (pathnumindex >= ref_pathatomnum[atomindex]) {
				ref_pathatomnum[atomindex] += MAXPATHATOMNUM;
				fprintf
					(stdout, "\nInfo: the number of the path atoms exceeds MAXPATHATOMNUM(%d) for atom[%d],extend the size and reallocate the memory automatically",
					 ref_pathatomnum[atomindex], atomindex);
				ref_pathscore[atomindex] =
					(double *) realloc(ref_pathscore[atomindex],
									   ref_pathatomnum[atomindex] *
									   sizeof(double));
				if (ref_pathscore[atomindex] == NULL) {
					fprintf(stdout,
							" reallocate memory for pathscore[%d] failed\n",
							atomindex);
					exit(1);
				}
			}
			return;
		}
		start = atom[startnum].con[i];
		for (k = 0; k < selectnum; k++)
			if (start == ref_selectelement[k]) {
				resetindex = 1;
				break;
			}
		if (resetindex == 1) {
			resetindex = -1;
			continue; 
		}
		if (start == -1)
			return;
		scorepath_ref(atom, selectnum, start, max_path_length, ngeom, geom);
		/* we have already visited this atom */
	}
}

void match_atompair(int input_atomnum, int ref_atomnum, ATOM *input_atom, ATOM *ref_atom, int *input_equ_atom_id, int *ref_equ_atom_id, int *input_atom_id, int *ref_atom_id, int input_bondnum, int ref_bondnum, BOND *input_bond, BOND *ref_bond, int max_path_length, int igeom) {
int i, j, k, m, n;
int  input_ngeom = 0;
int  ref_ngeom = 0;
int equindex;
int matchindex;
int suc_flag = 0;
int ncyc, nfail;
int atid, atid1, atid2;
double sum;
GEOM *input_geom = NULL;
GEOM *ref_geom = NULL;
/* allocation memory */
        if(input_atomnum > MAX_EQU_ATOM) {
                fprintf(stdout, "The number of atoms (%d) exceed the MAX_EQU_ATOM (%d) defined in equatom.c, extend MAX_EQU_ATOM and recompile the program\n",
             input_atomnum, MAX_EQU_ATOM);
                exit(1);
        }
        if(ref_atomnum > MAX_EQU_ATOM) {
                fprintf(stdout, "The number of atoms (%d) exceed the MAX_EQU_ATOM (%d) defined in equatom.c, extend MAX_EQU_ATOM and recompile the program\n",
             ref_atomnum, MAX_EQU_ATOM);
                exit(1);
        }

	input_selectindex = (int *) malloc(sizeof(int) * input_atomnum);
	if (input_selectindex == NULL) {
		fprintf(stdout, "memory allocation error for *input_selectindex\n");
		exit(1);
	}
	ref_selectindex = (int *) malloc(sizeof(int) * ref_atomnum);
	if (ref_selectindex == NULL) {
		fprintf(stdout, "memory allocation error for *ref_selectindex\n");
		exit(1);
	}

	input_pathnum = (int *) malloc(sizeof(int) * input_atomnum);
	if (input_pathnum == NULL) {
		fprintf(stdout, "memory allocation error for *input_pathnum\n");
		exit(1);
	}
	ref_pathnum = (int *) malloc(sizeof(int) * ref_atomnum);
	if (ref_pathnum == NULL) {
		fprintf(stdout, "memory allocation error for *ref_pathnum\n");
		exit(1);
	}

	input_selectelement = (int *) malloc(sizeof(int) * input_atomnum);
	if (input_selectelement == NULL) {
		fprintf(stdout, "memory allocation error for *input_selectelement\n");
		exit(1);
	}
	ref_selectelement = (int *) malloc(sizeof(int) * ref_atomnum);
	if (ref_selectelement == NULL) {
		fprintf(stdout, "memory allocation error for *ref_selectelement\n");
		exit(1);
	}

	input_pathatomnum = (int *) malloc(sizeof(int) * input_atomnum);
	if (input_pathatomnum == NULL) {
		fprintf(stdout, "memory allocation error for *input_pathatomnum\n");
		exit(1);
	}
	input_pathscore = (double**)malloc(MAX_EQU_ATOM *sizeof(double*));
	if (input_pathscore == NULL) {
        	fprintf(stderr, "memory allocation error for **input_pathscore\n");
        	exit(0);
	}
	for (i = 0; i < input_atomnum; i++) {
		input_pathatomnum[i] = MAXPATHATOMNUM;
		input_pathscore[i] = (double *) calloc(input_pathatomnum[i], sizeof(double));
		if (input_pathscore == NULL) {
			fprintf(stdout, "memory allocation error for *input_pathscore[%d]\n",
					i + 1);
			exit(1);
		}
	}
	ref_pathatomnum = (int *) malloc(sizeof(int) * ref_atomnum);
	if (ref_pathatomnum == NULL) {
		fprintf(stdout, "memory allocation error for *ref_pathatomnum\n");
		exit(1);
	}
	ref_pathscore = (double**)malloc(MAX_EQU_ATOM *sizeof(double*));
	if (ref_pathscore == NULL) {
        	fprintf(stderr, "memory allocation error for **ref_pathscore\n");
        	exit(0);
	}
	for (i = 0; i < ref_atomnum; i++) {
		ref_pathatomnum[i] = MAXPATHATOMNUM;
		ref_pathscore[i] = (double *) calloc(ref_pathatomnum[i], sizeof(double));
		if (ref_pathscore == NULL) {
			fprintf(stdout, "memory allocation error for *ref_pathscore[%d]\n",
					i + 1);
			exit(1);
		}
	}

	if(igeom == 1)  {
		read_geom_parm();
		input_ngeom = 0;
		ref_ngeom = 0;
		if(ngeomparm > 0) {
			input_geom = (GEOM *) malloc(sizeof(GEOM) * MAXGEOM);
			if (input_geom == NULL) {
				fprintf(stdout, "memory allocation error for *input_geom\n");
				exit(1);
			}
			geometry(input_atom, input_atomnum, input_bond, input_bondnum, input_geom, &input_ngeom);

			ref_geom = (GEOM *) malloc(sizeof(GEOM) * MAXGEOM);
			if (ref_geom == NULL) {
				fprintf(stdout, "memory allocation error for *ref_geom\n");
				exit(1);
			}
			geometry(ref_atom, ref_atomnum, ref_bond, ref_bondnum, ref_geom, &ref_ngeom);
		}
	}
/* for input */
        for (i = 0; i < input_atomnum; i++)
                input_equ_atom_id[i] = -1;
        for (i = 0; i < input_atomnum; i++) {
                selectnum = 0;
                pathnumindex = 0;
                atomindex = i;
                for (j = 0; j < input_atomnum; j++) {
                        input_selectindex[j] = -1;
                        input_selectelement[i] = -1;
                }
                scorepath_input(input_atom, 0, i, max_path_length, input_ngeom, input_geom);
                input_pathnum[i] = pathnumindex;
        }
        selectnum = 0;
        for (i = 0; i < input_atomnum; i++) {
                sum = 0.0;
                for (j = 0; j < input_pathnum[i]; j++) {
                        sum += input_pathscore[i][j];
                }
        }
        for (i = 0; i < input_atomnum; i++) {
                sort(input_pathscore[i], input_pathnum[i]);
/*		
		printf("INPUT %5d %5s %5d\n", i+1, input_atom[i].name, input_pathnum[i]);
*/
	}
/* 	for equivalent atoms within input */
        for (i = 0; i < input_atomnum; i++) {
                for (j = i + 1; j < input_atomnum; j++) {
                        if (input_equ_atom_id[j] >= 0)
                                continue;
                        equindex = 1;
                        if (input_pathnum[i] != input_pathnum[j])
                                continue;
                        for (k = 0; k < input_pathnum[i]; k++)
                                if (input_pathscore[i][k] != input_pathscore[j][k]) {
                                        equindex = -1;
                                        break;
                                }
                        if (equindex == 1) {
                                input_equ_atom_id[j] = i;
                                input_equ_atom_id[i] = -2;
			}
                }
        }
/* for ref*/
        for (i = 0; i < ref_atomnum; i++)
                ref_equ_atom_id[i] = -1;
        for (i = 0; i < ref_atomnum; i++) {
                selectnum = 0;
                pathnumindex = 0;
                atomindex = i;
                for (j = 0; j < ref_atomnum; j++) {
                        ref_selectindex[j] = -1;
                        ref_selectelement[i] = -1;
                }
                scorepath_ref(ref_atom, 0, i, max_path_length, ref_ngeom, ref_geom);
                ref_pathnum[i] = pathnumindex;
        }
        selectnum = 0;
        for (i = 0; i < ref_atomnum; i++) {
                sum = 0.0;
                for (j = 0; j < ref_pathnum[i]; j++) {
                        sum += ref_pathscore[i][j];
                }
        }
        for (i = 0; i < ref_atomnum; i++)  {
                sort(ref_pathscore[i], ref_pathnum[i]);
/*		
		printf("REF   %5d %5s %5d\n", i+1, ref_atom[i].name, ref_pathnum[i]);
*/
	}
/* 	for equivalent atoms within ref */
        for (i = 0; i < ref_atomnum; i++) {
                for (j = i + 1; j < ref_atomnum; j++) {
                        if (ref_equ_atom_id[j] >= 0)
                                continue;
                        equindex = 1;
                        if (ref_pathnum[i] != ref_pathnum[j])
                                continue;
                        for (k = 0; k < ref_pathnum[i]; k++)
                                if (ref_pathscore[i][k] != ref_pathscore[j][k]) {
                                        equindex = -1;
                                        break;
                                }
                        if (equindex == 1) {
                                ref_equ_atom_id[j] = i;
                                ref_equ_atom_id[i] = -2;
			}
                }
        }
/*
	for(i=0;i<input_atomnum;i++) {
		if(input_equ_atom_id[i] >= 0)
			printf("%5d %5s %5d %5s\n", i+1, input_atom[i].name, input_equ_atom_id[i], input_atom[input_equ_atom_id[i]].name);
		else
			printf("%5d %5s %5d\n", i+1, input_atom[i].name, input_equ_atom_id[i]);
	}
	printf("\n\n");
	for(i=0;i<ref_atomnum;i++) {
		if(ref_equ_atom_id[i] >= 0)
			printf("%5d %5s %5d %5s\n", i+1, ref_atom[i].name, ref_equ_atom_id[i], ref_atom[ref_equ_atom_id[i]].name);
		else
			printf("%5d %5s %5d\n", i+1, ref_atom[i].name, ref_equ_atom_id[i]);
	}
*/
/*	now make pairs */
        for (i = 0; i < input_atomnum; i++) 
		input_atom_id[i]=-1;
        for (i = 0; i < ref_atomnum; i++) 
		ref_atom_id[i]=-1;
/* 	first handle non-equvalent atoms*/
	nfail = 0;
        for (i = 0; i < input_atomnum; i++) {
		if(input_equ_atom_id[i] != -1) {
			nfail++;
			continue;
		}
                for (j = 0; j < ref_atomnum; j++) {
                        if (ref_atom_id[j] != -1)
                                continue;
			if(ref_equ_atom_id[j] != -1) 
				continue;
                       	if (input_pathnum[i] != ref_pathnum[j])
                               continue;
                        matchindex = 1;
                        for (k = 0; k < input_pathnum[i]; k++)
                                if (input_pathscore[i][k] != ref_pathscore[j][k]) {
                                        matchindex = -1;
                                        break;
                                }
                        if (matchindex == 1) {
                                input_atom_id[i] = j;
				ref_atom_id[j] = i;
				break;
			}
                }
		if(input_atom_id[i] < 0) nfail ++;
	}	
	if(nfail > 0) {
		if(nfail == input_atomnum) {
/* all atoms are equalvalent*/
			input_atom_id[0]=0;
			ref_atom_id[0] = 0;
		}
		ncyc = 0;
		while(ncyc < nfail) {
			ncyc++;
			atid1 = -1;
			atid2 = -1;
        		for (i = 0; i < input_atomnum; i++) {
				if(input_atom_id[i] >= 0) continue;
				suc_flag = 0;
				for(m=0; m<6; m++) {
					atid1 = input_atom[i].con[m];
					if(atid1 < 0) break;
					if(input_atom_id[atid1] >=0) {
						atid2 = input_atom_id[atid1];
						for(n=0; n<6; n++) {
							j = ref_atom[atid2].con[n]; 
							if(ref_atom_id[j] < 0) {
                       						if (input_pathnum[i] != ref_pathnum[j])
                               						continue;
                        					matchindex = 1;
                        					for (k = 0; k < input_pathnum[i]; k++)
                                					if (input_pathscore[i][k] != ref_pathscore[j][k]) {
                                        					matchindex = -1;
                                        					break;
                                					}
                        					if (matchindex == 1) {
                                					input_atom_id[i] = j;
									ref_atom_id[j] = i;
									suc_flag = 1;
								}
							} 
/*
printf("i=%5d %5s, j=%5d %5s, atid1=%5s, atid2=%5s, suc=%d\n", i+1, input_atom[i].name, j+1, ref_atom[j].name, input_atom[atid1].name, ref_atom[atid2].name, suc_flag);
*/
							if(suc_flag == 1) break;
						}
						if(suc_flag == 1) break;
					}
					if(suc_flag == 1) break;
				}
				if(suc_flag == 1) break;
			}
		}
    	}
}
