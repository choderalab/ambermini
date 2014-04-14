# include "rotate.c"
# define MAXPATHATOMNUM 5000
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
int *selectindex;
int *pathnum;
int *selectelement;
int *pathatomnum;
double *pathscore[MAX_EQU_ATOM] ;

int selectnum = 0;
int pathnumindex = 0;
int atomindex = 0;

void geometry(ATOM *atom, int atomnum, BOND *bond, int bondnum, GEOM *geom, int *ngeom) {

int i,j,m,n;
FILE *fp;
char filename[MAXCHAR];
char line[MAXCHAR];
int at1, at2, at3, at4;
int suc;
int count;
double tor;
/* first of all, reads in parameters */
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


void scorepath(ATOM atm[], int selectnum, int startnum, int max_path_length, int ngeom, GEOM geom[])
{
	int i, j, k, m;
	int start;
	double score;
	double coef;
	int resetindex;
	start = -1;
	resetindex = -1;
	selectindex[startnum] = selectnum;
	selectelement[selectnum++] = startnum;
	if(max_path_length != -1 && selectnum > max_path_length)
		return;
	for (i = 0; i < 6; i++) {
		if (atm[startnum].con[i] == -1) {
			score = 0.0;
			for (j = 0; j < selectnum; j++) {
				coef = 1.0;
				if(ngeom > 0) {
					for(m=0; m< ngeom; m++) {
						if(j>=3 && selectelement[j - 3] == geom[m].id1 && 
                                                           selectelement[j - 2] == geom[m].id2 &&
                                                           selectelement[j - 1] == geom[m].id3 &&
                                                           selectelement[j]     == geom[m].id4) {
						 	if(geom[m].type == 1) 
								coef = 1.01;		
							else 
								coef = 0.99;
							break;
						}
						if(j>=3 && selectelement[j - 3] == geom[m].id4 && 
                                                           selectelement[j - 2] == geom[m].id3 &&
                                                           selectelement[j - 1] == geom[m].id2 &&
                                                           selectelement[j]     == geom[m].id1) {
						 	if(geom[m].type == 1) 
								coef = 1.01;		
							else 
								coef = 0.99;
							break;
						}
					}

					for(m=0; m< ngeom; m++) {
						if((selectnum - j)>3 && selectelement[j + 3] == geom[m].id4 && 
                                                                        selectelement[j + 2] == geom[m].id3 &&
                                                                        selectelement[j + 1] == geom[m].id2 &&
                                                                        selectelement[j]     == geom[m].id1) {
						 	if(geom[m].type == 1) 
								coef = 1.01;		
							else 
								coef = 0.99;
							break;
						}
						if((selectnum - j)>3 && selectelement[j + 3] == geom[m].id1 && 
                                                                        selectelement[j + 2] == geom[m].id2 &&
                                                                        selectelement[j + 1] == geom[m].id3 &&
                                                                        selectelement[j]     == geom[m].id4) {
						 	if(geom[m].type == 1) 
								coef = 1.01;		
							else 
								coef = 0.99;
							break;
						}
					}
				}
				score += coef*((j + 1) * 0.11 +
					       atm[selectelement[j]].atomicnum * 0.08);
			}
			pathscore[atomindex][pathnumindex++] = score;
			if (pathnumindex >= pathatomnum[atomindex]) {
				pathatomnum[atomindex] += MAXPATHATOMNUM;
				fprintf
					(stdout, "\nInfo: the number of the path atoms exceeds MAXPATHATOMNUM(%d) for atom[%d],extend the size and reallocate the memory automatically",
					 pathatomnum[atomindex], atomindex);
				pathscore[atomindex] =
					(double *) realloc(pathscore[atomindex],
									   pathatomnum[atomindex] *
									   sizeof(double));
				if (pathscore[atomindex] == NULL) {
					fprintf(stdout,
							" reallocate memory for pathscore[%d] failed\n",
							atomindex);
					exit(1);
				}
			}
			return;
		}
		start = atm[startnum].con[i];
		for (k = 0; k < selectnum; k++)
			if (start == selectelement[k]) {
				resetindex = 1;
				break;
			}
		if (resetindex == 1) {
			resetindex = -1;
			continue; 
		}
		if (start == -1)
			return;
		scorepath(atm, selectnum, start, max_path_length, ngeom, geom);
		/* we have already visited this atom */
	}
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

void equatom(int atomnum, ATOM *atom, int *equ_atom_id, int max_path_length, int ngeom, GEOM geom[])
{
	int i, j, k;
	int equindex;
	double sum;
	for (i = 0; i < atomnum; i++)
		equ_atom_id[i] = -1;
	for (i = 0; i < atomnum; i++) {
		selectnum = 0;
		pathnumindex = 0;
		atomindex = i;
		for (j = 0; j < atomnum; j++) {
			selectindex[j] = -1;
			selectelement[i] = -1;
		}
		scorepath(atom, 0, i, max_path_length, ngeom, geom);
		pathnum[i] = pathnumindex;
	}
	selectnum = 0;
	for (i = 0; i < atomnum; i++) {
		sum = 0.0;
		for (j = 0; j < pathnum[i]; j++) {
			sum += pathscore[i][j];
		}
	}
	for (i = 0; i < atomnum; i++) {
		sort(pathscore[i], pathnum[i]);
	}
	for (i = 0; i < atomnum; i++) {
		for (j = i + 1; j < atomnum; j++) {
			if (equ_atom_id[j] != -1)
				continue;
			equindex = 1;
			if (pathnum[i] != pathnum[j])
				continue;
			for (k = 0; k < pathnum[i]; k++)
				if (pathscore[i][k] != pathscore[j][k]) {
					equindex = -1;
					break;
				}
			if (equindex == 1)
				equ_atom_id[j] = i;
		}
	}
}

void identify_equatom(int atomnum, ATOM *atom, int *equ_atom_id, int max_path_length, int bondnum, BOND *bond, int igeom) {
int i;
int ngeom;
GEOM *geom = NULL;
        if(atomnum > MAX_EQU_ATOM) {
                fprintf(stdout, "The number of atoms (%d) exceed the MAX_EQU_ATOM (%d) defined in equatom.c, extend MAX_EQU_ATOM and recompile the program\n",
             atomnum, MAX_EQU_ATOM);
                exit(1);
        }

	selectindex = (int *) malloc(sizeof(int) * atomnum);
	if (selectindex == NULL) {
		fprintf(stdout, "memory allocation error for *selectindex\n");
		exit(1);
	}
	pathnum = (int *) malloc(sizeof(int) * atomnum);
	if (pathnum == NULL) {
		fprintf(stdout, "memory allocation error for *pathnum\n");
		exit(1);
	}
	selectelement = (int *) malloc(sizeof(int) * atomnum);
	if (selectelement == NULL) {
		fprintf(stdout, "memory allocation error for *selectelement\n");
		exit(1);
	}
	pathatomnum = (int *) malloc(sizeof(int) * atomnum);
	if (pathatomnum == NULL) {
		fprintf(stdout, "memory allocation error for *pathatomnum\n");
		exit(1);
	}
	if(atomnum > MAX_EQU_ATOM) {
		fprintf(stdout, "The number of atoms (%d) exceed the MAX_EQU_ATOM (%d) defined in equatom.c, extend MAX_EQU_ATOM and recompile the program\n", 
             atomnum, MAX_EQU_ATOM);
		exit(1);
	}
	for (i = 0; i < atomnum; i++) {
		pathatomnum[i] = MAXPATHATOMNUM;
		pathscore[i] = (double *) calloc(pathatomnum[i], sizeof(double));
		if (pathscore == NULL) {
			fprintf(stdout, "memory allocation error for *pathscore[%d]\n",
					i + 1);
			exit(1);
		}
	}
	ngeom = 0;
	if(igeom == 1)  {
		geom = (GEOM *) malloc(sizeof(GEOM) * MAXGEOM);
		if (geom == NULL) {
			fprintf(stdout, "memory allocation error for *geom\n");
			exit(1);
		}
		geometry(atom, atomnum, bond, bondnum, geom, &ngeom);
	}
	equatom(atomnum, atom, equ_atom_id, max_path_length, ngeom, geom);
}
