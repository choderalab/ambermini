# include <stdio.h>
# include <math.h>
# include <ctype.h>
# include <stdlib.h>
# include <string.h>
char *amberhome;
# include "common.h"
# include "define.h"
# include "atom.h"
# include "utility.c"
# include "common.c"
# include "ring.c"
# include "ac.c"
# include "pdb.c"
# include "mol2.c"
# include "mdl.c"
# include "rst.c"
# define MAXATOMNAME 100

typedef struct {
	double x;
	double y;
	double z;
} FITATOM;
typedef struct {
	char str[20];
} DEF_ATOMNAME;

# define COLORTEXT "YES"
# define debug 0
char ifilename[MAXCHAR];
char rfilename[MAXCHAR];
char ofilename[MAXCHAR];
char mfilename[MAXCHAR] = "match.matrix";
char lfilename[MAXCHAR] = "match.log";
char line[MAXCHAR];
ATOM *atom;
ATOM *refatom;
BOND *bond;
BOND *refbond;
AROM *arom;
AROM *refarom;
RING *ring;
RING *refring;
FITATOM *fit_atom;
FITATOM *fit_refatom;

MOLINFO minfo;
CONTROLINFO cinfo;

int i;
int jobtype = 1;
int iinput = 0;
int ioutput = 0;
int iref = 0;
int imatrix = 0;
int atomnum = 0, refatomnum = 0;
int fit_atomnum = 0, fit_refatomnum = 0;
int bondnum, refbondnum;
int stype = 0;
int nsel = -1;

FILE *fpout, *fpmatrix, *fpdef, *fplog;
double rms;
double CG0[3] = {0,0,0}; 
double CG1[3]={0,0,0}; 
double ROT[3][3]={1,0,0,0,1,0,0,0,1};
char tmpchar1[MAXCHAR];
char tmpchar2[MAXCHAR];
char def_filename[MAXCHAR];
char def_atomname_str[MAXCHAR];
int ndef_atomname = 0;
DEF_ATOMNAME def_atomname[MAXATOMNAME];

int format = 1;

void memory(int flag, int maxatom, int maxbond, int maxring)
{
        if (flag == 0) {
                atom = (ATOM *) malloc(sizeof(ATOM) * maxatom);
                if (atom == NULL) {
                        fprintf(stdout, "memory allocation error for *atom\n");
                        exit(1);
                }
                arom = (AROM *) malloc(sizeof(AROM) * maxatom);
                if (arom == NULL) {
                        fprintf(stdout, "memory allocation error for *arom\n");
                        exit(1);
                }
                bond = (BOND *) malloc(sizeof(BOND) * maxbond);
                if (bond == NULL) {
                        fprintf(stdout, "memory allocation error for *bond\n");
                        exit(1);
                }
                int i;
                for (i = 0; i < maxbond; ++i) {
                        bond[i].jflag = -1; /* bond type has not been assigned */
                }
        }
        if (flag == 1 || flag == 4 || flag == 5 || flag == 7) {
                free(atom);
                atom = (ATOM *) malloc(sizeof(ATOM) * maxatom);
                if (atom == NULL) {
                        fprintf(stdout, "memory allocation error for *atom\n");
                        exit(1);
                }
        }
        if (flag == 2 || flag == 4 || flag == 6 || flag == 7) {
                free(bond);
                bond = (BOND *) malloc(sizeof(BOND) * maxbond);
                if (bond == NULL) {
                        fprintf(stdout, "memory allocation error for *bond\n");
                        exit(1);
                }
                int i;
                for (i = 0; i < maxbond; ++i) {
                        bond[i].jflag = -1; /* bond type has not been assigned */
                }
        }
        if (flag == 3 || flag == 5 || flag == 6 || flag == 7) {
                free(arom);
                arom = (AROM *) malloc(sizeof(AROM) * maxatom);
                if (arom == NULL) {
                        fprintf(stdout, "memory allocation error for *arom\n");
                        exit(1);
                }
        }
        if (flag == 8) {
                free(ring);
                ring = (RING *) malloc(sizeof(RING) * maxring);
                if (ring == NULL) {
                        fprintf(stdout, "memory allocation error for *ring\n");
                        exit(1);
                }
        }
}

void ref_memory(int flag, int maxatom, int maxbond, int maxring)
{
        if (flag == 0) {
                refatom = (ATOM *) malloc(sizeof(ATOM) * maxatom);
                if (refatom == NULL) {
                        fprintf(stdout, "memory allocation error for *refatom\n");
                        exit(1);
                }
                refarom = (AROM *) malloc(sizeof(AROM) * maxatom);
                if (refarom == NULL) {
                        fprintf(stdout, "memory allocation error for *refarom\n");
                        exit(1);
                }
                refbond = (BOND *) malloc(sizeof(BOND) * maxbond);
                if (refbond == NULL) {
                        fprintf(stdout, "memory allocation error for *refbond\n");
                        exit(1);
                }
                int i;
                for (i = 0; i < maxbond; ++i) {
                        refbond[i].jflag = -1; /* bond type has not been assigned */
                }
        }
        if (flag == 1 || flag == 4 || flag == 5 || flag == 7) {
                free(refatom);
                refatom = (ATOM *) malloc(sizeof(ATOM) * maxatom);
                if (refatom == NULL) {
                        fprintf(stdout, "memory allocation error for *refatom\n");
                        exit(1);
                }
        }
        if (flag == 2 || flag == 4 || flag == 6 || flag == 7) {
                free(refbond);
                refbond = (BOND *) malloc(sizeof(BOND) * maxbond);
                if (refbond == NULL) {
                        fprintf(stdout, "memory allocation error for *refbond\n");
                        exit(1);
                }
                int i;
                for (i = 0; i < maxbond; ++i) {
                        refbond[i].jflag = -1; /* bond type has not been assigned */
                }
        }
        if (flag == 3 || flag == 5 || flag == 6 || flag == 7) {
                free(refarom);
                refarom = (AROM *) malloc(sizeof(AROM) * maxatom);
                if (refarom == NULL) {
                        fprintf(stdout, "memory allocation error for *refarom\n");
                        exit(1);
                }
        }
        if (flag == 8) {
                free(refring);
                refring = (RING *) malloc(sizeof(RING) * maxring);
                if (refring == NULL) {
                        fprintf(stdout, "memory allocation error for *refring\n");
                        exit(1);
                }
        }
}

int rmatrix(char *filename) {
	char line[MAXCHAR];
	FILE *fpin;

	if ((fpin = fopen(filename, "r")) == NULL) {
		fprintf(stderr, "Cannot open file %s, exit\n", filename);
		exit(1);
	}
	for (;;) {
		if (fgets(line, MAXCHAR, fpin) == NULL) break;
		if (strncmp("CG0", line, 3) == 0) sscanf(&line[3], "%lf%lf%lf", &CG0[0], &CG0[1], &CG0[2]);
		if (strncmp("CG1", line, 3) == 0) sscanf(&line[3], "%lf%lf%lf", &CG1[0], &CG1[1], &CG1[2]);
		if (strncmp("ROT0", line, 4) == 0) sscanf(&line[4], "%lf%lf%lf", &ROT[0][0], &ROT[0][1], &ROT[0][2]);
		if (strncmp("ROT1", line, 4) == 0) sscanf(&line[4], "%lf%lf%lf", &ROT[1][0], &ROT[1][1], &ROT[1][2]);
		if (strncmp("ROT2", line, 4) == 0) sscanf(&line[4], "%lf%lf%lf", &ROT[2][0], &ROT[2][1], &ROT[2][2]);
	}
}

double lsfit(FITATOM atm1[], FITATOM atm2[], int n) {
//atm1 - reference
//atm2 - data
int i,j,k;
int ict;
int ix,iy,iz;
int iflag;
/* ict: counter for number of iterations. must do at least 3 iterations.
ix, iy,iz: pointer used in iterative least squares. may have the value 1,2,3
the value changes on each iteration. */
double tol=0.00001;
double sig,gam, sg, bb,cc;
double **x0;
double **x1;
double **xyzfit;
double aa[3][3];
double an;
double rms,ssum;

x0 = (double**)malloc(3*sizeof(double*));
for(i = 0 ; i<3 ; i++) {
        *(x0+i) = (double*) malloc(fit_atomnum*sizeof(double));
        if (*(x0+i) == NULL) {
                fprintf(stderr, "memory allocation error for *x0[%5d]\n", i);
                exit(0);
        }
}
if (x0 == NULL) {
        fprintf(stderr, "memory allocation error for **x0\n");
        exit(0);
}

x1 = (double**)malloc(3*sizeof(double*));
for(i = 0 ; i<3 ; i++) {
        *(x1+i) = (double*) malloc(fit_atomnum*sizeof(double));
        if (*(x1+i) == NULL) {
                fprintf(stderr, "memory allocation error for *x1[%5d]\n", i);
                exit(0);
        }
}
if (x1 == NULL) {
        fprintf(stderr, "memory allocation error for **x1\n");
        exit(0);
}

xyzfit = (double**)malloc(3*sizeof(double*));
for(i = 0 ; i<3 ; i++) {
        *(xyzfit+i) = (double*) malloc(fit_atomnum*sizeof(double));
        if (*(xyzfit+i) == NULL) {
                fprintf(stderr, "memory allocation error for *xyzfit[%5d]\n", i);
                exit(0);
        }
}
if (xyzfit == NULL) {
        fprintf(stderr, "memory allocation error for **xyzfit\n");
        exit(0);
}

  ssum=0.;
  for(i=0;i<n;i++)
   {
     x1[0][i]=atm2[i].x;
     x1[1][i]=atm2[i].y;
     x1[2][i]=atm2[i].z;
   }
  for(i=0;i<n;i++)
   {
    x0[0][i]=atm1[i].x;
    x0[1][i]=atm1[i].y;
    x0[2][i]=atm1[i].z;
   }
  an=1.0/(double)n;
  for(i=0;i<3;i++)
  {
  CG0[i]=0.;
  CG1[i]=0.;
  }
  for(i=0;i<n;i++)
   for(j=0;j<3;j++)
    {
    CG0[j]=CG0[j]+x0[j][i];
    CG1[j]=CG1[j]+x1[j][i];
   }
   for(i=0;i<3;i++)
   {
     CG0[i]*=an;
     CG1[i]*=an;
   }
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
     aa[i][j]=0.0;
  for(k=0;k<n;k++)
  {
   aa[0][0]=aa[0][0]+(x1[0][k]-CG1[0])*(x0[0][k]-CG0[0]);
   aa[1][0]=aa[1][0]+(x1[1][k]-CG1[1])*(x0[0][k]-CG0[0]);
   aa[2][0]=aa[2][0]+(x1[2][k]-CG1[2])*(x0[0][k]-CG0[0]);
   aa[0][1]=aa[0][1]+(x1[0][k]-CG1[0])*(x0[1][k]-CG0[1]);
   aa[1][1]=aa[1][1]+(x1[1][k]-CG1[1])*(x0[1][k]-CG0[1]);
   aa[2][1]=aa[2][1]+(x1[2][k]-CG1[2])*(x0[1][k]-CG0[1]);
   aa[0][2]=aa[0][2]+(x1[0][k]-CG1[0])*(x0[2][k]-CG0[2]);
   aa[1][2]=aa[1][2]+(x1[1][k]-CG1[1])*(x0[2][k]-CG0[2]);
   aa[2][2]=aa[2][2]+(x1[2][k]-CG1[2])*(x0[2][k]-CG0[2]);
  }
  for(i=0;i<3;i++)
    {
    for(j=0;j<3;j++)
      ROT[i][j]=0.;
    ROT[i][i]=1.0;
   }
  ict=0.;
  goto a51;
  a50:
  ix++;
  if(ix<4) goto a52;
  if(iflag==0) goto a70;
  a51:
  iflag=0;
  ix=1;
  a52:
  ict=ict+1;
  if(ict>1000) goto a70;
  iy=ix+1;
  if(iy==4)  iy=1;
  iz=6-ix-iy;
  sig=aa[iz-1][iy-1]-aa[iy-1][iz-1];
  gam=aa[iy-1][iy-1]+aa[iz-1][iz-1];
  sg=sqrt(sig*sig+gam*gam);
  if(sg==0) goto a50;
  sg=1.0/sg;
  if(fabs(sig)<tol*fabs(gam)) goto a50;
  for(k=0;k<3;k++)
  {
   bb=gam*aa[iy-1][k]+sig*aa[iz-1][k];
   cc=gam*aa[iz-1][k]-sig*aa[iy-1][k];
   aa[iy-1][k]=bb*sg;
   aa[iz-1][k]=cc*sg;
   bb=gam*ROT[iy-1][k]+sig*ROT[iz-1][k];
   cc=gam*ROT[iz-1][k]-sig*ROT[iy-1][k];
   ROT[iy-1][k]=bb*sg;
   ROT[iz-1][k]=cc*sg;
  }
  iflag=1;
  goto a50;
  a70:
/* the following code translates the ligand the center of mass of the
receptor site and rotates it according to the rotation matrxouneror ix
calculated by orient */
for(i=0;i<n;i++) {
 	xyzfit[0][i]=CG0[0]+ROT[0][0]*(x1[0][i]-CG1[0])+ \
             ROT[0][1]*(x1[1][i]-CG1[1])+
             ROT[0][2]*(x1[2][i]-CG1[2]);
 	xyzfit[1][i]=CG0[1]+ROT[1][0]*(x1[0][i]-CG1[0])+ \
             ROT[1][1]*(x1[1][i]-CG1[1])+
             ROT[1][2]*(x1[2][i]-CG1[2]);
 	xyzfit[2][i]=CG0[2]+ROT[2][0]*(x1[0][i]-CG1[0])+ \
             ROT[2][1]*(x1[1][i]-CG1[1])+
             ROT[2][2]*(x1[2][i]-CG1[2]);
}

for(j=0;j<n;j++)
  for(i=0;i<3;i++)
 	ssum+=(xyzfit[i][j]-x0[i][j])*(xyzfit[i][j]-x0[i][j]);
rms=ssum/n;
rms=sqrt(rms);
//	write out transition matrix

fprintf(fpmatrix, "CG0\t%16.10lf%16.10lf%16.10lf\n", CG0[0], CG0[1], CG0[2]);
fprintf(fpmatrix, "CG1\t%16.10lf%16.10lf%16.10lf\n", CG1[0], CG1[1], CG1[2]);
fprintf(fpmatrix, "ROT0\t%16.10lf%16.10lf%16.10lf\n", ROT[0][0], ROT[0][1], ROT[0][2]);
fprintf(fpmatrix, "ROT1\t%16.10lf%16.10lf%16.10lf\n", ROT[1][0], ROT[1][1], ROT[1][2]);
fprintf(fpmatrix, "ROT2\t%16.10lf%16.10lf%16.10lf\n", ROT[2][0], ROT[2][1], ROT[2][2]);

/*
 for(i=0;i<n;i++)
 {
  atm3[i].x=xyzfit[0][i];
  atm3[i].y=xyzfit[1][i];
  atm3[i].z=xyzfit[2][i];
 }
*/
return rms;
}

int main(int argc, char *argv[]) {
	int i,j,k;
	int overflow_flag = 0;
	int count;
	int atid1, atid2;
	int resid1, resid2;
	int oldresid;
	int suc;
	double sum;
	double x,y,z;
	char line[MAXCHAR];
	char defstr[MAXCHAR];
	char atomname[MAXCHAR];

    	amberhome = (char *) getenv("AMBERHOME");
   	if( amberhome == NULL ){
       		fprintf( stdout, "AMBERHOME is not set!\n" );
       		exit(1);
    	}

	if (strcmp(COLORTEXT, "YES") == 0 || strcmp(COLORTEXT, "yes") == 0) {
		if (argc == 2
			&& (strcmp(argv[1], "-h") == 0
				|| strcmp(argv[1], "-H") == 0)) {
			printf
				("[31mUsage: match  -i [0m input file name \n"
				 "[31m              -r [0m reference file name \n"
				 "[31m              -f [0m format: 1-pdb (the default), 2-ac, 3-mol2, 4-sdf, 5-crd/rst\n"
				 "[31m              -o [0m output file name\n"
				 "[31m              -l [0m run log file name, default is \"match.log\"\n"
				 "[31m              -s [0m selection mode\n"
				 "[34m                  0:[0m use all atoms (the default)\n"
				 "[34m                  1:[0m specify atom names\n"
				 "[34m                  2:[0m use atom defination file\n"
				 "[34m                  3:[0m use residue defination file - original residue IDs\n"
				 "[34m                  4:[0m use residue defination file - renumbered residue IDs\n"
				 "[31m              -ds [0mdefinition string if selection modes of '1' or '3' or '4'\n"
				 "[31m                  [0me.g. 'C,N,O,CA', or 'HET' which stands for heavy atoms for '-ds 1')\n"
				 "[31m              -df [0mdefinition file if selection mode of '2' or '3' or '4'\n"
				 "[31m                  [0mrecords take a form of 'ATOM atom_id_input atom_id_reference'\n"
				 "[31m                  [0mor 'RES res_id_input res_id_reference'\n"
				 "[31m              -n [0m number of atoms participating ls-fitting,\n"
				 "[31m                 [0m default is -1, which implies to use all the selected atoms\n"
				 "[31m              -m [0m matrix file, default is \"match.matrix\"\n"
				 "[31m              -t [0m job type:\n"
				 "[35m                  0:[0m calculate rms only, need -i and -r\n"    
				 "[35m                  1:[0m lsfit, need -i, -r and -o the default\n"    
				 "[35m                  2:[0m translation/rotation, need -i, -o and -m\n"); 
			exit(0);
		}
		if (argc != 23 && argc != 21 && argc != 19 && argc != 17 && argc != 15 && argc != 13 && argc != 11 && argc !=9 && argc !=7 && argc !=5) {
			printf
				("[31mUsage: match  -i [0m input file name \n"
				 "[31m              -r [0m reference file name \n"
				 "[31m              -f [0m format: 1-pdb (the default), 2-ac, 3-mol2, 4-sdf, 5-crd/rst\n"
				 "[31m              -o [0m output file name\n"
				 "[31m              -l [0m run log file name, default is \"match.log\"\n"
				 "[31m              -s [0m selection mode\n"
				 "[34m                  0:[0m use all atoms (the default)\n"
				 "[34m                  1:[0m specify atom names\n"
				 "[34m                  2:[0m use atom defination file\n"
				 "[34m                  3:[0m use residue defination file - original residue IDs\n"
				 "[34m                  4:[0m use residue defination file - renumbered residue IDs\n"
				 "[31m              -ds [0mdefinition string if selection modes of '1' or '3' or '4'\n"
				 "[31m                  [0me.g. 'C,N,O,CA', or 'HET' which stands for heavy atoms for '-ds 1')\n"
				 "[31m              -df [0mdefinition file if selection mode of '2' or '3' or '4'\n"
				 "[31m                  [0mrecords take a form of 'ATOM atom_id_input atom_id_reference'\n"
				 "[31m                  [0mor 'RES res_id_input res_id_reference'\n"
				 "[31m              -n [0m number of atoms participating ls-fitting,\n"
				 "[31m                 [0m default is -1, which implies to use all the selected atoms\n"
				 "[31m              -m [0m matrix file, default is \"match.matrix\"\n"
				 "[31m              -t [0m job type:\n"
				 "[35m                  0:[0m calculate rms only, need -i and -r\n"    
				 "[35m                  1:[0m lsfit, need -i, -r and -o the default\n"    
				 "[35m                  2:[0m translation/rotation, need -i, -o and -m\n"); 
			exit(0);
		}
	} else {
		if (argc == 2
			&& (strcmp(argv[1], "-h") == 0
				|| strcmp(argv[1], "-H") == 0)) {
			printf
				("Usage: match  -i  input file name \n"
				 "              -r  reference file name \n"
				 "              -f  format: 1-pdb (the default), 2-ac, 3-mol2, 4-sdf, 5-crd/rst\n"
				 "              -o  output file name\n"
				 "              -l  run log file name, default is \"match.log\"\n"
				 "              -s  selection mode\n"
				 "                  0: use all atoms (the default)\n"
				 "                  1: specify atom names\n"
				 "                  2: use atom defination file\n"
				 "                  3: use residue defination file - original residue IDs\n"
				 "                  4: use residue defination file - renumbered residue IDs\n"
				 "              -ds definition string if selection modes of '1' or '3' or '4'\n"
				 "                  e.g. 'C,N,O,CA', or 'HET' which stands for heavy atoms for '-ds 1')\n"
				 "              -df definition file if selection mode of '2' or '3' or '4'\n"
				 "                  records take a form of 'ATOM atom_id_input atom_id_reference'\n"
				 "                  or 'RES res_id_input res_id_reference'\n"
				 "              -n  number of atoms participating ls-fitting,\n"
				 "                  default is -1, which implies to use all the selected atoms\n"
				 "              -m  matrix file, default is \"match.matrix\"\n"
				 "              -t  job type:\n"
				 "                  0: calculate rms only, need -i and -r\n"    
				 "                  1: lsfit, need -i, -r and -o the default\n"    
				 "                  2: translation/rotation, need -i, -o and -m\n"); 
			exit(0);
		}
		if (argc != 23 && argc != 21 && argc != 19 && argc != 17 && argc != 15 && argc != 13 && argc != 11 && argc !=9 && argc !=7 && argc !=5) {
			printf
				("Usage: match  -i  input file name \n"
				 "              -r  reference file name \n"
				 "              -f  format: 1-pdb (the default), 2-ac, 3-mol2, 4-sdf, 5-crd/rst\n"
				 "              -o  output file name\n"
				 "              -l  run log file name, default is \"match.log\"\n"
				 "              -s  selection mode\n"
				 "                  0: use all atoms (the default)\n"
				 "                  1: specify atom names\n"
				 "                  2: use atom defination file\n"
				 "                  3: use residue defination file - original residue IDs\n"
				 "                  4: use residue defination file - renumbered residue IDs\n"
				 "              -ds definition string if selection modes of '1' or '3' or '4'\n"
				 "                  e.g. 'C,N,O,CA', or 'HET' which stands for heavy atoms for '-ds 1')\n"
				 "              -df definition file if selection mode of '2' or '3' or '4'\n"
				 "                  records take a form of 'ATOM atom_id_input atom_id_reference'\n"
				 "                  or 'RES res_id_input res_id_reference'\n"
				 "              -n  number of atoms participating ls-fitting,\n"
				 "                  default is -1, which implies to use all the selected atoms\n"
				 "              -m  matrix file, default is \"match.matrix\"\n"
				 "              -t  job type:\n"
				 "                  0: calculate rms only, need -i and -r\n"    
				 "                  1: lsfit, need -i, -r and -o the default\n"    
				 "                  2: translation/rotation, need -i, -o and -m\n"); 
			exit(0);
		}

	}

	for (i = 1; i < argc; i += 2) {
		if (strcmp(argv[i], "-i") == 0) {
			strcpy(ifilename, argv[i + 1]);
			iinput = 1;
		}
		if (strcmp(argv[i], "-f") == 0) 
			format=atoi(argv[i+1]);
		if (strcmp(argv[i], "-o") == 0) {
			strcpy(ofilename, argv[i + 1]);
			ioutput = 1;
		}
		if (strcmp(argv[i], "-r") == 0) {
			strcpy(rfilename, argv[i + 1]);
			iref = 1;
		}
		if (strcmp(argv[i], "-s") == 0) 
			stype = atoi(argv[i+1]);
		if (strcmp(argv[i], "-ds") == 0) 
			strcpy(def_atomname_str, argv[i + 1]);
		if (strcmp(argv[i], "-df") == 0) 
			strcpy(def_filename, argv[i + 1]);
		if (strcmp(argv[i], "-m") == 0) {
			strcpy(mfilename, argv[i + 1]);
			imatrix = 1;
		}
		if (strcmp(argv[i], "-t") == 0)
			jobtype = atoi(argv[i+1]);
		if (strcmp(argv[i], "-n") == 0)
			nsel = atoi(argv[i+1]);
		if (strcmp(argv[i], "-l") == 0) 
			strcpy(lfilename, argv[i + 1]);
	}

	if(jobtype !=0 && jobtype != 1 && jobtype != 2) jobtype = 1;
	if(format != 1 && format != 2 && format != 3 && format != 4 && format != 5) format = 1;
	if(jobtype == 0) 
		if(iinput == 0 || iref == 0) {
			fprintf(stderr, "RMS calculation needs an input and a reference files, exit\n");
			exit(0);
		}
	if(jobtype == 1) 
		if(iinput == 0 || ioutput == 0 || iref == 0) {
			fprintf(stderr, "The least-square fitting calculation needs an input, an output and a reference files, exit\n");
			exit(0);
		}

	if(jobtype == 2) 
		if(iinput == 0 || ioutput == 0) {
			fprintf(stderr, "The translation/rotation job needs an input and an output files, exit\n");
			exit(0);
		}
	if(jobtype == 0 || jobtype == 1) {
		if ((fplog = fopen(lfilename, "w")) == NULL) {
			fprintf(stderr, "Cannot open log file %s to write, exit\n", lfilename);
			exit(1);
		}
	}
/*	read in input files*/
        memory(0, MAXATOM, MAXBOND, MAXRING);
	if(format == 1) {
                overflow_flag = rpdb(ifilename, &atomnum, atom, cinfo, minfo, 0);
                if (overflow_flag) {
                        cinfo.maxatom = atomnum + 10;
                        cinfo.maxbond = bondnum + 10;
                        memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
                        overflow_flag = rpdb(ifilename, &atomnum, atom, cinfo, minfo, 0);
                }
	}
	if(format == 2) {
                overflow_flag = rac(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo);
                if (overflow_flag) {
                        cinfo.maxatom = atomnum + 10;
                        cinfo.maxbond = bondnum + 10;
                        memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
                	overflow_flag = rac(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo);
                }
	}
	if(format == 3) {
                overflow_flag = rmol2(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo, 0);
                if (overflow_flag) {
                        cinfo.maxatom = atomnum + 10;
                        cinfo.maxbond = bondnum + 10;
                        memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
                        overflow_flag =
                		rmol2(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo, 0);
                }
	}
	if(format == 4) {
                overflow_flag = rmdl(ifilename, &atomnum, atom, &bondnum, bond, cinfo, minfo);
                if (overflow_flag) {
                        cinfo.maxatom = atomnum + 10;
                        cinfo.maxbond = bondnum + 10;
                        memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
                        overflow_flag = rmdl(ifilename, &atomnum, atom, &bondnum, bond, cinfo, minfo);
                }
	}
	if(format == 5) {
                overflow_flag = rrst(ifilename, &atomnum, atom, cinfo);
                if (overflow_flag) {
                        cinfo.maxatom = atomnum + 10;
                        memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
                        overflow_flag = rrst(ifilename, &atomnum, atom, cinfo);
                }
	}
	free(bond);
	free(ring);
	free(arom);
	cinfo.maxatom = MAXATOM;
	cinfo.maxbond = MAXBOND;
	cinfo.maxring = MAXRING;
        ref_memory(0, MAXATOM, MAXBOND, MAXRING);
	if(format == 1) {
                overflow_flag = rpdb(rfilename, &refatomnum, refatom, cinfo, minfo, 0);
                if (overflow_flag) {
                        cinfo.maxatom = refatomnum + 10;
                        cinfo.maxbond = refbondnum + 10;
                        ref_memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
                        overflow_flag = rpdb(rfilename, &refatomnum, refatom, cinfo, minfo, 0);
                }
	}
	if(format == 2) {
               	overflow_flag = rac(rfilename, &refatomnum, refatom, &refbondnum, refbond, &cinfo, &minfo);
                if (overflow_flag) {
                        cinfo.maxatom = refatomnum + 10;
                        cinfo.maxbond = refbondnum + 10;
                        ref_memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
               		overflow_flag = rac(rfilename, &refatomnum, refatom, &refbondnum, refbond, &cinfo, &minfo);
                }
	}
	if(format == 3) {
                overflow_flag = rmol2(rfilename, &refatomnum, refatom, &refbondnum, refbond, &cinfo, &minfo, 0);
                if (overflow_flag) {
                        cinfo.maxatom = refatomnum + 10;
                        cinfo.maxbond = refbondnum + 10;
                        ref_memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
                        overflow_flag = rmol2(rfilename, &refatomnum, refatom, &refbondnum, refbond, &cinfo, &minfo, 0);
                }
	}
	if(format == 4) {
                overflow_flag = rmdl(rfilename, &refatomnum, refatom, &refbondnum, refbond, cinfo, minfo);
                if (overflow_flag) {
                        cinfo.maxatom = refatomnum + 10;
                        cinfo.maxbond = refbondnum + 10;
                        ref_memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
                        overflow_flag = rmdl(rfilename, &refatomnum, refatom, &refbondnum, refbond, cinfo, minfo);
                }
	}
	if(format == 5) {
                overflow_flag = rrst(rfilename, &refatomnum, refatom, cinfo);
                if (overflow_flag) {
                        cinfo.maxatom = refatomnum + 10;
                        ref_memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
                        overflow_flag = rrst(rfilename, &refatomnum, refatom, cinfo);
                }
	}
	free(refbond);
	free(refring);
	free(refarom);
	if(debug == 1) {
		for(i=0;i<atomnum;i++)
			printf("ATOM %5d %5s %5s %5d %8.3lf %8.3lf %8.3lf\n", i+1, atom[i].name, atom[i].aa, atom[i].resno, atom[i].x, atom[i].y, atom[i].z);
		for(i=0;i<refatomnum;i++)
			printf("ATOM %5d %5s %5s %5d %8.3lf %8.3lf %8.3lf\n", i+1, refatom[i].name, refatom[i].aa, refatom[i].resno, refatom[i].x, refatom[i].y, refatom[i].z);
	}
	
if(jobtype == 0 || jobtype == 1) {
        fit_atom = (FITATOM *) malloc(sizeof(FITATOM) * atomnum);
        if (fit_atom == NULL) {
        	fprintf(stderr, "memory allocation error for *fit_atom\n");
                exit(1);
        }
        fit_refatom = (FITATOM *) malloc(sizeof(FITATOM) * refatomnum);
        if (fit_refatom == NULL) {
        	fprintf(stderr, "memory allocation error for *fit_refatom\n");
                exit(1);
        }
	if(stype == 0) {
		for(i=0;i<atomnum;i++) {
			fit_atom[i].x=atom[i].x;
			fit_atom[i].y=atom[i].y;
			fit_atom[i].z=atom[i].z;
		}
		for(i=0;i<refatomnum;i++) {
			fit_refatom[i].x=refatom[i].x;
			fit_refatom[i].y=refatom[i].y;
			fit_refatom[i].z=refatom[i].z;
		}
		fit_atomnum = atomnum;
		fit_refatomnum = refatomnum;
	}
	if(stype == 1) {
		count = 0;
		for(i=0;i <= strlen(def_atomname_str);i++) {
			if(def_atomname_str[i]==' ') continue;
			if(def_atomname_str[i]==',') def_atomname_str[i] = ' ';
			defstr[count] = def_atomname_str[i];
			count++;
		}
		sscanf(defstr, "%s", atomname);	
		strcpy(def_atomname[ndef_atomname].str, atomname);
		ndef_atomname++;
		for(i=1;i<strlen(defstr);i++) 
			if(defstr[i] == ' ') {
				sscanf(&defstr[i], "%s", atomname);
				strcpy(def_atomname[ndef_atomname].str, atomname);
				ndef_atomname++;
			}
		for(i=0;i<atomnum;i++) {
			suc = 0;
			for(j=0;j<ndef_atomname;j++) {
				if(strcmp(def_atomname[j].str, atom[i].name) == 0) {
					suc = 1;
					break;
				}
				if(strcmp(def_atomname[j].str, "HET") == 0) 
					if(atom[i].name[0] != 'H') {
						suc = 1;
						break;
					}
			}
			if(suc == 1) {
				fit_atom[fit_atomnum].x=atom[i].x;
				fit_atom[fit_atomnum].y=atom[i].y;
				fit_atom[fit_atomnum].z=atom[i].z;
				fit_atomnum++;
			}
				
		}
		for(i=0;i<refatomnum;i++) {
			suc = 0;
			for(j=0;j<ndef_atomname;j++) {
				if(strcmp(def_atomname[j].str, refatom[i].name) == 0) {
					suc = 1;
					break;
				}
				if(strcmp(def_atomname[j].str, "HET") == 0) 
					if(refatom[i].name[0] != 'H') {
						suc = 1;
						break;
					}
			}
			if(suc == 1) {
				fit_refatom[fit_refatomnum].x=refatom[i].x;
				fit_refatom[fit_refatomnum].y=refatom[i].y;
				fit_refatom[fit_refatomnum].z=refatom[i].z;
				fit_refatomnum++;
			}
		}
	}
	if(stype == 2) {
		if ((fpdef = fopen(def_filename, "r")) == NULL) {
			fprintf(stderr, "Cannot open file %s, exit\n", def_filename);
			exit(1);
		}
		count = 0;
		for(;;) {
			if (fgets(line, MAXCHAR, fpdef) == NULL) break;
			if(strncmp(line, "ATOM", 4) == 0) {
				sscanf(&line[4], "%d%d", &atid1, &atid2);
				atid1--;
				atid2--;
				fit_atom[count].x=atom[atid1].x;
				fit_atom[count].y=atom[atid1].y;
				fit_atom[count].z=atom[atid1].z;
				fit_refatom[count].x=refatom[atid2].x;
				fit_refatom[count].y=refatom[atid2].y;
				fit_refatom[count].z=refatom[atid2].z;
				count++;
			}
		}
		fit_atomnum = count;
		fit_refatomnum = count;
		fclose(fpdef);
	}

	if(stype == 3 || stype == 4) {
/*first of all, read in atom name definition strings*/
		count = 0;
		for(i=0;i <= strlen(def_atomname_str);i++) {
			if(def_atomname_str[i]==' ') continue;
			if(def_atomname_str[i]==',') def_atomname_str[i] = ' ';
			defstr[count] = def_atomname_str[i];
			count++;
		}
		sscanf(defstr, "%s", atomname);	
		strcpy(def_atomname[ndef_atomname].str, atomname);
		ndef_atomname++;
		for(i=1;i<strlen(defstr);i++) 
			if(defstr[i] == ' ') {
				sscanf(&defstr[i], "%s", atomname);
				if(strcmp(atomname, "HET") == 0) {
					fprintf(stderr, "Warning: HET cannot be defined with -s of 3\n");
					continue;
				}
				strcpy(def_atomname[ndef_atomname].str, atomname);
				ndef_atomname++;
			}
		if(debug == 1) {
			for(i=0;i<ndef_atomname;i++)
				printf("DEF %5d %5s\n", i+1, def_atomname[i].str);

		}
		if(stype == 4) {
			count = 0;
			oldresid = -99999;
			for(i=0;i<atomnum;i++) {
				if(atom[i].resno != oldresid) {
					oldresid = atom[i].resno;
					count++;
				}
				atom[i].resno = count;
			}

			count = 0;
			oldresid = -99999;
			for(i=0;i<refatomnum;i++) {
				if(refatom[i].resno != oldresid) {
					oldresid = refatom[i].resno;
					count++;
				}
				refatom[i].resno = count;
			}
		}
/* now working on residue definition file */
		if ((fpdef = fopen(def_filename, "r")) == NULL) {
			fprintf(stderr, "Cannot open file %s, exit\n", def_filename);
			exit(1);
		}
		count = 0;
		for(;;) {
			if (fgets(line, MAXCHAR, fpdef) == NULL) break;
			if(strncmp(line, "RES", 3) == 0) {
				sscanf(&line[3], "%d%d", &resid1, &resid2);
				if(debug == 1) 
					printf("RES %5d %5d\n", resid1, resid2);	
				for(i=0; i< ndef_atomname; i++) {
					atid1 = -1;
					atid2 = -1;
					for(j=0;j<atomnum;j++) 
						if(strcmp(def_atomname[i].str, atom[j].name) == 0 && atom[j].resno == resid1) {
							atid1 = j;	
							break;
						}
					for(j=0;j<refatomnum;j++) 
						if(strcmp(def_atomname[i].str, refatom[j].name) == 0 && refatom[j].resno == resid2) {
							atid2 = j;	
							break;
						}
					
					if(atid1 >= 0 && atid2 >= 0) {	
						fit_atom[count].x=atom[atid1].x;
						fit_atom[count].y=atom[atid1].y;
						fit_atom[count].z=atom[atid1].z;
						fit_refatom[count].x=refatom[atid2].x;
						fit_refatom[count].y=refatom[atid2].y;
						fit_refatom[count].z=refatom[atid2].z;
						count++;
					}
				}
			}
		}

		fit_atomnum = count;
		fit_refatomnum = count;
		fclose(fpdef);
	}


	if(fit_atomnum != fit_refatomnum || fit_atomnum <=0) {
		fprintf(stderr, "\nThe numumber of fitting atoms in input file (%d) is different from the number of fitting atoms in ref file (%d)", fit_atomnum, fit_refatomnum);
		exit(0);
	}
	if(nsel > 0 && fit_atomnum > nsel) fit_atomnum = nsel;	
	if(fit_atomnum <= 0) {
		fprintf(stderr, "\nThe number of fitting atoms is smaller than or equal to zero, exit!");
		exit(1);
	}
	sum = 0;
	for(i=0;i<fit_atomnum;i++) {
		sum+=(fit_atom[i].x-fit_refatom[i].x) * (fit_atom[i].x-fit_refatom[i].x); 
		sum+=(fit_atom[i].y-fit_refatom[i].y) * (fit_atom[i].y-fit_refatom[i].y); 
		sum+=(fit_atom[i].z-fit_refatom[i].z) * (fit_atom[i].z-fit_refatom[i].z); 
	}
	rms = sqrt(sum/fit_atomnum);
	fprintf(fplog, "The rmsd before least-square fitting is %9.4lf\n", rms);
	fprintf(fplog, "Number of atoms participate least-square fitting is %9d\n", fit_atomnum);
	fprintf(stdout, "\nThe rmsd before least-square fitting is %9.4lf", rms);
	fprintf(stdout, "\nNumber of atoms participate least-square fitting is %9d", fit_atomnum);
	if(jobtype == 0) {
		fclose(fplog);
		return 0;
	}
	if ((fpmatrix = fopen(mfilename, "w")) == NULL) {
		fprintf(stderr, "Cannot open file %s, exit\n", mfilename);
		exit(1);
	}
	rms = lsfit(fit_refatom, fit_atom, fit_atomnum) ;
	for(i=0;i<fit_atomnum;i++) {
               	x=CG0[0]+ROT[0][0]*(fit_atom[i].x-CG1[0])+ ROT[0][1]*(fit_atom[i].y-CG1[1])+ ROT[0][2]*(fit_atom[i].z-CG1[2]);
               	y=CG0[1]+ROT[1][0]*(fit_atom[i].x-CG1[0])+ ROT[1][1]*(fit_atom[i].y-CG1[1])+ ROT[1][2]*(fit_atom[i].z-CG1[2]);
               	z=CG0[2]+ROT[2][0]*(fit_atom[i].x-CG1[0])+ ROT[2][1]*(fit_atom[i].y-CG1[1])+ ROT[2][2]*(fit_atom[i].z-CG1[2]);

		sum =(x-fit_refatom[i].x) * (x-fit_refatom[i].x);
		sum+=(y-fit_refatom[i].y) * (y-fit_refatom[i].y);
		sum+=(z-fit_refatom[i].z) * (z-fit_refatom[i].z);
		sum = sqrt(sum);
		fprintf(fplog, "ATOM    %5d %9.3lf %9.3lf %9.3lf ", i+1,  fit_atom[i].x, fit_atom[i].y, fit_atom[i].z);
		fprintf(fplog, "REFATOM %9.3lf %9.3lf %9.3lf %9.3lf\n", fit_refatom[i].x, fit_refatom[i].y, fit_refatom[i].z, sum);
	}
	fclose(fpmatrix);
	fprintf(fplog, "\nThe rmsd after least-square fitting is %9.4lf\n", rms);
	fprintf(stdout, "\nThe rmsd after least-square fitting is %9.4lf\n", rms);
	fclose(fplog);
}
if(jobtype == 2) 
	rmatrix(mfilename);
if(jobtype == 1 || jobtype == 2) {
	for(i=0;i<atomnum;i++) {
                x=CG0[0]+ROT[0][0]*(atom[i].x-CG1[0])+ ROT[0][1]*(atom[i].y-CG1[1])+ ROT[0][2]*(atom[i].z-CG1[2]);
                y=CG0[1]+ROT[1][0]*(atom[i].x-CG1[0])+ ROT[1][1]*(atom[i].y-CG1[1])+ ROT[1][2]*(atom[i].z-CG1[2]);
                z=CG0[2]+ROT[2][0]*(atom[i].x-CG1[0])+ ROT[2][1]*(atom[i].y-CG1[1])+ ROT[2][2]*(atom[i].z-CG1[2]);
		atom[i].x = x;
		atom[i].y = y;
		atom[i].z = z;
	}
	if(format == 1) wpdb(ofilename, atomnum, atom);
	if(format == 2) wac(ofilename, atomnum, atom, bondnum, bond, cinfo, minfo);
	if(format == 3) wmol2(ofilename, atomnum, atom, bondnum, bond, arom, cinfo, minfo);
	if(format == 4) wmdl(ofilename, atomnum, atom, bondnum, bond, cinfo);
	if(format == 5) wrst(ofilename, atomnum, atom);
}
return 0;
}
