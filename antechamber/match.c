# include <stdio.h>
# include <math.h>
# include <ctype.h>
# include <stdlib.h>
# include <string.h>
# define MAXATOM 10000
# define MAXCHAR 256
# define COLORTEXT "YES"
typedef struct {
        char name[5];
        double x;
        double y;
        double z;
} ATOM;

char ifilename[MAXCHAR];
char rfilename[MAXCHAR];
char ofilename[MAXCHAR];
char mfilename[MAXCHAR];
char line[MAXCHAR];
ATOM atom1[MAXATOM];
ATOM atom2[MAXATOM];
ATOM atom3[MAXATOM];
int i;
int numarg = 0;
int jobtype = 1;
int atomnum1, atomnum2, atomnum;
FILE *fpout;
double rms;
double CG0[3] = {0,0,0}; 
double CG1[3]={0,0,0}; 
double ROT[3][3]={1,0,0,0,1,0,0,0,1};
char tmpchar1[MAXCHAR];
char tmpchar2[MAXCHAR];

/* PDB */
int rpdb(char *filename, int *atomnum, ATOM * atom) {
	int numatom;
	int tmpint;
	char line[MAXCHAR];
	double x, y, z;
	FILE *fpin;

	if ((fpin = fopen(filename, "r")) == NULL) {
		fprintf(stdout, "Cannot open the pdb file %s in rpdb(), exit\n", filename);
		exit(1);
	}
	numatom = 0;
	for (;;) {
		if (fgets(line, MAXCHAR, fpin) == NULL) break;
		if (strncmp("ATOM", line, 4) == 0 || strncmp("HETATM", line, 6) == 0) {
			sscanf(&line[22], "%d%lf%lf%lf", &tmpint, &x, &y, &z);
			atom[numatom].name[0] = line[12];
			atom[numatom].name[1] = line[13];
			atom[numatom].name[2] = line[14];
			atom[numatom].name[3] = line[15];
			atom[numatom].x = x;
			atom[numatom].y = y;
			atom[numatom].z = z;
			numatom++;
		}
	}
	*atomnum = numatom;
	fclose(fpin);
}

int wpdb(char *filename1, char *filename2, ATOM *atom) {
int numatom;
int tmpint;
double x, y, z;
FILE *fpin, *fpout;

if ((fpin = fopen(filename1, "r")) == NULL) {
	fprintf(stdout, "Cannot open file %s to read in wpdb(), exit\n", filename1);
	exit(1);
}
if ((fpout = fopen(filename2, "w")) == NULL) {
	fprintf(stdout, "Cannot open file %s to write in wpdb(), exit\n", filename2);
	exit(1);
}
numatom = 0;
for (;;) {
	if (fgets(line, MAXCHAR, fpin) == NULL) break;
	if (strncmp("ATOM", line, 4) == 0 || strncmp("HETATM", line, 6) == 0) {
		strcpy(tmpchar1, "");
		tmpchar1[0]='\0';
		for(i=0;i<30;i++)
			tmpchar1[i] = line[i];			
		tmpchar1[30] = '\0'; 
		tmpint = strlen(line);
		strcpy(tmpchar2, "");
		tmpchar2[0]='\0';
		for(i=54;i<tmpint-1;i++)
			tmpchar2[i-54] = line[i];			
		tmpchar2[tmpint-1] = '\n';
		fprintf(fpout, "%s%8.3lf%8.3lf%8.3lf%s\n", tmpchar1, atom[numatom].x, atom[numatom].y, atom[numatom].z, tmpchar2);	
		numatom++;
	} 
 	else
	fprintf(fpout, "%s",line);
}
fclose(fpin);
fclose(fpout);
}

int rmatrix(char *filename) {
	char line[MAXCHAR];
	FILE *fpin;

	if ((fpin = fopen(filename, "r")) == NULL) {
		fprintf(stdout, "Cannot open the matrix file %s, exit\n", filename);
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

double lsfit(ATOM atm1[], ATOM atm2[], ATOM atm3[],int n) {
int i,j,k;
int ict;
int ix,iy,iz;
int iflag;
/* ict: counter for number of iterations. must do at least 3 iterations.
ix, iy,iz: pointer used in iterative least squares. may have the value 1,2,3
the value changes on each iteration. */
double tol=0.00001;
double sig,gam, sg, bb,cc;
double x0[3][MAXATOM];
double x1[3][MAXATOM];
double xyzfit[3][MAXATOM];
double aa[3][3];
double rot[3][3];
double cg0[3];
double cg1[3];
double an;
double rms,ssum;
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
  cg0[i]=0.;
  cg1[i]=0.;
  }
  for(i=0;i<n;i++)
   for(j=0;j<3;j++)
    {
    cg0[j]=cg0[j]+x0[j][i];
    cg1[j]=cg1[j]+x1[j][i];
   }
   for(i=0;i<3;i++)
   {
     cg0[i]*=an;
     cg1[i]*=an;
   }
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
     aa[i][j]=0.0;
  for(k=0;k<n;k++)
  {
   aa[0][0]=aa[0][0]+(x1[0][k]-cg1[0])*(x0[0][k]-cg0[0]);
   aa[1][0]=aa[1][0]+(x1[1][k]-cg1[1])*(x0[0][k]-cg0[0]);
   aa[2][0]=aa[2][0]+(x1[2][k]-cg1[2])*(x0[0][k]-cg0[0]);
   aa[0][1]=aa[0][1]+(x1[0][k]-cg1[0])*(x0[1][k]-cg0[1]);
   aa[1][1]=aa[1][1]+(x1[1][k]-cg1[1])*(x0[1][k]-cg0[1]);
   aa[2][1]=aa[2][1]+(x1[2][k]-cg1[2])*(x0[1][k]-cg0[1]);
   aa[0][2]=aa[0][2]+(x1[0][k]-cg1[0])*(x0[2][k]-cg0[2]);
   aa[1][2]=aa[1][2]+(x1[1][k]-cg1[1])*(x0[2][k]-cg0[2]);
   aa[2][2]=aa[2][2]+(x1[2][k]-cg1[2])*(x0[2][k]-cg0[2]);
  }
  for(i=0;i<3;i++)
    {
    for(j=0;j<3;j++)
      rot[i][j]=0.;
    rot[i][i]=1.0;
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
   bb=gam*rot[iy-1][k]+sig*rot[iz-1][k];
   cc=gam*rot[iz-1][k]-sig*rot[iy-1][k];
   rot[iy-1][k]=bb*sg;
   rot[iz-1][k]=cc*sg;
  }
  iflag=1;
  goto a50;
  a70:
 /* the following code translates the ligand the center of mass of the
 receptor site and rotates it according to the rotation matrxouneror ix
 calculated by orient */
 for(i=0;i<n;i++)
 {
 xyzfit[0][i]=cg0[0]+rot[0][0]*(x1[0][i]-cg1[0])+ \
             rot[0][1]*(x1[1][i]-cg1[1])+
             rot[0][2]*(x1[2][i]-cg1[2]);
 xyzfit[1][i]=cg0[1]+rot[1][0]*(x1[0][i]-cg1[0])+ \
             rot[1][1]*(x1[1][i]-cg1[1])+
             rot[1][2]*(x1[2][i]-cg1[2]);
 xyzfit[2][i]=cg0[2]+rot[2][0]*(x1[0][i]-cg1[0])+ \
             rot[2][1]*(x1[1][i]-cg1[1])+
             rot[2][2]*(x1[2][i]-cg1[2]);
 }
 fprintf(fpout, "CG0\t%16.10lf%16.10lf%16.10lf\n", cg0[0], cg0[1], cg0[2]);
 fprintf(fpout, "CG1\t%16.10lf%16.10lf%16.10lf\n", cg1[0], cg1[1], cg1[2]);
 fprintf(fpout, "ROT0\t%16.10lf%16.10lf%16.10lf\n", rot[0][0], rot[0][1], rot[0][2]);
 fprintf(fpout, "ROT1\t%16.10lf%16.10lf%16.10lf\n", rot[1][0], rot[1][1], rot[1][2]);
 fprintf(fpout, "ROT2\t%16.10lf%16.10lf%16.10lf\n", rot[2][0], rot[2][1], rot[2][2]);

for(j=0;j<n;j++)
  for(i=0;i<3;i++)
 ssum+=(xyzfit[i][j]-x0[i][j])*(xyzfit[i][j]-x0[i][j]);
 rms=ssum/n;
 rms=sqrt(rms);
 for(i=0;i<n;i++)
 {
  atm3[i].x=xyzfit[0][i];
  atm3[i].y=xyzfit[1][i];
  atm3[i].z=xyzfit[2][i];
 }
return rms;
}



int main(int argc, char *argv[]) {
	int i;
	int index;
	int format = 0;
	int judge_flag = 0;
	char command_at[2 * MAXCHAR];
	char command_bt[2 * MAXCHAR];
	double sum = 0.0;
	double coordx, coordy, coordz;
	if (strcmp(COLORTEXT, "YES") == 0 || strcmp(COLORTEXT, "yes") == 0) {
		if (argc == 2
			&& (strcmp(argv[1], "-h") == 0
				|| strcmp(argv[1], "-H") == 0)) {
			printf
				("[31mUsage: match  -i[0m input file name in pdb format \n"
				 "[31m              -r[0m ref file name in pdb format\n"
				 "[31m              -o[0m output file name in pdb format\n"
				 "[31m              -m[0m matrix file\n"
				 "[31m              -t[0m type:\n"
				 "[31m                [0m 	1: lsfit, need -i, -r, -o, -m;\n"    
				 "[31m                [0m	2: translation, need -i, -o, -m;\n"); 
			exit(1);
		}
		if (argc != 9 && argc != 11) {
			printf
				("[31mUsage: match  -i[0m input file name in pdb format \n"
				 "[31m              -r[0m ref file name in pdb format\n"
				 "[31m              -o[0m output file name in pdb format\n"
				 "[31m              -m[0m matrix file\n"
				 "[31m              -t[0m type:\n"
				 "[31m                [0m 	1: lsfit, need -i, -r, -o, -m;\n"    
				 "[31m                [0m	2: translation, need -i, -o, -m;\n"); 
			exit(1);
		}
	} else {
		if (argc == 2
			&& (strcmp(argv[1], "-h") == 0
				|| strcmp(argv[1], "-H") == 0)) {
			printf("Usage: match  -i input file name in pdb format \n"
				   "              -r ref file name in pdb format\n"
				   "              -o output file name in pdb format\n"
				   "              -m matrix file \n"
				   "              -t type \n"
				   "                 	1: lsfit, need -i, -r, -o, -m;\n"
				   "                 	2: translation, need -i, -o, -m;\n");
			exit(1);
		}
		if (argc != 9 && argc != 11) {
			printf("Usage: match  -i input file name in pdb format \n"
				   "              -r ref file name in pdb format\n"
				   "              -o output file name in pdb format\n"
				   "              -m matrix file \n"
				   "              -t type \n"
				   "                 	1: lsfit, need -i, -r, -o, -m;\n"
				   "                 	2: translation, need -i, -o, -m;\n");
			exit(1);
		}

	}

	for (i = 1; i < argc; i += 2) {
		if (strcmp(argv[i], "-i") == 0) {
			strcpy(ifilename, argv[i + 1]);
			numarg++;
		}
		if (strcmp(argv[i], "-o") == 0) {
			strcpy(ofilename, argv[i + 1]);
			numarg++;
		}
		if (strcmp(argv[i], "-r") == 0) {
			strcpy(rfilename, argv[i + 1]);
			numarg++;
		}
		if (strcmp(argv[i], "-m") == 0) {
			strcpy(mfilename, argv[i + 1]);
			numarg++;
		}
		if (strcmp(argv[i], "-t") == 0)
			jobtype = atoi(argv[i+1]);
	}
if(jobtype == 1 && numarg != 4) {
	printf("The number of arguments is wrong, exit\n");
	exit(1);
}
if(jobtype == 2 && numarg != 3) {
	printf("The number of arguments is wrong, exit\n");
	exit(1);
}
if(jobtype == 1) {
	rpdb(ifilename, &atomnum1, atom1);
	rpdb(rfilename, &atomnum2, atom2);
	if(atomnum1 != atomnum2 || atomnum1 <=0) {
		printf("\nThe numumber of atoms in input file (%d) is different from the number of atoms in ref file (%d)", atomnum1, atomnum2);
		exit(1);
	}
	atomnum = atomnum1;
	for(i=0;i<atomnum;i++) {
		sum+=(atom1[i].x-atom2[i].x) * (atom1[i].x-atom2[i].x); 
		sum+=(atom1[i].y-atom2[i].y) * (atom1[i].y-atom2[i].y); 
		sum+=(atom1[i].z-atom2[i].z) * (atom1[i].z-atom2[i].z); 
	}
	rms = sqrt(sum/atomnum);
	printf("\nThe rmsd before least-square fitting is %9.4lf", rms);
	if ((fpout = fopen(mfilename, "w")) == NULL) {
		fprintf(stdout, "Cannot open file %s to write in main(), exit\n", mfilename);
		exit(1);
	}
printf("\nwjm %5d\n", atomnum);
for(i=0;i<atomnum;i++) 
	printf("ATOM %5d %5s %9.2lf %9.2lf %9.2lf\n", i+1, atom1[i].name, atom1[i].x, atom1[i].y, atom1[i].z);

for(i=0;i<atomnum;i++) 
	printf("ATOM %5d %5s %9.2lf %9.2lf %9.2lf\n", i+1, atom2[i].name, atom2[i].x, atom2[i].y, atom2[i].z);

for(i=0;i<atomnum;i++) {
	sum = (atom1[i].x-atom2[i].x) * (atom1[i].x-atom2[i].x);
	sum+=(atom1[i].y-atom2[i].y) * (atom1[i].y-atom2[i].y);
	sum+=(atom1[i].z-atom2[i].z) * (atom1[i].z-atom2[i].z);
	sum = sqrt(sum);
	printf("ATOM %5d %5s %9.2lf\n", i+1, atom2[i].name, sum);
}
	rms = lsfit(atom2, atom1, atom3, atomnum) ;
	printf("\nThe rmsd after least-square fitting is %9.4lf\n", rms);
	wpdb(ifilename, ofilename, atom3);
}
if(jobtype == 2) {
	rpdb(ifilename, &atomnum1, atom1);
	rmatrix(mfilename);
	for(i=0;i<atomnum1;i++) {
		coordx = atom1[i].x;
		coordy = atom1[i].y;
		coordz = atom1[i].z;
		atom1[i].x = CG0[0] + ROT[0][0] * (coordx - CG1[0]) + 
				      ROT[0][1] * (coordy - CG1[1]) +
				      ROT[0][2] * (coordz - CG1[2]) ;
		atom1[i].y = CG0[1] + ROT[1][0] * (coordx - CG1[0]) + 
				      ROT[1][1] * (coordy - CG1[1]) +
				      ROT[1][2] * (coordz - CG1[2]) ;
		atom1[i].z = CG0[2] + ROT[2][0] * (coordx - CG1[0]) + 
				      ROT[2][1] * (coordy - CG1[1]) +
				      ROT[2][2] * (coordz - CG1[2]) ;
	}
	wpdb(ifilename, ofilename, atom1);
}
return 0;
}
