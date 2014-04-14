# include <stdio.h>
# include <math.h>
# include <ctype.h>
# include <stdlib.h>
# include <string.h>
# define MAX_AT_CORR 200
# define MAX_BT_CORR 1000
# define MAXCHAR 256 
# define COLORTEXT "YES" 
typedef struct {
        char str[10];
} STR_ARRAY;

char *amberhome;
long i,j,k=0;
long atomnum = 0;
long top_atomnum = 0;
long top_resnum = 0;
long top_bondnum = 0;
long top_bondnum_h = 0;
long top_bondnum_noh = 0;
long crd_atomnum = 0;
long bondnum = 0;
long bondnum_h = 0;
long bondnum_noh = 0;
long resnum = 0;
long resnum_water = 0;
int corratnum = 0;
int corrbtnum = 0;
int tmpint1, tmpint2;

long *bondat1;
long *bondat2;
long *bondat1_h;
long *bondat2_h;
long *bondat1_noh;
long *bondat2_noh;
long *resno;
long *reshead;

double *charge;
double *coord_x; 
double *coord_y; 
double *coord_z; 

STR_ARRAY *res;
STR_ARRAY *atomname;
STR_ARRAY *atomtype;
STR_ARRAY *atomtypebak;
STR_ARRAY *atomres;
STR_ARRAY *bondtype_h;
STR_ARRAY *bondtype_noh;

char atcorr[MAX_AT_CORR][10];
char atorg[MAX_AT_CORR][10];
char bt[MAX_BT_CORR][10];
char btan1[MAX_BT_CORR][10];
char btan2[MAX_BT_CORR][10];

int output_atomtype_flag = 2;
int output_bondtype_flag = 2;
int water_flag = 0;

char pfilename[MAXCHAR];
char cfilename[MAXCHAR];
char ofilename[MAXCHAR];
char atfilename[MAXCHAR];
char btfilename[MAXCHAR];

double tmpf1, tmpf2, tmpf3, tmpf4;
char line[MAXCHAR];
char tmpchar1[MAXCHAR], tmpchar2[MAXCHAR], tmpchar3[MAXCHAR];

FILE *fpinp, *fpinc, *fpout ;
FILE *fpat, *fpbt;

void ratcorr(void) {
	corratnum = 0;
        for (;;) {
                if (fgets(line, MAXCHAR, fpat) == NULL) break;
		sscanf(line, "%s%s", atorg[corratnum], atcorr[corratnum]);
		corratnum++;
		if(corratnum > MAX_AT_CORR) {
			fprintf(stdout, "Error: the nubmer of corresponding atom type pairs (%d) exceed MAX_AT_CORR (%d), extends MAX_AT_CORR defined in top2mol2.c and recompile the program", corratnum, MAX_AT_CORR);
			exit(1);
		}
	}
}
void rbtcorr(void) {
	corrbtnum = 0;
        for (;;) {
                if (fgets(line, MAXCHAR, fpbt) == NULL) break;
		sscanf(line, "%s%s%s", btan1[corrbtnum], btan2[corrbtnum], bt[corrbtnum]);
		corrbtnum++;
		if(corrbtnum > MAX_BT_CORR) {
			fprintf(stdout, "Error: the nubmer of corresponding bond type pairs (%d) exceed MAX_BT_CORR (%d), extends MAX_BT_CORR defined in top2mol2.c and recompile the program", corrbtnum, MAX_BT_CORR);
			exit(1);
		}
	}
}

void rcrd() {
long i, j;
i = 0;
j = 0;
for (;;) {
	if (fgets(line, MAXCHAR, fpinc) == NULL) break;
        j++;
        if (j == 2) {
        	sscanf(line, "%ld", &crd_atomnum);
		if(crd_atomnum != top_atomnum) {
			printf("\nError: the atom number from the topology (%ld) and the crd (%ld) files is not same\n", top_atomnum, crd_atomnum);
			exit(1);
		}
	}
        if (strncmp(".", &line[4], 1) == 0
                && strncmp(".", &line[16], 1) == 0
                && strncmp(".", &line[28], 1) == 0) {
        	sscanf(&line[1], "%lf%lf%lf%lf%lf%lf", 
		&coord_x[i], &coord_y[i], &coord_z[i], 
		&coord_x[i+1], &coord_y[i+1], &coord_z[i+1]);
                i = i + 2;
                if (i >= (crd_atomnum + 10)) return;
        }
}
}

void rtop(int mode) {
long tmpint;
int flag = -999;

for(;;) {
	line[7] = ' ';
	line[15] = ' ';
	line[23] = ' ';
	line[31] = ' ';
	line[39] = ' ';
	line[47] = ' ';
	line[55] = ' ';
	line[63] = ' ';
	line[71] = ' ';
	line[79] = ' ';
	if(fgets(line, MAXCHAR, fpinp)==NULL) break;
    	sscanf(line, "%s%s", tmpchar1, tmpchar2);
    	if (strcmp("%FLAG", tmpchar1) == 0 &&(flag == 1 || flag == 3 || flag == 5 || flag ==7 ||
					      flag == 9 || flag == 11|| flag == 13|| flag ==15)) 
		flag = -999;
    	if (strcmp("%FLAG", tmpchar1) == 0 &&strcmp("POINTERS", tmpchar2) == 0) {
		flag = 0; 
		continue; 
    	} 
    	if (strcmp("%FLAG", tmpchar1) == 0 &&strcmp("ATOM_NAME", tmpchar2) == 0) {
		flag = 2; 
		continue; 
    	} 
    	if (strcmp("%FLAG", tmpchar1) == 0 &&strcmp("CHARGE", tmpchar2) == 0) {
		flag = 4; 
		continue; 
    	} 
    	if (strcmp("%FLAG", tmpchar1) == 0 &&strcmp("RESIDUE_LABEL", tmpchar2) == 0) {
		flag = 6; 
		continue; 
    	} 
    	if (strcmp("%FLAG", tmpchar1) == 0 &&strcmp("RESIDUE_POINTER", tmpchar2) == 0) {
		flag = 8; 
		continue; 
    	} 
    	if (strcmp("%FLAG", tmpchar1) == 0 &&strcmp("BONDS_INC_HYDROGEN", tmpchar2) == 0) {
		flag = 10; 
		continue; 
    	} 
    	if (strcmp("%FLAG", tmpchar1) == 0 &&strcmp("BONDS_WITHOUT_HYDROGEN", tmpchar2) == 0) {
		flag = 12; 
		continue; 
    	} 
    	if (strcmp("%FLAG", tmpchar1) == 0 &&strcmp("AMBER_ATOM_TYPE", tmpchar2) == 0) {
		flag = 14; 
		continue; 
    	} 

    	if(flag==0) {flag = 1 ; continue; }
    	if(flag==2) {flag = 3; i =0 ; continue; }
    	if(flag==4) {flag = 5; i =0 ; continue; }
    	if(flag==6) {flag = 7; i =0 ; continue; }
    	if(flag==8) {flag = 9; i =0 ; continue; }
    	if(flag==10) {flag = 11; i =0 ; k = 0; continue; }
    	if(flag==12) {flag = 13; i =0 ; k = 0; continue; }
    	if(flag==14) {flag = 15; i =0 ; continue; }

    	if(mode == 0 && flag ==1) {
		sscanf(line, "%ld%ld%ld%ld", &top_atomnum, &tmpint, &top_bondnum_h, &top_bondnum_noh);
		top_bondnum = top_bondnum_h + top_bondnum_noh;
    		fgets(line, MAXCHAR, fpinp);
		sscanf(line, "%ld%ld", &tmpint, &top_resnum);
		return;
    	}	
    	if(flag ==3) {
		for(j=0;j<20;j++) {
          		if(line[j*4]==' '||line[j*4]=='\0'||j*4>=strlen(line)-2) break;
	  		atomname[i+j].str[0] = line[j*4]; 
	  		if(line[j*4+1]!=' ') 
				atomname[i+j].str[1] = line[j*4+1];
	  		else {
				atomname[i+j].str[1] = '\0';
				continue;
	  		}
	  		if(line[j*4+2]!=' ') 
				atomname[i+j].str[2] = line[j*4+2];
	  		else {
				atomname[i+j].str[2] = '\0';
				continue;
	  		}
	  		if(line[j*4+3]!=' ') 
				atomname[i+j].str[3] = line[j*4+3];
	  		else {
				atomname[i+j].str[3] = '\0';
				continue;
	  		}
       		}
        	i+=j;
        	atomnum = i;	
    	}
    	if(flag ==5) {
        	sscanf(line, "%lf%lf%lf%lf%lf", &charge[i], &charge[i+1], &charge[i+2], &charge[i+3],&charge[i+4]);
        	charge[i]/=18.2223;
        	charge[i+1]/=18.2223;
        	charge[i+2]/=18.2223;
        	charge[i+3]/=18.2223;
        	charge[i+4]/=18.2223;
		i+=5;
    	}
    	if(flag ==7) {
		for(j=0;j<20;j++) {
          		if(line[j*4]==' '||line[j*4]=='\0'||j*4>=strlen(line)-2) break;
	  		res[i+j].str[0] = line[j*4]; 
	  		if(line[j*4+1]!=' ') 
				res[i+j].str[1] = line[j*4+1];
	  		else {
				res[i+j].str[1] = '\0';
				continue;
	  		}
	  		if(line[j*4+2]!=' ') 
				res[i+j].str[2] = line[j*4+2];
	  		else {
				res[i+j].str[2] = '\0';
				continue;
	  		}
	  		if(line[j*4+3]!=' ') 
				res[i+j].str[3] = line[j*4+3];
	  		else {
				res[i+j].str[3] = '\0';
				continue;
	  		}
        	}
        	i+=j;
        	resnum = i;	
    	}
    	if(flag ==9) {
        	sscanf(line, "%ld%ld%ld%ld%ld%ld%ld%ld%ld%ld", &reshead[i], &reshead[i+1], &reshead[i+2],
		&reshead[i+3],&reshead[i+4], &reshead[i+5], &reshead[i+6], &reshead[i+7],
		&reshead[i+8],&reshead[i+9]);
		i+=10;
    	}
    	if(flag ==11) {
                j = 0;
                while(j <= 10) {
                        j++;
                        if(line[j*8-1] <='9' && line[j*8-1]>='0') 
                                sscanf(&line[8*j-8], "%ld", &tmpint);
                        else break;
                        if(k==0) {
                                bondat1_h[bondnum_h] = tmpint;
                                k=1;
                                continue;
                        }
                        if(k==1) {
                                bondat2_h[bondnum_h] = tmpint;
                                k=2;
                                continue;
                        }
                        if(k==2) {
                                k=0;
                                bondnum_h++;
                                continue;
                        }
                }
    	}
    	if(flag ==13) {
		j = 0;
		while(j <= 10) {
			j++;
                	if(line[j*8-1] <='9' && line[j*8-1]>='0') 
                        	sscanf(&line[8*j-8], "%ld", &tmpint);
			else break;
			if(k==0) {
				bondat1_noh[bondnum_noh] = tmpint;
				k=1;
				continue;
			}
			if(k==1) {
				bondat2_noh[bondnum_noh] = tmpint;
				k=2;
				continue;
			}
			if(k==2) {
				k=0;
				bondnum_noh++;
				continue;
			}
		}
    	}
    	if(flag ==15) {
		for(j=0;j<20;j++) {
          	if(line[j*4]==' '||line[j*4]=='\0'||j*4>=strlen(line)-2) break;
	  	atomtype[i+j].str[0] = line[j*4]; 
	  	if(line[j*4+1]!=' ') 
			atomtype[i+j].str[1] = line[j*4+1];
	  	else {
			atomtype[i+j].str[1] = '\0';
			continue;
	  	}
	  	if(line[j*4+2]!=' ') 
			atomtype[i+j].str[2] = line[j*4+2];
	  	else {
			atomtype[i+j].str[2] = '\0';
			continue;
	  	}
	  	if(line[j*4+3]!=' ') 
			atomtype[i+j].str[3] = line[j*4+3];
	  	else {
			atomtype[i+j].str[3] = '\0';
			continue;
	  	}
        }
        i+=j;
    	}
}
}

static void write(void) {
char tmpc1[MAXCHAR];	
char tmpc2[MAXCHAR];	
long i;
long tmpint;
fprintf(fpout, "@<TRIPOS>MOLECULE\n");
fprintf(fpout, "%s\n", "MOL");
if(water_flag == 0)
	fprintf(fpout, "%5ld%6ld%6ld%6d%6d\n", atomnum-3*resnum_water, bondnum-3*resnum_water, resnum-resnum_water, 0, 0);
else
	fprintf(fpout, "%5ld%6ld%6ld%6d%6d\n", atomnum, bondnum, resnum, 0, 0);
fprintf(fpout, "PROTEIN\n");
fprintf(fpout, "%s\n\n\n", "AMBER_CHARGES");
fprintf(fpout, "%s\n%s\n", "@<TRIPOS>", "PROTEIN protein");
fprintf(fpout, "@<TRIPOS>ATOM\n");

tmpint = 0;

for (i = 0; i < atomnum; i++) {
	if(water_flag == 0 && strcmp(atomres[i].str, "WAT") ==0) continue;
	tmpint++;
	strcpy(tmpc1, "");
	strcpy(tmpc2, "");
	sprintf(tmpc1, "%ld", resno[i]);  
	/*newitoa(resno[i], tmpc1);*/
	strcpy(tmpc2, atomres[i].str);
	strcat(tmpc2, tmpc1);
	fprintf(fpout, "%7ld %-8s%10.4lf%10.4lf%10.4lf %-6s%5ld %-8s%10.4lf\n", tmpint, atomname[i].str, 
		coord_x[i], coord_y[i], coord_z[i], atomtype[i].str, resno[i],tmpc2, charge[i]);
}
fprintf(fpout, "@<TRIPOS>BOND\n");
tmpint = 0;
for (i = 0; i < bondnum_noh; i++) {
	tmpint ++;
	fprintf(fpout, "%6ld %5ld %5ld %-4s\n", tmpint, bondat1_noh[i]+1, bondat2_noh[i]+1, bondtype_noh[i].str);
}

for (i = 0; i < bondnum_h; i++) {
	if(water_flag == 0) {
		if(strcmp(atomres[bondat1_h[i]].str, "WAT") == 0) continue;
		if(strcmp(atomres[bondat2_h[i]].str, "WAT") == 0) continue;
	}
	tmpint ++;
	fprintf(fpout, "%6ld %5ld %5ld %-4s\n", tmpint, bondat1_h[i]+1, bondat2_h[i]+1, bondtype_h[i].str);
}

fprintf(fpout, "@<TRIPOS>SUBSTRUCTURE\n");
if(water_flag == 0)
	for (i = 0; i < resnum-resnum_water; i++) {
		strcpy(tmpc2, "");
		sprintf(tmpc1, "%ld", i+1);  
		/* newitoa(i+1, tmpc1);	*/
		strcpy(tmpc2, res[i].str);
		strcat(tmpc2, tmpc1);
		fprintf(fpout, "%6ld %-7s %6ld RESIDUE           1 A\n", i + 1, tmpc2, reshead[i]);
	}
if(water_flag == 1)
	for (i = 0; i < resnum; i++) {
		strcpy(tmpc2, "");
		/* newitoa(i+1, tmpc1);	*/
		sprintf(tmpc1, "%ld", i+1);
		strcpy(tmpc2, res[i].str);
		strcat(tmpc2, tmpc1);
		fprintf(fpout, "%6ld %-7s %6ld RESIDUE           1 A\n", i + 1, tmpc2, reshead[i]);
	}
}


int main(int argc, char *argv[]) {
int i;
/*  readin a topology file and crd file,then save out a mol2 file */

amberhome = (char *) getenv("AMBERHOME");
if( amberhome == NULL ){
   fprintf( stdout, "AMBERHOME is not set!\n" );
   exit(1);
}

if (strcmp(COLORTEXT, "YES") == 0 || strcmp(COLORTEXT, "yes") == 0) {
	if (argc == 2
		&& (strcmp(argv[1], "-h") == 0
			|| strcmp(argv[1], "-H") == 0)) {
		printf("[31mUsage: top2mol2 -p [0m topology file name\n"
			   "[31m                -c [0m rst or crd file name\n"
			   "[31m                -o [0m output file name)\n"
			   "[31m                -ac[0m atom type corresponding file (optional)\n"
			   "[31m                -bc[0m bond type corresponding file (optional)\n"
			   "[31m                -at[0m atom type: sybyl (the default) or amber, optional\n"
			   "[31m                -bt[0m bond type: sybyl (the default) or amber (all set to 1), optional\n"
			   "[31m                -wt[0m keep water flag: 1 (inclduing water) or 0 (removing water, the default), optional\n");
		exit(1);
	}
	if (argc != 7 && argc != 9 && argc != 11 && argc != 13 && argc != 15 && argc != 17) {
		printf("[31mUsage: top2mol2 -p [0m topology file name\n"
			   "[31m                -c [0m rst or crd file name\n"
			   "[31m                -o [0m output file name)\n"
			   "[31m                -ac[0m atom type corresponding file (optional)\n"
			   "[31m                -bc[0m bond type corresponding file (optional)\n"
			   "[31m                -at[0m atom type: sybyl (the default) or amber, optional\n"
			   "[31m                -bt[0m bond type: sybyl (the default) or amber (all set to 1), optional\n"
			   "[31m                -wt[0m keep water flag: 1 (inclduing water) or 0 (removing water, the default), optional\n");
		exit(1);
	}
} else {
	if (argc == 2)
	if (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "-H") == 0) {
		printf("\n Usage top2mol2 -p  topology file name");
		printf("\n                -c  rst or crd file name");
		printf("\n                -o  output file name");
		printf("\n                -ac atom type corresponding file (optional)\n ");
		printf("\n                -bc bond type corresponding file (optional)\n ");
		printf("\n                -at atom type: sybyl (the default) or amber, optional");
		printf("\n                -bt bond type: sybyl (the default) or amber (all set to 1, optional)");
		printf("\n                -wt keep water flag: 1 (including water) or 0 (removing water, the default), optional)");
		exit(1);
	}
	if (argc != 7 && argc != 9 && argc != 11 && argc != 13 && argc != 15 && argc != 17) {
		printf("\n Usage top2mol2 -p  topology file name");
		printf("\n                -c  rst or crd file name");
		printf("\n                -o  output file name");
		printf("\n                -ac atom type corresponding file (optional)\n ");
		printf("\n                -bc bond type corresponding file (optional)\n ");
		printf("\n                -at atom type: sybyl (the default) or amber, optional");
		printf("\n                -bt bond type: sybyl (the default) or amber (all set to 1, optional)");
		printf("\n                -wt keep water flag: 1 (including water) or 0 (removing water, the default), optional)");
		exit(1);
	}
}

strcpy(atfilename, amberhome);
strcat(atfilename, "/dat/antechamber/ATOMTYPE_CHECK.TAB");
strcpy(btfilename, amberhome);
strcat(btfilename, "/dat/antechamber/BONDTYPE_CHECK.TAB");

for (i = 1; i < argc; i += 2) {
	if (strcmp(argv[i], "-p") == 0) 
		strcpy(pfilename, argv[i + 1]);
	if (strcmp(argv[i], "-c") == 0) 
		strcpy(cfilename, argv[i + 1]);
	if (strcmp(argv[i], "-o") == 0) 
		strcpy(ofilename, argv[i + 1]);
	if (strcmp(argv[i], "-ac") == 0) 
		strcpy(atfilename, argv[i + 1]);
	if (strcmp(argv[i], "-bc") == 0) 
		strcpy(btfilename, argv[i + 1]);
	if (strcmp(argv[i], "-at") == 0) {
		if(strcmp("Amber", argv[i+1]) == 0 || strcmp("AMBER", argv[i+1]) == 0 || 
		strcmp("amber", argv[i+1]) == 0 || strcmp("a", argv[i+1]) == 0) 
			output_atomtype_flag = 1;
		if(strcmp("Sybyl", argv[i+1]) == 0 || strcmp("SYBYL", argv[i+1]) == 0 || 
		strcmp("sybyl", argv[i+1]) == 0 || strcmp("s", argv[i+1]) == 0) 
			output_atomtype_flag = 2;
	}
	if (strcmp(argv[i], "-bt") == 0) {
		if(strcmp("Amber", argv[i+1]) == 0 || strcmp("AMBER", argv[i+1]) == 0 || 
		strcmp("amber", argv[i+1]) == 0 || strcmp("a", argv[i+1]) == 0) 
			output_bondtype_flag = 1;
		if(strcmp("Sybyl", argv[i+1]) == 0 || strcmp("SYBYL", argv[i+1]) == 0 || 
		strcmp("sybyl", argv[i+1]) == 0 || strcmp("s", argv[i+1]) == 0) 
			output_bondtype_flag = 2;
	}
	if (strcmp(argv[i], "-wt") == 0) 
		water_flag = atoi(argv[i+1]);
}
if((fpinp=fopen(pfilename,"r"))==NULL) {
    	fprintf(stdout, "\n Cannot open topology file: %s, exit", pfilename);
	exit(1);
}
if((fpinc=fopen(cfilename,"r"))==NULL) {
    	fprintf(stdout, "\n Cannot open coordinate file %s, exit", cfilename);
	exit(1);
}
if((fpout=fopen(ofilename,"w"))==NULL) {
	fprintf(stdout, "\n Cannot open file to write: %s, exit", ofilename);
	exit(1);
}
if(output_atomtype_flag == 2)
	if((fpat=fopen(atfilename,"r"))==NULL) {
    		printf("\n Cannot open atom type corresponding file %s, exit\n", atfilename);
		exit(1);
	}
if(output_bondtype_flag == 2)
	if((fpbt=fopen(btfilename,"r"))==NULL) {
    		printf("\n Cannot open bond type corresponding file %s, exit\n", btfilename);
		exit(1);
	}

rtop(0);
rewind(fpinp); 

bondat1_h = (long *) malloc(sizeof(long) * (top_bondnum_h + 10));
if (bondat1_h == NULL) {
	fprintf(stdout, "memory allocation error for *bondat1_h\n");
	exit(1);
 }

bondat2_h = (long *) malloc(sizeof(long) * (top_bondnum_h + 10));
if (bondat2_h == NULL) {
	fprintf(stdout, "memory allocation error for *bondat2_h\n");
	exit(1);
 }

bondat1_noh = (long *) malloc(sizeof(long) * (top_bondnum_noh + 10));
if (bondat1_noh == NULL) {
        fprintf(stdout, "memory allocation error for *bondat1_noh\n");
        exit(1);
 }
                                                                                                                                                                            
bondat2_noh = (long *) malloc(sizeof(long) * (top_bondnum_noh + 10));
if (bondat2_noh == NULL) {
        fprintf(stdout, "memory allocation error for *bondat2_noh\n");
        exit(1);
 }

resno = (long *) malloc(sizeof(long) * (top_atomnum + 10));
if (resno == NULL) {
	fprintf(stdout, "memory allocation error for *resno\n");
	exit(1);
 }

reshead = (long *) malloc(sizeof(long) * (top_resnum + 10));
if (reshead == NULL) {
	fprintf(stdout, "memory allocation error for *reshead\n");
	exit(1);
 }

charge = (double *) malloc(sizeof(double) * (top_atomnum + 10));
if (charge == NULL) {
	fprintf(stdout, "memory allocation error for *charge\n");
	exit(1);
 }

coord_x = (double *) malloc(sizeof(double) * (top_atomnum + 10));
if (coord_x == NULL) {
	fprintf(stdout, "memory allocation error for *coord_x\n");
	exit(1);
 }

coord_y = (double *) malloc(sizeof(double) * (top_atomnum + 10));
if (coord_y == NULL) {
	fprintf(stdout, "memory allocation error for *coord_y\n");
	exit(1);
 }

coord_z = (double *) malloc(sizeof(double) * (top_atomnum + 10));
if (coord_z == NULL) {
	fprintf(stdout, "memory allocation error for *coord_z\n");
	exit(1);
 }

res = (STR_ARRAY *) malloc(sizeof(STR_ARRAY) * (top_resnum + 10));
if (res == NULL) {
	fprintf(stdout, "memory allocation error for *res\n");
	exit(1);
 }

bondtype_h = (STR_ARRAY *) malloc(sizeof(STR_ARRAY) * (top_bondnum_h + 10));
if (bondtype_h == NULL) {
	fprintf(stdout, "memory allocation error for *bondtype_h\n");
	exit(1);
 }

bondtype_noh = (STR_ARRAY *) malloc(sizeof(STR_ARRAY) * (top_bondnum_noh + 10));
if (bondtype_noh == NULL) {
	fprintf(stdout, "memory allocation error for *bondtype_noh\n");
	exit(1);
 }

atomname = (STR_ARRAY *) malloc(sizeof(STR_ARRAY) * (top_atomnum + 10));
if (atomname == NULL) {
	fprintf(stdout, "memory allocation error for *atomname\n");
	exit(1);
 }

atomtype = (STR_ARRAY *) malloc(sizeof(STR_ARRAY) * (top_atomnum + 10));
if (atomtype == NULL) {
	fprintf(stdout, "memory allocation error for *atomtype\n");
	exit(1);
 }

atomtypebak = (STR_ARRAY *) malloc(sizeof(STR_ARRAY) * (top_atomnum + 10));
if (atomtypebak == NULL) {
	fprintf(stdout, "memory allocation error for *atomtypebak\n");
	exit(1);
 }

atomres = (STR_ARRAY *) malloc(sizeof(STR_ARRAY) * (top_atomnum + 10));
if (atomres == NULL) {
	fprintf(stdout, "memory allocation error for *atomres\n");
	exit(1);
 }

rcrd(); 
fclose(fpinc);

if(output_atomtype_flag == 2) {
	ratcorr();
	fclose(fpat);
}
if(output_atomtype_flag == 2) {
	rbtcorr();
	fclose(fpbt);
}

rtop(1);
bondnum = bondnum_h + bondnum_noh;
if(bondnum > top_bondnum)
	bondnum = top_bondnum;
fclose(fpinp);
j=0;
for(i=0;i<resnum-1;i++) {
	if(water_flag == 0)
		if(strcmp(res[i].str, "WAT") == 0)
			resnum_water++;
	for(j=reshead[i]-1;j<  reshead[i+1] -1;j++) {
     		resno[j] =i+1;
     		strcpy(atomres[j].str, res[i].str);
  	}
}
for(i = reshead[resnum-1]-1; i<atomnum; i++) {
     	resno[i] =resnum;
     	strcpy(atomres[i].str, res[resnum-1].str);
}

if(water_flag == 0)
	if(strcmp(res[resnum-1].str, "WAT") == 0)
		resnum_water++;

for(i=0;i<bondnum_h;i++) {
	strcpy(bondtype_h[i].str, "1");
	bondat1_h[i]=bondat1_h[i]/3 ;
	bondat2_h[i]=bondat2_h[i]/3 ;
}

for(i=0;i<bondnum_noh;i++) {
	strcpy(bondtype_noh[i].str, "1");
	bondat1_noh[i]=bondat1_noh[i]/3 ;
	bondat2_noh[i]=bondat2_noh[i]/3 ;
}
/*
for(i=0;i<corratnum;i++)
	printf("\n%5d%5s%5s", i+1, atorg[i], atcorr[i]);
for(i=0;i<corrbtnum;i++)
  	printf("\n%5d%5s%5s%5s", i+1, btan1[i],btan2[i], bt[i]);
*/
for(i=0;i<atomnum;i++)
	strcpy(atomtypebak[i].str, atomtype[i].str);

if(output_atomtype_flag == 2) 
	for(i=0;i<atomnum;i++) {
		if(water_flag == 0 && strcmp(atomres[i].str, "WAT") == 0) 
			break;
  		for(j=0;j<corratnum;j++) 
			if(strcmp(atomtype[i].str, atorg[j])==0) {
				strcpy(atomtype[i].str, atcorr[j]);	
				break;
			}
	}
for(i=0;i<bondnum_h;i++) {
	if(water_flag == 0) {
		if(strcmp(atomres[bondat1_h[i]].str, "WAT") == 0) break;
		if(strcmp(atomres[bondat2_h[i]].str, "WAT") == 0) break;
	}
	if(output_bondtype_flag == 1) 
		strcpy(bondtype_h[i].str, "1");
	if(output_bondtype_flag == 2) {
		for(j=0;j<corrbtnum;j++) {
			if(strcmp(atomtypebak[bondat1_h[i]].str, btan1[j])==0 && strcmp(atomtypebak[bondat2_h[i]].str, btan2[j])==0) {
				strcpy(bondtype_h[i].str, bt[j]);	
				break;
			}
			if(strcmp(atomtypebak[bondat1_h[i]].str, btan2[j])==0 && strcmp(atomtypebak[bondat2_h[i]].str, btan1[j])==0) {
				strcpy(bondtype_h[i].str, bt[j]);	
				break;
			}
		} 	
	}
}

for(i=0;i<bondnum_noh;i++) {
	if(output_bondtype_flag == 1) 
		strcpy(bondtype_noh[i].str, "1");
	if(output_bondtype_flag == 2) {
		for(j=0;j<corrbtnum;j++) {
			if(strcmp(atomtypebak[bondat1_noh[i]].str, btan1[j])==0 && strcmp(atomtypebak[bondat2_noh[i]].str, btan2[j])==0) {
				strcpy(bondtype_noh[i].str, bt[j]);	
				break;
			}
			if(strcmp(atomtypebak[bondat1_noh[i]].str, btan2[j])==0 && strcmp(atomtypebak[bondat2_noh[i]].str, btan1[j])==0) {
				strcpy(bondtype_noh[i].str, bt[j]);	
				break;
			}
		} 	
	}
}
write();
fclose(fpout);
return(0);
}

