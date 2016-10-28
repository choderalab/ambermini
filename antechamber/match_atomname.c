/*
************************************************************************
*           All Copyright Reserved!                                    *
*                                                                      *
*  Prog:    respgen                                                    *
*  Version: version 1.0                                                *
*  Author:  Junmei Wang                                                *
*                                                                      *
*  Department of Pharmaceutical Chemistry                              *
*  School of Pharmacy                                                  *
*  University of California                                            *
*  San Francisco   CA 94143                                            *
*  October, 2001                                                      *
************************************************************************
*/
char *amberhome;

# include "common.h"
# include "define.h"
# include "atom.h"
# include "utility.c"
# include "common.c"
# include "equatom2.c"
# include "pdb.c"
# include "ac.c"
# include "prep.c"
# include "mol2.c"
# define debug 0

ATOM *atom;
BOND *bond;
AROM *arom;
int atomnum = 0;
int bondnum = 0;
int *hcount; /*how many hydrogen atoms before the current atom*/
char iformat[MAXCHAR];
char rformat[MAXCHAR];
CONTROLINFO cinfo;
MOLINFO minfo;

ATOM *input_atom;
BOND *input_bond;
int input_atomnum = 0;
int input_bondnum = 0;
int *input_atom_id;
int *input_equ_atom_id;

ATOM *ref_atom;
BOND *ref_bond;
int ref_atomnum = 0;
int ref_bondnum = 0;
CONTROLINFO ref_cinfo;
MOLINFO ref_minfo;
int *ref_atom_id;
int *ref_equ_atom_id;

ATOM atm;
BOND bnd;
int max_path_length = -1;
int igeom = 0;
int ihydrogen = 0;

int i, j, k, l;
FILE *fpin;
FILE *fpout;
FILE *fpcharge;
FILE *fpaddinfo;
char line[MAXCHAR];
char ifilename[MAXCHAR];
char ofilename[MAXCHAR];
char rfilename[MAXCHAR];

int overflow_flag = 0;

void wprepi(char *input_filename, char *output_filename, int atomnum, ATOM atom[]) {
int i;
int len;
int flag;
int id;
int suc;
char atomname[MAXCHAR];
char atomname1[MAXCHAR];
char atomname2[MAXCHAR];
char atomname3[MAXCHAR];
char atomname4[MAXCHAR];
char tmpchar1[MAXCHAR];
char tmpchar2[MAXCHAR];
char tmpchar3[MAXCHAR];
char tmpchar4[MAXCHAR];
int  tmpint1, tmpint2, tmpint3;
double x,y,z,c;
FILE *fpin, *fpout;

if ((fpin = fopen(input_filename, "r")) == NULL) {
        fprintf(stdout, "Cannot open file %s to read in rprepi(), exit\n", input_filename);
        return;
}
if ((fpout = fopen(output_filename, "w")) == NULL) {
        fprintf(stdout, "Cannot open file %s to write in wprepi(), exit\n", output_filename);
        return;
}
flag = 0;
for(;;) {
	if (fgets(line, MAXCHAR, fpin) == NULL) break;
	len = strlen(line);
	strcpy(tmpchar1, "");
	strcpy(tmpchar2, "");
	strcpy(tmpchar3, "");
	strcpy(tmpchar4, "");
	tmpchar1[0] = '\0';
	tmpchar2[0] = '\0';
	tmpchar3[0] = '\0';
	tmpchar4[0] = '\0';
        sscanf(line, "%s%s%s%s", tmpchar1, tmpchar2, tmpchar3, tmpchar4);
	if(strcmp(tmpchar1, "CORRECT") == 0 || strcmp(tmpchar1, "CHANGE") == 0) {
		fprintf(fpout, "%s", line);
		fgets(line, MAXCHAR, fpin);
		fprintf(fpout, "%s", line);
		fgets(line, MAXCHAR, fpin);
		fprintf(fpout, "%s", line);
		fgets(line, MAXCHAR, fpin);
		fprintf(fpout, "%s", line);
		fgets(line, MAXCHAR, fpin);
		fprintf(fpout, "%s", line);
		flag = 1;
		continue;	
	}	
	if(strcmp(tmpchar1, "LOOP") == 0) {
		flag = 2;
		fprintf(fpout, "%s", line);
		continue;
	}
	if(strcmp(tmpchar1, "IMPROPER") == 0) {
		flag = 3;
		fprintf(fpout, "%s", line);
		continue;
	}
	if(strcmp(tmpchar1, "DONE") == 0) {
		flag = 0;
		fprintf(fpout, "%s", line);
		continue;
	}
	if(len <= 2) flag = 0;
	if(flag == 0) {
		fprintf(fpout, "%s", line);
		continue;
	}
	if(flag == 1) {
		strcpy(atomname, "");
		atomname[0]='\0';
		sscanf(line, "%d%s%s%s%d%d%d%lf%lf%lf%lf", &id, tmpchar1, tmpchar2, tmpchar3, &tmpint1, &tmpint2, &tmpint3, &x, &y, &z, &c);
		if(ihydrogen == 1 || (ihydrogen == 0 && tmpchar1[0]!='H')) {
			suc = 0;
			for(i=0;i<atomnum;i++) {
				if(strcmp(tmpchar1, atom[i].aa) == 0) {
					strcpy(atomname, atom[i].name);
					suc = 1;		
					break;
				}
			}
			if(suc == 0) {
				fprintf(stderr, "Cannot find out matched atom names for atom name %s\n", tmpchar1);			
				exit(1);
			}
		}
		else
			strcpy(atomname, tmpchar1);
                fprintf(fpout, "%4d  %-6s%-6s%-4s%-4d%-4d%-4d%8.3lf%10.3lf%10.3lf%10.6f\n",
                                id, atomname, tmpchar2, tmpchar3, tmpint1, tmpint2, tmpint3, x, y, z, c); 
	}
	if(flag == 2) {
		strcpy(atomname1, "");
		atomname1[0]='\0';
		strcpy(atomname2, "");
		atomname2[0]='\0';
		if(ihydrogen == 1 || (ihydrogen == 0 && tmpchar1[0]!='H')) {
			suc = 0;
			for(i=0;i<atomnum;i++) {
				if(strcmp(tmpchar1, atom[i].aa) == 0) {
					strcpy(atomname1, atom[i].name);
					suc = 1;		
					break;
				}
			}
			if(suc == 0) {
				fprintf(stderr, "Cannot find out matched atom names for atom name %s\n", tmpchar1);			
				exit(1);
			}
		}
		else 
			strcpy(atomname1, tmpchar1);

		if(ihydrogen == 1 || (ihydrogen == 0 && tmpchar2[0]!='H')) {
			suc = 0;
			for(i=0;i<atomnum;i++) {
				if(strcmp(tmpchar2, atom[i].aa) == 0) {
					strcpy(atomname2, atom[i].name);
					suc = 1;		
					break;
				}
			}
			if(suc == 0) {
				fprintf(stderr, "Cannot find out matched atom names for atom name %s\n", tmpchar2);			
				exit(1);
			}
		}
		else 
			strcpy(atomname2, tmpchar2);

		fprintf(fpout, "%5s%5s\n", atomname1, atomname2);
	}
	if(flag == 3) {
		strcpy(atomname1, "");
		strcpy(atomname2, "");
		strcpy(atomname3, "");
		strcpy(atomname4, "");
		atomname1[0]='\0';
		atomname2[0]='\0';
		atomname3[0]='\0';
		atomname4[0]='\0';
		if(strcmp(tmpchar1, "+M") == 0) {
			strcpy(atomname1, "+M");
			suc = 1;
		}
		else if(strcmp(tmpchar1, "-M") == 0) {
			strcpy(atomname1, "-M");
			suc = 1;
		}
		else {
			if(ihydrogen == 1 || (ihydrogen == 0 && tmpchar1[0]!='H')) {
                		suc = 0;
                		for(i=0;i<atomnum;i++) {
                        		if(strcmp(tmpchar1, atom[i].aa) == 0) {
                                		strcpy(atomname1, atom[i].name);
                                		suc = 1;  
                                		break;
                        		}
                		}
                		if(suc == 0) {
                        		fprintf(stderr, "Cannot find out matched atom names for atom name %s\n", tmpchar1);
                        		exit(1);
                		}
			}
			else
				strcpy(atomname1, tmpchar1);
		}

                if(strcmp(tmpchar2, "+M") == 0) {
                        strcpy(atomname2, "+M");
                        suc = 1;
                }
                else if(strcmp(tmpchar2, "-M") == 0) {
                        strcpy(atomname2, "-M");
                        suc = 1;
                }
                else {
			if(ihydrogen == 1 || (ihydrogen == 0 && tmpchar2[0]!='H')) {
                        	suc = 0;
                        	for(i=0;i<atomnum;i++) {
                                	if(strcmp(tmpchar2, atom[i].aa) == 0) {
                                        	strcpy(atomname2, atom[i].name);
                                        	suc = 1;
                                        	break;
                                	}
                        	}
                        	if(suc == 0) {
                                	fprintf(stderr, "Cannot find out matched atom names for atom name %s\n", tmpchar2);
                                	exit(1);
                        	}
			}
			else 
				strcpy(atomname2, tmpchar2);
                }

                if(strcmp(tmpchar3, "+M") == 0) {
                        strcpy(atomname3, "+M");
                        suc = 1;
                }
                else if(strcmp(tmpchar3, "-M") == 0) {
                        strcpy(atomname3, "-M");
                        suc = 1;
                }
                else {
			if(ihydrogen == 1 || (ihydrogen == 0 && tmpchar3[0]!='H')) {
                        	suc = 0;
                        	for(i=0;i<atomnum;i++) {
                                	if(strcmp(tmpchar3, atom[i].aa) == 0) {
                                        	strcpy(atomname3, atom[i].name);
                                        	suc = 1;
                                        	break;
                                	}
                        	}
                        	if(suc == 0) {
                                	fprintf(stderr, "Cannot find out matched atom names for atom name %s\n", tmpchar3);
                                	exit(1);
                        	}
			}
			else 
				strcpy(atomname3, tmpchar3);
                }

                if(strcmp(tmpchar4, "+M") == 0) {
                        strcpy(atomname4, "+M");
                        suc = 1;
                }
                else if(strcmp(tmpchar4, "-M") == 0) {
                        strcpy(atomname4, "-M");
                        suc = 1;
                }
                else {
			if(ihydrogen == 1 || (ihydrogen == 0 && tmpchar4[0]!='H')) {
                        	suc = 0;
                        	for(i=0;i<atomnum;i++) {
                                	if(strcmp(tmpchar4, atom[i].aa) == 0) {
                                        	strcpy(atomname4, atom[i].name);
                                        	suc = 1;
                                        	break;
                                	}
                        	}
                        	if(suc == 0) {
                                	fprintf(stderr, "Cannot find out matched atom names for atom name %s\n", tmpchar4);
                                	exit(1);
                        	}
			}
			else 
				strcpy(atomname4, tmpchar4);
                }

                fprintf(fpout, "%5s%5s%5s%5s\n", atomname1, atomname2, atomname3, atomname4);

	}
}
fclose(fpin);
fclose(fpout);
}


void wprepc(char *input_filename, char *output_filename, int atomnum, ATOM atom[]) {
int i;
int len;
int flag;
int id;
int suc;
char atomname[MAXCHAR];
char atomname1[MAXCHAR];
char atomname2[MAXCHAR];
char atomname3[MAXCHAR];
char atomname4[MAXCHAR];
char tmpchar1[MAXCHAR];
char tmpchar2[MAXCHAR];
char tmpchar3[MAXCHAR];
char tmpchar4[MAXCHAR];
int  tmpint1, tmpint2, tmpint3;
double x,y,z,c;
FILE *fpin, *fpout;

if ((fpin = fopen(input_filename, "r")) == NULL) {
        fprintf(stdout, "Cannot open file %s to read in rprepc(), exit\n", input_filename);
        return;
}
if ((fpout = fopen(output_filename, "w")) == NULL) {
        fprintf(stdout, "Cannot open file %s to write in wprepc(), exit\n", output_filename);
        return;
}
flag = 0;
for(;;) {
	if (fgets(line, MAXCHAR, fpin) == NULL) break;
	len = strlen(line);
	strcpy(tmpchar1, "");
	strcpy(tmpchar2, "");
	strcpy(tmpchar3, "");
	strcpy(tmpchar4, "");
	tmpchar1[0] = '\0';
	tmpchar2[0] = '\0';
	tmpchar3[0] = '\0';
	tmpchar4[0] = '\0';
        sscanf(line, "%s%s%s%s", tmpchar1, tmpchar2, tmpchar3, tmpchar4);
	if(strcmp(tmpchar1, "CORRECT") == 0 || strcmp(tmpchar1, "CHANGE") == 0) {
		fprintf(fpout, "%s", line);
		fgets(line, MAXCHAR, fpin);
		fprintf(fpout, "%s", line);
		fgets(line, MAXCHAR, fpin);
		fprintf(fpout, "%s", line);
		fgets(line, MAXCHAR, fpin);
		fprintf(fpout, "%s", line);
		fgets(line, MAXCHAR, fpin);
		fprintf(fpout, "%s", line);
		flag = 1;
		continue;	
	}	
	if(strcmp(tmpchar1, "LOOP") == 0) {
		flag = 2;
		fprintf(fpout, "%s", line);
		continue;
	}
	if(strcmp(tmpchar1, "IMPROPER") == 0) {
		flag = 3;
		fprintf(fpout, "%s", line);
		continue;
	}
	if(strcmp(tmpchar1, "DONE") == 0) {
		flag = 0;
		fprintf(fpout, "%s", line);
		continue;
	}
	if(len <= 2) flag = 0;
	if(flag == 0) {
		fprintf(fpout, "%s", line);
		continue;
	}
	if(flag == 1) {
		strcpy(atomname, "");
		atomname[0]='\0';
		sscanf(line, "%d%s%s%s%lf%lf%lf%lf", &id, tmpchar1, tmpchar2, tmpchar3, &x, &y, &z, &c);
		if(ihydrogen == 1 || (ihydrogen == 0 && tmpchar1[0]!='H')) {
			suc = 0;
			for(i=0;i<atomnum;i++) {
				if(strcmp(tmpchar1, atom[i].aa) == 0) {
					strcpy(atomname, atom[i].name);
					suc = 1;		
					break;
				}
			}
			if(suc == 0) {
				fprintf(stderr, "Cannot find out matched atom names for atom name %s\n", tmpchar1);			
				exit(1);
			}
		}
		else
			strcpy(atomname, tmpchar1);
                fprintf(fpout, "%4d  %-6s%-6s%-4s%15.6lf%12.6lf%12.6lf%12.6f",
                                id, atomname, tmpchar2, tmpchar3, x, y, z, c); 
	}
	if(flag == 2) {
		strcpy(atomname1, "");
		atomname1[0]='\0';
		strcpy(atomname2, "");
		atomname2[0]='\0';
		if(ihydrogen == 1 || (ihydrogen == 0 && tmpchar1[0]!='H')) {
			suc = 0;
			for(i=0;i<atomnum;i++) {
				if(strcmp(tmpchar1, atom[i].aa) == 0) {
					strcpy(atomname1, atom[i].name);
					suc = 1;		
					break;
				}
			}
			if(suc == 0) {
				fprintf(stderr, "Cannot find out matched atom names for atom name %s\n", tmpchar1);			
				exit(1);
			}
		}
		else 
			strcpy(atomname1, tmpchar1);

		if(ihydrogen == 1 || (ihydrogen == 0 && tmpchar2[0]!='H')) {
			suc = 0;
			for(i=0;i<atomnum;i++) {
				if(strcmp(tmpchar2, atom[i].aa) == 0) {
					strcpy(atomname2, atom[i].name);
					suc = 1;		
					break;
				}
			}
			if(suc == 0) {
				fprintf(stderr, "Cannot find out matched atom names for atom name %s\n", tmpchar2);			
				exit(1);
			}
		}
		else 
			strcpy(atomname2, tmpchar2);

		fprintf(fpout, "%5s%5s\n", atomname1, atomname2);
	}
	if(flag == 3) {
		strcpy(atomname1, "");
		strcpy(atomname2, "");
		strcpy(atomname3, "");
		strcpy(atomname4, "");
		atomname1[0]='\0';
		atomname2[0]='\0';
		atomname3[0]='\0';
		atomname4[0]='\0';
		if(strcmp(tmpchar1, "+M") == 0) {
			strcpy(atomname1, "+M");
			suc = 1;
		}
		else if(strcmp(tmpchar1, "-M") == 0) {
			strcpy(atomname1, "-M");
			suc = 1;
		}
		else {
			if(ihydrogen == 1 || (ihydrogen == 0 && tmpchar1[0]!='H')) {
                		suc = 0;
                		for(i=0;i<atomnum;i++) {
                        		if(strcmp(tmpchar1, atom[i].aa) == 0) {
                                		strcpy(atomname1, atom[i].name);
                                		suc = 1;  
                                		break;
                        		}
                		}
                		if(suc == 0) {
                        		fprintf(stderr, "Cannot find out matched atom names for atom name %s\n", tmpchar1);
                        		exit(1);
                		}
			}
			else
				strcpy(atomname1, tmpchar1);
		}

                if(strcmp(tmpchar2, "+M") == 0) {
                        strcpy(atomname2, "+M");
                        suc = 1;
                }
                else if(strcmp(tmpchar2, "-M") == 0) {
                        strcpy(atomname2, "-M");
                        suc = 1;
                }
                else {
			if(ihydrogen == 1 || (ihydrogen == 0 && tmpchar2[0]!='H')) {
                        	suc = 0;
                        	for(i=0;i<atomnum;i++) {
                                	if(strcmp(tmpchar2, atom[i].aa) == 0) {
                                        	strcpy(atomname2, atom[i].name);
                                        	suc = 1;
                                        	break;
                                	}
                        	}
                        	if(suc == 0) {
                                	fprintf(stderr, "Cannot find out matched atom names for atom name %s\n", tmpchar2);
                                	exit(1);
                        	}
			}
			else 
				strcpy(atomname2, tmpchar2);
                }

                if(strcmp(tmpchar3, "+M") == 0) {
                        strcpy(atomname3, "+M");
                        suc = 1;
                }
                else if(strcmp(tmpchar3, "-M") == 0) {
                        strcpy(atomname3, "-M");
                        suc = 1;
                }
                else {
			if(ihydrogen == 1 || (ihydrogen == 0 && tmpchar3[0]!='H')) {
                        	suc = 0;
                        	for(i=0;i<atomnum;i++) {
                                	if(strcmp(tmpchar3, atom[i].aa) == 0) {
                                        	strcpy(atomname3, atom[i].name);
                                        	suc = 1;
                                        	break;
                                	}
                        	}
                        	if(suc == 0) {
                                	fprintf(stderr, "Cannot find out matched atom names for atom name %s\n", tmpchar3);
                                	exit(1);
                        	}
			}
			else 
				strcpy(atomname3, tmpchar3);
                }

                if(strcmp(tmpchar4, "+M") == 0) {
                        strcpy(atomname4, "+M");
                        suc = 1;
                }
                else if(strcmp(tmpchar4, "-M") == 0) {
                        strcpy(atomname4, "-M");
                        suc = 1;
                }
                else {
			if(ihydrogen == 1 || (ihydrogen == 0 && tmpchar4[0]!='H')) {
                        	suc = 0;
                        	for(i=0;i<atomnum;i++) {
                                	if(strcmp(tmpchar4, atom[i].aa) == 0) {
                                        	strcpy(atomname4, atom[i].name);
                                        	suc = 1;
                                        	break;
                                	}
                        	}
                        	if(suc == 0) {
                                	fprintf(stderr, "Cannot find out matched atom names for atom name %s\n", tmpchar4);
                                	exit(1);
                        	}
			}
			else 
				strcpy(atomname4, tmpchar4);
                }

                fprintf(fpout, "%5s%5s%5s%5s\n", atomname1, atomname2, atomname3, atomname4);

	}
}
fclose(fpin);
fclose(fpout);
}

int main(int argc, char *argv[])
{
	int i;
	int atid;
	int bondi, bondj;
	int numatom;
	int numbond;
	int numcon;
	int nhydrogen;

    	amberhome = (char *) getenv("AMBERHOME");
    	if( amberhome == NULL ){
       		fprintf( stdout, "Error: AMBERHOME is not set!\n" );
       		exit(1);
    	}

	if (strcmp(COLORTEXT, "YES") == 0 || strcmp(COLORTEXT, "yes") == 0) {
		if (argc == 2
			&& (strcmp(argv[1], "-h") == 0
				|| strcmp(argv[1], "-H") == 0)) {
			printf("[31mUsage: match_atomname -i [0m input file name\n"
			       "[31m                      -fi[0m input format (pdb, ac, prepi, prepc, mol2)\n"
			       "[31m                      -r [0m reference file name\n"
			       "[31m                      -fr[0m reference format (pdb, ac, prepi, prepc, mol2)\n"
			       "[31m                      -o [0m output file name\n"
			       "[31m                      -h [0m include hydrogen atoms or not\n"
			       "[32m                          0[0m not, the default\n"
			       "[32m                          1[0m yes\n"
			       "[31m                      -g [0m geometric info (such as E/Z configuration) is considered to describe chemical environment\n"
			       "[32m                          0[0m no, the default\n"
			       "[32m                          1[0m yes\n"
			       "[31m                      -l [0m maximum path length, default is -1 (full length)\n"
			       "[31m                         [0m if it takes very long time and/or core dump occur, a value between 8 to 10 is recommended\n");
			exit(1);
		}
		if (argc != 11 && argc != 13 && argc !=15 && argc !=17) {
			printf("[31mUsage: match_atomname -i [0m input file name\n"
			       "[31m                      -fi[0m input format (pdb, ac, prepi, prepc, mol2)\n"
			       "[31m                      -r [0m reference file name\n"
			       "[31m                      -fr[0m reference format (pdb, ac, prepi, prepc, mol2)\n"
			       "[31m                      -o [0m output file name\n"
			       "[31m                      -h [0m include hydrogen atoms or not\n"
			       "[32m                          0[0m not, the default\n"
			       "[32m                          1[0m yes\n"
			       "[31m                      -g [0m geometric info (such as E/Z configuration) is considered to describe chemical environment\n"
			       "[32m                          0[0m no, the default\n"
			       "[32m                          1[0m yes\n"
			       "[31m                      -l [0m maximum path length, default is -1 (full length)\n"
			       "[31m                         [0m if it takes very long time and/or core dump occur, a value between 8 to 10 is recommended\n");
			exit(1);
		}
	} else {
		if (argc == 2
			&& (strcmp(argv[1], "-h") == 0
				|| strcmp(argv[1], "-H") == 0)) {
			printf("Usage: match_atomname -i  input file name\n"
			       "                      -fi input format (pdb, ac, prepi, prepc, mol2)\n"
			       "                      -r  reference file name\n"
			       "                      -fr reference format (pdb, ac, prepi, prepc, mol2)\n"
			       "                      -o  output file name\n"
			       "                      -h  include hydrogen atoms or not\n"
			       "                          0 not, the default\n"
			       "                          1 yes\n"
			       "                      -g  geometric info (such as E/Z configuration) is considered to describe chemical environment\n"
			       "                          0 no, the default\n"
			       "                          1 yes\n"
			       "                      -l  maximum path length, default is -1 (full length)\n"
			       "                          if it takes very long time and/or core dump occur, a value between 8 to 10 is recommended\n");
			exit(1);
		}
		if (argc != 11 && argc != 13 && argc !=15 && argc !=17) {
			printf("Usage: match_atomname -i  input file name\n"
			       "                      -fi input format (pdb, ac, prepi, prepc, mol2)\n"
			       "                      -r  reference file name in ac\n"
			       "                      -fr reference format (pdb, ac, prepi, prepc, mol2)\n"
			       "                      -o  output file name\n"
			       "                      -h  include hydrogen atoms or not\n"
			       "                          0 not, the default\n"
			       "                          1 yes\n"
			       "                      -g  geometric info (such as E/Z configuration) is considered to describe chemical environment\n"
			       "                          0 no, the default\n"
			       "                          1 yes\n"
			       "                      -l  maximum path length, default is -1 (full length)\n"
			       "                          if it takes very long time and/or core dump occur, a value between 8 to 10 is recommended\n");
			exit(1);
		}
	}

	max_path_length = -1;
	for (i = 1; i < argc; i += 2) {
		if (strcmp(argv[i], "-i") == 0)
			strcpy(ifilename, argv[i + 1]);
		if (strcmp(argv[i], "-fi") == 0)
			strcpy(iformat, argv[i + 1]);
		if (strcmp(argv[i], "-fr") == 0)
			strcpy(rformat, argv[i + 1]);
		if (strcmp(argv[i], "-o") == 0)
			strcpy(ofilename, argv[i + 1]);
		if (strcmp(argv[i], "-r") == 0)
			strcpy(rfilename, argv[i + 1]);
		if (strcmp(argv[i], "-g") == 0)
			igeom = atoi(argv[i+1]); 
		if (strcmp(argv[i], "-l") == 0)
			max_path_length = atoi(argv[i+1]); 
		if (strcmp(argv[i], "-h") == 0)
			ihydrogen = atoi(argv[i+1]); 
	}
	if(ihydrogen !=0 && ihydrogen !=1) 
		ihydrogen = 0;
	if(igeom !=0 && igeom !=1) 
		igeom = 0;
	if(strcmp(iformat, "pdb") != 0 && strcmp(iformat, "ac") != 0 &&
	   strcmp(iformat, "prepi") != 0 && strcmp(iformat, "prepc") != 0 && strcmp(iformat, "mol2") !=0) {
		fprintf(stderr, "The input format, %s, is not supported\n", iformat);    
		exit(1);
	}
	if(strcmp(rformat, "pdb") != 0 && strcmp(rformat, "ac") != 0 &&
	   strcmp(rformat, "prepi") != 0 && strcmp(rformat, "prepc") != 0 && strcmp(rformat, "mol2") !=0) {
		fprintf(stderr, "The input format, %s, is not supported\n", rformat);    
		exit(1);
	}

	default_cinfo(&cinfo);
	default_minfo(&minfo);
	default_cinfo(&ref_cinfo);
	default_minfo(&ref_minfo);

	atom = (ATOM *) malloc(sizeof(ATOM) * cinfo.maxatom);
	if (atom == NULL) {
		fprintf(stdout, "memory allocation error for *atom\n");
		exit(1);
	}
	ref_atom = (ATOM *) malloc(sizeof(ATOM) * ref_cinfo.maxatom);
	if (ref_atom == NULL) {
		fprintf(stdout, "memory allocation error for *ref_atom\n");
		exit(1);
	}

	bond = (BOND *) malloc(sizeof(BOND) * cinfo.maxbond);
	if (bond == NULL) {
		fprintf(stdout, "memory allocation error for *bond\n");
		exit(1);
	}
	ref_bond = (BOND *) malloc(sizeof(BOND) * ref_cinfo.maxbond);
	if (ref_bond == NULL) {
		fprintf(stdout, "memory allocation error for *ref_bond\n");
		exit(1);
	}


/*for input atom */
	if(strcmp(iformat, "ac") == 0) {
		overflow_flag =
		rac(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo);
	}
	if(strcmp(iformat, "pdb") == 0) {
        	overflow_flag =
                rpdb(ifilename, &atomnum, atom, cinfo, minfo, 0);
	}
	if(strcmp(iformat, "mol2") == 0) {
        	overflow_flag =
                rmol2(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo, 0);
	}
	if(strcmp(iformat, "prepi") == 0) {
        	overflow_flag =
                rprepi(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo);
	}
	if(strcmp(iformat, "prepc") == 0) {
        	overflow_flag =
                rprepc(ifilename, &atomnum, atom, &cinfo, &minfo);
	}

	if (overflow_flag) {
		cinfo.maxatom = atomnum + 10;
		cinfo.maxbond = bondnum + 10;
		free(atom);
		free(bond);
		atom = (ATOM *) malloc(sizeof(ATOM) * cinfo.maxatom);
		if (atom == NULL) {
			fprintf(stdout, "memory allocation error for *atom\n");
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

		if(strcmp(iformat, "ac") == 0) {
			overflow_flag =
			rac(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo);
		}
		if(strcmp(iformat, "pdb") == 0) {
        		overflow_flag =
                	rpdb(ifilename, &atomnum, atom, cinfo, minfo, 0);
		}
		if(strcmp(iformat, "mol2") == 0) {
        		overflow_flag =
                	rmol2(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo, 0);
		}
		if(strcmp(iformat, "prepi") == 0) {
        		overflow_flag =
                	rprepi(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo);
		}
		if(strcmp(iformat, "prepc") == 0) {
        		overflow_flag =
                	rprepc(ifilename, &atomnum, atom, &cinfo, &minfo);
		}
	}

/*for ref atom */
	if(strcmp(rformat, "ac") == 0) {
		overflow_flag =
		rac(rfilename, &ref_atomnum, ref_atom, &ref_bondnum, ref_bond, &ref_cinfo, &ref_minfo);
	}
	if(strcmp(rformat, "pdb") == 0) {
        	overflow_flag =
                rpdb(rfilename, &ref_atomnum, ref_atom, ref_cinfo, ref_minfo, 0);
	}
	if(strcmp(rformat, "mol2") == 0) {
        	overflow_flag =
                rmol2(rfilename, &ref_atomnum, ref_atom, &ref_bondnum, ref_bond, &ref_cinfo, &ref_minfo, 0);
	}
	if(strcmp(rformat, "prepi") == 0) {
        	overflow_flag =
                rprepi(rfilename, &ref_atomnum, ref_atom, &ref_bondnum, ref_bond, &ref_cinfo, &ref_minfo);
	}
	if(strcmp(rformat, "prepc") == 0) {
        	overflow_flag =
                rprepc(rfilename, &ref_atomnum, ref_atom, &ref_cinfo, &ref_minfo);
	}

	if (overflow_flag) {
		ref_cinfo.maxatom = ref_atomnum + 10;
		ref_cinfo.maxbond = ref_bondnum + 10;
		free(ref_atom);
		free(ref_bond);
		ref_atom = (ATOM *) malloc(sizeof(ATOM) * ref_cinfo.maxatom);
		if (ref_atom == NULL) {
			fprintf(stdout, "memory allocation error for *ref_atom\n");
			exit(1);
		}
		ref_bond = (BOND *) malloc(sizeof(BOND) * ref_cinfo.maxbond);
		if (ref_bond == NULL) {
			fprintf(stdout, "memory allocation error for *ref_bond\n");
			exit(1);
		}
		for (i = 0; i < ref_cinfo.maxbond; ++i) {
			ref_bond[i].jflag = -1; /* bond type has not been assigned */
		}
		if(strcmp(rformat, "ac") == 0) {
			overflow_flag =
			rac(rfilename, &ref_atomnum, ref_atom, &ref_bondnum, ref_bond, &ref_cinfo, &ref_minfo);
		}
		if(strcmp(rformat, "pdb") == 0) {
        		overflow_flag =
                	rpdb(rfilename, &ref_atomnum, ref_atom, ref_cinfo, ref_minfo, 0);
		}
		if(strcmp(rformat, "mol2") == 0) {
        		overflow_flag =
                	rmol2(rfilename, &ref_atomnum, ref_atom, &ref_bondnum, ref_bond, &ref_cinfo, &ref_minfo, 0);
		}	
		if(strcmp(rformat, "prepi") == 0) {
        		overflow_flag =
                	rprepi(rfilename, &ref_atomnum, ref_atom, &ref_bondnum, ref_bond, &ref_cinfo, &ref_minfo);
		}
		if(strcmp(rformat, "prepc") == 0) {
        		overflow_flag =
                	rprepc(rfilename, &ref_atomnum, ref_atom, &ref_cinfo, &ref_minfo);
		}
	}
	input_atom = (ATOM *) malloc(sizeof(ATOM) * cinfo.maxatom);
	if (input_atom == NULL) {
		fprintf(stdout, "memory allocation error for *input_atom\n");
		exit(1);
	}
	input_bond = (BOND *) malloc(sizeof(BOND) * cinfo.maxbond);
	if (input_bond == NULL) {
		fprintf(stdout, "memory allocation error for *input_bond\n");
		exit(1);
	}

	if(cinfo.maxatom > ref_cinfo.maxatom) 
		hcount = (int *) malloc(sizeof(int) * cinfo.maxatom);
	else
		hcount = (int *) malloc(sizeof(int) * ref_cinfo.maxatom);
	if (hcount == NULL) {
		fprintf(stdout, "memory allocation error for *hcount\n");
		exit(1);
	}

	for (i = 0; i < cinfo.maxbond; ++i) {
		input_bond[i].jflag = -1; /* bond type has not been assigned */
	}
	for (i = 0; i < ref_cinfo.maxbond; ++i) {
		ref_bond[i].jflag = -1; /* bond type has not been assigned */
	}

/*for input atom */
/*      for connect.tpl and radius parameter files */
        build_dat_path(minfo.connect_file, "CONNECT.TPL",
                sizeof minfo.connect_file, 0);
        build_dat_path(ref_minfo.connect_file, "CONNECT.TPL",
                sizeof ref_minfo.connect_file, 0);

	atomicnum(atomnum, atom);
	adjustatomname(atomnum, atom, 1);
	atomicnum(ref_atomnum, ref_atom);
	adjustatomname(ref_atomnum, ref_atom, 1);

	if(strcmp(iformat, "pdb")   == 0 ||
	   strcmp(iformat, "prepc") == 0 || 
	   strcmp(iformat, "prepi") == 0) {
        	connect(minfo.connect_file, atomnum, atom, &bondnum, bond, cinfo.maxbond);
	}
	if(strcmp(rformat, "pdb")   == 0 || 
	   strcmp(rformat, "prepc") == 0 ||
	   strcmp(rformat, "prepi") == 0) {
        	connect(ref_minfo.connect_file, ref_atomnum, ref_atom, &ref_bondnum, ref_bond, ref_cinfo.maxbond);
	}
/*	we make a copy of atomname to aa for prepi/prepc, which is not used anyway*/
	if(strcmp(iformat, "prepc") == 0 || strcmp(iformat, "prepi") == 0) {
        	for (i = 0; i < atomnum; i++) 
			strcpy(atom[i].aa, atom[i].name);
	}
/* 
	generate input_atom, input_bond, ref_atom and ref_bond
	need to adjust the atom.con and bond.bondi, bond.bondj accordingly
*/
	if(ihydrogen == 0) {
/*for input*/
		nhydrogen = 0;
		numatom = 0;
		for(i=0;i<atomnum;i++) {
			if(atom[i].atomicnum == 1) {
				hcount[i] = -1;
				nhydrogen++;
			}
			else {
				hcount[i] = nhydrogen;
				input_atom[numatom]=atom[i];
				numatom++;
			}
		}
		input_atomnum = numatom;	
		numbond = 0;
		for(i=0;i<bondnum;i++) {
			bondi = bond[i].bondi ;
			bondj = bond[i].bondj ;
			if(hcount[bondi] < 0 || hcount[bondj] < 0)	
				continue;
			input_bond[numbond] = bond[i];
			input_bond[numbond].bondi -= hcount[bondi];
			input_bond[numbond].bondj -= hcount[bondj];
			numbond++;
		}
		input_bondnum = numbond;	
		for(i=0; i< input_atomnum; i++) {
			input_atom[i].connum = 0;
			input_atom[i].con[0] = -1;
			input_atom[i].con[1] = -1;
			input_atom[i].con[2] = -1;
			input_atom[i].con[3] = -1;
			input_atom[i].con[4] = -1;
			input_atom[i].con[5] = -1;
		}
		for(i=0; i< input_bondnum; i++) {
			bondi = input_bond[i].bondi ;
			bondj = input_bond[i].bondj ;
			numcon = input_atom[bondi].connum;
			input_atom[bondi].con[numcon] = bondj;		
			input_atom[bondi].connum ++;
			numcon = input_atom[bondj].connum;
			input_atom[bondj].con[numcon] = bondi;		
			input_atom[bondj].connum ++;
		}
/* for ref*/
	        nhydrogen = 0;
                numatom = 0;
                for(i=0;i<ref_atomnum;i++) {
                        if(ref_atom[i].atomicnum == 1) {
                                hcount[i] = -1;
                                nhydrogen++;
                        }
                        else {
                                hcount[i] = nhydrogen;
                                ref_atom[numatom]=ref_atom[i];
                                numatom++;
                        }
                }
                ref_atomnum = numatom;        

                numbond = 0;
                for(i=0;i<ref_bondnum;i++) {
			bondi = ref_bond[i].bondi;
			bondj = ref_bond[i].bondj;
                        if(hcount[bondi] < 0 || hcount[bondj] < 0)      
                                continue;
			ref_bond[numbond] = ref_bond[i];
                        ref_bond[numbond].bondi -= hcount[bondi];
                        ref_bond[numbond].bondj -= hcount[bondj];
                        numbond++;
                }
                ref_bondnum = numbond;        
                for(i=0; i< ref_atomnum; i++) {
                        ref_atom[i].connum = 0;
                        ref_atom[i].con[0] = -1;
                        ref_atom[i].con[1] = -1;
                        ref_atom[i].con[2] = -1;
                        ref_atom[i].con[3] = -1;
                        ref_atom[i].con[4] = -1;
                        ref_atom[i].con[5] = -1;
                }
                for(i=0; i< ref_bondnum; i++) {
                        bondi = ref_bond[i].bondi ;
                        bondj = ref_bond[i].bondj ;
                        numcon = ref_atom[bondi].connum;
                        ref_atom[bondi].con[numcon] = bondj;
                        ref_atom[bondi].connum ++;
                        numcon = ref_atom[bondj].connum;
                        ref_atom[bondj].con[numcon] = bondi;
                        ref_atom[bondj].connum ++;
                }
	}
	if(ihydrogen == 1) {
		for(i=0;i<atomnum;i++)
			input_atom[i]=atom[i];
		input_atomnum = atomnum;
		for(i=0;i<bondnum;i++) {
			input_bond[i]=bond[i];
		}
		input_bondnum = bondnum;
	}

/* 	continue to allocate memory */
	input_equ_atom_id = (int *) malloc(sizeof(int) * input_atomnum);
	if (input_equ_atom_id == NULL) {
		fprintf(stdout, "memory allocation error for *input_equ_atom_id\n");
		exit(1);
	}
	input_atom_id = (int *) malloc(sizeof(int) * input_atomnum);
	if (input_atom_id == NULL) {
		fprintf(stdout, "memory allocation error for *input_atom_id\n");
		exit(1);
	}
	ref_equ_atom_id = (int *) malloc(sizeof(int) * ref_atomnum);
	if (ref_equ_atom_id == NULL) {
		fprintf(stdout, "memory allocation error for *ref_equ_atom_id\n");
		exit(1);
	}
	ref_atom_id = (int *) malloc(sizeof(int) * ref_atomnum);
	if (ref_atom_id == NULL) {
		fprintf(stdout, "memory allocation error for *ref_atom_id\n");
		exit(1);
	}
        for (i = 0; i < input_atomnum; i++) {
                input_equ_atom_id[i] = -1;
                input_atom_id[i] = -1;
	}
        for (i = 0; i < ref_atomnum; i++) {
                ref_equ_atom_id[i] = -1;
                ref_atom_id[i] = -1;
	}
	wac("MATCH_ATOMNAME_INPUT.AC", input_atomnum, input_atom, input_bondnum, input_bond, cinfo, minfo);
	wac("MATCH_ATOMNAME_REF.AC", ref_atomnum, ref_atom, ref_bondnum, ref_bond, ref_cinfo, ref_minfo);
/*	begin the real work */
	if(input_atomnum != ref_atomnum) {
		printf("Error: The numbers of atoms are different for the input and the reference\n"
		       "       molecules: %d versus %d \n", input_atomnum, ref_atomnum);
		exit(1);
	}
	if(max_path_length > 0) { 
		if(max_path_length > input_atomnum) 
			max_path_length = -1;
	}
	match_atompair(input_atomnum, ref_atomnum, input_atom, ref_atom, input_equ_atom_id, ref_equ_atom_id, input_atom_id, ref_atom_id, input_bondnum, ref_bondnum, input_bond, ref_bond, max_path_length, igeom); 
	if(debug == 1) {
		for(i=0;i<input_atomnum;i++) {
			if(input_atom_id[i] >=0)
				printf("Input %5d %5d %5s %5s\n", i+1, input_atom_id[i]+1, input_atom[i].name, ref_atom[input_atom_id[i]].name);
			else
				printf("Input %5d %5d %5s\n", i+1, input_atom_id[i]+1, input_atom[i].name);
		}
	}
/*	revise atom names */
	numatom = 0;
	for(i=0;i<atomnum;i++) {
		if(ihydrogen == 0) {
			if(atom[i].atomicnum == 1) continue;
			if(input_atom_id[numatom] >=0) 
				strcpy(atom[i].name, ref_atom[input_atom_id[numatom]].name);
			else
				fprintf(stderr, "Not successfully find the substituted atom name in the reference for ATOM %5d (%5s)\n", i+1, atom[i].name);
			numatom++;
		}
		else  {
			if(input_atom_id[i] >= 0) 
				strcpy(atom[i].name, ref_atom[input_atom_id[i]].name);
			else
				fprintf(stderr, "Not successfully find the substituted atom name in the reference for ATOM %5d (%5s)\n", i+1, atom[i].name);
		}
	}
	if(strcmp(iformat, "ac") == 0) 
		wac(ofilename, atomnum, atom, bondnum, bond, cinfo, minfo);
	if(strcmp(iformat, "pdb") == 0) 
		wpdb(ofilename, atomnum, atom); 
	if(strcmp(iformat, "mol2") == 0) 
	        wmol2(ofilename, atomnum, atom, bondnum, bond, arom, cinfo, minfo);
	if(strcmp(iformat, "prepi") == 0) 
	        wprepi(ifilename, ofilename, atomnum, atom); 
	if(strcmp(iformat, "prepc") == 0) 
	        wprepc(ifilename, ofilename, atomnum, atom); 
	return (0);

}
