# include <stdio.h>
# include <math.h>
# include <ctype.h>
# include <stdlib.h>
# include <string.h>
# define MAXCHAR 256 
# define COLORTEXT "YES" 
# define debug 0
# define MAXATOMTYPE 512
# define MAXPARM     512
# define FACTOR	 6.74834256
// This program reads in a topology file and write out a new topology file suitable for polarization calculation
typedef struct {
	char name[10];
	char type[10];
	double mass;
	double pol;
	double damp_factor;
} ATOM;

int atomnum = 0;
ATOM *atom;

// atom types
typedef struct {
	char name[10];
	char corr[10];
	int  pid;
} TYPE;
TYPE type[MAXPARM];
int typenum = 0;

// parameters of a polarizability model
typedef struct {
	char name[10];
	double pol;
	double damp_factor;
} PARM;

PARM parm[MAXPARM];
int parmnum = 0;

char ifilename[MAXCHAR];
char ofilename[MAXCHAR];
char pfilename[MAXCHAR];
char tfilename[MAXCHAR];
FILE *fpin,  *fpout ;
FILE *fpparm, *fptype;

char line[MAXCHAR];
char tmpchar1[MAXCHAR], tmpchar2[MAXCHAR], tmpchar3[MAXCHAR];
int i,j,k=0;
int tmpint1, tmpint2;
int ireplace = 1;
int ireadparm = 0;
int ireadat = 0;
int idamp_factor = 1;
char model[10];

double damp_factor = -999999;
int  pol_flag = 0; 
int  df_flag = 0;

void rparm() {
int i;
int read_flag;
char line[MAXCHAR];
char tmpchar[10];
if((fpparm=fopen(pfilename,"r"))==NULL) {
        fprintf(stdout, "\n Cannot open parm file to read: %s, exit", pfilename);
        exit(1);
}
read_flag = 0;
for(;;) {
        if(fgets(line, MAXCHAR, fpparm)==NULL) break;
	if(read_flag == 0 && strncmp(line, "MODEL_ID", 8) == 0) {
		sscanf(&line[9], "%s", tmpchar);
		if(strcmp(tmpchar, model) == 0) read_flag = 1;
		continue;
	}
	if(read_flag == 1 && strncmp(line, "MODEL_ID", 8) == 0) {
		read_flag = 0;
		break;
	} 
	if(read_flag == 1) {
		if(strncmp(line, "SCREEN_LENGTH", 13) == 0) 
        		sscanf(&line[14], "%lf", &damp_factor);
		if(strncmp(line, "PARM", 4) == 0) {
			parm[parmnum].damp_factor = damp_factor;
        		sscanf(&line[4], "%s%lf%lf", parm[parmnum].name, &parm[parmnum].pol, &parm[parmnum].damp_factor);
			parm[parmnum].pol /= FACTOR ;
			parmnum++;
		}
	}
}
if(debug == 1) {
	for(i=0;i<parmnum;i++)
		printf("PARM	%d	%5s	%8.4lf	%8.4lf\n", i+1, parm[i].name, parm[i].pol, parm[i].damp_factor);
}
fclose(fpparm);
}

void rparm2() {
int i;
char line[MAXCHAR];
if((fpparm=fopen(pfilename,"r"))==NULL) {
        fprintf(stdout, "\n Cannot open parm file to read: %s, exit", pfilename);
        exit(1);
}
for(;;) {
        if(fgets(line, MAXCHAR, fpparm)==NULL) break;
	if(strncmp(line, "SCREEN_LENGTH", 13) == 0) 
        	sscanf(&line[14], "%lf", &damp_factor);
	if(strncmp(line, "PARM", 4) == 0) {
		parm[parmnum].damp_factor = damp_factor;
        	sscanf(&line[4], "%s%lf%lf", parm[parmnum].name, &parm[parmnum].pol, &parm[parmnum].damp_factor);
		parm[parmnum].pol/= FACTOR;
		parmnum++;
	}
}
if(debug == 1) {
	for(i=0;i<parmnum;i++)
		printf("PARM	%d	%5s	%8.4lf	%8.4lf\n", i+1, parm[i].name, parm[i].pol, parm[i].damp_factor);
}
fclose(fpparm);
}

void rcorr() {
int i;
int suc;
char line[MAXCHAR];
if((fptype=fopen(tfilename,"r"))==NULL) {
        fprintf(stdout, "\n Cannot open atom type file to read: %s, exit", tfilename);
        exit(1);
}
for(;;) {
        if(fgets(line, MAXCHAR, fptype)==NULL) break;
	if(strncmp(line, "CORR", 4) == 0) {
        	sscanf(&line[4], "%s%s", type[typenum].name, type[typenum].corr);
		suc = 0;
                for(i=0;i<parmnum;i++) 
			if(strcmp(parm[i].name, type[typenum].corr) == 0) {
				type[typenum].pid = i;
				suc = 1;
				break;
			}
		if(suc == 0) {
			fprintf(stdout, "The atom type (%s),  which is corresponding to %s, does not show up in the parm file %s !\n", 
				type[typenum].name, type[typenum].corr, pfilename);
			exit(0);
		}
		typenum++;
	}
}
if(debug == 1) {
	for(i=0;i<typenum;i++)
		printf("TYPE	%d	%5s	%5s	%5d\n", i+1, type[i].name, type[i].corr, type[i].pid);
}
fclose(fptype);
}

void rtopparm() {
int tmpint1, tmpint2, tmpint3, tmpint4, tmpint5;
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
	if(fgets(line, MAXCHAR, fpin)==NULL) break;
    	sscanf(line, "%s%s", tmpchar1, tmpchar2);
    	if (strcmp("%FLAG", tmpchar1) == 0 &&strcmp("POINTERS", tmpchar2) == 0) {
    		fgets(line, MAXCHAR, fpin);
    		fgets(line, MAXCHAR, fpin);
		sscanf(line, "%d", &atomnum);
		return;
    	} 
}
}

void rtopatom() {
char str[10];
double v[5];
int numatom;
int tmpint;
int flag = -999;
int ipol = 0;

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
	if(fgets(line, MAXCHAR, fpin)==NULL) break;
    	sscanf(line, "%s%s", tmpchar1, tmpchar2);
    	if (strcmp("%FLAG", tmpchar1) == 0 &&strcmp("ATOM_NAME", tmpchar2) == 0) {
		flag = 0; 
		continue; 
    	} 
    	if (strcmp("%FLAG", tmpchar1) == 0 &&strcmp("AMBER_ATOM_TYPE", tmpchar2) == 0) {
		flag = 2; 
		continue; 
    	} 

    	if (strcmp("%FLAG", tmpchar1) == 0 &&strcmp("POLARIZABILITY", tmpchar2) == 0) {
		pol_flag = 1;
		ipol = 1;
		flag = 4; 
		continue; 
    	} 

    	if (strcmp("%FLAG", tmpchar1) == 0 &&strcmp("DIPOLE_DAMP_FACTOR", tmpchar2) == 0) {
		df_flag = 1;
		continue; 
    	} 

    	if(flag==0) {flag = 1; i= 0; continue; }
    	if(flag==2) {flag = 3; i =0 ; continue; }
    	if(flag==4) {flag = 5; i =0 ; continue; }

    	if(flag ==1) {
		for(j=0;j<20;j++) {
			strcpy(str, "");
			str[0] = '\0';
			str[0]=line[j*4];
			if(line[j*4+1] == ' ' || line[j*4+1] == '\0' || line[j*4+1] == '\n')
				str[1] = '\0';
			else {
				str[1]=line[j*4+1];
				if(line[j*4+2] == ' ' || line[j*4+2] == '\0' || line[j*4+2] == '\n')
					str[2] = '\0';
				else {
					str[2]=line[j*4+2];
					if(line[j*4+3] == ' ' || line[j*4+3] == '\0' || line[j*4+3] == '\n')
						str[3] = '\0';
					else {
						str[3]=line[j*4+3];
						str[4]= '\0';
					}
				}
			}
			strcpy(atom[i].name, str);
			i++;
			if(i>=atomnum) {
				flag = -999;
				break;
			}
		}
    	}
    	if(flag ==3) {
                for(j=0;j<20;j++) {
                        strcpy(str, "");
                        str[0] = '\0';
                        str[0]=line[j*4];
                        if(line[j*4+1] == ' ' || line[j*4+1] == '\0' || line[j*4+1] == '\n')
                                str[1] = '\0';
                        else {
                                str[1]=line[j*4+1];
                                if(line[j*4+2] == ' ' || line[j*4+2] == '\0' || line[j*4+2] == '\n')
                                        str[2] = '\0';
                                else {
                                        str[2]=line[j*4+2];
                                        if(line[j*4+3] == ' ' || line[j*4+3] == '\0' || line[j*4+3] == '\n')
                                                str[3] = '\0';
                                        else {
                                                str[3]=line[j*4+3];
                                                str[4]= '\0';
                                        }
                                }
                        }
                        strcpy(atom[i].type, str);
                        i++;
                        if(i>=atomnum) { 
                                flag = -999;
                                break;
                        }
                }

    	}
    	if(flag ==5) {
		sscanf(line, "%lf%lf%lf%lf%lf", &v[0], &v[1], &v[2], &v[3], &v[4]);
		for(j=0;j<5;j++) {
			atom[i].pol = v[j];	
			i++;
			if(i>=atomnum) {
				flag = -999;
				break;
			}
		} 
	}
}
if(ipol == 0) ireplace = 1;
if(debug == 1) {
	for(i=0;i<atomnum;i++) 
		printf("ATOM %5d %5s %5s %9.4lf\n", i+1, atom[i].name, atom[i].type, atom[i].pol);
}
}

void wtop() {
int i;
int count;
int write_flag = 1;

for(;;) {
        if(fgets(line, MAXCHAR, fpin)==NULL) break;
        sscanf(line, "%s%s", tmpchar1, tmpchar2);
        if (strcmp("%FLAG", tmpchar1) == 0 &&strcmp("POLARIZABILITY", tmpchar2) == 0) {
		fprintf(fpout, "%s\n%s\n", "%FLAG POLARIZABILITY","%FORMAT(5E16.8)");
		count = 0;
		for(i=0;i<atomnum;i++) {
			fprintf(fpout, "%16.8E", atom[i].pol);	
			count++;
			if(count == 5) {
				count = 0;
				fprintf(fpout, "\n");
			}
		}		
		if(count !=0) fprintf(fpout, "\n");

		if(idamp_factor == 1) {
			fprintf(fpout, "%s\n%s\n", "%FLAG DIPOLE_DAMP_FACTOR", "%FORMAT(5E16.8)");
			count = 0;
			for(i=0;i<atomnum;i++) {
				fprintf(fpout, "%16.8E", atom[i].damp_factor);	
				count++;
				if(count == 5) {
					count = 0;
					fprintf(fpout, "\n");
				}
			}		
			if(count !=0) fprintf(fpout, "\n");
		}
		strcpy(tmpchar1, "");
		tmpchar1[0]='\0';
		strcpy(tmpchar2, "");
		tmpchar2[0]='\0';
                write_flag = 0;
                continue;
        }
	if(write_flag == 0 && strcmp("%FLAG", tmpchar1) == 0) write_flag = 1;
        if(strcmp("%FLAG", tmpchar1) == 0 &&strcmp("DIPOLE_DAMP_FACTOR", tmpchar2) == 0) {
		strcpy(tmpchar1, "");
		tmpchar1[0]='\0';
		strcpy(tmpchar2, "");
		tmpchar2[0]='\0';
                write_flag = 0;
                continue;
	}
	if(write_flag == 0 && strcmp("%FLAG", tmpchar1) == 0) write_flag = 1;
	if(write_flag == 1) 
		fprintf(fpout, "%s", line);
}
if(pol_flag == 0) {
	fprintf(fpout, "%s\n%s\n", "%FLAG POLARIZABILITY","%FORMAT(5E16.8)");
	count = 0;
	for(i=0;i<atomnum;i++) {
		fprintf(fpout, "%16.8E", atom[i].pol);	
		count++;
		if(count == 5) {
			count = 0;
			fprintf(fpout, "\n");
		}
	}		
	if(count !=0) fprintf(fpout, "\n");
	if(idamp_factor == 1) {
		fprintf(fpout, "%s\n%s\n", "%FLAG DIPOLE_DAMP_FACTOR", "%FORMAT(5E16.8)");
		count = 0;
		for(i=0;i<atomnum;i++) {
			fprintf(fpout, "%16.8E", atom[i].damp_factor);	
			count++;
			if(count == 5) {
				count = 0;
				fprintf(fpout, "\n");
			}
		}		
		if(count !=0) fprintf(fpout, "\n");
	}
}
}

int main(int argc, char *argv[]) {
int i,j;
int suc;
int pid;
char *amberhome;
if (strcmp(COLORTEXT, "YES") == 0 || strcmp(COLORTEXT, "yes") == 0) {
	if (argc == 2
		&& (strcmp(argv[1], "-h") == 0
			|| strcmp(argv[1], "-H") == 0)) {
		printf("[31mUsage: poltop -i [0m topology file name\n"
		       "[31m              -o [0m output file name\n"
		       "[31m              -r [0m replace the existing polarizability in the topology? \n"
		       "[34m                  1:[0m yes (the default) \n"
		       "[34m                  0:[0m no  (use with care !)\n"
		       "[31m              -m [0m polarizability models, suppressed by argument '-p'\n"
		       "[34m                  A1:[0m Thole linear,      no 1-2 and 1-3, %83 of 1-4, variable screening lengths \n"
		       "[34m                  A2:[0m Thole exponential, no 1-2 and 1-3, %83 of 1-4, variable screening lengths \n"
		       "[34m                  A3:[0m Thole tinker-like, no 1-2 and 1-3, %83 of 1-4, variable screening lengths \n"
		       "[34m                  B0:[0m Applequist,        no 1-2 and 1-3, %50 of 1-4, also known as CA\n"
		       "[34m                  B1:[0m Thole linear,      no 1-2 and 1-3, %50 of 1-4, also known as CL\n"
		       "[34m                  B2:[0m Thole exponential, no 1-2 and 1-3, %50 of 1-4, also known as CE\n"
		       "[34m                  B3:[0m Thole tinker-like, no 1-2 and 1-3, %50 of 1-4, also known as CT\n"
		       "[34m                  C0:[0m Thole Applequist,  no 1-2 and 1-3, %83 of 1-4, also known as BA\n"
		       "[34m                  C1:[0m Thole linear,      no 1-2 and 1-3, %83 of 1-4, also known as BL\n"
		       "[34m                  C2:[0m Thole exponential, no 1-2 and 1-3, %83 of 1-4, also known as BE\n"
		       "[34m                  C3:[0m Thole tinker-like, no 1-2 and 1-3, %83 of 1-4, also known as BT\n"
		       "[34m                  D0:[0m Thole Applequist,  no 1-2 and 1-3, %100 of 1-4, also known as AA\n"
		       "[34m                  D1:[0m Thole linear,      no 1-2 and 1-3, %100 of 1-4, also known as AL\n"
		       "[34m                  D2:[0m Thole exponential, no 1-2 and 1-3, %100 of 1-4, also known as AE\n"
		       "[34m                  D3:[0m Thole tinker-like, no 1-2 and 1-3, %100 of 1-4, also known as AT\n"
		       "[31m              -t [0m atom type file, optional, the default is $AMBERHOME/dat/antechamber/ATCOR2.DAT\n"
		       "[31m              -p [0m parameter file, optional\n");
		exit(1);
	}
	if (argc != 7 && argc !=9 && argc != 11 && argc != 13 && argc != 15) {
		printf("[31mUsage: poltop -i [0m topology file name\n"
		       "[31m              -o [0m output file name\n"
		       "[31m              -r [0m replace the existing polarizability in the topology? \n"
		       "[34m                  1:[0m yes (the default) \n"
		       "[34m                  0:[0m no  (use with care !)\n"
		       "[31m              -m [0m polarizability models, suppressed by argument '-p'\n"
		       "[34m                  A1:[0m Thole linear,      no 1-2 and 1-3, %83 of 1-4, variable screening lengths \n"
		       "[34m                  A2:[0m Thole exponential, no 1-2 and 1-3, %83 of 1-4, variable screening lengths \n"
		       "[34m                  A3:[0m Thole tinker-like, no 1-2 and 1-3, %83 of 1-4, variable screening lengths \n"
		       "[34m                  B0:[0m Applequist,        no 1-2 and 1-3, %50 of 1-4, also known as CA\n"
		       "[34m                  B1:[0m Thole linear,      no 1-2 and 1-3, %50 of 1-4, also known as CL\n"
		       "[34m                  B2:[0m Thole exponential, no 1-2 and 1-3, %50 of 1-4, also known as CE\n"
		       "[34m                  B3:[0m Thole tinker-like, no 1-2 and 1-3, %50 of 1-4, also known as CT\n"
		       "[34m                  C0:[0m Thole Applequist,  no 1-2 and 1-3, %83 of 1-4, also known as BA\n"
		       "[34m                  C1:[0m Thole linear,      no 1-2 and 1-3, %83 of 1-4, also known as BL\n"
		       "[34m                  C2:[0m Thole exponential, no 1-2 and 1-3, %83 of 1-4, also known as BE\n"
		       "[34m                  C3:[0m Thole tinker-like, no 1-2 and 1-3, %83 of 1-4, also known as BT\n"
		       "[34m                  D0:[0m Thole Applequist,  no 1-2 and 1-3, %100 of 1-4, also known as AA\n"
		       "[34m                  D1:[0m Thole linear,      no 1-2 and 1-3, %100 of 1-4, also known as AL\n"
		       "[34m                  D2:[0m Thole exponential, no 1-2 and 1-3, %100 of 1-4, also known as AE\n"
		       "[34m                  D3:[0m Thole tinker-like, no 1-2 and 1-3, %100 of 1-4, also known as AT\n"
		       "[31m              -t [0m atom type file, optional, the default is $AMBERHOME/dat/antechamber/ATCOR2.DAT\n"
		       "[31m              -p [0m parameter file, optional\n");
		exit(1);
	}
} else {
	if (argc == 2)
	if (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "-H") == 0) {
		printf("Usage: poltop -i  topology file name\n"
		       "              -o  output file name\n"
		       "              -r  replace the existing polarizability in the topology? \n"
		       "                  1: yes (the default) \n"
		       "                  0: no  (use with care !)\n"
		       "              -m  polarizability models, suppressed by argument '-p'\n"
		       "                  A1: Thole linear,      no 1-2 and 1-3, %83 of 1-4, variable screening lengths \n"
		       "                  A2: Thole exponential, no 1-2 and 1-3, %83 of 1-4, variable screening lengths \n"
		       "                  A3: Thole tinker-like, no 1-2 and 1-3, %83 of 1-4, variable screening lengths \n"
		       "                  B0: Applequist,        no 1-2 and 1-3, %50 of 1-4, also known as CA\n"
		       "                  B1: Thole linear,      no 1-2 and 1-3, %50 of 1-4, also known as CL\n"
		       "                  B2: Thole exponential, no 1-2 and 1-3, %50 of 1-4, also known as CE\n"
		       "                  B3: Thole tinker-like, no 1-2 and 1-3, %50 of 1-4, also known as CT\n"
		       "                  C0: Thole Applequist,  no 1-2 and 1-3, %83 of 1-4, also known as BA\n"
		       "                  C1: Thole linear,      no 1-2 and 1-3, %83 of 1-4, also known as BL\n"
		       "                  C2: Thole exponential, no 1-2 and 1-3, %83 of 1-4, also known as BE\n"
		       "                  C3: Thole tinker-like, no 1-2 and 1-3, %83 of 1-4, also known as BT\n"
		       "                  D0: Thole Applequist,  no 1-2 and 1-3, %100 of 1-4, also known as AA\n"
		       "                  D1: Thole linear,      no 1-2 and 1-3, %100 of 1-4, also known as AL\n"
		       "                  D2: Thole exponential, no 1-2 and 1-3, %100 of 1-4, also known as AE\n"
		       "                  D3: Thole tinker-like, no 1-2 and 1-3, %100 of 1-4, also known as AT\n"
		       "              -t  atom type file, optional, the default is $AMBERHOME/dat/antechamber/ATCOR2.DAT\n"
		       "              -p  parameter file, optional\n");
		exit(1);
	}
	if (argc != 7 && argc !=9 && argc != 11 && argc != 13) {
		printf("Usage: poltop -i  topology file name\n"
		       "              -o  output file name\n"
		       "              -r  replace the existing polarizability in the topology? \n"
		       "                  1: yes (the default) \n"
		       "                  0: no  (use with care !)\n"
		       "              -m  polarizability models, suppressed by argument '-p'\n"
		       "                  A1: Thole linear,      no 1-2 and 1-3, %83 of 1-4, variable screening lengths \n"
		       "                  A2: Thole exponential, no 1-2 and 1-3, %83 of 1-4, variable screening lengths \n"
		       "                  A3: Thole tinker-like, no 1-2 and 1-3, %83 of 1-4, variable screening lengths \n"
		       "                  B0: Applequist,        no 1-2 and 1-3, %50 of 1-4, also known as CA\n"
		       "                  B1: Thole linear,      no 1-2 and 1-3, %50 of 1-4, also known as CL\n"
		       "                  B2: Thole exponential, no 1-2 and 1-3, %50 of 1-4, also known as CE\n"
		       "                  B3: Thole tinker-like, no 1-2 and 1-3, %50 of 1-4, also known as CT\n"
		       "                  C0: Thole Applequist,  no 1-2 and 1-3, %83 of 1-4, also known as BA\n"
		       "                  C1: Thole linear,      no 1-2 and 1-3, %83 of 1-4, also known as BL\n"
		       "                  C2: Thole exponential, no 1-2 and 1-3, %83 of 1-4, also known as BE\n"
		       "                  C3: Thole tinker-like, no 1-2 and 1-3, %83 of 1-4, also known as BT\n"
		       "                  D0: Thole Applequist,  no 1-2 and 1-3, %100 of 1-4, also known as AA\n"
		       "                  D1: Thole linear,      no 1-2 and 1-3, %100 of 1-4, also known as AL\n"
		       "                  D2: Thole exponential, no 1-2 and 1-3, %100 of 1-4, also known as AE\n"
		       "                  D3: Thole tinker-like, no 1-2 and 1-3, %100 of 1-4, also known as AT\n"
		       "              -t  atom type file, optional, the default is $AMBERHOME/dat/antechamber/ATCOR2.DAT\n"
		       "              -p  parameter file, optional\n");
		exit(1);
	}
}

ireplace = 1;
ireadat = 0;
ireadparm = 0;
strcpy(model, "A0");

for (i = 1; i < argc; i += 2) {
	if (strcmp(argv[i], "-i") == 0) 
		strcpy(ifilename, argv[i + 1]);
	if (strcmp(argv[i], "-o") == 0) 
		strcpy(ofilename, argv[i + 1]);
	if (strcmp(argv[i], "-p") == 0) {
		strcpy(pfilename, argv[i + 1]);
		ireadparm = 1;
	}
	if (strcmp(argv[i], "-t") == 0) {
		strcpy(tfilename, argv[i + 1]);
		ireadat = 1;
	}
	if (strcmp(argv[i], "-r") == 0) 
		ireplace = atoi(argv[i + 1]);
	if (strcmp(argv[i], "-m") == 0) 
		strcpy(model, argv[i + 1]);
}
if((fpin=fopen(ifilename,"r"))==NULL) {
    	fprintf(stdout, "\n Cannot open topology file: %s, exit", ifilename);
	exit(1);
}
if((fpout=fopen(ofilename,"w"))==NULL) {
	fprintf(stdout, "\n Cannot open file to write: %s, exit", ofilename);
	exit(1);
}
rtopparm();
atom = (ATOM *) malloc(sizeof(ATOM) * (atomnum + 10));
if (atom == NULL) {
        fprintf(stdout, "memory allocation error for *atom\n");
        exit(1);
}

if(ireadparm == 0) {
    amberhome = (char *) getenv("AMBERHOME");
    if( amberhome == NULL ){
       fprintf( stdout, "AMBERHOME is not set!\n" );
       exit(1);
    }
    strcpy(pfilename, amberhome);
    strcat(pfilename, "/dat/antechamber/POL.PARM"); 
    rparm();
}
else if(ireadparm == 1) rparm2();

if(ireadat == 0) {
    amberhome = (char *) getenv("AMBERHOME");
    if( amberhome == NULL ){
       fprintf( stdout, "AMBERHOME is not set!\n" );
       exit(1);
    }
    strcpy(tfilename, amberhome);
    strcat(tfilename, "/dat/antechamber/ATCOR2.DAT"); 
}
rcorr();
rewind(fpin);
rtopatom();

// find the parameters
for(i=0;i<atomnum;i++) {
	suc = 0;
	for(j=0;j<typenum;j++) {
		if(strcmp(atom[i].type, type[j].name) == 0) {
			pid = type[j].pid;
			if(ireplace == 1) atom[i].pol = parm[pid].pol;	
			atom[i].damp_factor = parm[pid].damp_factor;	
			suc = 1;
			break;
		}
	}
	if(suc == 0) {
		fprintf(stdout, "The atom type (%s) does not show up in the atom type file %s !\n", atom[i].type, tfilename);
		exit(0);
	}
	if(debug == 1) {
		printf("ATOM %d %5s %5s %9.5lf %9.5lf\n", i+1, atom[i].name, atom[i].type, atom[i].pol, atom[i].damp_factor);
	}
}

if(atom[0].damp_factor < -999990) idamp_factor = 0;

rewind(fpin);
wtop();
fclose(fpin);
fclose(fpout);
return(0);
}

