# include <stdio.h>
# include <math.h>
# include <ctype.h>
# include <stdlib.h>
# include <string.h>
# define MAXCHAR 256 
# define COLORTEXT "YES" 
# define PI 3.1415926
# define debug 0
typedef struct {
	int  flag;
	int  id;
	char name[10];
	char type[10];
	double mass;
	double pol;
} ATOM;
typedef struct {
	int  id;
	int flag;
	char at1[10];
	char at2[10];
	double len;
	double force;
} BOND;
typedef struct {
	int id;
	int flag;
	char at1[10];
	char at2[10];
	char at3[10];
	double angle;
	double force;
} ANGLE;
typedef struct {
	int id;
	int type;
	int flag;
	char at1[10];
	char at2[10];
	char at3[10];
	char at4[10];
	double periodicity;
	double phase;
	double force;
} TOR;

/* for vdw parameters */
typedef struct {
	char type[10];
	double e;
	double r;
} VDW;

/* for atom-atom pairs*/
typedef struct {
	int index;
	int ati;
	int atj;
	double a;
	double b;
} VDWPAIR;

/* for mixed information of ATOM and VDW*/
typedef struct {
        char type[10];
        double mass;
        double pol;
        double e;
        double r;
} VDW2;

int i,j,k=0;
int atomnum = 0;
int typenum = 0; /* vdw type number */
int pairnum = 0;
int btnum = 0;
int atnum = 0;
int ttnum = 0;

int btnum_count = 0;
int atnum_count = 0;
int ttnum_count = 0;
int improper_flag = 0;

ATOM *atom;
BOND *bond;
ANGLE *angle;
TOR *tor;
VDW *vdw;
VDWPAIR *vdwpair;

BOND *bond2;
ANGLE *angle2;
TOR *tor2;
VDW2 *vdw2;

int atomtypenum = 0; /*true atom type number */
int btnum2 = 0; /* ture bond type number */
int atnum2 = 0; /* true angle type number*/
int ttnum2 = 0; /* true torsional type number */
int nbonh = 0;
int nbona = 0;
int ntheth = 0;
int ntheta = 0;
int nphih = 0;
int nphia = 0;


char pfilename[MAXCHAR];
char ofilename[MAXCHAR];
FILE *fptop,  *fpout ;

char line[MAXCHAR];
char tmpchar1[MAXCHAR], tmpchar2[MAXCHAR], tmpchar3[MAXCHAR];
int tmpint1, tmpint2;


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
	if(fgets(line, MAXCHAR, fptop)==NULL) break;
    	sscanf(line, "%s%s", tmpchar1, tmpchar2);
    	if (strcmp("%FLAG", tmpchar1) == 0 &&strcmp("POINTERS", tmpchar2) == 0) {
    		fgets(line, MAXCHAR, fptop);
    		fgets(line, MAXCHAR, fptop);
		sscanf(line, "%d%d%d%d%d%d%d%d", &atomnum, &typenum, &nbonh, &nbona, &ntheth, &ntheta, &nphih, &nphia);
    		fgets(line, MAXCHAR, fptop);
		sscanf(line, "%d%d%d%d%d%d%d%d", &tmpint1, &tmpint2, &tmpint3, &tmpint4, &tmpint5, &btnum, &atnum, &ttnum);
    	} 
}
}

void rtopatom() {
char str[10];
double v[5];
int   index[10];
int numatom;
int tmpint;
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
	if(fgets(line, MAXCHAR, fptop)==NULL) break;
    	sscanf(line, "%s%s", tmpchar1, tmpchar2);
    	if (strcmp("%FLAG", tmpchar1) == 0 &&strcmp("ATOM_NAME", tmpchar2) == 0) {
		flag = 0; 
		continue; 
    	} 
    	if (strcmp("%FLAG", tmpchar1) == 0 &&strcmp("MASS", tmpchar2) == 0) {
		flag = 2; 
		continue; 
    	} 
    	if (strcmp("%FLAG", tmpchar1) == 0 &&strcmp("ATOM_TYPE_INDEX", tmpchar2) == 0) {
		flag = 4; 
		continue; 
    	} 
    	if (strcmp("%FLAG", tmpchar1) == 0 &&strcmp("AMBER_ATOM_TYPE", tmpchar2) == 0) {
		flag = 6; 
		continue; 
    	} 

    	if (strcmp("%FLAG", tmpchar1) == 0 &&strcmp("POLARIZABILITY", tmpchar2) == 0) {
		flag = 8; 
		continue; 
    	} 

    	if(flag==0) {flag = 1; i= 0; continue; }
    	if(flag==2) {flag = 3; i =0 ; continue; }
    	if(flag==4) {flag = 5; i =0 ; continue; }
    	if(flag==6) {flag = 7; i =0 ; continue; }
    	if(flag==8) {flag = 9; i =0 ; continue; }

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
		sscanf(line, "%lf%lf%lf%lf%lf", &v[0], &v[1], &v[2], &v[3], &v[4]);
		for(j=0;j<5;j++) {
			atom[i].mass = v[j];	
			i++;
			if(i>=atomnum) {
				flag = -999;
				break;
			}
		} 
	}
    	if(flag ==5) {
		sscanf(line, "%d%d%d%d%d%d%d%d%d%d", &index[0], &index[1], &index[2], &index[3], &index[4],
                                                     &index[5], &index[6], &index[7], &index[8], &index[9]);
		for(j=0;j<10;j++) {
			atom[i].id = index[j];	
			i++;
			if(i>=atomnum) {
				flag = -999;
				break;
			}
		} 
	}
    	if(flag ==7) {
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
    	if(flag ==9) {
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
if(debug == 1) {
	for(i=0;i<atomnum;i++) 
		printf("ATOM %5d %5s %5s %9.4lf %9.4lf\n", i+1, atom[i].name, atom[i].type, atom[i].mass, atom[i].pol);
}
}

void rtopff() {
int tmpint;
int flag = -999;
int tortype;
double v[5];
int at1, at2, at3, at4;

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
	if(fgets(line, MAXCHAR, fptop)==NULL) break;
    	sscanf(line, "%s%s", tmpchar1, tmpchar2);
    	if (strcmp("%FLAG", tmpchar1) == 0 &&strcmp("BOND_FORCE_CONSTANT", tmpchar2) == 0) {
		flag = 0; 
		continue; 
    	} 
    	if (strcmp("%FLAG", tmpchar1) == 0 &&strcmp("BOND_EQUIL_VALUE", tmpchar2) == 0) {
		flag = 2; 
		continue; 
    	} 
    	if (strcmp("%FLAG", tmpchar1) == 0 &&strcmp("ANGLE_FORCE_CONSTANT", tmpchar2) == 0) {
		flag = 4; 
		continue; 
    	} 
    	if (strcmp("%FLAG", tmpchar1) == 0 &&strcmp("ANGLE_EQUIL_VALUE", tmpchar2) == 0) {
		flag = 6; 
		continue; 
    	} 
    	if (strcmp("%FLAG", tmpchar1) == 0 &&strcmp("DIHEDRAL_FORCE_CONSTANT", tmpchar2) == 0) {
		flag = 8; 
		continue; 
    	} 
    	if (strcmp("%FLAG", tmpchar1) == 0 &&strcmp("DIHEDRAL_PERIODICITY", tmpchar2) == 0) {
		flag = 10; 
		continue; 
    	} 
    	if (strcmp("%FLAG", tmpchar1) == 0 &&strcmp("DIHEDRAL_PHASE", tmpchar2) == 0) {
		flag = 12; 
		continue; 
    	} 
    	if (strcmp("%FLAG", tmpchar1) == 0 &&strcmp("BONDS_INC_HYDROGEN", tmpchar2) == 0) {
		flag = 14; 
		continue; 
    	} 
    	if (strcmp("%FLAG", tmpchar1) == 0 &&strcmp("BONDS_WITHOUT_HYDROGEN", tmpchar2) == 0) {
		flag = 16; 
		continue; 
    	} 
    	if (strcmp("%FLAG", tmpchar1) == 0 &&strcmp("ANGLES_INC_HYDROGEN", tmpchar2) == 0) {
		flag = 18; 
		continue; 
    	} 
    	if (strcmp("%FLAG", tmpchar1) == 0 &&strcmp("ANGLES_WITHOUT_HYDROGEN", tmpchar2) == 0) {
		flag = 20; 
		continue; 
    	} 
    	if (strcmp("%FLAG", tmpchar1) == 0 &&strcmp("DIHEDRALS_INC_HYDROGEN", tmpchar2) == 0) {
		flag = 22; 
		continue; 
    	} 
    	if (strcmp("%FLAG", tmpchar1) == 0 &&strcmp("DIHEDRALS_WITHOUT_HYDROGEN", tmpchar2) == 0) {
		flag = 24; 
		continue; 
    	} 
    	if (strcmp("%FLAG", tmpchar1) == 0 &&strcmp("LENNARD_JONES_ACOEF", tmpchar2) == 0) {
		flag = 26; 
		continue; 
    	} 
    	if (strcmp("%FLAG", tmpchar1) == 0 &&strcmp("LENNARD_JONES_BCOEF", tmpchar2) == 0) {
		flag = 28; 
		continue; 
    	} 

    	if(flag==0) {flag = 1; i =0 ; continue; }
    	if(flag==2) {flag = 3; i =0 ; continue; }
    	if(flag==4) {flag = 5; i =0 ; continue; }
    	if(flag==6) {flag = 7; i =0 ; continue; }
    	if(flag==8) {flag = 9; i =0 ; continue; }
    	if(flag==10) {flag = 11; i =0 ; continue; }
    	if(flag==12) {flag = 13; i =0 ; continue; }
    	if(flag==14) {flag = 15; i =0 ; k=0;  continue; }
    	if(flag==16) {flag = 17; i =0 ; k=0;  continue; }
    	if(flag==18) {flag = 19; i =0 ; k=0;  continue; }
    	if(flag==20) {flag = 21; i =0 ; k=0;  continue; }
    	if(flag==22) {flag = 23; i =0 ; k=0;  continue; }
    	if(flag==24) {flag = 25; i =0 ; k=0;  continue; }
    	if(flag==26) {flag = 27; i =0 ; k=0;  continue; }
    	if(flag==28) {flag = 29; i =0 ; k=0;  continue; }

        if(flag ==1) {
                sscanf(line, "%lf%lf%lf%lf%lf", &v[0], &v[1], &v[2], &v[3], &v[4]);
                for(j=0;j<5;j++) {
                        bond[i].force = v[j];
                        i++;
                        if(i>=btnum) {
                                flag = -999;
                                break;
                        }
                }
        }

        if(flag ==3) {
                sscanf(line, "%lf%lf%lf%lf%lf", &v[0], &v[1], &v[2], &v[3], &v[4]);
                for(j=0;j<5;j++) {
                        bond[i].len = v[j];
			bond[i].flag = 0;
                        i++;
                        if(i>=btnum) {
                                flag = -999;
                                break;
                        }
                }
        }

        if(flag ==5) {
                sscanf(line, "%lf%lf%lf%lf%lf", &v[0], &v[1], &v[2], &v[3], &v[4]);
                for(j=0;j<5;j++) {
                        angle[i].force = v[j];
                        i++;
                        if(i>=atnum) {
                                flag = -999;
                                break;
                        }
                }
        }

        if(flag ==7) {
                sscanf(line, "%lf%lf%lf%lf%lf", &v[0], &v[1], &v[2], &v[3], &v[4]);
                for(j=0;j<5;j++) {
                        angle[i].angle = v[j];
			angle[i].flag = 0;
                        i++;
                        if(i>=atnum) {
                                flag = -999;
                                break;
                        }
                }
        }

        if(flag ==9) {
                sscanf(line, "%lf%lf%lf%lf%lf", &v[0], &v[1], &v[2], &v[3], &v[4]);
                for(j=0;j<5;j++) {
                        tor[i].force = v[j];
                        i++;
                        if(i>=ttnum) {
                                flag = -999;
                                break;
                        }
                }
        }

        if(flag ==11) {
                sscanf(line, "%lf%lf%lf%lf%lf", &v[0], &v[1], &v[2], &v[3], &v[4]);
                for(j=0;j<5;j++) {
                        tor[i].periodicity = v[j];
                        i++;
                        if(i>=ttnum) {
                                flag = -999;
                                break;
                        }
                }
        }
        if(flag ==13) {
                sscanf(line, "%lf%lf%lf%lf%lf", &v[0], &v[1], &v[2], &v[3], &v[4]);
                for(j=0;j<5;j++) {
                        tor[i].phase = v[j];
			tor[i].flag = -1;
                        i++;
                        if(i>=ttnum) {
                                flag = -999;
                                break;
                        }
                }
        }
        if(flag ==15 || flag == 17) {
                j = 0;
                while(j <= 10) {
                        j++;
                        if(line[j*8-1] <='9' && line[j*8-1]>='0')
                                sscanf(&line[8*j-8], "%d", &tmpint);
                        else break;
                        if(k==0) {
                                at1 = tmpint/3;
                                k=1;
                                continue;
                        }
                        if(k==1) {
                                at2 = tmpint/3;
                                k=2;
                                continue;
                        }
                        if(k==2) {
                                k=0;
				strcpy(bond2[btnum2].at1, atom[at1].type);	
				strcpy(bond2[btnum2].at2, atom[at2].type);	
				bond2[btnum2].id = tmpint - 1;
				bond2[btnum2].flag = 0;
				btnum2++;
                        	if(btnum2>=(nbonh + nbona)) {
                                	flag = -999;
                                	break;
                        	}
                                continue;
                        }
                } 
	}
	
        if(flag ==19 || flag == 21) {
                j = 0;
                while(j <= 10) {
                        j++;
                        if(line[j*8-1] <='9' && line[j*8-1]>='0')
                                sscanf(&line[8*j-8], "%d", &tmpint);
                        else break;
                        if(k==0) {
                                at1 = tmpint/3;
                                k=1;
                                continue;
                        }
                        if(k==1) {
                                at2 = tmpint/3;
                                k=2;
                                continue;
                        }
                        if(k==2) {
                                at3 = tmpint/3;
                                k=3;
                                continue;
                        }
                        if(k==3) {
                                k=0;
                                strcpy(angle2[atnum2].at1, atom[at1].type);
                                strcpy(angle2[atnum2].at2, atom[at2].type);
                                strcpy(angle2[atnum2].at3, atom[at3].type);
				angle2[atnum2].id = tmpint -1;
				angle2[atnum2].flag = 0;
				atnum2 ++;
                        	if(atnum2>=(ntheth + ntheta)) {
                                	flag = -999;
                                	break;
                        	}
                        }
                }
        }

        if(flag ==23 || flag == 25) {
                j = 0;
		tortype = 0;
                while(j <= 10) {
                        j++;
                        if(line[j*8-1] <='9' && line[j*8-1]>='0')
                                sscanf(&line[8*j-8], "%d", &tmpint);
                        else break;
                        if(k==0) {
                                at1 = fabs(tmpint)/3;
                                k=1;
                                continue;
                        }
                        if(k==1) {
                                at2 = fabs(tmpint)/3;
                                k=2;
                                continue;
                        }
                        if(k==2) {
                                at3 = fabs(tmpint)/3;
                                k=3;
                                continue;
                        }
                        if(k==3) {
				if(tmpint < 0) tortype = 1;
                                at4 = fabs(tmpint)/3;
                                k=4;
                                continue;
                        }
                        if(k==4) {
                                k=0;
                                strcpy(tor2[ttnum2].at1, atom[at1].type);
                                strcpy(tor2[ttnum2].at2, atom[at2].type);
                                strcpy(tor2[ttnum2].at3, atom[at3].type);
                                strcpy(tor2[ttnum2].at4, atom[at4].type);
				tor2[ttnum2].type = tortype;
				tor2[ttnum2].id = tmpint - 1;
				tor2[ttnum2].flag = 0;
				ttnum2 ++;
                        	if(ttnum2>=(nphih + nphia)) {
                                	flag = -999;
                                	break;
                        	}
                                continue;
                        }
                }
        }
        if(flag ==27) {
                sscanf(line, "%lf%lf%lf%lf%lf", &v[0], &v[1], &v[2], &v[3], &v[4]);
                for(j=0;j<5;j++) {
                        vdwpair[i].a = v[j];
                        i++;
                        if(i>=pairnum) {
                                flag = -999;
                                break;
                        }
                }
        }
        if(flag ==29) {
                sscanf(line, "%lf%lf%lf%lf%lf", &v[0], &v[1], &v[2], &v[3], &v[4]);
                for(j=0;j<5;j++) {
                        vdwpair[i].b = v[j];
                        i++;
                        if(i>=pairnum) {
                                flag = -999;
                                break;
                        }
                }
        }
}
}

void write(void) {
int i,j;
int suc;
fprintf(fpout, "remark goes here\n");
fprintf(fpout, "MASS\n");
for(i=0;i<atomtypenum;i++) 
	fprintf(fpout, "%-4s%10.3lf%10.3lf\n", vdw2[i].type, vdw2[i].mass, vdw2[i].pol);
fprintf(fpout, "\n");
fprintf(fpout, "BOND\n");
for(i=0;i<btnum2;i++){
	if(bond2[i].flag == 0) continue;
	if(strcmp(bond2[i].at1, bond2[i].at2) > 0) 
		fprintf(fpout, "%-2s-%-2s%10.3lf%10.3lf\n", bond2[i].at2, bond2[i].at1, bond2[i].force, bond2[i].len);
	else
		fprintf(fpout, "%-2s-%-2s%10.3lf%10.3lf\n", bond2[i].at1, bond2[i].at2, bond2[i].force, bond2[i].len);
}

fprintf(fpout, "\n");
fprintf(fpout, "ANGLE\n");
for(i=0;i<atnum2;i++){
	if(angle2[i].flag == 0) continue;
	if(strcmp(angle2[i].at1, angle2[i].at3) > 0) 
		fprintf(fpout, "%-2s-%-2s-%-2s%10.3lf%10.3lf\n", angle2[i].at3, angle2[i].at2, angle2[i].at1, 
				angle2[i].force, angle2[i].angle * 180.0/PI);
	else
		fprintf(fpout, "%-2s-%-2s-%-2s%10.3lf%10.3lf\n", angle2[i].at1, angle2[i].at2, angle2[i].at3, 
				angle2[i].force, angle2[i].angle * 180.0/PI);
}

fprintf(fpout, "\n");
fprintf(fpout, "DIHE\n");
for(i=0;i<ttnum2;i++){
	if(tor2[i].flag == 0) continue;
	if(tor2[i].type == 1) continue;
	suc = 0;
	if(i<ttnum2-1 && tor2[i+1].type != 1) {
		if(strcmp(tor2[i].at1, tor2[i+1].at1) == 0 &&
		   strcmp(tor2[i].at2, tor2[i+1].at2) == 0 &&
		   strcmp(tor2[i].at3, tor2[i+1].at3) == 0 &&
		   strcmp(tor2[i].at4, tor2[i+1].at4) == 0)
			suc = 1;
		if(suc == 0) {
			if(strcmp(tor2[i].at1, tor2[i+1].at4) == 0 &&
		   	   strcmp(tor2[i].at2, tor2[i+1].at3) == 0 &&
		           strcmp(tor2[i].at3, tor2[i+1].at2) == 0 &&
		           strcmp(tor2[i].at4, tor2[i+1].at1) == 0)
					suc = 1;
		}
		if(tor2[i+1].flag == 0) suc = 0;
		if(tor2[i+1].type == 1) suc = 0;
	}
	if(strcmp(tor2[i].at1, tor2[i].at4) > 0) {
		if(suc == 1) 
			fprintf(fpout, "%-2s-%-2s-%-2s-%-2s%4d%15.9lf%15.3lf%15.3lf\n", tor2[i].at4, tor2[i].at3, tor2[i].at2, tor2[i].at1, 1, 
						tor2[i].force, tor2[i].phase*180/PI, -tor2[i].periodicity);
		else
			fprintf(fpout, "%-2s-%-2s-%-2s-%-2s%4d%15.9lf%15.3lf%15.3lf\n", tor2[i].at4, tor2[i].at3, tor2[i].at2, tor2[i].at1, 1, 
						tor2[i].force, tor2[i].phase*180/PI, tor2[i].periodicity);
	}
	else {
		if(suc == 1)
			fprintf(fpout, "%-2s-%-2s-%-2s-%-2s%4d%15.9lf%15.3lf%15.3lf\n", tor2[i].at1, tor2[i].at2, tor2[i].at3, tor2[i].at4, 1, 
						tor2[i].force, tor2[i].phase*180/PI, -tor2[i].periodicity);
		else
			fprintf(fpout, "%-2s-%-2s-%-2s-%-2s%4d%15.9lf%15.3lf%15.3lf\n", tor2[i].at1, tor2[i].at2, tor2[i].at3, tor2[i].at4, 1, 
						tor2[i].force, tor2[i].phase*180/PI, tor2[i].periodicity);
	}
}
fprintf(fpout, "\n");
fprintf(fpout, "IMPROPER\n");
if(improper_flag == 1) {
	for(i=0;i<ttnum2;i++){
		if(tor2[i].flag == 0) continue;
		if(tor2[i].type == 0) continue;
		fprintf(fpout, "%-2s-%-2s-%-2s-%-2s%4d%15.3lf%15.3lf%15.3lf\n", tor2[i].at1, tor2[i].at2, tor2[i].at3, tor2[i].at4, 1, 
					tor2[i].force, tor2[i].phase*180/PI, tor2[i].periodicity);
	}
}
fprintf(fpout, "\n");
fprintf(fpout, "NONBON\n");
for(i=0;i<atomtypenum;i++)
	fprintf(fpout, "  %-2s%16.4lf%8.4lf\n", vdw2[i].type, vdw2[i].r, vdw2[i].e);
fprintf(fpout, "\n\n");

}


int main(int argc, char *argv[]) {
int i,j,m,n;
int suc;
int suc2;
int count;
int begin, end;
char atseq[3][MAXCHAR];
char tmpchar[MAXCHAR];
/*  readin a topology file and crd file,then save out a mol2 file */
if (strcmp(COLORTEXT, "YES") == 0 || strcmp(COLORTEXT, "yes") == 0) {
	if (argc == 2
		&& (strcmp(argv[1], "-h") == 0
			|| strcmp(argv[1], "-H") == 0)) {
		printf("[31mUsage: top2ff -i [0m topology file name (input)\n"
		       "[31m              -o [0m output file name)\n"
		       "[31m              -f [0m flag of output improper dihedral angle parameters\n"
		       "[32m                  1[0m Yes\n"
		       "[32m                  0[0m No, the default\n");
		exit(1);
	}
	if (argc != 5 && argc != 7) {
		printf("[31mUsage: top2ff -i [0m topology file name (input)\n"
		       "[31m              -o [0m output file name)\n"
		       "[31m              -f [0m flag of output improper dihedral angle parameters\n"
		       "[32m                  1[0m Yes\n"
		       "[32m                  0[0m No, the default\n");
		exit(1);
	}
} else {
	if (argc == 2)
	if (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "-H") == 0) {
		printf("Usage: top2ff -i  topology file name (input)\n"
		       "              -o  output file name)\n"
		       "              -f  flag of output improper dihedral angle parameters\n"
		       "                  1 Yes\n"
		       "                  0 No, the default\n");
		exit(1);
	}
	if (argc != 5 && argc != 7) {
		printf("Usage: top2ff -i  topology file name (input)\n"
		       "              -o  output file name)\n"
		       "              -f  flag of output improper dihedral angle parameters\n"
		       "                  1 Yes\n"
		       "                  0 No, the default\n");
		exit(1);
	}
}

for (i = 1; i < argc; i += 2) {
	if (strcmp(argv[i], "-i") == 0) 
		strcpy(pfilename, argv[i + 1]);
	if (strcmp(argv[i], "-o") == 0) 
		strcpy(ofilename, argv[i + 1]);
	if (strcmp(argv[i], "-f") == 0) 
		improper_flag = atoi (argv[i+1]); 
}
if(improper_flag != 0 && improper_flag != 1) 
	improper_flag = 0;
if((fptop=fopen(pfilename,"r"))==NULL) {
    	fprintf(stdout, "\n Cannot open topology file: %s, exit", pfilename);
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
vdw = (VDW *) malloc(sizeof(VDW) * (typenum + 10));
if (vdw == NULL) {
        fprintf(stdout, "memory allocation error for *vdw\n");
        exit(1);
}
vdwpair = (VDWPAIR *) malloc(sizeof(VDWPAIR) * (typenum*typenum + 10));
if (vdwpair == NULL) {
        fprintf(stdout, "memory allocation error for *vdwpair\n");
        exit(1);
}

bond = (BOND *) malloc(sizeof(BOND) * (btnum + 10));
if (bond == NULL) {
        fprintf(stdout, "memory allocation error for *bond\n");
        exit(1);
}
angle = (ANGLE *) malloc(sizeof(ANGLE) * (atnum + 10));
if (angle == NULL) {
        fprintf(stdout, "memory allocation error for *angle\n");
        exit(1);
}
tor = (TOR *) malloc(sizeof(TOR) * (ttnum + 10));
if (tor == NULL) {
        fprintf(stdout, "memory allocation error for *tor\n");
        exit(1);
}

btnum2 = nbonh + nbona; 
atnum2 = ntheth + ntheta; 
ttnum2 = nphih + nphia; 

bond2 = (BOND *) malloc(sizeof(BOND) * (btnum2 + 10));
if (bond2 == NULL) {
        fprintf(stdout, "memory allocation error for *bond\n");
        exit(1);
}
angle2 = (ANGLE *) malloc(sizeof(ANGLE) * (atnum2 + 10));
if (angle2 == NULL) {
        fprintf(stdout, "memory allocation error for *angle\n");
        exit(1);
}
tor2 = (TOR *) malloc(sizeof(TOR) * (ttnum2 + 10));
if (tor2 == NULL) {
        fprintf(stdout, "memory allocation error for *tor\n");
        exit(1);
}

btnum2 = 0;
atnum2 = 0;
ttnum2 = 0;

/* generate vdw pair index, this index is a function of type number*/
count = 0;
for(i=1; i<= typenum; i++)
	for(j=1; j<= i; j++) {
		vdwpair[count].index = count;
		vdwpair[count].ati   = i - 1;
		vdwpair[count].atj   = j - 1;
		count++;
	}
pairnum = count;

rewind(fptop);
rtopatom();
rewind(fptop);
rtopff();

fclose(fptop);

/* 
   determine vdw's radius and well depth parameters when the two atoms 
   of a pair have the same atom type
   there are 'count' distinct vdw parameters, ie cout = typenum
*/
count = 0;
for(i=0; i<pairnum;i++) {
	if(vdwpair[i].ati == vdwpair[i].atj) {
/*
        aij=vdw[i].e * pow(vdw[i].r*2, 12);
        bij=2*vdw[i].e * pow(vdw[i].r*2, 6);
*/		
		if(vdwpair[i].b <= 0) {
			vdw[count].r = 0;
			vdw[count].e = 0;
		}
		else {
			vdw[count].r = pow(2*vdwpair[i].a/vdwpair[i].b, 1/6.0);
			vdw[count].e = vdwpair[i].a / pow(vdw[count].r, 12);
			vdw[count].r *= 0.5;
		}
		count++;
	}
}

/*find out how many amber atom types */
count = 0;
for(i=0;i<atomnum;i++) {
	suc = 1;
	for(j=0;j<i;j++) {
		if(strcmp(atom[i].type, atom[j].type) == 0) {
			suc = 0;
			break;
		}
	}
	if(suc == 1) 
		count++;
}

vdw2 = (VDW2 *) malloc(sizeof(VDW2) * (count + 10));
if (vdw2 == NULL) {
        fprintf(stdout, "memory allocation error for *vdw2\n");
        exit(1);
}

/* assign values from arrays of 'atom' and 'vdw' to 'vdw2'*/
count = 0;
for(i=0;i<atomnum;i++) {
	suc = 1;
	for(j=0;j<i;j++) {
		if(strcmp(atom[i].type, atom[j].type) == 0) {
			suc = 0;
			break;
		}
	}
	if(suc == 1) {
		strcpy(vdw2[count].type, atom[i].type);
		vdw2[count].mass = atom[i].mass;
		vdw2[count].pol  = atom[i].pol;
		vdw2[count].r    = vdw[atom[i].id-1].r;
		vdw2[count].e    = vdw[atom[i].id-1].e;
		count++;
	}
}
atomtypenum = count;

/* find true bond types */ 
for(i=0;i<btnum2;i++) {
	suc = 1;
	for(j=0;j<i;j++) {
		if(strcmp(bond2[i].at1, bond2[j].at1) == 0 && strcmp(bond2[i].at2, bond2[j].at2) == 0) {
			suc = 0;
			break;
		}
		if(strcmp(bond2[i].at1, bond2[j].at2) == 0 && strcmp(bond2[i].at2, bond2[j].at1) == 0) {
			suc = 0;
			break;
		}
	}
	if(suc == 1) {
		bond2[i].flag  = 1;
		bond2[i].force = bond[bond2[i].id].force;
		bond2[i].len   = bond[bond2[i].id].len;
	}
}

/* find true angle types */ 
for(i=0;i<atnum2;i++) {
	suc = 1;
	for(j=0;j<i;j++) 
		if(strcmp(angle2[i].at2, angle2[j].at2) == 0) {
			if(strcmp(angle2[i].at1, angle2[j].at1) == 0 && strcmp(angle2[i].at3, angle2[j].at3) == 0) {
				suc = 0;
				break;
			}
			if(strcmp(angle2[i].at1, angle2[j].at3) == 0 && strcmp(angle2[i].at3, angle2[j].at1) == 0) {
				suc = 0;
				break;
			}
		}
	if(suc == 1) {
		angle2[i].flag  = 1;
		angle2[i].force = angle[angle2[i].id].force;
		angle2[i].angle = angle[angle2[i].id].angle;
	}
}
/* find true torsional angle types*/
for(i=0;i<ttnum2;i++) {
	if(tor2[i].type == 1) continue;
	suc = 1;
	for(j=0;j<i;j++) {
		if(strcmp(tor2[i].at2, tor2[j].at2) == 0 && strcmp(tor2[i].at3, tor2[j].at3) == 0 &&
		   strcmp(tor2[i].at1, tor2[j].at1) == 0 && strcmp(tor2[i].at4, tor2[j].at4) == 0) {
			suc = 0;
			break;
		}
		if(strcmp(tor2[i].at2, tor2[j].at3) == 0 && strcmp(tor2[i].at3, tor2[j].at2) == 0 &&
		   strcmp(tor2[i].at1, tor2[j].at4) == 0 && strcmp(tor2[i].at4, tor2[j].at1) == 0) {
			suc = 0;
			break;
		}
	}
	if(suc == 1) {
		tor2[i].flag  = 1;
		tor2[i].force = tor[tor2[i].id].force;
		tor2[i].phase = tor[tor2[i].id].phase;
		tor2[i].periodicity= tor[tor2[i].id].periodicity;
		begin = i;
		while(1) {
			if(i>= ttnum2) break;
			suc2 = 0;
			if(strcmp(tor2[i].at1, tor2[i+1].at1) == 0 &&
		   	   strcmp(tor2[i].at2, tor2[i+1].at2) == 0 &&
		           strcmp(tor2[i].at3, tor2[i+1].at3) == 0 &&
		           strcmp(tor2[i].at4, tor2[i+1].at4) == 0) 
				suc2 = 1;
			if(suc2 == 0) {
				if(strcmp(tor2[i].at1, tor2[i+1].at4) == 0 &&
		   	  	   strcmp(tor2[i].at2, tor2[i+1].at3) == 0 &&
		           	   strcmp(tor2[i].at3, tor2[i+1].at2) == 0 &&
		           	   strcmp(tor2[i].at4, tor2[i+1].at1) == 0) 
				suc2 = 1;
			}
/* as i and i+1 belong to the same torsional angle, their periodicity must be different*/
			for(j= begin; j<= i; j++)
				if(tor2[j].periodicity == tor[tor2[i+1].id].periodicity) {
					suc2 = 0;  
					break;
				}
			if(suc2 == 0) break;

			i++;
			tor2[i].flag  = 1;
			tor2[i].force = tor[tor2[i].id].force;
			tor2[i].phase = tor[tor2[i].id].phase;
			tor2[i].periodicity= tor[tor2[i].id].periodicity;
		}
	}
}
/*
for(i=0;i<ttnum2;i++) {
	printf("%5d %5d %5d %5s %5s %5s %5s %5.0lf %5.2lf %7.2lf\n", i+1, tor2[i].flag, tor2[i].id, tor2[i].at1, tor2[i].at2, tor2[i].at3, 
}
*/

/* find true improper torsional angle types*/
if(improper_flag == 1) {
	for(i=0;i<ttnum2;i++) {
		if(tor2[i].type == 0) continue;
		strcpy(atseq[0], tor2[i].at1); 
		strcpy(atseq[1], tor2[i].at2); 
		strcpy(atseq[2], tor2[i].at4); 
		for(m=0;m<2;m++)
			for(n=1;n<3;n++) 
				if(strcmp(atseq[m], atseq[n]) > 0) {
					strcpy(tmpchar, atseq[m]);
					strcpy(atseq[m], atseq[n]);
					strcpy(atseq[n], tmpchar);
				}
		strcpy(tor2[i].at1, atseq[0]);
		strcpy(tor2[i].at2, atseq[1]);
		strcpy(tor2[i].at4, atseq[2]);

		suc = 1;
		for(j=0;j<i;j++) 
			if(strcmp(tor2[i].at1, tor2[j].at1) == 0 && 
		   	strcmp(tor2[i].at2, tor2[j].at2) == 0 && 
		   	strcmp(tor2[i].at3, tor2[j].at3) == 0 && 
		   	strcmp(tor2[i].at4, tor2[j].at4) == 0) {
					suc = 0;
					break;
				}
		if(suc == 1) {
			tor2[i].flag  = 1;
			tor2[i].force = tor[tor2[i].id].force;
			tor2[i].phase = tor[tor2[i].id].phase;
			tor2[i].periodicity= tor[tor2[i].id].periodicity;
		}
	}
}
write();
fclose(fpout);
return(0);
}


