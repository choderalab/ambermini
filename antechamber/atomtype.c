/*
************************************************************************
*           All Copyright Reserved!                                    *
*                                                                      *
*  Prog:    atomtype                                                   *
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
# include "ring.c"
# include "rotate.c"
# include "ac.c"
# include "mol2.c"
# define debug 0

RING *ring;
int ringnum = 0;
MOLINFO minfo;
CONTROLINFO cinfo;

CHEMENV ces[MAXCES][MAXBEED];
int ceslen[MAXCES];
int cesindex[MAXBEED];
int ces_bond_index = 0;
int ces_bond_num = 0;
int select_chain_index[MAXCES];
char ces_bond_at1[MAX_CES_BOND][10];
char ces_bond_at2[MAX_CES_BOND][10];
char ces_bond_type[MAX_CES_BOND][10];
NAME *atomcesname;

ATOM *atom;
BOND *bond;
AROM *arty;
int atomnum = 0;
int bondnum = 0;

int *selectchain;
int *selectindex;
int *sb;
int *SB;
int *db;
int *DB;
int *DL;
int *tb;
int *TB;
int *AB;
int *nr;

int caid = 0;
FILE *fpin;
FILE *fpparm;
char ifilename[MAXCHAR];
char ofilename[MAXCHAR];
char atomtype_def[MAXCHAR] = "ATOMTYPR_GFF.DEF";
int i, j, k, l;
int selectnum = 0;
int chain = 0;
int pindex = 0;
int findex = 0;
int anindex = 0;
int suc_flag = 0;
int maxchain = 0;
int ipostadjustment = 0;


typedef struct {
	int num;
	int connum[30];
	int atomicnum[30];
	char wildname[10];
} WILDATOM;

WILDATOM wildatom[MAXWILDATOM];
int wenum;
int weindex = 1;

typedef struct {
	int num;
	int id;
	int at[MAXBEED];
} SCHAIN;

SCHAIN schain[MAXSCHAIN];
int ssindex[MAXSCHAIN];
int scnum;



void assignwildatom(int id1, int id2, char *str)
{
	int i;
	int num1, num2;
	char name[10];
	char connum[10];
	num1 = 0;
	num2 = 0;
	for (i = 0; i < strlen(str); i++) {
		if ((str[i] <= 'Z' && str[i] >= 'A')
			|| (str[i] <= 'z' && str[i] >= 'a'))
			name[num1++] = str[i];
		if ((str[i] <= '9' && str[i] >= '0'))
			connum[num2++] = str[i];
	}
	name[num1] = '\0';
	connum[num2] = '\0';
	wildatom[id1].connum[id2] = 0;
	wildatom[id1].connum[id2] = atoi(connum);

	switch (name[0]) {
		case 'A':
			if (name[1] == 'c') 
				wildatom[id1].atomicnum[id2] = 89;
			else if (name[1] == 'g') 
				wildatom[id1].atomicnum[id2] = 47;
			else if (name[1] == 'l')
				wildatom[id1].atomicnum[id2] = 13;
			else if (name[1] == 'm')
				wildatom[id1].atomicnum[id2] = 95;
			else if (name[1] == 'r') 
				wildatom[id1].atomicnum[id2] = 18;
			else if (name[1] == 's') 
				wildatom[id1].atomicnum[id2] = 33;
			else if (name[1] == 't') 
				wildatom[id1].atomicnum[id2] = 85;
			else if (name[1] == 'u') 
				wildatom[id1].atomicnum[id2] = 79;
			break;
		case 'B':
			if (name[1] == 'r' || name[1] == 'R') 
				wildatom[id1].atomicnum[id2] = 35;
			else if (name[1] == 'a') 
				wildatom[id1].atomicnum[id2] = 56;
			else if (name[1] == 'e') 
				wildatom[id1].atomicnum[id2] = 4;
			else if (name[1] == 'h') 
				wildatom[id1].atomicnum[id2] = 107;
			else if (name[1] == 'i') 
				wildatom[id1].atomicnum[id2] = 83;
			else if (name[1] == 'k') 
				wildatom[id1].atomicnum[id2] = 97;
			else 
				wildatom[id1].atomicnum[id2] = 5;
			break;
		case 'C':
			if (name[1] == 'l' || name[1] == 'L') 
				wildatom[id1].atomicnum[id2] = 17;
			else if (name[1] == 'a') 
				wildatom[id1].atomicnum[id2] = 20;
			else if (name[1] == 'd')  
				wildatom[id1].atomicnum[id2] = 48;
			else if (name[1] == 'e')  
				wildatom[id1].atomicnum[id2] = 58;
			else if (name[1] == 'f')  
				wildatom[id1].atomicnum[id2] = 98;
			else if (name[1] == 'm')  
				wildatom[id1].atomicnum[id2] = 96;
			else if (name[1] == 'o') 
				wildatom[id1].atomicnum[id2] = 27;
			else if (name[1] == 'r')  
				wildatom[id1].atomicnum[id2] = 24;
			else if (name[1] == 's')  
				wildatom[id1].atomicnum[id2] = 55;
			else if (name[1] == 'u')  
				wildatom[id1].atomicnum[id2] = 29;
			else 
				wildatom[id1].atomicnum[id2] = 6;
			break;
		case 'D':
			if (name[1] == 'b')  
				wildatom[id1].atomicnum[id2] = 105;
			else if (name[1] == 's')  
				wildatom[id1].atomicnum[id2] = 110;
			else if (name[1] == 'y')  
				wildatom[id1].atomicnum[id2] = 66;
			else 
				wildatom[id1].atomicnum[id2] = 1;
			break;
		case 'E':
			if (name[1] == 'P') 
				wildatom[id1].atomicnum[id2] = 0;
			else if (name[1] == 'r') 
				wildatom[id1].atomicnum[id2] = 68;
			else if (name[1] == 's') 
				wildatom[id1].atomicnum[id2] = 99;
			else if (name[1] == 'u') 
				wildatom[id1].atomicnum[id2] = 63;
			break;
		case 'F':
			if (name[1] == 'e') 
				wildatom[id1].atomicnum[id2] = 26;
			else if (name[1] == 'm') 
				wildatom[id1].atomicnum[id2] = 100;
			else if (name[1] == 'r') 
				wildatom[id1].atomicnum[id2] = 87;
			else 
				wildatom[id1].atomicnum[id2] = 9;
			break;
		case 'G':
			if (name[1] == 'a') 
				wildatom[id1].atomicnum[id2] = 31;
			else if (name[1] == 'd') 
				wildatom[id1].atomicnum[id2] = 64;
			else if (name[1] == 'e') 
				wildatom[id1].atomicnum[id2] = 32;
			break;
		case 'H':
			if (name[1] == 'e') 
				wildatom[id1].atomicnum[id2] = 2;
			else if (name[1] == 'f') 
				wildatom[id1].atomicnum[id2] = 72;
			else if (name[1] == 'g') 
				wildatom[id1].atomicnum[id2] = 80;
			else if (name[1] == 'o') 
				wildatom[id1].atomicnum[id2] = 67;
			else if (name[1] == 's') 
				wildatom[id1].atomicnum[id2] = 108;
			else 
				wildatom[id1].atomicnum[id2] = 1;
			break;
		case 'I':
			if (name[1] == 'n') 
				wildatom[id1].atomicnum[id2] = 49;
			if (name[1] == 'r') 
				wildatom[id1].atomicnum[id2] = 77;
			else 
				wildatom[id1].atomicnum[id2] = 53;
			break;
		case 'K':
			if (name[1] == 'r') 
				wildatom[id1].atomicnum[id2] = 36;
			else 
				wildatom[id1].atomicnum[id2] = 19;
			break;
		case 'l':
			if (name[1] == 'p')
				wildatom[id1].atomicnum[id2] = 0;
			break;
		case 'L':
			if (name[1] == 'i') 
				wildatom[id1].atomicnum[id2] = 3;
			else if (name[1] == 'a') 
				wildatom[id1].atomicnum[id2] = 57;
			else if (name[1] == 'r') 
				wildatom[id1].atomicnum[id2] = 103;
			else if (name[1] == 'u') 
				wildatom[id1].atomicnum[id2] = 71;
			else if (name[1] == 'P') 
				wildatom[id1].atomicnum[id2] = 0;
			break;
		case 'M':
			if (name[1] == 'n') 
				wildatom[id1].atomicnum[id2] = 25;
			else if (name[1] == 'g') 
				wildatom[id1].atomicnum[id2] = 12;
			else if (name[1] == 'd') 
				wildatom[id1].atomicnum[id2] = 101;
			else if (name[1] == 'o') 
				wildatom[id1].atomicnum[id2] = 42;
			else if (name[1] == 't') 
				wildatom[id1].atomicnum[id2] = 109;
			break;
		case 'N':
			if (name[1] == 'i') 
				wildatom[id1].atomicnum[id2] = 28;
			else if (name[1] == 'a') 
				wildatom[id1].atomicnum[id2] = 11;
			else if (name[1] == 'b') 
				wildatom[id1].atomicnum[id2] = 41;
			else if (name[1] == 'd') 
				wildatom[id1].atomicnum[id2] = 60;
			else if (name[1] == 'e') 
				wildatom[id1].atomicnum[id2] = 10;
			else if (name[1] == 'o') 
				wildatom[id1].atomicnum[id2] = 102;
			else if (name[1] == 'p') 
				wildatom[id1].atomicnum[id2] = 93;
			else  
				wildatom[id1].atomicnum[id2] = 7;
			break;
		case 'O':
			if (name[1] == 's') 
				wildatom[id1].atomicnum[id2] = 76;
			else 
				wildatom[id1].atomicnum[id2] = 8;
			break;
		case 'P':
			if (name[1] == 'd') 
				wildatom[id1].atomicnum[id2] = 46;
			else if (name[1] == 't') 
				wildatom[id1].atomicnum[id2] = 78;
			else if (name[1] == 'b') 
				wildatom[id1].atomicnum[id2] = 82;
			else if (name[1] == 'a') 
				wildatom[id1].atomicnum[id2] = 91;
			else if (name[1] == 'm') 
				wildatom[id1].atomicnum[id2] = 61;
			else if (name[1] == 'o') 
				wildatom[id1].atomicnum[id2] = 84;
			else if (name[1] == 'r') 
				wildatom[id1].atomicnum[id2] = 59;
			else if (name[1] == 'u') 
				wildatom[id1].atomicnum[id2] = 94;
			else 
				wildatom[id1].atomicnum[id2] = 15;
			break;
		case 'R':
			if (name[1] == 'u') 
				wildatom[id1].atomicnum[id2] = 44;
			else if (name[1] == 'h') 
				wildatom[id1].atomicnum[id2] = 45;
			else if (name[1] == 'a') 
				wildatom[id1].atomicnum[id2] = 88;
			else if (name[1] == 'b') 
				wildatom[id1].atomicnum[id2] = 37;
			else if (name[1] == 'e') 
				wildatom[id1].atomicnum[id2] = 75;
			else if (name[1] == 'f') 
				wildatom[id1].atomicnum[id2] = 104;
			else if (name[1] == 'n') 
				wildatom[id1].atomicnum[id2] = 86;
			break;
		case 'S':
			if (name[1] == 'i' || name[1] == 'I') 
				wildatom[id1].atomicnum[id2] = 14;
			else if (name[1] == 'c') 
				wildatom[id1].atomicnum[id2] = 21;
			else if (name[1] == 'e') 
				wildatom[id1].atomicnum[id2] = 34;
			else if (name[1] == 'r') 
				wildatom[id1].atomicnum[id2] = 38;
			else if (name[1] == 'b') 
				wildatom[id1].atomicnum[id2] = 51;
			else if (name[1] == 'g') 
				wildatom[id1].atomicnum[id2] = 106;
			else if (name[1] == 'm') 
				wildatom[id1].atomicnum[id2] = 62;
			else if (name[1] == 'n') 
				wildatom[id1].atomicnum[id2] = 50;
			else  
				wildatom[id1].atomicnum[id2] = 16;
			break;
		case 'T':
			if (name[1] == 'i') 
				wildatom[id1].atomicnum[id2] = 22;
			else if (name[1] == 'l') 
				wildatom[id1].atomicnum[id2] = 81;
			else if (name[1] == 'a') 
				wildatom[id1].atomicnum[id2] = 73;
			else if (name[1] == 'b') 
				wildatom[id1].atomicnum[id2] = 65;
			else if (name[1] == 'c') 
				wildatom[id1].atomicnum[id2] = 43;
			else if (name[1] == 'e') 
				wildatom[id1].atomicnum[id2] = 52;
			else if (name[1] == 'h') 
				wildatom[id1].atomicnum[id2] = 90;
			else if (name[1] == 'm') 
				wildatom[id1].atomicnum[id2] = 69;
			else 
				wildatom[id1].atomicnum[id2] = 1;
			break;
		case 'U':
			wildatom[id1].atomicnum[id2] = 92;
			break;
		case 'V':
			wildatom[id1].atomicnum[id2] = 23;
			break;
		case 'W':
			wildatom[id1].atomicnum[id2] = 74;
			break;
		case 'X':
			if (name[1] == 'e') 
				wildatom[id1].atomicnum[id2] = 54;
			break;
		case 'Y':
			if (name[1] == 'b') 
				wildatom[id1].atomicnum[id2] = 70;
			else 
				wildatom[id1].atomicnum[id2] = 39;
			break;
		case 'Z':
			if (name[1] == 'n') 
				wildatom[id1].atomicnum[id2] = 30;
			else if (name[1] == 'r') 
				wildatom[id1].atomicnum[id2] = 40;
			break;
		default:
                        printf("\n Unrecognized atomic name %5s, exit", name);
                        wildatom[id1].atomicnum[id2] = 0;
                }
		
}


void atname(void)
{
	int i, j, k;
	int countatom[150];
	int typenum;
	char tmpchar[MAXCHAR];

	if (anindex == 1) {
		for (i = 0; i < 150; i++)
			countatom[i] = 0;
		for (j = 0; j < atomnum; j++) {
			countatom[atom[j].atomicnum]++;
			/* newitoa(countatom[atom[j].atomicnum], tmpchar); */
			sprintf(tmpchar, "%d", countatom[atom[j].atomicnum]);
			strcpy(atom[j].name, atom[j].element);
			strcat(atom[j].name, tmpchar);
		}
	}
	if (anindex == 2) {
		for (j = 0; j < atomnum; j++) {
			typenum = 1;
			for (k = 0; k < j; k++)
				if (strcmp(atom[k].ambername, atom[j].ambername) == 0)
					typenum++;
			sprintf(tmpchar, "%d", typenum); 
			/* newitoa(typenum, tmpchar); */
			strcpy(atom[j].name, atom[j].ambername);
			strcat(atom[j].name, tmpchar);
		}
	}
}

int jbond(int id1, int id2, int type)
{

	int i;
	int tmpint = 0;
	for (i = 0; i < bondnum; i++)
		if ((bond[i].bondi == id1 && bond[i].bondj == id2)
			|| (bond[i].bondi == id2 && bond[i].bondj == id1)) {
			if (type == 1
				&& (bond[i].type == 1 || bond[i].type == 7
					|| bond[i].type == 9||bond[i].type == 10))
				tmpint = 1;
			if (type == -1 && (bond[i].type == 1 || bond[i].type == 9))
				tmpint = 1;
			if (type == 2
				&& (bond[i].type == 2 || bond[i].type == 8
					|| bond[i].type == 9 ||bond[i].type == 10))
				tmpint = 1;
			if (type == -2 && (bond[i].type == 2 || bond[i].type == 9))
				tmpint = 1;
			if (type == 3 && bond[i].type == 3)
				tmpint = 1;
			if (type == -3 && bond[i].type == 3)
				tmpint = 1;
			if (type == 4
				&& (bond[i].type == 7 || bond[i].type == 8
					|| bond[i].type == 10))
				tmpint = 1;
			if (type == 9 && bond[i].type == 9)
				tmpint = 1;
			break;
		}
	return tmpint;
}

void improper(void)
{
	int i, j;
	int index;
	int tmpint;
	for (i = 0; i < atomnum; i++) {
		if (atom[i].atomicnum == 6 && atom[i].connum == 3)
			atom[i].improper = 1;
		if (atom[i].atomicnum == 7 && atom[i].connum == 3 &&
			(arty[i].ar1 + arty[i].ar2 + arty[i].ar3) > 0)
			atom[i].improper = 1;
		if (atom[i].ambername[0] == 'c' && atom[i].ambername[1] == ' ')
			/* sp2 carbonyl carbon */
		{
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
			if (index == 1)
				atom[i].improper = 1;
		}
	}
}
int cesbondcheck(int id1, int id2, char *str)
{
	int i;
	int index = 0;
	for (i = 0; i < bondnum; i++)
		if ((bond[i].bondi == id1 && bond[i].bondj == id2) ||
			(bond[i].bondi == id2 && bond[i].bondj == id1)) {
			if (strcmp(str, "any") == 0 || strcmp(str, "ANY") == 0
				|| strcmp(str, "Any") == 0) {
				index = 1;
				break;
			}
			if (bond[i].type == 1 && strcmp(str, "sb") == 0) {
				index = 1;
				break;
			}
			if (bond[i].type == 2 && strcmp(str, "db") == 0) {
				index = 1;
				break;
			}
			if (bond[i].type == 3 && strcmp(str, "tb") == 0) {
				index = 1;
				break;
			}
           		if (bond[i].type == 7 &&
                           (strcmp(str, "sb") == 0 || strcmp(str, "SB") == 0 ||
                                    strcmp(str,"AB") == 0 || strcmp(str,"ab") == 0)) {
                                index = 1;
                                break;
                        }
                        if (bond[i].type == 8 &&
                           (strcmp(str, "db") == 0 || strcmp(str, "DB") == 0 ||
                                    strcmp(str,"AB") == 0 || strcmp(str,"ab") == 0)) {
                                index = 1;
                                break;
                        }

			if (bond[i].type == 9
				&& (strcmp(str, "sb") == 0 || strcmp(str, "SB") == 0)) {
				index = 1;
				break;
			}
			if (bond[i].type == 9
				&& (strcmp(str, "db") == 0 || strcmp(str, "DB") == 0)) {
				index = 1;
				break;
			}
           		if (bond[i].type == 10
                                && (strcmp(str, "ab") == 0 || strcmp(str, "AB") == 0)) {
                                index = 1;
                                break;
                        }
		}
	return index;
}


int apcheck(int atmid, int preatmid, char *str)
{
/*there are apconsnum constraint (separated by ,) and each constraint has several APUNIT (separated by .)*/
	typedef struct {
		char propstr[20];
	}APUNIT; 
	typedef struct {
        	int num;
		APUNIT unit[MAXAPUNIT];
	}AP; 
	AP ap[MAXAPCONS];
	int apconsnum = 0;

	int i, j, k;
	int suc;
	int propcount;
	int stopflag ;
	int numunit;
	int numchar ;
	
	char propname[20];
	char tmpchar[20];
	char tmpchar1[20];
	char tmpchar2[20];
	int tmpint1, tmpint2;
	int tmprg, tmprg3, tmprg4, tmprg5, tmprg6, tmprg7, tmprg8, tmprg9, tmprg10;
	int tmpar1, tmpar2, tmpar3, tmpar4, tmpar5;
	int tmpsb, tmpSB, tmpdb, tmpDB, tmpDL, tmptb, tmpTB, tmpAB, tmpnr;
	int sbindex = 0;
	int SBindex = 0;
	int dbindex = 0;
	int DBindex = 0;
	int DLindex = 0;
	int tbindex = 0;
	int TBindex = 0;
	int ABindex = 0;

	tmprg = arty[atmid].rg[0];
	tmprg3 = arty[atmid].rg[3];
	tmprg4 = arty[atmid].rg[4];
	tmprg5 = arty[atmid].rg[5];
	tmprg6 = arty[atmid].rg[6];
	tmprg7 = arty[atmid].rg[7];
	tmprg8 = arty[atmid].rg[8];
	tmprg9 = arty[atmid].rg[9];
	tmprg10 = arty[atmid].rg[10];
	tmpar1 = arty[atmid].ar1;
	tmpar2 = arty[atmid].ar2;
	tmpar3 = arty[atmid].ar3;
	tmpar4 = arty[atmid].ar4;
	tmpar5 = arty[atmid].ar5;
	tmpnr = arty[atmid].nr;
	tmpsb = sb[atmid];
	tmpSB = SB[atmid];
	tmpdb = db[atmid];
	tmpDB = DB[atmid];
	tmpDL = DL[atmid];
	tmptb = tb[atmid];
	tmpTB = TB[atmid];
	tmpAB = AB[atmid];

/* parse */
	apconsnum = 0;
	numchar = 0;
	numunit = 0;
	for (i = 0; i < strlen(str); i++) {
		if (str[i] == '[')
			continue;
		if (str[i] == ']') {
			ap[apconsnum].unit[numunit].propstr[numchar] = '\0';
			numunit++;
			if(numunit >= MAXAPUNIT) {
				fprintf(stdout, "Too many units, extend MAXAPUNIT and recompile, exit\n");
				exit(1);
			}
			ap[apconsnum].num = numunit;	
			apconsnum++;	
			if(apconsnum >= MAXAPCONS) {
				fprintf(stdout, "Too many ap constraints, extend MAXAPCONS and recompile, exit\n");
				exit(1);
			}
			break;
		}
		if (str[i] == '.') {
			ap[apconsnum].unit[numunit].propstr[numchar] = '\0';
			numunit ++;
			if(numunit >= MAXAPUNIT) {
				fprintf(stdout, "Too many units, extend MAXAPUNIT and recompile, exit\n");
				exit(1);
			}
			numchar = 0;
			continue;
		}
		if (str[i] == ',') {
			ap[apconsnum].unit[numunit].propstr[numchar] = '\0';
			numunit++;
			if(numunit >= MAXAPUNIT) {
				fprintf(stdout, "Too many units, extend MAXAPUNIT and recompile, exit\n");
				exit(1);
			}
			ap[apconsnum].num = numunit;	
			apconsnum++;	
			if(apconsnum >= MAXAPCONS) {
				fprintf(stdout, "Too many ap constraints, extend MAXAPCONS and recompile, exit\n");
				exit(1);
			}
			/* initialization*/
			ap[apconsnum].num = 0;
			numunit = 0;
			numchar = 0;
			continue;
		}
		ap[apconsnum].unit[numunit].propstr[numchar] = str[i];
		numchar ++;
	}

	for (i = 0; i< apconsnum; i++) {
		suc = 0;
		for(j=0;j<ap[i].num;j++) {
/* first of all, separate name and count from propstr*/
			strcpy(tmpchar, "");
			tmpchar[0]='\0';
			strcpy(tmpchar, ap[i].unit[j].propstr);
			strcpy(tmpchar1, "");
			tmpchar1[0]='\0';
			strcpy(tmpchar2, "");
			tmpchar2[0]='\0';
			tmpint1 = 0;
			tmpint2 = 0;
			stopflag = 0;
			for(k=0;k<strlen(tmpchar);k++) {
				if(stopflag == 0 && tmpchar[k] >='0' && tmpchar[k] <= '9') {
					tmpchar1[tmpint1] = tmpchar[k];
					tmpint1++;
				}
				else {
					tmpchar2[tmpint2] = tmpchar[k];
					stopflag = 1;
					tmpint2++;
				}
			}
			tmpchar1[tmpint1] = '\0';	
			tmpchar2[tmpint2] = '\0';	
			if(tmpint1 > 0) propcount = atoi(tmpchar1); 
			else 
				propcount = -1;
			strcpy(propname, tmpchar2);
			/* the following code judge if atom has this special property */
                        if (strcmp(propname, "RG") == 0) 
				if (propcount == tmprg || (propcount == -1 && tmprg > 0)) {
					suc = 1;
					break;
				}
                        if (strcmp(propname, "RG3") == 0) 
				if (propcount == tmprg3 || (propcount == -1 && tmprg3 > 0)) {
					suc = 1;
					break;
				}
                        if (strcmp(propname, "RG4") == 0) 
				if (propcount == tmprg4 || (propcount == -1 && tmprg4 > 0)) {
					suc = 1;
					break;
				}
                        if (strcmp(propname, "RG5") == 0) 
				if (propcount == tmprg5 || (propcount == -1 && tmprg5 > 0)) {
					suc = 1;
					break;
				}
                        if (strcmp(propname, "RG6") == 0) 
				if (propcount == tmprg6 || (propcount == -1 && tmprg6 > 0)) {
					suc = 1;
					break;
				}
                        if (strcmp(propname, "RG7") == 0) 
				if (propcount == tmprg7 || (propcount == -1 && tmprg7 > 0)) {
					suc = 1;
					break;
				}
                        if (strcmp(propname, "RG8") == 0) 
				if (propcount == tmprg8 || (propcount == -1 && tmprg8 > 0)) {
					suc = 1;
					break;
				}
                        if (strcmp(propname, "RG9") == 0) 
				if (propcount == tmprg9 || (propcount == -1 && tmprg9 > 0)) {
					suc = 1;
					break;
				}
                        if (strcmp(propname, "RG10") == 0) 
				if (propcount == tmprg10 || (propcount == -1 && tmprg10 > 0)) {
					suc = 1;
					break;
				}
                        if (strcmp(propname, "NR") == 0) 
				if (propcount == tmpnr || (propcount == -1 && tmpnr > 0)) {
					suc = 1;
					break;
				}
                        if (strcmp(propname, "AR1") == 0) 
				if (propcount == tmpar1 || (propcount == -1 && tmpar1 > 0)) {
					suc = 1;
					break;
				}
                        if (strcmp(propname, "AR2") == 0) 
				if (propcount == tmpar2 || (propcount == -1 && tmpar2 > 0)) {
					suc = 1;
					break;
				}
                        if (strcmp(propname, "AR3") == 0) 
				if (propcount == tmpar3 || (propcount == -1 && tmpar3 > 0)) {
					suc = 1;
					break;
				}
                        if (strcmp(propname, "AR4") == 0) 
				if (propcount == tmpar4 || (propcount == -1 && tmpar4 > 0)) {
					suc = 1;
					break;
				}
                        if (strcmp(propname, "AR5") == 0) 
				if (propcount == tmpar5 || (propcount == -1 && tmpar5 > 0)) {
					suc = 1;
					break;
				}

                       if (propname[0] == 'S' && propname[1] == 'B') {
				SBindex = 0;
                                if (propname[2] == '\'')
                                        SBindex = 1;
                                if (SBindex == 1 && propname[3] == '\'')
                                        SBindex = 2;
				if (propcount == tmpSB || (propcount == -1 && tmpSB > 0)) {
					if(SBindex == 0) suc = 1;
                                	if (SBindex == 1 && preatmid != -1 &&
                                            jbond(atmid, preatmid, -1) == 1)
						suc = 1;
                                	if (SBindex == 2 && preatmid != -1 &&
                                            jbond(atmid, preatmid, -1) == 0)
						suc = 1;
					if(suc == 1) break;
				}
                        }

                       if (propname[0] == 's' && propname[1] == 'b') {
				sbindex = 0;
                                if (propname[2] == '\'')
                                        sbindex = 1;
                                if (sbindex == 1 && propname[3] == '\'')
                                        sbindex = 2;
				if (propcount == tmpsb || (propcount == -1 && tmpsb > 0)) {
					if(sbindex == 0) suc = 1;
                                	if (sbindex == 1 && preatmid != -1 &&
                                            jbond(atmid, preatmid, 1) == 1)
						suc = 1;
                                	if (sbindex == 2 && preatmid != -1 &&
                                            jbond(atmid, preatmid, 1) == 0)
						suc = 1;
					if(suc == 1) break;
				}
                        }

                       if (propname[0] == 'D' && propname[1] == 'B') {
                                DBindex = 0;
                                if (propname[2] == '\'')
                                        DBindex = 1;
                                if (DBindex == 1 && propname[3] == '\'')
                                        DBindex = 2;
                                if (propcount == tmpDB || (propcount == -1 && tmpDB > 0)) {
                                        if(DBindex == 0) suc = 1;
                                        if (DBindex == 1 && preatmid != -1 &&
                                            jbond(atmid, preatmid, -2) == 1)
                                                suc = 1;
                                        if (DBindex == 2 && preatmid != -1 &&
                                            jbond(atmid, preatmid, -2) == 0)
                                                suc = 1;
                                        if(suc == 1) break;
                                }
                        }

                       if (propname[0] == 'd' && propname[1] == 'b') {
                                dbindex = 0;
                                if (propname[2] == '\'')
                                        dbindex = 1;
                                if (dbindex == 1 && propname[3] == '\'')
                                        dbindex = 2;
                                if (propcount == tmpdb || (propcount == -1 && tmpdb > 0)) {
                                        if(dbindex == 0) suc = 1;
                                        if (dbindex == 1 && preatmid != -1 &&
                                            jbond(atmid, preatmid, 2) == 1)
                                                suc = 1;
                                        if (dbindex == 2 && preatmid != -1 &&
                                            jbond(atmid, preatmid, 2) == 0)
                                                suc = 1;
                                        if(suc == 1) break;
                                }
                        }

                       if (propname[0] == 'T' && propname[1] == 'B') {
                                TBindex = 0;
                                if (propname[2] == '\'')
                                        TBindex = 1;
                                if (TBindex == 1 && propname[3] == '\'')
                                        TBindex = 2;
                                if (propcount == tmpTB || (propcount == -1 && tmpTB > 0)) {
                                        if(TBindex == 0) suc = 1;
                                        if (TBindex == 1 && preatmid != -1 &&
                                            jbond(atmid, preatmid, -3) == 1)
                                                suc = 1;
                                        if (TBindex == 2 && preatmid != -1 &&
                                            jbond(atmid, preatmid, -3) == 0)
                                                suc = 1;
                                        if(suc == 1) break;
                                }
                        }

                       if (propname[0] == 't' && propname[1] == 'b') {
                                tbindex = 0;
                                if (propname[2] == '\'')
                                        tbindex = 1;
                                if (tbindex == 1 && propname[3] == '\'')
                                        tbindex = 2;
                                if (propcount == tmptb || (propcount == -1 && tmptb > 0)) {
                                        if(tbindex == 0) suc = 1;
                                        if (tbindex == 1 && preatmid != -1 &&
                                            jbond(atmid, preatmid, 3) == 1)
                                                suc = 1;
                                        if (tbindex == 2 && preatmid != -1 &&
                                            jbond(atmid, preatmid, 3) == 0)
                                                suc = 1;
                                        if(suc == 1) break;
                                }
                        }

                       if (propname[0] == 'D' && propname[1] == 'L') {
                                DLindex = 0;
                                if (propname[2] == '\'')
                                        DLindex = 1;
                                if (DLindex == 1 && propname[3] == '\'')
                                        DLindex = 2;
                                if (propcount == tmpDL || (propcount == -1 && tmpDL > 0)) {
                                        if(DLindex == 0) suc = 1;
                                        if (DLindex == 1 && preatmid != -1 &&
                                            jbond(atmid, preatmid, 9) == 1)
                                                suc = 1;
                                        if (DLindex == 2 && preatmid != -1 &&
                                            jbond(atmid, preatmid, 9) == 0)
                                                suc = 1;
                                        if(suc == 1) break;
                                }
                        }
                       if (propname[0] == 'A' && propname[1] == 'B') {
                                ABindex = 0;
                                if (propname[2] == '\'')
                                        ABindex = 1;
                                if (ABindex == 1 && propname[3] == '\'')
                                        ABindex = 2;
                                if (propcount == tmpAB || (propcount == -1 && tmpAB > 0)) {
                                        if(ABindex == 0) suc = 1;
                                        if (ABindex == 1 && preatmid != -1 &&
                                            jbond(atmid, preatmid, 4) == 1)
                                                suc = 1;
                                        if (ABindex == 2 && preatmid != -1 &&
                                            jbond(atmid, preatmid, 4) == 0)
                                                suc = 1;
                                        if(suc == 1) break;
                                }
                        }
			if(suc == 1) break;
		}
		if(suc == 0) return 0; 
	}
	return suc;
}


void bondinfo(void)
{
	int i;
	for (i = 0; i < atomnum; i++) {
		sb[i] = 0;
		SB[i] = 0;
		db[i] = 0;
		DB[i] = 0;
		tb[i] = 0;
		TB[i] = 0;
		AB[i] = 0;
		DL[i] = 0;
	}
	for (i = 0; i < bondnum; i++) {
		switch (bond[i].type) {
		case 1:
			sb[bond[i].bondi]++;
			sb[bond[i].bondj]++;
			SB[bond[i].bondi]++;
			SB[bond[i].bondj]++;
			continue;
		case 2:
			db[bond[i].bondi]++;
			db[bond[i].bondj]++;
			DB[bond[i].bondi]++;
			DB[bond[i].bondj]++;
			continue;
		case 3:
			tb[bond[i].bondi]++;
			tb[bond[i].bondj]++;
			TB[bond[i].bondi]++;
			TB[bond[i].bondj]++;
			continue;
		case 7:
			AB[bond[i].bondi]++;
			AB[bond[i].bondj]++;
			sb[bond[i].bondi]++;
			sb[bond[i].bondj]++;
			continue;
		case 8:
			AB[bond[i].bondi]++;
			AB[bond[i].bondj]++;
			db[bond[i].bondi]++;
			db[bond[i].bondj]++;
			continue;
		case 9:
			sb[bond[i].bondi]++;
			sb[bond[i].bondj]++;
			SB[bond[i].bondi]++;
			SB[bond[i].bondj]++;
			DL[bond[i].bondi]++;
			DL[bond[i].bondj]++;
			continue;
		case 10:
			AB[bond[i].bondi]++;
			AB[bond[i].bondj]++;
			continue;
		default:
/*    fprintf(stdout, "unrecognized bond type %d\n", bond[i].type); */
			continue;
		}
	}
}

int chain_check()
{
	int i, j, k;
	int min;
	int flag;
	int tmpint1, tmpint2;
	int array1[MAXBEED];
	int array2[MAXBEED];

	for (i = 0; i < chain; i++)
		for (j = i + 1; j < chain; j++) {
			tmpint1 = select_chain_index[i];
			tmpint2 = select_chain_index[j];
			if (tmpint1 == tmpint2)
				return 0;
			for (k = 0; k < schain[tmpint1].num; k++)
				array1[k] = schain[tmpint1].at[k];
			for (k = 0; k < schain[tmpint1].num; k++)
				array2[k] = schain[tmpint2].at[k];
			if (schain[tmpint1].num < schain[tmpint2].num)
				min = schain[tmpint1].num;
			else
				min = schain[tmpint2].num;
			flag = 0;
			for (k = 0; k < min; k++)
				if (array1[k] != array2[k]) {
					flag = k + 1;
					break;
				}
			if (flag == 0)
				return 0;
/* New Code Since AMVBER10*/
			if (flag != 0)
                                for (k = 0; k < min; k++) {
                                        if (array1[k] == array2[k])
                                                if (strcmp(ces[i][k].cesname, ces[j][k].cesname) != 0) return 0;
                                        if (array1[k] != array2[k])
                                                if (strcmp(ces[i][k].cesname, ces[j][k].cesname) == 0) return 0;
                                 }
		}
	if (ces_bond_num != 0) {
		for (i = 0; i < atomnum; i++)
			strcpy(atomcesname[i].name, "");
		for (i = 0; i < chain; i++)
			for (j = 0; j < scnum; j++)
				if (select_chain_index[i] == j) {
					if (debug == 1)
						printf("\n%5s %5d %5d %5d", atom[caid].name, i, j,
							   ceslen[i]);
					for (k = 0; k < ceslen[i]; k++)
						strcpy(atomcesname[schain[j].at[k]].name,
							   ces[i][k].cesname);
				}
		strcpy(atomcesname[caid].name, "sa");
		for (i = 0; i < ces_bond_num; i++) {
			tmpint1 = -1;
			tmpint2 = -1;
			for (j = 0; j < atomnum; j++) {
				if (strcmp(atomcesname[j].name, ces_bond_at1[i]) == 0)
					tmpint1 = j;
				if (strcmp(atomcesname[j].name, ces_bond_at2[i]) == 0)
					tmpint2 = j;
			}
			if (tmpint1 == tmpint2 && tmpint1 != -1)
				return 0;
			if (cesbondcheck(tmpint1, tmpint2, ces_bond_type[i]) == 0)
				return 0;
		}
	}
	return 1;
}

void dccheck(int startnum)
{
	int i;
	if (suc_flag == 1)
		return;
	startnum++;
	for (i = 0; i < scnum; i++) {
		if (schain[i].id != startnum - 1)
			continue;
		if (schain[i].id == startnum - 1)
			select_chain_index[startnum - 1] = i;

		if (startnum == chain) {
			suc_flag = chain_check();
			if (debug == 1) {
				printf("\n--- begin ---");
				for (j = 0; j < chain; j++)
					printf("\n%5s %5d", atom[caid].name,
						   select_chain_index[j]);
				printf("\n--- end ---");
			}
			if (suc_flag == 1)
				return;
		}
		if (startnum < chain)
			dccheck(startnum);
	}
	return;
}

int wematch(char *name, int atmid)
{
	int i, j;
	int index = 0;
	int breakindex = 0;
	if (wenum == 0)
		return index;
	for (i = 0; i < wenum; i++)
		if (name[0] == wildatom[i].wildname[0])
			if (name[1] == wildatom[i].wildname[1]) {
				for (j = 0; j < wildatom[i].num; j++)
					if ((atom[atmid].atomicnum == wildatom[i].atomicnum[j])
						&& (atom[atmid].connum == wildatom[i].connum[j]
							|| wildatom[i].connum[j] == 0)) {
						index = 1;
						breakindex = 1;
						break;
					}
				if (breakindex == 1)
					break;
			}
	return index;
}

void cematch(ATOM atm[], int selectnum, int startnum)
{
	int i, j, k, m;
	int start;
	int index;
	int index1;
	int breakindex = 0;
	int tmpint;
	start = -1;
	selectchain[selectnum] = startnum;
	selectnum++;
	selectindex[startnum] = 1;
	for (k = 0; k < chain; k++) {
		index = 1;
/* 		if (cesindex[k] != 0) continue; */
		if (selectnum - 1 == ceslen[k]) {
/*check atom name, connum, property string */
			for (j = 1; j < selectnum; j++) {
				index1 = 1;
				tmpint = 0;
				if (atm[selectchain[j]].connum != ces[k][j - 1].atconnum
					&& ces[k][j - 1].atconnum != 0) {
					index = 0;
					break;
				}
				if (ces[k][j - 1].atname[0] == 'E'
					&& ces[k][j - 1].atname[1] == 'W')
					if (atm[selectchain[j]].ewd != 1) {
						index = 0;
						break;
					}
				if (atm[selectchain[j]].name[0] != ces[k][j - 1].atname[0])
					index1 = 0;
				if(index1 == 1) {
					if (strlen(ces[k][j - 1].atname) >= 2)
						if (atm[selectchain[j]].name[1] != ces[k][j - 1].atname[1])
							index1 = 0;
                                	if ((atm[selectchain[j]].element[1] >= 'a' && atm[selectchain[j]].element[1] <= 'z') ||
                                    	    (atm[selectchain[j]].element[1] >= 'A' && atm[selectchain[j]].element[1] <= 'Z'))
                                        	if(strlen(ces[k][j - 1].atname) == 1)
                                        		index1 = 0;
				}

				if (index1 == 0) {
					tmpint = wematch(ces[k][j - 1].atname, selectchain[j]);
					if (tmpint == 0) {
						index = 0;
						break;
					}
				}

				if (ces[k][j - 1].apindex == 1) {
					if (j == 1)
						if (apcheck(selectchain[j], caid, ces[k][j - 1].ap)
							== 0) {
							index = 0;
							break;
						}
					if (j != 1)
						if (apcheck
							(selectchain[j], selectchain[j - 1],
							 ces[k][j - 1].ap) == 0) {
							index = 0;
							break;
						}
				}
			}
			if (index == 1) {
				cesindex[k]++;
				schain[scnum].num = selectnum - 1;
				schain[scnum].id = k;
				for (m = 0; m < selectnum - 1; m++)
					schain[scnum].at[m] = selectchain[m + 1];
				scnum++;
				if (scnum >= MAXSCHAIN) {
					printf
						("\nnumber of schain exceeds MAXSCHAIN, please increase MAXSCHAIN, exit");
					exit(1);
				}
			}
		}
	}
	for (i = 0; i < 6; i++) {
		if (atm[startnum].con[i] == -1)
			return;
		start = atm[startnum].con[i];
		/* we have already visited this atom */
		for (j = 0; j < selectnum; j++)
			if (selectchain[j] == start) {
				breakindex = 1;
				break;
			}
		if (breakindex == 1) {
			breakindex = 0;
			continue;
		}
		if (selectnum > maxchain)
			return;
		cematch(atm, selectnum, start);
	}
}


int jatspecial(int atomno, char keyword[1024], char ces_bond_string[1024])
{
	int i, j;
	int tmpint1, tmpint2;
	int layer;
	int index;
	int index0;
	int apindex[MAXBEED];
	int apnum;
	int tmpapindex = 0;
	char tmpap[MAXCHAR] = "";
	char ap[MAXBEED][MAXCHAR];
	char tmpchar[10];
	char atname[MAXBEED][10];
	int atconnum[MAXBEED];

	/*for defined atom name in chemical environment string */
	int cesname_index = 0;		/*for atom connectivity in chemical environment string */
	int cesname_num;
	char tmpcesname[MAXCHAR] = "";
	char cesname[MAXBEED][MAXCHAR];

	/*for ces_bond_string */
	int ccnum;
	int ccindex;
	char ccs[MAXCHAR];

	int cea_id = 1;

	/* check first */
	for (i = 0; i < MAXCES; i++)
		for (j = 0; j < MAXBEED; j++) {
			ces[i][j].atconnum = 0;
			ces[i][j].atname[0] = '\0';
		}
	for (j = 0; j < MAXBEED; j++) {
		atconnum[j] = 0;
		apindex[j] = 0;
		atname[j][0] = '\0';
		ap[j][0] = '\0';
		cesname[j][0] = '\0';
	}
	index0 = 0;
	tmpint1 = 0;
	tmpint2 = 0;
	for (i = 0; i < strlen(keyword); i++) {
		if (keyword[i] == '(')
			tmpint1++;
		if (keyword[i] == ')')
			tmpint2++;
	}
	if (tmpint1 != tmpint2) {
		fprintf(stdout, "( and ) does not match for %s, exit\n", keyword);
		exit(1);
	}
	tmpint1 = 0;
	tmpint2 = 0;
	for (i = 0; i < strlen(keyword); i++) {
		if (keyword[i] == '[')
			tmpint1++;
		if (keyword[i] == ']')
			tmpint2++;
	}
	if (tmpint1 != tmpint2) {
		fprintf(stdout, "[ and ] does not match for %s, exit\n", keyword);
		exit(1);
	}
	tmpint1 = 0;
	tmpint2 = 0;
	for (i = 0; i < strlen(keyword); i++) {
		if (keyword[i] == '<')
			tmpint1++;
		if (keyword[i] == '>')
			tmpint2++;
	}
	if (tmpint1 != tmpint2) {
		fprintf(stdout, "< and > does not match for %s, exit\n", keyword);
		exit(1);
	}
	layer = 0;
	chain = 0;
	for (i = 0; i < strlen(keyword); i++) {
                if (tmpapindex ==0 && cesname_index == 0)
			if(((keyword[i]   >= 'a' && keyword[i]   <= 'z') || (keyword[i]   >= 'A' && keyword[i]   <= 'Z')) &&
                           ((keyword[i+1] >= 'a' && keyword[i+1] <= 'z') || (keyword[i+1] >= 'A' && keyword[i+1] <= 'Z')))
                        	continue;

		if (keyword[i] == '(') {
			layer++;
			if (layer >= MAXBEED) {
				printf
					("\nNumber of beed in ces chains exceed MAXBEED, please increase MAXBEED, exit\n");
				exit(1);
			}
		}
		if (keyword[i] == ')')
			layer--;

		if (tmpapindex == 0 && keyword[i] == '[') {
			tmpapindex = 1;
			strcpy(tmpap, "");
			tmpap[0] = '[';
			apnum = 1;
			continue;
		}
		if (tmpapindex == 1 && keyword[i] == ']') {
			apindex[layer] = 1;
			tmpap[apnum++] = ']';
			tmpap[apnum] = '\0';
			strcpy(ap[layer], tmpap);
			tmpapindex = 0;
			continue;
		}
		if (tmpapindex == 1) {
			tmpap[apnum++] = keyword[i];
			continue;
		}

/*for defined atom name in chemical environment string */
		if (cesname_index == 0 && keyword[i] == '<') {
			strcpy(tmpcesname, "");
			cesname_index = 1;
			cesname_num = 0;
			continue;
		}
		if (cesname_index == 1 && keyword[i] == '>') {
			tmpcesname[cesname_num] = '\0';
			strcpy(cesname[layer], tmpcesname);
			cesname_index = 0;
			continue;
		}
		if (cesname_index == 1) {
			tmpcesname[cesname_num++] = keyword[i];
			continue;
		}


		if (keyword[i] == ',' && keyword[i - 1] != ')') {
			if (layer >= MAXBEED) {
				printf
					("\nNumber of beed in ces chains exceed MAXBEED, please increase MAXBEED, exit\n");
				exit(1);
			}
			for (j = 0; j <= layer; j++) {
				strcpy(ces[chain][j].atname, atname[j + 1]);
				ces[chain][j].atconnum = atconnum[j + 1];
				ces[chain][j].apindex = apindex[j + 1];
				strcpy(ces[chain][j].ap, ap[j + 1]);
				strcpy(ces[chain][j].cesname, cesname[j + 1]);
			}
			chain++;
			ceslen[chain - 1] = layer;
			if (chain >= MAXCES) {
				printf
					("\nNumber of ces chains exceed MAXCES, please increase MAXCES, exit\n");
				exit(1);
			}
		}
		if (keyword[i] == ')' && keyword[i - 1] != ')') {
			for (j = 0; j <= layer + 1; j++) {
				strcpy(ces[chain][j].atname, atname[j + 1]);
				ces[chain][j].atconnum = atconnum[j + 1];
				ces[chain][j].apindex = apindex[j + 1];
				strcpy(ces[chain][j].ap, ap[j + 1]);
				strcpy(ces[chain][j].cesname, cesname[j + 1]);
			}
			chain++;
			ceslen[chain - 1] = layer + 1;
			if (chain >= MAXCES) {
				printf
					("\nNumber of ces chains exceed MAXCES, please increase MAXCES, exit\n");
				exit(1);
			}
		}


		if (((keyword[i]     >= 'a' && keyword[i]     <= 'z') || (keyword[i]     >= 'A' && keyword[i]     <= 'Z')) && 
                    ((keyword[i + 1] >= 'a' && keyword[i + 1] <= 'z') || (keyword[i + 1] >= 'A' && keyword[i + 1] <= 'Z')))
			continue;
		if ((keyword[i] >= 'a' && keyword[i] <= 'z') || 
                    (keyword[i] >= 'A' && keyword[i] <= 'Z')) {
			index0 = 1;
			if((keyword[i - 1] >= 'a' && keyword[i - 1] <= 'z') || 
                           (keyword[i - 1] >= 'A' && keyword[i - 1] <= 'Z')) {
				atname[layer][0] = keyword[i - 1];
				atname[layer][1] = keyword[i];

				strcpy(cesname[layer], "");
				strcpy(ap[layer], "");
				apindex[layer] = 0;

			} else {
				atname[layer][0] = keyword[i];
				atname[layer][1] = '\0';

				strcpy(cesname[layer], "");
				strcpy(ap[layer], "");
				apindex[layer] = 0;
			}
                        strcpy(tmpchar, "");
                        tmpchar[0] = '\0';
                        sprintf(tmpchar, "%d", cea_id);
                        cea_id++;
                        strcpy(cesname[layer], "cea_");
                        strcat(cesname[layer],tmpchar);
		}
		/* if(isdigit(keyword[i])==16) */
		if (keyword[i] <= '9' && keyword[i] >= '0') {
			tmpchar[0] = keyword[i];
			atconnum[layer] = atoi(tmpchar);
		} else if (index0 == 1) {
			atconnum[layer] = 0;
			index0 = 0;
		}
	}
	/* parse ces_bond_string: like c1:c2:db,c1:c2:sb */
	ces_bond_num = 0;
	if (ces_bond_index == 1) {
		ccnum = 0;
		ccindex = 0;
		strcpy(ccs, "");
		for (i = 0; i <= strlen(ces_bond_string); i++) {
			if (ces_bond_string[i] == ',' || i == strlen(ces_bond_string)) {
				ccs[ccnum] = '\0';
				strcpy(ces_bond_type[ces_bond_num], ccs);
				strcpy(ccs, "");
				ccindex = 0;
				ccnum = 0;
				ces_bond_num++;
				if (ces_bond_num >= MAX_CES_BOND) {
					printf
						("\nNumber of connected atom pairs in chemical environment string exceed MAX_CES_CON, please increase MAX_CES_CON, exit\n");
					exit(1);
				}

				continue;
			}
			if (ces_bond_string[i] == ':') {
				ccs[ccnum] = '\0';
				if (ccindex == 0)
					strcpy(ces_bond_at1[ces_bond_num], ccs);
				if (ccindex == 1)
					strcpy(ces_bond_at2[ces_bond_num], ccs);
				ccnum = 0;
				strcpy(ccs, "");
				ccindex++;
				continue;
			}
			ccs[ccnum] = ces_bond_string[i];
			ccnum++;
		}
	}

	if (debug == 1) {
		printf("\nChemical environment stringS are listed as the following");
		printf
			("\nNo. length  elem  connum  atom_property_index atom_property ces_atom_name");
		for (i = 0; i < chain; i++) {
			for (j = 0; j < ceslen[i]; j++)
				printf("\n%5d%5d %s %5d  %d %s %s", i, j, ces[i][j].atname,
					   ces[i][j].atconnum, ces[i][j].apindex, ces[i][j].ap,
					   ces[i][j].cesname);
		}
		printf("\nBond Connectivity in Chemical environment strings: %5d%5d",
			   ces_bond_num, ces_bond_index);
		for (i = 0; i < ces_bond_num; i++)
			printf("\n%5d%5s %5s %5s", i + 1, ces_bond_at1[i], ces_bond_at2[i],
				   ces_bond_type[i]);
	}
	maxchain = -1;
	scnum = 0;
	for (i = 0; i < chain; i++) {
		cesindex[i] = 0;
		if (ceslen[i] > maxchain)
			maxchain = ceslen[i];
	}
	for (i = 0; i < atomnum; i++) {
		selectindex[i] = -1;
		selectchain[i] = -1;
		strcpy(atomcesname[i].name, "");
	}

	cematch(atom, 0, atomno);
	index = 1;
	for (i = 0; i < chain; i++)
		if (cesindex[i] == 0) {
			index = 0;
			break;
		}
	if (index == 1)
		if (debug == 1) {
			printf("\nATOM : %s", atom[caid].name);
			printf("\n\nThere are %5d  schains ", scnum);
			printf("\nno   num   ce_id  cesindex[ce_id)]");
			printf("\nschain[i].at[j], atom[schain[i].at[j]].name");
			for (i = 0; i < scnum; i++) {
				printf("\n====================================== ");
				printf("\n%5d %5d %5d %5d", i + 1, schain[i].num,
					   schain[i].id + 1, cesindex[schain[i].id]);
				printf("\n-------------------------------------- ");
				for (j = 0; j < schain[i].num; j++)
					printf("\n%5d %5s", schain[i].at[j],
						   atom[schain[i].at[j]].name);
				printf("\n====================================== ");

			}
		}
	if (index == 1) {
		suc_flag = 0;
		for (i = 0; i < chain; i++)
			select_chain_index[i] = -1;
		dccheck(0);
		index = suc_flag;
	}

	if (index == 1)
		return 1;
	else
		return 0;
}


/* The following code offer new method to judge the amber atomic type */
void jat(void)
{
	int i, j;
	char line[3 * MAXCHAR];
	char keyword0[10];			/* atomic type */
	char keyword1[10];			/* residue name */
	char keyword2[10];			/* atomic number */
	char keyword3[10];			/* connected atomic number */
	char keyword4[10];			/* the number of directly connected hydrogen */
	char keyword5[10];			/* the number of electron withdraw atom of the 
								   immediately connected atom, meaningful only for hydrogen */
	char keyword6[MAXCHAR];		/* atomic property */
	char keyword7[MAXCHAR];		/* special chemical environment */
	char keyword8[MAXCHAR];		/* bond connectivity in chemical environment string */
	char tmpchar[10];
	char wildname[10];
	int tmpint, tmpint1, tmpint2;
	int resindex;
	int readindex;
	int namenum;
	FILE *fpin;

	if ((fpin = fopen(atomtype_def, "r")) == NULL) {
		fprintf(stdout, "Cannot open the atom type defination file: %s, exit\n", atomtype_def);
		exit(1);
	}
	
/* start to judge amber atomic type */
	wenum = 0;
	for (i = 0; i < atomnum; i++) {	/* 1st cycle */
		for (;;) {				/*2nd cycle */
			if (fgets(line, MAXCHAR * 2, fpin) == NULL) {
				fprintf(stdout, "Finished reading atomtype file.\n");
				break;
			}
			if (strncmp("WILDATOM", line, 8) == 0 && weindex == 1) {
				sscanf(&line[9], "%s", wildatom[wenum].wildname);
				namenum = 0;
				readindex = 0;

				for (j = 9; j < strlen(line); j++) {
					if (readindex == 0 && line[j] != ' ') {
						readindex = 2;
						continue;
					}
					if (readindex == 2 && line[j] != ' ')
						continue;
					if (readindex == 2 && line[j] == ' ') {
						readindex = 3;
						continue;
					}
					if (readindex == 3) {
						sscanf(&line[j], "%s", wildname);
						assignwildatom(wenum, namenum, wildname);
						namenum++;
						readindex = 4;
						continue;
					}
					if (readindex == 4 && line[j] != ' ')
						continue;
					if (readindex == 4 && line[j] == ' ') {
						readindex = 3;
						continue;
					}

				}
				wildatom[wenum].num = namenum;
				wenum++;
				if (wenum >= MAXWILDATOM) {
					printf
						("\nnumber of wild elements exceeds MAXWILDATOM, increase MAXWILDATOM");
					exit(1);
				}
			}
			if (strncmp("ATD", line, 3) == 0) {
				sscanf(line, "%s%s%s%s%s%s%s%s%s%s", tmpchar,
					   keyword0, keyword1, keyword2, keyword3, keyword4,
					   keyword5, keyword6, keyword7, keyword8);
/* resideu name */
				if (strcmp(keyword1, "&") == 0) {
					strcpy(atom[i].ambername, keyword0);
					break;
				}
				resindex = 0;
				if (strcmp(keyword1, "*") == 0
					|| strcmp(keyword2, atom[i].aa) == 0)
					resindex = 1;
				if (strcmp(keyword1, "AA") == 0
					|| strcmp(keyword2, "BIO") == 0) {
					if (strcmp(atom[i].aa, "ALA") == 0)
						resindex = 1;
					if (strcmp(atom[i].aa, "GLY") == 0)
						resindex = 1;
					if (strcmp(atom[i].aa, "SER") == 0)
						resindex = 1;
					if (strcmp(atom[i].aa, "THR") == 0)
						resindex = 1;
					if (strcmp(atom[i].aa, "LEU") == 0)
						resindex = 1;
					if (strcmp(atom[i].aa, "ILE") == 0)
						resindex = 1;
					if (strcmp(atom[i].aa, "VAL") == 0)
						resindex = 1;
					if (strcmp(atom[i].aa, "ASN") == 0)
						resindex = 1;
					if (strcmp(atom[i].aa, "GLN") == 0)
						resindex = 1;
					if (strcmp(atom[i].aa, "ARG") == 0)
						resindex = 1;
					if (strcmp(atom[i].aa, "HID") == 0)
						resindex = 1;
					if (strcmp(atom[i].aa, "HIE") == 0)
						resindex = 1;
					if (strcmp(atom[i].aa, "HIP") == 0)
						resindex = 1;
					if (strcmp(atom[i].aa, "TRP") == 0)
						resindex = 1;
					if (strcmp(atom[i].aa, "PHE") == 0)
						resindex = 1;
					if (strcmp(atom[i].aa, "TYR") == 0)
						resindex = 1;
					if (strcmp(atom[i].aa, "GLU") == 0)
						resindex = 1;
					if (strcmp(atom[i].aa, "ASP") == 0)
						resindex = 1;
					if (strcmp(atom[i].aa, "LYS") == 0)
						resindex = 1;
					if (strcmp(atom[i].aa, "PRO") == 0)
						resindex = 1;
					if (strcmp(atom[i].aa, "CYS") == 0)
						resindex = 1;
					if (strcmp(atom[i].aa, "CYX") == 0)
						resindex = 1;
					if (strcmp(atom[i].aa, "MET") == 0)
						resindex = 1;
					if (strcmp(atom[i].aa, "ASH") == 0)
						resindex = 1;
					if (strcmp(atom[i].aa, "GLH") == 0)
						resindex = 1;
					if (strcmp(atom[i].aa, "LYH") == 0)
						resindex = 1;
				}
				if (strcmp(keyword1, "NA") == 0
					|| strcmp(keyword2, "BIO") == 0) {
					if (strcmp(atom[i].aa, "DADE") == 0)
						resindex = 1;
					if (strcmp(atom[i].aa, "RADE") == 0)
						resindex = 1;
					if (strcmp(atom[i].aa, "DTHY") == 0)
						resindex = 1;
					if (strcmp(atom[i].aa, "RURA") == 0)
						resindex = 1;
					if (strcmp(atom[i].aa, "DGUA") == 0)
						resindex = 1;
					if (strcmp(atom[i].aa, "RGUA") == 0)
						resindex = 1;
					if (strcmp(atom[i].aa, "DCYT") == 0)
						resindex = 1;
					if (strcmp(atom[i].aa, "RCYT") == 0)
						resindex = 1;
				}
				if (resindex == 0)
					continue;

/* atomic number*/
				if (strcmp(keyword2, "&") == 0) {
					strcpy(atom[i].ambername, keyword0);
					break;
				}
				tmpint1 = atoi(keyword2);
				if (tmpint1 != atom[i].atomicnum)
					continue;
/* connected atomic number*/
				if (strcmp(keyword3, "&") == 0) {
					strcpy(atom[i].ambername, keyword0);
					break;
				}
				if (strcmp(keyword3, "*") != 0) {
					tmpint1 = atoi(keyword3);
					if (tmpint1 != atom[i].connum)
						continue;
				}


/* the number of connected hydrogen */
				if (strcmp(keyword4, "&") == 0) {
					strcpy(atom[i].ambername, keyword0);
					break;
				}
				if (strcmp(keyword4, "*") != 0) {
					tmpint1 = atoi(keyword4);
					tmpint2 = 0;
					for (j = 0; j < 6; j++) {
						tmpint = atom[i].con[j];
						if (tmpint == -1)
							break;
						if (atom[tmpint].atomicnum == 1)
							tmpint2++;
					}
					if (tmpint1 != tmpint2)
						continue;
				}

/* electron withdram number */
				if (strcmp(keyword5, "&") == 0) {
					strcpy(atom[i].ambername, keyword0);
					break;
				}
				if (strcmp(keyword5, "*") != 0) {
					tmpint = atom[i].con[0];
					tmpint1 = atoi(keyword5);
					tmpint2 = 0;
					for (j = 0; j < 6; j++) {
						if (atom[tmpint].con[j] == -1)
							break;
						if (atom[atom[tmpint].con[j]].ewd == 1)
							tmpint2++;
					}
					if (tmpint1 != tmpint2)
						continue;
				}


/*atomic properity */
				if (strcmp(keyword6, "&") == 0) {
					strcpy(atom[i].ambername, keyword0);
					break;
				}
				if (strcmp(keyword6, "*") != 0) {
					tmpint = apcheck(i, -1, keyword6);
					if (tmpint == 0)
						continue;
				}
/* special chemical environment*/
				if (strcmp(keyword7, "&") == 0) {
					strcpy(atom[i].ambername, keyword0);
					break;
				}
				ces_bond_index = 1;
				if (strcmp(keyword8, "&") == 0)
					ces_bond_index = 0;
				tmpint = 0;
				if (debug == 1)
					printf
						("\nBegin a new chemical environment defination for atom %d (%s) %s %s",
						 i + 1, atom[i].name, keyword7, keyword8);
				caid = i;
				tmpint = jatspecial(caid, keyword7, keyword8);
				if (tmpint == 1) {
					strcpy(atom[i].ambername, keyword0);
					break;
				} else
					continue;
			}
		}
		if (weindex == 1)
			weindex = 0;
		rewind(fpin);
	}
	fclose(fpin);
}

void atadjust(void)
{
	int i;
	int index = 0;
	int num = 0;
	int *index1;
	int *index2;
	int bondi, bondj;
	int flag;

        index1 = (int *) malloc(sizeof(int) * (atomnum + 10));
        if (index1 == NULL) {
                fprintf(stdout, "memory allocation error for *index1 in atadjust()\n");
                exit(1);
        }
        index2 = (int *) malloc(sizeof(int) * (atomnum + 10));
        if (index2 == NULL) {
                fprintf(stdout, "memory allocation error for *index2 in atadjust()\n");
                exit(1);
        }

	for (i = 0; i < atomnum; i++) {
		index1[i] = 0;
		index2[i] = 0;
	}
	for (i = 0; i < atomnum; i++)
		if (strcmp(atom[i].ambername, "cc") == 0 ||
			strcmp(atom[i].ambername, "ce") == 0 ||
			strcmp(atom[i].ambername, "cg") == 0 ||
			strcmp(atom[i].ambername, "pc") == 0 ||
			strcmp(atom[i].ambername, "pe") == 0 ||
			strcmp(atom[i].ambername, "nc") == 0 ||
			strcmp(atom[i].ambername, "ne") == 0) {
			index2[i] = 1;
			if (index == 0) {
				index1[i] = 1;
				index = 1;
			}
			num++;
		}
	if (num == 0)
		return;
	num--;
	while (num > 0) {
		num--;
		flag = 0;
		for (i = 0; i < bondnum; i++) {
			bondi=bond[i].bondi;
			bondj=bond[i].bondj;
			if ((index2[bondi] + index2[bondj]) != 2) continue;
			if (flag ==0 && index1[bondi] == 0 && index1[bondj] == 0 ) index1[bondi] =1;
			if (index1[bondi] == 0 && index1[bondj] != 0) {
				flag = 1;
				if (bond[i].type == 1 || bond[i].type == 7)
					index1[bondi] = index1[bondj];
				if (bond[i].type == 2 || bond[i].type == 8 || bond[i].type == 3)
					index1[bondi] = -index1[bondj];
			}
			if (index1[bondj] == 0 && index1[bondi] != 0 ) {
				flag = 1;
				if (bond[i].type == 1 || bond[i].type == 7)
					index1[bondj] = index1[bondi];
				if (bond[i].type == 2 || bond[i].type == 8 || bond[i].type == 3)
					index1[bondj] = -index1[bondi];
			}
		}
	}
	for (i = 0; i < atomnum; i++)
		if (index1[i] == -1) {
			if (strcmp(atom[i].ambername, "cc") == 0)
				strcpy(atom[i].ambername, "cd");
			if (strcmp(atom[i].ambername, "ce") == 0)
				strcpy(atom[i].ambername, "cf");
			if (strcmp(atom[i].ambername, "cg") == 0)
				strcpy(atom[i].ambername, "ch");
			if (strcmp(atom[i].ambername, "pc") == 0)
				strcpy(atom[i].ambername, "pd");
			if (strcmp(atom[i].ambername, "pe") == 0)
				strcpy(atom[i].ambername, "pf");
			if (strcmp(atom[i].ambername, "nc") == 0)
				strcpy(atom[i].ambername, "nd");
			if (strcmp(atom[i].ambername, "ne") == 0)
				strcpy(atom[i].ambername, "nf");
		}
	free(index1);
	free(index2);
}

void cpadjust(void)
{
	int i;
	int index = 0;
	int num = 0;
	int *index1;
	int *index2;
	int bondi, bondj;
        index1 = (int *) malloc(sizeof(int) * (atomnum + 10));
        if (index1 == NULL) {
                fprintf(stdout, "memory allocation error for *index1 in cpadjust()\n");
                exit(1);
        }
        index2 = (int *) malloc(sizeof(int) * (atomnum + 10));
        if (index2 == NULL) {
                fprintf(stdout, "memory allocation error for *index2 in cpadjust()\n");
                exit(1);
        }


	for (i = 0; i < atomnum; i++) {
		index1[i] = 0;
		index2[i] = 0;
	}
	for (i = 0; i < atomnum; i++)
		if (strcmp(atom[i].ambername, "cp") == 0) {
			index2[i] = 1;
			if (index == 0) {
				index1[i] = 1;
				index = 1;
			}
			num++;
		}
	if (num == 0)
		return;
	num--;

	while (num > 0) {
		num--;
		for (i = 0; i < bondnum; i++) {
			bondi=bond[i].bondi;
			bondj=bond[i].bondj;
			if ((index2[bondi] + index2[bondj]) != 2) continue;
			if (index1[bondi] == 0 && index1[bondj] != 0) {
				if (bond[i].type == 1) 
					index1[bondi] = index1[bondj];
				else			
					index1[bondi] = -index1[bondj];
			}
			if (index1[bondj] == 0 && index1[bondi] != 0 ) {
				if (bond[i].type == 1)
					index1[bondj] = index1[bondi];
				else	 	
					index1[bondj] = -index1[bondi];
			}
		}
	}
	for (i = 0; i < atomnum; i++)
		if (index1[i] == -1 && strcmp(atom[i].ambername, "cp") == 0)
				strcpy(atom[i].ambername, "cq");
	free(index1);
	free(index2);
}

void error() {
	int *index;
	int bondi, bondj;
	int i;
 
        index = (int *) malloc(sizeof(int) * (atomnum + 10));
        if (index == NULL) {
                fprintf(stdout, "memory allocation error for *index in error()\n");
                exit(1);
        }

	for (i = 0; i < atomnum; i++) {
		index[i] =0;
		if(strcmp(atom[i].ambername, "cc") == 0 ||
		   strcmp(atom[i].ambername, "ce") == 0 ||
		   strcmp(atom[i].ambername, "nc") == 0 ||
		   strcmp(atom[i].ambername, "ne") == 0 ||
		   strcmp(atom[i].ambername, "pc") == 0 ||
		   strcmp(atom[i].ambername, "pe") == 0 ||
		   strcmp(atom[i].ambername, "cg") == 0) 
			index[i] = 1;

		if(strcmp(atom[i].ambername, "cd") == 0 ||
		   strcmp(atom[i].ambername, "cf") == 0 ||
		   strcmp(atom[i].ambername, "nd") == 0 ||
		   strcmp(atom[i].ambername, "nf") == 0 ||
		   strcmp(atom[i].ambername, "pd") == 0 ||
		   strcmp(atom[i].ambername, "pf") == 0 ||
		   strcmp(atom[i].ambername, "ch") == 0) 
			index[i] = -1;

	}
	for (i = 0; i < bondnum; i++) {
		bondi=bond[i].bondi;
		bondj=bond[i].bondj;

		if((index[bondi]== 1 &&index[bondj]==1) || (index[bondi] == -1 && index[bondj]==-1))  	
			if(bond[i].type!=1 && bond[i].type!=7) 
				printf("\nWARNING: atom type of %5s (%s) and %5s (%s) may be wrong", atom[bondi].name, atom[bondi].ambername,
					atom[bondj].name, atom[bondj].ambername); 

		if(index[bondi]*index[bondj]==-1)  	
			if(bond[i].type==1 || bond[i].type==7) 
				printf("\nWARNING: atom type of %5s (%s) and %5s (%s) may be wrong", atom[bondi].name, atom[bondi].ambername,
					atom[bondj].name, atom[bondj].ambername); 

		if((strcmp(atom[bondi].ambername, "cp")==0 && strcmp(atom[bondj].ambername, "cp")==0) ||
		   (strcmp(atom[bondi].ambername, "cq")==0 && strcmp(atom[bondj].ambername, "cq")==0)) 
			if(bond[i].type!=1) 
				printf("\nWARNING: atom type of %5s (%s) and %5s (%s) may be wrong", atom[bondi].name, atom[bondi].ambername,
					atom[bondj].name, atom[bondj].ambername); 

		if((strcmp(atom[bondi].ambername, "cp")==0 && strcmp(atom[bondj].ambername, "cq")==0) ||
		   (strcmp(atom[bondi].ambername, "cq")==0 && strcmp(atom[bondj].ambername, "cp")==0)) 
			if(bond[i].type==1) 
				printf("\nWARNING: atom type of %5s (%s) and %5s (%s) may be wrong", atom[bondi].name, atom[bondi].ambername,
					atom[bondj].name, atom[bondj].ambername); 
	}
	free(index);
}

int main(int argc, char *argv[])
{
	int index1, index2;
	int i;
	int overflow_flag = 0;			/*if overflow_flag ==1, reallocate memory */
	char *fname;

    amberhome = (char *) getenv("AMBERHOME");
    if( amberhome == NULL ){
       fprintf( stdout, "AMBERHOME is not set!\n" );
       exit(1);
    }
	if (strcmp(COLORTEXT, "YES") == 0 || strcmp(COLORTEXT, "yes") == 0) {
		if (argc == 2
			&& (strcmp(argv[1], "-h") == 0
				|| strcmp(argv[1], "-H") == 0)) {
			printf("[31mUsage: atomtype -i[0m input file name\n"
				   "[31m                -o[0m output file name(ac)\n"
				   "[31m                -f[0m input file format(ac (the default) or mol2)\n"
				   "[31m                -p[0m amber or gaff or bcc or gas or sybyl, it is supressed by \"-d\" option\n"
				   "[31m                -d[0m atom type defination file, optional\n"
				   "[31m                -a[0m do post atom type adjustment (it is or sybyl applied with \"-d\" option)\n"
				   "                   1: yes, 0: no (the default)\n");
			exit(1);
		}
		if (argc != 15 && argc != 13 && argc != 11 && argc != 9 && argc != 7 && argc != 5
			&& argc != 3) {
			printf("[31mUsage: atomtype -i[0m input file name\n"
				   "[31m                -o[0m output file name (ac)\n"
				   "[31m                -f[0m input file format(ac (the default) or mol2)\n"
				   "[31m                -p[0m amber or gaff or bcc or gas or sybyl, it is supressed by \"-d\" option\n"
				   "[31m                -d[0m atom type defination file, optional\n"
				   "[31m                -a[0m do post atom type adjustment (it is applied with \"-d\" option)\n"
				   "                   1: yes, 0: no (the default)\n");
			exit(1);
		}
	} else {
		if (argc == 2)
			if (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "-H") == 0) {
				printf("\n Usage atomtype -i input file name");
				printf("\n                -o output file name (ac)");
				printf
					("\n                -f input file format (ac (default) or mol2)");
				printf
					("\n                -p amber or gaff or bcc or gas or sybyl, it is supressed by \"-d\" option");
				printf
					("\n                -d atom type defination file, optional\n ");
				printf
					("\n                -a do post atom type adjustment (it is applied with \"-d\" option\n"); 
				printf
					("\n                   1: yes, 0: no (the default)\n"); 
				/*
				   printf("\n Usage atomtype -i   inputfile ");   
				   printf("\n              -bcc yes_or_no, optional "); 
				   printf("\n              -def atom_type_defination_file, optional ");
				   printf("\n              -inf information_file, optional ");
				   printf("\n              -an  elem_or_type_or_orig (default is orig), optional \n");
				 */
				exit(1);
			}
		if (argc != 15 && argc != 13 && argc != 11 && argc != 9 && argc != 7 && argc != 5
			&& argc != 3) {
			printf("\n Usage atomtype -i input file name(ac)");
			printf("\n                -o output file name (ac)");
			printf
				("\n                -f input file format (ac (default) or mol2)");
			printf
				("\n                -p amber or gaff or bcc or gas or sybyl, it is supressed by \"-d\" option");
			printf
				("\n                -d atom type defination file, optional\n ");
			printf
				("\n                -a do post atom type adjustment (it is applied with \"-d\" option\n"); 
			printf
				("\n                   1: yes, 0: no (the default)\n"); 
			exit(1);
		}
	}
	index1 = 0;
	index2 = 0;
	ipostadjustment = 0;
	for (i = 1; i < argc; i += 2) {
		if (strcmp(argv[i], "-i") == 0) {
			strcpy(ifilename, argv[i + 1]);
			index1++;
		}
		if (strcmp(argv[i], "-o") == 0) {
			strcpy(ofilename, argv[i + 1]);
			index1++;
		}
		if (strcmp(argv[i], "-f") == 0)
			if (strcmp("mol2", argv[i + 1]) == 0
				|| strcmp("MOL2", argv[i + 1]) == 0)
				findex = 1;
		if (strcmp(argv[i], "-inf") == 0)
			strcpy(minfo.inf_filename, argv[i + 1]);
		if (strcmp(argv[i], "-p") == 0) {
			if (strcmp("amber", argv[i + 1]) == 0
				|| strcmp("AMBER", argv[i + 1]) == 0)
				pindex = 1;
			if (strcmp("gaff", argv[i + 1]) == 0
				|| strcmp("GAFF", argv[i + 1]) == 0)
				pindex = 0;
			if (strcmp("gff", argv[i + 1]) == 0
				|| strcmp("GFF", argv[i + 1]) == 0)
				pindex = 0;
			if (strcmp("bcc", argv[i + 1]) == 0
				|| strcmp("BCC", argv[i + 1]) == 0)
				pindex = 2;
			if (strcmp("am1bcc", argv[i + 1]) == 0
				|| strcmp("AM1BCC", argv[i + 1]) == 0)
				pindex = 2;
			if (strcmp("gas", argv[i + 1]) == 0
				|| strcmp("GAS", argv[i + 1]) == 0)
				pindex = 3;
			if (strcmp("sybyl", argv[i + 1]) == 0
				|| strcmp("SYBYL", argv[i + 1]) == 0)
				pindex = 4;
		}

		if (strcmp(argv[i], "-an") == 0) {
			if (strcmp("orig", argv[i + 1]) == 0
				|| strcmp("ORIG", argv[i + 1]) == 0)
				anindex = 0;
			if (strcmp("elem", argv[i + 1]) == 0
				|| strcmp("ELEM", argv[i + 1]) == 0)
				anindex = 1;
			if (strcmp("type", argv[i + 1]) == 0
				|| strcmp("TYPE", argv[i + 1]) == 0)
				anindex = 2;
		}
		if (strcmp(argv[i], "-d") == 0) {
			strcpy(atomtype_def, argv[i + 1]);
			index2 = 1;
		}
		if (strcmp(argv[i], "-a") == 0) {
			ipostadjustment = atoi(argv[i + 1]);
			if(ipostadjustment != 0 && ipostadjustment != 1) {
				fprintf(stdout, "The -a flag must be 0 or 1\n"); 
				exit(1);
			}
		}	
	}
	if (index1 != 2) {
		printf("\n Usage atomtype -i inputfile (ac)");
		printf("\n                -o outputfile (ac)");
		printf("\n                -f file format (ac (default) or mol2)");
		printf
			("\n                -p amber or gaff, it is supressed by \"-d\" option");
		printf
			("\n                -d atom_type_defination_file, optional\n ");
		printf
			("\n                -a do post atom type adjustment (it is applied with \"-d\" option\n"); 
		printf
			("\n                   1: yes, 0: no (the default)\n"); 
		exit(1);
	}


	if (index2 == 0) {
		if(pindex == 0 ) 
			ipostadjustment = 1;
		else
			ipostadjustment = 0;
		if (pindex == 0)
			fname = "ATOMTYPE_GFF.DEF";
		else if (pindex == 1)
			fname = "ATOMTYPE_AMBER.DEF";
		else if (pindex == 2)
			fname = "ATOMTYPE_BCC.DEF";
		else if (pindex == 3)
			fname = "ATOMTYPE_GAS.DEF";
		else if (pindex == 4)
			fname = "ATOMTYPE_SYBYL.DEF";
		else {
			fprintf (stderr,
		 	"Cannot find atom type defination file!, define with \"-def\" option\n");
			fname = NULL;
		}
		if (fname != NULL)
		    build_dat_path(atomtype_def, fname, sizeof atomtype_def, 0);
	}

	default_cinfo(&cinfo);
	default_minfo(&minfo);
/*	memory allocation */
	atom = (ATOM *) malloc(sizeof(ATOM) * cinfo.maxatom);
	if (atom == NULL) {
		fprintf(stdout, "memory allocation error for *atom\n");
		exit(1);
	}
	arty = (AROM *) malloc(sizeof(AROM) * cinfo.maxatom);
	if (arty == NULL) {
		fprintf(stdout, "memory allocation error for *arty\n");
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

	if (findex == 1)
		overflow_flag =
			rmol2(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo,
				  0);
	else
		overflow_flag =
			rac(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo);
	if (overflow_flag) {
		cinfo.maxatom = atomnum + 10;
		cinfo.maxbond = bondnum + 10;
		free(atom);
		free(arty);
		free(bond);
		atom = (ATOM *) malloc(sizeof(ATOM) * cinfo.maxatom);
		if (atom == NULL) {
			fprintf(stdout, "memory allocation error for *atom\n");
			exit(1);
		}
		arty = (AROM *) malloc(sizeof(AROM) * cinfo.maxatom);
		if (arty == NULL) {
			fprintf(stdout, "memory allocation error for *arty\n");
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
		if (findex == 1)
			rmol2(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo,
				  0);
		else
			rac(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo);
	}

	atomicnum(atomnum, atom);
	adjustatomname(atomnum, atom, 1);

	selectchain = (int *) malloc(sizeof(int) * (atomnum + 10));
	if (selectchain == NULL) {
		fprintf(stdout, "memory allocation error for *selectchain\n");
		exit(1);
	}

	selectindex = (int *) malloc(sizeof(int) * (atomnum + 10));
	if (selectindex == NULL) {
		fprintf(stdout, "memory allocation error for *selectindex\n");
		exit(1);
	}

	sb = (int *) malloc(sizeof(int) * (atomnum + 10));
	if (sb == NULL) {
		fprintf(stdout, "memory allocation error for *sb\n");
		exit(1);
	}

	SB = (int *) malloc(sizeof(int) * (atomnum + 10));
	if (SB == NULL) {
		fprintf(stdout, "memory allocation error for *SB\n");
		exit(1);
	}
	db = (int *) malloc(sizeof(int) * (atomnum + 10));
	if (db == NULL) {
		fprintf(stdout, "memory allocation error for *db\n");
		exit(1);
	}
	DB = (int *) malloc(sizeof(int) * (atomnum + 10));
	if (DB == NULL) {
		fprintf(stdout, "memory allocation error for *DB\n");
		exit(1);
	}
	DL = (int *) malloc(sizeof(int) * (atomnum + 10));
	if (DL == NULL) {
		fprintf(stdout, "memory allocation error for *DL\n");
		exit(1);
	}
	tb = (int *) malloc(sizeof(int) * (atomnum + 10));
	if (tb == NULL) {
		fprintf(stdout, "memory allocation error for *tb\n");
		exit(1);
	}
	TB = (int *) malloc(sizeof(int) * (atomnum + 10));
	if (TB == NULL) {
		fprintf(stdout, "memory allocation error for *TB\n");
		exit(1);
	}
	AB = (int *) malloc(sizeof(int) * (atomnum + 10));
	if (AB == NULL) {
		fprintf(stdout, "memory allocation error for *AB\n");
		exit(1);
	}
	nr = (int *) malloc(sizeof(int) * (atomnum + 10));
	if (nr == NULL) {
		fprintf(stdout, "memory allocation error for *nr\n");
		exit(1);
	}
	atomcesname = (NAME *) malloc(sizeof(NAME) * (atomnum + 10));
	if (atomcesname == NULL) {
		fprintf(stdout, "memory allocation error for *atomcesname\n");
		exit(1);
	}

	bondinfo();


	if (pindex != 2)
		overflow_flag = ringdetect(atomnum, atom, bondnum, bond, &ringnum, ring, arty,
				   cinfo.maxatom, cinfo.maxring, minfo.inf_filename, 0);
	else
		overflow_flag = ringdetect(atomnum, atom, bondnum, bond, &ringnum, ring, arty,
				   cinfo.maxatom, cinfo.maxring, minfo.inf_filename, 1);
	if(overflow_flag == 1) {
                cinfo.maxring = ringnum + 10;
                free(ring);
        	ring = (RING *) malloc(sizeof(RING) * cinfo.maxring);
        	if (ring == NULL) {
                	fprintf(stdout, "memory allocation error for *ring\n");
                	exit(1);
        	}
		if (pindex != 2)
			overflow_flag = ringdetect(atomnum, atom, bondnum, bond, &ringnum, ring, arty,
				   cinfo.maxatom, cinfo.maxring, minfo.inf_filename, 0);
		else
			overflow_flag = ringdetect(atomnum, atom, bondnum, bond, &ringnum, ring, arty,
				   cinfo.maxatom, cinfo.maxring, minfo.inf_filename, 1);
	}
	jat();
	improper();
	if (anindex != 0)
		atname();
	if(ipostadjustment == 1) {
		atadjust();
		cpadjust();
	}
	error();
	wac(ofilename, atomnum, atom, bondnum, bond, cinfo, minfo);
/*
	free(atom);
	free(arty);
	free(bond);
	free(selectchain);
	free(selectindex);
	free(sb);
	free(SB);
	free(db);
	free(DB);
	free(DL);
	free(tb);
	free(TB);
	free(AB);
	free(nr);
	free(wildatom);
	free(schain);
	free(ssindex);
*/
/*	printf("\n"); */
	return (0);
}
