char *amberhome;
# define MAXCORR 256 
# define debug 0
# include "common.h"
# include "define.h"
# include "atom.h"
# include "utility.c"
# include "common.c"
# include "ring.c"
# include "rotate.c"
# include "ac.c"
# include "charmm.c"
# include "mol2.c"
# include "mopcrt.c"
# include "divcrt.c"
# include "mopint.c"
# include "mopout.c"
# include "divout.c"
# include "sqmcrt.c"
# include "sqmout.c"
# include "gesp.c"
# include "gcrt.c"
# include "gzmat.c"
# include "gout.c"
# include "pdb.c"
# include "csd.c"
# include "mdl.c"
# include "alc.c"
# include "hin.c"
# include "prep.c"
# include "rst.c"
# include "jzmat.c"
# include "jcrt.c"
# include "jout.c"
MOLINFO minfo2;

static int overflow_flag = 0;			/*if overflow_flag ==1, reallocate memory */
int atomtype_flag;                          /*judge atom type? */
int bondtype_flag;                          /*judge bond type? */
int default_flag;                           /*assign default information? */
int atomname_flag;                          /*assign atom name? */
int atomicnum_flag;                         /*judge atomic number according to atom name ? */
int adjustatomname_flag;            /*adjust atom name? */
int duplicatedname_flag;            /*check atom name duplication? */
int cartcoord_flag;                         /*generate coordinate from internal coordinate ? */
int connect_flag;                           /*judge atom connectivity and generate bond ? */

int atomnum = 0;
int bondnum = 0;
int ringnum = 0;
ATOM *atom;
BOND *bond;
RING *ring;
AROM *arom;
MOLINFO minfo;
CONTROLINFO cinfo;

char line[MAXCHAR];
char ifilename[MAXCHAR];

FILE *fpin;
FILE *fpout;

typedef struct {
        char type[10];
        char element[10];
	int connum;
} CORR;

typedef struct {
        double low;
        double high;
	double bond;
	double offset;
	int valid;
} BONDLENGTH;

int iatomtype = 0;
int ibondtype = 0;

double bondradius[120];
BONDLENGTH bl;

void bondrange(int atomicnum1, int atomicnum2) {
double bond;
double offset;
bl.low = 0;
bl.high = 0;
if(bondradius[atomicnum1] > 0 &&  bondradius[atomicnum2] > 0) {
	bond = bondradius[atomicnum1] + bondradius[atomicnum2];
	if (bond <= 1.5) offset = bond * 0.15;
	if (bond > 1.5 && bond <= 1.90) offset = bond * 0.11;
	if (bond > 1.90 && bond <= 2.05) offset = bond * 0.09;
	if (bond > 2.05) offset = bond * 0.08;
	bl.offset = offset;
	bl.bond = bond;
/* definition in connect() */
/*
	bl.low = bond * 0.5;
	bl.high = bond + offset;
*/
/* now apply more stringent definition */
	bl.low = bond * 0.65;
	bl.high = bond + offset * 0.75;
	bl.valid = 1;
}
else
	bl.valid = 0;
fflush(stdout);
}



void formatchk(char *filename, char *format) {
int imol, iatom, ibond; /* for mol2 format */
char tmpchar[MAXCHAR];
char tmpchar1[MAXCHAR];
char tmpchar2[MAXCHAR];
char tmpchar3[MAXCHAR];
char FIVEBLANKS[]  = "     ";
char FOURBLANKS[]  = "    ";
char THREEBLANKS[] = "   ";
int ieps; /* for gout format */

FILE *fpin;
if ((fpin = fopen(filename, "r")) == NULL) {
        fprintf(stdout, "Cannot open the input file: %s, exit\n", filename);
        exit(1);
}

if(strcmp(format, "mol2") == 0) {
	imol = 0;
	iatom = 0;
	ibond = 0;
	printf("\n\n[34m-- Check Format for mol2 File --[0m\n");
	for (;;) {
                if (fgets(line, MAXCHAR, fpin) == NULL)  break;
		sscanf(line, "%s", tmpchar);
		if(strcmp(tmpchar, "@<TRIPOS>MOLECULE") == 0) imol = 1;
		if(strcmp(tmpchar, "@<TRIPOS>ATOM") == 0) iatom = 1;
		if(strcmp(tmpchar, "@<TRIPOS>BOND") == 0) ibond = 1;
		if(imol+iatom+ibond == 3) break;
	}
	if(imol == 0) {
			fprintf(stdout, "[31mError: no @@<TRIPOS>MOLECULE field, exit\n[0m");
			exit(1);
		}
	if(iatom == 0) {
			fprintf(stdout, "[31mError: no @@<TRIPOS>ATOM field, exit\n[0m");
			exit(1);
		}
	if(ibond == 0) {
			fprintf(stdout, "[31mError: no @@<TRIPOS>BOND field, exit\n[0m");
			exit(1);
		}

	printf("\n\n[32m-- Status: pass --[0m\n");
}
else if(strcmp(format, "mopint") == 0) {

}
else if(strcmp(format, "mopcrt") == 0) {

}
else if(strcmp(format, "mopout") == 0) {

}
else if(strcmp(format, "gcrt") == 0) {

}
else if(strcmp(format, "gzmat") == 0) {

}
else if(strcmp(format, "gout") == 0) {
	printf("\n\n[34m-- Check Format for Gaussian Output File --[0m\n");
	ieps = 0;
	for (;;) {
                if (fgets(line, MAXCHAR, fpin) == NULL)  break;
		sscanf(line, "%s%s%s", tmpchar1, tmpchar2, tmpchar3);
		if(strcmp(tmpchar1, "ESP") == 0 && strcmp(tmpchar2, "Fit") == 0 && strcmp(tmpchar3, "Center") == 0) {
			ieps = 1;
			break;
		}
	}
	if(ieps == 0) {
		fprintf(stdout, "[33mWarning: No ESP information in the gaussian output file, this file cannot be used to generate RESP charge\n[0m");
		fprintf(stdout, "         Using keywords like \"#HF/6-31G* SCF=tight Test Pop=MK iop(6/33=2) iop(6/42=6) opt\" in gaussian calculations\n");
	}
	printf("\n\n[32m-- Status: pass --[0m\n");
}
else if(strcmp(format, "jcrt") == 0) {

}
else if(strcmp(format, "jzmat") == 0) {

}
else if(strcmp(format, "jout") == 0) {

}
else if(strcmp(format, "pdb") == 0 || strcmp(format, "mpdb") == 0 || strcmp(format, "ac") == 0) {
	printf("\n\n[34m-- Check Format for pdb,mpdb or ac File --[0m\n");
	for (;;) {
                if (fgets(line, MAXCHAR, fpin) == NULL)  break;
		if(strncmp(line, "ATOM", 4) == 0) {
			if( strncmp( & line[6], FIVEBLANKS, sizeof(FIVEBLANKS) - 1 ) == 0 ) {
				fprintf(stdout, "[31mError: Atom IDs must in Columns 7-11\n[0m");
				exit(1);
			}
			if( strncmp( & line[12], FOURBLANKS, sizeof(FOURBLANKS) - 1 ) == 0 ) {
				fprintf(stdout, "[31mError: Atom Names must in Columns 13-16\n[0m");
				exit(1);
			}
			if( strncmp( & line[17], THREEBLANKS, sizeof(THREEBLANKS) - 1 ) == 0 ) {
				fprintf(stdout, "[31mError: Residue Names must in Columns 18-20\n[0m");
				exit(1);
			}
			if( strncmp( & line[22], FOURBLANKS, sizeof(FOURBLANKS) - 1 ) == 0 ) {
				fprintf(stdout, "[31mError: Residue IDs must in Columns 23-26\n[0m");
				exit(1);
			}
			if(line[34]!='.' || line[42] != '.' || line[50]!='.') {
				fprintf(stdout, "[31mError: Coordinates must be in Columns 31-38, 39-46 and 47-54 in %%8.3format\n[0m");
				exit(1);
			}

		}
	}
	printf("\n\n[32m-- Status: pass --[0m\n");
}
else if(strcmp(format, "csd") == 0) {

}
else if(strcmp(format, "mdl") == 0 || strcmp(format, "sdf") == 0) {
	printf("\n\n[34m-- Check Format for mdl or sdf File --[0m\n");
	fgets(line, MAXCHAR, fpin) ; 
	fgets(line, MAXCHAR, fpin) ; 
	if(line[20] == '2' && (line[21] == 'D' || line[21] == 'd')) {
		fprintf(stdout, "[31mError: antechamber cannot handle 2D-mdl or sdf format, try to fill in open valences\n[0m");
		exit(1);	
	}
	fgets(line, MAXCHAR, fpin) ; 
	fgets(line, MAXCHAR, fpin) ; 
	fgets(line, MAXCHAR, fpin) ; 
	if(line[6] >= '0' && line[6] <= '9') {	
		fprintf(stdout, "[31mError: Too many atoms and bonds, MDL sdf format can only hold molecules that have less than 999 atoms and 999 bonds, exit\n[0m");
        	exit(1);
	}
	printf("\n\n[32m-- Status: pass --[0m\n");
}
else if(strcmp(format, "alc") == 0) {

}
else if(strcmp(format, "hin") == 0) {

}
else if(strcmp(format, "prepi") == 0) {

}
else if(strcmp(format, "prepc") == 0) {

}
else if(strcmp(format, "divcrt") == 0) {

}
else if(strcmp(format, "divout") == 0) {

}
}

void molchk() {
int i,j,k;
int numh = 0;
int numn = 0;
int numo = 0;
int numc = 0;
int nums = 0;
int nump = 0;
int numx = 0; /* number of halogens */
FILE *fp;
char filename[MAXCHAR];
char line[MAXCHAR];
CORR *corr; 
int corrnum = 0;
int tvalence ;
int val_ar;
int connum;
int tmpint;
int selectflag; /*used by unit check*/
int errorflag;
int numselect;
int bondflag;
double dist;
double radius;

/*part1: element */
printf("\n\n[34m-- Check Unusual Elements --[0m\n");
for(i=0;i<atomnum;i++){
	if(atom[i].atomicnum == 1) numh++;
	else if (atom[i].atomicnum == 9 || atom[i].atomicnum == 17 || atom[i].atomicnum == 35 || atom[i].atomicnum == 53)
		numx ++;
	else if (atom[i].atomicnum == 6) numc ++; 
	else if (atom[i].atomicnum == 7) numn ++; 
	else if (atom[i].atomicnum == 8) numo ++; 
	else if (atom[i].atomicnum == 15) nump ++; 
	else if (atom[i].atomicnum == 16) nums ++; 
	else  
		fprintf(stdout, "[33mWarning: unusual element for ATOM %d (%s) : %s\n[0m", i+1, atom[i].name, atom[i].element);
}
if((atomnum - numh - numn - numo -numc -nums -nump - numx) > 0){
	fprintf(stdout, "[33m         GAFF does not have sufficient parameters for molecules having unusual [0m\n");
	fprintf(stdout, "[33m         elements (other than H,C,N,O,S,P and halogens[0m\n");
	fprintf(stdout, "[33m         To gauentee antechamber package working soundly, one may need to designate[0m\n");
	fprintf(stdout, "[33m         bond types for bonds involved in unusual elements[0m\n");
	fprintf(stdout, "[33m         To do so, simply frozen the bond types by appending \"F\" or \"f\" to the[0m\n");
	fprintf(stdout, "[33m         corresponding bond types in ac or mol2 files)[0m\n");
}
else if (atomnum == numh+numn+numo+numc+nums+nump+numx) {
	printf("\n\n[32m-- Status: pass --[0m\n");
}

/*part2: unfilled valences */
printf("\n\n[34m-- Check Open Valences --[0m\n");
/* how many hydrogen and halogen atoms */
if(numh+numx == 0) 
	fprintf(stdout, "[33mWarning: There are no hydrogen and halogens in this molecule, it is highly possible that this molecule has open valence\n[0m"); 

/* check if atom[].connum consistant with that in CORR_NAME_TYPE.DAT */
if(iatomtype == 1) {
	corr = (CORR *) malloc(sizeof(CORR) * MAXCORR);
	if (corr == NULL) {
        	fprintf(stdout, "memory allocation error for *corr\n");
        	exit(1);
	}
		build_dat_path(filename, "CORR_NAME_TYPE.DAT", sizeof filename, 0);
	if ((fp = fopen(filename, "r")) == NULL) {
        	fprintf(stdout, "Cannot open the CORR_NAME_TYPE.DAT file: %s, exit\n", filename);
        	exit(1);
	}
	for(;;) {
        	if (fgets(line, MAXCHAR, fp) == NULL)  break;
        	sscanf(line, "%s%s%d", corr[corrnum].type, corr[corrnum].element, &corr[corrnum].connum); 
		corrnum++;
		if(corrnum == MAXCORR) {
        		corr = (CORR *) realloc(corr, sizeof(CORR) * (MAXCORR + corrnum));
                	if (corr == NULL) {
                        	fprintf(stdout, "Error: memory reallocation with realloc for *corr\n");
                        	exit(1);
                	}
        	}
	}
	fclose(fp);

	if(debug)
		for(i=0;i<atomnum;i++)
			printf("\nATOM %d %s %8.3lf %8.3lf %8.3lf %8.3lf %s\n", i+1, atom[i].name, atom[i].x, atom[i].y, atom[i].z, atom[i].charge, atom[i].ambername);
	errorflag = 0;
	for(i=0;i<atomnum;i++) 
		for(j=0;j<corrnum;j++)
			if(strcmp(atom[i].ambername, corr[j].type) == 0) 
				if(corr[j].connum != -1  && corr[j].connum != atom[i].connum) { 
					fprintf(stdout, "[33mWarning: the number of bonded atoms (%d) for atom (ID: %d, Name: %s) does not match\n[0m", atom[i].connum, i+1, atom[i].name);
					fprintf(stdout, "[33m         the desired connectivity (%d) based on the atom type (%s) defined in CORR_NAME_TYPE.DAT.\n[0m", corr[j].connum, corr[j].type); 
					errorflag = 1;
/*		exit(1); */
/* Although error happens, it is ideal not to exit since many atom names (N2) same to AMBER atom types, such as N2 */
				}
	if(errorflag == 1) 
			fprintf(stdout, "[32mBUT, You may safely ingnor the warnings if your input molecule uses atom names or element as atom types\n[0m"); 
	free(corr);
}

/*part3: geometry */
printf("\n\n[34m-- Check Geometry --[0m\n");
printf("\n\n[34m-- For those bonded --[0m\n");

if ((fp = fopen(minfo.connect_file, "r")) == NULL) {
        fprintf(stdout, "Cannot open the minfo.connect_file: %s , exit\n", minfo.connect_file);
        exit(1);
}
for(i=0;i<120;i++)
	bondradius[i] = -1;
for (;;) {
        if (fgets(line, MAXCHAR, fp) == NULL)
                break;
        if (line[10] == '.') {
                sscanf(&line[3], "%d%lf", &tmpint, &radius);
		bondradius[tmpint] = radius;
	}
}

for(i=0;i<bondnum;i++) {
	bondrange(atom[bond[i].bondi].atomicnum, atom[bond[i].bondj].atomicnum);	
	if(bl.valid == 0) continue;
	dist = (atom[bond[i].bondi].x -  atom[bond[i].bondj].x)	* (atom[bond[i].bondi].x -  atom[bond[i].bondj].x);
	dist += (atom[bond[i].bondi].y -  atom[bond[i].bondj].y) * (atom[bond[i].bondi].y -  atom[bond[i].bondj].y);
	dist += (atom[bond[i].bondi].z -  atom[bond[i].bondj].z) * (atom[bond[i].bondi].z -  atom[bond[i].bondj].z);
	dist = sqrt(dist);
	if(dist < bl.low)
		fprintf(stdout, "[33m\nWarning: small distance for BOND\t%d\t%s\t%s\t%d\t%9.2lf  [%-5.2lf-%5.2lf][0m", i+1, atom[bond[i].bondi].name, atom[bond[i].bondj].name, bond[i].type, dist, bl.low, bl.high); 
	if(dist > bl.high)
		fprintf(stdout, "[33m\nWarning: large distance for BOND\t%d\t%s\t%s\t%d\t%9.2lf  [%-5.2lf-%5.2lf][0m", i+1, atom[bond[i].bondi].name, atom[bond[i].bondj].name, bond[i].type, dist,bl.low, bl.high); 
	if(debug)
		fprintf(stdout, "[33m\nBOND\t%d\t%s\t%s\t%d\t%9.2lf[0m", i+1, atom[bond[i].bondi].name, atom[bond[i].bondj].name, bond[i].type, dist); 
}

/* now check if those unbonded atoms come too close */
printf("\n\n[34m-- For those unbonded --[0m\n");
for(i=0;i<atomnum-1;i++)
	for(j=i+1;j<atomnum;j++) {
		bondflag = 0;	
		for(k=0;k<bondnum;k++) {
			if(bond[k].bondi == i && bond[k].bondj == j) {
				bondflag =1; 
				break; 
			}	
			if(bond[k].bondi == j && bond[k].bondj == i) {
				bondflag =1; 
				break; 
			}	
		}
		if(bondflag == 1) continue;
		bondrange(atom[i].atomicnum, atom[j].atomicnum);	
		if(bl.valid == 0) continue;
		dist =  (atom[i].x - atom[j].x) * (atom[i].x -  atom[j].x);
		dist += (atom[i].y - atom[j].y) * (atom[i].y -  atom[j].y);
		dist += (atom[i].z - atom[j].z) * (atom[i].z -  atom[j].z);
		dist = sqrt(dist);
		if(dist < (bl.bond + bl.offset * 1.25)) {
			fprintf(stdout, "[33m\nWarning: although no bond formed between ATOM\t%d\t%s and ATOM\t%d\t%s, the distance\t%7.2lf is very small[0m", i+1, atom[i].name, j+1, atom[j].name, dist); 

		}
	}
printf("\n\n[32m-- Status: pass --[0m\n");

/*part4  now check bonds */
printf("\n\n[34m-- check weird bonds --[0m\n");
if(ibondtype == 1) {
	if(debug)
		for(i=0;i<bondnum;i++)
			printf("\nBOND %d  %d  %s  %d  %s  %d\n", i+1, bond[i].bondi+1, atom[bond[i].bondi].name,  bond[i].bondj+1, atom[bond[i].bondj].name, bond[i].type);
	for(i=0;i<atomnum;i++) {  
		tvalence = 0;
		val_ar = 1;
		for(j=0;j<bondnum;j++) 
			if(bond[j].bondi == i || bond[j].bondj == i) {
				if(bond[j].type <= 3) 
					tvalence += bond[j].type;
				else if(bond[j].type == 9 || bond[j].type == 10) {
					tvalence += val_ar;
					if(val_ar == 1) 
						val_ar = 2;
					else 
						val_ar = 1;
				}
			}
		if(atom[i].atomicnum == 6) {
			if(tvalence > 4) {
				fprintf(stdout, "[31mError: weird atomic valence (%d) for atom (ID: %d, Name: %s), please check atomic connectivity\n[0m", tvalence, i+1, atom[i].name);	
				exit(1);
			}
			if(tvalence < 4) {
				fprintf(stdout, "[31mError: weird atomic valence (%d) for atom (ID: %d, Name: %s), possible have open valence\n[0m", tvalence, i+1, atom[i].name);	
				exit(1);
			}
		}

		if(atom[i].atomicnum == 1 || atom[i].atomicnum == 9 || atom[i].atomicnum == 17 || atom[i].atomicnum == 35 || atom[i].atomicnum == 53)  
			if(tvalence !=1) {
				fprintf(stdout, "[31mError: weird atomic valence (%d) for atom (ID: %d, Name: %s), please check atomic connectivity\n[0m", tvalence, i+1, atom[i].name);	
				exit(1);
			}

		if(atom[i].atomicnum == 7) {
			if(tvalence > 4) {
				fprintf(stdout, "[31mError: weird atomic valence (%d) for atom (ID: %d, Name: %s), please check atomic connectivity\n[0m", tvalence, i+1, atom[i].name);	
				exit(1);
			}
			if(tvalence < 3) {
				fprintf(stdout, "[31mError: weird atomic valence (%d) for atom (ID: %d, Name: %s), possible have open valence\n[0m", tvalence, i+1, atom[i].name);	
				exit(1);
			}
		}
		if(atom[i].atomicnum == 6) {
			if(tvalence > 4) {
				fprintf(stdout, "[31mError: weird atomic valence (%d) for atom (ID: %d, Name: %s), please check atomic connectivity\n[0m", tvalence, i+1, atom[i].name);	
				exit(1);
			}
			if(tvalence <= 3) {
				fprintf(stdout, "[31mError: weird atomic valence (%d) for atom (ID: %d, Name: %s), possible have open valence\n[0m", tvalence, i+1, atom[i].name);	
				exit(1);
			}
		}
		if(atom[i].atomicnum == 8) {
			if(tvalence > 2) {
				fprintf(stdout, "[31mError: weird atomic valence (%d) for atom (ID: %d, Name: %s), please check atomic connectivity\n[0m", tvalence, i+1, atom[i].name);	
				exit(1);
			}
			if(tvalence <  1) {
				fprintf(stdout, "[31mError: weird atomic valence (%d) for atom (ID: %d, Name: %s), possible have open valence\n[0m", tvalence, i+1, atom[i].name);	
				exit(1);
			}
		}
	}
}
/* check if bonded atoms exceed 6 */
if(ibondtype == 0) 
	for(i=0;i<atomnum;i++) {  
		connum = 0;
		for(j=0;j<bondnum;j++)  {
			if(bond[j].bondi == i || bond[j].bondj == i)
				connum ++;
			if(connum >6) break;
		}
		if(connum > 6) {
			fprintf(stdout, "[31mError: Number of bonded atoms (%d) for atom (ID: %d, Name: %s) exceed 6, antechamber cannot handle such kinds of molecules, exit\n[0m", connum, i+1, atom[i].name);	
			exit(1);
		}
	}

printf("\n\n[32m-- Status: pass --[0m\n");

/*part5: check if all atoms are linked together through paths*/
printf("\n\n[34m-- Check if the input is only one unit --[0m\n");
for(i=0;i<atomnum;i++)   
	atom[i].select = 0;
atom[0].select = 1;
selectflag = 1;
while(selectflag) {
	selectflag = 0;
	for(i=0;i<atomnum;i++) {  
		if(atom[i].select == 0) continue;
		for(j=0;j<atom[i].connum;j++)  
			if(atom[atom[i].con[j]].select == 0) {
				atom[atom[i].con[j]].select = 1;
				selectflag = 1;
			}
	}
}
numselect = 0;
for(i=0;i<atomnum;i++)   
	numselect+= atom[i].select;
if(numselect < atomnum) {
	fprintf(stdout, "[31mError: There are more than one unit, antechamber can only handle one unit molecule, exit\n[0m");	
	fprintf(stdout, "[31m       If you think your input is a single molecule, the geometry may be bad and the connectivity must be wrong.\n[0m");	
	fprintf(stdout, "[31m       Please convert your molecule to a mol2 file with antechamber using \"-j 5 -at sybyl\" flags\n[0m");	
	fprintf(stdout, "[31m       And then check your molecue with visual program and manually add new bonds or delete unwanted bonds\n[0m");	
	exit(0);
}

printf("\n\n[32m-- Status: pass --[0m\n");

printf("\n");

}


void usage() {
	if (strcmp(COLORTEXT, "YES") == 0 || strcmp(COLORTEXT, "yes") == 0) {
		printf("[31mUsage: acdoctor -i [0m input file name\n"
			   "[31m                -f [0m input file format\n");
		printf  ("\n               [31m List of the File Formats [0m \n");
		printf ("\n	 	file format type  abbre. index | file format type abbre. index");
		printf ("\n		--------------------------------------------------------------- ");
		printf ("\n		Antechamber        ac       1  | Sybyl Mol2         mol2    2 ");
		printf ("\n		PDB                pdb      3  | Modified PDB       mpdb    4 ");
		printf ("\n		AMBER PREP (int)   prepi    5  | AMBER PREP (car)   prepc   6 ");
		printf ("\n		Gaussian Z-Matrix  gzmat    7  | Gaussian Cartesian gcrt    8 ");
		printf ("\n		Mopac Internal     mopint   9  | Mopac Cartesian    mopcrt 10 ");
		printf ("\n		Gaussian Output    gout    11  | Mopac Output       mopout 12 ");
		printf ("\n		Alchemy            alc     13  | CSD                csd    14 ");
		printf ("\n		MDL                mdl     15  | Hyper              hin    16 ");
		printf ("\n		AMBER Restart      rst     17  | Jaguar Cartesian   jcrt   18 ");
		printf ("\n		Jaguar Z-Matrix    jzmat   19  | Jaguar Output      jout   20 ");
		printf ("\n		Divcon Input       divcrt  21  | Divcon Output      divout 22 ");
		printf ("\n		SQM Input          sqmcrt  23  | SQM Output         sqmout 24 ");
		printf ("\n		Charmm             charmm  25  | Gaussian ESP       gesp   26 ");
		printf ("\n		--------------------------------------------------------------\n");
	} else {
                printf("Usage: acdoctor -i input file name\n"
                       "                -f input file format\n");
                printf ("\n                List of the File Formats \n"); 
		printf ("\n            file format type  abbre. index | file format type abbre. index");
                printf ("\n            --------------------------------------------------------------- ");
                printf ("\n            Antechamber        ac       1  | Sybyl Mol2         mol2    2 ");
                printf ("\n            PDB                pdb      3  | Modified PDB       mpdb    4 ");
                printf ("\n            AMBER PREP (int)   prepi    5  | AMBER PREP (car)   prepc   6 ");
                printf ("\n            Gaussian Z-Matrix  gzmat    7  | Gaussian Cartesian gcrt    8 ");
                printf ("\n            Mopac Internal     mopint   9  | Mopac Cartesian    mopcrt 10 ");
                printf ("\n            Gaussian Output    gout    11  | Mopac Output       mopout 12 ");
                printf ("\n            Alchemy            alc     13  | CSD                csd    14 ");
                printf ("\n            MDL                mdl     15  | Hyper              hin    16 ");
                printf ("\n            AMBER Restart      rst     17  | Jaguar Cartesian   jcrt   18 ");
                printf ("\n            Jaguar Z-Matrix    jzmat   19  | Jaguar Output      jout   20 ");
                printf ("\n            Divcon Input       divcrt  21  | Divcon Output      divout 22 ");
		printf ("\n	       SQM Input          sqmcrt  23  | SQM Output         sqmout 24 ");
		printf ("\n	       Charmm             charmm  25  | Gaussian ESP       gesp   26 ");
                printf ("\n            --------------------------------------------------------------\n");
	}
}

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
/*flag = 1  <->atom
       = 2  <->bond
       = 3  <->arom
       = 4  <->atom + bond
       = 5  <->atom + arom 
       = 6  <->bond + arom
       = 7  <->atom + arom +bond
*/
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

void judgebondtype(int atomnum, ATOM * atom, int bondnum, BOND * bond,
				   CONTROLINFO cinfo, MOLINFO minfo, int bondtype_flag)
{
	char tmpchar[MAXCHAR];
	int status = 0;
	size_t copied_size;
	char *options;
	wac("ACDOCTOR_INIT.ac", atomnum, atom, bondnum, bond, cinfo,
		minfo);

        copied_size = build_exe_path(tmpchar, "bondtype", sizeof tmpchar, 1);
        if (bondtype_flag == 1)
                options = " -j part -i ACDOCTOR_INIT.ac"
                          " -o ACDOCTOR_BOND.ac -f ac" ;
        else
                options = " -j full -i ACDOCTOR_INIT.ac"
                          " -o ACDOCTOR_BOND.ac -f ac" ;
        strncat(tmpchar, options, MAXCHAR - copied_size );

	fprintf(stdout, "Running: %s\n", tmpchar);
	status = system(tmpchar);
        if(status != 0) {
                fprintf(stdout, "Error: cannot run \"%s\" in judgebondtype() of antechamber.c properly, exit\n", tmpchar);
                exit(1);
        }
	minfo2 = minfo;
	rac("ACDOCTOR_BOND.ac", &atomnum, atom, &bondnum, bond, &cinfo, &minfo2);
}

void preparation(void) {
        if (atomicnum_flag)
                atomicnum(atomnum, atom);
        if (atomname_flag)
                atomname(atomnum, atom);
        if (adjustatomname_flag) {
                if(strcmp(cinfo.intype, "mol2")==0 || strcmp(cinfo.intype, "2")==0 ||
                   strcmp(cinfo.intype, "ac")==0 || strcmp(cinfo.intype, "1")==0)
                        adjustatomname(atomnum, atom, 1);
                else
                        adjustatomname(atomnum, atom, 0);
        }
        if (default_flag)
                default_inf(atomnum, atom, default_flag);
        if (cartcoord_flag)
                cartcoord(atomnum, atom);
        if (connect_flag) {
                overflow_flag =
                        connect(minfo.connect_file, atomnum, atom, &bondnum, bond,
                                        cinfo.maxbond);
                if (overflow_flag) {
                        cinfo.maxbond = bondnum + 10;
                        memory(2, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
                        overflow_flag =
                                connect(minfo.connect_file, atomnum, atom, &bondnum, bond,
                                                cinfo.maxbond);
                }
        }
        if (duplicatedname_flag)
                duplicatedname(atomnum, atom);
}

void try(void) {
char tmpchar[MAXCHAR];
int status;
/*      bond type prediction */
fprintf(stdout, "[34m -- Now try to judge bond type -- [0m\n");
if (bondnum > 0) 
	judgebondtype(atomnum, atom, bondnum, bond, cinfo, minfo, bondtype_flag);

fprintf(stdout, "[34m -- Now try to assign atom type -- [0m\n");
/*      atom type assignment    */

build_exe_path(tmpchar, "atomtype", sizeof tmpchar, 1);
strcat(tmpchar, " -i ACDOCTOR_BOND.ac -o ACDOCTOR_ATOM.ac -p ");
strcat(tmpchar, minfo.atom_type_def);
fprintf(stdout, "\nRunning: %s\n", tmpchar);
status = system(tmpchar);
if(status != 0) {
        fprintf(stdout, "Error: cannot run \"%s\" in main() of antechamber.c properly, exit\n", tmpchar);
        exit(1);
}

fprintf(stdout, "[34m -- Now write out ACDOCTOR.mol2 -- [0m\n");

build_exe_path(tmpchar, "atomtype", sizeof tmpchar, 1);
strcat(tmpchar, " -i ACDOCTOR_BOND.ac -o ACDOCTOR_SYBYL.ac -p sybyl");
status = system(tmpchar);
if(status != 0) {
        fprintf(stdout, "Error: cannot run \"%s\" in main() of antechamber.c properly, exit\n", tmpchar);
        exit(1);
}
minfo2 = minfo;
rac("ACDOCTOR_SYBYL.ac", &atomnum, atom, &bondnum, bond, &cinfo, &minfo2);
wmol2("ACDOCTOR.mol2", atomnum, atom, bondnum, bond, arom, cinfo,minfo);
fprintf(stdout, "\n[32m -- IF no error occurs in bond and atom type assignments,[0m");
fprintf(stdout, "\n[32m    %s should have no problem for the antechamber package. -- \n[0m", ifilename);
fprintf(stdout, "\nACDOCTOR_INIT.ac  : ac file before bond type and atom type assignment"); 
fprintf(stdout, "\nACDOCTOR_BT.ac    : ac file after bond type assignment"); 
fprintf(stdout, "\nACDOCTOR_AT.ac    : ac file after gaff atom type assignment"); 
fprintf(stdout, "\nACDOCTOR_SYBYL.ac : ac file after sybyl atom type assignment"); 
fprintf(stdout, "\nACDOCTOR.mol2     : mol2 file after sybyl atom type assignment\n"); 

}

int main(int argc, char *argv[])
{
	int i,j,k;
	int index;

    amberhome = (char *) getenv("AMBERHOME");
    if( amberhome == NULL ){
       fprintf( stdout, "AMBERHOME is not set!\n" );
       exit(1);
    }
	if (argc == 2)
		if (strncmp(argv[1], "-h", 2) == 0
			|| strncmp(argv[1], "-H", 2) == 0) {
			usage();
			exit(1);
		}
	if (argc == 1) {
		usage();
		exit(1);
	}

/* 	set defaults information */
	default_cinfo(&cinfo);
	default_minfo(&minfo);
        atomtype_flag = 0;
        bondtype_flag = 0;
/* 	
	0 <-> not set 
	1 <-> set atom[].chain, atom[].ter and atom[].ambername	
	2 <-> set atom[].chain, atom[].ter
*/
        default_flag = 0; 
        atomname_flag = 0;
        atomicnum_flag = 0;
        adjustatomname_flag = 0;
        duplicatedname_flag = 1;
        cartcoord_flag = 0;
        connect_flag = 0;


	for (i = 1; i < argc - 1; i += 2) {
		if (strcmp(argv[i], "-i") == 0) {
			strcpy(ifilename, argv[i + 1]);
			continue;
		} else if (strcmp(argv[i], "-f") == 0) {
			strcpy(cinfo.intype, argv[i + 1]);
			continue;
		} 
	}

/* 	for connect.tpl and radius parameter files */
	build_dat_path(minfo.connect_file, "CONNECT.TPL", sizeof minfo.connect_file, 0);
	build_dat_path(minfo.radius_file, "RADIUS.DAT", sizeof minfo.radius_file, 0);


/*      allocate memory using default parameters MAXATOM and MAXBOND */
	memory(0, MAXATOM, MAXBOND, MAXRING);

/******************************************/
/* 	The following codes readin input file */
/******************************************/

	if (strcmp("ac", cinfo.intype) == 0 || strcmp("1", cinfo.intype) == 0) {
		formatchk(ifilename, "ac");
		overflow_flag =
			rac(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo);
		if (overflow_flag) {
			cinfo.maxatom = atomnum + 10;
			cinfo.maxbond = bondnum + 10;
			memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
			overflow_flag =
				rac(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo);
		}
		adjustatomname_flag = 1;
		atomicnum_flag = 1;
		iatomtype = 1;
		ibondtype = 1;
	}

	if (strcmp("mol2", cinfo.intype) == 0
		|| strcmp("2", cinfo.intype) == 0) {

		formatchk(ifilename, "mol2");
		overflow_flag =
			rmol2(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo, 1);
		if (overflow_flag) {
			cinfo.maxatom = atomnum + 10;
			cinfo.maxbond = bondnum + 10;
			memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
			overflow_flag =
				rmol2(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo, 1);
		}
		default_flag = 2;
		atomicnum_flag = 1;
		adjustatomname_flag = 1;
		iatomtype = 1;
		ibondtype = 1;
	}
	if (strcmp("mopint", cinfo.intype) == 0
		|| strcmp("9", cinfo.intype) == 0) {
		formatchk(ifilename, "mopint");
		overflow_flag = rmopint(ifilename, &atomnum, atom, cinfo, minfo);
		if (overflow_flag) {
			cinfo.maxatom = atomnum + 10;
			cinfo.maxbond = bondnum + 10;
			memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
			overflow_flag =
				rmopint(ifilename, &atomnum, atom, cinfo, minfo);
		}
		cartcoord_flag = 1;
		atomicnum_flag = 1;
		atomname_flag = 1;
		default_flag = 1;
		connect_flag = 1;
		bondtype_flag = 2;
	}

	if (strcmp("mopcrt", cinfo.intype) == 0
		|| strcmp("10", cinfo.intype) == 0) {
		formatchk(ifilename, "mopcrt");
		overflow_flag = rmopcrt(ifilename, &atomnum, atom, cinfo, minfo);
		if (overflow_flag) {
			cinfo.maxatom = atomnum + 10;
			cinfo.maxbond = bondnum + 10;
			memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
			overflow_flag =
				rmopcrt(ifilename, &atomnum, atom, cinfo, minfo);
		}
		atomicnum_flag = 1;
		atomname_flag = 1;
		default_flag = 1;
		connect_flag = 1;
		bondtype_flag = 2;
	}

	if (strcmp("mopout", cinfo.intype) == 0
		|| strcmp("12", cinfo.intype) == 0) {
		formatchk(ifilename, "mopout");
		overflow_flag = rmopout(ifilename, &atomnum, atom, &cinfo, &minfo);
		if (overflow_flag) {
			cinfo.maxatom = atomnum + 10;
			cinfo.maxbond = bondnum + 10;
			memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
			overflow_flag =
				rmopout(ifilename, &atomnum, atom, &cinfo, &minfo);
		}
		atomicnum_flag = 1;
		atomname_flag = 1;
		default_flag = 1;
		connect_flag = 1;
		bondtype_flag = 2;
	}

	if (strcmp("gcrt", cinfo.intype) == 0
		|| strcmp("8", cinfo.intype) == 0) {
		formatchk(ifilename, "gcrt");
		overflow_flag = rgcrt(ifilename, &atomnum, atom, cinfo, minfo);
		if (overflow_flag) {
			cinfo.maxatom = atomnum + 10;
			cinfo.maxbond = bondnum + 10;
			memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
			overflow_flag = rgcrt(ifilename, &atomnum, atom, cinfo, minfo);
		}
		atomname_flag = 1;
		default_flag = 1;
		connect_flag = 1;
		bondtype_flag = 2;
	}
	if (strcmp("gzmat", cinfo.intype) == 0
		|| strcmp("7", cinfo.intype) == 0) {
		formatchk(ifilename, "gzmat");
		overflow_flag = rgzmat(ifilename, &atomnum, atom, cinfo, minfo);
		if (overflow_flag) {
			cinfo.maxatom = atomnum + 10;
			cinfo.maxbond = bondnum + 10;
			memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
			overflow_flag =
				rgzmat(ifilename, &atomnum, atom, cinfo, minfo);
		}
		cartcoord_flag = 1;
		atomicnum_flag = 1;
		atomname_flag = 1;
		default_flag = 1;
		connect_flag = 1;
		bondtype_flag = 2;
	}

	if (strcmp("gout", cinfo.intype) == 0
		|| strcmp("11", cinfo.intype) == 0) {
		formatchk(ifilename, "gout");
		overflow_flag = rgout(ifilename, &atomnum, atom, cinfo, &minfo);
		if (overflow_flag) {
			cinfo.maxatom = atomnum + 10;
			cinfo.maxbond = bondnum + 10;
			memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
			overflow_flag = rgout(ifilename, &atomnum, atom, cinfo, &minfo);
		}
		atomname_flag = 1;
		default_flag = 1;
		connect_flag = 1;
		bondtype_flag = 2;
	}

	if (strcmp("jcrt", cinfo.intype) == 0
		|| strcmp("18", cinfo.intype) == 0) {
		formatchk(ifilename, "jcrt");
		overflow_flag = rjcrt(ifilename, &atomnum, atom, cinfo, minfo);
		if (overflow_flag) {
			cinfo.maxatom = atomnum + 10;
			cinfo.maxbond = bondnum + 10;
			memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
			overflow_flag = rjcrt(ifilename, &atomnum, atom, cinfo, minfo);
		}
		default_flag = 1;
		connect_flag = 1;
		bondtype_flag = 2;
	}

        if (strcmp("jzmat", cinfo.intype) == 0
                || strcmp("19", cinfo.intype) == 0) {
		formatchk(ifilename, "jzmat");
                overflow_flag = rjzmat(ifilename, &atomnum, atom, cinfo, minfo);
                if (overflow_flag) {
                        cinfo.maxatom = atomnum + 10;
                        cinfo.maxbond = bondnum + 10;
                        memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
                        overflow_flag =
                                rjzmat(ifilename, &atomnum, atom, cinfo, minfo);
                }
                cartcoord_flag = 1;
                atomicnum_flag = 1;
                atomname_flag = 0;
                default_flag = 1;
                connect_flag = 1;
                bondtype_flag = 2;
        }

        if (strcmp("jout", cinfo.intype) == 0
                || strcmp("20", cinfo.intype) == 0) {
		formatchk(ifilename, "jout");
                overflow_flag = rjout(ifilename, &atomnum, atom, cinfo, &minfo);
                if (overflow_flag) {
                        cinfo.maxatom = atomnum + 10;
                        cinfo.maxbond = bondnum + 10;
                        memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
                        overflow_flag = rjout(ifilename, &atomnum, atom, cinfo, &minfo);
                }
                default_flag = 1;
                connect_flag = 1;
                bondtype_flag = 2;
        }

	if (strcmp("pdb", cinfo.intype) == 0 || strcmp("3", cinfo.intype) == 0) {
		overflow_flag = rpdb(ifilename, &atomnum, atom, cinfo, minfo, 0);
		formatchk(ifilename, "pdb");
		if (overflow_flag) {
			cinfo.maxatom = atomnum + 10;
			cinfo.maxbond = bondnum + 10;
			memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
			overflow_flag =
				rpdb(ifilename, &atomnum, atom, cinfo, minfo, 0);
		}
		adjustatomname_flag = 1;
		atomicnum_flag = 1;
		default_flag = 1;
		bondtype_flag = 2;
		index = 0;
		for(i=0; i<atomnum; i++)
			if(atom[i].connum > 0) {
				index = 1;
				break;
			}
		if(index == 1) {
			bondnum = 0;
			for(i=0; i<atomnum-1; i++)
				for(j=i+1; j<atomnum; j++) 
					for(k=0; k<6; k++)  {
						if(atom[i].con[k] == -1) break;
						if(atom[i].con[k] == j) {
							bond[bondnum].bondi = i;
							bond[bondnum].bondj = j;
							bondnum++;
		  				}
					}
		}
		else 
			connect_flag = 1;
	}

	if (strcmp("mpdb", cinfo.intype) == 0
		|| strcmp("4", cinfo.intype) == 0) {
		formatchk(ifilename, "mpdb");
		overflow_flag = rpdb(ifilename, &atomnum, atom, cinfo, minfo, 1);
		if (overflow_flag) {
			cinfo.maxatom = atomnum + 10;
			cinfo.maxbond = bondnum + 10;
			memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
			overflow_flag =
				rpdb(ifilename, &atomnum, atom, cinfo, minfo, 1);
		}
		adjustatomname_flag = 1;
		atomicnum_flag = 1;
		default_flag = 0;
		bondtype_flag = 2;
		index = 0;
		for(i=0; i<atomnum; i++)
			if(atom[i].connum > 0) {
				index = 1;
				break;
			}
                if(index == 1) {
                        bondnum = 0;
                        for(i=0; i<atomnum-1; i++)
                                for(j=i+1; j<atomnum; j++) 
                                        for(k=0; k<6; k++)  {       
                                                if(atom[i].con[k] == -1) break;
                                                if(atom[i].con[k] == j) {
                                                        bond[bondnum].bondi = i;
                                                        bond[bondnum].bondj = j;
                                                        bondnum++;
                                                }
					}
                }
		else 
			connect_flag = 1;
	}
	if (strcmp("csd", cinfo.intype) == 0
		|| strcmp("14", cinfo.intype) == 0) {
		formatchk(ifilename, "csd");
		overflow_flag = rcsd(ifilename, &atomnum, atom, cinfo, minfo);
		if (overflow_flag) {
			cinfo.maxatom = atomnum + 10;
			cinfo.maxbond = bondnum + 10;
			memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
			overflow_flag = rcsd(ifilename, &atomnum, atom, cinfo, minfo);
		}
		atomicnum_flag = 1;
		default_flag = 1;
		connect_flag = 1;
		bondtype_flag = 2;
	}

	if (strcmp("mdl", cinfo.intype) == 0
		|| strcmp("15", cinfo.intype) == 0) {
		formatchk(ifilename, "mdl");
		overflow_flag =
			rmdl(ifilename, &atomnum, atom, &bondnum, bond, cinfo, minfo);
		if (overflow_flag) {
			cinfo.maxatom = atomnum + 10;
			cinfo.maxbond = bondnum + 10;
			memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
			overflow_flag =
				rmdl(ifilename, &atomnum, atom, &bondnum, bond, cinfo,
					 minfo);
		}
		atomicnum_flag = 1;
		atomname_flag = 1;
		default_flag = 1;
		bondtype_flag = 1;
		ibondtype = 1;
	}

	if (strcmp("alc", cinfo.intype) == 0
		|| strcmp("13", cinfo.intype) == 0) {
		formatchk(ifilename, "alc");
		overflow_flag =
			ralc(ifilename, &atomnum, atom, &bondnum, bond, cinfo, minfo);
		if (overflow_flag) {
			cinfo.maxatom = atomnum + 10;
			cinfo.maxbond = bondnum + 10;
			memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
			overflow_flag =
				ralc(ifilename, &atomnum, atom, &bondnum, bond, cinfo,
					 minfo);
		}
		atomicnum_flag = 1;
		atomname_flag = 1;
		bondtype_flag = 1;
		default_flag = 1;
		ibondtype = 1;
	}

	if (strcmp("hin", cinfo.intype) == 0
		|| strcmp("16", cinfo.intype) == 0) {
		formatchk(ifilename, "hin");
		overflow_flag =
			rhin(ifilename, &atomnum, atom, &bondnum, bond, cinfo, minfo);
		if (overflow_flag) {
			cinfo.maxatom = atomnum + 10;
			cinfo.maxbond = bondnum + 10;
			memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
			overflow_flag =
				rhin(ifilename, &atomnum, atom, &bondnum, bond, cinfo,
					 minfo);
		}
		atomicnum_flag = 1;
		atomname_flag = 1;
		bondtype_flag = 1;
		default_flag = 1;
		ibondtype = 1;
	}

	if (strcmp("prepi", cinfo.intype) == 0
		|| strcmp("5", cinfo.intype) == 0) {
		formatchk(ifilename, "prepi");
		overflow_flag = rprepi(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo);
		if (overflow_flag) {
			cinfo.maxatom = atomnum + 10;
			cinfo.maxbond = bondnum + 10;
			memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
			overflow_flag =
				rprepi(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo);
		}
		atomicnum_flag = 1;
		bondtype_flag = 2;
		connect_flag = 0;
		iatomtype = 1;
	}

	if (strcmp("prepc", cinfo.intype) == 0
		|| strcmp("5", cinfo.intype) == 0) {
		formatchk(ifilename, "prepc");
		overflow_flag = rprepc(ifilename, &atomnum, atom, &cinfo, &minfo);
		if (overflow_flag) {
			cinfo.maxatom = atomnum + 10;
			cinfo.maxbond = bondnum + 10;
			memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
			overflow_flag =
				rprepc(ifilename, &atomnum, atom, &cinfo, &minfo);
		}
		atomicnum_flag = 1;
		bondtype_flag = 2;
		connect_flag = 1;
		iatomtype = 1;
	}

	if (strcmp("rst", cinfo.intype) == 0
		|| strcmp("17", cinfo.intype) == 0) {
		fprintf(stdout,
				"RST (17) file format can only be additional file because it only has coordinate information\n");
		exit(1);
	}

        if (strcmp("divcrt", cinfo.intype) == 0
                || strcmp("21", cinfo.intype) == 0) {
		formatchk(ifilename, "divcrt");
                overflow_flag = rdivcrt(ifilename, &atomnum, atom, cinfo, minfo);
                if (overflow_flag) {
                        cinfo.maxatom = atomnum + 10;
                        cinfo.maxbond = bondnum + 10;
                        memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
                        overflow_flag =
                                rdivcrt(ifilename, &atomnum, atom, cinfo, minfo);
                }
                atomicnum_flag = 1;
                atomname_flag = 1;
                default_flag = 1;
                connect_flag = 1;
                bondtype_flag = 2;
        }
                                                                                                                                                                                                           
        if (strcmp("divout", cinfo.intype) == 0
                || strcmp("22", cinfo.intype) == 0) {
		formatchk(ifilename, "divout");
                overflow_flag = rdivout(ifilename, &atomnum, atom, &cinfo, &minfo);
                if (overflow_flag) {
                        cinfo.maxatom = atomnum + 10;
                        cinfo.maxbond = bondnum + 10;
                        memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
                        overflow_flag =
                                rdivout(ifilename, &atomnum, atom, &cinfo, &minfo);
                }
                atomicnum_flag = 1;
                atomname_flag = 1;
                default_flag = 1;
                connect_flag = 1;
                bondtype_flag = 2;
        }

        if (strcmp("sqmcrt", cinfo.intype) == 0
                || strcmp("23", cinfo.intype) == 0) {
                formatchk(ifilename, "sqmcrt");
                overflow_flag = rsqmcrt(ifilename, &atomnum, atom, cinfo, minfo);
                if (overflow_flag) {
                        cinfo.maxatom = atomnum + 10;
                        cinfo.maxbond = bondnum + 10;
                        memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
                        overflow_flag =
                                rsqmcrt(ifilename, &atomnum, atom, cinfo, minfo);
                }               
                atomicnum_flag = 1;
                atomname_flag = 1;
                default_flag = 1;
                connect_flag = 1;
                bondtype_flag = 2;
        }       
         
        if (strcmp("sqmout", cinfo.intype) == 0
                || strcmp("24", cinfo.intype) == 0) {
                formatchk(ifilename, "sqmout");
                overflow_flag = rsqmout(ifilename, &atomnum, atom, &cinfo, &minfo);
                if (overflow_flag) {
                        cinfo.maxatom = atomnum + 10;
                        cinfo.maxbond = bondnum + 10;
                        memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
                        overflow_flag =
                                rsqmout(ifilename, &atomnum, atom, &cinfo, &minfo);
                }
                atomicnum_flag = 1;
                atomname_flag = 1;
                default_flag = 1;
                connect_flag = 1;
                bondtype_flag = 2;
        }
        if (strcmp("charmm", cinfo.intype) == 0
                || strcmp("25", cinfo.intype) == 0) {

                overflow_flag =
                        rcharmm(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo);
                if (overflow_flag) {
                        cinfo.maxatom = atomnum + 10;
                        cinfo.maxbond = bondnum + 10;
                        memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
                        overflow_flag =
                                rcharmm(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo);
                }
                default_flag = 1;
                atomicnum_flag = 1;
                adjustatomname_flag = 1;
        }
        if (strcmp("gesp", cinfo.intype) == 0
                || strcmp("26", cinfo.intype) == 0) {
                overflow_flag = rgesp(ifilename, &atomnum, atom, cinfo, minfo);
                if (overflow_flag) {
                        cinfo.maxatom = atomnum + 10;
                        cinfo.maxbond = bondnum + 10;
                        memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
                        overflow_flag = rgesp(ifilename, &atomnum, atom, cinfo, minfo);
                }
                atomicnum_flag = 1;
                atomname_flag = 1;
                default_flag = 1;
                connect_flag = 1;
                bondtype_flag = 2;
        }

	if (adjustatomname_flag) {
		if(strcmp(cinfo.intype, "mol2")==0 || strcmp(cinfo.intype, "2")==0 ||
		   strcmp(cinfo.intype, "ac")==0 || strcmp(cinfo.intype, "1")==0)
			adjustatomname(atomnum, atom, 1);
		else
			adjustatomname(atomnum, atom, 0);
	}
	if (atomicnum_flag) 
		atomicnum(atomnum, atom);
	if (atomname_flag)
		atomname(atomnum, atom);
	if (default_flag)
		default_inf(atomnum, atom, default_flag);
	if (cartcoord_flag)
		cartcoord(atomnum, atom);
	if (connect_flag) {
		overflow_flag =
			connect(minfo.connect_file, atomnum, atom, &bondnum, bond,
					cinfo.maxbond);
		if (overflow_flag) {
			cinfo.maxbond = bondnum + 10;
			memory(2, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
			overflow_flag =
				connect(minfo.connect_file, atomnum, atom, &bondnum, bond,
						cinfo.maxbond);
		}
	}
/* 	print out information */
	if(debug) {
		for(i=0;i<atomnum;i++)
			printf("ATOM\t%d\t%s\t%d\t%s\t%8.3lf\t%8.3lf\t%8.3lf\t%8.3lf\t%d\t%d\t%d\t%d\t%d\t%d\n",
				i+1, atom[i].name, atom[i].resno, atom[i].aa, atom[i].x, atom[i].y, atom[i].z, atom[i].charge, atom[i].con[0]+1, atom[i].con[1]+1, atom[i].con[2]+1, atom[i].con[3]+1, atom[i].con[4]+1, atom[i].con[5]+1);
	for(i=0;i<bondnum;i++)
			printf("BOND\t%d\t%d\t%s\t%d\t%s\n", i+1, bond[i].bondi+1, atom[bond[i].bondi].name, bond[i].bondj+1, atom[bond[i].bondj].name); 
	}

	preparation();
	molchk();
	try();

	printf("\n");
	return (0);
}
