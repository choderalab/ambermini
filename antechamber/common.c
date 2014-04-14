#include <math.h>
# define NUMBER_OF_CHEMICAL_ELEMENTS 118

/*
   Build the full path of an executable/lib/data file.
   The executable/lib/file name is in the argument command.
   The argument sizeof_path is the total size of the memory
   pointed to by path, ie including space for the null terminator.
   Return the size of the path it tried to build; thus, overflow
   is indicated if the return value is greater than the sizeof_path;
   However, on overflow an error is emitted and the program exits.
   Assumptions:  strncat correctly handles a non-positive 3rd argument.
   'for_system' indicates if the path will be used in a 'system'
   command and therefore needs quotes if there are spaces in the
   AMBERHOME variable.
 */
size_t build_path(char *path, const char *subdir, const char *fname,
		size_t sizeof_path, int for_system)
{
    /* number of chars copied plus 1; on overflow > sizeof_path */
	char *quote = NULL;
    size_t c = 1;
    path[0] = '\0';
	if (for_system && index(amberhome, ' ') != NULL) {
		if (index(amberhome, '"') == NULL)
			quote = "\"";
		else
			quote = "'";
		strncat(path, quote, sizeof_path - c);
		c += 1;
	}
    strncat( path, amberhome, sizeof_path - c );
    c += strlen( amberhome );
    strncat( path, subdir, sizeof_path - c );
    c += strlen( subdir );
    strncat( path, fname, sizeof_path - c );
    c += strlen( fname );
	if (quote != NULL) {
		strncat(path, quote, sizeof_path - c);
		c += 1;
	}
    if (c > sizeof_path) {
        fprintf( stdout, "Error: insufficient string size !\n" );
        fprintf( stdout, "  Increase MAXCHAR in define.h and rebuild.\n" );
        fprintf( stdout, "  Truncated string: %s\n", path );
        exit( 1 );
    }
    return c;
}
size_t build_exe_path( char *path, const char *command, size_t sizeof_path, int for_system )
{
	return build_path(path, "/bin/", command, sizeof_path, for_system);
}
size_t build_dat_path( char *path, const char *command, size_t sizeof_path, int for_system )
{
	return build_path(path, "/dat/antechamber/", command, sizeof_path, for_system);
}

void default_minfo(MOLINFO * minfo) {
	strcpy((*minfo).dkeyword, "CARTESIAN AM1 STANDARD DIRECT OPT=BFGS XTEST=0.0001");
	strcpy((*minfo).mkeyword, "AM1 ANALYT MMOK GEO-OK PRECISE");
	strcpy((*minfo).skeyword, "  qm_theory='AM1', grms_tol=0.0002,\n  tight_p_conv=1, scfconv=1.d-10, ");
	strcpy((*minfo).gkeyword, "#HF/6-31G* SCF=tight Test Pop=MK iop(6/33=2) iop(6/42=6) opt");
	strcpy((*minfo).gesp,   "g09.gesp");
	strcpy((*minfo).ekeyword, "");
	strcpy((*minfo).gm, "");
	strcpy((*minfo).gn, "");
	strcpy((*minfo).resname, "MOL");
	strcpy((*minfo).atom_type_def, "gaff");
	strcpy((*minfo).resfilename, "molecule.res");
	strcpy((*minfo).chkfile, "molecule");
	strcpy((*minfo).gfilename, "GASPARM.DAT");
	strcpy((*minfo).connect_file, "CONNECT.TPL");
	strcpy((*minfo).radius_file, "RADIUS.DAT");
	strcpy((*minfo).inf_filename, "ATOMTYPE.INF");
	(*minfo).divcon = 2;
	(*minfo).longresname[0] = '\0';
	(*minfo).multiplicity = -9999;
	(*minfo).icharge = -9999;
	(*minfo).usercharge = -9999;
	(*minfo).dcharge = -9999.0;
	(*minfo).gv = 0;
	(*minfo).igkeyword = 0;
	(*minfo).eqcharge = -1;
}


/* Global variable containing a comment for the output file. */
char output_file_comment[MAXCHAR] = "";

/*
   Use the initial command line as the comment for the output file.
   Assumptions:  strncat correctly handles a non-positive 3rd argument.
 */
void create_output_file_comment(int argc, char *argv[])
{
    const char SEPARATOR[] = " ";
    int c; /* number of chars copied plus 1; may be incorrect on exit */
    int i;
    for (output_file_comment[0] = '\0', c = 1, i = 0;
         c < MAXCHAR && i < argc;
         c += strlen( argv[i] ), ++i) {
        strncat( output_file_comment, SEPARATOR, MAXCHAR - c );
        c += strlen( SEPARATOR );
        strncat( output_file_comment, argv[i], MAXCHAR - c );
    }
}


void default_cinfo(CONTROLINFO * cinfo)
{
	(*cinfo).intype[0] = '\0';
	(*cinfo).outtype[0] = '\0';
	(*cinfo).atype[0] = '\0';
	(*cinfo).chargetype[0] = '\0';
	(*cinfo).rnindex = 0;
	(*cinfo).intstatus = 1;
	(*cinfo).pfindex = 0;
	(*cinfo).prediction_index = 4;
	(*cinfo).bpindex = 1;
	(*cinfo).maxatom = MAXATOM;
	(*cinfo).maxbond = MAXBOND;
	(*cinfo).maxring = MAXRING;
	(*cinfo).max_path_length = -1;
	(*cinfo).verify_pdb_atomname = 1;
}



void default_inf(int atomnum, ATOM atom[], int index)
{
	int i;
	for (i = 0; i < atomnum; i++) {
		strcpy(atom[i].chain, " ");
		atom[i].ter = -1;
		if (index == 2)
			strcpy(atom[i].ambername, atom[i].name);
	}
}



int intcharge(int atomnum, ATOM atom[])
{
	int i;
	int icharge = 0;
	double fraction = 0.0;		/*decimal fraction */
	double dcharge = 0.0;		/*double precision charge */
	double tmpf;


	for (i = 0; i < atomnum; i++)
		dcharge += atom[i].charge;
	fraction = modf(dcharge, &tmpf);
	icharge = (int) tmpf;
	if (fabs(fraction) >= 0.50) {
		if (dcharge < 0)
			icharge--;
		if (dcharge > 0)
			icharge++;
	}
	return icharge;
}

void formula(int atomnum, ATOM atom[], char *form)
{
	int i, j;
	int countatom[NUMBER_OF_CHEMICAL_ELEMENTS];
	char tmpchar[MAXCHAR];
	char elemname[NUMBER_OF_CHEMICAL_ELEMENTS][5];

	for (i = 0; i < NUMBER_OF_CHEMICAL_ELEMENTS; i++) {
		countatom[i] = 0;
		strcpy(elemname[i], "X");
	}
        strcpy(elemname[1], "H");
        strcpy(elemname[2], "He");
        strcpy(elemname[3], "Li");
        strcpy(elemname[4], "Be");
        strcpy(elemname[5], "B");
        strcpy(elemname[6], "C");
        strcpy(elemname[7], "N");
        strcpy(elemname[8], "O");
        strcpy(elemname[9], "F");
        strcpy(elemname[10], "Ne");
        strcpy(elemname[11], "Na");
        strcpy(elemname[12], "Mg");
        strcpy(elemname[13], "Al");
        strcpy(elemname[14], "Si");
        strcpy(elemname[15], "P");
        strcpy(elemname[16], "S");
        strcpy(elemname[17], "Cl");
        strcpy(elemname[18], "Ar");
        strcpy(elemname[19], "K");
        strcpy(elemname[20], "Ca");
        strcpy(elemname[21], "Sc");
        strcpy(elemname[22], "Ti");
        strcpy(elemname[23], "V");
        strcpy(elemname[24], "Cr");
        strcpy(elemname[25], "Mn");
        strcpy(elemname[26], "Fe");
        strcpy(elemname[27], "Co");
        strcpy(elemname[28], "Ni");
        strcpy(elemname[29], "Cu");
        strcpy(elemname[30], "Zn");
        strcpy(elemname[31], "Ga");
        strcpy(elemname[32], "Ge");
        strcpy(elemname[33], "As");
        strcpy(elemname[34], "Se");
        strcpy(elemname[35], "Br");
        strcpy(elemname[36], "Kr");
        strcpy(elemname[37], "Rb");
        strcpy(elemname[38], "Sr");
        strcpy(elemname[39], "Y");
        strcpy(elemname[40], "Zr");
        strcpy(elemname[41], "Nb");
        strcpy(elemname[42], "Mo");
        strcpy(elemname[43], "Tc");
        strcpy(elemname[44], "Ru");
        strcpy(elemname[45], "Rh");
        strcpy(elemname[46], "Pd");
        strcpy(elemname[47], "Ag");
        strcpy(elemname[48], "Cd");
        strcpy(elemname[49], "In");
        strcpy(elemname[50], "Sn");
        strcpy(elemname[51], "Sb");
        strcpy(elemname[52], "Te");
        strcpy(elemname[53], "I");
        strcpy(elemname[54], "Xe");
        strcpy(elemname[55], "Cs");
        strcpy(elemname[56], "Ba");
        strcpy(elemname[57], "La");
/* Lanthanide Series*/
        strcpy(elemname[58], "Ce");
        strcpy(elemname[59], "Pr");
        strcpy(elemname[60], "Nd");
        strcpy(elemname[61], "Pm");
        strcpy(elemname[62], "Sm");
        strcpy(elemname[63], "Eu");
        strcpy(elemname[64], "Gd");
        strcpy(elemname[65], "Tb");
        strcpy(elemname[66], "Dy");
        strcpy(elemname[67], "Ho");
        strcpy(elemname[68], "Er");
        strcpy(elemname[69], "Tm");
        strcpy(elemname[70], "Yb");
        strcpy(elemname[71], "Lu");

        strcpy(elemname[72], "Hf");
        strcpy(elemname[73], "Ta");
        strcpy(elemname[74], "W");
        strcpy(elemname[75], "Re");
        strcpy(elemname[76], "Os");
        strcpy(elemname[77], "Ir");
        strcpy(elemname[78], "Pt");
        strcpy(elemname[79], "Au");
        strcpy(elemname[80], "Hg");
        strcpy(elemname[81], "Tl");
        strcpy(elemname[82], "Pb");
        strcpy(elemname[83], "Bi");
        strcpy(elemname[84], "Po");
        strcpy(elemname[85], "At");
        strcpy(elemname[86], "Rn");
        strcpy(elemname[87], "Fr");
        strcpy(elemname[88], "Ra");
        strcpy(elemname[89], "Ac");
/* Actinide Series*/
        strcpy(elemname[90], "Th");
        strcpy(elemname[91], "Pa");
        strcpy(elemname[92], "U");
        strcpy(elemname[93], "Np");
        strcpy(elemname[94], "Pu");
        strcpy(elemname[95], "Am");
        strcpy(elemname[96], "Cm");
        strcpy(elemname[97], "Bk");
        strcpy(elemname[98], "Cf");
        strcpy(elemname[99], "Es");
        strcpy(elemname[100], "Fm");
        strcpy(elemname[101], "Md");
        strcpy(elemname[102], "No");
        strcpy(elemname[103], "Lr");

        strcpy(elemname[104], "Rf");
        strcpy(elemname[105], "Db");
        strcpy(elemname[106], "Sg");
        strcpy(elemname[107], "Bh");
        strcpy(elemname[108], "Hs");
        strcpy(elemname[109], "Mt");
        strcpy(elemname[110], "Ds");
	strcpy(form, "");
	for (j = 0; j < atomnum; j++)
		countatom[atom[j].atomicnum]++;
	for (i = 0; i < NUMBER_OF_CHEMICAL_ELEMENTS; i++)
		if (countatom[i] >= 1) {
			strcat(form, elemname[i]);
			sprintf(tmpchar, "%d", countatom[i]);	
			/* newitoa(countatom[i], tmpchar); */
			strcat(form, tmpchar);
			strcat(form, " ");
		}
}

void initial(int num, ATOM * atom, char *resname)
{
	int i;
	for (i = 0; i < num; i++) {
		atom[i].connum = 0;
		atom[i].resno = 1;
		atom[i].con[0] = -1;
		atom[i].con[1] = -1;
		atom[i].con[2] = -1;
		atom[i].con[3] = -1;
		atom[i].con[4] = -1;
		atom[i].con[5] = -1;
		atom[i].charge = 0.0;
		atom[i].improper = 0;
		strcpy(atom[i].name, "NA");
		strcpy(atom[i].aa, resname);
	}
}


void read_radius(char *radius_parm, int atomnum, ATOM atom[])
{
	typedef struct {
		char name[5];
		double radius;
	} RADII;
	RADII radii[200];
	int i, j, k;
	int index;
	int num = 0;
	int parmnum;
	double tmpfloat;
	char line[MAXCHAR];
	char tmpchar[MAXCHAR];
	FILE *fpin;


	if ((fpin = fopen(radius_parm, "r")) == NULL) {
		fprintf(stdout, "Cannot open the radius parm file %s to read, exit\n", radius_parm);
		return;
	}
	for (;;) {
		if (fgets(line, LINELEN_MAX, fpin) == NULL)
			break;
		if (strncmp(line, "RADIUS", 6) == 0) {
			sscanf(&line[7], "%s%lf", tmpchar, &tmpfloat);
			strcpy(radii[num].name, tmpchar);
			radii[num++].radius = tmpfloat;
		}
	}
	fclose(fpin);
	parmnum = num;

	for (i = 0; i < atomnum; i++) {
		index = 0;
		for (j = 0; j < parmnum; j++)
			if (strcmp(radii[j].name, atom[i].ambername) == 0) {
				atom[i].radius = radii[j].radius;
				index = 1;
				break;
			}
		if (index == 0) {
			tmpchar[0] = '\0';
			num = 0;
			for (k = 0; k <= strlen(atom[i].ambername); k++)
				if ((atom[i].ambername[k] <= 'Z'
					 && atom[i].ambername[k] >= 'A')
					|| (atom[i].ambername[k] <= 'z'
						&& atom[i].ambername[k] >= 'a'))
					tmpchar[num++] = atom[i].ambername[k];
			tmpchar[num] = '\0';
			for (j = 0; j < parmnum; j++)
				if (strcmp(tmpchar, radii[j].name) == 0) {
					atom[i].radius = radii[j].radius;
					break;
				}
		}
	}
}



void adjustatomname(int atomnum, ATOM atom[], int check_atom_type)
{								/* just strips blanks out of the atomname */
	int i, j;
	int num, num1, num2;
	int tmpint;
	char tmpchar[MAXCHAR];
	char tmpchar1[MAXCHAR];
	char tmpchar2[MAXCHAR];
	char filename[MAXCHAR];
	char line[MAXCHAR];
	FILE *fp;
	int index;

	for (i = 0; i < atomnum; i++) {
		num = 0;
		tmpchar[0] = '\0';
		for (j = 0; j < strlen(atom[i].name); j++) {
			if (atom[i].name[j] == ' ')
				continue;
			tmpchar[num] = atom[i].name[j];
			if (tmpchar[num] == '\0')
				break;
			num++;
		}
		tmpchar[num] = '\0';
		strcpy(atom[i].name, tmpchar);
	}
	if(check_atom_type == 1) {
		strcpy(filename, amberhome);
		strcat(filename, "/dat/antechamber/CORR_NAME_TYPE.DAT");

		if ((fp = fopen(filename, "r")) == NULL) {
                	fprintf(stdout, "Cannot open the CORR_NAME_TYPE.DAT file %s to read, exit\n", filename);
                	exit(1);
        	}
		for (i = 0; i < atomnum; i++) {
			if(strcmp(atom[i].ambername, "Du") == 0 || strcmp(atom[i].ambername, "DU") == 0)
				continue;
			rewind(fp);	
			for(;;) {
			  	if (fgets(line, MAXCHAR, fp) == NULL)  break;
				sscanf(line, "%s%s", tmpchar1, tmpchar2);
				tmpint = strlen(tmpchar2);
				if(strcmp(tmpchar1, atom[i].ambername) == 0) {
					if(strncmp(atom[i].name, tmpchar2, tmpint)!=0)
						for(j=0;j<tmpint;j++)	
							atom[i].name[j] = tmpchar2[j];
					break;
				}
			}
		}
		fclose(fp);
	}
	for (i = 0; i < atomnum; i++) {
		num1 = 0;
		num2 = 0;
		index = 0;
		strcpy(tmpchar1, "");
		strcpy(tmpchar2, "");

		for (j = 0; j < strlen(atom[i].name); j++) {
			if (index == 0 && atom[i].name[0] >= '0'
				&& atom[i].name[0] <= '9')
				index = 1;
			if (index == 1
				&& (atom[i].name[j] >= '0' && atom[i].name[j] <= '9'))
				tmpchar1[num1++] = atom[i].name[j];
			if (index == 1
				&& (atom[i].name[j] < '0' || atom[i].name[j] > '9'))
				index = 2;
			if (index == 0 || index == 2)
				tmpchar2[num2++] = atom[i].name[j];
		}
		tmpchar2[num2] = '\0';
		if (index == 2) {
			tmpchar1[num1] = '\0';
			strcat(tmpchar2, tmpchar1);
		}
		strcpy(atom[i].name, tmpchar2);
	}
}


void atomname(int atomnum, ATOM * atom)
{
	int i, j;
	int countatom[NUMBER_OF_CHEMICAL_ELEMENTS];
	char tmpchar[MAXCHAR];

	for (i = 0; i < NUMBER_OF_CHEMICAL_ELEMENTS; i++)
		countatom[i] = 0;
	for (j = 0; j < atomnum; j++) {
		strcpy(atom[j].name, atom[j].element);
		countatom[atom[j].atomicnum]++;
		sprintf(tmpchar, "%d", countatom[atom[j].atomicnum]);
		/* newitoa(countatom[atom[j].atomicnum], tmpchar); */
		strcat(atom[j].name, tmpchar);
	}
}

void atomicnum(int atomnum, ATOM * atom)
{
/* type definition*/
/* 0 - non metal
   1 - metal
   2 - semi-metal
   3 - transition metal 
   4 - rare earth
   5 - nobel gas
*/
	int i;
	for (i = 0; i < atomnum; i++) {
		switch (atom[i].name[0]) {
		case 'A':
			if (atom[i].name[1] == 'c') {
				atom[i].atomicnum = 89;
				strcpy(atom[i].element, "Ac");
				atom[i].type = 4;
			}
			else if (atom[i].name[1] == 'g') {
				atom[i].atomicnum = 47;
				strcpy(atom[i].element, "Ag");
				atom[i].type = 3;
			}
			else if (atom[i].name[1] == 'l'){
				atom[i].atomicnum = 13;
				strcpy(atom[i].element, "Al");
				atom[i].type = 1;
			}
			else if (atom[i].name[1] == 'm'){
				atom[i].atomicnum = 95;
				strcpy(atom[i].element, "Am");
				atom[i].type = 4;
			}
			else if (atom[i].name[1] == 'r') {
				atom[i].atomicnum = 18;
				strcpy(atom[i].element, "Ar");
				atom[i].type = 5;
			}
			else if (atom[i].name[1] == 's') {
				atom[i].atomicnum = 33;
				strcpy(atom[i].element, "As");
				atom[i].type = 2;
			}
			else if (atom[i].name[1] == 't') {
				atom[i].atomicnum = 85;
				strcpy(atom[i].element, "At");
				atom[i].type = 0;
			}
			else if (atom[i].name[1] == 'u') {
				atom[i].atomicnum = 79;
				strcpy(atom[i].element, "Au");
				atom[i].type = 3;
			}
			break;
		case 'B':
			if (atom[i].name[1] == 'r' || atom[i].name[1] == 'R') {
				atom[i].atomicnum = 35;
				strcpy(atom[i].element, "Br");
				atom[i].type = 0;
			}
			else if (atom[i].name[1] == 'a') {
				atom[i].atomicnum = 56;
				strcpy(atom[i].element, "Ba");
				atom[i].type = 1;
			}
			else if (atom[i].name[1] == 'e') {
				atom[i].atomicnum = 4;
				strcpy(atom[i].element, "Be");
				atom[i].type = 1;
			}
			else if (atom[i].name[1] == 'h') {
				atom[i].atomicnum = 107;
				strcpy(atom[i].element, "Bh");
				atom[i].type = 3;
			}
			else if (atom[i].name[1] == 'i') {
				atom[i].atomicnum = 83;
				strcpy(atom[i].element, "Bi");
				atom[i].type = 1;
			}
			else if (atom[i].name[1] == 'k') {
				atom[i].atomicnum = 97;
				strcpy(atom[i].element, "Bk");
				atom[i].type = 1;
			}
			else {
				atom[i].atomicnum = 5;
				strcpy(atom[i].element, "B");
				atom[i].type = 2;
			}
			break;
		case 'C':
			if (atom[i].name[1] == 'l' || atom[i].name[1] == 'L') {
				atom[i].atomicnum = 17;
				strcpy(atom[i].element, "Cl");
				atom[i].type = 0;
			}
			else if (atom[i].name[1] == 'a') {
				atom[i].atomicnum = 20;
				strcpy(atom[i].element, "Ca");
				atom[i].type = 1;
			}
			else if (atom[i].name[1] == 'd')  {
				atom[i].atomicnum = 48;
				strcpy(atom[i].element, "Cd");
				atom[i].type = 3;
			}
			else if (atom[i].name[1] == 'e')  {
				atom[i].atomicnum = 58;
				strcpy(atom[i].element, "Ce");
				atom[i].type = 4;
			}
			else if (atom[i].name[1] == 'f')  {
				atom[i].atomicnum = 98;
				strcpy(atom[i].element, "Cf");
				atom[i].type = 4;
			}
			else if (atom[i].name[1] == 'm')  {
				atom[i].atomicnum = 96;
				strcpy(atom[i].element, "Cm");
				atom[i].type = 4;
			}
			else if (atom[i].name[1] == 'o') {
				atom[i].atomicnum = 27;
				strcpy(atom[i].element, "Co");
				atom[i].type = 3;
			}
			else if (atom[i].name[1] == 'r')  {
				atom[i].atomicnum = 24;
				strcpy(atom[i].element, "Cr");
				atom[i].type = 3;
			}
			else if (atom[i].name[1] == 's')  {
				atom[i].atomicnum = 55;
				strcpy(atom[i].element, "Cs");
				atom[i].type = 1;
			}
			else if (atom[i].name[1] == 'u')  {
				atom[i].atomicnum = 29;
				strcpy(atom[i].element, "Cu");
				atom[i].type = 3;
			}
			else {
				atom[i].atomicnum = 6;
				strcpy(atom[i].element, "C");
				atom[i].type = 0;
			}
			break;
		case 'D':
			if (atom[i].name[1] == 'b')  {
				atom[i].atomicnum = 105;
				strcpy(atom[i].element, "Db");
				atom[i].type = 3;
			}
			else if (atom[i].name[1] == 's')  {
				atom[i].atomicnum = 110;
				strcpy(atom[i].element, "Ds");
				atom[i].type = 3;
			}
			else if (atom[i].name[1] == 'y')  {
				atom[i].atomicnum = 66;
				strcpy(atom[i].element, "Dy");
				atom[i].type = 4;
			}
			else {
				atom[i].atomicnum = 1;
				strcpy(atom[i].element, "D");
				atom[i].type = 0;
			}
			break;
		case 'E':
			if (atom[i].name[1] == 'P') {
				atom[i].atomicnum = 0;
				strcpy(atom[i].element, "EP");
				atom[i].type = 0;
			}
			else if (atom[i].name[1] == 'r') {
				atom[i].atomicnum = 68;
				strcpy(atom[i].element, "Er");
				atom[i].type = 4;
			}
			else if (atom[i].name[1] == 's') {
				atom[i].atomicnum = 99;
				strcpy(atom[i].element, "Es");
				atom[i].type = 4;
			}
			else if (atom[i].name[1] == 'u') {
				atom[i].atomicnum = 63;
				strcpy(atom[i].element, "Eu");
				atom[i].type = 4;
			}
			break;
		case 'F':
			if (atom[i].name[1] == 'e') {
				atom[i].atomicnum = 26;
				strcpy(atom[i].element, "Fe");
				atom[i].type = 3;
			}
			else if (atom[i].name[1] == 'm') {
				atom[i].atomicnum = 100;
				strcpy(atom[i].element, "Fm");
				atom[i].type = 4;
			}
			else if (atom[i].name[1] == 'r') {
				atom[i].atomicnum = 87;
				strcpy(atom[i].element, "Fr");
				atom[i].type = 1;
			}
			else {
				atom[i].atomicnum = 9;
				strcpy(atom[i].element, "F");
				atom[i].type = 0;
			}
			break;
		case 'G':
			if (atom[i].name[1] == 'a') {
				atom[i].atomicnum = 31;
				strcpy(atom[i].element, "Ga");
				atom[i].type = 2;
			}
			else if (atom[i].name[1] == 'd') {
				atom[i].atomicnum = 64;
				strcpy(atom[i].element, "Gd");
				atom[i].type = 4;
			}
			else if (atom[i].name[1] == 'e') {
				atom[i].atomicnum = 32;
				strcpy(atom[i].element, "Ge");
				atom[i].type = 2;
			}
			break;
		case 'H':
			if (atom[i].name[1] == 'e') {
				atom[i].atomicnum = 2;
				strcpy(atom[i].element, "He");
				atom[i].type = 5;
			}
			else if (atom[i].name[1] == 'f') {
				atom[i].atomicnum = 72;
				strcpy(atom[i].element, "Hf");
				atom[i].type = 3;
			}
			else if (atom[i].name[1] == 'g') {
				atom[i].atomicnum = 80;
				strcpy(atom[i].element, "Hg");
				atom[i].type = 3;
			}
			else if (atom[i].name[1] == 'o') {
				atom[i].atomicnum = 67;
				strcpy(atom[i].element, "Ho");
				atom[i].type = 4;
			}
			else if (atom[i].name[1] == 's') {
				atom[i].atomicnum = 108;
				strcpy(atom[i].element, "Hs");
				atom[i].type = 3;
			}
			else {
				atom[i].atomicnum = 1;
				strcpy(atom[i].element, "H");
				atom[i].type = 0;
			}
			break;
		case 'I':
			if (atom[i].name[1] == 'n') {
				atom[i].atomicnum = 49;
				strcpy(atom[i].element, "In");
				atom[i].type = 1;
			}
			if (atom[i].name[1] == 'r') {
				atom[i].atomicnum = 77;
				strcpy(atom[i].element, "Ir");
				atom[i].type = 3;
			}
			else {
				atom[i].atomicnum = 53;
				strcpy(atom[i].element, "I");
				atom[i].type = 0;
			}
			break;
		case 'K':
			if (atom[i].name[1] == 'r') {
				atom[i].atomicnum = 36;
				strcpy(atom[i].element, "Kr");
				atom[i].type = 5;
			}
			else {
				atom[i].atomicnum = 19;
				strcpy(atom[i].element, "K");
				atom[i].type = 1;
			}
			break;
		case 'l':
			if (atom[i].name[1] == 'p')
				atom[i].atomicnum = 0;
				strcpy(atom[i].element, "lp");
				atom[i].type = 0;
			break;
		case 'L':
			if (atom[i].name[1] == 'i') {
				atom[i].atomicnum = 3;
				strcpy(atom[i].element, "Li");
				atom[i].type = 1;
			}
			else if (atom[i].name[1] == 'a') {
				atom[i].atomicnum = 57;
				strcpy(atom[i].element, "La");
				atom[i].type = 4;
			}
			else if (atom[i].name[1] == 'r') {
				atom[i].atomicnum = 103;
				strcpy(atom[i].element, "Lr");
				atom[i].type = 4;
			}
			else if (atom[i].name[1] == 'u') {
				atom[i].atomicnum = 71;
				strcpy(atom[i].element, "Lu");
				atom[i].type = 4;
			}
			else if (atom[i].name[1] == 'P') {
				atom[i].atomicnum = 0;
				strcpy(atom[i].element, "LP");
				atom[i].type = 0;
			}
			break;
		case 'M':
			if (atom[i].name[1] == 'n') {
				atom[i].atomicnum = 25;
				strcpy(atom[i].element, "Mn");
				atom[i].type = 3;
			}
			else if (atom[i].name[1] == 'g') {
				atom[i].atomicnum = 12;
				strcpy(atom[i].element, "Mg");
				atom[i].type = 1;
			}
			else if (atom[i].name[1] == 'd') {
				atom[i].atomicnum = 101;
				strcpy(atom[i].element, "Md");
				atom[i].type = 4;
			}
			else if (atom[i].name[1] == 'o') {
				atom[i].atomicnum = 42;
				strcpy(atom[i].element, "Mo");
				atom[i].type = 3;
			}
			else if (atom[i].name[1] == 't') {
				atom[i].atomicnum = 109;
				strcpy(atom[i].element, "Mt");
				atom[i].type = 3;
			}
			break;
		case 'N':
			if (atom[i].name[1] == 'i') {
				atom[i].atomicnum = 28;
				strcpy(atom[i].element, "Ni");
				atom[i].type = 3;
			}
			else if (atom[i].name[1] == 'a') {
				atom[i].atomicnum = 11;
				strcpy(atom[i].element, "Na");
				atom[i].type = 1;
			}
			else if (atom[i].name[1] == 'b') {
				atom[i].atomicnum = 41;
				strcpy(atom[i].element, "Nb");
				atom[i].type = 3;
			}
			else if (atom[i].name[1] == 'd') {
				atom[i].atomicnum = 60;
				strcpy(atom[i].element, "Nd");
				atom[i].type = 4;
			}
			else if (atom[i].name[1] == 'e') {
				atom[i].atomicnum = 10;
				strcpy(atom[i].element, "Ne");
				atom[i].type = 5;
			}
			else if (atom[i].name[1] == 'o') {
				atom[i].atomicnum = 102;
				strcpy(atom[i].element, "No");
				atom[i].type = 4;
			}
			else if (atom[i].name[1] == 'p') {
				atom[i].atomicnum = 93;
				strcpy(atom[i].element, "Np");
				atom[i].type = 4;
			}
			else  {
				atom[i].atomicnum = 7;
				strcpy(atom[i].element, "N");
				atom[i].type = 0;
			}
			break;
		case 'O':
			if (atom[i].name[1] == 's') {
				atom[i].atomicnum = 76;
				strcpy(atom[i].element, "Os");
				atom[i].type = 3;
			}
			else {
				atom[i].atomicnum = 8;
				strcpy(atom[i].element, "O");
				atom[i].type = 0;
			}
			break;
		case 'P':
			if (atom[i].name[1] == 'd') {
				atom[i].atomicnum = 46;
				strcpy(atom[i].element, "Pd");
				atom[i].type = 3;
			}
			else if (atom[i].name[1] == 't') {
				atom[i].atomicnum = 78;
				strcpy(atom[i].element, "Pt");
				atom[i].type = 3;
			}
			else if (atom[i].name[1] == 'b') {
				atom[i].atomicnum = 82;
				strcpy(atom[i].element, "Pb");
				atom[i].type = 3;
			}
			else if (atom[i].name[1] == 'a') {
				atom[i].atomicnum = 91;
				strcpy(atom[i].element, "Pa");
				atom[i].type = 4;
			}
			else if (atom[i].name[1] == 'm') {
				atom[i].atomicnum = 61;
				strcpy(atom[i].element, "Pm");
				atom[i].type = 4;
			}
			else if (atom[i].name[1] == 'o') {
				atom[i].atomicnum = 84;
				strcpy(atom[i].element, "Po");
				atom[i].type = 2;
			}
			else if (atom[i].name[1] == 'r') {
				atom[i].atomicnum = 59;
				strcpy(atom[i].element, "Pr");
				atom[i].type = 4;
			}
			else if (atom[i].name[1] == 'u') {
				atom[i].atomicnum = 94;
				strcpy(atom[i].element, "Pu");
				atom[i].type = 4;
			}
			else {
				atom[i].atomicnum = 15;
				strcpy(atom[i].element, "P");
				atom[i].type = 0;
			}
			break;
		case 'R':
			if (atom[i].name[1] == 'u') {
				atom[i].atomicnum = 44;
				strcpy(atom[i].element, "Ru");
				atom[i].type = 3;
			}
			else if (atom[i].name[1] == 'h') {
				atom[i].atomicnum = 45;
				strcpy(atom[i].element, "Rh");
				atom[i].type = 3;
			}
			else if (atom[i].name[1] == 'a') {
				atom[i].atomicnum = 88;
				strcpy(atom[i].element, "Ra");
				atom[i].type = 1;
			}
			else if (atom[i].name[1] == 'b') {
				atom[i].atomicnum = 37;
				strcpy(atom[i].element, "Rb");
				atom[i].type = 1;
			}
			else if (atom[i].name[1] == 'e') {
				atom[i].atomicnum = 75;
				strcpy(atom[i].element, "Re");
				atom[i].type = 3;
			}
			else if (atom[i].name[1] == 'f') {
				atom[i].atomicnum = 104;
				strcpy(atom[i].element, "Rf");
				atom[i].type = 3;
			}
			else if (atom[i].name[1] == 'n') {
				atom[i].atomicnum = 86;
				strcpy(atom[i].element, "Rn");
				atom[i].type = 5;
			}
			break;
		case 'S':
			if (atom[i].name[1] == 'i' || atom[i].name[1] == 'I') {
				atom[i].atomicnum = 14;
				strcpy(atom[i].element, "Si");
				atom[i].type = 2;
			}
			else if (atom[i].name[1] == 'c') {
				atom[i].atomicnum = 21;
				strcpy(atom[i].element, "Sc");
				atom[i].type = 3;
			}
			else if (atom[i].name[1] == 'e') {
				atom[i].atomicnum = 34;
				strcpy(atom[i].element, "Se");
				atom[i].type = 0;
			}
			else if (atom[i].name[1] == 'r') {
				atom[i].atomicnum = 38;
				strcpy(atom[i].element, "Sr");
				atom[i].type = 1;
			}
			else if (atom[i].name[1] == 'b') {
				atom[i].atomicnum = 51;
				strcpy(atom[i].element, "Sb");
				atom[i].type = 2;
			}
			else if (atom[i].name[1] == 'g') {
				atom[i].atomicnum = 106;
				strcpy(atom[i].element, "Sg");
				atom[i].type = 3;
			}
			else if (atom[i].name[1] == 'm') {
				atom[i].atomicnum = 62;
				strcpy(atom[i].element, "Sm");
				atom[i].type = 4;
			}
			else if (atom[i].name[1] == 'n') {
				atom[i].atomicnum = 50;
				strcpy(atom[i].element, "Sn");
				atom[i].type = 3;
			}
			else  {
				atom[i].atomicnum = 16;
				strcpy(atom[i].element, "S");
				atom[i].type = 0;
			}
			break;
		case 'T':
			if (atom[i].name[1] == 'i') {
				atom[i].atomicnum = 22;
				strcpy(atom[i].element, "Ti");
				atom[i].type = 3;
			}
			else if (atom[i].name[1] == 'l') {
				atom[i].atomicnum = 81;
				strcpy(atom[i].element, "Tl");
				atom[i].type = 1;
			}
			else if (atom[i].name[1] == 'a') {
				atom[i].atomicnum = 73;
				strcpy(atom[i].element, "Ta");
				atom[i].type = 3;
			}
			else if (atom[i].name[1] == 'b') {
				atom[i].atomicnum = 65;
				strcpy(atom[i].element, "Tb");
				atom[i].type = 4;
			}
			else if (atom[i].name[1] == 'c') {
				atom[i].atomicnum = 43;
				strcpy(atom[i].element, "Tc");
				atom[i].type = 3;
			}
			else if (atom[i].name[1] == 'e') {
				atom[i].atomicnum = 52;
				strcpy(atom[i].element, "Te");
				atom[i].type = 2;
			}
			else if (atom[i].name[1] == 'h') {
				atom[i].atomicnum = 90;
				strcpy(atom[i].element, "Th");
				atom[i].type = 4;
			}
			else if (atom[i].name[1] == 'm') {
				atom[i].atomicnum = 69;
				strcpy(atom[i].element, "Tm");
				atom[i].type = 4;
			}
			else {
				atom[i].atomicnum = 1;
				strcpy(atom[i].element, "T");
				atom[i].type = 0;
			}
			break;
		case 'U':
			atom[i].atomicnum = 92;
			strcpy(atom[i].element, "U");
			atom[i].type = 4;
			break;
		case 'V':
			atom[i].atomicnum = 23;
			strcpy(atom[i].element, "V");
			atom[i].type = 3;
			break;
		case 'W':
			atom[i].atomicnum = 74;
			strcpy(atom[i].element, "W");
			atom[i].type = 3;
			break;
		case 'X':
			if (atom[i].name[1] == 'e') {
				atom[i].atomicnum = 54;
				strcpy(atom[i].element, "Xe");
				atom[i].type = 5;
			}
			break;
		case 'Y':
			if (atom[i].name[1] == 'b') {
				atom[i].atomicnum = 70;
				strcpy(atom[i].element, "Yb");
				atom[i].type = 4;
			}
			else {
				atom[i].atomicnum = 39;
				strcpy(atom[i].element, "Y");
				atom[i].type = 3;
			}
			break;
		case 'Z':
			if (atom[i].name[1] == 'n') {
				strcpy(atom[i].element, "Zn");
				atom[i].type = 3;
				atom[i].atomicnum = 30;
			}
			else if (atom[i].name[1] == 'r') {
				atom[i].atomicnum = 40;
				strcpy(atom[i].element, "Zr");
				atom[i].type = 3;
			}
			break;
		default:
			printf("\n Unrecognized atomic name %5s, exit", atom[i].name);
		}
	}
}

void element(int atomnum, ATOM * atom)
{
	int i;
	for (i = 0; i < atomnum; i++)
		switch (atom[i].atomicnum) {
		case 1:
			strcpy(atom[i].element, "H");
			break;
		case 2:
			strcpy(atom[i].element, "He");
			break;
		case 3:
			strcpy(atom[i].element, "Li");
			break;
		case 4:
			strcpy(atom[i].element, "Be");
			break;
		case 5:
			strcpy(atom[i].element, "B");
			break;
		case 6:
			strcpy(atom[i].element, "C");
			break;
		case 7:
			strcpy(atom[i].element, "N");
			break;
		case 8:
			strcpy(atom[i].element, "O");
			break;
		case 9:
			strcpy(atom[i].element, "F");
			break;
		case 10:
			strcpy(atom[i].element, "Ne");
			break;
		case 11:
			strcpy(atom[i].element, "Na");
			break;
		case 12:
			strcpy(atom[i].element, "Mg");
			break;
		case 13:
			strcpy(atom[i].element, "Al");
			break;
		case 14:
			strcpy(atom[i].element, "Si");
			break;
		case 15:
			strcpy(atom[i].element, "P");
			break;
		case 16:
			strcpy(atom[i].element, "S");
			break;
		case 17:
			strcpy(atom[i].element, "Cl");
			break;
		case 18:
			strcpy(atom[i].element, "Ar");
			break;
		case 19:
			strcpy(atom[i].element, "K");
			break;
		case 20:
			strcpy(atom[i].element, "Ca");
			break;
		case 21:
			strcpy(atom[i].element, "Sc");
			break;
		case 22:
			strcpy(atom[i].element, "Ti");
			break;
		case 23:
			strcpy(atom[i].element, "V");
			break;
		case 24:
			strcpy(atom[i].element, "Cr");
			break;
		case 25:
			strcpy(atom[i].element, "Mn");
			break;
		case 26:
			strcpy(atom[i].element, "Fe");
			break;
		case 27:
			strcpy(atom[i].element, "Co");
			break;
		case 28:
			strcpy(atom[i].element, "Ni");
			break;
		case 29:
			strcpy(atom[i].element, "Cu");
			break;
		case 30:
			strcpy(atom[i].element, "Zn");
			break;
		case 31:
			strcpy(atom[i].element, "Ga");
			break;
		case 32:
			strcpy(atom[i].element, "Ge");
			break;
		case 33:
			strcpy(atom[i].element, "As");
			break;
		case 34:
			strcpy(atom[i].element, "Se");
			break;
		case 35:
			strcpy(atom[i].element, "Br");
			break;
		case 36:
			strcpy(atom[i].element, "Kr");
			break;
		case 37:
			strcpy(atom[i].element, "Rb");
			break;
		case 38:
			strcpy(atom[i].element, "Sr");
			break;
		case 39:
			strcpy(atom[i].element, "Y");
			break;
		case 40:
			strcpy(atom[i].element, "Zr");
			break;
		case 41:
			strcpy(atom[i].element, "Nb");
			break;
		case 42:
			strcpy(atom[i].element, "Mo");
			break;
		case 43:
			strcpy(atom[i].element, "Tc");
			break;
		case 44:
			strcpy(atom[i].element, "Ru");
			break;
		case 45:
			strcpy(atom[i].element, "Rh");
			break;
		case 46:
			strcpy(atom[i].element, "Pd");
			break;
		case 47:
			strcpy(atom[i].element, "Ag");
			break;
		case 48:
			strcpy(atom[i].element, "Cd");
			break;
		case 49:
			strcpy(atom[i].element, "In");
			break;
		case 50:
			strcpy(atom[i].element, "Sn");
			break;
		case 51:
			strcpy(atom[i].element, "Sb");
			break;
		case 52:
			strcpy(atom[i].element, "Te");
			break;
		case 53:
			strcpy(atom[i].element, "I");
			break;
		case 54:
			strcpy(atom[i].element, "Xe");
			break;
		case 55:
			strcpy(atom[i].element, "Cs");
			break;
		case 56:
			strcpy(atom[i].element, "Ba");
			break;
		case 57:
			strcpy(atom[i].element, "La");
			break;
		case 58:
			strcpy(atom[i].element, "Ce");
			break;
		case 59:
			strcpy(atom[i].element, "Pr");
			break;
		case 60:
			strcpy(atom[i].element, "Nd");
			break;
		case 61:
			strcpy(atom[i].element, "Pm");
			break;
		case 62:
			strcpy(atom[i].element, "Sm");
			break;
		case 63:
			strcpy(atom[i].element, "Eu");
			break;
		case 64:
			strcpy(atom[i].element, "Gd");
			break;
		case 65:
			strcpy(atom[i].element, "Tb");
			break;
		case 66:
			strcpy(atom[i].element, "Dy");
			break;
		case 67:
			strcpy(atom[i].element, "Ho");
			break;
		case 68:
			strcpy(atom[i].element, "Er");
			break;
		case 69:
			strcpy(atom[i].element, "Tm");
			break;
		case 70:
			strcpy(atom[i].element, "Yb");
			break;
		case 71:
			strcpy(atom[i].element, "Lu");
			break;
		case 72:
			strcpy(atom[i].element, "Hf");
			break;
		case 73:
			strcpy(atom[i].element, "Ta");
			break;
		case 74:
			strcpy(atom[i].element, "W");
			break;
		case 75:
			strcpy(atom[i].element, "Re");
			break;
		case 76:
			strcpy(atom[i].element, "Os");
			break;
		case 77:
			strcpy(atom[i].element, "Ir");
			break;
		case 78:
			strcpy(atom[i].element, "Pt");
			break;
		case 79:
			strcpy(atom[i].element, "Au");
			break;
		case 80:
			strcpy(atom[i].element, "Hg");
			break;
		case 81:
			strcpy(atom[i].element, "Tl");
			break;
		case 82:
			strcpy(atom[i].element, "Pb");
			break;
		case 83:
			strcpy(atom[i].element, "Bi");
			break;
		case 84:
			strcpy(atom[i].element, "Po");
			break;
		case 85:
			strcpy(atom[i].element, "At");
			break;
		case 86:
			strcpy(atom[i].element, "Rn");
			break;
		case 87:
			strcpy(atom[i].element, "Fr");
			break;
		case 88:
			strcpy(atom[i].element, "Ra");
			break;
		case 89:
			strcpy(atom[i].element, "Ac");
			break;
		case 90:
			strcpy(atom[i].element, "Th");
			break;
		case 91:
			strcpy(atom[i].element, "Pa");
			break;
		case 92:
			strcpy(atom[i].element, "U");
			break;
		case 93:
			strcpy(atom[i].element, "Np");
			break;
		case 94:
			strcpy(atom[i].element, "Pu");
			break;
		case 95:
			strcpy(atom[i].element, "Am");
			break;
		case 96:
			strcpy(atom[i].element, "Cm");
			break;
		case 97:
			strcpy(atom[i].element, "Bk");
			break;
		case 98:
			strcpy(atom[i].element, "Cf");
			break;
		case 99:
			strcpy(atom[i].element, "Es");
			break;
		case 100:
			strcpy(atom[i].element, "Fm");
			break;
		case 101:
			strcpy(atom[i].element, "Md");
			break;
		case 102:
			strcpy(atom[i].element, "No");
			break;
		case 103:
			strcpy(atom[i].element, "Lr");
			break;
		case 104:
			strcpy(atom[i].element, "Rf");
			break;
		case 105:
			strcpy(atom[i].element, "Db");
			break;
		case 106:
			strcpy(atom[i].element, "Sg");
			break;
		case 107:
			strcpy(atom[i].element, "Bh");
			break;
		case 108:
			strcpy(atom[i].element, "Hs");
			break;
		case 109:
			strcpy(atom[i].element, "Mt");
			break;
		case 110:
			strcpy(atom[i].element, "Ds");
			break;
		default:
			strcpy(atom[i].element, "du");
			break;
		}
}


void duplicatedname(int atomnum, ATOM * atom)
{
	int i, j, k;
	int index;
	int id;
	char tmpchar[MAXCHAR];
	char name[MAXCHAR];
	char elemname[NUMBER_OF_CHEMICAL_ELEMENTS][5];
	for (i = 0; i < NUMBER_OF_CHEMICAL_ELEMENTS; i++)
		strcpy(elemname[i], "X");
        strcpy(elemname[1], "H");
        strcpy(elemname[2], "He");
        strcpy(elemname[3], "Li");
        strcpy(elemname[4], "Be");
        strcpy(elemname[5], "B");
        strcpy(elemname[6], "C");
        strcpy(elemname[7], "N");
        strcpy(elemname[8], "O");
        strcpy(elemname[9], "F");
        strcpy(elemname[10], "Ne");
        strcpy(elemname[11], "Na");
        strcpy(elemname[12], "Mg");
        strcpy(elemname[13], "Al");
        strcpy(elemname[14], "Si");
        strcpy(elemname[15], "P");
        strcpy(elemname[16], "S");
        strcpy(elemname[17], "Cl");
        strcpy(elemname[18], "Ar");
        strcpy(elemname[19], "K");
        strcpy(elemname[20], "Ca");
        strcpy(elemname[21], "Sc");
        strcpy(elemname[22], "Ti");
        strcpy(elemname[23], "V");
        strcpy(elemname[24], "Cr");
        strcpy(elemname[25], "Mn");
        strcpy(elemname[26], "Fe");
        strcpy(elemname[27], "Co");
        strcpy(elemname[28], "Ni");
        strcpy(elemname[29], "Cu");
        strcpy(elemname[30], "Zn");
        strcpy(elemname[31], "Ga");
        strcpy(elemname[32], "Ge");
        strcpy(elemname[33], "As");
        strcpy(elemname[34], "Se");
        strcpy(elemname[35], "Br");
        strcpy(elemname[36], "Kr");
        strcpy(elemname[37], "Rb");
        strcpy(elemname[38], "Sr");
        strcpy(elemname[39], "Y");
        strcpy(elemname[40], "Zr");
        strcpy(elemname[41], "Nb");
        strcpy(elemname[42], "Mo");
        strcpy(elemname[43], "Tc");
        strcpy(elemname[44], "Ru");
        strcpy(elemname[45], "Rh");
        strcpy(elemname[46], "Pd");
        strcpy(elemname[47], "Ag");
        strcpy(elemname[48], "Cd");
        strcpy(elemname[49], "In");
        strcpy(elemname[50], "Sn");
        strcpy(elemname[51], "Sb");
        strcpy(elemname[52], "Te");
        strcpy(elemname[53], "I");
        strcpy(elemname[54], "Xe");
        strcpy(elemname[55], "Cs");
        strcpy(elemname[56], "Ba");
        strcpy(elemname[57], "La");
/* Lanthanide Series*/
        strcpy(elemname[58], "Ce");
        strcpy(elemname[59], "Pr");
        strcpy(elemname[60], "Nd");
        strcpy(elemname[61], "Pm");
        strcpy(elemname[62], "Sm");
        strcpy(elemname[63], "Eu");
        strcpy(elemname[64], "Gd");
        strcpy(elemname[65], "Tb");
        strcpy(elemname[66], "Dy");
        strcpy(elemname[67], "Ho");
        strcpy(elemname[68], "Er");
        strcpy(elemname[69], "Tm");
        strcpy(elemname[70], "Yb");
        strcpy(elemname[71], "Lu");

        strcpy(elemname[72], "Hf");
        strcpy(elemname[73], "Ta");
        strcpy(elemname[74], "W");
        strcpy(elemname[75], "Re");
        strcpy(elemname[76], "Os");
        strcpy(elemname[77], "Ir");
        strcpy(elemname[78], "Pt");
        strcpy(elemname[79], "Au");
        strcpy(elemname[80], "Hg");
        strcpy(elemname[81], "Tl");
        strcpy(elemname[82], "Pb");
        strcpy(elemname[83], "Bi");
        strcpy(elemname[84], "Po");
        strcpy(elemname[85], "At");
        strcpy(elemname[86], "Rn");
        strcpy(elemname[87], "Fr");
        strcpy(elemname[88], "Ra");
        strcpy(elemname[89], "Ac");
/* Actinide Series*/
        strcpy(elemname[90], "Th");
        strcpy(elemname[91], "Pa");
        strcpy(elemname[92], "U");
        strcpy(elemname[93], "Np");
        strcpy(elemname[94], "Pu");
        strcpy(elemname[95], "Am");
        strcpy(elemname[96], "Cm");
        strcpy(elemname[97], "Bk");
        strcpy(elemname[98], "Cf");
        strcpy(elemname[99], "Es");
        strcpy(elemname[100], "Fm");
        strcpy(elemname[101], "Md");
        strcpy(elemname[102], "No");
        strcpy(elemname[103], "Lr");

        strcpy(elemname[104], "Rf");
        strcpy(elemname[105], "Db");
        strcpy(elemname[106], "Sg");
        strcpy(elemname[107], "Bh");
        strcpy(elemname[108], "Hs");
        strcpy(elemname[109], "Mt");
        strcpy(elemname[110], "Ds");

	for (i = 0; i < atomnum; i++)
		for (j = i + 1; j < atomnum; j++) 
			if (strcmp(atom[i].name, atom[j].name) == 0 && atom[i].resno == atom[j].resno) {
				id = 1;
				index = 1;
				while (index) {
					sprintf(tmpchar, "%d", id);
					/* newitoa(id, tmpchar); */
					strcpy(name, elemname[atom[j].atomicnum]);
					strcat(name, tmpchar);
					for (k = 0; k < atomnum; k++)
						if (strncmp(name, atom[k].name,strlen(name)) == 0) {
							index = 0;
							break;
						}
					if (index == 1) 
						break;	
					if (index == 0) {
						id++;
					 	index = 1;	
					}
				}
				strcpy(atom[j].name, name);
			}
}

double translate(char *str) {
        int i;
        int pos = 0;
        int flag = 0;
        char tmpc1[MAXCHAR];
        char tmpc2[MAXCHAR];
        double f, e;

        for(i=0;i<strlen(str);i++) {
                if(str[i]=='D' || str[i]=='d'||str[i]=='E'||str[i]=='e') {
                        pos = i;
                        flag = 1;
                        continue;
                }
                if(flag == 0)
                        tmpc1[i]=str[i];
                if(flag == 1)
                        tmpc2[i-pos-1]=str[i];

        }

        f=atof(tmpc1);
        e=atof(tmpc2);
        return f*pow(10.0,e);
}

void info(int atomnum, ATOM atom[], int bondnum, BOND bond[], AROM arom[],
		  CONTROLINFO cinfo, MOLINFO minfo)
{
	int i;
	printf("\natomnum: %d", atomnum);
	printf("\nbondnum: %d", bondnum);

	printf("\n-----------------------------------------\n");
	for (i = 0; i < atomnum; i++)
		printf("ATOM%7d  %-4s%-4s%5d%5s%5s%12.3f%8.3f%8.3f%8d\n",
			   i + 1, atom[i].name, atom[i].aa, atom[i].resno,
			   atom[i].ambername, atom[i].element, atom[i].x, atom[i].y,
			   atom[i].z, atom[i].atomicnum);

	printf("\n-----------------------------------------\n");

	for (i = 0; i < bondnum; i++)
		printf("BOND  %5d  %5d  %5d  %5s  %5s  %5d\n", i + 1,
			   bond[i].bondi, bond[i].bondj, atom[bond[i].bondi].name,
			   atom[bond[i].bondj].name, bond[i].type);

	printf("\n-----------------------------------------\n");

	for (i = 0; i < atomnum; i++) {
		printf("\nFor Atom %d %s\n\n", i + 1, atom[i].name);
		printf("RING nr: %5d\n", arom[i].nr);
		printf("RING ar1=%d, ar2=%d, ar3=%d, ar4=%d, ar5=%d \n",
			   arom[i].ar1, arom[i].ar2, arom[i].ar3, arom[i].ar4,
			   arom[i].ar5);
		printf("RING rg3=%d, rg4=%d, rg5=%d, rg6=%d, rg7=%d \n",
			   arom[i].rg[3], arom[i].rg[4], arom[i].rg[5], arom[i].rg[6],
			   arom[i].rg[7]);
	}
	printf("\n-----------------------------------------\n");

	printf("\ndkeyword: %s", minfo.dkeyword);
	printf("\nmkeyword: %s", minfo.mkeyword);
	printf("\nekeyword: %s", minfo.ekeyword);
	printf("\ngkeyword: %s", minfo.gkeyword);
	printf("\nresname: %s", minfo.resname);
	printf("\natom_type_def: %s", minfo.atom_type_def);
	printf("\nresfilename: %s", minfo.resfilename);
	printf("\ngfilename: %s", minfo.gfilename);
	printf("\nconnect_file: %s", minfo.connect_file);
	printf("\nradius_file: %s", minfo.radius_file);
	printf("\ninf_filename: %s", minfo.inf_filename);
	printf("\nlongresname: %s", minfo.longresname);
	printf("\nmultiplicity: %d", minfo.multiplicity);
	printf("\nicharge: %d", minfo.icharge);
	printf("\nusercharge: %d", minfo.usercharge);
	printf("\ndcharge: %lf", minfo.dcharge);

	printf("\n-----------------------------------------\n");
	printf("\nintype: %s", cinfo.intype);
	printf("\nouttype: %s", cinfo.outtype);
	printf("\natype: %s", cinfo.atype);
	printf("\nchargetype: %s", cinfo.chargetype);
	printf("\nmaxatom: %d", cinfo.maxatom);
	printf("\nmaxbond: %d", cinfo.maxbond);
	printf("\nmaxring: %d", cinfo.maxring);
	printf("\nrnindex: %d", cinfo.rnindex);
	printf("\nintstatus: %d", cinfo.intstatus);
	printf("\npfindex: %d", cinfo.pfindex);
	printf("\nprediction_index %d", cinfo.prediction_index);
	printf("\nbpindex: %d", cinfo.bpindex);
	printf("\n-----------------------------------------\n\n");


}

void adjust_sequence_order(int atomnum, ATOM atom[], int bondnum, BOND bond[])
{
        int i, j;
	int flag;
        int suc_flag = 1;
	int *seq;
	int *select;
	ATOM *atm;
	int num = 0;


/* first of all check if the sequence order is OK */
        for (i = 1; i < atomnum; i++) {
		flag = 0;
		for(j=0;j<atom[i].connum;j++)
			if(atom[i].con[j] < i) {
				flag = 1;
				break;
			}
/* if flag == 0, the sequence order has some problem */
		if(flag == 0) {
			suc_flag = 0;
			break;
		}
	}
	if(suc_flag == 0) {
		fprintf(stdout, "WARNING: the sequence order of the atoms has some problem, readjust it\n");	

        	select = (int *) malloc(sizeof(int) * atomnum);
        	if (select == NULL) {
                	fprintf(stdout, "memory allocation error for *select in the adjust_sequence_order function of common.c\n");
                	exit(1);
        	}

        	seq = (int *) malloc(sizeof(int) * atomnum);
        	if (seq == NULL) {
                	fprintf(stdout, "memory allocation error for *select in the adjust_sequence_order function of common.c\n");
                	exit(1);
        	}

        	atm = (ATOM *) malloc(sizeof(ATOM) * atomnum);
        	if (seq == NULL) {
                	fprintf(stdout, "memory allocation error for *atm in the adjust_sequence_order function of common.c\n");
                	exit(1);
        	}

		for(i = 0; i< atomnum; i++) {
			select[i] = -1;
			seq[i] = -1;
		}
		seq[0] = 0;	
		select[0] = 0;
		num = 1;
		flag = 1;
		while(flag == 1) {
			flag = 0;
			for(i = 0; i< atomnum; i++) {
				if(select[i] == -1) continue;
				for(j=0; j< atom[i].connum; j++) 
					if(select[atom[i].con[j]] == -1) {
						select[atom[i].con[j]] = num;	
						seq[num] = atom[i].con[j];
						num++;
						flag = 1;
					}
			}
		}
/* rearrange the sequence order of the atoms */
		for(i=0;i<atomnum;i++) 
			atm[i] = atom[seq[i]];
		for(i=0;i<atomnum;i++) 
			for(j=0;j<atm[i].connum;j++) 
				atm[i].con[j] = select[atm[i].con[j]];
		for(i=0;i<atomnum;i++) 
			atom[i] = atm[i];	
		for(i=0;i<bondnum;i++)  {
			bond[i].bondi = select[bond[i].bondi];
			bond[i].bondj = select[bond[i].bondj];
		}
		free(seq);
		free(select);
		free(atm);
		return;
	}
}

