
/*
************************************************************************
*           All Copyright Reserved!                                    *
*                                                                      *
*  Prog:    delphigen                                                  *
*  Version: version 1.0                                                *
*  Author:  Junmei Wang                                                *
*                                                                      *
*  Department of Pharmaceutical Chemistry                              *
*  School of Pharmacy                                                  *
*  University of California                                            *
*  San Francisco   CA 94143                                            *
*  Octomber, 2001                                                      *
************************************************************************
*/

# include "common.h"
# include "define.h"
# include "atom.h"
# include "common.c"
# include "ac.c"
# include "pdb.c"

char *amberhome;
ATOM *atom;
BOND *bond;
int atomnum = 0;
int bondnum = 0;
int cindex = 0;
int radiusindex = 0;
int pindex = 0;
int mindex = 0;
int ch_index = 1;
int rh_index = 1;
CONTROLINFO cinfo;
MOLINFO minfo;

typedef struct {
	char name[5];
	double radius;
} RADII;
RADII radii[200];
int parmnum = 0;

char line[MAXCHAR];
char ifilename[MAXCHAR];
char cfilename[MAXCHAR];
char rfilename[MAXCHAR];
char mfilename[MAXCHAR];
char pfilename[MAXCHAR];

FILE *fp;
FILE *fpin;
FILE *fpoutc;
FILE *fpoutr;
FILE *fpoutp;

void radius(void)
{
	int i, j, k;
	int index;
	int num = 0;
	char tmpchar[6];
	double tmpfloat;
	FILE *fpin;
	FILE *fpout;


	if ((fpin = fopen(pfilename, "r")) == NULL) {
		printf("\n Cannot open parameter file (pfilename): %s in radius(), exit", pfilename);
		exit(1);
	}
	if (radiusindex == 1) {
		index = 1;
		if ((fpout = fopen(rfilename, "r")) == NULL) {
/*     printf("\n Cannot open file %s, exit", rfilename);*/
			index = 0;
		}
		if (index == 1)
			fclose(fpout);
		if ((fpout = fopen(rfilename, "a")) == NULL) {
			printf("\n Cannot open radius file (rfilename): %s in radius(), exit", rfilename);
			exit(1);
		}
		if (index == 0 && rh_index == 1)
			fprintf(fpout, "atom__res_radius_color_\n");
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
		if (radiusindex == 1)
			fprintf(fpout, "%-5s%-5s%8.5lf\n", atom[i].name, atom[i].aa,
					atom[i].radius);

	}
	if (radiusindex == 1)
		fclose(fpout);
}

void wdelc(char *filename)
{
	int i;
	int index;
	char line[180];
	FILE *fpin;
	FILE *fpout;

	index = 1;
	if ((fpin = fopen(filename, "r")) == NULL)
		index = -1;
	else {
		if (fgets(line, 180, fpin) == NULL)
			index = -1;
		else if (fgets(line, 180, fpin) == NULL)
			index = -1;
	}
	if (index == 1)
		fclose(fpin);
	/* the above code judge if the file exists and if this file has a remark line */

	if ((fpout = fopen(filename, "a")) == NULL) {
		printf("\n Cannot open a file (filename) to append: %s in wdelc(), exit", filename);
		exit(1);
	}
	if (index == -1 && ch_index == 1)
		fprintf(fpout, "atom__resnumbc_charge_\n");
	for (i = 0; i < atomnum; i++)
		fprintf(fpout, "%-6s%c%c%c%13.5lf\n", atom[i].name, atom[i].aa[0],
				atom[i].aa[1], atom[i].aa[2], atom[i].charge);

	fclose(fpout);
}

int main(int argc, char *argv[])
{
	int i;
	int overflow_flag = 0;			/*if overflow_flag ==1, reallocate memory */

    amberhome = (char *) getenv("AMBERHOME");
    if( amberhome == NULL ){
       fprintf( stdout, "AMBERHOME is not set!\n" );
       exit(1);
    }
	if (strcmp(COLORTEXT, "YES") == 0 || strcmp(COLORTEXT, "yes") == 0) {
		if (argc == 2
			&& (strcmp(argv[1], "-h") == 0
				|| strcmp(argv[1], "-H") == 0)) {
			printf("[31mUsage: delphigen -i[0m input file name (ac)\n"
				   "[31m                 -c[0m charge file name \n"
				   "[31m                 -r[0m radius file name \n"
				   "[31m                 -m[0m modified pdb file name (optional)\n"
				   "[31m                 -p[0m radius parameter file name (optional)\n"
				   "[31m                 -ch[0m head line in charge file (Yes (the default) or No, optional)\n"
				   "[31m                 -rh[0m head line in radius file (Yes (the default) or No, optional)\n");

			exit(1);
		}
		if (argc != 5 && argc != 7 && argc != 9 && argc != 11) {
			printf("[31mUsage: delphigen -i[0m input file name (ac)\n"
				   "[31m                 -c[0m charge file name \n"
				   "[31m                 -r[0m radius file name \n"
				   "[31m                 -m[0m modified pdb file name (optional)\n"
				   "[31m                 -p[0m radius parameter file name (optional)\n"
				   "[31m                 -ch[0m head line in charge file (Yes (the default) or No, optional)\n"
				   "[31m                 -rh[0m head line in radius file (Yes (the default) or No, optional)\n");
			exit(1);
		}
	} else {
		if (argc == 2)
			if (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "-H") == 0) {
				printf("Usage: delphigen -i  input file name (ac) \n");
				printf("                 -c  charge file name \n");
				printf("                 -r  radius file name \n");
				printf
					("                 -m  modified pdb file name (optional)\n");
				printf
					("                 -p  radius parameter file name (optional)\n");
				printf
					("                 -ch  head line in charge file (Yes (the default) or No, optional)\n");
				printf
					("                 -rh  head line in radius file (Yes (the default) or No, optional)\n");
				exit(1);
			}
		if (argc != 5 && argc != 7 && argc != 9 && argc != 11) {
			printf("Usage: delphigen -i  input file name (ac) \n");
			printf("                 -c  charge file name \n");
			printf("                 -r  radius file name \n");
			printf
				("                 -m  modified pdb file name (optional)\n");
			printf
				("                 -p  radius parameter file name (optional)\n");
			printf
				("                 -ch  head line in charge file (Yes (the default) or No, optional)\n");
			printf
				("                 -rh  head line in radius file (Yes (the default) or No, optional)\n");
			exit(1);
		}
	}
	for (i = 1; i < argc; i += 2) {
		if (strcmp(argv[i], "-i") == 0)
			strcpy(ifilename, argv[i + 1]);
		if (strcmp(argv[i], "-c") == 0) {
			strcpy(cfilename, argv[i + 1]);
			cindex = 1;
		}
		if (strcmp(argv[i], "-r") == 0) {
			strcpy(rfilename, argv[i + 1]);
			radiusindex = 1;
		}
		if (strcmp(argv[i], "-m") == 0) {
			strcpy(mfilename, argv[i + 1]);
			mindex = 1;
		}
		if (strcmp(argv[i], "-p") == 0) {
			strcpy(pfilename, argv[i + 1]);
			pindex = 1;
		}
		if (strcmp(argv[i], "-ch") == 0) {
			if(strcmp("Yes", argv[i + 1])==0 ||strcmp("YES", argv[i + 1])==0 ||strcmp("yes", argv[i + 1])==0);
				ch_index = 1;
			if(strcmp("Y", argv[i + 1])==0 ||strcmp("y", argv[i + 1])==0);
				ch_index = 1;
			if(strcmp("No", argv[i + 1])==0 ||strcmp("NO", argv[i + 1])==0 ||strcmp("no", argv[i + 1])==0);
				ch_index = 0;
			if(strcmp("N", argv[i + 1])==0 ||strcmp("n", argv[i + 1])==0);
				ch_index = 0;
		}
		if (strcmp(argv[i], "-rh") == 0) {
			if(strcmp("Yes", argv[i + 1])==0 ||strcmp("YES", argv[i + 1])==0 ||strcmp("yes", argv[i + 1])==0);
				rh_index = 1;
			if(strcmp("Y", argv[i + 1])==0 ||strcmp("y", argv[i + 1])==0);
				rh_index = 1;
			if(strcmp("No", argv[i + 1])==0 ||strcmp("NO", argv[i + 1])==0 ||strcmp("no", argv[i + 1])==0);
				rh_index = 0;
			if(strcmp("N", argv[i + 1])==0 ||strcmp("n", argv[i + 1])==0);
				rh_index = 0;
		}
	}

/* for radius file */

	if (pindex != 1) {
		build_dat_path(pfilename, "RADIUS.DAT", sizeof pfilename, 0);
	}

/* allocate memory*/
	default_cinfo(&cinfo);
	default_minfo(&minfo);
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
	int i;
	for (i = 0; i < cinfo.maxbond; ++i) {
		bond[i].jflag = -1; /* bond type has not been assigned */
	}

	overflow_flag =
		rac(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo);
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
		int i;
		for (i = 0; i < cinfo.maxbond; ++i) {
			bond[i].jflag = -1; /* bond type has not been assigned */
		}
		overflow_flag =
			rac(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo);
	}
	atomicnum(atomnum, atom);
	adjustatomname(atomnum, atom, 1);

/* The following codes readin input file */
	radius();
	if (cindex == 1)
		wdelc(cfilename);
	if (mindex == 1)
		wmpdb(mfilename, atomnum, atom);
/* The end */
	return (0);
}
