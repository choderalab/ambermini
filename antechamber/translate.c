/*
************************************************************************
*           All Copyright Reserved!                                    *
*                                                                      *
*  Prog:    translate                                                    *
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
char *amberhome;
# include "common.h"
# include "define.h"
# include "atom.h"
# include "utility.c"
# include "common.c"
# include "rotate.c"
# include "ac.c"
# include "pdb.c"
# include "mol2.c"
# include "prep.c"
# include "lsfit.c"

ATOM *atom;
AROM *arom;
BOND *bond;
ATOM *ref_atom;
AROM *ref_arom;
BOND *ref_bond;
ATOM tmpatom1;
ATOM tmpatom2;
RING ring[MAXRING];
MOLINFO minfo;
CONTROLINFO cinfo;
int atomnum = 0;
int bondnum = 0;
int ringnum;
char ifilename[MAXCHAR];
char rfilename[MAXCHAR];
char ofilename[MAXCHAR];
char line[MAXCHAR];
int i, j, k;
int at1 = -99999;
int at2 = -99999;
int at3 = -99999;
double vectx = -99999.0;
double vecty = -99999.0;
double vectz = -99999.0;
double coord_x1 = -99999.0;
double coord_y1 = -99999.0;
double coord_z1 = -99999.0;
double coord_x2 = -99999.0;
double coord_y2 = -99999.0;
double coord_z2 = -99999.0;
double degree = -99999.0;
double sum_coordx = 0.0;
double sum_coordy = 0.0;
double sum_coordz = 0.0;
double coordx;
double coordy;
double coordz;
double rmsd;
double w1, w2;
double minx = 99999, maxx = -99999;
double miny = 99999, maxy = -99999;
double minz = 99999, maxz = -99999;
int iflag = 0;
int oflag = 0;
int rflag = 0;
int cflag = 0;
char command[MAXCHAR];
FILE *fp, *fpout;

int main(int argc, char *argv[]) {
	int i;
	int format;
	default_cinfo(&cinfo);
	default_minfo(&minfo);
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
			printf("[31mUsage: translate -i [0m input file name (pdb, ac or mol2)\n"
				   "[31m                 -o [0m output file name\n"
				   "[31m                 -r [0m reference file name\n"
				   "[31m                 -f [0m file format\n"
				   "[31m                 -c [0m command (center, translate, rotate1, rotate2, match)\n"
				   "[34m                     check:	  [0m\n"
				   "[34m                     center:    [0m need -a1;\n"
				   "[34m                     translate: [0m need -vx, -vy and -vz;\n"
				   "[34m                     rotate1:   [0m need -a1, -a2 and -d;\n"
				   "[34m                     rotate2:   [0m need -x1, -y1, -z1, -x2, -y2, -z2 and -d;\n"
				   "[34m                     match:     [0m need -r;\n"
				   "[34m                     alignx:    [0m align to X-axis, need -x1, -y1, -z1, -x2, -y2, -z2;\n"
				   "[34m                     aligny:    [0m align to Y-axis, need -x1, -y1, -z1, -x2, -y2, -z2;\n"
				   "[34m                     alignz:    [0m align to Z-axis, need -x1, -y1, -z1, -x2, -y2, -z2;\n"
				   "[34m                     image:     [0m reflect by (0,0,0)\n"
				   "[34m                     imagex:    [0m reflect by Y-Z plane\n"
				   "[34m                     imagey:    [0m reflect by X-Z plane\n"
				   "[34m                     imagez:    [0m reflect by X-Y plane\n"
				   "[31m                 -d [0m degree to be rotated\n"
				   "[31m                 -vx[0m x vector\n"
				   "[31m                 -vy[0m y vector\n"
				   "[31m                 -vz[0m z vector\n"
				   "[31m                 -a1[0m id of atom 1 (0 coordinate center)\n"
				   "[31m                 -a2[0m id of atom 2\n"
				   "[31m                 -x1[0m coord x for point 1\n"
				   "[31m                 -y1[0m coord y for point 1\n"
				   "[31m                 -z1[0m coord z for point 1\n"
				   "[31m                 -x2[0m coord x for point 2\n"
				   "[31m                 -y2[0m coord y for point 2\n"
				   "[31m                 -z2[0m coord z for point 2\n");
			exit(1);
		}
		if (argc < 5) {
			printf("[31mUsage: translate -i [0m input file name (pdb, ac or mol2)\n"
				   "[31m                 -o [0m output file name\n"
				   "[31m                 -r [0m reference file name\n"
				   "[31m                 -f [0m file format\n"
				   "[31m                 -c [0m command (center, translate, rotate1, rotate2, match)\n"
				   "[34m                     check:	  [0m\n"
				   "[34m                     center:    [0m need -a1;\n"
				   "[34m                     translate: [0m need -vx, -vy and -vz;\n"
				   "[34m                     rotate1:   [0m need -a1, -a2 and -d;\n"
				   "[34m                     rotate2:   [0m need -x1, -y1, -z1, -x2, -y2, -z2 and -d;\n"
				   "[34m                     match:     [0m need -r;\n"
				   "[34m                     alignx:    [0m align to X-axis, need -x1, -y1, -z1, -x2, -y2, -z2;\n"
				   "[34m                     aligny:    [0m align to Y-axis, need -x1, -y1, -z1, -x2, -y2, -z2;\n"
				   "[34m                     alignz:    [0m align to Z-axis, need -x1, -y1, -z1, -x2, -y2, -z2;\n"
				   "[34m                     image:     [0m reflect by (0,0,0)\n"
				   "[34m                     imagex:    [0m reflect by Y-Z plane\n"
				   "[34m                     imagey:    [0m reflect by X-Z plane\n"
				   "[34m                     imagez:    [0m reflect by X-Y plane\n"
				   "[31m                 -d [0m degree to be rotated\n"
				   "[31m                 -vx[0m x vector\n"
				   "[31m                 -vy[0m y vector\n"
				   "[31m                 -vz[0m z vector\n"
				   "[31m                 -a1[0m id of atom 1 (0 coordinate center)\n"
				   "[31m                 -a2[0m id of atom 2\n"
				   "[31m                 -x1[0m coord x for point 1\n"
				   "[31m                 -y1[0m coord y for point 1\n"
				   "[31m                 -z1[0m coord z for point 1\n"
				   "[31m                 -x2[0m coord x for point 2\n"
				   "[31m                 -y2[0m coord y for point 2\n"
				   "[31m                 -z2[0m coord z for point 2\n");
			exit(1);
		}
	} else {
		if (argc == 2
			&& (strcmp(argv[1], "-h") == 0
				|| strcmp(argv[1], "-H") == 0)) {
			printf("Usage: translate -i  input file name (pdb, ac or mol2)\n"
				   "             -o  output file name\n"
				   "             -r  reference file name\n"
				   "             -f  file format\n"
				   "             -c  command (center, translate, rotate1, rotate2, rotate3, match)\n"
				   "		     check:\n"	
				   "		     center:     need -a1;\n"	
				   "		     translate:  need -vx, -vy and -vz;\n"	
				   "		     rotate1:    need -a1, -a2 and -d;\n"	
				   "		     rotate2:    need -x1, -y1, -z1, -x2, -y2, -z2 and -d;\n"	
				   "		     match:      need -r;\n"	
				   "                 alignx:     align to X-axis, need -x1, -y1, -z1, -x2, -y2, -z2;\n"
				   "                 aligny:     align to Y-axis, need -x1, -y1, -z1, -x2, -y2, -z2;\n"
				   "                 alignz:     align to Z-axis, need -x1, -y1, -z1, -x2, -y2, -z2;\n"
				   "                 image:      reflect by (0,0,0)\n"
				   "                 imagex:     reflect by Y-Z plane\n"
				   "                 imagey:     reflect by X-Z plane\n"
				   "                 imagez:     reflect by X-Y plane\n"
				   "             -d  degree to be rotated\n"
				   "             -vx x vector\n"
				   "             -vy y vector\n"
				   "             -vz z vector\n"
				   "             -a1 id of atom 1 (0 coordinate center)\n"
				   "             -a2 id of atom 2\n"
				   "             -x1 coord x for point 1\n"
				   "             -y1 coord y for point 1\n"
				   "             -z1 coord z for point 1\n"
				   "             -x2 coord x for point 2\n"
				   "             -y2 coord y for point 2\n"
				   "             -z2 coord z for point 2\n");
			exit(1);
		}
		if (argc < 5) {
                        printf("Usage: translate -i  input file name (pdb, ac or mol2)\n"
                                   "             -o  output file name\n"
                                   "             -r  reference file name\n"
                                   "             -f  file format\n"
                                   "             -c  command (center, translate, rotate1, rotate2, rotate3, match)\n"
				   "		     check:\n"	
                                   "                 center:     need -a1;\n"         
                                   "                 translate: need -vx, -vy and -vz;\n"
                                   "                 rotate1:    need -a1, -a2 and -d;\n"
                                   "                 rotate2:    need -x1, -y1, -z1, -x2, -y2, -z2 and -d;\n"          
                                   "                 match:      need -r;\n"  
				   "                 alignx:     align to X-axis, need -x1, -y1, -z1, -x2, -y2, -z2;\n"
				   "                 aligny:     align to Y-axis, need -x1, -y1, -z1, -x2, -y2, -z2;\n"
				   "                 alignz:     align to Z-axis, need -x1, -y1, -z1, -x2, -y2, -z2;\n"
				   "                 image:      reflect by (0,0,0)\n"
				   "                 imagex:     reflect by Y-Z plane\n"
				   "                 imagey:     reflect by X-Z plane\n"
				   "                 imagez:     reflect by X-Y plane\n"
				   "             -d  degree to be rotated\n"
                                   "             -vx x vector\n"
                                   "             -vy y vector\n"
                                   "             -vz z vector\n"
                                   "             -a1 id of atom 1 (0 coordinate center)\n"
                                   "             -a2 id of atom 2\n"
                                   "             -x1 coord x for point 1\n"
                                   "             -y1 coord y for point 1\n"
                                   "             -z1 coord z for point 1\n"
                                   "             -x2 coord x for point 2\n"
                                   "             -y2 coord y for point 2\n"
                                   "             -z2 coord z for point 2\n");
                        exit(1);
		}
	}
	format = -1;
	for (i = 1; i < argc; i += 2) {
		if (strcmp(argv[i], "-c")== 0) { 
			strcpy(command, argv[i + 1]);
			cflag = 1;
		}
		if (strcmp(argv[i], "-i") == 0) {
			strcpy(ifilename, argv[i + 1]);
			iflag = 1;
		}
		if (strcmp(argv[i], "-o") == 0) {
			strcpy(ofilename, argv[i + 1]);
			oflag = 1;
		}
		if (strcmp(argv[i], "-r") == 0) {
			strcpy(rfilename, argv[i + 1]);
			rflag = 1;
		}
		if (strcmp(argv[i], "-f") == 0) {
			if (strcmp(argv[i + 1], "ac") == 0)
				format = 0;
			if (strcmp(argv[i + 1], "pdb") == 0)
				format = 1;
			if (strcmp(argv[i + 1], "mol2") == 0)
				format = 2;
		}
		if (strcmp(argv[i], "-a1") == 0)  
			at1 = atoi(argv[i+1]);
		if (strcmp(argv[i], "-a2") == 0)  
			at2 = atoi(argv[i+1]);
		if (strcmp(argv[i], "-a3") == 0)  
			at3 = atoi(argv[i+1]);
		if (strcmp(argv[i], "-d") == 0)  
			degree = atof(argv[i+1]);
		if (strcmp(argv[i], "-vx") == 0)  
			vectx = atof(argv[i+1]);
		if (strcmp(argv[i], "-vy") == 0)  
			vecty = atof(argv[i+1]);
		if (strcmp(argv[i], "-vz") == 0)  
			vectz = atof(argv[i+1]);
		if (strcmp(argv[i], "-x1") == 0)  
			coord_x1 = atof(argv[i+1]);
		if (strcmp(argv[i], "-y1") == 0)  
			coord_y1 = atof(argv[i+1]);
		if (strcmp(argv[i], "-z1") == 0)  
			coord_z1 = atof(argv[i+1]);
		if (strcmp(argv[i], "-x2") == 0)  
			coord_x2 = atof(argv[i+1]);
		if (strcmp(argv[i], "-y2") == 0)  
			coord_y2 = atof(argv[i+1]);
		if (strcmp(argv[i], "-z2") == 0)  
			coord_z2 = atof(argv[i+1]);
	}

       	if (format != 0 && format != 1 && format != 2) {
		fprintf(stdout, "\nError, no or wrong format, must be 'ac', 'pdb' or 'mol2'\n");
		exit(0);
	}

	if(cflag == 0) {
		fprintf(stdout, "\nNo command specified, exit\n");
		exit(1);
	}
	if(iflag == 0) {
		fprintf(stdout, "\nNo input file specified, exit\n");
		exit(1);
	}
	/*read in prep or ac file */
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

	if (format == 0)
		overflow_flag =
			rac(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo);
	if (format == 1)
		overflow_flag =
			 rpdb(ifilename, &atomnum, atom, cinfo, minfo, 0); 	
	if (format == 2)
		overflow_flag =
			rmol2(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo, 1);
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
        	if (format == 0)
                	overflow_flag =
                        	rac(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo);
        	if (format == 1)
                	overflow_flag =
                         	rpdb(ifilename, &atomnum, atom, cinfo, minfo, 0); 
        	if (format == 2)
                	overflow_flag =
                        	rmol2(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo, 1);
	}
	if(strcmp(command, "check") == 0) {
		for(i=0;i<atomnum;i++) {
			if(atom[i].x < minx) minx = atom[i].x;
			if(atom[i].x > maxx) maxx = atom[i].x;
			if(atom[i].y < miny) miny = atom[i].y;
			if(atom[i].y > maxy) maxy = atom[i].y;
			if(atom[i].z < minz) minz = atom[i].z;
			if(atom[i].z > maxz) maxz = atom[i].z;
			sum_coordx += atom[i].x;
			sum_coordy += atom[i].y;
			sum_coordz += atom[i].z;
		}
		coordx = sum_coordx/atomnum;
		coordy = sum_coordy/atomnum;
		coordz = sum_coordz/atomnum;
		printf("\nThe molecule center is at: X= %9.4lf, Y= %9.4lf, Z = %9.4lf", coordx, coordy, coordz);
		printf("\nThe molecule dimensions are listed as follows");
		printf("\n MIN X= %9.4lf, MAX X= %9.4lf", minx, maxx); 
		printf("\n MIN Y= %9.4lf, MAX Y= %9.4lf", miny, maxy); 
		printf("\n MIN Z= %9.4lf, MAX Z= %9.4lf\n", minz, maxz); 
	}	

	if(strcmp(command, "center") == 0) {
		if(at1 <-99990) {
			printf("\nNo -a1 specified\n");
			exit(1);
		}
		if(at1 == 0) {
			for(i=0;i<atomnum;i++) {
				sum_coordx += atom[i].x;
				sum_coordy += atom[i].y;
				sum_coordz += atom[i].z;
			}
			coordx = sum_coordx/atomnum;
			coordy = sum_coordy/atomnum;
			coordz = sum_coordz/atomnum;
		}
		else {
			coordx = atom[at1-1].x;
			coordy = atom[at1-1].y;
			coordz = atom[at1-1].z;
		}
		for(i=0;i<atomnum;i++) {
			 atom[i].x-=coordx;
			 atom[i].y-=coordy;
			 atom[i].z-=coordz;
		}
		printf("\nThe molecule is translated: Vector X= %9.4lf, Vector Y= %9.4lf, Vector Z = %9.4lf\n", coordx, coordy, coordz);
	}	
	if(strcmp(command, "translate") == 0) {
		if(vectx <-99990) {
			printf("\nNo -vx specified\n");
			exit(1);
		}
		if(vecty <-99990) {
			printf("\nNo -vy specified\n");
			exit(1);
		}
		if(vectz <-99990) {
			printf("\nNo -vz specified\n");
			exit(1);
		}
		coordx =vectx; 
		coordy = vecty;
		coordz = vectz;
		for(i=0;i<atomnum;i++) {
			 atom[i].x-=coordx;
			 atom[i].y-=coordy;
			 atom[i].z-=coordz;
		}
		printf("\nThe molecule is translated: Vector X= %9.4lf, Vector Y= %9.4lf, Vector Z = %9.4lf\n", coordx, coordy, coordz);
	}	
	if(strcmp(command, "rotate1") == 0) {
		if(at1 <-99990) {
			printf("\nNo -a1 specified\n");
			exit(1);
		}
		if(at2 <-99990) {
			printf("\nNo -a2 specified\n");
			exit(1);
		}
		if(degree <-99990) {
			printf("\nNo -d specified\n");
			exit(1);
		}
		omegarotate(atom[at1-1], atom[at2-1], &w1, &w2);
		omegarotate2(atom, atomnum, atom[at2-1], degree, w1, w2);
	}	
	if(strcmp(command, "rotate2") == 0) {
		if(coord_x1 <-99990) {
			printf("\nNo -x1 specified\n");
			exit(1);
		}
		if(coord_y1 <-99990) {
			printf("\nNo -y1 specified\n");
			exit(1);
		}
		if(coord_z1 <-99990) {
			printf("\nNo -z1 specified\n");
			exit(1);
		}
		if(coord_x2 <-99990) {
			printf("\nNo -x2 specified\n");
			exit(1);
		}
		if(coord_y2 <-99990) {
			printf("\nNo -y2 specified\n");
			exit(1);
		}
		if(coord_z2 <-99990) {
			printf("\nNo -z2 specified\n");
			exit(1);
		}
		if(degree <-99990) {
			printf("\nNo -d specified\n");
			exit(1);
		}
		tmpatom1.x = coord_x1;
		tmpatom1.y = coord_y1;
		tmpatom1.z = coord_z1;
		tmpatom2.x = coord_x2;
		tmpatom2.y = coord_y2;
		tmpatom2.z = coord_z2;
		omegarotate(tmpatom1, tmpatom2, &w1, &w2);
		omegarotate2(atom, atomnum, tmpatom2, degree, w1, w2);
	}	
	if(strcmp(command, "match") == 0) {
		if(rflag == 0) {
			printf("\nNo -r specified\n");
			exit(1);
		}	
                ref_atom = (ATOM *) malloc(sizeof(ATOM) * cinfo.maxatom);
                if (ref_atom == NULL) {
                        fprintf(stdout, "memory allocation error for *refatom\n");
                        exit(1);
                }
                ref_bond = (BOND *) malloc(sizeof(BOND) * cinfo.maxbond);
                if (ref_bond == NULL) {
                        fprintf(stdout, "memory allocation error for *refbond\n");
                        exit(1);
                }
                int i;
                for (i = 0; i < cinfo.maxbond; ++i) {
                        ref_bond[i].jflag = -1; /* bond type has not been assigned */
                }
                if (format == 0)
                        overflow_flag =
                                rac(rfilename, &atomnum, ref_atom, &bondnum, ref_bond, &cinfo, &minfo);
                if (format == 1)
                        overflow_flag =
                                rpdb(rfilename, &atomnum, ref_atom, cinfo, minfo, 0);
                if (format == 2)
                        overflow_flag =
                                rmol2(rfilename, &atomnum, ref_atom, &bondnum, ref_bond, &cinfo, &minfo, 1);
		rmsd = lsfit(atom, ref_atom, atom, atomnum); 
		printf("\nThe rmsd is %9.4lf\n", rmsd);	
	}
        if(strcmp(command, "alignx") == 0) {
                if(coord_x1 <-99990) {
                        printf("\nNo -x1 specified\n");
                        exit(1);
                }
                if(coord_y1 <-99990) {
                        printf("\nNo -y1 specified\n");
                        exit(1);
                }
                if(coord_z1 <-99990) {
                        printf("\nNo -z1 specified\n");
                        exit(1);
                }
                if(coord_x2 <-99990) {
                        printf("\nNo -x2 specified\n");
                        exit(1);
                }
                if(coord_y2 <-99990) {
                        printf("\nNo -y2 specified\n");
                        exit(1);
                }
                if(coord_z2 <-99990) {
                        printf("\nNo -z2 specified\n");
                        exit(1);
                }
		
		if(coord_x1 > coord_x2) {
                	tmpatom1.x = coord_x1;
                	tmpatom1.y = coord_y1;
                	tmpatom1.z = coord_z1;
                	tmpatom2.x = coord_x2;
                	tmpatom2.y = coord_y2;
                	tmpatom2.z = coord_z2;
		}
		else {
                	tmpatom1.x = coord_x2;
                	tmpatom1.y = coord_y2;
                	tmpatom1.z = coord_z2;
                	tmpatom2.x = coord_x1;
                	tmpatom2.y = coord_y1;
                	tmpatom2.z = coord_z1;
		}
		alignx(tmpatom1, tmpatom2, atom, atomnum);
        }

        if(strcmp(command, "aligny") == 0) {
                if(coord_x1 <-99990) {
                        printf("\nNo -x1 specified\n");
                        exit(1);
                }
                if(coord_y1 <-99990) {
                        printf("\nNo -y1 specified\n");
                        exit(1);
                }
                if(coord_z1 <-99990) {
                        printf("\nNo -z1 specified\n");
                        exit(1);
                }
                if(coord_x2 <-99990) {
                        printf("\nNo -x2 specified\n");
                        exit(1);
                }
                if(coord_y2 <-99990) {
                        printf("\nNo -y2 specified\n");
                        exit(1);
                }
                if(coord_z2 <-99990) {
                        printf("\nNo -z2 specified\n");
                        exit(1);
                }

                if(coord_y1 > coord_y2) {
                        tmpatom1.x = coord_x1;
                        tmpatom1.y = coord_y1;
                        tmpatom1.z = coord_z1;
                        tmpatom2.x = coord_x2;
                        tmpatom2.y = coord_y2;
                        tmpatom2.z = coord_z2;
                }
                else {
                        tmpatom1.x = coord_x2;
                        tmpatom1.y = coord_y2;
                        tmpatom1.z = coord_z2;
                        tmpatom2.x = coord_x1;
                        tmpatom2.y = coord_y1;
                        tmpatom2.z = coord_z1;
                }
                aligny(tmpatom1, tmpatom2, atom, atomnum);
        }

        if(strcmp(command, "alignz") == 0) {
                if(coord_x1 <-99990) {
                        printf("\nNo -x1 specified\n");
                        exit(1);
                }
                if(coord_y1 <-99990) {
                        printf("\nNo -y1 specified\n");
                        exit(1);
                }
                if(coord_z1 <-99990) {
                        printf("\nNo -z1 specified\n");
                        exit(1);
                }
                if(coord_x2 <-99990) {
                        printf("\nNo -x2 specified\n");
                        exit(1);
                }
                if(coord_y2 <-99990) {
                        printf("\nNo -y2 specified\n");
                        exit(1);
                }
                if(coord_z2 <-99990) {
                        printf("\nNo -z2 specified\n");
                        exit(1);
                }

                if(coord_z1 > coord_z2) {
                        tmpatom1.x = coord_x1;
                        tmpatom1.y = coord_y1;
                        tmpatom1.z = coord_z1;
                        tmpatom2.x = coord_x2;
                        tmpatom2.y = coord_y2;
                        tmpatom2.z = coord_z2;
                }
                else {
                        tmpatom1.x = coord_x2;
                        tmpatom1.y = coord_y2;
                        tmpatom1.z = coord_z2;
                        tmpatom2.x = coord_x1;
                        tmpatom2.y = coord_y1;
                        tmpatom2.z = coord_z1;
                }
                alignz(tmpatom1, tmpatom2, atom, atomnum);
        }
        if(strcmp(command, "image") == 0) {
                for(i=0;i<atomnum;i++) {
                         atom[i].x*=-1;
                         atom[i].y*=-1;
                         atom[i].z*=-1;
                }
                printf("\nThe molecule is reflected by (0,0,0)\n");
        }
 
        if(strcmp(command, "imageX") == 0 || strcmp(command, "imagex") == 0) {
                for(i=0;i<atomnum;i++) {
                         atom[i].x*=-1;
                }
                printf("\nThe molecule is reflected by Y-Z plane\n");
        }

        if(strcmp(command, "imagey") == 0 || strcmp(command, "imagey") == 0) {
                for(i=0;i<atomnum;i++) {
                         atom[i].y*=-1;
                }
                printf("\nThe molecule is reflected by X-Z plane\n");
        }

        if(strcmp(command, "imagez") == 0 || strcmp(command, "imagez") == 0) {
                for(i=0;i<atomnum;i++) {
                         atom[i].z*=-1;
                }
                printf("\nThe molecule is reflected by X-Y plane\n");
        }

	if(oflag == 1) {
		if(format == 0) wac(ofilename, atomnum, atom, bondnum, bond, cinfo, minfo); 
 		if(format == 1) wpdb(ofilename, atomnum, atom);
               	if(format == 2) wmol2(ofilename, atomnum, atom, bondnum, bond, arom, cinfo, minfo); 
	}
	return (0);
}
