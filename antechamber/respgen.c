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
*  Octomber, 2001                                                      *
************************************************************************
*/
char *amberhome;

# include "common.h"
# include "define.h"
# include "atom.h"
# include "utility.c"
# include "common.c"
# include "ac.c"
# include "equatom.c"
# define MAX_RESP_PARM 100
# define debug 0

typedef struct {
        char atomtype[10];
        double atqwt;
        double refcharge;
} RESPPARM;
RESPPARM respparm[MAX_RESP_PARM];
int parmnum;

ATOM *atom;
BOND *bond;
double *atqwt;
double *refcharge;
int atomnum = 0;
int bondnum = 0;
CONTROLINFO cinfo;
MOLINFO minfo;

int *equ_atom_id;
int max_path_length = -1;

int i, j, k, l;
FILE *fpin;
FILE *fpout;
FILE *fpcharge;
FILE *fpaddinfo;
char line[MAXCHAR];
char ifilename[MAXCHAR];
char ofilename[MAXCHAR];
char afilename[MAXCHAR];
char dfilename[MAXCHAR];
char pfilename[MAXCHAR];

int overflow_flag = 0;
int method = -1;
double charge = 0.0;
double weight = 0.0;

int iaddinfo = 0;
double *pcharge; /* predefined partial charge */
int numcharge = 0;
int *icharge; /*flag of predefined partial charge */
int nconf = 1;

int iequ = 1;

void readparm(char *filename) {
char line[MAXCHAR];
int i;
FILE *fp;
	if ((fp = fopen(filename, "r")) == NULL) {
		printf("\n Cannot open file %s to read in readparm(), exit", filename);
		return;
	}
	parmnum = 0;
	for(;;) {
		if (fgets(line, MAXCHAR, fp) == NULL) break;
        	if (strncmp("ATD", line, 3) == 0) {
			sscanf(&line[3], "%s%lf%lf", respparm[parmnum].atomtype, &respparm[parmnum].atqwt, &respparm[parmnum].refcharge);
			parmnum++;
			if(parmnum >= MAX_RESP_PARM) {
				fprintf(stdout, "Too many resp parameters defined (%d), increase MAX_RESP_PARM (%d) and recompile respgen\n",
					parmnum, MAX_RESP_PARM);
				exit(1);
			}
		}
	}
	if(debug == 1)
		for(i=0; i<parmnum; i++)
			printf("PARM	%d	%s	%9.4lf	%9.4lf\n", i+1, respparm[i].atomtype, respparm[i].atqwt, respparm[i].refcharge);
	fclose(fp);
}

void assignparm(char *filename) {
int i;
int suc;
int count;
char atomtype[20];
char line[MAXCHAR];
FILE *fp;
	if ((fp = fopen(filename, "r")) == NULL) {
		printf("\n Cannot open file %s to read in assignparm(), exit", filename);
		return;
	}
	count = 0;
	for(;;) {
		if (fgets(line, MAXCHAR, fp) == NULL) break;
        	if (strncmp("ATOM", line, 4) == 0) {
			sscanf(&line[65], "%s", atomtype);
			suc = 0;
			for(i=0;i<parmnum;i++) 
				if(strcmp(respparm[i].atomtype, atomtype) == 0) { 
					atqwt[count] = respparm[i].atqwt;
					refcharge[count] = respparm[i].refcharge;
					suc = 1;
					break;
				}
			if(suc == 0) {
				atqwt[count] = 0;
				refcharge[count] = 0;
			}
			count++;
		}
	}
fclose(fp);
}

void readinfo() {
int i,j;
int tmpint;
int num;
double tmpf;
numcharge = 0;
for (;;) {
	if (fgets(line, MAXCHAR, fpaddinfo) == NULL) break;
        if (strncmp("CHARGE", line, 6) == 0) {
		sscanf(&line[6], "%lf%d", &tmpf, &tmpint); 
		if(tmpint < 1) {
			fprintf(stdout, "\nAtom ID should not smaller than 1, exit");
			exit(1);	
		}
		pcharge[tmpint-1] = tmpf;
		icharge[tmpint-1] = 1;
		numcharge ++;
        }
}
num = 0;
for (i = 0; i < nconf; i++)
	for (j = 0; j < atomnum; j++) {
        	fprintf(fpcharge, "%10.6lf", pcharge[j]);
        	num++;
        	if (num == 8) {
                	num = 0;
                	fprintf(fpcharge, "\n");
        	}
	}

}

void respin(int method)
{
	int i, j;
	int cindex;
	int *chindex;
	int times = 0;
	int tmpint;
	int count;
	int tcount;
	int tmpnum;

	double tmpcharge;
	FILE *fpout;

	chindex = (int *) malloc(sizeof(int) * (atomnum + 10));
        if (chindex == NULL) {
                fprintf(stdout, "memory allocation error for *chindex in respin()\n");
                exit(1);
        }
	
	if ((fpout = fopen(ofilename, "w")) == NULL) {
		printf("\n Cannot open file %s to write in respin(), exit", ofilename);
		return;
	}
	for (i = 0; i < atomnum; i++)
		chindex[i] = -1;
	if(method != 3) {
		for (i = 0; i < atomnum; i++) {
			cindex = 0;
			if (atom[i].name[0] == 'C')
				for (j = 0; j < 6; j++)
					if (atom[i].con[j] >= 0)
						if (atom[atom[i].con[j]].atomicnum == 1)
							cindex++;
			if (cindex >= 2) {
				chindex[i] = 1;
				for (j = 0; j < 6; j++)
					if (atom[i].con[j] >= 0)
						if (atom[atom[i].con[j]].atomicnum == 1) 
							chindex[atom[i].con[j]] = 1;
			}
		}
	}
/* for(i=0;i<atomnum;i++)
 printf("\n %5d%5s%5d%5d", i+1, atom[i].name, chindex[i], equ_atom_id[i]);
 */
	fprintf(fpout, "%s\n\n", "Resp charges for organic molecule");
	fprintf(fpout, "%s\n\n", " &cntrl");
	if (method == 3) 
		fprintf(fpout, "%s\n", " irun = 1,");
	fprintf(fpout, " nmol = %d,\n", nconf);
	fprintf(fpout, "%s\n", " ihfree = 1,");
	fprintf(fpout, "%s\n", " ioutopt = 1,");
	if (method == 0) {
		fprintf(fpout, "%s\n", " iqopt = 2,");
		fprintf(fpout, "%s\n", " irstrnt = 2,");
	}
	if (method == 1) {
		if(numcharge > 0) {
			fprintf(fpout, "%s\n", " iqopt = 2,");
			fprintf(fpout, " qwt =%8.5lf,\n", weight);
		}
		else 
			fprintf(fpout, " qwt =%8.5lf,\n", weight);
	}
	if (method == 2) {
		fprintf(fpout, "%s\n", " iqopt = 2,");
		fprintf(fpout, " qwt = %8.5lf,\n", weight);
	}
	if (method == 3) {
		if(numcharge > 0) {
			fprintf(fpout, "%s\n", " iqopt = 2,");
			fprintf(fpout, " qwt =%8.5lf,\n", weight);
		}
		else 
			fprintf(fpout, " qwt =%8.5lf,\n", weight);
	}
        if (method == 4) {
                fprintf(fpout, "%s\n", " qwt = 0.0005,");
                fprintf(fpout, "%s\n", " ipol = 4,");
                fprintf(fpout, "%s\n", " n12=0, n13=0, n14=1, scee=1.0,");
                fprintf(fpout, "%s\n", " iscreen=0,");
                fprintf(fpout, "%s\n", " ipol_iter=50,");
        }
        if (method == 5) {
                fprintf(fpout, "%s\n", " iqopt = 2,");
                fprintf(fpout, "%s\n", " qwt = 0.001,");
        }
	if (method == 6) {
		fprintf(fpout, " qwt = %8.5lf,\n", weight);
	}
        if (method == 7) {
		fprintf(fpout, "%s\n", " iqopt = 2,");
		fprintf(fpout, "%s\n", " irstrnt = 2,");
                fprintf(fpout, "%s\n", " ipol = 4,");
                fprintf(fpout, "%s\n", " n12=0, n13=0, n14=1, scee=1.0,");
                fprintf(fpout, "%s\n", " iscreen=0,");
                fprintf(fpout, "%s\n", " ipol_iter=50,");
        }
	if (method == -1) {
		fprintf(fpout, "%s\n", " idqrtrnt = 0,");
		fprintf(fpout, "%s\n", " dmscale = 1,");
	}
	fprintf(fpout, "\n%s\n", " &end");
	while(times < nconf) {
		fprintf(fpout, "%s\n", "    1.0");
		fprintf(fpout, "%s\n", "Resp charges for organic molecule");
		fprintf(fpout, "%5d%5d\n", (int) charge, atomnum);
		if (method == -1) {
			for (i = 0; i < atomnum; i++) {
				if(iaddinfo == 1 && icharge[i] == 1)
					fprintf(fpout, "%5d%5d", atom[i].atomicnum, -99);
				else {
					if (equ_atom_id[i] == -1)
						fprintf(fpout, "%5d%5d", atom[i].atomicnum, 0);
					else
						fprintf(fpout, "%5d%5d", atom[i].atomicnum,
								equ_atom_id[i] + 1);
				}
				if (chindex[i] == 1 && atom[i].atomicnum == 6)
					fprintf(fpout, "%8.4lf\n", 0.001);
				else
					fprintf(fpout, "\n");
			}
		}
		if (method == 0 || method == 7) { 
			for (i = 0; i < atomnum; i++) 
				fprintf(fpout, "%5d%5d\n", atom[i].atomicnum, -99);
		}
		if (method == 1) {
			for (i = 0; i < atomnum; i++) {
				if(iaddinfo == 1 && icharge[i] == 1)
					fprintf(fpout, "%5d%5d\n", atom[i].atomicnum, -99);
				else {
					if (equ_atom_id[i] == -1)
						fprintf(fpout, "%5d%5d\n", atom[i].atomicnum, 0);
					else if (chindex[i] == 1)
						fprintf(fpout, "%5d%5d\n", atom[i].atomicnum, 0);
					else
						fprintf(fpout, "%5d%5d\n", atom[i].atomicnum,
								equ_atom_id[i] + 1);
				}
			}
		}

		if (method == 2) {
			for (i = 0; i < atomnum; i++) {
				if(iaddinfo == 1 && icharge[i] == 1)
					fprintf(fpout, "%5d%5d\n", atom[i].atomicnum, -99);
				else {
					if (chindex[i] != 1)
						fprintf(fpout, "%5d%5d\n", atom[i].atomicnum, -99);
					else if (equ_atom_id[i] == -1)
						fprintf(fpout, "%5d%5d\n", atom[i].atomicnum, 0);
					else
						fprintf(fpout, "%5d%5d\n", atom[i].atomicnum,
								equ_atom_id[i] + 1);
				}
			}
		}

		if (method == 3) {
			for (i = 0; i < atomnum; i++) {
				if(iaddinfo == 1 && icharge[i] == 1)
					fprintf(fpout, "%5d%5d%12.4lf%12.4lf\n", atom[i].atomicnum, -99, atqwt[i],refcharge[i]);
				else {
					if (equ_atom_id[i] == -1)
						fprintf(fpout, "%5d%5d%12.4lf%12.4lf\n",atom[i].atomicnum, 0, atqwt[i], refcharge[i]);
					else if (chindex[i] == 1)
						fprintf(fpout, "%5d%5d%12.4lf%12.4lf\n",atom[i].atomicnum, 0, atqwt[i], refcharge[i]);
					else
						fprintf(fpout, "%5d%5d%12.4lf%12.4lf\n",atom[i].atomicnum, equ_atom_id[i] + 1, 
											atqwt[i], refcharge[i]);
				}
			}
		}
                if (method == 4) {
                        for (i = 0; i < atomnum; i++)
                                fprintf(fpout, "%5d%5d\n", atom[i].atomicnum, 0);
                }
                if (method == 5) {
                        for (i = 0; i < atomnum; i++) {
                                if(iaddinfo == 1 && icharge[i] == 1)
                                        fprintf(fpout, "%5d%5d\n", atom[i].atomicnum, -99);
                                else {
                                        if (equ_atom_id[i] == -1)
                                                fprintf(fpout, "%5d%5d\n", atom[i].atomicnum, 0);
                                        else
                                                fprintf(fpout, "%5d%5d\n", atom[i].atomicnum, equ_atom_id[i] + 1);
                                }
                        }
                }
/*
//	according to my study, it is good to apply no constrain (except the net charge of the whole molecule) for the first step
// 	and use a kind of one-step fitting in the second step
		if (method == 4) {
			for (i = 0; i < atomnum; i++) {
				if(iaddinfo == 1 && icharge[i] == 1)
					fprintf(fpout, "%5d%5d\n", atom[i].atomicnum, -99);
				else {
					if (equ_atom_id[i] == -1)
						fprintf(fpout, "%5d%5d\n", atom[i].atomicnum, 0);
					else if (chindex[i] == 1)
						fprintf(fpout, "%5d%5d\n", atom[i].atomicnum, 0);
					else
						fprintf(fpout, "%5d%5d\n", atom[i].atomicnum,
								equ_atom_id[i] + 1);
				}
			}
		}

		if (method == 5) {
			for (i = 0; i < atomnum; i++) {
				if(iaddinfo == 1 && icharge[i] == 1)
					fprintf(fpout, "%5d%5d\n", atom[i].atomicnum, -99);
				else {
					if (chindex[i] != 1)
						fprintf(fpout, "%5d%5d\n", atom[i].atomicnum, -99);
					else if (equ_atom_id[i] == -1)
						fprintf(fpout, "%5d%5d\n", atom[i].atomicnum, 0);
					else
						fprintf(fpout, "%5d%5d\n", atom[i].atomicnum,
								equ_atom_id[i] + 1);
				}
			}
		}
*/
                if (method == 6) {
                        for (i = 0; i < atomnum; i++) {
                                if(iaddinfo == 1 && icharge[i] == 1)
                                        fprintf(fpout, "%5d%5d\n", atom[i].atomicnum, -99);
                                else {
                                        if (equ_atom_id[i] == -1)
                                                fprintf(fpout, "%5d%5d\n", atom[i].atomicnum, 0);
					else
                                                fprintf(fpout, "%5d%5d\n", atom[i].atomicnum, equ_atom_id[i] + 1);
                                }
                        }
                }
		if(nconf > 1) fprintf(fpout, "\n");
		times++;
   	}
	if(iaddinfo == 1) {
/* group */
		rewind(fpaddinfo);
		tmpnum = -1;
		for (;;) {
        		if (fgets(line, MAXCHAR, fpaddinfo) == NULL) break;
        		if (strncmp("GROUP", line, 5) == 0) {
				if(tmpnum > 0) {
					if(tmpnum != tcount) {
						fprintf(stdout,"The number of ATOM lines does not equal to the number defined in the GROUP line\n"); 
						exit(1);
					}
					else 
						fprintf(fpout, "\n"); 
				} 
                		sscanf(&line[5], "%d%lf", &tmpnum, &tmpcharge);
				fprintf(fpout, "%5d%8.3lf\n", tmpnum, tmpcharge); 
				count = 0;
				tcount = 0;
				continue;
			}
        		if (strncmp("ATOM", line, 4) == 0) {
                		sscanf(&line[4], "%d", &tmpint);
				fprintf(fpout, "%5d%5d", 1, tmpint); 
				count ++;
				tcount ++;
				if(count >= 8 && tcount != tmpnum) {
					count = 0;
					fprintf(fpout, "\n");
				}	
			}
        	}

		if(tmpnum > 0) {
			if(tmpnum != tcount) {
				fprintf(stdout,"The number of ATOM lines does not equal to the number defined in the GROUP line\n"); 
				exit(1);
			}
			else
				fprintf(fpout, "\n\n"); 
		} 
/* equaliation*/
		for(i=1; i< nconf; i++)
			for(j=0; j< atomnum ; j++) {
				fprintf(fpout, "%5d\n", 2); 
				fprintf(fpout, "%5d%5d%5d%5d\n", 1, j+1, i+1, j+1); 
			}
	}
	if(iaddinfo == 0 && nconf > 1) {
		fprintf(fpout, "\n");
		for(i=1; i< nconf; i++)
			for(j=0; j< atomnum ; j++) {
				fprintf(fpout, "%5d\n", 2); 
				fprintf(fpout, "%5d%5d%5d%5d\n", 1, j+1, i+1, j+1); 
			}
	}
	free(chindex);
	fprintf(fpout, "\n");
	fclose(fpout);
}

int main(int argc, char *argv[])
{
	int i;
        int status = 0;
	char command[MAXCHAR];
	size_t copied_size;  

    	amberhome = (char *) getenv("AMBERHOME");
    	if( amberhome == NULL ){
       		fprintf( stdout, "AMBERHOME is not set!\n" );
       		exit(1);
    	}

	if (strcmp(COLORTEXT, "YES") == 0 || strcmp(COLORTEXT, "yes") == 0) {
		if (argc == 2
			&& (strcmp(argv[1], "-h") == 0
				|| strcmp(argv[1], "-H") == 0)) {
			printf("[31mUsage: respgen -i[0m input file name(ac)\n"
				   "[31m               -o[0m output file name\n"
				   "[31m               -l[0m maximum path length (default is -1, only recommand if the default setting\n"
				   "                  taking long time and/or cause core dump. If applied, a value of 8 to 10 should good)\n"
				   "[31m               -f[0m output file format (resp1 or resp2) \n"
				   "[32m                  resp0[0m - evaluation the current charges \n"
				   "[32m                  resp1[0m - first stage resp fitting \n"
				   "[32m                  resp2[0m - second stage resp fitting\n"
				   "[32m                  iresp0[0m - evaluation the current charges for polarizable model\n"
				   "[32m                  iresp1[0m- first stage of i_resp fitting \n"
				   "[32m                  iresp2[0m- second stage of i_resp fitting\n"
				   "[32m                  resp3[0m - one-stage resp fitting\n"
				   "[32m                  resp4[0m - calculating ESP from point charges\n"
				   "[31m               -e[0m equalizing atomic charge, default is 1\n"
				   "[32m                  0[0m not use \n"
				   "[32m                  1[0m by atomic paths\n"
				   "[32m                  2[0m by atomic paths and structural information, i.e. E/Z confirgurations\n"
				   "[31m               -a[0m additional input data (predefined charges, atom groups etc).)\n"
				   "[31m               -n[0m number of conformations (default is 1)\n"
				   "[31m               -w[0m weight of charge constraint, in default, 0.0005 for resp1 and 0.001 fore resp2\n");
			exit(1);
		}
		if (argc != 7 && argc != 9 && argc != 11 && argc != 13 && argc != 15 && argc != 17) {
			printf("[31mUsage: respgen -i[0m input file name(ac)\n"
				   "[31m               -o[0m output file name\n"
				   "[31m               -l[0m maximum path length (default is -1, only recommand if the default setting\n"
				   "                  taking long time and/or cause core dump. If applied, a value of 8 to 10 should good)\n"
				   "[31m               -f[0m output file format (resp1 or resp2) \n"
				   "[32m                  resp0[0m - evaluation the current charges \n"
				   "[32m                  resp1[0m - first stage resp fitting \n"
				   "[32m                  resp2[0m - second stage resp fitting\n"
				   "[32m                  iresp0[0m - evaluation the current charges for polarizable model\n"
				   "[32m                  iresp1[0m- first stage of i_resp fitting \n"
				   "[32m                  iresp2[0m- second stage of i_resp fitting\n"
				   "[32m                  resp3[0m - one-stage resp fitting\n"
				   "[32m                  resp4[0m - calculating ESP from point charges\n"
				   "[31m               -e[0m equalizing atomic charge, default is 1\n"
				   "[32m                  0[0m not use \n"
				   "[32m                  1[0m by atomic paths\n"
				   "[32m                  2[0m by atomic paths and structural information, i.e. E/Z confirgurations\n"
				   "[31m               -a[0m additional input data (predefined charges, atom groups etc).)\n"
				   "[31m               -n[0m number of conformations (default is 1)\n"
				   "[31m               -w[0m weight of charge constraint, in default, 0.0005 for resp1 and 0.001 fore resp2\n");
			exit(1);
		}
	} else {
		if (argc == 2
			&& (strcmp(argv[1], "-h") == 0
				|| strcmp(argv[1], "-H") == 0)) {
			printf("Usage respgen -i input file name(ac)\n");
			printf("              -o output file name\n");
			printf("	      -l maximum path length (default is -1, only recommand if the default setting\n");
			printf("                 taking long time and/or cause core dump. If applied, a value of 8 to 10 should good)\n");
			printf("              -f output file format (resp1 or resp2)\n");
		  	printf("	         resp0 - evaluation the current charges \n");
			printf("                 resp1 - first stage resp fitting\n");
			printf("                 resp2 - second stage resp fitting \n");
		  	printf("	         iresp0 - evaluation the current charges for polarizable model \n");
			printf("                 iresp1- first stage of i_resp fitting\n");
			printf("                 iresp2- second stage of i_resp fitting \n");
		        printf("                 resp3 - one-stage resp fitting\n");
			printf("                 resp4 - calculating ESP from point charges \n");
			printf("	      -e equalizing atomic charge, default is 1\n"); 
			printf("                  0 not use \n");
		 	printf("                  1 by atomic paths\n");
			printf("                  2 by atomic paths and structural information, i.e. E/Z confirgurations\n");
		        printf("              -a additional input data (predefined charges, atom groups etc).)\n");
		        printf("              -n number of conformations (default is 1)\n");
			printf("	      -w weight of charge constraint, in default, 0.0005 for resp1 and 0.001 fore resp2\n");
			exit(1);
		}
		if (argc != 7 && argc != 9 && argc != 11 && argc != 13 && argc != 15 && argc != 17) {
			printf("Usage respgen -i input file name(ac)\n");
			printf("              -o output file name\n");
			printf("	      -l maximum path length (default is -1, only recommand if the default setting\n");
			printf("                 taking long time and/or cause core dump. If applied, a value of 8 to 10 should good)\n");
			printf("              -f output file format (resp1 or resp2)\n");
		  	printf("	         resp0 - evaluation the current charges \n");
			printf("                 resp1 - first stage resp fitting\n");
			printf("                 resp2 - second stage resp fitting \n");
		  	printf("	         iresp0 - evaluation the current charges for polarizable model \n");
			printf("                 iresp1- first stage of i_resp fitting\n");
			printf("                 iresp2- second stage of i_resp fitting \n");
		        printf("                 resp3 - one-stage resp fitting\n");
			printf("                 resp4 - calculating ESP from point charges \n");
			printf("	      -e equalizing atomic charge, default is 1\n"); 
			printf("                  0 not use \n");
		 	printf("                  1 by atomic paths\n");
			printf("                  2 by atomic paths and structural information, i.e. E/Z confirgurations\n");
		        printf("              -a additional input data (predefined charges, atom groups etc).)\n");
			printf("	      -w weight of charge constraint, in default, 0.0005 for resp1 and 0.001 fore resp2\n");
			exit(1);
		}
	}

	method = -1;
	max_path_length = -1;
	for (i = 1; i < argc; i += 2) {
		if (strcmp(argv[i], "-i") == 0)
			strcpy(ifilename, argv[i + 1]);
		if (strcmp(argv[i], "-o") == 0)
			strcpy(ofilename, argv[i + 1]);
		if (strcmp(argv[i], "-a") == 0) {
			strcpy(afilename, argv[i + 1]);
			iaddinfo = 1;
		}
		if (strcmp(argv[i], "-n") == 0)
			nconf = atoi(argv[i+1]);
		if (strcmp(argv[i], "-l") == 0)
			max_path_length = atoi(argv[i+1]); 
		if (strcmp(argv[i], "-f") == 0) {
			if (strcmp("resp", argv[i + 1]) == 0) 
				method = -1;
			if (strcmp("resp0", argv[i + 1]) == 0) 
				method = 0;
			if (strcmp("resp1", argv[i + 1]) == 0) {
				method = 1;
				weight = 0.0005;
			}
			if (strcmp("resp2", argv[i + 1]) == 0) {
				method = 2;
				weight = 0.001;
			}
			if (strcmp("resp3", argv[i + 1]) == 0) 
				method = 3;
                        if (strcmp("iresp1", argv[i + 1]) == 0) {
                                method = 4;
				weight = 0.0005;
			}
                        if (strcmp("iresp2", argv[i + 1]) == 0) {
                                method = 5;
				weight = 0.001;
			}
			if (strcmp("resp4", argv[i + 1]) == 0) {
				method = 6;
				weight = 0.000;
			}
			if (strcmp("iresp0", argv[i + 1]) == 0) 
				method = 7;
		}
		if (strcmp(argv[i], "-w") == 0)
			weight = atof(argv[i+1]); 
		if (strcmp(argv[i], "-e") == 0)
			iequ = atoi(argv[i+1]); 
	}
	if(nconf < 1) {
		printf("\nNumber of conformations must be equal to or larger than 1"); 
		exit(1);
	}
	if(iequ != 0 && iequ != 1 && iequ !=2)
		iequ = 1;
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
	if(minfo.dcharge >= -9990) charge = minfo.dcharge;

	equ_atom_id = (int *) malloc(sizeof(int) * atomnum);
	if (equ_atom_id == NULL) {
		fprintf(stdout, "memory allocation error for *equ_atom_id\n");
		exit(1);
	}
        for (i = 0; i < atomnum; i++)
                equ_atom_id[i] = -1;

	if(method == 3) {
		atqwt = (double *) malloc(sizeof(double) * atomnum);
		if (atqwt == NULL) {
			fprintf(stdout, "memory allocation error for *atqwt\n");
			exit(1);
		}
		refcharge = (double *) malloc(sizeof(double) * atomnum);
		if (refcharge == NULL) {
			fprintf(stdout, "memory allocation error for *refcharge\n");
			exit(1);
		}
		for(i=0;i<atomnum;i++) {
			atqwt[i] = 0;
			refcharge[i] = 0;
		}
//	now read in atomic weights and reference charge parameters 
                pfilename[0] = '\0';
                strcpy(pfilename, amberhome);
                strcat(pfilename, "/dat/antechamber/RESPPARM.DAT");
		readparm(pfilename);

        	wac("ANTECHAMBER_RESP.AC", atomnum, atom, bondnum, bond, cinfo, minfo);
        	copied_size = build_exe_path(command, 
                "atomtype -i ANTECHAMBER_RESP.AC -o ANTECHAMBER_RESP_AT.AC -d ",
                 sizeof command, 1 );
        	dfilename[0] = '\0';
        	strcpy(dfilename, amberhome);
        	strcat(dfilename, "/dat/antechamber/ATOMTYPE_RESP.DEF");
        	strncat(command, dfilename, MAXCHAR - copied_size );

        	if (cinfo.intstatus == 2)
                	fprintf(stdout, "Running: %s\n", command);
        	status = system(command);
        	if(status != 0) {
                	fprintf(stdout, "Error: cannot run \"%s\" in respgen.c properly, exit\n", command);
                	exit(1);
        	}
		assignparm("ANTECHAMBER_RESP_AT.AC");
	}

	if(iequ == 1)
		identify_equatom(atomnum, atom, equ_atom_id, max_path_length, bondnum, bond, 0);
	if(iequ == 2)
		identify_equatom(atomnum, atom, equ_atom_id, max_path_length, bondnum, bond, 1);

	icharge = (int *) malloc(sizeof(int) * atomnum);
       	if (icharge == NULL) {
               	fprintf(stdout, "memory allocation error for *icharge in respin()\n");
               	exit(1);
       	}
	pcharge = (double *) malloc(sizeof(double) * atomnum);
       	if (pcharge == NULL) {
               	fprintf(stdout, "memory allocation error for *pcharge in respin()\n");
               	exit(1);
       	}
	for(i=0;i<atomnum;i++) {
		pcharge[i] = 0.0;	
		icharge[i] = 0;
	}

	if(iaddinfo == 1) {
		if ((fpaddinfo = fopen(afilename, "r")) == NULL) {
        		fprintf(stdout, "Cannot open the additional file %s to read in main(), exit\n", afilename);
        		exit(1);
		}
		if ((fpcharge = fopen("QIN", "w")) == NULL) {
        		fprintf(stdout, "Cannot open file QIN, exit\n");
        		exit(1);
		}
		readinfo();
	}
	respin(method);
	if(iaddinfo == 1) {
		fclose(fpcharge);
		fclose(fpaddinfo);
	}
	printf("\n");
/*
	 free(atom);
	 free(selectindex);
	 free(selectelement);
	 free(equ_atom_id);
	 free(pathnum);
	 free(pathatomnum);
	 for (i =0 ;i <atomnum; i++) free(pathscore[i]);
*/
	return (0);

}
