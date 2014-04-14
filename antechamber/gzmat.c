/* GZMAT */
int rgzmat(char *filename, int *atomnum, ATOM * atom, CONTROLINFO cinfo,
		   MOLINFO minfo)
{
        typedef struct {
                char name[MAXCHAR];
        } STRNAME;

	FILE *fpin;
	int i, j, index, index0;
	int overflow_flag = 0;
	int findindex;
	int numatom;
	int coordinate_flag = 0;
	STRNAME *bondstr;	
	STRNAME *anglestr;	
	STRNAME *twiststr;	
	char tmpchar1[MAXCHAR];
	char tmpchar2[MAXCHAR];
	char line[MAXCHAR];


	if ((fpin = fopen(filename, "r")) == NULL) {
		fprintf(stdout, "Cannot open the input file %s to read in rgzmat(), exit\n", filename);
		exit(1);
	}
        bondstr = (STRNAME *) malloc(sizeof(STRNAME) * (cinfo.maxatom +10));
        if (bondstr == NULL) {
                fprintf(stdout, "memory allocation error for *bondstr in rgzmat()\n");
                exit(1);
        }
        anglestr = (STRNAME *) malloc(sizeof(STRNAME) * (cinfo.maxatom +10));
        if (anglestr == NULL) {
                fprintf(stdout, "memory allocation error for *anglestr in rgzmat()\n");
                exit(1);
        }
        twiststr = (STRNAME *) malloc(sizeof(STRNAME) * (cinfo.maxatom +10));
        if (twiststr == NULL) {
                fprintf(stdout, "memory allocation error for *twiststr in rgzmat()\n");
                exit(1);
        }

	initial(cinfo.maxatom, atom, minfo.resname);
	index = 0;
	index0 = 1;
	numatom = 0;
	for (;;) {
		if (fgets(line, MAXCHAR, fpin) == NULL) {
/*       printf("\nFinished reading %s file.", cinfo.ifilename); */
			break;
		}
		if (spaceline(line) == 1) {
			if(coordinate_flag == 1) break;
			index++;
			continue;
		}
		if (index >= 2)
			index++;
		if (index <= 3)
			continue;
		if (index >= 4) {
			if (spaceline(line) == 1 || strncmp(line, "Vari", 4) == 0
				|| strncmp(line, "vari", 4) == 0)
				index0 = -1;
			if (strncmp(line, "Const", 5) == 0
				|| strncmp(line, "const", 5) == 0)
				index0 = -1;
		}
		if (index == 4) {
			if (overflow_flag == 0)
				sscanf(line, "%s", atom[numatom].name);
			numatom++;
			if (numatom >= cinfo.maxatom && overflow_flag == 0) {
				printf
					("\nInfo: the atom number exceeds the MAXATOM, reallocate memory automatically");
				overflow_flag = 1;
			}
			continue;
		}
		if (index == 5) {
			if (overflow_flag == 0) {
				sscanf(line, "%s%d%s", atom[numatom].name,
					   &atom[numatom].bondatom, bondstr[numatom].name);
				if(atom[numatom].bondatom > numatom) {
					printf("\nError: bond atom ID is larger than ID of current atom (%d,%s), exit", 
						numatom+1,  atom[numatom].name); 
					exit(1);
				}
			}
			numatom++;
			if (numatom >= cinfo.maxatom && overflow_flag == 0) {
				printf
					("\nInfo: the atom number exceeds the MAXATOM, reallocate memory automatically");
				overflow_flag = 1;
			}
			continue;
		}
		if (index == 6) {
			if (overflow_flag == 0) {
				sscanf(line, "%s%d%s%d%s", atom[numatom].name,
					   &atom[numatom].bondatom, bondstr[numatom].name,
					   &atom[numatom].angleatom, anglestr[numatom].name);
				if(atom[numatom].bondatom > numatom) {
					printf("\nError: bond atom ID is larger than ID of current atom (%d,%s), exit", 
						numatom+1,  atom[numatom].name); 
					exit(1);
				}
				if(atom[numatom].angleatom > numatom) {
					printf("\nError: angle atom ID is larger than ID of current atom (%d,%s), exit", 
						numatom+1,  atom[numatom].name); 
					exit(1);
				}
			}
			numatom++;
			if (numatom >= cinfo.maxatom && overflow_flag == 0) {
				printf
					("\nInfo: the atom number exceeds the MAXATOM, reallocate memory automatically");
				overflow_flag = 1;
			}
			continue;
		}
		if (index0 != -1) {
			if (overflow_flag == 0) {
				sscanf(line, "%s%d%s%d%s%d%s", atom[numatom].name,
					   &atom[numatom].bondatom, bondstr[numatom].name,
					   &atom[numatom].angleatom, anglestr[numatom].name,
					   &atom[numatom].twistatom, twiststr[numatom].name);
				if(atom[numatom].bondatom > numatom) {
					printf("\nError: bond atom ID is larger than ID of current atom (%d,%s), exit", 
						numatom+1,  atom[numatom].name); 
					exit(1);
				}
				if(atom[numatom].angleatom > numatom) {
					printf("\nError: angle atom ID is larger than ID of current atom (%d,%s), exit", 
						numatom+1,  atom[numatom].name); 
					exit(1);
				}
				if(atom[numatom].twistatom > numatom) {
					printf("\nError: torsional atom ID is larger than ID of current atom (%d,%s), exit", 
						numatom+1,  atom[numatom].name); 
					exit(1);
				}
			}
			numatom++;
			if (numatom >= cinfo.maxatom && overflow_flag == 0) {
				printf
					("\nInfo: the atom number exceeds the MAXATOM, reallocate memory automatically");
				overflow_flag = 1;
			}

			continue;
		}
		if (index0 == -1) {
         if (strchr(line, '=') == NULL)
             continue;
			coordinate_flag = 1;
			for(i=0;i<strlen(line);i++)
                       		if(line[i] == '=') line[i] = ' ';
			sscanf(line, "%s %s", tmpchar1, tmpchar2);
			for (i = 1; i < numatom; i++) {
				findindex = 1;
				for (j = 0; j < strlen(bondstr[i].name); j++)
					if (bondstr[i].name[j] != tmpchar1[j]) {
						findindex = 0;
						break;
					}
				if (findindex == 1) {
					strcpy(bondstr[i].name, tmpchar2);
					break;
				}

			}
			for (i = 2; i < numatom; i++) {
				findindex = 1;
				for (j = 0; j < strlen(anglestr[i].name); j++)
					if (anglestr[i].name[j] != tmpchar1[j]) {
						findindex = 0;
						break;
					}
				if (findindex == 1) {
					strcpy(anglestr[i].name, tmpchar2);
					break;
				}
			}
			for (i = 3; i < numatom; i++) {
				findindex = 1;
				for (j = 0; j < strlen(twiststr[i].name); j++)
					if (twiststr[i].name[j] != tmpchar1[j]) {
						findindex = 0;
						break;
					}
				if (findindex == 1) {
					strcpy(twiststr[i].name, tmpchar2);
					break;
				}
			}
		}
	}
	atom[1].bondatom--;
	atom[2].bondatom--;
	atom[2].angleatom--;
	for (i = 3; i < numatom; i++) {
		atom[i].bondatom--;
		atom[i].angleatom--;
		atom[i].twistatom--;
	}
	for (i = 1; i < numatom; i++)
		atom[i].bond = atof(bondstr[i].name);
	for (i = 2; i < numatom; i++)
		atom[i].angle = atof(anglestr[i].name);
	for (i = 3; i < numatom; i++)
		atom[i].twist = atof(twiststr[i].name);
	*atomnum = numatom;
/* printf("\n atom number is  %5d", *atomnum); */
	fclose(fpin);
	free(bondstr);
	free(anglestr);
	free(twiststr);
	return overflow_flag;

}

void wgzmat(char *filename, int atomnum, ATOM atom[], MOLINFO minfo)
{
        FILE *fpin;
        FILE *fpout;
        char *amberhome;
        char espparm_file[MAXCHAR];
        char line[MAXCHAR];
        char akeyword[MAXCHAR]="";
        char ckeyword[MAXCHAR];
	char tmpchar0[MAXCHAR];
	char tmpchar1[MAXCHAR];
	char tmpchar2[MAXCHAR];
	char tmpchar3[MAXCHAR];
        int i,j;
        int iradius_flag;
        int ibs= 0;
        int esp_flag;
        int nespparm = 0;
        int nbasisset = 0;
        double default_radius = 1.7;
        ESPPARM espparm[120];
        BASISSET basisset[100];


	if ((fpout = fopen(filename, "w")) == NULL) {
		fprintf(stdout, "Cannot open a file %s to write in wgzmat(), exit\n", filename);
		exit(1);
	}
	intercoord(atomnum, atom);
        fprintf(fpout, "%s\n", "--Link1--");
	if(strlen(minfo.gn) >= 5) 
		fprintf(fpout, "%s\n", minfo.gn);
	fprintf(fpout, "%s%s\n", "%chk=", minfo.chkfile);
	if(strlen(minfo.gm) >= 4) 
		fprintf(fpout, "%s\n", minfo.gm);

//      check ESP-related keyword
        esp_flag = 0;
	for(i=0;i<strlen(minfo.gkeyword);i++)
		ckeyword[i] = toupper(minfo.gkeyword[i]);
        if((strstr(ckeyword, "POP=") != 0 || strstr(ckeyword, "POP(")  != 0)  &&
           (strstr(ckeyword, "MK")   != 0 || strstr(ckeyword, "CHELP") != 0)) 
                esp_flag = 1;
//      when the default gaussian keyword is used, or esp_flag ==1, read ESP.PARM
        if(minfo.igkeyword == 0 || esp_flag == 1) {
                amberhome = (char *) getenv("AMBERHOME");
                if( amberhome == NULL ){
                        fprintf( stdout, "AMBERHOME is not set!\n" );
                        exit(1);
                }
                strcpy(espparm_file, amberhome);
                strcat(espparm_file, "/dat/antechamber/ESPPARM.DAT");
                if ((fpin = fopen(espparm_file, "r")) == NULL) {
                        fprintf(stdout, "Cannot open espparm_file %s in read_espparm(), exit\n", espparm_file);
                        exit(1);
                }
                for (;;) {
                        if (fgets(line, MAXCHAR, fpin) == NULL) break;
                        if(strncmp(line, "DEFAULT RADIUS", 14) == 0)
                                sscanf(&line[14], "%lf", &default_radius);
                        if(strncmp(line, "BASIS SET", 9) == 0) {
                                sscanf(&line[9], "%d%s", &basisset[nbasisset].id,
                                                           basisset[nbasisset].bs);
                                nbasisset++;
                        }
                        if(strncmp(line, "PARM", 4) == 0) {
                                sscanf(&line[4], "%d%s%lf%lf%d%d", &espparm[nespparm].atomicnum, espparm[nespparm].elem,
                                                                    &espparm[nespparm].vdw,      &espparm[nespparm].mk,
                                                                    &espparm[nespparm].flag,     &espparm[nespparm].bs);
                                nespparm++;
                        }
                }
                fclose(fpin);

                iradius_flag = 0;
                ibs = 0;
                for(i=0;i<atomnum;i++) {
                        for(j=0;j<nespparm;j++)
                                if(atom[i].atomicnum == espparm[j].atomicnum || strcmp(atom[i].element, espparm[j].elem) == 0) {
                                        if(ibs < espparm[j].bs) ibs=espparm[j].bs;
                                        if(espparm[j].flag != 0) {
                                                iradius_flag = 1;
                                                espparm[j].flag = 2;
                                        }
                                        break;
                                }
                }

                if(minfo.igkeyword == 0) {
                        strcpy(minfo.gkeyword, "#HF/");
                        strcat(minfo.gkeyword, basisset[ibs-1].bs);
                        strcat(minfo.gkeyword, " SCF=tight Test Pop=MK iop(6/33=2) iop(6/42=6) opt");
                }
        }
//      additional keywords
        if(esp_flag == 1) {
                if(iradius_flag == 1) {
                        if(strstr(minfo.gkeyword, "ReadRadii") == 0 &&
                           strstr(minfo.gkeyword, "READRADII") == 0 &&
                           strstr(minfo.gkeyword, "readradii") == 0) {
                                strcat(akeyword, "Pop=ReadRadii");
                        }
                }
                if(minfo.gv == 1) {
                        if(strstr(minfo.gkeyword, "6/50=1") == 0) {
                                strcat(akeyword, " iop(6/50=1)");
                        }
                }
        }

        if(strlen(akeyword) >= 1) {
                fprintf(fpout, "%s\n", minfo.gkeyword);
                fprintf(fpout, "#%s\n", akeyword);
                fprintf(fpout, "\n");
        }
        else
                fprintf(fpout, "%s\n\n", minfo.gkeyword);
	fprintf(fpout, "%s\n\n", "remark line goes here");
	fprintf(fpout, "%d%4d\n", minfo.icharge, minfo.multiplicity);
	element(atomnum, atom);
	for (i = 0; i < atomnum; i++) {
		/* newitoa(i + 1, tmpchar0); */
		sprintf(tmpchar0, "%d", i+1);
		if (i == 0) {
			fprintf(fpout, "%5s\n", atom[i].element);
			continue;
		}
		if (i == 1) {
			strcpy(tmpchar1, "b");
			strcat(tmpchar1, tmpchar0);
			fprintf(fpout, "%5s%5d%8s\n", atom[i].element,
					atom[i].bondatom + 1, tmpchar1);
			continue;
		}
		if (i == 2) {
			strcpy(tmpchar1, "b");
			strcat(tmpchar1, tmpchar0);
			fprintf(fpout, "%5s%5d%8s", atom[i].element,
					atom[i].bondatom + 1, tmpchar1);
			strcpy(tmpchar2, "a");
			strcat(tmpchar2, tmpchar0);
			fprintf(fpout, "%5d%8s\n", atom[i].angleatom + 1, tmpchar2);
			continue;
		}
		strcpy(tmpchar1, "b");
		strcat(tmpchar1, tmpchar0);
		fprintf(fpout, "%5s%5d%8s", atom[i].element, atom[i].bondatom + 1,
				tmpchar1);
		strcpy(tmpchar2, "a");
		strcat(tmpchar2, tmpchar0);
		fprintf(fpout, "%5d%8s", atom[i].angleatom + 1, tmpchar2);
		strcpy(tmpchar3, "t");
		strcat(tmpchar3, tmpchar0);
		fprintf(fpout, "%5d%8s\n", atom[i].twistatom + 1, tmpchar3);
	}

	fprintf(fpout, "Variables:\n");
	fprintf(fpout, "b2= %8.4lf\n", atom[1].bond);
	fprintf(fpout, "b3= %8.4lf\n", atom[2].bond);
	fprintf(fpout, "a3= %8.4lf\n", atom[2].angle);
	for (i = 3; i < atomnum; i++) {
		/* newitoa(i + 1, tmpchar0); */
		sprintf(tmpchar0, "%d", i+1);
		strcpy(tmpchar1, "b");
		strcat(tmpchar1, tmpchar0);
		strcpy(tmpchar2, "a");
		strcat(tmpchar2, tmpchar0);
		strcpy(tmpchar3, "t");
		strcat(tmpchar3, tmpchar0);
		fprintf(fpout, "%s= %8.4lf\n", tmpchar1, atom[i].bond);
		fprintf(fpout, "%s= %8.4lf\n", tmpchar2, atom[i].angle);
		fprintf(fpout, "%s= %8.4lf\n", tmpchar3, atom[i].twist);
	}
        if(esp_flag == 1) {
                if(minfo.gv == 1 && iradius_flag == 1) {
                        for(i=0;i<nespparm;i++)
                                if(espparm[i].flag == 2) {
                                        if(espparm[i].mk <= 0) espparm[i].mk=default_radius;
                                        fprintf(fpout, "\n%s     %8.2lf", espparm[i].elem, espparm[i].mk);
                                }
                        fprintf(fpout, "\n\n%s\n", minfo.gesp);
                        for(i=0;i<nespparm;i++)
                                if(espparm[i].flag == 2) {
                                        if(espparm[i].mk <= 0) espparm[i].mk=default_radius;
                                        fprintf(fpout, "\n%s     %8.2lf", espparm[i].elem, espparm[i].mk);
                                }
                        fprintf(fpout, "\n\n%s\n", minfo.gesp);
                }
                if(minfo.gv == 0 && iradius_flag == 1) {
                        for(i=0;i<nespparm;i++)
                                if(espparm[i].flag == 2) {
                                        if(espparm[i].mk <= 0) espparm[i].mk=default_radius;
                                        fprintf(fpout, "\n%s     %8.2lf", espparm[i].elem, espparm[i].mk);
                                }
                        fprintf(fpout, "\n");
                        for(i=0;i<nespparm;i++)
                                if(espparm[i].flag == 2) {
                                        if(espparm[i].mk <= 0) espparm[i].mk=default_radius;
                                        fprintf(fpout, "\n%s     %8.2lf", espparm[i].elem, espparm[i].mk);
                                }
                        fprintf(fpout, "\n");
                }
                if(minfo.gv == 1 && iradius_flag == 0) {
                        fprintf(fpout, "\n%s\n", minfo.gesp);
                        fprintf(fpout, "\n%s\n", minfo.gesp);
                }
        }
	fprintf(fpout, "\n\n");
	fclose(fpout);
}
