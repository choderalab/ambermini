/* GZMAT */
int rjzmat(char *filename, int *atomnum, ATOM * atom, CONTROLINFO cinfo,
		   MOLINFO minfo)
{
        typedef struct {
                char str[MAXCHAR];
        } STRNAME;

	FILE *fpin;
	int i, flag;
	int index = 0;
	int tmpint;
	int overflow_flag = 0;
	int findindex;
	int numatom;
	char tmpchar[MAXCHAR];
	char tmpchar1[MAXCHAR];
	char tmpchar2[MAXCHAR];
	char tmpchar3[MAXCHAR];
	char tmpchar4[MAXCHAR];
	char tmpchar5[MAXCHAR];
	char tmpchar6[MAXCHAR];
	STRNAME *bondstr;
	STRNAME *anglestr;
	STRNAME *twiststr;
	char line[MAXCHAR];

        bondstr = (STRNAME *) malloc(sizeof(STRNAME) * (cinfo.maxatom + 10));
        if (bondstr == NULL) {
                fprintf(stdout, "memory allocation error for *bondstr in rjzmat()\n");
                exit(1);
        }
        anglestr = (STRNAME *) malloc(sizeof(STRNAME) * (cinfo.maxatom + 10));
        if (anglestr == NULL) {
                fprintf(stdout, "memory allocation error for *anglestr in rjzmat()\n");
                exit(1);
        }
        twiststr = (STRNAME *) malloc(sizeof(STRNAME) * (cinfo.maxatom + 10));
        if (twiststr == NULL) {
                fprintf(stdout, "memory allocation error for *twiststr in rjzmat()\n");
                exit(1);
        }

	if ((fpin = fopen(filename, "r")) == NULL) {
		fprintf(stdout, "Cannot open the input file %s to read in rjzmat(), exit\n", filename);
		exit(1);
	}
	initial(cinfo.maxatom, atom, minfo.resname);
	numatom = 0;
        flag = 0;

        for (;;) {
                if (fgets(line, MAXCHAR, fpin) == NULL) {
/*     printf("\nFinished reading %s file.", cinfo.ifilename); */
                        break;
                }
		strcpy(tmpchar, "");
		tmpchar[0]='\0';
                sscanf(line, "%s", tmpchar);
                if(strncmp(tmpchar, "&zmat", 5)==0) {
                        flag =1;
                        continue;
                }
                if(strncmp(tmpchar, "&zvar", 5)==0) {
                        index =1;
			break;
                }
                if (flag == 1 && strcmp(tmpchar, "&") ==0) {
			flag = 0;
			continue;
		}
                if(strncmp(tmpchar, "molchg=", 7)==0) {
                        tmpint = 0;
                        for(i=7;i<= strlen(tmpchar);i++)
                                tmpchar2[tmpint++] = tmpchar[i];
                        minfo.icharge = atoi(tmpchar2);
                        continue;
                }
                if(strncmp(tmpchar, "multip=", 7)==0) {
                        tmpint = 0;
                        for(i=7;i<= strlen(tmpchar);i++)
                                tmpchar2[tmpint++] = tmpchar[i];
                        minfo.multiplicity = atoi(tmpchar2);
                        continue;
                }
                if (overflow_flag == 0 && flag == 1) {
                	sscanf(line, "%s", atom[numatom].name);
			if (numatom >= cinfo.maxatom && overflow_flag == 0) {
				printf ("\nInfo: the atom number exceeds the MAXATOM, reallocate memory automatically");
				overflow_flag = 1;
			}
			numatom++;
			continue;
		}
	}
	*atomnum = numatom;
	rewind(fpin);
	flag = 0;
	numatom = 0;	
        for (;;) {
                if (fgets(line, MAXCHAR, fpin) == NULL) break;
		strcpy(tmpchar, "");
		tmpchar[0]='\0';
                sscanf(line, "%s", tmpchar);
                if(strncmp(tmpchar, "&zmat", 5)==0) {
                        flag =1;
                        continue;
                }
                if(strncmp(tmpchar, "&zvar", 5)==0) {
                        flag =2;
			continue;
                }
		if(flag == 1 && strcmp(tmpchar, "&")==0) {
			flag = 0;
			if(index == 1) 
				continue;
			else
				break;
		}
		if(flag == 2 && strcmp(tmpchar, "&")==0) 
				break;

                if (overflow_flag == 0 && flag == 1) {
			strcpy(tmpchar1, "");
			strcpy(tmpchar2, "");
			strcpy(tmpchar3, "");
			strcpy(tmpchar4, "");
			tmpchar1[0]='\0';
			tmpchar2[0]='\0';
			tmpchar3[0]='\0';
			tmpchar4[0]='\0';
			sscanf(line, "%s%s%s%s%s%s%s", tmpchar1, tmpchar2, bondstr[numatom].str, tmpchar3, anglestr[numatom].str, tmpchar4, twiststr[numatom].str);
			if(numatom > 0) {
				findindex = 0;
				for(i=0;i<numatom;i++) {
					if(strcmp(atom[i].name, tmpchar2)==0){
						atom[numatom].bondatom = i;
						findindex = 1;
						break;	
					}
				}
				if(findindex == 0) {
					printf("\nThe bond atom (%s) does not appear before %s", tmpchar2, tmpchar1);
					exit(1);
				}
			}

			if(numatom > 1) {
				findindex = 0;
				for(i=0;i<numatom;i++) {
					if(strcmp(atom[i].name, tmpchar3)==0){
						atom[numatom].angleatom = i;
						findindex = 1;
						break;	
					}
				}
				if(findindex == 0) {
					printf("\nThe angle atom (%s) does not appear before %s", tmpchar3, tmpchar1);
					exit(1);
				}
			}

			if(numatom > 2) {
				findindex = 0;
				for(i=0;i<numatom;i++) {
					if(strcmp(atom[i].name, tmpchar4)==0){
						atom[numatom].twistatom = i;
						findindex = 1;
						break;	
					}
				}
				if(findindex == 0) {
					printf("\nThe torsional angle atom (%s) does not appear before %s", tmpchar4, tmpchar1);
					exit(1);
				}
			}
			numatom++;

		}

                if (flag == 2 && overflow_flag == 0) {
			for(i=0; i<strlen(line); i++)
				if(line[i]=='=') line[i]=' ';	
			strcpy(tmpchar1, "");
			tmpchar1[0]='\0';
			strcpy(tmpchar2, "");
			tmpchar2[0]='\0';
			strcpy(tmpchar3, "");
			tmpchar3[0]='\0';
			strcpy(tmpchar4, "");
			tmpchar4[0]='\0';
			strcpy(tmpchar5, "");
			tmpchar5[0]='\0';
			strcpy(tmpchar6, "");
			tmpchar6[0]='\0';
			sscanf(line, "%s%s%s%s%s%s", tmpchar1, tmpchar2, tmpchar3, tmpchar4, tmpchar5, tmpchar6); 
			for(i=0;i<*atomnum;i++) {
				if(strcmp(bondstr[i].str, tmpchar1)==0) 
					strcpy(bondstr[i].str, tmpchar2);
				if(strcmp(anglestr[i].str, tmpchar3)==0) 
					strcpy(anglestr[i].str, tmpchar4);
				if(strcmp(twiststr[i].str, tmpchar5)==0) 
					strcpy(twiststr[i].str, tmpchar6);
			}
		}
	}

	for (i = 1; i < *atomnum; i++)
		atom[i].bond = atof(bondstr[i].str);
	for (i = 2; i < numatom; i++)
		atom[i].angle = atof(anglestr[i].str);
	for (i = 3; i < numatom; i++)
		atom[i].twist = atof(twiststr[i].str);
/* printf("\n atom number is  %5d", *atomnum); */
	fclose(fpin);
	free(bondstr);
	free(anglestr);
	free(twiststr);
	return overflow_flag;

}
void wjzmat(char *filename, int atomnum, ATOM atom[], MOLINFO minfo)
{
	FILE *fpout;
	int i;
	/* int index; */
	char tmpchar0[MAXCHAR];
	char tmpchar1[MAXCHAR];
	char tmpchar2[MAXCHAR];
	char tmpchar3[MAXCHAR];

	if ((fpout = fopen(filename, "w")) == NULL) {
		fprintf(stdout, "Cannot open a file %s to write in wjzmat(), exit\n", filename);
		exit(1);
	}
	intercoord(atomnum, atom);
        fprintf(fpout, "%s\n", "This Jaguar input file generated by Antechamber");
	fprintf(fpout, "%s\n", "&gen");
        fprintf(fpout, "%s%d\n", "molchg=", minfo.icharge);
        fprintf(fpout, "%s%d\n", "multip=", minfo.multiplicity);
        fprintf(fpout, "igeopt=1\n");
        fprintf(fpout, "basis=6-31G*\n");
        fprintf(fpout, "dft=b3lyp\n");
        fprintf(fpout, "%s\n", "&");
        fprintf(fpout, "%s\n", "&zmat");

	for (i = 0; i < atomnum; i++) {
		/* newitoa(i + 1, tmpchar0); */
		sprintf(tmpchar0, "%d", i+1);
		if (i == 0) {
			fprintf(fpout, "%8s\n", atom[i].name);
			continue;
		}
		if (i == 1) {
			strcpy(tmpchar1, "r");
			strcat(tmpchar1, tmpchar0);
			fprintf(fpout, "%8s%8s%8s\n", atom[i].name,
					atom[atom[i].bondatom].name, tmpchar1);
			continue;
		}
		if (i == 2) {
			strcpy(tmpchar1, "r");
			strcat(tmpchar1, tmpchar0);
			fprintf(fpout, "%8s%8s%8s", atom[i].name,
					atom[atom[i].bondatom].name , tmpchar1);
			strcpy(tmpchar2, "a");
			strcat(tmpchar2, tmpchar0);
			fprintf(fpout, "%8s%8s\n", atom[atom[i].angleatom].name , tmpchar2);
			continue;
		}
		strcpy(tmpchar1, "r");
		strcat(tmpchar1, tmpchar0);
		fprintf(fpout, "%8s%8s%8s", atom[i].name, atom[atom[i].bondatom].name ,
				tmpchar1);
		strcpy(tmpchar2, "a");
		strcat(tmpchar2, tmpchar0);
		fprintf(fpout, "%8s%8s", atom[atom[i].angleatom].name , tmpchar2);
		strcpy(tmpchar3, "d");
		strcat(tmpchar3, tmpchar0);
		fprintf(fpout, "%8s%8s\n", atom[atom[i].twistatom].name, tmpchar3);
	}

	fprintf(fpout, "&\n&zvar\n");
	fprintf(fpout, "r2= %13.8lf\n", atom[1].bond);
	fprintf(fpout, "r3= %13.8lf  ", atom[2].bond);
	fprintf(fpout, "a3= %13.8lf\n", atom[2].angle);
	for (i = 3; i < atomnum; i++) {
		/* newitoa(i + 1, tmpchar0); */
		sprintf(tmpchar0, "%d", i+1);
		strcpy(tmpchar1, "r");
		strcat(tmpchar1, tmpchar0);
		strcpy(tmpchar2, "a");
		strcat(tmpchar2, tmpchar0);
		strcpy(tmpchar3, "d");
		strcat(tmpchar3, tmpchar0);
		fprintf(fpout, "%s= %13.8lf  ", tmpchar1, atom[i].bond);
		fprintf(fpout, "%s= %13.8lf  ", tmpchar2, atom[i].angle);
		fprintf(fpout, "%s= %13.8lf\n", tmpchar3, atom[i].twist);
	}
	fprintf(fpout, "&\n\n");
	fclose(fpout);
}
