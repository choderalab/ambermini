int rhin(char *filename, int *atomnum, ATOM atom[], int *bondnum,
		 BOND * bond, CONTROLINFO cinfo, MOLINFO minfo)
{
	int i;
	int numatom, numbond;
	int overflow_flag = 0;
	int tmpint1, tmpint2, tmpint3, tmpint4;
	int tmpint5, tmpint6, tmpint7, tmpint8;
	int index;
	char tmpchar[MAXCHAR];
	char tmpchar1[MAXCHAR], tmpchar2[MAXCHAR], tmpchar3[MAXCHAR],
		tmpchar4[MAXCHAR];
	char tmpchar5[MAXCHAR], tmpchar6[MAXCHAR], tmpchar7[MAXCHAR],
		tmpchar8[MAXCHAR];
	char tmpchar9[MAXCHAR], tmpchar10[MAXCHAR], tmpchar11[MAXCHAR];
	char line[MAXCHAR];
	double tmpfloat1, tmpfloat2, tmpfloat3, tmpfloat4;
	FILE *fpin;

	if ((fpin = fopen(filename, "r")) == NULL) {
		fprintf(stdout, "Cannot open the hin file %s to read, exit\n", filename);
		exit(1);
	}
	initial(cinfo.maxatom, atom, minfo.resname);
	numatom = 0;
	numbond = 0;
	index = 0;
	for (;;) {
		if (fgets(line, MAXCHAR, fpin) == NULL) {
			/*  printf("\nFinished reading %s file.", filename); */
			break;
		}
		index++;
		if (index == 1) {
			sscanf(line, "%s", tmpchar);
			continue;
		}
		tmpint1 = -1;
		tmpint2 = -1;
		tmpint3 = -1;
		tmpint4 = -1;
		tmpint5 = -1;
		tmpint6 = -1;
		tmpint7 = -1;
		tmpint8 = -1;

		if (strncmp("atom", line, 4) == 0) {
			strcpy(tmpchar1, "");
			tmpchar1[0] = '\0';
			strcpy(tmpchar2, "");
			tmpchar2[0] = '\0';
			strcpy(tmpchar3, "");
			tmpchar3[0] = '\0';
			strcpy(tmpchar4, "");
			tmpchar4[0] = '\0';
			strcpy(tmpchar5, "");
			tmpchar5[0] = '\0';
			strcpy(tmpchar6, "");
			tmpchar6[0] = '\0';
			strcpy(tmpchar7, "");
			tmpchar7[0] = '\0';
			strcpy(tmpchar8, "");
			tmpchar8[0] = '\0';
			strcpy(tmpchar9, "");
			tmpchar9[0] = '\0';
			strcpy(tmpchar10, "");
			tmpchar10[0] = '\0';
			strcpy(tmpchar11, "");
			tmpchar11[0] = '\0';
			sscanf(line,
				   "%s%d%s%s%s%s%lf%lf%lf%lf%d%d%s%d%s%d%s%d%s%d%s%d%s",
				   tmpchar1, &tmpint1, tmpchar2, tmpchar3, tmpchar4,
				   tmpchar5, &tmpfloat1, &tmpfloat2, &tmpfloat3,
				   &tmpfloat4, &tmpint2, &tmpint3, tmpchar6, &tmpint4,
				   tmpchar7, &tmpint5, tmpchar8, &tmpint6, tmpchar9,
				   &tmpint7, tmpchar10, &tmpint8, tmpchar11);
			if (overflow_flag == 0) {
				atom[numatom].charge = tmpfloat1;
				atom[numatom].x = tmpfloat2;
				atom[numatom].y = tmpfloat3;
				atom[numatom].z = tmpfloat4;
				strcpy(atom[numatom].name, tmpchar3);
			}
			if (tmpint3 > 0) {
				if (overflow_flag == 0) {
					atom[numatom].con[0] = tmpint3 - 1;
					atom[numatom].connum++;
				}
				if (numatom < tmpint3 - 1) {
					if (overflow_flag == 0) {
						bond[numbond].bondi = numatom;
						bond[numbond].bondj = tmpint3 - 1;
						bond[numbond].type = 1;
						if (tmpchar6[0] == 'd')
							bond[numbond].type = 2;
						if (tmpchar6[0] == 't')
							bond[numbond].type = 3;
					}
					numbond++;
					if (numbond >= cinfo.maxbond && overflow_flag == 0) {
						printf
							("\nInfo: the bond number exceeds the MAXBOND, reallocate memory automatically");
						overflow_flag = 1;
					}
				}
			}

			if (tmpint4 > 0) {
				if (overflow_flag == 0) {
					atom[numatom].con[1] = tmpint4 - 1;
					atom[numatom].connum++;
				}
				if (numatom < tmpint4 - 1) {
					if (overflow_flag == 0) {
						bond[numbond].bondi = numatom;
						bond[numbond].bondj = tmpint4 - 1;
						bond[numbond].type = 1;
						if (tmpchar7[0] == 'd')
							bond[numbond].type = 2;
						if (tmpchar7[0] == 't')
							bond[numbond].type = 3;
					}
					numbond++;
					if (numbond >= cinfo.maxbond && overflow_flag == 0) {
						printf
							("\nInfo: the bond number exceeds the MAXBOND, reallocate memory automatically");
						overflow_flag = 1;
					}
				}
			}


			if (tmpint5 > 0) {
				if (overflow_flag == 0) {
					atom[numatom].con[2] = tmpint5 - 1;
					atom[numatom].connum++;
				}
				if (numatom < tmpint5 - 1) {
					if (overflow_flag == 0) {
						bond[numbond].bondi = numatom;
						bond[numbond].bondj = tmpint5 - 1;
						bond[numbond].type = 1;
						if (tmpchar8[0] == 'd')
							bond[numbond].type = 2;
						if (tmpchar8[0] == 't')
							bond[numbond].type = 3;
					}
					numbond++;
					if (numbond >= cinfo.maxbond && overflow_flag == 0) {
						printf
							("\nInfo: the bond number exceeds the MAXBOND, reallocate memory automatically");
						overflow_flag = 1;
					}
				}
			}
			if (tmpint6 > 0) {
				if (overflow_flag == 0) {
					atom[numatom].con[3] = tmpint6 - 1;
					atom[numatom].connum++;
				}
				if (numatom < tmpint6 - 1) {
					if (overflow_flag == 0) {
						bond[numbond].bondi = numatom;
						bond[numbond].bondj = tmpint6 - 1;
						bond[numbond].type = 1;
						if (tmpchar9[0] == 'd')
							bond[numbond].type = 2;
						if (tmpchar9[0] == 't')
							bond[numbond].type = 3;
					}
					numbond++;
					if (numbond >= cinfo.maxbond && overflow_flag == 0) {
						printf
							("\nInfo: the bond number exceeds the MAXBOND, reallocate memory automatically");
						overflow_flag = 1;
					}
				}
			}
			if (tmpint7 > 0) {
				if (overflow_flag == 0) {
					atom[numatom].con[4] = tmpint7 - 1;
					atom[numatom].connum++;
				}
				if (numatom < tmpint7 - 1) {
					if (overflow_flag == 0) {
						bond[numbond].bondi = numatom;
						bond[numbond].bondj = tmpint7 - 1;
						bond[numbond].type = 1;
						if (tmpchar10[0] == 'd')
							bond[numbond].type = 2;
						if (tmpchar10[0] == 't')
							bond[numbond].type = 3;
					}
					numbond++;
					if (numbond >= cinfo.maxbond && overflow_flag == 0) {
						printf
							("\nInfo: the bond number exceeds the MAXBOND, reallocate memory automatically");
						overflow_flag = 1;
					}
				}
			}
			if (tmpint8 > 0) {
				if (overflow_flag == 0) {
					atom[numatom].con[5] = tmpint8 - 1;
					atom[numatom].connum++;
				}
				if (numatom < tmpint8 - 1) {
					if (overflow_flag == 0) {
						bond[numbond].bondi = numatom;
						bond[numbond].bondj = tmpint8 - 1;
						bond[numbond].type = 1;
						if (tmpchar11[0] == 'd')
							bond[numbond].type = 2;
						if (tmpchar11[0] == 't')
							bond[numbond].type = 3;
					}
					numbond++;
					if (numbond >= cinfo.maxbond && overflow_flag == 0) {
						printf
							("\nInfo: the bond number exceeds the MAXBOND, reallocate memory automatically");
						overflow_flag = 1;
					}
				}
			}
			numatom++;
			if (numatom >= cinfo.maxbond && overflow_flag == 0) {
				printf
					("\nInfo: the atom number exceeds the MAXATOM, reallocate memory automatically");
				overflow_flag = 1;
			}

		}
	}
	*atomnum = numatom;
	*bondnum = numbond;
	fclose(fpin);

	for (i = 0; i < *atomnum; i++)
		strcpy(atom[i].aa, tmpchar);
	return overflow_flag;
}
void whin(char *filename, int atomnum, ATOM * atom, int bondnum,
		  BOND * bond)
{
        typedef struct {
                char name[6];
        } ATOMNAME;
	
	int i, j, k;
	ATOMNAME *name;
/*	char tmpchar[MAXCHAR]; */
	char type[5];
	FILE *fpout;
	if ((fpout = fopen(filename, "w")) == NULL) {
		fprintf(stdout, "Cannot open a file %s to write in whin(), exit\n", filename);
		return;
	}
/*
        for(i=0;i<strlen(filename);i++) {
                if(filename[i]=='.') break;
                tmpchar[i]=filename[i];
        }
        tmpchar[i]='\0';
*/
        name = (ATOMNAME *) malloc(sizeof(ATOMNAME) * (atomnum +10));
        if (name == NULL) {
                fprintf(stdout, "memory allocation error for *name in whin()\n");
                exit(1);
        }
                                                                                                                                                                                                 
	element(atomnum, atom);
	for (i = 0; i < atomnum; i++)
		strcpy(name[i].name, atom[i].element);
	fprintf(fpout, "%s 1 \"%s\"\n", atom[i].aa, filename);
	for (i = 0; i < atomnum; i++) {
		fprintf(fpout,
				"atom %d - %s ** - %10.5lf%10.5lf%10.5lf%10.5lf %d ",
				i + 1, name[i].name, atom[i].charge, atom[i].x, atom[i].y,
				atom[i].z, atom[i].connum);
		for (j = 0; j < atom[i].connum; j++) {
			fprintf(fpout, "%d ", atom[i].con[j] + 1);
			for (k = 0; k < bondnum; k++)
				if ((bond[k].bondi == i && bond[k].bondj == atom[i].con[j])
					|| (bond[k].bondi == atom[i].con[j]
						&& bond[k].bondj == i)) {
					strcpy(type, "s");
					if (bond[k].type == 1)
						strcpy(type, "s");
					if (bond[k].type == 2)
						strcpy(type, "d");
					if (bond[k].type == 3)
						strcpy(type, "t");
					if (bond[k].type == 7)
						strcpy(type, "s");
					if (bond[k].type == 8)
						strcpy(type, "d");
					if (bond[k].type == 9)
						strcpy(type, "d");
					if (bond[k].type == 10)
						strcpy(type, "s");
					fprintf(fpout, "%s ", type);
					break;

				}
		}
		fprintf(fpout, "\n");
	}
    fprintf(fpout, "endmol 1\n");
	fclose(fpout);
	free(name);
}
