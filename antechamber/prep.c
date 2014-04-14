/* PREP */
int rprepi(char *filename, int *atomnum, ATOM * atom, int *bondnum, BOND *bond, CONTROLINFO *cinfo,
		   MOLINFO *minfo)
{
	FILE *fpin;
	int i,j,k;
	int overflow_flag = 0;
	int number = 0;
	int number0 = 0;
	int numbond = 0;
	int suc_flag = 0;
	
	ATOM *atm;
	int tmpint1, tmpint2, tmpint3, tmpint4;
	int readindex;
	char tmpchar1[MAXCHAR];
	char tmpchar2[MAXCHAR];
	char tmpchar3[MAXCHAR];
	char tmpchar[MAXCHAR];
/*	char resname[MAXCHAR]; */
	char line[MAXCHAR];
	double tmpfloat1;
	double tmpfloat2;
	double tmpfloat3;

	if ((fpin = fopen(filename, "r")) == NULL) {
		fprintf(stdout, "Cannot open file %s, exit\n", filename);
		exit(1);
	}
	initial((*cinfo).maxatom, atom, (*minfo).resname);
	i = 0;
	readindex = 0;

	atm = (ATOM *) malloc(sizeof(ATOM) * ((*cinfo).maxatom + 3));
	if (atm == NULL) {
		fprintf(stdout, "memory allocation error for *atm in rprep\n");
		exit(1);
	}

	for (;;) {
		if (fgets(line, MAXCHAR, fpin) == NULL) break;
/*
		if(strlen(line) >= 2 && iscntrl(line[strlen(line)-2])) {
			line[strlen(line)-2] = '\n'; 	
			line[strlen(line)-1] = '\0'; 	
		}
*/
		i++;
		if (i > (*cinfo).maxatom + 11) {
			printf
				("\nInfo: the atom number exceeds the MAXATOM, reallocate memory automatically");
			overflow_flag = 1;
		}
		if (i == 5) {
			sscanf(line, "%s", tmpchar);
/*
			(*minfo).resname[0] = tmpchar[0];
			(*minfo).resname[1] = tmpchar[1];
			(*minfo).resname[2] = tmpchar[2];
			(*minfo).resname[3] = '\0';
*/
			strcpy((*minfo).resname, tmpchar);
		}
		if (i == 8)
			readindex = 1;
		if (strncmp(line, "CHARGE", 6) == 0) {
			readindex = 3;
			continue;
		}
		if (strncmp(line, "LOOP", 4) == 0) {
			readindex = 4;
			continue;
		}
		if (spaceline(line) == 1 && readindex == 1)
			readindex = -2;
		if (spaceline(line) == 1 && readindex == 3)
			readindex = -3;
		if (spaceline(line) == 1 && readindex == 4)
			readindex = -4;
		if (readindex == 1) {
			sscanf(line, "%d%s%s%s%d%d%d%lf%lf%lf%lf", &tmpint1, tmpchar1,
				   tmpchar2, tmpchar3, &tmpint2, &tmpint3, &tmpint4,
				   &tmpfloat1, &tmpfloat2, &tmpfloat3, &atm[i - 8].charge);
			number++;
			if (i == 8) {
				atm[0].x = 0.0;
				atm[0].y = 0.0;
				atm[0].z = 0.0;
				if ((tmpchar2[0] == 'd' || tmpchar2[0] == 'D')
					&& (tmpchar2[1] == 'u' || tmpchar2[1] == 'U'))
					number0++;
				continue;
			}
			if (i == 9) {
				atm[1].x = tmpfloat1;
				atm[1].y = 0.0;
				atm[1].z = 0.0;
				if ((tmpchar2[0] == 'd' || tmpchar2[0] == 'D')
					&& (tmpchar2[1] == 'u' || tmpchar2[1] == 'U'))
					number0++;
				continue;
			}
			if (i == 10) {
				atm[2].x = atm[1].x - tmpfloat1 * cos(DEGRAD * tmpfloat2);
				atm[2].y = tmpfloat1 * sin(DEGRAD * tmpfloat2);
				atm[2].z = 0.0;
				if ((tmpchar2[0] == 'd' || tmpchar2[0] == 'D')
					&& (tmpchar2[1] == 'u' || tmpchar2[1] == 'U'))
					number0++;
				continue;
			}
			if (overflow_flag == 0) {
				strcpy(atm[i - 8].name, tmpchar1);
				strcpy(atm[i - 8].ambername, tmpchar2);
				if(i >= 12) {
					bond[numbond].bondi = numbond + 1 ;
					bond[numbond].bondj = tmpint2 - 4 ;
					numbond ++ ;
				}
				rotate0(atm[tmpint4 - 1], atm[tmpint3 - 1],
						atm[tmpint2 - 1], &atm[i - 8], tmpfloat1,
						tmpfloat2, tmpfloat3);
				continue;
			}
		}
		if (readindex == 3 && overflow_flag == 0) {
			sscanf(line, "%lf%lf%lf%lf%lf", &atm[number0].charge,
				   &atm[number0 + 1].charge, &atm[number0 + 2].charge,
				   &atm[number0 + 3].charge, &atm[number0 + 4].charge);
			number0 += 5;
			continue;
		}
		if (readindex == 4) {
			sscanf(line, "%s%s", tmpchar1, tmpchar2); 
			suc_flag = 0;
			for (j = 3; j < number; j++)  
				if(strcmp(atm[j].name, tmpchar1) == 0)	{
					suc_flag = 1;
					break;
				}
			if(suc_flag == 0) {
				fprintf(stdout, "%s in the LOOP section is not a valid atom name\n", tmpchar1); 	
				exit(1);
			}
			for (k = 3; k < number; k++)  
				if(strcmp(atm[k].name, tmpchar2) == 0)	{
					suc_flag = 1;
					break;
				}
			if(suc_flag == 0) {
				fprintf(stdout, "%s in the LOOP section is not a valid atom name\n", tmpchar2); 	
				exit(1);
			}
			bond[numbond].bondi = j - 3;
			bond[numbond].bondj = k - 3;
			numbond ++;
			continue;
		}
	}
	fclose(fpin);
	if (overflow_flag == 0)
		for (i = 3; i < number; i++) {
			atom[i - 3].x = atm[i].x;
			atom[i - 3].y = atm[i].y;
			atom[i - 3].z = atm[i].z;
			atom[i - 3].charge = atm[i].charge;
			strcpy(atom[i - 3].name, atm[i].name);
			strcpy(atom[i - 3].ambername, atm[i].ambername);
			if ((*cinfo).rnindex == 0)
				strcpy(atom[i - 3].aa, (*minfo).resname);
		}
	free(atm);
	if (overflow_flag == 0)
		for (i = 0; i < numbond; i++) {
			j = bond[i].bondi;
			k = bond[i].bondj;
			atom[j].con[atom[j].connum++] = k;
			atom[k].con[atom[k].connum++] = j;
		}
	*atomnum = number - 3;
	*bondnum = numbond;
	return overflow_flag;
}

int rprepc(char *filename, int *atomnum, ATOM * atom, CONTROLINFO *cinfo,
		   MOLINFO *minfo)
{
	FILE *fpin;
	int i;
	int overflow_flag = 0;
	int number = 0;
	int number0 = 0;
	ATOM *atm;
	int tmpint1;
	int readindex;
	char tmpchar1[MAXCHAR];
	char tmpchar2[MAXCHAR];
	char tmpchar3[MAXCHAR];
	char tmpchar[MAXCHAR];
/*	char resname[MAXCHAR]; */
	char line[MAXCHAR];
	double tmpfloat1;
	double tmpfloat2;
	double tmpfloat3;

	if ((fpin = fopen(filename, "r")) == NULL) {
		fprintf(stdout, "Cannot open file %s, exit\n", filename);
		exit(1);
	}
	atm = (ATOM *) malloc(sizeof(ATOM) * ((*cinfo).maxatom + 3));
	if (atm == NULL) {
		fprintf(stdout, "memory allocation error for *atm in rprep\n");
		exit(1);
	}

	initial((*cinfo).maxatom, atom, (*minfo).resname);
	i = 0;
	readindex = 0;
	for (;;) {
		if (fgets(line, 150, fpin) == NULL)
			break;
		i++;
		if (i > (*cinfo).maxatom + 11) {
			printf
				("\nInfo: the atom number exceeds the MAXATOM, reallocate memory automatically");
			overflow_flag = 1;
		}
		if (i == 5) {
			sscanf(line, "%s", tmpchar);
/*
			(*minfo).resname[0] = tmpchar[0];
			(*minfo).resname[1] = tmpchar[1];
			(*minfo).resname[2] = tmpchar[2];
			(*minfo).resname[3] = '\0';
*/
			strcpy((*minfo).resname, tmpchar) ;
		}
		if (i == 8)
			readindex = 1;
		if (strncmp(line, "CHARGE", 6) == 0) {
			readindex = 3;
			continue;
		}
		if (spaceline(line) == 1 && readindex == 1)
			readindex = 2;
		if (spaceline(line) == 1 && readindex == 3)
			readindex = 4;
		if (readindex == 1) {
			number++;
			if (overflow_flag == 0) {
				sscanf(line, "%d%s%s%s%lf%lf%lf%lf", &tmpint1, tmpchar1,
					   tmpchar2, tmpchar3, &tmpfloat1, &tmpfloat2,
					   &tmpfloat3, &atm[i - 8].charge);
				strcpy(atm[i - 8].name, tmpchar1);
				strcpy(atm[i - 8].ambername, tmpchar2);
				atm[i - 8].x = tmpfloat1;
				atm[i - 8].y = tmpfloat2;
				atm[i - 8].z = tmpfloat3;
				continue;
			}
		}
		if (readindex == 3 && overflow_flag == 0) {
			sscanf(line, "%lf%lf%lf%lf%lf", &atm[number0].charge,
				   &atm[number0 + 1].charge, &atm[number0 + 2].charge,
				   &atm[number0 + 3].charge, &atm[number0 + 4].charge);
			number0 += 5;
			continue;
		}
	}
	fclose(fpin);
	if (overflow_flag == 0)
		for (i = 3; i < number; i++) {
			atom[i - 3].x = atm[i].x;
			atom[i - 3].y = atm[i].y;
			atom[i - 3].z = atm[i].z;
			atom[i - 3].charge = atm[i].charge;
			strcpy(atom[i - 3].name, atm[i].name);
			strcpy(atom[i - 3].ambername, atm[i].ambername);
			if ((*cinfo).rnindex == 0)
				strcpy(atom[i - 3].aa, (*minfo).resname);
		}
	*atomnum = number - 3;
	return overflow_flag;
}

void wprep(char *filename, char *ifilename, int atomnum, ATOM * atom,
		   int bondnum, BOND * bond, CONTROLINFO cinfo, MOLINFO *minfo,
		   int index)
{
	char tmpchar[MAXCHAR];
	size_t copied_size;
	int status = 0;

	wac("ANTECHAMBER_PREP.AC0", atomnum, atom, bondnum, bond, cinfo,
		*minfo);
	status = system("cp -rf ANTECHAMBER_PREP.AC0 ANTECHAMBER_PREP.AC");
        if(status != 0) {
                fprintf(stdout, "Error: cannot run \"%s\" in wprep() of prep.c properly, exit\n", "cp -rf ANTECHAMBER_PREP.AC0 ANTECHAMBER_PREP.AC");
                exit(1);
        }

/*part1: if intype is not prepi, prepc and ac, judge atom type*/
/*
	if ((strcmp(cinfo.intype, "prepi") != 0 &&
		 strcmp(cinfo.intype, "prepc") != 0 &&
		 strcmp(cinfo.intype, "5") != 0 &&
		 strcmp(cinfo.intype, "6") != 0 &&
		 strcmp(cinfo.intype, "ac") != 0 &&
		 strcmp(cinfo.intype, "1") != 0) || cinfo.prediction_index == 1
		|| cinfo.prediction_index == 3) {
                copied_size = build_exe_path(tmpchar, "atomtype",
			sizeof tmpchar, 1);
                strncat(tmpchar, " -i ANTECHAMBER_PREP.AC0 -o ANTECHAMBER_PREP.AC"
			" -p ", sizeof tmpchar - copied_size);
                strncat(tmpchar, minfo->atom_type_def, sizeof tmpchar - strlen(tmpchar) );

		if (cinfo.intstatus == 2)
			fprintf(stdout, "\nRunning: %s\n", tmpchar);
		status = system(tmpchar);
        	if(status != 0) {
                	fprintf(stdout, "Error: cannot run \"%s\" in wprep() of prep.c properly, exit\n", tmpchar);
                	exit(1);
        	}
	}
*/
/*part2: if intype isn't prepi and the outype is prepi*/
/*
	if (strcmp(cinfo.intype, "prepi") != 0 &&
		strcmp(cinfo.intype, "5") != 0 && index == 1) {
*/
	if (index == 1) {
		copied_size = build_exe_path(tmpchar, "prepgen",
			sizeof tmpchar, 1);
		strncat(tmpchar, " -i ANTECHAMBER_PREP.AC -f int -o " ,
			sizeof tmpchar - copied_size);
		strncat(tmpchar, filename, sizeof(tmpchar) - strlen(tmpchar) -1);
		strcat(tmpchar, " -rn \"");
		strcat(tmpchar, atom[0].aa);
		strcat(tmpchar, "\" -rf ");
		strcat(tmpchar, minfo->resfilename);
		if (cinfo.intstatus == 2)
			fprintf(stdout, "\nRunning: %s\n", tmpchar);

		status = system(tmpchar);
                if(status != 0) {
                        fprintf(stdout, "Error: cannot run \"%s\" in wprep() of prep.c properly, exit\n", tmpchar);
                        exit(1);
                }
		return;

/*  system("rm -f ANTECHAMBER_PREP.AC");*/
	}

/*part3: if intype isn't prepc and the outype is prepc*/
/*
	if (strcmp(cinfo.intype, "prepc") != 0 &&
		strcmp(cinfo.intype, "6") != 0 && index == 0) {
*/
	if (index == 0) {
                copied_size = build_exe_path(tmpchar, "prepgen",
			sizeof tmpchar, 1);
		strncat(tmpchar, " -i ANTECHAMBER_PREP.AC -f car -o ",
			sizeof tmpchar - copied_size);
		strncat(tmpchar, filename, sizeof(tmpchar) - strlen(tmpchar) -1);
		strcat(tmpchar, " -rn \"");
		strcat(tmpchar, atom[0].aa);
		strcat(tmpchar, "\" -rf ");
		strcat(tmpchar, minfo->resfilename);
		if (cinfo.intstatus == 2)
			fprintf(stdout, "\nRunning: %s\n", tmpchar);
		status = system(tmpchar);
                if(status != 0) {
                        fprintf(stdout, "Error: cannot run \"%s\" in wprep() of prep.c properly, exit\n", tmpchar);
                        exit(1);
                }
		return;
/*  system("rm -f ANTECHAMBER_PREP.AC");*/
	}

/*part3: if intype is prepi and the outype is prepi*/
/*
	if ((strcmp(cinfo.intype, "prepi") == 0
		 || strcmp(cinfo.intype, "5") == 0)
		&& index == 1) {

		if ((fpin = fopen(ifilename, "r")) == NULL) {
			fprintf(stdout, "Cannot open file %s, exit\n", ifilename);
			exit(1);
		}
		if ((fpout = fopen(filename, "w")) == NULL) {
			fprintf(stdout, "Cannot open file %s, exit\n", filename);
			exit(1);
		}
		if (cinfo.prediction_index == 1 || cinfo.prediction_index == 3)
			rac("ANTECHAMBER_PREP.AC", &atomnum, atom, &bondnum, bond,
				&cinfo, minfo);
		chargeindex = 0;
		typeindex = 0;
		num = 0;
		for (;;) {
			if (fgets(line, MAXCHAR, fpin) == NULL)
				break;
			num++;
			if (num == 11
				&& (cinfo.prediction_index == 1
					|| cinfo.prediction_index == 3))
				typeindex = 1;
			if (strncmp(line, "CHARGE", 6) == 0) {
				fprintf(fpout, "%s", line);
				chargeindex = 1;
				continue;
			}
			if (chargeindex == 0 && typeindex == 0)
				fprintf(fpout, "%s", line);
			if (chargeindex == 2 && spaceline(line) == 1) {
				fprintf(fpout, "%s", line);
				chargeindex = 0;
				continue;
			}
			if (typeindex == 1 && spaceline(line) == 1) {
				fprintf(fpout, "%s", line);
				typeindex = 0;
				continue;
			}
			if (typeindex == 1) {
				sscanf(line, "%d%s%s%s%d%d%d%lf%lf%lf%lf", &tmpint1,
					   tmpchar1, tmpchar2, tmpchar3, &tmpint2, &tmpint3,
					   &tmpint4, &tmpfloat1, &tmpfloat2, &tmpfloat3,
					   &charge);
				fprintf(fpout, "%4d  %-5s %-5s %s", tmpint1, tmpchar1,
						atom[num - 11].ambername, tmpchar3);
				fprintf(fpout, "%5d%4d%4d%10.3f%10.3f%10.3f%10.6f\n",
						tmpint2, tmpint3, tmpint4, tmpfloat1, tmpfloat2,
						tmpfloat3, charge);
			}
			if (chargeindex == 1) {
				changelineindex = 0;
				for (i = 0; i < atomnum; i++) {
					fprintf(fpout, "%9.6f", atom[i].charge);
					changelineindex++;
					if (changelineindex == 5) {
						fprintf(fpout, "\n");
						changelineindex = 0;
					}
				}
				chargeindex = 2;
				continue;
			}
		}
		fclose(fpin);
		fclose(fpout);
	}
	if ((strcmp(cinfo.intype, "prepc") == 0
		 || strcmp(cinfo.intype, "6") == 0)
		&& index == 0) {
		if ((fpin = fopen(ifilename, "r")) == NULL) {
			fprintf(stdout, "Cannot open file %s, exit\n", ifilename);
			exit(1);
		}
		if ((fpout = fopen(filename, "w")) == NULL) {
			fprintf(stdout, "Cannot open file %s, exit\n", filename);
			exit(1);
		}
		chargeindex = 0;
		typeindex = 0;
		num = 0;
		if (cinfo.prediction_index == 1 || cinfo.prediction_index == 3)
			rac("ANTECHAMBER_PREP.AC", &atomnum, atom, &bondnum, bond,
				&cinfo, minfo);
		for (;;) {
			if (fgets(line, MAXCHAR, fpin) == NULL)
				break;
			num++;
			if (num == 11
				&& (cinfo.prediction_index == 1
					|| cinfo.prediction_index == 3))
				typeindex = 1;
			if (strncmp(line, "CHARGE", 6) == 0) {
				fprintf(fpout, "%s", line);
				chargeindex = 1;
				continue;
			}
			if (chargeindex == 0 && typeindex == 0)
				fprintf(fpout, "%s", line);
			if (chargeindex == 2 && spaceline(line) == 1) {
				fprintf(fpout, "%s", line);
				chargeindex = 0;
				continue;
			}
			if (typeindex == 1 && spaceline(line) == 1) {
				fprintf(fpout, "%s", line);
				typeindex = 0;
				continue;
			}
			if (typeindex == 1) {
				sscanf(line, "%d%s%s%s%lf%lf%lf%lf", &tmpint1, tmpchar1,
					   tmpchar2, tmpchar3, &tmpfloat1, &tmpfloat2,
					   &tmpfloat3, &charge);

				fprintf(fpout, "%4d  %-5s %-5s %s", tmpint1, tmpchar1,
						atom[num - 11].ambername, tmpchar3);
				fprintf(fpout, "        %10.6f%12.6f%12.6f%12.6f\n",
						tmpfloat1, tmpfloat2, tmpfloat3, charge);

			}
			if (chargeindex == 1) {
				changelineindex = 0;
				for (i = 0; i < atomnum; i++) {
					fprintf(fpout, "%9.6f", atom[i].charge);
					changelineindex++;
					if (changelineindex == 5) {
						fprintf(fpout, "\n");
						changelineindex = 0;
					}
				}
				chargeindex = 2;
			}
		}
		fclose(fpin);
		fclose(fpout);
	}
*/

}
