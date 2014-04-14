/* MOPAC OUTPUT */
int rsqmout(char *filename, int *atomnum, ATOM * atom, CONTROLINFO *cinfo,
			MOLINFO *minfo)
{
	/* int i,j; */
	int flag;
	int number = 0;
	int tmpint1, tmpint2;
	int overflow_flag = 0;
	FILE *fpout;
	char tmpchar1[MAXCHAR];
	char tmpchar2[MAXCHAR];
	char tmpchar3[MAXCHAR];
	char tmpchar4[MAXCHAR];
	char tmpchar5[MAXCHAR];
	char tmpchar6[MAXCHAR];
	char tmpchar7[MAXCHAR];
	char line[MAXCHAR];
	double tmpfloat1;
	double tmpfloat2;
	double tmpfloat3;

	if ((fpout = fopen(filename, "r")) == NULL) {
		fprintf(stdout, "Cannot open sqm output file: %s in rsqmout(), exit\n", filename);
		exit(1);
	}
	initial((*cinfo).maxatom, atom, (*minfo).resname);
	flag = 0;

	for (;;) {
		if (fgets(line, MAXCHAR, fpout) == NULL)
			break;

		sscanf(line, "%s%s%s%s%s%s%s", tmpchar1, tmpchar2, tmpchar3,
			   tmpchar4, tmpchar5, tmpchar6, tmpchar7);
		if (strcmp(tmpchar1, "Total") == 0 
			&& strcmp(tmpchar2, "Mulliken") == 0 
			&& strcmp(tmpchar3, "Charge") == 0 
			&&  strcmp(tmpchar4, "=") == 0) {
			sscanf(line, "%s%s%s%s%lf", tmpchar1, tmpchar2, tmpchar3,
				   tmpchar4, &(*minfo).dcharge);
			(*minfo).icharge = (int) (*minfo).dcharge;
			continue;
		}
		if (strcmp(tmpchar1, "Final") == 0 && strcmp(tmpchar2, "Structure") == 0) {
			flag = 1;	
			continue;
		}
		if (flag == 1 && strcmp("QM_NO.", tmpchar2) == 0
			      && strcmp("MM_NO.", tmpchar3) == 0
			      && strcmp("ATOM",   tmpchar4) == 0
			      && strcmp("X",      tmpchar5) == 0
			      && strcmp("Y",      tmpchar6) == 0
			      && strcmp("Z",      tmpchar7) == 0) {
			flag = 2;
			continue;
		}
		if (flag == 2 && spaceline(line) == 1) {
			flag = 3;
			*atomnum = number;
			break;
		}
		if (flag == 2) {
			if (overflow_flag == 0) {
				sscanf(line, "%s%d%d%s%lf%lf%lf", tmpchar1, &tmpint1, &tmpint2, tmpchar2, &tmpfloat1,
					   &tmpfloat2, &tmpfloat3);
				strcpy(atom[number].name, tmpchar2);
				atom[number].x = tmpfloat1;
				atom[number].y = tmpfloat2;
				atom[number].z = tmpfloat3;
				strcpy(atom[number].aa, (*minfo).resname);
			}
			number++;
			if (number >= (*cinfo).maxatom && overflow_flag == 0) {
				printf
					("\nInfo: the atom number exceeds the MAXATOM, reallocate memory automatically");
				overflow_flag = 1;
			}
		}
	}
	fclose(fpout);
	return overflow_flag;
}

void wsqmout(char *filename, int atomnum, ATOM * atom, CONTROLINFO cinfo,
			 MOLINFO minfo)
{
	char tmpchar[MAXCHAR];
	int status1 = 0;
	int status2 = 0;

        wsqmcrt("sqm.in", atomnum, atom, minfo);
        build_exe_path(tmpchar, "sqm", sizeof tmpchar, 1);
        strcat(tmpchar, " -O -i sqm.in -o sqm.out");

	if (cinfo.intstatus == 1 || cinfo.intstatus == 2)
		fprintf(stdout, "\nRunning: %s\n", tmpchar);
	status1 = system(tmpchar);
	if(status1 == 0) {
		strcpy(tmpchar, "mv -f sqm.out ");
		strcat(tmpchar, filename);
		status2 = system(tmpchar);
		if(status2 != 0) {
			fprintf(stdout, "Error: cannot run \"%s\" in wsqmout() of sqmout.c properly, exit\n", tmpchar);
			exit(1);
		}
	}
	if(status1 != 0) {
		fprintf(stdout, "Error: cannot run \"%s\" in wsqmout() of sqmout.c properly, exit\n", tmpchar);
		exit(1);
	}
}

void rsqmout_coord(char *filename, ATOM * atom) {
	/* int i,j; */
	int flag;
	int number = 0;
	int tmpint1, tmpint2;
	FILE *fpout;
	char tmpchar1[MAXCHAR];
	char tmpchar2[MAXCHAR];
	char tmpchar3[MAXCHAR];
	char tmpchar4[MAXCHAR];
	char tmpchar5[MAXCHAR];
	char tmpchar6[MAXCHAR];
	char tmpchar7[MAXCHAR];
	char line[MAXCHAR];
	double tmpfloat1;
	double tmpfloat2;
	double tmpfloat3;

	if ((fpout = fopen(filename, "r")) == NULL) {
		fprintf(stdout, "Cannot open sqm output file: %s in rsqmout(), exit\n", filename);
		exit(1);
	}
	flag = 0;

	for (;;) {
		if (fgets(line, MAXCHAR, fpout) == NULL)
			break;

		sscanf(line, "%s%s%s%s%s%s%s", tmpchar1, tmpchar2, tmpchar3,
			   tmpchar4, tmpchar5, tmpchar6, tmpchar7);
		if (strcmp(tmpchar1, "Final") == 0 && strcmp(tmpchar2, "Structure") == 0) {
			flag = 1;	
			continue;
		}
		if (flag == 1 && strcmp("QM_NO.", tmpchar2) == 0
			      && strcmp("MM_NO.", tmpchar3) == 0
			      && strcmp("ATOM",   tmpchar4) == 0
			      && strcmp("X",      tmpchar5) == 0
			      && strcmp("Y",      tmpchar6) == 0
			      && strcmp("Z",      tmpchar7) == 0) {
			flag = 2;
			continue;
		}
		if (flag == 2 && spaceline(line) == 1) break;
		if (flag == 2) {
			sscanf(line, "%s%d%d%s%lf%lf%lf", tmpchar1, &tmpint1, &tmpint2, tmpchar2, &tmpfloat1, &tmpfloat2, &tmpfloat3);
			atom[number].x = tmpfloat1;
			atom[number].y = tmpfloat2;
			atom[number].z = tmpfloat3;
			number++;
		}
	}
	fclose(fpout);
}
