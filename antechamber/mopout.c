/* MOPAC OUTPUT */
int rmopout(char *filename, int *atomnum, ATOM * atom, CONTROLINFO *cinfo,
			MOLINFO *minfo)
{
	/* int i,j; */
	int index;
	int number = 0;
	int tmpint;
	int overflow_flag = 0;
	FILE *fpout;
	char tmpchar[MAXCHAR];
	char tmpchar1[MAXCHAR];
	char tmpchar2[MAXCHAR];
	char tmpchar3[MAXCHAR];
	char tmpchar4[MAXCHAR];
	char tmpchar5[MAXCHAR];
	char tmpchar6[MAXCHAR];
	char line[MAXCHAR];
	double tmpfloat1;
	double tmpfloat2;
	double tmpfloat3;

	if ((fpout = fopen(filename, "r")) == NULL) {
		fprintf(stdout, "Cannot open mopac output file: %s in rmopout(), exit\n", filename);
		exit(1);
	}
	initial((*cinfo).maxatom, atom, (*minfo).resname);
	index = 0;

	for (;;) {
		if (fgets(line, MAXCHAR, fpout) == NULL)
			break;

		sscanf(line, "%s%s%s%s%s%s", tmpchar1, tmpchar2, tmpchar3,
			   tmpchar4, tmpchar5, tmpchar6);
		if (strcmp(tmpchar2, "CHARGE") == 0 && strcmp(tmpchar3, "ON") == 0
			&& strcmp("SYSTEM", tmpchar4) == 0) {
			sscanf(line, "%s%s%s%s%s%lf", tmpchar1, tmpchar2, tmpchar3,
				   tmpchar4, tmpchar5, &(*minfo).dcharge);
			(*minfo).icharge = (int) (*minfo).dcharge;
			continue;
		}
		if (index == 0 && strcmp("NET", tmpchar1) == 0
			&& strcmp("ATOMIC", tmpchar2) == 0
			&& strcmp("CHARGES", tmpchar3) == 0
			&& strcmp("AND", tmpchar4) == 0
			&& strcmp("DIPOLE", tmpchar5) == 0
			&& strcmp("CONTRIBUTIONS", tmpchar6) == 0) {
			index = 1;
			continue;
		}
		if (index == 1 && strcmp("ATOM", tmpchar1) == 0
			&& strcmp("NO.", tmpchar2) == 0) {
			index = 2;
			continue;
		}
		if (index == 2 && strcmp("CARTESIAN", tmpchar1) == 0
			&& strcmp("COORDINATES", tmpchar2) == 0) {
			index = 3;
			continue;
		}
		if (index == 3 && strcmp("NO.", tmpchar1) == 0
			&& strcmp("ATOM", tmpchar2) == 0) {
			index = 4;
			continue;
		}
		if (index == 4) {
			index = 5;
			number = 0;
			continue;
		}
		if (index == 5 && spaceline(line) == 1) {
			index = 6;
			*atomnum = number;
			break;
		}
		if (index == 5) {
			if (overflow_flag == 0) {
				sscanf(line, "%d%s%lf%lf%lf", &tmpint, tmpchar, &tmpfloat1,
					   &tmpfloat2, &tmpfloat3);
				strcpy(atom[number].name, tmpchar);
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

void wmopout(char *filename, int atomnum, ATOM * atom, CONTROLINFO cinfo,
			 MOLINFO minfo)
{
	char tmpchar[MAXCHAR];
	size_t copied_size;
	int status1 = 0;
	int status2 = 0;
        strcpy(minfo.ekeyword, "AM1 ANALYT MMOK GEO-OK PRECISE");
	wmopcrt("mopac.in", atomnum, atom, minfo);
	copied_size = build_exe_path(tmpchar, "mopac.sh", sizeof tmpchar, 1 );
	if (cinfo.intstatus == 1 || cinfo.intstatus == 2)
		fprintf(stdout, "\nRunning: %s\n", tmpchar);
	status1 = system(tmpchar);
	if(status1 == 0) {
		strcpy(tmpchar, "mv -f mopac.out ");
		strcat(tmpchar, filename);
		status2 = system(tmpchar);
		if(status2 != 0){
			fprintf(stdout, "Error: cannot run \"%s\" of wmopout() in mopout.c properly, exit\n", tmpchar);
			exit(1);
		}
	}
	if(status1 != 0){
		fprintf(stdout, "Error: cannot run \"%s\" of wmopout() in mopout.c properly, exit\n", tmpchar);
		exit(1);
	}
}


void rmopout_coord(char *filename, ATOM * atom) {
	int index;
	int number = 0;
	int tmpint;
	FILE *fpout;
	char tmpchar[MAXCHAR];
	char tmpchar1[MAXCHAR];
	char tmpchar2[MAXCHAR];
	char tmpchar3[MAXCHAR];
	char tmpchar4[MAXCHAR];
	char tmpchar5[MAXCHAR];
	char tmpchar6[MAXCHAR];
	char line[MAXCHAR];
	double tmpfloat1;
	double tmpfloat2;
	double tmpfloat3;

	if ((fpout = fopen(filename, "r")) == NULL) {
		fprintf(stdout, "Cannot open mopac output file: %s in rmopout(), exit\n", filename);
		exit(1);
	}
	index = 0;

	for (;;) {
		if (fgets(line, MAXCHAR, fpout) == NULL)
			break;

		sscanf(line, "%s%s%s%s%s%s", tmpchar1, tmpchar2, tmpchar3,
			   tmpchar4, tmpchar5, tmpchar6);
		if (index == 0 && strcmp("NET", tmpchar1) == 0
			&& strcmp("ATOMIC", tmpchar2) == 0
			&& strcmp("CHARGES", tmpchar3) == 0
			&& strcmp("AND", tmpchar4) == 0
			&& strcmp("DIPOLE", tmpchar5) == 0
			&& strcmp("CONTRIBUTIONS", tmpchar6) == 0) {
			index = 1;
			continue;
		}
		if (index == 1 && strcmp("ATOM", tmpchar1) == 0
			&& strcmp("NO.", tmpchar2) == 0) {
			index = 2;
			continue;
		}
		if (index == 2 && strcmp("CARTESIAN", tmpchar1) == 0
			&& strcmp("COORDINATES", tmpchar2) == 0) {
			index = 3;
			continue;
		}
		if (index == 3 && strcmp("NO.", tmpchar1) == 0
			&& strcmp("ATOM", tmpchar2) == 0) {
			index = 4;
			continue;
		}
		if (index == 4) {
			index = 5;
			number = 0;
			continue;
		}
		if (index == 5 && spaceline(line) == 1) break;
		if (index == 5) {
			sscanf(line, "%d%s%lf%lf%lf", &tmpint, tmpchar, &tmpfloat1, &tmpfloat2, &tmpfloat3);
			atom[number].x = tmpfloat1;
			atom[number].y = tmpfloat2;
			atom[number].z = tmpfloat3;
			number++;
		}
	}
	fclose(fpout);
}
