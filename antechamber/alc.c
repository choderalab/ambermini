/* ALC */
int ralc(char *filename, int *atomnum, ATOM * atom, int *bondnum,
		 BOND * bond, CONTROLINFO cinfo, MOLINFO minfo)
{
	int i, j;
	int index = 0;
	int overflow_flag = 0;
	int tmpint1;
	int tmpint2;
	int tmpint3;
	int tmpint4;
	int tmpint5;
	int itype;
	char type[20];
	char line[MAXCHAR];
/* double charge;*/
	char tmpchar1[30], tmpchar2[30], tmpchar3[30];
	FILE *fpin;

	if ((fpin = fopen(filename, "r")) == NULL) {
		fprintf(stdout, "Cannot open the alc input file: %s, exit\n", filename);
		exit(1);
	}
	initial(cinfo.maxatom, atom, minfo.resname);
	i = 0;
	j = 0;
	for (;;) {
		if (fgets(line, MAXCHAR, fpin) == NULL) {
/*       printf("\nFinished reading %s file.", cinfo.ifilename); */
			break;
		}
		index++;
		if (index == 1) {
			sscanf(line, "%d%s%d%s%d%s", &tmpint1, tmpchar1, &tmpint2,
				   tmpchar2, &tmpint3, tmpchar3);
			continue;
		}
		if (index > 1 && index < tmpint1 + 2) {
			if (overflow_flag == 0)
				sscanf(line, "%d%s%lf%lf%lf%lf", &tmpint3, atom[i].name,
					   &atom[i].x, &atom[i].y, &atom[i].z,
					   &atom[i].charge);
			i++;
			if (i >= cinfo.maxatom && overflow_flag == 0) {
				printf
					("\nInfo: the atom number exceeds the MAXATOM, reallocate memory automatically");
				overflow_flag = 1;
			}
		}
		if (index >= tmpint1 + 2 && index < tmpint1 + tmpint2 + 2) {
			sscanf(line, "%d%d%d%s", &tmpint3, &tmpint4, &tmpint5, type);
			if (overflow_flag == 0) {
				itype = 0;
				atom[tmpint4 - 1].con[atom[tmpint4 - 1].connum++] =
					tmpint5 - 1;
				atom[tmpint5 - 1].con[atom[tmpint5 - 1].connum++] =
					tmpint4 - 1;
				if (strcmp(type, "1") == 0)
					itype = 1;
				if (strcmp(type, "2") == 0)
					itype = 2;
				if (strcmp(type, "3") == 0)
					itype = 3;
				if (strcmp(type, "am") == 0)
					itype = 1;
				if (strcmp(type, "ar") == 0)
					itype = 10;
				if (strcmp(type, "SINGLE") == 0)
					itype = 1;
				if (strcmp(type, "DOUBLE") == 0)
					itype = 2;
				if (strcmp(type, "TRIPLE") == 0)
					itype = 3;
				bond[j].bondi = tmpint4 - 1;
				bond[j].bondj = tmpint5 - 1;
				bond[j].type = itype;
			}
			j++;
			if (j >= cinfo.maxbond && overflow_flag == 0) {
				printf
					("\nInfo: the bond number exceeds the MAXBOND, reallocate memory automatically");
				overflow_flag = 1;
			}
		}
	}
	*atomnum = i;
	*bondnum = j;
/*  printf("\n The atomic number is %5d", *atomnum); */
	fclose(fpin);
	return overflow_flag;
}
void walc(char *filename, int atomnum, ATOM * atom, int bondnum,
		  BOND * bond)
{
        typedef struct {
                char name[6];
        } ATOMNAME;
        ATOMNAME *name;
	int i;
	/* int breakindex; */
	char type[10];
	FILE *fpout;

	if ((fpout = fopen(filename, "w")) == NULL) {
		fprintf(stdout, "Cannot open a file (%s) to write in alc format, exit\n", filename);
		return;
	}
	name = (ATOMNAME *) malloc(sizeof(ATOMNAME) * (atomnum +10));
        if (name == NULL) {
                fprintf(stdout, "memory allocation error for *name in walc()\n");
                exit(1);
        }

	element(atomnum, atom);
	for (i = 0; i < atomnum; i++)
		strcpy(name[i].name, atom[i].element);

	fprintf(fpout, "%5d ATOMS, %5d BONDS, %5.1lf CHARGE\n", atomnum,
			bondnum, 0.0);
	for (i = 0; i < atomnum; i++)
		fprintf(fpout, "%5d %-5s%9.4lf%9.4lf%9.4lf%11.4lf\n",
				i + 1, name[i].name, atom[i].x, atom[i].y, atom[i].z,
				atom[i].charge);
/*
	for (i = 0; i < atomnum; i++)
		for (j = i + 1; j < atomnum; j++)
			for (k = 0; k < bondnum; k++)
				if ((bond[k].bondi == i && bond[k].bondj == j) ||
					(bond[k].bondi == j && bond[k].bondj == i)) {
					fprintf(fpout, "%5d%6d%6d", k + 1, bond[k].bondi + 1,
							bond[k].bondj + 1);
					strcpy(type, "0");
					if (bond[k].type == 1)
						strcpy(type, "SINGLE");
					if (bond[k].type == 2)
						strcpy(type, "DOUBLE");
					if (bond[k].type == 3)
						strcpy(type, "TRIPLE");
					if (bond[k].type == 7)
						strcpy(type, "SINGLE");
					if (bond[k].type == 8)
						strcpy(type, "DOUBLE");
					if (bond[k].type == 9)
						strcpy(type, "DOUBLE");
					if (bond[k].type == 10)
						strcpy(type, "SINGLE");
					fprintf(fpout, "%8s\n", type);
				}
*/
	for (i = 0; i < bondnum; i++) {
		strcpy(type, "0");
		if (bond[i].type == 1)
			strcpy(type, "SINGLE");
		if (bond[i].type == 2)
			strcpy(type, "DOUBLE");
		if (bond[i].type == 3)
			strcpy(type, "TRIPLE");
		if (bond[i].type == 7)
			strcpy(type, "SINGLE");
		if (bond[i].type == 8)
			strcpy(type, "DOUBLE");
		if (bond[i].type == 9)
			strcpy(type, "DOUBLE");
		if (bond[i].type == 10) 
			strcpy(type, "SINGLE");
		fprintf(fpout, "%5d%6d%6d%8s\n", i + 1, bond[i].bondi + 1,
				bond[i].bondj + 1, type);
	}
	fclose(fpout);
	free(name);
}
