/* AMBER RST */
int rrst(char *filename, int *atomnum, ATOM atom[], CONTROLINFO cinfo)
{
	int i, j;
	int overflow_flag = 0;
	FILE *fpin;
	char line[MAXCHAR];

/* since rst file only has coordinates information, it can only be the additional file */
	if ((fpin = fopen(filename, "r")) == NULL) {
		fprintf(stdout, "Cannot open rst file: %s, exit\n", filename);
		exit(1);
	}
/* readin the amber rst file */
	i = 0;
	j = 0;
	for (;;) {
		if (fgets(line, MAXCHAR, fpin) == NULL) {
/*       printf("\nFinished reading %s .", filename); */
			break;
		}
		j++;
		if (j == 2)
			sscanf(line, "%d", atomnum);
		if (strncmp(".", &line[4], 1) == 0
			&& strncmp(".", &line[16], 1) == 0
			&& strncmp(".", &line[28], 1) == 0) {
			if (overflow_flag == 0)
				sscanf(&line[1], "%lf%lf%lf%lf%lf%lf", &atom[i].x,
					   &atom[i].y, &atom[i].z, &atom[i + 1].x,
					   &atom[i + 1].y, &atom[i + 1].z);
			i = i + 2;
			if (i >= *atomnum) break;
			if (i >= cinfo.maxatom && overflow_flag == 0) {
				printf
					("\nInfo: the atom number exceeds the MAXATOM, need reallocate memory automatically");
				overflow_flag = 1;
			}

		}
	}
	fclose(fpin);
	return overflow_flag;

}
void wrst(char *filename, int atomnum, ATOM atom[])
{
	int i, j;
	FILE *fpout;
	if ((fpout = fopen(filename, "w")) == NULL) {
		fprintf(stdout, "Cannot open a file %s to write, exit\n", filename);
		exit(1);
	}
	fprintf(fpout, "MOLECULE\n");
	fprintf(fpout, "%5d\n", atomnum);
	i = 0;
	j = 0;
	while (i < atomnum) {
		fprintf(fpout, "%12.7lf%12.7lf%12.7lf", atom[i].x, atom[i].y,
				atom[i].z);
		i++;
		j++;
		if (j == 2) {
			fprintf(fpout, "\n");
			j = 0;
		}
	}
	fclose(fpout);
}
