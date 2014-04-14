/* CSD */
int rcsd(char *filename, int *atomnum, ATOM atom[], CONTROLINFO cinfo,
		 MOLINFO minfo)
{
	int i;
	int index = 0;
	int overflow_flag = 0;
	int tmpint;
	int tmpint0;
	int tmpint1;
	int tmpint2;
	int tmpint3;
	int tmpint4;
	int tmpint5;
	int tmpint6;
	char line[MAXCHAR];

	FILE *fpin;

	if ((fpin = fopen(filename, "r")) == NULL) {
		fprintf(stdout, "Cannot open the csd file %s to read, exit\n", filename);
		exit(1);
	}
	initial(cinfo.maxatom, atom, minfo.resname);
	i = 0;
	for (;;) {
		if (fgets(line, MAXCHAR, fpin) == NULL) {
/*       printf("\nFinished reading %s file.", cinfo.ifilename); */
			break;
		}
		index++;
		if (index == 3) {
			sscanf(line, "%d", &tmpint);
			continue;
		}
		tmpint1 = -1;
		tmpint2 = -1;
		tmpint3 = -1;
		tmpint4 = -1;
		tmpint5 = -1;
		tmpint6 = -1;

		if (index >= 4) {
			if (overflow_flag == 0) {
				sscanf(line, "%d%s%lf%lf%lf%d%d%d%d%d%d", &tmpint0,
					   atom[i].name, &atom[i].x, &atom[i].y, &atom[i].z,
					   &tmpint1, &tmpint2, &tmpint3, &tmpint4, &tmpint5,
					   &tmpint6);
				if (tmpint1 > 0) {
					atom[i].con[0] = tmpint1 - 1;
					atom[i].connum++;
				}
				if (tmpint2 > 0) {
					atom[i].con[1] = tmpint2 - 1;
					atom[i].connum++;
				}
				if (tmpint3 > 0) {
					atom[i].con[2] = tmpint3 - 1;
					atom[i].connum++;
				}
				if (tmpint4 > 0) {
					atom[i].con[3] = tmpint4 - 1;
					atom[i].connum++;
				}
				if (tmpint5 > 0) {
					atom[i].con[4] = tmpint5 - 1;
					atom[i].connum++;
				}
				if (tmpint6 > 0) {
					atom[i].con[5] = tmpint6 - 1;
					atom[i].connum++;
				}
			}
			i++;
			if (i >= cinfo.maxatom && overflow_flag == 0) {
				printf
					("\nInfo: the atom number exceeds the MAXATOM, reallocate memory automatically");
				overflow_flag = 1;
			}
		}
	}
	*atomnum = i;
/*  printf("\n The atomic number is %5d", *atomnum); */
	fclose(fpin);
/*   for(i=0;i<*atomnum;i++)
   printf("\n%5d%5s%9.3lf%9.3lf%9.3lf%5d%5d%5d%5d%5d%5d", i+1, atom[i].name,
   atom[i].x, atom[i].y, atom[i].z,atom[i].con[0], atom[i].con[1], atom[i].con[2], atom[i].con[3], atom[i].con[4],atom[i].con[5]);
   */
	return overflow_flag;
}
void wcsd(char *filename, int atomnum, ATOM atom[])
{
	int i;
	FILE *fpout;
	if ((fpout = fopen(filename, "w")) == NULL) {
		fprintf(stdout, "Cannot open file %s to write in wcsd(), exit\n", filename);
		exit(1);
	}
	fprintf(fpout,
			"REFRENCE STRUCTURE = 00000    A,B,C, =   1.000   1.000   1.000\n");
	fprintf(fpout,
			"ALPHA,BETA,GAMMA = 90.000   90.000   90.000   SPGR =   P1\n");
	fprintf(fpout, "%4d\n", atomnum);
	for (i = 0; i < atomnum; i++) {
		fprintf(fpout, "%4d %-4s%10.5lf%10.5lf%10.5lf", i + 1,
				atom[i].name, atom[i].x, atom[i].y, atom[i].z);
		if (atom[i].con[0] > -1)
			fprintf(fpout, "%5d", atom[i].con[0] + 1);
		if (atom[i].con[1] > -1)
			fprintf(fpout, "%4d", atom[i].con[1] + 1);
		if (atom[i].con[2] > -1)
			fprintf(fpout, "%4d", atom[i].con[2] + 1);
		if (atom[i].con[3] > -1)
			fprintf(fpout, "%4d", atom[i].con[3] + 1);
		if (atom[i].con[4] > -1)
			fprintf(fpout, "%4d", atom[i].con[4] + 1);
		if (atom[i].con[5] > -1)
			fprintf(fpout, "%4d", atom[i].con[5] + 1);
		fprintf(fpout, "\n");
	}
	fclose(fpout);
}
