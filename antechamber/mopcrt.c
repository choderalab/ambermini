/* MOPAC CRT */
int rmopcrt(char *filename, int *atomnum, ATOM atom[], CONTROLINFO cinfo,
			MOLINFO minfo)
{
	FILE *fpin;
	int index;
	int numatom;
	int tmpint1, tmpint2, tmpint3;
	int overflow_flag = 0;
	char line[MAXCHAR];
	char tmpchar[20];
	double tmpfloat1, tmpfloat2, tmpfloat3;

	if ((fpin = fopen(filename, "r")) == NULL) {
		fprintf(stdout, "Cannot open file %s to read in rmopcrt(), exit\n", filename);
		exit(1);
	}
	initial(cinfo.maxatom, atom, minfo.resname);
	index = 0;
	numatom = 0;
	for (;;) {
		if (fgets(line, MAXCHAR, fpin) == NULL) {
/*       printf("\nFinished reading %s file.", cinfo.ifilename); */
			break;
		}
		index++;
		if (index < 4)
			continue;
		sscanf(line, "%s%lf%d%lf%d%lf%d", tmpchar, &tmpfloat1, &tmpint1,
			   &tmpfloat2, &tmpint2, &tmpfloat3, &tmpint3);
		if (spaceline(line) == 1)
			break;
		if (overflow_flag == 0) {
			strcpy(atom[numatom].name, tmpchar);
			atom[numatom].x = tmpfloat1;
			atom[numatom].y = tmpfloat2;
			atom[numatom].z = tmpfloat3;
		}
		numatom++;
		if (numatom >= cinfo.maxatom && overflow_flag == 0) {
			printf
				("\nInfo: the atom number exceeds the MAXATOM, reallocate memory automatically");
			overflow_flag = 1;
		}
	}
	fclose(fpin);
	*atomnum = numatom;
/* printf("\n atom number is  %5d", *atomnum); */
	return overflow_flag;
}
void wmopcrt(char *filename, int atomnum, ATOM atom[], MOLINFO minfo)
{
	FILE *fpout;
	int i, nelectrons; 
	char tmpkeyword[MAXCHAR];

	if ((fpout = fopen(filename, "w")) == NULL) {
		fprintf(stdout, "Cannot open file %s to write in wmopcrt(), exit\n", filename);
		exit(1);
	}

	for (i = 0; i < strlen(minfo.ekeyword); i++)
		tmpkeyword[i] = toupper(minfo.ekeyword[i]);
	if (strstr(tmpkeyword, "CHAR") == NULL) {
		fprintf(fpout, "%s", minfo.ekeyword);
		fprintf(fpout, " CHARGE=%d", minfo.icharge);
	}
	if (strstr(tmpkeyword, "DOUBLET") == NULL && minfo.multiplicity == 2)
		fprintf(fpout, " DOUBLET");
	if (strstr(tmpkeyword, "TRIPLET") == NULL && minfo.multiplicity == 3)
		fprintf(fpout, " TRIPLET");
/*
#if 0
	fprintf(fpout, "\ncreated by wmopcrt()\n\n" );
#else
	fprintf(fpout, "\ncreated by wmopcrt()\n");
#endif
*/
	fprintf(fpout, "\ncreated by wmopcrt() for mopac\n\n");
	element(atomnum, atom);
	nelectrons = 0;
	for (i = 0; i < atomnum; i++){
		fprintf(fpout, "%5s%12.4lf  1  %12.4lf  1  %12.4lf  1   \n",
			atom[i].element, atom[i].x, atom[i].y, atom[i].z);
/*  check that the number of electrons is even:   */
		nelectrons += atom[i].atomicnum;
	}
	fprintf(fpout, "\n");
	fclose(fpout);
	nelectrons -= minfo.icharge;
	fprintf( stdout, "Total number of electrons: %d; net charge: %d\n",
		nelectrons,minfo.icharge );
	if( nelectrons%2 != 0 ){
		fprintf( stdout, "INFO: Number of electrons is odd: %d\n", nelectrons) ;
		fprintf( stdout, "      Please check the total charge (-nc flag) and spin multiplicity (-m flag)\n"); 
	}
}
