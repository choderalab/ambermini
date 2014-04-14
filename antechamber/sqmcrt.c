/* divcon CRT */
int rsqmcrt(char *filename, int *atomnum, ATOM atom[], CONTROLINFO cinfo,
                        MOLINFO minfo)
{
        FILE *fpin;
	int read_flag = 0;
        int numatom;
        int overflow_flag = 0;
	int i,j;
	int icharge;
        char line[MAXCHAR];
        char tmpchar1[MAXCHAR];
        char tmpchar2[MAXCHAR];
        double tmpfloat1, tmpfloat2, tmpfloat3;
                
        if ((fpin = fopen(filename, "r")) == NULL) {
                fprintf(stdout, "Cannot open the divcrt file %s, exit\n", filename);
                exit(1);
        }
        initial(cinfo.maxatom, atom, minfo.resname);
        numatom = 0;
        for (;;) {
                if (fgets(line, MAXCHAR, fpin) == NULL) {
/*       printf("\nFinished reading %s file.", cinfo.ifilename); */
                        break;
                }
		sscanf(line, "%s%s", tmpchar1, tmpchar2);
		if(read_flag ==0 && strstr(line, "qmcharge") != 0) 
			for(i=7;i<strlen(line);i++) {
				if(line[i-7] == 'q' && line[i-6] == 'm' && line[i-5]== 'c' &&
				   line[i-4] == 'h' && line[i-3] == 'a' && line[i-2]== 'r' &&
				   line[i-1] == 'g' && line[i]   == 'e') {
						for(j=i;j<strlen(line);j++) 
							if(line[j] == '=') {
								line[j] = ' ';
								sscanf(&line[j], "%d", &icharge);
								minfo.icharge = icharge;
								minfo.dcharge = icharge;
								break;
							}
						break;
					} 
			}

		if(strcmp(tmpchar1, "/") == 0) {
			read_flag = 1;
                        continue;
		}
		if(read_flag == 1 && strlen(line) < 5)
			break;
		if(read_flag == 1) {
                	sscanf(line, "%s%s%lf%lf%lf", tmpchar1, tmpchar2, &tmpfloat1, &tmpfloat2, &tmpfloat3);
                	if (overflow_flag == 0) {
                        	strcpy(atom[numatom].name, tmpchar2);
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
        }
        fclose(fpin);
        *atomnum = numatom;
/* printf("\n atom number is  %5d", *atomnum); */
        return overflow_flag;
}

void wsqmcrt(char *filename, int atomnum, ATOM atom[], MOLINFO minfo)
{
	FILE *fpout;
	int i, nelectrons; 

	if ((fpout = fopen(filename, "w")) == NULL) {
		fprintf(stdout, "Cannot open file %s to write in wsqmcrt(), exit\n", filename);
		exit(1);
	}

    /*  write initial keywords and control parameters:  */
	fprintf(fpout, "Run semi-empirical minimization\n");
	fprintf(fpout, " &qmmm\n" );
	fprintf(fpout, "  %s  qmcharge=%d,\n /\n",
       minfo.ekeyword, minfo.icharge );

	/* element(atomnum, atom); */

	nelectrons = 0;
	for (i = 0; i < atomnum; i++){
		fprintf(fpout, "%4d %5s  %12.4lf  %12.4lf  %12.4lf \n",
			atom[i].atomicnum, atom[i].name, atom[i].x, atom[i].y, atom[i].z);
		nelectrons += atom[i].atomicnum;
	}
	fprintf(fpout, "\n");
	fclose(fpout);
	nelectrons -= minfo.icharge;
	fprintf( stdout, "Total number of electrons: %d; net charge: %d\n",
		nelectrons,minfo.icharge );
    /*  check that the number of electrons is even:   */
	if( nelectrons%2 != 0 ){
		fprintf( stdout, "INFO: Number of electrons is odd: %d\n", nelectrons) ;
		fprintf( stdout, "      Please check the total charge (-nc flag) and spin multiplicity (-m flag)\n"); 
	}
}
