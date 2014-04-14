/* MOL2 */
int rmol2(char *filename, int *atomnum, ATOM atom[], int *bondnum,
		  BOND bond[], CONTROLINFO *cinfo, MOLINFO *minfo, int flag)
{
/*if flag =1, read in atom type, if flag ==0, do not read in atom type */
	int i;
	int numatom;
	int numbond;
	int index = 0;
	int atomrecord = 0;
	int bondrecord = 0;
	int tmpint1;
	int tmpint2;
	int tmpint3;
	int tmpint4;
	int tmpint5;
	int mf1 = 1;
	int itype = 1;
	int overflow_flag = 0;
	int read_atomnum;
	int read_bondnum;
	int len;
	double tmpf1, tmpf2, tmpf3, tmpf4;
	char type[MAXCHAR];
	char tmpchar1[MAXCHAR], tmpchar2[MAXCHAR], tmpchar3[MAXCHAR],
		tmpchar4[MAXCHAR], tmpchar5[MAXCHAR];
	char tmpc1[MAXCHAR];
	char line[MAXCHAR];
	FILE *fpin;

	if ((fpin = fopen(filename, "r")) == NULL) {
		fprintf(stdout, "Cannot open file %s to read in rmol2(), exit\n", filename);
		exit(1);
	}
	
	initial((*cinfo).maxatom, atom, (*minfo).resname);
	numatom = 0;
	numbond = 0;
	for (;;) {
		if (fgets(line, MAXCHAR, fpin) == NULL) break;

		sscanf(line, "%s%s%s", tmpchar1, tmpchar2, tmpchar3);
		if (strcmp("@<TRIPOS>MOLECULE", tmpchar1) == 0) {
			index = 1;
			continue;
		}
		if (strcmp("@<TRIPOS>ATOM", tmpchar1) == 0) {
			atomrecord = 1;
			continue;
		}
		if (strcmp("@<TRIPOS>BOND", tmpchar1) == 0) {
			bondrecord = 1;
			atomrecord = 0;
			continue;
		}
		if (bondrecord == 1 && strncmp("@<TRIPOS>", tmpchar1, 9) == 0) { /* end of bonds */
			bondrecord = 0;
			continue;
		}
		if (index == 1) {
			if((*cinfo).rnindex != 1) 
				strcpy((*minfo).longresname, tmpchar1);
			index = 2;
			continue;
		}
		if (index == 2) {
			sscanf(line, "%d%d", &read_atomnum, &read_bondnum);
 			if (read_atomnum > (*cinfo).maxatom || read_bondnum > (*cinfo).maxbond) {
			        *atomnum = read_atomnum;
       				*bondnum = read_bondnum;
				return 1;
			}
			index = -1;
			continue;
		}

		/*  process the ATOM section:    */

		if (atomrecord) {
			strcpy(tmpchar5, "MOL");
			sscanf(line, "%d%s%lf%lf%lf%s%d%s%lf", &tmpint1, tmpc1,
				   &tmpf1, &tmpf2, &tmpf3, tmpchar4,
				   &tmpint2, tmpchar5, &tmpf4);
			if( numatom==0 ) mf1 = tmpint1;  /* atom number of first atom */
			if (overflow_flag == 0) {
				if (strlen(tmpc1) > 5) {
					atom[numatom].name[0] = tmpc1[0];
					atom[numatom].name[1] = tmpc1[1];
					atom[numatom].name[2] = tmpc1[2];
					atom[numatom].name[3] = tmpc1[3];
					atom[numatom].name[4] = tmpc1[4];
					atom[numatom].name[5] = '\0';
				} else
					strcpy(atom[numatom].name, tmpc1);
				if (tmpc1[0] == '*') {
					for (i = 0; i < strlen(tmpchar4); i++) {
						if (tmpchar4[i] == '.')
							break;
						if (i > 3)
							break;
						atom[numatom].name[i] = tmpchar4[i];
					}
					atom[numatom].name[i] = '\0';
				}
				if (strlen(tmpc1) > 5) {
					atom[numatom].ambername[0] = tmpchar4[0];
					atom[numatom].ambername[1] = tmpchar4[1];
					atom[numatom].ambername[2] = tmpchar4[2];
					atom[numatom].ambername[3] = tmpchar4[3];
					atom[numatom].ambername[4] = tmpchar4[4];
					atom[numatom].ambername[5] = '\0';
				} else
					strcpy(atom[numatom].ambername, tmpchar4);

				atom[numatom].x = tmpf1;
				atom[numatom].y = tmpf2;
				atom[numatom].z = tmpf3;
				atom[numatom].charge = tmpf4;
				atom[numatom].resno = tmpint2;
				if (flag == 1){
					if (strlen(tmpchar4) > 5) {
						atom[numatom].ambername[0] = tmpchar4[0];
						atom[numatom].ambername[1] = tmpchar4[1];
						atom[numatom].ambername[2] = tmpchar4[2];
						atom[numatom].ambername[3] = tmpchar4[3];
						atom[numatom].ambername[4] = tmpchar4[4];
						atom[numatom].ambername[5] = '\0';
					} else {
						strcpy(atom[numatom].ambername, tmpchar4);
					}
				}
				if ((*cinfo).rnindex == 1)
					strcpy(atom[numatom].aa, (*minfo).resname);
				else {
					if (strlen(tmpchar5) > 3) {
						atom[numatom].aa[0] = tmpchar5[0];
						atom[numatom].aa[1] = tmpchar5[1];
						atom[numatom].aa[2] = tmpchar5[2];
						atom[numatom].aa[3] = '\0';
					} else
						strcpy(atom[numatom].aa, tmpchar5);
				}
			}
			numatom++;
			if (numatom >= (*cinfo).maxatom && overflow_flag == 0) {
				printf
					("\nThe atom number exceeds the MAXATOM, reallocate memory");
				overflow_flag = 1;
			}
		}
		/*  process the BOND section:    */

		if (bondrecord) {
			if (numbond >= read_bondnum) continue;	
			sscanf(line, "%d%d%d%s", &tmpint3, &tmpint4, &tmpint5, type);
			if (overflow_flag == 0) {
				bond[numbond].jflag = -1;
				len = strlen(type);
                                if(len >=2 && (type[len-1] == 'f' || type[len-1] == 'F')) {
/* bond type is frozen during bond type assignment */
                                        type[len-1] = '\0'; 
                                        bond[numbond].jflag = 1;
					if(strncmp(type, "ar", 2) ==0 || strncmp(type, "AR", 2) == 0) {
						fprintf(stdout, "You cannot freeze this bond (%d-%d,%s) since it is unclear whether it is a single or double bond, you may change it to '7' (aromatic single) or '8' (aromatic double), exit\n", tmpint4-mf1 + 1, tmpint5-mf1 + 1, type);
						exit(1);
					}
					if(strncmp(type, "un", 2) ==0 || strncmp(type, "UN", 2) == 0) {
						fprintf(stdout, "You cannot freeze this bond (%d-%d,%s) since the bond type is unclear, exit\n", tmpint4-mf1 + 1, tmpint5-mf1 + 1, type);
						exit(1);
					}
                                }

				atom[tmpint4 - mf1].con[atom[tmpint4 - mf1].connum++] =
					tmpint5 - mf1;
				atom[tmpint5 - mf1].con[atom[tmpint5 - mf1].connum++] =
					tmpint4 - mf1;
				if (strcmp(type, "1") == 0)
					itype = 1;
				else if (strcmp(type, "2") == 0)
					itype = 2;
				else if (strcmp(type, "3") == 0)
					itype = 3;
				else if (strcmp(type, "am") == 0)
					itype = 1;
				else if (strcmp(type, "7") == 0)
					itype = 1;
				else if (strcmp(type, "8") == 0)
					itype = 2;
				else if (strcmp(type, "ar") == 0)
					itype = 10;
				else if (strcmp(type, "SINGLE") == 0)
					itype = 1;
				else if (strcmp(type, "DOUBLE") == 0)
					itype = 2;
				else if (strcmp(type, "TRIPLE") == 0)
					itype = 3;
				else if (strcmp(type, "-1") == 0)
					itype = -1;
				else if (strcmp(type, "-2") == 0)
					itype = -2;
				else if (strcmp(type, "-3") == 0)
					itype = -3;
				else if (strcmp(type, "-am") == 0)
					itype = -1;
				else if (strcmp(type, "-ar") == 0)
					itype = -10;
				else if (strcmp(type, "-SINGLE") == 0)
					itype = -1;
				else if (strcmp(type, "-DOUBLE") == 0)
					itype = -2;
				else if (strcmp(type, "-TRIPLE") == 0)
					itype = -3;
				else {
					printf("\nUnknown bond type for BOND(%s-%s:%s), assuming it is a single bond", atom[tmpint4-mf1].name, atom[tmpint5-mf1].name, type);
					itype = 1;
				}
				bond[numbond].bondi = tmpint4 - mf1;
				bond[numbond].bondj = tmpint5 - mf1;
				bond[numbond].type = itype;
			}
			numbond++;
			if (numbond >= (*cinfo).maxbond && overflow_flag == 0) {
				printf
					("\nThe bond number exceeds the MAXBOND, reallocate memory");
				overflow_flag = 1;
			}
		}
	}
	fclose(fpin);
	*atomnum = numatom;
	*bondnum = numbond;
/*
	if(overflow_flag == 0) {
		tmpf1 = 0;
		tmpf2 = 0;
		for(i=0;i<numatom;i++) {
			if(atom[i].charge > 0)
				tmpf1+=atom[i].charge;
			else
				tmpf2+=atom[i].charge;
		}
		printf("\nchargep is %9.4lf", tmpf1);
		printf("\nchargen is %9.4lf", tmpf2);
	}
*/
	return overflow_flag;
}



void wmol2(char *filename, int atomnum, ATOM atom[], int bondnum,
		   BOND bond[], AROM arom[], CONTROLINFO cinfo, MOLINFO minfo) {
	int i,j;
	int flag; 
	int suc;
	int bondi, bondj;
	FILE *fpout;
	if ((fpout = fopen(filename, "w")) == NULL) {
		fprintf(stdout, "Cannot open file %s to write in wmol2(), exit\n", filename);
		return;
	}
	fprintf(fpout, "@<TRIPOS>MOLECULE\n");
	if (strlen(minfo.longresname) <= 1)
		strcpy(minfo.longresname, atom[0].aa);
	fprintf(fpout, "%s\n", minfo.longresname);
	fprintf(fpout, "%5d%6d%6d%6d%6d\n", atomnum, bondnum, 1, 0, 0);
	fprintf(fpout, "SMALL\n");
	if(strlen(cinfo.chargetype)<=1) strcpy(cinfo.chargetype, "No Charge or Current Charge"); 
	fprintf(fpout, "%s\n\n\n", cinfo.chargetype);
	fprintf(fpout, "@<TRIPOS>ATOM\n");
	for (i = 0; i < atomnum; i++) 
		fprintf(fpout,
			"%7d %-8s%10.4lf%10.4lf%10.4lf %-6s%5d %-8s%9.6lf\n",
			i + 1, atom[i].name, atom[i].x, atom[i].y, atom[i].z,
			atom[i].ambername, atom[i].resno, atom[i].aa, atom[i].charge);
	fprintf(fpout, "@<TRIPOS>BOND\n");
	for (i = 0; i < bondnum; i++) {
		fprintf(fpout, "%6d%5d%5d", i + 1, bond[i].bondi + 1, bond[i].bondj + 1);
		bondi = bond[i].bondi;
		bondj = bond[i].bondj;
		if (bond[i].type == 1) {
			if(strcmp(atom[bondi].ambername, "N.am") == 0 && strcmp(atom[bondj].ambername, "C.2") == 0) {
				flag = 0;
                              	for (j = 0; j < 3; j++) 
                                       	if(atom[atom[bondj].con[j]].atomicnum == 8 && atom[atom[bondj].con[j]].connum == 1) {
						fprintf(fpout, " %-4s\n", "am");
						flag = 1;
						break;
					}
				if(flag == 1) continue; 
			}
			if (strcmp(atom[bondj].ambername, "N.am") == 0 && strcmp(atom[bondi].ambername, "C.2") == 0)  {
				flag = 0;
                               	for (j = 0; j < 3; j++) 
                                       	if(atom[atom[bondi].con[j]].atomicnum == 8 && atom[atom[bondi].con[j]].connum == 1) {
						fprintf(fpout, " %-4s\n", "am");
						flag = 1;
						break;
					}
				if(flag == 1) continue;
			}
			fprintf(fpout, " %-4s\n", "1");
			continue;
		}
		if (bond[i].type == 2) {
			if(strcmp(atom[bondi].ambername, "C.2") == 0 && strcmp(atom[bondj].ambername, "O.co2") == 0) 
				fprintf(fpout, " %-4s\n", "1");
			else if(strcmp(atom[bondj].ambername, "C.2") == 0 && strcmp(atom[bondi].ambername, "O.co2") == 0) 
				fprintf(fpout, " %-4s\n", "1");
			else
				fprintf(fpout, " %-4s\n", "2");
			continue;
		}
		if (bond[i].type == 3) {
			fprintf(fpout, " %-4s\n", "3");
			continue;
		}
		if (bond[i].type == 7) {
			suc = 1;
			if(strcmp(atom[bondi].ambername, "C.ar") != 0 && 
                           strcmp(atom[bondi].ambername, "N.ar") != 0 &&
			   strcmp(atom[bondi].ambername, "ca")   != 0 && 
			   strcmp(atom[bondi].ambername, "cp")   != 0 && 
			   strcmp(atom[bondi].ambername, "pb")   != 0 && 
                           strcmp(atom[bondi].ambername, "nb")   != 0) 
				suc = 0;
			if(strcmp(atom[bondj].ambername, "C.ar") != 0 && 
                           strcmp(atom[bondj].ambername, "N.ar") != 0 &&
			   strcmp(atom[bondj].ambername, "ca")   != 0 && 
			   strcmp(atom[bondj].ambername, "cp")   != 0 && 
			   strcmp(atom[bondj].ambername, "pb")   != 0 && 
                           strcmp(atom[bondj].ambername, "nb")   != 0) 
				suc = 0;
			if(suc == 1) 
				fprintf(fpout, " %-4s\n", "ar");
			else 
				fprintf(fpout, " %-4s\n", "1");
			continue;
		}
		if (bond[i].type == 8) {
			suc = 1;
			if(strcmp(atom[bondi].ambername, "C.ar") != 0 && 
                           strcmp(atom[bondi].ambername, "N.ar") != 0 &&
			   strcmp(atom[bondi].ambername, "ca")   != 0 && 
			   strcmp(atom[bondi].ambername, "cp")   != 0 && 
			   strcmp(atom[bondi].ambername, "pb")   != 0 && 
                           strcmp(atom[bondi].ambername, "nb")   != 0) 
				suc = 0;
			if(strcmp(atom[bondj].ambername, "C.ar") != 0 && 
                           strcmp(atom[bondj].ambername, "N.ar") != 0 &&
			   strcmp(atom[bondj].ambername, "ca")   != 0 && 
			   strcmp(atom[bondj].ambername, "cp")   != 0 && 
			   strcmp(atom[bondj].ambername, "pb")   != 0 && 
                           strcmp(atom[bondj].ambername, "nb")   != 0) 
				suc = 0;
			if(suc == 1)
				fprintf(fpout, " %-4s\n", "ar");
			else
				fprintf(fpout, " %-4s\n", "2");
			continue;
		}
		if (bond[i].type == 10) {
			fprintf(fpout, " %-4s\n", "ar");
			continue;
		}
		if (bond[i].type == 9) {
			fprintf(fpout, " %-4s\n", "1");
			continue;
		}
		if (bond[i].type == 6) {
			fprintf(fpout, " %-4s\n", "2");
			continue;
		}
		fprintf(fpout, " %-4s\n", "1");
	}
	fprintf(fpout, "@<TRIPOS>SUBSTRUCTURE\n");
	fprintf(fpout, "%s%-12s%s", "     1 ", atom[0].aa,
			"1 TEMP              0 ****  ****    0 ROOT\n");
	fclose(fpout);
}

