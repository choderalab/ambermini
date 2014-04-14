int rac(char *filename, int *atomnum, ATOM * atom, int *bondnum,
		BOND * bond, CONTROLINFO *cinfo, MOLINFO *minfo)
{
	int index;
	int tmpint, tmpint1, tmpint2, tmpint3;
	double tmpf1, tmpf2, tmpf3, tmpf4;
	int numatom;
	int numbond;
	int terindex;
	int overflow_flag = 0;
	FILE *fpin;
	char line[MAXCHAR];
	char tmpchar[MAXCHAR];
	if ((fpin = fopen(filename, "r")) == NULL) {
		fprintf(stdout, "Cannot open the ac input file: %s, exit\n", filename);
		exit(1);
	}
	initial((*cinfo).maxatom, atom, (*minfo).resname);
	numatom = 0;
	numbond = 0;
	(*cinfo).bpindex = 0;
	terindex = -1;
	for (;;) {
		if (fgets(line, MAXCHAR, fpin) == NULL) {
/*       printf("\nFinished reading %s file.", filename); */
			break;
		}
		if (strncmp("TER", line, 3) == 0) {
			terindex = 1;
			continue;
		}
		if (strncmp("ATOM", line, 4) == 0) {
			if (overflow_flag == 0) {
                                line[26] = ' ';
				if (line[12] != ' ')
					index = 12;
				else
					index = 13;
				atom[numatom].name[0] = line[index];
				atom[numatom].name[1] = line[index + 1];
				atom[numatom].name[2] = line[index + 2];
				atom[numatom].name[3] = line[index + 3];

                                if(atom[numatom].name[1] == ' ')
                                        atom[numatom].name[1] = '\0';
                                else if(atom[numatom].name[2] == ' ')
                                        atom[numatom].name[2] = '\0';
                                else if(atom[numatom].name[3] == ' ')
                                        atom[numatom].name[3] = '\0';
                                else
                                        atom[numatom].name[4] = '\0';

				index = 17;
				if ((*cinfo).rnindex == 0) {
					atom[numatom].aa[0] = line[index];
					atom[numatom].aa[1] = line[index + 1];
					atom[numatom].aa[2] = line[index + 2];
					atom[numatom].aa[3] = line[index + 3];
				}

                                if(atom[numatom].aa[1] == ' ')
                                        atom[numatom].aa[1] = '\0';
                                else if(atom[numatom].aa[2] == ' ')
                                        atom[numatom].aa[2] = '\0';
                                else if(atom[numatom].aa[3] == ' ')
                                        atom[numatom].aa[3] = '\0';
                                else
                                        atom[numatom].aa[4] = '\0';

				if (line[21] != ' ')
					atom[numatom].chain[0] = line[21];
				if (terindex == 1) {
					atom[numatom].ter = terindex;
/*the first atom followed after TER line has terindex of 1*/
					terindex = -1;
				}
				sscanf(&line[22], "%d%lf%lf%lf%lf%s",
					   &tmpint, &tmpf1, &tmpf2, &tmpf3, &tmpf4, tmpchar);

				atom[numatom].resno = tmpint;
				atom[numatom].x = tmpf1;
				atom[numatom].y = tmpf2;
				atom[numatom].z = tmpf3;
				atom[numatom].charge = tmpf4;
				strcpy(atom[numatom].ambername, tmpchar);
			}
			numatom++;

			if (numatom >= (*cinfo).maxatom && overflow_flag == 0) {
				printf
					("\nInfo: the atom number exceeds the MAXATOM, reallocate memory automatically");
				overflow_flag = 1;
			}
		}
		if (strncmp("BOND", line, 4) == 0) {
			if (overflow_flag == 0) {
				sscanf(&line[4], "%d%d%d%s", &tmpint, &tmpint1,
					   &tmpint2, tmpchar);
				tmpint1--;
				tmpint2--;
				atom[tmpint1].con[atom[tmpint1].connum++] = tmpint2;
				atom[tmpint2].con[atom[tmpint2].connum++] = tmpint1;
				bond[numbond].bondi = tmpint1;
				bond[numbond].bondj = tmpint2;
				bond[numbond].jflag = -1;
				tmpint3 = strlen(tmpchar);
				if(tmpint3 >=2 && (tmpchar[tmpint3-1] == 'f' || tmpchar[tmpint3-1] == 'F')) {
/* bond type is frozen during bond type assignment */
					tmpchar[tmpint3-1] = '\0'; 
					bond[numbond].jflag = 1;
				}
				bond[numbond].type = atoi(tmpchar);
				if(bond[numbond].jflag == 1 && bond[numbond].type == 10) {
					fprintf(stdout, "You cannot freeze this bond (%s-%s,%d) since it is unclear whether it is a single or double bond\n",atom[bond[numbond].bondi].name, atom[bond[numbond].bondj].name, bond[numbond].type); 
					fprintf(stdout, "Try to use '7' for aromatic single and '8' for aromatic double bond, exit\n"); 
					exit(1);
				}
				if (tmpint3 == 0)
					(*cinfo).bpindex = 1;
			}
			numbond++;
			if (numbond >= (*cinfo).maxbond && overflow_flag == 0) {
				printf
					("\nInfo: the bond number exceeds the MAXBOND, reallocate memory, automatically");
				overflow_flag = 1;
			}
		}
		if (strncmp("CHARGE", line, 6) == 0) {
			sscanf(&line[6], "%lf", &(*minfo).dcharge);
			sscanf(&line[18], "%d", &(*minfo).icharge);
		}

	}
	*atomnum = numatom;
	*bondnum = numbond;
/* printf("\n The atomic number is %5d\n", atomnum); */
	fclose(fpin);
	return overflow_flag;

}


/* Antechamber */
void wac(char *filename, int atomnum, ATOM * atom, int bondnum,
		 BOND * bond, CONTROLINFO cinfo, MOLINFO minfo)
{
	int i;
	char form[5 * MAXCHAR];
	double fraction;
	double tmpf;
	char resname[10], atomname[10];

	FILE *fpout;
	if ((fpout = fopen(filename, "w")) == NULL) {
		fprintf(stdout, "Cannot open a file (%s) to write in ac format, exit\n", filename);
		exit(1);
	}
        if (minfo.usercharge > -9999) { /*charge read in with -nc flag, it is unlikely users input a charge smaller than -9999 */
                minfo.icharge = minfo.usercharge;
                minfo.dcharge = minfo.usercharge;
        }

        if(minfo.dcharge < -9990) {
        	minfo.dcharge = 0.0;
                for (i = 0; i < atomnum; i++)
                	minfo.dcharge += atom[i].charge;
                fraction = modf(minfo.dcharge, &tmpf);
                minfo.icharge = (int) tmpf;
                if (fabs(fraction) >= 0.50) {
                	if (minfo.dcharge < 0)
                        	minfo.icharge--;
                        if (minfo.dcharge > 0)
                        	minfo.icharge++;
                }
	}

/*zero weird charges */
        if(minfo.usercharge < -9990 && (minfo.icharge <= -100 || minfo.icharge >= 100)) {
                fprintf(stdout, "Warning: Weird total charge: %d!"
            "The net charge is assumed to be 0.\n"
                      "         If the weird charge was correct, "
            "specify it via the -nc net charge flag.\n", minfo.icharge);
                minfo.icharge = 0;
        }
/*
	if (minfo.dcharge < -9990) {
		minfo.dcharge = 0.0;
		for (i = 0; i < atomnum; i++)
			minfo.dcharge += atom[i].charge;
	}
	if (minfo.usercharge < -9990.)
		minfo.icharge = intcharge(atomnum, atom);
        if (minfo.icharge < -9990.)
                minfo.icharge = intcharge(atomnum, atom);
	if (minfo.icharge < -100.0 || minfo.icharge > 100.0) 
		fprintf(stdout, "Warning: weird total charge (%d) !\n", minfo.icharge); 
*/
	fprintf(fpout, "CHARGE %9.2lf ( %d )\n", minfo.dcharge, minfo.icharge);
	formula(atomnum, atom, form);
	fprintf(fpout, "Formula: %s\n", form);
	for (i = 0; i < atomnum; i++) {
		if(strlen(atom[i].name) > 4) {
			atomname[0] = atom[i].name[0];
			atomname[1] = atom[i].name[1];
			atomname[2] = atom[i].name[2];
			atomname[3] = atom[i].name[3];
			atomname[4] = '\0';
		} 
		else                  
			strcpy(atomname, atom[i].name);
		if(strlen(atom[i].aa) > 3) {
			resname[0] = atom[i].aa[0];
			resname[1] = atom[i].aa[1];
			resname[2] = atom[i].aa[2];
			resname[3] = '\0';
		} 
		else                  
			strcpy(resname, atom[i].aa);
		fprintf(fpout,
				"ATOM%7d  %-4s%-4s%5d%12.3f%8.3f%8.3f%10.6lf%10s\n",
				i + 1, atomname, resname, atom[i].resno, atom[i].x,
				atom[i].y, atom[i].z, atom[i].charge, atom[i].ambername);
	}

/*
	for (i = 0; i < atomnum; i++)
		for (j = i + 1; j < atomnum; j++)
			for (k = 0; k < bondnum; k++)
				if ((bond[k].bondi == i && bond[k].bondj == j) ||
					(bond[k].bondi == j && bond[k].bondj == i))
					fprintf(fpout, "BOND%5d%5d%5d%5d  %5s%5s\n", 1 + num++,
							i + 1, j + 1, bond[k].type, atom[i].name,
							atom[j].name);
*/
	for (i = 0; i < bondnum; i++) 
		if(bond[i].jflag == 1)
			fprintf(fpout, "BOND%5d%5d%5d%5d%-s  %5s%5s\n", i + 1, bond[i].bondi + 1, bond[i].bondj + 1, 
				bond[i].type, "F", atom[bond[i].bondi].name, atom[bond[i].bondj].name);
		else
			fprintf(fpout, "BOND%5d%5d%5d%5d  %5s%5s\n", i + 1, bond[i].bondi + 1, bond[i].bondj + 1, 
				bond[i].type, atom[bond[i].bondi].name, atom[bond[i].bondj].name);
	fclose(fpout);
}
