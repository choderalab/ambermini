/* CHARMM */
int rcharmm(char *filename, int *atomnum, ATOM atom[], int *bondnum,
		  BOND bond[], CONTROLINFO *cinfo, MOLINFO *minfo)
{
int i;
int numatom;
int numatom2;
int numbond;
char tmpchar[MAXCHAR], tmpchar1[MAXCHAR], tmpchar2[MAXCHAR], tmpchar3[MAXCHAR], tmpchar4[MAXCHAR], tmpchar5[MAXCHAR];
char tmpchar6[MAXCHAR], tmpchar7[MAXCHAR], tmpchar8[MAXCHAR], tmpchar9[MAXCHAR], tmpchar10[MAXCHAR];
char line[MAXCHAR];
char rtf[MAXCHAR];
int irtf = 0;
int rflag = 0;
int flag1, flag2;
int bondi, bondj;
int overflow_flag = 0;
FILE *fpin;

if ((fpin = fopen(filename, "r")) == NULL) {
	fprintf(stdout, "Cannot open the charmm input file: %s, exit\n", filename);
	exit(1);
}
initial((*cinfo).maxatom, atom, (*minfo).resname);
strcpy(rtf, filename);
numatom = 0;
numbond = 0;
for (;;) {
	if (fgets(line, MAXCHAR, fpin) == NULL) break;
	if(line[0] == '*') continue;
	strcpy(tmpchar1, "");
	strcpy(tmpchar2, "");
	strcpy(tmpchar3, "");
	strcpy(tmpchar4, "");
	strcpy(tmpchar5, "");
	tmpchar1[0]='\0';
	tmpchar2[0]='\0';
	tmpchar3[0]='\0';
	tmpchar4[0]='\0';
	tmpchar5[0]='\0';
	sscanf(line, "%s%s%s%s%s", tmpchar1, tmpchar2, tmpchar3, tmpchar4, tmpchar5);
	 if((irtf == 0) && 
	   (strcmp(tmpchar1, "Read") == 0 || strcmp(tmpchar1, "READ") == 0 || strcmp(tmpchar1, "read") == 0) && 
	   (strcmp(tmpchar2, "Rtf") == 0  || strcmp(tmpchar2, "RTF") == 0  || strcmp(tmpchar2, "rtf") == 0) &&
	   (strcmp(tmpchar3, "Card") == 0 || strcmp(tmpchar3, "CARD") == 0 || strcmp(tmpchar3, "card") == 0) && 
           (strcmp(tmpchar4, "Name") == 0 || strcmp(tmpchar4, "NAME") == 0 || strcmp(tmpchar4, "name") == 0)) {
		strcpy(rtf, tmpchar5);
		irtf = 1;
		continue;
	}
	 if((rflag == 0) && 
	   (strcmp(tmpchar1, "Read") == 0 || strcmp(tmpchar1, "READ") == 0 || strcmp(tmpchar1, "read") == 0) && 
           (strcmp(tmpchar2, "Coor") == 0 || strcmp(tmpchar2, "COOR") == 0 || strcmp(tmpchar2, "coor") == 0) &&
	   (strcmp(tmpchar3, "Card") == 0 || strcmp(tmpchar3, "CARD") == 0 || strcmp(tmpchar3, "card") == 0)) { 
			rflag = 1;
			continue;
	}
	if(rflag == 1) {
		sscanf(line, "%d", &numatom2);
		if (numatom2 >= (*cinfo).maxatom && overflow_flag == 0) {
			printf ("\nThe atom number exceeds the MAXATOM, reallocate memory");
			overflow_flag = 1;
			return overflow_flag;
		}
		rflag = 2;
		numatom = 0;
		continue;	
	}
	if(rflag == 2) {
		sscanf(line, "%d%s%s%s%lf%lf%lf", &atom[numatom].id, atom[numatom].resid, atom[numatom].aa, atom[numatom].name,
						  &atom[numatom].x, &atom[numatom].y, &atom[numatom].z);
		numatom++;
		if(numatom == numatom2) break;
	}

}
fclose(fpin);

if ((fpin = fopen(rtf, "r")) == NULL) {
        fprintf(stdout, "Cannot open the rtf file in rcharmm() : %s, exit\n", rtf);
        exit(1);
}
numatom = 0;
numbond = 0;
for (;;) {
        if (fgets(line, MAXCHAR, fpin) == NULL) break;
        if(line[0] == '*') continue;
	strcpy(tmpchar, "");
	strcpy(tmpchar1, "");
	strcpy(tmpchar2, "");
	strcpy(tmpchar3, "");
	strcpy(tmpchar4, "");
	strcpy(tmpchar5, "");
	strcpy(tmpchar6, "");
	strcpy(tmpchar7, "");
	strcpy(tmpchar8, "");
	strcpy(tmpchar9, "");
	strcpy(tmpchar10, "");
	tmpchar[0]='\0';
	tmpchar1[0]='\0';
	tmpchar2[0]='\0';
	tmpchar3[0]='\0';
	tmpchar4[0]='\0';
	tmpchar5[0]='\0';
	tmpchar6[0]='\0';
	tmpchar7[0]='\0';
	tmpchar8[0]='\0';
	tmpchar9[0]='\0';
	tmpchar10[0]='\0';
        sscanf(line, "%s%s%s%s%s%s%s%s%s%s%s", tmpchar, tmpchar1, tmpchar2, tmpchar3, tmpchar4, tmpchar5, tmpchar6,
						 tmpchar7, tmpchar8, tmpchar9, tmpchar10);
	if(strcmp(tmpchar, "RESI") == 0)
		(*minfo).dcharge = atof(tmpchar2);
	if(strcmp(tmpchar, "ATOM") == 0) {
		strcpy(atom[numatom].ambername, tmpchar2);
		atom[numatom].charge = atof(tmpchar3);
		numatom++;
	}
	if(strcmp(tmpchar, "BOND") == 0) {
		if(strlen(tmpchar1) > 0 && strlen(tmpchar2) > 0) {
			flag1 = 0;
			flag2 = 0;
			for(i=0;i<numatom2;i++) {
				if(flag1 == 1 && flag2 == 1) break;
				if(strcmp(tmpchar1, atom[i].name) == 0) {
					bondi = i;
					flag1 = 1;
				}
				if(strcmp(tmpchar2, atom[i].name) == 0) {
					bondj = i;
					flag2 = 1;
				}
			}
		
			if(flag1 ==0 || flag2 == 0) continue;
			if(overflow_flag == 0) {
				bond[numbond].bondi = bondi;	
				bond[numbond].bondj = bondj;	
			}
			numbond++;
                	if (numbond >= (*cinfo).maxbond && overflow_flag == 0) {
                		printf ("\nThe bond number exceeds the MAXBOND, reallocate memory");
                        	overflow_flag = 1;
                	}
		}
                if(strlen(tmpchar3) > 0 && strlen(tmpchar4) > 0) {
                        flag1 = 0;
                        flag2 = 0;
                        for(i=0;i<numatom2;i++) {
                                if(flag1 == 1 && flag2 == 1) break;
                                if(strcmp(tmpchar3, atom[i].name) == 0) {
                                        bondi = i;
                                        flag1 = 1;
                                }
                                if(strcmp(tmpchar4, atom[i].name) == 0) {
                                        bondj = i;
                                        flag2 = 1;
                                }
                        }
			if(flag1 ==0 || flag2 == 0) continue;
                        if(overflow_flag == 0) {
                                bond[numbond].bondi = bondi;
                                bond[numbond].bondj = bondj;
                        }
                        numbond++;
                        if (numbond >= (*cinfo).maxbond && overflow_flag == 0) {
                                printf ("\nThe bond number exceeds the MAXBOND, reallocate memory");
                                overflow_flag = 1;
                        }
                }
                if(strlen(tmpchar5) > 0 && strlen(tmpchar6) > 0) {
                        flag1 = 0;
                        flag2 = 0;
                        for(i=0;i<numatom2;i++) {
                                if(flag1 == 1 && flag2 == 1) break;
                                if(strcmp(tmpchar5, atom[i].name) == 0) {
                                        bondi = i;
                                        flag1 = 1;
                                }
                                if(strcmp(tmpchar6, atom[i].name) == 0) {
                                        bondj = i;
                                        flag2 = 1;
                                }
                        }

			if(flag1 ==0 || flag2 == 0) continue;
                        if(overflow_flag == 0) {
                                bond[numbond].bondi = bondi;
                                bond[numbond].bondj = bondj;
                        }
                        numbond++;
                        if (numbond >= (*cinfo).maxbond && overflow_flag == 0) {
                                printf ("\nThe bond number exceeds the MAXBOND, reallocate memory");
                                overflow_flag = 1;
                        }
                }

                if(strlen(tmpchar7) > 0 && strlen(tmpchar8) > 0) {
                        flag1 = 0;
                        flag2 = 0;
                        for(i=0;i<numatom2;i++) {
                                if(flag1 == 1 && flag2 == 1) break;
                                if(strcmp(tmpchar7, atom[i].name) == 0) {
                                        bondi = i;
                                        flag1 = 1;
                                }
                                if(strcmp(tmpchar8, atom[i].name) == 0) {
                                        bondj = i;
                                        flag2 = 1;
                                }
                        }

			if(flag1 ==0 || flag2 == 0) continue;
                        if(overflow_flag == 0) {
                                bond[numbond].bondi = bondi;
                                bond[numbond].bondj = bondj;
                        }
                        numbond++;
                        if (numbond >= (*cinfo).maxbond && overflow_flag == 0) {
                                printf ("\nThe bond number exceeds the MAXBOND, reallocate memory");
                                overflow_flag = 1;
                        }
                }

                if(strlen(tmpchar9) > 0 && strlen(tmpchar10) > 0) {
                        flag1 = 0;
                        flag2 = 0;
                        for(i=0;i<numatom2;i++) {
                                if(flag1 == 1 && flag2 == 1) break;
                                if(strcmp(tmpchar9, atom[i].name) == 0) {
                                        bondi = i;
                                        flag1 = 1;
                                }
                                if(strcmp(tmpchar10, atom[i].name) == 0) {
                                        bondj = i;
                                        flag2 = 1;
                                }
                        }

			if(flag1 ==0 || flag2 == 0) continue;
                        if(overflow_flag == 0) {
                                bond[numbond].bondi = bondi;
                                bond[numbond].bondj = bondj;
                        }
                        numbond++;
                        if (numbond >= (*cinfo).maxbond && overflow_flag == 0) {
                                printf ("\nThe bond number exceeds the MAXBOND, reallocate memory");
                                overflow_flag = 1;
                        }
                }
	}

}
*bondnum = numbond;
*atomnum = numatom2;
return overflow_flag;
}

/* CHARMM */
void wcharmm(char *filename, char *ifilename, int atomnum, ATOM * atom,
                  int bondnum, BOND * bond, CONTROLINFO cinfo, MOLINFO minfo)
{
	char tmpchar[MAXCHAR];
	size_t copied_size;
	wac("ANTECHAMBER_PREP.AC0", atomnum, atom, bondnum, bond, cinfo,
		minfo);
	system("cp -rf ANTECHAMBER_PREP.AC0 ANTECHAMBER_PREP.AC");
/*part1: if intype is not prepi, prepc and ac, judge atom type*/
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
                strncat(tmpchar, minfo.atom_type_def,
			sizeof(tmpchar) - strlen(tmpchar) -1 );
                if(cinfo.intstatus == 2)
			fprintf(stdout, "\nRunning: %s\n", tmpchar);
		system(tmpchar);
	}

        copied_size = build_exe_path(tmpchar, "charmmgen", sizeof tmpchar, 1);
	strncat(tmpchar, " -i ANTECHAMBER_PREP.AC -f ac -o ",
		sizeof tmpchar - copied_size);
        strncat(tmpchar, filename, sizeof(tmpchar) - strlen(tmpchar) -1);
        strcat(tmpchar, " -r ");
        strcat(tmpchar, minfo.resname);
	if(cinfo.intstatus == 2)
       	     	fprintf(stdout, "\nRunning: %s\n", tmpchar);
        system(tmpchar);
/*  system("rm -f ANTECHAMBER_PREP.AC");*/
}

