/* GCRT */
int rgcrt(char *filename, int *atomnum, ATOM * atom, CONTROLINFO cinfo,
		  MOLINFO minfo)
{
	FILE *fpin;
	int i, index;
	int nameindex;
	int numatom;
	int read_charge = 1;
	int overflow_flag = 0;
	char tmpchar[MAXCHAR];
	char line[MAXCHAR];


	if ((fpin = fopen(filename, "r")) == NULL) {
		fprintf(stdout, "Cannot open the gcrt file %s to read, exit\n", filename);
		exit(1);
	}
	initial(cinfo.maxatom, atom, minfo.resname);
	index = 0;
	numatom = 0;
	nameindex = -1;
	for (;;) {
		if (fgets(line, MAXCHAR, fpin) == NULL) {
			break;
		}
		if (spaceline(line) == 1 || strlen(line) <= 1) {
			index++;
			continue;
		}
		if (index == 0 || index == 1) continue;
		if (read_charge ==1 && index == 2) {
			sscanf(line, "%d%d", &(minfo.icharge), &(minfo.multiplicity));
			read_charge = 0;
			continue;
		}
		if (index == 3) {
			break;
		}
		sscanf(line, "%s", tmpchar);
		if (nameindex == -1) {
			if (tmpchar[0] >= '0' && tmpchar[0] <= '9')
				nameindex = 0;
			else
				nameindex = 1;
		}
		if (overflow_flag == 0) {
			if (nameindex == 0)
				sscanf(line, "%d%lf%lf%lf", &atom[numatom].atomicnum,
					   &atom[numatom].x, &atom[numatom].y,
					   &atom[numatom].z);
			else
				sscanf(line, "%s%lf%lf%lf", atom[numatom].name,
					   &atom[numatom].x, &atom[numatom].y,
					   &atom[numatom].z);
		}
		numatom++;
		if (numatom >= cinfo.maxatom && overflow_flag == 0) {
			printf
				("\nInfo: the atom number exceeds the MAXATOM, reallocate memory automatically");
			overflow_flag = 1;
		}
	}
	*atomnum = numatom;
	if (nameindex == 0) {
		element(*atomnum, atom);
		for (i = 0; i < *atomnum; i++)
			strcpy(atom[i].name, atom[i].element);
	}
	if (nameindex == 1)
		atomicnum(*atomnum, atom);
/* printf("\n atom number is  %5d", *atomnum); */
	fclose(fpin);
	return overflow_flag;
}

void wgcrt(char *filename, int atomnum, ATOM atom[], MOLINFO minfo)
{
	FILE *fpin;
	FILE *fpout;
	char *amberhome;
	char espparm_file[MAXCHAR];
	char line[MAXCHAR];
	char akeyword[MAXCHAR]="";
	char ckeyword[MAXCHAR];
	int i,j;
	int iradius_flag;
	int ibs= 0;
	int esp_flag; 
	int nespparm = 0;
	int nbasisset = 0;
	double default_radius = 1.7;
	ESPPARM espparm[120];
	BASISSET basisset[100];

	/* int index; */
	if ((fpout = fopen(filename, "w")) == NULL) {
		fprintf(stdout, "Cannot open a file %s to write in wgcrt(), exit\n", filename);
		exit(1);
	}
	fprintf(fpout, "%s\n", "--Link1--");
	if(strlen(minfo.gn) >= 5) 
		fprintf(fpout, "%s\n", minfo.gn);
	fprintf(fpout, "%s%s\n", "%chk=", minfo.chkfile);
	if(strlen(minfo.gm) >= 4) 
		fprintf(fpout, "%s\n", minfo.gm);

//	check ESP-related keyword
	esp_flag = 0;
        for(i=0;i<strlen(minfo.gkeyword);i++)
                ckeyword[i] = toupper(minfo.gkeyword[i]);
        if((strstr(ckeyword, "POP=") != 0 || strstr(ckeyword, "POP(")  != 0)  &&
           (strstr(ckeyword, "MK")   != 0 || strstr(ckeyword, "CHELP") != 0))
                esp_flag = 1;

//	when the default gaussian keyword is used, or esp_flag ==1, read ESP.PARM
	if(minfo.igkeyword == 0 || esp_flag == 1) {
	    	amberhome = (char *) getenv("AMBERHOME");
    	    	if( amberhome == NULL ){
       			fprintf( stdout, "AMBERHOME is not set!\n" );
       			exit(1);
    		}
		strcpy(espparm_file, amberhome);	
		strcat(espparm_file, "/dat/antechamber/ESPPARM.DAT");
        	if ((fpin = fopen(espparm_file, "r")) == NULL) {
                	fprintf(stdout, "Cannot open espparm_file %s in read_espparm(), exit\n", espparm_file);
                	exit(1);
        	}
        	for (;;) {
                	if (fgets(line, MAXCHAR, fpin) == NULL) break;
			if(strncmp(line, "DEFAULT RADIUS", 14) == 0) 
				sscanf(&line[14], "%lf", &default_radius);
			if(strncmp(line, "BASIS SET", 9) == 0) {
				sscanf(&line[9], "%d%s", &basisset[nbasisset].id, 
                                                           basisset[nbasisset].bs);
				nbasisset++;
			}
			if(strncmp(line, "PARM", 4) == 0) {
				sscanf(&line[4], "%d%s%lf%lf%d%d", &espparm[nespparm].atomicnum, espparm[nespparm].elem,
							            &espparm[nespparm].vdw,      &espparm[nespparm].mk,	
                                                                    &espparm[nespparm].flag,     &espparm[nespparm].bs);
				nespparm++;
			}
        	}
		fclose(fpin);

		iradius_flag = 0;
		ibs = 0;
		for(i=0;i<atomnum;i++) {
			for(j=0;j<nespparm;j++) 
				if(atom[i].atomicnum == espparm[j].atomicnum || strcmp(atom[i].element, espparm[j].elem) == 0) {
					if(ibs < espparm[j].bs) ibs=espparm[j].bs;
					if(espparm[j].flag != 0) {
						iradius_flag = 1;
						espparm[j].flag = 2;
					}	
					break;
				}
		}

		if(minfo.igkeyword == 0) {
			strcpy(minfo.gkeyword, "#HF/");
			strcat(minfo.gkeyword, basisset[ibs-1].bs);
			strcat(minfo.gkeyword, " SCF=tight Test Pop=MK iop(6/33=2) iop(6/42=6) opt");
		}
	}
//	additional keywords
	if(esp_flag == 1) {
		if(iradius_flag == 1) {
			if(strstr(minfo.gkeyword, "ReadRadii") == 0 &&
                           strstr(minfo.gkeyword, "READRADII") == 0 && 
                           strstr(minfo.gkeyword, "readradii") == 0) {
				strcat(akeyword, "Pop=ReadRadii");	
			}
		}
		if(minfo.gv == 1) {
			if(strstr(minfo.gkeyword, "6/50=1") == 0) {
				strcat(akeyword, " iop(6/50=1)");
			}
		}
	}	
	
	if(strlen(akeyword) >= 1) {
		fprintf(fpout, "%s\n", minfo.gkeyword);
		fprintf(fpout, "#%s\n", akeyword);	
		fprintf(fpout, "\n");
	} else {
		fprintf(fpout, "%s\n\n", minfo.gkeyword);
	}
	fprintf(fpout, "%s\n\n", "remark line goes here");
	fprintf(fpout, "%d%4d\n", minfo.icharge, minfo.multiplicity);
	element(atomnum, atom);

	for (i = 0; i < atomnum; i++)
//		fprintf(fpout, "%5s%12.4lf    %12.4lf    %12.4lf     \n",
//				atom[i].element, atom[i].x, atom[i].y, atom[i].z);
		fprintf(fpout, "%5s%16.10lf    %16.10lf    %16.10lf     \n",
				atom[i].element, atom[i].x, atom[i].y, atom[i].z);

	if(esp_flag == 1) {
		if(minfo.gv == 1 && iradius_flag == 1) {
			for(i=0;i<nespparm;i++)
				if(espparm[i].flag == 2) {
					if(espparm[i].mk <= 0) espparm[i].mk=default_radius;
					fprintf(fpout, "\n%s     %8.2lf", espparm[i].elem, espparm[i].mk);
				}
			fprintf(fpout, "\n\n%s\n", minfo.gesp);
			for(i=0;i<nespparm;i++)
				if(espparm[i].flag == 2) {
					if(espparm[i].mk <= 0) espparm[i].mk=default_radius;
					fprintf(fpout, "\n%s     %8.2lf", espparm[i].elem, espparm[i].mk);
				}
			fprintf(fpout, "\n\n%s\n", minfo.gesp);
		}
		if(minfo.gv == 0 && iradius_flag == 1) {
			for(i=0;i<nespparm;i++)
				if(espparm[i].flag == 2) {
					if(espparm[i].mk <= 0) espparm[i].mk=default_radius;
					fprintf(fpout, "\n%s     %8.2lf", espparm[i].elem, espparm[i].mk);
				}
			fprintf(fpout, "\n");
			for(i=0;i<nespparm;i++)
				if(espparm[i].flag == 2) {
					if(espparm[i].mk <= 0) espparm[i].mk=default_radius;
					fprintf(fpout, "\n%s     %8.2lf", espparm[i].elem, espparm[i].mk);
				}
			fprintf(fpout, "\n");
		}
		if(minfo.gv == 1 && iradius_flag == 0) {
			fprintf(fpout, "\n%s\n", minfo.gesp);
			fprintf(fpout, "\n%s\n", minfo.gesp);
		}
	}
	fprintf(fpout, "\n\n");
	fclose(fpout);
}
