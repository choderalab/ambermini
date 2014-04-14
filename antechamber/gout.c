/* Gaussian OUT */
int rgout(char *filename, int *atomnum, ATOM * atom, CONTROLINFO cinfo,
		  MOLINFO *minfo)
{
/*
The order of coordinates to be used 
(1) coordinates of standard orientation immediately after stationary is found
(2) coordinates of the last standard orientation
(3) coordinates of the last input orientation
*/
	int i;
	int numatom = 0;
	int Found_Stationary = 0;
	int Standard = 0;
	int read_flag = 0;
	int suc_flag = 0;	
	int scan_flag = 1;
	int charge_flag = 1;
	int overflow_flag = 0;
	int format = 0;
	long input_pos = -9999;
	long stand_pos = -9999;
	char tmpchar1[MAXCHAR];
	char tmpchar2[MAXCHAR];
	char tmpchar3[MAXCHAR];
	char tmpchar4[MAXCHAR];
	char tmpchar5[MAXCHAR];
	char line[MAXCHAR];
	FILE *fpin;

	if ((fpin = fopen(filename, "r")) == NULL) {
		fprintf(stdout, "Cannot open the input file %s in rgout(), exit\n", filename);
		exit(1);
	}
	initial(cinfo.maxatom, atom, (*minfo).resname);
	i = 0;
	for (;;) {
		sscanf(&line[0], "%s%s%s%s%s", tmpchar1, tmpchar2, tmpchar3,
			   tmpchar4, tmpchar5);
		if (scan_flag == 1 && charge_flag == 1 && strcmp(tmpchar1, "Charge") == 0) {
			sscanf(&line[9], "%lf", &(*minfo).dcharge);
			sscanf(&line[27], "%d", &(*minfo).multiplicity);
			(*minfo).icharge = (int) (*minfo).dcharge;
			charge_flag = 0;
			continue;
		}
		if (fgets(line, MAXCHAR, fpin) == NULL) {
			/* now go to standard_pos and read coordinates */
			scan_flag = 0;
			if(read_flag == 0 && suc_flag == 0) {
				read_flag = 1;
				if(stand_pos >= 0)
					fseek(fpin,stand_pos,0);
				continue;
			}
			/* now go to input_pos and read coordinates */
			if(read_flag == 1 && suc_flag == 0) {
				read_flag = 2;
				if(input_pos >= 0)
					fseek(fpin,input_pos,0);
				continue;
			}
			if(read_flag == 2 && suc_flag == 0) {
				fprintf(stdout, "No atom read in, the gaussian output file may not complete, exit\n");
				exit(1);
			}
			break;
		}
		if (scan_flag == 1 && strcmp("Input", tmpchar1) == 0
			&& strcmp("orientation:", tmpchar2) == 0) {
			input_pos = ftell(fpin);
			continue;
		}
		if (scan_flag == 1 && strcmp("Standard", tmpchar1) == 0
			&& strcmp("orientation:", tmpchar2) == 0) {
			if(Standard == 0) {
				stand_pos = ftell(fpin);
				if(Found_Stationary == 1)
					Standard = 1;
			}
/* we prefer to use the first standard orientation after Found_Stationary == 1 */
			continue;
		}
		if (scan_flag == 1 && Found_Stationary == 0 
			&& strcmp("--", tmpchar1) == 0
			&& strcmp("Stationary", tmpchar2) == 0
			&& strcmp("point", tmpchar3) == 0
			&& strcmp("found.", tmpchar4) == 0) {
			Found_Stationary = 1;
			continue;
		}
		if(scan_flag == 0 && (read_flag == 1 || read_flag == 2)) {
			if (line[28] == '.' && line[40] == '.' && line[52] == '.') {
				suc_flag = 1;
				format = 1;
				if (overflow_flag == 0) {
					sscanf(&line[10], "%d%lf%lf%lf", &atom[numatom].atomicnum,
					   	&atom[numatom].x, &atom[numatom].y, &atom[numatom].z);
					atom[numatom].charge = 0.0;
				}
				numatom++;
				if (numatom >= cinfo.maxatom && overflow_flag == 0) {
					printf
						("\nInfo: the atom number exceeds the MAXATOM, reallocate memory automatically");
					overflow_flag = 1;
				}
			}
			if (line[39] == '.' && line[51] == '.' && line[63] == '.') {
				suc_flag = 1;
				format = 2;
                        	if (overflow_flag == 0) {
                                	sscanf(&line[10], "%d", &atom[numatom].atomicnum);
                                	sscanf(&line[31], "%lf%lf%lf", &atom[numatom].x, &atom[numatom].y,
                                           	&atom[numatom].z);
                                	atom[numatom].charge = 0.0;
                        	}
                        	numatom++;
                        	if (numatom >= cinfo.maxatom && overflow_flag == 0) {
                                	printf
                                        	("\nInfo: the atom number exceeds the MAXATOM, reallocate memory automatically");
                                	overflow_flag = 1;
                        	}
			}
			if (suc_flag == 1 && format == 1 && line[28] != '.' && line[40] != '.' && line[52] != '.') 
				break;
			if (suc_flag == 1 && format == 2 && line[39] != '.' && line[51] != '.' && line[63] != '.') 
				break;
		}
	}
	*atomnum = numatom;
	fclose(fpin);
	element(*atomnum, atom);
	for (i = 0; i < *atomnum; i++)
		strcpy(atom[i].name, atom[i].element);
	return overflow_flag;
}


void wgout()
{
	printf
		("\n Sorry, you may get the gaussian output by running a Gaussian program (g98, g03, etc.)");
}
