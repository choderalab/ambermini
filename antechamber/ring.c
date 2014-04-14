int ring_detect_cycle(ATOM atm[], int selectnum, int startnum,
					  int *ringnum, RING ring[], int *selectarom,
					  int *selectaromindex, int run_flag)
{
/* determain the aromatic atom through an iterative way */
	int i, j, k, m;
	int start;
	int breakindex = 0;
	start = -1;
	selectarom[selectnum++] = startnum;
	selectaromindex[startnum] = 1;
	for (i = 0; i < 4; i++) {
		if (atm[startnum].con[i] == -1)
			return 0;
		start = atm[startnum].con[i];
		if (atm[start].atomicnum != 6 && atm[start].atomicnum != 7 &&
			atm[start].atomicnum != 8 && atm[start].atomicnum != 16 &&
			atm[start].atomicnum != 15)
			continue;
		if (atm[start].atomicnum == 6 && atm[start].connum <= 2)
			continue;
		if (atm[start].atomicnum == 8 && atm[start].connum == 1)
			continue;
		if (atm[start].atomicnum == 16 && atm[start].connum == 1)
			continue;
		for (j = 0; j < selectnum; j++)
			if (selectarom[j] == start) {
				breakindex = 1;
				break;
			}
		if (breakindex == 1) {
			breakindex = 0;
			continue;
		}
		if (selectnum > 10)
			return 0;
		/* we have already visited this atom */
		for (j = 2; j < 10; j++)	/* deal with 3-10 member ring */
			if (selectnum == j)
				for (k = 0; k < 4; k++) {
					if (atm[selectarom[0]].con[k] == -1)
						break;
					if (atm[selectarom[0]].con[k] == start) {
						if(run_flag == 1) {
							for (m = 0; m < j; m++)
								ring[*ringnum].atomno[m] = selectarom[m];
							ring[*ringnum].atomno[j] = start;
							ring[*ringnum].num = j + 1;
						}
						(*ringnum)++;
						break;
					}
				}
		if (start == -1)
			return 0;
		ring_detect_cycle(atm, selectnum, start, ringnum, ring, selectarom, selectaromindex, run_flag);
	}
	return 0;
}

void purify(ATOM atm[], int ringnum, RING ring[])
{
	int i, j, k;
	RING *ringbak; 
	int ringbaknum = 0;
	int tmpint;
	int index;

	if(ringnum != 0) {
        	ringbak = (RING *) malloc(sizeof(RING) * ringnum );
        	if (ringbak == NULL) {
                	fprintf(stdout, "memory allocation error for *ringbak\n");
                	exit(1);
        	}
	}
	else 
		return;

	for (i = 0; i < ringnum; i++)
		for (j = 0; j < ring[i].num; j++)
			for (k = j + 1; k < ring[i].num; k++) {
				if (ring[i].atomno[j] > ring[i].atomno[k]) {
					tmpint = ring[i].atomno[k];
					ring[i].atomno[k] = ring[i].atomno[j];
					ring[i].atomno[j] = tmpint;
				}
			}

	for (i = 0; i < ringnum; i++)
		for (j = i + 1; j < ringnum; j++)
			if (ring[i].num == ring[j].num && ring[i].num != 0) {
				index = 1;
				for (k = 0; k < ring[i].num; k++)
					if (ring[i].atomno[k] != ring[j].atomno[k]) {
						index = 0;
						break;
					}
				if (index == 1)
					ring[j].num = 0;
			}

	for (i = 0; i < ringnum; i++)
		for (j = 0; j < ring[i].num; j++) {
			index = 0;
			for (k = 0; k < ring[i].num; k++) {
				if (atm[ring[i].atomno[j]].con[0] == ring[i].atomno[k])
					index++;
				if (atm[ring[i].atomno[j]].con[1] == ring[i].atomno[k])
					index++;
				if (atm[ring[i].atomno[j]].con[2] == ring[i].atomno[k])
					index++;
				if (atm[ring[i].atomno[j]].con[3] == ring[i].atomno[k])
					index++;
				if (atm[ring[i].atomno[j]].con[4] == ring[i].atomno[k])
					index++;
				if (atm[ring[i].atomno[j]].con[5] == ring[i].atomno[k])
					index++;
			}
			if (index == 3) {
				ring[i].num = 0;
				break;
			}
		}

	for (i = 0; i < ringnum; i++)
		if (ring[i].num != 0) {
			for (j = 0; j < ring[i].num; j++)
				ringbak[ringbaknum].atomno[j] = ring[i].atomno[j];
			ringbak[ringbaknum].num = ring[i].num;
			ringbaknum++;

		}
	ringnum = ringbaknum;
	for (i = 0; i < ringnum; i++) {
		ring[i].num = ringbak[i].num;
		for (j = 0; j < ring[i].num; j++)
			ring[i].atomno[j] = ringbak[i].atomno[j];
	}
	free(ringbak);
}
void ringproperty(int ringnum, RING ring[], AROM arom[])
{
	int i, j;
	int tmpint;

	for (i = 0; i < ringnum; i++) {
		for (j = 0; j < ring[i].num; j++) {
			tmpint = ring[i].atomno[j];
			arom[tmpint].rg[0]++;
			arom[tmpint].rg[ring[i].num]++;
		}
	}
}



void aromatic(int atnum, ATOM atm[], int bondnum, BOND bond[], int ringnum,
			  int maxatom, RING ring[], AROM arom[])
{
	int i, j, k;
	int *initarom;
	int tmpint;
	int tmpint1;
/*
  int tmpint2;
  int tmpint3;
  int tmpint4;
  int tmpint5;
*/
	int index;
	int index0;
	int continue_index = 0;

	initarom = (int *) malloc(sizeof(int) * maxatom);
	if (initarom == NULL) {
		fprintf(stdout, "memory allocation error for *initarom\n");
		exit(1);
	}
	for (i = 0; i < atnum; i++) {
		initarom[i] = 0;
		arom[i].ar1 = 0;
		arom[i].ar2 = 0;
		arom[i].ar3 = 0;
		arom[i].ar4 = 0;
		arom[i].ar5 = 0;
		atm[i].saturate = -1;
		switch (atm[i].atomicnum) {
		case 6:
			atm[i].ewd = 0;
			if (atm[i].connum == 3)
				initarom[i] = 2;
			if (atm[i].connum == 4) {
				initarom[i] = -2;
				atm[i].saturate = 1;
			}
			break;
		case 7:
			atm[i].ewd = 1;
			if (atm[i].connum <= 3)
				initarom[i] = 2;
			if (atm[i].connum >= 3)
				atm[i].saturate = 1;
			break;
		case 8:
			atm[i].ewd = 1;
			if (atm[i].connum == 2) {
				initarom[i] = 1;
				atm[i].saturate = 1;
			}
			break;
		case 15:
			atm[i].ewd = 0;
			if (atm[i].connum == 2)
				initarom[i] = 2;
//WJM 2010
//			if (atm[i].connum > 2)
			if (atm[i].connum == 3)
				initarom[i] = 1;
			if (atm[i].connum >= 4)
				initarom[i] = -1;
			if (atm[i].connum >= 3)
				atm[i].saturate = 1;
			break;
		case 16:
			atm[i].ewd = 1;
//WJM 2010
//			if (atm[i].connum >= 2)
			if (atm[i].connum == 2)
				initarom[i] = 1;
			if (atm[i].connum >= 3)
				initarom[i] = -1;
			if (atm[i].connum >= 3)
				atm[i].saturate = 1;
			break;
		case 9:
			atm[i].ewd = 1;
			atm[i].saturate = 1;
			break;
		case 17:
			atm[i].ewd = 1;
			atm[i].saturate = 1;
			break;
		case 35:
			atm[i].ewd = 1;
			atm[i].saturate = 1;
			break;
		case 53:
			atm[i].ewd = 1;
			atm[i].saturate = 1;
			break;
		case  1:
			atm[i].ewd = 0;
			atm[i].saturate = 1;
			break;
		default:
			atm[i].ewd = 0;
		}
	}

	for (i = 0; i < ringnum; i++) {
		tmpint = 0;
		for (j = 0; j < ring[i].num; j++)
			tmpint += initarom[ring[i].atomno[j]];
/* for pure aliphatic rings*/
		if (tmpint == -2 * ring[i].num) {
			for (j = 0; j < ring[i].num; j++)
				arom[ring[i].atomno[j]].ar5++;
			continue;
		}

/* for rings contains sp3 carbon */
		for (k = 0; k < ring[i].num; k++)
			if (initarom[ring[i].atomno[k]] < 0) {
				for (j = 0; j < ring[i].num; j++)
					arom[ring[i].atomno[j]].ar4++;
				continue_index = 1;
				break;
			}
		if (continue_index == 1) {
			continue_index = 0;
			continue;
		}

/* for planar rings formed "outside" double bonds */
		if (tmpint >= ring[i].num && tmpint <= 2 * ring[i].num)
			for (j = 0; j < bondnum; j++) {
				index = 0;
				for (k = 0; k < ring[i].num; k++)
					if (bond[j].bondi == ring[i].atomno[k])
						if (arom[bond[j].bondj].rg[0] == 0)
							index++;
				for (k = 0; k < ring[i].num; k++)
					if (bond[j].bondj == ring[i].atomno[k])
						if (arom[bond[j].bondi].rg[0] == 0)
							index++;
				if (index == 1 && (bond[j].type == 2 || bond[j].type == 8)) {
					for (j = 0; j < ring[i].num; j++)
						arom[ring[i].atomno[j]].ar3++;
					continue_index = 1;
					break;
				}
			}
		if (continue_index == 1) {
			continue_index = 0;
			continue;
		}

/* for pure aromatic rings*/
		if (tmpint == 12 && ring[i].num == 6) {

			index = 0;
			for (j = 0; j < ring[i].num; j++)
				if (atm[ring[i].atomno[j]].atomicnum == 7
					|| atm[ring[i].atomno[j]].atomicnum == 15) {
					tmpint1 = ring[i].atomno[j];
					index0 = 0;
					for (k = 0; k < bondnum; k++) {
						if (bond[k].bondi == tmpint1
							&& (bond[k].type == 8 || bond[k].type == 2 || bond[k].type == 10))
							index0 = 1;
						if (bond[k].bondj == tmpint1
							&& (bond[k].type == 8 || bond[k].type == 2 || bond[k].type == 10))
							index0 = 1;
					}
					if (index0 == 0)
						index = 1;
				}
                        if(index == 0) {
                                for (j = 0; j < ring[i].num; j++)
                                        arom[ring[i].atomno[j]].ar1++;
                                continue;
                        }

		}

/* for other planar rings*/
		if (tmpint >= ring[i].num + 3) {
			for (j = 0; j < ring[i].num; j++)
				arom[ring[i].atomno[j]].ar2++;
			continue;
		}
/* for the other rings */
		for (j = 0; j < ring[i].num; j++)
			arom[ring[i].atomno[j]].ar4++;
	}
	for (i = 0; i < atnum; i++) {
		if (arom[i].ar1 > 0) {
			arom[i].nr = 0;
			continue;
		}
		if (arom[i].ar2 > 0) {
			arom[i].nr = 0;
			continue;
		}
		if (arom[i].ar3 > 0) {
			arom[i].nr = 0;
			continue;
		}
		if (arom[i].ar4 > 0) {
			arom[i].nr = 0;
			continue;
		}
		if (arom[i].ar5 > 0) {
			arom[i].nr = 0;
			continue;
		}
	}
	free(initarom);
}

int ringdetect(int atnum, ATOM atm[], int bondnum, BOND bond[],
			   int *ringnum, RING ring[], AROM arom[], int maxatom, int maxring, 
			   char *filename, int index)
{
	int i, j, k;
	int num = 1;
	int *selectaromindex;
	int *selectarom;
	int id, flag;
	FILE *fptest;
	selectaromindex = (int *) malloc(sizeof(int) * maxatom);
	if (selectaromindex == NULL) {
		fprintf(stdout, "memory allocation error for *selectaromindex\n");
		exit(1);
	}
	selectarom = (int *) malloc(sizeof(int) * maxatom);
	if (selectarom == NULL) {
		fprintf(stdout, "memory allocation error for *selectarom\n");
		exit(1);
	}
	if ((fptest = fopen(filename, "w")) == NULL) {
		fprintf(stdout, "Cannot open a file %s to write in ringdetect(), exit\n", filename);
		exit(1);
	}
	for (i = 0; i < atnum; i++) {
		arom[i].nr = 1;
		for (j = 0; j <= 10; j++)
			arom[i].rg[j] = 0;
	}
	for (i = 0; i < maxring; i++) {
		for (j = 0; j < 10; j++)
			ring[i].atomno[j] = -1;
		ring[i].num = 0;
	}
	*ringnum = 0;
	for (i = 0; i < atnum; i++) {
		if (atm[i].atomicnum != 6 && atm[i].atomicnum != 7 &&
			atm[i].atomicnum != 8 && atm[i].atomicnum != 16 &&
			atm[i].atomicnum != 15)
			continue;
		if (atm[i].atomicnum == 6 && atm[i].connum <= 2)
			continue;
		if (atm[i].atomicnum == 8 && atm[i].connum == 1)
			continue;
		if (atm[i].atomicnum == 16 && atm[i].connum == 1)
			continue;
		for (j = 0; j < atnum; j++)
			selectaromindex[j] = -1;
		for (j = 0; j < 6; j++)
			selectarom[j] = -1;
		ring_detect_cycle(atm, 0, i, ringnum, ring, selectarom, selectaromindex,0);
	}
	if (*ringnum > maxring)	 {
		printf("\nInfo: the actual number of rings (%d) exceeds the defaut ring size (%d), reallocate memory automatically", *ringnum, maxring);
		return 1;
	}
	/*reallocate memory*/

	*ringnum = 0;
        for (i = 0; i < atnum; i++) {
                if (atm[i].atomicnum != 6 && atm[i].atomicnum != 7 &&
                        atm[i].atomicnum != 8 && atm[i].atomicnum != 16 &&
                        atm[i].atomicnum != 15)
                        continue;
                if (atm[i].atomicnum == 6 && atm[i].connum <= 2)
                        continue;
                if (atm[i].atomicnum == 8 && atm[i].connum == 1)
                        continue;
                if (atm[i].atomicnum == 16 && atm[i].connum == 1)
                        continue;
                for (j = 0; j < atnum; j++)
                        selectaromindex[j] = -1;
                for (j = 0; j < 6; j++)
                        selectarom[j] = -1;
                ring_detect_cycle(atm, 0, i, ringnum, ring, selectarom, selectaromindex,1);
                purify(atm, *ringnum, ring);
        }
	purify(atm, *ringnum, ring);
	ringproperty(*ringnum, ring, arom);
	aromatic(atnum, atm, bondnum, bond, *ringnum, maxatom, ring, arom);

/*for am1-bcc, five-memberred ring in indole-like system is not considered as aromatic ring */
	if (index == 1) {
		id = -1;
		for (i = 0; i < atnum; i++)
			if (arom[i].rg[5] >= 1 && arom[i].rg[6] >= 1
				&& ((arom[i].ar1 + arom[i].ar2) >= 2)) {
				id = i;
				for (j = 0; j < *ringnum; j++)
					if (ring[j].num == 5) {
						flag = 0;
						for (k = 0; k < 5; k++)
							if (id == ring[j].atomno[k]) {
								flag = 1;
								break;
							}
						if (flag == 1)
							for (k = 0; k < 5; k++)
								arom[ring[j].atomno[k]].ar2--;
					}
				id = -1;
			}
	}
	fprintf(fptest,
			"====================================================================================");
	fprintf(fptest,
			"\n-------------------------------ring property (I)------------------------------------\n");

	for (i = 0; i < *ringnum; i++) {
		if (ring[i].num == 0)
			continue;
		fprintf(fptest, "\n               Ring %d               \n\n",
				num++);
		for (j = 0; j < ring[i].num; j++)
			fprintf(fptest, "%5d %5d %5d %5s\n", i + 1, j + 1,
					ring[i].atomno[j] + 1, atm[ring[i].atomno[j]].name);
	}

	fprintf(fptest,
			"\n====================================================================================");
	fprintf(fptest,
			"\n-------------------------------ring property (II)------------------------------------\n");
	for (i = 0; i < *ringnum; i++) {
		if (ring[i].num == 0)
			continue;
		for (j = 0; j < ring[i].num; j++)
			fprintf(fptest,
					"atom[%2d] (%-4s) belongs to one member of %d-membered  ring (No %2d)\n",
					ring[i].atomno[j] + 1, atm[ring[i].atomno[j]].name,
					ring[i].num, i + 1);
	}
	fprintf(fptest,
			"\n====================================================================================");
	fprintf(fptest,
			"\n-------------------------------ring property (III)-----------------------------------");
	for (i = 0; i < atnum; i++) {
		for (j = 1; j <= 10; j++)
			if (arom[i].rg[j] > 0)
				fprintf(fptest,
						"\natom[%2d] (%-4s) involves in %d %d-member of ring(s)",
						i + 1, atm[i].name, arom[i].rg[j], j);

		if (arom[i].nr > 0)
			fprintf(fptest,
					"\natom[%2d] (%-4s) is not in any ring (nr[%d]=%d)",
					i + 1, atm[i].name, i + 1, arom[i].nr);
	}

	fprintf(fptest,
			"\n\n====================================================================================");
	fprintf(fptest,
			"\n-------------------------------aromatic property------------------------------------");

	for (i = 0; i < atnum; i++) {
		if (arom[i].ar1 >= 1)
			fprintf(fptest,
					"\natom[%2d] (%-4s) is in %d pure aromatic ring(s) (AR1)",
					i + 1, atm[i].name, arom[i].ar1);
		if (arom[i].ar2 >= 1)
			fprintf(fptest,
					"\natom[%2d] (%-4s) is in %d planar ring(s) (AR2)",
					i + 1, atm[i].name, arom[i].ar2);
		if (arom[i].ar3 >= 1)
			fprintf(fptest,
					"\natom[%2d] (%-4s) is in %d planar ring(s), which has/have \"outside\" bonds (AR3)",
					i + 1, atm[i].name, arom[i].ar3);
		if (arom[i].ar4 >= 1)
			fprintf(fptest,
					"\natom[%2d] (%-4s) is in %d non-planar ring(s) (AR4)",
					i + 1, atm[i].name, arom[i].ar4);
		if (arom[i].ar5 >= 1)
			fprintf(fptest,
					"\natom[%2d] (%-4s) is in %d pure aliphatic ring, which is/are made of sp3 carbons (AR5)",
					i + 1, atm[i].name, arom[i].ar5);
	}
	fprintf(fptest,
			"\n\n====================================================================================");
	fprintf(fptest,
			"\n-------------------------------electronic property----------------------------------");

	for (i = 0; i < atnum; i++)
		if (atm[i].ewd == 1)
			fprintf(fptest,
					"\natom [%2d] (%-4s) is an electron-withdrew atom",
					i + 1, atm[i].name);
		else
			fprintf(fptest,
					"\natom [%2d] (%-4s) is not an electron-withdrew atom",
					i + 1, atm[i].name);
	fprintf(fptest,
			"\n\n====================================================================================");
	fprintf(fptest,
			"\n--------------------------------connectivity property-------------------------------");

	for (i = 0; i < atnum; i++) {
		if (atm[i].con[0] != -1)
			fprintf(fptest, "\natom[%2d] (%-4s) %5d %5s", i + 1,
					atm[i].name, atm[i].con[0] + 1,
					atm[atm[i].con[0]].name);
		if (atm[i].con[1] != -1)
			fprintf(fptest, "\natom[%2d] (%-4s) %5d %5s", i + 1,
					atm[i].name, atm[i].con[1] + 1,
					atm[atm[i].con[1]].name);
		if (atm[i].con[2] != -1)
			fprintf(fptest, "\natom[%2d] (%-4s) %5d %5s", i + 1,
					atm[i].name, atm[i].con[2] + 1,
					atm[atm[i].con[2]].name);
		if (atm[i].con[3] != -1)
			fprintf(fptest, "\natom[%2d] (%-4s) %5d %5s", i + 1,
					atm[i].name, atm[i].con[3] + 1,
					atm[atm[i].con[3]].name);
		if (atm[i].con[4] != -1)
			fprintf(fptest, "\natom[%2d] (%-4s) %5d %5s", i + 1,
					atm[i].name, atm[i].con[4] + 1,
					atm[atm[i].con[4]].name);
		if (atm[i].con[5] != -1)
			fprintf(fptest, "\natom[%2d] (%-4s) %5d %5s", i + 1,
					atm[i].name, atm[i].con[5] + 1,
					atm[atm[i].con[5]].name);
	}
	fprintf(fptest,
			"\n----------------------------------------END-----------------------------------------\n\n");
	fclose(fptest);
	return 0;
}

int ringdetect2(int atnum, ATOM atm[], int bondnum, BOND bond[],
			   int *ringnum, RING ring[], AROM arom[], int maxatom, int maxring, 
			   char *filename, int index)
{
	int i, j, k;
	int num = 1;
	int *selectaromindex;
	int *selectarom;
	int id, flag;
	int bondi, bondj, type;
	int flag1, flag2;
	FILE *fptest;
	selectaromindex = (int *) malloc(sizeof(int) * maxatom);
	if (selectaromindex == NULL) {
		fprintf(stdout, "memory allocation error for *selectaromindex\n");
		exit(1);
	}
	selectarom = (int *) malloc(sizeof(int) * maxatom);
	if (selectarom == NULL) {
		fprintf(stdout, "memory allocation error for *selectarom\n");
		exit(1);
	}
	if ((fptest = fopen(filename, "w")) == NULL) {
		fprintf(stdout, "Cannot open file %s to write in ringdetect2(), exit\n", filename);
		exit(1);
	}
	for (i = 0; i < atnum; i++) {
		arom[i].nr = 1;
		for (j = 0; j <= 10; j++)
			arom[i].rg[j] = 0;
	}
	for (i = 0; i < maxring; i++) {
		for (j = 0; j < 10; j++)
			ring[i].atomno[j] = -1;
		ring[i].num = 0;
	}
	*ringnum = 0;
	for (i = 0; i < atnum; i++) {
		if (atm[i].atomicnum != 6 && atm[i].atomicnum != 7 &&
			atm[i].atomicnum != 8 && atm[i].atomicnum != 16 &&
			atm[i].atomicnum != 15)
			continue;
		if (atm[i].atomicnum == 6 && atm[i].connum <= 2)
			continue;
		if (atm[i].atomicnum == 8 && atm[i].connum == 1)
			continue;
		if (atm[i].atomicnum == 16 && atm[i].connum == 1)
			continue;
		for (j = 0; j < atnum; j++)
			selectaromindex[j] = -1;
		for (j = 0; j < 6; j++)
			selectarom[j] = -1;
		ring_detect_cycle(atm, 0, i, ringnum, ring, selectarom, selectaromindex,0);
	}
	if (*ringnum > maxring)	 {
		printf("\nInfo: the actual number of rings (%d) exceeds the defaut ring size (%d), reallocate memory automatically", *ringnum, maxring);
		return 1;
	}
	/*reallocate memory*/

	*ringnum = 0;
        for (i = 0; i < atnum; i++) {
                if (atm[i].atomicnum != 6 && atm[i].atomicnum != 7 &&
                        atm[i].atomicnum != 8 && atm[i].atomicnum != 16 &&
                        atm[i].atomicnum != 15)
                        continue;
                if (atm[i].atomicnum == 6 && atm[i].connum <= 2)
                        continue;
                if (atm[i].atomicnum == 8 && atm[i].connum == 1)
                        continue;
                if (atm[i].atomicnum == 16 && atm[i].connum == 1)
                        continue;
                for (j = 0; j < atnum; j++)
                        selectaromindex[j] = -1;
                for (j = 0; j < 6; j++)
                        selectarom[j] = -1;
                ring_detect_cycle(atm, 0, i, ringnum, ring, selectarom, selectaromindex,1);
                purify(atm, *ringnum, ring);
        }
	purify(atm, *ringnum, ring);
	ringproperty(*ringnum, ring, arom);
	aromatic(atnum, atm, bondnum, bond, *ringnum, maxatom, ring, arom);

/*for am1-bcc, five-memberred ring in indole-like system is not considered as aromatic ring */
	if (index == 1) {
		id = -1;
		for (i = 0; i < atnum; i++)
			if (arom[i].rg[5] >= 1 && arom[i].rg[6] >= 1
				&& ((arom[i].ar1 + arom[i].ar2) >= 2)) {
				id = i;
				for (j = 0; j < *ringnum; j++)
					if (ring[j].num == 5) {
						flag = 0;
						for (k = 0; k < 5; k++)
							if (id == ring[j].atomno[k]) {
								flag = 1;
								break;
							}
						if (flag == 1)
							for (k = 0; k < 5; k++)
								arom[ring[j].atomno[k]].ar2--;
					}
				id = -1;
			}
	}
	
/* output atom and bond types */
	for (i = 0; i < atnum; i++) {
		atm[i].type2 = 0;
		atm[i].type3 = 0;

	}
	for (i = 0; i < *ringnum; i++) {
		if (ring[i].num == 0)
			continue;
		for (j = 0; j < ring[i].num; j++)
			atm[ring[i].atomno[j]].type3 = 1;
	}

	for (i = 0; i< bondnum; i++) {
		bondi = bond[i].bondi;
		bondj = bond[i].bondj;
		type = bond[i].type ;
		bond[i].type2 = 0;

/* No H, F, Cl, Br, I etc */
		if(atm[bondi].connum == 1 || atm[bondj].connum == 1) continue; 
/* No double triple and aromatic bonds */
		if(type == 2 || type == 3 || type == 8 || type == 9 || type == 10) continue; 
		

/* exclude some single bonds in a functional group */
/* O=C-N , O=C-O*/
	if(atm[bondi].connum == 3 && atm[bondi].atomicnum == 6) {
		flag1 = 0;
		flag2 = 0;	
		if(flag1==0 && atm[atm[bondi].con[0]].atomicnum == 8 && atm[atm[bondi].con[0]].connum == 1) flag1 = 1;
		if(flag1==0 && atm[atm[bondi].con[1]].atomicnum == 8 && atm[atm[bondi].con[1]].connum == 1) flag1 = 1;
		if(flag1==0 && atm[atm[bondi].con[2]].atomicnum == 8 && atm[atm[bondi].con[2]].connum == 1) flag1 = 1;

		if(flag2==0 && atm[atm[bondi].con[0]].atomicnum == 7 && atm[atm[bondi].con[0]].connum == 3) flag2 = 1;
		if(flag2==0 && atm[atm[bondi].con[1]].atomicnum == 7 && atm[atm[bondi].con[1]].connum == 3) flag2 = 1;
		if(flag2==0 && atm[atm[bondi].con[2]].atomicnum == 7 && atm[atm[bondi].con[2]].connum == 3) flag2 = 1;
		if(flag2==0 && atm[atm[bondi].con[0]].atomicnum == 8 && atm[atm[bondi].con[0]].connum == 2) flag2 = 1;
		if(flag2==0 && atm[atm[bondi].con[1]].atomicnum == 8 && atm[atm[bondi].con[1]].connum == 2) flag2 = 1;
		if(flag2==0 && atm[atm[bondi].con[2]].atomicnum == 8 && atm[atm[bondi].con[2]].connum == 2) flag2 = 1;
		if(flag1 + flag2 == 2) continue;
	}
	if(atm[bondj].connum == 3 && atm[bondj].atomicnum == 6) {
		flag1 = 0;
		flag2 = 0;	
		if(flag1==0 && atm[atm[bondj].con[0]].atomicnum == 8 && atm[atm[bondj].con[0]].connum == 1) flag1 = 1;
		if(flag1==0 && atm[atm[bondj].con[1]].atomicnum == 8 && atm[atm[bondj].con[1]].connum == 1) flag1 = 1;
		if(flag1==0 && atm[atm[bondj].con[2]].atomicnum == 8 && atm[atm[bondj].con[2]].connum == 1) flag1 = 1;

		if(flag2==0 && atm[atm[bondj].con[0]].atomicnum == 7 && atm[atm[bondj].con[0]].connum == 3) flag2 = 1;
		if(flag2==0 && atm[atm[bondj].con[1]].atomicnum == 7 && atm[atm[bondj].con[1]].connum == 3) flag2 = 1;
		if(flag2==0 && atm[atm[bondj].con[2]].atomicnum == 7 && atm[atm[bondj].con[2]].connum == 3) flag2 = 1;
		if(flag2==0 && atm[atm[bondj].con[0]].atomicnum == 8 && atm[atm[bondj].con[0]].connum == 2) flag2 = 1;
		if(flag2==0 && atm[atm[bondj].con[1]].atomicnum == 8 && atm[atm[bondj].con[1]].connum == 2) flag2 = 1;
		if(flag2==0 && atm[atm[bondj].con[2]].atomicnum == 8 && atm[atm[bondj].con[2]].connum == 2) flag2 = 1;
		if(flag1 + flag2 == 2) continue;
	}
/* O-C#N, S-C#N */
        if(atm[bondi].connum == 2 && atm[bondi].atomicnum == 6) {
                flag1 = 0;
                flag2 = 0;
                if(flag1==0 && atm[atm[bondi].con[0]].atomicnum == 7 && atm[atm[bondi].con[0]].connum == 1) flag1 = 1;
                if(flag1==0 && atm[atm[bondi].con[1]].atomicnum == 7 && atm[atm[bondi].con[1]].connum == 1) flag1 = 1;
                if(flag1==0 && atm[atm[bondi].con[2]].atomicnum == 7 && atm[atm[bondi].con[2]].connum == 1) flag1 = 1;

                if(flag2==0 && atm[atm[bondi].con[0]].atomicnum == 8 && atm[atm[bondi].con[0]].connum == 2) flag2 = 1;
                if(flag2==0 && atm[atm[bondi].con[1]].atomicnum == 8 && atm[atm[bondi].con[1]].connum == 2) flag2 = 1;
                if(flag2==0 && atm[atm[bondi].con[2]].atomicnum == 8 && atm[atm[bondi].con[2]].connum == 2) flag2 = 1;
                if(flag2==0 && atm[atm[bondi].con[0]].atomicnum == 16 && atm[atm[bondi].con[0]].connum == 2) flag2 = 1;
                if(flag2==0 && atm[atm[bondi].con[1]].atomicnum == 16 && atm[atm[bondi].con[1]].connum == 2) flag2 = 1;
                if(flag2==0 && atm[atm[bondi].con[2]].atomicnum == 16 && atm[atm[bondi].con[2]].connum == 2) flag2 = 1;
                if(flag1 + flag2 == 2) continue;
        }
        if(atm[bondj].connum == 2 && atm[bondj].atomicnum == 6) {
                flag1 = 0;
                flag2 = 0;
                if(flag1==0 && atm[atm[bondj].con[0]].atomicnum == 7 && atm[atm[bondj].con[0]].connum == 1) flag1 = 1;
                if(flag1==0 && atm[atm[bondj].con[1]].atomicnum == 7 && atm[atm[bondj].con[1]].connum == 1) flag1 = 1;
                if(flag1==0 && atm[atm[bondj].con[2]].atomicnum == 7 && atm[atm[bondj].con[2]].connum == 1) flag1 = 1;

                if(flag2==0 && atm[atm[bondj].con[0]].atomicnum == 8 && atm[atm[bondj].con[0]].connum == 2) flag2 = 1;
                if(flag2==0 && atm[atm[bondj].con[1]].atomicnum == 8 && atm[atm[bondj].con[1]].connum == 2) flag2 = 1;
                if(flag2==0 && atm[atm[bondj].con[2]].atomicnum == 8 && atm[atm[bondj].con[2]].connum == 2) flag2 = 1;
                if(flag2==0 && atm[atm[bondj].con[0]].atomicnum == 16 && atm[atm[bondj].con[0]].connum == 2) flag2 = 1;
                if(flag2==0 && atm[atm[bondj].con[1]].atomicnum == 16 && atm[atm[bondj].con[1]].connum == 2) flag2 = 1;
                if(flag2==0 && atm[atm[bondj].con[2]].atomicnum == 16 && atm[atm[bondj].con[2]].connum == 2) flag2 = 1;
                if(flag1 + flag2 == 2) continue;
        }
/*S-S or O-O */
        if(atm[bondi].connum == 2 && atm[bondi].atomicnum == 16 && atm[bondj].connum == 2 && atm[bondj].atomicnum == 16)
		continue; 
        if(atm[bondi].connum == 2 && atm[bondi].atomicnum == 8 && atm[bondj].connum == 2 && atm[bondj].atomicnum == 8)
		continue; 

/* O=N-O*/
	if(atm[bondi].connum == 2 && atm[bondi].atomicnum == 7) {
		flag1 = 0;
		flag2 = 0;	
		if(flag1==0 && atm[atm[bondi].con[0]].atomicnum == 8 && atm[atm[bondi].con[0]].connum == 1) flag1 = 1;
		if(flag1==0 && atm[atm[bondi].con[1]].atomicnum == 8 && atm[atm[bondi].con[1]].connum == 1) flag1 = 1;

		if(flag2==0 && atm[atm[bondi].con[0]].atomicnum == 8 && atm[atm[bondi].con[0]].connum == 2) flag2 = 1;
		if(flag2==0 && atm[atm[bondi].con[1]].atomicnum == 8 && atm[atm[bondi].con[1]].connum == 2) flag2 = 1;
		if(flag1 + flag2 == 2) continue;
	}
	if(atm[bondj].connum == 2 && atm[bondj].atomicnum == 7) {
		flag1 = 0;
		flag2 = 0;	
		if(flag1==0 && atm[atm[bondj].con[0]].atomicnum == 8 && atm[atm[bondj].con[0]].connum == 1) flag1 = 1;
		if(flag1==0 && atm[atm[bondj].con[1]].atomicnum == 8 && atm[atm[bondj].con[1]].connum == 1) flag1 = 1;

		if(flag2==0 && atm[atm[bondj].con[0]].atomicnum == 8 && atm[atm[bondj].con[0]].connum == 2) flag2 = 1;
		if(flag2==0 && atm[atm[bondj].con[1]].atomicnum == 8 && atm[atm[bondj].con[1]].connum == 2) flag2 = 1;
		if(flag1 + flag2 == 2) continue;
	}

	if(atm[bondi].connum == 3 && atm[bondi].atomicnum == 7) {
		flag1 = 0;
		flag2 = 0;	
		if(flag1==0 && atm[atm[bondi].con[0]].atomicnum == 8 && atm[atm[bondi].con[0]].connum == 1) flag1 = 1;
		if(flag1==0 && atm[atm[bondi].con[1]].atomicnum == 8 && atm[atm[bondi].con[1]].connum == 1) flag1 = 1;
		if(flag1==0 && atm[atm[bondi].con[2]].atomicnum == 8 && atm[atm[bondi].con[2]].connum == 1) flag1 = 1;

		if(flag2==0 && atm[atm[bondi].con[0]].atomicnum == 8 && atm[atm[bondi].con[0]].connum == 2) flag2 = 1;
		if(flag2==0 && atm[atm[bondi].con[1]].atomicnum == 8 && atm[atm[bondi].con[1]].connum == 2) flag2 = 1;
		if(flag2==0 && atm[atm[bondi].con[2]].atomicnum == 8 && atm[atm[bondi].con[2]].connum == 2) flag2 = 1;
		if(flag1 + flag2 == 2) continue;
	}
	if(atm[bondj].connum == 3 && atm[bondj].atomicnum == 7) {
		flag1 = 0;
		flag2 = 0;	
		if(flag1==0 && atm[atm[bondj].con[0]].atomicnum == 8 && atm[atm[bondj].con[0]].connum == 1) flag1 = 1;
		if(flag1==0 && atm[atm[bondj].con[1]].atomicnum == 8 && atm[atm[bondj].con[1]].connum == 1) flag1 = 1;
		if(flag1==0 && atm[atm[bondj].con[2]].atomicnum == 8 && atm[atm[bondj].con[2]].connum == 1) flag1 = 1;

		if(flag2==0 && atm[atm[bondj].con[0]].atomicnum == 8 && atm[atm[bondj].con[0]].connum == 2) flag2 = 1;
		if(flag2==0 && atm[atm[bondj].con[1]].atomicnum == 8 && atm[atm[bondj].con[1]].connum == 2) flag2 = 1;
		if(flag2==0 && atm[atm[bondj].con[2]].atomicnum == 8 && atm[atm[bondj].con[2]].connum == 2) flag2 = 1;
		if(flag1 + flag2 == 2) continue;
	}
/* O=P-O*/
	if(atm[bondi].connum == 5 && atm[bondi].atomicnum == 15) {
		flag1 = 0;
		flag2 = 0;	
		if(flag1==0 && atm[atm[bondi].con[0]].atomicnum == 8 && atm[atm[bondi].con[0]].connum == 1) flag1 = 1;
		if(flag1==0 && atm[atm[bondi].con[1]].atomicnum == 8 && atm[atm[bondi].con[1]].connum == 1) flag1 = 1;
		if(flag1==0 && atm[atm[bondi].con[2]].atomicnum == 8 && atm[atm[bondi].con[2]].connum == 1) flag1 = 1;
		if(flag1==0 && atm[atm[bondi].con[3]].atomicnum == 8 && atm[atm[bondi].con[3]].connum == 1) flag1 = 1;
		if(flag1==0 && atm[atm[bondi].con[4]].atomicnum == 8 && atm[atm[bondi].con[4]].connum == 1) flag1 = 1;

		if(flag2==0 && atm[atm[bondi].con[0]].atomicnum == 8 && atm[atm[bondi].con[0]].connum == 2) flag2 = 1;
		if(flag2==0 && atm[atm[bondi].con[1]].atomicnum == 8 && atm[atm[bondi].con[1]].connum == 2) flag2 = 1;
		if(flag2==0 && atm[atm[bondi].con[2]].atomicnum == 8 && atm[atm[bondi].con[2]].connum == 2) flag2 = 1;
		if(flag2==0 && atm[atm[bondi].con[3]].atomicnum == 8 && atm[atm[bondi].con[3]].connum == 2) flag2 = 1;
		if(flag2==0 && atm[atm[bondi].con[4]].atomicnum == 8 && atm[atm[bondi].con[4]].connum == 2) flag2 = 1;
		if(flag1 + flag2 == 2) continue;
	}
	if(atm[bondj].connum == 5 && atm[bondj].atomicnum == 7) {
		flag1 = 0;
		flag2 = 0;	
		if(flag1==0 && atm[atm[bondj].con[0]].atomicnum == 8 && atm[atm[bondj].con[0]].connum == 1) flag1 = 1;
		if(flag1==0 && atm[atm[bondj].con[1]].atomicnum == 8 && atm[atm[bondj].con[1]].connum == 1) flag1 = 1;
		if(flag1==0 && atm[atm[bondj].con[2]].atomicnum == 8 && atm[atm[bondj].con[2]].connum == 1) flag1 = 1;
		if(flag1==0 && atm[atm[bondj].con[3]].atomicnum == 8 && atm[atm[bondj].con[3]].connum == 1) flag1 = 1;
		if(flag1==0 && atm[atm[bondj].con[4]].atomicnum == 8 && atm[atm[bondj].con[4]].connum == 1) flag1 = 1;

		if(flag2==0 && atm[atm[bondj].con[0]].atomicnum == 8 && atm[atm[bondj].con[0]].connum == 2) flag2 = 1;
		if(flag2==0 && atm[atm[bondj].con[1]].atomicnum == 8 && atm[atm[bondj].con[1]].connum == 2) flag2 = 1;
		if(flag2==0 && atm[atm[bondj].con[2]].atomicnum == 8 && atm[atm[bondj].con[2]].connum == 2) flag2 = 1;
		if(flag2==0 && atm[atm[bondj].con[3]].atomicnum == 8 && atm[atm[bondj].con[3]].connum == 2) flag2 = 1;
		if(flag2==0 && atm[atm[bondj].con[4]].atomicnum == 8 && atm[atm[bondj].con[4]].connum == 2) flag2 = 1;
		if(flag1 + flag2 == 2) continue;
	}
/* bondi and bondj are not both ring atoms */
	if((atm[bondi].type3 + atm[bondj].type3) <= 1) {
		bond[i].type2 = 1;
		atm[bondi].type2 = 1;
		atm[bondj].type2 = 1;
		continue;
	}
	num = 0;
	for (j = 0; j < *ringnum; j++) {
		if (ring[j].num == 0) continue;
		num = 0;	
		for (k = 0; k < ring[j].num; k++) {
			if(ring[j].atomno[k] == bondi || ring[j].atomno[k] == bondj) 
				num++;
		}
		if(num == 2) break;
	}
	if(num != 2) {  /* bond that links two rings*/
		bond[i].type2 = 1;
		atm[bondi].type2 = 1;
		atm[bondj].type2 = 1;
	}
}
for (i = 0; i < atnum; i++) 
	fprintf(fptest, "\nAtom[%2d](%-4s)\t%5d\t%5d\t%5d", i + 1, atm[i].name, atm[i].type, atm[i].type2, atm[i].type3);
fprintf(fptest,"\n");
for (i = 0; i < bondnum; i++) 
	fprintf(fptest, "\nBond[%2d]: %5d(%-4s)-%5d(%-4s)\t%5d\t%5d\n", i + 1, bond[i].bondi + 1, atm[bond[i].bondi].name, 
			bond[i].bondj+1, atm[bond[i].bondj].name, bond[i].type, bond[i].type2);
fclose(fptest);
return 0;
}
