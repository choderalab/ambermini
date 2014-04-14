void cycle(ATOM atm[], int selectnum, int startnum)
{
/* determain the aromatic atom through an iterative way */
	int i, j, k, l, m;
	int aromindex1;
	int aromindex2;
	int start;
	int breakindex = 0;
	start = -1;
	aromindex1 = 0;
	selectarom[selectnum++] = startnum;
	selectaromindex[startnum] = 1;
	for (i = 0; i < 4; i++) {
		if (atm[startnum].con[i] == -1)
			return;
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
		if (selectnum > 6)
			return;
		/* we have already visited this atom */
		for (j = 2; j < 6; j++) {	/* deal with 3-7 member ring */
			if (selectnum == j) {
				for (k = 0; k < 4; k++) {
					if (atm[selectarom[0]].con[k] == -1)
						break;
					if (atm[selectarom[0]].con[k] == start) {
						for (m = 0; m < j; m++)
							ring[start][ringindex[start]++] =
								selectarom[m];
						ring[start][ringindex[start]++] = -99;
						aromindex2 = 0;
						for (l = 0; l < selectnum; l++)
							if (atm[selectarom[l]].atomicnum == 6
								&& atm[selectarom[l]].connum == 4)
								aromindex2 = 1;
						if (atm[start].atomicnum == 6
							&& atm[start].connum == 4)
							aromindex2 = 1;
						if (aromindex2 == 0)
							atm[selectarom[0]].arom += j + 1;
						if (aromindex2 == 1)
							atm[selectarom[0]].aliph += (j + 1);
						aromindex1 = 1;
						break;
					}
				}
				if (aromindex1 == 1)
					break;
				if (aromindex2 == 1)
					break;
			}
		}
		if (aromindex1 == 1)
			return;
		if (start == -1)
			return;
		cycle(atm, selectnum, start);
	}
}

void aromatic(void)
{
	int i, j;
/* The following code judge weather the atom is electron withdrawing one or not */
	for (i = 0; i < atomnum; i++) {
		switch (atom[i].atomicnum) {
		case 6:
			atom[i].ewd = 0;
			if (atom[i].connum < 4)
				atom[i].saturate = -1;
			else
				atom[i].saturate = 1;
			break;
		case 7:
			atom[i].ewd = 1;
			if (atom[i].connum < 3)
				atom[i].saturate = -1;
			else
				atom[i].saturate = 1;
			break;
		case 8:
			atom[i].ewd = 1;
			if (atom[i].connum < 2)
				atom[i].saturate = -1;
			else
				atom[i].saturate = 1;
			break;
		case 16:
			atom[i].ewd = 1;	/* S is considered electron withdraw group */
			break;
		case 9:
			atom[i].ewd = 1;
			atom[i].saturate = 1;
			break;
		case 17:
			atom[i].ewd = 1;
			atom[i].saturate = 1;
			break;
		case 35:
			atom[i].ewd = 1;
			atom[i].saturate = 1;
			break;
		case 53:
			atom[i].ewd = 1;
			atom[i].saturate = 1;
			break;
		case  1:
			atom[i].ewd = 0;
			atom[i].saturate = 1;
			break;
		default:
			atom[i].ewd = 0;
		}
	}
/* The following code judge whether an atom is aromatic atom or not */
	for (i = 0; i < atomnum; i++) {
		atom[i].arom = 0;
		atom[i].aliph = 0;
	}
	for (i = 0; i < atomnum; i++) {
		if (atom[i].atomicnum != 6 && atom[i].atomicnum != 7 &&
			atom[i].atomicnum != 8 && atom[i].atomicnum != 16 &&
			atom[i].atomicnum != 15)
			continue;
		if (atom[i].atomicnum == 6 && atom[i].connum <= 2)
			continue;
		if (atom[i].atomicnum == 8 && atom[i].connum == 1)
			continue;
		if (atom[i].atomicnum == 16 && atom[i].connum == 1)
			continue;
		for (j = 0; j < atomnum; j++)
			selectaromindex[j] = -1;
		for (j = 0; j < 6; j++)
			selectarom[j] = -1;
		cycle(atom, 0, i);
	}

	fprintf(fptest, "              AROMATIC RING\n\n");
	for (i = 0; i < atomnum; i++) {
		if (atom[i].arom == 0)
			continue;
		else {
			atom[i].arom /= 2;
			fprintf(fptest,
					"atom[%2d] belongs to one member of %2d aromatic ring\n",
					i + 1, atom[i].arom);
		}
	}
	fprintf(fptest, "\n              ALIPHATIC RING \n\n");
	for (i = 0; i < atomnum; i++) {
		if (atom[i].aliph == 0)
			continue;
		else {
			atom[i].aliph /= 2;
			fprintf(fptest,
					"atom[%2d] belongs to one member of %2d aliphatic ring\n",
					i + 1, atom[i].aliph);
		}
	}
}
