void memory(int flag, int maxatom, int maxbond, int maxring)
{
	if (flag == 0) {
		atom = (ATOM *) malloc(sizeof(ATOM) * maxatom);
		if (atom == NULL) {
			fprintf(stdout, "meory allocation error for *atom\n");
			exit(1);
		}
		arom = (AROM *) malloc(sizeof(AROM) * maxatom);
		if (arom == NULL) {
			fprintf(stdout, "meory allocation error for *arom\n");
			exit(1);
		}
		bond = (BOND *) malloc(sizeof(BOND) * maxbond);
		if (bond == NULL) {
			fprintf(stdout, "meory allocation error for *bond\n");
			exit(1);
		}
		int i;
		for (i = 0; i < maxbond; ++i) {
			bond[i].jflag = -1; /* bond type has not been assigned */
		}
		ring = (RING *) malloc(sizeof(RING) * maxring);
		if (ring == NULL) {
			fprintf(stdout, "meory allocation error for *ring\n");
			exit(1);
		}
	}
	if (flag == 3) {
		free(atom);
		free(arom);
		free(bond);
		free(ring);
		atom = (ATOM *) malloc(sizeof(ATOM) * maxatom);
		if (atom == NULL) {
			fprintf(stdout, "meory allocation error for *atom\n");
			exit(1);
		}
		arom = (AROM *) malloc(sizeof(AROM) * maxatom);
		if (arom == NULL) {
			fprintf(stdout, "meory allocation error for *arom\n");
			exit(1);
		}
		bond = (BOND *) malloc(sizeof(BOND) * maxbond);
		if (bond == NULL) {
			fprintf(stdout, "meory allocation error for *bond\n");
			exit(1);
		}
		int i;
		for (i = 0; i < maxbond; ++i) {
			bond[i].jflag = -1; /* bond type has not been assigned */
		}
		ring = (RING *) malloc(sizeof(RING) * maxring);
		if (ring == NULL) {
			fprintf(stdout, "meory allocation error for *ring\n");
			exit(1);
		}
	}
	if (flag == 1) {
		free(atom);
		free(arom);
		free(bond);
		atom = (ATOM *) malloc(sizeof(ATOM) * maxatom);
		if (atom == NULL) {
			fprintf(stdout, "meory allocation error for *atom\n");
			exit(1);
		}
		arom = (AROM *) malloc(sizeof(AROM) * maxatom);
		if (arom == NULL) {
			fprintf(stdout, "meory allocation error for *arom\n");
			exit(1);
		}
		bond = (BOND *) malloc(sizeof(BOND) * maxbond);
		if (bond == NULL) {
			fprintf(stdout, "meory allocation error for *bond\n");
			exit(1);
		}
		int i;
		for (i = 0; i < maxbond; ++i) {
			bond[i].jflag = -1; /* bond type has not been assigned */
		}
	}
	if (flag == 2) {
		free(ring);
		ring = (RING *) malloc(sizeof(RING) * maxring);
		if (ring == NULL) {
			fprintf(stdout, "meory allocation error for *ring\n");
			exit(1);
		}
	}
}
