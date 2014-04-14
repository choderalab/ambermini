/*
************************************************************************
*           All Copyright Reserved!                                    *
*                                                                      *
*  Prog:    parmcal                                                    *
*  Version: version 1.0                                                *
*  Author:  Junmei Wang                                                *
*                                                                      *
*  Department of Pharmaceutical Chemistry                              *
*  School of Pharmacy                                                  *
*  University of California                                            *
*  San Francisco   CA 94143                                            *
*  Octomber, 2001                                                      *
************************************************************************
*/

#include <math.h>
#include "common.h"
#include "define.h"
#include "utility.c"

char *amberhome;
char filename[MAXCHAR];
char line[MAXCHAR];
int i, j, k, l;

int id1, id2, id3;
char name1[5], name2[5], name3[5];

/*H-1, C-2, N-3, O-4, F-5, Cl-6, Br-7, I-8, S-9, P-10*/
double radius[11];
double elecneg[11];
double refbondlength[11][11];
double bfkai[11][11];
double anglec[11];
double anglez[11];

double bondlength;
double blforce;
double bondlength1, bondlength2;
double bondangle;
double baforce;
long angleparmnum = 0;


FILE *fp, *fpout;

typedef struct {
	char name1[5];
	char name2[5];
	char name3[5];
	double angle;
	double force;
} ANGLE;
ANGLE angleparm[10 * MAXATOMTYPE * MAXATOMTYPE];

void readangleparm(char *filename)
{
	int mindex = -1;
	int bindex = 0;
	int aindex = 0;
	int num = 0;
	FILE *fp;
	char line[MAXCHAR];

	if ((fp = fopen(filename, "r")) == NULL) {
		printf("\n Cannot open the angle parameter file %s in readangleparm(), exit", filename);
		return;
	}
	for (;;) {
		if (fgets(line, MAXCHAR, fp) == NULL)
			break;
		num++;
		if (mindex == -1 && num == 1) {
			mindex = 1;
			continue;
		}
		if (mindex == 1 && spaceline(line) == 1) {
			mindex = 0;
			num = 0;
			bindex = -1;
			continue;
		}
		if (bindex == -1 && num == 1) {
			bindex = 1;
			continue;
		}
		if (bindex == 1 && spaceline(line) == 1) {
			bindex = 0;
			aindex = 1;
			continue;
		}
		if (aindex == 1 && spaceline(line) == 1) {
			aindex = 0;
			break;
		}
		if (aindex == 1) {
			angleparm[angleparmnum].name1[0] = line[0];
			angleparm[angleparmnum].name1[1] = line[1];
			angleparm[angleparmnum].name2[0] = line[3];
			angleparm[angleparmnum].name2[1] = line[4];
			angleparm[angleparmnum].name3[0] = line[6];
			angleparm[angleparmnum].name3[1] = line[7];
			sscanf(&line[8], "%lf%lf", &angleparm[angleparmnum].force,
				   &angleparm[angleparmnum].angle);
			angleparmnum++;
			if(angleparmnum >= 10 * MAXATOMTYPE * MAXATOMTYPE) {
				fprintf(stdout, "Error, the number of angle parameters (%ld) exceed 10 * MAXATOMTYPE * MAXATOMTYPE,\nincrease MAXATOMTYPE (%d) in define.h\n", angleparmnum, MAXATOMTYPE);
			}
		}
	}
	fclose(fp);
}

void bondlengthcal(void)
{
	int boid = 0;
	int atid1 = 0;
	int atid2 = 0;
	bondlength =
		radius[id1] + radius[id2] -
		0.05 * pow(fabs(elecneg[id1] - elecneg[id2]), 1.4) - 0.02;
	do {
		printf("Do you think the bond A-B is \n");
		printf("1. single bond\n");
		printf("2. pure double bond\n");
		printf("3. aromatic double  bond\n");
		printf("4. conjugated double  bond\n");
		printf("5. triple bond\n");
		scanf("%d", &boid);
	} while (boid != 1 && boid != 2 && boid != 3 && boid != 4
			 && boid != 5);
	if (boid == 2)
		bondlength -= 0.2;
	if (boid == 3)
		bondlength -= 0.15;
	if (boid == 4)
		bondlength -= 0.17;
	if (boid == 5)
		bondlength -= 0.34;
	if (boid == 1) {
		do {
			printf("For single bond, we need more information\n");
			printf("For Atom A, which is the best type to describe it?\n");
			printf(" 1. sp2 carbon \n");
			printf(" 2. sp2 carbon in carbonyl group (C=O)\n");
			printf(" 3. sp  carbon\n");
			printf(" 4. sp3 nitrogen connected to three atoms\n");
			printf(" 5. sp3 nitrogen connected to four atoms\n");
			printf(" 6. sp2 nitrogen in a double bond group (N=X)\n");
			printf(" 7. sp2 nitrogen in amind group\n");
			printf
				(" 8. sp2 nitrogen in an aromatic ring connected to three atoms\n");
			printf(" 9. sp2 nitrogen connected to an aromatic atom\n");
			printf("10. sp  nitrogen\n");
			printf("11. sp2 sulfur in a double bond group (S=X)\n");
			printf("12. sp3 sulfur\n");
			printf("13. OH oxygen in hydroxyl group\n");
			printf("14. OS oxygen in ester or ether\n");
			printf("15. None of the above\n");
			scanf("%d", &atid1);
		} while (atid1 != 1 && atid1 != 2 && atid1 != 3 && atid1 != 4
				 && atid1 != 5 && atid1 != 6 && atid1 != 7 && atid1 != 8
				 && atid1 != 9 && atid1 != 10 && atid1 != 11 && atid1 != 12
				 && atid1 != 13 && atid1 != 14 && atid1 != 15);
		do {
			printf("For single bond, we need more information\n");
			printf("For Atom B, which is the best type to describe it?\n");
			printf(" 1. sp2 carbon \n");
			printf(" 2. sp2 carbon in carbonyl group (C=O)\n");
			printf(" 3. sp  carbon\n");
			printf(" 4. sp3 nitrogen connected to three atoms\n");
			printf(" 5. sp3 nitrogen connected to four atoms\n");
			printf(" 6. sp2 nitrogen in a double bond group (N=X)\n");
			printf(" 7. sp2 nitrogen in amind group\n");
			printf
				(" 8. sp2 nitrogen in an aromatic ring connected to three atoms\n");
			printf(" 9. sp2 nitrogen connected to an aromatic atom\n");
			printf("10. sp  nitrogen\n");
			printf("11. sp2 sulfur in a double bond group (S=X)\n");
			printf("12. sp3 sulfur\n");
			printf("13. OH oxygen in hydroxyl group\n");
			printf("14. OS oxygen in ester or ether\n");
			printf("15. None of the above\n");
			scanf("%d", &atid2);
		} while (atid2 != 1 && atid2 != 2 && atid2 != 3 && atid2 != 4
				 && atid2 != 5 && atid2 != 6 && atid2 != 7 && atid2 != 8
				 && atid2 != 9 && atid2 != 10 && atid2 != 11 && atid2 != 12
				 && atid2 != 13 && atid2 != 14 && atid1 != 15);

		if (atid1 == 1)
			bondlength -= 0.03;
		if (atid1 == 2)
			bondlength -= 0.01;
		if (atid1 == 3)
			bondlength -= 0.06;
		if (atid1 == 6)
			bondlength -= 0.06;
		if (atid1 == 10)
			bondlength -= 0.08;
		if (atid1 == 11)
			bondlength -= 0.03;
		if (atid2 == 1)
			bondlength -= 0.03;
		if (atid2 == 2)
			bondlength -= 0.01;
		if (atid2 == 3)
			bondlength -= 0.06;
		if (atid2 == 6)
			bondlength -= 0.06;
		if (atid2 == 10)
			bondlength -= 0.08;
		if (atid2 == 11)
			bondlength -= 0.03;

		if (atid1 == 1 || atid1 == 2 || atid1 == 3 || atid1 == 6
			|| atid1 == 10 || atid1 == 11) {
			if (atid2 == 4)
				bondlength -= 0.04;
			if (atid2 == 5)
				bondlength -= 0.03;
			if (atid2 == 7)
				bondlength -= 0.05;
			if (atid2 == 8)
				bondlength -= 0.04;
			if (atid2 == 9)
				bondlength -= 0.05;
			if (atid2 == 12)
				bondlength -= 0.03;
			if (atid2 == 13)
				bondlength -= 0.06;
			if (atid2 == 14)
				bondlength -= 0.05;
		}

		if (atid2 == 1 || atid2 == 2 || atid2 == 3 || atid2 == 6
			|| atid2 == 10 || atid2 == 11) {
			if (atid1 == 4)
				bondlength -= 0.04;
			if (atid1 == 5)
				bondlength -= 0.03;
			if (atid1 == 7)
				bondlength -= 0.05;
			if (atid1 == 8)
				bondlength -= 0.04;
			if (atid1 == 9)
				bondlength -= 0.05;
			if (atid1 == 12)
				bondlength -= 0.03;
			if (atid1 == 13)
				bondlength -= 0.06;
			if (atid1 == 14)
				bondlength -= 0.05;
		}

	}
}
void bondlengthfc(void)
{
	blforce = exp(bfkai[id1][id2] - 4.5 * log(bondlength));
}

int bondanglecal(void)
{
	int i;
	int id;
	int num1 = -1;
	int num2 = -1;
	int index;
	char typename1[5];
	char typename2[5];
	char typename3[5];
	char parmname1[5];
	char parmname2[5];
	char parmname3[5];

	if (id1 == 1)
		strcpy(typename1, "h");
	if (id2 == 1)
		strcpy(typename2, "h");
	if (id3 == 1)
		strcpy(typename3, "h");
	if (id1 == 5)
		strcpy(typename1, "f");
	if (id2 == 5)
		strcpy(typename2, "f");
	if (id3 == 5)
		strcpy(typename3, "f");
	if (id1 == 6)
		strcpy(typename1, "cl");
	if (id2 == 6)
		strcpy(typename2, "cl");
	if (id3 == 6)
		strcpy(typename3, "cl");
	if (id1 == 7)
		strcpy(typename1, "br");
	if (id2 == 7)
		strcpy(typename2, "br");
	if (id3 == 7)
		strcpy(typename3, "br");
	if (id1 == 8)
		strcpy(typename1, "i");
	if (id2 == 8)
		strcpy(typename2, "i");
	if (id3 == 8)
		strcpy(typename3, "i");

	printf("For Atom A, which is the best type to describe it?\n");
	if (id1 == 2) {
		id = 0;
		do {
			printf(" 1. sp3 carbon \n");
			printf(" 2. sp2 carbon in carbonyl group (C=O)\n");
			printf(" 3. sp2 carbon (aromatic)\n");
			printf(" 4. sp2 carbon (non-aromatic)\n");
			printf(" 5. sp  carbon\n");
			scanf("%d", &id);
		}
		while (id != 1 && id != 2 && id != 3 && id != 4 && id != 5);
		if (id == 1)
			strcpy(typename1, "c3");
		if (id == 2)
			strcpy(typename1, "c");
		if (id == 3)
			strcpy(typename1, "ca");
		if (id == 4)
			strcpy(typename1, "c2");
		if (id == 5)
			strcpy(typename1, "c1");
	}
	if (id1 == 3) {
		id = 0;
		do {
			printf(" 1. sp3 nitrogen in amine (X3)\n");
			printf(" 2. sp3 nitrogen in amine (X4)\n");
			printf(" 3. sp2 nitrogen in double bond group (N=X)\n");
			printf(" 4. sp2 nitrogen nitro group (O=N=O)\n");
			printf(" 5. sp2 nitrogen in amind group\n");
			printf
				(" 6. sp2 nitrogen in an aromatic ring connected to three atoms\n");
			printf(" 7. sp2 nitrogen connected to an aromatic atom\n");
			printf(" 8. sp  nitrogen\n");
			scanf("%d", &id);
		}
		while (id != 1 && id != 2 && id != 3 && id != 4 && id != 5
			   && id != 6 && id != 7 && id != 7);
		if (id == 1)
			strcpy(typename1, "n3");
		if (id == 2)
			strcpy(typename1, "n4");
		if (id == 3)
			strcpy(typename1, "n2");
		if (id == 4)
			strcpy(typename1, "no");
		if (id == 5)
			strcpy(typename1, "n");
		if (id == 6)
			strcpy(typename1, "na");
		if (id == 7)
			strcpy(typename1, "nh");
	}

	if (id1 == 4) {
		id = 0;
		do {
			printf(" 1. OH oxygen in hydroxyl group\n");
			printf(" 2. OS oxygen in ester or ether\n");
			printf(" 3. sp2 oxygen (O=X)\n");
			scanf("%d", &id);
		}
		while (id != 1 && id != 2 && id != 3);
		if (id == 1)
			strcpy(typename1, "oh");
		if (id == 2)
			strcpy(typename1, "os");
		if (id == 3)
			strcpy(typename1, "o ");
	}

	if (id1 == 9) {
		id = 0;
		do {
			printf(" 1. SH sulfur in S-H\n");
			printf(" 2. SS sulfur in thio-ester or thio-ether or S-S\n");
			printf(" 3. sp2 sulfur in sulfone (-S=O)\n");
			printf(" 4. sp2 sulfur in sulfone (O=S=O)\n");
			printf(" 5. sp2 sulfur (S=X)\n");
			scanf("%d", &id);
		}
		while (id != 1 && id != 2 && id != 3 && id != 4 && id != 5);
		if (id == 1)
			strcpy(typename1, "sh");
		if (id == 2)
			strcpy(typename1, "ss");
		if (id == 3)
			strcpy(typename1, "s4");
		if (id == 4)
			strcpy(typename1, "s6");
		if (id == 5)
			strcpy(typename1, "s2");
	}
	if (id1 == 10) {
		id = 0;
		do {
			printf(" 1. P3 phosphate with three connected atoms\n");
			printf(" 2. P4 phosphate with five connected atoms\n");
			scanf("%d", &id);
		}
		while (id != 1 && id != 2 && id != 3);
		if (id == 1)
			strcpy(typename1, "p3");
		if (id == 2)
			strcpy(typename1, "p4");
	}


	printf("For Atom B, which is the best type to describe it?\n");
	if (id2 == 2) {
		id = 0;
		do {
			printf(" 1. sp3 carbon \n");
			printf(" 2. sp2 carbon in carbonyl group (C=O)\n");
			printf(" 3. sp2 carbon (aromatic)\n");
			printf(" 4. sp2 carbon (non-aromatic)\n");
			printf(" 5. sp  carbon\n");
			scanf("%d", &id);
		}
		while (id != 1 && id != 2 && id != 3 && id != 4 && id != 5);
		if (id == 1)
			strcpy(typename2, "c3");
		if (id == 2)
			strcpy(typename2, "c");
		if (id == 3)
			strcpy(typename2, "ca");
		if (id == 4)
			strcpy(typename2, "c2");
		if (id == 5)
			strcpy(typename2, "c1");
	}
	if (id2 == 3) {
		id = 0;
		do {
			printf(" 1. sp3 nitrogen in amine (X3)\n");
			printf(" 2. sp3 nitrogen in amine (X4)\n");
			printf(" 3. sp2 nitrogen in double bond group (N=X)\n");
			printf(" 4. sp2 nitrogen nitro group (O=N=O)\n");
			printf(" 5. sp2 nitrogen in amind group\n");
			printf
				(" 6. sp2 nitrogen in an aromatic ring connected to three atoms\n");
			printf(" 7. sp2 nitrogen connected to an aromatic atom\n");
			printf(" 8. sp  nitrogen\n");
			scanf("%d", &id);
		}
		while (id != 1 && id != 2 && id != 3 && id != 4 && id != 5
			   && id != 6 && id != 7 && id != 7);
		if (id == 1)
			strcpy(typename2, "n3");
		if (id == 2)
			strcpy(typename2, "n4");
		if (id == 3)
			strcpy(typename2, "n2");
		if (id == 4)
			strcpy(typename2, "no");
		if (id == 5)
			strcpy(typename2, "n");
		if (id == 6)
			strcpy(typename2, "na");
		if (id == 7)
			strcpy(typename2, "nh");
	}

	if (id2 == 4) {
		id = 0;
		do {
			printf(" 1. OH oxygen in hydroxyl group\n");
			printf(" 2. OS oxygen in ester or ether\n");
			printf(" 3. sp2 oxygen (O=X)\n");
			scanf("%d", &id);
		}
		while (id != 1 && id != 2 && id != 3);
		if (id == 1)
			strcpy(typename2, "oh");
		if (id == 2)
			strcpy(typename2, "os");
		if (id == 3)
			strcpy(typename2, "o ");
	}

	if (id2 == 9) {
		id = 0;
		do {
			printf(" 1. SH sulfur in S-H\n");
			printf(" 2. SS sulfur in thio-ester or thio-ether or S-S\n");
			printf(" 3. sp2 sulfur in sulfone (-S=O)\n");
			printf(" 4. sp2 sulfur in sulfone (O=S=O)\n");
			printf(" 5. sp2 sulfur (S=X)\n");
			scanf("%d", &id);
		}
		while (id != 1 && id != 2 && id != 3 && id != 4 && id != 5);
		if (id == 1)
			strcpy(typename2, "sh");
		if (id == 2)
			strcpy(typename2, "ss");
		if (id == 3)
			strcpy(typename2, "s4");
		if (id == 4)
			strcpy(typename2, "s6");
		if (id == 5)
			strcpy(typename2, "s2");
	}
	if (id2 == 10) {
		id = 0;
		do {
			printf(" 1. P3 phosphate with three connected atoms\n");
			printf(" 2. P4 phosphate with five connected atoms\n");
			scanf("%d", &id);
		}
		while (id != 1 && id != 2 && id != 3);
		if (id == 1)
			strcpy(typename2, "p3");
		if (id == 2)
			strcpy(typename2, "p4");
	}



	printf("For Atom C, which is the best type to describe it?\n");
	if (id3 == 2) {
		id = 0;
		do {
			printf(" 1. sp3 carbon \n");
			printf(" 2. sp2 carbon in carbonyl group (C=O)\n");
			printf(" 3. sp2 carbon (aromatic)\n");
			printf(" 4. sp2 carbon (non-aromatic)\n");
			printf(" 5. sp  carbon\n");
			scanf("%d", &id);
		}
		while (id != 1 && id != 2 && id != 3 && id != 4 && id != 5);
		if (id == 1)
			strcpy(typename3, "c3");
		if (id == 2)
			strcpy(typename3, "c");
		if (id == 3)
			strcpy(typename3, "ca");
		if (id == 4)
			strcpy(typename3, "c2");
		if (id == 5)
			strcpy(typename3, "c1");
	}
	if (id3 == 3) {
		id = 0;
		do {
			printf(" 1. sp3 nitrogen in amine (X3)\n");
			printf(" 2. sp3 nitrogen in amine (X4)\n");
			printf(" 3. sp2 nitrogen in double bond group (N=X)\n");
			printf(" 4. sp2 nitrogen nitro group (O=N=O)\n");
			printf(" 5. sp2 nitrogen in amind group\n");
			printf
				(" 6. sp2 nitrogen in an aromatic ring connected to three atoms\n");
			printf(" 7. sp2 nitrogen connected to an aromatic atom\n");
			printf(" 8. sp  nitrogen\n");
			scanf("%d", &id);
		}
		while (id != 1 && id != 2 && id != 3 && id != 4 && id != 5
			   && id != 6 && id != 7 && id != 7);
		if (id == 1)
			strcpy(typename3, "n3");
		if (id == 2)
			strcpy(typename3, "n4");
		if (id == 3)
			strcpy(typename3, "n2");
		if (id == 4)
			strcpy(typename3, "no");
		if (id == 5)
			strcpy(typename3, "n");
		if (id == 6)
			strcpy(typename3, "na");
		if (id == 7)
			strcpy(typename3, "nh");
	}

	if (id3 == 4) {
		id = 0;
		do {
			printf(" 1. OH oxygen in hydroxyl group\n");
			printf(" 2. OS oxygen in ester or ether\n");
			printf(" 3. sp2 oxygen (O=X)\n");
			scanf("%d", &id);
		}
		while (id != 1 && id != 2 && id != 3);
		if (id == 1)
			strcpy(typename3, "oh");
		if (id == 2)
			strcpy(typename3, "os");
		if (id == 3)
			strcpy(typename3, "o ");
	}

	if (id3 == 9) {
		id = 0;
		do {
			printf(" 1. SH sulfur in S-H\n");
			printf(" 2. SS sulfur in thio-ester or thio-ether or S-S\n");
			printf(" 3. sp2 sulfur in sulfone (-S=O)\n");
			printf(" 4. sp2 sulfur in sulfone (O=S=O)\n");
			printf(" 5. sp2 sulfur (S=X)\n");
			scanf("%d", &id);
		}
		while (id != 1 && id != 2 && id != 3 && id != 4 && id != 5);
		if (id == 1)
			strcpy(typename3, "sh");
		if (id == 2)
			strcpy(typename3, "ss");
		if (id == 3)
			strcpy(typename3, "s4");
		if (id == 4)
			strcpy(typename3, "s6");
		if (id == 5)
			strcpy(typename3, "s2");
	}
	if (id3 == 10) {
		id = 0;
		do {
			printf(" 1. P3 phosphate with three connected atoms\n");
			printf(" 2. P4 phosphate with five connected atoms\n");
			scanf("%d", &id);
		}
		while (id != 1 && id != 2 && id != 3);
		if (id == 1)
			strcpy(typename3, "p3");
		if (id == 2)
			strcpy(typename3, "p4");
	}

	if (typename1[1] == '\0') {
		typename1[1] = ' ';
		typename1[2] = '\0';
	}
	if (typename2[1] == '\0') {
		typename2[1] = ' ';
		typename2[2] = '\0';
	}
	if (typename3[1] == '\0') {
		typename3[1] = ' ';
		typename3[2] = '\0';
	}
	for (i = 0; i < angleparmnum; i++) {
		index = 0;
		parmname1[0] = angleparm[i].name1[0];
		parmname1[1] = angleparm[i].name1[1];
		if (parmname1[0] == 'h')
			parmname1[1] = ' ';
/*for hydrogen, we assume hc, ha, h1, h2, h3 ... are equivalent*/
		parmname2[0] = angleparm[i].name2[0];
		parmname2[1] = angleparm[i].name2[1];
		if (parmname2[0] == 'h')
			parmname2[1] = ' ';
		parmname3[0] = angleparm[i].name3[0];
		parmname3[1] = angleparm[i].name3[1];
		if (parmname3[0] == 'h')
			parmname3[1] = ' ';

		if (parmname1[0] == typename1[0] && parmname1[1] == typename1[1] &&
			parmname2[0] == typename2[0] && parmname2[1] == typename2[1] &&
			parmname3[0] == typename3[0] && parmname3[1] == typename3[1]) {
			num1 = i;
			index = 1;
			baforce = angleparm[i].force;
			bondangle = angleparm[i].angle;
			break;
		}

		if (parmname1[0] == typename1[0] && parmname1[1] == typename1[1] &&
			parmname2[0] == typename2[0] && parmname2[1] == typename2[1] &&
			parmname3[0] == typename1[0] && parmname3[1] == typename1[1]) {
			num1 = i;
			break;
		}
	}
	if (index == 1)
		return 1;
	if (num1 == -1)
		return 0;				/*unsuccessful */
	for (i = 0; i < angleparmnum; i++) {
		index = 0;
		parmname1[0] = angleparm[i].name1[0];
		parmname1[1] = angleparm[i].name1[1];
		if (parmname1[0] == 'h')
			parmname1[1] = ' ';
		parmname2[0] = angleparm[i].name2[0];
		parmname2[1] = angleparm[i].name2[1];
		if (parmname2[0] == 'h')
			parmname2[1] = ' ';
		parmname3[0] = angleparm[i].name3[0];
		parmname3[1] = angleparm[i].name3[1];
		if (parmname3[0] == 'h')
			parmname3[1] = ' ';

		if (parmname1[0] == typename3[0] && parmname1[1] == typename3[1] &&
			parmname2[0] == typename2[0] && parmname2[1] == typename2[1] &&
			parmname3[0] == typename1[0] && parmname3[1] == typename1[1]) {
			num1 = i;
			index = 1;
			baforce = angleparm[i].force;
			bondangle = angleparm[i].angle;
			break;
		}

		if (parmname1[0] == typename3[0] && parmname1[1] == typename3[1] &&
			parmname2[0] == typename2[0] && parmname2[1] == typename2[1] &&
			parmname3[0] == typename3[0] && parmname3[1] == typename3[1]) {
			num2 = i;
			break;
		}
	}

	if (index == 1)
		return 1;
	if (num2 == -1)
		return 0;				/*unsuccessful */
	bondangle = 0.5 * (angleparm[num1].angle + angleparm[num2].angle);
	return 1;					/*successful */
}

void bondanglefc(void)
{
	double d;
	if (baforce != 0.0)
		return;
	d = (bondlength1 - bondlength2) * (bondlength1 - bondlength2);
	d /= (bondlength1 + bondlength2) * (bondlength1 + bondlength2);
	baforce =
		143.9 * anglez[id1] * anglec[id2] * anglez[id3] * exp(-2 * d) /
		(bondlength1 + bondlength2);
	baforce /= sqrt(bondangle * 3.1415926 / 180.0);
}

void assignparm(void)
{
	int i, j;
	/* H<->1, C<->2, N<->3, O<->4, F<->5,P<->6, S<->7, 
	   Cl<->8, Br<->9, I<->10 */
	for (i = 0; i < 11; i++)
		for (j = 0; j < 11; j++) {
			refbondlength[i][j] = 0.0;
			bfkai[i][j] = 0.0;
		}
	refbondlength[1][1] = 0.7383;
	bfkai[1][1] = 4.661;
	refbondlength[1][2] = 1.090;
	bfkai[1][2] = 6.217;
	refbondlength[1][3] = 1.010;
	bfkai[1][3] = 6.057;
	refbondlength[1][4] = 0.960;
	bfkai[1][4] = 5.794;
	refbondlength[1][5] = 0.920;
	bfkai[1][5] = 5.600;
	refbondlength[1][6] = 1.280;
	bfkai[1][6] = 6.937;
	refbondlength[1][7] = 1.410;
	bfkai[1][7] = 7.301;
	refbondlength[1][8] = 1.600;
	bfkai[1][8] = 7.802;
	refbondlength[1][9] = 1.41;
	bfkai[1][9] = 7.257;
	refbondlength[1][10] = 1.34;
	bfkai[1][10] = 7.018;

	refbondlength[2][2] = 1.526;
	bfkai[2][2] = 7.6425;
	refbondlength[2][3] = 1.470;
	bfkai[2][3] = 7.5039;
	refbondlength[2][4] = 1.440;
	bfkai[2][4] = 7.3465;
	refbondlength[2][5] = 1.370;
	bfkai[2][5] = 7.227;
	refbondlength[2][6] = 1.800;
	bfkai[2][6] = 8.241;
	refbondlength[2][7] = 1.940;
	bfkai[2][7] = 8.478;
	refbondlength[2][8] = 2.160;
	bfkai[2][8] = 8.859;
	refbondlength[2][9] = 1.83;
	bfkai[2][9] = 8.237;
	refbondlength[2][10] = 1.82;
	bfkai[2][10] = 8.117;

	refbondlength[3][3] = 1.4406;
	bfkai[3][3] = 7.634;
	refbondlength[3][4] = 1.42;
	bfkai[3][4] = 7.526;
	refbondlength[3][5] = 1.420;
	bfkai[3][5] = 7.475;
	refbondlength[3][6] = 1.750;
	bfkai[3][6] = 8.266;
	refbondlength[3][7] = 1.930;
	bfkai[3][7] = 8.593;
	refbondlength[3][8] = 2.120;
	bfkai[3][8] = 8.963;
	refbondlength[3][9] = 1.72;
	bfkai[3][9] = 8.212;
	refbondlength[3][10] = 1.69;
	bfkai[3][10] = 8.073;

	refbondlength[4][4] = 1.46;
	bfkai[4][4] = 7.561;
	refbondlength[4][5] = 1.410;
	bfkai[4][5] = 7.375;
	refbondlength[4][6] = 1.700;
	bfkai[4][6] = 8.097;
	refbondlength[4][7] = 1.790;
	bfkai[4][7] = 8.276;
	refbondlength[4][8] = 2.110;
	bfkai[4][8] = 8.854;
	refbondlength[4][9] = 1.64;
	bfkai[4][9] = 7.957;
	refbondlength[4][10] = 1.65;
	bfkai[4][10] = 7.922;

	refbondlength[5][5] = 1.4057;
	bfkai[5][5] = 7.358;
	refbondlength[5][9] = 1.500;
	bfkai[5][9] = 7.592;
	refbondlength[5][10] = 1.580;
	bfkai[5][10] = 7.733;

	refbondlength[6][6] = 2.0308;
	bfkai[6][6] = 8.648;
	refbondlength[6][9] = 2.040;
	bfkai[6][9] = 8.656;
	refbondlength[6][10] = 2.030;
	bfkai[6][10] = 8.619;

	refbondlength[7][7] = 2.3365;
	bfkai[7][7] = 9.012;
	refbondlength[7][9] = 2.24;
	bfkai[7][9] = 8.729;
	refbondlength[7][10] = 2.21;
	bfkai[7][10] = 8.728;

	refbondlength[8][8] = 2.8357;
	bfkai[8][8] = 9.511;
	refbondlength[8][9] = 2.49;
	bfkai[8][9] = 9.058;
	refbondlength[8][10] = 2.56;
	bfkai[8][10] = 9.161;

	refbondlength[9][9] = 2.3239;
	bfkai[9][9] = 8.805;
	refbondlength[9][10] = 2.12;
	bfkai[9][10] = 8.465;

	refbondlength[10][10] = 2.038;
	bfkai[10][10] = 8.316;

        for(i=0;i<=10;i++)
                for(j=0;j<=10;j++) {
                        if(bfkai[i][j] == 0.0) bfkai[i][j] = bfkai[j][i];
                        if(refbondlength[i][j] == 0.0) refbondlength[i][j] = refbondlength[j][i];
                }

	radius[1] = 0.32;
	elecneg[1] = 2.1;

	radius[2] = 0.77;
	elecneg[2] = 2.5;

	radius[3] = 0.75;
	elecneg[3] = 3.0;

	radius[4] = 0.73;
	elecneg[4] = 3.5;

	radius[5] = 0.72;
	elecneg[5] = 4.0;

	radius[6] = 0.99;
	elecneg[6] = 3.0;

	radius[7] = 1.14;
	elecneg[7] = 2.8;

	radius[8] = 1.40;
	elecneg[8] = 2.5;

	radius[9] = 1.02;
	elecneg[9] = 2.5;

	radius[10] = 1.06;
	elecneg[10] = 2.1;

	anglec[1] = 0.0;
	anglez[1] = 0.784;

	anglec[2] = 1.339;
	anglez[2] = 1.183;

	anglec[3] = 1.3;
	anglez[3] = 1.212;

	anglec[4] = 1.249;
	anglez[4] = 1.219;

	anglec[5] = 0.0;
	anglez[5] = 1.166;

	anglec[6] = 0.0;
	anglez[6] = 1.272;

	anglec[7] = 0.0;
	anglez[7] = 1.378;

	anglec[8] = 0.0;
	anglez[8] = 1.398;

	anglec[9] = 0.906;
	anglez[9] = 1.620;

	anglec[10] = 1.448;
	anglez[10] = 1.280;
}

void assignbl(char *name1, char *name2)
{
	if (tolower(name1[0]) == 'h')
		id1 = 1;
	if (tolower(name1[0]) == 'c' && tolower(name1[1]) == 'l')
		id1 = 6;
	if (tolower(name1[0]) == 'c' && tolower(name1[1]) != 'l')
		id1 = 2;
	if (tolower(name1[0]) == 'n')
		id1 = 3;
	if (tolower(name1[0]) == 'o')
		id1 = 4;
	if (tolower(name1[0]) == 'f')
		id1 = 5;
	if (tolower(name1[0]) == 'b')
		id1 = 7;
	if (tolower(name1[0]) == 'i')
		id1 = 8;
	if (tolower(name1[0]) == 'p')
		id1 = 9;
	if (tolower(name1[0]) == 's')
		id1 = 10;
	if (tolower(name2[0]) == 'h')
		id2 = 1;
	if (tolower(name2[0]) == 'c' && tolower(name2[1]) == 'l')
		id2 = 6;
	if (tolower(name2[0]) == 'c' && tolower(name2[1]) != 'l')
		id2 = 2;
	if (tolower(name2[0]) == 'n')
		id2 = 3;
	if (tolower(name2[0]) == 'o')
		id2 = 4;
	if (tolower(name2[0]) == 'f')
		id2 = 5;
	if (tolower(name2[0]) == 'b')
		id2 = 7;
	if (tolower(name2[0]) == 'i')
		id2 = 8;
	if (tolower(name2[0]) == 'p')
		id2 = 9;
	if (tolower(name2[0]) == 's')
		id2 = 10;
}

void assignba(char *name1, char *name2, char *name3)
{
	if (tolower(name1[0]) == 'h')
		id1 = 1;
	if (tolower(name1[0]) == 'c' && tolower(name1[1]) == 'l')
		id1 = 6;
	if (tolower(name1[0]) == 'c' && tolower(name1[1]) != 'l')
		id1 = 2;
	if (tolower(name1[0]) == 'n')
		id1 = 3;
	if (tolower(name1[0]) == 'o')
		id1 = 4;
	if (tolower(name1[0]) == 'f')
		id1 = 5;
	if (tolower(name1[0]) == 'b')
		id1 = 7;
	if (tolower(name1[0]) == 'i')
		id1 = 8;
	if (tolower(name1[0]) == 'p')
		id1 = 9;
	if (tolower(name1[0]) == 's')
		id1 = 10;
	if (tolower(name2[0]) == 'h')
		id2 = 1;
	if (tolower(name2[0]) == 'c' && tolower(name2[1]) == 'l')
		id2 = 6;
	if (tolower(name2[0]) == 'c' && tolower(name2[1]) != 'l')
		id2 = 2;
	if (tolower(name2[0]) == 'n')
		id2 = 3;
	if (tolower(name2[0]) == 'o')
		id2 = 4;
	if (tolower(name2[0]) == 'f')
		id2 = 5;
	if (tolower(name2[0]) == 'b')
		id2 = 7;
	if (tolower(name2[0]) == 'i')
		id2 = 8;
	if (tolower(name2[0]) == 'p')
		id2 = 9;
	if (tolower(name2[0]) == 's')
		id2 = 10;
	if (tolower(name3[0]) == 'h')
		id3 = 1;
	if (tolower(name3[0]) == 'c' && tolower(name3[1]) == 'l')
		id3 = 6;
	if (tolower(name3[0]) == 'c' && tolower(name3[1]) != 'l')
		id3 = 2;
	if (tolower(name3[0]) == 'n')
		id3 = 3;
	if (tolower(name3[0]) == 'o')
		id3 = 4;
	if (tolower(name3[0]) == 'f')
		id3 = 5;
	if (tolower(name3[0]) == 'b')
		id3 = 7;
	if (tolower(name3[0]) == 'i')
		id3 = 8;
	if (tolower(name3[0]) == 'p')
		id3 = 9;
	if (tolower(name3[0]) == 's')
		id3 = 10;
}



int main()
{
	int cmdid = 1;

    amberhome = (char *) getenv("AMBERHOME");
    if( amberhome == NULL ){
       fprintf( stdout, "AMBERHOME is not set!\n" );
       exit(1);
    }
	strcpy(filename, amberhome);
	strcat(filename, "/dat/leap/parm/gaff.dat");
	assignparm();
	readangleparm(filename);	/*principle parameter file */
	cmdid = 1;
	while (cmdid == 1 || cmdid == 2) {
		blforce = 0.0;
		bondlength = 0.0;
		bondlength1 = 0.0;
		bondlength2 = 0.0;
		bondangle = 0.0;
		baforce = 0.0;
		printf("Please select:\n");
		printf("1. calculate the bond length parameter: A-B\n");
		printf("2. calculate the bond angle parameter: A-B-C\n");
		printf("3. exit\n");
		scanf("%d", &cmdid);
		if (cmdid == 3)
			exit(1);
		if (cmdid == 1) {
			printf("Please input the element name of atom A in A-B\n");
			scanf("%s", name1);
			printf("Please input the element name of atom B in A-B\n");
			scanf("%s", name2);
			assignbl(name1, name2);
			printf
				("Please input the bond length (a non-positive number\n");
			printf("means to calculate it according to empirical rules\n");
			scanf("%lf", &bondlength);
			if (bondlength <= 0)
				bondlengthcal();
			bondlengthfc();
			printf("\n\nBOND %2s-%2s%9.3lf%12.3lf\n\n\n\n", name1, name2,
				   blforce, bondlength);
		}
		if (cmdid == 2) {
			printf("Please input the element name of atom A in A-B-C\n");
			scanf("%s", name1);
			printf("Please input the element name of atom B in A-B-C\n");
			scanf("%s", name2);
			printf("Please input the element name of atom C in A-B-C\n");
			scanf("%s", name3);
			assignba(name1, name2, name3);
			printf("Please input the bond length between A-B\n");
			scanf("%lf", &bondlength1);
			printf("Please input the bond length between B-C\n");
			scanf("%lf", &bondlength2);
			printf("Please input the bond angle (a non-positive number\n");
			printf("means to calculate it according to empirical rules\n");
			scanf("%lf", &bondangle);
			if (bondangle <= 0)
				bondanglecal();
			bondanglefc();
			printf("\n\nANGLE %2s-%2s-%2s%9.3lf%12.3lf\n\n\n", name1,
				   name2, name3, baforce, bondangle);
		}
	}
	return (0);
}
