#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include <stdarg.h>

enum _parseOperators {
	PARSE_NOOP,
	PARSE_CONTINUATION,
	PARSE_ATOM,
	PARSE_RESIDUE
};

typedef struct _parseEntry {
	char *token;
	int operator;
	int isnumber;
} parseEntry;


#define NAME_SIZE 5
#define NAME_DEFAULT "    "
typedef char Name[NAME_SIZE];



typedef struct _stackType {
	void *entry;
	struct _stackType *next;
} stackType;

typedef void (*fxnPrintStackEntry) (void *);

#ifndef skipWhitespace
#define skipWhitespace( xxx ) \
  while( ((xxx[0] == '\n') || isspace(xxx[0])) && strlen(xxx) ) xxx++;
#endif


void *safe_malloc(int size)
{
	void *mem;
	if (size) {
		if ((mem = (void *) malloc((size_t) size)) == NULL)
			fprintf(stderr, "safe_malloc: Error in alloc of %i bytes",
					size);
	} else
		mem = NULL;
	return (mem);
}


parseEntry *parseToken(char **textp, int operator)
{
	char *text;
	int i, j;
	parseEntry *p;

	text = *textp;

	for (i = 0;; i++) {
		if (text[i] == (char) 0 ||
			text[i] == ',' ||
			text[i] == ':' ||
			text[i] == '@' ||
			(i > 0 && text[i] == '-' && !isalpha(text[i - 1])) ||
			isspace(text[i]))
			break;
	}

	p = (parseEntry *) safe_malloc(sizeof(parseEntry));
	p->isnumber = 0;

	if (i > 0)
		p->token = (char *) safe_malloc(sizeof(char) * (i + 1));

	for (j = 0; j < i; j++) {
		p->token[j] = text[j];
	}
	p->token[i] = (char) 0;


	p->operator = PARSE_NOOP;
	switch (text[i]) {

	case '-':
		p->operator = PARSE_CONTINUATION;

	case ',':

		text++;

	}

	if (isdigit(p->token[0])) {
		p->isnumber = 1;
		for (j = 1; j < i; j++) {
			if (isdigit(p->token[j]) == 0) {
				p->isnumber = 0;
			}
		}
	}

	text = text + i;
	skipWhitespace(text);
	/*
	 *  extra special check to handle extra spacing between continuation
	 *  operators
	 */
	if (text[0] == '-') {
		p->operator = PARSE_CONTINUATION;
		text += 1;
		skipWhitespace(text);
	} else if (text[0] == ',') {
		text += 1;
		skipWhitespace(text);
	}

	*textp = text;
	return (p);


}


int isMatch(char *s1, char *s2)
{
	int i;

	/*
	 *  straight match
	 */

	if (strcmp(s1, s2) == 0)
		return 1;

	/*
	 *  fuzzy match: NOTE this will break if s1 > s2
	 */

	if (strchr(s1, '?') != NULL || strchr(s1, '*') != NULL) {

		/*
		 *  look over the minimal map between the two strings
		 */
		for (i = 0; i < strlen(s1) && i < strlen(s2); i++) {

			if (s1[i] != s2[i]) {
				switch (s1[i]) {

				case '*':		/* wild card, multiple characters */
					return 1;
				case '?':		/* wild card, single character    */
					if (isspace(s2[i]))
						return 0;
					break;
				default:		/* mismatch                       */
					return 0;
				}
			}
		}
		return 1;
	}
	return 0;
}


void safe_free(void *pointer)
{
	if (pointer != NULL)
		free(pointer);
}


void pushStack(stackType ** stackp, void *entry)
{
	stackType *sp;

	sp = safe_malloc(sizeof(stackType));
	sp->entry = entry;
	sp->next = *stackp;
	*stackp = sp;
}


void pushBottomStack(stackType ** stackp, void *entry)
{
	stackType *sp;

	if (*stackp == NULL)
		pushStack(stackp, entry);
	else {
		sp = *stackp;
		while (sp->next != NULL)
			sp = sp->next;
		sp->next = safe_malloc(sizeof(stackType));
		sp->next->entry = entry;
		sp->next->next = NULL;
	}
}


void *popStack(stackType ** stackp)
{
	void *entry;
	stackType *sp;

	if (*stackp == NULL) {
		return ((char *) NULL);
	}

	sp = *stackp;
	entry = sp->entry;

	*stackp = sp->next;
	sp->next = NULL;
	free(sp);
	return (entry);
}

void clearStack(stackType ** stackp)
{
	stackType *sp, *tmp;

	if (stackp == NULL) {
		return;
	}

	tmp = NULL;
	for (sp = *stackp; sp != NULL;) {
		if (tmp != NULL)
			safe_free((void *) tmp);
		safe_free((void *) sp->entry);
		sp->entry = NULL;
		tmp = sp;
		sp = sp->next;
	}
	*stackp = NULL;
}


void printStack(stackType ** stackp, fxnPrintStackEntry fxn, char *babble)
{
	stackType *p;
	int i;
	for (p = *stackp, i = 1; p != NULL; p = p->next, i++) {
		if (babble != NULL)
			fprintf(stdout, "%s (%i)\n", babble, i);
		fxn(p->entry);
	}

}




int *processAtomMaskDetailed(char *maskString, int atoms, int residues,
							 Name * atomName, Name * residueName,
							 int *ipres)
{
	int actualAtoms;
	int not = 0;
	char *maskp;
	char *tmp;

	int *residueMask;
	int residueMaskActive = 0;
	int *atomMask;
	int atomMaskActive = 0;
	int res1, res2;
	int atom1, atom2;
	int continuation = 0;
	int i, j;
	stackType *residueStack = NULL;
	stackType *atomStack = NULL;
	parseEntry *pp;

	Name name;

	maskp = maskString;
	skipWhitespace(maskp);

	/*
	 *  allocate mask strings
	 */
	atomMask = (int *) safe_malloc(sizeof(int) * atoms);
	residueMask = (int *) safe_malloc(sizeof(int) * residues);
	memset(atomMask, 0, sizeof(int) * atoms);
	memset(residueMask, 0, sizeof(int) * residues);

	/*
	 *  special case, choose ALL atoms
	 */
	if (maskp[0] == (char) 0 || maskp[0] == '*') {
		for (i = 0; i < atoms; i++)
			atomMask[i] = 1;
		goto clean_return;
	}


	/*
	 *  get rid of all NOT characters "~"; only one is significant
	 *  and set NOT status if one is found...
	 */
	while ((tmp = strchr(maskString, '~')) != NULL) {
		not = 1;
		tmp[0] = ' ';
	}

	/*
	 *  check for error
	 */
	if (strchr(maskp, ':') == NULL && strchr(maskp, '@') == NULL) {
		fprintf(stdout,
				"WARNING: Error in mask string, no \"@\" or \":\" present (%s)\n",
				maskString);
		safe_free(atomMask);
		safe_free(residueMask);
		return NULL;
	}

	/*
	 *  do the main "parsing"
	 */
	skipWhitespace(maskp);
	while (maskp[0] != (char) 0) {

		if (maskp[0] == ':') {
			maskp++;
			for (;;) {

				skipWhitespace(maskp);
				pp = parseToken(&maskp, PARSE_RESIDUE);
				pushBottomStack(&residueStack, (void *) pp);
				if (maskp[0] == (char) 0 ||
					maskp[0] == '@' || maskp[0] == ':')
					break;
			}
		}

		if (maskp[0] == '@') {
			maskp++;
			for (;;) {

				skipWhitespace(maskp);
				pp = parseToken(&maskp, PARSE_ATOM);
				pushBottomStack(&atomStack, (void *) pp);
				if (maskp[0] == (char) 0 || maskp[0] == ':')
					break;
				if (maskp[0] == '@')
					maskp++;
			}
		}

		/*
		 *  now process the atomStack and residueStack
		 */


		if (not) {
			for (i = 0; i < atoms; i++)
				atomMask[i] = 1;
			for (i = 0; i < residues; i++)
				residueMask[i] = 1;
		}

		while (residueStack != NULL) {

			if (continuation) {
				res1 = res2;
				res2 = -1;
			}

			pp = (parseEntry *) popStack(&residueStack);
			if (pp->isnumber) {
				if (sscanf(pp->token, "%i", &res2) != 1) {
					fprintf(stdout, "WARNING: error parsing atom mask\n");
					safe_free(atomMask);
					safe_free(residueMask);
					return NULL;
				}
				res2--;

				if (continuation) {
					continuation = 0;
					if (res1 < 0)
						res1 = 0;
					if (res2 >= residues)
						res2 = residues - 1;
					if (res2 < res1)
						res2 = res1;

					for (i = res1; i <= res2; i++) {
						residueMask[i] = (not ? 0 : 1);
						residueMaskActive = 1;
					}
				} else {
					residueMask[res2] = (not ? 0 : 1);
					residueMaskActive = 1;
				}

				if (pp->operator == PARSE_CONTINUATION)
					continuation = 1;

			} else {

				continuation = 0;
				strcpy(name, "    ");
				for (i = 0; i < strlen(pp->token) && i < NAME_SIZE; i++) {
					name[i] = pp->token[i];
				}
				name[NAME_SIZE - 1] = (char) 0;

				for (i = 0; i < residues; i++)
					if (isMatch(name, residueName[i])) {
						residueMask[i] = (not ? 0 : 1);
						residueMaskActive = 1;
					}
			}
			safe_free(pp->token);
			pp->token = NULL;
			safe_free(pp);
		}

		while (atomStack != NULL) {

			if (continuation) {
				atom1 = atom2;
				atom2 = -1;
			}

			pp = (parseEntry *) popStack(&atomStack);
			if (pp->isnumber) {
				if (sscanf(pp->token, "%i", &atom2) != 1) {
					fprintf(stdout, "WARNING: error parsing atom mask\n");
					safe_free(atomMask);
					safe_free(residueMask);
					return NULL;
				}
				atom2--;

				if (continuation) {
					continuation = 0;
					if (atom1 < 0)
						atom1 = 0;
					if (atom2 > atoms)
						atom2 = atoms;
					if (atom2 < atom1)
						atom2 = atom1;

					for (i = atom1; i < atom2; i++) {
						atomMask[i] = (not ? 0 : 1);
						atomMaskActive = 1;
					}
				} else {
					atomMask[atom2] = (not ? 0 : 1);
					atomMaskActive = 1;
				}

				if (pp->operator == PARSE_CONTINUATION)
					continuation = 1;

			} else {

				continuation = 0;
				strcpy(name, "    ");
				for (i = 0; i < strlen(pp->token) && i < NAME_SIZE; i++) {
					name[i] = pp->token[i];
				}
				name[NAME_SIZE - 1] = (char) 0;

				for (i = 0; i < atoms; i++)
					if (isMatch(name, atomName[i])) {
						atomMask[i] = (not ? 0 : 1);
						atomMaskActive = 1;
					}
			}
			safe_free(pp->token);
			pp->token = NULL;
			safe_free(pp);
		}

		if (atomMaskActive && residueMaskActive) {
			for (i = 0; i < residues; i++) {

				if (residueMask[i] == 0) {
					for (j = ipres[i] - 1; j < ipres[i + 1] - 1; j++) {
						atomMask[j] = 0;
					}
				}
			}
		} else if (residueMaskActive) {
			for (i = 0; i < residues; i++) {

				for (j = ipres[i] - 1; j < ipres[i + 1] - 1; j++) {
					if (residueMask[i])
						atomMask[j] = 1;
					else if (not)
						atomMask[j] = 0;
				}
			}
		}

		atomMaskActive = 0;
		residueMaskActive = 0;
	}

	clean_return:

	actualAtoms = 0;
	for (i = 0; i < atoms; i++)
		if (atomMask[i])
			actualAtoms++;


	if (tmp = strchr(maskString, '\n'))
		tmp[0] = (char) 0;

/*  ---no printout inside these routines any more!   */

	if (actualAtoms > 0) {
/*		fprintf(stdout, "Mask %s%s] represents %i atoms\n",
				(not ? "[~" : "["), maskString, actualAtoms);   */
	} else {
/*		fprintf(stdout, "Mask %s%s] represents %i atoms ",
				(not ? "[~" : "["), maskString, actualAtoms);
		fprintf(stdout, "!!!NO ATOMS DETECTED!!!\n");          */
		safe_free(atomMask);
		atomMask = NULL;
	}

	safe_free(residueMask);
	return (atomMask);

}

/*  Fortran-callable wrapper, for use for calling from sander  */

#ifdef CLINK_CAPS
#  define parsemask_ PARSEMASK
#endif
#ifdef CLINK_PLAIN
#  define parsemask_ parsemask
#endif

void
parsemask_(char *maskString, int *atoms, int *residues,
		   char *atomName, char *residueName, int *ipres, int *mask,
		   int *natc)
{
	int i, j, Natoms, Nresidues;
	int *atomMask;
	Name *paddedAtoms, *paddedResidues;

	Natoms = *atoms;

	paddedAtoms = (Name *) safe_malloc(sizeof(Name) * Natoms);
	j = 0;
	for (i = 0; i < Natoms; i++) {
		strcpy(paddedAtoms[i], "    \0");
		strncpy(paddedAtoms[i], atomName + j, 4);
		j += 4;
	}

	Nresidues = *residues;
	paddedResidues = (Name *) safe_malloc(sizeof(Name) * Nresidues);
	j = 0;
	for (i = 0; i < Nresidues; i++) {
		strcpy(paddedResidues[i], "    \0");
		strncpy(paddedResidues[i], residueName + j, 4);
		j += 4;
	}

	atomMask = processAtomMaskDetailed(maskString, Natoms, Nresidues,
									   paddedAtoms, paddedResidues, ipres);

	if (!atomMask) {
		*natc = 0;
		for (i = 0; i < Natoms; i++) {
			mask[i] = 0;
		}
	} else {
		*natc = 0;
		for (i = 0; i < Natoms; i++) {
			mask[i] = atomMask[i];
			*natc += mask[i];
		}
	}
	safe_free(atomMask);
	safe_free(paddedAtoms);
	safe_free(paddedResidues);

}
