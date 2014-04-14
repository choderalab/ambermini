/* SMILES */

#include "common.h"
#include "define.h"
#include "atom.h"
#include "globals.h"

void rsmiles(char *filename)
{
	char tmpchar[MAXCHAR];
	int status = 0;
	strcpy(tmpchar, "trigo concord4.02 < ");
	strcat(tmpchar, ifilename);
	if (intstatus == 2)
		fprintf(stdout, "\nRunning: %s\n", tmpchar);
	status = system(tmpchar);
	if(status != 0) {
                fprintf(stdout, "Error: cannot run \"%s\" in rsmiles() of smiles.c properly, exit\n", tmpchar);
                exit(1);
        }
	rmol2("coords.mol2");
}
void wsmiles(char *filename)
{
	printf("\n Sorry, cannot write smiles format at present");
}
