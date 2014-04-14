
/*      File:   mol3File.c
 ______________________________________________________________________________
|                                                                              |
|       Author: Mason Louchart (2011)                                          |
|       http://q4md-forcefieldtools.org                                        |
|       Universite de Picardie - Jules Verne, Amiens                           |
|                                                                              |
|       Tutorial available at                                                  |
|       http://q4md-forcefieldtools.org/Tutorial/leap-mol3.php                 |
|_____________________________________________________________________________*/

#include        <limits.h>

#include        "basics.h"
#include        "classes.h"
#include        "dictionary.h"
#include        "database.h"
#include        "pdb.h"
#include        "varArray.h"
#include        "avl.h"
#include        "sort.h"
#include        "matrix.h"
#include        "model.h"
#include        "minimizer.h"
#include        "leap.h"
#include        "defaults.h"
#include        "atom.h"
#include        "tripos.h"
#include        "stdio.h"
#include        "parmLib.h"
#include        "unitio.h"

#include	"mol2File.h"


/* Main function, write the Mol3 file */
void 
Mol3Write( FILE *fOut, UNIT uUnit, int choice ){
    int i,j,iResNumber = 0;
    VARARRAY  vaRes;
    STRING sName;
    char* cPTemp;

    Mol2Write( fOut, uUnit, choice );

    fprintf(fOut, "@<TRIPOS>HEADTAIL\n") ;

    /* Write the HEAD in the file */
    if ( aUnitHead(uUnit) != NULL ){
        iResNumber=(((CONTAINER)aUnitHead(uUnit))->cContainedBy)->iSequence;

        cPTemp = ((CONTAINER)aUnitHead(uUnit))->sName;
        strncpy( sName, cPTemp, 4 );
      
        if ( strlen(cPTemp) > 4 ) {
            VP0(( " Truncating residue name for PDB format: %s -> %s\n", ((CONTAINER)aUnitHead(uUnit))->sName, sName ));
            fprintf(fOut,"%s %d\n",sName,iResNumber);
        }
        else {
            fprintf(fOut,"%s %d\n",((CONTAINER)aUnitHead(uUnit))->sName,iResNumber);
        }
    }
    else {
        fprintf(fOut,"0 0\n");
    }

    /* Write the TAIL in the file */
    if ( aUnitTail(uUnit) != NULL ){
        iResNumber=(((CONTAINER)aUnitTail(uUnit))->cContainedBy)->iSequence;

        cPTemp = ((CONTAINER)aUnitTail(uUnit))->sName;
        strncpy( sName, cPTemp, 4 );

        if ( strlen(cPTemp) > 4 ) {
            VP0(( " Truncating residue name for PDB format: %s -> %s\n", ((CONTAINER)aUnitTail(uUnit))->sName, sName ));
            fprintf(fOut,"%s %d\n",sName,iResNumber);
        }
        else{
            fprintf(fOut,"%s %d\n",((CONTAINER)aUnitTail(uUnit))->sName,iResNumber);
        }
    }else{
        fprintf(fOut,"0 0\n");
    }

    fprintf(fOut,"@<TRIPOS>RESIDUECONNECT\n");

    /* Write Atoms Connect for every Residue */
    vaRes=uUnit->vaResidues;
    for (i = 0; i < iVarArrayElementCount((VARARRAY)vaRes); i++) {
        SAVERESIDUEt prRes = *PVAI(vaRes, SAVERESIDUEt, i);
        fprintf(fOut,"%d",prRes.iSequenceNumber);
        for (j = 0; j <= 5; j++) {
            if ( prRes.rResidue->aaConnect[j] != NULL ){
                cPTemp = sContainerName(prRes.rResidue->aaConnect[j]);
                strncpy( sName, cPTemp, 4 );
                if ( strlen(cPTemp) > 4 ) {
                    VP0(( " Truncating residue name for PDB format: %s -> %s\n", sContainerName(prRes.rResidue->aaConnect[j]), sName ));
                    fprintf(fOut," %s",sName);
                }else{
                    fprintf(fOut," %s",sContainerName(prRes.rResidue->aaConnect[j]) );
                }
            }else{
                fprintf(fOut," 0");
            }
        }
        fprintf(fOut,"\n");
    }
}
