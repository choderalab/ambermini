/*
 *      File:   amber.c
 *
 ************************************************************************
 *                            LEAP                                      *
 *                                                                      *
 *                   Copyright (c) 1992, 1995                           *
 *           Regents of the University of California                    *
 *                     All Rights Reserved.                             *
 *                                                                      *
 *  This software provided pursuant to a license agreement containing   *
 *  restrictions on its disclosure, duplication, and use. This software *
 *  contains confidential and proprietary information, and may not be   *
 *  extracted or distributed, in whole or in part, for any purpose      *
 *  whatsoever, without the express written permission of the authors.  *
 *  This notice, and the associated author list, must be attached to    *
 *  all copies, or extracts, of this software. Any additional           *
 *  restrictions set forth in the license agreement also apply to this  *
 *  software.                                                           *
 ************************************************************************
 *                                                                      *
 *     Designed by:    Christian Schafmeister                           *
 *     Author:         Christian Schafmeister                           *
 *                                                                      *
 *     VERSION: 1.0                                                     *
 *     Programmers:                                                     *
 *             Christian Schafmeister                                   *
 *             David Rivkin                                             *
 *                                                                      *
 *     Principal Investigator: Peter A. Kollman                         *
 *                                                                      *
 ************************************************************************
 *
 *      Description:
 *              Read AMBER prep and parameter database files into
 *              UNITs and PARMSETs.
 */

#include        "basics.h"
#include        "vector.h"
#include        "classes.h"
#include        "leap.h"
#include        "tools.h"
#include        "amber.h"
#include        "defaults.h"
#include        "sort.h"
#include        "zMatrix.h"
#include        "cmap.h"

#define FGETS(s,f)      { strcpy(s,""); fgets(s,sizeof(s),f); }

/*
 *      Mapping of AMBER types to ELEMENT names used by iAtomTypesToElement
 */

typedef struct  {
        char    sType[5];
        int     iElement;
        int     iHybridization;
} TYPEELEMENTt;

static VARARRAY         vaAtomTypes = NULL;

typedef struct  {
        char            sNonBondType[NAMELEN];
        double          dRStar;
        double          dDepth;
        double          dRStar14;
        double          dDepth14;
} NONBONDt;

typedef struct  {
        char            sType[NAMELEN];
        double          dMass;
        double          dPolar;
        double          dScreenF;
        int             iElement;
        int             iHybridization;
} MASSt;        

#define MAXEQUIV        20
typedef char    SHORTt[NAMELEN];
typedef struct  {
        char            sName[NAMELEN];
        int             iEquivs;
        SHORTt          saEquivs[MAXEQUIV];
} EQUIVt;

CMAP *cmap=NULL;
CMAPLST *cmaplst=NULL;
int mapnum=0;
/*
 *----------------------------------------------------
 *
 *      Private routines
 */

#define NODASHES(s)     {int z;for(z=0;z<9;z++) if(s[z]=='-')s[z]=' ';}
#define NODASHESWL(s)     {int z;for(z=0;z<strlen(s);z++) if(s[z]=='-')s[z]=' ';}


/*
 *      zAmberConvertWildCard
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      If the atom type is an "X", an old AMBER wild card
 *      character then convert it to the new wild card character.
 */
static void
zAmberConvertWildCard( char *sType )
{
    if ( strcmp( sType, AMBER_WILD_CARD ) == 0 ) {
        strcpy( sType, WILD_CARD_TYPE );
    }
}



/*
 *      zAmberReadMasses
 *
 *      Read the atomic masses into a VARARRAY
 */
static void     
zAmberReadParmSetMasses( VARARRAY *vaPMasses, FILE *fIn )
{
STRING          sLine;
int             iRead;
VARARRAY        vaMasses;
MASSt           mMass;
double          dPolar;
double          dScreenF;
int             iElement;
STRING          sHybridization;

                /* Read the atom masses into a VARARRAY until we */
                /* read the non-bond parameters */

    memset(&mMass, 0, sizeof(mMass));                   /* for Purify */
    MESSAGE(( "Reading masses\n" ));
    vaMasses = vaVarArrayCreate( sizeof(MASSt) );

    while (1) {
        FGETS( sLine, fIn );
        MESSAGE(( "Read: %s\n", sLine ));
        iRead = sscanf( sLine, "%s %lf %lf %lf %d %s", 
                        mMass.sType, &mMass.dMass, &dPolar, &dScreenF,
                        &iElement, sHybridization );
        if ( iRead <= 0 ) 
                break;
        if ( iRead == 1 ) {
                VP0(( "ERROR: %s: mass not read - omitting\n", mMass.sType ));
                continue;
        }
        if ( mMass.dMass < 0.0 ) {
                VP0(( "ERROR: %s: mass %lf - omitting\n", 
                                                mMass.sType, mMass.dMass ));
                continue;
        }
        mMass.dPolar = -1.0;    /* flag that default=0.0 will be used */
        mMass.dScreenF = 0.0;
        mMass.iElement=-10;   /* NOELEMENT==EP==-1 ; -10 seems safe */
        mMass.iHybridization=HUNDEFINED;
        switch( iRead ){
            default:
            case 6:
                StringLower( sHybridization );
                if ( !strcmp( sHybridization, "sp2") )
                       mMass.iHybridization = HSP2;
                else if ( !strcmp( sHybridization, "sp3") )
                        mMass.iHybridization = HSP3;
                else if ( !strcmp( sHybridization, "sp") )
                        mMass.iHybridization = HSP1;
                else if ( !strcmp( sHybridization, "sp0") )
                        mMass.iHybridization = 0;
                else {
                        mMass.iHybridization = HUNDEFINED;
                        VP0(( "atom type %s - unknown hybridization %s\n",
                                   mMass.sType, sHybridization ));
                }
            case 5:
                mMass.iElement = iElement;
            case 4:
                mMass.dScreenF = dScreenF;
            case 3:
                if ( dPolar < 0.0 ) {
                    if ( dPolar < -0.99 && dPolar > -1.01 )
                        break;
                    VP0(( "ERROR: %s: polarization %lf - omitting\n",
                                                mMass.sType, dPolar ));
                    continue;
                }
                if ( dPolar > 15.0 )
                    VP0(( "WARNING: %s: polarization %lf\n",
                                                mMass.sType, dPolar ));
                mMass.dPolar = dPolar;
            case 2:
                break;
        };

        VarArrayAdd( vaMasses, (GENP)&mMass );
    }
    *vaPMasses = vaMasses;
}





/*
 *      zAmberReadParmSetBonds
 *
 *      Read the bond parameter terms.
 */
static void
zAmberReadParmSetBonds( PARMSET psParms, FILE *fIn )
{
STRING          sLine;
int             iRead;
STRING          saStr[10];
double          dKb, dR0;

    memset(saStr, 0, sizeof(saStr));                    /* for Purify */
    while (1) {
        FGETS( sLine, fIn );
        NODASHES(sLine);
        iRead = sscanf( sLine, "%s %s %lf %lf", saStr[0], saStr[1], 
                                &dKb, &dR0 );
        if ( iRead<=0 ) break;
        MESSAGE(( "Read: %s\n", sLine ));
        iParmSetAddBond( psParms, saStr[0], saStr[1], dKb, dR0, "" );
    } 
       
} 



/*
 *      zAmberReadParmSetAngles
 *
 *      Read the angle parameter terms.
 */
static void
zAmberReadParmSetAngles( PARMSET psParms, FILE *fIn )
{
STRING          sLine;
int             iRead;
STRING          saStr[10];
double          dKt, dT0, dTkub, dRkub, zero;

    zero = 0.0;
    memset(saStr, 0, sizeof(saStr));                /* for Purify */
    while (1) {
        FGETS( sLine, fIn );
        NODASHES(sLine);
        if( GDefaults.iCharmm ){
            iRead = sscanf( sLine, "%s %s %s %lf %lf %lf %lf", 
                            saStr[0], saStr[1], saStr[2], &dKt, &dT0, &dTkub, &dRkub );
            if ( iRead != 7 ) break;
            MESSAGE(( "Read: %s\n", sLine ));
            iParmSetAddAngle( psParms, saStr[0], saStr[1], saStr[2],
                              dKt, dT0*DEGTORAD, dTkub, dRkub, "" );
        } else {
            iRead = sscanf( sLine, "%s %s %s %lf %lf", 
                            saStr[0], saStr[1], saStr[2], &dKt, &dT0 );
            if ( iRead != 5 ) break;
            MESSAGE(( "Read: %s\n", sLine ));
            iParmSetAddAngle( psParms, saStr[0], saStr[1], saStr[2],
                              dKt, dT0*DEGTORAD, zero, zero, "" );
        }
    } 
}

static void
zAmberReadParmSetIPOL( FILE *fIn )
{
STRING          sLine;
int             iRead,myiPOL;
    while (1) {
        FGETS( sLine, fIn );
        NODASHES(sLine);
        iRead = sscanf( sLine, "%i", &myiPOL );
        if ( iRead != 1 ) break;
        MESSAGE(( "Read: %s\n", sLine ));
        if (GDefaults.iIPOLset == 1 && GDefaults.iIPOL != myiPOL ) {
           VP0(( "Error: Conflicting IPOL. %i %i\n", myiPOL, GDefaults.iIPOL));
        } else {
           GDefaults.iIPOL = myiPOL;
           GDefaults.iIPOLset = 1;
        }
    } 
}
/*
 *      zAmberReadParmSetCMAP
 *
 */
static void
zAmberReadParmSetCMAP( VARARRAY *vaFoo, FILE *fIn )
{
    int    iRead, flag;
    STRING sLine,sLine0;
    STRING saStr[10];
    STRING tmpchar1, tmpchar2;
    int imap;
    int tmpint1;
    int i;
    CMNT *cmntt, *cmnt, *cmnt0;
    CMAPLST *cmaplstt;
    CMAP *cmapt;
    static int called = 0;

    if (!called) {
// initialize cmap record
          cmaplst = (CMAPLST *) malloc(sizeof( CMAPLST ));
          cmaplst->cmap = NULL;
          cmaplst->next = NULL;
          cmap = NULL;
          mapnum = 0;
          called ++;
      }

    cmaplstt = cmaplst;
    while (cmaplstt->next != NULL) {
         cmaplstt=cmaplstt->next;
      }

    cmnt=NULL;
    cmnt0=NULL;
    cmntt=NULL;
    flag = 0;
    imap = 0;
//    mapnum = 0;
    cmnt0=(CMNT *)malloc(sizeof(CMNT));
    cmntt=cmnt0;
    cmntt->next=NULL;
    
    //This is a pseudo reader
    while (1) {
        
        FGETS(sLine0, fIn);
        strcpy(sLine,sLine0);

        NODASHESWL(sLine);
        if ( sLine[0] == '%' ) {
            sscanf(sLine, "%s%s", tmpchar1, tmpchar2);
            if (strcmp("%FLAG", tmpchar1) == 0) {
                flag=0;
                if (strcmp("CMAP_COUNT", tmpchar2)==0){
                    sscanf(sLine, "%s%s%d", tmpchar1, tmpchar2, &tmpint1);

//                    if (mapnum <= 1) mapnum = 1;
                    // initialize storage for maps
//                    cmap = (CMAP *)malloc(sizeof(CMAP) * (mapnum));
                    //printf("mjhsieh: mapnum =%d\n",mapnum);
//                    mapnum += tmpint1;
//                    VP0(("Read %i cmaps\n",tmpint1));

                }else if(strcmp("CMAP_TITLE", tmpchar2)==0){
                    FGETS(sLine, fIn);
                    cmap = (CMAP *)malloc(sizeof(CMAP) );
                    while (cmaplstt->cmap != NULL ) {
                          cmaplstt = cmaplstt->next;
                       }
                    cmaplstt->cmap = cmap;
                    cmaplstt->next = (CMAPLST *) malloc(sizeof(CMAPLST));
                    cmaplstt->next->cmap = NULL;
                    cmaplstt->next->next = NULL;

                    for (i=0;i<80 && sLine[i] != '\n'; i++) cmap->title[i]=sLine[i];
                    cmap->title[i] = '\0';
//                    VP0(("Read cmap %s\n",cmap->title));
                    cmnt=(CMNT *) malloc(sizeof(CMNT));
                    cmnt->next=NULL;
                    cmap->cmnt=cmnt;
                    mapnum ++;

// set default to the main chain
                    cmap->residx[0]=-1;
                    cmap->residx[1]= 0;
                    cmap->residx[2]= 0;
                    cmap->residx[3]= 0;
                    cmap->residx[4]= 1;

                    cmap->nresidx[0]= 0;
                    cmap->nresidx[1]= 0;
                    cmap->nresidx[2]= 0;
                    cmap->nresidx[3]= 0;
                    cmap->nresidx[4]= 1;

                    cmap->cresidx[0]=-1;
                    cmap->cresidx[1]= 0;
                    cmap->cresidx[2]= 0;
                    cmap->cresidx[3]= 0;
                    cmap->cresidx[4]= 0;

                    strcpy(cmap->atmname[0],"C");
                    strcpy(cmap->atmname[1],"N");
                    strcpy(cmap->atmname[2],"CA");
                    strcpy(cmap->atmname[3],"C");
                    strcpy(cmap->atmname[4],"N");

                    strcpy(cmap->natmname[0],"H1");
                    strcpy(cmap->natmname[1],"N");
                    strcpy(cmap->natmname[2],"CA");
                    strcpy(cmap->natmname[3],"C");
                    strcpy(cmap->natmname[4],"N");

                    strcpy(cmap->catmname[0],"C");
                    strcpy(cmap->catmname[1],"N");
                    strcpy(cmap->catmname[2],"CA");
                    strcpy(cmap->catmname[3],"C");
                    strcpy(cmap->catmname[4],"OXT");

                    cmap->termmap = 1; // assume applicable to terminal residues by default

                } else if(strcmp("CMAP_RESLIST",tmpchar2)==0){
                    char tmpchar[8];
                    sscanf(sLine,"%s%s%d",tmpchar1,tmpchar2,&tmpint1);
                    cmap->nres = tmpint1;
                    cmap->reslist=(WRD *) malloc(sizeof(WRD)*tmpint1 );
                    cmap->creslist=(WRD *) malloc(sizeof(WRD)*tmpint1 );
                    cmap->nreslist=(WRD *) malloc(sizeof(WRD)*tmpint1 );
                    for (i=0; i<tmpint1; i++) {
                        fscanf(fIn, "%s", tmpchar);
                        //printf("%s\n",tmpchar);
                        strcpy(cmap->reslist[i], tmpchar);
                        strcpy(cmap->creslist[i], "C"); strcat(cmap->creslist[i], tmpchar);
                        strcpy(cmap->nreslist[i], "N"); strcat(cmap->nreslist[i], tmpchar);
                    }
                } else if (strcmp("CMAP_RESOLUTION", tmpchar2) == 0) {
                    sscanf(sLine, "%s%s%d", tmpchar1, tmpchar2, &tmpint1);
                    cmap->resolution = tmpint1;
                    cmap->map=(double *) malloc(sizeof(double) * tmpint1 * tmpint1 );
                    //printf("%s\n",sLine);
                } else if (strcmp("CMAP_PARAMETER", tmpchar2) == 0) {
                    double tmpdbl;
                    tmpint1 = cmap->resolution * cmap->resolution;
                    for (i=0; i<tmpint1; i++) {
                        fscanf(fIn, "%lf", &tmpdbl);
                        cmap->map[i]=tmpdbl;
                    }
                    //printf("%s\n",sLine);
                } else if(strcmp("CMAP_ATMLIST",tmpchar2)==0){
                    sscanf(sLine,"%s%s%s%s%s%s%s",tmpchar1,tmpchar2,
                           cmap->atmname[0], cmap->atmname[1], cmap->atmname[2],
                           cmap->atmname[3], cmap->atmname[4]);
                           int l;
                           for (l=0; l<5; l++) {
                                strcpy(cmap->catmname[l], cmap->atmname[l]);
                                strcpy(cmap->natmname[l], cmap->atmname[l]);
                              }
                           if (strcmp("C",cmap->atmname[0]) == 0) strcpy(cmap->natmname[0], "H1");
                           if (strcmp("C",cmap->atmname[4]) == 0) strcpy(cmap->natmname[4], "H1");
                           if (strcmp("N",cmap->atmname[0]) == 0) strcpy(cmap->catmname[0], "OXT");
                           if (strcmp("N",cmap->atmname[4]) == 0) strcpy(cmap->catmname[4], "OXT");

                } else if(strcmp("CMAP_RESIDX",tmpchar2)==0){
                    sscanf(sLine0,"%s%s%d%d%d%d%d",tmpchar1,tmpchar2,
                           &cmap->residx[0], &cmap->residx[1], &cmap->residx[2],
                           &cmap->residx[3], &cmap->residx[4]);
                           int l;
                           for (l=0; l<5; l++) {
                                cmap->cresidx[l] = cmap->residx[l];
                                cmap->nresidx[l] = cmap->residx[l];
                              }
                              cmap->nresidx[0] = 0;
                              cmap->cresidx[4] = 0;

                } else if(strcmp("CMAP_TERMMAP",tmpchar2)==0){
                    sscanf(sLine0,"%s%s%d",tmpchar1,tmpchar2,
                           &cmap->termmap);

                } else 
                {
                    VP0(("Unknown Flag : %s\n",sLine));
                }
            }
            if (strcmp("%COMMENT",tmpchar1) == 0 ) {
                if (cmnt != NULL) {
                    for (i=0;i<80 && sLine[i] != '\n'; i++) {
                        continue;
                    }
                    sLine[i]='\0';
                    cmnt->record=(char *) malloc(sizeof(char) * 256);
                    strcpy(cmnt->record, &sLine[9]);
                    cmnt->next=(CMNT *) malloc(sizeof(CMNT));
                    cmnt=cmnt->next;
                    cmnt->next=NULL;
                } else {
                    // Well, this may be a comment line not specific for a particular map.
                    // simply copy it.
                    cmntt->record=(char *) malloc(sizeof(char) * 256);
                    strcpy(cmntt->record, &sLine[9]);
                    cmntt->next=(CMNT *) malloc(sizeof(CMNT));
                    cmntt=cmntt->next;
                    cmntt->next=NULL;
                }
            }
        } else if (sLine[0]=='\0'){
            break;
        } else {
            //printf("%s\n",sLine);
            continue;
        }
    }
}

/*
 *      zAmberReadParmSetPropers
 *
 *      Read the proper torsion parameter terms.
 */
static void
zAmberReadParmSetPropers( PARMSET psParms, FILE *fIn )
{
STRING          sLine;
int             iRead, iN, nLine;
STRING          saStr[10], taStr[10];
double          dDivisions, dKp, dP0, dN, dScEE, dScNB, tScEE, tScNB;
char            *cScEE, *cScNB;

    memset(saStr, 0, sizeof(saStr));                    /* for Purify */
    memset(taStr, 0, sizeof(taStr));
    while (1) {
        FGETS( sLine, fIn );
        NODASHES(sLine);

// beginning of a group.
/* Read in parmxx.dat format
 *
 *        DI-DJ-DK-DL division kp PHI periodicity SCEE SCNB
 */
        iRead = sscanf( sLine, "%s %s %s %s %lf %lf %lf %lf %lf %lf", 
                                saStr[0], saStr[1], saStr[2], saStr[3], 
                                &dDivisions, &dKp, &dP0, &dN,
                                &dScEE, &dScNB );

        if ( iRead <= 0 ) 
                break;
        if ( sLine[0] == ' ' && sLine[1] == ' ') {
                VP0(( "bad proper torsion definition (skipping):\n (%s)\n", 
                                                                sLine ));
                continue;
        }
/* Try GLYCAM format
 *
 *        DI-DJ-DK-DL division kp PHI periodicity SCEE=scee SCNB=scnb
 */
        if (iRead < 9) {
        cScEE = strstr(sLine, "SCEE");
        if(cScEE!=NULL) {
                iRead += sscanf( cScEE, "SCEE=%lf", &dScEE);
                }

        cScNB = strstr(sLine, "SCNB");
        if(cScNB!=NULL) {
                iRead += sscanf( cScNB, "SCNB=%lf", &dScNB);
                }
         }

/*
 *  Neither SCEE nor SCNB has been set. Mark dScEE = 0.0, dScNB = 0.0 "not read".
 */
        if (iRead < 9) {
                dScEE = -1.0;
                dScNB = -1.0;
        } else if (iRead < 10) {
                dScNB = -1.0;
        }

/* For some reason, we must divide dKp by dDivisions */

        if ( dDivisions == 0.0 ) 
                dDivisions = 1.0;
        dKp /= dDivisions;
        iN = (int)floor(dN+0.5);
        zAmberConvertWildCard( saStr[0] );
        zAmberConvertWildCard( saStr[1] );
        zAmberConvertWildCard( saStr[2] );
        zAmberConvertWildCard( saStr[3] );

//        tScEE = 0.0;
//        tScNB = 0.0;
//        if ( iN > 0 ) {
//             tScEE = dScEE;
//             tScNB = dScNB;
//           }
        iParmSetAddProperTerm( psParms, 
                                saStr[0], saStr[1], saStr[2], saStr[3],
                                abs(iN), dKp, dP0*DEGTORAD, dScEE,
                                dScNB, "" );
        while( iN < 0 ) {
                FGETS( sLine, fIn );
                NODASHES(sLine);

// The same group of torsions share identical SCEE, SCNB

        iRead = sscanf( &sLine[11], "%lf %lf %lf %lf",
                &dDivisions, &dKp, &dP0, &dN);

        if ( iRead<=0 ) break;

        if ( dDivisions == 0.0 ) 
             dDivisions = 1.0;
        dKp /= dDivisions;
        iN = (int)floor(dN+0.5);
//        tScEE = 0.0;
//        tScNB = 0.0;
//        if ( iN > 0 ) {
//             tScEE = dScEE;
//             tScNB = dScNB;
//           }
        iParmSetAddProperTerm( psParms,
             saStr[0], saStr[1], saStr[2], saStr[3],
             abs(iN), dKp, dP0*DEGTORAD, dScEE, dScNB, "" );
        MESSAGE(( "Read extra term: %s\n", sLine ));
        }
        if ( iRead <= 0 )
                break;
    }
}


/*
 *      zAmberReadParmSetImpropers
 *
 *      Read the improper torsion parameter terms.
 */
static void
zAmberReadParmSetImpropers( PARMSET psParms, FILE *fIn )
{
STRING          sLine;
int             iRead, iN;
STRING          saStr[10];
double          dKp, dP0, dN, dScEE, dScNB;
BOOL            bPrintLine;

    memset(saStr, 0, sizeof(saStr));                    /* for Purify */
    while (1) {
        FGETS( sLine, fIn );
        NODASHES(sLine);
        /*
         *  get atoms & values, skipping possible 
         *      extraneous idivf value allowed by
         *      the input spec for PARM
         */
        iRead = sscanf( sLine, "%s %s %s %s",
                                saStr[0], saStr[1], saStr[2], saStr[3] );
        if ( iRead != 4 ) break;
        iRead = sscanf( sLine+15, "%lf %lf %lf", &dKp, &dP0, &dN );
        if ( iRead != 3 ) break;
        MESSAGE(( "Read: %s\n", sLine ));
        dN = (int)floor(dN+0.5);
        zAmberConvertWildCard( saStr[0] );
        zAmberConvertWildCard( saStr[1] );
        zAmberConvertWildCard( saStr[2] );
        zAmberConvertWildCard( saStr[3] );
        iN = (int)dN;
        dScEE = 0.0;
        dScNB = 0.0;
        /*
         *  check everything in case a format or other user error
         *      led to wrong values (e.g. IDIVF offset)
         */
        bPrintLine = FALSE;
        if ( dKp <= 0.0 ) {
            VP0((" WARNING: expected Improper Torsion PK>0 (%f)\n",
                dKp ));
            bPrintLine = TRUE;
        }
        if ( dP0 < 179.999  ||  dP0 > 180.001 ) {
            VP0((" WARNING: expected Improper Torsion PHASE=180 (%f)\n",
                dP0 ));
            bPrintLine = TRUE;
        }
        if ( iN < 1  ||  iN > 6 ) {
            VP0((" WARNING: unexpected Improper Torsion PN term (%d)\n",
                iN ));
            bPrintLine = TRUE;
        }
        if ( bPrintLine )
            VP0(( " Here is the Improper Torsion line in question:\n%s", 
                                                                sLine ));

        iParmSetAddImproperTerm( psParms, 
                                saStr[0], saStr[1], saStr[2], saStr[3],
                                iN, dKp, dP0*DEGTORAD, dScEE, dScNB, "" );
    }
    if ( iRead > 0 )
        VP0(( "WARNING: incomplete Improper Torsion line:\n%s", sLine ));
}

        


/*
 *      zAmberReadParmSetHBonds
 *
 *      Read hbond parameters.
 */
static void
zAmberReadParmSetHBonds( PARMSET psParms, FILE *fIn )
{
STRING          sLine;
int             iRead;
STRING          saStr[10];
double          dA, dB;

    memset(saStr, 0, sizeof(saStr));                    /* for Purify */
    while (1) {
        FGETS( sLine, fIn );
        iRead = sscanf( sLine, "%s %s %lf %lf", saStr[0], saStr[1], 
                                &dA, &dB );
        if ( iRead<=0 ) break;
        MESSAGE(( "Read: %s\n", sLine ));
        iParmSetAddHBond( psParms, saStr[0], saStr[1], dA, dB, "" );
    } 
}
 

/*
 *      zAmberReadParmSetNonBonds
 *
 *      Read non-bond parameters.
 */
static void
zAmberReadParmSetNonBonds( VARARRAY *vaPNonBonds, FILE *fIn )
{
STRING          sLine;
int             iRead;
STRING          saStr[10];
double          dRStar, dDepth, dRStar14, dDepth14;
NONBONDt        nbNonBond;
VARARRAY        vaNonBonds;

        memset(saStr, 0, sizeof(saStr));                        /* for Purify */
        memset(&nbNonBond, 0, sizeof(nbNonBond));               /* for Purify */
        vaNonBonds = vaVarArrayCreate(sizeof(NONBONDt));
        while (1) {
                FGETS( sLine, fIn );
                if( GDefaults.iCharmm ){
                        iRead = sscanf( sLine, "%s %lf %lf %lf %lf", 
                                saStr[0], &dRStar, &dDepth, &dRStar14, &dDepth14 );
                        if ( iRead <= 0 ) break;
                        MESSAGE(( "Read: %s\n", sLine ));
                        strcpy( nbNonBond.sNonBondType, saStr[0] );
                        nbNonBond.dRStar = dRStar;
                        nbNonBond.dDepth = dDepth;
                        nbNonBond.dRStar14 = dRStar14;
                        nbNonBond.dDepth14 = dDepth14;
                        VarArrayAdd( vaNonBonds, (GENP)&nbNonBond );
                } else {
                        iRead = sscanf( sLine, "%s %lf %lf", saStr[0], &dRStar, &dDepth );
                        if ( iRead <= 0 ) break;
                        MESSAGE(( "Read: %s\n", sLine ));
                        strcpy( nbNonBond.sNonBondType, saStr[0] );
                        nbNonBond.dRStar = dRStar;
                        nbNonBond.dDepth = dDepth;
                        nbNonBond.dRStar14 = dRStar;
                        nbNonBond.dDepth14 = dDepth;
                        VarArrayAdd( vaNonBonds, (GENP)&nbNonBond );
                }
        }
        *vaPNonBonds = vaNonBonds;
}


/*
 *      zAmberReadParmSetNBPairEdits
 *
 *      Author: David Cerutti (Case group, 2013)
 *
 *      Read non-bonded parameter pair changes specific to a force field
 *      which break the standard Lennard-Jones combining rules.
 */
static void
zAmberReadParmSetNBPairEdits( PARMSET psParms, FILE *fIn, int segfound )
{
  STRING          sLine;
  int             iRead;
  STRING          saStr[10];
  double          dEI, dEJ, dRI, dRJ;

  memset(saStr, 0, sizeof(saStr));                        /* for Purify */

  /* Scan the file for the card */
  while ( segfound == 0) {
    if (fgets( sLine, MAXSTRINGLENGTH, fIn ) == NULL) {
      break;
    }
    iRead = sscanf( sLine, "%s", saStr[0]);
    if (iRead > 0) {
      if (strcmp(saStr[0], "LJEDIT") == 0) {
	segfound = 1;
      }
    }
  }

  /* Scan the file for Lennard-Jones pair edits */
  while (segfound == 1) {
    if (fgets( sLine, MAXSTRINGLENGTH, fIn ) == NULL) {
      break;
    }
    iRead = sscanf( sLine, "%s %s %lf %lf %lf %lf", saStr[0], saStr[1],
		    &dRI, &dEI, &dRJ, &dEJ );
    if ( iRead <= 0 ) break;
    MESSAGE(( "Read: %s\n", sLine ));
    iParmSetAddNBEdit( psParms, saStr[0], saStr[1], dEI, dEJ, dRI, dRJ, "" );
  }
}


/*
 *      zbAmberDetermineParmSetFrcModType
 *
 *      Determine whether this is an old AMBER parameter set
 *      or if it is one of the new FRCMOD parameter sets.
 */
static BOOL
zbAmberDetermineParmSetFrcModType( FILE *fIn, BOOL *bPMass, BOOL *bPNonBond )
{
STRING          sLine;
BOOL            bNew;

    bNew = FALSE;
    *bPMass = FALSE;
    *bPNonBond = FALSE;
    do {
        FGETS( sLine, fIn );
        if ( strncmp( sLine, "MASS", 4 ) == 0 ) {
            bNew = TRUE;
            *bPMass = TRUE;
        } else if ( strncmp( sLine, "BOND", 4 ) == 0 ) {
            bNew = TRUE;
        } else if ( strncmp( sLine, "ANGL", 4 ) == 0 ) {
            bNew = TRUE;
        } else if ( strncmp( sLine, "DIHE", 4 ) == 0 ) {
            bNew = TRUE;
        } else if ( strncmp( sLine, "IMPR", 4 ) == 0 ) {
            bNew = TRUE;
        } else if ( strncmp( sLine, "HBON", 4 ) == 0 ) {
            bNew = TRUE;
        } else if ( strncmp( sLine, "NONB", 4 ) == 0 ) {
            bNew = TRUE;
            *bPNonBond = TRUE;
        }
    } while ( !feof(fIn) );
    fseek( fIn, 0, 0 );
    return(bNew);
}

/*
 *      iAtomTypeToElement
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Convert a type to an element number.
 *      If the type is unknown then return NOELEMENT.
 */
static int
iAtomTypeToElement( char *sType )
{
int             iElement, i, iCount;
STRING          sLowerType, sLowerTable;
TYPEELEMENTt    *tP;

    iCount = iVarArrayElementCount( vaAtomTypes );
    if ( iCount == 0 ) {
        VP1(( " (UNKNOWN ATOM TYPE: %s [no types loaded])\n", sType ));
        return( NOELEMENT );
    }

    iElement = NOELEMENT;
    strcpy( sLowerType, sType );
    StringLower( sLowerType );
    tP = PVAI( vaAtomTypes, TYPEELEMENTt, 0 );
    for (i=0; i<iCount; i++, tP++) {
        strcpy( sLowerTable, tP->sType );
        StringLower( sLowerTable ); 
        if ( strcmp( sLowerTable, sLowerType ) == 0 ) {
            iElement = tP->iElement;
            break;
        }
    }
    if ( iElement == NOELEMENT ) {
        VP1(( "(UNKNOWN ATOM TYPE: %s)\n", sType ));
    }
    return(iElement);
}




/*
 *      zpsAmberReadParmSetOld
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Read a PARMSET from an AMBER parameter database file.
 */
static PARMSET
zpsAmberReadParmSetOld( FILE *fIn )
{
STRING          sLine;
STRING          saStr[10];
VARARRAY        vaMasses, vaEquivs, vaNonBonds;
EQUIVt          eEquiv;
PARMSET         psParms;
NONBONDt        *nbNonBondP;
int             i, j, k, iRead, iSet;

    memset(saStr, 0, sizeof(saStr));                    /* for Purify */
    psParms = (PARMSET)oCreate(PARMSETid);
    
                /* Read and write the title */

    FGETS( sLine, fIn );
    VP0(("Reading title:\n%s", sLine ));

                /* Read the masses */

    zAmberReadParmSetMasses( &vaMasses, fIn );


                /* Skip hydrophilicity stuff */
    while (1) {
        FGETS( sLine, fIn );
        iRead = sscanf( sLine, 
                "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s",
                        saStr[0], saStr[0], saStr[0], saStr[0], saStr[0],
                        saStr[0], saStr[0], saStr[0], saStr[0], saStr[0],
                        saStr[0], saStr[0], saStr[0], saStr[0], saStr[0],
                        saStr[0], saStr[0], saStr[0], saStr[0], saStr[0] );
        if ( iRead < 20 ) break;
    } 

                /* Read the bond parameters */

    MESSAGE(( "Reading bond parameters\n" ));
    zAmberReadParmSetBonds( psParms, fIn );



                /* Read the angle parameters */
                
    MESSAGE(( "Reading angle parameters\n" ));
    zAmberReadParmSetAngles( psParms, fIn );



                /* Read the torsion parameters */
                
    MESSAGE(( "Reading torsion parameters\n" ));
    zAmberReadParmSetPropers( psParms, fIn );



                /* Read the improper parameters */

    MESSAGE(( "Reading improper parameters\n" ));
    zAmberReadParmSetImpropers( psParms, fIn );




                /* Read the H-bond parameters */
                
    MESSAGE(( "Reading H-Bond parameters\n" ));
    zAmberReadParmSetHBonds( psParms, fIn );






                /* Read the equivalences for non-bond parameters */

    MESSAGE(( "Reading non-bond equivalences\n" ));
    vaEquivs = vaVarArrayCreate(sizeof(EQUIVt));
    memset(&eEquiv, 0, sizeof(eEquiv));                 /* for Purify */
    while (1) {
        FGETS( sLine, fIn );
                        /* Remove the \n at the end */
        sLine[strlen(sLine)-1] = '\0';
        sRemoveLeadingSpaces( sLine );
        if ( strlen(sLine)==0 ) 
                break;
        MESSAGE(( "Read: %s\n", sLine ));
        sRemoveFirstString( sLine, eEquiv.sName );
        for (i=0; i<MAXEQUIV; i++) {
            sRemoveLeadingSpaces(sLine);
            if ( strlen(sLine)==0 ) 
                break;
            sRemoveFirstString( sLine, eEquiv.saEquivs[i] );
        }
        if ( strlen(sLine)!=0 ) 
                VP0(("Parmset: equiv %s: Unrecognized equivs (> %d)!!\n",
                        eEquiv.sName, MAXEQUIV ));
        eEquiv.iEquivs = i;
        VarArrayAdd( vaEquivs, (GENP)&eEquiv );
    }

                /* Read in the Non-Bond parameters, look them up in the */
                /* Masses table and then add them and their equivalences */
                /* to the PARMSET */

                /* Skip a line */
    FGETS( sLine, fIn );
    MESSAGE(( "Reading non-bond parameters\n" ));
    zAmberReadParmSetNonBonds( &vaNonBonds, fIn );

    /*
     *  correlate nb's & masses
     */
    nbNonBondP = PVAI( vaNonBonds, NONBONDt, 0 );
    for ( i=0; i<iVarArrayElementCount(vaNonBonds); i++, nbNonBondP++ ) {
        double          dMass, dPolar, dDepth, dRStar, dDepth14, dRStar14;
        double          dScreenF;
        int             iElement, iHybridization;
        MASSt           *MassP;
        EQUIVt          *EquivP;

        strcpy( saStr[0], nbNonBondP->sNonBondType );
        dRStar = nbNonBondP->dRStar;
        dDepth = nbNonBondP->dDepth;
        dRStar14 = nbNonBondP->dRStar14;
        dDepth14 = nbNonBondP->dDepth14;

                        /* Lookup the MASS if it is defined */
        /*
         *  look up mass/polar param TODO - use avl instead of vararray
         *      AND CHECK FOR DUPLICATE DEFN'S!
         */
        iSet = 0;
        MassP = PVAI(vaMasses,MASSt,0);
        for ( j=0; j<iVarArrayElementCount(vaMasses); j++, MassP++ ) {
            if ( strcmp( saStr[0], MassP->sType )==0 ) {
                dMass = MassP->dMass;
                dPolar = MassP->dPolar;
                dScreenF = MassP->dScreenF;
                iElement = MassP->iElement;
                iHybridization = MassP->iHybridization;
                iSet = 1;
            }
        }
        if ( iSet == 0 ) {
            VP0(( "No mass was defined for non-bond atom type: %s - ignoring\n",
                        saStr[0] ));
        } else {
                MESSAGE(( "Adding %d atom type: %s  %lf %lf %lf %lf\n",
                        i, saStr[0], dMass, dPolar, dRStar, dDepth ));
                iParmSetAddAtom( psParms, saStr[0], dMass, dPolar, 
                                dDepth, dRStar, dDepth14, dRStar14,
                                dScreenF,
                                iElement!=-10?iElement
                                              :iAtomTypeToElement(saStr[0]),
                                iHybridization!=HUNDEFINED?iHybridization
                                              :iAtomTypeHybridization(saStr[0]),
                                "" );
        }
        EquivP = PVAI(vaEquivs,EQUIVt,0);
        for ( j=0; j<iVarArrayElementCount(vaEquivs); j++, EquivP++ ) {
            if ( strcmp( saStr[0], EquivP->sName )==0 ) {
                for ( k=0; k<EquivP->iEquivs; k++ ) {
                    if ( iSet == 0 ) {
                        VP0(( "No mass equivalenced type: %s - skipping\n",
                                    EquivP->saEquivs[k] ));
                        continue;
                    }
                    iParmSetAddAtom( psParms, EquivP->saEquivs[k],
                       dMass, dPolar, dDepth, dRStar, dDepth14, dRStar14,
                       dScreenF,
                       iElement!=-10?iElement
                                    :iAtomTypeToElement(saStr[0]),
                       iHybridization!=HUNDEFINED?iHybridization
                                    :iAtomTypeHybridization(saStr[0]),
                       "" );
                }
            }
        }
    }

    /* Read in non-bonded pair interaction edits, check their */
    /* validity against existing atom types, and add them to  */
    /* the parmset.                                           */
    MESSAGE(( "Seeking non-bond pair edits\n" ));
    zAmberReadParmSetNBPairEdits( psParms, fIn, 0 );

    VarArrayDestroy( &vaMasses );
    VarArrayDestroy( &vaEquivs );
    VarArrayDestroy( &vaNonBonds );
    fclose(fIn);

    return(psParms);
}




/*
 *      zpsAmberReadParmSetFrcMod
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Read a PARMSET from an AMBER parameter database file.
 */
static PARMSET
zpsAmberReadParmSetFrcMod( FILE *fIn )
{
STRING          sLine, saStr[10];
VARARRAY        vaMasses, vaNonBonds;
    VARARRAY        vaFoo;
PARMSET         psParms;
int             lastsegment;
//NONBONDt        nbNonBond;
int             i, j;

    memset(saStr, 0, sizeof(saStr));                    /* for Purify */
    psParms = (PARMSET)oCreate(PARMSETid);
    
                /* Read and write the title */

    FGETS( sLine, fIn );
    VP0(("Reading title:\n%s", sLine ));

    vaMasses = NULL;
    vaNonBonds = NULL;
    lastsegment = 0;

    do {
        FGETS( sLine, fIn );
        if ( strncmp( sLine, "MASS", 4 ) == 0 ) {
            zAmberReadParmSetMasses( &vaMasses, fIn );
            lastsegment = 1;
        } else if ( strncmp( sLine, "BOND", 4 ) == 0 ) {
            zAmberReadParmSetBonds( psParms, fIn );
            lastsegment = 2;
        } else if ( strncmp( sLine, "ANGL", 4 ) == 0 ) {
            zAmberReadParmSetAngles( psParms, fIn );
            lastsegment = 3;
        } else if ( strncmp( sLine, "DIHE", 4 ) == 0 ) {
            zAmberReadParmSetPropers( psParms, fIn );
            lastsegment = 4;
        } else if ( strncmp( sLine, "IMPR", 4 ) == 0 ) {
            zAmberReadParmSetImpropers( psParms, fIn );
            lastsegment = 5;
        } else if ( strncmp( sLine, "HBON", 4 ) == 0 ) {
            zAmberReadParmSetHBonds( psParms, fIn );
            lastsegment = 6;
        } else if ( strncmp( sLine, "NONB", 4 ) == 0 ) {
            zAmberReadParmSetNonBonds( &vaNonBonds, fIn );
            lastsegment = 7;
        } else if ( strncmp( sLine, "CMAP", 4 ) == 0 ) {
            zAmberReadParmSetCMAP( &vaFoo, fIn );
            lastsegment = 8;
        } else if ( strncmp( sLine, "IPOL", 4 ) == 0 ) {
            zAmberReadParmSetIPOL( fIn );
            lastsegment = 9;
	} else if ( strncmp( sLine, "LJEDIT", 6 ) == 0 ) {
	  zAmberReadParmSetNBPairEdits( psParms, fIn, 1 );
	  lastsegment = 10;
        } else if ( strncmp( sLine, "END", 3 ) == 0 ) {
                continue;
        } else {
            sLine[10] = '\0';
            if ( strlen(sLine) > 1 ) {
                VP0(( "Unknown keyword: %s in parameter file. Perhaps a format error?\n", sLine ));
            }
        }
    } while ( !feof(fIn) );

                /* If there are NonBond parameters then use them */
                /* this implies that there are MASSES */

    if ( vaNonBonds ) {
        NONBONDt        *nbP = PVAI( vaNonBonds, NONBONDt, 0);

        for ( i=0; i<iVarArrayElementCount(vaNonBonds); i++, nbP++ ) {
            double      dMass, dPolar, dRStar, dDepth, dRStar14, dDepth14;
            double      dScreenF;
            int         iElement, iHybridization;
            int         iSet = 0;
            MASSt       *mP;

            strcpy( saStr[0], nbP->sNonBondType );
            dRStar = nbP->dRStar;
            dDepth = nbP->dDepth;
            dRStar14 = nbP->dRStar14;
            dDepth14 = nbP->dDepth14;

                            /* Lookup the MASS if it is defined */
            mP = PVAI(vaMasses,MASSt,0);
            for ( j=0; j<iVarArrayElementCount(vaMasses); j++, mP++ ) {
                if ( strcmp( saStr[0], mP->sType )==0 ) {
                        dMass = mP->dMass;
                        dPolar = mP->dPolar;
                        dScreenF = mP->dScreenF;
                        iElement = mP->iElement;
                        iHybridization = mP->iHybridization;
                        iSet = 1;
                        break;
                }
            }
            if ( iSet == 0 ) {
                VP0(( "No mass defined for non-bond atom type: %s - skipping\n",
                            saStr[0] ));
                continue;
            }
            MESSAGE(( "Adding %d atom type: %s  %lf %lf %lf %lf\n",
                            i, saStr[0], dMass, dPolar, dRStar, dDepth ));
            iParmSetAddAtom( psParms, saStr[0], dMass, dPolar, dDepth, dRStar,
                    dDepth14, dRStar14, dScreenF,
                                iElement!=-10?iElement
                                             :iAtomTypeToElement(saStr[0]),
                                iHybridization!=HUNDEFINED?iHybridization
                                             :iAtomTypeHybridization(saStr[0]),
                                    "" );
        }
    }

    if ( vaMasses )
        VarArrayDestroy( &vaMasses );
    if ( vaNonBonds )
        VarArrayDestroy( &vaNonBonds );
    fclose(fIn);

    return(psParms);
}



/*
 *      uAmberReadUnitFromPrep
 *
 *      Author: Christian Schafmeister (1991)
 *
 *TODO: Currently this routine does not do a distance search to
 *TODO: create bonds between atoms a specified distance apart.
 *TODO: The PREP file contains the distance.
 *
 *      Read a UNIT from a PREP file.
 *      If the first line in the file is 'STOP' then return NULL.
 */
typedef struct {
        VECTOR  vPos;
        ATOM    aAtom;
        int     iUnfilled;
        int     iBack;
} TREENODEt;

                /* Define a special type for sorting prep format lines */
typedef struct  {
        int     iLineNumber;
        STRING  sLine;
} INPUTLINEt;

static UNIT
uAmberReadUnitFromPrep( FILE *fIn )
{
TREENODEt       *cPCoor;
TREENODEt       cCoor;
STRING          sLine, sChain;
STRING          saStr[15];
int             iRead, iIndex;
UNIT            uUnit = NULL;
RESIDUE         rRes;
ATOM            aAtom, aAtom2;
VARARRAY        vaCoor = NULL;
VARARRAY        vaLines = NULL;
STRING          sName, sType, sSave;
double          dX, dY, dZ, dCharge;
int             iBond, iAngle, iTorsion;
double          dBond, dAngle, dTorsion;
VECTOR          vPos, vPos1, vPos2, vPos3;
ATOM            aMain0, aMain1;
BOOL            bFirstTime;
double          da[10];
int             i, iLen, iErr, iCharge, iCharges = 0, iChgWarn = 0, iAtoms = 0;
double          dCutoff;
INPUTLINEt      ilLine;
int             iaTreeStack[500];
int             iTreeStackTop;
BOOL            bUseFirstColumn;
#define         TREEPOP()       ( iTreeStackTop-- )
#define         TREEPUSH( i )   ( iaTreeStack[++iTreeStackTop] = i )
#define         iTREETOP()      ( iaTreeStack[iTreeStackTop] )
#define         iTREESTACKPOS() ( iTreeStackTop )
#define         bTREEATTOP()    ( iTreeStackTop == 0 )

                /* Read the long name of the UNIT and null-terminate it */
    FGETS( sLine, fIn );
    if ( strncmp( sLine, "STOP", 4 ) == 0 ) {
        return(NULL);
    }
    if ( feof( fIn ) ) {
        VP0(( "Unexpected EOF: Assuming STOP\n" ));
        return(NULL);
    }
    sLine[strlen(sLine)-1] = '\0';

    uUnit = (UNIT)oCreate(UNITid);
    rRes  = (RESIDUE)oCreate(RESIDUEid);
    ContainerAdd( (CONTAINER)uUnit, (OBJEKT)rRes );

    ResidueSetDescription( rRes, sLine );

                /* Skip a line */
    FGETS( sLine, fIn );
    if ( feof( fIn ) ) {
        VP0(( "Unexpected EOF: discarding residue (%s)\n", 
                                                sResidueDescription(rRes) ));
        goto ERROR;
    }

                /* Define an array to store the coordinates generated */
                /* from the tree-matrix */

    vaCoor = vaVarArrayCreate( sizeof(TREENODEt) );

                /* Read in the short name, INT/XYZ, 0/1 */

    FGETS( sLine, fIn );
    if ( feof( fIn ) ) {
        VP0(( "Unexpected EOF: discarding residue (%s)\n",
                                                sResidueDescription(rRes) ));
        goto ERROR;
    }
    sscanf( sLine, "%s %s %s", saStr[1], saStr[2], saStr[3] );
    ContainerSetName( uUnit, saStr[1] );
    ContainerSetName( rRes, saStr[1] );

                /* Read the line: CORR/CHANGE, OMIT/NOMIT, DU, ALL/BEG */

    FGETS( sLine, fIn );
    if ( feof( fIn ) ) {
        VP0(( "Unexpected EOF: discarding residue (%s)\n",
                                                sResidueDescription(rRes) ));
        goto ERROR;
    }
    sscanf( sLine, "%s", saStr[1] );
    bUseFirstColumn = TRUE;
    if ( strncmp( saStr[1], "CORR", 4 ) == 0 ) {
        bUseFirstColumn = FALSE;
    }

                /* Read the CUT line */

    FGETS( sLine, fIn );
    if ( feof( fIn ) ) {
        VP0(( "Unexpected EOF: discarding residue (%s)\n",
                                                sResidueDescription(rRes) ));
        goto ERROR;
    }
    sscanf( sLine, "%lf", &dCutoff );


                /* Read all of the lines defining the atoms so that */
                /* we can sort then by the first number */

    vaLines = vaVarArrayCreate(sizeof(INPUTLINEt));
    memset(ilLine.sLine, 0, sizeof(ilLine.sLine));              /* for Purify */
    for (;;) {
        FGETS( sLine, fIn );
        if ( feof( fIn ) ) {
                VP0(( "Unexpected EOF: discarding residue (%s)\n",
                                                sResidueDescription(rRes) ));
                return(NULL);
        }
        iRead = sscanf( sLine, "%s %s %s %s %s %s %s %s %s %s %s",
                        saStr[1], saStr[2], saStr[3], saStr[4], saStr[5],
                        saStr[6], saStr[7], saStr[8], saStr[9], saStr[10],
                        saStr[11] );
        if ( iRead <= 0 ) 
                break;
        ilLine.iLineNumber = atoi(saStr[1]);
        strcpy( ilLine.sLine, sLine );
        MESSAGE(( "iRead = %d   Adding: |%s|\n", iRead, sLine ));
        VarArrayAdd( vaLines, (GENP)&ilLine );
    }

        /* Save the last line read, because it will be a command */
        /* that has to be parsed after all the atoms are defined */

    strcpy( sSave, sLine );

#ifdef  DEBUG
MESSAGE(( "------------\n" ));
for ( i=0; i<iVarArrayElementCount(vaLines); i++ ) {
    MESSAGE(( "%s", PVAI(vaLines,INPUTLINEt,i)->sLine ));
}
MESSAGE(( "============\n" ));
#endif
    if ( bUseFirstColumn ) {
        INPUTLINEt      *Pil = PVAI(vaLines,INPUTLINEt,0);

        SortByInteger( (GENP) Pil,
                        iVarArrayElementCount(vaLines),
                        sizeof(INPUTLINEt),
                        (GENP) &Pil->iLineNumber, TRUE );
    }
#ifdef  DEBUG
MESSAGE(( "------------\n" ));
for ( i=0; i<iVarArrayElementCount(vaLines); i++ ) {
    MESSAGE(( "%s", PVAI(vaLines,INPUTLINEt,i)->sLine ));
}
MESSAGE(( "============\n" ));
#endif


                /* Read in ATOM lines until the first string is 
                 * the keyword LOOP, IMPROPER, DONE, CHARGE 
                 * Create atoms for every entry except for those of 
                 * type DU 
                 * Read the tree type records and build a tree using 
                 * the back indices (iBack) in the vaCoor array 
                 * and the index stack iaTreeStack/iTreeStackTop 
                 * Keep adding back indices to a node until the 
                 * iUnfilled field is zero, then pop the stack */

    iTreeStackTop = -1;
    TREEPUSH(-1);
    MESSAGE(( "Reading atom records\n" ));
    aMain0 = NULL;
    aMain1 = NULL;
    iIndex = 0;
    memset(&cCoor, 0, sizeof(cCoor));                   /* for Purify */
    for ( i=0; i<iVarArrayElementCount(vaLines); i++ ) {
        strcpy( sLine, PVAI( vaLines, INPUTLINEt, i )->sLine );
        iRead = sscanf( sLine, "%s %s %s %s %s %s %s %s %s %s %s",
                        saStr[1], saStr[2], saStr[3], saStr[4], saStr[5],
                        saStr[6], saStr[7], saStr[8], saStr[9], saStr[10],
                        saStr[11] );
        if ( iRead <= 0 ) {
            if ( feof(fIn) ) goto DONE;
            continue;
        }
        MESSAGE(( "Read atom: %s\n", saStr[2] ));
        iAtoms++;

                /* Create an entry in the vaCoor array */

        cCoor.iBack = 0;
        cCoor.iUnfilled = 0;
        VarArrayAdd( vaCoor, (GENP)&cCoor );
        cPCoor = PVAI( vaCoor, TREENODEt, iIndex );

        iErr = 0;
        strcpy( sName, saStr[2] );
        iLen = strlen( sName );
        if ( iLen > 4 ) {
            VP0(( "residue %s atom %s: name must be <= 4 characters\n",
                                sResidueDescription(rRes), sName ));
            iErr++;
        }

        strcpy( sType, saStr[3] );
        iLen = strlen( sType );
        if ( iLen > 2 ) {
            VP0(( 
                "residue %s atom %s type %s: type must be 1 or 2 characters\n",
                sResidueDescription(rRes), sName, sType ));
            iErr++;
        }

        strcpy( sChain, saStr[4] );
        iLen = strlen( sChain );
        if ( iLen > 1 ) {
            VP0(( "residue %s atom %s: chain type [%s] must be 1 character\n",
                                sResidueDescription(rRes), sName, sChain ));
            iErr++;
        }

        if ( iErr )
            goto ERROR;


                /* Define the node type */

        switch ( sChain[0] ) {
            case 'M':
                cPCoor->iUnfilled = 10000;
                break;
            case 'B':
                cPCoor->iUnfilled = 2;
                break;
            case '3':
                cPCoor->iUnfilled = 3;
                break;
            case 'S':
                cPCoor->iUnfilled = 1;
                break;
            case 'E':
                cPCoor->iUnfilled = 0;
                break;
            case '4':
                cPCoor->iUnfilled = 4;
                break;
            case '5':
                cPCoor->iUnfilled = 5;
                break;
            case '6':
                cPCoor->iUnfilled = 6;
                break;
            default:
                VP0(( "Atom %s: Illegal chain specifier [%c] in PREP file.\n", 
                        saStr[2], sChain[0] ));
                VP0(( "%s\n", sLine ));
                goto ERROR;
        }

                /* Point the new node back to the top node on the stack */
                /* If the top node is not accepting any more */
                /* connections then pop it. */

        cPCoor->iBack = iTREETOP();
        if ( iTREETOP() >= 0 ) {
            PVAI( vaCoor, TREENODEt, iTREETOP() )->iUnfilled--;
            if ( PVAI( vaCoor, TREENODEt, iTREETOP() )->iUnfilled == 0 ) {
                TREEPOP();
                if ( iTREESTACKPOS() <= 0 ) DFATAL(( "Stack is empty" ));
            }
        }

                /* Now if the current node is unfilled then push it on */
                /* the tree stack */

        if ( cPCoor->iUnfilled != 0 ) TREEPUSH( iIndex );

                /* Get the coordinates for the atom */

        iBond = cPCoor->iBack;
        if ( iRead == 7 ) {
            bStringToDouble( saStr[5], &dX );
            bStringToDouble( saStr[6], &dY );
            bStringToDouble( saStr[7], &dZ );
            dCharge = 0.0;
        } else if ( iRead == 8 ) {
            bStringToDouble( saStr[5], &dX );
            bStringToDouble( saStr[6], &dY );
            bStringToDouble( saStr[7], &dZ );
            bStringToDouble( saStr[8], &dCharge );
            iCharges++;
        } else {
            if ( iBond >= 0 ) {
                iAngle= PVAI( vaCoor, TREENODEt, iBond )->iBack;
                if ( iAngle >= 0 ) {
                    iTorsion = PVAI( vaCoor, TREENODEt, iAngle )->iBack;
                }
            }
            bStringToDouble( saStr[8], &dBond );
            bStringToDouble( saStr[9], &dAngle );
            dAngle *= DEGTORAD;
            bStringToDouble( saStr[10], &dTorsion );
            dTorsion *= DEGTORAD;
            if ( iBond == -1 ) {
                ZMatrixNothing( &vPos );
            } else if ( iAngle == -1 ) {
                vPos1 = PVAI( vaCoor, TREENODEt, iBond )->vPos;
                ZMatrixBond( &vPos, &vPos1, dBond );
            } else if ( iTorsion == -1 ) {
                vPos1 = PVAI( vaCoor, TREENODEt, iBond )->vPos;
                vPos2 = PVAI( vaCoor, TREENODEt, iAngle )->vPos;
                ZMatrixBondAngle( &vPos, &vPos1, &vPos2, dBond, dAngle );
            } else {
                vPos1 = PVAI( vaCoor, TREENODEt, iBond )->vPos;
                vPos2 = PVAI( vaCoor, TREENODEt, iAngle )->vPos;
                vPos3 = PVAI( vaCoor, TREENODEt, iTorsion )->vPos;
                ZMatrixBondAngleTorsion( &vPos, &vPos1, &vPos2, &vPos3,
                                dBond, dAngle, dTorsion );
            }
            dX = dVX(&vPos);
            dY = dVY(&vPos);
            dZ = dVZ(&vPos);
            if (iRead < 11)
                dCharge = 0.0;
            else {
                bStringToDouble( saStr[11], &dCharge );
                iCharges++;
            }
/* printf("--- %f [%d]\n", dCharge, iRead);*/
/* %%% */
        }
        if ( iIndex >= iVarArrayElementCount(vaCoor) ) {
            VarArraySetSize( vaCoor, iIndex+1 );
        }
        cPCoor = PVAI( vaCoor, TREENODEt, iIndex );
        VectorDef( &(cPCoor->vPos), dX, dY, dZ );
        aAtom = NULL;
        if ( strcmp( sType, "DU" ) != 0 ) {
            MESSAGE(( "Creating atom: %s\n", sName ));
            aAtom = (ATOM)oCreate(ATOMid);
            ContainerSetName( aAtom, sName );
            AtomSetElement( aAtom, iAtomTypeToElement(sType) );
            AtomSetType( aAtom, sType );
            VectorDef( &vPos, dX, dY, dZ );
            AtomSetPosition( aAtom, vPos );
            AtomSetCharge( aAtom, dCharge );
            ContainerAdd( (CONTAINER)rRes, (OBJEKT)aAtom );
            if ( iBond != -1 ) {
                aAtom2 = (PVAI( vaCoor, TREENODEt, iBond))->aAtom;
                if ( aAtom2 != NULL ) {
                    MESSAGE(( "Created bond between: %s - %s\n",
                                sContainerName(aAtom),
                                sContainerName(aAtom2) ));
                    AtomBondTo( aAtom, aAtom2 );
                }
            }
            if ( (strcmp( sChain, "M" )==0) && (strcmp( sType, "DU" )!=0) ) {
                if ( aMain0 == NULL ) aMain0 = aAtom;
                aMain1 = aAtom;         
            }
        }
        else {
            if ( strncmp( sName, "DUMM", 4 ) != 0 ) {
                VP0(( "Entry of type DUmmy; atom %s not created\n", sName ));
            }
        }
        cPCoor->aAtom = aAtom;
        iIndex++;
    }

    if ( iCharges == 0 ) {
        VP0(("(no charges read on atoms lines in %s)\n",
                                                sResidueDescription(rRes) ));
    } else {
        if ( iCharges != iAtoms )
            VP0((
              " WARNING - %d of %d atoms missing charges on atoms lines: %s\n",
                                        iAtoms-iCharges, iAtoms,
                                        sResidueDescription(rRes) ));
        iChgWarn = 1;
        iCharges = 0;
    }

        /* Get the last line read which is the first command */

    strcpy( sLine, sSave );

    MESSAGE(( "Parsing commands after atom records\n" ));

                /* Parse the commands that follow the atom records */
    bFirstTime = TRUE;
    while (1) {
        if ( bFirstTime ) bFirstTime = FALSE;
        else            FGETS( sLine, fIn );
        if ( feof( fIn ) ) {
                VP0(( "Unexpected EOF: discarding residue (%s)\n",
                                                sResidueDescription(rRes) ));
                return(NULL);
        }
        iRead = sscanf( sLine, "%s", saStr[1] );
        if ( iRead <= 0 ) continue;
        MESSAGE(( "Parsed command: %s\n", saStr[1] ));
        TESTMEMORY();
        if ( strcmp( saStr[1], "DONE" ) == 0 ) break;
        else if ( strcmp( saStr[1], "LOOP" ) == 0 ) {
            int problem = 0;
            while (1) {
                FGETS( sLine, fIn );
                if ( feof( fIn ) ) {
                        VP0(("Unexpected EOF: discarding (%s)\n",
                                                sResidueDescription(rRes) ));
                        goto ERROR;
                }
                iRead = sscanf( sLine, "%s %s", saStr[1], saStr[2] );
                if ( iRead <= 0 ) break;
                if ( iRead != 2 ) {
                    VP0(( "** only got 1 atom (%s) from LOOP line - skipping\n",
                                        sLine ));
                    problem++;
                    continue;
                }
                aAtom = (ATOM)cContainerFindName( (CONTAINER)rRes, 
                                                        ATOMid, saStr[1] );
                aAtom2= (ATOM)cContainerFindName( (CONTAINER)rRes, 
                                                        ATOMid, saStr[2] );
                if ( aAtom == NULL ) {
                    VP0(( "** LOOP atom %s not found - bond not formed\n",
                            saStr[1] ));
                    problem++;
                    continue;
                }
                if ( aAtom2 == NULL ) {
                    VP0(( "** LOOP atom %s not found - bond not formed\n",
                            saStr[2] ));
                    problem++;
                    continue;
                }
                if (bAtomBondedTo(aAtom, aAtom2)) {
                    VP0(( "%s:  LOOP: redundant bond  %s--%s  ignored\n",
                        sResidueDescription(rRes), saStr[1], saStr[2] ));
                } else
                        AtomBondTo( aAtom, aAtom2 );
            }
            if (problem) {
                VP0(( "Discarding residue (%s) to EOF\n",
                                        sResidueDescription(rRes) ));
                goto ERROR;
            }
        } else if ( strcmp( saStr[1], "CHARGE" )==0 ) {
            if ( iChgWarn )
                VP0((
                "Warning: per-line charges being overridden by CHARGE block in %s\n",
                                                sResidueDescription(rRes) ));
            iCharge = 3;
            while (1) {
                FGETS( sLine, fIn );
                if ( feof( fIn ) ) {
                        VP0(( "Unexpected EOF: discarding residue (%s)\n",
                                sResidueDescription(rRes) ));
                        goto ERROR;
                }
                iRead = sscanf( sLine, "%lf %lf %lf %lf %lf", 
                                &da[0], &da[1], &da[2],
                                &da[3], &da[4] );
                if ( iRead<=0 ) break;
                for ( i=0; i<iRead; i++ ) {
                    if ( iCharge >= iIndex ) {
                        VP0((" WARNING - skipping extra charge in: %s\n",
                                                sResidueDescription(rRes) ));
                        continue;
                    }
                    AtomSetCharge( PVAI( vaCoor, TREENODEt, iCharge )->aAtom,
                                        da[i] );
                    iCharge++;
                }
            }
            if ( iCharge < iAtoms ) {
                    VP0(( " ERROR - incomplete CHARGE block in: %s: %d/%d\n",
                                                sResidueDescription(rRes),
                                                iCharge, iAtoms ));
                    goto ERROR;
            }
        } else if ( strcmp( saStr[1], "IMPROPER" )==0 ) {
            int problem = 0;
            IMPROPERt   Imp;

                /* maybe will use this info someday; for now this is
                   the 1st step towards saving it
                        (debate continues over whether it has value)  */

            rRes->vaImpropers = vaVarArrayCreate( sizeof(IMPROPERt) );
            while (1) {
                FGETS( sLine, fIn );
                if ( feof( fIn ) ) {
                        VP0(( "Unexpected EOF: discarding residue (%s)\n",
                                                sResidueDescription(rRes) ));
                        goto ERROR;
                }
                iRead = sscanf( sLine, "%s %s %s %s", 
                                saStr[1], saStr[2], saStr[3], saStr[4] );
                if ( iRead <= 0 ) break;
                if ( iRead != 4 ) {
                    VP0(( "** got %d atom%s from IMPROPER line (%s) %s\n",
                                iRead, (iRead>1 ? "s":""), sLine,
                                " - expected 4 (skipping line)" ));
                    problem++;
                    continue;
                }
                for ( i=1; i<=4; i++ ) {
                    if ( strcmp(saStr[i], "-M") && strcmp(saStr[i], "+M") &&
                             !cContainerFindName( (CONTAINER)rRes, ATOMid, saStr[i] ) ) {
                        VP0(( "** %s: 'IMPROPER' atom %s not found\n", 
                                sResidueDescription(rRes), saStr[i] ));
                        problem++;
                    }
                }
                if ( !problem ) {
                    memset(&Imp, ' ', sizeof(Imp));
                    memcpy(Imp.sName1, saStr[1], strlen(saStr[1]));
                    memcpy(Imp.sName2, saStr[2], strlen(saStr[2]));
                    memcpy(Imp.sName3, saStr[3], strlen(saStr[3]));
                    memcpy(Imp.sName4, saStr[4], strlen(saStr[4]));
                    VarArrayAdd( rRes->vaImpropers, (GENP)&Imp );
                }
            }
            if (problem && 0) {
                VP0(( "Discarding residue (%s) to EOF\n",
                                        sResidueDescription(rRes) ));
                goto ERROR;
            }
        }
    }

    TESTMEMORY();

    UnitSetHead( uUnit, aMain0 );
    ResidueSetConnectAtom( rRes, NEND, aMain0 );
    UnitSetTail( uUnit, aMain1 );
    ResidueSetConnectAtom( rRes, CEND, aMain1 );

                /* Check if there was a cutoff value, if there was */
                /* use it in a Distance Search to create bonds */

    if ( bUseFirstColumn == TRUE  &&  dCutoff > 0.01 ) {
        VP1(( "Distance search to create bonds for: %s  distance: %10.5lf\n",
                sContainerName(uUnit), dCutoff ));
        iToolDistanceSearch( (CONTAINER)uUnit, dCutoff, TRUE, 
                                DISTANCE_SEARCH_CREATE_BONDS );
    }

DONE:
    VarArrayDestroy( &vaCoor );
    VarArrayDestroy( &vaLines );

    return(uUnit);

ERROR:
    VarArrayDestroy( &vaCoor );
    VarArrayDestroy( &vaLines );

    Destroy( (OBJEKT *)&uUnit );

    return(NULL);
}


/*
 *------------------------------- 
 *       Public routines
 */


/*
 *      AmberAddAtomTypes
 *
 *      Author: Bill Ross (1996)
 *
 *      Add to array of atom types / hybridizations
 */
void    
AmberAddAtomTypes( LIST lEntries )
{
LISTLOOP        llEntries, llEntry;
ASSOC           aEntry1, aEntry2;
LIST            lEntry;
OBJEKT          oEntry2;
TYPEELEMENTt    tTypeInfo;

    if ( iObjectType(lEntries) != LISTid ) {
        DFATAL(( "AmberAddAtomTypes: Need LIST" ));
    }
    if ( vaAtomTypes == NULL )
        vaAtomTypes = vaVarArrayCreate( sizeof( TYPEELEMENTt ) );

    llEntries = llListLoop(lEntries);
    while ( (aEntry1 = (ASSOC)oListNext(&llEntries)) ) {
        /*
         *  prepare to go over the Type, Element, Hybridization
         *      items - list items w/in the overall list
         */
        lEntry = (LIST)oAssocObject(aEntry1);   /* view current entry as list */
        llEntry = llListLoop(lEntry);           /* set up to loop over it */
        aEntry2 = (ASSOC)oListNext(&llEntry);   /* get 1st thing */

        if ( aEntry2 == NULL ) {
                VP0(( "addAtomTypes: null atom type entry\n" ));
                continue;
        }
        oEntry2 = (OBJEKT)oAssocObject(aEntry2); /* get the objekt it is */
        if ( iObjectType(oEntry2) != OSTRINGid ) {
                VP0(( "addAtomTypes: atom type - expected string\n" ));
                continue;
        }
#if 0
        /*  Joung-Cheatham use ion types with three characters,
            so skip this warning  */
        if ( strlen( sOString( oEntry2 ) ) > 2 ) {
                VP0(( "addAtomTypes: type %s - max length is 2\n", 
                                        sOString( oEntry2 ) ));
                continue;
        }
#endif
        memset(&tTypeInfo, 0, sizeof(tTypeInfo));       /* for Purify */
        strcpy( tTypeInfo.sType, sOString( oEntry2 ) );

        aEntry2 = (ASSOC)oListNext(&llEntry);   /* get next thing */
        if ( aEntry2 == NULL ) {
                VP0(( "addAtomTypes: incomplete type entry (%s)\n", 
                                        tTypeInfo.sType ));
                continue;
        }
        oEntry2 = (OBJEKT)oAssocObject(aEntry2); /* get the objekt it is */
        if ( iObjectType(oEntry2) != OSTRINGid ) {
                VP0(( "addAtomTypes: %s - expected string for element\n",
                                                tTypeInfo.sType ));
                continue;
        }
        if ( strlen( sOString( oEntry2 ) ) == 0 ) {
                /* NULL ("") entry for dummy atoms */
                tTypeInfo.iElement = NOELEMENT;
        } else {
                tTypeInfo.iElement = iElementNumber( sOString( oEntry2 ) );
                if ( tTypeInfo.iElement == NOELEMENT ) {
                        VP0(( "atom type %s - unknown element %s\n",
                                                tTypeInfo.sType,
                                                sOString( oEntry2 ) ));
                        continue;
                }
        }
        aEntry2 = (ASSOC)oListNext(&llEntry);   /* get next thing */
        if ( aEntry2 == NULL ) {
                VP0(( "incomplete atomtype entry (%s)\n", tTypeInfo.sType ));
                continue;
        }
        oEntry2 = (OBJEKT)oAssocObject(aEntry2); /* get the objekt it is */
        if ( iObjectType(oEntry2) == OSTRINGid ) {
                StringLower( sOString( oEntry2 ) );
                if ( !strcmp( sOString( oEntry2 ), "sp2") )
                        tTypeInfo.iHybridization = HSP2;
                else if ( !strcmp( sOString( oEntry2 ), "sp3") )
                        tTypeInfo.iHybridization = HSP3;
                else {
                        VP0(( "atom type %s - unknown hybridization %s\n",
                                                tTypeInfo.sType, sOString( oEntry2 ) ));
                        continue;
                }
        } else if ( iObjectType(oEntry2) == OINTEGERid ) {
                tTypeInfo.iHybridization = iOInteger( oEntry2 );
                if (tTypeInfo.iHybridization != HSP2  &&
                    tTypeInfo.iHybridization != HSP3  &&
                    tTypeInfo.iHybridization != 0) {
                        VP0(( "atom type %s - unknown hybridization %d\n",
                                                tTypeInfo.iHybridization ));
                        continue;
                }
        } else {
          VP0(( "atom type %s - expected string or integer for hybridization\n",
                                                tTypeInfo.sType ));
          continue;
        }
        VarArrayAdd( vaAtomTypes, (GENP)&tTypeInfo );
    }

    /*
     *  now delete it all.. except the top list, which parser.y does,
     *  but this needs to be DEREF'd for some reason.. the overall
     *  structure:
     *
     *          assoc -> list(lEntries) 
     *                  [commands.c passes lEntries to amber.c;
     *                   both the assoc & lEntries destroyed in parser.y]
     *
     *          lEntries: list of assoc, each containing a list of strings.
     */
    llEntries = llListLoop(lEntries);
    while ( (aEntry1 = (ASSOC)oListNext(&llEntries)) ) {
        lEntry = (LIST)oAssocObject(aEntry1);
        llEntry = llListLoop(lEntry);
        while ( (aEntry2 = (ASSOC)oListNext(&llEntry)) ) {
                Destroy( (OBJEKT *)&aEntry2 );
        }
        Destroy( (OBJEKT *)&aEntry1 );
        Destroy( (OBJEKT *)&lEntry );
    }
    DEREF( lEntries );
}

/*
 *      dAmberReadPrepFile
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Read an AMBER PREP file into a dictionary.
 */
DICTIONARY
dAmberReadPrepFile( char *sFilename )
{
DICTIONARY      dUnits;
UNIT            uUnit;
FILE            *fIn;
STRING          sTemp;

    fIn = FOPENCOMPLAIN( sFilename, "r" );
    if ( fIn == NULL ) 
        return(NULL);
    VP0(( "Loading Prep file: %s\n", GsBasicsFullName ));

    dUnits = dDictionaryCreate();

                /* Skip the two header lines */
    FGETS( sTemp, fIn );
    FGETS( sTemp, fIn );

                /* Read each UNIT in the PREP file */

    while ( !feof(fIn) && (uUnit = uAmberReadUnitFromPrep(fIn)) != NULL ) {
        DictionaryAdd( dUnits, sContainerName(uUnit), (GENP)uUnit );
    }
   
    fclose( fIn ); 
    return(dUnits);
}





/*
 *      psAmberReadParmSet
 *
 *      Read either an old AMBER parm set or a new frcmod parmset.
 *      Determine which type.
 */
PARMSET 
psAmberReadParmSet( FILE *fIn, char *sFilename )
{
BOOL            bFrcMod;
BOOL            bMass, bNonBonds;
PARMSET         psParms;

    bFrcMod = zbAmberDetermineParmSetFrcModType( fIn, &bMass, &bNonBonds );
    if ( !bFrcMod ) {
        psParms = zpsAmberReadParmSetOld( fIn );
    } else {
        /*
         *  frcmod
         */
        if ( bMass != bNonBonds ) {
                VP0(( "Modified force field files must contain %s\n",
                        "both a MASS and NONB entry, or neither" ));
                VP0(( "Could not load parameter set.\n" ));
                fclose(fIn);
                return(NULL);
        }
        VP0(( "Reading force field modification type file (frcmod)\n" ));
        psParms = zpsAmberReadParmSetFrcMod( fIn );
    }

    strcpy( sParmName(psParms), sFilename );
    return(psParms);
}





/*
 *      iAtomTypeHybridization
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Return the element type of an AMBER atom type.
 *      Return HUNKNOWN if it is undefined.
 */
int
iAtomTypeHybridization( char *sType )
{
int             i, iCount;
TYPEELEMENTt    *tP;

    iCount = iVarArrayElementCount( vaAtomTypes );
    if ( iCount == 0 )
        return( HUNKNOWN );

    tP = PVAI( vaAtomTypes, TYPEELEMENTt, 0 );
    for (i=0; i<iCount; i++, tP++) {
        if ( strcmp( sType, tP->sType ) == 0 ) 
            return(tP->iHybridization);
    }
    return(HUNKNOWN);
}




