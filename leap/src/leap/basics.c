/*
 *      File: basics.c
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
 *              Declare global variables, and some routines which
 *              are always required.
 */



#include	"basics.h"
#include	"defaults.h"
#include	"zMatrix.h"
#include        "cmap.h"
defaultstruct GDefaults;

#include	<pwd.h>
#include	<signal.h>
#include	<stdarg.h>
#include	<stdlib.h>



           /* Globals for doing quick integer powers */
double  zzdS, zzdS1, zzdS6; 

STRING		GsBasicsFullName;

extern	double	strtod();


BOOL	GbInterrupt = FALSE;

int	GiUnitEditors = 0;

/*
 * ------------------------------------------------------------------
 *
 *      typedefs
 *
 *	For Fail-safe memory management.
 *
 */

#define FILENAMELEN     40
#define TRAILERLEN      100
                        /* MUST not be longer than TRAILERLEN+1 */
#define CHECKSTR        \
"usedUSEDThis is a buffer at the end of the object to check for overwrites"
                        /* MUST not be longer than TRAILERLEN+1 */
#define FREESTR         \
"FREEfreeNow this is free also and will be free until it is malloced again" 


typedef struct  MEMHEADERs {
	char                    sFile[FILENAMELEN];
	int                     iLineNumber;
	long                    lSize;
	struct MEMHEADERs	*mPNext;
	char                    sCheck[TRAILERLEN];
} MEMHEADER;

static  MEMHEADER	*SmPMallocList = NULL;

#include "varArray.h"
#define TSIZE	1000000
static VARARRAY	vtest = NULL;

void 
IMem()
{
	int   i, *ip;
	return;
	vtest = vaVarArrayCreate( sizeof(int) );
	VarArraySetSize( vtest, TSIZE );
	ip = PVAI( vtest, int, 0);
	for (i=0; i<TSIZE; i++) ip[i] = i;
}

void 
TMem()
{
	int	i, *ip;
	char	*p = NULL;
	return;
	ip = PVAI( vtest, int, 0);
	for (i=0; i<TSIZE; i++)
	if (ip[i] != i) {
		fprintf(stderr, "element %d=%d\n", i, ip[i]);
		*p = ' ';
	}
}


#ifdef DEBUG
/*
 *	zBasicsTrapBUS
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Trap SIGBUS exceptions.
 */
static void	
zBasicsTrapBUS( int iSignal, int iCode )
{
    DFATAL(( "Bus error.\n" ));
}

/*
 *	zBasicsTrapSEGV
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Trap SIGSEGV exceptions.
 */
static void	
zBasicsTrapSEGV( int iSignal, int iCode )
{
TMem();

    DFATAL(( "Segmentation violation.\n" ));
}
#endif




/*
 *      myAcos
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return the acos of a number, but don't crap
 *      out on domain errors.
 */
double  
myAcos( double d )
{
	if ( d >= 1.0 ) 
		return(0.0);
	if ( d <= -1.0 )
		return(PI);
	return(acos(d));
}



/*
 *	myPow
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Return the power of a number, but don't crap out
 *	on possible underflow errors.
 */
double	
myPow( double x, double y )
{
	if ( fabs(x) < VERYSMALL ) 
		return(0.0);
	return(pow(x,y));
}




/*
 *      iDoubleCompare
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Compare two double precision values.
 *
 *      The values are considered to be the same if
 *      the ranges    [dA-TOLERANCE,dA+TOLERANCE] contains dB.
 *
 *      Return 0 if they are the same, a value <0 if dA<dB and
 *      a value >0 if dA>dB.  Just like strcmp.
 */
int     
iDoubleCompare( double dA, double dB )
{
	if ( (dA-TOLERANCE<dB) && (dA+TOLERANCE>dB) )
		return(0);
	if ( dA < dB )
		return(-1);
	return(1);
}



/*
 *	bStringToDouble
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Return TRUE if the string passed to this routine is
 *	could be completely converted into a double precision
 *	value.  Only change the value of *dPData if the string
 *	completely represents a double precision value.
 */
BOOL	
bStringToDouble( char *cPData, double *dPData )
{
	char	*cPEnd;
	double	dValue;

	dValue = (double)strtod( cPData, &cPEnd );
	if ( cPEnd - cPData == strlen(cPData) ) {
		*dPData = dValue;
		return(TRUE);
	}
	return(FALSE);
}

/*
 *	bStringToInt - By David Rivkin
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Return TRUE if the string passed to this routine is
 *	could be completely converted into a integer
 *	value.  Only change the value of *dPData if the string
 *	completely represents a integer value.
 */
BOOL	
bStringToInt( char *cPData, int *iPData )
{
	char	*cp;

	for (cp=cPData; *cp; cp++) {
		if (*cp == '-') {
			if (cp != cPData) {
				VP0(( "'-' embedded in number %s\n", cPData ));
				return(FALSE);
			}
			continue;
		}
		if (!isdigit(*cp)) {
			VP0(( "non-digit in %s\n", cPData ));
			return(FALSE);
		}
		
	}
	*iPData = atoi( cPData );
	return(TRUE);
}



/*
 *      StringLower
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Convert the string to lower case.
 */
void    
StringLower( char *sStr )
{
    while ( *sStr != '\0' ) {
        *sStr = cLower((int)*sStr);
        sStr++;
    }
}






/*
 *--------------------------------------------------------------------
 *
 *      Memory management with consistcncy checking
 *
 */

STRING  GsLastMemOpFile;
int     GiLastMemOpLine;
int	GiMemoryAllocated;
BOOL	GbTestMemory = FALSE;

/*
 *      DebugCheckMemoryBlock
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Check if a memory block is a MALLOC'd memory block.
 */
static void    
DebugCheckMemoryBlock( MEMHEADER *mPMem )
{
    if ( strcmp( mPMem->sCheck, CHECKSTR ) != 0 ) {
        DFATAL(( "In FREE or REALLOC with a BAD memory block!\n" ));
    }
}



/*
 *      DebugMemoryTest
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Test the linked list of malloc'd memory blocks for
 *      consistancy.
 *      Display an error message and exit if the memory has been
 *      corrupted.
 *	If bReport is TRUE then write the memory header for each
 *	block as it is checked to the LOG file.
 */
static void
DebugMemoryTest( char *sFile, int iLine, BOOL bReport )
{
	MEMHEADER	*mPMem, *mPPreviousMem;
	char		*sTrailer;

	if ( !GbTestMemory ) 
		return;

	if ( bReport ) {
		VPLOG(( "Dumping unreleased memory\n" ));
	}

	MESSAGE(( "---TESTMEMORY\n" ));
    
	mPPreviousMem = NULL;
	mPMem = SmPMallocList;
	while ( mPMem != NULL ) {
		sTrailer = ((char*)mPMem)+sizeof(MEMHEADER)+mPMem->lSize;
		if ( strcmp( mPMem->sCheck, CHECKSTR ) != 0 ||
		     strcmp( mPMem->sCheck, sTrailer ) != 0 ) {
			PRINTF(( "MEMORY ERROR!!!!\n" ));
			PRINTF(( "A memory block has been corrupted!\n\n" ));
			PRINTF(( 
			  "Last memory operation before error in: %15s:%5d\n",
                          GsLastMemOpFile, GiLastMemOpLine ));
			PRINTF(( 
			  "Memory op. that discovered error in  : %15s:%5d\n",
			  sFile, iLine ));
			PRINTF(( "--------------------\n" ));
			PRINTF(( 
			  "Corrupted block malloc'd in: %s:%d   size: %ld\n",
                          mPMem->sFile, mPMem->iLineNumber,
                          mPMem->lSize ));
			PRINTF(( "Header= |%s|\n", mPMem->sCheck ));
			PRINTF(( "Trailer=|%s|\n", sTrailer ));
			PRINTF(( 
			  "NOTE: Header must be the same as Trailer!!!!!\n" ));
			if ( mPPreviousMem != NULL ) {
				PRINTF(( "------------------------\n" ));
				PRINTF(( "The preceding memory block is:\n" ));
				PRINTF(( 
				  "malloc'd in: %s  line: %d   size: %ld\n",
				  mPPreviousMem->sFile, 
				  mPPreviousMem->iLineNumber,
				  mPPreviousMem->lSize ));
				PRINTF(( "Header= |%s|\n",
				  mPPreviousMem->sCheck ));
			} else 
				PRINTF(( "This is the first memory block.\n" ));
			abort();
		} else {
			if ( bReport ) {
				VPLOG(( 
				  "Block malloc'd in: %s:%d   size: %ld\n",
				  mPMem->sFile, mPMem->iLineNumber,
				  mPMem->lSize ));
			}
		}
		mPPreviousMem = mPMem;
		mPMem = mPMem->mPNext;
	}
    
	MESSAGE(( "---TESTMEMORY GOOD!!!!!!!\n" ));
}




    
                        

/*
 *      DebugMalloc
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Malloc the amount of memory required with some extra
 *      information at the head and tail of the memory.
 *
 *      In the head place the name of the file and the line number
 *      where the memory was malloced, also place the memory in a 
 *      linked list.  Then place a 1 byte value just before the block of
 *      memory that will be returned to the user and the same value
 *      immediatly after.  This will be used to check if the block
 *      is being overwritten.
 */
char *
DebugMalloc( long lSize, char *sFile, int iLine )
{
MEMHEADER	*mPMem;
char		*sTrailer;
char		*cPBlock;

                /* Check for consistancy */

    MESSAGE(( "---MALLOC\n" ));

    DebugMemoryTest( sFile, iLine, FALSE );

                /* Malloc what the user wants and a header and trailer */
                
    mPMem = (MEMHEADER*)malloc( lSize + sizeof(MEMHEADER) + TRAILERLEN );
    FILEMESSAGE( "MEMORY", ( "In: %s:%d Malloc at: %lX,  %ld bytes\n", 
				sFile, iLine, mPMem, lSize ));
    if ( mPMem == NULL ) {
        DFATAL(( "Could not malloc: %s\n", strerror(errno) ));
    }
    
    cPBlock = ((char*)mPMem)+sizeof(MEMHEADER);

    if ( strcmp( mPMem->sCheck, "" ) != 0 ) {
        MESSAGE(( "---In MALLOC, previous contents=%s\n", mPMem->sCheck ));
    }
    
                /* Fill the header */

    strcpy( mPMem->sFile, sFile );
    mPMem->iLineNumber = iLine;
    mPMem->lSize = lSize;
    strcpy( mPMem->sCheck, CHECKSTR );
    mPMem->mPNext = SmPMallocList;
    SmPMallocList = mPMem;
    
    sTrailer = ((char*)mPMem)+sizeof(MEMHEADER)+lSize;
    memmove( sTrailer, mPMem->sCheck, TRAILERLEN );
    
    return ((char*)mPMem) + sizeof(MEMHEADER);
}





/*
 *      DebugRealloc
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Check ALL of memory, and if everything is OK then
 *      perform the realloc.  After the realloc is performed, update
 *      the trailer of the memory block.
 */
char *
DebugRealloc( char *cPBlock, long lSize, char *sFile, int iLine )
{
MEMHEADER	*mPMem, *mPPrevious, *mPCur;
char		*sTrailer;

                /* Check for consistancy */

    MESSAGE(( "---REALLOC\n" ));
    DebugMemoryTest( sFile, iLine, FALSE );

    
                /* Find where the actual memory block begins */

    mPMem = (MEMHEADER*)( cPBlock - sizeof(MEMHEADER) );
    FILEMESSAGE( "MEMORY", ( "<<<In %s:%d  Realloc at: %lX,  %ld bytes\n", 
			sFile, iLine, cPBlock, mPMem->lSize ));
    DebugCheckMemoryBlock(mPMem);

                /* Find the memory block in the linked list */
                /* to remove it, invalidate the current block */

    mPPrevious = NULL;
    mPCur = SmPMallocList;
    while ( mPCur != NULL && mPCur != mPMem ) {
        mPPrevious = mPCur;
        mPCur = mPCur->mPNext;
    }
    if ( mPCur == NULL ) {
        DFATAL(( "In DebugRealloc, MEMORY CORRUPTION!, block not found!" ));
    }
                /* Remove the block from the linked list */
    if ( mPPrevious == NULL ) {
        SmPMallocList = mPCur->mPNext;
    } else {
        mPPrevious->mPNext = mPCur->mPNext;
    }

                /* Allocate memory for the new block, copy the old stuff */
                /* into it, and free the old stuff */

    mPPrevious = mPMem;
    mPMem = (MEMHEADER*)malloc( lSize+sizeof(MEMHEADER)+TRAILERLEN );
    FILEMESSAGE( "MEMORY", ( "In %s:%d  >>>Realloc at: %lX,  %ld bytes\n",
				sFile, iLine, mPMem, lSize ));
    if ( mPMem == NULL ) {
        DFATAL(( "Could not malloc in REALLOC: %s\n", strerror(errno) ));
    }
    memmove( mPMem, mPPrevious, 
           MIN( lSize, mPPrevious->lSize) + sizeof(MEMHEADER)+TRAILERLEN );
                     
    strcpy( mPPrevious->sCheck, FREESTR );
    free(mPPrevious);
    mPMem->lSize = lSize;
    mPMem->mPNext = SmPMallocList;
    SmPMallocList = mPMem;

                /* Fill the trailer */
    
    sTrailer = ((char*)mPMem)+sizeof(MEMHEADER)+lSize;
    memmove( sTrailer, mPMem->sCheck, strlen(mPMem->sCheck) );
    sTrailer[strlen(mPMem->sCheck)]='\0';
    
    return ((char*)mPMem) + sizeof(MEMHEADER);
}






/*
 *      DebugFree
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Check ALL of memory, and if everything is OK then
 *      perform the free.
 */
void    
DebugFree( char *cPBlock, char *sFile, int iLine )
{
MEMHEADER	*mPMem, *mPCur, *mPPrevious;

                /* Check for consistancy */

    MESSAGE(( "---FREE\n" ));
    DebugMemoryTest( sFile, iLine, FALSE );


                /* Find where the actual memory block begins */

    mPMem = (MEMHEADER*)( cPBlock - sizeof(MEMHEADER) );
    FILEMESSAGE( "MEMORY", ( "In %s:%d Free at: %lX  %ld bytes\n", 
			sFile, iLine, cPBlock, mPMem->lSize));
    DebugCheckMemoryBlock(mPMem);

    strcpy( mPMem->sCheck, FREESTR );

                /* Find the memory block in the linked list */
                /* to remove it */

    mPPrevious = NULL;
    mPCur = SmPMallocList;
    while ( mPCur != NULL && mPCur != mPMem ) {
        mPPrevious = mPCur;
        mPCur = mPCur->mPNext;
    }
    if ( mPCur == NULL ) {
        DFATAL(( "In DebugFree, MEMORY CORRUPTION!, Memory block could not be found!!!!!" ));
    }
                /* Remove the block from the linked list */
    if ( mPPrevious == NULL ) {
        SmPMallocList = mPCur->mPNext;
    } else {
        mPPrevious->mPNext = mPCur->mPNext;
    }
    
                /* Perform the free */

    free( mPMem );
}




/*
 *---------------------------------------------------------------
 *
 *      Message printing.
 *
 *      Maintain and look up source files that currently display
 *      MESSAGE statements.
 */
 
static  int     SiMessageFiles = 0;
static  STRING  SsaMessageFiles[MAXMESSAGEFILES];



 

/*
 *      bMessageCheck
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Check if a message should be printed.
 */
BOOL    
bMessageCheck( char *sFile )
{
int     i;
    if ( SiMessageFiles == 0 ) 
	return FALSE;

    for ( i=0; i<SiMessageFiles; i++ ) {
        if ( SsaMessageFiles[i][0] == '*' ) 
		return TRUE;
        if ( strncmp( sFile, SsaMessageFiles[i], 
                        strlen(SsaMessageFiles[i]) ) == 0 ) 
		return TRUE;
    }
    return FALSE;
}


/*
 *      MessageAddFile
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Add a file to the list of files which should produce
 *      messages.
 */
void    
MessageAddFile( char *sFile )
{

    if ( SiMessageFiles >= MAXMESSAGEFILES ) 
	return;
    strcpy( SsaMessageFiles[SiMessageFiles], sFile );
    SiMessageFiles++;

    PRINTF(( "TURNING ON MESSAGES FROM: <%s>\n", sFile ));
}


/*
 *      MessageRemoveFile
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Remove a file from the list of files that should produce
 *      messages.
 */
void    
MessageRemoveFile( char *sFile )
{
int             i, j;

                /* Find the file */
    for ( i=0; i<SiMessageFiles; i++ ) {
        if ( strcmp( SsaMessageFiles[i], sFile )== 0 ) {
            for ( j=i;j<SiMessageFiles-1;j++ ) {
                strcpy( SsaMessageFiles[j], SsaMessageFiles[j+1] );
            }
            SiMessageFiles--;
            break;
        }
    }
    PRINTF(( "TURNING OFF MESSAGES FROM: <%s>\n", sFile ));
}


/*
 *      MessageFileList
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Print a list of the files which will generate messages.
 */
void    
MessageFileList()
{
	int	i;


    for ( i=0; i<SiMessageFiles; i++ ) {
        myPrintf( "%s\n", SsaMessageFiles[i] );
    }
    myPrintf( "------\n" );
}    



#ifdef DEBUG
/*
 *	MessageInitialize
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Initialize message printing from 
 *	those files mentioned in the MESSAGEON environment variable.
 *
 *	The MESSAGEON variable contains a sequence of names
 *	seperated by spaces.
 */
static void	
MessageInitialize()
{
STRING	sName;
char	*cPGet, *cPPut;

    PRINTF(( "Reading MESSAGEON environment variable.\n" ));
    cPGet = (char*)getenv("MESSAGEON");

    if ( cPGet == NULL ) 
	return;

    while ( *cPGet ) {
        cPPut = sName;
	while ( *cPGet && *cPGet != ' ' ) {
	    *cPPut = *cPGet;
	    cPPut++;
	    cPGet++;
	}
	*cPPut = '\0';
	if ( *cPGet ) cPGet++;
	MESSAGEON(sName);
    }
}
#endif



/*
 *=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
 *
 *	File handling routines.
 *
 */
#define	MAXDIRECTORIES	10

static STRING		SsaDirectories[MAXDIRECTORIES];
static int		SiCurDirectory = 0;

static int	
iExpandDir( char *sExpanded, char *sOriginal )
{
char		user[100];
int		i;
struct passwd 	*pw, *getpwnam();

/* 
TODO: add string size protection
*/
    if ( sOriginal[0] != '~' ) {
	/*
	**	nothing to expand
	*/
    	strcpy( sExpanded, sOriginal );
    	return(0);
    }

    /*
     *  there is a leading ~ which needs to be expanded
     */

    if ( sOriginal[1] == '\0'  ||  sOriginal[1] == '/' ) {

	/*
	 *  relative to user's home dir
	 */

	char	*home;

	home = (char*)getenv("HOME");
	if ( home == (char*)NULL ) {
		VP0(( "~: Could not get HOME from environment\n" ));
		return(1);
	}
	strcpy( sExpanded, home);
	if ( sOriginal[1] == '/' )
		strcat( sExpanded, &sOriginal[1] );
	return(0);
    }

    /*
     *  relative to a specified user's home dir - get user name &
     *	look it up in pwd file
     */

    for (i=0; i<99; i++) {
	if ( sOriginal[i+1] == '\0' || sOriginal[i+1] == '/' ) {
		user[i] = '\0';
		break;
	}
	user[i] = sOriginal[i+1];
    }
    i++;
    pw = getpwnam( user );
    if ( pw == NULL ) {
	VP0(( "%s: Could not get from password file: %s\n", user,
						strerror(errno) ));
	return(1);
    }

    strcpy( sExpanded, pw->pw_dir );
    if ( sOriginal[i] == '/' )
	strcat( sExpanded, &sOriginal[i] );
	
    return(0);
}

/*
 *	BasicsAddDirectory
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Add a directory to the list of directories to search
 *	for a file when it cannot be opened in the current directory.
 */
int		
BasicsAddDirectory( STRING sDirectory, int bomb )
{
STRING		sTmpDir;
FILESTATUSt	fsStatus;

    if ( SiCurDirectory + 1 >= MAXDIRECTORIES ) {
	VP0(( "Limit reached on include directories - can't add dir\n" ));
	if ( bomb )
		exit(1);
	return(0);
    }
    if ( iExpandDir( sTmpDir, sDirectory ) ) {
	if ( bomb )
		exit(1);
	return(0);
    }
    fsStatus = fsSysdependFileStatus(sTmpDir);
    if ( !(fsStatus.fMode & FILEDIRECTORY) ) {
	if (errno)
		VP0(( "%s: %s\n", sTmpDir, strerror(errno) ));
	else
		VP0(( "%s: not a directory %s\n", sTmpDir ));
	if ( bomb )
		exit(1);
	return(0);
    }
    strcpy( SsaDirectories[SiCurDirectory], sTmpDir );
    SiCurDirectory++;
    return(1);
}







/*
 *	fBasicsMyFopen
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	This routine is a wrapper to the library fopen routine.
 *	If the path is absolute (starts w/ root or '.'), try it. 
 *	If not, try the current directory then search through the 
 *	list of directories in SvaDirectories.
 */
FILE *
fBasicsMyFopen( char *sFilename, char *sAttributes, BOOL bComplain )
{
FILE		*fFile;
int		i, iExistErr;
FILESTATUSt	fsStatus;
STRING		sExpanded;

    if ( strlen(sFilename) == 0 ) 
	return(NULL);

    if ( strlen(sFilename) == 1 && sFilename[0] == '-' ) {
	VP2(( "Reading from standard input.\n" ));
	return(stdin);
    }

    if ( iExpandDir( sExpanded, sFilename ) )
	return(NULL);
    fsStatus = fsSysdependFileStatus(sExpanded);

    if ( fsStatus.fMode & FILEDIRECTORY ) {
	VP0(( "%s is a directory\n", sExpanded ));
	return(NULL);
    }

    if ( sExpanded[0] == '/' || sExpanded[0] == '.' ) {
	/* (this is where an sExpanded ~ dir would fall) */
        strcpy ( GsBasicsFullName, sExpanded );
        fFile = fopen( sExpanded, sAttributes );
	if ( fFile == NULL  &&  bComplain ) {
		VP0(( "Could not open file %s: %s\n", 
						sExpanded, strerror(errno) ));
	}
	return(fFile);
    }
    sprintf ( GsBasicsFullName, "./%s", sExpanded );

    fFile = fopen( sExpanded, sAttributes );
    if ( fFile != NULL ) 
	return(fFile);
    iExistErr = 1;
    for ( i=0; i<SiCurDirectory; i++ ) {
	strcpy( GsBasicsFullName, SsaDirectories[i] );
	strcat( GsBasicsFullName, "/" );
	strcat( GsBasicsFullName, sExpanded );
	fFile = fopen( GsBasicsFullName, sAttributes );
	if ( fFile != NULL ) 
	    return(fFile);
	if ( errno != ENOENT ) {
	    VP0(( "Opening %s: %s\n", GsBasicsFullName, strerror(errno) ));
	    iExistErr = 0;
	}
    }
    if ( bComplain ) {
	VP0(( "Could not open file %s: %s\n", sExpanded,
			( iExistErr ? "not found" : "system error" )  ));
    }
    return(NULL);
}




/*
 *------------------------------------------------------------------
 *
 *	Handle multiple output sinks.
 *
 *	The current output sink is the one where all output from
 *	the VPx commands is sent.
 *
 */


typedef	struct	{
	BOOL		bSinkUsed;
	BOOL		bPrintPrefix;
	STRING		sPrefix;
	VFUNCTION	fCallback;
	GENP		PData;
} SINKINFOt;

static	int		SiNumberOfSinks = 0;
static	SINKINFOt	*SsiPSinks = NULL;

#define	MAX_SINK_STACK		30

static int		SiaSinkStack[MAX_SINK_STACK];
static int		SiNextSink = 0;
static int		iCurrentPrintSink = -1;

/*
 *	iCreatePrintSink
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Create a new print sink and return an integer
 *	handle for the sink that can be used by subsequent
 *	calls to DefineCurrentPrintSink or DestroyPrintSink.
 *
 *	The PData field can be used to pass additional info to
 *	the output callback, like a Widget.
 *
 *	The Callback function has the following form.
 *
 *	callback( cPString, PData )
 *
 */
int
iCreatePrintSink( VFUNCTION fOutputCallback, char *sPrefix, GENP PData )
{
int		i;
BOOL		bFoundOne;

	/* First check if there isn't a free sink */

    bFoundOne = FALSE;
    for ( i=0; i<SiNumberOfSinks; i++ ) {
	if ( !SsiPSinks[i].bSinkUsed ) {
	    bFoundOne = TRUE;
	    break;
	}
    }

		/* If one is not found then allocate space for one */

    if ( !bFoundOne ) {
	if ( SsiPSinks == NULL ) {
	    i = 0;
	    MALLOC( SsiPSinks, SINKINFOt*, sizeof(SINKINFOt) );
	    SiNumberOfSinks = 1;
	} else {
	    i = SiNumberOfSinks;
	    REALLOC( SsiPSinks, SINKINFOt*, SsiPSinks,
				sizeof(SINKINFOt)*(SiNumberOfSinks+1) );
	    SiNumberOfSinks++;
	    /*
	     *  realloc can mean that globals pointing to the
	     *	current print sink now point to freed memory
	     */
	    if (iCurrentPrintSink == -1)
		DFATAL(("iCurrentPrintSink == -1\n"));
	    GcPPrefix = SsiPSinks[iCurrentPrintSink].sPrefix;
	    GbPrintPrefix = SsiPSinks[iCurrentPrintSink].bPrintPrefix;
	    GfPrintStringCallback = SsiPSinks[iCurrentPrintSink].fCallback;
	    GPData = SsiPSinks[iCurrentPrintSink].PData;
	}
    }

		/* Define the information */

    SsiPSinks[i].bSinkUsed = TRUE;
    SsiPSinks[i].bPrintPrefix = TRUE;
    SsiPSinks[i].fCallback = fOutputCallback;
    strcpy( SsiPSinks[i].sPrefix, sPrefix );
    SsiPSinks[i].PData = PData;
    MESSAGE(( "Created print sink: %d\n", i ));

    return(i);
}



GENP	GPData;

/*
 *	zDefineCurrentPrintSink
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Using the handle passed in iHandle, define the
 *	current print sink.
 */
static void
zDefineCurrentPrintSink( int iHandle )
{
    iCurrentPrintSink = iHandle;

    if ( iHandle >= SiNumberOfSinks || iHandle < 0 ) {
	DFATAL(( "Illegal print sink handle: %d, must be in [0..%d]\n",
			iHandle, SiNumberOfSinks-1 ));
    }

    MESSAGE(( "Selecting print sink: %d\n", iHandle ));

    if ( SsiPSinks[iHandle].bSinkUsed ) {
	GcPPrefix = SsiPSinks[iHandle].sPrefix;
	GbPrintPrefix = SsiPSinks[iHandle].bPrintPrefix;
	GfPrintStringCallback = SsiPSinks[iHandle].fCallback;
	GPData = SsiPSinks[iHandle].PData;
    } else {
	DFATAL(( "Unused print sink: %d\n", iHandle ));
    }
}


/*
 *	DestroyPrintSink
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Destroy the print sink indicated by iHandle.
 */
void	
DestroyPrintSink( int iHandle )
{
    if ( iHandle >= SiNumberOfSinks || iHandle < 0 ) {
	DFATAL(( "Illegal print sink handle: %d, must be in [0..%d]\n",
			iHandle, SiNumberOfSinks-1 ));
    }

    MESSAGE(( "Destroying print sink: %d\n", iHandle ));

    if ( SsiPSinks[iHandle].bSinkUsed ) {
	SsiPSinks[iHandle].bSinkUsed = FALSE;
    } else {
	DFATAL(( "Unused print sink: %d\n", iHandle ));
    }
}



/*
 *	PushCurrentPrintSink
 *
 *	Push the sink in (iSink) onto the stack as the
 *	current print sink.
 */
void
PushCurrentPrintSink( int iSink )
{
    if ( SiNextSink >= MAX_SINK_STACK ) {
	DFATAL(( "Exhausted Print Sink Stack\n" ));
    }
    SiaSinkStack[SiNextSink] = iSink;
    SiNextSink++;
    zDefineCurrentPrintSink(iSink);
}



/*
 *	PopCurrentPrintSink
 *
 *	Pop the top print sink from the stack.
 */
void
PopCurrentPrintSink()
{
	if ( SiNextSink <= 0 ) {
		DFATAL(( "Underflow in Print Sink Stack\n" ));
	}
	SiNextSink--;
	if ( SiNextSink == 0 ) {
	    GfPrintStringCallback = NULL;
	} else {
	    zDefineCurrentPrintSink(SiaSinkStack[SiNextSink-1]);
	}
}





/*
 *------------------------------------------------------------------
 *
 *      My own printf which checks verbosity levels
 *      and writes output to an optional log file.
 */

	/* Define variables used by myPrintf which maintain what */
	/* function to call to print a string on the display */
	/* and a buffer that stores the output to be printed */
	/* By default myPrintf will write its non-log output to */
	/* stdout using 'puts'.  But this can be re-directed using */
	/* DefinePrintStringCallback(func(char*)). */
	/* This is useful for output to Widgets */


FILE	*GfLog = NULL;
int     GiTraceIndentationLevel = 0;
int     GiVerbosityLevel = 0;
int     GiVerbosity;                    /* This changes for every P# */
BOOL	GbPrintPrefix = TRUE;
char	*GcPPrefix = NULL;


/*
 *	myPuts
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Put a string to the stdout without appending a newline to the end.
 *	PData is not used.
 */
static void	
myPuts( char *sLine, GENP PData )
{
    printf( "%s", sLine );
}

void	(*GfPrintStringCallback)() = myPuts;


#define	MAXCHARSPERPRINTF	5000		/* 5000 characters max */
static	char	SsTempBuffer[MAXCHARSPERPRINTF];

/*
 *	myPrintString
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Print the string to the appropriate devices.
 *	Break up the line at \n characters.
 *	Prefix the string with GcPPrefix if the last character
 *	printed was a \n.
 */
int itest=0;
void	
myPrintString( char *cPString )
{
char	*cPStart, *cPStop, *cPPrint;
BOOL	bPrintPrefix;

    cPStart = cPString;
    cPStop  = cPString;

    while ( *cPStop ) {

			/* Seperate out the \n, \0 terminated pieces */
			/* and print them possibly prefixed with GcPPrefix */
	/*
	 *  set cPStop to next newline or end of string
	 */
	while ( *cPStop ) {
	    if ( *cPStop == '\n' ) break;
	    cPStop++;
	}

	if ( *cPStop == '\n' ) {
	    /*
	     *  will print this line: copy into buffer to null-terminate
	     */
	    strncpy( SsTempBuffer, cPStart, (cPStop-cPStart)+1 );
	    SsTempBuffer[(cPStop-cPStart)+1] = '\0';
	    cPPrint = SsTempBuffer;
	    bPrintPrefix = TRUE;
	    /*
	     *  prepare for next line
	     */
	    cPStart = ++cPStop;
	} else {
	    /*
	     *  end of string: print as-is
	     */
	    cPPrint = cPStart;
	    bPrintPrefix = FALSE;
	}

	if ( GbPrintPrefix  &&  GcPPrefix != NULL ) {

	    if ( GfPrintStringCallback == NULL ) {
		DFATAL(( "Illegal print string callback!" ));
	    }
	    if ( GiVerbosityLevel >= GiVerbosity ) {
		GfPrintStringCallback( GcPPrefix, GPData );
	    }
	    if ( GfLog != NULL  &&  GiVerbosity >= 0 ) {
		fputs( GcPPrefix, GfLog );
	    }
	}
	GbPrintPrefix = bPrintPrefix;
	if ( GiVerbosityLevel >= GiVerbosity ) {
	    GfPrintStringCallback( cPPrint, GPData );
	}
	if ( (GfLog != NULL) && ( GiVerbosity >= 0 ) ) {
	    fputs( cPPrint, GfLog );
	}
    }
}




/*
 *	BasicsInitialize
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Initialize basic routines.
 */
void	
BasicsInitialize()
{
	int i;
	(void)zMatrixInit();

/*   initialize some of the defaults   */

	GDefaults.pdbwritecharges = 0;
	GDefaults.dGridSpace = 1.0;
	GDefaults.dShellExtent = 4.0;
	GDefaults.dDipoleDampFactor = 0.0;
	GDefaults.dSceeScaleFactor = 1.2;
	GDefaults.dScnbScaleFactor = 2.0;
#define DIEL_R2		1
	GDefaults.iDielectricFlag = DIEL_R2;
	GDefaults.iGBparm = 2;
	GDefaults.iOldPrmtopFormat = 0;
	GDefaults.iGibbs = 0;
	GDefaults.iCharmm = 0;
	GDefaults.iResidueImpropers = 0;
	GDefaults.iDeleteExtraPointAngles = 1;
	GDefaults.iPdbLoadSequential = 1;
	GDefaults.iHaveResIds = 0;
	GDefaults.iUseResIds = 0;
	GDefaults.iFlexibleWater = 0;
	for( i=0; i<MAXRESID; i++ ){
		GDefaults.sResidueId[i][0] = '\0';
	}
	GDefaults.iCMAP = 0;
	GDefaults.iIPOL = 0;
	GDefaults.iIPOLset = 0;
	GDefaults.nocenter = 0;

/*
    signal( SIGSEGV, zBasicsTrapSEGV );
    signal( SIGBUS, zBasicsTrapBUS );
*/
#ifdef	DEBUG
    MessageInitialize();
#endif

#ifdef DEBUG
    ALWAYS(( "Memory testing on = %s\n", sBOOL(bTEST_MEMORY_ON()) ));
#endif

} 

/*
 *	BasicsKlassMisMatchPanic
 *
 *	A Klass mismatch has occured so stop.
 */
void	
BasicsKlassMismatchPanic( GENP PObj, char *cPFile, int iLine )
{
    printf( "ERROR: Klass mismatch in file: %s  line: %d\n",
		cPFile, iLine );
    if ( !PObj ) {
	printf( "Object is NULL!\n" );
    } else {
	printf( "Object is not NULL but of unknown type.\n" );
    }

    abort();
}
