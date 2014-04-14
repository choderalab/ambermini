#include <stdio.h>
#include <stdlib.h>

/*
**  C routines for linking w/ ftn
*/

#ifdef CLINK_PLAIN
#define openbinpos_ openbinpos
#define readbinpos_ readbinpos
#define startbinpos_ startbinpos
#define writebinpos_ writebinpos
#endif
#ifdef CLINK_CAPS
#define openbinpos_ OPENBINPOS
#define readbinpos_ READBINPOS
#define startbinpos_ STARTBINPOS
#define writebinpos_ WRITEBINPOS
#endif

void
openbinpos_()
{
	char	magic[ 10 ];

	if( fread( magic, 1, 4, stdin ) != 4 ){
		fprintf( stderr, "Couldn't read magic number from BINPOS\n" );
		exit( -1 );
	}

	magic[ 4 ] = '\0';
	if( strcmp( magic, "fxyz" ) != 0 ){
		fprintf( stderr, "bad magic number \"%s\"\n", magic );
		exit( -1 );
	}
}

void
readbinpos_( int *n_atom, float apos[], int *eoflag )
{
	int	count;
	char buffer[4];
	*eoflag = 0 ;

	if( fread( buffer, sizeof( int ), 1, stdin ) != 1 ) {
		*eoflag = 1;
		return; 
	}
	if( 0==strncmp(buffer, "fxyz", 4)) { 
		fread( n_atom,sizeof( int ), 1, stdin );
	} else {
		n_atom = (int *)buffer;
	}
	if( ( count = fread( apos, sizeof( float ), 3 * *n_atom, stdin ) )
		!= 3 * *n_atom ){
		fprintf( stderr, "Could only read %d of %d atoms requested\n",
			count / 3, *n_atom );
		exit( -1 );
	}
}

void
startbinpos_()
{
#ifndef ABSOFT_WINDOWS
    if(isatty(1)) {
        /* don't send binary junk to terminal */
        fprintf(stderr,"please redirect output to file or pipe\n");
	}
#endif
    /* write magic number */
    fwrite("fxyz",4,1,stdout);
}

void
writebinpos_( int *n_atom, float apos[] )
{
	fwrite( n_atom, sizeof( int ), 1, stdout ) ;
	fwrite( apos, sizeof( float ), 3 * *n_atom, stdout );
}

