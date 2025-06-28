#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mafft.h"

int main( void )
{
	int i;
	int argc;
	char **argv;
	char **seq;
	char **name;
	char *message;
	int res;
	int n, l, mlen;

//	printf( "This is interface.\n" );


	l = 10000;	
	n = 13;
	seq = (char **)calloc( n, sizeof( char * ) );
	name = (char **)calloc( n, sizeof( char * ) );
	for( i=0; i<n; i++ ) seq[i] = calloc( l+1, sizeof( char ) );
	for( i=0; i<n; i++ ) name[i] = calloc( 100, sizeof( char ) );

	for( i=0; i<1; i++ )
	{
		strcpy( name[i*13+0], "name0" );
		strcpy( name[i*13+1], "name1" );
		strcpy( name[i*13+2], "name2" );
		strcpy( name[i*13+3], "name3" );
		strcpy( name[i*13+4], "name4" );
		strcpy( name[i*13+5], "name5" );
		strcpy( name[i*13+6], "name6" );
		strcpy( name[i*13+7], "name7" );
		strcpy( name[i*13+8], "name8" );
		strcpy( name[i*13+9], "name9" );
		strcpy( name[i*13+10], "name10" );
		strcpy( name[i*13+11], "name11" );
		strcpy( name[i*13+12], "name12" );

		strcpy( seq[i*13+0], "miyata");
		strcpy( seq[i*13+1], "takashi" );
		strcpy( seq[i*13+2], "miyata" );
		strcpy( seq[i*13+3], "takashi" );
		strcpy( seq[i*13+4], "miyata" );
		strcpy( seq[i*13+5], "miyata--" );
		strcpy( seq[i*13+6], "miyata--" );
		strcpy( seq[i*13+7], "miyata--" );
		strcpy( seq[i*13+8], "miyata--" );
		strcpy( seq[i*13+9], "miyata-" );
		strcpy( seq[i*13+10], "-" );
		strcpy( seq[i*13+11], "-" );
		strcpy( seq[i*13+12], "-" );

	}

	argc = 18;
	argv = (char **)calloc( argc, sizeof( char * ) );
	for( i=0; i<argc; i++ ) argv[i] = calloc( 100, sizeof( char ) );
	strcpy( argv[0], "disttbfast" );
	strcpy( argv[1], "-W" );
	strcpy( argv[2], "6" );
	strcpy( argv[3], "-b" );
	strcpy( argv[4], "62" );
	strcpy( argv[5], "-Q" );
	strcpy( argv[6], "100" );
	strcpy( argv[7], "-h" );
	strcpy( argv[8], "0" );
	strcpy( argv[9], "-F" );
	strcpy( argv[10], "-X" );
	strcpy( argv[11], "-s" );
	strcpy( argv[12], "0.0" );
	strcpy( argv[13], "-f" );
	strcpy( argv[14], "-1.53" );
	strcpy( argv[15], "-C" );
	strcpy( argv[16], "0" );
	strcpy( argv[17], "-P" );  // Necessary. DNA -> -D; Protein -> -P

	mlen = 5000;
	message = (char *)calloc( mlen+1, sizeof( char ) );

	fprintf( stderr, "first run\n" );
	res = disttbfast( n, l, mlen, name, seq, &message, argc, argv );
	fprintf( stderr, "second run\n" );
	res = disttbfast( n, l, mlen, name, seq, &message, argc, argv );
	fprintf( stderr, "third run\n" );
	res = disttbfast( n, l, mlen, name, seq, &message, argc, argv );

	fprintf( stderr, "\n\n\nmessage in interface = :%s:\n", message );
	free( message );

	if( res == GUI_LENGTHOVER )
	{
		fprintf( stderr, "length over!" );
	}
	else
	{
		fprintf( stderr, "res = %d\n", res );
		fprintf( stdout, "Output:\n" );
		for( i=0; i<n; i++ ) 
			fprintf( stdout, "%s\n", seq[i] );
	}
	fprintf( stderr, "argv = \n" );
	for( i=0; i<argc; i++ )
		fprintf( stderr, "%s ", argv[i] );
	fprintf( stderr, "\n" );
	
	for( i=0; i<n; i++ ) free( seq[i] );
	free( seq );
	for( i=0; i<n; i++ ) free( name[i] );
	free( name );
	for( i=0; i<argc; i++ ) free( argv[i] );
	free( argv );

}
