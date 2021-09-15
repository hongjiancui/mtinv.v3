#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif
#ifndef SWAP
#define SWAP 1
#endif
#ifndef NOSWAP
#define NOSWAP 0
#endif


int main(int ac, char **av)
{
	int isBigEndian(void);
	int TestByteOrder(void);
	int TestFileByteOrder(char *);
	int byteorder;

/*** one way ***/

	if( isBigEndian() )
	{
		fprintf( stdout, "I am Big Endian\n" );
	}
	else
	{
		fprintf( stdout, "I am Little Endian\n" );
	}

/*** another way ***/

	fprintf( stdout, "TestByteOrder=%d\n", TestByteOrder() );

	if( TestByteOrder() == LITTLE_ENDIAN )
		fprintf( stdout, "I am Little Endian\n" );

	if( TestByteOrder() ==  BIG_ENDIAN )
		fprintf( stdout, "I am Big Endian\n" );

	if( ac <= 1 ) exit(1);
	byteorder = TestFileByteOrder(av[1]);

	if( byteorder == SWAP )
		fprintf( stdout, "I need to be swapped\n" );
	else if( byteorder == NOSWAP )
		fprintf( stdout, "I do not need to be swapped\n" );
	else
		fprintf( stdout, "I cannot figure out endian-ness\n" );
}

int TestFileByteOrder(char *filename)
{
	FILE *fp;
	char cbuf4[4];

	union {
		char cval[4];
		float fval;
	} u1;

	union {
                char cval[4];
                float fval;
        } u2;

	fp = fopen(filename, "rb");
	fread( &cbuf4, 4, 1, fp );
	fclose(fp);

/*** NOSWAP ***/
	u1.cval[0] = cbuf4[0];
	u1.cval[1] = cbuf4[1];
	u1.cval[2] = cbuf4[2];
	u1.cval[3] = cbuf4[3];
	fprintf( stdout, "%g\n", u1.fval );

/*** SWAP ***/
	u2.cval[0] = cbuf4[3];
        u2.cval[1] = cbuf4[2];
        u2.cval[2] = cbuf4[1];
        u2.cval[3] = cbuf4[0];
	fprintf( stdout, "%g\n", u2.fval );

	if( fabs(u1.fval) > 1.0E+4 ) return SWAP;
	if( fabs(u2.fval) > 1.0E+4 ) return NOSWAP;
}

int TestByteOrder(void)
{
        int word = 0x0001;
        char *byte = (char *) &word;
        return(byte[0] ? LITTLE_ENDIAN : BIG_ENDIAN);
}

int isBigEndian(void)
{	
	char filename[64];
	char inBytes[4];
	int outInt = 1;
	FILE *fp;
	int bigendian;

/*******************************************/
/*** write number one to file as integer ***/
/*******************************************/
	sprintf( filename, "/tmp/IndianCheck_%s_%d.ThrowMeAway", 
		getenv( "USER" ),
		getpid() );

	if( (fp = fopen( filename, "wb" )) == NULL )
	{
		fprintf( stderr, "cannot open file %s\n", filename );
		exit(-1);
	}	
	fwrite( &outInt, 4, 1, fp );
	fclose(fp);

/*******************************************/
/*** read number one from file as integer ***/
/*******************************************/
	if( (fp=fopen(filename, "rb")) == NULL )
	{
		fprintf( stderr, "cannot open file %s\n", filename );
		exit(-1);
	}
	fread( inBytes, 1, 4, fp );
	fclose(fp);

/***********************/
/*** clean up        ***/
/***********************/
	remove(filename);

/***********************/
/*** check the value ***/
/***********************/
	if( inBytes[0] == 0 &&
	    inBytes[1] == 0 &&
	    inBytes[2] == 0 &&
	    inBytes[3] == 1 ) 
	{
		bigendian = TRUE;
	}
	else if( 
	    inBytes[0] == 1 &&
	    inBytes[1] == 0 &&
 	    inBytes[2] == 0 &&
	    inBytes[3] == 0 ) 
	{
		bigendian = FALSE;
	}
	else
	{
	    fprintf( stderr, "cannot test for endian-ness, Assumes bigendian\n" );
	    bigendian = TRUE;
	}

	return(bigendian);
}
