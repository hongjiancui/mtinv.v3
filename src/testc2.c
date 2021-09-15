/*** this program is for prototype testing 
	1. Num Rec util matrix() memory allocation
	2. pthreads and fopen for reading and writting data to disk
***/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <pthread.h>

#include "../include/nrutil.h"     /** numerical recipes **/

typedef struct {
	long id;
	char description[128];
	int nt;
	float *x;
	int nrow, ncol;
	float **z;
} Data;

typedef struct {
	char filename[128];
	int i;
	Data c;
} ThreadData;

char progname[128];

int main( int ac, char **av )
{
	int i;
	int ndata;
	Data *a;
	Data *b;
	char filename[128];
	char outfilename[128];

/*** pthreads stuff ***/

	pthread_t *thread;
	ThreadData *td;
	pthread_attr_t attr;

	void *readData_parallel( void *ptr );

/*** functional prototypes ***/

	void readData( char *filename, Data *a );
	void makeData( Data *a, int id, char *description, int nt, int nrow, int ncol );
	void writeData( char *filename, Data *a );
	void printData( Data *a );
	int setpar( int, char **), mstpar(), getpar();
	void endpar();
	float *vector( int, int );
        float **matrix( int, int, int, int );

/*** start program, get command line args and parameters args ***/

	strcpy( progname, av[0] );
	setpar( ac, av );
	mstpar( "f", "s", filename );
	mstpar( "ndata", "d", &ndata );
	endpar();

	a = malloc( ndata * sizeof(Data) );

	for( i = 0; i < ndata; i++ )
	{
		makeData( &a[i], i, "This is a test structure", 100, 10, 10 );
		/* printData( &a[i]  ); */
	}

	for( i = 0; i < ndata; i++ )
	{
		sprintf( outfilename, "test.%02d.out", i );
		writeData( outfilename, &a[i] );
	}

	b = (Data *) malloc( ndata * sizeof(Data) );

/*** multithreaded reader ***/

	fprintf( stderr,"allocated memory for threads\n" );
        fflush(stderr);

	thread = (pthread_t *)malloc( ndata * sizeof(pthread_t) );
	td = malloc( ndata * sizeof(ThreadData) );


	/*
	for( i = 0; i < ndata; i++ )
		printData( &a[i] );
	*/

	for( i = 0; i < ndata; i++ )
	{
		sprintf( filename, "test.%02d.out", i );
	
		fprintf( stderr, "filename=%s i=%d\n", filename, i );
		fflush( stderr );

		strcpy( td[i].filename, filename );

		fprintf( stderr, "td.filename=%s\n", td[i].filename );
		fflush( stderr );
		td[i].i = i;
		td[i].c = b[i];
	}

	pthread_attr_init( &attr );
	/* pthread_attr_setdetachstate( &attr, PTHREAD_CREATE_JOINABLE ); */
	pthread_attr_setdetachstate( &attr, PTHREAD_CREATE_DETACHED );

	for( i = 0; i < ndata; i++ )
	{
		fprintf( stderr, "launching thread i=%d filename=%s\n", 
			i, 
			td[i].filename ); 
		fflush( stderr );

		pthread_create( &thread[i], &attr, readData_parallel, &td[i] );
	}

/************
	for( i = 0; i < ndata; i++ )
	{
		fprintf( stderr, "joining thread i=%d\n", i );
		fflush(stderr);
		pthread_join( thread[i], NULL );
	}
************/

	pthread_attr_destroy( &attr );

/************/
/****
	for( i = 0; i < ndata; i++ )
	{
		fprintf( stderr, "transfer i=%d\n", i );

		b[i].id          = td[i].c.id;

		strcpy( b[i].description, td[i].c.description );

		b[i].nt          = td[i].c.nt;

		fprintf( stdout, "%ld %s %d\n", b[i].id, b[i].description, b[i].nt );

		b[i].x           = vector( 0, b[i].nt );
		b[i].x           = (float *)td[i].c.x;
		b[i].nrow        = td[i].c.nrow;
		b[i].ncol        = td[i].c.ncol;
		b[i].z           = matrix( 0, b[i].nrow, 0, b[i].ncol );
        	b[i].z           = (float **)td[i].c.z;
	}
	
	for( i = 0; i < ndata; i++ )
	{
		sprintf( outfilename, "test.%02d.out", i );
		readData( outfilename, &b[i]  );
	}
	for( i = 0; i < ndata; i++ )
		printData( &b[i]  );
****/

}

void *readData_parallel( void *ptr )
{
	FILE *fp;
	char *filename;
	Data *a;
	ThreadData *td;
	int i;

	void printData( Data *a );

/*** start ***/
	for( i = 0 ; i < 100000000; i++ )
	{
		fprintf( stdout, "hello\n" );
	}

	a = malloc( sizeof(Data) );

	td = (ThreadData *) ptr;

	strcpy( filename, td->filename );
	fprintf( stderr, "thread %d opening file %s %s\n", i, filename, td->filename );
	fflush(stderr);

/***
	a->id   = (long)     td->c.id;
	a->nt   = (int)      td->c.nt;
	a->x    = (float *)  td->c.x;
	a->nrow = (int)      td->c.nrow;
	a->ncol = (int)      td->c.ncol;
	a->z    = (float **) td->c.z;
***/

        if( (fp = fopen( filename, "rb" )) == NULL )
	{
		fprintf( stderr, "error cannot read file %s\n",
			filename );
		exit(-1);
	}

        fread( &(a->id), sizeof(long), 1, fp );
        fread( a->description, 128*sizeof(char), 1, fp );
        fread( &(a->nt), sizeof(int), 1, fp );
        a->x = vector( 0, a->nt );
        fread( &(a->x[0]), a->nt*sizeof(float), 1, fp );
        fread( &(a->nrow), sizeof(int), 1, fp );
        fread( &(a->ncol), sizeof(int), 1, fp );
        a->z = matrix( 0, a->nrow, 0, a->ncol );
        fread( &(a->z), a->nrow*a->ncol*sizeof(float), 1, fp );
        fclose(fp);

/*
	fprintf( stderr, "done reading %s\n", filename );
	fflush(stderr);
*/
	printData( a );
	pthread_exit((void *)0);
}

void readData( char *filename, Data *a )
{
	FILE *fp;
	fp = fopen( filename, "rb" );
	fread( &(a->id), sizeof(long), 1, fp );
	fread( a->description, 128*sizeof(char), 1, fp );
	fread( &(a->nt), sizeof(int), 1, fp );
	a->x = vector( 0, a->nt );
	fread( &(a->x[0]), a->nt*sizeof(float), 1, fp );
	fread( &(a->nrow), sizeof(int), 1, fp );
	fread( &(a->ncol), sizeof(int), 1, fp );
	a->z = matrix( 0, a->nrow, 0, a->ncol );
	fread( &(a->z), a->nrow*a->ncol*sizeof(float), 1, fp );
	fclose(fp);
}

void writeData( char *filename, Data *a )
{
	FILE *fp;

	fp = fopen( filename, "wb" );
	fwrite( &(a->id), sizeof(long), 1, fp );
	fwrite( a->description, 128*sizeof(char), 1, fp );
	fwrite( &(a->nt), sizeof(int), 1, fp );
	fwrite( &(a->x[0]), a->nt*sizeof(float), 1, fp );
	fwrite( &(a->nrow), sizeof(int), 1, fp );
	fwrite( &(a->ncol), sizeof(int), 1, fp );
	fwrite( &(a->z), a->nrow*a->ncol*sizeof(float), 1, fp );
	fclose(fp);
}

void printData( Data *a )
{
	int i, j;
	fprintf( stdout, "id=%ld %s\n", a->id, a->description );
	fprintf( stdout, "nt=%d\n\tx = ", a->nt );

	for( i = 0; i < a->nt; i++ )
	{
		fprintf( stdout, " %g ", a->x[i] );
	}

	fprintf( stdout, "\n" );

	fprintf( stdout, "nrow = %d ncol = %d\n", a->ncol, a->nrow );

	for( i = 0; i < a->nrow; i++ )
	{
		for( j = 0; j < a->ncol; j++ )
		{
			fprintf( stdout, "%03.0f ", a->z[i][j] );
		}
		fprintf( stdout, "\n" );
	}
	fprintf( stdout, "\n" );
}

void makeData( Data *a, int ndata, char *description, int nt, int nrow, int ncol )
{
	int i, j;

	float *vector( int, int );
	float **matrix( int, int, int, int );

	a->id = (long)ndata + 10000000;
	strcpy( a->description, description );
	a->nt = nt;
	a->nrow = nrow;
	a->ncol = ncol;

	a->x = vector( 0, nt );

	for( i = 0; i < nt; i++ ) a->x[i] = (float)i;

	a->z = matrix( 0, nrow, 0, ncol );
	
	for( i = 0; i < nrow; i++ )
	{
		for( j = 0; j < ncol; j++ )
		{
			a->z[i][j] = (float)i*(float)j;
		}
	}
}
