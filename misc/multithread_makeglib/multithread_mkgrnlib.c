#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <string.h>
#include <unistd.h>
#include <signal.h>
#include <errno.h>
#include <sys/wait.h>

typedef struct {
	float dt;
	char modfilename[64];
	char staname[8];
	char network[8];
	pid_t pid;
} PARFILE;

#define MAXTHREADS 50  /*** this is the maximum number of threads (stations/rows) allowed to be launched at once ***/

char progname[128];

int main( int ac, char *av[] )
{
	char args[][64] = { "/Users/ichinose/Work/mtinv.v3.0.3/bin/mkgrnlib", "par=", "stnm=", "net=", "dt=", "\0" };
	char executable_pathname[256];

	/* char args[128][128]; */

	PARFILE pf[MAXTHREADS];
	int i, j, n;
	int retv, wait_status, wait_pid, last_pid;
	char parfilename[128];
	char executable_path[256];
	
	PARFILE *read_parfile( char *, PARFILE *, int * );

	int setpar(int, char **), mstpar(), getpar();
	void endpar();
	int verbose = 0;

	strcpy( progname, av[0] );

	setpar(ac,av);
	mstpar("parfile", "s", parfilename);
	mstpar("executable_pathname", "s", executable_pathname );
	getpar("verbose", "b", &verbose );
	endpar();

	fprintf(stderr, "%s: STDERR: start: parfile=%s executable_pathname=%s\n",
		progname, parfilename, executable_pathname );
	if(verbose)
	{
	  fprintf( stdout, "%s: STDOUT executable_pathname=%s parfilename=%s\n", 
		progname, executable_pathname, parfilename );
	  fflush(stdout);
	}

/*******************************************************/
/*** do a safety check for file existance            ***/
/*** and avoid some problems with runaway fork bombs ***/
/*******************************************************/

	if( access( executable_pathname, F_OK|X_OK ) != 0 ) 
	{
		fprintf( stderr, "%s: STDERR: cannot find executable_pathname=%s. Ensure it's set correctly\n",
			 progname, executable_pathname );
		fprintf( stdout, "%s: STDOUT: cannot find executable_pathname=%s. Ensure it's set correctly \n", 
                         progname, executable_pathname );
		exit(-1);
	}

	if( access( parfilename, F_OK ) != 0 ) 
	{
                fprintf( stderr, "%s: STDERR: cannot find parfile=%s. Ensure it's set correctly\n",
                         progname, parfilename );
                fprintf( stdout, "%s: STDOUT: cannot find parfile=%s. Ensure it's set correctly \n", 
                         progname, parfilename );
                exit(-1);
        }


/*** do a safety check for file existance ***/

	read_parfile( parfilename, pf, &n );

	if( n > MAXTHREADS )
	{
	  fprintf( stderr, "%s: STDERR: n = %d MAXTHREADS %d exceeded. EXIT!\n", progname, n, MAXTHREADS );
	  fprintf( stdout, "%s: STDOUT: n = %d MAXTHREADS %d exceeded. EXIT!\n", progname, n, MAXTHREADS );
	  fflush(stderr);
	  fflush(stdout);
	  exit(-1);
	}

	if(verbose)
	{
		fprintf( stdout, "nrows read = %d\n", n );
		fflush(stdout);
	}

	if(verbose)
	{
	  for( i = 0; i < n; i++ )
	  {
		strcpy(  args[0], executable_pathname );
                sprintf( args[1], "par=%s",  pf[i].modfilename );
                sprintf( args[2], "stnm=%s", pf[i].staname );
                sprintf( args[3], "net=%s",  pf[i].network );
                sprintf( args[4], "dt=%g",   pf[i].dt );

		for( j = 0; j <= 4; j++ )
			fprintf( stdout, "(%d)%s ", j, args[j] );

		fprintf( stdout, "\n" );
                fflush(stdout);
	  }
	}
	sleep(3);

/***********************************************************************************/
	for( i = 0; i < n; i++ )
	{
		pf[i].pid = fork();

		if( pf[i].pid < 0 )
		{
			fprintf(stdout, "fork error");
		}
		else if( pf[i].pid == 0 )
		{
			strcpy( args[0], executable_pathname );
			sprintf( args[1], "par=%s",  pf[i].modfilename );
			sprintf( args[2], "stnm=%s", pf[i].staname );
			sprintf( args[3], "net=%s",  pf[i].network );
			sprintf( args[4], "dt=%g",   pf[i].dt );

			retv = execlp( args[0], args[0], args[1], args[2], args[3], args[4], (char *)0 );

		} /* if fork pid */

	} /* loop over n */

	while( (wait_pid = waitpid( -1, NULL, 0 ))  ) 
	{
		if( errno == ECHILD )
		{
			break;
		}
	}
	if(verbose)fprintf( stdout, "%s: STDOUT finished.\n", progname );
	fprintf(stderr, "%s: STDERR: finished.\n\n", progname );
}

PARFILE *read_parfile( char *parfilename, PARFILE *pf, int *nrows )
{
	FILE *fp;
	int i = 0;
	char rec[256];
	extern char progname[128];

	if( (fp = fopen(parfilename, "r")) == NULL )
	{
		fprintf( stderr, "%s: cannot open parfilename=%s\n", 
			progname, parfilename );
		exit(-1);
	}
	
	while( fgets(rec,256,fp) != NULL )
	{
		/**** skip comment ****/
		if( rec[0] == '#' ) continue;

		sscanf(rec,"%s %s %s %f",
			pf[i].staname,
			pf[i].network,
			pf[i].modfilename,
			&(pf[i].dt) 
		);
		i++;
	}
	fclose(fp);
	*nrows = i;
	return (PARFILE *)pf;
}
