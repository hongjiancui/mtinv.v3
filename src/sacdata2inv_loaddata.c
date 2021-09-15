#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#include "../include/mt.h"         /** global datatype and structure declarations **/
#include "../include/nrutil.h"     /** numerical recipes **/

char progname[128];

#define MAXDATAPTS 1000000

void sacdata2inv_loaddata( 
	EventInfo *ev,
        int nsta,
        char *sacdatadir,
        char *respdir,
        char *currentdir,
	int verbose )
{

/*** local variable declarations ***/
	int    ista;
	int    ifiles, number_files, n3cmp, found;
	char   sacfiles[256][2048];
	int    errno;

/*** local function prototypes ***/
	
	/*** sacdata2inv_subs.c ***/
	void getsacfiles(
                char *stnm,
                char *net,
                char *pathname,
                char sacfiles[256][2048],
                int *nfiles,
                int verbose );

        /*** sacdata2inv_subs.c ***/
        float *getsacdata(
                char *cmp,
                float *x,
                Sac_Header *sp,
                char *sacfile,
                char *filename,
                int *ifound,
                int verbose );

        /*** getrespfile_sub.c ***/
        void  getrespfile(
                char *pathname,
                char *sta,
                char *net,
                char *cmp,
                char *khole,
                int verbose,
                char respfile[256] );

        /*** sacdata2inv_subs.c ***/
        void  fix_component_names( EventInfo *ev );

	/*** sacdata2inv_subs.c ***/
        void  writesacfile( EventInfo *ev );

	/*** transfer/sactransfer.c ***/
        void  transfer_response(
                float *data,
                int npts,
                float delta,
                char *sacpzfile,
                int verbose );

/*** start subroutine ***/

	for( ista = 0; ista < nsta; ista++ )
	{

/*********************************************************/
/*** given a directory path, get a file list of SAC    ***/
/*** files written from RDSEED                         ***/
/*********************************************************/
		if( verbose )
		{
       			fprintf( stdout,
			  "%s: sacdata2inv_loaddata(): calling getsacfiles()\n", 
				progname );
		}

		getsacfiles(
               		ev[ista].stnm,
               		ev[ista].net,
               		sacdatadir,
               		sacfiles,
               		&number_files,
               		verbose );

       		if( verbose )
       		{
                	fprintf(stdout,
                   	"%s: sacdata2inv_loadata(): number of SAC files=%d for station=%s network=%s\n",
                        	progname, number_files, ev[ista].stnm, ev[ista].net );

                	for( ifiles=0; ifiles<number_files; ifiles++ )
                	{
                        	fprintf( stdout, "%s: sacdata2inv_loadata(): ifiles=%d sacfile=%s\n",
                                	progname, ifiles, sacfiles[ifiles] );
                	}

        	}

/************************************************/
/*** go to the working directory              ***/
/************************************************/
        	if( (errno = chdir( sacdatadir )) != 0 )
        	{
                	fprintf( stderr, "%s: chdir errno=%d\n", progname, errno );
                	if( errno==-1 )
                	{
                        	fprintf(stderr, "%s: No Directory %s\n\n",
                                	progname, sacdatadir );
                	}
                	fprintf(stderr, "%s: Error in chdir %s\n",
                        	progname, sacdatadir );
        	}
        	if( verbose )
          		fprintf( stdout, "%s: sacdata2inv_loadata(): calling directory %s\n",
                		progname, sacdatadir );

/************************************************/
/*** load sac files into memory, figure out   ***/
/*** the component orientations               ***/
/************************************************/

        	if(verbose)
          		fprintf( stdout, "%s: sacdata2inv_loadata(): allocating memory to load data\n", progname );

        	ev[ista].z.data  = (float *)calloc(MAXDATAPTS,sizeof(float));
        	ev[ista].ns.data = (float *)calloc(MAXDATAPTS,sizeof(float));
        	ev[ista].ew.data = (float *)calloc(MAXDATAPTS,sizeof(float));

		n3cmp = 0;
		for( ifiles=0; ifiles < number_files; ifiles++ )
		{
			found = 0;
			if( n3cmp >= 3 ) break; /*** if more than 3 SAC files do not read the rest ***/
			ev[ista].ew.data = (float *)getsacdata( "EW",  ev[ista].ew.data, &(ev[ista].ew.s),
                                sacfiles[ifiles], ev[ista].ew.filename, &found, verbose );
			if(found)
                	{
                        	n3cmp++;
                        	continue;
                	}
			ev[ista].ns.data = (float *)getsacdata( "NS",  ev[ista].ns.data, &(ev[ista].ns.s),
                                        sacfiles[ifiles], ev[ista].ns.filename, &found, verbose );
			if(found)
                	{
                        	n3cmp++;
                        	continue;
                	}
			ev[ista].z.data  = (float *)getsacdata( "VER", ev[ista].z.data,  &(ev[ista].z.s),
                                        sacfiles[ifiles], ev[ista].z.filename, &found, verbose );
			if( found )
                	{
                        	n3cmp++;
                        	continue;
                	}

			if(verbose)
                	{
                        	fprintf(stdout,
                         	"%s: sacdata2inv_loadata(): done with getsacdata only n3cmp=%d will be read\n",
                                progname, n3cmp );
                	}

                	if( n3cmp < 3 )
                	{
                        	fprintf( stderr,
                          	"%s: STDERR sacdata2inv_loadata(): FATAL ERROR only n3cmp=%d read 3 needed. Quitting.\n",
                                progname, n3cmp );

                        	fprintf( stdout,
                          	"%s: STDOUT sacdata2inv_loadata(): FATAL ERROR only n3cmp=%d read 3 needed. Quitting.\n",
                                progname, n3cmp );

                        	exit(-1);
                	}

        	} /*** loop over ifiles ***/


/********************************************************************************/
/*** give a warning that more than 3 sac files exist for this station.network ***/
/*** the files are alphabetically sorted so typical the first 3 are data      ***/
/*** while the others are later arriving packets (see mergesac.c in misc)     ***/
/********************************************************************************/

        	if( n3cmp > 3 )
        	{
          	fprintf(stdout,
            	"%s: sacdata2inv_loadata(): warning more than 3 components found for station=%s.%s(ista=%d).\n",
                	progname, ev[ista].stnm, ev[ista].net, ista );
        	}


/*****************************************************************/
/*** fix broken kcmpnm  and khole character fields set by IRIS ***/
/*****************************************************************/

        	fix_component_names( &ev[ista] );

        	if(verbose)
        	{
                	fprintf( stdout,  "%s: sacdata2inv_loadata(): %s khole= ns=(%s) ew=(%s) z=(%s)\n",
                        	progname,
                        	ev[ista].z.filename,
                        	ev[ista].ns.s.khole,
                        	ev[ista].ew.s.khole,
                        	ev[ista].z.s.khole );
        	}

/********************************************************/
/*** find the SAC_PZs files in the Response directory ***/
/********************************************************/

		chdir( currentdir );

		if(verbose)
        	{
                	fprintf( stdout,
                  	"%s: sacdata2inv_loadata(): cwd=%s\n\t looking for SAC_PZs files in dir=%s\n",
                        progname, currentdir, respdir );
        	}

		sprintf( ev[ista].ew.sacpzfile, " " );
        	sprintf( ev[ista].ns.sacpzfile, " " );
        	sprintf( ev[ista].z.sacpzfile, " " );

		getrespfile(
                	respdir,
                	ev[ista].stnm,
                	ev[ista].net,
                	ev[ista].ew.s.kcmpnm,
                	ev[ista].ew.s.khole,
                	verbose,
                	ev[ista].ew.sacpzfile );
	
		getrespfile(
                	respdir,
                	ev[ista].stnm,
                	ev[ista].net,
                	ev[ista].ns.s.kcmpnm,
                	ev[ista].ns.s.khole,
                	verbose,
                	ev[ista].ns.sacpzfile );
	
		getrespfile(
                	respdir,
                	ev[ista].stnm,
                	ev[ista].net,
                	ev[ista].z.s.kcmpnm,
                	ev[ista].z.s.khole,
                	verbose,
                	ev[ista].z.sacpzfile );

		fprintf( stdout,
          	"%s: sacdata2inv_loadata(): Response files =\n \tew=(%s)\n \tns=(%s)\n \tz=(%s)\n",
                	progname,
                	ev[ista].ew.sacpzfile,
                	ev[ista].ns.sacpzfile,
                	ev[ista].z.sacpzfile );
	
	/*** reset the current directory for writting output files ***/

        	chdir( currentdir );

	} /*** loop over ista ***/

	fprintf( stdout,
		"%s: sacdata2inv_loadata(): done loading\n\n", progname );
}
