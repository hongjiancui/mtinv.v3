#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <sys/stat.h>
#include <errno.h>

#include "../include/sacfile.h"
#include "../include/mt_version.h" /*** get the version label ***/

char progname[128];

int main( int ac, char **av )
{	
	Sac_Header *s;
	float *data;
	char sacfilename[256];
	FILE *fp;
	MyTime ot;
	char timestring[25];
	char comment[256];
	char loc[8];
	float depth = 0, evla=-999, evlo=-999;
	double drdistkm, daz, dbaz;

	int verbose = 0;
	int igmt5 = 0;
	int i = 1, nsta = 0;
	int irt = 0;

/*** sacio.c ***/

	float *readsac( Sac_Header *s, char *filename, int verbose );

	int distaz( double olat, double olon, double tlat, double tlon, 
		double *rdistkm, double *az, double *baz );

/*** local ***/

	void make_gmt5_map( Sac_Header *s, int nsta, char *script_filename, int irt, int verbose );

	void make_gmt4_map( Sac_Header *s, int nsta, char *script_filename, int irt, int verbose );

	void make_glib_parfile( MyTime *ot, Sac_Header *s, int nsta, float depth,
		char *comment, int igmt5, int irt, char *script_filename, int verbose );

	void Usage( int ac, char **av );

/*** timesubs.o ***/
        void parsestring( MyTime *ot, char *str );
	void WriteMyTime2STDERR( MyTime *t );

/*** shorten_path.c ***/

        char *shorten_path( char *pathname, char *filename );
	char pathname[128];

/*** start progname ***/
	strcpy( pathname, av[0] );
        shorten_path( pathname, progname );

	if( verbose )
        {
          fprintf( stdout, "%s: STDOUT: Version=%s ReleaseDate=%s exec full path=%s\n",
                progname, Version_Label, Version_Date, pathname );
        }

        fprintf( stderr, "%s: STDERR: Version=%s ReleaseDate=%s exec full path=%s\n",
                progname, Version_Label, Version_Date, pathname );

	sprintf( timestring, "2001/01/01,01:01:59.99" );
	sprintf( comment, "Comment Here" );

	if( ac <= 1 )
	{
	  fprintf( stderr,
	  "%s: Usage %s -z [depth km] -ev [lat lon] -ot [\"YYYY/MM/DD,hh:mm:ss.ss\"] -com [\"Comment\"] -gmt5 -help -verbose [sac files]\n", 
		progname, progname );
	  exit(-1);
	}

	s = (Sac_Header *)malloc( ac * sizeof(Sac_Header) );

/***  go through list of sac files and read headers for station information ***/
/** i = 0, skip because is av[0]=progname **/

	nsta = 0;
	while( i < ac )
	{
		strcpy( sacfilename, av[i] );

		if( strcmp( sacfilename, "-h" ) == 0 || strcmp( sacfilename, "-help" ) == 0 ) 
		{
			Usage(ac,av);
		}

		if( strcmp( sacfilename, "-realtime" ) == 0 || strcmp( sacfilename, "-rt" ) == 0 )
		{
			fprintf( stderr, "%s: -realtime or -rt stop user interactive map display\n", progname );

			irt = 1;
			i++;
			continue;
			
		} 

		if( strcmp( sacfilename, "-ev" ) == 0 || strcmp( sacfilename, "-event" ) == 0 )
		{
			i++;
			evla = atof( av[i] );
			i++;
			evlo = atof( av[i] );
			fprintf( stderr, "%s: -ev or -event : evla=%g evlo=%g\n", progname, evla, evlo );
			i++;
			continue;
		}

		if( strcmp( sacfilename, "-z" ) == 0 || strcmp( sacfilename, "-depth" ) == 0  )
		{
			i++;
			depth = atof( av[i] );
			fprintf( stderr, "%s: -z or -depth = %g\n", progname, depth );
			i++;
			continue;
		}

		if( strcmp( sacfilename, "-com" ) == 0 || strcmp( sacfilename, "-comment" ) == 0  )
		{
			i++;
			strcpy( comment, av[i] );
			fprintf( stderr, "%s: -com or -comment : %s\n", progname, comment );
			i++;
			continue;
		}

		if( strcmp( sacfilename, "-ot" ) == 0 || strcmp( sacfilename, "-origintime" ) == 0  ) 
		{
			i++;
			strcpy( timestring, av[i] );
			parsestring( &ot, timestring );
			fprintf( stderr, "%s: -ot or -origintime : ", progname ); 
				WriteMyTime2STDERR( &ot );
			i++;
			continue;
		}

		if( strcmp( sacfilename, "-verbose" ) == 0 || strcmp( sacfilename, "-v" ) == 0  )
		{
			verbose = 1;
			fprintf( stderr, "%s: Switch %s found for verbosity, continuing...\n", 
                                progname, sacfilename );
			i++;
                        continue;
		}

		if( strcmp( sacfilename, "-gmt5" ) == 0 ) 
		{
			fprintf( stderr, "%s: Switch %s found for GMT version 5.x.x, continuing...\n", 
				progname, sacfilename );
			igmt5 = 1;
			i++;
			continue;
		}

		if( (fp=fopen(sacfilename,"rb")) == NULL )
		{
			fprintf( stderr, "%s: cannot fopen file %s, continuing...\n", 
				progname, sacfilename );
			i++;
			continue;
		}
		else
		{
			fclose(fp);
			if(verbose)
			{
			  fprintf( stderr, "%s: i=%d nsta=%d %s\n",
				progname, i, nsta, sacfilename );
			}
			if( (data = readsac( &s[nsta++], sacfilename, verbose )) == NULL )
			{
				nsta--;
				fprintf( stderr, "%s: cannot fopen file %s, continuing...\n",
                                	progname, sacfilename );
				i++;
				continue;
			}
                	i++;
		}
        }

	if( nsta > 0 )
	{
		fprintf( stderr, "%s: processing %d number of stations.\n", progname, nsta );
	}
	else
	{
		fprintf( stderr, "%s: no stations found. EXITING.\n", progname );
		exit(0);
	}

	if( evla != -999 && evlo != -999 )
	{
		for( i = 0; i < nsta; i++ )
		{
			s[i].evla = evla;
			s[i].evlo = evlo;
			distaz( (double)evla, (double)evlo, (double)s[i].stla, (double)s[i].stlo,
				&drdistkm, &daz, &dbaz );
			s[i].dist = (float)drdistkm;
			s[i].az   = (float)daz;
			s[i].baz  = (float)dbaz;
		}
	}

/*** fix khole problem with downloading from webservice***/

	for( i = 0; i < nsta; i++ )
	{
		/* fprintf( stdout, "khole=(%s)\n", s[i].khole ); */
		if( strncmp( s[i].khole, "-12345", 6 ) == 0 ) 
		{
			strcpy( s[i].khole, "" );
		}
		/* fprintf( stdout, "khole=(%s)\n", s[i].khole ); */
	}

	if( igmt5 )
	{
		if(verbose) fprintf(stderr, "%s: make gmt5 map\n", progname );
		make_gmt5_map( s, nsta, "plotmap.csh", irt, verbose );
	}
	else
	{
		if(verbose)fprintf(stderr, "%s: make gmt4 map\n", progname );
		make_gmt4_map( s, nsta, "plotmap.csh", irt, verbose );
	}

	if(verbose)
	{
		fprintf( stderr, "%s: main(): calling make_glib_parfile(): \n", progname );
	}
	make_glib_parfile( &ot, s, nsta, depth, comment, igmt5, irt, "makeglib.csh", verbose );

	exit(0);
}

void make_glib_parfile( 
	MyTime *ot, 
	Sac_Header *s, 
	int nsta, 
	float depth,
	char *comment, 
	int igmt5,
	int irt,
	char *script_filename, 
	int verbose )
{
	FILE *fp;
	int i = 0, j = 0;
	float dt = 0.35;
	char model_name[128];
	char gmt_version_string[12];
	char realtime_string[12];
	char com[3];

/*** sort by distance ***/
	int *indx;
	float *dist;
	void sort_sac_headers_by_dist( int n, Sac_Header *s, float *arrin, int *indx );

/************************/
/*** start subroutine ***/
/************************/

	if(verbose)
	{
	  fprintf( stderr,
		"%s: make_glib_parfile(): calling sort_sac_headers_by_dist(): nsta=%d\n",
			progname, nsta );
	}

	dist = calloc( nsta+1, sizeof(float) );
	indx = calloc( nsta+1, sizeof(int) );

	sort_sac_headers_by_dist( nsta, s, dist, indx );

	if(verbose)
	  fprintf( stderr, "%s: make_glib_parfile(): done with sort_sac_headers_by_dist():\n", progname );

	sprintf( model_name, "wus" );

	if( (fp=fopen( script_filename,"w" )) == NULL )
        {               
                fprintf( stderr, "%s: cannot open file %s for writting\n",
                        progname, script_filename );
                exit(-1);
        }
	
	fprintf( fp, "#!/bin/csh\n" );
	fprintf( fp, "###\n" );
	fprintf( fp, "### SET the executable path \n" );
	fprintf( fp, "### set MTINV_PATH=/Users/ichinose1/Work/mtinv.v3.0.5\n" );
	fprintf( fp, "###\n" );
	fprintf( fp, "set MTINV_PATH = ( ` find ${HOME} -name mtinv.v%s -print` ) \n", Version_Label );
	fprintf( fp, "echo \" MTINV_PATH = ${MTINV_PATH} \" \n" );
	fprintf( fp, "\n" );

	fprintf( fp, "### \n" );
	fprintf( fp, "### %4d/%02d/%02dT%02dh%02dm%05.2f %g %g %s\n",
		ot->year, ot->month, ot->mday, ot->hour, ot->min, ot->fsec,
		s[0].evla, s[0].evlo, comment );
	fprintf( fp, "### \n" );
	
	fprintf( fp, "\n" );
	fprintf( fp, "cat >! %s.par << EOF\n", model_name );
	fprintf( fp, "velmod=%s\n", model_name );

	if( depth >= 0 && depth < 33 )
	{
		fprintf( fp, "zrange=3,3,33\n" );
	}
	else if( depth >= 33 && depth < 100 )
	{
		fprintf( fp, "zrange=10,10,100\n" );
	}
	else if( depth >= 100 && depth < 300 )
	{
		fprintf( fp, "zrange=10,20,300\n" );
	}
	else if( depth >= 300 && depth < 700 )
	{
		fprintf( fp, "zrange=100,50,700\n" );
	}
	else
	{
		fprintf( fp, "zrange=3,3,33\n" );
	}

	fprintf( fp, "evla=%g\n", s[0].evla );
	fprintf( fp, "evlo=%g\n", s[0].evlo );
	fprintf( fp, "dt=0.15\n" );
	fprintf( fp, "nt=2048\n" );
	fprintf( fp, "fmax=0.4\n" );
	fprintf( fp, "t0=0.\n" );
	fprintf( fp, "redv=18.\n" );
	fprintf( fp, "damp=1.\n" );
	fprintf( fp, "kmax=20000\n" );
	fprintf( fp, "eps=0.0005\n" );
	fprintf( fp, "smin=0.0005\n" );
	fprintf( fp, "modeldb=${MTINV_PATH}/data/modeldb/\n" );
	fprintf( fp, "stadb=../Data/rdseed.stations\n" );
	fprintf( fp, "noverbose\n" );
	fprintf( fp, "nodump\n" );
	fprintf( fp, "EOF\n" );
	fprintf( fp, "\n" );

	fprintf( fp, "cat >! mkgrnlib.par << EOF\n" );

	for( j = 0; j < nsta; j++ )
	{
		i = indx[j+1]-1;

		if(verbose)
                {
                  fprintf( stderr,  "%s: j=%d indx[j+1]=%d i=%d %s\n",
			progname, j, indx[j+1], i, s[i].kstnm );
                }

		if( strcmp( s[i].khole, "" ) == 0 || strcmp( s[i].khole, "00" ) == 0 || strcmp( s[i].khole, "-12345" ) == 0 )
                  strcpy( com, "" );
                else
                  strcpy( com, "#" );

		dt = 0.35;
		if( s[i].dist < 2000 ) dt = 0.30;
		if( s[i].dist < 1000 ) dt = 0.25;
		if( s[i].dist < 800 )  dt = 0.20;
		if( s[i].dist < 500 )  dt = 0.15;
		if( s[i].dist < 200 )  dt = 0.10;
		if( s[i].dist < 100 )  dt = 0.05;
		if( s[i].dist < 50 )   dt = 0.02;

		if( s[i].dist  > 0 && s[i].dist < 2000 )
		{
		  fprintf( fp, "%s%-8s %-2s %s.par %.2f ### R=%6.0f Az=%03.0f %s.%s.%s\n", 
			com,
			s[i].kstnm, 
			s[i].knetwk, 
			model_name, 
			dt, 
			s[i].dist, 
			s[i].az,
			s[i].knetwk, s[i].kstnm, s[i].khole );
		}
		else
		{
		  fprintf( fp, "#%s%-8s %-2s %s.par %.2f ### R=%6.0f Az=%03.0f %s.%s.%s\n",
			com,
                        s[i].kstnm, 
                        s[i].knetwk,
                        model_name,
                        dt, 
                        s[i].dist, 
                        s[i].az,
			s[i].knetwk, s[i].kstnm, s[i].khole );
		}
	}

	fprintf( fp, "EOF\n" );
	
	fprintf( fp, "\n" );

	fprintf( fp, "multithread_mkgrnlib \\\n" );
	fprintf( fp, "     parfile=mkgrnlib.par \\\n" );
	fprintf( fp, "     executable_pathname=${MTINV_PATH}/bin/mkgrnlib > multithread_mkgrnlib.out\n" );
	
	fprintf( fp, "\n" );

	if( igmt5 ) strcpy( gmt_version_string, "gmt5" );
	else strcpy( gmt_version_string, "nogmt5" );
	if( irt )
	{
		strcpy( realtime_string, "realtime" );
	}
	else
	{
		strcpy( realtime_string, "norealtime" );
	}

	fprintf( fp, "makepar com=\"%s\" \\\n", comment );
        fprintf( fp, "    date=\"%4d/%02d/%02d,%02d:%02d:%05.2f\" \\\n", 
		ot->year, ot->month, ot->mday, 
		ot->hour, ot->min, ot->fsec );
        fprintf( fp, "    DataDir=../Data \\\n" );
        fprintf( fp, "    RespDir=../Resp \\\n" );
        fprintf( fp, "    %s nooracle nomysql sqlite\\\n", gmt_version_string );
	fprintf( fp, "    lf=0.02 hf=0.05 \\\n" );
	fprintf( fp, "    minsnr=3.0 ctol=0.85 maxshift=10 %s *.glib\n", realtime_string );

	fclose(fp);

	free(dist);
	free(indx);

	chmod( script_filename, S_IRWXU|S_IRWXG|S_IRWXO );
}

void sort_sac_headers_by_dist( int nsta, Sac_Header *s, float *dist, int *indx )
{
	int i;
	void indexx( int n, float *arrin, int *indx );
	int verbose = 0;

	for( i = 0; i < nsta; i++ )
	{
		dist[i+1] = s[i].dist;
		indx[i+1] = i;
	}

	if( nsta == 1 ) {
		dist[1] = s[0].dist;
		indx[1] = 1;
		return;
	}

	indexx( nsta, dist, indx );

	if(verbose)
	{
	  for( i = 1; i <= nsta; i++ )
	  {
		fprintf( stdout, "nsta=%d i=%d dist[i]=%g indx[i]=%d dist[indx]=%g %g %s.%s.%s\n", 
			nsta, i, dist[i], indx[i], dist[indx[i]], 
			s[indx[i]-1].dist, s[indx[i]-1].knetwk, s[indx[i]-1].kstnm, s[indx[i]-1].khole );
	  }
	}
}

void make_gmt4_map( Sac_Header *s, int nsta, char *script_filename, int irt, int verbose )
{
	FILE *fp;
        int i = 1;
        float maxlat = -999, maxlon = -999;
        float minlat = +999, minlon = +999;
	float aspect_ratio;
        float sumlat = 0, sumlon = 0;
        float ticklength = 0.5, annotation = 1.0;
        char PS_filename[128], JPG_filename[128];
        char command_line[512];

        sprintf( PS_filename, "sac2gmtmap.ps" );
        sprintf( JPG_filename, "sac2gmtmap.jpg" );

        if( (fp=fopen( script_filename,"w" )) == NULL )
        {
                fprintf( stderr, "%s: cannot open file %s for writting\n",
                        progname, script_filename );
                exit(-1);
        }

	maxlat = s[0].evla;
        minlat = s[0].evla;
        maxlon = s[0].evlo;
        minlon = s[0].evlo;

        for( i = 0; i < nsta; i++ )
        {
                if(  s[i].stla > maxlat ) maxlat = s[i].stla;
                if(  s[i].stla < minlat ) minlat = s[i].stla;
                if(  s[i].stlo > maxlon ) maxlon = s[i].stlo;
                if(  s[i].stlo < minlon ) minlon = s[i].stlo;

        }
        maxlat =  ceil( maxlat ) + 1;
        minlat = floor( minlat ) - 1;
        maxlon =  ceil( maxlon ) + 1;
        minlon = floor( minlon ) - 1;

	aspect_ratio = fabs(maxlat - minlat) / fabs( maxlon - minlon );
        if(verbose) fprintf( stderr, "%s: aspect_ratio = %g\n", progname, aspect_ratio );

        if( fabs( maxlon - minlon ) < 3 )
        {
                maxlon += 2;
                minlon -= 2;
        }

/*** make GMT version 4.5.x map ***/

	fprintf( fp, "#!/bin/csh\n");
        fprintf( fp, "###################################################################################\n");
        fprintf( fp, "## This C-shell script was automatically generated by mtinv                      ##\n");
        fprintf( fp, "## and requires Generic Mapping Tools (GMT) http://gmt.soest.hawaii.edu/         ##\n");
        fprintf( fp, "## The script plots the mechanism at the event location and station locations    ##\n");
        fprintf( fp, "## on a map. This script uses GMT version 4.x.x only see usage for ver 5.x.x     ##\n");
        fprintf( fp, "###################################################################################\n");
        fprintf( fp, "\n" );

        fprintf( fp, "gmtset BASEMAP_TYPE plain\n" );
        fprintf( fp, "\n");

        fprintf( fp, "pscoast -R%g/%g/%g/%g -JM6i -Di ", minlon, maxlon, minlat, maxlat );
        fprintf( fp, " -N1/1.2p,black,5_2:0p -N2/0.8p,black,5_2:0p " );
        fprintf( fp, " -A1000 -W1p,black -Glightgray -P -K >! %s\n", PS_filename );
        fprintf( fp, "\n");

        fprintf( fp, "psxy -R -JM -m -O -K >> %s << EOF\n", PS_filename );
        for( i = 0; i < nsta; i++ )
        {
                fprintf( fp, "> -W1p,black %s.%s.%s%s\n",
                        s[i].kstnm, s[i].knetwk, s[i].kcmpnm, s[i].khole );
                fprintf( fp, "%g %g\n", s[i].evlo, s[i].evla );
                fprintf( fp, "%g %g\n", s[i].stlo, s[i].stla );
        }
        fprintf( fp, "EOF\n" );
        fprintf( fp, "\n");

        fprintf( fp, "psxy -R -JM -St0.15i -W1p,black -Gred -O -K >> %s << EOF\n", PS_filename );
        for( i = 0; i < nsta; i++ )
        {
                fprintf( fp, "%g %g %s.%s.%s%s\n",
                        s[i].stlo, s[i].stla,
                        s[i].kstnm, s[i].knetwk,
                        s[i].kcmpnm, s[i].khole );
        }
        fprintf( fp, "EOF\n" );
        fprintf( fp, "\n");

	fprintf( fp, "#Reads (x,y[,fontinfo,angle,justify],text) from <table> [or stdin].\n" );

        fprintf( fp, "pstext -R -JM -N -D0.05i/0.1i -Wwhite -O -K >> %s << EOF\n", PS_filename );
        for( i = 0; i < nsta; i++ )
        {
		if( strncmp( s[i].knetwk, "-12345", 6 ) == 0 ) 
		{
			fprintf( fp, "%g %g 6 0 0 0 %s\n", 
				s[i].stlo, s[i].stla, s[i].kstnm );
		}
		else
		{
                  fprintf( fp, "%g %g 6 0 0 0 %s.%s.%s%s\n",
                        s[i].stlo, s[i].stla,
                        s[i].kstnm, s[i].knetwk,
                        s[i].kcmpnm, s[i].khole );	
		}
        }
        fprintf( fp, "EOF\n" );
        fprintf( fp, "\n");

/*** check event location ***/

        if( ( s[0].evlo == -12345 || s[0].evlo == 0 ) &&
            ( s[0].evla == -12345 || s[0].evla == 0 ) )
        {
                fprintf( fp, "## %s: evla evlo undefined in SAC file\n", progname );
                fprintf( stderr, "%s: evla evlo undefined in SAC file\n",
                        progname );
                for( i = 0; i < nsta; i++ )
                {
                        sumlat += s[i].stla;
                        sumlon += s[i].stlo;
                }
                sumlat /= nsta;
                sumlon /= nsta;
          fprintf( fp, "echo %g %g | ", sumlon, sumlat );
          fprintf( fp, "psxy -R -JM -Sa0.2i -W1p,black -Gmagenta -O -K >> %s\n", PS_filename );
        }
        else
        {
          fprintf( fp, "echo %g %g | ", s[0].evlo, s[0].evla );
          fprintf( fp, "psxy -R -JM -Sa0.2i -W1p,black -Gmagenta -O -K >> %s\n", PS_filename );
        }
        fprintf( fp, "\n");

        if( fabs( maxlon - minlon ) < 5 )
        {
                ticklength = 0.5;
                annotation = 2;
        }
        else
        {
                ticklength = 1;
                annotation = 5;
        }

        fprintf( fp, "psbasemap -R -JM -Bf%ga%g/f%ga%gNSEW -Lx1/1/%g/100 -O >> %s\n",
                ticklength, annotation, ticklength, annotation, s[0].evla, PS_filename );
        fprintf( fp, "\n");

        fprintf( fp, "ps2raster -A  -Tj -E300 %s\n", PS_filename );
        fprintf( fp, "/bin/rm -f %s gmt.conf gmt.history .gmtdefaults4 .gmtcommands4\n", PS_filename );

	if(irt)
		fprintf( fp, "#echo use eog xv open to view jpg\n" );
	else
		fprintf( fp, "eog %s\n", JPG_filename );

        fclose(fp);

        chmod( script_filename, S_IRWXU|S_IRWXG|S_IRWXO );
	
	sprintf( command_line, "./%s", script_filename );
	fprintf( stderr, "%s: running command %s\n", progname, command_line );
	system( command_line );
}

/*** GMT Version 5.x.x ****/

void make_gmt5_map( Sac_Header *s, int nsta, char *script_filename, int irt, int verbose )
{
	FILE *fp;
	int i = 1;
	float maxlat = -999, maxlon = -999;
	float minlat = +999, minlon = +999;
	float aspect_ratio;
	float sumlat = 0, sumlon = 0;
        float ticklength = 0.5, annotation = 1.0;
	char PS_filename[128], JPG_filename[128];
	char command_line[512];

	sprintf( PS_filename, "sac2gmtmap.ps" );
	sprintf( JPG_filename, "sac2gmtmap.jpg" );

	if( (fp=fopen( script_filename,"w" )) == NULL )
	{
		fprintf( stderr, "%s: cannot open file %s for writting\n",
			progname, script_filename );	
		exit(-1);
	}

	maxlat = s[0].evla;
	minlat = s[0].evla;
	maxlon = s[0].evlo;
	minlon = s[0].evlo;

	for( i = 0; i < nsta; i++ )
	{
		if(  s[i].stla > maxlat ) maxlat = s[i].stla;
                if(  s[i].stla < minlat ) minlat = s[i].stla;
                if(  s[i].stlo > maxlon ) maxlon = s[i].stlo;
                if(  s[i].stlo < minlon ) minlon = s[i].stlo;

        }
        maxlat =  ceil( maxlat ) + 1;
        minlat = floor( minlat ) - 1;
        maxlon =  ceil( maxlon ) + 1;
        minlon = floor( minlon ) - 1;

	aspect_ratio = fabs(maxlat - minlat) / fabs( maxlon - minlon );
	if(verbose) fprintf( stderr, "%s: aspect_ratio = %g\n", progname, aspect_ratio );

        if( fabs( maxlon - minlon ) < 3 )
        {
                maxlon += 2;
                minlon -= 2;
        }

	if( fabs(minlon) > 180 ) minlon = -180;
	if( fabs(maxlon) > 180 ) maxlon = 180;
	if( fabs(minlat) > 89 ) minlat = -89;
	if( fabs(maxlat) > 89 ) maxlat = 89;

/*** make GMT version 5.x.x map ***/

	fprintf( fp, "#!/bin/csh\n");
        fprintf( fp, "###################################################################################\n");
        fprintf( fp, "## This C-shell script was automatically generated by mtinv                      ##\n");
        fprintf( fp, "## and requires Generic Mapping Tools (GMT) http://gmt.soest.hawaii.edu/         ##\n");
        fprintf( fp, "## The script plots the mechanism at the event location and station locations    ##\n");
        fprintf( fp, "## on a map. This script uses GMT version 5.x.x only see usage for ver 4.x.x     ##\n");
        fprintf( fp, "###################################################################################\n");
	fprintf( fp, "\n" );

	fprintf( fp, "gmt set MAP_FRAME_TYPE plain\n" );
        fprintf( fp, "gmt set FORMAT_GEO_OUT DG\n" );
        fprintf( fp, "gmt set FORMAT_GEO_MAP DG\n" );
        fprintf( fp, "\n");

	fprintf( fp, "pscoast -R%g/%g/%g/%g -JM6i -Di ", minlon, maxlon, minlat, maxlat );
	fprintf( fp, " -N1/1.2p,black,5_2:0p -N2/0.8p,black,5_2:0p " );
	fprintf( fp, " -A1000 -W1p,black -Glightgray -P -K >! %s\n", PS_filename );
	fprintf( fp, "\n");

        fprintf( fp, "psxy -R -JM -O -K >> %s << EOF\n", PS_filename );
        for( i = 0; i < nsta; i++ )
        {
                fprintf( fp, "> -W1p,black %s.%s.%s%s\n", 
			s[i].kstnm, s[i].knetwk, s[i].kcmpnm, s[i].khole );
                fprintf( fp, "%g %g\n", s[i].evlo, s[i].evla );
                fprintf( fp, "%g %g\n", s[i].stlo, s[i].stla );
        }
        fprintf( fp, "EOF\n" );
	fprintf( fp, "\n");

        fprintf( fp, "psxy -R -JM -St0.15i -W1p,black -Gred -O -K >> %s << EOF\n", PS_filename );
        for( i = 0; i < nsta; i++ )
        {
                fprintf( fp, "%g %g %s.%s.%s%s\n", 
			s[i].stlo, s[i].stla, 
			s[i].kstnm, s[i].knetwk, 
			s[i].kcmpnm, s[i].khole );
        }
        fprintf( fp, "EOF\n" );
	fprintf( fp, "\n");

	fprintf( fp, "#Reads (x,y[,fontinfo,angle,justify],text) from <table> [or stdin].\n" );

        fprintf( fp, "pstext -R -JM -C0.01i/0.01i -N -D0.05i/0.1i -W0p,white -Tc -Gwhite " );
	fprintf( fp, "  -F+f8p,Times-Roman,blue+jBL -O -K >> %s << EOF\n", PS_filename );
        for( i = 0; i < nsta; i++ )
        {
                fprintf( fp, "%g %g %s.%s.%s%s\n", 
			s[i].stlo, s[i].stla, 
			s[i].kstnm, s[i].knetwk,
			s[i].kcmpnm, s[i].khole );
        }
        fprintf( fp, "EOF\n" );
	fprintf( fp, "\n");

/*** check event location ***/
	if( ( s[0].evlo == -12345 || s[0].evlo == 0 ) &&
	    ( s[0].evla == -12345 || s[0].evla == 0 ) )
	{
		fprintf( fp, "## %s: evla evlo undefined in SAC file\n", progname );
		fprintf( stderr, "%s: evla evlo undefined in SAC file\n",
			progname );
		for( i = 0; i < nsta; i++ )
		{
			sumlat += s[i].stla;
			sumlon += s[i].stlo;
		}
		sumlat /= nsta;
		sumlon /= nsta;
	  fprintf( fp, "echo %g %g | ", sumlon, sumlat );
          fprintf( fp, "psxy -R -JM -Sa0.2i -W1p,black -Gmagenta -O -K >> %s\n", PS_filename );
	}
	else
	{
          fprintf( fp, "echo %g %g | ", s[0].evlo, s[0].evla );
	  fprintf( fp, "psxy -R -JM -Sa0.2i -W1p,black -Gmagenta -O -K >> %s\n", PS_filename );
	}
	fprintf( fp, "\n");

        if( fabs( maxlon - minlon ) < 5 )
        {
                ticklength = 0.5;
                annotation = 2;
        }
        else
        {
                ticklength = 1;
                annotation = 5;
        }

        fprintf( fp, "psbasemap -R -JM -Bxf%ga%g -Byf%ga%g -BNSEW -Lx1/1/%g/100 -O >> %s\n",
                ticklength, annotation, ticklength, annotation, s[0].evla, PS_filename );
	fprintf( fp, "\n");

        fprintf( fp, "psconvert -A -Z -Tj -E300 %s\n", PS_filename );
	fprintf( fp, "/bin/rm -f gmt.conf gmt.history\n" );

	if(irt)
		fprintf( fp, "#echo use eog xv open to view jpg\n" );
	else
		fprintf( fp, "open %s\n", JPG_filename );

	fclose(fp);

	chmod( script_filename, S_IRWXU|S_IRWXG|S_IRWXO );

	sprintf( command_line, "./%s", script_filename ); 
	fprintf( stderr, "%s: running command %s\n", progname, command_line );
       	system( command_line );
}

void Usage( int ac, char **av )
{
	fprintf( stderr, "\n\n" );
	fprintf( stderr, "%s: \n", av[0] );
	fprintf( stderr,
          "\t Usage %s -z [depth km] -ev [lat lon] -ot [\"YYYY/MM/DD,hh:mm:ss.ss\"] -com [\"Comment\"] -gmt5 -help -realtime -verbose [sac files]\n",
                progname );
	fprintf( stderr, "\n\n" );
	fprintf( stderr, "\t -z  (float) source depth in kilometers (does not overwrite SAC file)\n" );
	fprintf( stderr, "\t -ev (float float) source latitude and longitude in decimal degrees (does overwrite SAC files)\n" );
	fprintf( stderr, "\t -ot (string) source origin time format: YYYY/MM/DD,hh:mm:ss.ss or YYYY-MM-DDThh:mm:ss.ss\n" );
	fprintf( stderr, "\t -com (string) a comment that gets printed on the moment tensor results plot (typical region name and mb evid)\n" );
	fprintf( stderr, "\t -gmt5 (boolean flag) creates map with event and station locations using GMT version 5.x.x [default GMT v4.5.x]\n" );
	fprintf( stderr, "\t -h or -help (boolean flag) prints this page\n" );
	fprintf( stderr, "\t -v or -verbose (boolean flag) verbose output\n" );
	fprintf( stderr, "\t list of sac files can use wildcards\n" );
	fprintf( stderr, "\t -realtime or -rt (boolean flag) interactive run script and display map on screen [default off]\n" );
	fprintf( stderr, "\n\n" );
	fprintf( stderr, "\t This program creates a makeglib.csh file that computes Green functions.\n" );
	fprintf( stderr, "\t Add SAC files \"../Data/*.BHZ.?.SAC\" to tell setupMT which stations to use\n" );
	fprintf( stderr, "\t in the moment tensor inversion\n" );
	fprintf( stderr, "\n\n" );

	exit(0);
}
