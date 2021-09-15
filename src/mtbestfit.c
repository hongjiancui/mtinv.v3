#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <sys/stat.h>
#include <errno.h>

#include <sys/utsname.h>
#include <netdb.h>
#include <unistd.h>

#include "../include/nrutil.h"
#include "../include/mt.h"

char progname[128];

typedef struct {
	float ot_shift;
	float z;
	float fsec;
	float vred;
	float fit;
	float pdc;
	float piso;
	float pclvd;
	char mt_type[6];
	int mtdegfree;
	char author[32];
	char comment[128];
	char email_file[256];
	double originTimeEpoch;
	MyTime ot;
	float lat, lon;

	int ishift;
	float cortol;
	float maxtimeshift;
	int iuse_snr;
	float minsnr;
	int igmt5;
	int sqlite3_db_write;
	int mysql_db_write;
	int oracle_db_write;

	float vred_max, fit_max;
	float fit_max_threshold, vred_diff_threshold;
	float vred_diff;

	int npages;
	char **psplotfiles;
} bestFitMT;

int main( int ac, char **av )
{
	FILE *fp;
	int i, j, i_best_vred, i_best_fit, nsol;
	bestFitMT *a;
	char tag[32], pstag[32];
	char rec[512];
	float vred_max, fit_max;

	char filename[] = { "automt.txt" };
	float fit_max_threshold = 20;
	float vred_diff_threshold = 30;
	int force_best_vred = 1;

/*** functional prototypes ***/

	void write_mtbest( bestFitMT *a );
	void make_run( bestFitMT *a, char *fitType );
	void createResultsWebpage( bestFitMT *a );

	void WriteMyTime2STDOUT( MyTime *t );
	MyTime *epoch2time( MyTime *t, double epoch );

	int setpar(int ac, char **av);
	int getpar();
	void endpar();

/*** begin main(), open input file ***/

	strcpy( progname, av[0] );
	
	setpar( ac, av );
	getpar( "force_best_vred", "b", &force_best_vred );
	endpar();

	if( (fp = fopen( filename, "r" )) == NULL )
	{
		fprintf( stderr, "%s: Cannot open file %s\n", 
			progname, filename );
		exit(-1);
	}

/*** read file and parse ***/

	a =calloc(1,sizeof(bestFitMT));
	i = 0;
	while( fgets( rec, 512, fp ) != NULL )
	{
		sscanf( rec, "%s", tag );

		if( strcmp( tag, "BEGIN" ) == 0 )
		{
			a = realloc( a, (i+1)*sizeof(bestFitMT));
		
			sscanf( rec, "%s %f %f %f %f %f %f %s %d %s",
				tag,
				&(a[i].ot_shift),
				&(a[i].z),
				&(a[i].fsec),
				&(a[i].vred),
				&(a[i].fit),
				&(a[i].pdc),
				a[i].mt_type,
				&(a[i].mtdegfree),
				a[i].author );
		}

		if( strcmp( tag, "CMDLINE" ) == 0 )
		{
			/*   CMDLINE   1 0.85 10 1 3 1 0 0 0 **/
			sscanf( rec, "%s %d %g %g %d %g %d %d %d %d",
				tag,
				&(a[i].ishift),
				&(a[i].cortol),
				&(a[i].maxtimeshift),
				&(a[i].iuse_snr),
				&(a[i].minsnr),
				&(a[i].igmt5),
				&(a[i].sqlite3_db_write),
				&(a[i].mysql_db_write),
				&(a[i].oracle_db_write) );	
		}

		if( strcmp( tag, "COM" ) == 0 )
			sscanf( rec, "%s %[^\n]", tag, a[i].comment );

		if( strcmp( tag, "EMAIL" ) == 0 )
			sscanf( rec, "%s %s", tag, a[i].email_file );

		if( strcmp( tag, "NPAGES" ) == 0 )
		{
			sscanf( rec, "%s %d", tag, &(a[i].npages) );

		/*** allocate memory for PS/JPEG plot files ***/

			a[i].psplotfiles = calloc( a[i].npages, sizeof(char *) );
			for( j = 0; j < a[i].npages; j++ )
			{
				a[i].psplotfiles[j] = calloc( 128, sizeof(char) );
			}

			for( j = 0; j < a[i].npages; j++ )
			{
				fgets( rec, 512, fp );
				sscanf( rec, "%s", tag );

				sprintf( pstag, "PSPLOT%02d", j+1 );

				/*** Debug ***/
				/* fprintf( stdout, "tag=%s pstag=%s rec=%s\n", tag, pstag, rec ); */

				if( strcmp( tag, pstag ) == 0 ) 
				{
					sscanf( rec, "%s %s", pstag, a[i].psplotfiles[j] );
				}
			}

		} /*** reading NPAGES PS/JPEG plot files ***/

		if( strcmp( tag, "PISO" )  == 0 )
			sscanf( rec, "%s %f", tag, &(a[i].piso) );

		if( strcmp( tag, "PCLVD" ) == 0 )
			sscanf( rec, "%s %f", tag, &(a[i].pclvd) );

		if( strcmp( tag, "LAT" ) == 0 )
			sscanf( rec, "%s %f", tag, &(a[i].lat) );

		if( strcmp( tag, "LON" ) == 0 )
			sscanf( rec, "%s %f", tag, &(a[i].lon) );

		if( strcmp( tag, "TIME" ) == 0 )
		{	
			sscanf( rec, "%s %lf", tag, &(a[i].originTimeEpoch) );
			epoch2time( &(a[i].ot), a[i].originTimeEpoch );
			/* WriteMyTime2STDOUT( &(a[i].ot) ); */
		}
	
		
		if( strcmp( tag, "END" ) == 0 )
		{
			i++;
		}
	}
	nsol = i;
	/* rewind(fp); */
	fclose(fp);

	i_best_fit  = 0;
	i_best_vred = 0;
	vred_max    = 0;
	fit_max     = 0;

	for( i = 0; i < nsol; i++ )
	{
		a[i].fit_max_threshold = fit_max_threshold;
		a[i].vred_diff_threshold = vred_diff_threshold;
		if( a[i].fit > fit_max )
		{
			fit_max = a[i].fit;
			i_best_fit = i;
		}

		if( a[i].vred > vred_max ) 
		{
			vred_max = a[i].vred;
			i_best_vred = i;
		}

		/* write_mtbest( &(a[i]) ); */
	}
	/* write_mtbest( &(a[i_best_vred]) ); */

/*** determine which is better fit or vred ? ***/

	a[i_best_vred].vred_max = vred_max;
	a[i_best_fit].vred_max  = vred_max;
	a[i_best_vred].fit_max  = fit_max;
	a[i_best_fit].fit_max   = fit_max;
        a[i_best_vred].vred_diff= fabs( a[i_best_vred].vred - a[i_best_fit].vred );
	a[i_best_fit].vred_diff = fabs( a[i_best_vred].vred - a[i_best_fit].vred );

	if( force_best_vred )
	{
		make_run( &(a[i_best_vred]), "FORCE BEST VARIANCE-REDUCTION" );
                createResultsWebpage( &(a[i_best_vred]) );
	}
	else
	{
	  if(	( fit_max > fit_max_threshold ) && 
		( fabs( a[i_best_vred].vred - a[i_best_fit].vred ) > vred_diff_threshold ) )
	  {
		make_run( &(a[i_best_fit]), "BEST FIT" );
		createResultsWebpage( &(a[i_best_fit]) );
	  }
	  else
	  {
		make_run( &(a[i_best_vred]), "BEST VARIANCE-REDUCTION" );
		createResultsWebpage( &(a[i_best_vred]) );
	  }
	}


} /*** END OF MAIN() ***/


/******************************/
/**** void make_run() *********/
/******************************/

void make_run( bestFitMT *a, char *fitType )
{
	int i;
	FILE *fp;
	char gmtstring[8];
	char dbprog[8];

/*** set some flags ***/

	a->igmt5 = 1;
        strcpy( gmtstring, "gmt5" ); /*** default ***/

        if(a->igmt5)
          strcpy( gmtstring, "gmt5" );
        else
          strcpy( gmtstring, "nogmt5" );

/*** use this one for automt ***/

	a->sqlite3_db_write = 1;
        strcpy( dbprog, "sqlite3" ); /*** default ***/

        if(a->oracle_db_write)
          strcpy( dbprog, "oracle" );
        else if(a->mysql_db_write)
          strcpy( dbprog, "mysql" );
        else if(a->sqlite3_db_write)
          strcpy( dbprog, "sqlite" );

/*** begin write run2.csh file ***/

	fp = fopen( "run2.csh", "w" );

	fprintf( fp, "#!/bin/csh \n");

	fprintf( fp, "###\n" );
	fprintf( fp, "### Created by %s fitType=%s\n", progname, fitType );
	fprintf( fp, "### ot=%g fsec=%g z=%g vred=%g fit=%g pdc=%g piso=%g pclvd=%g\n",
		a->ot_shift, a->fsec, a->z, a->vred, a->fit,
		roundf(a->pdc), roundf(a->piso), roundf(a->pclvd) );

	fprintf( fp, "### fit_max=%g fit_max_threshold=%g vred_diff=%g vred_diff_threshold=%g\n", 
			a->fit_max, 
			a->fit_max_threshold, 
			a->vred_diff, 
			a->vred_diff_threshold );

	fprintf( fp, "###\n" );

        fprintf( fp, "setenv MTINV_GMT_GRID_FILE /Users/ichinose1/Work/topogmt/etopo5.grd\n" );
        fprintf( fp, "setenv MTINV_GMT_INT_FILE  /Users/ichinose1/Work/topogmt/etopo5.int\n" );
        fprintf( fp, "setenv MTINV_GMT_CPT_FILE  /Users/ichinose1/Work/topogmt/etopo5.cpt\n" );
        fprintf( fp, "setenv MT_DATABASE_FILE  /Users/ichinose1/Work/mtinv.v3.0.5/data/mt.db\n" );

        if( strcmp( a->mt_type, "FULL" ) == 0 )
        {
           fprintf( fp, "set DEGFREE=6  # 1-isotropic_mt 5-deviatoric_mt 6-full_mt\n" );
        }
        else if( strcmp( a->mt_type, "DEV" ) == 0 )
        {
           fprintf( fp, "set DEGFREE=5  # 1-isotropic_mt 5-deviatoric_mt 6-full_mt\n" );
        }
	else if( strcmp( a->mt_type, "ISO" ) == 0 )
	{
		fprintf( fp, "set DEGFREE=1  # 1-isotropic_mt 5-deviatoric_mt 6-full_mt\n" );
	}
	else
	{
		fprintf( fp, "set DEGFREE=5  # 1-isotropic_mt 5-deviatoric_mt 6-full_mt\n" );
	}
        fprintf( fp, "\n" );

	fprintf( fp, "###\n" );
        fprintf( fp, "### Clean Up\n" );
        fprintf( fp, "###\n" );

        fprintf( fp, "/bin/rm -f email_T???.?sec_Z???.?km_.txt plot_T???.?sec_Z???.?km_.p??.jpg *.ps\n" );
        fprintf( fp, "\n" );

        fprintf( fp, "### MAKE GMT PLOT WITH LOCATION/STATION/SOLUTION ###\n" );
        fprintf( fp, "###\n" );

        fprintf( fp, "mtinv \\\n" );
        fprintf( fp, "\t ts0=%g         \\\n", -1 * a->ot_shift );
        fprintf( fp, "\t par=mtinv.par  \\\n" );
        fprintf( fp, "\t %s             \\\n", gmtstring );
        fprintf( fp, "\t mtdegfree=${DEGFREE} \\\n" );
        fprintf( fp, "\t gmtmap         \\\n" );
        fprintf( fp, "\t nodumpxy       \\\n" );
        fprintf( fp, "\t %s             \\\n", dbprog );
        fprintf( fp, "\t PltXcorLabel   \\\n" );
        fprintf( fp, "\t use_snr        \\\n" );
        fprintf( fp, "\t minsnr=%g      \\\n", a->minsnr );
        fprintf( fp, "\t ctol=%g        \\\n", a->cortol );
        fprintf( fp, "\t maxshift=%g >> mtinv.out\n", a->maxtimeshift );
        fprintf( fp, "\n" );


	for( i = 0; i < a->npages; i++ )
	{
	  if(a->igmt5)
          {
		fprintf( fp, "psconvert -Z -Tj -E300 plot_T???.?sec_Z???.?km_.p??.ps\n" );
                /* fprintf( fp, "psconvert -Z -Tj -E300 %s\n", a->psplotfiles[i] ); */
          }
          else
          {
		fprintf( fp, "ps2raster -Tj -E300 plot_T???.?sec_Z???.?km_.p??.ps\n" );
                /* fprintf( fp, "ps2raster -Tj -E300 %s\n", a->psplotfiles[i] ); */
          }
	}

        fprintf( fp, "csh gmtmap.csh\n" );
	fprintf( fp, "#open gmtmap.jpg\n" );
	fprintf( fp, "\n" );

	fprintf( fp, "### dumpxy option for publication quality GMT plots\n" );
        fprintf( fp, "#gmtwf.csh\n" );
        fprintf( fp, "\n" );

	if( a->sqlite3_db_write )
        {
                fprintf( fp, "\n" );
                fprintf( fp, "sqlite3 ${MT_DATABASE_FILE} << EOF\n" );
                fprintf( fp, ".read insert.sql\n" );
                fprintf( fp, ".quit\n" );
                fprintf( fp, "EOF\n" );

                fprintf( fp, "updateMTdb\n" );
                fprintf( fp, "\n" );
                fprintf( fp, "# list_MTdb.csh ${MT_DATABASE_FILE}\n" );
                fprintf( fp, "print_MTdb.csh > db.txt\n" );
                fprintf( fp, "# remove_MTdb.csh\n" );
        }
        fclose(fp);

        chmod( "run2.csh", S_IRWXU|S_IRGRP|S_IXGRP|S_IROTH|S_IXOTH );
}

void write_mtbest( bestFitMT *a )
{
	int j;
	fprintf( stdout, "%f %f %f %f %f %f %s %d %s %s %s %d : \n",
		a->ot_shift,
		a->z,
		a->fsec,
		a->vred,
		a->fit,
		a->pdc,
		a->mt_type,
		a->mtdegfree,
		a->author,
		a->comment,
		a->email_file,
		a->npages );
	for( j = 0; j < a->npages; j++ )
		fprintf( stdout, "\t(%s)\n", a->psplotfiles[j] );
}

void createResultsWebpage( bestFitMT *a )
{
        FILE *fp;
	int i = 0;

        fp = fopen( "index.html", "w" );
        fprintf( fp, "<!DOCTYPE html>\n" );
        fprintf( fp, "<HTML>\n" );
        fprintf( fp, "<HEAD>\n" );
        fprintf( fp, "<TITLE>\n" );

	fprintf( fp, "%4d-%02d-%02d (%03d) %02d:%02d:%05.2f %.4f %.4f %s\n", 
		a->ot.year,
		a->ot.month,
		a->ot.mday,
		a->ot.jday,
		a->ot.hour,
		a->ot.min,
		a->ot.fsec,
		a->lat,
		a->lon,
		a->comment );

        fprintf( fp, "</TITLE>\n" );
        fprintf( fp, "</HEAD>\n" );

/*** email file ***/

	fprintf( fp, "<H1>%4d-%02d-%02d (%03d) %02d:%02d:%05.2f %.4f %.4f %s</H1>\n", 
                a->ot.year,
                a->ot.month,
                a->ot.mday,
                a->ot.jday,
                a->ot.hour,
                a->ot.min,
                a->ot.fsec,
                a->lat,
                a->lon,
                a->comment );

        fprintf( fp, "<P><A HREF=\"%s\">%s</A>\n", a->email_file, a->email_file );
        fprintf( fp, "<DIV ID=\"list\">\n" );
        fprintf( fp, "<P>\n" );
	fprintf( fp, "<IFRAME SRC=\"%s\" WIDTH=1200 HEIGHT=300 FRAMEBORDER=1></IFRAME>\n",
			a->email_file );
        fprintf( fp, "</P></DIV>\n" );
	fprintf( fp, "<P><IMG src=\"gmtmap.jpg\" style=\"width: 60%%; height: auto\" ></P>\n" );
	fprintf( fp, "<P>\n" );
	for( i = 0; i < a->npages; i++ )
	{
		fprintf( fp, "<IMG SRC=\"%s\"style=\"width: 60%%; height: auto\" >\n", 
				a->psplotfiles[i] );
	}
	fprintf( fp, "</P>\n" );
	fprintf( fp, "<P><IMG SRC=\"results.5.jpg\" style=\"width: 60%%; height: auto\"/></P>\n" );
	fprintf( fp, "<P><A HREF=\"plotmech.jpg\">plotmech.jpg</A></P>\n" );
	fprintf( fp, "<P><A HREF=\"plotz.jpg\">plotz.jpg</A></P>\n" );
	fprintf( fp, "<P><A HREF=\"db.txt\">db.txt</A>\n" );
	fprintf( fp, "<DIV ID=\"list\">\n" );
	fprintf( fp, "<P>\n" );
	fprintf( fp, "<IFRAME SRC=\"db.txt\" width=1200 height=300 frameborder=2 ></IFRAME>\n" );
	fprintf( fp, "</P>\n" );
	fprintf( fp, "</DIV>\n" );
	fprintf( fp, "</BODY>\n" );
        fprintf( fp, "</HTML>\n" );
        fclose(fp);
}

