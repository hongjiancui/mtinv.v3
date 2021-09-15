#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <errno.h>
#include "../include/nrutil.h"
#include "../include/mt.h"

extern char progname[128];

static float half=0.5;
static float third=0.333333333333;
static float sixth=0.166666666667;

void write_sac_file( char *filename, Sac_File *s, int verbose )
{
	FILE *fp;
	if( (fp = fopen(filename,"w")) == NULL )
		fprintf(stderr, "%s: %s: %s: cannot open file for writting\n",
			progname, __FILE__, __func__ );

	fwrite(&(s->s),sizeof(Sac_Header),1,fp);
	fwrite(&(s->data[0]),s->s.npts*sizeof(float),1,fp);
	fclose(fp);
}

void compute_synthetics( int is, int iz, EventInfo *ev, Greens **grn, Solution *sol, int mtdegfree )
{
	int it, nt;
	float d2r, fi;
	float *w, *u, *v;
	float *rss, *rds, *rdd, *rep, *zss, *zds, *zdd, *zep, *tss, *tds;
	float Mxx, Myy, Mzz, Mxy, Mxz, Myz;

/* see include/mt.h Mo = math.pow( 10.0, 1.5*(Mw+10.73) ) = 1.2445146117713818e+16; Reference Mw = 0.0 */

	d2r = M_PI/180;
	fi = grn[is][iz].az * d2r;
	nt = grn[is][iz].nt;

	Mxx = sol[iz].moment_tensor[1][1] / base_moment;
	Myy = sol[iz].moment_tensor[2][2] / base_moment;
	Mzz = sol[iz].moment_tensor[3][3] / base_moment;
	Mxy = sol[iz].moment_tensor[1][2] / base_moment;
	Mxz = sol[iz].moment_tensor[1][3] / base_moment;
	Myz = sol[iz].moment_tensor[2][3] / base_moment;

	w = ev[is].syn_t.data;
	u = ev[is].syn_r.data;
	v = ev[is].syn_z.data;

	rss = grn[is][iz].g.rss;
	rds = grn[is][iz].g.rds;
	rdd = grn[is][iz].g.rdd;
	rep = grn[is][iz].g.rep;
	zss = grn[is][iz].g.zss;
	zds = grn[is][iz].g.zds;
	zdd = grn[is][iz].g.zdd;
	zep = grn[is][iz].g.zep;
	tss = grn[is][iz].g.tss;
	tds = grn[is][iz].g.tds;

	for( it = 0; it < nt; it++ )
	{
	  if( mtdegfree == FULL_MOMENT_TENSOR )
	  {
	    w[it]=( Mxx * (  half*sin(2*fi) * tss[it] )
	          + Myy * ( -half*sin(2*fi) * tss[it] )
		  + Mxy * (      -cos(2*fi) * tss[it] )
		  + Mxz * (       sin(fi) * tds[it] )
		  + Myz * (      -cos(fi) * tds[it] ) );

	    u[it]=( Mxx * (  half*cos(2*fi) * rss[it] - sixth*rdd[it] + third*rep[it] )
	          + Myy * ( -half*cos(2*fi) * rss[it] - sixth*rdd[it] + third*rep[it] )
		  + Mxy * (   sin(2*fi) * rss[it] ) 
		  + Mxz * (   cos(fi)   * rds[it] ) 
		  + Myz * (   sin(fi)   * rds[it] )
		  + Mzz * ( third * rdd[it] + third * rep[it] ) );

	    v[it]=( Mxx * (  half*cos(2*fi) * zss[it] - sixth*zdd[it] + third*zep[it] )
                  + Myy * ( -half*cos(2*fi) * zss[it] - sixth*zdd[it] + third*zep[it] )
                  + Mxy * (   sin(2*fi) * zss[it] )
                  + Mxz * (   cos(fi)   * zds[it] )
                  + Myz * (   sin(fi)   * zds[it] )
                  + Mzz * ( third * zdd[it] + third * zep[it] ) );
	  }
	  else if( mtdegfree == DEVIATORIC_MOMENT_TENSOR )
	  {
            w[it]=( Mxx * (  half*sin(2*fi) * tss[it] )
                  + Myy * ( -half*sin(2*fi) * tss[it] )
                  + Mxy * (      -cos(2*fi) * tss[it] )
                  + Mxz * (       sin(fi) * tds[it] )
                  + Myz * (      -cos(fi) * tds[it] ) );
                                                                                                                                    
            u[it]=( Mxx * (  half*cos(2*fi) * rss[it] - half*rdd[it] )
                  + Myy * ( -half*cos(2*fi) * rss[it] - half*rdd[it] )
                  + Mxy * (   sin(2*fi) * rss[it] )
                  + Mxz * (   cos(fi)   * rds[it] )
                  + Myz * (   sin(fi)   * rds[it] ) );
                                                                                                                                    
            v[it]=( Mxx * (  half*cos(2*fi) * zss[it] - half*zdd[it] )
                  + Myy * ( -half*cos(2*fi) * zss[it] - half*zdd[it] )
                  + Mxy * (   sin(2*fi) * zss[it] )
                  + Mxz * (   cos(fi)   * zds[it] )
                  + Myz * (   sin(fi)   * zds[it] ) );
	  }
	  else if( mtdegfree == FORCE_EXPLOSION )
	  {
	    w[it]=0;
	    u[it]=(Mxx*(third*rep[it]) + Myy*(third*rep[it]) + Mzz*(third*rep[it]));
	    v[it]=(Mxx*(third*zep[it]) + Myy*(third*zep[it]) + Mzz*(third*zep[it]));
	  }
	  else
	  {
		fprintf( stderr, "%s: %s: %s: unknown mtdegfree %d\n",
			progname, __FILE__, __func__, mtdegfree );
	  }
	}
}

float compute_l2norm_error( float *a, float *b, int n )
{
	int i;
	float err=0, tmp;
	for( i=0; i<n; i++ )
	{
		tmp = a[i] - b[i]; /*** data - syn ***/
		err += tmp*tmp;
	}
	err /= ( n - 1 );
	return err;
}

float mean( float *d, int n )
{
	int i;
	float sum = 0;
	for( i=1; i<=n; i++ )
	{
		sum += d[i];
	}
	return sum/n;
}

float variance( float *d, int n, float mean )
{
	int i;
	float sum = 0, tmp;
	for( i=1; i<=n; i++ )
	{
		tmp = d[i] - mean;
		sum += tmp*tmp;
	}
	return sum;
}
	
float root_mean_square_variance( float *d, float *s, int n )
{
	int i;
	float sum = 0, tmp;
	for( i = 1; i <=n; i++ )
	{
		tmp = d[i] - s[i];
		sum += tmp * tmp;
	}
	return sqrt(sum/(float)n);
}

float variance_reduction( float *d, float *s, int n0, int n1 )
{
        int i;
        float sum1 = 0, sum2 = 0;
        for( i = n0; i < n1; i++ )
        {	
		sum1 += ( (d[i]-s[i]) * (d[i]-s[i]) );
		sum2 += ( d[i] * d[i] );
        }
        return 100.*(1.0-(sum1/sum2));
}

/******************************************************************************/
/*** write out plotmech files New version appends output to plotmech.txt    ***/
/*** rather than writing out to separate files                              ***/
/******************************************************************************/

void write_gmt_xy_values( Solution *sol, EventInfo *ev, Greens **grn, int iz, int nsta, int verbose )
{
	char gmtfilename[256];
	FILE *fp;
	float ot_tmp;

	ot_tmp = ev[0].ts0 * -1; /*** this is correction so that OT shift times  ***/
			/*** are correctly plotted, also see write_sol_rec(): results.out  ***/
	
	sprintf( gmtfilename, "plotmech.txt" );

	if(verbose) fprintf(stdout, "%s: %s: %s: file=%s\n",
		progname, __FILE__, __func__, gmtfilename );

	if( (fp=fopen(gmtfilename, "a+")) == NULL )
	{
		fprintf(stderr, "%s: %s: %s: cannot write file %s\n",
			progname, __FILE__, __func__, gmtfilename );
		exit(-1);
	}

	if( verbose )
	{
		fprintf( stdout,
		 "%s: %s: %s: s->ot.fsec=%g ev.ot.fsec=%g ot_tmp=%g ev.ts0=%g grn.ts0=%g\n",
			progname, __FILE__, __func__, 
			ev[0].ew.ot.fsec, ev[0].ot.fsec, ot_tmp, 
			ev[0].ts0, grn[0][iz].ts0 );
	}

	/*** psmeca -S<format=m> X, Y, depth, mrr, mtt, mff, mrt, mrf, mtf, exp, newX, newY, event_title ***/

			/* 1   2     3     4     5     6     7     8     9     10    11   12  13   14   15   16  17 */
	fprintf( fp, "%6.1f %6.1f %6.1f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %.1f %.4f %.4f %.1f %.0f %.2f %g %g\n",
		grn[0][iz].evdp,  /* 1 evdp */
		ot_tmp,           /* 2 ot */
		grn[0][iz].evdp,  /* 3 evdp */
		sol[iz].mrr,      /* 4 mrr */
		sol[iz].mtt,      /* 5 mtt */
		sol[iz].mff,      /* 6 mff */
		sol[iz].mrt,      /* 7 mrt */
		sol[iz].mrf,      /* 8 mrf */
		sol[iz].mtf,      /* 9 mtf */
		sol[iz].total_fitness2,     /* 10 total fitness */
                grn[0][iz].evlo,   /* 11 evlo */
                grn[0][iz].evla,   /* 12 evla */
                sol[iz].var_red,   /* 13 var_red */
                sol[iz].PDC,       /* 14 pdc */
                sol[iz].mw,        /* 15 mw */
		sol[iz].ot,         /* 16 ot */
		grn[0][iz].evdp );   /* 17 depth */

	fclose(fp);

} /*** END OF SUBROUTINE : write_gmt_xy_values() ***/


/***********************************************************************************/
/*** BEGIN SUBROUTINE : write_email_message()  make an ASCII plot for fast email ***/
/***********************************************************************************/

void write_email_message( FILE *fp, int nsta, int iz, Solution *sol, EventInfo *ev, Greens **grn, int ifwd )
{
	int ista, nsta_used, kk;
	char grd_mo_type, defining;
	char *username;

	MyTime *gmtnow, *localnow;
	MyTime *myGMTtime( MyTime * );
	MyTime *mylocaltime( MyTime * );

	void ascii_mech( FILE *, int, int, int, float, float, float, float, float, float, float, float, float, float );

/*** get the username from shell environment ***/
	username = getenv( "USER" );

/*** get the system time now ***/

	gmtnow = (MyTime *)calloc(1,sizeof(MyTime));
	gmtnow = myGMTtime( gmtnow );

	localnow = (MyTime *)calloc(1,sizeof(MyTime));
	localnow = mylocaltime( localnow );

/*** is this a inversion result or forward calculation ***/
	if( ifwd )
		fprintf( fp, "Seismic Moment Tensor (Forward Calculaiton)\n\n" );
	else
		fprintf( fp, "Seismic Moment Tensor Solution\n\n" );

	fprintf(fp, "%4d/%02d/%02d (%03d) %02d:%02d:%05.2f %.4f %.4f %s\n",
		ev[0].ot.year, ev[0].ot.month, ev[0].ot.mday, ev[0].ot.jday,
		ev[0].ot.hour, ev[0].ot.min, sol[iz].ot,
		grn[0][iz].evla, grn[0][iz].evlo, ev[0].comment );

	fprintf(fp, "\tDepth = %5.1f (km)\n", grn[0][iz].evdp );
	fprintf(fp, "\tMw    = %5.2f\n", sol[iz].mw );
	fprintf(fp, "\tMo    = %5.2fx10^%2d (dyne x cm)\n", sol[iz].abcassa, sol[iz].exponent );
	fprintf(fp, "\n" );

	if( ifwd )
	{
		fprintf( fp, "\tPercent Double Couple=%3.0f %%\n", sol[iz].PDC );
		fprintf( fp, "\tno CLVD calculated\n" );
		fprintf( fp, "\tno ISO calculated\n" );
		fprintf( fp, "\tno Epsilon, k, or F-factor calculated\n" );
	}
	else if( ifwd == 0 && sol[iz].mt_type == FULL_MT )
	{
		fprintf( fp, "\tPercent Double Couple = %3.0f %%\n", sol[iz].PDC );
		fprintf( fp, "\tPercent CLVD          = %3.0f %%\n", sol[iz].PCLVD );
		fprintf( fp, "\tPercent Isotropic     = %3.0f %%\n", sol[iz].PISO );
		fprintf( fp, "\tLune latitude=%.3f longitude=%.3f\n", 
			sol[iz].lune_lat, sol[iz].lune_lon );
		fprintf( fp, "\tEpsilon=%0.2f k=%0.2f F-factor=%0.2f\n",
			sol[iz].epsilon, sol[iz].k, sol[iz].f_factor );
	}
	else
	{
		fprintf( fp, "\tPercent Double Couple = %3.0f %%\n", sol[iz].PDC );
		fprintf( fp, "\tPercent CLVD          = %3.0f %%\n", sol[iz].PCLVD );
		fprintf( fp, "\tno ISO calculated\n" );
		fprintf( fp, "\tEpsilon=%0.2f\n", sol[iz].epsilon );
	}

	fprintf(fp, "\t Percent Variance Reduction = %6.2f %%\n", sol[iz].var_red );
	fprintf(fp, "\t Total Fit                  =  %.2f \n", sol[iz].total_fitness1 );

	fprintf(fp, "\tMajor Double Couple\n" );
	fprintf(fp, "\t\t            strike dip   rake\n");
	fprintf(fp, "\t\tNodal Plane 1: %3.0f  %2.0f %4.0f\n", 
		sol[iz].Maj.P1.s, sol[iz].Maj.P1.d, sol[iz].Maj.P1.r );
	fprintf(fp, "\t\tNodal Plane 2: %3.0f  %2.0f %4.0f\n", 
		sol[iz].Maj.P2.s, sol[iz].Maj.P2.d, sol[iz].Maj.P2.r );
	fprintf(fp, "\n" );

	if( sol[iz].mt_type == DEVIATORIC ) fprintf( fp, "\tDEVIATORIC MOMENT TENSOR\n" );
	if( sol[iz].mt_type == FULL_MT    ) fprintf( fp, "\tFULL MOMENT TENSOR\n" );
	if( sol[iz].mt_type == EXPLOSION  ) fprintf( fp, "\tEXPLOSION ONLY MOMENT TENSOR\n" );

	fprintf(fp, "\n\tMoment Tensor Elements: Spherical Coordinates\n" );
	fprintf(fp, "\t\tMrr=%6.2f Mtt=%6.2f Mff=%6.2f\n", sol[iz].mrr, sol[iz].mtt, sol[iz].mff );
	fprintf(fp, "\t\tMrt=%6.2f Mrf=%6.2f Mtf=%6.2f EXP=%2d\n", 
		sol[iz].mrt, sol[iz].mrf, sol[iz].mtf, sol[iz].exponent );
	fprintf(fp, "\n" );

	fprintf(fp, "\n\tMoment Tensor Elements: Cartesian Coordinates\n" );
	fprintf(fp, "\t\t%5.2f %5.2f %5.2f\n", sol[iz].mxx, sol[iz].mxy, sol[iz].mxz );
	fprintf(fp, "\t\t%5.2f %5.2f %5.2f\n", sol[iz].mxy, sol[iz].myy, sol[iz].myz );
	fprintf(fp, "\t\t%5.2f %5.2f %5.2f\n", sol[iz].mxz, sol[iz].myz, sol[iz].mzz );

	fprintf(fp, "\n\tEigenvalues:\n");
	fprintf(fp, "\t\tT-axis eigenvalue=%6.2f\n", sol[iz].FullMT.T.ev );
	fprintf(fp, "\t\tN-axis eigenvalue=%6.2f\n", sol[iz].FullMT.B.ev );
	fprintf(fp, "\t\tP-axis eigenvalue=%6.2f\n", sol[iz].FullMT.P.ev );

	fprintf(fp, "\n\tEigenvalues and eigenvectors of the Major Double Couple:\n");
	fprintf(fp, "\t\tT-axis ev=%5.2f trend=%.0f plunge=%.0f\n", sol[iz].Maj.T.ev, sol[iz].Maj.T.az, sol[iz].Maj.T.pl );
	fprintf(fp, "\t\tN-axis ev=%5.2f trend=%.0f plunge=%.0f\n", sol[iz].Maj.B.ev, sol[iz].Maj.B.az, sol[iz].Maj.B.pl );
	fprintf(fp, "\t\tP-axis ev=%5.2f trend=%.0f plunge=%.0f\n", sol[iz].Maj.P.ev, sol[iz].Maj.P.az, sol[iz].Maj.P.pl );

	nsta_used = 0;
	for( ista=0; ista<nsta; ista++ ) 
	{
		if( ev[ista].iused == 1 ) nsta_used++;
	}

	fprintf( fp, "\n\tMaximum Azmuithal Gap=%.0f Distance to Nearest Station=%5.1f (km)\n",
		sol[iz].maxgap, sol[iz].dminkm );

	fprintf(fp, 
	  "\n\tNumber of Stations (D=Displacement/V=Velocity) Used=%d (defining only)\n", 
		nsta_used );


	fprintf(fp,"\t\t");
	kk=0;
	for( ista=0; ista<nsta; ista++ )
	{
		if( ev[ista].iused == 1 )
		{
			if( kk%4 == 0 ) fprintf(fp,"\n\t");
			grd_mo_type = ' ';
			if( ev[ista].grd_mo_type == DISPLACEMENT ) grd_mo_type = 'D';
			if( ev[ista].grd_mo_type == VELOCITY ) grd_mo_type = 'V';
			fprintf(fp, " %s.%s.%c", 
				grn[ista][iz].stnm, grn[ista][iz].net, grd_mo_type );
			kk++;
		}
	}	
	fprintf(fp, "\n\n\n");
                                                                                                                     
	ascii_mech( fp, 12, 12, 12,	
		sol[iz].Maj.P1.s, sol[iz].Maj.P1.d, sol[iz].Maj.P1.r,
		sol[iz].Maj.P2.s, sol[iz].Maj.P2.d, sol[iz].Maj.P2.r,
		sol[iz].Maj.P.pl, sol[iz].Maj.P.az,
		sol[iz].Maj.T.pl, sol[iz].Maj.T.az
	);
                                                                                                                     
	fprintf(fp, "\n\n");
	fprintf(fp, "All Stations defining and nondefining: \n" );
	fprintf(fp, "Station.Net \tDef \tDistance \tAzi    \tBazi  \tlo-f \thi-f vmodel\n" );
	fprintf(fp, "            \t    \t(km)     \t(deg)  \t(deg) \t(Hz) \t(Hz)    \n" );
	for( ista=0; ista<nsta; ista++ )
	{
		grd_mo_type = ' ';
		defining = ' ';

		if( ev[ista].iused == 1 ) 
			defining = 'Y';
		else 
			defining = 'N';

		if( ev[ista].grd_mo_type == DISPLACEMENT ) grd_mo_type = 'D';
		if( ev[ista].grd_mo_type == VELOCITY ) grd_mo_type = 'V';

		fprintf(fp, "%s.%s (%c) \t%c \t%8.1f  \t%3.0f  \t%3.0f  \t%.3f \t%.3f %s\n",
			grn[ista][iz].stnm, grn[ista][iz].net, grd_mo_type, defining,
			grn[ista][iz].rdist, grn[ista][iz].az, grn[ista][iz].baz,
			ev[ista].lf, ev[ista].hf, ev[ista].glib_filename );
	}

	fprintf( fp, "\n (V)-velocity (D)-Displacement\n\n" );
	fprintf( fp, "Author: %s\n", username );


	fprintf( fp, "Date: %4d/%02d/%02d %02d:%02d:%02.0f %s\n\n",
		gmtnow->year, gmtnow->month, gmtnow->mday, 
		gmtnow->hour, gmtnow->min, gmtnow->fsec, gmtnow->tzone );


	fprintf( fp, "Date: %4d/%02d/%02d %02d:%02d:%02.0f %s\n\n",
                localnow->year, localnow->month, localnow->mday,
                localnow->hour, localnow->min, localnow->fsec, localnow->tzone );

	fprintf( fp, "%s Version %s %s\n",
		progname, Version_Label, Version_Date );

	free(gmtnow);
	free(localnow);

}  /*** END OF SUBROUTINE     write_email_message() ****/


/*************************************************************************************************************/
/*** BEGIN SUBROUTINE : calc_max_azi_gap() - calculate the azimuthal gap from an unsorted arrary of floats ***/
/*************************************************************************************************************/

float calc_max_azi_gap( float *azi, int nsta )
{	
	int ista;
	float maxgap;
	void floatsort( float *, int );
	if( nsta <= 2 ) return 360;
	floatsort( azi, nsta );
	maxgap = azi[1] + azi[nsta] - azi[nsta-1];
	for( ista=1; ista<nsta; ista++ )
	{
		if( azi[ista+1] - azi[ista] > maxgap )
		{
			maxgap = azi[ista+1] - azi[ista];
		}
	}
	return (float)maxgap;
}

/***************************************************************************/
/**** BEGIN SUBROUTINE : calc_azi_gap_and_dmin() assumes iz == iz_best ***/
/***************************************************************************/

void calc_azi_gap_and_dmin( int iz, int nsta, EventInfo *ev, Solution *sol, Greens **g )
{
	int ista, nsta_on=0, itmp;
	float *azi, maxgap;
	float dmin;
	void floatsort( float *, int );
	
/*** compute dmin from only defining stations ***/

	dmin = 99999;
	for( ista=0; ista<nsta; ista++ )
	{
		if( ev[ista].iused == 1 )
		{
			if( g[ista][iz].rdist < dmin ) dmin = g[ista][iz].rdist;
			nsta_on++;
		}
	}

/*** compute azimuthal gap from only defining stations ***/

	azi = (float *)calloc( 2+nsta, sizeof(float) );
	itmp = 1;

	/* fprintf( stdout, "nsta=%d\n", nsta ); */

	for( ista=0; ista<nsta; ista++ )
	{
		/* fprintf( stdout, "ista=%d iused=%d az=%g\n", 
			ista, ev[ista].iused, g[ista][iz].az );
		*/

		if( ev[ista].iused == 1 )
		{
			azi[itmp] = g[ista][iz].az;
			itmp++;
		}
	}

	floatsort( azi, nsta_on );

	/* fprintf( stdout, "nsta_on=%d\n", nsta_on ); */
	/*
	for( ista=1; ista<=nsta_on; ista++ )
	{
		fprintf( stdout, "ista=%d azi=%g\n", ista, azi[ista] );
	}
	*/

	maxgap = 0;
	for( ista=1; ista<=nsta_on-1; ista++ )
	{
		if( azi[ista+1]-azi[ista] > maxgap ) maxgap = azi[ista+1]-azi[ista];
	}

	if( ((360-azi[nsta_on]) + azi[1] ) > maxgap ) 
		maxgap = ((360-azi[nsta_on]) + azi[1] );

	/* fprintf( stdout, "maxgap=%g\n", maxgap ); */

	sol[iz].nsta = nsta_on;
	sol[iz].maxgap = maxgap;
	sol[iz].dminkm = dmin;
	free(azi);

}  /*** END SUBROUTINE : calc_azi_gap_and_dmin() **/

/***************************************************************************/
/**** BEGIN SUBROUTINE : write_mterror() -  write out an ascii record *****/
/***************************************************************************/

void write_mterror( FILE *fp, int nz, Solution *sol, EventInfo *ev, Greens **g )
{
	int iz, ista=0;

	for( iz=0; iz<nz; iz++ )
	{
	  fprintf( fp, 
	    "%6.2f %6.2f Mij %7.3e %7.3e %7.3e %7.3e %7.3e %7.3e ERR %7.3e %7.3e %7.3e %7.3e %7.3e %7.3e %4d/%02d/%02d %02d:%02d:%05.2f\n",
			ev[ista].ts0,
			g[ista][iz].evdp,
			sol[iz].moment_tensor[1][1],
			sol[iz].moment_tensor[2][2],
			sol[iz].moment_tensor[1][2],
			sol[iz].moment_tensor[1][3],
			sol[iz].moment_tensor[2][3],
			sol[iz].moment_tensor[3][3],
			sol[iz].smxx,
			sol[iz].smyy,
			sol[iz].smxy,
			sol[iz].smxz,
			sol[iz].smyz,
			sol[iz].smzz,
			ev[ista].ot.year,
			ev[ista].ot.month,
			ev[ista].ot.mday,
			ev[ista].ot.hour,
			ev[ista].ot.min,
			sol[iz].ot
		);
	}

} /*** END SUBROUTINE write_mterror() ***/

/***************************************************************************/
/**** BEGIN SUBROUTINE : write_sol_rec() -  write out an ascii record *****/
/***************************************************************************/

void write_sol_rec( FILE *fp, int iz, int nsta, EventInfo *ev, Solution *sol, Greens **grn )
{
	int ista;
	int mtdegfree;
	float ot_tmp;

	ot_tmp = ev[0].ts0 * -1;  /*** this corrects for OT shift problem,  ***/
			/*** see similar correction to write_gmt_xy_values() : plotmech.out ***/

	if( sol[iz].mt_type == EXPLOSION  ) mtdegfree = 1;
	if( sol[iz].mt_type == DEVIATORIC ) mtdegfree = 5;
	if( sol[iz].mt_type == FULL_MT    ) mtdegfree = 6;

/***            1   2     3   4     5           6                  7      8     9   10    11    12   13    14    15    16    ***/
/***            x   y     z  ts0    yr/mo/da     hr:mn:sec           Mw    Mo   PDC PCLVD PISO       VRED  fit1  fit2   l2   ***/ 
fprintf(fp, "%.4f %.4f %5.1f %5.2f %4d/%02d/%02d %02d:%02d:%05.2f %5.2f %5.2e %3.0f %3.0f %3.0f %02d %6.2f %6.2f %6.2f %5.2e ",
                grn[0][iz].evlo,
                grn[0][iz].evla,
                grn[0][iz].evdp,
		ot_tmp, 
                ev[0].ot.year,
                ev[0].ot.month,
                ev[0].ot.mday,
                ev[0].ot.hour,
                ev[0].ot.min,
                sol[iz].ot,
                sol[iz].mw,
                sol[iz].dmoment,
                sol[iz].PDC,
		sol[iz].PCLVD,
		sol[iz].PISO,
		mtdegfree,
                sol[iz].var_red,
                sol[iz].total_fitness1,
                sol[iz].total_fitness2,
                sol[iz].l2norm_error
        );

                /*
                sol[iz].mrr,
                sol[iz].mtt,
                sol[iz].mff,
                sol[iz].mrt,
                sol[iz].mrf,
                sol[iz].mtf,
                */

/**
               sol[iz].moment_tensor[3][3],
                sol[iz].moment_tensor[1][1],
                sol[iz].moment_tensor[2][2],
                sol[iz].moment_tensor[1][3],
                -sol[iz].moment_tensor[2][3],
                -sol[iz].moment_tensor[1][2],
**/
	/***          17    18   19    20    21    22   23   24    25    26     ***/
	/***         mzz   mxx   myy   mxz  -myz  -mxy  exp  Teig  Peig  Beig   ***/
        fprintf(fp, "%6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %2d %6.2f %6.2f %6.2f ",
		sol[iz].mrr,
                sol[iz].mtt,
                sol[iz].mff,
                sol[iz].mrt,
                sol[iz].mrf,
                sol[iz].mtf,
                sol[iz].exponent,
		sol[iz].FullMT.T.ev,
		sol[iz].FullMT.P.ev,
		sol[iz].FullMT.B.ev
        );
 
	/***                     27           28                    29  ***/
	/***    double couple moment, clvd moment, and isotropic moment ***/
	fprintf( fp, "%6.3e %6.3e %6.3e ", sol[iz].Mdc, sol[iz].Mclvd, sol[iz].Miso );

/***          30   31    32    33     34   35    36    37     38    39    40    41     42   43    44    45   46   47     ***/
/***          eigt stkt  plgt  eigp  stkp  plgp  eign  stkn  plgn  str0  dip0   rak0  str1  dip1 rak1   eps   k   f_fac  ***/
fprintf(fp, "%6.2f %3.0f %3.0f %6.2f %3.0f %3.0f %6.2f %3.0f %3.0f %3.0f %2.0f %4.0f %3.0f %2.0f %4.0f %5.2f %5.2f %8.2f ",
                sol[iz].Maj.T.ev,
                sol[iz].Maj.T.az,
                sol[iz].Maj.T.pl,
                sol[iz].Maj.P.ev,
                sol[iz].Maj.P.az,
                sol[iz].Maj.P.pl,
		sol[iz].Maj.B.ev,
		sol[iz].Maj.B.az,
		sol[iz].Maj.B.pl,
                sol[iz].stk0, sol[iz].dip0, sol[iz].rak0,
                sol[iz].stk1, sol[iz].dip1, sol[iz].rak1,
                sol[iz].epsilon, sol[iz].k, sol[iz].f_factor
        );

        for( ista=0; ista<nsta; ista++ )
	{
		if( ev[ista].iused == 1 )
		{
			fprintf(fp, " %-4.4s", ev[ista].z.s.kstnm );
		}
	}
        fprintf(fp, "\n");

}  /*** write_sol_rec() ***/


/*****************************************************************/
/*** BEGIN SUBROUTINE : sdr_to_mt()                            ***/
/*****************************************************************/

void sdr_to_mt( float **mt, float str, float dip, float rak, float Mw, int verbose )
{
	float Mo, Mo0, pi, d2r, strr, dipr, rakr;
	Tensor M;
	pi = M_PI;
	d2r = pi/180.;
	Mo  = pow( 10.0, (1.5*(Mw+10.73)) );
	Mo0 = pow( 10.0, (1.5*10.73) );
	strr = str * d2r;
	dipr = dip * d2r;
	rakr = rak * d2r;

	M.xx=-(sin(dipr)*cos(rakr)*sin(2*strr)+sin(2*dipr)*sin(rakr)*sin(strr)*sin(strr));
	M.yy= (sin(dipr)*cos(rakr)*sin(2*strr)-sin(2*dipr)*sin(rakr)*cos(strr)*cos(strr));
	M.zz= (sin(2*dipr)*sin(rakr));
	M.xy= (sin(dipr)*cos(rakr)*cos(2*strr)+0.5*sin(2*dipr)*sin(rakr)*sin(2*strr));
	M.xz=-(cos(dipr)*cos(rakr)*cos(strr)+cos(2*dipr)*sin(rakr)*sin(strr));
	M.yz=-(cos(dipr)*cos(rakr)*sin(strr)-cos(2*dipr)*sin(rakr)*cos(strr));
	M.yx = M.xy;
	M.zx = M.xz;
	M.zy = M.yz;
	if(verbose)
	{
		printf("%s:sdr_to_mt(): coef Mxx=%g Myy=%g Mzz=%g Mxy=%g Mxz=%g Myz=%g\n",
			progname, M.xx, M.yy, M.zz, M.xy, M.xz, M.yz );
	}

        mt[1][1] = M.xx * Mo/Mo0;
        mt[1][2] = M.xy * Mo/Mo0;
        mt[1][3] = M.xz * Mo/Mo0;
        mt[2][1] = M.yx * Mo/Mo0;
        mt[2][2] = M.yy * Mo/Mo0;
        mt[2][3] = M.yz * Mo/Mo0;
        mt[3][1] = M.zx * Mo/Mo0;
        mt[3][2] = M.zy * Mo/Mo0;
        mt[3][3] = M.zz * Mo/Mo0;

	if(verbose)
	{
		printf( "%s: sdr_to_mt: moment tensor\n", progname );
		printf( "%s: sdr_to_mt:\t %6.2e %6.2e %6.2e\n", progname, mt[1][1], mt[1][2], mt[1][3] );
		printf( "%s: sdr_to_mt:\t %6.2e %6.2e %6.2e\n", progname, mt[2][1], mt[2][2], mt[2][3] );
		printf( "%s: sdr_to_mt:\t %6.2e %6.2e %6.2e\n", progname, mt[3][1], mt[3][2], mt[3][3] );
	}
}

float fsign( float x )
{
	float y = 1;
	if( x >= 0 ) return y;
	if( x <  0 ) return -y;
	return 1;
}

void floatsort( float *x, int n )
{
	int k,j,m;
	float hold;
	n=n+1;
	for( k=1; k <=n-2; k++ )
	{
		m = k;
		for( j = k + 1; j <= n - 1; j++ )
		{
			if( x[j] < x[m] ) m = j;
		}
		hold = x[m];
		x[m] = x[k];
		x[k] = hold;
	}
	return;
}

int roundflt2int( float x )
{
	if( (x-floor(x)) > 0.5 ) return ceil(x);
	if( (x-floor(x)) < 0.5 ) return floor(x);
	return ceil(x);
}

void normalize_moment_tensor( MomentTensor *Ma, MomentTensor *Mn, float Mo, int verbose )
{
	int i, j;

	if(verbose)
	{
		fprintf( stdout,
		  "%s: %s: %s: calling normalize_moment_tensor: \n",
			progname, __FILE__, __func__ );
	}

/*** compute the total scaler seismic moment ***/

/*********************************************************************************************/
/*** Bowers and Hudson (1999) shows that Silver and Jordan (1982) method only              ***/
/*** works for pure DC not pure CLVD 0.86 too small and pure EX 1.22 too big.              ***/
/*** Bowers and Hudson (1999-Oct) BSSA v89n5pp1390-1394 suggest Mtotal = Miso + Mdev where ***/
/***	eigenvalues: m1, m2, m3                                                            ***/
/***	Miso = fabs| ( m1 + m2 + m3 ) / 3 |                                                ***/
/***	Mdev = | max(m) - Miso |                                                           ***/
/*********************************************************************************************/
/***
	for( i = 1; i <= 3; i++ )
	{
		for( j = 1; j <= 3; j++ )
		{
			Ma->moment += (Ma->mt[i][j] * Ma->mt[i][j]);
		}
	}
	Ma->moment  = sqrt( Ma->moment / 2.0 ) * base_moment;
***/
	Ma->moment = Mo;
	Ma->Mw      = log10( Ma->moment )/1.5 - 10.73;
	Ma->expon    = (int)log10( Ma->moment );
	Ma->abcassa = Ma->moment/pow(10.0, Ma->expon );

/*** normalized moment tensor ***/

	Mn->xx = Ma->xx / pow( 10.0, Ma->expon );
        Mn->xy = Ma->xy / pow( 10.0, Ma->expon );
        Mn->xz = Ma->xz / pow( 10.0, Ma->expon );
        Mn->yx = Ma->yx / pow( 10.0, Ma->expon );
        Mn->yy = Ma->yy / pow( 10.0, Ma->expon );
        Mn->yz = Ma->yz / pow( 10.0, Ma->expon );
        Mn->zx = Ma->zx / pow( 10.0, Ma->expon );
        Mn->zy = Ma->zy / pow( 10.0, Ma->expon );
        Mn->zz = Ma->zz / pow( 10.0, Ma->expon );

	Mn->mt[1][1] = Mn->xx;
        Mn->mt[1][2] = Mn->xy;
        Mn->mt[1][3] = Mn->xz;
        Mn->mt[2][1] = Mn->yx;
        Mn->mt[2][2] = Mn->yy;
        Mn->mt[2][3] = Mn->yz;
        Mn->mt[3][1] = Mn->zx;
        Mn->mt[3][2] = Mn->zy;
        Mn->mt[3][3] = Mn->zz;

/************************************************************/
/*** Moment tensor components  T=theta F=phi              ***/
/***    MXX  MXY  MXZ or if using Harvard  MTT -MTF  MRT  ***/
/***    MYX  MYY  MYZ     notation        -MTF  MFF -MRF  ***/
/***    MZX  MZY  MZZ                      MRT -MRF  MRR  ***/
/************************************************************/
/**********************************************************************************/
/***               mzz, mxx, myy, mxz, -myz, -mxy, exp                          ***/
/***  X, Y, depth, mrr, mtt, mff, mrt,  mrf,  mtf, exp, newX, newY, event_title ***/
/**********************************************************************************/
	Mn->rr =  Mn->zz;
	Mn->tt =  Mn->xx;
	Mn->ff =  Mn->yy;
	Mn->rt =  Mn->xz;
	Mn->rf = -Mn->yz;
	Mn->tf = -Mn->xy;

	Mn->moment = Ma->moment;
	Mn->Mw     = Ma->Mw;
	Mn->expon  = Ma->expon;
	Mn->abcassa= Ma->abcassa;
	
	if(verbose)
	{
		fprintf( stdout, "%s: %s: %s: Normalized Moment Tensor:\n", progname, __FILE__, __func__ );
		fprintf( stdout, "\t %5.2f %5.2f %5.2f\n", Mn->xx, Mn->xy, Mn->xz );
		fprintf( stdout, "\t %5.2f %5.2f %5.2f\n", Mn->yx, Mn->yy, Mn->yz );
		fprintf( stdout, "\t %5.2f %5.2f %5.2f\n", Mn->zx, Mn->zy, Mn->zz );
		fprintf( stdout, "\t Mo=%5.2e Mw=%5.2f Mo=%5.2f X 10^%2d Dyne cm\n",
			Ma->moment, Ma->Mw, Ma->abcassa, Ma->expon );
		fprintf( stdout, "\n" );
	}
}

void set_moment_tensor( MomentTensor *Ma, float *x, int idf, int verbose )
{
	int i, j;

	if(verbose)
	  fprintf(stdout, "%s: mtinv_subs.c: set_moment_tensor():\n", progname );

/*** set the ABS moment tensor from x vector ***/
	Ma->xx = x[1];
	Ma->yy = x[2];
	Ma->xy = x[3];
	Ma->xz = x[4];
	Ma->yz = x[5];
	Ma->zz = x[6];

	if( idf == 5 ) Ma->zz = -( Ma->xx + Ma->yy );
	Ma->yx = Ma->xy;
	Ma->zx = Ma->xz;
	Ma->zy = Ma->yz;
	
	Ma->mt[1][1] = Ma->xx;
	Ma->mt[1][2] = Ma->xy;
	Ma->mt[1][3] = Ma->xz;
	Ma->mt[2][1] = Ma->yx;
	Ma->mt[2][2] = Ma->yy;
	Ma->mt[2][3] = Ma->yz;
	Ma->mt[3][1] = Ma->zx;
	Ma->mt[3][2] = Ma->zy;
	Ma->mt[3][3] = Ma->zz;
}

/*****************************************************************************/
/*** calculate the eigenvalues and eigenvectors of the full moment tensors ***/
/*** also calculates the isotropic component, kiso and epsilon values      ***/
/*****************************************************************************/

void mt2eig( MomentTensor Ma, Solution *sol, int iz, int verbose )
{
	float **z1, **z2;
	float *eval1, *eval2, *evp, *evd;
	int i, j, nx=3;
	void tred2( float **, int, float *, float * );
	void tqli( float *, float *, int, float ** );

	if(verbose)
	  fprintf(stdout, "%s: mtinv_subs.c: mt2eig(): \n", progname );

/*** allocate memory ***/
	eval1 = vector( 0, 4 );
	eval2 = vector( 0, 4 );
	evp = vector( 0, 4 );
	evd = vector( 0, 4 );
	z1 = matrix( 0, 4, 0, 4 );
	z2 = matrix( 0, 4, 0, 4 );

/*** initalize arrays ***/
	for( i=1; i<=nx; i++ )
	{
		for( j=1; j<=nx; j++)
		{
			z1[i][j] = Ma.mt[i][j];
			z2[i][j] = z1[i][j];
		}
	}

/*** do the eigenvalue/eigenvector calculations ***/

	tred2( z2, nx, eval1, eval2 );
	tqli( eval1, eval2, nx, z2 );

/*** write out ***/

	if(verbose)
	{
		for( i=1; i<=nx; i++ )
		{
			/*** eigenvalues ***/
			fprintf( stdout, "%s: mtinv_subs.c: mt2eig(): eval1=%6.2e z2=", 
				progname, eval1[i]);
			/*** eigenvectors ***/
			for( j=1; j<=nx; j++)
			{
				printf("%6.2f ", z2[i][j] );
			}
			fprintf(stdout, "\n" );
		}
	}

	sol[iz].FullMT.eig[1].val = eval1[1];
	sol[iz].FullMT.eig[2].val = eval1[2];
	sol[iz].FullMT.eig[3].val = eval1[3];

	sol[iz].FullMT.eig[1].vec[1] = z2[1][1];
	sol[iz].FullMT.eig[1].vec[2] = z2[2][1];
	sol[iz].FullMT.eig[1].vec[3] = z2[3][1];

	sol[iz].FullMT.eig[2].vec[1] = z2[1][2];
	sol[iz].FullMT.eig[2].vec[2] = z2[2][2];
	sol[iz].FullMT.eig[2].vec[3] = z2[3][2];

	sol[iz].FullMT.eig[3].vec[1] = z2[1][3];
	sol[iz].FullMT.eig[3].vec[2] = z2[2][3];
	sol[iz].FullMT.eig[3].vec[3] = z2[3][3];

}  /*** end of mt2eig() ***/


/*******************************************************************************/
/*** void eig2iso( EventInfo *ev, Solution *sol, int iz, int verbose ) ***/
/*******************************************************************************/
void eig2iso( Solution *sol, int iz, int verbose )
{
	float etmp[5];
	float sign[5];
	int nx=3, i, j;
	int index[5];
	float iso;
	float x_vec[8];
	MomentTensor Ma, Mn;
	static float zero = 1.0E-09;

	void set_moment_tensor( MomentTensor *, float *, int, int );
	void normalize_moment_tensor( MomentTensor *Ma, MomentTensor *Mn, float Mo, int verbose );
	void floatsort( float *, int );
	float fsign( float );
	void indexx( int, float *, int * );
	void sac_absmax( float *, int, float *);

/** set the eigenvalues of the full moment tensor ***/
	
	if(verbose) 
	{
	  fprintf(stdout, "%s: mtinv_subs.c: eig2iso(): Full moment tensor:\n", progname );
	}

	for( i=1; i<=nx; i++ )
	{
		if(verbose)	
		  fprintf( stdout, "%s: mtinv_subs.c: eig2iso(): val=%6.2e vec=", 
			progname, 
			sol[iz].FullMT.eig[i].val );

		etmp[i] = sol[iz].FullMT.eig[i].val;

		if( verbose )
		{
			for( j=1; j<=nx; j++ )
			{
				fprintf( stdout, "%6.2e ", sol[iz].FullMT.eig[i].vec[j] );
			}
			fprintf( stdout, "\n");
		}
	}

/*************************************************************************/
/*** sort the princple moments index[1] points to the largest negative ***/
/***                           index[3] points to the largest positive ***/
/***                           index[2] points to the smallest         ***/
/*************************************************************************/
	indexx( 3, &etmp[0], &index[0] );
	sol[iz].FullMT.B.ev = sol[iz].FullMT.eig[index[2]].val;
	sol[iz].FullMT.P.ev = sol[iz].FullMT.eig[index[1]].val;
	sol[iz].FullMT.T.ev = sol[iz].FullMT.eig[index[3]].val;

	if(verbose)
	{
		fprintf( stdout, "%s: mtinv_subs.c: eig2iso(): FULL sorted eigenvalues = ", progname );
		fprintf( stdout, "%e <= %e <= %e.\n",
			sol[iz].FullMT.eig[index[1]].val,
			sol[iz].FullMT.eig[index[2]].val,
			sol[iz].FullMT.eig[index[3]].val );
	}

/*** isotropic moment ***/

	iso = ( sol[iz].FullMT.eig[1].val + 
		sol[iz].FullMT.eig[2].val +
		sol[iz].FullMT.eig[3].val )/3.0;

	sol[iz].Miso = fabs( iso );

	if(verbose)
	{
	  fprintf( stdout, "%s: %s: %s: raw isotropic component = %e Miso = %e\n", 
		progname, __FILE__, __func__, iso, sol[iz].Miso );

	  fprintf( stdout, "%s: %s: %s: M xx=%e yy=%e xy=%e xz=%e yz=%e zz=%e\n",
		progname, 
		__FILE__, __func__,
		sol[iz].moment_tensor[1][1],
		sol[iz].moment_tensor[2][2],
		sol[iz].moment_tensor[1][2],
		sol[iz].moment_tensor[1][3],
		sol[iz].moment_tensor[2][3],
		sol[iz].moment_tensor[3][3] );
	}

/******************************************************************/
/*** form the moment tensor from solution vector x              ***/
/*** col1-Mxx, col2-Myy, col3-Mxy, col4-Mxz, col5-Myz, col6-Mzz ***/
/******************************************************************/
	if( fabs(iso) > 1.0E+12 ) 
	{
		if(verbose)
		{
        		fprintf( stdout,
		  "%s: mtinv_subs.c: eig2iso(): fabs(iso)=%e > 1.0E+12\n", 
				progname, 
				fabs(iso) ); 
		}

		x_vec[1] = sol[iz].moment_tensor[1][1] - iso;
		x_vec[2] = sol[iz].moment_tensor[2][2] - iso;
		x_vec[3] = sol[iz].moment_tensor[1][2];
		x_vec[4] = sol[iz].moment_tensor[1][3];
		x_vec[5] = sol[iz].moment_tensor[2][3];
		x_vec[6] = sol[iz].moment_tensor[3][3] - iso;
	}
	else
	{
		x_vec[1] = sol[iz].moment_tensor[1][1];
		x_vec[2] = sol[iz].moment_tensor[2][2];
		x_vec[3] = sol[iz].moment_tensor[1][2];
		x_vec[4] = sol[iz].moment_tensor[1][3];
		x_vec[5] = sol[iz].moment_tensor[2][3];
		x_vec[6] = sol[iz].moment_tensor[3][3];
	}

/*** calculate deviatoric moment ***/
/***********************************************************************************/
/*** remove the isotropic part from the principle moments (eigenvalues) from the ***/
/*** from the full moment tenosr to the the deviatoric principle moments         ***/
/*** also copy over the eigenvalues                                              ***/
/***********************************************************************************/

	if(verbose)
	  fprintf( stdout,
		"%s: mtinv_subs.c: eig2iso(): deviatoric moment tensor\n",
			progname );

	for( i=1; i<=3; i++ )
	{
		sol[iz].Dev.eig[i].val = sol[iz].FullMT.eig[i].val - iso;
		etmp[i] = sol[iz].Dev.eig[i].val;

		if(verbose)
		{
			fprintf( stdout,
			  "%s: mtinv_subs.c: eig2iso():\t val=%6.2e vec=",
				progname, sol[iz].Dev.eig[i].val );
		}

		for( j=1; j<=3; j++ )
		{
			sol[iz].Dev.eig[i].vec[j] = sol[iz].FullMT.eig[i].vec[j];
			if(verbose)
			  fprintf( stdout, "%6.2f ", sol[iz].Dev.eig[i].vec[j] );
		}
		if(verbose) fprintf( stdout, "\n");
	}
	
/***********************************************************************************/
/*** DEVIATORIC (isotropic removed                                               ***/
/*** sort the principle moments index[1] points to the largest negative eig      ***/
/***                            index[2] points to the smallest near zero        ***/
/***                            index[3] points to the largest positive eig      ***/
/***********************************************************************************/

	indexx(3,&etmp[0],&index[0]);
	if(verbose)
	{
		fprintf( stdout, "%s: mtinv_subs.c: eig2iso(): deviatoric sorted eigenvalues: ", 
			progname );
		fprintf( stdout, "%s: mtinv_subs.c: eig2iso(): %g %g %g\n",
			progname, 
			sol[iz].Dev.eig[index[1]].val,
			sol[iz].Dev.eig[index[2]].val,
			sol[iz].Dev.eig[index[3]].val );
	}

/*************************************************/
/*** sort the absolute value principle moments ***/
/*** etmp[1] is the smallest                   ***/
/*** etmp[3] is the largest                    ***/
/*** pure isotropic source will have 0 deviatoric     ***/
/*** which weaks havoc when dividing by zero          ***/
/********************************************************/

	for( i=1; i<=3; i++ )
	{
		if( fabs( sol[iz].Dev.eig[i].val ) <= zero ) sol[iz].Dev.eig[i].val = zero;

		etmp[i] = fabs( sol[iz].Dev.eig[i].val );
	}
	floatsort( &etmp[0], 3 );

/********************************************************/
/*** copy the sign of the eigenvalue/principle moment ***/
/********************************************************/

	if( fabs(sol[iz].Dev.eig[1].val ) == etmp[1] ) sign[1] = fsign( sol[iz].Dev.eig[1].val );
	if( fabs(sol[iz].Dev.eig[1].val ) == etmp[2] ) sign[2] = fsign( sol[iz].Dev.eig[1].val );
	if( fabs(sol[iz].Dev.eig[1].val ) == etmp[3] ) sign[3] = fsign( sol[iz].Dev.eig[1].val );
	if( fabs(sol[iz].Dev.eig[2].val ) == etmp[1] ) sign[1] = fsign( sol[iz].Dev.eig[2].val );
	if( fabs(sol[iz].Dev.eig[2].val ) == etmp[2] ) sign[2] = fsign( sol[iz].Dev.eig[2].val );
	if( fabs(sol[iz].Dev.eig[2].val ) == etmp[3] ) sign[3] = fsign( sol[iz].Dev.eig[2].val );
	if( fabs(sol[iz].Dev.eig[3].val ) == etmp[1] ) sign[1] = fsign( sol[iz].Dev.eig[3].val );
	if( fabs(sol[iz].Dev.eig[3].val ) == etmp[2] ) sign[2] = fsign( sol[iz].Dev.eig[3].val );
	if( fabs(sol[iz].Dev.eig[3].val ) == etmp[3] ) sign[3] = fsign( sol[iz].Dev.eig[3].val );

	if(verbose) 
	{
	 fprintf( stdout, "%s: mtinv_subs.c: eig2iso(): Absolute Value Eigenvalue Sorted : E1 > E2 > E3 \n", progname );
	 fprintf( stdout, "%s: mtinv_subs.c: eig2iso(): |Eigval|\t %e %e %e\n", progname, etmp[1], etmp[2], etmp[3] );
	 fprintf( stdout, "%s: mtinv_subs.c: eig2iso(): signs \t %g %g %g\n", progname, sign[1], sign[2], sign[3] );
	}

/********************************************************/
/*** k value source type Hudson et al 1989            ***/
/********************************************************/

/* ************ New way, as per Gordon Kraft (JIRN, 19 Feb 2014) ********* */
/* Sort the absolute values of the deviatoric components and sort          */
/* a copy of the deviatoric components in the corresponding manner.        */
/* Because the first arguments is ignored, pass the address of the         */
/* first index, not the zeroth index.   The logic below is from lam2Tk.m   */
/* lam2Tk.m from Carl Tape (12/2012 matlab code):                          */
/* http://compearth.googlecode.com/svn/trunk/momenttensor/matlab/lam2Tk.m  */
/* Using the lam2Tk.m convenction, T should then be defined as -2*epsilon. */
/*****************************************************************************/

	void flt_sort_asc(int npts, float *fvec_a, float *fvec_b);

	float  *tmp_abs_dev_eigen   = (float *)calloc( 4, sizeof(float) );
	float  *tmp_dev_eigen       = (float *)calloc( 4, sizeof(float) );

	float   max_abs_dev_eigen   = 0.0f;
	float   min_abs_dev_eigen   = 0.0f;
	float   min_dev_eigen_corsp = 0.0f;

	tmp_dev_eigen[0]     = 0.0f;
	tmp_dev_eigen[1]     = sol[iz].Dev.eig[1].val;
	tmp_dev_eigen[2]     = sol[iz].Dev.eig[2].val;
	tmp_dev_eigen[3]     = sol[iz].Dev.eig[3].val;

	tmp_abs_dev_eigen[0] = 0.0f;
	tmp_abs_dev_eigen[1] = fabs( sol[iz].Dev.eig[1].val );
	tmp_abs_dev_eigen[2] = fabs( sol[iz].Dev.eig[2].val );
	tmp_abs_dev_eigen[3] = fabs( sol[iz].Dev.eig[3].val );

/**********************************************************************************/
/* Sort the abs(eigen) values, but keep a one-to-one match with the signed values */
/* Need to pass index one, since flt_sort_asc uses standard C array indexing.     */
/**********************************************************************************/

	flt_sort_asc( 3, &(tmp_abs_dev_eigen[1]), &(tmp_dev_eigen[1]) );


/* Print to screen the sorted (ascending) abs( eigen[] ) values */
	if(verbose)
	{
	  fprintf( stdout,
		"%s: %s: %s(): Dev eigen values, in ascending order w.r.t their abs values: %6.2e %6.2e %6.2e\n",
		progname,
		__FILE__,
		__func__,
		tmp_dev_eigen[1], tmp_dev_eigen[2], tmp_dev_eigen[3]);
	}

	min_abs_dev_eigen   = tmp_abs_dev_eigen[1]; /* min( abs(eigen[]) )                   */
	min_dev_eigen_corsp = tmp_dev_eigen[1];     /* Signed version of the min(abs(eigen)) */
	max_abs_dev_eigen   = tmp_abs_dev_eigen[3]; /* max( abs(eigen[]) )                   */

	if(verbose)
	{
	  fprintf( stdout, "%s: %s: %s: min_abs_dev_eigen=%g min_dev_eigen_corsp=%g max_abs_dev_eigen=%g\n",
		progname, __FILE__, __func__, min_abs_dev_eigen, min_dev_eigen_corsp, max_abs_dev_eigen );
	}

	sol[iz].Mdev = max_abs_dev_eigen;

/* Get the values for k and epsilon */
	sol[iz].k           = iso / ( fabs(iso) + max_abs_dev_eigen );

	if( fabs(min_dev_eigen_corsp) <= zero || fabs(max_abs_dev_eigen) <= zero )
	{
		sol[iz].epsilon  = 0;
	}
	else
	{
		sol[iz].epsilon  = (-1.0) * min_dev_eigen_corsp / max_abs_dev_eigen;
	}
	
/* Free the allocated memory */
	if( tmp_abs_dev_eigen )  free( tmp_abs_dev_eigen );
	if( tmp_dev_eigen )      free( tmp_dev_eigen );

/********************************************************/
/*** epsilon value source type ( Hudson et al 1989 ) ***/
/********************************************************/

/**********************************************************************************/
/* Gene, It was much cleaner to keep the computation of the k and epsilon   */
/* together, since they new way reuse the same sorted arrays.  For viewing  */
/* the new of computing epsilon, go up ~30 lines of code JIRN, 19 Feb 2014  */
/**********************************************************************************/

	sol[iz].Mclvd = 2.0 * fabs( sol[iz].epsilon ) * sol[iz].Mdev;

	sol[iz].Mdc   = sol[iz].Mdev - sol[iz].Mclvd;

	sol[iz].Mtotal = sol[iz].Mdev + sol[iz].Miso;

/***********************************************************/
/*** tectonic release f-factor (Toksoz and Kehrner 1972) ***/
/***********************************************************/

	if( fabs( sol[iz].Miso ) != 0 )
	{
	  sol[iz].f_factor = fabs( sol[iz].Mdev ) / fabs( sol[iz].Miso );
	}
	else
	{
          sol[iz].f_factor = fabs( sol[iz].Mdev ) / zero;
	}

	if( sol[iz].mt_type != FULL_MT )
	{
		sol[iz].f_factor = 0;
		sol[iz].k = 0;
	}

	if( verbose )
	{
	  fprintf( stdout, 
"%s: %s: %s: iz=%d z=%g k=%5.3f epsilon=%5.3f Mdev=%g F-factor=%g (Mdc=%g Mclvd=%g Miso=%g sum=%g) Mtotal=%g(Miso+Mdev)\n",
		progname,
		__FILE__,
		__func__,
		iz,
		sol[iz].evdp,
		sol[iz].k, 
		sol[iz].epsilon,
		sol[iz].Mdev,
		sol[iz].f_factor,
		sol[iz].Mdc,
		sol[iz].Mclvd,
		sol[iz].Miso,
		(sol[iz].Mdc+sol[iz].Mclvd+sol[iz].Miso),
		sol[iz].Mtotal );
	}

	sol[iz].bev = sign[1] * etmp[1];
	if( sign[3] >= 0 )
	{
		sol[iz].tev = sign[3] * etmp[3];
		sol[iz].pev = sign[2] * etmp[2];
	}
	if( sign[3] < 0 )
	{
		sol[i].pev = sign[3] * etmp[3];
		sol[i].tev = sign[2] * etmp[2];
	}

/*** use the k=%iso and 2*epsilon=%clvd ***/

	sol[iz].fIso  = sol[iz].k;
	sol[iz].fclvd = 2 * fabs( sol[iz].epsilon );
	sol[iz].PDEV  = 100*( sol[iz].Mdev  / sol[iz].Mtotal );
	sol[iz].PISO  = 100*( sol[iz].Miso  / sol[iz].Mtotal );
	sol[iz].PCLVD = 100*( sol[iz].Mclvd / sol[iz].Mtotal );
	sol[iz].PDC   = 100*( sol[iz].Mdc   / sol[iz].Mtotal );
	
	if(verbose) 
	{
	  fprintf( stdout,
		"%s: %s: %s: PDEV=%4.1f PISO=%4.1f PCLVD=%4.1f PDC=%4.1f\n",
			progname,
			__FILE__,
			__func__,
			sol[iz].PDEV,
			sol[iz].PISO,
			sol[iz].PCLVD,
			sol[iz].PDC );
	}

} /*** end of eig2iso ***/

/*******************************************************************/
/*** calculate the eigenvalues and eigenvectors for the Major DC ***/
/*** convert eig to P,T,N axes and nodal planes str,dip,rak      ***/
/*******************************************************************/

void Eig2MajorDC( Solution *sol, int iz, int verbose )
{
	int index[4], i, j;
	float evd[4];
	void indexx( int, float *, int * );
	void trans1( Solution *, int, int *, int );
	void tpdss( Solution *, int );

	if(verbose) fprintf(stdout, "%s: calling Eig2MajorDC\n", progname );

/*** transfer the deviatoric eigenvalues to majordc ***/
	for( i=1; i<=3; i++ )
	{
		sol[iz].Maj.eig[i].val = sol[iz].Dev.eig[i].val;
		for( j=1; j<=3; j++ )
			sol[iz].Maj.eig[i].vec[j] = sol[iz].Dev.eig[i].vec[j];
	}

/*** sort and make majordc ***/
	for(i=1; i<=3; i++) evd[i]=sol[iz].Dev.eig[i].val;
	indexx( 3, &evd[0], &index[0] );
	sol[iz].Maj.eig[index[2]].val = +0.0;
	sol[iz].Maj.eig[index[1]].val = -1.0 * sol[iz].Maj.eig[index[3]].val;
	sol[iz].Maj.eig[index[3]].val = +1.0 * sol[iz].Maj.eig[index[3]].val;

/*** find ptn axis and nodal planes ***/
	trans1( sol, iz, index, 0 );
	tpdss( sol, iz );

	sprintf( sol[iz].Maj.T.ktype, "T " );
	sol[iz].Maj.T.az = sol[iz].stkt;
	sol[iz].Maj.T.pl = sol[iz].plnt;
	sol[iz].Maj.T.ev = sol[iz].Maj.eig[index[3]].val;

	sprintf( sol[iz].Maj.P.ktype, "P " );
	sol[iz].Maj.P.az = sol[iz].stkp;
	sol[iz].Maj.P.pl = sol[iz].plnp;
	sol[iz].Maj.P.ev = sol[iz].Maj.eig[index[1]].val;

	sprintf( sol[iz].Maj.B.ktype, "B " );
	sol[iz].Maj.B.az = sol[iz].stkn;
	sol[iz].Maj.B.pl = sol[iz].plnn;
	sol[iz].Maj.B.ev = sol[iz].Maj.eig[index[2]].val;

	sol[iz].Maj.P1.s = sol[iz].stk0;
	sol[iz].Maj.P1.d = sol[iz].dip0;
	sol[iz].Maj.P1.r = sol[iz].rak0;
	sol[iz].Maj.P2.s = sol[iz].stk1;
	sol[iz].Maj.P2.d = sol[iz].dip1;
	sol[iz].Maj.P2.r = sol[iz].rak1;

	if(verbose)
	{
	  fprintf(stdout, "Eig2MajorDC(): iz=%d Major Double Couple:\n", iz );
	  for( i=1; i<=3; i++ )
	  {
		printf( "Eig2MajorDC: \t val=%6.2e vec=", sol[iz].Maj.eig[index[i]].val );
		for( j=1; j<=3; j++ )
			printf( "%5.2f ", sol[iz].Maj.eig[index[i]].vec[j]);
		printf("\n");
	  }
	  fprintf(stdout, "Eig2MajorDC:\t T-axis trend=%3.0f pl=%3.0f Ev=%6.2e\n", sol[iz].Maj.T.az, sol[iz].Maj.T.pl, sol[iz].Maj.T.ev );
	  fprintf(stdout, "Eig2MajorDC:\t P-axis trend=%3.0f pl=%3.0f Ev=%6.2e\n", sol[iz].Maj.P.az, sol[iz].Maj.P.pl, sol[iz].Maj.P.ev );
	  fprintf(stdout, "Eig2MajorDC:\t B-axis trend=%3.0f pl=%3.0f Ev=%6.2e\n", sol[iz].Maj.B.az, sol[iz].Maj.B.pl, sol[iz].Maj.B.ev );
	  fprintf(stdout, "Eig2MajorDC:\t NP1: stk=%3.0f dip=%3.0f rak=%3.0f\n", sol[iz].Maj.P1.s, sol[iz].Maj.P1.d, sol[iz].Maj.P1.r );
	  fprintf(stdout, "Eig2MajorDC:\t NP2: stk=%3.0f dip=%3.0f rak=%3.0f\n", sol[iz].Maj.P2.s, sol[iz].Maj.P2.d, sol[iz].Maj.P2.r );
	  fprintf(stdout, "\n");
	}

} /*** End of Eig2MajorDC ***/

/******************************************************************/
/*** calculate the eigenvalues and eigenvectors of the MINOR DC ***/
/******************************************************************/
void Eig2MinorDC( Solution *sol, int iz, int verbose )
{
	int index[4], i, j;
	float evd[4];
	void indexx( int, float *, int * );
	void trans1( Solution *, int, int *, int );
	void tpdss( Solution *, int );

	if(verbose) fprintf(stdout, "%s: calling Eig2MinorDC\n", progname );

	for( i=1; i<=3; i++ )
	{
		sol[iz].Min.eig[i].val = sol[iz].Dev.eig[i].val;
		for( j=1; j<=3; j++ )
			sol[iz].Min.eig[i].vec[j] = sol[iz].Dev.eig[i].vec[j];
	}

	for(i=1; i<=3; i++) evd[i] = fabs(sol[iz].Dev.eig[i].val);
	indexx( 3, &evd[0], &index[0] );

	sol[iz].Min.eig[index[1]].val = +1.0 * sol[iz].Dev.eig[index[1]].val;
	sol[iz].Min.eig[index[2]].val = -1.0 * sol[iz].Dev.eig[index[1]].val;
	sol[iz].Min.eig[index[3]].val = +0.0;

	trans1( sol, iz, &index[0], 1 );
	tpdss( sol, iz );

	sprintf( sol[iz].Min.T.ktype, "T " );
	sol[iz].Min.T.az = sol[iz].stkt;
	sol[iz].Min.T.pl = sol[iz].plnt;
	sol[iz].Min.T.ev = sol[iz].Min.eig[index[1]].val;

	sprintf( sol[iz].Min.P.ktype, "P " );
	sol[iz].Min.P.az = sol[iz].stkp;
	sol[iz].Min.P.pl = sol[iz].plnp;
	sol[iz].Min.P.ev = sol[iz].Min.eig[index[2]].val;

	sprintf( sol[iz].Min.B.ktype, "B " );
	sol[iz].Min.B.az = sol[iz].stkn;
	sol[iz].Min.B.pl = sol[iz].plnn;
	sol[iz].Min.B.ev = sol[iz].Min.eig[index[3]].val;

	sol[iz].Min.P1.s = sol[iz].stk0;
	sol[iz].Min.P1.d = sol[iz].dip0;
	sol[iz].Min.P1.r = sol[iz].rak0;
	sol[iz].Min.P2.s = sol[iz].stk1;
	sol[iz].Min.P2.d = sol[iz].dip1;
	sol[iz].Min.P2.r = sol[iz].rak1;

	if(verbose)
	{
	  fprintf(stdout, "Eig2MinorDC(): Minor Double Couple:\n");
	  for( i=1; i<=3; i++ )
	  {
		printf( "Eig2MinorDC:\t val=%6.2e vec=", sol[iz].Min.eig[index[i]].val );
		for( j=1; j<=3; j++ )
			printf( "%6.3f ", sol[iz].Min.eig[i].vec[j]);
		printf("\n");
	  }
	  fprintf(stdout, "Eig2MinorDC:\t T-axis trend=%3.0f pl=%3.0f Ev=%6.2e\n", sol[iz].Min.T.az, sol[iz].Min.T.pl, sol[iz].Min.T.ev );
	  fprintf(stdout, "Eig2MinorDC:\t P-axis trend=%3.0f pl=%3.0f Ev=%6.2e\n", sol[iz].Min.P.az, sol[iz].Min.P.pl, sol[iz].Min.P.ev );
	  fprintf(stdout, "Eig2MinorDC:\t B-axis trend=%3.0f pl=%3.0f Ev=%6.2e\n", sol[iz].Min.B.az, sol[iz].Min.B.pl, sol[iz].Min.B.ev );
	  fprintf(stdout, "Eig2MinorDC:\t NP1: stk=%3.0f dip=%3.0f rak=%3.0f\n", sol[iz].Min.P1.s, sol[iz].Min.P1.d, sol[iz].Min.P1.r );
	  fprintf(stdout, "Eig2MinorDC:\t NP2: stk=%3.0f dip=%3.0f rak=%3.0f\n", sol[iz].Min.P2.s, sol[iz].Min.P2.d, sol[iz].Min.P2.r );
	  fprintf(stdout, "\n");
	}

} /*** end of Eig2MinorDC ***/


/*********************************************************/
/*** shift each station by its own specific time shift ***/
/*********************************************************/

void time_shift( EventInfo *ev, int nsta, int verbose )
{
	float dt, *ver, *rad, *tra, ts;
	int ista, nt, it, its;
	void shiftx( int, float *, float *, int );
	
	for( ista=0; ista<nsta; ista++ )
	{
		ts = ev[ista].time_shift_all;
		if(verbose)
		{
		  fprintf( stdout, "%s: time_shift(): ista=%d ts=%g\n", 
			progname, ista, ts );
		}
		if( ts == 0 ) continue;

	/*** transverse component ***/	
		nt = ev[ista].ew.s.npts;
		dt = ev[ista].ew.s.delta;
		its = ts/dt;
		tra = (float *)calloc( nt, sizeof(float) );
		shiftx( nt, ev[ista].ew.data, tra, its );

	/*** radial component ***/
		nt = ev[ista].ns.s.npts;
		dt = ev[ista].ns.s.delta;
		its = ts/dt;
		rad = (float *)calloc( nt, sizeof(float) );
		shiftx( nt, ev[ista].ns.data, rad, its );

	/*** vertical component ***/
		nt = ev[ista].z.s.npts;
		dt = ev[ista].z.s.delta;
		its = ts/dt;
		ver = (float *)calloc( nt, sizeof(float) );
		shiftx( nt, ev[ista].z.data, ver, its );

		if(verbose) 
		{
		  fprintf(stdout,"%s: time_shift(): %s.%s ista=%d nt=%d dt=%g its=%d\n",
			progname, 
			ev[ista].stnm, 
			ev[ista].net, 
			ista, 
			nt, 
			dt, 
			its );
		}

	/*** put back shift data into its original place in memory ***/
		for( it=0; it<nt; it++ )
		{
			ev[ista].ew.data[it] = tra[it];
			ev[ista].ns.data[it] = rad[it];
			ev[ista].z.data[it]  = ver[it];
		}
		free(rad);
		free(ver);
		free(tra);
	}
}

/******************************/
/*** load green's functions ***/
/******************************/

float *load_greens( EventInfo *ev, Greens **grn, int nsta, int *nz_tmp, int verbose )
{
        FILE *fpg;
        int ista, iz, nz;
        float *z;
        int it, nt;
        float min,max;

        for( ista=0; ista<nsta; ista++ )
        {
                if( (fpg=fopen( ev[ista].ginv_filename, "rb" )) == NULL )
                {
                        fprintf( stdout, "%s: Fatal Error cannot open file %s\n",
                                progname, ev[ista].ginv_filename );
                        exit(-1);
                }

                if(verbose)
		{
                  fprintf(stdout, 
			"%s: getting depth increment information \n", progname );
		}

                fread( &nz, sizeof(int), 1, fpg );
                if(verbose) 
		{
			fprintf( stdout, "%s: ista=%d nz=%d from %s\n", 
				progname, ista, nz, ev[ista].ginv_filename );
		}

                z = malloc( nz*sizeof(*z) );
                fread( &z[0], nz*sizeof(float), 1, fpg );

                if(verbose) 
		{
			fprintf( stdout, 
			  "%s: allocating memory for grn[nsta=%d][nz=%d]\n", 
				progname, ista, nz );
		}

                if(verbose) 
		{
			fprintf( stdout, "%s: opened file=%s nz=%d\n", 
				progname, ev[ista].ginv_filename, nz );
		}

                grn[ista] = (Greens *)malloc( nz*sizeof(Greens));
                for( iz=0; iz<nz; iz++ )
                {
                        if(verbose)
			{
				fprintf( stdout, 
				  "%s: reading grn[ista=%d][iz=%d]\n",
                                        progname, ista, iz );
			}
                        fread( &grn[ista][iz], sizeof(Greens), 1, fpg );

                        nt = grn[ista][iz].nt;
                        min = grn[ista][iz].g.rss[0];
                        max = grn[ista][iz].g.rss[0];

			for( it=1; it<nt; it++ )
                        {
                          if( grn[ista][iz].g.rss[it] < min ) min = grn[ista][iz].g.rss[0];
                          if( grn[ista][iz].g.rss[it] > max ) max = grn[ista][iz].g.rss[0];
                        }
                        if(verbose)
                        {
                          fprintf( stdout, 
			    "%s: iz=%d z=%g %s.%s rdist=%g az=%g t0=%g dt=%g nt=%d min=%g max=%g fn=%s modf=%s\n",
                                progname, iz, grn[ista][iz].evdp, grn[ista][iz].stnm, grn[ista][iz].net,
                                grn[ista][iz].rdist, grn[ista][iz].az, grn[ista][iz].t0, grn[ista][iz].dt,
                                grn[ista][iz].nt, min, max, grn[ista][iz].filename, grn[ista][iz].v.modfile );
                        }
                        grn[ista][iz].ts0 = ev[ista].ts0;

			grn[ista][iz].tra = (float *)calloc(grn[ista][iz].nt, sizeof(float) );
			grn[ista][iz].rad = (float *)calloc(grn[ista][iz].nt, sizeof(float) );
			grn[ista][iz].ver = (float *)calloc(grn[ista][iz].nt, sizeof(float) );

                }  /*** depth loop ****/

                fclose(fpg);

        } /*** station loop ****/

        *nz_tmp = nz;
        return (float *)z;
}

/************************************************/
/*** loop over the stations and load the data ***/
/************************************************/

void load_the_data( EventInfo *ev, int nsta, float ts0, int verbose )
{
	float max, min, mean;

        int ista;
        FILE *fpd;
        void shift_data( EventInfo *, int, float, int );
        void time_shift( EventInfo *, int, int );
	void sac_minmax( float *, int, float *, float *, float * );

        if(verbose) 
	{
		fprintf( stdout, 
		  "%s: load_the_data() nsta=%d\n", progname, nsta );
	}

        for( ista=0; ista<nsta; ista++ )
        {
                if(verbose)
		{
                  fprintf( stdout, 
			"%s: open file %s\n", progname, ev[ista].data_filename );
		}

                if( (fpd=fopen(ev[ista].data_filename, "rb" ))==NULL )
                {
                        fprintf(stderr, "%s: Fatal Error, cannot open file %s\n",
                                progname, ev[ista].data_filename );
                        exit(-1);
                }

                if(verbose)
		{
			fprintf( stdout,
			  "%s: reading data header...\n", progname );
		}

                fread( &ev[ista], sizeof(EventInfo), 1, fpd );

                if(verbose)
                {
                  fprintf( stdout, 
                    "%s: ista=%3d data=%s glib=%s ginv=%s np=%d pass=%d lf=%g hf=%g nt=%d dt=%g tr=%g tt=%g\n",
                        progname, ista, ev[ista].data_filename, ev[ista].glib_filename,
                        ev[ista].ginv_filename, ev[ista].npole, ev[ista].npass, ev[ista].lf,
                        ev[ista].hf, ev[ista].nt, ev[ista].dt, ev[ista].tr, ev[ista].tt );

                  fprintf( stdout, "%s: origin time=\t\t", progname );
                        WriteMyTime2STDOUT( &(ev[ista].ot ) );
                }

                if(verbose)
		{
			fprintf( stdout, "%s: Allocating memory for data...", progname );
		}

                ev[ista].ew.data = calloc( ev[ista].ew.s.npts, sizeof(float) );
                ev[ista].ns.data = calloc( ev[ista].ns.s.npts, sizeof(float) );
                ev[ista].z.data  = calloc( ev[ista].z.s.npts, sizeof(float) );

                if(verbose)
			fprintf( stdout, "%s: Reading the data...", progname );

                fread( &(ev[ista].ew.data[0]), ev[ista].ew.s.npts * sizeof(float), 1, fpd );
                fread( &(ev[ista].ns.data[0]), ev[ista].ns.s.npts * sizeof(float), 1, fpd );
                fread( &(ev[ista].z.data[0]),  ev[ista].z.s.npts  * sizeof(float), 1, fpd );
                fclose(fpd);

                if(verbose)
		{
			fprintf( stdout, "%s: DONE!\n", progname );
		}
		

/**********************************************************************************************/
/*** check to see if sacqc has flagged the data to user2=999 and kuser="SACQC"              ***/
/*** if any of the 3-channels is flagged as being out of QC statistical tolerance           ***/
/*** then the data are not iused=0 but still displayed or iused=-1 not displayed or used    ***/
/**********************************************************************************************/

		if(  ev[ista].z.s.user2 == 999 && strcmp( ev[ista].z.s.kuser2, "SACQC" ) == 0 )
			ev[ista].iused = -1;
		if( ev[ista].ns.s.user2 == 999 && strcmp( ev[ista].ns.s.kuser2, "SACQC" ) == 0 )
                        ev[ista].iused = -1;   
		if( ev[ista].ew.s.user2 == 999 && strcmp( ev[ista].ew.s.kuser2, "SACQC" ) == 0 )
                        ev[ista].iused = -1;   

                if(verbose)
                {
                  fprintf( stdout, 
			"%s: Read the file %s\n", progname, ev[ista].data_filename );

                  fprintf( stdout, 
			"%s: cmp=%-3.3s nt=%d dt=%g sta=%s net=%s min=%.3e max=%.3e\n",
                        progname, ev[ista].ew.s.kcmpnm, ev[ista].ew.s.npts,
                        ev[ista].ew.s.delta, ev[ista].ew.s.kstnm,
                        ev[ista].ew.s.knetwk, ev[ista].ew.s.depmin, ev[ista].ew.s.depmax );

                  fprintf( stdout, 
			"%s: cmp=%-3.3s nt=%d dt=%g sta=%s net=%s min=%.3e max=%.3e\n",
                        progname, ev[ista].ns.s.kcmpnm, ev[ista].ns.s.npts,
                        ev[ista].ns.s.delta, ev[ista].ns.s.kstnm,
                        ev[ista].ns.s.knetwk, ev[ista].ns.s.depmin, ev[ista].ns.s.depmax );

                  fprintf( stdout, 
			"%s: cmp=%-3.3s nt=%d dt=%g sta=%s net=%s min=%.3e max=%.3e\n",
                        progname, ev[ista].z.s.kcmpnm, ev[ista].z.s.npts,
                        ev[ista].z.s.delta, ev[ista].z.s.kstnm,
                        ev[ista].z.s.knetwk, ev[ista].z.s.depmin, ev[ista].z.s.depmax );
                }
        }

/*********************************************************/
/*** shift the data by time_shift_all value set in the ***/
/*** glib2inv.par file for just this station           ***/
/*********************************************************/
        if(verbose)
                fprintf( stdout, "%s: load_the_data(): calling time_shift()\n", progname );

        time_shift( ev, nsta, verbose );

/*********************************************************/
/*** shift the data by ts0 seconds                     ***/
/*********************************************************/
	if(verbose)
		fprintf( stdout, "%s: load_the_data(): calling shift_data()\n", progname );

        shift_data( ev, nsta, ts0, verbose );
}

/*********************************************************/
/*** shift the data by ts0 seconds                     ***/
/*** origin time shift                                 ***/
/*********************************************************/

void shift_data( EventInfo *ev, int nsta, float ts0, int verbose )
{
        float dt, *ver, *rad, *tra;
        int it, nt, ista, its0;
        void shiftx( int, float *, float *, int );

        if(verbose)
        {
          fprintf(stdout, "%s: shift_data(): nsta=%d ts0=%g \n",
                progname, nsta, ts0 );
        }

        for(ista=0; ista<nsta; ista++ )
        {

        /*** fill in some values ***/
                ev[ista].ts0 = ts0;

        /*** fix the time structure ***/
                epoch2time( &(ev[ista].ot), (ev[ista].ot.epoch - ts0 ) );
                ev[ista].ot_shift = ev[ista].ot.fsec;
                if( ts0 == 0 ) break;

                if(verbose)
                {
                        fprintf( stdout, "%s: shift_data(): shifting ts0=%g \t", 
				progname, ts0 );
                        WriteMyTime2STDOUT( &(ev[ista].ot) );
			fprintf( stdout, "\n" );
                }

        /*** transverse compoment ***/
                nt = ev[ista].ew.s.npts;
                dt = ev[ista].ew.s.delta;
                its0 = ts0/dt;
		/* fprintf( stdout, "shift tra ts0=%g its0=%d\n", ts0, its0 ); */
                tra = (float *)calloc( nt, sizeof(float) );
                shiftx( nt, ev[ista].ew.data, tra, its0 );

        /*** radial component ***/
                nt = ev[ista].ns.s.npts;
                dt = ev[ista].ns.s.delta;
                its0 = ts0/dt;
                rad = (float *)calloc( nt, sizeof(float) );
                shiftx( nt, ev[ista].ns.data, rad, its0 );

        /*** vertical component ***/
                nt = ev[ista].z.s.npts;
                dt = ev[ista].z.s.delta;
                its0 = ts0/dt;
                ver = (float *)calloc( nt, sizeof(float) );
                shiftx( nt, ev[ista].z.data, ver, its0 );

                if(verbose)
                {
                  fprintf( stdout, "%s: shift_data(): %s.%s. ista=%d nt=%d dt=%g its0=%d\n",
			progname, 
                        ev[ista].stnm, 
			ev[ista].net, 
			ista, 
			nt, 
			dt, 
			its0 );
                }

                for( it=0; it<nt; it++ )
                {
                        ev[ista].ew.data[it] = tra[it];
                        ev[ista].ns.data[it] = rad[it];
                        ev[ista].z.data[it]  = ver[it];
                }
                free(tra);
                free(rad);
                free(ver);
        }
}

/*******************************************/
/*** make output directory for GMT plots ***/
/*******************************************/

void mk_gmtplot_dir( char *output_pathname, int verbose )
{
        mode_t mode;
	extern char progname[128];

        if( mkdir( output_pathname, mode ) != 0 )
        {
                if( errno == EEXIST )
                {
                  if(verbose)
		  {
                    fprintf( stdout, "%s: mk_gmtplot_dir():  directory %s already exists\n", 
				progname, output_pathname );
		  }
                }
                else
                {
                  fprintf( stderr, "%s: mk_gmtplot_dir(): unknown mkdir() errno=%d %s\n", 
			progname, errno, output_pathname );
                  exit(-1);
                }
        }
        chmod( output_pathname, S_IRWXU|S_IRGRP|S_IXGRP|S_IROTH|S_IXOTH );
}

void dumpxy( EventInfo *ev, Solution *sol, Greens **grn, int nsta, int iz, int verbose )
{
	FILE *fp;
	char filename[256], synt_filename[256], data_filename[256];
	char grd_mo_type[5];
	float cm2mm = 10.;
	float cm2microns = 10000.;
	float cm2nanometers = 10000000;
	float scale;
	float *x, *syn_z, *syn_r, *syn_t, *dat_z, *dat_r, *dat_t;
	int ista, it, npts;
	float dt, beg;

/******************************/
/*** write the output files ***/
/******************************/

	scale = cm2nanometers;

	for( ista=0; ista<nsta; ista++ )
	{
		npts = ev[ista].ew.s.npts;
		dt   = ev[ista].ew.s.delta;
                beg  = ev[ista].ew.s.b;

                syn_z = (float *) malloc( npts * sizeof(float) );
                syn_r = (float *) malloc( npts * sizeof(float) );
                syn_t = (float *) malloc( npts * sizeof(float) );
                dat_z = (float *) malloc( npts * sizeof(float) );
                dat_r = (float *) malloc( npts * sizeof(float) );
                dat_t = (float *) malloc( npts * sizeof(float) );
                x     = (float *) malloc( npts * sizeof(float) );

		for( it=0; it<npts; it++ )
                {
                        dat_t[it] = ev[ista].ew.data[it]    * scale;
                        dat_r[it] = ev[ista].ns.data[it]    * scale;
                        dat_z[it] = ev[ista].z.data[it]     * scale;
                        syn_t[it] = ev[ista].syn_t.data[it] * scale;
                        syn_r[it] = ev[ista].syn_r.data[it] * scale;
                        syn_z[it] = ev[ista].syn_z.data[it] * scale;
                        x[it] = beg + (float)it * dt;
                }

                sprintf( grd_mo_type, "D" );
                if( ev[ista].grd_mo_type == VELOCITY ) sprintf( grd_mo_type, "V" );

	/*** DATA observations ***/

	/*** vertical component ***/

                sprintf( filename, "%s.%s.%s.%03d.z.dat.xy",
                        grn[ista][iz].stnm, grn[ista][iz].net, grd_mo_type, ista );
                fp = fopen( filename, "w" );
                fprintf( fp, "> %s %s %s %d VERTICAL DATA dist=%g az=%g tshift=%g\n",
                        grn[ista][iz].stnm, grn[ista][iz].net, grd_mo_type, ista,
                        grn[ista][iz].rdist, grn[ista][iz].az, ev[ista].time_shift_all );
                for( it=0; it<npts; it++ ) fprintf( fp, "%g %g\n", x[it], dat_z[it] );
                fclose(fp);

	/*** radial component ***/

                sprintf( filename, "%s.%s.%s.%03d.r.dat.xy",
                        grn[ista][iz].stnm, grn[ista][iz].net, grd_mo_type, ista );
                fp = fopen( filename, "w" );
                fprintf( fp, "> %s %s %s %d RADIAL DATA dist=%g az=%g tshift=%g\n",
                        grn[ista][iz].stnm, grn[ista][iz].net, grd_mo_type, ista,
                        grn[ista][iz].rdist, grn[ista][iz].az, ev[ista].time_shift_all );
                for( it=0; it<npts; it++ ) fprintf( fp, "%g %g\n", x[it], dat_r[it] );
                fclose(fp);

	/*** transverse component ***/

		sprintf( filename, "%s.%s.%s.%03d.t.dat.xy",
                        grn[ista][iz].stnm, grn[ista][iz].net, grd_mo_type, ista );
                fp = fopen( filename, "w" );
                fprintf( fp, "> %s %s %s %d TRANSVERSE DATA dist=%g az=%g tshift=%g\n",
                        grn[ista][iz].stnm, grn[ista][iz].net, grd_mo_type, ista,
                        grn[ista][iz].rdist, grn[ista][iz].az, ev[ista].time_shift_all );
                for( it=0; it<npts; it++ ) fprintf( fp, "%g %g\n", x[it], dat_t[it] );
                fclose(fp);

        /*** SYNTHETICS ***/

	/*** vertical component ***/

                sprintf( filename, "%s.%s.%s.%03d.z.syn.xy",
                        grn[ista][iz].stnm, grn[ista][iz].net, grd_mo_type, ista );
                fp = fopen( filename, "w" );
                fprintf( fp, "> %s %s %s %d VERTICAL SYN dist=%g az=%g tshift=%g\n",
                        grn[ista][iz].stnm, grn[ista][iz].net, grd_mo_type, ista,
                        grn[ista][iz].rdist, grn[ista][iz].az, ev[ista].time_shift_all );
                for( it=0; it<npts; it++ ) fprintf( fp, "%g %g\n", x[it], syn_z[it] );
                fclose(fp);
	
	/*** radial component ***/

                sprintf( filename, "%s.%s.%s.%03d.r.syn.xy",
                        grn[ista][iz].stnm, grn[ista][iz].net, grd_mo_type, ista );
                fp = fopen( filename, "w" );
                fprintf( fp, "> %s %s %s %d RADIAL SYN dist=%g az=%g tshift=%g\n",
                        grn[ista][iz].stnm, grn[ista][iz].net, grd_mo_type, ista,
                        grn[ista][iz].rdist, grn[ista][iz].az, ev[ista].time_shift_all );
                for( it=0; it<npts; it++ ) fprintf( fp, "%g %g\n", x[it], syn_r[it] );
                fclose(fp);

	/*** transverse component ***/

                sprintf( filename, "%s.%s.%s.%03d.t.syn.xy",
                        grn[ista][iz].stnm, grn[ista][iz].net, grd_mo_type, ista );
                fp = fopen( filename, "w" );
                fprintf( fp, "> %s %s %s %d TRANSVERSE SYN dist=%g az=%g tshift=%g\n",
                        grn[ista][iz].stnm, grn[ista][iz].net, grd_mo_type, ista,
                        grn[ista][iz].rdist, grn[ista][iz].az, ev[ista].time_shift_all );
                for( it=0; it<npts; it++ ) fprintf( fp, "%g %g\n", x[it], syn_t[it] );
                fclose(fp);

		free(x);
                free(syn_z);
                free(syn_r);
                free(syn_t);
                free(dat_z);
                free(dat_r);
                free(dat_t);

        } /*** loop and stations ***/
}
/*** finish writing out GMT input files ***/



void flt_sort_asc(int npts, float *fvec_a, float *fvec_b)
{
	/* Begin declaration(s) of function prototype(s) */
	void flt_swap(float *, float *);

	/* Begin declaration(s) of local variable(s)     */
	int i;
	int  j  = 1;

	/* Sort the a and b arrays; b is sorted with respect of a */
	if( npts < 1 )  {
		fprintf(stderr, "*** %s(): Numer of points less than one.  Exiting!\n", __func__);
		return;
	}

	for(i=0; i<npts; i++)  {
		j = i+1;
		while( j < npts )  {
			if( fvec_a[i] > fvec_a[j])  {

				flt_swap( &(fvec_a[i]), &(fvec_a[j]));

				if( fvec_b != NULL)  {
					flt_swap( &(fvec_b[i]), &(fvec_b[j]));
				}
			}
			j++;
		}
	}

	return;
}


void flt_swap(float *a, float *b)
{
	/* Begin declaration(s) of local variable(s) */
	float tempflt = 0.0f;

	/* Verify that the arguments are not NULL.  Otherwise, swap the values */
	if( (a == NULL) || (b == NULL) )  {
		fprintf(stderr, "*** %s(): WARNING:  One or both of the arguments is/are NULL.\n", __func__);
	}
	else  {
		tempflt = *a;
		*a      = *b;
		*b      = tempflt;
	}

	return;
}
