#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../include/nrutil.h"     /** numerical recipes **/
#include "../include/mt.h"         /** global datatype and structure declarations **/

struct event { 
	float str,dip,rak,Mw,Mo,tr,tt; 
	float Miso,Mdev,Mclvd,Mdc,Mtotal;
	float Piso,Pdev,Pclvd,Pdc;
	float stla,stlo,evla,evlo,evdp;
	float Mxx,Myy,Mzz,Mxy,Myz,Mxz;
	MyTime ot;
	MyTime beg;
};

char progname[128];

#define MT_INPUT	0
#define SDR_MO_INPUT	1
#define SDR_MISO_INPUT	2
#define ISO_CLVD_INPUT	3

int main( int ac, char **av )
{
        Greens grn;
	int nz,iz;
	int ista;
	FILE *fp;
	char filename[256];
	float *z, my_z, ztol;
	struct event ev;

	int verbose=0;
	int mtdegfree=5;
	int ifoundz=0;
	int input_type;
	int idumpgrn=0;
	int inoise;
	int iseed;
	float noise_Mw;

/******************************************************************************/
/*** timesubs.o ***/
/******************************************************************************/
	char dateid[256];
	MyTime ot;
	void parsestring( MyTime *, char * );
        void clone_mytime(  MyTime *, MyTime * );
        void WriteMyTime2STDOUT( MyTime * );

/******************************************************************************/
/*** function prototypes ***/
/******************************************************************************/
	void wrtgrn2sac( Greens *, int );
	void grn2disp2( Greens *, struct event, int, int, int, float, int, int );
	int setpar(int,char **), mstpar(), getpar();
	void endpar();
	void Usage_Print();

/**********************/
/**** begin program ***/
/**********************/
	strcpy( progname, av[0] );
	fprintf( stderr, "\n\n%s: Version=%s Release Date=%s\n",
                progname, Version_Label, Version_Date );

	if( ac <= 1 ) Usage_Print();

	setpar(ac,av);
	mstpar("glib",	  "s", filename );
	mstpar("z",	  "f", &my_z );
	getpar("verbose", "b", &verbose );
	getpar("dumpgrn", "b", &idumpgrn );

	strcpy( dateid, "2008/01/01,00:00:00\0" );
	getpar( "date", "s", dateid );

	getpar("noise",   "b", &inoise );
	if( inoise )
	{
		mstpar("nMw", "f", &noise_Mw );
		getpar("seed", "d", &iseed );
	}
	
	if( !idumpgrn )
	{
	  	if(verbose)
		  fprintf(stderr, "%s: idumpgrn is off\n", progname );

		input_type = SDR_MO_INPUT;  /*** default ***/
		mstpar("type", "d", &input_type );

		if( input_type == MT_INPUT )
		{
			mtdegfree = 6;
			mstpar("Mxx","f",&ev.Mxx );
			mstpar("Myy","f",&ev.Myy );
			mstpar("Mzz","f",&ev.Mzz );
			mstpar("Mxy","f",&ev.Mxy );
			mstpar("Mxz","f",&ev.Mxz );
			mstpar("Myz","f",&ev.Myz );
			mstpar("Mo", "f",&ev.Mo );
		}
		else if( input_type == SDR_MO_INPUT )
		{
			mtdegfree = 5;
			mstpar("str",	"f", &ev.str );
			mstpar("dip",	"f", &ev.dip );
			mstpar("rak",	"f", &ev.rak );
			mstpar("Mw",	"f", &ev.Mw );
			ev.Piso = 0;
			ev.Pdev = 1;
			ev.Pdc  = 1;
			ev.Pclvd = 0;
		}
		else if( input_type == SDR_MISO_INPUT )
		{
			mtdegfree = 5;
			mstpar("str",   "f", &ev.str );
			mstpar("dip",   "f", &ev.dip );
			mstpar("rak",   "f", &ev.rak );
			mstpar("Mw",    "f", &ev.Mw );

			mstpar("Piso",  "f", &ev.Piso );
			mstpar("Pclvd", "f", &ev.Pclvd );
			mstpar("Pdc",   "f", &ev.Pdc );
		}

	/*** check the ratio of isotropic to deviatoric moment ***/

		ev.Pdev = 1. - fabs(ev.Piso);

		if( ( fabs(ev.Pclvd) + fabs(ev.Pdc) ) == ev.Pdev )
		{
			fprintf(stdout, "%s: Pdev= %g Pclvd+Pdc=%g\n",
				progname, ev.Pdev, fabs(ev.Pclvd) + fabs(ev.Pdc)  );
		}
		else
		{
			fprintf(stderr,
				"%s: ERROR Percent Deviatoric does not equal 100%% - Isotropic: Pdev= %g Pclvd+Pdc=%g\n",
				progname, ev.Pdev, fabs(ev.Pclvd) + fabs(ev.Pdc)  );
			/* exit(-1); */
		}

	/*** check that all the percent moments adds to 100 percent ***/

		if( ( fabs(ev.Piso) + fabs(ev.Pclvd) + fabs(ev.Pdc) ) == 1.0 )
		{
			fprintf(stdout, "%s: Piso=%g Pclvd=%g Pdc=%g\n", progname,
				ev.Piso, ev.Pclvd, ev.Pdc );
		}
		else
		{
			fprintf(stderr,
				"%s: ERROR Percent Iso+Clvd+DC greater than 100%% \n",
				progname );
			/* exit(-1); */
		}

		ev.tr = 0; 
		ev.tt = 0;
		getpar("tr",    "f", &ev.tr );
		getpar("tt",    "f", &ev.tt );
	}
	else
	{
		if(verbose)
		 fprintf(stderr, "%s: dumpgrn is on\n", progname );
	}
	endpar();


/***********************************************************************************/
/*** END OF GETPAR                                                               ***/
/***********************************************************************************/
	
/*********************************************************/
/*** parse dateid string  and put into event structure ***/
/*********************************************************/
	if(verbose) 
		fprintf( stdout, "%s: dateid=%s\n", progname, dateid ); 

	parsestring( &ot, dateid );
	clone_mytime( &ot, &(ev.ot) );
	if(verbose) { fprintf( stdout, "%s: origin time=", progname ); WriteMyTime2STDOUT( &(ev.ot) ); }

	clone_mytime( &ot, &(ev.beg) );
	if(verbose) { fprintf( stdout, "%s: beginn time=", progname ); WriteMyTime2STDOUT( &(ev.beg) ); }

	if( verbose )
		fprintf(stderr, "%s: verbose is on\n", progname );

/**********************************************/
/*** open green's function file and read in ***/
/**********************************************/
	if( (fp = fopen( filename, "rb" ) ) == NULL )
	{
		fprintf(stderr, "%s: Fatal Error, cannot open file %s\n", 
			progname, filename );
		exit(-1);
	}

	fread(&nz,sizeof(int),1,fp);
	z = malloc(nz*sizeof(float));
	fread(&z[0],nz*sizeof(float),1,fp);

	if(verbose)
	{
	  fprintf( stdout, "%s: nz=%d ", progname, nz );
	  for(iz=0;iz<nz;iz++)printf(" z[%d]=%g ", iz, z[iz] );
	  printf("\n");
	}

	ifoundz = -1;
	ztol = 1.0E-05;
	for( iz=0; iz<nz; iz++ )
	{
		if( (my_z > z[iz]-ztol) && (my_z < z[iz]+ztol) )
		{
			ifoundz = iz;
			if(verbose)
			  fprintf( stdout, "%s: found iz=%d z=%g\n", 
			    progname,  ifoundz, z[ifoundz] );
			break;
		}
	}
	if( ifoundz < 0 )
	{
		fprintf( stderr, "%s: Fatal ERROR My_z=%g not found in z=",
			progname, my_z );
	}
	else
	{
		if(verbose)
		  fprintf(stderr, "%s: found my_z=%g ifoundz=%d iz=%d z=%g\n",
			progname, my_z, ifoundz, iz, z[ifoundz] );
	}

	for( iz=0; iz<nz; iz++ )
	{
		fread(&grn,sizeof(Greens),1,fp);
		if(verbose)
		{
		  fprintf( stdout, "%s: iz=%d z=%g rdist=%g az=%g evdp=%g t0=%g dt=%g nt=%d %s\n",
			progname, 
			iz, z[iz], grn.rdist, grn.az, grn.evdp, grn.t0, 
			grn.dt, grn.nt, grn.filename );
		}

		if( iz == ifoundz )
		{
		  ista = 0;
		  if( idumpgrn == 1 ) 
		  {
			if(verbose)
			  fprintf( stdout, "%s: calling wrtgrn2sac: \n", progname );
			wrtgrn2sac( &grn, ista );
		  }
		  if( idumpgrn == 0 ) 
		  {
		    if(verbose)
			fprintf( stdout, "%s: calling grn2disp2: \n", progname );
		    grn2disp2( &grn, ev, verbose, mtdegfree, inoise, noise_Mw, iseed, input_type );
		  }
	          break;
		}
	}
	fclose(fp);
	free(z);

} /*** end of main.c ***/

/**************************************************************************************************/
/*** GRN2DISP2.c         ***/
/**************************************************************************************************/

void grn2disp2( Greens *g, struct event ev, int verbose, int mtdegfree, 
  int inoise, float noise_Mw, int iseed, int input_type )
{
	float *tra, *rad, *ver, *ns, *ew;
	float *rss, *rds, *rdd, *rep, *zss, *zds, *zdd, *zep, *tss, *tds;

	float a1,a2,a3,a4,a5;
	Tensor M;
	float strr,dipr,rakr,Miso,tmpaz,angle,cmpinc,az;
	Sac_Header sp;
	char sacfile[256];
	FILE *fp;
	int it;
	float pi, d2r, dt, t0, e, tt, tr, fi, twin, r, dist, azimuth;
	float slip4pi_iso, slip4pi_dev, slip4pi_tot;
	float third, half, sixth;
	int nt;
	float cm2m = 0.01;
/* see include/mt.h Mo = math.pow( 10.0, 1.5*(Mw+10.73) ) = 1.2445146117713818e+16; Reference Mw = 0.0 */

	int new_nt;
	float new_dt;
	float *xtra, *xrad, *xver, *xns, *xew;

	float fac, noise_moment, peak;
	float gasdev( int * );
	void source_time_function( float *, int, float, float, float );
	void  rotate( float *, int, float *, float, float *, int, float *, float, float, int );
	void interpolate_wiggins( float *, int, float, float, float *, int, float );

	void write_SACPZ_file( char *, char *, char *,  int, int, int, int );

	half = 0.5;
	third = 0.33333333;
	sixth = 0.166666667;
	pi = acos(-1.0);
	d2r = pi/180;
	nt = g->nt;
	dt = g->dt;
	t0 = g->t0;
	twin = (float)nt*dt;
	e = t0 + (nt*dt);
	tt = ev.tt;
	tr = ev.tr;

/**************************************************************************************************/
/*** allocate memory for displacement vector ***/
/**************************************************************************************************/

	tra = malloc( nt*sizeof(float) );
	rad = malloc( nt*sizeof(float) );
	ver = malloc( nt*sizeof(float) );

/**************************************************************************************************/
/*** set t0 from reduction velocity redv ***/
/**************************************************************************************************/

	if( g->redv > 0 ) g->t0 = g->rdist/g->redv;
	t0 = g->t0;

/**************************************************************************************************/
/*** simplify by assigning local pointers to structure ***/
/**************************************************************************************************/

	rss = g->g.rss;
	rds = g->g.rds;
        rdd = g->g.rdd;
        rep = g->g.rep;
        zss = g->g.zss;
        zds = g->g.zds;
        zdd = g->g.zdd;
        zep = g->g.zep;
        tss = g->g.tss;
        tds = g->g.tds;

/**************************************************************************************************/
/*** find the peak value of the GFs for noise scaling ***/
/**************************************************************************************************/

	peak = 1E-29;
	for( it = 0; it < nt; it++ )
	{
          if(fabs(rss[it]) > peak) peak=fabs(rss[it]);
          if(fabs(rds[it]) > peak) peak=fabs(rds[it]);
          if(fabs(rdd[it]) > peak) peak=fabs(rdd[it]);
          if(fabs(rep[it]) > peak) peak=fabs(rep[it]);
          if(fabs(zss[it]) > peak) peak=fabs(zss[it]);
          if(fabs(zds[it]) > peak) peak=fabs(zds[it]);
          if(fabs(zdd[it]) > peak) peak=fabs(zdd[it]);
          if(fabs(zep[it]) > peak) peak=fabs(zep[it]);
          if(fabs(tss[it]) > peak) peak=fabs(tss[it]);
          if(fabs(tds[it]) > peak) peak=fabs(tds[it]);
	}

/**************************************************************************************************/
/*** convolve source time funciton ***/
/**************************************************************************************************/

	source_time_function( rss, nt, dt, tr, tt );
	source_time_function( rds, nt, dt, tr, tt );
	source_time_function( rdd, nt, dt, tr, tt );
	source_time_function( rep, nt, dt, tr, tt );
	source_time_function( zss, nt, dt, tr, tt );
	source_time_function( zds, nt, dt, tr, tt );
	source_time_function( zdd, nt, dt, tr, tt );
	source_time_function( zep, nt, dt, tr, tt );
	source_time_function( tss, nt, dt, tr, tt );
	source_time_function( tds, nt, dt, tr, tt );

/**************************************************************************************************/
/*** compute the 3 component displacements ***/
/**************************************************************************************************/

	fi = g->az * d2r;

	if( input_type == SDR_MO_INPUT || input_type == SDR_MISO_INPUT )
	{
		if(verbose)
			fprintf( stdout, "%s: SDR_MO_INPUT or SDR_MISO_INPUT format\n", progname );

		ev.Mo = pow( 10., (1.5*(ev.Mw + 10.73)));

		ev.Miso = ev.Piso * ev.Mo/base_moment;

		ev.Mdc  = ev.Pdc * ev.Mo/base_moment;

		ev.Mclvd = ev.Pclvd * ev.Mo/base_moment;

		ev.Mtotal = ev.Miso + ev.Mclvd + ev.Mdc;

		if(verbose) 
		{
		  fprintf( stdout, "%s: Mw=%g Mo=%g Miso=%g(%0.f%%) Mclvd=%g(%0.f%%) Mdc=%g(%0.f%%) Mtotal=%g\n",
			progname, 
			ev.Mw, 
			ev.Mo, 
			ev.Miso   * base_moment, ev.Piso  * 100, 
			ev.Mclvd  * base_moment, ev.Pclvd * 100, 
			ev.Mdc    * base_moment, ev.Pdc   * 100,
			ev.Mtotal * base_moment );
		}

		strr = ev.str * d2r;
		dipr = ev.dip * d2r;
		rakr = ev.rak * d2r;

		M.xx = -( sin(dipr) * cos(rakr) * sin(2*strr) + sin(2*dipr) * sin(rakr) * sin(strr)*sin(strr) );
		M.yy =  ( sin(dipr) * cos(rakr) * sin(2*strr) - sin(2*dipr) * sin(rakr) * cos(strr)*cos(strr) );
		M.zz =  ( sin(2*dipr) * sin( rakr ) );
		M.xy =  ( sin(dipr) * cos(rakr) * cos(2*strr) + 0.5*sin(2*dipr) * sin(rakr) * sin(2*strr) );
		M.xz = -( cos(dipr) * cos(rakr) * cos(strr) + cos(2*dipr) * sin(rakr) * sin(strr) );
		M.yz = -( cos(dipr) * cos(rakr) * sin(strr) - cos(2*dipr) * sin(rakr) * cos(strr) );

	/***
		if( fabs(M.xx) < 1.0E-6 ) M.xx=0;
		if( fabs(M.yy) < 1.0E-6 ) M.yy=0;
		if( fabs(M.zz) < 1.0E-6 ) M.zz=0;
		if( fabs(M.xy) < 1.0E-6 ) M.xy=0;
		if( fabs(M.xz) < 1.0E-6 ) M.xz=0;
		if( fabs(M.yz) < 1.0E-6 ) M.yz=0;
	****/

		if(verbose)
		{
		  fprintf( stdout, "%s: str=%g dip=%g rak=%g Mw=%g Mo=%g\n",
		    progname, ev.str, ev.dip, ev.rak, ev.Mw, ev.Mo );
		  fprintf( stdout, "%s: fi=%g strr=%g dipr=%g rakr=%g Mxx=%g Myy=%g Mzz=%g Mxy=%g Mxz=%g Myz=%g\n",
		    progname, fi, strr, dipr, rakr, M.xx, M.yy, M.zz, M.xy, M.xz, M.yz);
		}

	/*********************************/
	/*** compute the coefficients ****/
	/*********************************/
		a1 = half * ( M.xx - M.yy ) * cos( 2 * fi ) + M.xy * sin( 2 * fi );
		a2 = M.xz * cos( fi ) + M.yz * sin( fi );

		a3 = -half*( M.xx + M.yy );
		/* a3 = -sixth*( M.xx + M.yy - 2 * M.zz ); */

		a4 = half * ( M.xx - M.yy ) * sin( 2 * fi ) - M.xy * cos( 2 * fi );
		a5 = -M.yz * cos( fi ) + M.xz * sin( fi );
 
		if( verbose )
	 	  fprintf( stdout, "%s: a1=%g a2=%g a3=%g a4=%g a5=%g\n", 
			progname, a1, a2, a3, a4, a5 );

		for( it=0; it<nt; it++)
		{
		  ver[it] = ( a1 * zss[it] + a2 * zds[it] + a3 * zdd[it] ) * ev.Mdc - zdd[it] * ev.Mclvd + zep[it] * ev.Miso;
		  rad[it] = ( a1 * rss[it] + a2 * rds[it] + a3 * rdd[it] ) * ev.Mdc - rdd[it] * ev.Mclvd + rep[it] * ev.Miso;
		  tra[it] = ( a4 * tss[it] + a5 * tds[it] ) * ev.Mdc;           
		}
	}

	if( input_type == MT_INPUT )
	{
		ev.Mtotal = ev.Mo/base_moment;
		if( verbose ) 
		  fprintf( stdout, "%s: Mtotal=%e Mo=%e\n", progname, ev.Mtotal, ev.Mo );

		M.xx = ev.Mxx;
		M.yy = ev.Myy;
		M.zz = ev.Mzz;
		M.xy = ev.Mxy;
		M.xz = ev.Mxz;
		M.yz = ev.Myz;

		if( verbose )
		  fprintf( stdout, "%s: MT_INPUT format Mxx=%g Mxx=%g Mzz=%g Mxy=%g Mxz=%g Myz=%g\n",
			progname, 
			M.xx, M.yy, M.zz, M.xy, M.xz, M.yz );

		for( it = 0; it < nt; it++ )
		{
		  tra[it]=( M.xx * (  half * sin(2*fi) * tss[it] )
			+ M.yy * ( -half * sin(2*fi) * tss[it] )
			+ M.xy * (       -cos(2*fi) * tss[it] )
			+ M.xz * (        sin(fi)   * tds[it] )
			+ M.yz * (       -cos(fi)   * tds[it] ) ) * ev.Mtotal;

		  rad[it]=( M.xx * ( -sixth*rdd[it] + half*cos(2*fi)*rss[it] + third*rep[it] )
			+ M.yy * ( -sixth*rdd[it] - half*cos(2*fi)*rss[it] + third*rep[it] )
			+ M.xy * (  sin(2*fi) * rss[it] )
			+ M.xz * (  cos(fi)   * rds[it] )
			+ M.yz * (  sin(fi)   * rds[it] )
			+ M.zz * (  third*rdd[it] + third*rep[it] ) ) * ev.Mtotal;

		  ver[it]=( M.xx * ( -sixth*zdd[it] + half*cos(2*fi)*zss[it] + third*zep[it] )
			+ M.yy * ( -sixth*zdd[it] - half*cos(2*fi)*zss[it] + third*zep[it] )
			+ M.xy * (  sin(2*fi) * zss[it] )
			+ M.xz * (  cos(fi)   * zds[it] )
			+ M.yz * (  sin(fi)   * zds[it] )
			+ M.zz * (  third*zdd[it] + third*zep[it] ) ) * ev.Mtotal;
		}
	}

	if( inoise )
	{
		noise_moment = pow( 10, (1.5*(noise_Mw + 10.73)));
		fac = 0.1 * peak * ( noise_moment/base_moment );
		if(verbose) 
		  fprintf( stdout, "%s: fac=%e peak=%e noise_moment=%e noise_Mw=%g iseed=%d\n",
			progname, fac, peak, noise_moment, noise_Mw, iseed );
		for( it=0; it<nt; it++ )
		{
			tra[it] += fac*gasdev(&iseed);
			rad[it] += fac*gasdev(&iseed);
			ver[it] += fac*gasdev(&iseed);
		}
	}

/**************************************************************************************************/
/*** convert displacement from meters to cm ***/
/**************************************************************************************************/

	for( it=0; it<nt; it++ )
	{
		tra[it] *= cm2m;
		rad[it] *= cm2m;
		ver[it] *= cm2m;
	}
	if(verbose) 
		fprintf( stdout, "%s: converting amplitudes cm to meters\n", progname );

/**************************************************************************************************/
/*** rotate ***/
/**************************************************************************************************/

	ns = calloc( nt, sizeof(float) );
	ew = calloc( nt, sizeof(float) );
	for( it=0; it<nt; it++ )
	{
		ew[it] = tra[it];
		ns[it] = rad[it];
	}
	angle = 360 - g->az;
	az = g->az;
	tmpaz = az + 90;
	if( tmpaz > 360 ) tmpaz -= 360;
	cmpinc = 90;
	rotate( ns, nt, &az, cmpinc, ew, nt, &tmpaz, cmpinc, angle, verbose );
	if(verbose) 
	  fprintf( stdout, "%s: rotated to RADIAL(%g)->NS TRANSVERSE(%g)->EW by angle = %g\n", progname, g->az, tmpaz, angle );

/**************************************************************************************************/
/*** interpolate to dt = 0.05 sec/sample using wiggins interpolation scheme ***/
/**************************************************************************************************/

	new_dt = 0.05;
	new_nt = twin/new_dt;

	xtra = calloc( new_nt, sizeof(float) );
	xrad = calloc( new_nt, sizeof(float) );
	xver = calloc( new_nt, sizeof(float) );
	xew  = calloc( new_nt, sizeof(float) );
	xns  = calloc( new_nt, sizeof(float) );

	interpolate_wiggins( tra, nt, dt, sp.b, xtra, new_nt, new_dt );
	interpolate_wiggins( rad, nt, dt, sp.b, xrad, new_nt, new_dt );
	interpolate_wiggins( ver, nt, dt, sp.b, xver, new_nt, new_dt );
	interpolate_wiggins( ns,  nt, dt, sp.b, xns,  new_nt, new_dt );
        interpolate_wiggins( ew,  nt, dt, sp.b, xew,  new_nt, new_dt );

	nt = new_nt;
	dt = new_dt;
	if(verbose) 
		fprintf( stdout, "%s: interpolated to nt=%d dt=%g\n", 
			progname, nt, dt );

/*****************************************/
/**** SAC HEADER                     *****/
/*****************************************/
	sp = sac_null;

	sp.nvhdr = 6;
	sp.norid = 0;
	sp.nevid = 0;
	sp.iftype = ITIME;
	sp.idep = IUNKN;
	sp.iztype = IO;
	sp.ievtyp = IQUAKE;
	sp.leven = TRUE;
	sp.lpspol = FALSE;
	sp.lcalda = TRUE;

	sp.npts = nt;
	sp.delta = dt;
	sp.b = t0;
	sp.e = e;

/**************************************************************************************************/
/*** set origin time ***/
/**************************************************************************************************/

	sp.nzyear = ev.ot.year;
	sp.nzjday = ev.ot.jday;
	sp.nzhour = ev.ot.hour;
	sp.nzmin  = ev.ot.min;
	sp.nzsec  = ev.ot.isec;
	sp.nzmsec = ev.ot.msec;
	sp.o      = 0;
	strcpy( sp.ko, "OT" );

/**************************************************************************************************/
/*** set origin hypocenter and receiver locations ***/
/**************************************************************************************************/

	sp.evla = g->evla;
	sp.evlo = g->evlo;
	sp.evdp = g->evdp;
	sp.stla = g->stla;
	sp.stlo = g->stlo;
	sp.stel = g->stel;

	if( g->redv > 0 )
	{
	 sprintf( sp.kt1, "%gkm/s", g->redv ); 
	 sp.t1 = g->rdist/g->redv;  /*** Reduction Velocity km/km/sec ***/
	}
	
	sp.user0 = ev.str; strcpy( sp.kuser0, "Strike" );
	sp.user1 = ev.dip; strcpy( sp.kuser0, "Dip" );
	sp.user2 = ev.rak; strcpy( sp.kuser0, "Rake" );
	sp.user3 = ev.Mo;
	sp.user4 = ev.Mw;
	sp.user5 = ev.tr;
	sp.user6 = ev.tt;
	sp.mag  = ev.Mw;

	strcpy( sp.kstnm, g->stnm );
	strcpy( sp.knetwk, g->net );

/**************************************************************************************************/
/*** EAST-WEST component ***/
/**************************************************************************************************/

	strcpy( sp.kcmpnm, "BHE" );
	sp.cmpinc = 90;
	sp.cmpaz  = 90;

	write_SACPZ_file( sp.kstnm, sp.knetwk, sp.kcmpnm, sp.nzyear, sp.nzjday, sp.nzhour, sp.nzmin );

	sprintf( sacfile, "%4d.%03d.%02d.%02d.%02d.0000..%s.%s.%s.SAC",
		sp.nzyear, sp.nzjday, sp.nzhour, sp.nzmin, sp.nzsec, sp.knetwk, sp.kstnm, sp.kcmpnm );

	printf( "%s: writting file %s\n", progname, sacfile );

	if( (fp=fopen(sacfile,"w")) == NULL ) printf("ERROR cannot write to file\n" );
	fwrite( &sp, sizeof(Sac_Header), 1, fp );
	fwrite( &xew[0], sp.npts*sizeof(float), 1, fp );
	fclose(fp);

/**************************************************************************************************/
/*** NORTH-SOUTH Component ***/
/**************************************************************************************************/

	strcpy( sp.kcmpnm, "BHN" );
	sp.cmpinc = 90;
	sp.cmpaz  = 0;

	write_SACPZ_file( sp.kstnm, sp.knetwk, sp.kcmpnm, sp.nzyear, sp.nzjday, sp.nzhour, sp.nzmin );

	sprintf( sacfile, "%4d.%03d.%02d.%02d.%02d.0000..%s.%s.%s.SAC",
                sp.nzyear, sp.nzjday, sp.nzhour, sp.nzmin, sp.nzsec, sp.knetwk, sp.kstnm, sp.kcmpnm );
	printf( "%s: writting file %s\n", progname, sacfile );

	if( (fp=fopen(sacfile,"w")) == NULL ) printf("ERROR cannot write to file\n" );
	fwrite( &sp, sizeof(Sac_Header), 1, fp );
	fwrite( &xns[0], sp.npts*sizeof(float), 1, fp );
        fclose(fp);

/**************************************************************************************************/
/*** radial ***/
/**************************************************************************************************/
/****
	sprintf( sp.kcmpnm, "R%3.0f", g->az );
	sp.cmpinc = 90;
	sp.cmpaz = g->az;
	sprintf( sacfile, "%s.rad", g->filename );
	printf( "%s: writting file %s\n", progname, sacfile );
	if( (fp=fopen(sacfile,"w")) == NULL ) printf("ERROR cannot write to file\n" );
	fwrite( &sp, sizeof(Sac_Header), 1, fp );
	fwrite( &xrad[0], sp.npts*sizeof(float), 1, fp );
	fclose(fp);
****/

/**************************************************************************************************/
/*** vertical ***/
/**************************************************************************************************/

	sprintf( sp.kcmpnm, "BHZ" );
	sp.cmpinc = 0;
	sp.cmpaz = 0;

	write_SACPZ_file( sp.kstnm, sp.knetwk, sp.kcmpnm, sp.nzyear, sp.nzjday, sp.nzhour, sp.nzmin );

	sprintf( sacfile, "%4d.%03d.%02d.%02d.%02d.0000..%s.%s.%s.SAC",
                sp.nzyear, sp.nzjday, sp.nzhour, sp.nzmin, sp.nzsec, sp.knetwk, sp.kstnm, sp.kcmpnm );

	printf( "%s: writting file %s\n", progname, sacfile );

	if( (fp=fopen(sacfile,"w")) == NULL ) printf("ERROR cannot write to file\n" );
	fwrite( &sp, sizeof(Sac_Header), 1, fp );
	fwrite( &xver[0], sp.npts*sizeof(float), 1, fp );
	fclose(fp);

/**************************************************************************************************/
/*** transverse ***/
/**************************************************************************************************/
/****
	tmpaz = g->az+90;
	if( tmpaz>360 ) tmpaz-=360;
	sp.cmpaz = tmpaz;
	sp.cmpinc = 90;
	sprintf( sp.kcmpnm, "T%3.0f", tmpaz );
	sprintf( sacfile, "%s.tra", g->filename );
        printf( "%s: writting file %s\n", progname, sacfile );
        if( (fp=fopen(sacfile,"w")) == NULL ) printf("ERROR cannot write to file\n" );
        fwrite( &sp, sizeof(Sac_Header), 1, fp );
        fwrite( &xtra[0], sp.npts*sizeof(float), 1, fp );
        fclose(fp);
***/

	free( tra );
	free( rad );
	free( ver );
	free( xtra );
	free( xrad );
	free( xver );
	free( ns );
	free( ew );
	free( xns );
	free( xew );

} /*** end of grn2disp2.c ***/

void Usage_Print()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "USAGE: %s: glib= z= [no]verbose [no]dumpgrn tr=[0] tt=[0]", 
		progname );
	fprintf(stderr, "  type=[no default see below for additional options]\n\n" );

	fprintf(stderr, "\t REQUIRED ARGUMENTS:\n");
	fprintf(stderr, "\t\t glib= Green's function library (GLIB) name computed from mkgrnlib.c\n");
	fprintf(stderr, "\t\t z= the depth of the source to look up in the GLIB file\n");
	fprintf(stderr, "\n");

	fprintf(stderr, "\t OPTIONAL PARAMETERS:\n");
	fprintf(stderr, "\t\t [no]dumpgrn write out only Green's functions\n");
	fprintf(stderr, "\t\t [no]verbose verbosy output for diagnosis\n");
	fprintf(stderr, "\t\t date=YYYY/MM/DD,HH24:mm:ss.ss  origin time Default = [2008/01/01:00:00:00.00]\n");
	fprintf(stderr, "\n");

	fprintf(stderr, "\t\t tr=0.0      Rise time in seconds (triangle or trapazoid src time func)\n");
	fprintf(stderr, "\t\t tt=0.0      source duration in sec.  =0 if triangle, or >0 if trapazoid\n");
	fprintf(stderr, "\n");

	fprintf(stderr, "\t\t [no]noise   add white Gaussian noise\n");
	fprintf(stderr, "\t\t nMw=3.8     level of the noise in units of Mw for freq band of interest\n");
	fprintf(stderr, "\t\t seed=1      random seed\n");
	fprintf(stderr, "\n");

	fprintf(stderr, "\t\t type=       Source Input Mode Type\n");
	fprintf(stderr, "\t\t\t =0   Input Moment Tensor\n");
	fprintf(stderr, "\t\t\t =1   Input Pure Deviatoric Source (Strike/Dip/Rake)\n");
	fprintf(stderr, "\t\t\t =2   Input Pure Deviatoric source plus isotropic component\n");
	fprintf(stderr, "\n");

	fprintf(stderr, "\t\t if type=0 \n");
	fprintf(stderr, "\t\t\t Mxx=Normalized Element\n");
	fprintf(stderr, "\t\t\t Myy=Normalized Element\n");
	fprintf(stderr, "\t\t\t Mzz=Normalized Element\n");
	fprintf(stderr, "\t\t\t Mxy=Normalized Element\n");
	fprintf(stderr, "\t\t\t Mxz=Normalized Element\n");
	fprintf(stderr, "\t\t\t Myz=Normalized Element\n");
	fprintf(stderr, "\t\t\t Mo=Total Moment\n");
	fprintf(stderr, "\t\t if type=1 \n");
	fprintf(stderr, "\t\t\t str= fault strike\n");
	fprintf(stderr, "\t\t\t dip= fault dip\n");
	fprintf(stderr, "\t\t\t rak= fault rake\n");
	fprintf(stderr, "\t\t\t Mw= Total Moment Magnitude\n");
	fprintf(stderr, "\t\t if type=2 \n");
	fprintf(stderr, "\t\t\t str= fault strike\n");
	fprintf(stderr, "\t\t\t dip= fault dip\n");
	fprintf(stderr, "\t\t\t rak= fault rake\n");
	fprintf(stderr, "\t\t\t Mw= Total Moment Magnitude\n");
	fprintf(stderr, "\t\t\t Piso= The percent of the total moment allocated to isotropic\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Below mstpar will guide you:\n");

} /*** end of Print_Usage.c  ***/

void write_SACPZ_file( char *sta, char *net, char *cmp, int year, int jday, int hour, int min )
{
	FILE *fp;
	char sacpzfilename[256];
	sprintf( sacpzfilename, "SAC_PZs_%s_%s_%s__%4d.%03d.%02d.%02d.00.0000",
		net, sta, cmp, year, jday, hour, min );
	fp = fopen( sacpzfilename, "w" );
	fprintf( fp, "ZEROS 0\n" );
	fprintf( fp, "POLES 0\n" );
	fprintf( fp, "CONSTANT 1.00E+00\n" );
	fclose(fp);
}
