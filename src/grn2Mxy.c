#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../include/nrutil.h"     /** numerical recipes **/
#include "../include/mt.h"         /** global datatype and structure declarations **/

char progname[128];

#define MT_INPUT        0
#define SDR_MO_INPUT    1
#define SDR_MISO_INPUT  2
#define ISO_CLVD_INPUT  3

int main( int ac, char **av )
{
	Greens grn;
	int nz, iz;
	int ista;
	FILE *fp;
	char filename[256];
	float *z, my_z, ztol;
	MyTime ot;
	MyTime beg;
	
	int verbose = 0;
	int mtdegfree = 5;
	int ifoundz = 0;
	int input_type; 

/******************************************************************************/
/*** timesubs.o ***/
/******************************************************************************/
	char dateid[256];
	void parsestring( MyTime *t, char *str );
	void clone_mytime(  MyTime *t2, MyTime *t1 );
	void WriteMyTime2STDOUT( MyTime * );

	int setpar(int ac, char **av), mstpar(), getpar();
	void endpar();

	void conv2mxy( Greens *g );

/**********************/
/**** begin program ***/
/**********************/
	strcpy( progname, av[0] );
        fprintf( stderr, "\n\n%s: Version=%s Release Date=%s\n",
                progname, Version_Label, Version_Date );

	setpar(ac,av);
	mstpar("glib",    "s", filename );
	mstpar("z",       "f", &my_z );
	getpar("verbose", "b", &verbose );
	strcpy( dateid, "2008/01/01,00:00:00\0" );
	getpar( "date", "s", dateid );
	endpar();

/*** set the time ***/
/*
	parsestring( &ot, dateid );
	clone_mytime( &ot, &(ot) );
	clone_mytime( &ot, &(beg) );
*/

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

/**************************************************************************/
/*** search Green function library for the requested depth              ***/
/**************************************************************************/
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

		if( iz == ifoundz )
		{
			ista = 0;
			/*** do something ***/
			conv2mxy( &grn );
		}
		break;
	}
	fclose(fp);
	free(z);

} /*** end of main.c ***/

void conv2mxy( Greens *grn )
{
	Sac_Header sp;
	char sacfile[256];
	FILE *fp;
	int it;

	float *rss, *rds, *rdd, *rep, *zss, *zds, *zdd, *zep, *tss, *tds;
	float *txx, *txy, *txz, *tyy, *tyz;
        float *rxx, *rxy, *rxz, *ryy, *ryz, *rzz;
        float *zxx, *zxy, *zxz, *zyy, *zyz, *zzz;

	float sc;
	float third = 0.333333333333333333;
	float sixth = 0.166666666666666667;
	float half  = 0.500000000000000000;
	int nt;
	float dt, t0, twin, fi, e;

/*** begin subroutine ***/

	fi = grn->az * (M_PI/180);
	nt = grn->nt;
	dt = grn->dt;
	t0 = grn->t0;
	twin = (float)nt * dt;
	e = t0 + (nt*dt);
	
/**************************************************************************************************/
/*** allocate memory for green func ***/
/**************************************************************************************************/
	txx = calloc( nt, sizeof(float) );
	txy = calloc( nt, sizeof(float) );
	txz = calloc( nt, sizeof(float) );
	tyy = calloc( nt, sizeof(float) );
	tyz = calloc( nt, sizeof(float) );
	
	rxx = calloc( nt, sizeof(float) );
	rxy = calloc( nt, sizeof(float) );
	rxz = calloc( nt, sizeof(float) );
	ryy = calloc( nt, sizeof(float) );
	ryz = calloc( nt, sizeof(float) );
	rzz = calloc( nt, sizeof(float) );

	zxx = calloc( nt, sizeof(float) );
	zxy = calloc( nt, sizeof(float) );
	zxz = calloc( nt, sizeof(float) );
	zyy = calloc( nt, sizeof(float) );
	zyz = calloc( nt, sizeof(float) );
	zzz = calloc( nt, sizeof(float) );

/**************************************************************************************************/
/*** simplify by assigning local pointers to structure ***/
/**************************************************************************************************/
        rss = grn->g.rss;
        rds = grn->g.rds;
        rdd = grn->g.rdd;
        rep = grn->g.rep;
        zss = grn->g.zss;
        zds = grn->g.zds;
        zdd = grn->g.zdd;
        zep = grn->g.zep;
        tss = grn->g.tss;
        tds = grn->g.tds;

	sc = 1.0E+20 / base_moment;

	for( it = 0; it < nt; it++ )
	{
	  txx[it] = sc*( half * sin(2*fi) * tss[it] );
	  txy[it] = sc*(       -cos(2*fi) * tss[it] );
	  txz[it] = sc*(        sin(fi)   * tds[it] );
	  tyy[it] = sc*(-half * sin(2*fi) * tss[it] );
	  tyz[it] = sc*(       -cos(fi)   * tds[it] );
  
	  rxx[it] = sc*( +half*cos(2*fi) * rss[it] - sixth*rdd[it] + third*rep[it] );
	  rxy[it] = sc*( sin(2*fi) * rss[it] );
	  rxz[it] = sc*(   cos(fi) * rds[it] );
	  ryy[it] = sc*( -half*cos(2*fi) * rss[it] - sixth*rdd[it] + third*rep[it] );
	  ryz[it] = sc*(   sin(fi) * rds[it] );
	  rzz[it] = sc*( third * rdd[it] + third * rep[it] );

	  zxx[it] = sc*( +half*cos(2*fi) * zss[it] - sixth*zdd[it] + third*zep[it] );
	  zxy[it] = sc*(   sin(2*fi) * zss[it] );
	  zxz[it] = sc*(   cos(fi)   * zds[it] );
	  zyy[it] = sc*( -half*cos(2*fi) * zss[it] - sixth*zdd[it] + third*zep[it] );
	  zyz[it] = sc*(   sin(fi)   * zds[it] );
	  zzz[it] = sc*( third * zdd[it] + third * zep[it] );
	}

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
/*
        sp.nzyear = ot.year;
        sp.nzjday = ot.jday;
        sp.nzhour = ot.hour;
        sp.nzmin  = ot.min;
        sp.nzsec  = ot.isec;
        sp.nzmsec = ot.msec;
*/
        sp.o      = 0;
        strcpy( sp.ko, "OT" );
/**************************************************************************************************/
/*** set origin hypocenter and receiver locations ***/
/**************************************************************************************************/

        sp.evla = grn->evla;
        sp.evlo = grn->evlo;
        sp.evdp = grn->evdp;
        sp.stla = grn->stla;
        sp.stlo = grn->stlo;
        sp.stel = grn->stel;

        if( grn->redv > 0 )
        {
         sprintf( sp.kt1, "%gkm/s", grn->redv );
         sp.t1 = grn->rdist/grn->redv;  /*** Reduction Velocity km/km/sec ***/
        }

	strcpy( sp.kcmpnm, "TXX" );
	sprintf( sacfile, "%s.txx.grn", grn->filename );
	fp = fopen(sacfile,"w");
	fwrite( &sp, sizeof(Sac_Header), 1, fp );
        fwrite( &txx[0], sp.npts*sizeof(float), 1, fp );
        fclose(fp);

	strcpy( sp.kcmpnm, "TXY" );
        sprintf( sacfile, "%s.txy.grn", grn->filename );
        fp = fopen(sacfile,"w");
        fwrite( &sp, sizeof(Sac_Header), 1, fp );
        fwrite( &txy[0], sp.npts*sizeof(float), 1, fp );
        fclose(fp);

	strcpy( sp.kcmpnm, "TXZ" );
        sprintf( sacfile, "%s.txz.grn", grn->filename );
        fp = fopen(sacfile,"w");
        fwrite( &sp, sizeof(Sac_Header), 1, fp );
        fwrite( &txz[0], sp.npts*sizeof(float), 1, fp );
        fclose(fp);

	strcpy( sp.kcmpnm, "TYY" );
        sprintf( sacfile, "%s.tyy.grn", grn->filename );
        fp = fopen(sacfile,"w");
        fwrite( &sp, sizeof(Sac_Header), 1, fp );
        fwrite( &tyy[0], sp.npts*sizeof(float), 1, fp );
        fclose(fp);

	strcpy( sp.kcmpnm, "TYZ" );
        sprintf( sacfile, "%s.tyz.grn", grn->filename );
        fp = fopen(sacfile,"w");
        fwrite( &sp, sizeof(Sac_Header), 1, fp );
        fwrite( &tyz[0], sp.npts*sizeof(float), 1, fp );
        fclose(fp);

	strcpy( sp.kcmpnm, "RXX" );
        sprintf( sacfile, "%s.rxx.grn", grn->filename );
        fp = fopen(sacfile,"w");
        fwrite( &sp, sizeof(Sac_Header), 1, fp );
        fwrite( &rxx[0], sp.npts*sizeof(float), 1, fp );
        fclose(fp);

	strcpy( sp.kcmpnm, "RXY" );
        sprintf( sacfile, "%s.rxy.grn", grn->filename );
        fp = fopen(sacfile,"w");
        fwrite( &sp, sizeof(Sac_Header), 1, fp );
        fwrite( &rxy[0], sp.npts*sizeof(float), 1, fp );
        fclose(fp);

	strcpy( sp.kcmpnm, "RXZ" );
        sprintf( sacfile, "%s.rxz.grn", grn->filename );
        fp = fopen(sacfile,"w");
        fwrite( &sp, sizeof(Sac_Header), 1, fp );
        fwrite( &rxz[0], sp.npts*sizeof(float), 1, fp );
        fclose(fp);

	strcpy( sp.kcmpnm, "RYY" );
        sprintf( sacfile, "%s.ryy.grn", grn->filename );
        fp = fopen(sacfile,"w");
        fwrite( &sp, sizeof(Sac_Header), 1, fp );
        fwrite( &ryy[0], sp.npts*sizeof(float), 1, fp );
        fclose(fp);

	strcpy( sp.kcmpnm, "RYZ" );
        sprintf( sacfile, "%s.ryz.grn", grn->filename );
        fp = fopen(sacfile,"w");
        fwrite( &sp, sizeof(Sac_Header), 1, fp );
        fwrite( &ryz[0], sp.npts*sizeof(float), 1, fp );
        fclose(fp);

	strcpy( sp.kcmpnm, "RZZ" );
        sprintf( sacfile, "%s.rzz.grn", grn->filename );
        fp = fopen(sacfile,"w");
        fwrite( &sp, sizeof(Sac_Header), 1, fp );
        fwrite( &rzz[0], sp.npts*sizeof(float), 1, fp );
        fclose(fp);

	strcpy( sp.kcmpnm, "ZXX" );
        sprintf( sacfile, "%s.zxx.grn", grn->filename );
        fp = fopen(sacfile,"w");
        fwrite( &sp, sizeof(Sac_Header), 1, fp );
        fwrite( &zxx[0], sp.npts*sizeof(float), 1, fp );
        fclose(fp);

	strcpy( sp.kcmpnm, "ZXY" );
        sprintf( sacfile, "%s.zxy.grn", grn->filename );
        fp = fopen(sacfile,"w");
        fwrite( &sp, sizeof(Sac_Header), 1, fp );
        fwrite( &zxy[0], sp.npts*sizeof(float), 1, fp );
        fclose(fp);

	strcpy( sp.kcmpnm, "ZXZ" );
        sprintf( sacfile, "%s.zxz.grn", grn->filename );
        fp = fopen(sacfile,"w");
        fwrite( &sp, sizeof(Sac_Header), 1, fp );
        fwrite( &zxz[0], sp.npts*sizeof(float), 1, fp );
        fclose(fp);

	strcpy( sp.kcmpnm, "ZYY" );
        sprintf( sacfile, "%s.zyy.grn", grn->filename );
        fp = fopen(sacfile,"w");
        fwrite( &sp, sizeof(Sac_Header), 1, fp );
        fwrite( &zyy[0], sp.npts*sizeof(float), 1, fp );
        fclose(fp);

	strcpy( sp.kcmpnm, "ZYZ" );
        sprintf( sacfile, "%s.zyz.grn", grn->filename );
        fp = fopen(sacfile,"w");
        fwrite( &sp, sizeof(Sac_Header), 1, fp );
        fwrite( &zyz[0], sp.npts*sizeof(float), 1, fp );
        fclose(fp);

	strcpy( sp.kcmpnm, "ZZZ" );
        sprintf( sacfile, "%s.zzz.grn", grn->filename );
        fp = fopen(sacfile,"w");
        fwrite( &sp, sizeof(Sac_Header), 1, fp );
        fwrite( &zzz[0], sp.npts*sizeof(float), 1, fp );
        fclose(fp);

	free(txx);
	free(txy);
	free(txz);
	free(tyy);
	free(tyz);
	
	free(rxx);
        free(rxy);
        free(rxz);
        free(ryy);
        free(ryz);
	free(rzz);

	free(zxx);
        free(zxy);
        free(zxz);
        free(zyy);
        free(zyz);
	free(zzz);
}	
