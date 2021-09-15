#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <errno.h>
#include "../include/nrutil.h"     /** numerical recipes **/
#include "../include/mt.h"         /** global datatype and structure declarations **/

char progname[128];

void demultiplex( Greens **grn, EventInfo *ev, int iz, float *s_vector, int nsta )
{
        int ista, it, nt, ir;
        ir = 1;

        for( ista=0; ista<nsta; ista++ )
        {
		if( ev[ista].iused != 1 ) continue;

                nt = grn[ista][iz].nt;

                for( it=0; it<nt; it++)
                {
                        grn[ista][iz].tra[it] = s_vector[ir++];
                        ev[ista].syn_t.data[it] = grn[ista][iz].tra[it];
                }
                for( it=0; it<nt; it++)
                {
                        grn[ista][iz].rad[it] = s_vector[ir++];
                        ev[ista].syn_r.data[it] = grn[ista][iz].rad[it];
                }
                for( it=0; it<nt; it++)
                {
                        grn[ista][iz].ver[it] = s_vector[ir++];
                        ev[ista].syn_z.data[it] = grn[ista][iz].ver[it];
                }
        }
}

/****************************************************************************************************/
/****************************************************************************************************/
/**** call after main(): invert(): compute_synthetics(): with iz = izbest                       ***/
/*** since compute_synthetics works on both ev[ista].iused = 0 and 1 then we can calc           ***/
/*** on predictions too.  ***/
/****************************************************************************************************/
/****************************************************************************************************/

void computeStationChannel_VarianceReduction_v2(
	int nsta,
	int iz, /* iz = iz_best */
	EventInfo *ev,
	Greens **grn,
	Solution *sol,
	int verbose )
{
	int it, ix, nt, ista;
	float *tdat, *rdat, *zdat;
        float *tsyn, *rsyn, *zsyn;
        float *dat, *syn;

	FILE *fp;
	static char output_filename[] = { "var_red.out" };

	float variance_reduction( float *d, float *s, int n0, int n1 );

/*** start ***/
	fp = fopen( output_filename, "a+" );

	fprintf( fp,
"#   sta     net       evdp otsec ts0    dist  az  flo    fhi  used vr_ev  vr_sta  vr_z   vr_r   vr_t        snr_z        snr_r        snr_t\n" );
	fflush( fp );

	fprintf( stdout,
	  "%s: %s: %s: Event.Station.Channel Variance Reduction for best depth z= %g and ts0= %g otsec= %g var_red= %6.1f\n",
                progname,
		__FILE__,
		__func__,
		grn[0][iz].evdp,
		ev[0].ts0, 
		ev[0].ot_shift, 
		sol[iz].var_red );
	fflush( stdout );

/*** mtinv_subs.c: compute_synthetics()  start at ( it = 0; it < nt; it++ )***/
/*** mtinv_subs.c: load_the_data()  fread(&(ev[ista].z.data[0]),  ev[ista].z.s.npts  * sizeof(float), 1, fpd ); ***/

	for( ista = 0; ista < nsta; ista++ )
	{
		/* if( ev[ista].iused != 1 ) continue; */

		ix = 0;
		nt = grn[ista][iz].nt;

		tdat = calloc( nt+1, sizeof(float) );
                rdat = calloc( nt+1, sizeof(float) );
                zdat = calloc( nt+1, sizeof(float) );
                tsyn = calloc( nt+1, sizeof(float) );
                rsyn = calloc( nt+1, sizeof(float) );
                zsyn = calloc( nt+1, sizeof(float) );
		dat = calloc( 3*(nt+1), sizeof(float) );
		syn = calloc( 3*(nt+1), sizeof(float) );

		for( it = 0; it < nt; it++ )
		{
        		tsyn[it] = ev[ista].syn_t.data[it];
        		tdat[it] = ev[ista].ew.data[it];
        		syn[ix]  = ev[ista].syn_t.data[it];
        		dat[ix]  = ev[ista].ew.data[it];
        		ix++;
		}
		ev[ista].vred_tcmp = variance_reduction( tdat, tsyn, 0, nt );

		for( it = 0; it < nt; it++ )
                {
                        rsyn[it] = ev[ista].syn_r.data[it];
                        rdat[it] = ev[ista].ns.data[it];
                        syn[ix]  = ev[ista].syn_r.data[it];
                        dat[ix]  = ev[ista].ns.data[it];
                        ix++;
                }
                ev[ista].vred_rcmp = variance_reduction( rdat, rsyn, 0, nt );

		for( it = 0; it < nt; it++ )
                {
                        zsyn[it] = ev[ista].syn_z.data[it];
                        zdat[it] = ev[ista].z.data[it];
                        syn[ix]  = ev[ista].syn_z.data[it];
                        dat[ix]  = ev[ista].z.data[it];
                        ix++;
                }
                ev[ista].vred_zcmp = variance_reduction( zdat, zsyn, 0, nt );

		ev[ista].vred_sta =  variance_reduction( dat, syn, 0, (3*nt) );

                fprintf( stdout,
  "%s: %s: %s: ista= %2d %8s %8s iz= %2d z= %g otsec= %g iused= %2d vred_ev = %6.1f vred_sta= %6.1f z= %6.1f r= %6.1f t= %6.1f\n",
                        progname,
                        __FILE__,
                        __func__,
			ista, 
			ev[ista].z.s.kstnm,
			ev[ista].z.s.knetwk,
			iz,
                        grn[ista][iz].evdp,
			ev[0].ot_shift, 
			ev[ista].iused,
			sol[iz].var_red,
                        ev[ista].vred_sta,
                        ev[ista].vred_zcmp,
                        ev[ista].vred_rcmp,
                        ev[ista].vred_tcmp );

		/*** sta net evdp otsec ts0 dist az flo fhi iused vred_ev vred_sta vred_cmpz vred_cmpr vred_cmpt snr_z snr_r snr_t **/
		fprintf( fp, 
			"%03d %-8s %-8s %4g %3g %3g %8.1f %3.0f %6.3f %6.3f %2d %6.1f %6.1f %6.1f %6.1f %6.1f %12.5f %12.5f %12.5f\n",
			ista,
			ev[ista].z.s.kstnm,
			ev[ista].z.s.knetwk,
			grn[ista][iz].evdp,
			ev[0].ot_shift,
			ev[0].ts0,
			ev[ista].rdistkm,
			ev[ista].az,
			ev[ista].lf,
			ev[ista].hf,
			ev[ista].iused,
                        sol[iz].var_red,
                        ev[ista].vred_sta,
                        ev[ista].vred_zcmp,
                        ev[ista].vred_rcmp,
                        ev[ista].vred_tcmp,
			ev[ista].z.P2P_snr,
			ev[ista].ns.P2P_snr,
			ev[ista].ew.P2P_snr );
/*
                free(dat);
                free(syn);
                free(rdat);
                free(rsyn);
                free(zdat);
                free(zsyn);
                free(tdat);
                free(tsyn);
	
*/
	} /*** loop over each ista ***/

	fclose(fp);

} /*** end of subroutine: computeStationChannel_VarianceReduction_v2()  ***/



/****************************************************************************************************/
/****************************************************************************************************/
/*** this only works with s_vector which is only available after invert() call from main         ***/
/*** and EventInfo ev structure is overwritten at each depth so useless                         ****/
/*** create new version 2 which can be called after call to invert() and compute_synthetics()    ***/
/*** to do the data and syn from the best solution iz=iz_best                                    ***/
/****************************************************************************************************/
/****************************************************************************************************/

void computeStationChannel_VarianceReduction( 
	EventInfo *ev, Greens **grn, int nsta, int iz, float *s_vector, 
	float var_red, int verbose )
{
	int it, nt, ista;
	int ir, ix;

	float *tdat, *rdat, *zdat;
	float *tsyn, *rsyn, *zsyn;
	float *dat, *syn;

	float variance_reduction( float *d, float *s, int n0, int n1 );

	fprintf( stdout, "%s: %s: %s: event variance reduction iz=%d var_red=%.1f\n",
		progname, __FILE__, __func__, iz, var_red );

	ir = 1;
	for( ista = 0; ista < nsta; ista++ )
	{
		if( ev[ista].iused != 1 ) continue;

		ix = 1;
		nt = grn[ista][iz].nt;

		tdat = calloc( nt+1, sizeof(float) );
		rdat = calloc( nt+1, sizeof(float) );
		zdat = calloc( nt+1, sizeof(float) );
		dat = calloc( 3*(nt+1), sizeof(float) );

		tsyn = calloc( nt+1, sizeof(float) );
		rsyn = calloc( nt+1, sizeof(float) );
		zsyn = calloc( nt+1, sizeof(float) );
		syn  = calloc( 3*(nt+1), sizeof(float) );

		for( it=1; it<=nt; it++)
		{
			tsyn[it] = s_vector[ir];
			tdat[it] = ev[ista].ew.data[it];

			syn[ix]  = s_vector[ir];
			dat[ix]  = ev[ista].ew.data[it];
			ir++;
			ix++;
		}
/* fprintf( stdout, "%s: %s: %s: ista=%d ir=%d ix=%d nt=%d\n", progname, __FILE__, __func__, ista, ir, ix, nt ); */

		ev[ista].vred_tcmp = variance_reduction( tdat, tsyn, 0, nt );
		
		for( it=1; it<=nt; it++ )
		{
			rsyn[it] = s_vector[ir];
			rdat[it] = ev[ista].ns.data[it];

                        syn[ix]  = s_vector[ir];
                        dat[ix]  = ev[ista].ns.data[it];
                        ir++;
			ix++;
		}
/* fprintf( stdout, "%s: %s: %s: ista=%d ir=%d ix=%d nt=%d\n", progname, __FILE__, __func__, ista, ir, ix, nt ); */

		ev[ista].vred_rcmp = variance_reduction( rdat, rsyn, 0, nt );

		for( it=1; it<=nt; it++ )
                {
                        zsyn[it] = s_vector[ir];
                        syn[ix]  = s_vector[ir];
                        zdat[it] = ev[ista].z.data[it];
                        dat[ix]  = ev[ista].z.data[it];
                        ir++;
			ix++;
                }
/* fprintf( stdout, "%s: %s: %s: ista=%d ir=%d ix=%d nt=%d\n", progname, __FILE__, __func__, ista, ir, ix, nt ); */

                ev[ista].vred_zcmp = variance_reduction( zdat, zsyn, 0, nt );

		ev[ista].vred_sta =  variance_reduction( dat, syn, 0, (3*nt) );

		fprintf( stdout, "%s: %s: %s: iz=%02d ista=%02d sta=%4s vred_sta=%.1f z=%.1f r=%.1f t=%.1f\n",
			progname,
			__FILE__,
			__func__,
			iz,
			ista,
			ev[ista].z.s.kstnm,
			ev[ista].vred_sta,
			ev[ista].vred_zcmp,
			ev[ista].vred_rcmp,
			ev[ista].vred_tcmp );

		free(dat);
		free(syn);
		free(rdat);
		free(rsyn);
		free(zdat);
		free(zsyn);
		free(tdat);	
		free(tsyn);
	}
}

void mtinv_dumpSAC( EventInfo *ev, Greens **grn, int nsta, int iz, float *s_vector, int verbose )
{
	int ista;
	char sac_file_name[128];

	void demultiplex( Greens **grn, EventInfo *ev, int iz, float *s_vector, int nsta );
	void sac_minmax( float *, int, float *, float *, float * );
	void write_sac_file( char *, Sac_File *, int );

	for( ista = 0; ista < nsta; ista++ )
	{
		if( ev[ista].iused != 1 ) continue;

		fprintf( stdout, "%s: %s: %s: ev.sta=%s grn.sta=%s ista=%d iz=%d ev.nt=%d grn.nt=%d\n",
			progname, __FILE__, __func__,
			ev[ista].z.s.kstnm, grn[ista][iz].stnm, ista, iz,
			ev[ista].z.s.npts, grn[ista][iz].nt );

		grn[ista][iz].tra = calloc( grn[ista][iz].nt, sizeof(float) );
		grn[ista][iz].rad = calloc( grn[ista][iz].nt, sizeof(float) );
		grn[ista][iz].ver = calloc( grn[ista][iz].nt, sizeof(float) );
		ev[ista].syn_r.data = calloc( ev[ista].ns.s.npts, sizeof(float) );
		ev[ista].syn_t.data = calloc( ev[ista].ew.s.npts, sizeof(float) );
		ev[ista].syn_z.data = calloc( ev[ista].z.s.npts,  sizeof(float) );

	}

	demultiplex( grn, ev, iz, s_vector, nsta );

	for( ista = 0; ista < nsta; ista++ )
	{

	  if( ev[ista].iused != 1 ) continue;

/*** synthetic radial component ***/
	  sprintf( sac_file_name, "%s.%s.%02d.%02d.syn.r.sac", ev[ista].stnm, ev[ista].net, iz, ista );

	  ev[ista].syn_r.s = ev[ista].ns.s;

	  sac_minmax( ev[ista].syn_r.data, ev[ista].syn_r.s.npts,
                    &(ev[ista].syn_r.s.depmax), &(ev[ista].syn_r.s.depmin),
                    &(ev[ista].syn_r.s.depmen) );

	  write_sac_file( sac_file_name, &(ev[ista].syn_r), verbose );

/*** synthetic vertical component ***/
	  sprintf( sac_file_name, "%s.%s.%02d.%02d.syn.z.sac", ev[ista].stnm, ev[ista].net, iz, ista );

	  ev[ista].syn_z.s = ev[ista].z.s;

	  sac_minmax( ev[ista].syn_z.data, ev[ista].syn_z.s.npts,
                    &(ev[ista].syn_z.s.depmax), &(ev[ista].syn_z.s.depmin),
                    &(ev[ista].syn_z.s.depmen) );
	  write_sac_file( sac_file_name, &(ev[ista].syn_z), verbose );

/*** synthetic transverse component ***/
	  sprintf( sac_file_name, "%s.%s.%02d.%02d.syn.t.sac", ev[ista].stnm, ev[ista].net, iz, ista ); 

	  ev[ista].syn_t.s = ev[ista].ew.s;

	  sac_minmax( ev[ista].syn_t.data, ev[ista].syn_t.s.npts,
                    &(ev[ista].syn_t.s.depmax), &(ev[ista].syn_t.s.depmin),
                    &(ev[ista].syn_t.s.depmen) );

	  write_sac_file( sac_file_name, &(ev[ista].syn_t), verbose );

/*** data radial componet ***/
	  sprintf( sac_file_name, "%s.%s.%02d.%02d.dat.r.sac", ev[ista].stnm, ev[ista].net, iz, ista );
	  write_sac_file( sac_file_name, &(ev[ista].ns), verbose );
	  
/*** data vertical component ***/
	  sprintf( sac_file_name, "%s.%s.%02d.%02d.dat.z.sac", ev[ista].stnm, ev[ista].net, iz, ista ); 
	  write_sac_file( sac_file_name, &(ev[ista].z), verbose );

/*** data transverse component ***/
	  sprintf( sac_file_name, "%s.%s.%02d.%02d.dat.t.sac", ev[ista].stnm, ev[ista].net, iz, ista ); 
	  write_sac_file( sac_file_name, &(ev[ista].ew), verbose );
	}

/***
	for( ista = 0; ista < nsta; ista++ )
	{
		if( ev[ista].iused != 1 ) continue;

		free( ev[ista].syn_r.data );
		free( ev[ista].syn_r.data );
		free( ev[ista].syn_r.data );
		free( grn[ista][iz].rad );
		free( grn[ista][iz].tra );
		free( grn[ista][iz].ver );
	}
***/

}
