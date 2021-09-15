#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "../include/cgraph.h"
#include "../include/mt.h"

/**************************/
/*** default line types ***/
/**************************/

#define SOLID   0
#define DASHED  1
#define DOTTED  3
#define DASHDOT 4

/***************************/
/*** default line widths ***/
/***************************/

#define DATA_LINE_WIDTH 0.75
#define SYN_LINE_WIDTH 1.25

#define CENTIMETERS 0  /*** do nothing original amplitudes in cm ***/
#define MILIMETERS 1   /*** multiply by 10 to convert from cm to mm ***/
#define MICRONS  2     /*** multiply by ***/

/**************************************************************/
/*** line color structure and default red and black predefs ***/
/**************************************************************/

typedef struct {
	float gray,red,grn,blu;
} Color ;

		/**   Gray   R   G     B   **/
static Color Black = { 0.0, 0.0,  0.0,  0.0  };
static Color Red   = { 0.3, 1.0,  0.0,  0.0  };
static Color Green = { 0.3, 0.0,  1.0,  0.0  };
static Color Blue  = { 0.3, 0.0,  0.0,  1.0  };
static Color White = { 0.3, 1.0,  1.0,  1.0  };
static Color Gray  = { 0.3, 0.3,  0.3,  0.3  };
static Color Brown = { 0.3, 0.64, 0.16, 0.16 };
static Color Purple= { 0.3, 0.63, 0.13, 0.94 };
static Color Orange= { 0.3, 1.0,  0.64, 0.0  };
static Color Cyan  = { 0.3, 0.0,  1.0,  1.0  };
static Color Yellow= { 0.3, 1.0,  1.0,  0.0  };

char progname[128];

	/***************************/
	/*** the main subroutine ***/
	/***************************/

int psplot(
	int nsta,
	int iz,
	char *filenbase,
	EventInfo *ev, 
	Solution *sol,
	Greens **grn,
	int units,
	int verbose,
	int forward, 
	int PltXcorLabel,
	int PS_Output_Screen )
{

/*****************************/
/*** Cgraph plotting stuff ***/
/*****************************/
	
	extern char progname[128];

	float scale=1.0;   /*** inches=1.0  0.39370=cm ***/
	extern Color Red, Black, Green, Blue, White, Gray, Brown;
	Color line_col;
	float line_width, x0_shift, y0_shift;
	int ltype;
	char filename[256], xlabel[128], ylabel[128], label[128], grd_mo_type[5];
	float xaxsize, yaxsize;
	float cm2mm = 10.;
	float cm2microns = 10000.;
	float *syn_z, *syn_r, *syn_t, *dat_z, *dat_r, *dat_t;
	FILE *fp;
/* see include/mt.h Mo = math.pow( 10.0, 1.5*(Mw+10.73) ) = 1.2445146117713818e+16; Reference Mw = 0.0 */

	char ns, ew;

/*****************************/
/*** other variables       ***/
/*****************************/

	int jcol, irow, it, npts, ista;
	float new_col, x0, y0, x0_new, y0_new, xsize, ysize;
	float dt, beg, xmin, xmax, ymin, ymax;
	float *x, *y;
	int ipage;
	float new_page;

	float xc, yc, rad, stkt, plnt, min_ev, stkb, plnb, med_ev, stkp, plnp, max_ev;
	int planes;
	int ifill;

/******************************************/
/*** initialize some plotting variables ***/
/******************************************/

	void ps_line( float *, float *, int, float, float, int, Color, char * );
	void plot_dc_clvd_mech( float, float, float, Solution *, Greens **, int, int, int, int, int, int );
	void ps_axes( float, float, float, float, float, float, Color, char *, char * );
	void set_color( Color );
	void sac_minmax( float *, int, float *, float *, float * );

/********************************/
/*** loop over stations       ***/
/********************************/

	x0 = 1.25;
	y0 = 1.25;
	xsize = 3.5;
	ysize = 2.5;
	irow = -1;
	jcol = -1;
	ipage = 0;

/************************************/
/*** loop over pages and stations ***/
/************************************/

	for( ista=0; ista<=nsta; ista++ )
	{
		new_page = ista % 9;

		if(verbose)
		{
		  fprintf(stdout, 
		    "%s: psplot.c: new_page=%g ipage=%d ista=%d irow=%d jcol=%d new_col=%g\n",
			progname, new_page, ipage, ista, irow, jcol, new_col );
		}

	/***********************************************************/
	/*** if this is the first or 9th panel, start a new page ***/
	/***********************************************************/

		if( new_page == 0 )
		{
			if(verbose) fprintf( stdout, "%s: psplot.c: newpage\n", progname );

		/****************************************/
		/*** close old plot if not first page ***/
		/****************************************/

			if( ipage != 0 ) cg_showpage();

			ipage++;
			sprintf( filename, "%s.p%02d.ps", filenbase, ipage );
			if(verbose)
			{
			  fprintf( stdout, "%s: psplot.c: Plotting Postscript Plot %s\n", 
				progname, filename );
			}
			remove( filename );

		/*****************************/
		/*** initializing plotting ***/
		/*****************************/

			cg_launch_preview(0);
			cg_set_output_filename( filename );
			cg_setcreator( progname );
			cg_settitle( "Seismic Moment Tensor Inversion" );

			if( PS_Output_Screen == SCREEN_PRINT )
			{
				cg_setboundingbox( "10 10 810 610" );
				cg_init( PORTRAIT, 1, scale );
			}
			else 
			{
				cg_setboundingbox( "0 0 792 612" );
				cg_init( LANDSCAPE, 1, scale );
			}

			cg_coord_select( 1, 1 );
			irow = -1;
			jcol = -1;
		}

		new_col = ista % 3;

		if( new_col == 0 )
		{
			irow++; /*** increment row ***/
			jcol=0; /*** reset column to 0 ***/
		}
		else
		{
			jcol++;
		}

	/**************************/
	/*** set the new origin ***/
	/**************************/

		x0_new = x0 + ( xsize * (float)irow );
		y0_new = y0 + ( ysize * (float)jcol );
		cg_aorigin( x0_new, y0_new );

	/******************************************************/
	/*** plot the axis and set the workspace dimensions ***/
	/******************************************************/

		xaxsize = 2.3;
		yaxsize = 0.5 * ysize;

	/***************************/
	/*** if this is the last ***/
	/*** panel to plot then  ***/
	/*** plot the mechanism  ***/
	/***************************/
		if(verbose)
		{
		  fprintf( stdout, 
			"%s: psplot.c: psplot(): Last panel, mechanism plot next\n", 
				progname );
		  fflush(stdout);
		}

		if( ista == nsta )
		{
			cg_aorigin( x0_new-0.75, y0_new-0.75 );
			cg_axis_enable( 0, 0 );
			cg_linax_style( 0, 0, 7, 1, 6 );
			cg_fontbyname( 12, "Palatino-Roman" );
			cg_xaxis( 2.1, 0, 10, -0.15, 1, 2 );
			cg_yaxis( 2.1, 0, 10, -0.15, 1, 2 );
			cg_fontbyname( 10, "Palatino-Roman" );

		/*****************************************/
		/*** plot the full moment tensor       ***/
		/*****************************************/

			xc = 3.5;
			yc = 9.0;
			rad = 2.7;
			planes = 1;
			ifill = 1;
			plot_dc_clvd_mech( xc, yc, rad, sol, grn, iz, nsta, planes, ifill, FULLMT, verbose );

		/********************************************************/
		/*** plot the major double couple if not an explosion ***/
		/********************************************************/

			if( sol[iz].mt_type != EXPLOSION && forward == 0 )
			{
				xc = 8;
				yc = 8.5;
				rad = 1.5;
				planes = 1;
				ifill = 1;
				/*** nsta = 0 do not plot ***/
				plot_dc_clvd_mech( xc, yc, rad, sol, grn, iz, 0, planes, ifill, MAJORDC, verbose );
				sprintf( label, "Major" );
				cg_text( 7, 10.5, 0, label );
			}
			
			if( forward )
			{
				sprintf( label, "FORWARD CALCULATION" );
				cg_text( 7, 10.5, 0, label );
			}

		/********************************************************/
		/*** check for minor DC and then plot if necessary ***/
		/********************************************************/
			if( fabs( sol[iz].Min.P.ev ) + fabs( sol[iz].Min.B.ev ) + 
				fabs( sol[iz].Min.T.ev ) > 0.0001 )
			{
				xc = 12;
				yc = 8.5;
				rad = 1.5 * ( fabs( sol[iz].epsilon) + 0.3 );
				planes = 1;
				ifill = 1;
				/*** nsta = 0 do not plot takeoff angles ***/
				plot_dc_clvd_mech( xc, yc, rad, sol, grn, iz, 0, planes, ifill, MINORDC, verbose );
				sprintf( label, "Minor" );
				cg_text( 11, 10, 0, label );
			}

		/*************************/
		/*** comment           ***/
		/*************************/

			sprintf( label, "%s", ev[0].comment );
			cg_text( 0, 12, 0, label );
			
			if( grn[0][iz].evla < 0 ) ns = 'S'; else ns='N';
			if( grn[0][iz].evlo < 0 ) ew = 'W'; else ew='E';
			
			sprintf( label, "%4d/%02d/%02d/%02d:%02d:%04.1f %.3f%c %.3f%c %.1f km\n",
				ev[0].ot.year, ev[0].ot.month, ev[0].ot.mday,
				ev[0].ot.hour, ev[0].ot.min, sol[iz].ot,
				grn[0][iz].evla, ns, grn[0][iz].evlo, ew, grn[0][iz].evdp );

			cg_text( 0, 5.2, 0, label );
			
			sprintf( label, "NP1: strike=%.0f Dip=%.0f Rake=%.0f",
				sol[iz].Maj.P1.s, sol[iz].Maj.P1.d, sol[iz].Maj.P1.r );
			cg_text( 0, 4.5, 0, label );

			sprintf( label, "NP2: strike=%.0f Dip=%.0f Rake=%.0f",
				sol[iz].Maj.P2.s, sol[iz].Maj.P2.d, sol[iz].Maj.P2.r );
			cg_text( 0, 3.8, 0, label );

		/**************************/
		/*** moment tensor type ***/
		/**************************/

			if( sol[iz].mt_type == DEVIATORIC )
			{
				sprintf( label, "DEVIATORIC MT DC=%.1f%% CLVD=%.1f%%",
					sol[iz].PDC, sol[iz].PCLVD );	
			}
			else if( sol[iz].mt_type == FULL_MT )	
			{
				sprintf( label, "FULL MT ISO=%.1f%% CLVD=%.1f%% DC=%.1f%%",
					sol[iz].PISO, sol[iz].PCLVD, sol[iz].PDC );
			}
			else if( sol[iz].mt_type == EXPLOSION )
			{
				sprintf( label, "EXPLOSION ONLY" );
			}
			cg_text( 0, 3.1, 0, label );

			/*** reduce clutter ***/
			
			sprintf( label, "VarRed=%.1f%% TOTFIT=%.1f Mw=%.2f",
				sol[iz].var_red, sol[iz].total_fitness1, sol[iz].mw );
			/*
			sprintf( label, "VarRed=%.1f%% Mo=%0.2fE+%02d Mw=%.2f",
				sol[iz].var_red,
				sol[iz].abcassa,
                                sol[iz].exponent,
				sol[iz].mw );
			*/

			cg_text( 0, 2.4, 0, label );

		/*** small fonts ***/
			cg_fontbyname( 8, "Palatino-Roman" );
			sprintf( label, "Mo=%0.2fE+%02d xx=%.2f xy=%.2f xz=%.2f yy=%.2f yz=%.2f zz=%.2f",
		 	  sol[iz].abcassa,
			  sol[iz].exponent,
			  sol[iz].moment_tensor[1][1]/pow(10,sol[iz].exponent),
			  sol[iz].moment_tensor[1][2]/pow(10,sol[iz].exponent),
			  sol[iz].moment_tensor[1][3]/pow(10,sol[iz].exponent),
			  sol[iz].moment_tensor[2][2]/pow(10,sol[iz].exponent),
			  sol[iz].moment_tensor[2][3]/pow(10,sol[iz].exponent),
			  sol[iz].moment_tensor[3][3]/pow(10,sol[iz].exponent) );

			cg_text( 0, 1.5, 0, label );

			if( sol[iz].mt_type == DEVIATORIC )
			{
			  sprintf( label, "T ev=%6.2f az=%3.0f pl=%2.0f Epsilon=%.2f",
			    sol[iz].FullMT.T.ev, sol[iz].Maj.T.az, sol[iz].Maj.T.pl, sol[iz].epsilon );
			  cg_text( 0, 1.0, 0, label );

			  sprintf( label, "B ev=%6.2f az=%3.0f pl=%2.0f",
			   sol[iz].FullMT.B.ev, sol[iz].Maj.B.az, sol[iz].Maj.B.pl );
			  cg_text( 0, 0.5, 0, label );

			  sprintf( label, "P ev=%6.2f az=%3.0f pl=%2.0f",
				sol[iz].FullMT.P.ev, sol[iz].Maj.P.az, sol[iz].Maj.P.pl );
			  cg_text( 0, 0.0, 0, label );

			}
			else if( sol[iz].mt_type == FULL_MT )
			{
			  sprintf( label, "T ev=%6.2f az=%3.0f pl=%2.0f k=%.2f Epsilon=%.2f",
			    sol[iz].FullMT.T.ev, sol[iz].Maj.T.az, sol[iz].Maj.T.pl, 
				 sol[iz].k, sol[iz].epsilon );
			  cg_text( 0, 1.0, 0, label );

			  sprintf( label, "B ev=%6.2f az=%3.0f pl=%2.0f F-factor=%.2f", 
			    sol[iz].FullMT.B.ev, sol[iz].Maj.B.az, sol[iz].Maj.B.pl,
				sol[iz].f_factor );
			  cg_text( 0, 0.5, 0, label ); 

			  sprintf( label, "P ev=%6.2f az=%3.0f pl=%2.0f Lune lat=%.2f lon=%.2f",
				sol[iz].FullMT.P.ev, sol[iz].Maj.P.az, sol[iz].Maj.P.pl,
				sol[iz].lune_lat, sol[iz].lune_lon );
			  cg_text( 0, 0.0, 0, label );

		 	}
			else if( sol[iz].mt_type == EXPLOSION )
			{
			  sprintf( label, "T ev=%6.2f az=%3.0f pl=%2.0f",
                            sol[iz].FullMT.T.ev, sol[iz].Maj.T.az, sol[iz].Maj.T.pl );
                          cg_text( 0, 1.0, 0, label );
                                                                                                                     
                          sprintf( label, "B ev=%6.2f az=%3.0f pl=%2.0f",
                           sol[iz].FullMT.B.ev, sol[iz].Maj.B.az, sol[iz].Maj.B.pl );
                          cg_text( 0, 0.5, 0, label );

			  sprintf( label, "P ev=%6.2f az=%3.0f pl=%2.0f",
				sol[iz].FullMT.P.ev, sol[iz].Maj.P.az, sol[iz].Maj.P.pl );
			  cg_text( 0, 0.0, 0, label );
			}
			
			cg_stroke();

			break;
		} /** done with plotting the psmech panel ***/

	/***************************/
	/*** set the time values ***/
	/***************************/
		npts = ev[ista].ew.s.npts;
		dt   = ev[ista].ew.s.delta;
		beg  = ev[ista].ew.s.b;
		beg  = 0;

	/**************************/
	/*** convert amplitudes ***/
	/**************************/
		if(verbose)
		{
		  fprintf( stdout, 
			"%s: psplot.c: allocating local memory for plotting \n", progname );
		}
		syn_z = (float *) malloc( npts * sizeof(float) );
		syn_r = (float *) malloc( npts * sizeof(float) );
		syn_t = (float *) malloc( npts * sizeof(float) );
		dat_z = (float *) malloc( npts * sizeof(float) );
		dat_r = (float *) malloc( npts * sizeof(float) );
		dat_t = (float *) malloc( npts * sizeof(float) );
		x     = (float *) malloc( npts * sizeof(float) );

		if(verbose)
		{
		  fprintf( stdout, "%s: psplot.c: DONE allocating local memory for plotting\n", progname );
		}

		for( it=0; it<npts; it++ )
		{
			dat_t[it] = ev[ista].ew.data[it]    * cm2microns;
			dat_r[it] = ev[ista].ns.data[it]    * cm2microns;
			dat_z[it] = ev[ista].z.data[it]     * cm2microns;
			syn_t[it] = ev[ista].syn_t.data[it] * cm2microns;
			syn_r[it] = ev[ista].syn_r.data[it] * cm2microns;
			syn_z[it] = ev[ista].syn_z.data[it] * cm2microns;
			x[it] = beg + (float)it * dt;
		}

	/*** set the syn minmax seperate in the SAC header ***/

		sac_minmax( dat_t, npts, &(ev[ista].ew.s.depmax), &(ev[ista].ew.s.depmin), &(ev[ista].ew.s.depmen) );
		sac_minmax( dat_r, npts, &(ev[ista].ns.s.depmax), &(ev[ista].ns.s.depmin), &(ev[ista].ns.s.depmen) );
		sac_minmax( dat_z, npts, &(ev[ista].z.s.depmax), &(ev[ista].z.s.depmin), &(ev[ista].z.s.depmen) );
                sac_minmax( syn_t, npts, &(ev[ista].syn_t.s.depmax), &(ev[ista].syn_t.s.depmin), &(ev[ista].syn_t.s.depmen) );
                sac_minmax( syn_r, npts, &(ev[ista].syn_r.s.depmax), &(ev[ista].syn_r.s.depmin), &(ev[ista].syn_r.s.depmen) );
                sac_minmax( syn_z, npts, &(ev[ista].syn_z.s.depmax), &(ev[ista].syn_z.s.depmin), &(ev[ista].syn_z.s.depmen) );
		xmin = beg;
		xmax = beg+npts*dt;

		if(verbose)
		{
		 fprintf( stdout, "%s: ista=%d data t max=%g min=%g mean=%g\n", progname, ista, ev[ista].ew.s.depmax, ev[ista].ew.s.depmin, ev[ista].ew.s.depmen );
		 fprintf( stdout, "%s: ista=%d data r max=%g min=%g mean=%g\n", progname, ista, ev[ista].ns.s.depmax, ev[ista].ns.s.depmin, ev[ista].ns.s.depmen );
		 fprintf( stdout, "%s: ista=%d data z max=%g min=%g mean=%g\n", progname, ista, ev[ista].z.s.depmax,  ev[ista].z.s.depmin,  ev[ista].z.s.depmen );

		 fprintf( stdout, "%s: ista=%d syn  t max=%g min=%g mean=%g\n", progname, ista, ev[ista].syn_t.s.depmax, ev[ista].syn_t.s.depmin, ev[ista].syn_t.s.depmen );
		 fprintf( stdout, "%s: ista=%d syn  r max=%g min=%g mean=%g\n", progname, ista, ev[ista].syn_r.s.depmax, ev[ista].syn_r.s.depmin, ev[ista].syn_r.s.depmen );
		 fprintf( stdout, "%s: ista=%d syn  z max=%g min=%g mean=%g\n", progname, ista, ev[ista].syn_z.s.depmax, ev[ista].syn_z.s.depmin, ev[ista].syn_z.s.depmen );
		}

	/********************************************************/
	/*** plot the axis and set the workspace dimensions   ***/
	/*** estimate the total ymax size from all components ***/
	/********************************************************/

		line_width = 1.5;
		ymin = 0.0;
		ymax =  ( fabs( ev[ista].ew.s.depmin ) + fabs( ev[ista].ew.s.depmax ) ) +
			( fabs( ev[ista].ns.s.depmin ) + fabs( ev[ista].ns.s.depmax ) ) +
			( fabs( ev[ista].z.s.depmin ) + fabs( ev[ista].z.s.depmax ) );

	/*******************************************/
	/*** X and Y axis labels                 ***/
	/*******************************************/

		if( ev[ista].grd_mo_type == DISPLACEMENT ) 
			sprintf( ylabel, "DISP (microns)" );
		else if( ev[ista].grd_mo_type == VELOCITY ) 
			sprintf( ylabel, "VEL (microns/sec)" );
		else	
			sprintf( ylabel, "unknown" );

		if( ev[ista].redv > 0 )
			sprintf( xlabel, "reduced travel time [T-R/%4.1f] (sec)",
				ev[ista].redv );
		else
			sprintf( xlabel, "travel time (sec)" );
	
		ps_axes( xmin, xmax, ymax, xaxsize, yaxsize, line_width, Black, xlabel, ylabel );
	
	/*************************************************************/
	/*** plot data and syn on workspace give above predef axes ***/
	/*************************************************************/

	/*** blue = station is not used in inversion Red station is used ***/
		if( ev[ista].iused == 0 ) line_col = Blue;
		if( ev[ista].iused == 1 ) line_col = Red;
		if( ev[ista].iused == -1 ) line_col = Yellow;

	/*** Transverse component synthetics and data ***/
		y0_shift = fabs( ev[ista].ew.s.depmin );
		ps_line(x, dat_t, npts, y0_shift, DATA_LINE_WIDTH, SOLID,  Black,    "T"  );
		ps_line(x, syn_t, npts, y0_shift, SYN_LINE_WIDTH,  DOTTED, line_col, "\0" );

		cg_fontbyname( 10, "NewCenturySchlbk-Roman" );
		set_color( Brown );
		/**** too much clutter and font size too small to read ***/
		/* sprintf( label, "%g/%.2f/%.1f", ev[ista].ttlag, ev[ista].txcor, ev[ista].ew.P2P_snr ); */
		sprintf( label, "%g/%.2f", ev[ista].ttlag, ev[ista].txcor );
		if(PltXcorLabel) cg_text( x[0], y0_shift + fabs(ev[ista].ew.s.depmax)/2, 0, label );

	/*** radial component ***/
		y0_shift = y0_shift + 1.1*( fabs( ev[ista].ns.s.depmin ) + fabs( ev[ista].ew.s.depmax ) );
		ps_line( x, dat_r, npts, y0_shift, DATA_LINE_WIDTH,  SOLID, Black,    " R" ); /*Dat*/
                ps_line( x, syn_r, npts, y0_shift, SYN_LINE_WIDTH,  DOTTED, line_col, "\0" ); /*Syn*/

		cg_fontbyname( 10, "NewCenturySchlbk-Roman" );
		set_color( Brown );
		/** sprintf( label, "%g/%.2f/%.1f", ev[ista].rtlag, ev[ista].rxcor, ev[ista].ns.P2P_snr ); **/
		sprintf( label, "%g/%.2f", ev[ista].rtlag, ev[ista].rxcor );
		if(PltXcorLabel) cg_text( x[0], y0_shift + fabs(ev[ista].ns.s.depmax)/2, 0, label );

	/*** vertical component ***/
		y0_shift = y0_shift + 1.1*( fabs( ev[ista].z.s.depmin ) + fabs( ev[ista].ns.s.depmax ) );
		ps_line( x, dat_z, npts, y0_shift, DATA_LINE_WIDTH,  SOLID, Black,    " Z" ); /*Dat*/
		ps_line( x, syn_z, npts, y0_shift, SYN_LINE_WIDTH,  DOTTED, line_col, "\0" ); /*Syn*/

		cg_fontbyname( 10, "NewCenturySchlbk-Roman" );
		set_color( Brown );
		/** sprintf( label, "%g/%.2f/%.1f", ev[ista].ztlag, ev[ista].zxcor, ev[ista].z.P2P_snr ); **/
		sprintf( label, "%g/%.2f", ev[ista].ztlag, ev[ista].zxcor );
		if(PltXcorLabel) cg_text( x[0], y0_shift + fabs(ev[ista].z.s.depmax)/2, 0, label );

	/*** plot the header label ***/
		y0_shift = ymax + ymax*0.2;
		if( y0_shift < 0.1 ) y0_shift = ymax + 0.6*ymax;

		x0_shift = -0.25 * (x[npts-1] - x[0]);

		sprintf( grd_mo_type, "D" );
		if( ev[ista].grd_mo_type == VELOCITY ) sprintf( grd_mo_type, "V" );

		/*** reduce clutter ***/
		/*
		sprintf(label, "%s.%s(%s) R=%.0f km Az=%.0f ts=%g wt=%.0f",
			grn[ista][iz].stnm, grn[ista][iz].net, grd_mo_type,
			grn[ista][iz].rdist, grn[ista][iz].az,
			ev[ista].time_shift_all, ev[ista].weight );
		*/
		sprintf(label, "%s.%s R=%.0f km Az=%.0f ts=%g",
			grn[ista][iz].stnm, grn[ista][iz].net, 
                        grn[ista][iz].rdist, grn[ista][iz].az,
			ev[ista].time_shift_all
		 );

		set_color( Black );
		cg_fontbyname( 12, "Palatino-Roman" );
		cg_text( x0_shift, y0_shift, 0, label );
		cg_stroke();

		if(verbose)
		{
		  fprintf(stdout, "%s: psplot.c: freeing local memory in psplot... ", progname );
		  fflush(stdout);
		}

		free(x);
		free(syn_z);
		free(syn_r);
		free(syn_t);
		free(dat_z);
		free(dat_r);
		free(dat_t);

		if(verbose)
		{
		  fprintf( stdout, " DONE.\n" );
		  fflush(stdout);
		}

		if(verbose)
		{
		  fprintf( stdout, "%s: psplot.c: psplot(): looper->ista=%d new_page=%g ipage=%d irow=%d jcol=%d new_col=%g\n",
                       	progname, ista, new_page, ipage, irow, jcol, new_col );
		  fflush(stdout);
		}

	} /*** loop over pages and stations ***/

	if(verbose)
	{
	  fprintf( stdout, "%s: psplot.c: done loop\n", progname );
	  fflush(stdout);
	}

	cg_showpage();

	if(verbose)
	{
		fprintf( stdout, "%s: psplot.c: psplot(): done with psplot().\n", progname );
	}

	return ipage;
}

void ps_axes( float xmin, float xmax, float ymax, float xax_size, float yax_size, 
	float line_width, Color col, char *xlabel, char *ylabel )
{
	int i;
	int numstep;
	float xlen, ylen, xinc, yinc;
	float y_smallest=0.01;
	float Y_MAX;
	void set_color( Color );
	
/*************************************************/
/*** text marker 0=inch 1=user axis            ***/
/*************************************************/
	cg_axis_enable( 1, 1 );

/*************************************************/
/*** numdist = (fontsize/10) * numdist         ***/
/*** tnpos = 1 ticks outside numbers outside   ***/
/*** ticlen = units of 1/72 inch               ***/
/*** precision, numoff, numdist, tnpos, ticlen ***/
/*************************************************/
	cg_linax_style( 1, 0, 2, 1, 3 );

/*************************************************/
/*** size min max offset ticsep numsep(int)    ***/
/*************************************************/
	cg_linewidth( line_width );
	set_color( col );
	cg_fontbyname( 11, "Palatino-Roman" );

/*** x axes ***/
	xlen = 10;
	if( xmax <= 150 )
		xlen = roundf(xmax/10.0)*10.0;
	else
		xlen = roundf(xmax/100.0)*100.0;
/*
	for( i=0; i<150; i++ )
	{
		xlen += 10;
		if( xlen > xmax ) break;
	}
*/
	xinc = 5.;
	if( xlen >= 100 ) xinc = 25;
	if( xlen >= 200 ) xinc = 50;
	if( xlen >= 400 ) xinc = 100;

/*** y axes ***/

	Y_MAX = 0.0;
	while( Y_MAX <= ymax )
	{
		Y_MAX += y_smallest;
	}
	ylen = Y_MAX + y_smallest;
	ylen = roundf( ylen*100.0 )/100.0;
	yinc = ylen/2;

/*******************************************************/
/***     size    min   max    offset ticstep  numsep ***/
/*******************************************************/
	numstep = 2;
	cg_xaxis( xax_size, xmin, xlen, -0.15,  xinc, numstep );
	cg_yaxis( yax_size,    0, ylen, -0.15,  yinc, numstep );
	cg_xlabel( xlabel );
	cg_ylabel( ylabel );

	/* no grid lines */
	/* cg_mesh(); */
}

void ps_line( float *x, float *y, int npts, float y0_shift, float line_width, 
	int ltype, Color lcolor, char *label )
{
	int i;
	void set_color( Color );

/*************************************/
/*** set line type                 ***/
/*************************************/
	if( ltype==SOLID )	cg_dash( 0, 1.0 );
	if( ltype==DASHED )	cg_dash( 1, 0.5 );
	if( ltype==DOTTED )	cg_dash( 4, 0.1 );
	if( ltype==DASHDOT )	cg_dash( 2, 0.5 );

/*************************************/
/*** set line width                ***/
/*************************************/
	cg_linewidth( line_width );

/*************************************/
/*** set line color                ***/
/*************************************/
	set_color( lcolor );

/*************************************/
/*** draw the line                 ***/
/*************************************/
	cg_move( x[0], y[0]+y0_shift );
	for( i=1; i<npts; i++ )
	{
		cg_line( x[i], y[i]+y0_shift );
	}
	cg_stroke();

/*************************************/
/*** check if label exists         ***/
/*************************************/
	/* if( label == NULL ) return; */

/*************************************/
/*** make a label and reset color  ***/
/*************************************/
	if( label[0] != '\0' )
	{ 
		cg_fontbyname( 12, "Times-Bold" );
		cg_text( x[npts-1], y[npts-1]+y0_shift, 0, label );
		cg_stroke();
	}
}

void set_color( Color col )
{
	cg_grayrgbcolor( col.gray, col.red, col.grn, col.blu );
}

/*** font names 
	Helvetica
	Helvetica-Bold
	Helvetica-Oblique 
	Helvetica-BoldOblqiue
	Times-Roman
	Times-Bold
	Times-Italic
	Times-BoldItalic
	Courier
	Courier-Bold
	Courier-Oblqiue
	Courier-BoldOblique
	Symbol
	Helvetica-Narrow-Bold
	Helvetica-Narrow-Oblique
	Helvetica-Narrow-BoldOblique
	NewCentrurySchlbk-Roman
	NewCentrurySchlbk-Italic
	NewCentrurySchlbk-Bold
	NewCentrurySchlbk-BoldItalic
	Palatino-Roman
	Palatino-Italic
	Palatino-Bold
	Palatino-BoldItalic
***/
