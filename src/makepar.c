#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/utsname.h>

#include "../include/mt.h"

char progname[128];

#define MAX_STA 128

int main( int ac, char **av )
{
	Greens **grn;
	int nz, iz;
	float *z;
	int i, ix, ista;
	char comment[256];
	char dateid[256];
	FILE *fp;
	int it, nt, nsta;
	float min,max;
	int verbose=0;
	float *rdistkm;
	int *index;
	int indexx( int, float *, int * );
	float my_dt, lf=0.01, hf=0.1;
	int my_nt;
	float minsnr = 3.0;
	float ctol = 0.85;
	float maxshift = 10.0;
	int igmt5 = 1;
	char gmtstring[8]; /* nogmt5 gmt5 */
	int ifullmt = 0;
	float LOWEST_VELOCITY = 2.35;

	int irealtime = 0;
	int ioracle = 0;
	int imysql = 0;
	int isqlite = 1;
	char dbprog[8]; /* oracle mysql sqlite */

	char RespDir[512], DataDir[512];
	char use;
	char com[3];

	void select_dt_and_nt( float dist, float vel, int grn_npts, float grn_dt, float *dat_dt, int *dat_nt );

	int setpar(int,char**),getpar(),mstpar();
	void endpar();

	void Print_Usage(void);
	extern char progname[128];

	struct utsname my_uname;
	int ihelp = 0;

/*** start main ***/

	strcpy( progname, av[0] );
	fprintf( stderr, "%s: version=%s release date=%s\n",
		progname, Version_Label, Version_Date );

	if( uname(&my_uname) == -1 )
	{
		fprintf( stderr, "%s: uname() failed\n", progname );
	}
	else
	{
	/* Darwin Linux SunOS */

		fprintf( stdout, "%s: System name: %s\n", 
			progname, my_uname.sysname );

	/* ppc64 i386 i686 x86_64 amd64 ia64 sun4u */

		fprintf( stdout, "%s: Machine name: %s  Node Name: %s\n",
                        progname, my_uname.machine, my_uname.nodename );
	}

/*** defaults ***/
	strcpy( comment, "New Region" );
	strcpy( dateid, "2008/01/01,00:00:00.00" );
	strcpy( DataDir, "../Data" );
	strcpy( RespDir, "../Resp" );

	if( ac <= 1 ) 
	{
		Print_Usage();
	}

	setpar(ac,av);

	getpar( "com", "s", comment );
	getpar( "date", "s", dateid );
	getpar( "DataDir", "s", DataDir );
	getpar( "RespDir", "s", RespDir );

	getpar( "verbose", "b", &verbose );
	getpar( "minsnr", "f", &minsnr );
	getpar( "ctol", "f", &ctol );
	getpar( "maxshift", "f", &maxshift );
	getpar( "lf", "f", &lf );
	getpar( "hf", "f", &hf );
	getpar( "gmt5", "b", &igmt5 );
	getpar( "fullmt", "b", &ifullmt );
	getpar( "realtime", "b", &irealtime );

	getpar( "oracle", "b", &ioracle );
	getpar( "mysql", "b", &imysql );
	getpar( "sqlite", "b", &isqlite );
	getpar( "help", "b", &ihelp );
	endpar();

	if(ihelp)Print_Usage();

/*** set some flags ***/

	strcpy( gmtstring, "nogmt5" ); /*** default ***/
	if(igmt5) 
          strcpy( gmtstring, "gmt5" );
        else
          strcpy( gmtstring, "nogmt5" );

	strcpy( dbprog, "sqlite" ); /*** default ***/
        if(ioracle)
          strcpy( dbprog, "oracle" );
        else if(imysql)
          strcpy( dbprog, "mysql" );
        else if(isqlite)
          strcpy( dbprog, "sqlite" );

/***************************/
/*** load the glib files ***/
/***************************/

	grn = (Greens **)malloc( MAX_STA*sizeof(Greens *) );
	for( ista=0, i=1; i<ac; i++ )
	{
		if( (fp=fopen(av[i],"r")) == NULL )
		{
			continue;
		}
		fread( &nz, sizeof(int), 1, fp );
		if(verbose) fprintf( stderr, "%s: i=%d ista=%d file=%s nz=%d\n", 
				progname, i, ista, av[i], nz );
		z = malloc( nz*sizeof(float) );
		fread( &z[0], nz*sizeof(float), 1, fp );
		grn[ista] = (Greens *)malloc( nz*sizeof(Greens) );
		for( iz=0; iz<nz; iz++ )
		{
			fread( &grn[ista][iz], sizeof(Greens), 1, fp );
			nt = grn[ista][iz].nt;
			min = grn[ista][iz].g.rss[0];
			max = grn[ista][iz].g.rss[0];
			for( it=1; it<nt; it++ )
			{
			  if( grn[ista][iz].g.rss[it] < min ) min =  grn[ista][iz].g.rss[it];
			  if( grn[ista][iz].g.rss[it] > max ) max =  grn[ista][iz].g.rss[it];
			}
			if(verbose)
			{
			  fprintf(stderr, 
			    "%s: iz=%d z=%g %s.%s rdist=%.1f az=%.0f t0=%.2f dt=%.2f nt=%d min=%.2e max=%.2e fn=%s modf=%s\n",
				progname, iz, grn[ista][iz].evdp, grn[ista][iz].stnm,
				grn[ista][iz].net, grn[ista][iz].rdist, grn[ista][iz].az,
				grn[ista][iz].t0, grn[ista][iz].dt, grn[ista][iz].nt, min, max,
				grn[ista][iz].filename, grn[ista][iz].v.modfile );
			}
		}
		fclose(fp);
		ista++;
	}
	nsta= ista;
	if(verbose) fprintf(stderr, "nsta=%d\n", nsta );

/******************************************************************/
/*** sort the src-rec distance and print in increasing distance ***/
/******************************************************************/
	
	rdistkm = (float *)calloc( nsta+1, sizeof(float) );
	index = (int *)calloc( nsta+1, sizeof(int) );
	if(verbose) fprintf(stderr, "done allocating memory for rdistkm and index\n" );

	if( nsta == 1 )
	{
		rdistkm[1] = grn[0][0].rdist;
		index[1] = 1;
	}
	else
	{
		for( ista=0; ista<nsta; ista++ )
		{
			rdistkm[ista+1] = grn[ista][0].rdist;
		}
		indexx( nsta, rdistkm, index );
	}

/******************************************************************/
/*** begin writting run.csh file *****/
/******************************************************************/

	if( fopen( "run.csh", "r" ) != NULL )
	{
		fprintf( stderr, "%s: file run.csh already exists. ", progname );  
		fprintf( stderr, "I cannot overwrite. Please rename or delete the file run.csh\n\n" );
		exit(-1);
	}

	fprintf( stderr, "%s: writting output to file run.csh\n\n", progname );

	fp = fopen( "run.csh", "w" );

	fprintf( fp, "#!/bin/csh \n");
	fprintf( fp, "setenv MTINV_GMT_GRID_FILE /root/china.grd\n" );
	fprintf( fp, "setenv MTINV_GMT_INT_FILE  /root/china_i.int\n" );
	fprintf( fp, "setenv MTINV_GMT_CPT_FILE  /root/MYTOPO.cpt\n" );
      fprintf( fp, "\n" );
	if(irealtime) fprintf( fp, "#### Realtime version\n" );

	if( ifullmt )
	{
	   fprintf( fp, "\n### uncomment the one needed 1-isotropic mt  5-deviatoric mt  6-full mt\n" );
	   fprintf( fp, "set DEGFREE=6\n" );
	   fprintf( fp, "#set DEGFREE=1\n" );
	}
	else
	{
	   fprintf( fp, "set DEGFREE=5  # 1-isotropic_mt 5-deviatoric_mt 6-full_mt\n" );
	}
	
	fprintf( fp, "set ts0Arr=\"-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8\"\n");
	fprintf( fp, "\n" );

	fprintf( fp, "cat >! mtinv.par << EOF\n" );
	fprintf( fp, "#### REGION COMMENT ############################\n" );
	fprintf( fp, "CM %s\n", comment );
	fprintf( fp, "#### Date and Origin Time ######################\n" );
	fprintf( fp, "OT %s\n", dateid );
	fprintf( fp, "#### Forward Calculations ######################\n" );
	fprintf( fp, "## stk  dip rak   Mw  evlo  evla  Z   ##########\n" );
	fprintf( fp, "EV   0    0   0   0.0 %g    %g    15  \n", grn[0][0].evlo, grn[0][0].evla );
	fprintf( fp, "#####################################################################################################\n" );
	fprintf( fp, "# sta net model  np pas lf  hf  nt  dt   tr  tt v/d  mulfac used(Y/N)  ts0  weight ###              #\n" );

	
	for( ista=0; ista<nsta; ista++ )
	{
		ix = index[ista+1]-1;
		select_dt_and_nt( grn[ix][0].rdist, LOWEST_VELOCITY, grn[ix][0].nt, grn[ix][0].dt, &my_dt, &my_nt );

		if( my_dt < grn[ix][0].dt ) my_dt = grn[ix][0].dt;

		if( ista < 50 )
			use = 'y';
		else
			use = 'n';

		strcpy( com, "" );
		if( grn[ix][0].rdist > 1000 ) use = 'n';
/***		if( grn[ix][0].rdist > 1300 || ista > 7 ) strcpy( com, "#" );***/

		fprintf( fp, "%s%s\t%s\t%s 3 2 %.3f %.3f %5d %5.2f 0.0 0.0 d  1.0 %c  +0.00 +1.0 Surf/Pnl ### R=%.1f Az=%.0f\n",
			com,
			grn[ix][0].stnm,
			grn[ix][0].net,
			grn[ix][0].v.modfile,
			lf, 
			hf,
			my_nt, 
			my_dt,
			use, 
			grn[ix][0].rdist,
			grn[ix][0].az 
		);
	}
	fprintf( fp, "EOF\n" );
	fprintf( fp, "\n" );

	fprintf( fp, "### CLEAN UP ###\n" );
	fprintf( fp, "/bin/rm -f *.ginv *.data \n" );
	fprintf( fp, "/bin/rm -f plot_T???.?sec_Z???.?km_.p??.ps email_T???.?sec_Z???.?km_.txt *.?.dat.xy *.?.syn.xy\n" );
	fprintf( fp, "/bin/rm -f plot_T???.?sec_Z???.?km_.p??.jpg plot_T???.?sec_Z???.?km_.p??.pdf automt.out *.sql\n" );
	fprintf( fp, "/bin/rm -f results.?.??? plotmech.??? plotz.??? gmtmap.??? mtinv.out multithread_mkgrnlib.out snr.out var_red.out\n" );
	if( irealtime ) fprintf( fp, "/bin/rm -f automt.txt\n" );
	fprintf( fp, "\n" );

	fprintf( fp, "### PROCESS GREENS FUNCTIONS ###\n" );
	fprintf( fp, "glib2inv par=mtinv.par noverbose parallel\n" );
	fprintf( fp, "\n" );

	fprintf( fp, "### PROCESS DATA ### \n" );
	fprintf( fp, "sacdata2inv par=mtinv.par path=%s respdir=%s noverbose nodumpsac parallel\n", DataDir, RespDir );
	fprintf( fp, "\n" );

/*********************************
	fprintf( fp, "### FORWARD CALCULATION (Set Parameters in EV line of mtinv.par) ###\n" );
	fprintf( fp, "# mtinv ts0=0 %s par=mtinv.par mtdegfree=${DEGFREE} fwd  ### FixISOZ=1.0\n", gmtstring );
	fprintf( fp, "\n" );
*************************************/
	
	
	if( irealtime )
	{
		fprintf( fp, "foreach ts0 ( -8 -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 8 )\n" );
		fprintf(fp, " mtinv AutoAuth ts0=${ts0} par=mtinv.par %s mtdegfree=${DEGFREE} use_snr minsnr=%g shift ctol=%g maxshift=%g >> mtinv.out\n", 
				gmtstring, minsnr, ctol, maxshift );
		fprintf( fp, "end\n" );
		fprintf( fp, "\n" );
	}
	else
	{
		fprintf(fp, "multi_mtinv ${ts0Arr} par=mtinv.par gmt5 mtdegfree=${DEGFREE} use_snr minsnr=%g shift ctol=%g maxshift=%g >> mtinv.out\n",
				gmtstring, minsnr, ctol, maxshift);
		fprintf(fp, "\n");
	}
	

/***********************************/
/*** check origin time shift    ****/
/***********************************/
	fprintf( fp, "### CHECK ORIGIN TIME SHIFT ###\n" );
	fprintf( fp, "csh results.${DEGFREE}.csh\n" );

	if( !irealtime )
	{
	  if( strcmp( my_uname.sysname, "Darwin" ) == 0 )
	  {
		fprintf( fp, "# open results.${DEGFREE}.jpg\n" );
	  }
	  else if( strcmp( my_uname.sysname, "Linux" ) == 0 )
	  {
		fprintf( fp, "# eog results.${DEGFREE}.jpg\n" );
	  }
	  else if( strcmp( my_uname.sysname, "SunOS" ) == 0 )
	  {
		fprintf( fp, "# xv results.${DEGFREE}.jpg\n" );
	  }
	  else
	  {
		fprintf( fp, "# gs results.${DEGFREE}.ps\n" );
	  }
	}
	else
	{
	  fprintf( fp, "### uncomment for review\n" );
	  fprintf( fp, "# open results.?.jpg\n" );
	}
	fprintf( fp, "\n" );

/***********************************/
/*** view waveform fits          ***/
/***********************************/
/*
	fprintf( fp, "### Use Ghostview to view PS files ###\n" );
	fprintf( fp, "#gs -dEPSCrop plot_T???.?sec_Z???.?km_.p??.ps\n" );
*/
	fprintf( fp, "### convert each postscript file to jpg \n" );
	fprintf( fp, "foreach ps ( plot_T???.?sec_Z???.?km_.p??.ps ) \n" );
	if(igmt5)
	{
		fprintf( fp, "  gmt psconvert -Z -Tj -E300 ${ps}\n" ); /*** GMT 5.x.x ***/
	}
	else
	{
		fprintf( fp, "  gmt ps2raster -Tj -E120 -A ${ps}\n" ); /*** GMT 4.x.x ***/
	}
	fprintf( fp, "end\n" );

	if( !irealtime )
	{
	  if( strcmp( my_uname.sysname, "Darwin" ) == 0 )
          {
		fprintf( fp, "# open plot_T???.?sec_Z???.?km_.p??.jpg\n" );
          }
          else if( strcmp( my_uname.sysname, "Linux" ) == 0 )
          {
		fprintf( fp, "# eog plot_T???.?sec_Z???.?km_.p??.jpg\n" );
          }
          else if( strcmp( my_uname.sysname, "SunOS" ) == 0 )
          {
		fprintf( fp, "# xv plot_T???.?sec_Z???.?km_.p??.jpg\n" );
                fprintf( fp, "# acroread plot_T???.?sec_Z???.?km_.p??.pdf\n" );
		fprintf( fp, "# mergepdf.csh plot_T???.?sec_Z???.?km_.p??.pdf \n" );
		fprintf( fp, "# acroread merged.pdf\n" );
          }
          else
          {
                fprintf( fp, "# gs -dEPSCrop plot_T???.?sec_Z???.?km_.p??.ps\n" );
          }
	}
	else
	{
		fprintf( fp, "### uncomment \n" );
		fprintf( fp, "#open plot_T???.?sec_Z???.?km_.p??.jpg\n" );
	}
        fprintf( fp, "\n" );

/***********************************/
/*** make depth sensitivity plot ***/
/***********************************/

	fprintf( fp, "### MAKE DEPTH SENSITIVITY PLOT ###\n" );
	fprintf( fp, "csh plotz.csh\n" );

	if( strcmp( my_uname.sysname, "Darwin" ) == 0 )
        {
                fprintf( fp, "#open plotz.jpg\n" );
        }
        else if( strcmp( my_uname.sysname, "Linux" ) == 0 )
        {
                fprintf( fp, "#eog plotz.jpg\n" );
        }
        else if( strcmp( my_uname.sysname, "SunOS" ) == 0 )
        {
                fprintf( fp, "#xv plotz.jpg\n" );
        }
        else
        {
                fprintf( fp, "#gs plotz.ps\n" );
        }
        fprintf( fp, "\n" );

/**********************************************/
/*** MAKE DEPTH / OT-SHIFT SENSITIVITY PLOT ***/
/**********************************************/

	fprintf( fp, "### MAKE DEPTH / OT-SHIFT SENSITIVITY PLOT ###\n" );
	fprintf( fp, "csh plotmech.csh\n" );

	if( strcmp( my_uname.sysname, "Darwin" ) == 0 )
        {
                fprintf( fp, "#open plotmech.jpg\n" );
        }
        else if( strcmp( my_uname.sysname, "Linux" ) == 0 )
        {
                fprintf( fp, "#eog plotmech.jpg\n" );
        }
        else if( strcmp( my_uname.sysname, "SunOS" ) == 0 )
        {
                fprintf( fp, "#xv plotmech.jpg\n" );
        }
        else
        {
                fprintf( fp, "#gs plotmech.ps\n" );
        }
        fprintf( fp, "\n" );

/****************************************************/
/*** MAKE GMT PLOT WITH LOCATION/STATION/SOLUTION ***/
/*** this moves into the run2.csh script          ***/
/****************************************************/

	if(!irealtime)
	{
          fprintf( fp, "### MAKE GMT PLOT WITH LOCATION/STATION/SOLUTION ###\n" );
	  fprintf( fp, "\n" );
	  fprintf( fp, "mtinv \\\n" );
	  fprintf( fp, "\t ts0=0          \\\n" );
	  fprintf( fp, "\t par=mtinv.par  \\\n" );
	  fprintf( fp, "\t %s             \\\n", gmtstring );
	  fprintf( fp, "\t mtdegfree=${DEGFREE} \\\n" );
	  fprintf( fp, "\t gmtmap         \\\n" );
	  fprintf( fp, "\t nodumpxy       \\\n" );
	  fprintf( fp, "\t %s             \\\n", dbprog );
	  fprintf( fp, "\t PltXcorLabel   \\\n" );
	  fprintf( fp, "\t use_snr        \\\n" );
	  fprintf( fp, "\t minsnr=%g      \\\n", minsnr );
	  fprintf( fp, "\t ctol=%g        \\\n", ctol );
	  fprintf( fp, "\t maxshift=%g >> mtinv.out\n", maxshift );
	  fprintf( fp, "\n" );

          fprintf( fp, "csh gmtmap.csh\n" );

          if( strcmp( my_uname.sysname, "Darwin" ) == 0 )
          {
                fprintf( fp, "#open gmtmap.jpg\n" );
          }
          else if( strcmp( my_uname.sysname, "Linux" ) == 0 )
          {
                fprintf( fp, "#eog gmtmap.jpg\n" );
          }
          else if( strcmp( my_uname.sysname, "SunOS" ) == 0 )
          {
                fprintf( fp, "#xv gmtmap.jpg\n" );
          }
          else
          {
                fprintf( fp, "#gs gmtmap.ps\n" );
          }
	}

	fprintf( fp, "\n" );

	if( irealtime )
	{
		fprintf( fp, "mtbestfit\n" );
	}

	fclose(fp);

	chmod( "run.csh", S_IRWXU|S_IRGRP|S_IXGRP|S_IROTH|S_IXOTH );

/*************************************************************************/
/**** Write the CSHELL script to finalize the moment tensor solution ****/
/**** assumes the run.csh script was run and OT set to 0 offset      ****/
/*************************************************************************/

	if( fopen( "run2.csh", "r" ) != NULL )
        {
                fprintf( stderr, "%s: file run2.csh already exists. ", progname );
                fprintf( stderr, "I cannot overwrite. Please rename or delete the file run2.csh\n\n" );
                exit(-1);
        }
	else
	{
        	fprintf( stderr, "%s: writting output to file run2.csh\n\n", progname );
	}

        fp = fopen( "run2.csh", "w" );

        fprintf( fp, "#!/bin/csh \n");

	fprintf( fp, "setenv MTINV_GMT_GRID_FILE /root/china.grd\n" );
	fprintf( fp, "setenv MTINV_GMT_INT_FILE  /root/china_i.int\n" );
	fprintf( fp, "setenv MTINV_GMT_CPT_FILE  /root/MYTOPO.cpt\n" );
/***	fprintf( fp, "setenv MT_DATABASE_FILE  /Users/ichinose1/Work/mtinv.v3.0.5/data/mt.db\n" );***/
	fprintf( fp, "\n" );

        if( ifullmt )
        {
	   fprintf( fp, "###\n" );
           fprintf( fp, "### uncomment the one needed 1-isotropic mt  5-deviatoric mt  6-full mt\n" );
	   fprintf( fp, "###\n" );
           fprintf( fp, "set DEGFREE=6\n" );
           fprintf( fp, "#set DEGFREE=1\n" );
        }
        else
        {
	   fprintf( fp, "### 1-isotropic_mt  5-deviatoric_mt  6-full_mt\n" );
	   fprintf( fp, "###\n" );
           fprintf( fp, "set DEGFREE=5  # 1-isotropic_mt 5-deviatoric_mt 6-full_mt\n" );
        }
        fprintf( fp, "\n" );

	fprintf( fp, "###\n" );
	fprintf( fp, "### Clean Up\n" );
	fprintf( fp, "###\n" );

	fprintf( fp, "/bin/rm -f email_T???.?sec_Z???.?km_.txt plot_T???.?sec_Z???.?km_.p??.jpg *.ps\n" );
	fprintf( fp, "\n" );

	fprintf( fp, "###\n" );
	fprintf( fp, "### MAKE GMT PLOT WITH LOCATION/STATION/SOLUTION ###\n" );
	fprintf( fp, "###\n" );

        fprintf( fp, "mtinv \\\n" );
        fprintf( fp, "\t ts0=0          \\\n" );
        fprintf( fp, "\t par=mtinv.par  \\\n" );
        fprintf( fp, "\t %s             \\\n", gmtstring );
        fprintf( fp, "\t mtdegfree=${DEGFREE} \\\n" );
        fprintf( fp, "\t gmtmap         \\\n" );
        fprintf( fp, "\t nodumpxy       \\\n" );
        fprintf( fp, "\t %s             \\\n", dbprog );
        fprintf( fp, "\t PltXcorLabel   \\\n" );
        fprintf( fp, "\t use_snr        \\\n" );
        fprintf( fp, "\t minsnr=%g      \\\n", minsnr );
        fprintf( fp, "\t ctol=%g        \\\n", ctol );
        fprintf( fp, "\t maxshift=%g >> mtinv.out\n", maxshift );
        fprintf( fp, "\n" );

	if(igmt5)
        {
                fprintf( fp, "psconvert -Z -Tj -E300 plot_T???.?sec_Z???.?km_.p??.ps\n" );
        }
        else
        {
                fprintf( fp, "  ps2raster -Tj -E300 plot_T???.?sec_Z???.?km_.p??.ps\n" );
        }

        fprintf( fp, "csh gmtmap.csh\n" );

        if( strcmp( my_uname.sysname, "Darwin" ) == 0 )
        {
                fprintf( fp, "#open gmtmap.jpg\n" );
        }
        else if( strcmp( my_uname.sysname, "Linux" ) == 0 )
        {
                fprintf( fp, "#eog gmtmap.jpg\n" );
        }
        else if( strcmp( my_uname.sysname, "SunOS" ) == 0 )
        {
                fprintf( fp, "#xv gmtmap.jpg\n" );
        }
        else
        {
                fprintf( fp, "#gs gmtmap.ps\n" );
        }
	fprintf( fp, "\n" );

	fprintf( fp, "### dumpxy option for publication quality GMT plots\n" );
	fprintf( fp, "#gmtwf.csh\n" );
	fprintf( fp, "\n" );

	if( isqlite )
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
	fprintf( fp, "\n" );
	fclose(fp);

	chmod( "run2.csh", S_IRWXU|S_IRGRP|S_IXGRP|S_IROTH|S_IXOTH );

} /**** end of main() ***/

/**************************************************************************************************/
/*** Print Help Usage ***/
/**************************************************************************************************/

void Print_Usage()
{
	fprintf( stderr, "\n\t This program auto-generates the run.csh and run2.csh C-shell scripts that allow\n" );
	fprintf( stderr, "\t\t for processing and inversion for Regional longperiod surface wave data\n" );
	fprintf( stderr, "\t NOTE! This program will not overwrite old run.csh and run2.csh scripts.  User must delete manually for safety.\n\n" );

	fprintf( stderr, "Usage: \n" );

	fprintf( stderr,
	  "\t %s com=\"New Region Comment\" date=\"yyy/mm/dd,hh:mm:ss.ss\" lf=(float) hf=(float)\n", progname );
	fprintf( stderr,
	  "\t\t minsnr=(float) ctol=(float) maxshift=(float) [no]gmt5 [no]fullmt DataDir=(string) RespDir=(string)\n" );
	fprintf( stderr,
	  "\t\t [no]oracle [no]mysql [no]sqlite [no]help [no]verbose *.glib\n" );

	fprintf( stderr, "\n\t REQUIRED PARAMETERS:\n" );
	fprintf( stderr, "\t\t date=(string) format: yyyy/mm/dd,hh:mm:ss.ss - Origin Time in exact format\n" );
	fprintf( stderr, "\n\t OPTIONAL PARAMETERS:\n" );
	fprintf( stderr, "\t\t comment=(string) A short comment that goes in the CM tag [default \"New Region\"]\n" );
	fprintf( stderr, "\t\t DataDir=(string) Directory where the *.SAC files are [default ../Data]\n" );
	fprintf( stderr, "\t\t RespDir=(string) Directory where the SAC_PZs_* files are [default ../Resp]\n" );
	fprintf( stderr, "\t\t lf=(float) low frequency corner of the filter  [default 0.01 Hz]\n" );
	fprintf( stderr, "\t\t hf=(float) high frequency corner of the filter  [default 0.1 Hz]\n" );
	fprintf( stderr, "\t\t minsnr=(float) Minimum SNR allowed Peak-to-Peak amplitude in filtered band [default 3.0]\n" );
	fprintf( stderr, "\t\t ctol=(float)   Minimum Cross-Corrrelation needed to shift the waveforms [default 0.85]\n" );
	fprintf( stderr, "\t\t maxshift=(float) Maximum lag-time shift allowed when above ctol [default 10.0 sec]\n" ); 
	fprintf( stderr, "\t\t [no]gmt5       Make C-shell scripts compatible with GMT 5.x.x for plotting [default on]\n" );
	fprintf( stderr, "\t\t [no]fullmt     Do full moment tensor instead of deviatoric [default off]\n" );
	fprintf( stderr, "\t\t\t\t This can be easily changed in the run.csh and run2.csh scripts\n " );
	fprintf( stderr, "\t\t [no]oracle=(boolean) write SQL create & insert statements for ORACLE Database [default off]\n" );
	fprintf( stderr, "\t\t [no]mysql=(boolean) write SQL create & insert statements for MySQL Database [default off]\n" );
	fprintf( stderr, "\t\t [no]sqlite=(boolean) write SQL create & insert statements for Sqlite3 Database [default on]\n" );
	fprintf( stderr, "\t\t [no]help=(boolean)     prints this page [default off]\n" );
        fprintf( stderr, "\t\t [no]verbose=(boolean)  turn on verbosity [default off]\n" );
	fprintf( stderr, "\n\n" );

	exit(0);
}

void select_dt_and_nt( float dist, float vel, int grn_npts, float grn_dt, float *dat_dt, int *dat_nt )
{
	float dt, slow_time;
	int nt;
	float roundoff( float, float );
	slow_time = dist/vel;

	nt = 1024;
	dt = 0.45;

	if( dist < 2000 )
	{
		nt = 1024;
		dt = slow_time/nt;
		dt = roundoff( dt, 100 );
	}

	if( dist < 120 )
	{
		nt = 512;
		dt = slow_time/nt;
		dt = roundoff( dt, 100 );
	}

	if( dist <= 60 )
	{
		nt = 256;
		dt = grn_dt + 0.02;
		dt = roundoff( dt, 100 );
	}
	*dat_dt = dt;
	*dat_nt = nt;
}

float roundoff( float x, float base )
{
	return ((float)rint(x*base)/base);
}
