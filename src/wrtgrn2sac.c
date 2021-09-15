#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/mt.h"         /** global datatype and structure declarations **/

char progname[128];

void wrtgrn2sac( Greens *g, int ista )
{
        Sac_Header sp;
        FILE *fp;
        char filename[256];
        char cmp[][3] = { "rss", "rds", "rdd", "rep", "zss", "zds", "zdd", "zep", "tss", "tds" };
        int i;

        for( i=0; i<10; i++ )
        {
                sprintf( filename, "%s.%03d.%-3.3s.%g.grns", 
			g->filename, ista, cmp[i], g->evdp );

                if( (fp = fopen(filename,"wb")) == NULL )
                {
                        printf("%s: cannot open output file %s\n", progname, filename );
                        exit(-1);
                }
		
                sp              = sac_null;
                sp.b            = g->t0;
                sp.delta        = g->dt;
                sp.e            = g->t0 + (g->nt * g->dt);
                sp.stla         = g->stla;
                sp.stlo         = g->stlo;
                sp.evla         = g->evla;
                sp.evlo         = g->evlo;
                sp.evdp         = g->evdp;

                sp.nvhdr = 6;
                sp.norid = 0;
                sp.nevid = 0;
                sp.iftype = ITIME; /* data type: IRLIM spec re im | IAMPH amp ph | IXY general x,y  */
                sp.idep   = IUNKN; /* not disp vel acc or volts */
                sp.iztype = IB;    /* types: IUNKN,IB,IDAY,IO,IA,ITn  */
                sp.ievtyp = IUNKN; /* type of event IQUAKE - earthquake */
                sp.npts = g->nt;

                sp.leven  = TRUE;  /* is data evenly sampled  */
                sp.lpspol = FALSE; /* LLL sets it */
                sp.lcalda = TRUE;  /* should az,baz,gcarc,dist be calculated?  */
                sp.lhdr5  = FALSE; /* LLL sets it */

		sprintf( sp.kstnm,	"%-4.4s",	g->stnm );
		/*** this causes MacOS-X BSD unix to abort ***/
                /* sprintf( sp.kevnm,	"%-16.16s",	filename ); */
                sprintf( sp.kcmpnm,	"%-3.3s",	cmp[i] );

                fwrite( &sp, sizeof(Sac_Header), 1, fp );
                if( i==0 ) fwrite( &(g->g.rss[0]), g->nt * sizeof(float), 1, fp );
                if( i==1 ) fwrite( &(g->g.rds[0]), g->nt * sizeof(float), 1, fp );
                if( i==2 ) fwrite( &(g->g.rdd[0]), g->nt * sizeof(float), 1, fp );
                if( i==3 ) fwrite( &(g->g.rep[0]), g->nt * sizeof(float), 1, fp );
                if( i==4 ) fwrite( &(g->g.zss[0]), g->nt * sizeof(float), 1, fp );
                if( i==5 ) fwrite( &(g->g.zds[0]), g->nt * sizeof(float), 1, fp );
                if( i==6 ) fwrite( &(g->g.zdd[0]), g->nt * sizeof(float), 1, fp );
                if( i==7 ) fwrite( &(g->g.zep[0]), g->nt * sizeof(float), 1, fp );
                if( i==8 ) fwrite( &(g->g.tss[0]), g->nt * sizeof(float), 1, fp );
                if( i==9 ) fwrite( &(g->g.tds[0]), g->nt * sizeof(float), 1, fp );
                fclose(fp);
        }
}
