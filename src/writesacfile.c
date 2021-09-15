#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
                                                                                                                                        
#include "../include/mt.h"
                                                                                                                                        
char progname[128];

void writesacfile( EventInfo *ev )
{
        extern char progname[128];
        char outsacfile[256];
        FILE *fpsac;
                                                                                                                                        
        sprintf( outsacfile, "%s.%s.%-3.3s.cor", ev->stnm, ev->net, ev->z.s.kcmpnm );
        printf( "%s: insacfile=%s outsacfile=%s sta=%-4.4s net=%-2.2s cmp=%-3.3s n=%d dt=%6.4f caz=%03.0f cinc=%03.0f r=%.1f az=%.1f\n",
                progname,
                ev->z.filename, outsacfile, ev->z.s.kstnm, ev->z.s.knetwk, ev->z.s.kcmpnm,
                ev->z.s.npts, ev->z.s.delta, ev->z.s.cmpaz, ev->z.s.cmpinc, ev->z.s.dist, ev->z.s.az );
        fpsac = fopen( outsacfile, "wb" );
        fwrite( &ev->z.s, sizeof(Sac_Header), 1, fpsac );
        fwrite( &ev->z.data[0], ev->z.s.npts*sizeof(float), 1, fpsac );
        fclose(fpsac);
                                                                                                                                        
        sprintf( outsacfile, "%s.%s.%-3.3s.cor", ev->stnm, ev->net, ev->ns.s.kcmpnm );
        printf( "%s: insacfile=%s outsacfile=%s sta=%-4.4s net=%-2.2s cmp=%-3.3s n=%d dt=%6.4f caz=%03.0f cinc=%03.0f r=%.1f az=%.1f\n",
                progname,
                ev->ns.filename, outsacfile, ev->ns.s.kstnm, ev->ns.s.knetwk, ev->ns.s.kcmpnm,
                ev->ns.s.npts, ev->ns.s.delta, ev->ns.s.cmpaz, ev->ns.s.cmpinc, ev->ns.s.dist, ev->ns.s.az );
        fpsac = fopen( outsacfile, "wb" );
        fwrite( &ev->ns.s, sizeof(Sac_Header), 1, fpsac );
        fwrite( &ev->ns.data[0], ev->ns.s.npts*sizeof(float), 1, fpsac );
        fclose(fpsac);
                                                                                                                                        
        sprintf( outsacfile, "%s.%s.%-3.3s.cor", ev->stnm, ev->net, ev->ew.s.kcmpnm );
        printf( "%s: insacfile=%s outsacfile=%s sta=%-4.4s net=%-2.2s cmp=%-3.3s n=%d dt=%6.4f caz=%03.0f cinc=%03.0f r=%.1f az=%.1f\n",
                progname,
                ev->ew.filename, outsacfile, ev->ew.s.kstnm, ev->ew.s.knetwk, ev->ew.s.kcmpnm,
                ev->ew.s.npts, ev->ew.s.delta, ev->ew.s.cmpaz, ev->ew.s.cmpinc, ev->ew.s.dist, ev->ew.s.az );
                                                                                                                                        
        fpsac = fopen( outsacfile, "wb" );
        fwrite( &ev->ew.s, sizeof(Sac_Header), 1, fpsac );
        fwrite( &ev->ew.data[0], ev->ew.s.npts*sizeof(float), 1, fpsac );
        fclose(fpsac);
}
