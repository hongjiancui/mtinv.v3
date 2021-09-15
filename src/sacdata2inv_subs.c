#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <errno.h>

#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>

#include "../include/mt.h"

char progname[128];

/*** my own array of file names for myscandir() and getsacfiles()  ***/

typedef struct { char fn[256]; } FileNameList;

/*************************************************************************************/
/*** this subroutine has been retooled to use readdir() instead of scandir()       ***/
/*** differences between standard libc and BSD Unix is too great between platforms ***/
/*************************************************************************************/

int my_file_name_cmp( const void *a, const void *b )
{
	FileNameList *c1 = (FileNameList *)a;
	FileNameList *c2 = (FileNameList *)b;
	return strcmp( c1->fn, c2->fn );
}

void getsacfiles( char *stnm, char *net, char *pathname, char sacfiles[256][2048], int *nfiles, int verbose )
{
	int i, count=0, iword, ifile=0;
	FileNameList *filen;
	FileNameList *myscandir( char *, int *, FileNameList * );
	extern char progname[128];

	/** dont need, defined in stdlib.h **/
/***
     void qsort(void *base, size_t nel, size_t width, int (*compar)(const void *, const void *));
***/
	int my_file_name_cmp( const void *, const void * );

	char tmp[1024], *ptr, tok[2] = { '.', '\0' };
	void fatal_error( char *, char * );

	if( verbose )
	{
	  printf("%s: getsacfiles(): sta=%s net=%s calling directory %s\n", 
		progname, stnm, net, pathname );
	}

/*** allocate space for filename list ***/
	filen = (FileNameList *)malloc( sizeof(FileNameList) );

/*** get only sac file name list from remote directory ***/
	filen = (FileNameList *)myscandir( pathname, &count, filen );

/*** sort the file names ***/
	qsort( filen, count, sizeof(FileNameList), my_file_name_cmp );

	if( verbose )
	  fprintf( stdout, "%s: getsacfiles(): number of files = %d\n", progname, count );

	if( count <= 0 )
		fatal_error( "No files in this directory", pathname );

/*** find the sac files from the list that have the station name and network in them ***/
	for( i=0; i < count; i++ )
	{
		strcpy( tmp, filen[i].fn );
		iword = 1;
		ptr = strtok( tmp, tok );
		while( (ptr = strtok( NULL, tok )) != NULL )
		{
			iword++;
			if( strcmp( ptr, net ) == 0 )
			{
				ptr = strtok( NULL, tok );
				iword++;
				if( strcmp( ptr, stnm ) == 0 )
				{
					strcpy( sacfiles[ifile], filen[i].fn );
					if( verbose )
					{
					  fprintf( stdout, 
					    "%s: getsacfiles(): Found station=%s net=%s ifile=%d file=%s\n",
						progname, stnm, net, ifile, sacfiles[ifile] );
					}
					ifile++;
				}
			}
		}
	}
	*nfiles = ifile;
	free(filen);
}

FileNameList *myscandir( char *pathname, int *count, FileNameList *filen )
{
	struct stat f;
	DIR *dip;
	struct dirent *dit;
	int i=0;
	char *eptr;
	if( ( dip = opendir(pathname) ) == NULL )
	{
		perror("opendir");
		fprintf(stderr, "\n\nDirectory %s does not exist\n", pathname );
		exit(2);
	}

	while( (dit = readdir(dip)) != NULL )
	{
		if( (strcmp( dit->d_name, "." ) == 0) ||
			(strcmp(dit->d_name, ".." ) == 0 ) ) continue;
		stat( dit->d_name, &f );

		/** does not work on suns or linux **/
		/** gichinose Apr 2007 ***/
		/** if( S_ISDIR( f.st_mode ) ) continue; **/

		eptr = rindex( dit->d_name, '.' );
		if( eptr == NULL ) continue;
		if( strncmp(eptr, ".SAC", 4) == 0 || strncmp(eptr, ".sac", 4) == 0 )
		{
			strcpy( filen[i].fn, dit->d_name );
			i++;
			filen = (FileNameList *)realloc( filen, (i+1)*sizeof(FileNameList) );
		}
	}
	
	if( closedir(dip) == -1 )
	{
		perror("closedir");
		exit(3);
	}
	*count = i;
	return filen;
}

void fatal_error( char *msg1, char *msg2 )
{
	extern char progname[128];
        fprintf(stderr, "%s: FATAL ERROR msg=%s %s\n", progname, msg1, msg2 );
        fprintf(stdout, "%s: FATAL ERROR msg=%s %s\n", progname, msg1, msg2 );
        exit(-1);
}


float *getsacdata( 
	char *cmp,
	float *x,
	Sac_Header *sp,
	char *sacfile,
	char *filename,
	int *ifound,
	int verbose )
{
	extern char progname[128];
	int kk;
	FILE *fp;
	Sac_Header *a;

	float *readsac( Sac_Header *s, char *filename, int verbose );

	a = (Sac_Header *)malloc(sizeof(Sac_Header));

/*** readsac() - autodetects byteorder and loads sac file does swapbytes ***/

	readsac( a, sacfile, verbose );

/*** due to readsac, no longer needed ***/
/***
        if( (fp = fopen( sacfile, "rb" )) == NULL )
	{
		printf("error cannot open file %s\n", sacfile );
		exit(-1);
	}
        fread( a, sizeof(Sac_Header), 1, fp );
	rewind(fp);
****/

	*ifound = 0;

	if( strncmp( cmp, "VER", 3 ) == 0 )
	{
		if(
		  ( a->cmpinc ==   0 && a->cmpaz == 0 && a->kcmpnm[2] == 'Z' ) ||
		  ( a->cmpinc == 180 && a->cmpaz == 0 && a->kcmpnm[2] == 'Z' ))
		{
			strcpy( filename, sacfile );

		/*** readsac() - autodetects byteorder and loads sac file does swapbytes ***/
			x = readsac( sp, sacfile, verbose );

		/*** due to readsac, no longer needed ***/
		/***
			fread( sp, sizeof( Sac_Header), 1, fp );
			x = (float *)realloc( x, sp->npts*sizeof(float) );
			fread( &x[0], sp->npts * sizeof(float), 1, fp );
		***/
			for( kk=0; kk<8; kk++ )
			{
				if( sp->kstnm[kk]  == ' ' ) sp->kstnm[kk]='\0';
				if( sp->kcmpnm[kk] == ' ' ) sp->kcmpnm[kk]='\0';
				if( sp->knetwk[kk] == ' ' ) sp->knetwk[kk]='\0';
			}

			printf("%s: getsacdata(): VER sta=%s net=%s cmp=%s az=%7.2f inc=%7.2f %s\n",
				progname, 
				sp->kstnm, 
				sp->knetwk, 
				sp->kcmpnm,
				sp->cmpaz, 
				sp->cmpinc, 
				filename );
			*ifound = 1;
		}
	}
	else if( strncmp( cmp, "NS", 2 ) == 0 )
	{
		if(
		  ( a->cmpinc == 90 && a->kcmpnm[2] == 'N' ) ||
                  ( a->cmpinc == 90 && a->kcmpnm[2] == '1' ))
		{
			strcpy( filename, sacfile );

		/*** readsac() - autodetects byteorder and loads sac file does swapbytes ***/
                        x = readsac( sp, sacfile, verbose );
                                                                                                                                                         
                /*** due to readsac, no longer needed ***/
		/***
			fread( sp, sizeof(Sac_Header), 1, fp );
			x = (float *)realloc( x, sp->npts*sizeof(float) );
			fread( &x[0], sp->npts * sizeof(float), 1, fp );
		***/

			for( kk=0; kk<8; kk++ )
			{
				if( sp->kstnm[kk]  == ' ' ) sp->kstnm[kk]='\0';
				if( sp->kcmpnm[kk] == ' ' ) sp->kcmpnm[kk]='\0';
				if( sp->knetwk[kk] == ' ' ) sp->knetwk[kk]='\0';
			}

			printf("%s: getsacdata(): NS  sta=%s net=%s cmp=%s az=%7.2f inc=%7.2f %s\n",
                                progname,
                                sp->kstnm,
                                sp->knetwk,
                                sp->kcmpnm,
                                sp->cmpaz,
                                sp->cmpinc,
                                filename );
			*ifound = 1;
		}
	}
	else if( strncmp( cmp, "EW", 2 ) == 0 )
	{
		if(
                  ( a->cmpinc == 90 && a->kcmpnm[2] == 'E' ) ||
                  ( a->cmpinc == 90 && a->kcmpnm[2] == '2' ))
                {
			strcpy( filename, sacfile );

		/*** readsac() - autodetects byteorder and loads sac file does swapbytes ***/
			x = readsac( sp, sacfile, verbose );

		/*** due to readsac, no longer needed ***/
		/***
			fread( sp, sizeof(Sac_Header), 1, fp );
			x = (float *)realloc( x, sp->npts*sizeof(float) );
			fread( &x[0], sp->npts * sizeof(float), 1, fp );
		***/
			for( kk=0; kk<8; kk++ )
			{
				if( sp->kstnm[kk]  == ' ' ) sp->kstnm[kk]='\0';
				if( sp->kcmpnm[kk] == ' ' ) sp->kcmpnm[kk]='\0';
				if( sp->knetwk[kk] == ' ' ) sp->knetwk[kk]='\0';
			}

			printf("%s: getsacdata(): EW  sta=%s net=%s cmp=%s az=%7.2f inc=%7.2f %s\n",
                                progname,
                                sp->kstnm,
                                sp->knetwk,
                                sp->kcmpnm,
                                sp->cmpaz,
                                sp->cmpinc,
                                filename );
                        *ifound = 1;
		}
	}
	else
	{
		x = (float *)NULL;
		*ifound = 0;
	}
	free(a);
	/* fclose(fp); */
        return (float *)x;
}

void fix_component_names( EventInfo *ev )
{
      if( strncmp( ev->z.s.kcmpnm,  "LHZ", 3 ) == 0 ) sprintf( ev->z.s.kcmpnm,  "LHZ" );
      if( strncmp( ev->ns.s.kcmpnm, "LHN", 3 ) == 0 ) sprintf( ev->ns.s.kcmpnm, "LHN" );
      if( strncmp( ev->ew.s.kcmpnm, "LHE", 3 ) == 0 ) sprintf( ev->ew.s.kcmpnm, "LHE" );
      if( strncmp( ev->ns.s.kcmpnm, "LH1", 3 ) == 0 ) sprintf( ev->ns.s.kcmpnm, "LH1" );
      if( strncmp( ev->ew.s.kcmpnm, "LH2", 3 ) == 0 ) sprintf( ev->ew.s.kcmpnm, "LH2" );
                
      if( strncmp( ev->z.s.kcmpnm,  "BHZ", 3 ) == 0 ) sprintf( ev->z.s.kcmpnm,  "BHZ" );
      if( strncmp( ev->ns.s.kcmpnm, "BHN", 3 ) == 0 ) sprintf( ev->ns.s.kcmpnm, "BHN" );
      if( strncmp( ev->ew.s.kcmpnm, "BHE", 3 ) == 0 ) sprintf( ev->ew.s.kcmpnm, "BHE" );
      if( strncmp( ev->ns.s.kcmpnm, "BH1", 3 ) == 0 ) sprintf( ev->ns.s.kcmpnm, "BH1" );
      if( strncmp( ev->ew.s.kcmpnm, "BH2", 3 ) == 0 ) sprintf( ev->ew.s.kcmpnm, "BH2" );
               
      if( strncmp( ev->z.s.kcmpnm,  "BLZ", 3 ) == 0 ) sprintf( ev->z.s.kcmpnm,  "BLZ" );
      if( strncmp( ev->ns.s.kcmpnm, "BLN", 3 ) == 0 ) sprintf( ev->ns.s.kcmpnm, "BLN" );
      if( strncmp( ev->ew.s.kcmpnm, "BLE", 3 ) == 0 ) sprintf( ev->ew.s.kcmpnm, "BLE" );
      if( strncmp( ev->ns.s.kcmpnm, "BL1", 3 ) == 0 ) sprintf( ev->ns.s.kcmpnm, "BL1" );
      if( strncmp( ev->ew.s.kcmpnm, "BL2", 3 ) == 0 ) sprintf( ev->ew.s.kcmpnm, "BL2" );
              
      if( strncmp( ev->z.s.kcmpnm,  "SHZ", 3 ) == 0 ) sprintf( ev->z.s.kcmpnm,  "SHZ" );
      if( strncmp( ev->ns.s.kcmpnm, "SHN", 3 ) == 0 ) sprintf( ev->ns.s.kcmpnm, "SHN" );
      if( strncmp( ev->ew.s.kcmpnm, "SHE", 3 ) == 0 ) sprintf( ev->ew.s.kcmpnm, "SHE" );
      if( strncmp( ev->ns.s.kcmpnm, "SH1", 3 ) == 0 ) sprintf( ev->ns.s.kcmpnm, "SH1" );
      if( strncmp( ev->ew.s.kcmpnm, "SH2", 3 ) == 0 ) sprintf( ev->ew.s.kcmpnm, "SH2" );
             
      if( strncmp( ev->z.s.kcmpnm,  "HHZ", 3 ) == 0 ) sprintf( ev->z.s.kcmpnm,  "HHZ" );
      if( strncmp( ev->ns.s.kcmpnm, "HHN", 3 ) == 0 ) sprintf( ev->ns.s.kcmpnm, "HHN" );
      if( strncmp( ev->ew.s.kcmpnm, "HHE", 3 ) == 0 ) sprintf( ev->ew.s.kcmpnm, "HHE" );
      if( strncmp( ev->ns.s.kcmpnm, "HH1", 3 ) == 0 ) sprintf( ev->ns.s.kcmpnm, "HH1" );
      if( strncmp( ev->ew.s.kcmpnm, "HH2", 3 ) == 0 ) sprintf( ev->ew.s.kcmpnm, "HH2" );

      if( strncmp( ev->z.s.kcmpnm,  "HLZ", 3 ) == 0 ) sprintf( ev->z.s.kcmpnm,  "HLZ" );
      if( strncmp( ev->ns.s.kcmpnm, "HLN", 3 ) == 0 ) sprintf( ev->ns.s.kcmpnm, "HLN" );
      if( strncmp( ev->ew.s.kcmpnm, "HLE", 3 ) == 0 ) sprintf( ev->ew.s.kcmpnm, "HLE" );
      if( strncmp( ev->ns.s.kcmpnm, "HL1", 3 ) == 0 ) sprintf( ev->ns.s.kcmpnm, "HL1" );
      if( strncmp( ev->ew.s.kcmpnm, "HL2", 3 ) == 0 ) sprintf( ev->ew.s.kcmpnm, "HL2" );

      if( strncmp( ev->ns.s.khole, "-12345", 6 ) == 0 ) ev->ns.s.khole[0]='\0';
      if( strncmp( ev->ew.s.khole, "-12345", 6 ) == 0 ) ev->ew.s.khole[0]='\0';
      if( strncmp( ev->z.s.khole,  "-12345", 6 ) == 0 ) ev->z.s.khole[0]='\0';

      if( strncmp( ev->ns.s.khole, "  ", 2 ) == 0 ) ev->ns.s.khole[0]='\0';
      if( strncmp( ev->ew.s.khole, "  ", 2 ) == 0 ) ev->ew.s.khole[0]='\0';
      if( strncmp( ev->z.s.khole,  "  ", 2 ) == 0 ) ev->z.s.khole[0]= '\0';

      if( strncmp( ev->ns.s.khole, "00", 2 ) == 0 ) sprintf( ev->ns.s.khole, "00" );
      if( strncmp( ev->ew.s.khole, "00", 2 ) == 0 ) sprintf( ev->ew.s.khole, "00" );
      if( strncmp( ev->z.s.khole,  "00", 2 ) == 0 ) sprintf( ev->z.s.khole,  "00" );

      if( strncmp( ev->ns.s.khole, "10", 2 ) == 0 ) sprintf( ev->ns.s.khole, "10" );
      if( strncmp( ev->ew.s.khole, "10", 2 ) == 0 ) sprintf( ev->ew.s.khole, "10" );
      if( strncmp( ev->z.s.khole,  "10", 2 ) == 0 ) sprintf( ev->z.s.khole,  "10" );
}
