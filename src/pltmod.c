#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../include/mt.h"

char progname[128];

int main( int ac, char **av )
{
	VelMod vm;
	int init=0;
	int gmt=0;

	void create_mod( VelMod * );
	void reinit_mod( VelMod * );
	void print_mod1( FILE *, VelMod * );
	void plot_gmt_mod( VelMod *, int );
	void print_mod0( VelMod * );
	int setpar(),mstpar(),getpar();
	void endpar();

	strcpy( progname, av[0] );

	if( ac == 1 )
	{
		fprintf( stderr, "\t Usage: %s modeldb= velmod= [no]init [no]gmt\n",
			progname );
		fprintf( stderr, "\t\t or \n" );
		fprintf( stderr, "\t\t %s par=mkgrnlib.par [no]init [no]gmt\n", progname );
		exit(-1);
	}

	setpar(ac,av);
	mstpar( "modeldb", "s", &(vm.modpath) );
	mstpar( "velmod", "s", &(vm.modfile) );
	getpar( "init", "b", &init );
	getpar( "gmt", "b", &gmt );
	endpar();

	create_mod( &vm );

	if( init )
	{
		FILE *fp;
		reinit_mod( &vm );
		// print_mod0( &vm );
		fp = fopen("new.mod", "w");
		print_mod1( fp, &vm );
		fclose(fp);
	}

	if( gmt )
	{
		plot_gmt_mod( &vm, 1 );
	}

	if( init == 0 && gmt == 0 ) print_mod0( &vm );

	exit(-1);
}
