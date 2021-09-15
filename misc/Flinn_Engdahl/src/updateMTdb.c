#include <stdio.h>
#include <stdlib.h>
#include <strings.h>

char DataBaseFile[] = { "/Users/ichinose1/Work/mtinv.v3.0.5/data/mt.db" };
char TableName[] = { "MT_ORIGIN_STAGE" };
char FLINN_ENGDAHL_PROGRAM[] = { "/Users/ichinose1/Work/mtinv.v3.0.5/bin/FlinnEngdahl" };

int main( int ac, char **av )
{
	FILE *fp;
	char database[256];
	char table[256];
	char program[256];
	char commandline[1024];

	int imax = 1;  /*** true just do the last load ***/
	int cleanup = 1; /** delete tmp files ***/

	int setpar(int,char **),getpar(),mstpar();
	void endpar();

	strcpy( database, DataBaseFile );
	strcpy( table, TableName );
	strcpy( program, FLINN_ENGDAHL_PROGRAM );

	setpar(ac,av);
	getpar("db","s", database );
	getpar("tb","s", table );
	getpar("prog","s", program );
	getpar("max","b", &imax );
	getpar("clean","b", &cleanup );
	endpar();
	
	fp = fopen( "query.csh","w");
	fprintf( fp, "#!/bin/csh\n" );
	fprintf( fp, "sqlite3 %s << EOF\n", database );
	fprintf( fp, ".output dump.out\n" );
	fprintf( fp, ".mode column\n" );
	fprintf( fp, ".header off\n" );

	if(imax)
	  fprintf( fp, "select lon, lat, max(orid) from %s;\n", table );
	else
	  fprintf( fp, "select lon, lat, orid from %s;\n", table );

	fprintf( fp, ".quit\n" );
	fprintf( fp, "EOF\n" );
	fclose(fp);
	system( "/bin/csh ./query.csh;" );

	sprintf( commandline, "%s < dump.out > update.sql", program );
	system( commandline );
	/* system( "/Users/ichinose1/bin/FlinnEngdahl < dump.out > update.sql" ); */
	
	fp = fopen( "update.csh", "w" );
	fprintf( fp, "#!/bin/csh\n" );
	fprintf( fp, "sqlite3 %s << EOF\n", database );
	fprintf( fp, ".read update.sql\n" );
	fprintf( fp, ".quit\n" );
	fprintf( fp, "EOF\n" );
	fclose(fp);

	system( "/bin/csh ./update.csh;" );
	system( "/bin/rm -f dump.out update.csh update.sql query.csh" );
}
