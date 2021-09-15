#include <stdio.h>
#include <stdlib.h>
#include <math.h>


int main(int ac, char **av)
{
	float a,b,x,y,phi;
	float x0,y0,x1,y1,theta;
	int i;
	float D2R;
	int hudson = 0;
	float t,k,xt,yk;
	int hudson_transform(float,float,float*,float*);
	int setpar(int,char **), mstpar(), getpar();
	void endpar();

	D2R = M_PI/180;

	setpar(ac,av);
	mstpar( "a", "f", &a);
	mstpar( "b", "f", &b);
	mstpar( "x0",  "f", &x0 );
	mstpar( "y0",  "f", &y0 );
	mstpar( "theta", "f", &theta );
	getpar( "hudson", "b", &hudson );
	endpar();

/*** this makes it right with GMT ***/
	theta = (90-theta);
	/* theta = theta + 90; */

	fprintf( stdout, ">\n" );
	for( phi = 0; phi <= 360; phi++ )
	{
	  x = a*cos( phi * D2R );
	  y = b*sin( phi * D2R );

	  /* fprintf( stdout, "%g %g\n",  x, y ); */
	
	  x1 = x * cos(theta * D2R) - y * sin(theta * D2R);
	  y1 = x * sin(theta * D2R) + y * cos(theta * D2R);

	  if( hudson )
	  {
	     t = x1 + x0;
	     k = y1 + y0;
	     hudson_transform(t,k,&xt,&yk);
	     fprintf( stdout, "%g %g %g %g\n", xt, yk, t, k );
	  }
	  else
	  {
	     fprintf( stdout, "%g %g\n",  x1+x0, y1+y0 );
	  }
	}
}

int hudson_transform( float t, float k, float *xc, float *yc )
{
        float x,y;
        float sign(float,float);
        if( fabs(t) > 1 || fabs(k) > 1 ) return 1;
        if( fabs(t) == 1 && fabs(k) == 1 ) t = 0;
        if( ( k >= 0 && t <= 0 ) || ( k <= 0 && t >= 0 ) )
        {
                x = t * ( 1 - fabs(k) );
                y = k;
        }
        else
        {
                x = sign(1,k) / ( 1/(fabs(t) * (1-fabs(k))) - 0.5 );
                y = k * ( 1 + fabs(x)/2 );
                if( (y/x) < 0.25 )
                {
                        if( k != 0 )
                                y = ( 1/(1/k - (2*sign(1,k))));
                        else
                                y = 0;
                        x = t * ( 1 + fabs(y));
                }
        }
        *xc = x;
        *yc = y;
        return 0;
}
                                                                                                                                                         
float sign(float a,float b)
{
        if( b>= 0 ) return fabs(a);
        else if( b < 0 ) return -fabs(a);
	else return fabs(a);
}
